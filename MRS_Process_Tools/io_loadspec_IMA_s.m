%io_loadspec_IMA_s.m
%
% USAGE:
% [out, out_ref] = io_loadspec_IMA_s(dirString);
% 
% DESCRIPTION:
% Loads a directory with siemens .IMA files into matlab structure format.
% Each file is one average. Reference Scans for 'svs_slaser_dkd' sequences
% are determined by checking the number of averages expected with files in
% the directory.
% 
% INPUTS:
% dirString      = Directory containing several Siemens .IMA files to load.
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.
% out_ref    = Reference dataset in FID-A structure format. Empty if no
%              reference data was found/expected.
%
% Jamie Near, McGill University 2014
% Modified by Ivo Opitz, Charite Universitätsmedizin Berlin, 2022 to add
% option to load a whole directory of .IMA files into one structure.
% Modified by Ralf Mekle, Charite Universitätsmedizin Berlin, 2022

function [out, out_ref] = io_loadspec_IMA_s(dirString)

% %% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'io_loadspec_IMA_s.m';
sMsg_newLines		= sprintf('\n\n');
% sMsg_newLine		= sprintf('\n');
% disp(sMsg_newLines);


%% Modify for reading a whole directory of IMA files

% Load all DICOM files without changing working directory
%oldFolder	= cd(dirString);
%dirData		= dir('*.IMA');
dirData		= dir([dirString, '*.IMA']);
files		= char({dirData.name});
Nfiles		= size(files, 1);

% Indicate progress of loading all DICOM files
h		= waitbar (0, 'Loading DICOM Files...', 'Name', 'Loading DICOM Files');

% Load Dicom Info using Chris Rogers' "SiemensCsaParse.m" function:
fullFileName	= fullfile(dirString, files(1, :));
info			= SiemensCsaParse(fullFileName);
%info			= SiemensCsaParse(files(1, :));

% Read in first Dicom file using Chris Rogers' "SiemensCsaReadFid.m" function:
%[fids(:,1), info] = SiemensCsaReadFid(SiemensCsaParse(files(1, :)), 0, 'conj');
%[fids(:,1), info] = SiemensCsaReadFid(SiemensCsaParse(fullFileName), 0, 'conj');
[fids_shot, info] = SiemensCsaReadFid(SiemensCsaParse(fullFileName), 0, 'conj');

% Determine most important scan parameters from the DICOM info object
Naverages		= info.csa.NumberOfAverages;
% Bo = info.csa.MagneticFieldStrength;
    % Bo can be read from the metadata but does not seem to be accurate enough.
    % In the raw data Bo is lower than 3, but here it is exactly 3, not fitting
    % to the transmitting frequency given. So the latter is used and Bo
    % manually calculated
txfrq			= info.csa.ImagingFrequency * 1000000;
te				= info.csa.EchoTime;
tr				= info.csa.RepetitionTime;
dwelltime		= info.csa.RealDwellTime / 1000000000;
spectralwidth	= 1 / dwelltime;

% Typically the coils are already combined for MRS DICOM data
Ncoils = 1;

% Determine # of subspectra from DICOM info of first shot/average
% NOTE: This is an 'educated guess', since it is not clear how exactly the # of
% subspectra would be encoded in a real sequence; e.g., for SPECIAL spectra, the 
% subsprectra are obtained from averaging (subtraction), i.e. subspectra are encoded as
% averages
NsubSpectra		= info.csa.SpectroscopyAcquisitionOutofplanePhaseSteps;
if NsubSpectra > 1
	disp(sMsg_newLines);
	warning('%s: \tNsubSpectra = %d > 1!', sFunctionName, NsubSpectra);
	disp(sMsg_newLines);
end

% Determine # of dimensions of fids data struct for a single shot (average)
% (help ndims - The number of dimensions in an array is always greater than or 
% equal to 2. => ndims(fids) always > 1!
noDims_fids_shot	= ndims(fids_shot);
switch noDims_fids_shot
	case 4
	case 3
	case 2
		fids(:,1)	= fids_shot;
		
	otherwise
		error('%s: Unknown # of dimensions of fids_shot = ndims(fids_shot) = %d!', sFunctionName, noDims_fids_shot);
end		% End of switch ndims(fids_tmp)

% Read in all remaining DICOM files (one file for each average) and indicate progress of
% file loading
for i = 2:Nfiles
    progress		= (i-1.0) / Nfiles;
    progressStr		= sprintf('%d of %d files loaded', i-1, Nfiles);
    waitbar (progress, h, progressStr);
    %[fids(:,i), ~]	= SiemensCsaReadFid(SiemensCsaParse(files(i, :)), 0, 'conj');
	fullFileName	= fullfile(dirString, files(i, :));
	[fids(:,i), ~]	= SiemensCsaReadFid(SiemensCsaParse(fullFileName), 0, 'conj');
end

close(h);
%cd(oldFolder);

% Determine size of array of FIDs
sz = size(fids);


%% Determine which sequence protocol was used
% Specific entries/fields of the DICOM header might exist or might not exist, if MRS DICOM
% data were anonymized or pseudonomized
if(isfield(info, 'ProtocolName'))
	sequence		= info.ProtocolName;
else if(isfield(info.csa, 'SequenceName'))
		sequence		= info.csa.SequenceName;
	else
		sequence		= '';
	end
end

% In particular, determine if Dinesh's sLASER sequence was used
isSVSdkdseq		= contains(sequence,'svs_slaser_dkd');


%% Check on number of averages and determined numer of possible (water) refrence scans
if(Nfiles < Naverages)
    Naverages = Nfiles;
elseif(Nfiles == Naverages)
    noRefScans = 0;
else
    noRefScans = Nfiles - Naverages;
end


%% Determine dimensions
if ndims(fids)==4  %Default config when 4 dims are acquired
    dims.t=1;
    dims.coils=2;
    dims.averages=3;
    dims.subSpecs=4;
elseif ndims(fids)<4  %Too many permutations...ask user for dims.
    if Naverages == 1 && Ncoils == 1
		% Bug fix
		% (help ndims - The number of dimensions in an array is always greater than or 
		% equal to 2. => ndims(fids) always > 1!
		%if ndims(fids)>1
		if ndims(fids)>2
			error('%s: ndims(fids) = %d with Naverages = %d and Ncoils = %d! Unknown additional data dimension!', ...
				sFunctionName, ndims(fids), Naverages, Ncoils);
		end
		% ndims(fids) = 2 here, since <4 and not <2
		% With Naverages = 1, ndims(fids) = 2 the same for both cases, with our without 
		% subspectra; ndims cannot distinguish these two cases
		if NsubSpectra > 1
            dims.t=1;
            dims.coils=0;
            dims.averages=0;
            dims.subSpecs=2;
        else
            dims.t=1;
            dims.coils=0;
            dims.averages=0;
            dims.subSpecs=0;
        end
    elseif Naverages>1 && Ncoils==1
        if ndims(fids)>2
            dims.t=1;
            dims.coils=0;
            dims.averages=2;
            dims.subSpecs=3;
        else
            dims.t=1;
            dims.coils=0;
            dims.averages=2;
            dims.subSpecs=0;
        end
    elseif Naverages==1 && Ncoils>1
        if ndims(fids)>2
            dims.t=1;
            dims.coils=2;
            dims.averages=0;
            dims.subSpecs=3;
        else
            dims.t=1;
            dims.coils=2;
            dims.averages=0;
            dims.subSpecs=0;
        end
    elseif Naverages>1 && Ncoils>1
        if ndims(fids)>3
            dims.t=1;
            dims.coils=2;
            dims.averages=3;
            dims.subSpecs=4;
        else
            dims.t=1;
            dims.coils=2;
            dims.averages=3;
            dims.subSpecs=0;
        end
    end
%     dims.t=1;
%     dims.coils=input('Enter the coils Dimension (0 for none):  ');
%     dims.averages=input('Enter the averages Dimension (0 for none):  ');
%     dims.subSpecs=input('Enter the subSpecs Dimension (0 for none);  ');
end

% Necessary for some processing operations
dims.extras = 0;


%% Now get relevant scan parameters:*****************************

%Calculate Bo
Bo = txfrq / 42577000;

%Get Date
%date=info(:,1).InstanceCreationDate;
date='';

averages = Naverages;
rawAverages = Naverages;

subspecs = 1;
rawSubspecs = 1;


%% Calculate t and ppm arrays using the calculated parameters:
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(dims.t))):spectralwidth/(sz(dims.t)):(spectralwidth/2)-(spectralwidth/(2*sz(dims.t)))];
ppm = -f / (Bo*42.577);
ppm= ppm + 4.65;

t=[0:dwelltime:(sz(dims.t)-1)*dwelltime];


%% Sort data into averages and reference scans
if(noRefScans > 0)
    indicesRefScans	= [1:(noRefScans/2) (noRefScans/2+Naverages+1):Nfiles];
    indicesAverages	= [(noRefScans/2+1):(noRefScans/2+Naverages)];
    
    fidsAverages = fids(:, indicesAverages);
    fidsRefScans = fids(:, indicesRefScans);
    
    specsAverages = fftshift(ifft(fidsAverages,[],dims.t),dims.t);
    specsRefScans = fftshift(ifft(fidsRefScans,[],dims.t),dims.t);
else
    fidsAverages = fids;
    specsAverages = fftshift(ifft(fidsAverages,[],dims.t),dims.t);
end


%% FILLING IN DATA STRUCTURE
out.fids=fidsAverages;
out.specs=specsAverages;
out.sz=size(fidsAverages);
out.ppm=ppm;  
out.t=t;    
out.spectralwidth=spectralwidth;
out.dwelltime=dwelltime;
out.txfrq=txfrq;
out.date=date;
out.dims=dims;
out.Bo=Bo;
out.averages=averages;
out.rawAverages=rawAverages;
out.subspecs=subspecs;
out.rawSubspecs=rawSubspecs;
out.seq=sequence;
out.te=te;
out.tr=tr;
out.pointsToLeftshift=0;


%FILLING IN THE FLAGS
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=0;
out.flags.addedrcvrs=1;
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
if out.dims.subSpecs==0
    out.flags.isFourSteps=0;
else
    out.flags.isFourSteps=(out.sz(out.dims.subSpecs)==4);
end


%% Create a data structure for (water) reference scans
if ~isSVSdkdseq || noRefScans == 0
    % will be logical 0 if isempty is called
    out_ref = struct([]);
else
    out_ref				= out;
    out_ref.fids		= fidsRefScans;
    out_ref.specs		= specsRefScans;
    out_ref.sz			= size(fidsRefScans);
    out_ref.averages	= noRefScans;
    out_ref.rawAverages = noRefScans;
    out_ref.subspecs	= 1;
    out_ref.rawSubspecs = 1;
end

%DONE
