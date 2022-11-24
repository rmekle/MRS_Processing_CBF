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

function [out, out_ref] = io_loadspec_IMA_s(dirString)
%% Modify for reading a whole directory of IMA files
%Load Dicom Info using Chris Rogers' "SiemensCsaParse.m" function:

oldFolder = cd(dirString);
dirData = dir('*.IMA');
files = char({dirData.name});
Nfiles = size(files, 1);

h = waitbar (0, 'Loading DICOM Files...', 'Name', 'Loading DICOM Files');
[fids(:,1), info] = SiemensCsaReadFid(SiemensCsaParse(files(1, :)), 0, 'conj');
for i = 2:Nfiles
    progress = (i-1.0) / Nfiles;
    progressStr = sprintf('%d of %d files loaded', i-1, Nfiles);
    waitbar (progress, h, progressStr);
    [fids(:,i), ~] = SiemensCsaReadFid(SiemensCsaParse(files(i, :)), 0, 'conj');
end

close(h);
cd(oldFolder);

sz = size(fids);

%% Read the regular parameters of the scan
Naverages = info.csa.NumberOfAverages;
% Bo = info.csa.MagneticFieldStrength;
    % Bo can be read from the metadata but does not seem to be accurate enough.
    % In the raw data Bo is lower than 3, but here it is exactly 3, not fitting
    % to the transmitting frequency given. So the latter is used and Bo
    % manually calculated
txfrq = info.csa.ImagingFrequency * 1000000;
te = info.csa.EchoTime;
tr = info.csa.RepetitionTime;
dwelltime = info.csa.RealDwellTime / 1000000000;
spectralwidth = 1 / dwelltime;
sequence = info.ProtocolName;

%% Determine if sLaser was used
isSVSdkdseq = contains(sequence,'svs_slaser_dkd');

%% Determine the number of averages
if(Nfiles < Naverages)
    Naverages = Nfiles;
elseif(Nfiles == Naverages)
    noRefScans = 0;
else
    noRefScans = Nfiles - Naverages;
end

% Typically the coils are already combined
Ncoils = 1;

%% Determine dimensions
if ndims(fids)==4  %Default config when 4 dims are acquired
    dims.t=1;
    dims.coils=2;
    dims.averages=3;
    dims.subSpecs=4;
elseif ndims(fids)<4  %To many permutations...ask user for dims.
    if Naverages == 1 && Ncoils == 1
        if ndims(fids)>1
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

%% Create the reference structure
if ~isSVSdkdseq || noRefScans == 0
    % will be logical 0 if isempty is called
    out_ref = struct([]);
else
    out_ref = out;
    out_ref.fids		= fidsRefScans;
    out_ref.specs		= specsRefScans;
    out_ref.sz			= size(fidsRefScans);
    out_ref.averages	= noRefScans;
    out_ref.rawAverages = noRefScans;
    out_ref.subspecs	= 1;
    out_ref.rawSubspecs = 1;
end

%DONE
