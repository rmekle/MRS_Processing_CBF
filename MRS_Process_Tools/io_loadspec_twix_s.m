%io_loadspec_twix_s.m
%Jamie Near, McGill University 2014.
%Edits from Franck Lamberton, 2017, Ralf Mekle (RM), Charite, 2021.
%
% USAGE:
% [out, out_ref] = io_loadspec_twix_s(filename);
% 
% DESCRIPTION:
% Reads in siemens twix raw data (.dat file) using the mapVBVD.m and 
% twix_map_obj.m functions from Philipp Ehses (philipp.ehses@tuebingen.mpg.de).
% 
% op_loadspec_twix outputs the data in structure format, with fields corresponding to time
% scale, fids, frequency scale, spectra, and header fields containing
% information about the acquisition.  The resulting matlab structure can be
% operated on by the other functions in this MRS toolbox.
% 
% INPUTS:
% filename   = filename of Siemens twix data to load.
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.
% out_ref    = Data from MRS refrence scans in FID-A structure format.

function [out, out_ref] = io_loadspec_twix_s(filename)

%% CBF: Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'io_loadspec_twix_s.m';
%sMsg_newLines		= sprintf('\n\n');
%sMsg_newLine		= sprintf('\n');
%disp(sMsg_newLines);
fprintf('\n\n');
fprintf('%s: Loading Siemens twix data for filename = %s\n\n', sFunctionName, filename);    


%% Read in twix data
%read in the data using the new mapVBVD.  This code has been adapted to 
%handle both single RAID files and multi-RAID files.  The vast majority of
%Siemens twix data comes as a single RAID file, but I've encoundered a few 
%multi-RAID files, particularly when using VD13D.  The way to distinguish
%them here is that a for a single RAID file, mapVBVD will output a struct, 
%whereas for a multi-RAID file, mapVBVD will output a cell array of structs.
%This code assumes that the data of interest is in the last element of the 
%cell array (possibly a bad assumption under some circumstances):
twix_obj=mapVBVD(filename);
if isstruct(twix_obj)
    disp('single RAID file detected.');
    RaidLength=1;
elseif iscell(twix_obj)
    disp('multi RAID file detected.');
    RaidLength=length(twix_obj);
    %assume that the data of interest is in the last element of the cell.
    twix_obj=twix_obj{RaidLength};
end
dOut.data=twix_obj.image();
version=twix_obj.image.softwareVersion;
sqzSize=twix_obj.image.sqzSize; 
sqzDims=twix_obj.image.sqzDims;


%find out what sequence, the data were acquired with.  If this is a
%multi-raid file, then the header may contain multiple instances of
%'tSequenceFileName' for different scans (including a pre-scan).
%Therefore, if multi-raid file, we will need to do a bit of extra digging 
%to find the correct sequence name.  
sequence=twix_obj.hdr.Config.SequenceFileName;  

%Try to find out what sequence this is:
isSpecial=~isempty(strfind(sequence,'rm_special')) ||...  %Is this Ralf Mekle's SPECIAL sequence?
            ~isempty(strfind(sequence,'vq_special'));  %or the CIBM SPECIAL sequence?
isjnSpecial=~isempty(strfind(sequence,'jn_svs_special')) ||...  %or Jamie Near's SPECIAL sequence?
            ~isempty(strfind(sequence,'md_Adiab_Special')) ||... %or Masoumeh Dehghani's Adiabatic SPECIAL sequence?
            ~isempty(strfind(sequence,'md_Special')) ||... %or another version of Masoumeh Dehghani's SPECIAL sequence?
            ~isempty(strfind(sequence,'md_Inv_special')); %or Masoumeh Dehghani's Inversion Recovery SPECIAL sequence?
ishdSPECIAL=~isempty(strfind(sequence,'md_dvox_special')); %Is this Masoumeh Dehghani's hadamard-encoded dual-SPECIAL sequence?
isjnMP=~isempty(strfind(sequence,'jn_MEGA_GABA')); %Is this Jamie Near's MEGA-PRESS sequence?
isjnseq=~isempty(strfind(sequence,'jn_')) ||... %Is this another one of Jamie Near's sequences 
        ~isempty(strfind(sequence,'md_'));      %or a sequence derived from Jamie Near's sequences (by Masoumeh Dehghani)?
isWIP529=~isempty(strfind(sequence,'edit_529')); %Is this WIP 529 (MEGA-PRESS)?
isWIP859=~isempty(strfind(sequence,'edit_859')); %Is this WIP 859 (MEGA-PRESS)?
isMinn=~isempty(strfind(sequence,'eja_svs_')); %Is this one of Eddie Auerbach's (CMRR, U Minnesota) sequences?
isSiemens=(~isempty(strfind(sequence,'svs_se')) ||... %Is this the Siemens PRESS seqeunce?
            ~isempty(strfind(sequence,'svs_st'))) && ... % or the Siemens STEAM sequence?
            isempty(strfind(sequence,'eja_svs'));    %And make sure it's not 'eja_svs_steam'.

		
%% CBF: Adjust routine for data acquired using 'svs_slaser_dkd' and 'svs_eja_slaser' sequences from CMRR
% CBF: Is this Dinesh K. Deelchand's single voxel (sLaser) sequence from CMRR, U Minnesota
% and init # of reference scans
isSVSdkdseq = contains(sequence,'svs_slaser_dkd');
noRefScans	= 0;
% For svs_eja sequences, init variable to indicate extraction of real/relevant data
% points, which includes any leftshift
bLeftshifted	= 0;

%If this is the SPECIAL sequence, it probably contains both inversion-on
%and inversion-off subspectra on a single dimension, unless it is the VB
%version of Jamie Near's SPECIAL sequence, in which case the subspecs are
%already stored on separate dimensions.  
%Both Ralf Mekle's SPECIAL and the VD-VE version of Jamie Near's SPECIAL sequence 
%do not store the subspectra along a separate dimension of the data array, 
%so we will separate them artifically:
%25 Oct 2018: Due to a recent change, the VE version of Jamie Near's MEGA-PRESS 
%sequence also falls into this category. 
if isSpecial ||... %Catches Ralf Mekle's and CIBM version of the SPECIAL sequence
		(strcmp(version,'vd') && isjnSpecial) ||... %and the VD/VE versions of Jamie Near's SPECIAL sequence
		(strcmp(version,'vd') && isjnMP);  %and the VD/VE versions of Jamie Near's MEGA-PRESS sequence
	squeezedData=squeeze(dOut.data);
	if twix_obj.image.NCol>1 && twix_obj.image.NCha>1
		data(:,:,:,1)=squeezedData(:,:,[1:2:end-1]);
		data(:,:,:,2)=squeezedData(:,:,[2:2:end]);
		sqzSize=[sqzSize(1) sqzSize(2) sqzSize(3)/2 2];
	elseif twix_obj.NCol>1 && twixObj.image.NCha==1
		data(:,:,1)=squeezedData(:,[1:2:end-1]);
		data(:,:,2)=squeezedData(:,[2:2:end]);
		sqzSize=[sqzSize(1) sqzSize(2)/2 2];
	end
	if isjnseq
		sqzDims{end+1}='Set';
	else
		sqzDims{end+1}='Ida';
	end
elseif ishdSPECIAL %For Masoumeh Dehghani's hadamard-encoded dual-voxel SPECIAL sequence:
	squeezedData=squeeze(dOut.data);
	if twix_obj.image.NCol>1 && twix_obj.image.NCha>1
		data(:,:,:,1)=squeezedData(:,:,[1:4:end-3]);
		data(:,:,:,2)=squeezedData(:,:,[2:4:end-2]);
		data(:,:,:,3)=squeezedData(:,:,[3:4:end-1]);
		data(:,:,:,4)=squeezedData(:,:,[4:4:end]);
		sqzSize=[sqzSize(1) sqzSize(2) sqzSize(3)/4 4];
	elseif twix_obj.image.NCol>1 && twix_obj.image.NCha==1
		data(:,:,1)=squeezedData(:,[1:4:end-3]);
		data(:,:,2)=squeezedData(:,[2:4:end-2]);
		data(:,:,3)=squeezedData(:,[3:4:end-1]);
		data(:,:,4)=squeezedData(:,[4:4:end]);
		sqzSize=[sqzSize(1) sqzSize(2)/4 4];
	end
	if isjnseq
		sqzDims{end+1}='Set';
	else
		sqzDims{end+1}='Ida';
	end
	
	% CBF: If this is a dkd MRS sequence, check whether reference scans are present
	% in the data and extract these into a separate data object
	% else
	% 	data=dOut.data;
	% end
else if isSVSdkdseq
		% All shots, i.e. reference scans and averages are stored in 'Set' dimension of
		% twix object
		% Find index of set dimension within cell array of dimensions and use this index
		% into array of sizes of dimensions assuming that only sizes of non-singleton
		% dimensions are stored in array of sizes of dimensions
		indSet	=  contains(sqzDims,'Set');
		% Obtain true # of averages and compute # of reference scans
		noAverages	= twix_obj.hdr.Meas.Averages;
		noRefScans	= sqzSize(indSet) - twix_obj.hdr.Meas.Averages;
		
		% Squeeze data in twix object for easier processing
		squeezedData	= squeeze(dOut.data);
		%sqz_ndims		= ndims(squeezeData);
		if noRefScans > 0
			% Data includes reference scans and averages
			% Extract reference scans and averages into separate data objects and update
			% information about data objects
			% In this case here, half of the reference scans are stored at beginning and
			% the other half at the end of the 'Set' dimension
			indicesRefScans	= [1:(noRefScans/2) (noRefScans/2+noAverages+1):sqzSize(indSet)];
			indicesAverages	= [(noRefScans/2+1):(noRefScans/2+noAverages)];
			
			% Use substruct indexing to extract selected data independent of # of
			% dimensions of squeezed data object
			% NOTE: For this to work, 'indSet' must also refer to 'Set' dimension in 
			% squeezed data object (array)
			% Define indexing structure for squeezed data object (array)
			% S.type is character vector or string scalar containing (), {}, or .,
			% specifying the subscript type; here it is '()' that is used for indexing
			% S.subs is cell array, character vector, or string scalar containing the
			% actual subscripts; here it is cell array of {':'} for each data dimension
			S.type			= '()';
			S.subs			= repmat({':'}, 1, ndims(squeezedData));
			
			% Select subscripts (indices) in 'Set' dimension of squezzed data object
			% for reference data
			S.subs{indSet}	= indicesRefScans;
			refData			= subsref(squeezedData, S);
			
			% Select subscripts (indices) in 'Set' dimension of squezzed data object
			% for averages
			S.subs{indSet}	= indicesAverages;
			data			= subsref(squeezedData, S);
			
			% 			% Check whether "Set' dimension is last dimension of data array and extract
			% 			% reference scans and averages according to # of non-singleton data dimensions
			% 			if indSet == sqz_ndims
			% 				% Use substruct indexing to extract selected data
			% 				refData		= squeezedData(:, :, indicesRefScans);
			% 				data		= squeezedData(:, :, indicesAverages);
			% 			else
			% 				error('Set dimension in data array not last data dimension! indSet = %d \t sqz_ndims = %d', indSet, sqz_ndims);
			% 			end		% End of if indSet == sqz_ndims
			
			% Update information about data objects
			% (correct, if indSet has NOT changed in squeezed data object)
			sqzSize(indSet)		= sqzSize(indSet) - noRefScans;
		else
			% No reference scans, data only includes averages
			data	= squeezedData;
		end		% End of if noRefScans > 0
		% CBF: If this is a eja MRS sequence, i.e. from Minnesota (isMinn) determine
		% relevant # of columns in the data and only pass on these columns as data object
		% (In eja_svs sequences some additional samples before the echo (peak of FID) and
		% at end of FID are acquired to avoid digital filtering effects; these samples can
		% be discarded, their actual number depends on sequence parameter settings and TE)
	else if isMinn
			% Find index of col dimension within cell array of dimensions and use this
			% index into array of sizes of dimensions assuming that only sizes of
			% non-singleton dimensions are stored in array of sizes of dimensions
			indCol	=  contains(sqzDims,'Col');
			
			% Squeeze data in twix object for easier processing
			squeezedData	= squeeze(dOut.data);
			%sqz_ndims		= ndims(squeezeData);
			
			% According to Eddie Auerbach, the author of eja_svs sequences,
			% "MDH header parameters, but the one you need to look for is 
			% KSpaceCenterColumn or ushKSpaceCentreColumn. This will tell you the position
			% of the first true measured point of the FID. ...
			% There will also be a few points left over at the end, but that number is not
			% explicitly recorded anywhere. It is simply what is left over once you
			% discard the initial dummy points and keep the number of expected measured
			% points (there will typically be 8 left over). Those are meant to pad the end
			% of the FID and absorb glitches that can appear when removing the 
			% oversampling (actually it doesn’t really work)"
			%
			% In the twix_obj returned by mapVBVD, parameters kspaceCenterColumn are all 
			% empty ([]), but the parameter array twix_obj.image.centerCol seems to store 
			% the first real/relevant point of each FID/average (same value for all 
			% averages)
			%
			%indFID_first		= twix_obj.image.centerCol(1);
			%indFID_first		= (twix_obj.image.centerCol(1)*twix_obj.hdr.Meas.ReadoutOSFactor) - 1;
			%indFID_first		= (twix_obj.image.centerCol(1)*twix_obj.hdr.Meas.ReadoutOSFactor);
			%indFID_first		= 1;
			%
			% Actually, using twix_obj.image.centerCol(...) did not seem to result in
			% proper reading of these data; following the Python-based suspect package 
			% from Benny Rowland, the # of points to left shift (in iceParam) is used as 
			% offset for reading the FID, and an additional # of dummy points that might
			% have been acquired before the first point of the FID is obtained from the 
			% parameter .freeParam from struct twix_obj.image
			% (the latter is inferred from the suspect code, but not yet fully verified!)
			% +1 is added for the index, since Matlab indexing starts at 1
			FID_offset			= twix_obj.image.iceParam(5,1);
			num_dummy_points	= twix_obj.image.freeParam(1,1);
			indFID_first		= FID_offset + num_dummy_points + 1;
			
			%% If # of dummy points is non-zero, indicate this
			%if num_dummy_points ~= 0
				% Display info
				disp(sMsg_newLine);
				fprintf('%s: Sequence isMinn = %d \t FID_offset = %d \t num_dummy_points = %d \t indFID_first = %d', sFunctionName, isMinn, FID_offset, num_dummy_points, indFID_first);
				%disp(sMsg_newLines);
				fprintf('\n\n');
			%end
			
			% Calculate index of last real/relevant point of FID using the index of the 
			% first real/relavant point of FID and the expected vector size including 
			% oversampling (subtract 1 to obtain correct vector size (length) of FID)
			if isfield(twix_obj.hdr.Meas, 'ReadoutOSFactor')
				indFID_last			= indFID_first + ...
					(twix_obj.hdr.Meas.lVectorSize*twix_obj.hdr.Meas.ReadoutOSFactor) - 1;
			else
				% No oversampling (sometimes OS factor then not a field in header)
				indFID_last			= indFID_first + twix_obj.hdr.Meas.lVectorSize - 1;
			end		% End of if exist(twix_obj.hdr.Meas.ReadoutOSFactor)
			
			% Extract real/relevant samples of all FIDs into data object
			% Determine indices of samples to be extracted
			indicesFID	= [indFID_first:indFID_last];
			
			% Use substruct indexing to extract selected data independent of # of
			% dimensions of squeezed data object
			% NOTE: For this to work, 'indCol' must also refer to 'Col' dimension in 
			% squeezed data object (array)
			% Define indexing structure for squeezed data object (array)
			% S.type is character vector or string scalar containing (), {}, or .,
			% specifying the subscript type; here it is '()' that is used for indexing
			% S.subs is cell array, character vector, or string scalar containing the
			% actual subscripts; here it is cell array of {':'} for each data dimension
			S.type			= '()';
			S.subs			= repmat({':'}, 1, ndims(squeezedData));
			
			% Select subscripts (indices) in 'Col' dimension of squezzed data object
			% for real/relevant points of all FIDs
			S.subs{indCol}	= indicesFID;
			data			= subsref(squeezedData, S);
			
			% Update information about data objects
			% (correct, if indCol has NOT changed in squeezed data object)
			% (here it could also be set to (Vector Size * OversamplingFactor))
			new_data_size		= size(data);
			sqzSize(indCol)		= new_data_size(indCol);
			
			% Indicate that data (FIDs) have already been left shifted
			bLeftshifted	= 1;
		else
			data=dOut.data;
		end		% End of else if isMinn
	end		% End of else if isSVSdkdseq
end		% End of if isSpecial ||... %Catches Ralf Mekle's and CIBM version of the SPECIAL sequence


%Squeeze the data to remove singleton dims
fids=squeeze(data);

%noticed that in the Siemens PRESS and STEAM sequences, there is sometimes
%an extra dimension containing unwanted reference scans or something.  Remove them here.
if isSiemens && (strcmp(version,'vd') || strcmp(version,'ve')) && strcmp(sqzDims{end},'Phs')
    sqzDims=sqzDims(1:end-1);
    sqzSize=sqzSize(1:end-1);
    if ndims(fids)==4
        fids=fids(:,:,:,2);
        fids=squeeze(fids);
    elseif ndims(fids)==3
        fids=fids(:,:,2);
        fids=squeeze(fids);
    elseif ndims(fids)==2
        fids=fids(:,2);
        fids=squeeze(fids);
    end
end

%Make a pulse sequence identifier for the header (out.seq);
seq=sequence;

%Find the magnetic field strength:
Bo=twix_obj.hdr.Dicom.flMagneticFieldStrength;

%Find the number of averages:
Naverages=twix_obj.hdr.Meas.Averages;

%Find out if multiple coil elements were used:
Ncoils=twix_obj.hdr.Meas.iMaxNoOfRxChannels;  

% CBF: For sequence svs_slaser_dkd, values for partial TEs (TE1, TE2, etc.) are stored in 
% array of TE values of header; TE = TE1 + TE2 + TE3, i.e. sum of partial TEs
if isSVSdkdseq
	TE	= sum([twix_obj.hdr.MeasYaps.alTE{:}]);
else
	%Find the TE:
	TE = twix_obj.hdr.MeasYaps.alTE{1};  %Franck Lamberton
end

%Find the TR:
TR = twix_obj.hdr.MeasYaps.alTR{1};  %Franck Lamberton

%Now begin indexing the dimensions of the data array. ie. create the dims
%structure, which specifies which dimensions of the data array are being
%used to hold the time-domain data, the multiple coil channels, the
%average, the sub-spectra, and any additional dimensions.
sqzDims_update=sqzDims;
dimsToIndex=[1:length(sqzDims)];

%First index the dimension of the time-domain data
dims.t=find(strcmp(sqzDims,'Col'));
if ~isempty(dims.t)
    %remove the time dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.t);
else
    dims.t=0;
    error('%s:  ERROR: Spectrum contains no time domain information!!', sFunctionName);
end

%Now index the dimension of the coil channels
dims.coils=find(strcmp(sqzDims,'Cha'));
if ~isempty(dims.coils)
    %remove the coils dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.coils);
else
    dims.coils=0;
end

% CBF: Extract # of averages from 'Set' dimension also for 'svs_sLaser_dkd' sequence
%Now index the dimension of the averages
if strcmp(version,'vd') || strcmp(version,'ve')
    %if isMinn
	if isMinn || isSVSdkdseq
        dims.averages=find(strcmp(sqzDims,'Set'));
    else
        dims.averages=find(strcmp(sqzDims,'Ave'));
    end
else
    dims.averages=find(strcmp(sqzDims,'Set'));
end
if ~isempty(dims.averages)
    %remove the averages dimension from the dimsToIndex vector
    dimsToIndex=dimsToIndex(dimsToIndex~=dims.averages);
else
    %If no Averages dimension was found, then check for a "Repetitions"
    %dimension.  If that is found, store it under "averages".  If both
    %"Averages" and "Repetitions" dimensions are found, "Repetitions" will
    %be indexed under "Extras", since "Repetitions is not currently an
    %option in FID-A.
    dims.averages=find(strcmp(sqzDims,'Rep'));
    if ~isempty(dims.averages)
        dimsToIndex=dimsToIndex(dimsToIndex~=dims.averages);
    else
        %If neither an "Averages" or a "Repetitions" dimension is found,
        %then set the FID-A "Averages" dimension to zero.
        dims.averages=0;
    end
end

%Now we have indexed the dimensions containing the timepoints, the coil
%channels, and the averages.  As we indexed each dimension, we removed the
%corresponding index from the dimsToIndex vector.  At this point, if there
%are any values left in the dimsToIndex vector, then there must be some
%additional dimensions that need indexing.  We assume that if sub-spectra exist,
%then these must be indexed in either the 'Ida' dimension (for all Jamie
%Near's VB-version pulse sequences), the 'Set' dimension (for all Jamie 
%Near's VD/VE-version pulse sequences), the 'Eco' dimension (for the WIP
%529 MEGA-PRESS sequence or the Minnesota MEGA-PRESS sequence), or the 'Ide' 
% dimension (for the WIP 859 MEGA-PRESS sequence). 
if ~isempty(dimsToIndex)
    %Now index the dimension of the sub-spectra
    if isjnseq  || isSpecial
        if strcmp(version,'vd') || strcmp(version,'ve')
            dims.subSpecs=find(strcmp(sqzDims,'Set'));
        else
            dims.subSpecs=find(strcmp(sqzDims,'Ida'));
        end
    elseif isWIP529 || isMinn
        dims.subSpecs=find(strcmp(sqzDims,'Eco'));
    elseif isWIP859
        dims.subSpecs=find(strcmp(sqzDims,'Ide'));
    else
        dims.subSpecs=dimsToIndex(1);
    end
    if ~isempty(dims.subSpecs)
        %remove the sub-spectra dimension from the dimsToIndex vector
        dimsToIndex=dimsToIndex(dimsToIndex~=dims.subSpecs);
    else
        dims.subSpecs=0;
    end
else
    dims.subSpecs=0;
end

%And if any further dimensions exist after indexing the sub-spectra, call
%these the 'extras' dimension.  
if ~isempty(dimsToIndex)
    %Now index the 'extras' dimension
    dims.extras=dimsToIndex(1);
    if ~isempty(dims.extras)
        %remove the extras dimension from the dimsToIndex vector
        dimsToIndex=dimsToIndex(dimsToIndex~=dims.extras);
    else
        dims.extras=0;
    end
else
    dims.extras=0;
end

%Now that we've indexed the dimensions of the data array, we now need to
%permute it so that the order of the dimensions is standardized:  we want
%the order to be as follows:  
%   1) time domain data.  
%   2) coils.
%   3) averages.
%   4) subSpecs.
%   5) extras.
% CBF: Init and store dimorder, with which data dimensions are permuted for reuse
% Init dimorder with zeros to detect error, i.e. in case dimorder is not correctly set
dimorder	= zeros(1, length(sqzDims));		% [1:length(sqzDims)];
if length(sqzDims)==5
	dimorder = [dims.t dims.coils dims.averages dims.subSpecs dims.extras];
    fids=permute(fids,[dims.t dims.coils dims.averages dims.subSpecs dims.extras]);
    dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=4;dims.extras=5;	
elseif length(sqzDims)==4
    if dims.extras==0
		dimorder = [dims.t dims.coils dims.averages dims.subSpecs];
        fids=permute(fids,[dims.t dims.coils dims.averages dims.subSpecs]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=4;dims.extras=0;
    elseif dims.subSpecs==0
		dimorder = [dims.t dims.coils dims.averages dims.extras];
        fids=permute(fids,[dims.t dims.coils dims.averages dims.extras]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=0;dims.extras=4;
    elseif dims.averages==0
		dimorder = [dims.t dims.coils dims.subSpecs dims.extras];
        fids=permute(fids,[dims.t dims.coils dims.subSpecs dims.extras]);
        dims.t=1;dims.coils=2;dims;averages=0;dims.subSpecs=3;dims.extras=4;twix_obj.image.centerCol
    elseif dims.coils==0
		dimorder = [dims.t dims.averages dims.subSpecs dims.extras];
        fids=permute(fids,[dims.t dims.averages dims.subSpecs dims.extras]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=3;dims.extras=4;
    end
elseif length(sqzDims)==3
    if dims.extras==0 && dims.subSpecs==0
		dimorder = [dims.t dims.coils dims.averages];
        fids=permute(fids,[dims.t dims.coils dims.averages]);
        dims.t=1;dims.coils=2;dims.averages=3;dims.subSpecs=0;dims.extras=0;
    elseif dims.extras==0 && dims.averages==0
		dimorder = [dims.t dims.coils dims.subSpecs];
        fids=permute(fids,[dims.t dims.coils dims.subSpecs]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=3;dims.extras=0;
    elseif dims.extras==0 && dims.coils==0
		dimorder = [dims.t dims.averages dims.subSpecs];
        fids=permute(fids,[dims.t dims.averages dims.subSpecs]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=3;dims.extras=0;
    end
elseif length(sqzDims)==2
    if dims.extras==0 && dims.subSpecs==0 && dims.averages==0
		dimorder = [dims.t dims.coils];
        fids=permute(fids,[dims.t dims.coils]);
        dims.t=1;dims.coils=2;dims.averages=0;dims.subSpecs=0;dims.extras=0;
    elseif dims.extras==0 && dims.subSpecs==0 && dims.coils==0
		dimorder = [dims.t dims.averages];
        fids=permute(fids,[dims.t dims.averages]);
        dims.t=1;dims.coils=0;dims.averages=2;dims.subSpecs=0;dims.extras=0;
    elseif dims.extras==0 && dims.averages==0 && dims.coils==0
		dimorder = [dims.t dims.subSpecs];
        fids=permute(fids,[dims.t dims.subSpecs]);
        dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=2;dims.extras=0;
    end
elseif length(sqzDims)==1
	dimorder = [dims.t];
    fids=permute(fids,[dims.t]);
    dims.t=1;dims.coils=0;dims.averages=0;dims.subSpecs=0;dims.extras=0;
end

%Now get the size of the data array:
sz=size(fids);

%Now take fft of time domain to get fid:
specs=fftshift(ifft(fids,[],dims.t),dims.t);
    

%Now get relevant scan parameters:*****************************

%Get Spectral width and Dwell Time
dwelltime = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1}*1e-9;  %Franck Lamberton
spectralwidth=1/dwelltime;
    
%Get TxFrq
%txfrq=twix_obj.hdr.Meas.Frequency;
txfrq=twix_obj.hdr.Config.Frequency;

%Get Date
%date = getfield(regexp(twix_obj.hdr.MeasYaps.tReferenceImage0, ...
%'^".*\.(?<DATE>\d{8})\d*"$', 'names'), 'DATE');  %Franck Lamberton

date=''; %The above code for extracting the date from the header 
         %was causing problems.  Since date is not critical
         %for almost any applications, removing it now to be fixed at a
         %later date.

%Find the number of averages.  'averages' will specify the current number
%of averages in the dataset as it is processed, which may be subject to
%change.  'rawAverages' will specify the original number of acquired 
%averages in the dataset, which is unchangeable.
if dims.subSpecs ~=0
    if dims.averages~=0
        averages=sz(dims.averages)*sz(dims.subSpecs);
        rawAverages=averages;
    else
        averages=sz(dims.subSpecs);
        rawAverages=1;
    end
else
    if dims.averages~=0
        averages=sz(dims.averages);
        rawAverages=averages;
    else
        averages=1;
        rawAverages=1;
    end
end

%Find the number of subspecs.  'subspecs' will specify the current number
%of subspectra in the dataset as it is processed, which may be subject to
%change.  'rawSubspecs' will specify the original number of acquired 
%subspectra in the dataset, which is unchangeable.
if dims.subSpecs ~=0
    subspecs=sz(dims.subSpecs);
    rawSubspecs=subspecs;
else
    subspecs=1;
    rawSubspecs=subspecs;
end

% CBF: Include option of sequence svs_slaser_dkd
%Find the number of points acquired before the echo so that this
%information can be stored in the .pointsToLeftshfit field of the data
%structure.  Depending on the pulse sequence used to acquire the data, the
%header location of this parameter is different.  For product PRESS
%sequences, the value is located in twix_obj.image.freeParam(1).  For WIP
%sequences, the value is located in twix_obj.image.cutOff(1,1).  For CMRR
%sequences, the value is located in twix_obj.image.iceParam(5,1).  Special
%thanks to Georg Oeltzschner for decoding all of this and sharing the
%information with me:

if isWIP529 || isWIP859
    leftshift = twix_obj.image.cutOff(1,1);
elseif isSiemens
    leftshift = twix_obj.image.freeParam(1);
%elseif isMinn 
elseif isMinn || isSVSdkdseq
    leftshift = twix_obj.image.iceParam(5,1);
	if isSVSdkdseq
		%disp(sprintf('\n'));
		disp(sMsg_newLine);
		% Not clear whether this option holds for svs_slaser_dkd
		%warning('isSVSdkdseq = %d, check on parameter leftshift = twix_obj.image.iceParam(5,1) = %d\n', isSVSdkdseq, leftshift);
		fprintf('%s: isSVSdkdseq = %d, check on parameter leftshift = twix_obj.image.iceParam(5,1) = %d\n', sFunctionName, isSVSdkdseq, leftshift);
	end
else
    leftshift = twix_obj.image.freeParam(1);
end

% CBF Indicate for eja_svs sequences (isMinn) that left shift of data was already
% performed to extract real/relevant points of all FIDs
% (to avoid that a left shift is incorrectly applied once more outside of this routine)
if isMinn && bLeftshifted
	leftshift = 0;
end

%****************************************************************


%Calculate t and ppm arrays using the calculated parameters:
%f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
%ppm=-f/(Bo*42.577);
%ppm=ppm+4.65;
%Switch between different Nuclei - PT,2021
f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
nucleus=twix_obj.hdr.Config.Nucleus;
switch nucleus
    case '1H'
        gamma=42.576;
        ppm=-f/(Bo*gamma);
        ppm=ppm+4.65;
    case '31P'
        gamma=17.235;
        ppm=-f/(Bo*gamma);
    case '13C'
        gamma=10.7084;
        ppm=-f/(Bo*gamma);
end
t=[0:dwelltime:(sz(1)-1)*dwelltime];


%FILLING IN DATA STRUCTURE
out.fids=fids;
out.specs=specs;
out.sz=sz;
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
out.seq=seq;
out.te=TE/1000;
out.tr=TR/1000;
out.pointsToLeftshift=leftshift;
out.nucleus=nucleus;
out.gamma=gamma;

%FILLING IN THE FLAGS
out.flags.writtentostruct=1;
out.flags.gotparams=1;
% CBF: Indicate whether data were already leftshifted, e.g. for eja sequences from Minn
%out.flags.leftshifted=0;
out.flags.leftshifted=bLeftshifted;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=0;
out.flags.addedrcvrs=0;
out.flags.subtracted=0;
out.flags.writtentotext=0;
out.flags.downsampled=0;
if out.dims.subSpecs==0
    %out.flags.isISIS=0;
	out.flags.isFourSteps=0;
else
    %out.flags.isISIS=(out.sz(out.dims.subSpecs)==4);
	out.flags.isFourSteps=(out.sz(out.dims.subSpecs)==4);
end




%DONE


%% CBF: Create FID-A output data structure for possibly existing reference data
% So far, these reference scans have only been used by Dinesh Deelchand's sequence from
% the CMRR in Minnesota, svs_slaser_dkd; these can also set to zero
if ~isSVSdkdseq || noRefScans == 0
	% Display info
	disp(sMsg_newLine);
	fprintf('%s: Sequence isSVSdkdseq = %d \t # of MRS reference scans = noRefScans = %d, not extracted!', sFunctionName, isSVSdkdseq, noRefScans);    
	%disp(sMsg_newLines);
	fprintf('\n\n');
	
	% Create empyt output data structure for reference data
    %out_ref		= '';
	out_ref		= struct([]);	% Creates a 0x0 struct yielding isempty(out_ref) = 1
else
	% Display info
	disp(sMsg_newLine);
	fprintf('%s: Sequence isSVSdkdseq = %d \t # of MRS reference scans = noRefScans = %d', sFunctionName, isSVSdkdseq, noRefScans);    
	%disp(sMsg_newLines);
	fprintf('\n\n');
	
	% Init output data structure for reference scans by copying output data structure for 
	% regular scans and then replace data and adjust information about data accordingly;
	% most fields of the structure and data flags are the same for the reference scans and
	% the regular MRS scans
	out_ref		= out;
	
	% Squeeze the reference data to remove singleton dimensions
	fids_ref	= squeeze(refData);
	
	% Reorder dimensions for reference scans as performed above for the regular MRS data
	% Reuse stored order of dimensions for that
	% Now that we've indexed the dimensions of the data array, we now need to
	% permute it so that the order of the dimensions is standardized:  we want
	% the order to be as follows:
	%   1) time domain data.
	%   2) coils.
	%   3) averages.
	%   4) subSpecs.
	%   5) extras.
	fids_ref	= permute(fids_ref, dimorder);
	
	% Now get the size of the reference data array
	sz_ref		= size(fids_ref);
	
	% Now take fft of time domain data to get spectra
	specs_ref	= fftshift(ifft(fids_ref,[],dims.t),dims.t);
	
	% Find the number of averages for the reference scans.  
	% 'averages_ref' will specify the current number of averages
	%  in the reference scan dataset as it is processed, which may be subject to
	% change.  'rawAverages_ref' will specify the original number of acquired
	% averages in the reference scan dataset, which is unchangeable.
	% Reference scans are acquired as averages, i.e. # of reference scans = # of
	% averages_ref
	% Information in dims about subSpecs also correct for reference scans
	% (sz_ref(dims.averages) should be equal to noRefScans after permutation of data
	% dimensions)
	if dims.subSpecs ~= 0
		%if dims.averages~=0
		if noRefScans > 0
			averages_ref		= sz_ref(dims.averages)*sz_ref(dims.subSpecs);
			rawAverages_ref		= averages_ref;
		else
			averages_ref		= sz_ref(dims.subSpecs);
			rawAverages_ref		= 1;
		end
	else
		%if dims.averages~=0
		if noRefScans > 0
			averages_ref		= sz_ref(dims.averages);
			rawAverages_ref		= averages_ref;
		else
			averages_ref		= 1;
			rawAverages_ref		= 1;
		end
	end
	
	% Find the number of subspecs.  'subspecs_ref' will specify the current number of 
	% subspectra in the reference scan dataset as it is processed, which may be subject to
	% change.  'rawSubspecs_ref' will specify the original number of acquired
	% subspectra in the reference scan dataset, which is unchangeable.
	% Information in dims about subSpecs also correct for reference scans
	if dims.subSpecs ~= 0
		subspecs_ref		= sz_ref(dims.subSpecs);
		rawSubspecs_ref		= subspecs_ref;
	else
		subspecs_ref		= 1;
		rawSubspecs_ref		= subspecs_ref;
	end
	
	
	% Only adjust specific parameters and flags of output data structure for reference 
	% scans, since reference scans are acquired under same conditions as regular MRS scans
	out_ref.fids		= fids_ref;
	out_ref.specs		= specs_ref;
	out_ref.sz			= sz_ref;
	out_ref.averages	= averages_ref;
	out_ref.rawAverages = rawAverages_ref;
	out_ref.subspecs	= subspecs_ref;
	out_ref.rawSubspecs = rawSubspecs_ref;
	
	% Only adjust specific flags of output data structure for reference scans
	% Flags should all be the same after adjusting parameters in output data structure for
	% reference scans assuming that ISIS acquisition is also used for reference scans
% 	if out_ref.dims.subSpecs == 0
% 		out_ref.flags.isISIS = 0;
% 	else
% 		out_ref.flags.isISIS = (out_ref.sz(out_ref.dims.subSpecs)==4);
% 	end
	
end		% End of if ~isSVSdkdseq
	
	
	
