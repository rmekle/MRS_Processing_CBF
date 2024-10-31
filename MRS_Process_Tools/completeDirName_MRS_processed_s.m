function [strDir_Out] = completeDirName_MRS_processed_s(strDir_Out_Base, fileExt, strVOI, dataType, signals, leftshift, avgBlockSize, rmbadav, noSD, strSpecReg, driftCorr, bECC, strProcessTool)
%UNTITLED2 Summary of this function goes here

% Select directory for output data depending on voxel location, data type,
% # of SDs, and other options used for pre-processing of MR spectra or acquired
% macromolecules (MMs)
% Use variable 'dirSting_out_AddOn1' to include information about the type of
% signals (spectra or MMs), selected voxel location, data type (.dat or .IMA),
% and processing software (e.g. FID-A)
% Use variable 'dirSting_out_AddOn2' to include information about most important
% processing options, preferrably in the order of application

% Make output directories for acquired macromolecules (MMs) distinguishable from
% those for spectra
%if strcmp(signals, 'MMs')
%	strDir_Out_AddOn1	= sprintf('%s_%s_%s_FID-A_SD_%d_%d', signals, strVOI, fileExt, digits(1), digits(2));
%else
%	strDir_Out_AddOn1	= sprintf('%s_%s_FID-A_SD_%d_%d', strVOI, fileExt, digits(1), digits(2));
%end		% End of if strcmp(signals, 'MMs')
if strcmp(signals, 'MMs')
	strDir_Out_AddOn1	= sprintf('%s_%s_%s_%s', signals, strVOI, fileExt, strProcessTool);
else
	strDir_Out_AddOn1	= sprintf('%s_%s_%s', strVOI, fileExt, strProcessTool);
end		% End of if strcmp(signals, 'MMs')

% Init information about processing of MRS data
strDir_Out_AddOn2	= '';
% Leftshifting of data (to cut off points before true first point of FID)
% Include info about Leftshifting of data, only if applied
if leftshift > 0
	strDir_Out_AddOn2	= [strDir_Out_AddOn2, sprintf('_ls%d', leftshift)];
end		% End of if leftshift > 0

% Block averaging prior to processing to improve SNR
% Indicate in output directory name which size of block averaging (Bavg) was used,
% only if averaging of blocks was indeed performed
if avgBlockSize > 0
	strDir_Out_AddOn2	= [strDir_Out_AddOn2, sprintf('_Bavg%d', avgBlockSize)];
end

% Removal of bad averages
% Always include info about removal of bad averages independent of whether it was
% performed or not
digits_noSD			= [fix(noSD) round(abs(noSD-fix(noSD))*10)];
if strcmpi(rmbadav, 'y')	% Case-insensitive strcmp
	strDir_Out_AddOn2	= [strDir_Out_AddOn2, sprintf('_SD%d_%d', digits_noSD(1), digits_noSD(2))];
else
	strDir_Out_AddOn2	= [strDir_Out_AddOn2, '_NoRM'];
end		% End of if strcmpi(rmbadav, 'y')

% Spectral registreation / drift correction
% Always include info about spectral registration independent of whether it was
% performed or not
if strcmpi(driftCorr, 'y')	% Case-insensitive strcmp
	strDir_Out_AddOn2	= [strDir_Out_AddOn2, '_', strSpecReg];
else
	strDir_Out_AddOn2	= [strDir_Out_AddOn2, '_NoSR'];
end		% if strcmpi(driftCorr, 'y')

% Eddy current correction (ECC)
% Include info about ECC, only if applied
if bECC
	% Use reference (water) signals for ECC, if acquired
	% If not, then use an unsuppressed water signal, if acquired
	% If no reference and no water signals are acquired, check whether MR spectrum
	% is water signal itself; and if it is, use it for ECC
	% Indicate different options for ECC in the corresponding directory name; for
	% that, search for strings 'ref', 'w', and 'water' in string for MRS data type
	%refInd		= strfind(dataType, '_ref');
	%wInd		= strfind(dataType, '_w');
	%waterInd	= strfind(dataType, 'water');
	%if ~isempty(refInd)
	if contains(dataType, '_ref')
		strDir_Out_AddOn2	= [strDir_Out_AddOn2, '_ECCref'];
	else
		%if ~isempty(wInd)
		if contains(dataType, '_w')
			strDir_Out_AddOn2	= [strDir_Out_AddOn2, '_ECCw'];
		else
			%if ~isempty(waterInd)
			if contains(dataType, 'water')
				strDir_Out_AddOn2	= [strDir_Out_AddOn2, '_ECCwater'];
			else
				% No reference and no water signals and MR spectrum is not water
				% signal itself => ECC not possible
				error('%s: No reference and no water signals and MR spectrum is not water signal itself (dataType_MRS = %s) => ECC not possible!', sFunctionName, dataType_MRS);
			end		% End of if ~isempty(waterInd)
		end		% End of if ~isempty(wInd)
	end		% End of if ~isempty(refInd)
end		% End of if bECC

% Complete output directory name for processed MRS data
strDir_Out			= [strDir_Out_Base, strDir_Out_AddOn1, strDir_Out_AddOn2, filesep];

end		% End of function