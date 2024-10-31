%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% completeDirName_MRS_processed_s.m
%
%% Function to preprocess single volume magnetic resonance spectroscopy (MRS) data
%  using functions from the MRS processing toolkit FID-A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USAGE
% [strDir_Out] = completeDirName_MRS_processed_s(strDir_Out_Base, fileExt, strVOI, dataType, signals, leftshift, avgBlockSize, rmbadav, noSD, strSpecReg, driftCorr, bECC, strProcessTool)
% 
% DESCRIPTION:
% Function for processing Siemens MRS data in .dat format (twix raw data) or in .IMA
% format (DICOM data)
% using functions from the MRS processing toolkit FID-A
% Includes combination of receiver coil channels (if required), removal of bad averages, 
% frequency drift correction, averaging, eddy current correction (optional), phase and 
% frequency correction.
% 
% INPUTS:
% dirString    = String variable for the name of the directory containing
%                   the water suppressed .dat file or .IMA files
% dirString_w  = (Optional) ['DirectoryWater'] String variable for the name of the directory containing
%                   the water unsuppressed .dat file or .IMA files,
%                   optional because water unsupressed data is optional and
%                   dat file can also be in dirString
% outDirString = String variable for the name of the output directory, i.e. the directory,
%					where all output files are saved to
% filename     = (Optional) ['Filename'] String variable for the name of the water suppressed .dat file,
%					e.g. 3T_20170510_PetraO_meas_MID00091_FID153257_rm_special_RF_ACC.dat
%					or specialDLPFC.dat - Optional because not necessary
%					for .IMA data
% filename_w   = (Optional) ['UnsupressedFilename'] String variable for the name of the water unsuppressed .dat file,
%					e.g. 3T_20170510_PetraO_meas_MID00091_FID153259_rm_special_RF_ACC.dat
%					or specialDLPFC_w.dat - Optional because water
%					unsupressed signal is optional and water unsupressed
%					can be a directory containing .IMA files
% filename_r   = (Optional) ['ReportFilename'] String variable for the name
%                   of the file - Optional because the name will just be
%                   automatically generated otherwise
% seqType	   = String specifying the MRS sequence type used for data acquisition, e.g.
%					'PRESS', 'STEAM', 'SPECIAL', 'sLASER', 'MEGA-PRESS'
% dataType	   = String describing the type of MRS data:
%					'mrs'		= MR spectrum without water signal, 
%					'mrs_w'		= MR spectrum with unsuppressed water signal
%					'mrs_w_ref'	= MR spectrum with unsuppressed water signal and reference
%									(water) scans
%					'mrs_ref'	= MR spectrum with reference (water) scans
%					'water'		= MR spectrum is unsuppressed water signal itself
%					'water_ref' = MR spectrum is unsuppressed water signal itself with
%									reference (water) scans (should be very rare!)
% strOVS	   = (optional) ['OVS'] Character array that specifies whether MR spectrum acquired 
%					 with OVS ('wOVS') or withoutOVS ('woutOVS') is used for processing. 
%					 Default is 'woutOVS', which means that OVS was not used.
% strOVS_w	   = (optional) ['OVSwater'] Character array that specifies whether water acquired with OVS ('wOVS')
%					 or withoutOVS ('woutOVS') is used for processing.
%					 Default is 'woutOVS', which means that OVS was not used.
% leftshift	   = (optional) ['Leftshift'] # of points to leftshift all FIDs, i.e. to remove leading 
%					 datapoints from the FID  to get rid of 1st order phase. Default is 0.
% leftshift_w  = (optional) ['WaterLeftshift'] # of points to leftshift all water unsupressed FIDs, i.e. to remove leading 
%					 datapoints from the FID  to get rid of 1st order phase. Default is the same as leftshift.
% avgBlockSize = (Optional) ['avgBlockSize'] Block size of averages used for block
%					averaging prior to processing for noisy data. If equal to 0, block
%					averaging is not applied. Default value is 0.
% rmbadav	   = (Optional) ['RemoveBadAverages'] Character array that specifies whether 
%					removal of bad averages should be performed or not. 
%					Default is 'y'.
% nSD		   = (Optional) ['StandardDeviation'] # of standard deviations for bad average removal. Default
%					value is 3.2.
% strSpecReg   = (Optional) ['SpectralRegistrationID'] Character array that specifies ID,
%					i.e. name, of spectral registration (drift correction) that might be 
%					performed to distinguish results. Default is 'SR00'.
% driftCorr	   = (Optional) ['DriftCorrection'] Character array that specifies whether 
%					spectral registration (drift correction) should be performed or not. 
%					Default is 'y'.
% iterin       = (Optional) ['Iterations']  Maximum number of allowed iterations for the spectral
%                   registration to converge. Default is 20.
% aaDomain     = (Optional) ['aaDomain'] Perform the spectral registration (drift correction) using
%                   the full spectrum ('t'), or only a limited frequency range ('f').  Default is 'f'.
% tmaxin       = (Optional) ['DriftCorrectionDuration'] Duration (in sec.) of the time domain signal
%                   used in the spectral registration (drift correction).
%                   Default is 0.2 sec.
% bTmaxSet	   = (Optional) ['MaxTimeAlignmentSet'] Boolean, if 1, the given value for
%					tmaxin is used; if 0, the routine that aligns all averages determines
%					tmax from the data (currently only implemented for the time domain).
%					Default is 1.
% medin		   = (Optional) ['medianALignment] Selects reference, to which all averages
%					are aligned to; 'y', 'n', 'a' or 'ref' are possible; see
%					op_alignAverages_fd.m and op_alignAverages.m from FID-A for details.
%					Default is 'y'.
% ppmmin_fix   = (Optional) ['ppmMinimum_fix'] Initial minimum ppm value, for which
%					spectral registration is applied in the frequency domain. 
%					Default is 1.6.
% ppmmaxarreay_fix = (Optional) ['ppmMaximumArray_fix'] Initial array of maximum ppm
%					values, for which spectral registration is applied in the frequency
%					domain. Default is [3.5; 4.0; 5.5].
% bECC 		   = (optional) ['ECC'] Boolean that specifies whether eddy current correction (ECC) 
%					 should be performed or not. Default is 0.
% bPhaseCorrFreqShift = (Optional) ['PhaseFrequencyCorrection'] Boolean that specifies whether phase correction and
%					frequency shifting should be performed or not. Default is 0.
% strMinUserIn = (Optional) ['MinimizeUserInput'] String that specifies whether user input/interaction should be
%					 minimized or not; 'y' or 'Y' lead to minimization, 'n' or 'N' do not
%					Default is 'y'.
% plotSwitch   = (Optional)	['ShowPlots'] Switch for displaying plots: 1 = ON, 0 = OFF. Default is 0
% reportSwitch = (Optional) ['GenerateReport'] Switch for generating an html report with corresponding  
%					figures and a readme file: 1 = ON, 0 = OFF. Default is 1. 
% 
% OUTPUTS:
% strDir_Out		= Fully preprocessed, water suppressed output spectrum
%
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2024;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [strDir_Out] = completeDirName_MRS_processed_s(strDir_Out_Base, fileExt, strVOI, dataType, signals, leftshift, avgBlockSize, rmbadav, noSD, strSpecReg, driftCorr, bECC, strProcessTool)

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
				error('%s: No reference and no water signals and MR spectrum is not water signal itself (dataType = %s) => ECC not possible!', sFunctionName, dataType);
			end		% End of if contains(dataType, 'water')	%if ~isempty(waterInd)
		end		% End of if contains(dataType, '_w')	%if ~isempty(wInd)
	end		% End of if contains(dataType, '_ref')	%if ~isempty(refInd)
end		% End of if bECC

% Complete output directory name for processed MRS data
strDir_Out			= [strDir_Out_Base, strDir_Out_AddOn1, strDir_Out_AddOn2, filesep];

end		% End of function