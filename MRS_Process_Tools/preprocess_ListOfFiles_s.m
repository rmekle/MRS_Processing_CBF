%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% preprocess_ListOfFiles_s.m
%
%% Script to preprocess a list of files of magnetic resonance spectroscopy (MRS) data
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2018, 2019, 2020, 2021, 2022,
% 2023, 2024, 2025; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'preprocess_ListOfFiles_s';
%sMsg_newLines		= sprintf('\n\n');	%disp(sMsg_newLines);
fprintf('\n\n');


%% Init input parameters for preprocessing
% Obtain parameter settings for preprocessing from initialization routine
%configSel				= 'config_Study_sLASER_VOI_dat_MRS_lsN_SDx_y_SC1_ECC';
%configSel				= 'config_Study_sLASER_VOI_IMA_MMs_lsN_SDx_y_SR8_ECC';
configSel				= 'config_Study_sLASER_VOI_IMA_MMs_lsN_SDx_y_SC3_ECC';
[sParamsMRS_struct]		= initParams_MRS_s(configSel);

% Set addition to starting index into list of files to be preprocessed; 
% indexAdd				= 0;	% To start with first file in list of files
% indexAdd				= n;	% To start with (n+1)th file in list of files
indexAdd				= 0;

% Extract parameter settings from parameter struct
filename_In				= sParamsMRS_struct.filename;
filename_w_In			= sParamsMRS_struct.filename_w;
strStudy_MRS			= sParamsMRS_struct.strStudy_MRS;
seqType_MRS				= sParamsMRS_struct.seqType_MRS;
strVOI_MRS				= sParamsMRS_struct.strVOI_MRS;
fileExt_MRS				= sParamsMRS_struct.fileExt_MRS;
dataType_MRS			= sParamsMRS_struct.dataType_MRS;
signals_MRS				= sParamsMRS_struct.signals_MRS;
strOVS_In				= sParamsMRS_struct.strOVS;
strOVS_w_In				= sParamsMRS_struct.strOVS_w;
leftshift_In			= sParamsMRS_struct.leftshift;
avgBlockSize_In			= sParamsMRS_struct.avgBlockSize;

% Info about processing tool(s) mainly used
strProcessTool_In		= sParamsMRS_struct.strProcessTool;

% Parameters for removal of bad averages
rmbadav_In				= sParamsMRS_struct.rmbadav;
noSD_In					= sParamsMRS_struct.noSD;
%digits_noSD_In			= [fix(noSD_In) round(abs(noSD_In-fix(noSD_In))*10)];

% Parameters for aligning of averages/frequency and phase drift correction using
% one of the following techniques:
%	Spectral registration	performed in either frequency or time domain or
%	Cross-Correlation		performed in frequency domain
% Parameter to select frequency and phase drift correction method
driftCorr_In			= sParamsMRS_struct.driftCorr;
strFreqPhaseCorr_In		= sParamsMRS_struct.strFreqPhaseCorr;

% Extract setting to align subspectra independent of selected correction technique
alignSS_In				= sParamsMRS_struct.alignSS;	% For aligning subspectra (e.g. in SPECIAL)

% Extract settings for frequency and phase drift correction techniques into separate 
% structs for easier use
% Parameters for spectral registration (SR)
structSR_In				= sParamsMRS_struct.structSR;
% Parameters for spectral cross-correlation (SC)
structSC_In				= sParamsMRS_struct.structSC;


% driftCorr_In			= sParamsMRS_struct.driftCorr;
% iterin_In				= sParamsMRS_struct.iterin;
% aaDomain_In				= sParamsMRS_struct.aaDomain;
% tmaxin_In				= sParamsMRS_struct.tmaxin;
% bTmaxset_In				= sParamsMRS_struct.bTmaxset;
% ppmOption				= sParamsMRS_struct.ppmOption;
% medin_In				= sParamsMRS_struct.medin;
% Obtain parameters for drift correction depending on type of data, i.e. whether MRS
% data is spectrum or water signal
% NOTE: Check whether aligning of averages in frequency domain works, if the MR
% spectrum is water signal itself; if not, simply align averages in time domain
%ppmmin_fix_In			= sParamsMRS_struct.ppmmin_fix;
%ppmmaxarray_fix_In		= sParamsMRS_struct.ppmmaxarray_fix;

% Additional parameter settings
bECC_In					= sParamsMRS_struct.bECC;
bPhaseCorrFreqShift_In	= sParamsMRS_struct.bPhaseCorrFreqShift;
strMinUserIn_In			= sParamsMRS_struct.strMinUserIn;
plotSwitch_In			= sParamsMRS_struct.plotSwitch;
reportSwitch_In			= sParamsMRS_struct.reportSwitch;
bPrep_MetabQuant		= sParamsMRS_struct.bPrep_MetabQuant;


%% Additional input parameters specific to preparation of metabolite quantification
% After preprocessing of MRS data files, specific groups of these processed MRS data files
% are then copied into output directory for metabolite quantification using LCM analysis
% and corresponding list of filenames are written into text files
% 'bCopyFiles' used to turn on/off any copying of files (mostly used for debugging)
bWriteFilenames_In			= 1;
bCopyFiles_In				= 1;
%bCopyFiles_MRS_In			= 1;
%bCopyFiles_ref_Quant_In	= 0;
bCopyFiles_ref_ECC_In		= 0;
%bCopyFiles_w_In			= 1;
switch dataType_MRS
	case {'mrs_w_ref', 'mrs_ref'}
		% MR spectrum was acquired with reference scans that are typically used for
		% quantification
		bCopyFiles_MRS_In			= 1;
		bCopyFiles_ref_Quant_In		= 1;
		bCopyFiles_w_In				= 0;
	case {'mrs_w'}
		% MR spectrum was acquired without reference scans, but with water signals
		bCopyFiles_MRS_In			= 1;
		bCopyFiles_ref_Quant_In		= 0;
		bCopyFiles_w_In				= 1;
	case {'mrs'}
		% MR spectrum was acquired without reference scans and without water signals
		bCopyFiles_MRS_In			= 1;
		bCopyFiles_ref_Quant_In		= 0;
		bCopyFiles_w_In				= 0;
	case {'water', 'water_ref'}
		% MR spectrum is water signal itself without or with reference scans
		% (not really used for quantification)
		bCopyFiles_MRS_In			= 0;
		bCopyFiles_ref_Quant_In		= 0;
		bCopyFiles_w_In				= 1;

	otherwise
		error('%s: Unknown MRS dataType_MRS = %s!', sFunctionName, dataType_MRS);
end		% End of switch dataType_MRS


%% Set (additional) parameters depending on sequence type
switch seqType_MRS
	case 'SPECIAL'
		dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_RawData_dat_Files/';
		dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Z_Pain_Tmp/';
		%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Preprocessed_forLCModel_SD_4_0/';
		%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Preprocessed_forLCModel_SD_3_2/';
		%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Preprocessed_forLCModel_SD_2_6/';
		% Output directory for metabolite quantification using LCM analysis
		dirString_Out_LCM		= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Z_Pain_Tmp/';
	case 'MEGA-PRESS'
		% Data (input) directories
		dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopa/MRS_Dopamin_00_All_RawData_dat_Files_MRS/';
		%dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopa/MRS_Dopamin_00_All_RawData_dat_Files_MRS_New/';
		
		% Select directories for output data depending on # of SDs used for pre-processing
		% of MR spectra
		switch(noSD_In)
			case(2.6)
				dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_2_6/';
			case(3.2)
				%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_3_2/';
				%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_3_2_II/';
				dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_3_2_III/';
				%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_3_2_New/';
				%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2_New/';
			case(4.0)
				dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_4_0/';
				%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_4_0/';
			
			otherwise
				error('%s: ERROR: No directory/data for noSD_In =  %f!', sFunctionName, noSD_In);
		end
		% Output directory for metabolite quantification using LCM analysis
		dirString_Out_LCM			= [dirString_Out, 'DOPA_LCModel_Analysis_Data/'];
	case 'sLASER'
		% Select data input and output directories depending on study, MRS data type, 
		% i.e. file extension, study, and other parameters
		%digits_noSD_In		= [fix(noSD_In) round(abs(noSD_In-fix(noSD_In))*10)];
		switch strStudy_MRS
			case '3T_Trauma'
				% Data (input) directories
				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/';
				%dirString_In_AddOn1		= sprintf('MRS_Trauma_00_All_RawData_dat_Files_MRS_%s', strVOI_MRS);
				
				% Output data directory
				dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Trauma_Analysis/';
				
				% Directories depending on MRS data type
				switch fileExt_MRS
					case 'dat'
						% Select directories specific to MRS raw data (.dat)
						dirString_In_AddOn1		= sprintf('MRS_Trauma_00_All_RawData_dat_Files_MRS_%s', strVOI_MRS);
					case 'IMA'
						% Select directories specific to MRS DICOM data (.IMA)
						dirString_In_AddOn1		= sprintf('MRS_Trauma_00_All_DICOM_IMA_Files_MRS_%s', strVOI_MRS);
						
					otherwise
						error('%s: ERROR: Unknown file extension (data type) %s!', sFunctionName, fileExt_MRS);
				end			% End of switch fileExt_MRS
			case '3T_SBAM'
				% Data (input) directories
				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/SBAM/';
				
				% Output data directory
				dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_SBAM_Analysis/';
				
				% Directories depending on MRS data type
				switch fileExt_MRS
					case 'dat'
						% Select directories specific to MRS raw data (.dat)
						dirString_In_AddOn1		= sprintf('MRS_SBAM_00_All_RawData_dat_Files_MRS_%s', strVOI_MRS);
					case 'IMA'
						% Select directories specific to MRS DICOM data (.IMA)
						dirString_In_AddOn1		= sprintf('MRS_SBAM_00_All_DICOM_IMA_Files_MRS_%s', strVOI_MRS);
						
					otherwise
						error('%s: ERROR: Unknown file extension (data type) %s!', sFunctionName, fileExt_MRS);
				end			% End of switch fileExt_MRS
			case '3T_MMs_SBA'
				% Data (input) directories
				%dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/';
				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MMs/';
				%dirString_In_Base		= '/home/mekler/CSB_NeuroRad/destiana/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MMs/';
				%dirString_In_Base		= '/home/destiana/CSB_NeuroRad/destiana/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MMs/';

				% Output data directory
				%dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Trauma/Trauma_Z_Analysis/';
				%dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/ZZZZ_Test/';
				%dirString_Out_Base		= '/home/destiana/CSB_NeuroRad/destiana/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MMs/CodeResults/';
				dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Trauma_MMs_Analysis/';

				% Directories depending on MRS data type
				switch fileExt_MRS
					case 'dat'
						% Select directories specific to MRS raw data (.dat)
						% either for acquired macromolecules (MMs) or for spectra
						if strcmp(signals_MRS, 'MMs')
							dirString_In_AddOn1		= 'MMs_dat';
						else
							dirString_In_AddOn1		= 'Spectra_dat';
						end		% End of if strcmp(signals_MRS, 'MMs')
					case 'IMA'
						% Select directories specific to MRS DICOM data (.IMA)
						% either for acquired macromolecules (MMs) or for spectra
						if strcmp(signals_MRS, 'MMs')
							dirString_In_AddOn1		= 'MMs_IMA';
						else
							dirString_In_AddOn1		= 'Spectra_IMA';
						end		% End of if strcmp(signals_MRS, 'MMs')

					otherwise
						error('%s: ERROR: Unknown file extension (data type) %s!', sFunctionName, fileExt_MRS);
				end			% End of switch fileExt_MRS
			case '3T_MMs_SBAM'
				% Data (input) directory
				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MMs_SBAM/';
				
				% Output data directory
				dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_SBAM_MMs_Analysis/';
				
				% Directories depending on MRS data type
				switch fileExt_MRS
					case 'dat'
						% Select directories specific to MRS raw data (.dat)
						% for acquired macromolecules (MMs)
						dirString_In_AddOn1		= 'MMs_SBAM_dat';
					case 'IMA'
						% Select directories specific to MRS DICOM data (.IMA)
						% for acquired macromolecules (MMs)
						dirString_In_AddOn1		= 'MMs_SBAM_IMA';

					otherwise
						error('%s: ERROR: Unknown file extension (data type) %s!', sFunctionName, fileExt_MRS);
				end			% End of switch fileExt_MRS
			case '3T_TGA'
				% Data (input) directories
				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_MRS_TGA/';

				% Output data directory
				dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_MRS_TGA_Analysis/';

				% Directories depending on MRS data type
				switch fileExt_MRS
					case 'dat'
						% Select directories specific to MRS raw data (.dat)
						dirString_In_AddOn1		= sprintf('MRS_TGA_00_All_RawData_dat_Files_MRS_%s', strVOI_MRS);
					case 'IMA'
						% Select directories specific to MRS DICOM data (.IMA)
						dirString_In_AddOn1		= sprintf('MRS_TGA_00_All_DICOM_IMA_Files_MRS_%s', strVOI_MRS);

					otherwise
						error('%s: ERROR: Unknown file extension (data type) %s!', sFunctionName, fileExt_MRS);
						
				end			% End of switch fileExt_MRS
			case '3T_BPAPS'
				% Data (input) directories
				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_BPAPS/BPAPS_BCAN/';

				% Output data directory
				dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_BPAPS_Analysis/';

				% Directories depending on MRS data type
				switch fileExt_MRS
					case 'dat'
						% Select directories specific to MRS raw data (.dat)
						dirString_In_AddOn1		= sprintf('MRS_BPAPS_00_All_RawData_dat_Files_MRS_%s', strVOI_MRS);
					case 'IMA'
						% Select directories specific to MRS DICOM data (.IMA)
						dirString_In_AddOn1		= sprintf('MRS_BPAPS_00_All_DICOM_IMA_Files_MRS_%s', strVOI_MRS);

					otherwise
						error('%s: ERROR: Unknown file extension (data type) %s!', sFunctionName, fileExt_MRS);
						
				end			% End of switch fileExt_MRS				
			case '7T_KCL'
				% Data (input) directories
				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/7T_KCL/';
				%dirString_In_AddOn1		= sprintf('7T_KCL_00_ALL_RawData_dat_Files_MRS_eja_%s', strVOI_MRS);
				
				% Output data directory
				dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/7T_KCL_Analysis/';
				
				% Directories depending on MRS data type
				switch fileExt_MRS
					case 'dat'
						% Select directories specific to MRS raw data (.dat)
						dirString_In_AddOn1		= sprintf('7T_KCL_00_ALL_RawData_dat_Files_MRS_eja_%s', strVOI_MRS);
					case 'IMA'
						% Select directories specific to MRS DICOM data (.IMA)
						dirString_In_AddOn1		= sprintf('7T_KCL_00_ALL_DICOM_IMA_Files_MRS_eja_%s', strVOI_MRS);
						
					otherwise
						error('%s: ERROR: Unknown file extension (data type) %s!', sFunctionName, fileExt_MRS);
				end			% End of switch fileExt_MRS
			case '3T_Test'
				% Data (input) directories
				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/Z_Test_Data/3T_Test_MRS/';

				% Output data directory
				dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/Z_Test_Data_Analysis/3T_Test_MRS_Analysis/';

				% Directories depending on MRS data type
				switch fileExt_MRS
					case 'dat'
						% Select directories specific to MRS raw data (.dat)
						dirString_In_AddOn1		= sprintf('MRS_TGA_00_All_RawData_dat_Files_MRS_%s', strVOI_MRS);
					case 'IMA'
						% Select directories specific to MRS DICOM data (.IMA)
						dirString_In_AddOn1		= sprintf('MRS_TGA_00_All_DICOM_IMA_Files_MRS_%s', strVOI_MRS);

					otherwise
						error('%s: ERROR: Unknown file extension (data type) %s!', sFunctionName, fileExt_MRS);
				end			% End of switch fileExt_MRS

			otherwise
				error('%s: ERROR: Unknown study %s!', sFunctionName, strStudy_MRS);
		end				% End of switch strStudy_MRS
		
		% Complete names of data (input) directories
		%dirString_In_AddOn1		= [dirString_In_AddOn1, '_Test'];
		dirString_In			= [dirString_In_Base, dirString_In_AddOn1, filesep];
			
		% Select directory for output data depending on voxel location, data type,
		% # of SDs, and other options used for pre-processing of MR spectra or acquired
		% macromolecules (MMs)
		% Use variable 'dirSting_out_AddOn1' to include information about the type of
		% signals (spectra or MMs), selected voxel location, data type (.dat or .IMA),
		% and processing software (e.g. FID-A)
		% Use variable 'dirSting_out_AddOn2' to include information about most important
		% processing options, preferrably in the order of application
		
		% Complete output data directory name for preprocessed MRS data
		%dirString_Out			= [dirString_Out_Base, dirString_Out_AddOn1, dirString_Out_AddOn2, filesep];
		%dirString_Out			= [dirString_Out_Base, dirString_Out_AddOn1, dirString_Out_AddOn2, '_Test', filesep];
		dirString_Out			= completeDirName_MRS_processed_s(dirString_Out_Base, strVOI_MRS, fileExt_MRS, dataType_MRS, signals_MRS, leftshift_In, avgBlockSize_In, ...
			strProcessTool_In, rmbadav_In, noSD_In, strFreqPhaseCorr_In, driftCorr_In, bECC_In);

		% If directory for results from preprocessing does not exist, create it
		% else, if it exists, check whether it can be overwritten
		%if ~exist( dirString_Out, 'dir' )
		if not(isfolder(dirString_Out))	% Preferred over if ~exist(...) according to MATLAB
			fprintf('%s: Creating new output directory %s ...\n\n\n', sFunctionName, dirString_Out);
			if ~mkdir(dirString_Out)
				error('%s: Could not create (mkdir) new output directory %s!\n', sFunctionName, dirString_Out);
			end
		else
			% Output directory already exists, ask user whether to overwrite or not
			% (should help to avoid accidentally overwriting previously processed data)
			prompt			= sprintf('\n\nOutput directoy = %s\nDo you want to overwrite the existing output directory (y/n)?  ', dirString_Out);
			strOverwrite	= input(prompt, 's');
			if strcmp(strOverwrite, 'n') || strcmp(strOverwrite, 'N')
				fprintf('\n%s: Already existing output directory is not overwritten! Preprocessing aborted!\n\n\n', sFunctionName)
				return;
			else
				fprintf('\n%s: Already existing output directory is overwritten! Preprocessing is continued!\n\n\n', sFunctionName)
			end		% End of if strOverwrite == 'n' || strOverwrite  == 'N'
		end		% End of iif not(isfolder(dirString_Out))
		% Output directory for metabolite quantification using LCM analysis
		dirString_Out_LCM			= [dirString_Out, strVOI_MRS, '_LCModel_Data/'];

	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS);
end		% End of switch seqType_MRS



%% Obtain information about the list of files (.dat) or list of directories (.IMA)
% (assuming that all .dat files are included in the same directory)
% (On Linux, file list in Matlab also includes the two directories "." and "..", which
% means that the actual # of files in the directory is (# of entries in list - 2;
% however, if dir is used to list specific files, e.g. using a file extension, these two 
% directories are not included in the resulting list)
%cd(dirString_In);

% FLAG: CHANGE
% Added an option for loading a set of directories containing scans
switch fileExt_MRS
    case 'dat'
        structFileListing		= dir([dirString_In, '*.dat']);
        noEntriesListing		= length( structFileListing );
        %noDataFiles				= noEntriesListing - 2
    case 'IMA'
        structFileListingAll	= dir(dirString_In);
        subDir					= [structFileListingAll(:).isdir];
        structFileListing		= structFileListingAll(subDir);
        % Remove the two directories '.' and '..'
        structFileListing		= structFileListing(~ismember({structFileListing(:).name},{'.','..'}));
        noEntriesListing		= length(structFileListing);
    otherwise
        error('%s: ERROR: Unknown file extension %s!', sFunctionName, fileExt_MRS);
end



%% Preprocess all data files depending on sequence type

% FLAG: TODO: Update for all 3 cases the function to the new version with
% name-value pair parameters, depending on if imaDataSwitch is set or not


switch seqType_MRS
	case 'SPECIAL'
		% Pre-processing of .IMA files is not implemented for run_specialproc_CBF, only
        % for the sLASER sequence using preProcess_MRS_s
		if strcmp(fileExt_MRS, 'IMA')
			error('%s: ERROR: File extension %s incompatible with sequence type %s!', sFunctionName, fileExt_MRS, seqType_MRS);
		end

		% Here preprocessing of SPECIAL MR spectra together with the corresponding water 
		% files is performed
		% Assumptions:
		% - All data files are consecutively sorted, e.g. by date
		% - Data files come in group of three files:
		% - SPECIAL MR spectrum, SPECIAL water withOVS (wOVS), SPECIAL water withoutOVS (woutOVS)
		% - Order of files is same for all cases, i.e. spectrum; water_wOVS, water_woutOVS
		
		% Preprocess each case (spectrum) once using water_wOVS and water_woutOVS
		% Water signal is used for coil combination of the spectrum and is preprocessed 
		% as well
		%indexStart		= 3;	% To skip entries for directories "." and ".."
		indexStart		= 1;
		indexStep		= 3;	% For SPECIAL, since two different water signals exist
		for ind=indexStart : indexStep : 3		% noEntriesListing	% 6		% 3
			% Preprocess MR spectrum and water_withOVS
			strOVS_In		= 'wOVS';
			filename_In		= structFileListing(ind).name;
			filename_w_In	= structFileListing(ind+1).name;
			%disp(sMsg_newLines);
			fprintf('\n\n');
			disp([sprintf('ind = %d\t', ind), strOVS_In, sprintf('\t'), filename_In, sprintf('\t'), filename_w_In, sprintf('\n\n')]);
			[out,out_w,out_noproc,out_w_noproc]=run_specialproc_CBF(dirString_In,dirString_Out,filename_In,filename_w_In,noSD_In,strOVS_In,strMinUserIn_In,structSR_In.aaDomain,structSR_In.tmaxin,structSR_In.iterin);
			
			% Close all figures
			close all;
			
			% Preprocess MR spectrum and water_woutOVS
			% (spectral data remain the same and only processed differently now)
			strOVS_In		= 'woutOVS';
			filename_w_In	= structFileListing(ind+2).name;
			%disp(sMsg_newLines);
			fprintf('\n\n');
			disp([sprintf('ind = %d\t', ind), strOVS_In, sprintf('\t'), filename_In, sprintf('\t'), filename_w_In, sprintf('\n\n')]);
			[out,out_w,out_noproc,out_w_noproc]=run_specialproc_CBF(dirString_In,dirString_Out,filename_In,filename_w_In,noSD_In,strOVS_In,strMinUserIn_In,structSR_In.aaDomain,structSR_In.tmaxin,structSR_In.iterin);
		end
	case 'MEGA-PRESS'
		% Pre-processing of .IMA files is not implemented for run_megapressproc_CBF, only
        % for the sLASER sequence using preProcess_MRS_s
		if strcmp(fileExt_MRS, 'IMA')
			error('%s: ERROR: File extension IMA incompatible with sequence type %s!', sFunctionName, seqType_MRS);
		end
		
		% Here preprocessing of MEGA-PRESS MR spectraFSN-FZ-CSB-08 together with the corresponding 
		% water file is performed
		% Assumptions:
		% - All data files are consecutively sorted, e.g. by date
		% - Data files come in group of two files:
		% - MEGA-PRESS MR spectrum, MEGA-PRESS water signal
		% - Order of files is same for all cases, i.e. spectrum; water
		
		% Preprocess each case (spectrum)
		% Water signal is used for coil combination of the spectrum and is preprocessed 
		% as well
		indexStart		= 1;
		indexStep		= 2;	% For MEGA-PRESS, since only one water signals exists
		for ind=indexStart : indexStep : noEntriesListing		% noEntriesListing	% 4		% 2
			% Preprocess MR spectrum and water
			filename_In		= structFileListing(ind).name;
			filename_w_In	= structFileListing(ind+1).name;
			%disp(sMsg_newLines);
			fprintf('\n\n');
			disp([sprintf('ind = %d\t', ind), sprintf('\t'), filename_In, sprintf('\t'), filename_w_In, sprintf('\n\n')]);
			[diffSpecOut,sumSpecOut,subSpec1Out,subSpec2Out,outwOut,outw_subSpec1Out,outw_subSpec2Out,coilcombosOut]=run_megapressproc_CBF(dirString_In,dirString_Out,filename_In,filename_w_In,noSD_In,strMinUserIn_In,structSR_In.aaDomain,structSR_In.tmaxin,structSR_In.iterin,alignSS_In);
			
			% Close all figures
			%close all;
		end
	case 'sLASER'
		% Here preprocessing of sLASER MR spectra acquired with or without water reference
		% signals together with the corresponding water file is performed
		% Assumptions for MRS raw data (.dat) files:
		% - All data files are consecutively sorted, e.g. by date
		% - Data files come in group of two files:
		% - sLASER MR spectrum with or without reference scans, sLASER water signal
		% - Order of files is same for all cases, i.e. spectrum; water
		
		% Assumptions for MRS DICOM data (.IMA) files:
		% - All data directories are consecutively sorted, e.g. by date
		% - Data directories come in group of two directories:
		% - sLASER MR spectrum with or without reference scans, sLASER water signal
		% - Order of directories is same for all cases, i.e. spectrum; water
		
		% Preprocess each case (spectrum)
		% If separate water scans and reference scans were acquired, get coil phases from
		% both types of signals, and for coil combination use 
		% if reference scans exist,
		%	coil phases from reference scans for reference scans and MR spectra, since
		%	reference scans were acquired together with MR spectra
		%	coil phases from water scans for water scans
		%
		% if only water scans exist, 
		%	coil phases from water scans for water scans and for MR spectra
		%
		% if neither reference scans nor water scans exist,
		%	coil phases from MR spectra for MR spectra
		%
		% Select size for stepping through indices, i.e. list of files (.dat) or list of
		% directoried (.IMA), depending on data type, i.e. how many different signals 
		% (spectra and/or water signals) are included
		indexStart		= 1 + indexAdd;
		indexStep		= 2;	% Default for sLASER spectrum with one water signal
		switch dataType_MRS
			case {'mrs_w', 'mrs_w_ref'}
				% Spectra and water signals in list of files/directories
				indexStep		= 2;	
			case {'mrs', 'mrs_ref', 'water', 'water_ref'}
				% Only spectra or only water signals in list of files/directories
				indexStep		= 1;

			otherwise
				error('%s: Unknown MRS dataType_MRS = %s!', sFunctionName, dataType_MRS);
		end		% End of switch dataType_MRS
		for ind=indexStart : indexStep : noEntriesListing	% noEntriesListing	% 2  % 1
			% FLAG: CHANGE
			% Preprocess MR spectrum and water
            % with parameters set accoprding to data type (file extension)
			switch fileExt_MRS
				case 'dat'
					% MRS raw data (.dat)
					filename_In			= structFileListing(ind).name;
					% Select file for raw data water signals (.dat) depending on data 
					% type, i.e. how many different signals  (spectra and/or water 
					% signals) are included
					switch dataType_MRS
						case {'mrs_w', 'mrs_w_ref'}
							% Spectra and water signals in list of files/directories
							filename_w_In		= structFileListing(ind+1).name;
						case {'mrs', 'mrs_ref', 'water', 'water_ref'}
							% Only spectra or only water signals in list of files/directories
							filename_w_In		= '';

						otherwise
							error('%s: Unknown MRS dataType_MRS = %s!', sFunctionName, dataType_MRS);
					end		% End of switch dataType_MRS
					%disp(sMsg_newLines);
					fprintf('\n\n');
					disp([sprintf('ind = %d\t', ind), sprintf('\t'), filename_In, sprintf('\t'), filename_w_In, sprintf('\n\n')]);
					%[out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_RawData_s(dirString_In,dirString_Out,filename_In,filename_w_In,seqType_MRS,dataType_MRS,strOVS_In,strOVS_w_In,leftshift_In,noSD_In,structSR_In.aaDomain,structSR_In.tmaxin,structSR_In.iterin,bECC_In,bPhaseCorrFreqShift_In,plotSwitch_In,strMinUserIn_In,reportSwitch_In);
					% [out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_s(...
					% 	dirString_In,...
					% 	dirString_Out,...
					% 	seqType_MRS,...
					% 	dataType_MRS,...
					% 	fileExt_MRS,...
					% 	'Filename', filename_In,...
					% 	'WaterDirectory', dirString_In,...
					% 	'WaterFilename', filename_w_In,...
					% 	'OVS', strOVS_In,...
					% 	'WaterOVS', strOVS_w_In,...
					% 	'Leftshift', leftshift_In,...
					% 	'avgBlockSize', avgBlockSize_In,...
					% 	'RemoveBadAverages', rmbadav_In,...
					% 	'noStandardDeviation', noSD_In,...
					% 	'FreqPhaseCorrectionID', strFreqPhaseCorr_In,...
					% 	'DriftCorrection', driftCorr_In,...
					% 	'Iterations', structSR_In.iterin,...
					% 	'aaDomain', structSR_In.aaDomain,...
					% 	'MaxTimeAlignment', structSR_In.tmaxin,...
					% 	'MaxTimeAlignmentSet', structSR_In.bTmaxset,...
					% 	'medianAlignment', structSR_In.medin,...
					% 	'ppmMinimum_fix', structSR_In.ppmmin_fix,...
					% 	'ppmMaximumArray_fix', structSR_In.ppmmaxarray_fix,...
					% 	'ECC', bECC_In,...
					% 	'PhaseFrequencyCorrection', bPhaseCorrFreqShift_In,...
					% 	'MinimizeUserInput', strMinUserIn_In,...
					% 	'ShowPlots', plotSwitch_In,...						
					% 	'GenerateReport', reportSwitch_In);

					[out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_s(...
						dirString_In,...
						dirString_Out,...
						seqType_MRS,...
						dataType_MRS,...
						fileExt_MRS,...
						structSR_In,...
						structSC_In,...
						'Filename', filename_In,...
						'WaterDirectory', dirString_In,...
						'WaterFilename', filename_w_In,...
						'OVS', strOVS_In,...
						'WaterOVS', strOVS_w_In,...
						'Leftshift', leftshift_In,...
						'avgBlockSize', avgBlockSize_In,...
						'RemoveBadAverages', rmbadav_In,...
						'noStandardDeviation', noSD_In,...
						'FreqPhaseCorrectionID', strFreqPhaseCorr_In,...
						'DriftCorrection', driftCorr_In,...
						'ECC', bECC_In,...
						'PhaseFrequencyCorrection', bPhaseCorrFreqShift_In,...
						'MinimizeUserInput', strMinUserIn_In,...
						'ShowPlots', plotSwitch_In,...
						'GenerateReport', reportSwitch_In);

				case 'IMA'
					% MRS DICOM data (.IMA)
					dirString_In_IMA	= [dirString_In structFileListing(ind).name];
					% Select directory for DICOM water signals (.IMA) depending on data 
					% type, i.e. how many different signals  (spectra and/or water 
					% signals) are included
					switch dataType_MRS
						case {'mrs_w', 'mrs_w_ref'}
							% Spectra and water signals in list of files/directories
							dirString_w_In_IMA	= [dirString_In structFileListing(ind+1).name];
						case {'mrs', 'mrs_ref', 'water', 'water_ref'}
							% Only spectra or only water signals in list of files/directories
							dirString_w_In_IMA	= '';

						otherwise
							error('%s: Unknown MRS dataType_MRS = %s!', sFunctionName, dataType_MRS);
					end		% End of switch dataType_MRS
					% Display some info
					dirParts_In_IMA		= strsplit(dirString_In_IMA, filesep);
					dirParts_w_In_IMA	= strsplit(dirString_w_In_IMA, filesep);
					%disp(sMsg_newLines);
					fprintf('\n\n');
					if ~isempty(dirString_w_In_IMA)
						disp([sprintf('ind = %d\t', ind), sprintf('\t'), dirParts_In_IMA{end-1}, sprintf('\t'), dirParts_w_In_IMA{end-1}, sprintf('\n\n')]);
					else
						disp([sprintf('ind = %d\t', ind), sprintf('\t'), dirParts_In_IMA{end-1}, sprintf('\n\n')]);
					end		% End of if ~isempty(dirString_w_In_IMA)
					
					% [out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_s(...
					% 	dirString_In_IMA,...
					% 	dirString_Out,...
					% 	seqType_MRS,...
					% 	dataType_MRS,...
					% 	fileExt_MRS,...
					% 	'WaterDirectory', dirString_w_In_IMA,...
					% 	'OVS', strOVS_In,...
					% 	'WaterOVS', strOVS_w_In,...
					% 	'Leftshift', leftshift_In,...
					% 	'avgBlockSize', avgBlockSize_In,...
					% 	'RemoveBadAverages', rmbadav_In,...
					% 	'noStandardDeviation', noSD_In,...
					% 	'FreqPhaseCorrectionID', strFreqPhaseCorr_In,...
					% 	'DriftCorrection', driftCorr_In,...
					% 	'Iterations', structSR_In.iterin,...
					% 	'aaDomain', structSR_In.aaDomain,...
					% 	'MaxTimeAlignment', structSR_In.tmaxin,...
					% 	'MaxTimeAlignmentSet', structSR_In.bTmaxset,...
					% 	'medianAlignment', structSR_In.medin,...
					% 	'ppmMinimum_fix', structSR_In.ppmmin_fix,...
					% 	'ppmMaximumArray_fix', structSR_In.ppmmaxarray_fix,...
					% 	'ECC', bECC_In,...
					% 	'PhaseFrequencyCorrection', bPhaseCorrFreqShift_In,...
					% 	'MinimizeUserInput', strMinUserIn_In,...
					% 	'ShowPlots', plotSwitch_In,...						
					%	'GenerateReport', reportSwitch_In);

					[out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_s(...
						dirString_In_IMA,...
						dirString_Out,...
						seqType_MRS,...
						dataType_MRS,...
						fileExt_MRS,...
						structSR_In,...
						structSC_In,...
						'WaterDirectory', dirString_w_In_IMA,...
						'OVS', strOVS_In,...
						'WaterOVS', strOVS_w_In,...
						'Leftshift', leftshift_In,...
						'avgBlockSize', avgBlockSize_In,...
						'RemoveBadAverages', rmbadav_In,...
						'noStandardDeviation', noSD_In,...
						'FreqPhaseCorrectionID', strFreqPhaseCorr_In,...
						'DriftCorrection', driftCorr_In,...
						'ECC', bECC_In,...
						'PhaseFrequencyCorrection', bPhaseCorrFreqShift_In,...
						'MinimizeUserInput', strMinUserIn_In,...
						'ShowPlots', plotSwitch_In,...
						'GenerateReport', reportSwitch_In);
			end			% End of switch fileExt_MRS

			% Close all figures
			%close all;
		end		% End of or ind=indexStart : indexStep : noEntriesListing
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS);
end		% End of switch seqType_MRS


%% Prepare metabolite quantification using LCM analysis, if selected 
% by copying specific groups of processed MRS data files into output directory for LCM
% analysis and write corresponding list of filenames into text files
fprintf('\n\n');
if bPrep_MetabQuant == 1
	fprintf('Preparing metabolite quantification using LCM analysis ...\n\n');
	[status_prep, msg_prep] = prep_MetabQuant_s(dirString_Out, dirString_Out_LCM, seqType_MRS, ...
		'CopyFiles', bCopyFiles_In, 'CopyFiles_MRS', bCopyFiles_MRS_In, ...
		'CopyFiles_ref_Quant', bCopyFiles_ref_Quant_In, 'CopyFiles_ref_ECC', bCopyFiles_ref_ECC_In, ...
		'CopyFiles_w', bCopyFiles_w_In, 'WriteFilenames', bWriteFilenames_In);
	if ~status_prep
		disp(msg_prep);
		error('%s: Preparing metabolite quantification using LCM analysis for study %s, VOI %s, and sequence type %s failed!\n', sFunctionName, strStudy_MRS, strVOI_MRS, seqType_MRS);
	end
else
	fprintf('NO preparation of metabolite quantification using LCM analysis!\n\n');
end		% End of if bPrep_MetabQuant == 1


%% Save variables of workspace to file
% Obtain current date and time in specific format
% Since use of datestr was no longer recommended, code was changed to use datetime instead
%dt		= datestr(now,'yyyymmdd_HH_MM_SS');
%tTmp	= datetime('now', 'Format', 'yyyyMMdd_HH_mm_ss');
%dt		= char(tTmp);
dt		= char(datetime('now', 'Format', 'yyyyMMdd_HH_mm_ss'));

% Save workspace into output directory (optional with user input)
% (Extension".mat" in filename explicitly required, so that Matlab can correctly load 
% workspace file with a "." in its filename)
%strSavedWorkspaceFileName		= 'workspace_run_specialproc_CBF';
strSavedWorkspaceFileName		= ['workspace_', sFunctionName, '_', seqType_MRS, '_', dataType_MRS, '_', dt];
strSavedWorkspaceFileNameFull	= [dirString_Out, strSavedWorkspaceFileName, sprintf('_SD_%.1f.mat', noSD_In)];
%strSaveWorkspace	= input('Would you like to save all variables of the workspace to file?  ', 's');
strSaveWorkspace	= 'y';
if strcmp(strSaveWorkspace,'y') || strcmp(strSaveWorkspace,'Y')
	save(strSavedWorkspaceFileNameFull);
end

