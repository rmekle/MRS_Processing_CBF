%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% preprocess_ListOfFiles_s.m
%
%% Script to preprocess a list of files of magnetic resonance spectroscopy (MRS) data
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2018, 2019, 2020, 2021, 2022,
% 2023; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'preprocess_ListOfFiles_s';
sMsg_newLines		= sprintf('\n\n');
sMsg_newLine		= sprintf('\n');
disp(sMsg_newLines);


%% Init input parameters for preprocessing
%dirString_In			= '';
%dirString_Out			= '';
fileExtension           = 'IMA';		% Currently: 'dat' (raw data) or 'IMA' (DICOM)
filename_In				= '';
filename_w_In			= '';
strStudy				= '3T_MMs';		% '3T_Trauma';	'7T_KCL';	'3T_MMs';
strVOI					= 'PCG';			% 'PCG';	% 'HC'; % 'Pons'; % 'CB'; % 'PFC'; % 'PCC';
seqType_MRS				= 'sLASER';		% 'SPECIAL';	% 'MEGA-PRESS'; % 'sLASER';
dataType_MRS			= 'mrs_ref';		% 'mrs_w_ref';		'mrs_w';	% 'mrs_ref';	
signals_MRS				= 'MMs';		% 'MMs';	% 'Spectra';
strOVS_In				= 'wOVS';		% 'wOVS';	% 'woutOVS';
strOVS_w_In				= 'wOVS';		% 'wOVS';	% 'woutOVS';
leftshift_In			= 1;		% 3;	% 2;	% 0;	% 1;
noSD_In					= 3.2;			% 3.2;		2.6;		4.0;
digits					= [fix(noSD_In) round(abs(noSD_In-fix(noSD_In))*10)];

% Parameters for spectral registration (aligning of averages/frequency and phase drift
% correction) performed in either frequency or time domain
driftCorr_In			= 'y';		% 'y';		'n';
iterin_In				= 20;
aaDomain_In				= 'f';		% 'f';		't';
tmaxin_In				= 0.2;		% 0.2;		0.1;
bTmaxset_In				= 1;
medin_In				= 'y';		% 'y';	'n';	'a';	'ref';
% Set parameters for drift correction depending on type of data, i.e. whether MRS
% data is spectrum or water signal
% NOTE: Check whether aligning of averages in frequency domain works, if the MR
% spectrum is water signal itself; if not, simply align averages in time domain
switch dataType_MRS
	case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
		% MR spectrum is provided together without or with unsuppressed water
		% signal and/or with reference scans
		ppmmin_fix_In		= 1.6;		% 1.6;		1.8;
		%ppmmaxarray_fix_In	= [3.5; 4.0; 5.5];
		ppmmaxarray_fix_In	= [2.4,2.85,3.35,4.2,4.4,5.2];
	case {'water', 'water_ref'}
		% MR spectrum is water signal itself without or with reference scans
		ppmmin_fix_In		= 4.2;
		ppmmaxarray_fix_In	= [5.5 5.5 5.2];

	otherwise
		error('%s: Unknown MRS dataType_MRS = %s!', sFunctionName, dataType_MRS);
end		% End of switch dataType_MRS
alignSS_In				= 2;		% For aligning subspectra (e.g. in SPECIAL)
strSpecReg				= 'SR1';	% To distinguish settings for spectral registration

% Additional parameter settings
bECC_In					= 1;
bPhaseCorrFreqShift_In	= 0;
strMinUserIn_In			= 'y';
plotSwitch_In			= 0;
reportSwitch_In			= 1;


%% Set (additional) parameters depending on sequence type
switch seqType_MRS
	case 'SPECIAL'
		dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_RawData_dat_Files/';
		dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Z_Pain_Tmp/';
		%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Preprocessed_forLCModel_SD_4_0/';
		%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Preprocessed_forLCModel_SD_3_2/';
		%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Preprocessed_forLCModel_SD_2_6/';
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
	case 'sLASER'
		% Select data input and output directories depending on study, MRS data type, 
		% i.e. file extension, study, and other parameters
		%digits = [fix(noSD_In) round(abs(noSD_In-fix(noSD_In))*10)];
		switch strStudy
			case '3T_Trauma'
				% Data (input) directories
				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/';
				%dirString_In_AddOn1		= sprintf('MRS_Trauma_00_All_RawData_dat_Files_MRS_%s', strVOI);
				
				% Output data directory
				dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Trauma/Trauma_Z_Analysis/';
				
				% Directories depending on MRS data type
				switch fileExtension
					case 'dat'
						% Select directories specific to MRS raw data (.dat)
						dirString_In_AddOn1		= sprintf('MRS_Trauma_00_All_RawData_dat_Files_MRS_%s', strVOI);
					case 'IMA'
						% Select directories specific to MRS DICOM data (.IMA)
						dirString_In_AddOn1		= sprintf('MRS_Trauma_00_All_DICOM_IMA_Files_MRS_%s', strVOI);
						
					otherwise
						error('%s: ERROR: Unknown file extension (data type) %s!', sFunctionName, fileExtension);
				end			% End of switch fileExtension
			case '3T_MMs'
				% Data (input) directories
				%dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/';
				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/destiana/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MMs/';
				%dirString_In_Base		= '/home/destiana/CSB_NeuroRad/destiana/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MMs/';
				
				
				% Output data directory
				%dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Trauma/Trauma_Z_Analysis/';
				dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/ZZZZ_Test/';
				%dirString_Out_Base		= '/home/destiana/CSB_NeuroRad/destiana/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MMs/CodeResults/';

				% Directories depending on MRS data type
				switch fileExtension
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
						error('%s: ERROR: Unknown file extension (data type) %s!', sFunctionName, fileExtension);
				end			% End of switch fileExtension
			case '7T_KCL'
				% Data (input) directories
				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/7T_KCL/';
				%dirString_In_AddOn1		= sprintf('7T_KCL_00_ALL_RawData_dat_Files_MRS_eja_%s', strVOI);
				
				% Output data directory
				dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/7T_KCL_Analysis/';
				
				% Directories depending on MRS data type
				switch fileExtension
					case 'dat'
						% Select directories specific to MRS raw data (.dat)
						dirString_In_AddOn1		= sprintf('7T_KCL_00_ALL_RawData_dat_Files_MRS_eja_%s', strVOI);
					case 'IMA'
						% Select directories specific to MRS DICOM data (.IMA)
						dirString_In_AddOn1		= sprintf('7T_KCL_00_ALL_DICOM_IMA_Files_MRS_eja_%s', strVOI);
						
					otherwise
						error('%s: ERROR: Unknown file extension (data type) %s!', sFunctionName, fileExtension);
				end			% End of switch fileExtension
				
			otherwise
				error('%s: ERROR: Unknown study %s!', sFunctionName, strStudy);
		end				% End of switch strStudy
		
		% Complete names of data (input) directories
		dirString_In			= [dirString_In_Base, dirString_In_AddOn1, filesep];
			
		% Select directory for output data depending on voxel location, data type,
		% # of SDs, and other options used for pre-processing of MR spectra
		%digits = [fix(noSD_In) round(abs(noSD_In-fix(noSD_In))*10)];
		%dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Trauma/Trauma_Z_Analysis/';
		%dirString_Out_AddOn1	= sprintf('%s_FID-A_SD_%d_%d', strVOI, digits(1), digits(2));
		% Make output directories for acquired macromolecules (MMs) distinguishable from 
		% those for spectra
		if strcmp(signals_MRS, 'MMs')
			dirString_Out_AddOn1	= sprintf('%s_%s_%s_FID-A_SD_%d_%d', signals_MRS, strVOI, fileExtension, digits(1), digits(2));
		else
			dirString_Out_AddOn1	= sprintf('%s_%s_FID-A_SD_%d_%d', strVOI, fileExtension, digits(1), digits(2));
		end		% End of if strcmp(signals_MRS, 'MMs')
		dirString_Out_AddOn2	= '';
		if bECC_In
			% Use reference (water) signals for ECC, if acquired
			% If not, then use an unsuppressed water signal, if acquired
			% If no reference and no water signals are acquired, check whether MR spectrum
			% is water signal itself; and if it is, use it for ECC
			% Indicate different options for ECC in the corresponding directory name; for
			% that, search for strings 'ref', 'w', and 'water' in string for MRS data type
			refInd		= strfind(dataType_MRS, '_ref');
			wInd		= strfind(dataType_MRS, '_w');
			waterInd	= strfind(dataType_MRS, 'water');
			if ~isempty(refInd)
				dirString_Out_AddOn2	= '_ECCref';
			else
				if ~isempty(wInd)
					dirString_Out_AddOn2	= '_ECCw';
				else
					if ~isempty(waterInd)
						dirString_Out_AddOn2	= '_ECCwater';
					else
						% No reference and no water signals and MR spectrum is not water
						% signal itself => ECC not possible
						error('%s: No reference and no water signals and MR spectrum is not water signal itself (dataType_MRS = %s) => ECC not possible!', sFunctionName, dataType_MRS);
					end		% End of if ~isempty(waterInd)
				end		% End of if ~isempty(wInd)
			end		% if ~isempty(refInd)
		end		% End of if bECC_In
		if leftshift_In > 0
			dirString_Out_AddOn2	= [dirString_Out_AddOn2, sprintf('_ls%d', leftshift_In)];
		end		% End of if leftshift_In > 0
		% Indicate in output directory name which type of spectral registration was used
		if driftCorr_In == 'y' || driftCorr_In == 'Y'
			dirString_Out_AddOn2	= [dirString_Out_AddOn2, '_', strSpecReg];
		else
			dirString_Out_AddOn2	= [dirString_Out_AddOn2, '_NoSR'];
		end
		dirString_Out			= [dirString_Out_Base, dirString_Out_AddOn1, dirString_Out_AddOn2, filesep];
		%dirString_Out			= [dirString_Out_Base, dirString_Out_AddOn1, dirString_Out_AddOn2, '_Test', filesep];
		
		% If output directory is non-existent, create it
		if ~exist( dirString_Out, 'dir' )
			mkdir(dirString_Out);
		else
			% Output directory already exists, ask user whether to overwrite or not
			% (should help to avoid accidentally overwriting previosuly processed data)
			strOverwrite	= input('\n\nDo you want to overwrite the existing output directory (y/n)?  ','s');
			if strOverwrite == 'n' || strOverwrite  == 'N'
				fprintf('\n%s: Already existing output directory is not overwritten! Processing is aborted!\n\n', sFunctionName)
				return;
			else
				fprintf('\n%s: Already existing output directory is overwritten! Processing is continued!\n\n', sFunctionName)
			end		% End of if strOverwrite == 'n' || strOverwrite  == 'N'
		end		% End of if ~exist( dirString_Out, 'dir' )
		
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
switch fileExtension
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
        error('%s: ERROR: Unknown file extension %s!', sFunctionName, fileExtension);
end



%% Preprocess all data files depending on sequence type

% FLAG: TODO: Update for all 3 cases the function to the new version with
% name-value pair parameters, depending on if imaDataSwitch is set or not

switch seqType_MRS
	case 'SPECIAL'
		% Pre-processing of .IMA files is not implemented for run_specialproc_CBF, only
        % for the sLASER sequence using preProcess_MRS_s
        if strcmp(fileExtension, 'IMA')
            error('%s: ERROR: File extension IMA incompatible with sequence type %s!', sFunctionName, seqType_MRS);
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
			disp(sMsg_newLines);
			disp([sprintf('ind = %d\t', ind), strOVS_In, sprintf('\t'), filename_In, sprintf('\t'), filename_w_In, sprintf('\n\n')]);
			[out,out_w,out_noproc,out_w_noproc]=run_specialproc_CBF(dirString_In,dirString_Out,filename_In,filename_w_In,noSD_In,strOVS_In,strMinUserIn_In,aaDomain_In,tmaxin_In,iterin_In);
			
			% Close all figures
			close all;
			
			% Preprocess MR spectrum and water_woutOVS
			% (spectral data remain the same and only processed differently now)
			strOVS_In		= 'woutOVS';
			filename_w_In	= structFileListing(ind+2).name;
			disp(sMsg_newLines);
			disp([sprintf('ind = %d\t', ind), strOVS_In, sprintf('\t'), filename_In, sprintf('\t'), filename_w_In, sprintf('\n\n')]);
			[out,out_w,out_noproc,out_w_noproc]=run_specialproc_CBF(dirString_In,dirString_Out,filename_In,filename_w_In,noSD_In,strOVS_In,strMinUserIn_In,aaDomain_In,tmaxin_In,iterin_In);
		end
	case 'MEGA-PRESS'
		% Pre-processing of .IMA files is not implemented for run_megapressproc_CBF, only
        % for the sLASER sequence using preProcess_MRS_s
        if strcmp(fileExtension, 'IMA')
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
			disp(sMsg_newLines);
			disp([sprintf('ind = %d\t', ind), sprintf('\t'), filename_In, sprintf('\t'), filename_w_In, sprintf('\n\n')]);
			[diffSpecOut,sumSpecOut,subSpec1Out,subSpec2Out,outwOut,outw_subSpec1Out,outw_subSpec2Out,coilcombosOut]=run_megapressproc_CBF(dirString_In,dirString_Out,filename_In,filename_w_In,noSD_In,strMinUserIn_In,aaDomain_In,tmaxin_In,iterin_In,alignSS_In);
			
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
		indexStart		= 1;
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
			switch fileExtension
				case 'dat'
					% MRS raw data (.dat)
					filename_In			= structFileListing(ind).name;
					filename_w_In		= structFileListing(ind+1).name;
					disp(sMsg_newLines);
					disp([sprintf('ind = %d\t', ind), sprintf('\t'), filename_In, sprintf('\t'), filename_w_In, sprintf('\n\n')]);
					%[out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_RawData_s(dirString_In,dirString_Out,filename_In,filename_w_In,seqType_MRS,dataType_MRS,strOVS_In,strOVS_w_In,leftshift_In,noSD_In,aaDomain_In,tmaxin_In,iterin_In,bECC_In,bPhaseCorrFreqShift_In,plotSwitch_In,strMinUserIn_In,reportSwitch_In);
					[out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_s(...
						dirString_In,...
						dirString_Out,...
						seqType_MRS,...
						dataType_MRS,...
						'Filename', filename_In,...
						'WaterDirectory', dirString_In,...
						'WaterFilename', filename_w_In,...
						'OVS', strOVS_In,...
						'WaterOVS', strOVS_w_In,...
						'Leftshift', leftshift_In,...
						'noStandardDeviation', noSD_In,...
						'DriftCorrection', driftCorr_In,...
						'Iterations', iterin_In,...
						'aaDomain', aaDomain_In,...
						'MaxTimeAlignment', tmaxin_In,...
						'MaxTimeAlignmentSet', bTmaxset_In,...
						'medianAlignment', medin_In,...
						'ppmMinimum_fix', ppmmin_fix_In,...
						'ppmMaximumArray_fix', ppmmaxarray_fix_In,...
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
					disp(sMsg_newLines);
					if ~isempty(dirString_w_In_IMA)
						disp([sprintf('ind = %d\t', ind), sprintf('\t'), dirParts_In_IMA{end-1}, sprintf('\t'), dirParts_w_In_IMA{end-1}, sprintf('\n\n')]);
					else
						disp([sprintf('ind = %d\t', ind), sprintf('\t'), dirParts_In_IMA{end-1}, sprintf('\n\n')]);
					end		% End of if ~isempty(dirString_w_In_IMA)
					
					[out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_s(...
						dirString_In_IMA,...
						dirString_Out,...
						seqType_MRS,...
						dataType_MRS,...
						'WaterDirectory', dirString_w_In_IMA,...
						'OVS', strOVS_In,...
						'WaterOVS', strOVS_w_In,...
						'Leftshift', leftshift_In,...
						'noStandardDeviation', noSD_In,...
						'DriftCorrection', driftCorr_In,...
						'Iterations', iterin_In,...
						'aaDomain', aaDomain_In,...
						'MaxTimeAlignment', tmaxin_In,...
						'MaxTimeAlignmentSet', bTmaxset_In,...
						'medianAlignment', medin_In,...
						'ppmMinimum_fix', ppmmin_fix_In,...
						'ppmMaximumArray_fix', ppmmaxarray_fix_In,...
						'ECC', bECC_In,...
						'PhaseFrequencyCorrection', bPhaseCorrFreqShift_In,...
						'MinimizeUserInput', strMinUserIn_In,...
						'ShowPlots', plotSwitch_In,...						
						'GenerateReport', reportSwitch_In);
			end			% End of switch fileExtension
			
			% Close all figures
			%close all;
		end		% End of or ind=indexStart : indexStep : noEntriesListing
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS);
end		% End of switch seqType_MRS


%% Save variables of workspace to file
% Obtain current date and time in specific format
dt		= datestr(now,'yyyymmdd_HH_MM_SS');

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

