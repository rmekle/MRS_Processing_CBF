%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% measure_MRS_ListOfFiles_SNR_LW_s.m
%
%% Script to measure SNR and Linewidth (LW) in list of files of MR spectroscopy (MRS) data
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;



%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'measure_MRS_ListOfFiles_SNR_LW_s';
% sMsg_newLines		= sprintf('\n\n');
% sMsg_newLine		= sprintf('\n');
% disp(sMsg_newLines);
fprintf('\n\n');


%% Init parameters for measuring SNR and Linewidth (LW) in list of files of MRS data
% Parameters to select sequence, study, volume-of-interest (VOI)/voxel, and file extension
% of original MRS data
seqType_MRS_In			= 'sLASER';
strStudy				= '3T_Trauma';  	% 'Test'; '3T_Trauma'; 
strVOI					= 'PCG'; 	% 'PCG'; 'HC'; 'Pons'; 'CB'; 'PFC'; 'PCC';
fileExtension           = 'dat';		% Currently: 'dat' (raw data) or 'IMA' (DICOM)

% Parameters for saving of results to file
% Select range (here first field) in Excel to write cell array of results to
bSaveResults			= 1;
outNamingOption			= 1;
%outputFileName_Add_1	= '_SNR_FWHM';
acOutFileType			= '.xlsx';		% '.xlsx';	'.txt';
strRangeSel				= 'A4';		% 'A4';

% Input directory options
% Use base directory and directory AddOns (e.g. subfolder names) to allow flexible choice
% of output filename, if results are saved to file
switch seqType_MRS_In
	case 'sLASER'
		switch strStudy
			case 'Test'
				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/';
				dirString_In_AddOn_1	= 'Z_Test_Data';
				dirString_In_AddOn_2	= 'Test_MRS_SNR_LW';
				outNamingOption			= 2;
			case '3T_Trauma'
				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Trauma_Analysis/';
				switch fileExtension
					case 'dat'
						dirString_In_AddOn_1	= [strVOI, '_', fileExtension, '_FID-A_SD_3_2_ECCref_ls3_SR1'];
					case 'IMA'
						dirString_In_AddOn_1	= [strVOI, '_', fileExtension, '_FID-A_SD_3_2_ECCref_ls1_SR1'];

					otherwise
						error('%s: ERROR: Unknown file extension (data type) %s!', sFunctionName, fileExtension);
				end			% End of switch fileExtension
				dirString_In_AddOn_2	= [strVOI, '_LCModel_Data_MRS_only'];
				outNamingOption			= 1;

			otherwise
				error('%s: ERROR: Unknown study %s!', sFunctionName, strStudy);
		end				% End of switch strStudy

	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS_In);
end		% End of switch seqType_MRS_In
% Complete name of input directory
dirString_In			= [dirString_In_Base, dirString_In_AddOn_1, filesep, dirString_In_AddOn_2, filesep];

% Additional parameters
filename_MRS_In			= '';
filename_w_In			= '';
dataFormat_MRS_In		= 'lcmRAW';
signal_ppmRange_In		= [1.8, 2.2];
noise_ppmRange_In		= [-3.0, -1.0];
LWpeak_ppmRange_In		= [2.9, 3.1];
zp_factor_In			= 8;
dataType_MRS_In			= 'mrs';
bAutoPhase_In			= 0;
bOutFile_In				= 0;
plotswitch_In			= 0;
procParams_In			= struct([]);
Bo_field_In				= [];
spectralWidth_In		= [];
TE_In					= [];
TR_In					= [];

% Output directory and filenmae options
outDirString_AddOn_1	= '';
outDirString_In			= [dirString_In, outDirString_AddOn_1];
outputFileName_Add_1	= sprintf('_SNR_%.1f_%.1f_%.1f_%.1f_FWHM_%.2f_%.2f', ...
	signal_ppmRange_In(1), signal_ppmRange_In(2), noise_ppmRange_In(1), noise_ppmRange_In(2), ...
	LWpeak_ppmRange_In(1), LWpeak_ppmRange_In(2));
outputFileName_Add_2	= '';
if bAutoPhase_In
	outputFileName_Add_2	= '_phasedZO';
end


%% Obtain information about the list of files for (preprocessed) MR spectra
% (assuming that all files are included in the same directory)
% (On Linux, file list in Matlab also includes the two directories "." and "..", which
% means that the actual # of files in the directory is (# of entries in list - 2;
% however, if dir is used to list specific files, e.g. using a file extension, these two 
% directories are not included in the resulting list)
%cd(dirString_In);


% Select file extension to search for depending on format of MRS data
% (again, here assuming that all MRS data files are in same directory)
% Then determine filenames and # of files for all corresponding MRS data files
switch dataFormat_MRS_In
    case 'dat'
		acSearchString			= '*.dat';
        %structFileListing		= dir([dirString_In, '*.dat']);
        %noEntriesListing		= length( structFileListing );
        %noDataFiles				= noEntriesListing - 2
	case 'DICOM'
		acSearchString			= '*.dcm';
    case 'IMA'
		acSearchString			= '*.IMA';
        %structFileListingAll	= dir(dirString_In);
        %subDir					= [structFileListingAll(:).isdir];
        %structFileListing		= structFileListingAll(subDir);
        % Remove the two directories '.' and '..'
        %structFileListing		= structFileListing(~ismember({structFileListing(:).name},{'.','..'}));
        %noEntriesListing		= length(structFileListing);
	case 'rda'
		acSearchString			= '*.rda';
	case 'lcmRAW'
		acSearchString			= '*.RAW';

    otherwise
        error('%s: ERROR: Unknown dataFormat_MRS_In %s!', sFunctionName, dataFormat_MRS_In);
end
structFileListing		= dir([dirString_In, acSearchString]);
noEntriesListing		= length( structFileListing );


%% Measure SNR and LW (FWHM) for all MRS data files
% Allocate arrays to store MRS data, measurement values, and zero order phases applied
% prior to measurements
data_MRS		= cell(noEntriesListing, 1);
SNR				= zeros(noEntriesListing, 1);
FWHM			= zeros(noEntriesListing, 1);
phase0			= zeros(noEntriesListing, 1);

% Measure desired quantities for each case (spectrum)
% Select size for stepping through indices, i.e. list of files 
% depending on data type, i.e. how many different signals
% (spectra and/or water signals) are included
indexStart		= 1;
indexStep		= 1;
switch dataType_MRS_In
	case {'mrs_w', 'mrs_w_ref'}
		% Spectra and water signals in list of files/directories
		indexStep		= 2;
	case {'mrs', 'mrs_ref', 'water', 'water_ref'}
		% Only spectra or only water signals in list of files/directories
		indexStep		= 1;

	otherwise
		error('%s: Unknown MRS dataType_MRS_In = %s!', sFunctionName, dataType_MRS_In);
end		% End of switch dataType_MRS_In
for ind=indexStart : indexStep : noEntriesListing	% noEntriesListing	% 2  % 1
	filename_MRS_In			= structFileListing(ind).name;
	fprintf('\n\n');
	if indexStep == 2
		filename_w_In		= structFileListing(ind+1).name;
		fprintf('ind = %d\t\t%s\t%s\n\n', ind, filename_MRS_In, filename_w_In);
	else	% No water file
		fprintf('ind = %d\t\t%s\n\n', ind, filename_MRS_In);
	end
	[data_MRS{ind}, SNR(ind), FWHM(ind), phase0(ind), info]	= measure_MRS_SNR_LW_FIDA_s(dirString_In, filename_MRS_In, filename_w_In, dataFormat_MRS_In, signal_ppmRange_In, noise_ppmRange_In, LWpeak_ppmRange_In, zp_factor_In, outDirString_In, dataType_MRS_In, bAutoPhase_In, bOutFile_In, plotswitch_In, seqType_MRS_In, procParams_In, Bo_field_In, spectralWidth_In, TE_In, TR_In);
end		% End of or ind=indexStart : indexStep : noEntriesListing
fprintf('\n\n');


%% Include info about data files and results into one cell array
% Create cell arrays with info line (header), all MRS data filenames, and combine them
% with results into new cell array
% Dimensions of cell arrays have to match for that				% N = noEntriesListing
cellInfoLine		= {'MRS Data File' 'SNR' 'FWHM / Hz' 'Phase0 LW_Peak / deg'};	% Yields 1x4 cell array
cellDataFileNames	= {structFileListing(:).name}';						% Yields Nx1 cell array	
cellData			= [cellDataFileNames num2cell([SNR FWHM phase0])];	% Yields Nx4 cell array
cellInfoAndData		= [cellInfoLine; cellData];							% Yields (N+1)x4 cell array


%% Save results from SNR and LW measurements to file, if selected
if bSaveResults
	% Create output filename for results depending on selected naming option:
	%		Name of (input) subfolder chosen to best describe the MRS data
	%		User input
	switch outNamingOption
		case 1
			outFileName_Base	= dirString_In_AddOn_1;
		case 2
			outFileName_Base	= dirString_In_AddOn_2;
		case 9
			outFileName_Base	= input('\n\nPlease enter the base of the output filename: ', "s");

		otherwise
			error('%s: Unknown outNamingOption = %d!', sFunctionName, outNamingOption);
	end		% End of switch namingOption
	outFileName		= [outFileName_Base, outputFileName_Add_1, outputFileName_Add_2, acOutFileType];

	% Write cell array with info and results from SNR and LW measurements to file 
	% File type depends on chosen file extension: .xls is speradsheet and .txt is textfile
	fprintf('Saving results for SNR and LW measurenents to file ...\n\n');
	writecell( cellInfoAndData, fullfile(outDirString_In, outFileName), ...
		'WriteMode', 'inplace', 'AutoFitWidth', 1, 'Range', strRangeSel)
end		% End of if bSaveResults


%% Save variables of workspace to file
% Obtain current date and time in specific format
dt		= datestr(now,'yyyymmdd_HH_MM_SS');

% Save workspace into output directory (optional with user input)
% (Extension".mat" in filename explicitly required, so that Matlab can correctly load 
% workspace file with a "." in its filename)
strSavedWorkspaceFileName		= ['workspace_', sFunctionName, '_', seqType_MRS_In, '_', dataType_MRS_In, '_', dt];
strSavedWorkspaceFileNameFull	= [outDirString_In, strSavedWorkspaceFileName, '.mat'];
%strSaveWorkspace	= input('Would you like to save all variables of the workspace to file?  ', 's');
strSaveWorkspace	= 'y';
if strcmp(strSaveWorkspace,'y') || strcmp(strSaveWorkspace,'Y')
	save(strSavedWorkspaceFileNameFull);
end

