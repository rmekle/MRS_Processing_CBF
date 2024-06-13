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


%% Init input parameters for easuring SNR and Linewidth (LW) in list of files of MRS data
dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/Z_Test_Data/Test_MRS_SNR_LW/';
%dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Trauma_Analysis/PCG_dat_FID-A_SD_3_2_ECCref_ls3_SR1/PCG_LCModel_Data_MRS_only/';
filename_MRS_In			= '3T_SBA_C_0008_20210325_meas_MID00333_FID95325_svs_slaser_dkd_HC_TE23_WS256_wOVS_3.2_processed_lcm.RAW';
filename_w_In			= '';
dataFormat_MRS_In		= 'lcmRAW';
signal_ppmRange_In		= [1.8, 2.2];
noise_ppmRange_In		= [-2.0, 0];
LWpeak_ppmRange_In		= [2.85, 3.15];
zp_factor_In			= 8;
outDirString_AddOn_1	= '';
outDirString_In			= [dirString_In, outDirString_AddOn_1];
dataType_MRS_In			= 'mrs';
bOutFile_In				= 0;
plotswitch_In			= 0;
seqType_MRS_In			= 'sLASER';
procParams_In			= struct([]);
Bo_field_In				= [];
spectralWidth_In		= [];
TE_In					= [];
TR_In					= [];


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


%% Measure SNR and LW for all MRS data files

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
	filename_In			= structFileListing(ind).name;
	fprintf('\n\n');
	if indexStep == 2
		filename_w_In		= structFileListing(ind+1).name;
		fprintf('ind = %d\t\t%s\t%s\n\n', ind, filename_In, filename_w_In);
	else	% No water file
		fprintf('ind = %d\t\t%s\n\n', ind, filename_In);
	end
	[data_MRS, SNR, FWHM, info]	= measure_MRS_SNR_LW_FIDA_s(dirString_In, filename_MRS_In, filename_w_In, dataFormat_MRS_In, signal_ppmRange_In, noise_ppmRange_In, LWpeak_ppmRange_In, zp_factor_In, outDirString_In, dataType_MRS_In, bOutFile_In, plotswitch_In, seqType_MRS_In, procParams_In, Bo_field_In, spectralWidth_In, TE_In, TR_In);
end		% End of or ind=indexStart : indexStep : noEntriesListing



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

