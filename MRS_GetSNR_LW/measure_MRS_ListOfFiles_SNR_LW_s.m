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
filename_MRS_In			= '';
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



%[data_MRS, SNR, FWHM, info]	= measure_MRS_SNR_LW_FIDA_s(dirString, filename_MRS, filename_w, dataFormat_MRS, signal_ppmRange, noise_ppmRange, LWpeak_ppmRange, zp_factor, outDirString, dataType_MRS, bOutFile, plotswitch, seqType_MRS, procParams, Bo_field, spectralWidth, TE, TR);








%% Obtain information about the list of files (.dat) or list of directories (.IMA)
% (assuming that all .dat files are included in the same directory)
% (On Linux, file list in Matlab also includes the two directories "." and "..", which
% means that the actual # of files in the directory is (# of entries in list - 2;
% however, if dir is used to list specific files, e.g. using a file extension, these two 
% directories are not included in the resulting list)
%cd(dirString_In);

% % FLAG: CHANGE
% % Added an option for loading a set of directories containing scans
% switch fileExtension
%     case 'dat'
%         structFileListing		= dir([dirString_In, '*.dat']);
%         noEntriesListing		= length( structFileListing );
%         %noDataFiles				= noEntriesListing - 2
%     case 'IMA'
%         structFileListingAll	= dir(dirString_In);
%         subDir					= [structFileListingAll(:).isdir];
%         structFileListing		= structFileListingAll(subDir);
%         % Remove the two directories '.' and '..'
%         structFileListing		= structFileListing(~ismember({structFileListing(:).name},{'.','..'}));
%         noEntriesListing		= length(structFileListing);
%     otherwise
%         error('%s: ERROR: Unknown file extension %s!', sFunctionName, fileExtension);
% end



%% Preprocess all data files depending on sequence type



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

