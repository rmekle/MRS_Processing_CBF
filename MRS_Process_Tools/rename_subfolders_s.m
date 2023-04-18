%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% rename_subfolders_s.m
%
%% Script to rename subfolders of several directories
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2023; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'rename_subfolders_s';
sMsg_newLines		= sprintf('\n\n');
sMsg_newLine		= newline;
disp(sMsg_newLines);


%% Init input parameters for renaming subfolders
parentDir_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/';
parentDir_AddOn		= 'MRS_Trauma_00_All_DICOM_Test/';	% 'MRS_Trauma_00_All_DICOM/';
parentDir			= [parentDir_Base, parentDir_AddOn];

% Select start and end paatern for substring extraction
startPat		= 'MRS_Trauma_';
endPat			= 'DICOM';


%% Obtain information about subfolders in parent directory and select extraction pattern for substring
% Check whether parent directory exists
if ~isfolder(parentDir)
	error('%s: Parent directory %s is not a folder!\n', sFunctionName, parentDir);
end

% List subfolders in parent directory and only include folders (directories) and exclude
% Linux directories '.' and '..'
list_subfolders			= dir(parentDir);
list_subfolders			= list_subfolders([list_subfolders.isdir]);
list_subfolders_name	= {list_subfolders.name};
list_subfolders_name(strncmp(list_subfolders_name, '.', 1)) = [];

% Loop over all subfolders
% For each subfolder, extract desired substring (e.g. pseudonym for subject) from name of
% subfolder assuming that such information is included
% Then, add extracted substring to all folder names within each subfolder of parent
% directory
for iFolder = 1:1:numel(list_subfolders_name)
	% Extract desired substring from name of subfolder
	newChr		= extractBetween(list_subfolders_name{iFolder}, startPat, endPat);


end