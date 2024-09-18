%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% rename_files_MRS_Folders_s.m
%
%% Script to rename files distributed in several directories
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2024; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'rename_files_MRS_Folders_s';
% fprintf('\n\n');


%% Init input parameters for renaming subfolders
parentDir_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/';
parentDir_AddOn		= 'SBAM_Test/';	
					% 'SBAM/';
parentDir			= [parentDir_Base, parentDir_AddOn];

% Select start and end pattern for substring extraction, if required
% End pattern can be empty, if substring to be extracted is at end of name of subdirectory
% Select string additions that have to be added at start and/or end of new filenames
startPat		= 'MRS_Trauma_';
endPat			= '';		% '_DICOM';
aCh_AddStart	= '3T_';
aCh_AddEnd		= '';

% Select whether all files and folders in each subfolder are renamed using the same
% modification or whether individual modifcations are used
bSameRenaming	= 1;


%% Obtain information about subfolders in parent directory and select extraction pattern for substring
% Check whether parent directory exists
if ~isfolder(parentDir)
	error('%s: Parent directory %s is not a folder!\n', sFunctionName, parentDir);
end

% List subfolders in parent directory and only include folders (directories) and exclude
% Linux directories '.' and '..'
list_subfolders			= dir_s(parentDir);
list_subfolders			= list_subfolders([list_subfolders.isdir]);
list_subfolders_name	= {list_subfolders.name};
%list_subfolders_name(strncmp(list_subfolders_name, '.', 1)) = [];

% Loop over all subfolders
% For each subfolder, rename all files in subfolder with info composed of acquisition
% parameters, subject desgination, and date of data acquisition
%newChr		= cell(1);
newName		= '';
for iFolder = 1:1:numel(list_subfolders_name)
% 	% Extract desired substring from name of subfolder
% 	newChr		= extractBetween(list_subfolders_name{iFolder}, startPat, endPat);
% 	if isempty(newChr)
% 		error('%s: Error! Extracted substring newChr = %s is empty!\n', sFunctionName, newChr{1,1})
% 	end

	% Determine list of directories and files in current subdirectory
	subDir						= fullfile(parentDir, list_subfolders_name{iFolder});
	list_subfolders_sub			= dir_s(subDir);
	list_subfolders_sub_name	= {list_subfolders_sub.name};

	% Create new name for each file and subfolder
	% Check whether renaming of all files and folders is the same
	if bSameRenaming == 1

	% Assume that only one MRS DICOM file is in each subdirectory
	newName						= [list_subfolders_name{iFolder}, '.IMA'];
	% [status,msg]	= movefile( fullfile(subDir, list_subfolders_sub_name{1}), fullfile(subDir, newName) );
	% if status ~= 1
	% 	error('%s: Error moving/renaming file %s!\n\n%s', sFunctionName, list_subfolders_sub_name{iSub}, msg);
	% end
	else
		error('%s: Error! Renaming files and folders individually not yet implemented!\n\n', sFunctionName);
	end		% End of if bSameRenaming == 1

end		% End of for iFolder = 1:1:numel(list_subfolders_name)
