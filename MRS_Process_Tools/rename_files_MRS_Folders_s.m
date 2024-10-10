%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% rename_files_MRS_Folders_s.m
%
%% Script to rename files distributed in several directories following a specific pattern
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
					% 'SBAM/';		% SBAM_Test/';
parentDir			= [parentDir_Base, parentDir_AddOn];

% Select start and end pattern for substring extraction, if required
% End pattern can be empty, if substring to be extracted is at end of name of subdirectory
% Select string additions that have to be added at start and/or end of new filenames
startPat		= 'MRS_Trauma_';
endPat			= '';		% '_DICOM';
strAddStart		= '3T_';
strAddEnd		= '_';
lenStartPat		= length(startPat);

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
newNames	= '';
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
		% Add specific string to all filenames and foldernames of subfolder
		% Create string to be added by extracting information from the name of each
		% subfolder and concatenate this with selected additional patterns
		% Find selected patterns in name of subfolder and extract string between these
		% patterns
		strTmp				= list_subfolders_name{iFolder};
		lenSubfolderName	= length(strTmp);
		kStart				= strfind(strTmp, startPat);
		if isempty(kStart)
			% Starting pattern was not found in name of subfolder
			error('%s: Error! Starting pattern startPat = %s not found in name of subfolder = %s!\n\n%s', sFunctionName, startPat, list_subfolders_name{iFolder});
		end		% End of if isempty(kStart)
		kEnd				= strfind(strTmp, endPat);
		if isempty(kEnd)
			% End pattern was not found in name of subfolder
			% Set to length of name of subfolder plus one, so that subsequent code yields
			% end index equal to end of name of subfolder
			kEnd	= lenSubfolderName + 1;			
		end		% End of if isempty(kEnd)
		% Determine indices of start and end of string to be extracted from first
		% occurrences of start and end pattern in name of subfolder
		indStart		= kStart(1) + lenStartPat;
		indEnd			= kEnd(1) - 1;		
		strAddInfo		= strTmp(indStart:indEnd);

		% Create specific string  and add it to all filenames and foldernames of subfolder
		% Strcat(...) works with all elements of cell array
		strAddSame		= [strAddStart, strAddInfo, strAddEnd];
		newNames		= strcat(strAddSame	, list_subfolders_sub_name);
		% Rename all files and folders in subfolder
		for iFile = 1 : 1 : numel(list_subfolders_sub_name)
			[status,msg]	= movefile( fullfile(subDir, list_subfolders_sub_name{iFile}), fullfile(subDir, newNames{iFile}) );
			if status ~= 1
				error('%s: Error moving/renaming file %s!\n\n%s', sFunctionName, list_subfolders_sub_name{iFile}, msg);
			end		% End of if status ~= 1
		end		% End of for iFile = 1 : 1 : numel(list_subfolders_sub_name)
	else
		error('%s: Error! Renaming files and folders individually not yet implemented!\n\n', sFunctionName);
	end		% End of if bSameRenaming == 1

end		% End of for iFolder = 1:1:numel(list_subfolders_name)
