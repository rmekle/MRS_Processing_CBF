%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% prepare_MetabQuant_s.m
%
%% Script to prepare metabolite quantification of magnetic resonance spectroscopy (MRS) data
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2020; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= ' prepare_MetabQuant_s';
sMsg_newLines		= sprintf('\n\n');
%disp(sMsg_newLines);


%% Init input parameters for preprocessing routine
%dirString_In			= '';
%outDirString_In			= '';
filename_In				= '';
filenamew_In			= '';
noSD_In					= 4.0;			% 2.6;		3.2;		4.0;
%strOVS_In				='wOVS';
%strMinUserIn_In			= 'y';
%aaDomain_In				= 'f';
seqType					= 'MEGA-PRESS';		% 'SPECIAL';	% 'MEGA-PRESS';
bCopyFiles				= 1;
bWriteFilenames			= 1;

% Set (additional) parameters depending on sequence type
switch seqType
	case 'SPECIAL'
		dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_RawData_dat_Files/';
		outDirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Z_Pain_Tmp/';
	case 'MEGA-PRESS'
		% Select directory for (input) data files depending on # of SDs used for
		% pre-processing of MR spectra
		switch(noSD_In)
			case(2.6)
				dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_2_6/';
			case(3.2)
				%dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2/';
				dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2_New/';
			case(4.0)
				dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_4_0/';
			
			otherwise
				error('%s: ERROR: No directory/data for noSD_In =  %f!', sFunctionName, noSD_In);
		end
		outDirString			= [dirString_In, 'DOPA_LCModel_Analysis_Data/'];
		%outDirString_w			= [outDirString, 'Water_Signals/'];
		textFileName_diff_MRS	= 'list_filenames_MRS_Diff_Spectra.txt';
		textFileName_OFF		= 'list_filenames_MRS_editOFF_All.txt';
		textFileName_OFF_MRS	= 'list_filenames_MRS_editOFF_Spectra.txt';
		textFileName_OFF_w	= 'list_filenames_MRS_editOFF_Water.txt';
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType);
end

% If output directories for file copying do not exist, create them
if not(isfolder(outDirString))
	sMsg = sprintf('%s: Creating output directory %s ...\n', sFunctionName, outDirString);
    disp(sMsg);
    if ~mkdir(outDirString)
		error('%s: Could not create (mkdir) output directory %s!\n', sFunctionName, outDirString);
	end
end
% % Additional output directory for water signals
% if not(isfolder(outDirString_w))
% 	sMsg = sprintf('%s: Creating output directory %s ...\n', sFunctionName, outDirString_w);
%     disp(sMsg);
%     if ~mkdir(outDirString_w)
% 		error('%s: Could not create (mkdir) output directory %s!\n', sFunctionName, outDirString_w);
% 	end
% end


%% Obtain information about the list of different groups of files
% (assuming that all data files are included in the same directory)
% (On Linux, file list in Matlab also includes the two directories "." and "..", which
% means that the actual # of files in the directory is (# of entries in list - 2; 
% however, if dir is used to list specific files, e.g. using a file extension, these two 
% directories are not included in the resulting list)
%cd(dirString_In);
% MEGA-PRESS difference files (MRS + water)
structFileListing_diff		= dir([dirString_In, '*_diff_*.RAW']);
noEntriesListing_diff		= length( structFileListing_diff );
%noDataFiles_diff			= noEntriesListing_diff - 2

% MEGA-PRESS edit_OFF files (MRS + water)
structFileListing_OFF		= dir([dirString_In, '*_editOFF_*.RAW']);
noEntriesListing_OFF		= length( structFileListing_OFF );

% MEGA-PRESS water difference files (water only)
structFileListing_diff_w	= dir([dirString_In, '*_w_*_diff_*.RAW']);
noEntriesListing_diff_w		= length( structFileListing_diff_w );

% MEGA-PRESS water edit_OFF files (water only)
structFileListing_OFF_w		= dir([dirString_In, '*_w_*_editOFF_*.RAW']);
noEntriesListing_OFF_w		= length( structFileListing_OFF_w );

% Determine MEGA-PRESS MRS difference files (MRS only)
structFileListing_diff_MRS	= structFileListing_diff(~ismember({structFileListing_diff.name}, {structFileListing_diff_w.name}));
noEntriesListing_diff_MRS	= length( structFileListing_diff_MRS );

% Determine MEGA-PRESS MRS edit_OFF files (MRS only)
structFileListing_OFF_MRS	= structFileListing_OFF(~ismember({structFileListing_OFF.name}, {structFileListing_OFF_w.name}));
noEntriesListing_OFF_MRS	= length( structFileListing_OFF_MRS );


%% Copy data files depending on sequence type and create list of filenames, if selected
switch seqType
	case 'SPECIAL'
		% Nothing yet
	case 'MEGA-PRESS'
		% Copy only MEGA-PRESS MRS difference files and all edit_OFF files with .RAW 
		% extension to destination folder, if selected
		indexStart		= 1;
		indexStep		= 1;
		if bCopyFiles
			% Copy MEGA-PRESS MRS difference files
			for ind=indexStart : indexStep : noEntriesListing_diff_MRS
				% Copy MEGA-PRESS difference file
				[status, msg] = copyfile(fullfile(dirString_In, structFileListing_diff_MRS(ind).name), outDirString);
				if ~status
					disp(msg);
					error('%s: Could not copy file %s to output directory %s!\n', sFunctionName, structFileListing_diff_MRS(ind).name, outDirString);
				end
			end
			
			% Copy MEGA-PRESS edit_OFF files (MRS + water)
			for ind=indexStart : indexStep : noEntriesListing_OFF
				% Copy MEGA-PRESS edit_OFF file
				[status, msg] = copyfile(fullfile(dirString_In, structFileListing_OFF(ind).name), outDirString);
				if ~status
					disp(msg);
					error('%s: Could not copy file %s to output directory %s!\n', sFunctionName, structFileListing_OFF(ind).name, outDirString);
				end
			end
			
		end		% End of if bCopyFiles
			
		
		%% Write filenames for specific groups of data files into corresponding textfiles
		% (additional for loops are used to allow for greater flexibility), if selected
		istart		= 1;
		istep		= 1;
		
		if bWriteFilenames
			% Filenames of all MEGA-PRESS edit_OFF files (MRS + water)
			fid_OFF			= fopen(fullfile(outDirString, textFileName_OFF), 'wt+');
			if fid_OFF == -1
				error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF, outDirString);
			end
			for i=istart : istep : (noEntriesListing_OFF-1)
				nbytes	= fprintf(fid_OFF, '%s\n', structFileListing_OFF(i).name);
			end
			% Write last filename without new line character at end of line to textfile
			nbytes	= fprintf(fid_OFF, '%s', structFileListing_OFF(noEntriesListing_OFF).name);
			if fclose(fid_OFF) == -1
				error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF, outDirString);
			end
			
			% Filenames of MEGA-PRESS water edit_OFF files (water only)
			fid_OFF_w		= fopen(fullfile(outDirString, textFileName_OFF_w), 'wt+');
			if fid_OFF_w == -1
				error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF_w, outDirString);
			end
			for i=istart : istep : (noEntriesListing_OFF_w-1)
				nbytes	= fprintf(fid_OFF_w, '%s\n', structFileListing_OFF_w(i).name);
			end
			% Write last filename without new line character at end of line to textfile
			nbytes	= fprintf(fid_OFF_w, '%s', structFileListing_OFF_w(noEntriesListing_OFF_w).name);
			if fclose(fid_OFF_w) == -1
				error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF_w, outDirString);
			end
			
			% Filenames of MEGA-PRESS MRS edit_OFF files (MRS only)
			fid_OFF_MRS		= fopen(fullfile(outDirString, textFileName_OFF_MRS), 'wt+');
			if fid_OFF_MRS == -1
				error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF_MRS, outDirString);
			end
			for i=istart : istep : (noEntriesListing_OFF_MRS-1)
				nbytes	= fprintf(fid_OFF_MRS, '%s\n', structFileListing_OFF_MRS(i).name);
			end
			% Write last filename without new line character at end of line to textfile
			nbytes	= fprintf(fid_OFF_MRS, '%s', structFileListing_OFF_MRS(noEntriesListing_OFF_MRS).name);
			if fclose(fid_OFF_MRS) == -1
				error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF_MRS, outDirString);
			end
			
			% Filenames of MEGA-PRESS MRS difference files (MRS only)
			fid_diff_MRS		= fopen(fullfile(outDirString, textFileName_diff_MRS), 'wt+');
			if fid_diff_MRS == -1
				error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_diff_MRS, outDirString);
			end
			for i=istart : istep : (noEntriesListing_diff_MRS-1)
				nbytes	= fprintf(fid_diff_MRS, '%s\n', structFileListing_diff_MRS(i).name);
			end
			% Write last filename without new line character at end of line to textfile
			nbytes	= fprintf(fid_diff_MRS, '%s', structFileListing_diff_MRS(noEntriesListing_diff_MRS).name);
			if fclose(fid_diff_MRS) == -1
				error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_diff_MRS, outDirString);
			end
			
		end		% End of if bWriteFilenames
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType);
end


%% Save variables of workspace to file
% Obtain current date and time in specific format
dt		= datestr(now,'yyyymmdd_HH_MM_SS');

% Save workspace into output directory (optional with user input)
% (Extension".mat" in filename explicitly required, so that Matlab can correctly load 
% workspace file with a "." in its filename)
%strSavedWorkspaceFileName		= 'workspace_run_specialproc_CBF';
strSavedWorkspaceFileName		= ['workspace_', sFunctionName, '_', seqType, '_', dt];
strSavedWorkspaceFileNameFull	= [outDirString, strSavedWorkspaceFileName, sprintf('_SD_%.1f.mat', noSD_In)];
%strSaveWorkspace	= input('Would you like to save all variables of the workspace to file?  ', 's');
strSaveWorkspace	= 'n';
if strcmp(strSaveWorkspace,'y') || strcmp(strSaveWorkspace,'Y')
	save(strSavedWorkspaceFileNameFull);
end

