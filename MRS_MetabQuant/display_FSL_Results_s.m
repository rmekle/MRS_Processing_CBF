%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% display_FSL_Results_s.m
%
%% Script to display results from FSL brain extraction (bet) and segmentation (segment)
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2018, 2021; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'display_FSL_Results_s';
sMsg_newLines		= sprintf('\n\n');
sMsg_newLine		= sprintf('\n');
disp(sMsg_newLines);;


%% Init input parameters
%inputDirImages			= '';
command					= '';
status					= 0;
bShowOverlays			= 'Yes';			% 'Yes';		% 'No';
seqType					= 'sLASER';		% 'SPECIAL';	% 'MEGA-PRESS';		% 'sLASER';

% Set (additional) parameters depending on sequence type
switch seqType
	case 'SPECIAL'
		inputDirImages			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_NIfTI_Files/';
		%inputDirOverlays		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_Segmented/bet_coordCenterOfBrain_87_110_160_fractThresh_0_25/';
		%inputDirOverlays		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_Segmented/bet_coordCenterOfBrain_87_110_160_fractThresh_0_3/';
		inputDirOverlays		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_Segmented/bet_coordCenterOfBrain_87_115_150_fractThresh_0_3/';
		%inputDirOverlays		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_Segmented/bet_coordCenterOfBrain_87_115_180_fractThresh_0_3/';
	case 'MEGA-PRESS'
		%inputDirImages			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Files/';
		%inputDirImages			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Files_Copy/';
		inputDirImages			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Files_DistCorr/';
		%inputDirOverlays		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Segmented/DOPA_bet_CenterOfBrain_87_115_180_fractThresh_0_3/';
		inputDirOverlays		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Segmented/DOPA_bet_CenterOfBrain_87_115_180_fractThresh_0_3_DistCorr/';
		%inputDirOverlays		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Segmented/DOPA_bet_CenterOfBrain_87_115_180_fractThresh_0_3_ND/';
		%inputDirOverlays		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Segmented/DOPA_bet_CenterOfBrain_87_115_170_fractThresh_0_3/';
		%inputDirOverlays		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Segmented/DOPA_bet_CenterOfBrain_87_115_175_fractThresh_0_3/';
	case 'sLASER'
		inputDirImages 		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MPRAGE_NIfTI/';
		inputDirOverlays	= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MPRAGE_NIfTI_Segmented/Trauma_bet_CenterOfBrain_87_115_180_fractThresh_0_3/';
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType);
end
outputDir				= inputDirOverlays;


%% Obtain information about the list of image (e.g. NIfTI) files and overlays
% (assuming that all image files are included in the same directory)
% Directory of images
structFileListing_Images	= dir([inputDirImages, '*.nii']);
noEntriesListing_Images		= length( structFileListing_Images );

% Directory of overlays
structFileListing_Overlays	= dir([inputDirOverlays, '*_bet.nii.gz']);
noEntriesListing_Overlays	= length( structFileListing_Overlays );

% Check whether the # of files in these two file lists is the same
if( noEntriesListing_Images ~= noEntriesListing_Overlays )
	error('%s: # of image = %d and overlay files = %d are not the same!\n\n%s', sFunctionName, noEntriesListing_Images, noEntriesListing_Overlays);
end

% If '*_ND_bet.nii.gz' and '*_bet.nii.gz' overlays (e.g. from MPRAGE and MPRAGE_ND images)
% exist, the default order from "dir" is that '*_ND_bet.nii.gz' comes before
% '*_bet.nii.gz' (probably due to capital versus small letters), which is opposite to the
% order of the corresponding images; To amend this, the corresponding filenames are sorted
% in UPPER case, and the sorted indices used to reorder entries in the original struct of
% overlays
namesOverlays				= {structFileListing_Overlays(:).name};
[sortedOverlays, sortInd]	= sort(upper(namesOverlays));
structFileListing_Overlays	= structFileListing_Overlays(sortInd);


%% Display overlays (from brain extraction) onto (NIfTI) source images, if desired
% Assumptions:	
% - All image files are consecutively sorted, e.g. by date
% - Each image file has one corresponding overlay file, and the sorting order of these two
%	list of files has to match

% Display overlays (from brain extraction) onto (NIfTI) source images, if desired 
% using a system call that invokes the fsl utility "fsleyes"
indexStart		= 1;	
indexStep		= 4;	% Optionally adjustable step size
disp(sMsg_newLines);
if(strcmp(bShowOverlays, 'Yes'))
	fprintf('%s: Displaying overlays onto source images ...\n', sFunctionName);
	for ind=indexStart : indexStep : noEntriesListing_Images	% noEntriesListing_Images	% 2		% 1
		% Select pair of files from list of image and overlay files
		inputFileNameImage					= structFileListing_Images(ind).name;
		%[filepathImage,nameImage,extImage]	= fileparts(inputFileNameImage);
		inputFileNameOverlay					= structFileListing_Overlays(ind).name;
		disp(sMsg_newLines);
		disp([sprintf('ind = %d\t', ind), sprintf('\t'), inputFileNameImage, newline, sprintf('ind = %d\t', ind), sprintf('\t'), inputFileNameOverlay, sprintf('\n\n')]);
 		
		% Create command for display of overlays onto source images using "fsleyes" and
		% invoke system call for display
		% Set colormaps for display of both types of data
 		command				= sprintf( 'fsleyes %s --cmap greyscale %s --cmap brain_colours_5redyell', ...
			fullfile(inputDirImages, inputFileNameImage), fullfile(inputDirOverlays, inputFileNameOverlay) );
		[status,cmdout]		= system(command);
		if status ~= 0
			error('%s: Error for display of image %s and overlay %s!\n\n%s', sFunctionName, inputFileNameImage, inputFileNameOverlay, cmdout);
		end
	end		% End of for ind=indexStart : indexStep : noDataEntries_Images
else
	fprintf('%s: No display of overlays onto source images!\n', sFunctionName);
end		% End of if(strcmp(bShowOverlays, 'Yes'))


%% Save variables of workspace to file
% Save workspace into output directory (optional with user input)
% (Extension".mat" in filename explicitly required, so that Matlab can correctly load 
% workspace file with a "." in its filename)
strSavedWorkspaceFileName		= ['workspace_', sFunctionName];
strSavedWorkspaceFileNameFull	= [outputDir, strSavedWorkspaceFileName];
%strSaveWorkspace	= input('Would you like to save all variables of the workspace to file?  ', 's');
strSaveWorkspace	= 'n';
if strcmp(strSaveWorkspace,'y') || strcmp(strSaveWorkspace,'Y')
	disp(sMsg_newLines);
	fprintf('%s: Saving variables of workspace to file ...\n', sFunctionName);
	save(strSavedWorkspaceFileNameFull);
end

