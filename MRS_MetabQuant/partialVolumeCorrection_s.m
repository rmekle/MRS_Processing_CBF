%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% partialVolumeCorrection_s.m
%
%% Script to calculate partial volume correction tissue coefficients in MRS for brain
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2018, 2020, 2021, 2022; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'partialVolumeCorrection_s';
sMsg_newLines		= sprintf('\n\n');
sMsg_newLine		= sprintf('\n');
disp(sMsg_newLines);;


%% Init input parameters
command					= '';
status					= 0;
noTissues				= 3;
bCalcPartialVolCoeffs	= 'Yes';			% 'Yes';		% 'No';
winnerFileName			= 'winner.nii';
seqType					= 'MEGA-PRESS';		% 'SPECIAL';	% 'MEGA-PRESS';		% 'sLASER';
strVOI					= 'PCG';			% 'HC';		% 'PCG';
strPVCorr 				= 'PVCorr_bet_87_115_180_fractThresh_0_3';		% '';
strSeg					= 'Trauma_bet_CenterOfBrain_87_115_180_fractThresh_0_3';
strDistCorr				= 'DistCorr';		% 'DistCorr';		% 'ND';

% Set (additional) parameters depending on sequence type 
switch seqType
	case 'SPECIAL'
		outFileName_PVCorr		= '3T_Potsdam_Pain_TissueVolCoeffs.txt';
		dirData_MRS_rda			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_rda_Files_Spectra/';
		dirData_NIfTI			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_NIfTI_Files_DistCorr/';
		dirData_Seg				= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_Segmented/bet_coordCenterOfBrain_mixed_87_115_180_150_fractThresh_0_3_DistCorr/';
		outputDir_PVCorr		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_PartialVolumeCorrection/PVCorr_bet_mixed_87_115_180_150_fractThresh_0_3_DistCorr/';
	case 'MEGA-PRESS'
		outFileName_PVCorr		= '3T_MRS_Dopamin_TissueVolCoeffs.txt';
		%dirData_MRS_rda			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_rda_Files_MRS/';
		dirData_MRS_rda			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_rda_Files_MRS_editOFF/';
		dirData_NIfTI			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Files_DistCorr/';
		dirData_Seg				= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Segmented/DOPA_bet_CenterOfBrain_87_115_180_fractThresh_0_3_DistCorr/';
		outputDir_PVCorr		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_PartialVolumeCorrection/PVCorr_bet_87_115_180_fractThresh_0_3_DistCorr/';
	case 'sLASER'
		outFileName_PVCorr		= ['3T_MRS_Trauma_TissueVolCoeffs_', strVOI, '.txt'];
		dirData_MRS_rda			= ['/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_rda_Files_MRS_', strVOI, filesep];
		dirData_NIfTI 			= ['/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MPRAGE_NIfTI_', strVOI, filesep];
		dirData_Seg				= ['/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MPRAGE_NIfTI_Segmented_', strVOI, filesep, strSeg, filesep];
		outputDir_PVCorr_Base	= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_PartialVolumeCorrection/';
		if ~isempty(strPVCorr)
			% Append name of subdirectory 
			outputDir_PVCorr 	= [outputDir_PVCorr_Base, strPVCorr, '_', strVOI, filesep];
		else
			% Empty string strPVCorr, use base directorey 
			outputDir_PVCorr 	= [outputDir_PVCorr_Base, strVOI, filessep];
		end
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType);
end
fullOutFileName_PVCorr	= fullfile(outputDir_PVCorr, outFileName_PVCorr);


%% Add path to PTB tools used for calculation of partial volume coefficients (e.g. 
% pvoxel3.py) to Linux/Ubuntu environment PATH variable via system call
% DOES NOT WORK YET! Syntax incorrect? NEEDS TO BE DONE in .bashrc to work here
% command				= sprintf('export PATH="$PATH:/data01/MRI_MRS_Software/PTB_Tools/bin/"');
% [status,cmdout]		= system(command);
% if status ~= 0
% 	error('%s: Error when setting PATH variable to include tools for partial volume calculations!\n\n%s', sFunctionName, cmdout);
% end
%
% % If adding path to PTB tolld does not work, test system call that uses absolute path to
% % routine that should be executed
% command = sprintf('/data01/MRI_MRS_Software/PTB_Tools/bin/pvoxel3.py');
% system(command);


%% Obtain information about the list of rda files for MR spectra
% (assuming that all data directories are included in the same directory)
% (On Linux, file list in Matlab also includes the two directories "." and "..", which
% means that the actual # of files in the directory is (# of entries in list - 2;
% however, if dir is used to list specific files, e.g. using a file extension, these two 
% directories are not included in the resulting list)
structFileListing_rda		= dir([dirData_MRS_rda, '*.rda']);
noEntriesListing_rda		= length( structFileListing_rda );
%noDataEntries_rda			= noEntriesListing_rda - 2;

% Obtain information about the list of NIfTI files for partial volume correction
% (assuming that all data files are included in the same directory)
% (On Linux, file list in Matlab also includes the two directories "." and "..", which
% means that the actual # of files in the directory is (# of entries in list - 2)
%disp(sMsg_newLines);
structFileListing_NIfTI		= dir([dirData_NIfTI, '*.nii']);
noEntriesListing_NIfTI		= length( structFileListing_NIfTI );

% Obtain information about the list of files of segmented tissue volumes
% (for tissues WM. GM, and CSF, there should be 3 tissue volumes per case)
structFileListing_Seg		= dir([dirData_Seg, '*_pve_*.gz']);
noEntriesListing_Seg		= length( structFileListing_Seg );


%% Calculate tissue volume coefficients for partial volume correction in MRS, if desired
% Assumptions:	
% - All data files are consecutively sorted, e.g. by date
% - List of rda files, NIFTI files, and segmented tissue volumes are named and sorted in a
% coherent fashion, i.e. the first rda file in the list of rda files corresponds to the
% first NIfTI file in the list of NIfTI files, etc.

% Each calculation of partial tissue volume coefficients requires one rda file for
% obtaining the coordinates of the volume of interest (VOI)/voxel, a NIfTI file of
% T1-weighted images that were used for placement of the VOI and for segmentation, and for
% each segmented tissue a tissue volume file

% Here, for actual calculation the python routine pvoxel3.py is used
% Usage: pvoxel3.py [-d] [-w] [-e T1.nii] [-m maskSource.nii] [-r offsetRead]  [-p offsetPhase]  [-s offsetSlice]  spect.rda compartment1.nii ...
%       -w          winner mode, for each pixel the winning compartment takes it all
%       -d          debug mode: print names, create additional files called highLight...nii, for testing purpose
%       -m          copies maskSource.nii to a file with prefix "mask_" which is 1.0 inside and 0.0 outside of voxel
%       -l          display first highLight file, for testing purpose
%       -e T1.nii   geo info of T1.ni is applied to all compartment.nii
% Output: overlap of spect.rda's voxel with each compartment.nii.
% Options -w and -l exclude each other in the sense that the highlight files to be created
% for displying using option '-l' are not created when using option '-w'
indexStart		= 1;	
indexStep		= 1;	% Optionally adjustable step size
disp(sMsg_newLines);
if(strcmp(bCalcPartialVolCoeffs, 'Yes'))
	fprintf('%s: Calculation of tissue volume coefficients for partial volume correction ...\n', sFunctionName);
	% If output directory(ies) for file output does (do) not exist, create it (them)
	if not(isfolder(outputDir_PVCorr))
		sMsg = sprintf('%s: Creating output directory %s ...\n', sFunctionName, outputDir_PVCorr);
		disp(sMsg);
		if ~mkdir(outputDir_PVCorr)
			error('%s: Could not create (mkdir) output directory %s!\n', sFunctionName, outputDir_PVCorr);
		end
	end
	% If file with resulting tissue volume coefficients already exists, rename it
	if exist(fullOutFileName_PVCorr, 'file')
		% File exists, display warning message and rename file
		% (successful execution of 'movefile' returns status = 1)
		warningMessage	= sprintf('Warning: file %s \n already exists and will be renamed to %s\n\n', fullOutFileName_PVCorr, [outFileName_PVCorr, '_previous']);
		disp(sMsg_newLines);
		disp(warningMessage);
		[status,msg]	= movefile(fullOutFileName_PVCorr, [fullOutFileName_PVCorr, '_previous']);
		if status ~= 1
			error('%s: Error moving/renaming file %s!\n\n%s', sFunctionName, outFileName_PVCorr, msg);
		end
	end		% End of if exist(fullOutFileName_PVCorr, 'file')
	for ind=indexStart : indexStep : noEntriesListing_rda		% noEntriesListing_rda	% 1		% 2
		% Select file from list of rda and from NIfTI files 
		inFileName_rda			= structFileListing_rda(ind).name;
		inFileName_NIfTI		= structFileListing_NIfTI(ind).name;
		
		% Select tissue volume files assuming that numbering of volumes 2, 1, 0 
		% corresponds to tissues WM. GM, and CSF
		indTissues				= noTissues * (ind-1) + 1;
		inFileName_WM			= structFileListing_Seg(indTissues+2).name;
		inFileName_GM			= structFileListing_Seg(indTissues+1).name;
		inFileName_CSF			= structFileListing_Seg(indTissues).name;
		
		% Obtain parts of filename of rda file for renaming some resulting files
		[filepath_rda,name_rda,ext_rda]		= fileparts(inFileName_rda);
		
 		disp(sMsg_newLines);
 		disp([sprintf('ind = %d\t', ind), sprintf('\t'), inFileName_rda, sprintf('\t'), inFileName_NIfTI]);
		disp([sprintf('indTissues = %d\t', indTissues), newline, inFileName_WM, newline, inFileName_GM, newline, inFileName_CSF, sprintf('\n\n')]);

 		% Create command for calculation of of partial tissue volume coefficients using
 		% absolute pathnames for all files and
		% invoke system call for executing calculation
		fullInFileName_rda		= fullfile(dirData_MRS_rda, inFileName_rda);
		fullInFileName_NIfTI	= fullfile(dirData_NIfTI, inFileName_NIfTI);
		fullInFileName_WM		= fullfile(dirData_Seg, inFileName_WM);
		fullInFileName_GM		= fullfile(dirData_Seg, inFileName_GM);
		fullInFileName_CSF		= fullfile(dirData_Seg, inFileName_CSF);
		
		% If routines to be executed via system call cannot be found in Matlab, check
		% .bashrc file, where the path to these routines should be added to PATH
		% environment variable of system
		% NOTE that using absolute paths of routines to be executed via system call does
		% NOT work, since routine pvoxel3.py calls routine rpnND, which is not found if
		% PATH variable of system environment doe not include the path to all these tools
% 		command				= sprintf('/data01/MRI_MRS_Software/PTB_Tools/bin/pvoxel3.py -w -d -e %s %s %s %s %s | /data01/MRI_MRS_Software/PTB_Tools/bin/transpose.py | tee -a %s', ...
% 			fullInFileName_NIfTI, fullInFileName_rda, fullInFileName_WM, fullInFileName_GM, fullInFileName_CSF, fullOutFileName_PVCorr);
		command				= sprintf('pvoxel3.py -w -d -e %s %s %s %s %s | transpose.py | tee -a %s', ...
			fullInFileName_NIfTI, fullInFileName_rda, fullInFileName_WM, fullInFileName_GM, fullInFileName_CSF, fullOutFileName_PVCorr);

		[status,cmdout]		= system(command);
		if status ~= 0
			error('%s: Error in calculation of tissue volume coefficients for case %s!\n\n%s', sFunctionName, inFileName_rda, cmdout);
		end
		
		% When using pvoxel3.py with the option '-w', for each case, a file "winner.nii" 
		% is created in the directory, where MATLAB is executed in;
		% Rename "winner.nii" to be case-specific and move renamed file into output
		% directory for tissue volume coefficients
		% (successful execution of 'movefile' returns status = 1)
		[status,msg]	= movefile( winnerFileName, fullfile(outputDir_PVCorr, [name_rda, '_', winnerFileName]) );
		if status ~= 1
			error('%s: Error moving file %s!\n\n%s', sFunctionName, winnerFileName, msg);
		end
	end		% End of for ind=indexStart : indexStep : noDataEntries_rda
else
	fprintf('%s: No calculation of tissue volume coefficients for partial volume correction!\n', sFunctionName);
end		% End of if(strcmp(bCalcPartialVolCoeffs, 'Yes'))


%% Save variables of workspace to file
% Obtain current date and time in specific format
dt		= datestr(now,'yyyymmdd_HH_MM_SS');

% Save workspace into output directory (optional with user input)
% (Extension".mat" in filename explicitly required, so that Matlab can correctly load 
% workspace file with a "." in its filename)
strSavedWorkspaceFileName		= ['workspace_', sFunctionName, '_', dt];
strSavedWorkspaceFileNameFull	= [outputDir_PVCorr, strSavedWorkspaceFileName];
%strSaveWorkspace	= input('Would you like to save all variables of the workspace to file?  ', 's');
strSaveWorkspace	= 'y';
if strcmp(strSaveWorkspace,'y') || strcmp(strSaveWorkspace,'Y')
	disp(sMsg_newLines);
	fprintf('%s: Saving variables of workspace to file ...\n', sFunctionName);
	save(strSavedWorkspaceFileNameFull);
end

