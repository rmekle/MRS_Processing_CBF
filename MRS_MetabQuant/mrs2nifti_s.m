%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mrs2nifti_s.m
%
%% Script to convert a list of MRS data files of into corresponding NIfTI format
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2025; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'mrs2nifti_s';
% sMsg_newLines		= sprintf('\n\n');
% sMsg_newLine		= sprintf('\n');
% disp(sMsg_newLines);
fprintf('\n\n');


%% Init input parameters
inputDir				= '';
command					= '';
status					= 0;
bProcessNewFiles		= 0;
bConvert_mrs2nii		= 1;			% 1;		% 0;
seqType_MRS				= 'sLASER';		% 'SPECIAL';	% 'MEGA-PRESS'; % 'sLASER';
strStudy				= '3T_ENIGMA';	% '3T_Trauma';	'7T_KCL';	'3T_MMs';	'3T_SBAM';

% Set (additional) parameters depending on sequence type
switch seqType_MRS
	case 'SPECIAL'
		% % 3T CBF Potsdam_Pain study
		dirData_MRS				= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_DICOM_Files/';
		outputDir_NIfTI			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_NIfTI_Files/';
	case 'MEGA-PRESS'
		% 3T BCAN MRS_and_Dopamin study
		dirData_DICOM			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_DICOM_Files/';
		outputDir_NIfTI			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Files/';
	case 'sLASER'
		% Select data input and output directories and other parameters depending on study
		switch strStudy
			case '3T_Trauma'
				% 3T BCAN MRS Trauma study
				% SBA
				dirData_DICOM			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MPRAGE_DICOM/';
				outputDir_NIfTI			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MPRAGE_NIfTI/';
			case '3T_SBAM'
				% 3T BCAN MRS Trauma study
				% SBAM
				dirData_DICOM			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/SBAM/MRS_SBAM_00_All_MPRAGE_DICOM/';
				outputDir_NIfTI			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/SBAM/MRS_SBAM_00_All_MPRAGE_NIfTI/';

			otherwise
				error('%s: ERROR: Unknown study %s!', sFunctionName, strStudy);
		end				% End of switch strStudy
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS);
end

% Adjust directory names, if only new (newly acquired) data should be processed
if bProcessNewFiles
	dirData_DICOM		= [dirData_DICOM, '_New'];
	outputDir_NIfTI		= [outputDir_NIfTI, '_New'];
end
dirData_DICOM			= [dirData_DICOM, filesep];
outputDir_NIfTI			= [outputDir_NIfTI, filesep];

% If output directories for specific processing opitions do not exist, create them
if bConvert_mrs2nii
	if not(isfolder(outputDir_NIfTI))
		sMsg = sprintf('%s: Creating output directory %s ...\n', sFunctionName, outputDir_NIfTI);
		disp(sMsg);
		if ~mkdir(outputDir_NIfTI)
			error('%s: Could not create (mkdir) output directory %s!\n', sFunctionName, outputDir_NIfTI);
		end
	end
end		% End of if bConvert_dcm2ni


%% Obtain information about the list of DICOM files or list of directories, respectively
% (assuming that all data directories are included in the same directory)
% (On Linux, file list in Matlab also includes the two directories "." and "..", which
% means that the actual # of files in the directory is (# of entries in list - 2)
%cd(dirData_In);
structFileListing_DICOM		= dir(dirData_DICOM);
noEntriesListing_DICOM		= length( structFileListing_DICOM );
noDataEntries_DICOM			= noEntriesListing_DICOM - 2;


%% Convert all (MPRAGE) DICOM data files into NIfTI format, if desired
% Assumptions:	
% - All data files are consecutively sorted, e.g. by date
% - Each set of DICOM images is in one directory of the list of directories just obtained

% Convert DICOM images of each directory into NIfTI format, if desired 
% using a system call that invokes the utility "dcm2niix"
% and place all output NIfTI files into same directory
% (Note that here the counter for the for loop has to include all entries up to the last
% index!)
indexStart		= 3;	% To skip entries for directories "." and ".."
indexStep		= 1;	% Optionally adjustable step size
%disp(sMsg_newLines);
fprintf('\n\n');
if bConvert_dcm2nii 
	fprintf('%s: Conversion of DICOM data into NIfTI format ...\n', sFunctionName);
	for ind=indexStart : indexStep : noEntriesListing_DICOM		% noEntriesListing_DICOM	% 3		% 4
		subDirData_DICOM	= structFileListing_DICOM(ind).name;
		inputDir			= fullfile(dirData_DICOM, subDirData_DICOM, filesep);
		disp(sMsg_newLines);
		disp([sprintf('ind = %d\t', ind), sprintf('\t'), subDirData_DICOM, sprintf('\n\n')]);
		
		% Create command for conversion to NIfTI for each set of DICOM input images and
		% invoke system call for NIfTI conversion
		% (last argument on command line for dcm2niix is the input directory)
		command				= sprintf('dcm2niix -z n -f %s -v 1 -o %s %s', subDirData_DICOM, outputDir_NIfTI, inputDir);
		[status,cmdout]		= system(command);
		if status ~= 0
			error('%s: Error in DICOM to NIfTI conversion for data in %s!\n\n%s', sFunctionName, subDirData_DICOM, cmdout);
		end
	end		% End of for ind=indexStart : indexStep : noDataEntries_DICOM
else
	fprintf('%s: No conversion of DICOM data into NIfTI format!\n', sFunctionName);
end		% End of if bConvert_dcm2nii 


%% Save variables of workspace to file


% Obtain current date and time in specific format
dt		= datestr(now,'yyyymmdd_HH_MM_SS');

% Save workspace into desired output directory (optional with user input)
% (Extension".mat" in filename explicitly required, so that Matlab can correctly load 
% workspace file with a "." in its filename)
strSavedWorkspaceFileName		= ['workspace_', sFunctionName, '_', dt];
% Select output directory depending on processing options; if segmentation was included,
% save workspace to output directory for segmentation; if
%strSaveWorkspace	= input('Would you like to save all variables of the workspace to file?  ', 's');
strSaveWorkspace				= 'y';
strSavedWorkspaceFileNameFull	= [outputDir_NIfTI, strSavedWorkspaceFileName];
if(strcmp(bSegmentImages, 'Yes'))
	strSavedWorkspaceFileNameFull	= [outputDir_Seg, strSavedWorkspaceFileName];
else if(~strcmp(bConvert_dcm2nii, 'Yes'))
		% Do not save workspace, since no processing was done
		strSaveWorkspace	= 'n';
	end
end
if strcmp(strSaveWorkspace,'y') || strcmp(strSaveWorkspace,'Y')
	%disp(sMsg_newLines);
	fprintf('\n\n');
	fprintf('%s: Saving variables of workspace to file ...\n', sFunctionName);
	save(strSavedWorkspaceFileNameFull);
end

