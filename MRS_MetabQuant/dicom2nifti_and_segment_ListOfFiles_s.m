%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% dicom2nifti_and_segment_ListOfFiles_s.m
%
%% Script to convert a list of DICOM files of into NIfTI format and segment them
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2018, 2019, 2021; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'dicom2nifti_and_segment_ListOfFiles_s';
sMsg_newLines		= sprintf('\n\n');
sMsg_newLine		= sprintf('\n');
disp(sMsg_newLines);


%% Init input parameters
inputDir				= '';
command					= '';
status					= 0;
bProcessNewFiles		= 0;
bConvert_dcm2nii		= 'Yes';		% 'Yes';		% 'No';
bSegmentImages			= 'No';			% 'Yes';		% 'No';

% Set (additional) parameters
% % 3T Potsdam_Pain study
% dirData_DICOM			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_DICOM_Files/';
% outputDir_NIfTI			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_NIfTI_Files/';
% dirData_NIfTI			= outputDir_NIfTI;
% outputDir_Seg			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_Segmented/';

% % 3T BCAN MRS_and_Dopamin study
% %dirData_DICOM			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_DICOM_Files/';
% dirData_DICOM			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_DICOM_Files_New/';
% %outputDir_NIfTI			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Files/';
% outputDir_NIfTI			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Files_New/';
% dirData_NIfTI			= outputDir_NIfTI;
% %outputDir_Seg			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Segmented/';
% outputDir_Seg			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Segmented_New/';

% 3T BCAN MRS_and_Trauma study
dirData_DICOM			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MPRAGE_DICOM';
outputDir_NIfTI			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MPRAGE_NIfTI';
outputDir_Seg			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MPRAGE_NIfTI_Segmented';
% Adjust directory names, if only new (newly acquired) data should be processed
if bProcessNewFiles
	dirData_DICOM		= [dirData_DICOM, '_New'];
	outputDir_NIfTI		= [outputDir_NIfTI, '_New'];
	outputDir_Seg		= [outputDir_Seg, '_New'];
end
dirData_DICOM			= [dirData_DICOM, filesep];
outputDir_NIfTI			= [outputDir_NIfTI, filesep];
dirData_NIfTI			= outputDir_NIfTI;
outputDir_Seg			= [outputDir_Seg, filesep];


% Set parameters for brain extraction and segmentation
% FSL bet routine: options
%	-c <x y z>  centre-of-gravity (voxels not mm) of initial mesh surface
%	-f <f>      fractional intensity threshold (0->1); default=0.5; smaller values give larger brain outline estimates
%	-R          robust brain centre estimation (iterates BET several times)
% [85 107 175];		[89 113 142];	[87 110 160];	[87 115 180];	[87 115 150];	[87 115 170];
coordCenterOfBrain		= [87 115 180];	
fractIntensThresh		= 0.3;			% 0.5;		0.4;	0.3;	0.25;	0.2;

% FSL fast routine: options
% 	-t,--type       type of image 1=T1, 2=T2, 3=PD; default=T1
% 	-n,--class      number of tissue-type classes; default=3
%	-H,--Hyper      segmentation spatial smoothness; default=0.1
%	-I,--iter       number of main-loop iterations during bias-field removal; default=4
%	-l,--lowpass    bias field smoothing extent (FWHM) in mm; default=20
%	-o,--out        output basename


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
disp(sMsg_newLines);
if(strcmp(bConvert_dcm2nii, 'Yes'))
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
end		% End of if(strcmp(bConvert_dcm2nii, 'Yes'))


%% Segment NIfTI brain imaging data files into tissues WM, GM, and CSF. if desired
%% Obtain information about the list of NIfTI files to be segmented
% (assuming that all data files are included in the same directory)
% (On Linux, file list in Matlab also includes the two directories "." and "..", which
% means that the actual # of files in the directory is (# of entries in list - 2)
%cd(dirData_In);
disp(sMsg_newLines);
structFileListing_NIfTI		= dir([dirData_NIfTI, '*.nii']);
noEntriesListing_NIfTI		= length( structFileListing_NIfTI );
%noDataEntries_NIfTI			= noEntriesListing_NIfTI - 2

% Segment each NIfTI brain imaging dataset into tissues WM, GM, and CSF, if desired
% using routines from FSL as system calls 
indexStart		= 1;	% To possibly skip entries for specific files
indexStep		= 1;	% Optionally adjustable step size
disp(sMsg_newLines);
if(strcmp(bSegmentImages, 'Yes'))
	disp('Segmenting imaging datasets into tissues WM, GM, and CSF ...');
	for ind=indexStart : indexStep : noEntriesListing_NIfTI		% noEntriesListing_NIfTI	% 1		% 2
		% Select file from list of NIfTI files and obtain parts of filename
		inputFileNameSeg				= structFileListing_NIfTI(ind).name;
		[filepathSeg,nameSeg,extSeg]	= fileparts(inputFileNameSeg);
		disp(sMsg_newLines);
		disp([sprintf('ind = %d\t', ind), sprintf('\t'), inputFileNameSeg, sprintf('\n\n')]);
		
		% Use FSL brain extraction routine 'bet' to extract images of only brain tissue
		% via a system call
		% For MPRAGE data that includes substantial parts of the neck, the bet routine
		% yields improved results, when supplied with the approximate coordinates of the 
		% center of the brain in voxels (check e.g. with 'fsleyes' for that);
		% same coordinates can be generally used for a set of imaging cases
		betFileName			= [nameSeg, '_bet'];
		betFileNameFull		= [betFileName, extSeg, '.gz'];
		command				= sprintf('bet %s %s -v -R -c %d %d %d -f %.1f', ...
			fullfile(dirData_NIfTI, inputFileNameSeg), fullfile(outputDir_Seg, betFileName), coordCenterOfBrain(1), coordCenterOfBrain(2), coordCenterOfBrain(3), fractIntensThresh);		
			%fullfile(dirData_NIfTI, inputFileNameSeg), betFileName, coordCenterOfBrain(1), coordCenterOfBrain(2), coordCenterOfBrain(3), fractIntensThresh);
		[status,cmdout]		= system(command);
		if status ~= 0
			error('%s: Error in brain tissue extraction for file %s!\n\n%s', sFunctionName, inputFileNameSeg, cmdout);
		end
		
% 		% Move results from brain extraction (_bet) into output directory for segmentation
% 		command			= sprintf('mv %s %s', fullfile(dirData_NIfTI, betFileNameFull), fullfile(outputDir_Seg, betFileNameFull));
% 		[status,cmdout]	= system(command);
% 		if status ~= 0
% 			error('%s: Error moving results from brain extraction for file %s!\n\n%s', sFunctionName, betFileNameFull, cmdout);
% 		end
		
		% Use FSL routine 'fast' to segment T1-w MPRAGE images, from which all brain
		% tissue was extracted, into tissues WM, GM, and CSF
		% via s system call
		% fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -o /home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_Segmented/3T_20170510_PetraO_MPRAGE_1_0MM_SAG_0012_bet /home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_Segmented/3T_20170510_PetraO_MPRAGE_1_0MM_SAG_0012_bet.nii.gz
		command			= sprintf('fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -o %s %s', fullfile(outputDir_Seg, betFileName), fullfile(outputDir_Seg, betFileNameFull));
		[status,cmdout]	= system(command);
		if status ~= 0
			error('%s: Error segmenting file %s!\n\n%s', sFunctionName, betFileNameFull, cmdout);
		end		
		
	end		% End of for ind=indexStart : indexStep : noEntriesListing_NIfTI
else
	fprintf('%s: No segmentation of images!\n', sFunctionName);
end		% End of if(strcmp(bSegmentImages, 'Yes'))


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
	disp(sMsg_newLines);
	fprintf('%s: Saving variables of workspace to file ...\n', sFunctionName);
	save(strSavedWorkspaceFileNameFull);
end

