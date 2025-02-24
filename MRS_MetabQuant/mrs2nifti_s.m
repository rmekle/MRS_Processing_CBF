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
fprintf('\n\n');


%% Init input parameters
inputDir				= '';
command					= '';
status					= 0;
bProcessNewFiles		= 0;
bConvert_mrs2nii		= 1;			% 1;		% 0;
strStudy_MRS			= 'ENIGMA_3T_SBA';	% '3T_Trauma';	'7T_KCL';	'3T_MMs';	'3T_SBAM';
fileExt_MRS				= 'IMA';		% Currently: 'dat' (raw data) or 'IMA' (DICOM) or '.dcm' (enhanced DICOM)
seqType_MRS				= 'sLASER';		% 'SPECIAL';	% 'MEGA-PRESS'; % 'sLASER';
dataType_MRS			= 'mrs_ref';		% 'mrs_w_ref';		'mrs_w';	% 'mrs_ref';

% Set (additional) parameters depending on sequence type
switch seqType_MRS
	case 'SPECIAL'
		% % 3T CBF Potsdam_Pain study
		dirData_MRS				= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_DICOM_Files/';
		outputDir_NIfTI			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_MPRAGE_NIfTI_Files/';
	case 'MEGA-PRESS'
		% 3T BCAN MRS_and_Dopamin study
		dirData_MRS				= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_DICOM_Files/';
		outputDir_NIfTI			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_MPRAGE_NIfTI_Files/';
	case 'sLASER'
		% Select data input and output directories and other parameters depending on study
		switch strStudy_MRS
			case '3T_Trauma'
				% 3T BCAN MRS Trauma study
				% SBA
				dirData_MRS				= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MPRAGE_DICOM/';
				outputDir_NIfTI			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MPRAGE_NIfTI/';
			case '3T_SBAM'
				% 3T BCAN MRS Trauma study
				% SBAM
				dirData_MRS				= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/SBAM/MRS_SBAM_00_All_MPRAGE_DICOM/';
				outputDir_NIfTI			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/SBAM/MRS_SBAM_00_All_MPRAGE_NIfTI/';
			case 'ENIGMA_3T_SBA'
				% ENIGMA dataset selected from 3T BCAN MRS SBA study 
				dirData_MRS_base			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/ENIGMA_MRS/ENIGMA_3T_SBA/ENIGMA_3T_SBA_00_All_MRS_Files';
			case 'ENIGMA_3T_SBAM'
				% ENIGMA dataset selected from 3T BCAN MRS SBAM study
				dirData_MRS_base			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/ENIGMA_MRS/ENIGMA_3T_SBAM/ENIGMA_3T_SBAM_00_All_MRS_Files';				

			otherwise
				error('%s: ERROR: Unknown study %s!', sFunctionName, strStudy);
		end				% End of switch strStudy_MRS
		% Complete directory names based on extension/format of MRS data files 
		switch fileExt_MRS
			case 'dat'
				% MRS raw data .dat files
				dirData_MRS			= [dirData_MRS_base, '_RawData_dat/'];
				outputDir_NIfTI		= [dirData_MRS_base, '_RawData_NIfTI/'];
			case 'IMA'
				% MRS DICOM .IMA or files
				dirData_MRS			= [dirData_MRS_base, '_DICOM_IMA/'];
				outputDir_NIfTI		= [dirData_MRS_base, '_DICOM_NIfTI/'];
			case 'dcm'
				% MRS enhanced DICOM .dcm files
				error('%s: ERROR: MRS enhanced DICOM to NIfTI conversion for file extension %s NOT yet implemented!', sFunctionName, fileExt_MRS);

			otherwise
				error('%s: ERROR: Unknown file extension %s!', sFunctionName, fileExt_MRS);
		end		% End of switch fileExt_MRS

	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS);
end		% End of switch seqType_MRS

% Adjust directory names, if only new (newly acquired) data should be processed
if bProcessNewFiles
	dirData_MRS			= [dirData_MRS, '_New/'];
	outputDir_NIfTI		= [outputDir_NIfTI, '_New/'];
end
%dirData_MRS				= [dirData_MRS, filesep];
%outputDir_NIfTI			= [outputDir_NIfTI, filesep];

% If output directories for specific processing opitions do not exist, create them
if bConvert_mrs2nii
	if not(isfolder(outputDir_NIfTI))
		sMsg = sprintf('%s: Creating output directory %s for MRS NIfTI files ...\n', sFunctionName, outputDir_NIfTI);
		disp(sMsg);
		if ~mkdir(outputDir_NIfTI)
			error('%s: Could not create (mkdir) output directory %s for MRS NIfTI files!\n', sFunctionName, outputDir_NIfTI);
		end
	end
end		% End of if bConvert_mrs2nii


%% Obtain information about the list of MRS files or list of directories, respectively
% depending on MRS data type
% (assuming that all .dat files or all MRS DICOM data directories are included in the same
% directory)
% (On Linux, file list in Matlab also includes the two directories "." and "..", which
% means that the actual # of files in the directory is (# of entries in list - 2;
% however, if dir is used to list specific files, e.g. using a file extension, these two 
% directories are not included in the resulting list)
switch fileExt_MRS
    case 'dat'
		% MRS raw data .dat files
        structFileListing_MRS		= dir([dirData_MRS, '*.dat']);
        noEntriesListing_MRS		= length( structFileListing_MRS );
        %noDataFiles				= noEntriesListing_MRS - 2
    case 'IMA'
		% MRS DICOM .IMA or files
        structFileListingAll		= dir(dirData_MRS);
        subDir_MRS					= [structFileListingAll(:).isdir];
        structFileListing_MRS		= structFileListingAll(subDir_MRS);
        % Remove the two directories '.' and '..'
        structFileListing_MRS		= structFileListing_MRS(~ismember({structFileListing_MRS(:).name},{'.','..'}));
        noEntriesListing_MRS		= length(structFileListing_MRS);
	case 'dcm'
		% MRS enhanced DICOM .dcm files
		error('%s: ERROR: Case for MRS enhanced DICOM with file extension %s NOT yet implemented!', sFunctionName, fileExt_MRS);

    otherwise
        error('%s: ERROR: Unknown file extension %s!', sFunctionName, fileExt_MRS);
end		% End of switch fileExt_MRS


%% Convert all MRS data files into NIfTI format, if desired
% Assumptions:	
% - All data files are consecutively sorted, e.g. by date
% - Each set of DICOM MRS files is in one directory of the list of directories just
%	obtained

% Convert MRS data files of each directory into NIfTI format, if desired 
% using a system call that invokes the utility "spec2nii"
% and place all output NIfTI files into same directory
% (Note that here the counter for the for loop has to include all entries up to the last
% index!)
indexStart		= 1;	% To skip entries, if desired or needed
indexStep		= 1;	% Optionally adjustable step size
fprintf('\n\n');
if bConvert_mrs2nii 
	fprintf('%s: Conversion of MRS data files into NIfTI format ...\n', sFunctionName);
	for ind=indexStart : indexStep : noEntriesListing_MRS		% noEntriesListing_MRS	% 3		% 4
		dataMRS_In			= structFileListing_MRS(ind).name;
		dataMRS_InPath		= fullfile(dirData_MRS, dataMRS_In, filesep);
		fprintf('\n\n');
		fprintf('ind = %d\t .%s\t dataMRS_In = %s\n\n', ind, fileExt_MRS, dataMRS_In);
		%disp([sprintf('ind = %d\t', ind), sprintf('\t'), dataMRS_In, sprintf('\n\n')]);
		
		% Create command for conversion to NIfTI for each set of MRS data files 
		% depending on data format and type
		switch fileExt_MRS
			case 'dat'
				% MRS raw data .dat files
				% (-m can be used to specify which multi-raid file to convert if used on VE data;
				%  -m 2 refers then to the second RAID file within the .dat file that usually
				%  contains the MRS data; RAID file 1 seems to correspond to noise scans;
				%  option -m 2 does not seem to be required for the conversion command)
				%command				= sprintf('spec2nii twix -m 2 -e image -j %s -o %s', dataMRS_InPath, outputDir_NIfTI);
				command				= sprintf('spec2nii twix -e image -j %s -o %s', dataMRS_InPath, outputDir_NIfTI);
			case 'IMA'
				% MRS DICOM .IMA or files
				command				= sprintf('spec2nii dicom -j -f %s -o %s %s', dataMRS_In, outputDir_NIfTI, dataMRS_InPath);
			case 'dcm'
				% MRS enhanced DICOM .dcm files
				error('%s: ERROR: MRS enhanced DICOM to NIfTI conversion for file extension %s NOT yet implemented!', sFunctionName, fileExt_MRS);

			otherwise
				error('%s: ERROR: Unknown file extension %s!', sFunctionName, fileExt_MRS);
		end		% End of switch fileExt_MRS

		% Invoke system call for NIfTI conversion 
		[status,cmdout]		= system(command);
		if status ~= 0
			error('%s: Error in conversion of MRS data into NIfTI format for data in %s!\n\n%s', sFunctionName, dataMRS_In, cmdout);
		end
	end		% End of for ind=indexStart : indexStep : noEntriesListing_MRS
else
	fprintf('%s: No conversion of MRS data into NIfTI format!\n', sFunctionName);
end		% End of if bConvert_mrs2nii 



