%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% prepare_MetabQuant_s.m
%
%% Script to prepare metabolite quantification of magnetic resonance spectroscopy (MRS) data
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2020, 2021, 2022, 2023; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'prepare_MetabQuant_s';
sMsg_newLines		= sprintf('\n\n');
%disp(sMsg_newLines);


%% Init input parameters for preparing metabolite quantification
%dirString_Out			= '';
%dirString_Out			= '';
fileExtension           = 'dat';		% Currently: 'dat' (raw data) or 'IMA' (DICOM)
filename_In				= '';
filename_w_In			= '';
strStudy				= '3T_Trauma';		% '3T_Trauma';	'7T_KCL';	'3T_MMs';
strVOI					= 'HC'; 			% 'PCG';	% 'HC'; % 'Pons'; % 'CB'; % 'PFC'; % 'PCC';
seqType_MRS				= 'sLASER';		% 'SPECIAL';	% 'MEGA-PRESS'; % 'sLASER';
dataType_MRS			= 'mrs_w_ref';			% 'mrs_w_ref';		% 'mrs_w';
signals_MRS				= 'Spectra';		% 'MMs';	% 'Spectra';
strOVS_In				= 'wOVS';		% 'wOVS';	% 'woutOVS';
strOVS_w_In				= 'woutOVS';		% 'wOVS';	% 'woutOVS';
leftshift_In			= 2;
noSD_In					= 3.2;			% 3.2;		2.6;		4.0;
digits					= [fix(noSD_In) round(abs(noSD_In-fix(noSD_In))*10)];

% Parameters for spectral registration (aligning of averages/frequency and phase drift
% correction) performed in either frequency or time domain
driftCorr_In			= 'y';		% 'y';		'n';
iterin_In				= 20;
aaDomain_In				= 'f';		% 'f';		't';
tmaxin_In				= 0.2;		% 0.2;		0.1;
bTmaxset_In				= 1;
ppmOption				= 1;
medin_In				= 'y';		% 'y';	'n';	'a';	'ref';
% Set parameters for drift correction depending on type of data, i.e. whether MRS
% data is spectrum or water signal
% NOTE: Check whether aligning of averages in frequency domain works, if the MR
% spectrum is water signal itself; if not, simply align averages in time domain
switch dataType_MRS
	case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
		% MR spectrum is provided together without or with unsuppressed water
		% signal and/or with reference scans
		%ppmmin_fix_In		= 1.6;		% 1.6;		1.8;
		%ppmmaxarray_fix_In	= [3.5; 4.0; 5.5];
		%ppmmaxarray_fix_In	= [2.4,2.85,3.35,4.2,4.4,5.2];
		switch ppmOption
			case 1
				% For MR spectra
				ppmmin_fix_In			= 1.6;		% 1.6;		1.8;
				ppmmaxarray_fix_In		= [2.4,2.85,3.35,4.2,4.4,5.2];
			case 2
				% For MR spectra
				ppmmin_fix_In			= 1.6;
				ppmmaxarray_fix_In		= [3.5; 4.0; 5.5];
			case 3
				% For MR spectra using settings for water signals
				ppmmin_fix_In			= 4.2;
				ppmmaxarray_fix_In		= [5.5 5.5 5.2];
			case 4
				% Wide range to always inlcude water resonance
				ppmmin_fix_In			= 1.6;
				ppmmaxarray_fix_In		= [5.5 5.5 5.2];
			case 5
				% For MMs signals
				ppmmin_fix_In			= 0.2;
				ppmmaxarray_fix_In		= [3.35,4.2,4.4];
			case 6
				% For MMs signals
				ppmmin_fix_In			= 0.2;
				ppmmaxarray_fix_In		= [3.35,4.0,4.1];

			otherwise
				error('%s: Unknown ppmOption = %d!', sFunctionName, ppmOption);
		end			% End of switch ppmOption
	case {'water', 'water_ref'}
		% MR spectrum is water signal itself without or with reference scans
		ppmmin_fix_In		= 4.2;
		ppmmaxarray_fix_In	= [5.5 5.5 5.2];

	otherwise
		error('%s: Unknown MRS dataType_MRS = %s!', sFunctionName, dataType_MRS);
end		% End of switch dataType_MRS
alignSS_In				= 2;		% For aligning subspectra (e.g. in SPECIAL)
strSpecReg				= 'SR1';	% To distinguish settings for spectral registration

% Additional parameter settings
bECC_In					= 1;
bPhaseCorrFreqShift_In	= 0;
strMinUserIn_In			= 'y';
plotSwitch_In			= 0;
reportSwitch_In			= 1;


%% Additional input parameters specific to this routine
% 'bCopyFiles' used to turn on/off any copying of files (mostly used for debugging)
bCopyFiles				= 1;
bCopyFiles_MRS			= 1;
bCopyFiles_ref_Quant	= 0;
bCopyFiles_ref_ECC		= 0;
bCopyFiles_w			= 1;
bWriteFilenames			= 1;


%% Set (additional) parameters depending on sequence type
% 'dirString_Out' inidcates output directories from preprocessing of MR spectra
% Using 'dirString_Out' instead of 'dirString_In' allows to use same code for generating
% directory names as in routine preprocess_ListOfFiles_s.m
% 'outDirString_LCM' then holds information about directory, into which files are copied
% for use with metabolite quantification software (e.g. LCModel)
switch seqType_MRS
	case 'SPECIAL'
		dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_RawData_dat_Files/';
		outDirString_LCM		= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Z_Pain_Tmp/';
	case 'MEGA-PRESS'
		% Select directory for (input) data files depending on # of SDs used for
		% pre-processing of MR spectra
		switch(noSD_In)
			case(2.6)
				dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_2_6/';
			case(3.2)
				dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_3_2/';
				%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_3_2_II/';
				%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_3_2_III/';
				%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_3_2_New/';
			case(4.0)
				dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_4_0/';
			
			otherwise
				error('%s: ERROR: No directory/data for noSD_In =  %f!', sFunctionName, noSD_In);
		end			% End of switch(noSD_In)
		outDirString_LCM		= [dirString_Out, 'DOPA_LCModel_Analysis_Data/'];
		%outDirString_LCM_w		= [outDirString_LCM, 'Water_Signals/'];
		textFileName_diff_MRS	= 'list_filenames_MRS_Diff_Spectra.txt';
		textFileName_OFF		= 'list_filenames_MRS_editOFF_All.txt';
		textFileName_OFF_MRS	= 'list_filenames_MRS_editOFF_Spectra.txt';
		textFileName_OFF_w		= 'list_filenames_MRS_editOFF_Water.txt';
	case 'sLASER'
		% Select directory for input data depending on study, on # of SDs, and other 
		% options used for pre-processing of MR spectra
		%digits = [fix(noSD_In) round(abs(noSD_In-fix(noSD_In))*10)];	
		switch strStudy
			case '3T_Trauma'
				% Output data directory for preprocessing
				dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Trauma_Analysis/';
			case '3T_MMs'								
				% Output data directory
				%dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Trauma/Trauma_Z_Analysis/';
				%dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/ZZZZ_Test/';
				%dirString_Out_Base		= '/home/destiana/CSB_NeuroRad/destiana/Data_II/3T_BCAN_MRS_Trauma/MRS_Trauma_00_All_MMs/CodeResults/';
				dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Trauma_MMs_Analysis/';
			case '7T_KCL'
				% Output data directory for preprocessing
				dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/7T_KCL_Analysis/';
				
			otherwise
				error('%s: ERROR: Unknown study %s!', sFunctionName, strStudy);
		end			% End of switch strStudy

		% Complete output data directory name for preprocessing
		dirString_Out_AddOn1		= sprintf('%s_FID-A_SD_%d_%d', strVOI, digits(1), digits(2));
		dirString_Out_AddOn2		= '';
		if bECC_In
			% Use reference (water) signals for ECC, if acquired
			% If not, then use an unsuppressed water signal, if acquired
			% If no reference and no water signals are acquired, check whether MR spectrum
			% is water signal itself; and if it is, use it for ECC
			% Indicate different options for ECC in the corresponding directory name; for
			% that, search for strings 'ref', 'w', and 'water' in string for MRS data type
			refInd		= strfind(dataType_MRS, '_ref');
			wInd		= strfind(dataType_MRS, '_w');
			waterInd	= strfind(dataType_MRS, 'water');
			if ~isempty(refInd)
				dirString_Out_AddOn2	= '_ECCref';
			else
				if ~isempty(wInd)
					dirString_Out_AddOn2	= '_ECCw';
				else
					if ~isempty(waterInd)
						dirString_Out_AddOn2	= '_ECCwater';
					else
						% No reference and no water signals and MR spectrum is not water
						% signal itself => ECC not possible
						error('%s: No reference and no water signals and MR spectrum is not water signal itself (dataType_MRS = %s) => ECC not possible!', sFunctionName, dataType_MRS);
					end		% End of if ~isempty(waterInd)
				end		% End of if ~isempty(wInd)
			end		% if ~isempty(refInd)
			%dirString_Out_AddOn2	= '_ECC_Test';
		end		% End of if bECC_In
		dirString_Out			= [dirString_Out_Base, dirString_Out_AddOn1, dirString_Out_AddOn2, filesep];
		
		% Create strings for output directory and list of filenames of MR spectra and 
		% unsuppressed water signals
		outDirString_LCM		= [dirString_Out, strVOI, '_LCModel_Data/'];
		textFileName_MRS 		= 'list_filenames_MRS_Spectra.txt';
		textFileName_ref_ECC 	= 'list_filenames_MRS_Ref_ECC.txt';
		textFileName_ref_Quant 	= 'list_filenames_MRS_Ref_Quant.txt';
		textFileName_w 			= 'list_filenames_MRS_Water.txt';
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS);
end		% End of switch seqType_MRS

% If output directories for file copying do not exist, create them
if not(isfolder(outDirString_LCM))
	sMsg = sprintf('%s: Creating output directory %s ...\n', sFunctionName, outDirString_LCM);
    disp(sMsg);
    if ~mkdir(outDirString_LCM)
		error('%s: Could not create (mkdir) output directory %s!\n', sFunctionName, outDirString_LCM);
	end
end
% % Additional output directory for water signals
% if not(isfolder(outDirString_LCM_w))
% 	sMsg = sprintf('%s: Creating output directory %s ...\n', sFunctionName, outDirString_LCM_w);
%     disp(sMsg);
%     if ~mkdir(outDirString_LCM_w)
% 		error('%s: Could not create (mkdir) output directory %s!\n', sFunctionName, outDirString_LCM_w);
% 	end
% end


%% Copy data files depending on sequence type and create list of filenames, if selected
switch seqType_MRS
	case 'SPECIAL'
		% Nothing yet
	case 'MEGA-PRESS'		
		% Obtain information about the list of different groups of files
		% (assuming that all data files are included in the same directory)
		% (On Linux, file list in Matlab also includes the two directories "." and "..", which
		% means that the actual # of files in the directory is (# of entries in list - 2;
		% however, if dir is used to list specific files, e.g. using a file extension, these two
		% directories are not included in the resulting list)
		%cd(dirString_Out);
		% MEGA-PRESS difference files (MRS + water)
		structFileListing_diff		= dir([dirString_Out, '*_diff_*.RAW']);
		noEntriesListing_diff		= length( structFileListing_diff );
		%noDataFiles_diff			= noEntriesListing_diff - 2
		
		% MEGA-PRESS edit_OFF files (MRS + water)
		structFileListing_OFF		= dir([dirString_Out, '*_editOFF_*.RAW']);
		noEntriesListing_OFF		= length( structFileListing_OFF );
		
		% MEGA-PRESS water difference files (water only)
		structFileListing_diff_w	= dir([dirString_Out, '*_w_*_diff_*.RAW']);
		noEntriesListing_diff_w		= length( structFileListing_diff_w );
		
		% MEGA-PRESS water edit_OFF files (water only)
		structFileListing_OFF_w		= dir([dirString_Out, '*_w_*_editOFF_*.RAW']);
		noEntriesListing_OFF_w		= length( structFileListing_OFF_w );
		
		% Determine MEGA-PRESS MRS difference files (MRS only)
		structFileListing_diff_MRS	= structFileListing_diff( ~ismember({structFileListing_diff.name}, {structFileListing_diff_w.name}) );
		noEntriesListing_diff_MRS	= length( structFileListing_diff_MRS );
		
		% Determine MEGA-PRESS MRS edit_OFF files (MRS only)
		structFileListing_OFF_MRS	= structFileListing_OFF( ~ismember({structFileListing_OFF.name}, {structFileListing_OFF_w.name}) );
		noEntriesListing_OFF_MRS	= length( structFileListing_OFF_MRS );

		
		% Copy only MEGA-PRESS MRS difference files and all edit_OFF files with .RAW 
		% extension to destination folder, if selected
		indexStart		= 1;
		indexStep		= 1;
		if bCopyFiles
			% Copy MEGA-PRESS MRS difference files
			for ind=indexStart : indexStep : noEntriesListing_diff_MRS
				% Copy MEGA-PRESS difference file
				[status, msg] = copyfile(fullfile(dirString_Out, structFileListing_diff_MRS(ind).name), outDirString_LCM);
				if ~status
					disp(msg);
					error('%s: Could not copy file %s to output directory %s!\n', sFunctionName, structFileListing_diff_MRS(ind).name, outDirString_LCM);
				end
			end
			
			% Copy MEGA-PRESS edit_OFF files (MRS + water)
			for ind=indexStart : indexStep : noEntriesListing_OFF
				% Copy MEGA-PRESS edit_OFF file
				[status, msg] = copyfile(fullfile(dirString_Out, structFileListing_OFF(ind).name), outDirString_LCM);
				if ~status
					disp(msg);
					error('%s: Could not copy file %s to output directory %s!\n', sFunctionName, structFileListing_OFF(ind).name, outDirString_LCM);
				end
			end
			
		end		% End of if bCopyFiles
			
		
		%% Write filenames for specific groups of data files into corresponding textfiles
		% (additional for loops are used to allow for greater flexibility), if selected
		istart		= 1;
		istep		= 1;
		
		if bWriteFilenames
			% Filenames of all MEGA-PRESS edit_OFF files (MRS + water)
			fid_OFF			= fopen(fullfile(outDirString_LCM, textFileName_OFF), 'wt+');
			if fid_OFF == -1
				error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF, outDirString_LCM);
			end
			for i=istart : istep : (noEntriesListing_OFF-1)
				nbytes	= fprintf(fid_OFF, '%s\n', structFileListing_OFF(i).name);
			end
			% Write last filename without new line character at end of line to textfile
			nbytes	= fprintf(fid_OFF, '%s', structFileListing_OFF(noEntriesListing_OFF).name);
			if fclose(fid_OFF) == -1
				error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF, outDirString_LCM);
			end
			
			% Filenames of MEGA-PRESS water edit_OFF files (water only)
			fid_OFF_w		= fopen(fullfile(outDirString_LCM, textFileName_OFF_w), 'wt+');
			if fid_OFF_w == -1
				error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF_w, outDirString_LCM);
			end
			for i=istart : istep : (noEntriesListing_OFF_w-1)
				nbytes	= fprintf(fid_OFF_w, '%s\n', structFileListing_OFF_w(i).name);
			end
			% Write last filename without new line character at end of line to textfile
			nbytes	= fprintf(fid_OFF_w, '%s', structFileListing_OFF_w(noEntriesListing_OFF_w).name);
			if fclose(fid_OFF_w) == -1
				error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF_w, outDirString_LCM);
			end
			
			% Filenames of MEGA-PRESS MRS edit_OFF files (MRS only)
			fid_OFF_MRS		= fopen(fullfile(outDirString_LCM, textFileName_OFF_MRS), 'wt+');
			if fid_OFF_MRS == -1
				error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF_MRS, outDirString_LCM);
			end
			for i=istart : istep : (noEntriesListing_OFF_MRS-1)
				nbytes	= fprintf(fid_OFF_MRS, '%s\n', structFileListing_OFF_MRS(i).name);
			end
			% Write last filename without new line character at end of line to textfile
			nbytes	= fprintf(fid_OFF_MRS, '%s', structFileListing_OFF_MRS(noEntriesListing_OFF_MRS).name);
			if fclose(fid_OFF_MRS) == -1
				error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF_MRS, outDirString_LCM);
			end
			
			% Filenames of MEGA-PRESS MRS difference files (MRS only)
			fid_diff_MRS		= fopen(fullfile(outDirString_LCM, textFileName_diff_MRS), 'wt+');
			if fid_diff_MRS == -1
				error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_diff_MRS, outDirString_LCM);
			end
			for i=istart : istep : (noEntriesListing_diff_MRS-1)
				nbytes	= fprintf(fid_diff_MRS, '%s\n', structFileListing_diff_MRS(i).name);
			end
			% Write last filename without new line character at end of line to textfile
			nbytes	= fprintf(fid_diff_MRS, '%s', structFileListing_diff_MRS(noEntriesListing_diff_MRS).name);
			if fclose(fid_diff_MRS) == -1
				error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_diff_MRS, outDirString_LCM);
			end
			
		end		% End of if bWriteFilenames
	case 'sLASER'
		% Obtain information about the list of different groups of files
		% (assuming that all data files are included in the same directory)
		
		% sLASER all processed .RAW files
		structFileListing_all		= dir([dirString_Out, '*_processed_*.RAW']);
		noEntriesListing_all		= length( structFileListing_all );
		%noDataFiles_all			= noEntriesListing_all - 2
		
		% sLASER processed .RAW water reference files for ECC
		structFileListing_ref_ECC	= dir([dirString_Out, '*_ref_ECC*_processed_*.RAW']);
		noEntriesListing_ref_ECC	= length( structFileListing_ref_ECC );
		
		% sLASER processed .RAW water reference files for quantification (Quant)
		structFileListing_ref_Quant	= dir([dirString_Out, '*_ref_Quant*_processed_*.RAW']);
		noEntriesListing_ref_Quant	= length( structFileListing_ref_Quant );
				
		% sLASER processed .RAW unsuppressed water files
		structFileListing_w		= dir([dirString_Out, '*_w_*_processed_*.RAW']);
		noEntriesListing_w		= length( structFileListing_w );
		
 		% sLASER processed .RAW MR spectra files 
		% (find .RAW files that do not belong to any other group)
		structFileListing_Tmp1		= structFileListing_all( ~ismember({structFileListing_all.name}, {structFileListing_ref_Quant.name}) );
		structFileListing_Tmp2		= structFileListing_Tmp1( ~ismember({structFileListing_Tmp1.name}, {structFileListing_ref_ECC.name}) );
 		structFileListing_MRS		= structFileListing_Tmp2( ~ismember({structFileListing_Tmp2.name}, {structFileListing_w.name}) );
 		noEntriesListing_MRS		= length( structFileListing_MRS );

		% Copy (processed) MR Spectra files, water reference signals for quantification,
		% and unsuppresed water signals
		% with . RAW extension to destination folder, if selected
		indexStart		= 1;
		indexStep		= 1;
		if bCopyFiles
			% Copy (processed) MR spectra files, if selected
			if bCopyFiles_MRS
				for ind=indexStart : indexStep : noEntriesListing_MRS
					% Copy MR spectra file
					[status, msg] = copyfile(fullfile(dirString_Out, structFileListing_MRS(ind).name), outDirString_LCM);
					if ~status
						disp(msg);
						error('%s: Could not copy file %s to output directory %s!\n', sFunctionName, structFileListing_MRS(ind).name, outDirString_LCM);
					end
				end
			end		% End of if bCopyFiles_MRS
			
			% Copy (processed) water reference files quantification (Quant), if selected
			if bCopyFiles_ref_Quant
				for ind=indexStart : indexStep : noEntriesListing_ref_Quant
					% Copy water reference file quantification (Quant)
					[status, msg] = copyfile(fullfile(dirString_Out, structFileListing_ref_Quant(ind).name), outDirString_LCM);
					if ~status
						disp(msg);
						error('%s: Could not copy file %s to output directory %s!\n', sFunctionName, structFileListing_ref_Quant(ind).name, outDirString_LCM);
					end
				end
			end		% End of if bCopyFiles_ref_Quant
			
			% Copy (processed) water reference files for ECC, if selected
			if bCopyFiles_ref_ECC
				for ind=indexStart : indexStep : noEntriesListing_ref_ECC
					% Copy water reference file for ECC
					[status, msg] = copyfile(fullfile(dirString_Out, structFileListing_ref_ECC(ind).name), outDirString_LCM);
					if ~status
						disp(msg);
						error('%s: Could not copy file %s to output directory %s!\n', sFunctionName, structFileListing_ref_ECC(ind).name, outDirString_LCM);
					end
				end
			end		% End of if bCopyFiles_ref_ECC
			
			% Copy (processed) unsuppressed water files, if selected
			if bCopyFiles_w
				for ind=indexStart : indexStep : noEntriesListing_w
					% Copy unsuppressed water file
					[status, msg] = copyfile(fullfile(dirString_Out, structFileListing_w(ind).name), outDirString_LCM);
					if ~status
						disp(msg);
						error('%s: Could not copy file %s to output directory %s!\n', sFunctionName, structFileListing_w(ind).name, outDirString_LCM);
					end
				end
			end		% End of if bCopyFiles_w
			
		end		% End of if bCopyFiles
		
		
		%% Write filenames for specific groups of data files into corresponding textfiles
		% (additional for loops are used to allow for greater flexibility), if selected
		istart		= 1;
		istep		= 1;
		if bWriteFilenames
			% Filenames of (processed) MR spectra files
			if bCopyFiles_MRS
				fid_MRS			= fopen(fullfile(outDirString_LCM, textFileName_MRS), 'wt+');
				if fid_MRS == -1
					error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_MRS, outDirString_LCM);
				end
				for i=istart : istep : (noEntriesListing_MRS-1)
					nbytes	= fprintf(fid_MRS, '%s\n', structFileListing_MRS(i).name);
				end
				% Write last filename without new line character at end of line to textfile
				nbytes	= fprintf(fid_MRS, '%s', structFileListing_MRS(noEntriesListing_MRS).name);
				if fclose(fid_MRS) == -1
					error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_MRS, outDirString_LCM);
				end
			end		% End of if bCopyFiles_MRS
			
			% Filenames of (processed) water reference files for quantification (Quant)
			if bCopyFiles_ref_Quant
				fid_ref_Quant	= fopen(fullfile(outDirString_LCM, textFileName_ref_Quant), 'wt+');
				if fid_ref_Quant == -1
					error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_ref_Quant, outDirString_LCM);
				end
				for i=istart : istep : (noEntriesListing_ref_Quant-1)
					nbytes	= fprintf(fid_ref_Quant, '%s\n', structFileListing_ref_Quant(i).name);
				end
				% Write last filename without new line character at end of line to textfile
				nbytes	= fprintf(fid_ref_Quant, '%s', structFileListing_ref_Quant(noEntriesListing_ref_Quant).name);
				if fclose(fid_ref_Quant) == -1
					error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_ref_Quant, outDirString_LCM);
				end
			end		% End of if bCopyFiles_ref_Quant
			
			% Filenames of (processed) water reference files for ECC
			if bCopyFiles_ref_ECC
				fid_ref_ECC		= fopen(fullfile(outDirString_LCM, textFileName_ref_ECC), 'wt+');
				if fid_ref_ECC == -1
					error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_ref_ECC, outDirString_LCM);
				end
				for i=istart : istep : (noEntriesListing_ref_ECC-1)
					nbytes	= fprintf(fid_ref_ECC, '%s\n', structFileListing_ref_ECC(i).name);
				end
				% Write last filename without new line character at end of line to textfile
				nbytes	= fprintf(fid_ref_ECC, '%s', structFileListing_ref_ECC(noEntriesListing_ref_ECC).name);
				if fclose(fid_ref_ECC) == -1
					error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_ref_ECC, outDirString_LCM);
				end
			end		% End of if bCopyFiles_ref_ECC
			
			% Filenames of (processed) unsuppressed water files
			if bCopyFiles_w
				fid_w			= fopen(fullfile(outDirString_LCM, textFileName_w), 'wt+');
				if fid_w == -1
					error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_w, outDirString_LCM);
				end
				for i=istart : istep : (noEntriesListing_w-1)
					nbytes	= fprintf(fid_w, '%s\n', structFileListing_w(i).name);
				end
				% Write last filename without new line character at end of line to textfile
				nbytes	= fprintf(fid_w, '%s', structFileListing_w(noEntriesListing_w).name);
				if fclose(fid_w) == -1
					error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_w, outDirString_LCM);
				end
			end		% End of if bCopyFiles_w
			
		end		% End of if bWriteFilenames
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS);
end		% End of switch seqType_MRS


%% Save variables of workspace to file
% Obtain current date and time in specific format
dt		= datestr(now,'yyyymmdd_HH_MM_SS');

% Save workspace into output directory (optional with user input)
% (Extension".mat" in filename explicitly required, so that Matlab can correctly load 
% workspace file with a "." in its filename)
%strSavedWorkspaceFileName		= 'workspace_run_specialproc_CBF';
strSavedWorkspaceFileName		= ['workspace_', sFunctionName, '_', seqType_MRS, '_', dt];
strSavedWorkspaceFileNameFull	= [outDirString_LCM, strSavedWorkspaceFileName, sprintf('_SD_%.1f.mat', noSD_In)];
%strSaveWorkspace	= input('Would you like to save all variables of the workspace to file?  ', 's');
strSaveWorkspace	= 'y';
if strcmp(strSaveWorkspace,'y') || strcmp(strSaveWorkspace,'Y')
	save(strSavedWorkspaceFileNameFull);
end

