%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% preprocess_ListOfFiles_s.m
%
%% Script to preprocess a list of files of magnetic resonance spectroscopy (MRS) data
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2018, 2019, 2020, 2021; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'preprocess_ListOfFiles_s';
sMsg_newLines		= sprintf('\n\n');
sMsg_newLine		= sprintf('\n');
disp(sMsg_newLines);


%% Init input parameters for preprocessing routine
%dirString_In			= '';
%dirString_Out			= '';
filename_In				= '';
filename_w_In			= '';
strVOI					= 'PCG';			% 'HC';		% 'PCG';
seqType_MRS				= 'sLASER';		% 'SPECIAL';	% 'MEGA-PRESS'; % 'sLASER';
dataType_MRS			= 'mrs_w_ref';
strOVS_In				= 'wOVS';
strOVS_w_In				= 'wOVS';
leftshift_In			= 3;
noSD_In					= 3.2;			% 3.2;		2.6;		4.0;
strMinUserIn_In			= 'y';
aaDomain_In				= 'f';
tmaxin_In				= 0.2;
iterin_In				= 20;
alignSS_In				= 2;
bECC_In					= 1;
bPhaseCorrFreqShift_In	= 0;
plotSwitch_In			= 0;
reportSwitch_In			= 1;

% Set (additional) parameters depending on sequence type
switch seqType_MRS
	case 'SPECIAL'
		dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_RawData_dat_Files/';
		dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Z_Pain_Tmp/';
		%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Preprocessed_forLCModel_SD_4_0/';
		%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Preprocessed_forLCModel_SD_3_2/';
		%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Preprocessed_forLCModel_SD_2_6/';
	case 'MEGA-PRESS'
		% Data (input) directories
		dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_RawData_dat_Files_MRS/';
		%dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_RawData_dat_Files_MRS_New/';
		
		% Select directories for output data depending on # of SDs used for pre-processing
		% of MR spectra
		switch(noSD_In)
			case(2.6)
				dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_2_6/';
			case(3.2)
				%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2/';
				dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2_New/';
				%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_3_2/';
			case(4.0)
				dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_4_0/';
				%dirString_Out			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_4_0/';
			
			otherwise
				error('%s: ERROR: No directory/data for noSD_In =  %f!', sFunctionName, noSD_In);
		end
	case 'sLASER'
		% Data (input) directories
		dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma/';
		dirString_In_AddOn1		= sprintf('MRS_Trauma_00_All_RawData_dat_Files_MRS_%s', strVOI);
		dirString_In			= [dirString_In_Base, dirString_In_AddOn1, filesep];
		
		% Select directory for output data depending on # of SDs and other options used 
		% for pre-processing of MR spectra
		digits = [fix(noSD_In) round(abs(noSD_In-fix(noSD_In))*10)];
		dirString_Out_Base		= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Trauma/Trauma_Z_Analysis/';
		dirString_Out_AddOn1	= sprintf('%s_FID-A_SD_%d_%d', strVOI, digits(1), digits(2));
		dirString_Out_AddOn2	= '';
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
		end		% End of if bECC_In
		dirString_Out			= [dirString_Out_Base, dirString_Out_AddOn1, dirString_Out_AddOn2, filesep];
		%dirString_Out			= [dirString_Out_Base, dirString_Out_AddOn1, dirString_Out_AddOn2, '_Test', filesep];
		
		% If outout directory is non-existent, create it
		if ~exist( dirString_Out, 'dir' )
			mkdir(dirString_Out);
		end
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS);
end



%% Obtain information about the list of files
% (assuming that all data files are included in the same directory)
% (On Linux, file list in Matlab also includes the two directories "." and "..", which
% means that the actual # of files in the directory is (# of entries in list - 2;
% however, if dir is used to list specific files, e.g. using a file extension, these two 
% directories are not included in the resulting list)
%cd(dirString_In);
structFileListing		= dir([dirString_In, '*.dat']);
noEntriesListing		= length( structFileListing );
%noDataFiles				= noEntriesListing - 2


%% Preprocess all data files depending on sequence type
switch seqType_MRS
	case 'SPECIAL'
		% Here preprocessing of SPECIAL MR spectra together with the corresponding water 
		% files is performed
		% Assumptions:
		% - All data files are consecutively sorted, e.g. by date
		% - Data files come in group of three files:
		% - SPECIAL MR spectrum, SPECIAL water withOVS (wOVS), SPECIAL water withoutOVS (woutOVS)
		% - Order of files is same for all cases, i.e. spectrum; water_wOVS, water_woutOVS
		
		% Preprocess each case (spectrum) once using water_wOVS and water_woutOVS
		% Water signal is used for coil combination of the spectrum and is preprocessed 
		% as well
		%indexStart		= 3;	% To skip entries for directories "." and ".."
		indexStart		= 1;
		indexStep		= 3;	% For SPECIAL, since two different water signals exist
		for ind=indexStart : indexStep : 3		% noEntriesListing	% 6		% 3
			% Preprocess MR spectrum and water_withOVS
			strOVS_In		= 'wOVS';
			filename_In		= structFileListing(ind).name;
			filename_w_In	= structFileListing(ind+1).name;
			disp(sMsg_newLines);
			disp([sprintf('ind = %d\t', ind), strOVS_In, sprintf('\t'), filename_In, sprintf('\t'), filename_w_In, sprintf('\n\n')]);
			[out,out_w,out_noproc,out_w_noproc]=run_specialproc_CBF(dirString_In,dirString_Out,filename_In,filename_w_In,noSD_In,strOVS_In,strMinUserIn_In,aaDomain_In,tmaxin_In,iterin_In);
			
			% Close all figures
			close all;
			
			% Preprocess MR spectrum and water_woutOVS
			% (spectral data remain the same and only processed differently now)
			strOVS_In		= 'woutOVS';
			filename_w_In	= structFileListing(ind+2).name;
			disp(sMsg_newLines);
			disp([sprintf('ind = %d\t', ind), strOVS_In, sprintf('\t'), filename_In, sprintf('\t'), filename_w_In, sprintf('\n\n')]);
			[out,out_w,out_noproc,out_w_noproc]=run_specialproc_CBF(dirString_In,dirString_Out,filename_In,filename_w_In,noSD_In,strOVS_In,strMinUserIn_In,aaDomain_In,tmaxin_In,iterin_In);
		end
	case 'MEGA-PRESS'
		% Here preprocessing of MEGA-PRESS MR spectraFSN-FZ-CSB-08 together with the corresponding 
		% water file is performed
		% Assumptions:
		% - All data files are consecutively sorted, e.g. by date
		% - Data files come in group of two files:
		% - MEGA-PRESS MR spectrum, MEGA-PRESS water signal
		% - Order of files is same for all cases, i.e. spectrum; water
		
		% Preprocess each case (spectrum)
		% Water signal is used for coil combination of the spectrum and is preprocessed 
		% as well
		indexStart		= 1;
		indexStep		= 2;	% For MEGA-PRESS, since only one water signals exists
		for ind=indexStart : indexStep : noEntriesListing		% noEntriesListing	% 4		% 2
			% Preprocess MR spectrum and water
			filename_In		= structFileListing(ind).name;
			filename_w_In	= structFileListing(ind+1).name;
			disp(sMsg_newLines);
			disp([sprintf('ind = %d\t', ind), sprintf('\t'), filename_In, sprintf('\t'), filename_w_In, sprintf('\n\n')]);
			[diffSpecOut,sumSpecOut,subSpec1Out,subSpec2Out,outwOut,outw_subSpec1Out,outw_subSpec2Out,coilcombosOut]=run_megapressproc_CBF(dirString_In,dirString_Out,filename_In,filename_w_In,noSD_In,strMinUserIn_In,aaDomain_In,tmaxin_In,iterin_In,alignSS_In);
			
			% Close all figures
			%close all;
		end
	case 'sLASER'
		% Here preprocessing of sLASER MR spectra acquired with or without water reference
		% signals together with the corresponding water file is performed
		% Assumptions:
		% - All data files are consecutively sorted, e.g. by date
		% - Data files come in group of two files:
		% - sLASER MR spectrum with or without reference scans, sLASER water signal
		% - Order of files is same for all cases, i.e. spectrum; water
		
		% Preprocess each case (spectrum)
		% If separate water scans and reference scans were acquired, get coil phases from
		% both types of signals, and for coil combination use 
		% if reference scans exist,
		%	coil phases from reference scans for reference scans and MR spectra, since
		%	reference scans were acquired together with MR spectra
		%	coil phases from water scans for water scans
		%
		% if only water scans exist, 
		%	coil phases from water scans for water scans and for MR spectra
		%
		% if neither reference scans nor water scans eFSN-FZ-CSB-08xist,
		%	coil phases from MR spectra for MR spectra
		indexStart		= 1;
		indexStep		= 2;	% For sLASER, since only one water signals exists
		for ind=indexStart : indexStep : noEntriesListing		% noEntriesListing	% 4		% 2
			% Preprocess MR spectrum and water
			filename_In		= structFileListing(ind).name;
			filename_w_In	= structFileListing(ind+1).name;
			disp(sMsg_newLines);
			disp([sprintf('ind = %d\t', ind), sprintf('\t'), filename_In, sprintf('\t'), filename_w_In, sprintf('\n\n')]);
			[out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_RawData_s(dirString_In,dirString_Out,filename_In,filename_w_In,seqType_MRS,dataType_MRS,strOVS_In,strOVS_w_In,leftshift_In,noSD_In,aaDomain_In,tmaxin_In,iterin_In,bECC_In,bPhaseCorrFreqShift_In,plotSwitch_In,strMinUserIn_In,reportSwitch_In);
			
			% Close all figures
			%close all;
		end		
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS);
end


%% Save variables of workspace to file
% Obtain current date and time in specific format
dt		= datestr(now,'yyyymmdd_HH_MM_SS');

% Save workspace into output directory (optional with user input)
% (Extension".mat" in filename explicitly required, so that Matlab can correctly load 
% workspace file with a "." in its filename)
%strSavedWorkspaceFileName		= 'workspace_run_specialproc_CBF';
strSavedWorkspaceFileName		= ['workspace_', sFunctionName, '_', seqType_MRS, '_', dataType_MRS, '_', dt];
strSavedWorkspaceFileNameFull	= [dirString_Out, strSavedWorkspaceFileName, sprintf('_SD_%.1f.mat', noSD_In)];
%strSaveWorkspace	= input('Would you like to save all variables of the workspace to file?  ', 's');
strSaveWorkspace	= 'y';
if strcmp(strSaveWorkspace,'y') || strcmp(strSaveWorkspace,'Y')
	save(strSavedWorkspaceFileNameFull);
end

