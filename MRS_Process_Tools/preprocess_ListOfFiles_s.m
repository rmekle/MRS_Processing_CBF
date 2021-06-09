%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% preprocess_ListOfFiles_s.m
%
%% Script to preprocess a list of files of magnetic resonance spectroscopy (MRS) data
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2018, 2019, 2020; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'preprocess_ListOfFiles_s';
sMsg_newLines		= sprintf('\n\n');
%disp(sMsg_newLines);


%% Init input parameters for preprocessing routine
%dirString_In			= '';
%outDirString_In			= '';
filename_In				= '';
filenamew_In			= '';
noSD_In					= 4.0;			% 3.2;		2.6;		4.0;
strOVS_In				='wOVS';
strMinUserIn_In			= 'y';
aaDomain_In				= 'f';
tmaxin_In				= 0.2;
iterin_In				= 20;
alignSS_In				= 2;
seqType					= 'MEGA-PRESS';		% 'SPECIAL';	% 'MEGA-PRESS';

% Set (additional) parameters depending on sequence type
switch seqType
	case 'SPECIAL'
		dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_RawData_dat_Files/';
		outDirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Z_Pain_Tmp/';
		%outDirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Preprocessed_forLCModel_SD_4_0/';
		%outDirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Preprocessed_forLCModel_SD_3_2/';
		%outDirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Preprocessed_forLCModel_SD_2_6/';
	case 'MEGA-PRESS'
		% Data (input) directories
		dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_RawData_dat_Files_MRS/';
		%dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Dopamin/MRS_Dopamin_00_All_RawData_dat_Files_MRS_New/';
		
		% Select directories for output data depending on # of SDs used for pre-processing
		% of MR spectra
		switch(noSD_In)
			case(2.6)
				outDirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_2_6/';
			case(3.2)
				%outDirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2/';
				outDirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2_New/';
				%outDirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_3_2/';
			case(4.0)
				outDirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_4_0/';
				%outDirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_4_0/';
			
			otherwise
				error('%s: ERROR: No directory/data for noSD_In =  %f!', sFunctionName, noSD_In);
		end
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType);
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
switch seqType
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
			filenamew_In	= structFileListing(ind+1).name;
			disp(sMsg_newLines);
			disp([sprintf('ind = %d\t', ind), strOVS_In, sprintf('\t'), filename_In, sprintf('\t'), filenamew_In, sprintf('\n\n')]);
			[out,out_w,out_noproc,out_w_noproc]=run_specialproc_CBF(dirString_In,outDirString_In,filename_In,filenamew_In,noSD_In,strOVS_In,strMinUserIn_In,aaDomain_In,tmaxin_In,iterin_In);
			
			% Close all figures
			close all;
			
			% Preprocess MR spectrum and water_woutOVS
			% (spectral data remain the same and only processed differently now)
			strOVS_In		= 'woutOVS';
			filenamew_In	= structFileListing(ind+2).name;
			disp(sMsg_newLines);
			disp([sprintf('ind = %d\t', ind), strOVS_In, sprintf('\t'), filename_In, sprintf('\t'), filenamew_In, sprintf('\n\n')]);
			[out,out_w,out_noproc,out_w_noproc]=run_specialproc_CBF(dirString_In,outDirString_In,filename_In,filenamew_In,noSD_In,strOVS_In,strMinUserIn_In,aaDomain_In,tmaxin_In,iterin_In);
		end
	case 'MEGA-PRESS'
		% Here preprocessing of MEGA-PRESS MR spectra together with the corresponding 
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
			filenamew_In	= structFileListing(ind+1).name;
			disp(sMsg_newLines);
			disp([sprintf('ind = %d\t', ind), sprintf('\t'), filename_In, sprintf('\t'), filenamew_In, sprintf('\n\n')]);
			[diffSpecOut,sumSpecOut,subSpec1Out,subSpec2Out,outwOut,outw_subSpec1Out,outw_subSpec2Out,coilcombosOut]=run_megapressproc_CBF(dirString_In,outDirString_In,filename_In,filenamew_In,noSD_In,strMinUserIn_In,aaDomain_In,tmaxin_In,iterin_In,alignSS_In);
			
			% Close all figures
			%close all;
		end		
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
strSavedWorkspaceFileNameFull	= [outDirString_In, strSavedWorkspaceFileName, sprintf('_SD_%.1f.mat', noSD_In)];
%strSaveWorkspace	= input('Would you like to save all variables of the workspace to file?  ', 's');
strSaveWorkspace	= 'y';
if strcmp(strSaveWorkspace,'y') || strcmp(strSaveWorkspace,'Y')
	save(strSavedWorkspaceFileNameFull);
end

