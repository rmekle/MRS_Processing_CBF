%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% preprocess_ListOfFiles_s.m
%
%% Script to preprocess a list of files of magnetic resonance spectroscopy (MRS) data
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2018; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Init input parameters for preprocessing routine
%dirString_In			= '';
%outDirString_In			= '';
filename_In				= '';
filenamew_In			= '';
noSD_In					= 3.2;			% 3.2;		2.6;		4.0;
strOVS_In				='wOVS';
strMinUserIn_In			= 'y';
aaDomain_In				= 'f';
tmaxin_In				= 0.2;
iterin_In				= 20;

% Set (additional) parameters
dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_Potsdam_Pain/Potsdam_Pain_00_All_RawData_dat_Files/';
outDirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Z_Pain_Tmp/';
%outDirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Preprocessed_forLCModel_SD_4_0/';
%outDirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Preprocessed_forLCModel_SD_3_2/';
%outDirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Preprocessed_forLCModel_SD_2_6/';


%% Obtain information about the list of files
% (assuming that all data files are included in the same directory)
% (On Linux, file list in Matlab also includes the two directories "." and "..", which
% means that the actual # of files in the directory is (# of entries in list - 2)
%cd(dirString_In);
structFileListing		= dir(dirString_In);
noEntriesListing		= length( structFileListing );
noDataFiles				= noEntriesListing - 2


%% Preprocess all data files
% Here preprocessing of SPECIAL MR spectra together with the corresponding water files is
% performed
% Assumptions:	
% - All data files are consecutively sorted, e.g. by date
% - Data files come in group of three files:
% - SPECIAL MR spectrum, SPECIAL water withOVS (wOVS), SPECIAL water withoutOVS (woutOVS)
% - Order of files is the same for all cases, i.e. spectrum; water_wOVS, water_woutOVS

% Preprocess each case (spectrum) once using water_wOVS and water_woutOVS
% Water signal is used for coil combination of the spectrum and is preprocessed as well
indexStart		= 3;	% To skip entries for directories "." and ".."
indexStep		= 3;	% For SPECIAL spectra, since two different water signals exist
for ind=indexStart : indexStep : 3		% noDataFiles	% 6		% 3
	% Preprocess MR spectrum and water_withOVS
	strOVS_In		= 'wOVS';
	filename_In		= structFileListing(ind).name;
	filenamew_In	= structFileListing(ind+1).name;
	disp(sprintf('\n\n'));
	disp([sprintf('ind = %d\t', ind), strOVS_In, sprintf('\t'), filename_In, sprintf('\t'), filenamew_In, sprintf('\n\n')]);
	[out,out_w,out_noproc,out_w_noproc]=run_specialproc_CBF(dirString_In,outDirString_In,filename_In,filenamew_In,noSD_In,strOVS_In,strMinUserIn_In,aaDomain_In,tmaxin_In,iterin_In);
	
	% Close all figures
	close all;
	
	% Preprocess MR spectrum and water_woutOVS
	% (spectral data remain the same and only processed differently now)
	strOVS_In		= 'woutOVS';
	filenamew_In	= structFileListing(ind+2).name;
	disp(sprintf('\n\n'));
	disp([sprintf('ind = %d\t', ind), strOVS_In, sprintf('\t'), filename_In, sprintf('\t'), filenamew_In, sprintf('\n\n')]);
	[out,out_w,out_noproc,out_w_noproc]=run_specialproc_CBF(dirString_In,outDirString_In,filename_In,filenamew_In,noSD_In,strOVS_In,strMinUserIn_In,aaDomain_In,tmaxin_In,iterin_In);
end


%% Save variables of workspace to file
% Save workspace into output directory (optional with user input)
% (Extension".mat" in filename explicitly required, so that Matlab can correctly load 
% workspace file with a "." in its filename)
strSavedWorkspaceFileName		= 'workspace_run_specialproc_CBF';
strSavedWorkspaceFileNameFull	= [outDirString_In, strSavedWorkspaceFileName, sprintf('_SD_%.1f.mat', noSD_In)];
%strSaveWorkspace	= input('Would you like to save all variables of the workspace to file?  ', 's');
strSaveWorkspace	= 'y';
if strcmp(strSaveWorkspace,'y') || strcmp(strSaveWorkspace,'Y')
	save(strSavedWorkspaceFileNameFull);
end

