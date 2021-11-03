%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% run_LCModel_ListOfFiles_s.m
%
%% Script to run LCModel quantification on a list of files of magnetic resonance 
%% spectroscopy (MRS) data
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2018, 2019, 2020; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'run_LCModel_ListOfFiles_s';
sMsg_newLines		= sprintf('\n\n');
sMsg_newLine		= sprintf('\n');
disp(sMsg_newLines);


%% Init parameter settings
strVOI					= 'HC';			% 'HC';		% 'PCG';
seqType_MRS				= 'sLASER';		% 'SPECIAL';	% 'MEGA-PRESS'; % 'sLASER';
dataType_MRS			= 'mrs_w_ref';
% strOVS_In				= 'wOVS';
% strOVS_w_In				= 'wOVS';
noSD_In					= 3.2;			% 3.2;		2.6;		4.0;
% strMinUserIn_In			= 'y';
% aaDomain_In				= 'f';
% tmaxin_In				= 0.2;
% iterin_In				= 20;
% alignSS_In				= 2;
bECC_In					= 0;
strTissue				= 'GM';		% 'GM';	% 'WM';	% 'HC';	% 'PCG';
strWaterQuant			= '_ref_Quant';		% '_ref_Quant'; % '_ref_ECC';	% '_w'; %'';
b0nratio				= 1;
strAnalysisData			= 'MRS_editOFF';	% 'MRS_diff';	'MRS_editOFF';	'MRS_reg';
charVecSaveWorkspace	= 'y';


%% Set directories for data, results, and filenames for lists of data files 
% Sequence independent settings
filename_MRS					= '';
filename_MRS_w					= '';
filename_MRS_noExt				= '';
filename_MRS_w_noExt			= '';

% Settings depending on sequence type
switch seqType_MRS
	case 'SPECIAL'
		dirData							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Processed_forLCModel_SD_3_2/LCModel_Analysis_Data/';
		outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Processed_forLCModel_SD_3_2/LCModel_Results_nratio0/';
		%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/Potsdam_Pain/PotsdamPain_DataAnalysis/Processed_forLCModel_SD_3_2/LCModel_Results/';
		fullFilename_listOfFiles_MRS	= [dirData, 'list_filenames_MRS_Spectra.txt'];
		fullFilename_listOfFiles_MRS_w	= [dirData, 'list_filenames_MRS_Water.txt'];
		
	case 'MEGA-PRESS'
		% Select directories for (input) data files and for output data depending on # of 
		% SDs used for pre-processing of MR spectra and on analysis settings
		switch(noSD_In)
			case(2.6)
				dirData							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_2_6/DOPA_LCModel_Analysis_Data/';
				
				% Select output directory depending on type of data (spectra) used for
				% analysis
				if strcmp(strAnalysisData, 'MRS_diff')
					% MRS_diff
					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_2_6/DOPA_LCM_Out_noECC/';
					outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_2_6/DOPA_LCM_Out_noECC_0nratio/';
				else
					% MRS_editOFF
					outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_2_6/DOPA_LCM_Out_editOFF_GM_w_0nratio/';
				end
			case(3.2)
				dirData							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2/DOPA_LCModel_Analysis_Data/';
				%dirData							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2_New/DOPA_LCModel_Analysis_Data/';
				%dirData							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_3_2/DOPA_LCModel_Analysis_Data/';
				
				% Select output directory depending on type of data (spectra) used for
				% analysis
				if strcmp(strAnalysisData, 'MRS_diff')
					% MRS_diff
					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2/DOPA_LCM_Out_noECC/';
					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2/DOPA_LCM_Out_noECC_0nratio/';
					outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2_New/DOPA_LCM_Out_noECC/';
					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2_New/DOPA_LCM_Out_noECC_0nratio/';
					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_forLCModel_SD_3_2/DOPA_LCModel_Results_noECC/';
					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_3_2/DOPA_LCModel_Results_noECC/';
					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_3_2/DOPA_LCModel_Results_noECC_0nratio/';
					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_3_2/DOPA_LCModel_Results_noECC_noWS/';
				else
					% MRS_editOFF
					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2/DOPA_LCM_Out_editOFF_GM_w_0nratio_noECC/';
					outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2/DOPA_LCM_Out_editOFF_GM_w_0nratio/';
					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2/DOPA_LCM_Out_editOFF_GM_w/';
				end
			case(4.0)
				dirData							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_4_0/DOPA_LCModel_Analysis_Data/';
				%dirData							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_4_0/DOPA_LCModel_Analysis_Data/';

				% Select output directory depending on type of data (spectra) used for
				% analysis
				if strcmp(strAnalysisData, 'MRS_diff')
					% MRS_diff
					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_4_0/DOPA_LCM_Out_noECC/';
					outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_4_0/DOPA_LCM_Out_noECC_0nratio/';
					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_4_0/DOPA_LCModel_Results/';
					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_4_0/DOPA_LCModel_Results_noECC/';
					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_4_0/DOPA_LCModel_Results_noECC_noWS/';
				else
					% MRS_editOFF
					outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_4_0/DOPA_LCM_Out_editOFF_GM_w_0nratio/';
				end
				
			otherwise
				error('%s: ERROR: No directory/data for noSD_In =  %f!', sFunctionName, noSD_In);
		end
		% Lists of filenames for MRS spectra and water signals depending on type of data 
		% (spectra) used for analysis
		if strcmp(strAnalysisData, 'MRS_diff')
			% MRS_diff
			fullFilename_listOfFiles_MRS	= [dirData, 'list_filenames_MRS_Diff_Spectra.txt'];
			fullFilename_listOfFiles_MRS_w	= [dirData, 'list_filenames_MRS_editOFF_Water.txt'];
		else
			% MRS_editOFF
			fullFilename_listOfFiles_MRS	= [dirData, 'list_filenames_MRS_editOFF_Spectra.txt'];
			fullFilename_listOfFiles_MRS_w	= [dirData, 'list_filenames_MRS_editOFF_Water.txt'];
		end
	case 'sLASER'
		% Select directories for (input) data files and for output data depending on # of
		% SDs and other options used for pre-processing of MR spectra and on settings for 
		% LCModel analysis
		digits = [fix(noSD_In) round(abs(noSD_In-fix(noSD_In))*10)];
		dirData_Base		= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Trauma/Trauma_Z_Analysis/';
		dirData_AddOn1		= sprintf('%s_FID-A_SD_%d_%d', strVOI, digits(1), digits(2));
		dirData_AddOn2		= '';
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
				dirData_AddOn2	= '_ECCref';
			else
				if ~isempty(wInd)
					dirData_AddOn2	= '_ECCw';
				else
					if ~isempty(waterInd)
						dirData_AddOn2	= '_ECCwater';
					else
						% No reference and no water signals and MR spectrum is not water
						% signal itself => ECC not possible
						error('%s: No reference and no water signals and MR spectrum is not water signal itself (dataType_MRS = %s) => ECC not possible!', sFunctionName, dataType_MRS);
					end		% End of if ~isempty(waterInd)
				end		% End of if ~isempty(wInd)
			end		% if ~isempty(refInd)
			%dirData_AddOn2	= '_ECC_Test';
		end		% End of if bECC_In
		dirData_Processed	= [dirData_Base, dirData_AddOn1, dirData_AddOn2, filesep];
		% Add elements for voxel location and quantification analysis to input directory
		% name
		dirData				= [dirData_Processed, strVOI, '_LCModel_Data/'];
		
		% Lists of filenames for MRS spectra and (unsuppresed) water signals depending on 
		% type of data (spectra) used for analysis
		textFileName_MRS 		= 'list_filenames_MRS_Spectra.txt';
		textFileName_ref_ECC 	= 'list_filenames_MRS_Ref_ECC.txt';
		textFileName_ref_Quant 	= 'list_filenames_MRS_Ref_Quant.txt';
		textFileName_w 			= 'list_filenames_MRS_Water.txt';
		fullFilename_listOfFiles_MRS		= [dirData, textFileName_MRS];
		fullFilename_listOfFiles_ref_ECC	= [dirData, textFileName_ref_ECC];
		fullFilename_listOfFiles_ref_Quant	= [dirData, textFileName_ref_Quant];
		fullFilename_listOfFiles_w			= [dirData, textFileName_w];
		
		% Select output directory based on type of data and type of analysis being used
		outDir		= [dirData_Processed, strVOI, '_LCM_Out_', strTissue, strWaterQuant];
		if b0nratio
			outDir		= [outDir, '_0nratio'];
		end
		outDir		= [outDir, filesep];
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS);
end

% Check whether directories for data or results exist
% If directory for data does not exist or is empty, return with error
% if numel(dir(dirData)) <= 2, then folder is empty
if (not(isfolder(dirData))) && (numel(dir(dirData)) <= 2)
    error('%s: Data directory %s does NOT exist or is empty!\n', sFunctionName, dirData);
end

% If directory for results does not exist, create it
if not(isfolder(outDir))
	sMsg = sprintf('%s: Creating output directory %s ...\n', sFunctionName, outDir);
    disp(sMsg);
    if ~mkdir(outDir)
		error('%s: Could not create (mkdir) output directory %s!\n', sFunctionName, outDir);
	end
end



%% Select basis set and control file for LCModel analysis depending on sequence type
switch seqType_MRS
	case 'SPECIAL'
		% For 3T Potsdam pain study and SPECIAL
		dirBasis						= '/home/mekler/.lcmodel/basis-sets/Basis_Sets_Ralf/SPECIAL_SE/3TBasis_new_withAcquired_MM_Verio/';
		LCM_Basis						= '3T_sim_TE8-5_mac_ac.basis';
		LCM_Control						= '3T_RAW_mac_shortTE_ACC_water_nratio0';
		%LCM_Control						= '3T_RAW_mac_shortTE_ACC_water';
		
	case 'MEGA-PRESS'
		% For 3T BCAN Dopamin study and MEGA-PRESS depending on type of data (spectra) 
		% used for analysis
		if strcmp(strAnalysisData, 'MRS_diff')
			% MRS_diff
			dirBasis						= '/home/mekler/.lcmodel/basis-sets/Basis_Sets_Ralf/MEGA-PRESS/3TBasis_PurdueU/';
			LCM_Basis						= '3t_IU_MP_te68_diff_yesNAAG_noLac_Kaiser.basis';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_JM_Method4';
			LCM_Control						= '3T_RAW_MEGA-PRESS_JM_Method4_noECC';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_JM_Method4_noECC_nratio0';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_JM_Method4_noECC_noWScale';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_JM_Method4_noECC_noWScale_nratio0';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_JM_Method2';
		else
			% MRS_editOFF
			dirBasis						= '/home/mekler/.lcmodel/basis-sets/Basis_Sets_Ralf/MEGA-PRESS/3TBasis_PurdueU/';
			LCM_Basis						= '3t_IU_MP_te68_748_ppm_inv_Edit-Off.basis';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_editOFF_TE68_GM_water_nratio0_noECC';
			LCM_Control						= '3T_RAW_MEGA-PRESS_editOFF_TE68_GM_water_nratio0';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_editOFF_TE68_GM_water';
		end
	case 'sLASER'
		% TE = 23 ms
		dirBasis						= '/home/mekler/.lcmodel/basis-sets/Basis_Sets_DineshKD/sLASER/';
		LCM_Basis						= 'sead_3T_23ms_02Nov2017.BASIS';
		LCM_Control						= '3T_RAW_sLASER_TE23_GM_water_nratio0';
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS);
end
% Sequence independent settings
dirControl						= '/home/mekler/.lcmodel/profiles/ralf/control-defaults/';
fullFileName_LCM_Basis			= [dirBasis, LCM_Basis];
fullFileName_LCM_Control		= [dirControl, LCM_Control];
fullFileName_LCM_Control_New	= [outDir, LCM_Control, '_template.control'];
fullFileName_LCM_Control_case	= [outDir, LCM_Control, '_case.control'];


%% Select additional parameters for LCModel analysis or processing options 
% Indicate whether water scaling is used
% (in later version, this should be determined from loaded control file)
charWaterScaling				= 'Yes';		% 'Yes';	'No';
% Select water signals used for water scaling
if strcmp(charWaterScaling, 'Yes')
	switch strWaterQuant
		case '_ref_Quant'
			fullFilename_listOfFiles_MRS_water = fullFilename_listOfFiles_ref_Quant;
		case '_ref_ECC'
			fullFilename_listOfFiles_MRS_water = fullFilename_listOfFiles_ref_ECC;
		case '_w'
			fullFilename_listOfFiles_MRS_water = fullFilename_listOfFiles_w;
			
		otherwise
			error('%s: ERROR: water quantification option strWaterQuant = %s!', sFunctionName, strWaterQuant);
	end		% End of switch strWaterQuant
end		% End of if strcmp(charWaterScaling, 'Yes')


%% Check whether computer is system with LCModel installed
% Note: LCModel is only installed on virtual server s-csb-mrs01
cmd					= 'uname -n';
[status, cmdout]	= system( cmd );
%disp(cmdout);
cmd					= '';
if contains(cmdout,'s-csb-mrs01')
	%cmd = strcat( '$HOME/.lcmodel/bin/lcmodel < ', controlname );
    sMsg = sprintf('%s: LCModel installed on system %s\n', sFunctionName, cmdout);
    disp(sMsg);
else
	%cmd = strcat( 'ssh s-csb-mrs01 $HOME/.lcmodel/bin/lcmodel < ', controlname );
	%sMsg = sprintf('%s: No calculation of CEST data for ROIs!\n', sFunctionName);
	error('%s: LCModel not installed on system %s\n', sFunctionName, cmdout);
end


%% Create template for control file for running LCModel analysis from the command line
% Start creating new control file for LCModel command line analyis
% Open selected original control file, create header line in template control file, and
% copy contents of original into template/new control file, sot that template for new
% control file contains all the information that remains the same for all cases
fid_control		= fopen(fullFileName_LCM_Control, 'r');
if( fid_control == -1 )
	error('%s: Could not open control file %s!\n', sFunctionName, fullFileName_LCM_Control);
else
	fid_control_new		= fopen(fullFileName_LCM_Control_New, 'w');
	if( fid_control_new == -1 )
		status_control	= fclose(fid_control);
		error('%s: Could not open control file %s!\n', sFunctionName, fullFileName_LCM_Control_New);
	else
		nbytes = fprintf(fid_control_new,'$LCMODL\n');
		
		% Copy contents of old into new control file
		while ~feof(fid_control)
			line		= fgetl(fid_control);
			nbytes		= fprintf(fid_control_new,'%s\n', line);
		end
		% Close original control file
		status_control	= fclose(fid_control);
	end
end

% Specify basis file as input in new control file
% (this file is the same for the set of MR spectra to be analyzed)
charAddition		= strcat('FILBAS=''', fullfile(dirBasis, LCM_Basis), '''');
nbytes				= fprintf(fid_control_new, '%s\n', charAddition);    

% Close new control file
status_control_new	= fclose(fid_control_new);


%% Load list of filenames from text files into corresponding data structures
% Load list of filenames for MR spectra
fileID_1						= fopen(fullFilename_listOfFiles_MRS);
cell_listOfFiles_MRS_Spectra	= textscan(fileID_1, '%s');
status_1						= fclose(fileID_1);

% Load list of filenames for water, if selected
if strcmp(charWaterScaling, 'Yes')
	fileID_2						= fopen(fullFilename_listOfFiles_MRS_water);
	cell_listOfFiles_MRS_Water		= textscan(fileID_2, '%s');
	status_2						= fclose(fileID_2);
	sMsg = sprintf('%s: Water scaling is used ...\n', sFunctionName);
else
	sMsg = sprintf('%s: No water scaling is used ...\n', sFunctionName);
end
disp(sMsg);


%% Perform LCModel analysis for all MR spectra specified in the list of files
% For that, extract a MR spectrum from the corrresponding list of files, 
% and, if selected, also extract water file from corresponding list of files; 
% complete new control file for LCModel command line analysis for extracted file(s),
% and invoke LCModel anlysis via command line
noFiles_MRS			= length(cell_listOfFiles_MRS_Spectra{1, 1});
%disp(['List of MRS Spectra Files', sprintf('\t\t\t\t\t\t'), ' List of MRS Water Files']);
for Ind=1 :1 : noFiles_MRS			% noFiles_MRS	% 2		% 0
	% Obtain filename for MR spectrum with and without extension
	filename_MRS						= cell_listOfFiles_MRS_Spectra{1}{Ind};
	[pathstr, filename_MRS_noExt, ext]	= fileparts(filename_MRS);
	
	% Obtain filename for water signal, if water scaling is selected
	if strcmp(charWaterScaling, 'Yes')
		filename_MRS_w	= cell_listOfFiles_MRS_Water{1}{Ind};
		[pathstr_w, filename_MRS_w_noExt, ext_w]	= fileparts(filename_MRS_w);
		disp([filename_MRS, sprintf('\t\t'), filename_MRS_w]);
	else
		disp(filename_MRS);
	end
	
	% Copy template for new control file and complete copy with case-specific information
	% to generate a case-specific control file for each MR spectrum (case)
	[status, sMsg] = copyfile(fullFileName_LCM_Control_New, fullFileName_LCM_Control_case);
	if( status == 0 )
		disp(sMsg);
		error('%s: Could not copy control file %s!\n', sFunctionName, fullFileName_LCM_Control_New);
	end
	
	% Include .table, .ps, .csv, .coord, .print, and .coraw files as LCModel output and
	% append corresponding case-specific filenames to case-specific control file
	fid_control_case		= fopen(fullFileName_LCM_Control_case, 'a');
	if( fid_control_case == -1 )
		error('%s: Could not open case-specific control file %s!\n', sFunctionName, fullFileName_LCM_Control_case);
	end
	% Use ''' to include ' within final character array and thus within new control file,
	% where this syntax is required by LCModel
	charAddition	= strcat('LTABLE=7, FILTAB=''', fullfile(outDir, strcat(filename_MRS_noExt, '.table')), '''');
	nbytes			= fprintf(fid_control_case, '%s\n', charAddition);    
	charAddition	= strcat('FILPS=''', fullfile(outDir, strcat(filename_MRS_noExt, '.ps')), '''');
	nbytes			= fprintf(fid_control_case, '%s\n', charAddition);  
	charAddition	= strcat('LCSV=11, FILCSV=''', fullfile(outDir, strcat(filename_MRS_noExt, '.csv')), '''');
	nbytes			= fprintf(fid_control_case, '%s\n', charAddition);
	charAddition	= strcat('LCOORD=9, FILCOO=''', fullfile(outDir, strcat(filename_MRS_noExt, '.coord')), '''');
	nbytes			= fprintf(fid_control_case, '%s\n', charAddition); 
	charAddition	= strcat('LPRINT=6, FILPRI=''', fullfile(outDir, strcat(filename_MRS_noExt, '.print')), '''');
	nbytes			= fprintf(fid_control_case, '%s\n', charAddition); 
	charAddition	= strcat('LCORAW=10, FILCOR=''', fullfile(outDir, strcat(filename_MRS_noExt, '.coraw')), '''');
	nbytes			= fprintf(fid_control_case, '%s\n', charAddition); 
	
	% Specify files for water (if selected) and MR spectrum and and add to case-specific
	% control file
	if strcmp(charWaterScaling, 'Yes')
% 		% Replace the '.' in extension '.RAW' with '_' and use new extension '.h2o
% 		charAddition	= strcat('FILH2O=''', fullfile(dirData, strcat(filename_MRS_w_noExt, '_', ext_w(2:end),'.h2o')), '''');
		% Use .RAW extension and original filename of water signal
		charAddition	= strcat('FILH2O=''', fullfile(dirData, filename_MRS_w), '''');
		nbytes			= fprintf(fid_control_case, '%s\n', charAddition);  
	end
	charAddition	= strcat('FILRAW=''', fullfile(dirData, filename_MRS), '''');
	nbytes			= fprintf(fid_control_case, '%s\n', charAddition); 
	
	% Indicate end of case-specific control file for LCModel command line analysis
	nbytes			= fprintf(fid_control_case, '$END\n');
	
	% Close the case-specific control file
	status_control_case	= fclose(fid_control_case);
	
	
	%% Perform LCModel command line analyis using the case-specific control file
	% LCModel command line analysis is realized via a (Linux) system call
	sMsg	= sprintf('\n LCModel analysis ...\n');
	disp(sMsg);
	cmd		= strcat( '$HOME/.lcmodel/bin/lcmodel < ', fullFileName_LCM_Control_case ); 
	[status, cmdout]	= system( cmd );
    disp(['   LCModel command line output: ', cmdout]);
	
	% Convert .ps file from LCModel output into pdf file using corresponding linux routine
	% via a system call
	cmd_pdf						= ['ps2pdf ', fullfile(outDir, strcat(filename_MRS_noExt, '.ps')), ' ', fullfile(outDir, strcat(filename_MRS_noExt, '.pdf'))];
	[status_pdf, cmdout_pdf]	= system( cmd_pdf );

	% Delete case-specific control file, if not last case (last MR spectrum in list)
	% This is an approach to control any possible accumulation of (control) files, so that 
	% a "fresh" new file can be generated for next case while always 
	% using the same filename for the case-specifc control file; alternatively, different
	% filenames could be used for each case-specific control file
	% (Actually, this step is not required, rather it would work without it as well)
	%if( Ind ~= 2)
	if( Ind ~= noFiles_MRS)
		delete(fullFileName_LCM_Control_case);
	end
	disp(sMsg_newLines);
end		% End of for Ind=1 :1 : noFiles_MRS			% noFiles_MRS	% 2		% 0


%% Process output from LCModel analysis
% NOTE: SLIGHTLY INCORRECT FOR MEGA-PRESS FOR NOW!! (Check header lines, etc.)
% Set parameters, allocate variables, obtain list of .table files generated by LCModel 
% analysis, and select types of metabolite concentrations
% NOTE: In LCModel output there are three columns of values: Conc., %SD, and /Cr+PCr
% Here these values are termed "types of metabolite concentrations" and are selected by
% indices 1, 2, or 3 that indicate the position of their respective column
mydata					= {};
noHeaderLines			= 7;
charTmp					= fullfile(outDir, '*.table');
listFiles_table			= dir(charTmp);
%fullFileName_all		= fullfile(outDir, strcat('3T_', seqType_MRS, '_LCModel_', LCM_Basis, '_', LCM_Control, '_', dt, '_all.csv'));
%fullFileName_all		= fullfile(outDir, strcat('3T_', seqType_MRS, '_LCModel_', LCM_Basis, '_', LCM_Control, '_all.csv'));
fullFileName_all		= fullfile(outDir, strcat('3T_', 'LCModel_', LCM_Basis, '_', LCM_Control, '_all.csv'));
noFiles_table			= numel(listFiles_table);
selectedTypesMetabConc	= [1; 2; 3;];
noTypesMetabConc		= numel(selectedTypesMetabConc);

% From SC_GUI.m -> function finish_lcm:
%selectedTypesMetabConc	= handles.batch_selected_types_of_metabolite_concentrations;
%noTypesMetabConc		= numel(handles.batch_selected_types_of_metabolite_concentrations);

% Collect metabolite quantification results from all individual LCModel analyses into one
% new .csv file and compute mean value and standard deviation for all quantities
if( noFiles_table > 0 )
	
	% Determine # of metabolites from contents of first .table file
	charTmp		= fullfile(outDir, listFiles_table(1).name);
	
	txt			= fileread( charTmp );
	tokens		= strfind( txt, '$$CONC' );
	noMet		= sscanf(txt(tokens+6:tokens+9),'%d') - 1;    % Number of metabolites
	
	% Collect information about metabolites, concentrations, and Cramer-Rao Lower Bounds
	% (CRLBs)/%SD from LCModel output from all .table files into one data structure
	for inr=1 : 1 : noFiles_table
		charTmp		= fullfile(outDir, listFiles_table(inr).name);	
		fid_table	= fopen(charTmp);
		if( fid_table == -1 )
			error('%s: Could not open .table file %s!\n', sFunctionName, charTmp);
		end
		% Skip first few header lines to get to concentration values in .table file
		for lnr=1 : 1 : noHeaderLines
			line = fgetl(fid_table);
		end
		% Write values for concentrations, CRLBs, relative concentration and metabolite
		% names into data structure
		% (The textscan function reapplies formatSpec throughout the entire file and stops
		%  when it cannot match formatSpec to the data.)
		mydata{inr}		= textscan(fid_table, '%f %d%% %f %s');
		
		% Close .table file
		status_table	= fclose(fid_table);
	end

	% Write metabolite quantification results from all .table files to new .csv file in
	% text mode
	fid_all		= fopen(fullFileName_all, 'wt');
	if( fid_all == -1 )
		error('%s: Could not open _all.csv file %s!\n', sFunctionName, fullFileName_all);
	end
	
	% Header line in new .csv file
	nbytes		= fprintf(fid_all, '             \t');
	% Write name of each metabolite as many times as there are selected types of
	% metabolite concentrations into new .csv file as a heading for the subsequent values
	for im=1: 1 : noMet
		for ic=1 : 1 : noTypesMetabConc
			nbytes	= fprintf(fid_all, '  %s   \t',mydata{1}{4}{im});
		end
	end
	nbytes		= fprintf(fid_all, '\n');
	
	% Write concentration values and CRLBs/%SD from all .table files into new .csv file
	% Variable 'faverage is used to be able to select only specific (not all) types of
	% metabolite concentrations to be written to file
	faverage	= zeros(3, noMet, noFiles_table, 'double');
	for inr=1 : 1 : noFiles_table
		% Write filename of .table file into new .csv file
		nbytes		= fprintf(fid_all, '%s  \t', listFiles_table(inr).name);
		for im = 1: 1 : noMet
			for ic = 1: 1 : noTypesMetabConc
				nbytes	= fprintf(fid_all, '  %0.5e   \t', mydata{inr}{selectedTypesMetabConc(ic)}(im));			
				faverage(selectedTypesMetabConc(ic),(im),(inr)) = ...
					mydata{inr}{selectedTypesMetabConc(ic)}((im));
			end
		end
		nbytes		= fprintf(fid_all, '\n');
	end
	
	% Insert blank line into new .csv file using the same tab format as for the metabolite
	% concentrations, which makes it later easy to read in the new .csv file into Excel
	% using the tab as a delimiter
	nbytes		= fprintf(fid_all, '     \t');
	for im=1 : 1 : noMet
		for ic=1 : 1: noTypesMetabConc
			nbytes	= fprintf(fid_all, '  \t');
		end
	end
	nbytes		= fprintf(fid_all, '\n');
	
	% Compute mean values for all quantities and write results to new .csv file
	nbytes		=  fprintf(fid_all, 'Mean     \t');
	for im =1 : 1 : noMet
		for ic=1: 1 : noTypesMetabConc
			nbytes	= fprintf(fid_all, '  %0.5e  \t', mean(faverage(selectedTypesMetabConc(ic),(im),(1:noFiles_table))));
		end
	end
	nbytes		= fprintf(fid_all, '\n');
	
	% Compute standard deviation for all quantities and write results to new .csv file
	nbytes	=	fprintf(fid_all, 'Standard Deviation\t');
	for im=1: 1 : noMet
		for ic=1 : 1 : noTypesMetabConc	
			nbytes	= fprintf(fid_all, '  %0.5e  \t', std(faverage(selectedTypesMetabConc(ic),(im),(1:noFiles_table))));
		end
	end
	nbytes		= fprintf(fid_all, '\n');
	
	% Close new .csv file
	status_all	= fclose(fid_all);	
end		% End of if( noFiles_table > 0 )



%% Save variables of workspace to file
% Obtain current date and time in specific format
dt		= datestr(now,'yyyymmdd_HH_MM_SS');

% Save workspace into output directory (optional with user input)
% (Extension".mat" in filename explicitly required, so that Matlab can correctly load 
% workspace file with a "." in its filename)
%fullFileName_SavedWorkspace		= fullfile(outDir, strcat('3T_LCModel_', LCM_Basis, '_', LCM_Control, '.mat'));
%fullFileName_SavedWorkspace		= fullfile(outDir, strcat('3T_', seqType_MRS, '_LCModel_', LCM_Basis, '_', LCM_Control, '_', dt, '.mat'));
fullFileName_SavedWorkspace		= fullfile(outDir, strcat('3T_', 'LCModel_', LCM_Basis, '_', LCM_Control, '_', dt, '.mat'));
%charVecSaveWorkspace	= input('Would you like to save all variables of the workspace to file?  ', 's');
if strcmp(charVecSaveWorkspace,'y') || strcmp(charVecSaveWorkspace,'Y')
	save(fullFileName_SavedWorkspace);
end

