%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% run_LCModel_ListOfFiles_s.m
%
%% Script to run LCModel quantification on a list of files of magnetic resonance 
%% spectroscopy (MRS) data
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2018, 2019, 2020, 2021, 2022; 
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
charVecSaveWorkspace	= 'y';
strVOI					= 'PCG';			% 'HC';		% 'PCG';
seqType_MRS				= 'MEGA-PRESS';		% 'SPECIAL';	% 'MEGA-PRESS'; % 'sLASER';
dataType_MRS			= 'mrs_w';		% 'mrs';		% 'mrs_w';	% 'mrs_w_ref;		
%										% 'mrs_ref';	% 'water';	% 'water_ref';
% strOVS_In				= 'wOVS';
% strOVS_w_In				= 'wOVS';
noSD_In					= 3.2;			% 3.2;		2.6;		4.0;
digits					= [fix(noSD_In) round(abs(noSD_In-fix(noSD_In))*10)];
str_noSD_In				= sprintf('%d_%d', digits(1), digits(2));
% strMinUserIn_In			= 'y';
% aaDomain_In				= 'f';
% tmaxin_In				= 0.2;
% iterin_In				= 20;
% alignSS_In				= 2;
bECC_In					= 1;
strTissue				= 'PCG';		% 'GM';	% 'WM';	% 'HC';	% 'PCG';
strAnalysisData			= 'MRS_diff';	% 'MRS_diff';	'MRS_editOFF';	'MRS_reg';
b0nratio				= 1;
% Indicate whether water scaling is used
% (in later version, this should be determined from loaded control file)
charWaterScaling		= 'Yes';		% 'Yes';	'No';
strWaterQuant			= '_w';			% '_ref_Quant'; % '_ref_ECC';	% '_w'; %'';
% Set parameters for copying results from .csv file into formatted Excel file
bCopyIntoExcel			= 1;



%% Select basis set and control file for LCModel analysis depending on sequence type
% Init directory for control files and string(s) for output directory addition(s)
dirControl						= '/home/mekler/.lcmodel/profiles/ralf/control-defaults/';
outDir_AddControl				= '';
switch seqType_MRS
	case 'SPECIAL'
		% For 3T Potsdam pain study and SPECIAL
		dirBasis						= '/home/mekler/.lcmodel/basis-sets/Basis_Sets_Ralf/Basis_SPECIAL_SE/3TBasis_new_withAcquired_MM_Verio/';
		LCM_Basis						= '3T_sim_TE8-5_mac_ac.basis';
		dirControl						= '/home/mekler/.lcmodel/profiles/ralf/control-defaults/Other/';
		LCM_Control						= '3T_RAW_mac_shortTE_ACC_water_nratio0';
		%LCM_Control						= '3T_RAW_mac_shortTE_ACC_water';	
	case 'MEGA-PRESS'
		% Control files for MEGA-PRESS
		dirControl						= '/home/mekler/.lcmodel/profiles/ralf/control-defaults/LCM_Control_MEGA-PRESS/';
		% For 3T BCAN Dopamin study and MEGA-PRESS depending on type of data (spectra) 
		% used for analysis
		if strcmp(strAnalysisData, 'MRS_diff')
			% MRS_diff
			% Basis sets and control files for water-symmetric editing
			dirBasis						= '/home/mekler/.lcmodel/basis-sets/Basis_Sets_Ralf/Basis_MEGA-PRESS/3TBasis_PurdueU/';
			LCM_Basis						= '3t_IU_MP_te68_diff_yesNAAG_noLac_Kaiser.basis';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_JM_Method4_noECC';
			LCM_Control						= '3T_RAW_MEGA-PRESS_JM_Method4_noECC_nratio0';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_JM_Method4_noECC_noWScale';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_JM_Method4_noECC_noWScale_nratio0';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_JM_Method2';
			
			% Basis sets and control files for MM-suppressed (symmetric) editing
			%dirBasis						= '/home/mekler/.lcmodel/basis-sets/Basis_Sets_Ralf/Basis_MEGA-PRESS/3TBasis_JimMurdoch/';
			%LCM_Basis						= '3t_IU_MEGAPRESS_1915_te68_Kaiser_diff.basis';
			%LCM_Basis						= '3t_IU_MEGAPRESS_1915_te68_Kaiser_diff.basis';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_JM_Config1_noECC';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_JM_Config1_noECC_nratio0';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_MM-Symm_Diff_Config1_1_85-4_1ppm';
			
			% Set pattern to be searched for in filename of control file to shorten final
			% name of output directory
			strPatt							= 'Diff';
		else
			% MRS_editOFF
			% Basis sets and control files for water-symmetric editing
			dirBasis						= '/home/mekler/.lcmodel/basis-sets/Basis_Sets_Ralf/Basis_MEGA-PRESS/3TBasis_PurdueU/';
			LCM_Basis						= '3t_IU_MP_te68_748_ppm_inv_Edit-Off.basis';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_editOFF_TE68_GM_water_nratio0_noECC';
			LCM_Control						= '3T_RAW_MEGA-PRESS_editOFF_TE68_GM_water_nratio0';
			%LCM_Control						= '3T_RAW_MEGA-PRESS_editOFF_TE68_GM_water';
			
			% Set pattern to be searched for in filename of control file to shorten final
			% name of output directory
			strPatt							= 'editOFF';
		end
		
		% Create addition to output directory name according to (parts of) the
		% filename of the control file
		% (To keep resulting output directory paths short, only part of the filename
		% of the control file is at times included in the output directory name)
		len_LCM_Control					= length(LCM_Control);
		IndPatt							= strfind(LCM_Control, strPatt);
		if ~isempty(IndPatt)
			outDir_AddControl			= ['DOPA_LCM_', LCM_Control(IndPatt:len_LCM_Control)];
		else
			% String pattern not found in filename of control file, so add entire
			% filename
			outDir_AddControl			= ['DOPA_LCM_', LCM_Control];
		end
	case 'sLASER'
		% TE = 23 ms
		dirBasis						= '/home/mekler/.lcmodel/basis-sets/Basis_Sets_DineshKD/Basis_sLASER/';
		LCM_Basis						= 'sead_3T_23ms_02Nov2017.BASIS';
		dirControl						= '/home/mekler/.lcmodel/profiles/ralf/control-defaults/LCM_Control_sLASER_dkd_TE23/';
		switch strTissue
			case 'GM'
				LCM_Control						= '3T_RAW_sLASER_TE23_GM_water_nratio0';
			case 'HC'
				%LCM_Control						= '3T_RAW4093_sLASER_TE23_HC_water_nratio0_noECC_43772';
				%LCM_Control						= '3T_RAW4094_sLASER_TE23_HC_water_nratio0_noECC_43772';
				%LCM_Control						= '3T_RAW4094_sLASER_TE23_HC_water_nratio0_noECC_40592';
				LCM_Control						= '3T_RAW4094_sLASER_TE23_HC_water_nratio0_noECC_40592_mac';
			case 'PCG'
				%LCM_Control						= '3T_RAW4093_sLASER_TE23_PCG_water_nratio0_noECC_45322';
				%LCM_Control						= '3T_RAW4094_sLASER_TE23_PCG_water_nratio0_noECC_45322';
				%LCM_Control						= '3T_RAW4094_sLASER_TE23_PCG_water_nratio0_noECC_42708';
				LCM_Control						= '3T_RAW4094_sLASER_TE23_PCG_water_nratio0_noECC_42708_mac';
				%LCM_Control						= '3T_RAW4094_sLASER_TE23_PCG_water_noECC_42708';
				%LCM_Control						= '3T_RAW4094_sLASER_TE23_HC_water_nratio0_noECC_40592';
				
				
			otherwise
				error('%s: ERROR: No LCM control file found for strTissue = %s!', sFunctionName, strTissue);
		end		% End of switch strTissue
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS);
end		% End of switch seqType_MRS

% Sequence independent settings
%dirControl						= '/home/mekler/.lcmodel/profiles/ralf/control-defaults/';
fullFileName_LCM_Basis			= [dirBasis, LCM_Basis];
fullFileName_LCM_Control		= [dirControl, LCM_Control];
%fullFileName_LCM_Control_New	= [outDir, LCM_Control, '_template.control'];
%fullFileName_LCM_Control_case	= [outDir, LCM_Control, '_case.control'];


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
		fullFilename_listOfFiles_w		= [dirData, 'list_filenames_MRS_Water.txt'];
		%fullFilename_listOfFiles_MRS_w	= [dirData, 'list_filenames_MRS_Water.txt'];	
	case 'MEGA-PRESS'
		% Select directories for (input) data files and for output data depending on # of 
		% SDs used for pre-processing of MR spectra and on analysis settings
		dirDataAnalysis					= ['/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/'];
		dirData							= [dirDataAnalysis 'Z_DOPA_FID-A_SD_', str_noSD_In, filesep, 'DOPA_LCModel_Analysis_Data/'];
		outDir_Base						= [dirDataAnalysis 'Z_DOPA_FID-A_SD_', str_noSD_In, filesep];
		if isempty(outDir_AddControl)
			error('%s: ERROR: Addition to output directory name outDir_AddControl =  %s is empty!', sFunctionName, outDir_AddControl);
		end
		outDir							= [outDir_Base, outDir_AddControl, filesep];

		% Lists of filenames for MRS spectra and water signals depending on type of data 
		% (spectra) used for analysis
		if strcmp(strAnalysisData, 'MRS_diff')
			% MRS_diff
			fullFilename_listOfFiles_MRS	= [dirData, 'list_filenames_MRS_Diff_Spectra.txt'];
			fullFilename_listOfFiles_w		= [dirData, 'list_filenames_MRS_editOFF_Water.txt'];
			%fullFilename_listOfFiles_MRS_w	= [dirData, 'list_filenames_MRS_editOFF_Water.txt'];
		else
			% MRS_editOFF
			fullFilename_listOfFiles_MRS	= [dirData, 'list_filenames_MRS_editOFF_Spectra.txt'];
			fullFilename_listOfFiles_w		= [dirData, 'list_filenames_MRS_editOFF_Water.txt'];
			%fullFilename_listOfFiles_MRS_w	= [dirData, 'list_filenames_MRS_editOFF_Water.txt'];
		end		% End of if strcmp(strAnalysisData, 'MRS_diff')
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


%% Select additional parameters for LCModel analysis or processing options 
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
fullFileName_LCM_Control_New	= [outDir, LCM_Control, '_template.control'];
fullFileName_LCM_Control_case	= [outDir, LCM_Control, '_case.control'];
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
for Ind=1 : 1 : noFiles_MRS			% noFiles_MRS	% 2		% 0
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
%fullFileName_all		= fullfile(outDir, strcat('3T_', 'LCModel_', LCM_Basis, '_', LCM_Control, '_all.csv'));
FileName_all			= strcat(str_noSD_In, '_LCModel_', LCM_Basis, '_', LCM_Control);
fullFileName_all		= fullfile(outDir, [FileName_all '.csv']);
noFiles_table			= numel(listFiles_table);
selectedTypesMetabConc	= [1; 2; 3;];
noTypesMetabConc		= numel(selectedTypesMetabConc);
noMet					= 0;

% From SC_GUI.m -> function finish_lcm:
%selectedTypesMetabConc	= handles.batch_selected_types_of_metabolite_concentrations;
%noTypesMetabConc		= numel(handles.batch_selected_types_of_metabolite_concentrations);

% Collect metabolite quantification results from all individual LCModel analyses into one
% new .csv file and compute mean value and standard deviation for all quantities
if( noFiles_table > 0 )
	
	% Determine # of metabolites from contents of first .table file
	charTmp		= fullfile(outDir, listFiles_table(1).name);
	
	txt			= fileread( charTmp );
	
	% Check whether .table file contains any information about concentrations of 
	% metabolites
	% If LCModel fit fails, which is rare, but can happen, e.g. if data are corrupted and
	% restrictions on the degrees of freedom for fitting, e.g. by setting the parameter
	% dknmtn to a specific value, are imposed, no metabolite concentrations are reported
	% in the corresponding .table file; the index 'tokens' is then empty
	tokens		= strfind( txt, '$$CONC' );
	if isempty(tokens)
		% No metabolite concentrations found, exit with error
		error('%s: tokens = %d (empty) => No metabolite concentrations found in %s!\n Check on other .table files!\n', sFunctionName, tokens, charTmp);
	else
		noMet		= sscanf(txt(tokens+6:tokens+9),'%d') - 1;    % Number of metabolites
	end
	
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
	% Use term for case/subject to facilitate reading in of .csv file as a Matlab table
	%nbytes		= fprintf(fid_all, '             \t');
	nbytes		= fprintf(fid_all, '  Case   \t');
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
		
		% If concentration values exist, write them into new .csv file; if not, i.e. if 
		% fit for metabolite quantification failed for some reason,
		% insert zeros into new .csv file and into array of average values, so that 
		% routine continues to write information from subsequent .table files to file and 
		% does not issue an error
		% Concentration values exist, if last (fourth) cell array of data structure is 
		% non-empty and if correct # of metabolites that was determined from # of 
		% metabolites in first .table file was read into first cell array of data
		% structure
		listMet_tableFile	= mydata{1,inr}{1,4};
		noMet_tableFile		= length(mydata{1,inr}{1,1});
		if ~isempty(listMet_tableFile) && noMet_tableFile == noMet
			% Write concentration and other values into new .csv file
			for im = 1: 1 : noMet
				for ic = 1: 1 : noTypesMetabConc
					nbytes	= fprintf(fid_all, '  %0.5e   \t', mydata{inr}{selectedTypesMetabConc(ic)}(im));
					faverage(selectedTypesMetabConc(ic),(im),(inr)) = ...
						mydata{inr}{selectedTypesMetabConc(ic)}((im));
				end
			end
		else	% Concentration information not found or incomplete (fit probably failed)
			% Issue warning and fFor all values, insert zeros into new .csv files and 
			% into array for averaging
			warning('%s: .table file = %s\n\nEither empty cell array of concentrations isempty(listMet_tableFile) = %d and/or # of metabolite concentrations in this table file = noMet_tableFile = %d   ~=    noMet = %d = # of expected metabolite concentrations!\n\n=> Inserting zeros into file for all values\n\n', ...
				sFunctionName, listFiles_table(inr).name, isempty(listMet_tableFile), noMet_tableFile, noMet);
			for im = 1: 1 : noMet
				for ic = 1: 1 : noTypesMetabConc
					nbytes	= fprintf(fid_all, '  %0.5e   \t', 0);
					faverage(selectedTypesMetabConc(ic),(im),(inr)) = 0;
				end
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
	% Use only one word for standard deviation to facilitate reading in of .csv file as a
	% Matlab table
	%nbytes	=	fprintf(fid_all, 'Standard Deviation\t');
	nbytes	=	fprintf(fid_all, 'STD\t');
	for im=1: 1 : noMet
		for ic=1 : 1 : noTypesMetabConc	
			nbytes	= fprintf(fid_all, '  %0.5e  \t', std(faverage(selectedTypesMetabConc(ic),(im),(1:noFiles_table))));
		end
	end
	nbytes		= fprintf(fid_all, '\n');
	
	% Close new .csv file
	status_all	= fclose(fid_all);
	
	% Set parameters for modifying table of results from metabolite quantification and
	% saving/copying final table of results into (formatted) Excel file
	switch seqType_MRS
		case 'SPECIAL'
			disp('NOT YET COMPLETED ...');
			bCopyIntoExcel			= 0;
		case 'MEGA-PRESS'
			% Insert selected # of columns for additional parameters, e.g. water linewidth (LW_H2O),
			% into table of results
			% Create two variables of empty strings and zeros with fitting variable names and
			% same # of rows as table of results; create table from these two variables,
			% use # of colums to be moved for reordering of columns within table,
			% concatenate this table with table of results; then reorder columns in resulting
			% table to obtain desired order of columns
			L_DOPA					= strings(noFiles_MRS, 1);
			LW_H2O					= zeros(noFiles_MRS, 1);
			newColsTable_Add		= table(L_DOPA, LW_H2O);
			nLastRowsToDel			= 3;
			nColsToMove				= 1;
			bUseTemplateFile		= 1;
			fullPath_TemplateFile	= [dirDataAnalysis '3T_MRS_DOPA_Analysis_Template.xltx'];
			
			% Select sheet and range of Excel file, where table of results is to be saved
			% depending on type of spectra being analyzed
			if strcmp(strAnalysisData, 'MRS_diff')
				% MRS_diff
				strSheetSel				= 'MRS_Diff_All';
				strRangeSel				= 'A4';
			else
				% MRS_editOFF
				trSheetSel				= 'MRS_editOFF_All';
				strRangeSel				= 'A4';
			end		% End of if strcmp(strAnalysisData, 'MRS_diff')		
		case 'sLASER'
			disp('Coming soon ...');
			bCopyIntoExcel			= 0;
			
		otherwise
			error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS);
	end

	% Save/copy results from .csv file also into (formatted) Excel file, if selected
	if bCopyIntoExcel
		% Insert selected # of columns for additional parameters, e.g. water linewidth (LW_H2O),
		% into table of results, reorder final table of results and save it into Excel
		% file; Excel file will be formatted, if Excel template file is used
		[resultsTable_Out] = copyMetabResults_IntoExcel_s(fullFileName_all,nLastRowsToDel,newColsTable_Add,nColsToMove,outDir,strSheetSel,strRangeSel,bUseTemplateFile,fullPath_TemplateFile);
	end		% End of if bCopyIntoExcel
	
end		% End of if( noFiles_table > 0 )



%% Save variables of workspace to file
% Obtain current date and time in specific format
dt		= datestr(now,'yyyymmdd_HH_MM_SS');

% Save workspace into output directory (optional with user input)
% (Extension".mat" in filename explicitly required, so that Matlab can correctly load 
% workspace file with a "." in its filename)
%fullFileName_SavedWorkspace		= fullfile(outDir, strcat('3T_LCModel_', LCM_Basis, '_', LCM_Control, '.mat'));
%fullFileName_SavedWorkspace		= fullfile(outDir, strcat('3T_', seqType_MRS, '_LCModel_', LCM_Basis, '_', LCM_Control, '_', dt, '.mat'));
%fullFileName_SavedWorkspace		= fullfile(outDir, strcat('3T_', 'LCModel_', LCM_Basis, '_', LCM_Control, '_', dt, '.mat'));
fullFileName_SavedWorkspace		= fullfile(outDir, [FileName_all '_' dt '.mat']);
%charVecSaveWorkspace	= input('Would you like to save all variables of the workspace to file?  ', 's');
if strcmp(charVecSaveWorkspace,'y') || strcmp(charVecSaveWorkspace,'Y')
	save(fullFileName_SavedWorkspace);
end




%% Code no longer used
% 		switch(noSD_In)
% 			case(2.6)
% 				dirData							= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopamin_Analysis/Z_DOPA_FID-A_SD_2_6/DOPA_LCModel_Analysis_Data/';
% 				
% 				% Select output directory depending on type of data (spectra) used for
% 				% analysis
% 				if strcmp(strAnalysisData, 'MRS_diff')
% 					% MRS_diff
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_2_6/DOPA_LCM_Out_noECC/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_2_6/DOPA_LCM_Out_noECC_0nratio/';
% 					outDir							= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_2_6/DOPA_LCM_Out_noECC/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_2_6/DOPA_LCM_Out_noECC_0nratio/';
% 				else
% 					% MRS_editOFF
% 					outDir							= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_2_6/DOPA_LCM_Out_editOFF_GM_w_0nratio/';
% 				end
% 			case(3.2)
% 				dirData							= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_3_2/DOPA_LCModel_Analysis_Data/';
% 				%dirData							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2_New/DOPA_LCModel_Analysis_Data/';
% 				%dirData							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_3_2/DOPA_LCModel_Analysis_Data/';
% 				
% 				% Select output directory depending on type of data (spectra) used for
% 				% analysis
% 				if strcmp(strAnalysisData, 'MRS_diff')
% 					% MRS_diff
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2/DOPA_LCM_Out_noECC/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2/DOPA_LCM_Out_noECC_0nratio/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2_New/DOPA_LCM_Out_noECC/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2_New/DOPA_LCM_Out_noECC_0nratio/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_forLCModel_SD_3_2/DOPA_LCModel_Results_noECC/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_3_2/DOPA_LCModel_Results_noECC/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_3_2/DOPA_LCModel_Results_noECC_0nratio/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_3_2/DOPA_LCModel_Results_noECC_noWS/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopamin_Analysis/Z_DOPA_FID-A_SD_3_2/DOPA_LCM_Out_noECC/';
% 					outDir							= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_3_2/DOPA_LCM_Out_noECC_0nratio/';
% 				else
% 					% MRS_editOFF
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2/DOPA_LCM_Out_editOFF_GM_w_0nratio_noECC/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2/DOPA_LCM_Out_editOFF_GM_w_0nratio/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_3_2/DOPA_LCM_Out_editOFF_GM_w/';
% 					outDir							= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_3_2/DOPA_LCM_Out_editOFF_GM_w_0nratio/';
% 				end
% 			case(4.0)
% 				dirData							= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_4_0/DOPA_LCModel_Analysis_Data/';
% 				%dirData							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_4_0/DOPA_LCModel_Analysis_Data/';
% 				%dirData							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_4_0/DOPA_LCModel_Analysis_Data/';
% 
% 				% Select output directory depending on type of data (spectra) used for
% 				% analysis
% 				if strcmp(strAnalysisData, 'MRS_diff')
% 					% MRS_diff
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_4_0/DOPA_LCM_Out_noECC/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_Z_Analysis/DOPA_FID-A_SD_4_0/DOPA_LCM_Out_noECC_0nratio/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_4_0/DOPA_LCModel_Results/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_4_0/DOPA_LCModel_Results_noECC/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/DOPA_FID-A_forLCModel_SD_4_0/DOPA_LCModel_Results_noECC_noWS/';
% 					outDir							= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_4_0/DOPA_LCM_Out_noECC/';
% 					%outDir							= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_4_0/DOPA_LCM_Out_noECC_0nratio/';
% 				else
% 					% MRS_editOFF
% 					outDir							= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_4_0/DOPA_LCM_Out_editOFF_GM_w_0nratio/';
% 				end
% 				
% 			otherwise
% 				error('%s: ERROR: No directory/data for noSD_In =  %f!', sFunctionName, noSD_In);
% 		end
