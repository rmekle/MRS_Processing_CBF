%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% loadlcmdetail_ListOfFiles_s.m
%
%% Script to load detailed output from LCModel metabolite quantification for list of files
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% clearvars;
% close all;



%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'loadlcmdetail_ListOfFiles_s';
fprintf('\n\n');


%% Init input parameters for loading detailed output from LCModel metabolite quantification for list of files
% Obtain parameter settings for preprocessing from initialization routine
%configSel				= 'config_Study_sLASER_VOI_IMA_MRS_lsN_SDx_y_SR1_ECC';
configSel				= 'config_Study_sLASER_VOI_dat_MRS_lsN_SDx_y_SR1_ECC';
[sParamsMRS_struct]		= initParams_MRS_s(configSel);

% Extract parameter settings from parameter struct
filename_In				= sParamsMRS_struct.filename;
filename_w_In			= sParamsMRS_struct.filename_w;
strStudy_MRS			= sParamsMRS_struct.strStudy_MRS;
seqType_MRS				= sParamsMRS_struct.seqType_MRS;
strVOI_MRS				= sParamsMRS_struct.strVOI_MRS;
fileExt_MRS				= sParamsMRS_struct.fileExt_MRS;
dataType_MRS			= sParamsMRS_struct.dataType_MRS;
signals_MRS				= sParamsMRS_struct.signals_MRS;
strOVS_In				= sParamsMRS_struct.strOVS;
strOVS_w_In				= sParamsMRS_struct.strOVS_w;
leftshift_In			= sParamsMRS_struct.leftshift;
avgBlockSize_In			= sParamsMRS_struct.avgBlockSize;

% Info about processing tool(s) mainly used
strProcessTool_In		= sParamsMRS_struct.strProcessTool;

% Parameters for removal of bad averages
rmbadav_In				= sParamsMRS_struct.rmbadav;
noSD_In					= sParamsMRS_struct.noSD;
%digits_noSD_In			= [fix(noSD_In) round(abs(noSD_In-fix(noSD_In))*10)];

% Parameters for spectral registration (aligning of averages/frequency and phase drift
% correction) performed in either frequency or time domain
strSpecReg_In			= sParamsMRS_struct.strSpecReg;	% To distinguish settings for spectral registration
driftCorr_In			= sParamsMRS_struct.driftCorr;
iterin_In				= sParamsMRS_struct.iterin;
aaDomain_In				= sParamsMRS_struct.aaDomain;
tmaxin_In				= sParamsMRS_struct.tmaxin;
bTmaxset_In				= sParamsMRS_struct.bTmaxset;
ppmOption				= sParamsMRS_struct.ppmOption;
medin_In				= sParamsMRS_struct.medin;
alignSS_In				= sParamsMRS_struct.alignSS;	% For aligning subspectra (e.g. in SPECIAL)
% Obtain parameters for drift correction depending on type of data, i.e. whether MRS
% data is spectrum or water signal
% NOTE: Check whether aligning of averages in frequency domain works, if the MR
% spectrum is water signal itself; if not, simply align averages in time domain
ppmmin_fix_In			= sParamsMRS_struct.ppmmin_fix;
ppmmaxarray_fix_In		= sParamsMRS_struct.ppmmaxarray_fix;


% Additional parameter settings
bECC_In					= sParamsMRS_struct.bECC;
bPhaseCorrFreqShift_In	= sParamsMRS_struct.bPhaseCorrFreqShift;
strMinUserIn_In			= sParamsMRS_struct.strMinUserIn;
plotSwitch_In			= sParamsMRS_struct.plotSwitch;
reportSwitch_In			= sParamsMRS_struct.reportSwitch;
bPrep_MetabQuant		= sParamsMRS_struct.bPrep_MetabQuant;


%% Additional parameter settings
bSaveResults			= 0;


%% Select input and output directories and filename options
% outDirString_AddOn_1	= '';
% %outDirString_In			= [dirString_In, outDirString_AddOn_1];
% outputFileName_Add_1	= '';
% outputFileName_Add_2	= '';
% 
% % Select directory and other options depending on sequence type and study
% % Use base directory and directory AddOns (e.g. subfolder names) to allow flexible choice
% % of output filename, if results are saved to file
% switch seqType_MRS_In
% 	case 'sLASER'
% 		switch strStudy
% 			case 'Test'
% 				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/';
% 				dirString_In_AddOn_1	= 'Z_Test_Data';
% 				dirString_In_AddOn_2	= 'Test_MRS_SNR_LW';
% 				outNamingOption			= 2;
% 			case '3T_Trauma'
% 				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Trauma_Analysis/';
% 				switch fileExtension
% 					case 'dat'
% 						dirString_In_AddOn_1	= [strVOI, '_', fileExtension, '_FID-A_SD_3_2_ECCref_ls3_SR1'];
% 					case 'IMA'
% 						dirString_In_AddOn_1	= [strVOI, '_', fileExtension, '_FID-A_SD_3_2_ECCref_ls1_SR1'];
% 
% 					otherwise
% 						error('%s: ERROR: Unknown file extension (data type) %s!', sFunctionName, fileExtension);
% 				end			% End of switch fileExtension
% 				dirString_In_AddOn_2	= [strVOI, '_LCModel_Data_MRS_only'];
% 				outNamingOption			= 1;
% 			case '3T_TGA'
% 				dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_MRS_TGA';
% 				switch fileExtension
% 					case 'dat'
% 						%dirString_In_AddOn_1	= [strVOI, '_', fileExtension, '_FID-A_SD_3_2_ECCref_ls3_SR1'];
% 						% .dat svs_sLaser_dkd_LW scans
% 						dirString_In_AddOn_1	= ['MRS_TGA_00_All_RawData_dat_Files_LW_HC'];
% 						dirString_In_AddOn_2	= '';
% 						outNamingOption			= 1;
% 						outDirString_AddOn_1	= [dirString_In_AddOn_1, '_SNR_LW'];
% 						% Init MRS parameters for specific signals chosen for SNR and FWHM
% 						% measurements that are needed to process these signals
% 						[procParams_In]			= initParams_MRS_s('config_Study_sLASER_VOI_dat_water_lsN_SDx_y_SR1_ECC');
% 					case 'IMA'
% 						%dirString_In_AddOn_1	= [strVOI, '_', fileExtension, '_FID-A_SD_3_2_ECCref_ls1_SR1'];
% 						% .IMA svs_sLaser_dkd_LW scans
% 						dirString_In_AddOn_1	= ['MRS_TGA_00_All_DICOM_IMA_Files_LW_HC'];
% 						dirString_In_AddOn_2	= '';
% 						outNamingOption			= 1;
% 						outDirString_AddOn_1	= [dirString_In_AddOn_1, '_SNR_LW'];
% 
% 					otherwise
% 						error('%s: ERROR: Unknown file extension (data type) %s!', sFunctionName, fileExtension);
% 				end			% End of switch fileExtension
% 				%outNamingOption			= 1;
% 
% 			otherwise
% 				error('%s: ERROR: Unknown study %s!', sFunctionName, strStudy);
% 		end				% End of switch strStudy
% 
% 	otherwise
% 		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS_In);
% end		% End of switch seqType_MRS_In
% % Complete name of input and output directory
% %dirString_In			= [dirString_In_Base, dirString_In_AddOn_1, filesep, dirString_In_AddOn_2, filesep];
% dirString_In			= fullfile(dirString_In_Base, dirString_In_AddOn_1, dirString_In_AddOn_2);
% %outDirString_In			= [dirString_In, outDirString_AddOn_1];
% outDirString_In			= fullfile(dirString_In, outDirString_AddOn_1);

% Select input directory based on sequence, study, VOI (voxel), and other information
switch seqType_MRS
	case 'sLASER'
		switch strStudy_MRS
			case '3T_Trauma'
				switch strVOI_MRS
					case 'HC'
						switch fileExt_MRS
							case 'dat'
								dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Trauma_Analysis/HC_dat_FID-A_SD_3_2_ECCref_ls3_SR1/HC_LCM_Out_PCG_ref_Quant_Con8/LCM_print_Sel/';
							case 'IMA'
								dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Trauma_Analysis/HC_IMA_FID-A_SD_3_2_ECCref_ls1_SR1/HC_LCM_Out_HC_ref_Quant_Con8/LCM_print_Sel/';

							otherwise
								error('%s: ERROR: Unknown fileExt_MRS %s!', sFunctionName, fileExt_MRS);
						end			% End of switch fileExt_MRS
					case 'PCG'
						switch fileExt_MRS
							case 'dat'
								dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Trauma_Analysis/PCG_dat_FID-A_SD_3_2_ECCref_ls3_SR1/PCG_LCM_Out_PCG_ref_Quant_Con8/LCM_print_Sel/';
							case 'IMA'
								dirString_In			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Trauma_Analysis/PCG_IMA_FID-A_SD_3_2_ECCref_ls1_SR1/PCG_LCM_Out_PCG_ref_Quant_Con8/LCM_print_Sel/';

							otherwise
								error('%s: ERROR: Unknown fileExt_MRS %s!', sFunctionName, fileExt_MRS);
						end			% End of switch fileExt_MRS

					otherwise
						error('%s: ERROR: Unknown VOI %s for study %s for analyzing detailed LCM ooutput (correlation)!\n', sFunctionName, strVOI_MRS, strStudy_MRS);
				end			% End of switch strVOI_MRS
			case '3T_SBAM'
				error('%s: Not yet for seqType_MRS %s, strStudy_MRS %s, strVOI_MRS = %s, and fileExt_MRS %s!\n\n', seqType_MRS, strStudy_MRS, strVOI_MRS, fileExt_MRS);
			case '3T_TGA'
				error('%s: Not yet for seqType_MRS %s, strStudy_MRS %s, strVOI_MRS = %s, and fileExt_MRS %s!\n\n', seqType_MRS, strStudy_MRS, strVOI_MRS, fileExt_MRS);
			case '7T_KCL'
				error('%s: Not yet for seqType_MRS %s, strStudy_MRS %s, strVOI_MRS = %s, and fileExt_MRS %s!\n\n', seqType_MRS, strStudy_MRS, strVOI_MRS, fileExt_MRS);

			otherwise
				error('%s: ERROR: Unknown study %s!', sFunctionName, strStudy);
		end				% End of switch strStudy
	case 'MEGA-PRESS'
		error('%s: Not yet for seqType_MRS %s, strStudy_MRS %s, strVOI_MRS = %s, and fileExt_MRS %s!\n\n', seqType_MRS, strStudy_MRS, strVOI_MRS, fileExt_MRS);

	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType_MRS);
end		% End of switch seqType_MRS
% Select output directory
outDirString_In			= dirString_In;


%% Obtain information about the list of files for detailed LCModel output
% (assuming that all files are included in the same directory)
% (On Linux, file list in Matlab also includes the two directories "." and "..", which
% means that the actual # of files in the directory is (# of entries in list - 2;
% however, if dir is used to list specific files, e.g. using a file extension, these two 
% directories are not included in the resulting list)
%cd(dirString_In);


% Select file extension to search for depending on selected LCModel output
% (again, here assuming that all LCModel outputfiles are in same directory)
% Then determine filenames and # of files for all corresponding LCModel output files
acSearchString			= '*.print';

% NOTE: File separator at end of directory name is required so that searching for 
%		specific files works properly when using the subsequent notation; multiple file
%		separators at end do not cause an error
%structFileListing		= dir([dirString_In, acSearchString]);
structFileListing		= dir([dirString_In, filesep, acSearchString]);
noEntriesListing		= length( structFileListing );


%% Load detailed LCM output for all MRS data files and extract selected information
% Load detailed LCM output (correlation matrix) for first MRS data file, if existent, to
% obtain size information
if noEntriesListing == 0
	error('%s: No files for acSearchString = %s found in directory %s!\n\n', sFunctionName, acSearchString, dirString_In);
else
	fprintf('%s: %d files for acSearchString = %s found in directory %s!\n\n', sFunctionName, noEntriesListing, acSearchString, dirString_In);
end % End of if noEntriesListing == 0
filename_detailedLCM_In		= structFileListing(1).name;
[metabs, corrMatrix]		= io_loadlcmdetail(fullfile(dirString_In, filename_detailedLCM_In));

% Allocate arrays for detailed outputLCM output (correlation coefficients) for all MRS 
% data files
sz_corrMatrix	= size(corrMatrix);
noMetabs		= length(metabs);
metabs_all		= cell(noMetabs, noEntriesListing);
corrMatrix_all	= zeros([sz_corrMatrix noEntriesListing]);

% Insert detailed LCM output for first MRS data file into corresponding arrays
if ~ismatrix(corrMatrix)
	error('%s: Ndims of corrMatrix of first detailed LCM output file = %d ~= 2!!\n\n', sFunctionName, ndims(corrMatrix));
end		% End of if ndims(corrMatrix) ~= 2
metabs_all(:, 1)			= metabs;
corrMatrix_all(:, :, 1)		= corrMatrix;

% Load detailed LCM output (correlation coefficients) for remaining MRS data files into
% corresponding arrays
indexStart		= 2;
indexStep		= 1;
for ind=indexStart : indexStep : noEntriesListing	% noEntriesListing	% 2  % 1
	filename_detailedLCM_In		= structFileListing(ind).name;
	%fprintf('\n\n');
	fprintf('ind = %3d\t\t%s\n\n', ind, filename_detailedLCM_In);
	[metabs_all(:, ind), corrMatrix_all(:, :, ind)]		= io_loadlcmdetail(fullfile(dirString_In, filename_detailedLCM_In));
end		% End of or ind=indexStart : indexStep : noEntrieAspsListing
fprintf('\n\n');


%% Calculate statistics of aggregate detailed LCM output (correlation coefficients)
nDims_corrMatrix_all		= ndims(corrMatrix_all);
corrMatrix_mean				= mean(corrMatrix_all, nDims_corrMatrix_all);
corrMatrix_std				= std(corrMatrix_all, 0, nDims_corrMatrix_all);
corrMatrix_abs_mean			= mean(abs(corrMatrix_all), nDims_corrMatrix_all);
corrMatrix_max				= max(corrMatrix_all, [], nDims_corrMatrix_all);
corrMatrix_min				= min(corrMatrix_all, [], nDims_corrMatrix_all);
corrMatrix_all_Below1		= corrMatrix_all.*(corrMatrix_all < 1);
corrMatrix_Below1_max		= max(corrMatrix_all_Below1, [], nDims_corrMatrix_all);

% Extract specific inforamtion from all correlation coefficients
% Minimum correlation coefficient
% Maximum correlation coefficient below 1.0
corrMatrix_all_min			= min_array_s(corrMatrix_all);
corrMatrix_all_Below1_max	= max_array_s(corrMatrix_all_Below1);

% Asp and Gln frequently are moderately correlated with other metabolites for short TE MRS
% at 3T
% Asp
metabSel1							= 'Asp';
indMetabSel1						= find(strcmp(metabs, metabSel1));
corrMatrix_mean_metabSel1			= corrMatrix_mean(indMetabSel1, :)
corrMatrix_max_metabSel1			= corrMatrix_max(indMetabSel1, :)
corrMatrix_min_metabSel1			= corrMatrix_min(indMetabSel1, :)
corrMatrix_Below1_max_metabSel1		= corrMatrix_Below1_max(indMetabSel1, :);
[corrMatrix_metabSel1_min, corr_min_metabSel1_ind]			= min(corrMatrix_min_metabSel1);
[corrMatrix_Below1_metabSel1_max, corr_max_metabSel1_ind]	= max(corrMatrix_Below1_max_metabSel1);

% Gln
metabSel2							= 'Gln';
indMetabSel2						= find(strcmp(metabs, metabSel2));
corrMatrix_mean_metabSel2			= corrMatrix_mean(indMetabSel2, :)
corrMatrix_max_metabSel2			= corrMatrix_max(indMetabSel2, :)
corrMatrix_min_metabSel2			= corrMatrix_min(indMetabSel2, :)
corrMatrix_min_metabSel1			= corrMatrix_min(indMetabSel1, :)
corrMatrix_Below1_max_metabSel2		= corrMatrix_Below1_max(indMetabSel2, :);
[corrMatrix_min_metabSel2_min, corr_min_metabSel2_ind]		= min(corrMatrix_min_metabSel2);
[corrMatrix_Below1_metabSel2_max, corr_max_metabSel2_ind]	= max(corrMatrix_Below1_max_metabSel2);

% Display info
fprintf('\n\nCorrelation Coefficients from Detailed LCM Ouput for\n');
fprintf('seqType_MRS %s, strStudy_MRS %s, strVOI_MRS = %s, and fileExt_MRS %s:\n\n', seqType_MRS, strStudy_MRS, strVOI_MRS, fileExt_MRS);
fprintf('corrMatrix_all_min = %.3f\tcorrMatrix_all_Below1_max  = %.3f\n\n', corrMatrix_all_min, corrMatrix_all_Below1_max);


%% Display statistics of all correlation matrices for all metabolites
% fullFilename = '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Trauma_Analysis/PCG_dat_FID-A_SD_3_2_ECCref_ls3_SR1/PCG_LCM_Out_PCG_ref_Quant_Con8/3T_SBA_C_0012_20210401_meas_MID00174_FID96489_svs_slaser_dkd_PCG_TE23_WS128_wOVS_3.2_processed_lcm.print';
% [metabs,corrMatrix]=io_loadlcmdetail(fullFilename);
% figure, image(corrMatrix,'CDataMapping','scaled'), colorbar; set(gca, 'XTick', [1:length(metabs)], 'XTickLabel', metabs, 'YTick', [1:length(metabs)], 'YTickLabel', metabs);
%figure, image(corrMatrix_mean,'CDataMapping','scaled'), colorbar;
%set(gca, 'XTick', [1:length(metabs)], 'XTickLabel', metabs, 'YTick', [1:length(metabs)], 'YTickLabel', metabs);
resolution	= 600;
noFigures	= 5;
h_figs		= gobjects(noFigures, 1);
indFig		= 1;

% Mean of correlation matrices
h_figs(indFig)	= figure; imagesc(corrMatrix_mean), colorbar; 
title(sprintf('%s %s, Correlation Matrix Mean', strStudy_MRS, strVOI_MRS), 'Interpreter', 'none')
set(gca, 'XTick', [1:length(metabs)], 'XTickLabel', metabs, 'YTick', [1:length(metabs)], 'YTickLabel', metabs);
indFig			= indFig+1;

% Standard deviation (std) of correlation matrices
h_figs(indFig)		= figure; imagesc(corrMatrix_std), colorbar;
title(sprintf('%s %s, Correlation Matrix Std', strStudy_MRS, strVOI_MRS), 'Interpreter', 'none')
set(gca, 'XTick', [1:length(metabs)], 'XTickLabel', metabs, 'YTick', [1:length(metabs)], 'YTickLabel', metabs);
indFig			= indFig+1;

% Mean of absolute values of correlation matrices
h_figs(indFig)		= figure; imagesc(corrMatrix_abs_mean), colorbar;
title(sprintf('%s %s, Abs(Correlation Matrix) Mean', strStudy_MRS, strVOI_MRS), 'Interpreter', 'none')
set(gca, 'XTick', [1:length(metabs)], 'XTickLabel', metabs, 'YTick', [1:length(metabs)], 'YTickLabel', metabs);
indFig			= indFig+1;

% Maximum of correlation matrices
h_figs(indFig)		= figure; imagesc(corrMatrix_max), colorbar;
title(sprintf('%s %s, Correlation Matrix Max', strStudy_MRS, strVOI_MRS), 'Interpreter', 'none')
set(gca, 'XTick', [1:length(metabs)], 'XTickLabel', metabs, 'YTick', [1:length(metabs)], 'YTickLabel', metabs);
indFig			= indFig+1;

% Minimum of correlation matrices
h_figs(indFig)		= figure; imagesc(corrMatrix_min), colorbar;
title(sprintf('%s %s, Correlation Matrix Min', strStudy_MRS, strVOI_MRS), 'Interpreter', 'none')
set(gca, 'XTick', [1:length(metabs)], 'XTickLabel', metabs, 'YTick', [1:length(metabs)], 'YTickLabel', metabs);

% Save figures as .fig and .png files, if selected
% Select common output filename for figures and cell array of specific additions for each
% figure
outFileName			= sprintf('%s_%s_%s_LCM_CorrelationCoeffs_', strStudy_MRS, seqType_MRS, strVOI_MRS);
cellFigName_add		= {'Mean', 'Std', 'Mean_Abs', 'Max', 'Min'};
if bSaveResults
	for indFig=1 : 1 : noFigures
		strFigName_add	= cellFigName_add{indFig};
		figureName_fig	= [outFileName, strFigName_add, '.fig'];
		figureName_png	= [outFileName, strFigName_add, '.png'];
		saveFigure_s(h_figs(indFig), outDirString_In, figureName_fig, 'fig', resolution);
		saveFigure_s(h_figs(indFig), outDirString_In, figureName_png, 'png', resolution);
	end		% End of for i=1 : 1 : noFigures
end		% End of if bSaveResults
% Clear figure handles/graphics objects from workspace to avoid warning, when saving
% workspace
clear h_figs;


%% Save variables of workspace to file
% Obtain current date and time in specific format
% Since use of datestr was no longer recommended, code was changed to use datetime instead
%dt		= datestr(now,'yyyymmdd_HH_MM_SS');
dt		= char(datetime('now', 'Format', 'yyyyMMdd_HH_mm_ss'));

% Save workspace into output directory (optional with user input)
% (Extension".mat" in filename explicitly required, so that Matlab can correctly load 
% workspace file with a "." in its filename)
%strSavedWorkspaceFileName		= ['workspace_', sFunctionName, '_', seqType_MRS_In, '_', dataType_MRS, '_', dt];
%strSavedWorkspaceFileNameFull	= [outDirString_In, strSavedWorkspaceFileName, '.mat'];
strSavedWorkspaceFileName		= ['workspace_', sFunctionName, '_', seqType_MRS, '_', dataType_MRS, '_', dt, '.mat'];
strSavedWorkspaceFileNameFull	= fullfile(outDirString_In, strSavedWorkspaceFileName);
%strSaveWorkspace	= input('Would you like to save all variables of the workspace to file?  ', 's');
strSaveWorkspace	= 'n';
if strcmp(strSaveWorkspace,'y') || strcmp(strSaveWorkspace,'Y')
	save(strSavedWorkspaceFileNameFull);
end

