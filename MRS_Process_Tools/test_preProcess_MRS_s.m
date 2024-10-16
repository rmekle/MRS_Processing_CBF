clear all;
close all;

% MRS data in DICOM format (.IMA)
dicomDirString_HC		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma_Test/MRS_Trauma_SBA_C_0142_20211104_DICOM/SVS_SLASER_DKD_HC_TE23_WS256_0025/';
dicomDirString_w_HC		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma_Test/MRS_Trauma_SBA_C_0142_20211104_DICOM/SVS_SLASER_DKD_HC_TE23_W8_0027/';
dicomOutDirString_HC	= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma_Test/Test_Processed_DICOM_HC/';

dicomDirString_PCG		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma_Test/MRS_Trauma_SBA_C_0142_20211104_DICOM/SVS_SLASER_DKD_PCG_TE23_WS128_0032/';
dicomDirString_w_PCG	= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma_Test/MRS_Trauma_SBA_C_0142_20211104_DICOM/SVS_SLASER_DKD_PCG_TE23_W8_0034/';
dicomOutDirString_PCG	= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma_Test/Test_Processed_DICOM_PCG/';

% MRS data in raw data format (.dat)
rawDirString_HC			= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma_Test/MRS_Trauma_SBA_C_0142_20211104/';
rawOutDirString_HC		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma_Test/Test_Processed_Raw_HC/';
rawFilename_HC			= '3T_SBA_C_0142_20211104_meas_MID00117_FID135450_svs_slaser_dkd_HC_TE23_WS256.dat';
rawFilename_w_HC		= '3T_SBA_C_0142_20211104_meas_MID00118_FID135451_svs_slaser_dkd_HC_TE23_w8.dat';

rawDirString_PCG		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma_Test/MRS_Trauma_SBA_C_0142_20211104/';
rawOutDirString_PCG		= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma_Test/Test_Processed_Raw_PCG/';
rawFilename_PCG			= '3T_SBA_C_0142_20211104_meas_MID00128_FID135461_svs_slaser_dkd_PCG_TE23_WS128.dat';
rawFilename_w_PCG		= '3T_SBA_C_0142_20211104_meas_MID00129_FID135462_svs_slaser_dkd_PCG_TE23_w8.dat';

dicomRawOutDirString_HC	= '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma_Test/Test_Processed_DICOM_Raw_HC/';
rawDicomOutDirString_HC = '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_BCAN_MRS_Trauma_Test/Test_Processed_Raw_DICOM_HC/';


%% Test for DICOM HC
%{
[out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_s(...
    dicomDirString_HC,...
    dicomOutDirString_HC,...
	'sLASER',...
    'mrs_w_ref',...
    'Filename','',...
	'WaterDirectory', dicomDirString_w_HC,...
    'WaterFilename','',...
	'ReportFilename', 'FirstTests',...
	'OVS','wOVS',...
	'WaterOVS','wOVS',...
    'Leftshift',1,...
    'noStandardDeviation',3.2,...
	'aaDomain','f',...
	'MaxTimeAlignment',0.2,...
	'Iterations',20,...
    'ECC',1,...
    'PhaseFrequencyCorrection',0,...
    'ShowPlots', 0,...
	'MinimizeUserInput','y',...
	'GenerateReport',1);
%}

%% Test for DICOM PCG
%%{
[out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_s(...
    dicomDirString_PCG,...
	dicomOutDirString_PCG,...
	'sLASER',...
    'mrs_w_ref',...
    'Filename','',...
	'WaterDirectory', dicomDirString_w_PCG,...
    'WaterFilename','',...
	'ReportFilename', 'FirstTests',...
	'OVS','wOVS',...
	'WaterOVS','woutOVS',...
    'Leftshift',1,...
    'noStandardDeviation',3.2,...
	'aaDomain','f',...
	'MaxTimeAlignment',0.2,...
	'Iterations',20,...
    'ECC',1,...
    'PhaseFrequencyCorrection',0,...
    'ShowPlots', 0,...
	'MinimizeUserInput','y',...
	'GenerateReport',1);
%}

%% Test for RAW HC
%{
[out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_s(...
    rawDirString_HC,...
    rawOutDirString_HC,...
    'sLASER',...
    'mrs_w_ref',...
    'Filename',rawFilename_HC,...
	'WaterDirectory','',...
    'WaterFilename',rawFilename_w_HC,...
	'ReportFilename', 'FirstTests',...
	'OVS','wOVS',...
	'WaterOVS','wOVS',...
    'Leftshift',3,...
    'noStandardDeviation',3.2,...
	'aaDomain','f',...
	'MaxTimeAlignment',0.2,...
	'Iterations',20,...
    'ECC',1,...
    'PhaseFrequencyCorrection',0,...
    'ShowPlots', 0,...
	'MinimizeUserInput','y',...
	'GenerateReport',1);
%}

%% Test for RAW PCG
%{
[out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_s(...
    rawDirString_PCG,...
    rawOutDirString_PCG,...
    'sLASER',...
    'mrs_w_ref',...
    'Filename',rawFilename_PCG,...
	'WaterDirectory','',...
    'WaterFilename',rawFilename_w_PCG,...
	'ReportFilename', 'FirstTests',...
	'OVS','wOVS',...
	'WaterOVS','woutOVS',...
    'Leftshift',3,...
    'noStandardDeviation',3.2,...
	'aaDomain','f',...
	'MaxTimeAlignment',0.2,...
	'Iterations',20,...
    'ECC',1,...
    'PhaseFrequencyCorrection',0,...
    'ShowPlots', 0,...
	'MinimizeUserInput','y',...
	'GenerateReport',1);
%}

%% Test for Combination (WS: RAW, W: DICOM) HC
%{
[out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_s(...
    rawDirString_HC,...
    rawDicomOutDirString_HC,...
    'sLASER',...
    'mrs_w_ref',...
    'WaterDirectory', dicomDirString_w_HC,...
    'WaterLeftshift', 1,...
    'Filename',rawFilename_HC,...
    'Leftshift', 3,...
    'ReportFilename', 'firstTests',...
    'ECC',1,...
    'PhaseFrequencyCorrection',1,...
    'ShowPlots', 1);
%}

%% Test for Combination (WS: DICOM, W: RAW) HC
%{
[out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_s(...
    dicomDirString_HC,...
    dicomRawOutDirString_HC,...
    'sLASER',...
    'mrs_w_ref',...
    'Leftshift', 1,...
    'WaterDirectory', rawDirString_HC,...
    'WaterFilename',rawFilename_w_HC,...
    'WaterLeftshift', 3,...
    'ReportFilename', 'firstTests',...
    'ECC',1,...
    'PhaseFrequencyCorrection',1,...
    'ShowPlots', 1);
%}