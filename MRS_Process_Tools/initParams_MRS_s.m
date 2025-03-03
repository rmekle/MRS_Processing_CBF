%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% initParams_MRS_s.m
%
%% Function to initialize settings for processing of magnetic resonance spectroscopy (MRS) data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USAGE
% [paramsMRS_struct] = initParams_MRS_s()
%
% DESCRIPTION:
% Function to initialize parameter values for processing of magnetic resonance 
% spectroscopy (MRS) including preprocessing and metabolite quantification.
%
% INPUTS:
% 	= 
%						
%
% OUTPUTS:
% paramsMRS_struct	= Struct of parameters for processing of MRS data
%
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2025;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [paramsMRS_struct] = initParams_MRS_s()

%% Init input parameters for preprocessing of magnetic resonance spectroscopy (MRS) data
%paramsMRS_struct.dirString_In			= '';
%paramsMRS_struct.dirString_Out			= '';
paramsMRS_struct.fileExt_MRS			= 'dat';		% Currently: 'dat' (raw data) or 'IMA' (DICOM)
paramsMRS_struct.filename_In			= '';
paramsMRS_struct.filename_w_In			= '';
paramsMRS_struct.strStudy_MRS			= '3T_SBAM';		% '3T_Trauma';	'7T_KCL';	'3T_MMs'; '3T_SBAM';
paramsMRS_struct.strVOI_MRS				= 'PCG';			% 'PCG';	% 'HC'; % 'Pons'; % 'CB'; % 'PFC'; % 'PCC';
paramsMRS_struct.seqType_MRS			= 'sLASER';		% 'SPECIAL';	% 'MEGA-PRESS'; % 'sLASER';
paramsMRS_struct.dataType_MRS			= 'mrs_w_ref';		% 'mrs_w_ref';		'mrs_w';	% 'mrs_ref';
paramsMRS_struct.signals_MRS			= 'Spectra';		% 'MMs';	% 'Spectra';
paramsMRS_struct.strOVS_In				= 'wOVS';		% 'wOVS';	% 'woutOVS';
paramsMRS_struct.strOVS_w_In			= 'woutOVS';		% 'wOVS';	% 'woutOVS';
paramsMRS_struct.leftshift_In			= 3;		% 3;	% 2;	% 0;	% 1;
paramsMRS_struct.avgBlockSize_In		= 0;		% 0;	2;		4;		8;		16;

% Parameters for removal of bad averages
paramsMRS_struct.rmbadav_In				= 'y';		% 'y';		'n';
paramsMRS_struct.noSD_In				= 3.2;		% 3.2;	2.6;	5.0;	4.0;	3.0;	2.0;	1.8;
%paramsMRS_struct.digits_noSD_In		= [fix(noSD_In) round(abs(noSD_In-fix(noSD_In))*10)];

% Parameters for spectral registration (aligning of averages/frequency and phase drift
% correction) performed in either frequency or time domain
paramsMRS_struct.strSpecReg_In			= 'SR1';	% To distinguish settings for spectral registration
paramsMRS_struct.driftCorr_In			= 'y';		% 'y';		'n';
paramsMRS_struct.iterin_In				= 20;
paramsMRS_struct.aaDomain_In			= 'f';		% 'f';		't';
paramsMRS_struct.tmaxin_In				= 0.2;		% 0.2;		0.1;
paramsMRS_struct.bTmaxset_In			= 1;
paramsMRS_struct.ppmOption				= 1;
paramsMRS_struct.medin_In				= 'y';		% 'y';	'n';	'a';	'ref';
paramsMRS_struct.alignSS_In				= 2;		% For aligning subspectra (e.g. in SPECIAL)
% Set parameters for drift correction depending on type of data, i.e. whether MRS
% data is spectrum or water signal
% NOTE: Check whether aligning of averages in frequency domain works, if the MR
% spectrum is water signal itself; if not, simply align averages in time domain
switch paramsMRS_struct.dataType_MRS
	case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
		% MR spectrum is provided together without or with unsuppressed water
		% signal and/or with reference scans
		%ppmmin_fix_In		= 1.6;		% 1.6;		1.8;
		%ppmmaxarray_fix_In	= [3.5; 4.0; 5.5];
		%ppmmaxarray_fix_In	= [2.4,2.85,3.35,4.2,4.4,5.2];
		switch paramsMRS_struct.ppmOption
			case 1
				% For MR spectra
				paramsMRS_struct.ppmmin_fix_In			= 1.6;		% 1.6;		1.8;
				paramsMRS_struct.ppmmaxarray_fix_In		= [2.4,2.85,3.35,4.2,4.4,5.2];
			case 2
				% For MR spectra
				paramsMRS_struct.ppmmin_fix_In			= 1.6;
				paramsMRS_struct.ppmmaxarray_fix_In		= [3.5; 4.0; 5.5];
			case 3
				% For MR spectra using settings for water signals
				paramsMRS_struct.ppmmin_fix_In			= 4.2;
				paramsMRS_struct.ppmmaxarray_fix_In		= [5.5 5.5 5.2];
			case 4
				% Wide range to always include water resonance
				paramsMRS_struct.ppmmin_fix_In			= 1.6;
				paramsMRS_struct.ppmmaxarray_fix_In		= [5.5 5.5 5.2];
			case 5
				% For MMs signals
				paramsMRS_struct.ppmmin_fix_In			= 0.2;
				paramsMRS_struct.ppmmaxarray_fix_In		= [3.35,4.2,4.4];
			case 6
				% For MMs signals
				paramsMRS_struct.ppmmin_fix_In			= 0.2;
				pparamsMRS_struct.pmmaxarray_fix_In		= [3.35,4.0,4.1];

			otherwise
				error('%s: Unknown ppmOption = %d!', sFunctionName, paramsMRS_struct.ppmOption);
		end			% End of switch paramsMRS_struct.ppmOption
	case {'water', 'water_ref'}
		% MR spectrum is water signal itself without or with reference scans
		paramsMRS_struct.ppmmin_fix_In		= 4.2;
		paramsMRS_struct.ppmmaxarray_fix_In	= [5.5 5.5 5.2];

	otherwise
		error('%s: Unknown MRS dataType_MRS = %s!', sFunctionName, paramsMRS_struct.dataType_MRS);
end		% End of switch paramsMRS_struct.dataType_MRS

% Additional parameter settings
paramsMRS_struct.bECC_In					= 1;
paramsMRS_struct.bPhaseCorrFreqShift_In		= 0;
paramsMRS_struct.strMinUserIn_In			= 'y';
paramsMRS_struct.plotSwitch_In				= 0;
paramsMRS_struct.reportSwitch_In			= 1;
paramsMRS_struct.strProcessTool_In			= 'FID-A';
paramsMRS_struct.bPrep_MetabQuant			= 1;

end		% End of function