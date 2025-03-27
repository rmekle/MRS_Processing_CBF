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
% [paramsMRS_struct] = initParams_MRS_s(config)
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
% paramsMRS_struct	= Struct of parameter settings for processing of MRS data
%
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2025;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [paramsMRS_struct] = initParams_MRS_s(config)

%% Set string for name of routine and display blank lines for enhanced output visibility
sFunctionName		= 'initParams_MRS_s';
fprintf('\n\n');


%% Init input parameters for (pre)processing of magnetic resonance spectroscopy (MRS) data
% Parameter settings that are the same for all acquisitions
paramsMRS_struct.filename				= '';
paramsMRS_struct.filename_w				= '';

% Parameter settings depending on selected configuration
switch config
	case 'config_Study_sLASER_VOI_dat_MRS_lsN_SDx_y_SR1_ECC'
		% MR spectra in raw data (.dat) format processed using spectral registration (SR1)
		%paramsMRS_struct.dirString_In			= '';
		%paramsMRS_struct.dirString_Out			= '';
		%paramsMRS_struct.filename				= '';
		%paramsMRS_struct.filename_w				= '';
		paramsMRS_struct.strStudy_MRS			= '3T_TGA';		% '3T_Trauma';	'7T_KCL';	'3T_MMs'; '3T_SBAM';	'3T_TGA';
		paramsMRS_struct.seqType_MRS			= 'sLASER';		% 'SPECIAL';	% 'MEGA-PRESS'; % 'sLASER';
		paramsMRS_struct.strVOI_MRS				= 'HC';		% 'PCG';	% 'HC'; % 'Pons'; % 'CB'; % 'PFC'; % 'PCC';
		paramsMRS_struct.fileExt_MRS			= 'dat';		% Currently: 'dat' (raw data) or 'IMA' (DICOM)
		paramsMRS_struct.dataType_MRS			= 'mrs_ref';	% 'mrs_w_ref';		'mrs_w';	% 'mrs_ref';
		paramsMRS_struct.signals_MRS			= 'Spectra';	% 'MMs';	% 'Spectra';
		paramsMRS_struct.strOVS					= 'wOVS';		% 'wOVS';	% 'woutOVS';
		paramsMRS_struct.strOVS_w				= 'wOVS';	% 'wOVS';	% 'woutOVS';
		paramsMRS_struct.leftshift				= 3;			% 3;	% 2;	% 0;	% 1;
		paramsMRS_struct.avgBlockSize			= 0;			% 0;	2;		4;		8;		16;
		
		% Info about processing tool(s) mainly used
		paramsMRS_struct.strProcessTool			= 'FID-A';

		% Parameters for removal of bad averages
		paramsMRS_struct.rmbadav				= 'y';			% 'y';		'n';
		paramsMRS_struct.noSD					= 3.0;			% 3.2;	2.6;	5.0;	4.0;	3.0;	2.0;	1.8;
		%paramsMRS_struct.digits_noSD_In		= [fix(noSD_In) round(abs(noSD_In-fix(noSD_In))*10)];

		% Parameters for spectral registration (aligning of averages/frequency and phase 
		% drift correction) performed in either frequency or time domain
		paramsMRS_struct.strSpecReg				= 'SR1';	% To distinguish settings for spectral registration
		paramsMRS_struct.driftCorr				= 'y';		% 'y';		'n';
		paramsMRS_struct.iterin					= 20;
		paramsMRS_struct.aaDomain				= 'f';		% 'f';		't';
		paramsMRS_struct.tmaxin					= 0.2;		% 0.2;		0.1;
		paramsMRS_struct.bTmaxset				= 1;
		paramsMRS_struct.ppmOption				= 1;
		paramsMRS_struct.medin					= 'y';		% 'y';	'n';	'a';	'ref';
		paramsMRS_struct.alignSS				= 2;		% For aligning subspectra (e.g. in SPECIAL)
		% Set parameters for drift correction depending on type of data, i.e. whether MRS
		% data is spectrum or water signal
		% NOTE: Check whether aligning of averages in frequency domain works, if the MR
		% spectrum is water signal itself; if not, simply align averages in time domain
		switch paramsMRS_struct.dataType_MRS
			case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
				% MR spectrum is provided together without or with unsuppressed water
				% signal and/or with reference scans
				%ppmmin_fix		= 1.6;		% 1.6;		1.8;
				%ppmmaxarray_fix	= [3.5; 4.0; 5.5];
				%ppmmaxarray_fix	= [2.4,2.85,3.35,4.2,4.4,5.2];
				switch paramsMRS_struct.ppmOption
					case 1
						% For MR spectra
						paramsMRS_struct.ppmmin_fix			= 1.6;		% 1.6;		1.8;
						paramsMRS_struct.ppmmaxarray_fix	= [2.4,2.85,3.35,4.2,4.4,5.2];
					case 2
						% For MR spectra
						paramsMRS_struct.ppmmin_fix			= 1.6;
						paramsMRS_struct.ppmmaxarray_fix	= [3.5; 4.0; 5.5];
					case 3
						% For MR spectra using settings for water signals
						paramsMRS_struct.ppmmin_fix			= 4.2;
						paramsMRS_struct.ppmmaxarray_fix	= [5.5 5.5 5.2];
					case 4
						% Wide range to always include water resonance
						paramsMRS_struct.ppmmin_fix			= 1.6;
						paramsMRS_struct.ppmmaxarray_fix	= [5.5 5.5 5.2];
					case 5
						% For MMs signals
						paramsMRS_struct.ppmmin_fix			= 0.2;
						paramsMRS_struct.ppmmaxarray_fix	= [3.35,4.2,4.4];
					case 6
						% For MMs signals
						paramsMRS_struct.ppmmin_fix			= 0.2;
						paramsMRS_struct.pmmaxarray_fix_In	= [3.35,4.0,4.1];

					otherwise
						error('%s: Unknown ppmOption = %d!', sFunctionName, paramsMRS_struct.ppmOption);
				end			% End of switch paramsMRS_struct.ppmOption
			case {'water', 'water_ref'}
				% MR spectrum is water signal itself without or with reference scans
				paramsMRS_struct.ppmmin_fix			= 4.2;
				paramsMRS_struct.ppmmaxarray_fix	= [5.5 5.5 5.2];

			otherwise
				error('%s: Unknown MRS dataType_MRS = %s!', sFunctionName, paramsMRS_struct.dataType_MRS);
		end		% End of switch paramsMRS_struct.dataType_MRS

		% Additional parameter settings
		paramsMRS_struct.bECC					= 1;
		paramsMRS_struct.bPhaseCorrFreqShift	= 0;
		paramsMRS_struct.strMinUserIn			= 'y';
		paramsMRS_struct.plotSwitch				= 0;
		paramsMRS_struct.reportSwitch			= 1;
		paramsMRS_struct.bPrep_MetabQuant		= 1;

	case 'config_Study_sLASER_VOI_dat_water_lsN_SDx_y_SR1_ECC'
		% Unsuppressed MRS water signals in raw data (.dat) format processed using 
		% spectral registration (SR)
		paramsMRS_struct.strStudy_MRS			= '3T_TGA';
		paramsMRS_struct.seqType_MRS			= 'sLASER';
		paramsMRS_struct.strVOI_MRS				= 'HC';
		paramsMRS_struct.fileExt_MRS			= 'dat';
		paramsMRS_struct.dataType_MRS			= 'water';
		paramsMRS_struct.signals_MRS			= 'Spectra';
		paramsMRS_struct.strOVS					= 'woutOVS';
		paramsMRS_struct.strOVS_w				= 'woutOVS';
		paramsMRS_struct.leftshift				= 3;
		paramsMRS_struct.avgBlockSize			= 0;

		% Info about processing tool(s) mainly used
		paramsMRS_struct.strProcessTool			= 'FID-A';

		% Parameters for removal of bad averages
		paramsMRS_struct.rmbadav				= 'y';
		paramsMRS_struct.noSD					= 3.0;

		% Parameters for spectral registration (aligning of averages/frequency and phase 
		% drift correction) performed in either frequency or time domain
		paramsMRS_struct.strSpecReg				= 'SR1';
		paramsMRS_struct.driftCorr				= 'y';
		paramsMRS_struct.iterin					= 20;
		paramsMRS_struct.aaDomain				= 'f';
		paramsMRS_struct.tmaxin					= 0.2;
		paramsMRS_struct.bTmaxset				= 1;
		paramsMRS_struct.ppmOption				= 1;
		paramsMRS_struct.medin					= 'y';
		paramsMRS_struct.alignSS				= 2;
		% Set parameters for drift correction depending on type of data, i.e. whether MRS
		% data is spectrum or water signal
		% NOTE: Check whether aligning of averages in frequency domain works, if the MR
		% spectrum is water signal itself; if not, simply align averages in time domain
		switch paramsMRS_struct.dataType_MRS
			case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
				% MR spectrum is provided together without or with unsuppressed water
				% signal and/or with reference scans
				%ppmmin_fix		= 1.6;		% 1.6;		1.8;
				%ppmmaxarray_fix	= [3.5; 4.0; 5.5];
				%ppmmaxarray_fix	= [2.4,2.85,3.35,4.2,4.4,5.2];
				switch paramsMRS_struct.ppmOption
					case 1
						% For MR spectra
						paramsMRS_struct.ppmmin_fix			= 1.6;		% 1.6;		1.8;
						paramsMRS_struct.ppmmaxarray_fix	= [2.4,2.85,3.35,4.2,4.4,5.2];
					case 2
						% For MR spectra
						paramsMRS_struct.ppmmin_fix			= 1.6;
						paramsMRS_struct.ppmmaxarray_fix	= [3.5; 4.0; 5.5];
					case 3
						% For MR spectra using settings for water signals
						paramsMRS_struct.ppmmin_fix			= 4.2;
						paramsMRS_struct.ppmmaxarray_fix	= [5.5 5.5 5.2];
					case 4
						% Wide range to always include water resonance
						paramsMRS_struct.ppmmin_fix			= 1.6;
						paramsMRS_struct.ppmmaxarray_fix	= [5.5 5.5 5.2];
					case 5
						% For MMs signals
						paramsMRS_struct.ppmmin_fix			= 0.2;
						paramsMRS_struct.ppmmaxarray_fix	= [3.35,4.2,4.4];
					case 6
						% For MMs signals
						paramsMRS_struct.ppmmin_fix			= 0.2;
						paramsMRS_struct.pmmaxarray_fix_In	= [3.35,4.0,4.1];

					otherwise
						error('%s: Unknown ppmOption = %d!', sFunctionName, paramsMRS_struct.ppmOption);
				end			% End of switch paramsMRS_struct.ppmOption
			case {'water', 'water_ref'}
				% MR spectrum is water signal itself without or with reference scans
				paramsMRS_struct.ppmmin_fix			= 4.2;
				paramsMRS_struct.ppmmaxarray_fix	= [5.5 5.5 5.2];

			otherwise
				error('%s: Unknown MRS dataType_MRS = %s!', sFunctionName, paramsMRS_struct.dataType_MRS);
		end		% End of switch paramsMRS_struct.dataType_MRS

		% Additional parameter settings
		paramsMRS_struct.bECC						= 1;
		paramsMRS_struct.bPhaseCorrFreqShift		= 0;
		paramsMRS_struct.strMinUserIn				= 'y';
		paramsMRS_struct.plotSwitch					= 0;
		paramsMRS_struct.reportSwitch				= 0;
		paramsMRS_struct.bPrep_MetabQuant			= 1;

		case 'config_Study_sLASER_VOI_IMA_MRS_lsN_SDx_y_SR1_ECC'
		% MR spectra in DICOM (.IMA) format processed using spectral registration (SR1)
		paramsMRS_struct.strStudy_MRS			= '3T_TGA';
		paramsMRS_struct.seqType_MRS			= 'sLASER';
		paramsMRS_struct.strVOI_MRS				= 'HC';
		paramsMRS_struct.fileExt_MRS			= 'IMA';
		paramsMRS_struct.dataType_MRS			= 'mrs_ref';
		paramsMRS_struct.signals_MRS			= 'Spectra';
		paramsMRS_struct.strOVS					= 'wOVS';
		paramsMRS_struct.strOVS_w				= 'wOVS';
		paramsMRS_struct.leftshift				= 1;
		paramsMRS_struct.avgBlockSize			= 0;
		
		% Info about processing tool(s) mainly used
		paramsMRS_struct.strProcessTool			= 'FID-A';

		% Parameters for removal of bad averages
		paramsMRS_struct.rmbadav				= 'y';
		paramsMRS_struct.noSD					= 3.0;

		% Parameters for spectral registration (aligning of averages/frequency and phase 
		% drift correction) performed in either frequency or time domain
		paramsMRS_struct.strSpecReg				= 'SR1';
		paramsMRS_struct.driftCorr				= 'y';
		paramsMRS_struct.iterin					= 20;
		paramsMRS_struct.aaDomain				= 'f';
		paramsMRS_struct.tmaxin					= 0.2;
		paramsMRS_struct.bTmaxset				= 1;
		paramsMRS_struct.ppmOption				= 1;
		paramsMRS_struct.medin					= 'y';
		paramsMRS_struct.alignSS				= 2;
		% Set parameters for drift correction depending on type of data, i.e. whether MRS
		% data is spectrum or water signal
		% NOTE: Check whether aligning of averages in frequency domain works, if the MR
		% spectrum is water signal itself; if not, simply align averages in time domain
		switch paramsMRS_struct.dataType_MRS
			case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
				% MR spectrum is provided together without or with unsuppressed water
				% signal and/or with reference scans
				%ppmmin_fix		= 1.6;		% 1.6;		1.8;
				%ppmmaxarray_fix	= [3.5; 4.0; 5.5];
				%ppmmaxarray_fix	= [2.4,2.85,3.35,4.2,4.4,5.2];
				switch paramsMRS_struct.ppmOption
					case 1
						% For MR spectra
						paramsMRS_struct.ppmmin_fix			= 1.6;		% 1.6;		1.8;
						paramsMRS_struct.ppmmaxarray_fix	= [2.4,2.85,3.35,4.2,4.4,5.2];
					case 2
						% For MR spectra
						paramsMRS_struct.ppmmin_fix			= 1.6;
						paramsMRS_struct.ppmmaxarray_fix	= [3.5; 4.0; 5.5];
					case 3
						% For MR spectra using settings for water signals
						paramsMRS_struct.ppmmin_fix			= 4.2;
						paramsMRS_struct.ppmmaxarray_fix	= [5.5 5.5 5.2];
					case 4
						% Wide range to always include water resonance
						paramsMRS_struct.ppmmin_fix			= 1.6;
						paramsMRS_struct.ppmmaxarray_fix	= [5.5 5.5 5.2];
					case 5
						% For MMs signals
						paramsMRS_struct.ppmmin_fix			= 0.2;
						paramsMRS_struct.ppmmaxarray_fix	= [3.35,4.2,4.4];
					case 6
						% For MMs signals
						paramsMRS_struct.ppmmin_fix			= 0.2;
						paramsMRS_struct.pmmaxarray_fix_In	= [3.35,4.0,4.1];

					otherwise
						error('%s: Unknown ppmOption = %d!', sFunctionName, paramsMRS_struct.ppmOption);
				end			% End of switch paramsMRS_struct.ppmOption
			case {'water', 'water_ref'}
				% MR spectrum is water signal itself without or with reference scans
				paramsMRS_struct.ppmmin_fix			= 4.2;
				paramsMRS_struct.ppmmaxarray_fix	= [5.5 5.5 5.2];

			otherwise
				error('%s: Unknown MRS dataType_MRS = %s!', sFunctionName, paramsMRS_struct.dataType_MRS);
		end		% End of switch paramsMRS_struct.dataType_MRS

		% Additional parameter settings
		paramsMRS_struct.bECC					= 1;
		paramsMRS_struct.bPhaseCorrFreqShift	= 0;
		paramsMRS_struct.strMinUserIn			= 'y';
		paramsMRS_struct.plotSwitch				= 0;
		paramsMRS_struct.reportSwitch			= 1;
		paramsMRS_struct.bPrep_MetabQuant		= 1;

	otherwise
		error('%s: ERROR: Unknown configuration %s!', sFunctionName, config);
end		% End of switch config

end		% End of function