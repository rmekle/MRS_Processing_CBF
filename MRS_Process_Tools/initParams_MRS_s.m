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

% Set string for name of routine and display blank lines for enhanced output visibility
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
		paramsMRS_struct.strStudy_MRS			= '3T_Trauma';		% '3T_Trauma';	'7T_KCL';	'3T_MMs'; '3T_SBAM';	'3T_TGA';
		paramsMRS_struct.seqType_MRS			= 'sLASER';		% 'SPECIAL';	% 'MEGA-PRESS'; % 'sLASER';
		paramsMRS_struct.strVOI_MRS				= 'HC';		% 'PCG';	% 'HC'; % 'Pons'; % 'CB'; % 'PFC'; % 'PCC';
		paramsMRS_struct.fileExt_MRS			= 'dat';		% Currently: 'dat' (raw data) or 'IMA' (DICOM)
		paramsMRS_struct.dataType_MRS			= 'mrs_w_ref';	% 'mrs_w_ref';		'mrs_w';	% 'mrs_ref';
		paramsMRS_struct.signals_MRS			= 'Spectra';	% 'MMs';	% 'Spectra';
		paramsMRS_struct.strOVS					= 'wOVS';		% 'wOVS';	% 'woutOVS';
		paramsMRS_struct.strOVS_w				= 'woutOVS';	% 'wOVS';	% 'woutOVS';
		paramsMRS_struct.leftshift				= 3;			% 3;	% 2;	% 0;	% 1;
		paramsMRS_struct.avgBlockSize			= 0;			% 0;	2;		4;		8;		16;
		
		% Info about processing tool(s) mainly used
		paramsMRS_struct.strProcessTool			= 'FID-A';

		% Parameters for removal of bad averages
		paramsMRS_struct.rmbadav				= 'y';			% 'y';		'n';
		paramsMRS_struct.noSD					= 3.2;			% 3.2;	2.6;	5.0;	4.0;	3.0;	2.0;	1.8;
		%paramsMRS_struct.digits_noSD_In		= [fix(noSD_In) round(abs(noSD_In-fix(noSD_In))*10)];

		% Parameters for aligning of averages/frequency and phase drift correction using
		% one of the following techniques:
		%	Spectral registration	performed in either frequency or time domain or
		%	Cross-Correlation		performed in frequency domain
		%paramsMRS_struct.strSpecReg				= 'SR1';	% To distinguish settings for spectral registration
		paramsMRS_struct.strFreqPhaseCorr		= 'SR1';	
		paramsMRS_struct.driftCorr				= 'y';		% 'y';		'n';
		paramsMRS_struct.iterin					= 20;
		paramsMRS_struct.aaDomain				= 'f';		% 'f';		't';
		paramsMRS_struct.tmaxin					= 0.2;		% 0.2;		0.1;
		paramsMRS_struct.bTmaxset				= 1;
		paramsMRS_struct.ppmOption				= 1;
		paramsMRS_struct.medin					= 'y';		% 'y';	'n';	'a';	'ref';
		paramsMRS_struct.alignSS				= 2;		% For aligning subspectra (e.g. in SPECIAL)
		% Set parameters for drift correction using spectral registration (SR) 
		% depending on type of data, i.e. whether MRS data are spectra or water signals
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

		% Parameters for aligning of averages/frequency and phase drift correction using
		% one of the following techniques:
		%	Spectral registration	performed in either frequency or time domain or
		%	Cross-Correlation		performed in frequency domain
		paramsMRS_struct.strFreqPhaseCorr		= 'SR1';	
		paramsMRS_struct.driftCorr				= 'y';
		paramsMRS_struct.iterin					= 20;
		paramsMRS_struct.aaDomain				= 'f';
		paramsMRS_struct.tmaxin					= 0.2;
		paramsMRS_struct.bTmaxset				= 1;
		paramsMRS_struct.ppmOption				= 1;
		paramsMRS_struct.medin					= 'y';
		paramsMRS_struct.alignSS				= 2;
		% Set parameters for drift correction using spectral registration (SR) 
		% depending on type of data, i.e. whether MRS data are spectra or water signals
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

		% Parameters for aligning of averages/frequency and phase drift correction using
		% one of the following techniques:
		%	Spectral registration	performed in either frequency or time domain or
		%	Cross-Correlation		performed in frequency domain
		paramsMRS_struct.strFreqPhaseCorr		= 'SR1';	
		paramsMRS_struct.driftCorr				= 'y';
		paramsMRS_struct.iterin					= 20;
		paramsMRS_struct.aaDomain				= 'f';
		paramsMRS_struct.tmaxin					= 0.2;
		paramsMRS_struct.bTmaxset				= 1;
		paramsMRS_struct.ppmOption				= 1;
		paramsMRS_struct.medin					= 'y';
		paramsMRS_struct.alignSS				= 2;
		% Set parameters for drift correction using spectral registration (SR) 
		% depending on type of data, i.e. whether MRS data are spectra or water signals
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

	case 'config_Study_sLASER_VOI_IMA_MRS_lsN_SDx_y_SC1_ECC'
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

		% Parameters for aligning of averages/frequency and phase drift correction using
		% one of the following techniques:
		%	Spectral registration (SR)		performed in either frequency or time domain 
		%	or
		%	Spectral cross-correlation (SC)	performed in frequency domain
		paramsMRS_struct.strFreqPhaseCorr		= 'SC1';

		% Parameters for spectral registration (SR)
		paramsMRS_struct.driftCorr				= 'y';
		paramsMRS_struct.iterin					= 20;
		paramsMRS_struct.aaDomain				= 'f';
		paramsMRS_struct.tmaxin					= 0.2;
		paramsMRS_struct.bTmaxset				= 1;
		paramsMRS_struct.ppmOption				= 1;
		paramsMRS_struct.medin					= 'y';
		paramsMRS_struct.alignSS				= 2;
		% Set parameters for drift correction using spectral registration (SR) 
		% depending on type of data, i.e. whether MRS data are spectra or water signals
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
		
		% Parameters for spectral cross-correlation (SC)
		paramsMRS_struct.structSC.dataFlag		= 'conj';
		paramsMRS_struct.structSC.refSC			= 'f';
		paramsMRS_struct.structSC.filterFlagSC	= 0;
		paramsMRS_struct.structSC.plotFlagSC	= 0;
		paramsMRS_struct.structSC.XnuclOffsetSC	= 0;
		% Set parameters for drift correction using spectral cross-correlation (SC) 
		% depending on type of data, i.e. whether MRS data are spectra or water signals
		switch paramsMRS_struct.dataType_MRS
			case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
				% MR spectrum is provided together without or with unsuppressed water
				% signal and/or with reference scans
				paramsMRS_struct.structSC.minppmSC		= 1.8;
				paramsMRS_struct.structSC.maxppmSC		= 3.6;
			case {'water', 'water_ref'}
				% MR spectrum is water signal itself without or with reference scans
				pparamsMRS_struct.structSC.minppmSC		= 3.75;
				paramsMRS_struct.structSC.maxppmSC		= 5.55;

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

	case 'config_Test_seqType_VOI_ext_MRS_lsN_SDx_y_FreqPhaseCorr_bECC'
		% Test configuration
		% MR spectra in format fileExt processed using frequency and phase correction
		%paramsMRS_struct.dirString_In			= '';
		%paramsMRS_struct.dirString_Out			= '';
		%paramsMRS_struct.filename				= '';
		%paramsMRS_struct.filename_w				= '';
		paramsMRS_struct.strStudy_MRS			= '3T_Test';		% '3T_Trauma';	'7T_KCL';	'3T_MMs'; '3T_SBAM';	'3T_TGA';
		paramsMRS_struct.seqType_MRS			= 'sLASER';		% 'SPECIAL';	% 'MEGA-PRESS'; % 'sLASER';
		paramsMRS_struct.strVOI_MRS				= 'HC';		% 'PCG';	% 'HC'; % 'Pons'; % 'CB'; % 'PFC'; % 'PCC';
		paramsMRS_struct.fileExt_MRS			= 'IMA';		% Currently: 'dat' (raw data) or 'IMA' (DICOM)
		paramsMRS_struct.dataType_MRS			= 'mrs_ref';	% 'mrs_w_ref';		'mrs_w';	% 'mrs_ref';
		paramsMRS_struct.signals_MRS			= 'Spectra';	% 'MMs';	% 'Spectra';
		paramsMRS_struct.strOVS					= 'wOVS';		% 'wOVS';	% 'woutOVS';
		paramsMRS_struct.strOVS_w				= 'wOVS';	% 'wOVS';	% 'woutOVS';
		paramsMRS_struct.leftshift				= 1;			% 3;	% 2;	% 0;	% 1;
		paramsMRS_struct.avgBlockSize			= 0;			% 0;	2;		4;		8;		16;
		
		% Info about processing tool(s) mainly used
		paramsMRS_struct.strProcessTool			= 'FID-A';

		% Parameters for removal of bad averages
		paramsMRS_struct.rmbadav				= 'y';			% 'y';		'n';
		paramsMRS_struct.noSD					= 3.0;			% 3.2;	2.6;	5.0;	4.0;	3.0;	2.0;	1.8;
		%paramsMRS_struct.digits_noSD_In		= [fix(noSD_In) round(abs(noSD_In-fix(noSD_In))*10)];

		% Parameters for aligning of averages/frequency and phase drift correction using
		% one of the following techniques:
		%	Spectral registration (SR)		performed in either frequency or time domain 
		%	or
		%	Spectral cross-correlation (SC)	performed in frequency domain
		paramsMRS_struct.strFreqPhaseCorr		= 'SR1';

		% Parameters for spectral registration (SR)
		paramsMRS_struct.driftCorr				= 'y';		% 'y';		'n';
		paramsMRS_struct.iterin					= 20;
		paramsMRS_struct.aaDomain				= 'f';		% 'f';		't';
		paramsMRS_struct.tmaxin					= 0.2;		% 0.2;		0.1;
		paramsMRS_struct.bTmaxset				= 1;
		paramsMRS_struct.ppmOption				= 1;
		paramsMRS_struct.medin					= 'y';		% 'y';	'n';	'a';	'ref';
		paramsMRS_struct.alignSS				= 2;		% For aligning subspectra (e.g. in SPECIAL)
		% Set parameters for drift correction using spectral registration (SR) 
		% depending on type of data, i.e. whether MRS data are spectra or water signals
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


%% Check(s) on parameter settings
% For studies '3T_Trauma' and ''3T_SBAM', water signals for the HC and the PCG have been
% acquired with different OVS settings (by accident), so adjust for this automatically
%fprintf('\n\n');
if strcmp(paramsMRS_struct.seqType_MRS, 'sLASER')
	switch paramsMRS_struct.strStudy_MRS
		case {'3T_Trauma', '3T_SBAM'}
			switch paramsMRS_struct.strVOI_MRS
				case 'HC'
					paramsMRS_struct.strOVS_w		= 'wOVS';
				case 'PCG'
					paramsMRS_struct.strOVS_w		= 'woutOVS';

				otherwise
					error('%s: ERROR: Unknown VOI %s for study %s!\n', sFunctionName, paramsMRS_struct.strVOI_MRS, paramsMRS_struct.strStudy_MRS);
			end			% End of switch strVOI_MRS

		otherwise
			fprintf('%s: Settings (in (paramsMRS_struct) for OVS not adjusted for strStudy_MRS = %s and strVOI_MRS = %s!\n\n', sFunctionName, paramsMRS_struct.strStudy_MRS, paramsMRS_struct.strVOI_MRS);
	end				% End of switch strVOI_MRS
end		% End of if strcmp(paramsMRS_struct.seqType_MRS, 'sLASER')
fprintf('%s: strStudy_MRS = %s\t strVOI_MRS = %s\n\n\tSettings (in paramsMRS_struct) for OVS are strOVS = %s\t and\t strOVS_w = %s\n\n', sFunctionName, paramsMRS_struct.strStudy_MRS, paramsMRS_struct.strVOI_MRS, paramsMRS_struct.strOVS, paramsMRS_struct.strOVS_w);


end		% End of function [paramsMRS_struct] = initParams_MRS_s(config)


% Introduce two new functions to init struct for spectral registration and spectral
% cross-correlation
% Input parameters are paramsMRS_struct.strFreqPhaseCorr and paramsMRS_struct.dataType_MRS
% Init one struct for SR1, SR2, SR3, and SR4, and use settings for SR1 for all SC options
% Init one struct for SC1, SC2, SC3, and SC4, and use settings for SC1 for all SR options
% Return these two sturcts with larger parameter struct
% Init parameters from structs in preprocess_ListOfFiles_s(...)
% Pass on all parameters into preProcess_MRS_s(...)
% Call options for SC in preProcess_MRS_s(...)
% Check on categories in initParams_MRS_s(...)

function [paramsSpecReg_struct] = initParams_SpecReg_s(strFreqPhaseCorr, dataType_MRS)

% Set string for name of routine and display blank lines for enhanced output visibility
sFunctionName		= 'initParams_specReg_s';
%fprintf('\n\n');


%% Init parameters for frequency and phase correction using spectral registration (SR)
% Initialize (independent) settings as for 'SR1' for all options of frequency and phase
% correction
paramsSpecReg_struct.driftCorr			= 'y';		% 'y';		'n';
paramsSpecReg_struct.iterin				= 20;
paramsSpecReg_struct.aaDomain			= 'f';		% 'f';		't';
paramsSpecReg_struct.tmaxin				= 0.2;		% 0.2;		0.1;
paramsSpecReg_struct.bTmaxset			= 1;
paramsSpecReg_struct.ppmOption			= 1;
paramsSpecReg_struct.medin				= 'y';		% 'y';	'n';	'a';	'ref';
paramsSpecReg_struct.alignSS			= 2;		% For aligning subspectra (e.g. in SPECIAL)

% Modify settings for specifc options of spectral registration (SR), i.e. differences with
% respect to 'SR1'
switch strFreqPhaseCorr
	case {'SR1', 'SC1', 'SC2', 'SC3', 'SC4'}
		% SR in frequency domain
		% No modifications required
		% Use same settings as for 'SR1' for all options of spectral cross-correlation
		% (i.e. 'SC1', 'SC2', etc.) to have these parameters still available in calling
		% routines
	case 'SR2'
		% SR in frequency domain
		% Different range of ppm values
		paramsSpecReg_struct.ppmOption			= 2;
	case 'SR3'
		% SR in time domain
		% tmaxin is set to fixed value
		paramsSpecReg_struct.aaDomain			= 't';
		paramsSpecReg_struct.tmaxin				= 0.2;
		paramsSpecReg_struct.bTmaxset			= 1;
		paramsSpecReg_struct.ppmOption			= 2;	% Not relevant here, since SR in time domain
	case 'SR4'
		% SR in time domain
		% tmaxin is determined from the data in corresponding alignment routine
		paramsSpecReg_struct.aaDomain			= 't';
		paramsSpecReg_struct.tmaxin				= 0.2;
		paramsSpecReg_struct.bTmaxset			= 0;
		paramsSpecReg_struct.ppmOption			= 2;	% Not relevant here, since SR in time domain

	otherwise
		error('%s: ERROR: Unknown strFreqPhaseCorr = %s!', sFunctionName, strFreqPhaseCorr);
end		% End of switch strFreqPhaseCorr

% Set parameters for drift correction using spectral registration (SR)
% depending on type of data, i.e. whether MRS data are spectra or water signals and
% other settings for selected option for spectral registration (i.e. 'SR1', 'SR2', etc.)
% NOTE: Check whether aligning of averages in frequency domain works, if the MR
% spectrum is water signal itself; if not, simply align averages in time domain
switch dataType_MRS
	case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
		% MR spectrum is provided together without or with unsuppressed water
		% signal and/or with reference scans
		%ppmmin_fix		= 1.6;		% 1.6;		1.8;
		%ppmmaxarray_fix	= [3.5; 4.0; 5.5];
		%ppmmaxarray_fix	= [2.4,2.85,3.35,4.2,4.4,5.2];
		switch paramsSpecReg_struct.ppmOption
			case 1
				% For MR spectra
				paramsSpecReg_struct.ppmmin_fix			= 1.6;		% 1.6;		1.8;
				paramsSpecReg_struct.ppmmaxarray_fix	= [2.4,2.85,3.35,4.2,4.4,5.2];
			case 2
				% For MR spectra
				paramsSpecReg_struct.ppmmin_fix			= 1.6;
				paramsSpecReg_struct.ppmmaxarray_fix	= [3.5; 4.0; 5.5];
			case 3
				% For MR spectra using settings for water signals
				paramsSpecReg_struct.ppmmin_fix			= 4.2;
				paramsSpecReg_struct.ppmmaxarray_fix	= [5.5 5.5 5.2];
			case 4
				% Wide range to always include water resonance
				paramsSpecReg_struct.ppmmin_fix			= 1.6;
				paramsSpecReg_struct.ppmmaxarray_fix	= [5.5 5.5 5.2];
			case 5
				% For MMs signals
				paramsSpecReg_struct.ppmmin_fix			= 0.2;
				paramsSpecReg_struct.ppmmaxarray_fix	= [3.35,4.2,4.4];
			case 6
				% For MMs signals
				paramsSpecReg_struct.ppmmin_fix			= 0.2;
				paramsSpecReg_struct.pmmaxarray_fix		= [3.35,4.0,4.1];

			otherwise
				error('%s: Unknown ppmOption = %d!', sFunctionName, paramsSpecReg_struct.ppmOption);
		end			% End of switch paramsSpecReg_struct.ppmOption
	case {'water', 'water_ref'}
		% MR spectrum is water signal itself without or with reference scans
		paramsSpecReg_struct.ppmmin_fix			= 4.2;
		paramsSpecReg_struct.ppmmaxarray_fix	= [5.5 5.5 5.2];

	otherwise
		error('%s: Unknown MRS dataType_MRS = %s!', sFunctionName, paramsSpecReg_struct.dataType_MRS);
end		% End of switch dataType_MRS


end		%End of function [paramsSpecReg_struct] = initParams_SpecReg_s(strFreqPhaseCorr, dataType_MRS)





