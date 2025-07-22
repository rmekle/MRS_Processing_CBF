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
		% Parameter to select frequency and phase drift correction method
		% Setting to align subspectra independent of selected correction technique
		paramsMRS_struct.driftCorr				= 'y';		% 'y';		'n';
		paramsMRS_struct.strFreqPhaseCorr		= 'SR1';
		paramsMRS_struct.alignSS				= 2;		% For aligning subspectra (e.g. in SPECIAL)

		% Note that both structs for SR and SC will be initialized to have these settings
		% available in calling routines (possible use of both methods)
		% Parameters for spectral registration (SR)
		[paramsMRS_struct.structSR]				= initParams_SpecReg_s(paramsMRS_struct.strFreqPhaseCorr, ...
																		paramsMRS_struct.dataType_MRS);
		% Parameters for spectral cross-correlation (SC)
		plotFlagSC_init		= 0;
		[paramsMRS_struct.structSC]				= initParams_SpecCrossCorrel_s(paramsMRS_struct.strFreqPhaseCorr, ...
																		paramsMRS_struct.dataType_MRS, plotFlagSC_init);
		
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
		% Parameter to select frequency and phase drift correction method
		% Setting to align subspectra independent of selected correction technique
		paramsMRS_struct.driftCorr				= 'y';
		paramsMRS_struct.strFreqPhaseCorr		= 'SR1';
		paramsMRS_struct.alignSS				= 2;		% For aligning subspectra (e.g. in SPECIAL)

		% Note that both structs for SR and SC will be initialized to have these settings
		% available in calling routines (possible use of both methods)
		% Parameters for spectral registration (SR)
		[paramsMRS_struct.structSR]				= initParams_SpecReg_s(paramsMRS_struct.strFreqPhaseCorr, ...
																		paramsMRS_struct.dataType_MRS);
		% Parameters for spectral cross-correlation (SC)
		plotFlagSC_init		= 0;
		[paramsMRS_struct.structSC]				= initParams_SpecCrossCorrel_s(paramsMRS_struct.strFreqPhaseCorr, ...
																		paramsMRS_struct.dataType_MRS, plotFlagSC_init);

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
		% Parameter to select frequency and phase drift correction method
		% Setting to align subspectra independent of selected correction technique
		paramsMRS_struct.driftCorr				= 'y';
		paramsMRS_struct.strFreqPhaseCorr		= 'SR1';
		paramsMRS_struct.alignSS				= 2;		% For aligning subspectra (e.g. in SPECIAL)

		% Note that both structs for SR and SC will be initialized to have these settings
		% available in calling routines (possible use of both methods)
		% Parameters for spectral registration (SR)
		[paramsMRS_struct.structSR]				= initParams_SpecReg_s(paramsMRS_struct.strFreqPhaseCorr, ...
																		paramsMRS_struct.dataType_MRS);
		% Parameters for spectral cross-correlation (SC)
		plotFlagSC_init		= 0;
		[paramsMRS_struct.structSC]				= initParams_SpecCrossCorrel_s(paramsMRS_struct.strFreqPhaseCorr, ...
																		paramsMRS_struct.dataType_MRS, plotFlagSC_init);

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
		%	Spectral registration	performed in either frequency or time domain or
		%	Cross-Correlation		performed in frequency domain
		% Parameter to select frequency and phase drift correction method
		% Setting to align subspectra independent of selected correction technique
		paramsMRS_struct.driftCorr				= 'y';
		paramsMRS_struct.strFreqPhaseCorr		= 'SC1';
		paramsMRS_struct.alignSS				= 2;		% For aligning subspectra (e.g. in SPECIAL)

		% Note that both structs for SR and SC will be initialized to have these settings
		% available in calling routines (possible use of both methods)
		% Parameters for spectral registration (SR)
		[paramsMRS_struct.structSR]				= initParams_SpecReg_s(paramsMRS_struct.strFreqPhaseCorr, ...
																		paramsMRS_struct.dataType_MRS);
		% Parameters for spectral cross-correlation (SC)
		plotFlagSC_init		= 0;
		[paramsMRS_struct.structSC]				= initParams_SpecCrossCorrel_s(paramsMRS_struct.strFreqPhaseCorr, ...
																		paramsMRS_struct.dataType_MRS, plotFlagSC_init);
		
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
		%	Spectral registration	performed in either frequency or time domain or
		%	Cross-Correlation		performed in frequency domain
		% Parameter to select frequency and phase drift correction method
		% Setting to align subspectra independent of selected correction technique
		paramsMRS_struct.driftCorr				= 'y';
		paramsMRS_struct.strFreqPhaseCorr		= 'SR1';
		paramsMRS_struct.alignSS				= 2;		% For aligning subspectra (e.g. in SPECIAL)
		
		% Note that both structs for SR and SC will be initialized to have these settings
		% available in calling routines (possible use of both methods)
		% Parameters for spectral registration (SR)
		[paramsMRS_struct.structSR]				= initParams_SpecReg_s(paramsMRS_struct.strFreqPhaseCorr, ...
																		paramsMRS_struct.dataType_MRS);
		% Parameters for spectral cross-correlation (SC)
		plotFlagSC_init		= 0;
		[paramsMRS_struct.structSC]				= initParams_SpecCrossCorrel_s(paramsMRS_struct.strFreqPhaseCorr, ...
																		paramsMRS_struct.dataType_MRS, plotFlagSC_init);
		
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





function [paramsSpecReg_struct] = initParams_SpecReg_s(strFreqPhaseCorr, dataType_MRS)

% Set string for name of routine and display blank lines for enhanced output visibility
sFunctionName		= 'initParams_specReg_s';
%fprintf('\n\n');


%% Init parameters for frequency and phase correction using spectral registration (SR)
% Initialize (independent) settings as for 'SR1' for all options of frequency and phase
% correction
%paramsSpecReg_struct.driftCorr			= 'y';		% 'y';		'n';
paramsSpecReg_struct.iterin				= 20;
paramsSpecReg_struct.aaDomain			= 'f';		% 'f';		't';
paramsSpecReg_struct.tmaxin				= 0.2;		% 0.2;		0.1;
paramsSpecReg_struct.bTmaxset			= 1;
paramsSpecReg_struct.ppmOption			= 1;
paramsSpecReg_struct.medin				= 'y';		% 'y';	'n';	'a';	'ref';
%paramsSpecReg_struct.alignSS			= 2;		% For aligning subspectra (e.g. in SPECIAL)

% Modify settings for specifc options of spectral registration (SR), i.e. differences 
% with respect to 'SR1'
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
				error('%s: Unknown paramsSpecReg_struct.ppmOption = %d!', sFunctionName, paramsSpecReg_struct.ppmOption);
		end			% End of switch paramsSpecReg_struct.ppmOption
	case {'water', 'water_ref'}
		% MR spectrum is water signal itself without or with reference scans
		paramsSpecReg_struct.ppmmin_fix			= 4.2;
		paramsSpecReg_struct.ppmmaxarray_fix	= [5.5 5.5 5.2];

	otherwise
		error('%s: Unknown MRS dataType_MRS = %s!', sFunctionName, paramsSpecReg_struct.dataType_MRS);
end		% End of switch dataType_MRS


end		% End of function [paramsSpecReg_struct] = initParams_SpecReg_s(strFreqPhaseCorr, dataType_MRS)






function [paramsSC_struct] = initParams_SpecCrossCorrel_s(strFreqPhaseCorr, dataType_MRS, plotFlagSC)

% Set string for name of routine and display blank lines for enhanced output visibility
sFunctionName		= 'initParams_specReg_s';
%fprintf('\n\n');


%% Init parameters for frequency and phase correction using spectral cross-correlation (SC)
% Initialize (independent) settings as for 'SC1' for all options of frequency and phase
% correction
paramsSC_struct.dataFlag			= 'conj';	% FIDs were read in as conjugate complex
paramsSC_struct.refSC				= 'f';		% Use first FID as reference
paramsSC_struct.filterFlagSC		= 0;		% Default is no filtering (linebroadening)
paramsSC_struct.plotFlagSC			= plotFlagSC;	% Default is no plotting, (=0)
paramsSC_struct.XnuclOffsetSC		= 4.65;		% Default ppm offset for 1H
paramsSC_struct.ppmOption			= 1;		% Default setting for ppm values
paramsSC_struct.LB					= 0;		% Default linebroadening factor, if filterFlag = 0
paramsSC_struct.GF					= 1000;		% Default (Gaussian) apodization factor, if filterFlag = 0

% Modify settings for specifc options of spectral cross-correlation (SC), i.e. differences
% with respect to 'SC1'
% SC in frequency domain
switch strFreqPhaseCorr
	case {'SC1', 'SR1', 'SR2', 'SR3', 'SR4'}
		% For MR spectra wihtout baseline issue(s)
		% No modifications required
		% Use same settings as for 'SC1' for all options of spectral registration
		% (i.e. 'SR1', 'SR2', etc.) to have these parameters still available in calling
		% routines
	case 'SC2'
		% For MR spectra with baseline issue(s)
		% Different range of ppm values
		paramsSC_struct.ppmOption			= 2;
	case 'SC3'
		% Apply filtering of MRS data (linebroadening and Gaussian apodization)
		% Use default values for linebroadening factor and Gaussian apodization factor
		paramsSC_struct.filterFlagSC		= 1;
		paramsSC_struct.LB					= 5;
		paramsSC_struct.GF					= 0.15;		
	case 'SC4'
		% NOT YET!
		error('%s: ERROR: Option for strFreqPhaseCorr = %s not yet implemented!', sFunctionName, strFreqPhaseCorr);
		
	otherwise
		error('%s: ERROR: Unknown strFreqPhaseCorr = %s!', sFunctionName, strFreqPhaseCorr);
end		% End of switch strFreqPhaseCorr

% Set parameters for drift correction using spectral cross-correlation (SC)
% depending on type of data, i.e. whether MRS data are spectra or water signals and
% other settings for selected option for spectral cross-correlation (i.e. 'SC1', 'SC2', etc.)
switch dataType_MRS
	case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
		% MR spectrum is provided together without or with unsuppressed water
		% signal and/or with reference scans
		switch paramsSC_struct.ppmOption
			case 1
				% For MR spectra without baseline issue(s)
				paramsSC_struct.minppmSC		= 1.8;
				paramsSC_struct.maxppmSC		= 3.6;
			case 2
				% For MR spectra with baseline issue(s)
				paramsSC_struct.minppmSC		= 1.8;
				paramsSC_struct.maxppmSC		= 2.2;

			otherwise
				error('%s: Unknown paramsSC_struct.ppmOption = %d!', sFunctionName, paramsSC_struct.ppmOption);
		end			% End of switch paramsSC_struct.ppmOption
	case {'water', 'water_ref'}
		% MR spectrum is water signal itself without or with reference scans
		paramsSC_struct.minppmSC		= 3.75;
		paramsSC_struct.maxppmSC		= 5.55;

	otherwise
		error('%s: Unknown MRS dataType_MRS = %s!', sFunctionName, dataType_MRS);
end		% End of switch dataType_MRS


end		% End of function [paramsSC_struct] = initParams_SpecCrossCorrel_s(strFreqPhaseCorr, dataType_MRS, plotFlagSC)

