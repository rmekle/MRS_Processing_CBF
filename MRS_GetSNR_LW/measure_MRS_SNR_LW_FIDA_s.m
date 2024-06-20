%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% measure_MRS_SNR_LW_FID-A_s.m
%
%% Function to measure the signal-to-noise ratio (SNR) and the linewidth (LW)/FWHM of a 
%% selected resonance/peak in single volume magnetic resonance spectroscopy (MRS) data
%  using functions from the MRS processing toolkit FID-A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  USAGE
%	[data_MRS, SNR, FWHM, ph0, info]	= measure_MRS_SNR_LW_FIDA_s(dirString, filename_MRS, filename_w, dataFormat_MRS, signal_ppmRange, noise_ppmRange, LWpeak_ppmRange, zp_factor, outDirString, dataType_MRS, bAutoPhase, bOutFile, plotswitch, seqType_MRS, procParams, Bo_field, spectralWidth, TE, TR)
%
%  INPUTS
%	dirString		(required)	String variable for the name of the directory containing
%								the MRS data file(s)
%	filename_MRS	(required)	Input filename of MRS data							 
%	filename_w		(required)	Input filename of additional unsuppressed water signal
%								 that has to be in the same directory as the MRS data
%								 (can be empty, if not needed, e.g. if MRS data is
%								 already unsuppressed water signal to be analyzed)
%	dataFormat_MRS	(required)	String specifying data format of MRS data:
%								 'DICOM' = DICOM, 'IMA' = Siemens DICOM, 
%								 'rda' = Siemens rda, 'dat' = Siemens raw data, 
%								 'lcmRAW" = .RAW file for LCModel
%	signal_ppmRange (required)	ppm range (min, max) for signal measurement in MRS data
%	noise_ppmRange	(required)	ppm range (min, max) for noise measurement in MRS data
%	LWpeak_ppmRange	(required)	ppm range (min, max) for FWHM measurement in MRS data
%	zp_factor		(required)	zero-padding factor (used for method 1 of the linewidth
%								measurement); default value there is = 8, but here
%								zp_factor has to be specified, since linewidth measurement
%								routine also uses plotswitch as input argument
%	outDirString	(required)	Directory path for saving output & workspace
%   dataType_MRS    (required)  String describing the type of MRS data:
%					'mrs'		= MR spectrum without water signal, 
%					'mrs_w'		= MR spectrum with unsuppressed water signal
%					'mrs_w_ref'	= MR spectrum with unsuppressed water signal and reference
%									(water) scans
%					'mrs_ref'	= MR spectrum with reference (water) scans
%					'water'		= MR spectrum is unsuppressed water signal itself
%					'water_ref' = MR spectrum is unsuppressed water signal itself with
%									reference (water) scans (should be very rare!)
%   bAutoPhase		(optional)	Boolean (0 or 1) to select whether automatic zero order
%								 phasing with respect to LWpeak_ppmRange is applied prior
%								 to SNR and linewidth measurements or not; default = 0
%	bOutToFile		(optional)	Boolean (0 or 1) to select whether output values are
%									written to file or not; default = 1
%	plotswitch		(optional)	Switch for displaying plots: 1 = ON, 0 = OFF; default = 1
%   seqType_MRS     (optional)  String to specify sequence type; essential for processing
%                                raw data/dat files
%   procParams     (optional)  Struct with parameters for processing raw data/dat files:
%                               nSD         =  Number of standard deviations as threshold     
%									for removal of bad averages. Default is 3.2.
%                               aaDomain    = Perform the spectral registration (drift 
%                                   correction) using the full spectrum ('t'), or only a 
%                                   limited frequency range ('f').  Default is 'f'.
%                               tmaxin      = Duration (in sec.) of the time domain 
%                                   signal used in the spectral registration (drift 
%                                   correction). Default is 0.2 sec.
%                               iterin      = Maximum number of allowed iterations for 
%                                   the spectral registration to converge. Default is 20.
%	Bo_field		(optional)	Field strength of acquisition for MRS data in tesla (T);
%								if provided, it overwrites the calculated value
%	spectralWidth	(optional)	Spectral width (bandwidth) of MRS data acquisition Hz; 
%								if provided, it overwrites the calculated value
%	TE				(optional)	Echo time TE in ms
%	TR				(optional)	Repetition time TR in ms
%		
%  OUTPUTS
%	data_MRS		MRS dataset in FID-A structure format
%	SNR				Signal-to-noise ratio (SNR) of MRS data using a selected peak height
%	FWHM			Linewidth (LW)/FWHM of selected peak
%   ph0				Zero order phase reference (linewidth) peak of spectrum was phased to;
%					 equal to zero, if no phasing was applied
%	info			If dataFormat_MRS = 'DICOM' or 'IMA':
%					 struct that contains the metadata from the compliant DICOM or DICOS 
%					 file specified as MRS input data file obtained using the Matlab
%					 routine dicominfo plus a .csa field that is a struct that contains 
%					 values from the 'SIEMENS CSA HEADER' private group in a DICOM file 
%					 produced by a Siemens MR scanner;
%					If dataFormat_MRS = 'rda' or 'dat' or 'lcmRAW': Empty
%
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2020, 2021, 2024; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data_MRS, SNR, FWHM, ph0, info]	= measure_MRS_SNR_LW_FIDA_s(dirString, filename_MRS, filename_w, dataFormat_MRS, signal_ppmRange, noise_ppmRange, LWpeak_ppmRange, zp_factor, outDirString, dataType_MRS, bAutoPhase, bOutToFile, plotswitch, seqType_MRS, procParams, Bo_field, spectralWidth, TE, TR)

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'measure_MRS_SNR_LW_FIDA_s.m';
% sMsg_newLines		= sprintf('\n\n');
% disp(sMsg_newLines);
fprintf('\n\n');


%% Check on input arguments
% NOTE: From Matlab R2019b on, declare function argument validation 'arguments' can be
% used to validate input arguments to a function

% Set total # and # of required input arguments
noInputs_max	= 19;
noInputs_req	= 10;

% First input arguments have to be provided
if nargin < noInputs_req
    error('%s: Error: %d required input arguments have to be provided!\n\n', sFunctionName, noInputs_req);
elseif nargin == noInputs_req
	disp([sFunctionName, ': Parameter plotswitch will be set to 1!']);
	plotswitch = 1;
% elseif nargin > noInputs_req 
%     if isnumeric(plotswitch) == 0 | plotswitch < 0 | plotswitch > 1
%         disp([sFunctionName, ': Input parameter plotswitch not correctly identified, will be set to 1!']);
%         plotswitch = 1;
% 	end
end

% Check on required input arguments
if( isempty(filename_MRS) )
    error('%s: Error: Empty input filername_MRS = "%s!"\n\n', sFunctionName, filename_MRS);
end

% Check on optional input arguments
if nargin < noInputs_max
	TR	= [];
	if nargin < (noInputs_max-1)
		TE	= [];
		if nargin < (noInputs_max-2)
			spectralWidth	= [];
			if nargin < (noInputs_max-3)
				Bo_field	= [];
				if nargin < (noInputs_max-4)
					procParams.nSD		= 3.2;
					procParams.aaDomain	= 'f';
					procParams.tmaxin	= 0.2;
					procParams.iterin	= 20;
					if nargin < (noInputs_max-5)
						seqType_MRS	= '';
						if nargin < (noInputs_max-6)
							plotswitch = 1;
							if nargin < (noInputs_max-7)
								bOutToFile = 1;
								if nargin < (noInputs_max-8)
									bAutoPhase = 0;
								end
							end
						end
					end
				end
			end
		end
	end
end
% % Assign output directory to current directory, if not specified
% if( nargin < (noInputs_req+2))
% 	outDirString	= pwd;
% end

% Ensure that input and output directory string end with file separator, i.e. '/' or '\'
% (Windows handles the Unix '/' just fine)
if( ~strcmp( dirString(end), filesep ) )
	dirString	= [dirString, filesep];
end

if( ~strcmp( outDirString(end), filesep ) )
	outDirString	= [outDirString, filesep];
end


%% Set parameters for measurements of SNR and LW of MRS data
%dirString_In			= '';
%outDirString_In		= '';


% Set (additional) parameters
% Gyromagnetic ratio/larmor constant for 1H in MHz/Tesla
gammabar_FIDA		= 42.577;		% FID-A toolkit
gammabar_IDEA		= 42.575575;	% Siemens IDEA software


%% Load MRS data according to its data format into data structure used in FID-A toolkit
% Create full filename for MRS data to be analyzed and create empty info struct
fullFilename_MRS		= [dirString, filename_MRS];
info					= [];
switch dataFormat_MRS
	case {'DICOM', 'IMA'}
		if strcmp(dataFormat_MRS, 'DICOM')
			warning('%s: Warning: MRS data format %s currently loaded as Siemens DICOM = IMA!', sFunctionName, dataFormat_MRS);
		end
		% Read Siemens DICOM header to access information about scan parameters and
		% determine scan parameters from header that are required for data loading
		info				= SiemensCsaParse(fullFilename_MRS);
		Bo_field_csa		= info.csa.ImagingFrequency/gammabar_FIDA;
		spectralWidth_csa	= 1/(info.csa.RealDwellTime*1E-9);	% Dwell time is in ns
		TE_csa				= info.csa.EchoTime;
		TR_csa				= info.csa.RepetitionTime;
		% Use scan parameters from header for data loading, if not provided as input
		% arguments
		if( ~exist('Bo_field', 'var') || isempty(Bo_field) )
			Bo_field		= Bo_field_csa;
		end
		if( ~exist('spectralWidth', 'var') || isempty(spectralWidth) )
			spectralWidth	= spectralWidth_csa;
		end
		if( ~exist('TE', 'var') || isempty(TE) )
			TE				= TE_csa;
		end
		if( ~exist('TR', 'var') || isempty(TR) )
			TR				= TR_csa;
		end
		data_MRS	= io_loadspec_IMA(fullFilename_MRS, Bo_field, spectralWidth, TE, TR);
        % DICOM/IMA data is loaded using routine SiemensCsaParse.m by Chris
        % Rodgers, where data is converted into single precision;
        % To make calculations more coherent for all MRS data formats, data
        % in the struct data_MRS is converted back to double precision here
        data_MRS.fids       = double(data_MRS.fids);
        data_MRS.specs      = double(data_MRS.specs);
	case 'rda'
		data_MRS	= io_loadspec_rda(fullFilename_MRS);
	case 'dat'
		%data_MRS	= io_loadspec_twix(fullFilename_MRS);
		
        % Since MRS data stored as .dat files has not undergone coil
        % combination nor averaging these steps need to be performed before
        % subsequent SNR and linewidth measurements
		
		% Check whether sequence type is already provided; if not, obtain sequence type
		% from user input
		if( ~exist('seqType_MRS', 'var') || isempty(seqType_MRS) )
			seqType_MRS			= input('Please enter MRS sequence type: ','s');
		end
        
		% Load and preprocess MRS data according to sequence and data type
		% Set remaining preprocessing options
		strOVS			= 'woutOVS';
		strMinUserIn	= 'y';
		report_switch	= 1;		% 1;		% 0;
		[data_MRS,data_MRSw,data_MRS_noproc,data_MRSw_noproc] = preProcess_MRS_RawData_s(dirString,outDirString, ...
			filename_MRS,filename_w,seqType_MRS,dataType_MRS,strOVS,procParams.nSD,procParams.aaDomain, ...
			procParams.tmaxin,procParams.iterin,plotswitch,strMinUserIn,report_switch);
	case 'lcmRAW'
		% Load MRS data that was potentially preprocessed and save in .RAW format for
		% LCModel (using the FID-A routine io_writelcm(...))
		% (choosing type = 'dat' - .raw file generated by FID-A using io_writelcm)
		data_MRS		= io_readlcmraw(fullFilename_MRS, 'dat'); 
		
	otherwise
		error('%s: ERROR: Unknown MRS data format %s!', sFunctionName, dataFormat_MRS);
end


% %% If plotswitch == 0, turn off figure plotting
% if( plotswitch == 0)
% 	set(gcf,'Visible','off');
% 	set(0,'DefaultFigureVisible','off');
% end


%% Apply automatic zero order phasing with respect to LWpeak_ppmRange prior to SNR and 
% linewidth measurements, if selected
% Init desired zero order phase in degrees to phase reference peak to
ph				= 0;
%ph0				= 0:
%data_MRS_ph		= data_MRS;
if bAutoPhase
	[data_MRS_ph, ph0]		= op_autophase(data_MRS, LWpeak_ppmRange(1), LWpeak_ppmRange(2), ph);
else
	% No phasing applied
	ph0				= 0:
	data_MRS_ph		= data_MRS;
end		% End of if bAutoPhase


%% Compute the SNR of the MRS data for a selected peak and noise range
%[SNR,MRS_signal,MRS_noisesd]	= op_getSNR(data_MRS, signal_ppmRange(1), signal_ppmRange(2), noise_ppmRange(1), noise_ppmRange(2));
[SNR,MRS_signal,MRS_noiseSD]	= op_getSNR_s(data_MRS_ph, signal_ppmRange(1), signal_ppmRange(2), noise_ppmRange(1), noise_ppmRange(2), plotswitch);
%SNR	= 100;
if MRS_signal < 1 || MRS_noiseSD < 1
	fprintf('\nMRS_signal = %.2e\tMRS_noiseSD = %.2e\t=>\tSNR = %.2f \t = %.1f \n', MRS_signal, MRS_noiseSD, SNR, SNR);
else
	fprintf('\nMRS_signal = %.2f\tMRS_noiseSD = %.2f\t=>\tSNR = %.2f \t = %.1f \n', MRS_signal, MRS_noiseSD, SNR, SNR);
end


%% Determine the linewidth (LW) of the MRS data
%figure;		% Since first plot of op_getLW(---) uses plot wout figure command
%[FWHM]	= op_getLW(data_MRS, LWpeak_ppmRange(1), LWpeak_ppmRange(2), zp_factor);
[FWHM]	= op_getLW_s(data_MRS_ph, LWpeak_ppmRange(1), LWpeak_ppmRange(2), zp_factor, plotswitch);
%FWHM	= 6;


% %% Turn figure plotting back on for subsequent figures
% if( plotswitch == 0)
% 	%set(gcf,'Visible','on');
% 	set(0,'DefaultFigureVisible','on');
% end


%% Write results of measurements to a textfile, if selected
% Use same filename as for input filename for MRS data, but replace '.' with '_' to avoid
% two '.'s in new filename
[filepath_MRS, name_MRS, ext_MRS]    = fileparts(filename_MRS);
new_ext_MRS             = ['_', ext_MRS(2:end)];
filename_MRS_Results    = fullfile(outDirString, [name_MRS, new_ext_MRS, '_Results_SNR_LW.txt']);

% Open textfile and write results to textfile, if selected
if bOutToFile
	fprintf('\n%s: \tbOutToFile = %d => Writing values for SNR and linewidth of spectra to file!\n\n', sFunctionName, bOutToFile);
	fid_txt             = fopen(filename_MRS_Results, 'wt');
	if( fid_txt == -1 )
		error('%s: Could not open textfile %s!\n', sFunctionName, filename_MRS_Results);
	else
		nbytes = fprintf(fid_txt, 'Results from SNR and Linewidth Measurements for MRS Data\n');
		nbytes = fprintf(fid_txt, '=========================================================================================\n');
		nbytes = fprintf(fid_txt, '=========================================================================================\n\n');
		nbytes = fprintf(fid_txt, 'MRS data path \t= %s\n', dirString);
		nbytes = fprintf(fid_txt, 'MRS data file \t= %s\n', [name_MRS, ext_MRS]);
		nbytes = fprintf(fid_txt, 'Data format \t= %s\n', dataFormat_MRS);
		nbytes = fprintf(fid_txt, '\n\n');
		nbytes = fprintf(fid_txt, 'Signal-to-Noise Measurement\n');
		nbytes = fprintf(fid_txt, '-----------------------------------------------------------------------------------------\n\n');
		nbytes = fprintf(fid_txt, 'Signal as maximum peak height of magnitude spectrum in \t\tppm range = [%.2f, %.2f]\n', signal_ppmRange(1), signal_ppmRange(2));
		nbytes = fprintf(fid_txt, 'Noise as standard deviation of real part of signal in \t\tppm range = [%.2f, %.2f]\n\n', noise_ppmRange(1), noise_ppmRange(2));

		% Select precision for printing floating point numbers based on magnitude of values
		% (If values < 1, a precision string of %.2f will print 0.00 into textfile)
		if MRS_signal < 1 || MRS_noiseSD < 1
			nbytes = fprintf(fid_txt, 'MRS Signal \t=\t%.2e\n', MRS_signal);
			nbytes = fprintf(fid_txt, 'MRS Noise \t=\t%.2e\n', MRS_noiseSD);
		else
			nbytes = fprintf(fid_txt, 'MRS Signal \t=\t%.2f\n', MRS_signal);
			nbytes = fprintf(fid_txt, 'MRS Noise \t=\t%.2f\n', MRS_noiseSD);
		end
		nbytes = fprintf(fid_txt, 'MRS SNR \t=\t%.2f \t~=\t%.1f\n', SNR, SNR);
		nbytes = fprintf(fid_txt, '\n\n');
		nbytes = fprintf(fid_txt, 'Linewidth / Full Width Half Maximum (FWHM) Measurement\n');
		nbytes = fprintf(fid_txt, '-----------------------------------------------------------------------------------------\n\n');
		nbytes = fprintf(fid_txt, 'Linewidth of a reference peak in the spectrum in \t\tppm range = [%.2f, %.2f]\n\n', LWpeak_ppmRange(1), LWpeak_ppmRange(2));
		nbytes = fprintf(fid_txt, 'FWHM \t\t=\t%.2f Hz \t~=\t%.1f Hz\n', FWHM, FWHM);
		nbytes = fprintf(fid_txt, '\n\n');

		% Close textfile
		status_txt	= fclose(fid_txt);
	end			% End of if( fid_txt == -1 )
%else
	%fprintf('\n%s: \tbOutToFile = %d => Values for SNR and linewidth of spectra NOT written to file!\n\n', sFunctionName, bOutToFile);
end		% End of if bOutToFile


%% Save variables of workspace to file
% Save workspace into output directory (optional with user input)
% (Extension".mat" in filename explicitly required, so that Matlab can correctly load 
% workspace file with a "." in its filename)
strSavedWorkspaceFileName		= ['workspace_', sFunctionName];
strSavedWorkspaceFileNameFull	= [outDirString, strSavedWorkspaceFileName];
%strSaveWorkspace	= input('Would you like to save all variables of the workspace to file?  ', 's');
strSaveWorkspace	= 'n';		% strSaveWorkspace	= 'y';
if strcmp(strSaveWorkspace,'y') || strcmp(strSaveWorkspace,'Y')
	save(strSavedWorkspaceFileNameFull);
end

