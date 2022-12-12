%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% preProcess_MRS_s.m
%
%% Function to preprocess single volume magnetic resonance spectroscopy (MRS) data
%  using functions from the MRS processing toolkit FID-A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USAGE
% [out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_s(dirString,outDirString,seqType,dataType,options)
% 
% DESCRIPTION:
% Function for processing Siemens MRS data in .dat format (twix raw data) or in .IMA
% format (DICOM data)
% using functions from the MRS processing toolkit FID-A
% Includes combination of receiver coil channels (if required), removal of bad averages, 
% frequency drift correction, averaging, eddy current correction (optional), phase and 
% frequency correction.
% 
% INPUTS:
% dirString    = String variable for the name of the directory containing
%                   the water suppressed .dat file or .IMA files
% dirString_w  = (Optional) ['DirectoryWater'] String variable for the name of the directory containing
%                   the water unsuppressed .dat file or .IMA files,
%                   optional because water unsupressed data is optional and
%                   dat file can also be in dirString
% outDirString = String variable for the name of the output directory, i.e. the directory,
%					where all output files are saved to
% filename     = (Optional) ['Filename'] String variable for the name of the water suppressed .dat file,
%					e.g. 3T_20170510_PetraO_meas_MID00091_FID153257_rm_special_RF_ACC.dat
%					or specialDLPFC.dat - Optional because not necessary
%					for .IMA data
% filename_w   = (Optional) ['UnsupressedFilename'] String variable for the name of the water unsuppressed .dat file,
%					e.g. 3T_20170510_PetraO_meas_MID00091_FID153259_rm_special_RF_ACC.dat
%					or specialDLPFC_w.dat - Optional because water
%					unsupressed signal is optional and water unsupressed
%					can be a directory containing .IMA files
% filename_r   = (Optional) ['ReportFilename'] String variable for the name
%                   of the file - Optional because the name will just be
%                   automatically generated otherwise
% seqType	   = String specifying the MRS sequence type used for data acquisition, e.g.
%					'PRESS', 'STEAM', 'SPECIAL', 'sLASER', 'MEGA-PRESS'
% dataType	   = String describing the type of MRS data:
%					'mrs'		= MR spectrum without water signal, 
%					'mrs_w'		= MR spectrum with unsuppressed water signal
%					'mrs_w_ref'	= MR spectrum with unsuppressed water signal and reference
%									(water) scans
%					'mrs_ref'	= MR spectrum with reference (water) scans
%					'water'		= MR spectrum is unsuppressed water signal itself
%					'water_ref' = MR spectrum is unsuppressed water signal itself with
%									reference (water) scans (should be very rare!)
% strOVS	   = (optional) ['OVS'] Character array that specifies whether MR spectrum acquired 
%					 with OVS ('wOVS') or withoutOVS ('woutOVS') is used for processing. 
%					 Default is 'woutOVS', which means that OVS was not used.
% strOVS_w	   = (optional) ['OVSwater'] Character array that specifies whether water acquired with OVS ('wOVS')
%					 or withoutOVS ('woutOVS') is used for processing.
%					 Default is 'woutOVS', which means that OVS was not used.
% leftshift	   = (optional) ['Leftshift'] # of points to leftshift all FIDs, i.e. to remove leading 
%					 datapoints from the FID  to get rid of 1st order phase. Default is 0.
% leftshift_w  = (optional) ['WaterLeftshift'] # of points to leftshift all water unsupressed FIDs, i.e. to remove leading 
%					 datapoints from the FID  to get rid of 1st order phase. Default is the same as leftshift.
% nSD		   = (Optional) ['StandardDeviation'] # of standard deviations for bad average removal. Default
%					value is 3.2.
% aaDomain     = (Optional) ['aaDomain'] Perform the spectral registration (drift correction) using
%                   the full spectrum ('t'), or only a limited frequency
%                   range ('f').  Default is 'f'.
% tmaxin       = (Optional) ['DriftCorrectionDuration'] Duration (in sec.) of the time domain signal
%                   used in the spectral registration (drift correction).
%                   Default is 0.2 sec.
% iterin       = (Optional) ['Iterations']  Maximum number of allowed iterations for the spectral
%                   registration to converge. Default is 20.
% bECC 		   = (optional) ['ECC'] Boolean that specifies whether eddy current correction (ECC) 
%					 should be performed or not. Default is 0.
% bPhaseCorrFreqShift = (Optional) ['PhaseFrequencyCorrection'] Boolean that specifies whether phase correction and
%					frequency shifting should be performed or not. Default is 0.
% plotSwitch   = (Optional)	['ShowPlots'] Switch for displaying plots: 1 = ON, 0 = OFF. Default is 0
% strMinUserIn = (Optional) ['MinimizeUserInput'] String that specifies whether user input/interaction should be
%					 minimized or not; 'y' or 'Y' lead to minimization, 'n' or 'N' do not
%					Default is 'y'.
% reportSwitch = (Optional) ['GenerateReport'] Switch for generating an html report with corresponding  
%					figures and a readme file: 1 = ON, 0 = OFF. Default is 1. 
% 
% OUTPUTS:
% out          = Fully preprocessed, water suppressed output spectrum
% out_w        = Fully preprocessed, water unsuppressed output spectrum
% out_noproc   = Water suppressed output spectrum without preprocessing 
%                   (No bad-averages removal, no frequency drift correction)
% out_w_noproc = Water unsuppressed output spectrum without preprocessing
%
% If reference sigbals are present in input data:
%
% out_ref_ECC			= Fully preprocessed, water reference signal(s) for eddy current
%							compensation (ECC)
% out_ref_Quant			= Fully preprocessed, water reference signal(s) for metabolite
%							quantification (Quant)
% out_ref_ECC_noproc	= Water reference signal(s) for ECC without preprocessing 
% out_ref_Quant_noproc	= Water reference signal(s) for Quant without preprocessing
%
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2021;
% Ivo Opitz, Charite Universitätsmedizin Berlin, Germany, 2022;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out,out_w,out_noproc,out_w_noproc,out_ref_ECC,out_ref_Quant,out_ref_ECC_noproc,out_ref_Quant_noproc] = preProcess_MRS_s(dirString,outDirString,seqType,dataType,options)

%% Parse arguments
% FLAG: Modified
arguments
    dirString       {mustBeText}
    outDirString    {mustBeText}
    seqType         {mustBeText}
    dataType        {mustBeText}
    options.Filename                {mustBeText} = ''
    options.WaterDirectory          {mustBeText} = ''
    options.WaterFilename           {mustBeText} = ''
    options.ReportFilename          {mustBeText} = ''
    options.OVS			            {mustBeMember(options.OVS,{'wOVS', 'woutOVS'})} = 'woutOVS'
    options.WaterOVS				{mustBeMember(options.WaterOVS,{'wOVS', 'woutOVS'})} = 'woutOVS'
    options.Leftshift               (1,1) {mustBeNumeric}   = 0
    options.WaterLeftshift          (1,1) {mustBeNumeric}   = NaN
    options.noStandardDeviation     (1,1) double            = 3.2
    options.aaDomain                {mustBeMember(options.aaDomain,{'t', 'f'})} = 'f'
    options.MaxTimeAlignment		(1,1) double            = 0.2
    options.Iterations				(1,1) {mustBeNumeric}   = 20
    options.ECC						(1,1) {islogical}       = 0
    options.PhaseFrequencyCorrection(1,1) {islogical}		= 0
    options.ShowPlots				(1,1) {islogical}       = 0
    options.MinimizeUserInput		{mustBeMember(options.MinimizeUserInput,{'y', 'Y', 'n', 'N'})} = 'y'
    options.GenerateReport          (1,1) {islogical}       = 1
end

% Convert optional inputs to regular variables
    dirString_w     = options.WaterDirectory;
    filename        = options.Filename;
    filename_w      = options.WaterFilename;
    filename_r      = options.ReportFilename;
    strOVS          = options.OVS;
    strOVS_w        = options.WaterOVS;
    leftshift       = options.Leftshift;
    if isnan(options.WaterLeftshift)
        leftshift_w = leftshift;
    else
        leftshift_w = options.WaterLeftshift;
    end
    nSD             = options.noStandardDeviation;
    aaDomain        = options.aaDomain;
    tmaxin          = options.MaxTimeAlignment;
    iterin          = options.Iterations;
    bECC            = options.ECC;
    bPhaseCorrFreqShift = options.PhaseFrequencyCorrection;
    plotSwitch      = options.ShowPlots;
    strMinUserIn    = options.MinimizeUserInput;
    reportSwitch    = options.GenerateReport;

%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'preProcess_MRS_s.m';
sMsg_newLines		= sprintf('\n\n');
sMsg_newLine		= sprintf('\n');
disp(sMsg_newLines);


%% Process names of directories and input filenames and create directory for html report
% Use filenames for the MRS datasets provided as input arguments
% NOTE: It is assumed that water suppressed and water unsuppressed data are in same
% directory! Addendum: Not anymore! 
% For MRS raw data (.dat) assumption should hold
% For MRS DICOM data (.IMA) each spectrum or water signal is in its own directory

% Ensure that input and output directory names are not empty
% Note that the directory name for the water signal can be empty, e.g. when preprocessing
% is to be applied to MRS raw data (.dat), since then MR spectrum and water signal can be
% in same directory
if isempty(dirString)
	error('%s: Error: Name of directory for MRS signal %s is empty!\n\n', sFunctionName, dirString);
end

if isempty(outDirString)
	error('%s: Error: Name of output directory for processed MRS signals %s is empty!\n\n', sFunctionName, outDirString);
end

% Ensure that input and output directory strings end with file separator, i.e. '/' or '\'
% (Windows handles the Unix '/' just fine)
if( ~strcmp( dirString(end), filesep ) )
	dirString	= [dirString, filesep];
end

% FLAG: Modified
if(~isempty(dirString_w))
    if( ~strcmp( dirString_w(end), filesep ) )
        dirString_w	= [dirString_w, filesep];
    end
end

if( ~strcmp( outDirString(end), filesep ) )
	outDirString	= [outDirString, filesep];
end

% Obtain different parts of input filenames
% FLAG: Modified
% Set flags to indicate whether MRS data is in DICOM format (.IMA) or not
[sPathStrSpec,nameSpec,extSpec] 	= fileparts(filename);
if strcmp(extSpec, '.dat')
    isIMA = 0;
else
    isIMA = 1;
end

[sPathStr_w,name_w, ext_w] 			= fileparts(filename_w);
if strcmp(ext_w, '.dat')
    isIMA_w = 0;
else
    isIMA_w = 1;
end

% FLAG: Modified
% % Derive filenames for spectrum and water signal from either input filename of report
% % file or from current date and time
% % (for MRS DICOM data, since filenames are empty, since directories are provided instead)
% if filename_r ~= ""
% 	nameSpec	= filename_r;
% 	name_w		= [filename_r '_w'];
% else
% 	if nameSpec == ""
% 		strDate		= datestr(datetime(now, 'ConvertFrom', 'datenum'));
% 		strDate		= strrep(strDate, ' ', '_');
% 		strDate		= strrep(strDate, ':', '-');
% 		nameSpec	= strDate;
% 	end
% 	if name_w == ""
% 		strDate		= datestr(datetime(now, 'ConvertFrom', 'datenum'));
% 		strDate		= strrep(strDate, ' ', '_');
% 		strDate		= strrep(strDate, ':', '-');
% 		name_w		= [strDate '_w'];
% 	end
% end

% For MRS DICOM data (.IMA)), derive filenames for spectrum and water signal from
% corresponding sub-directory names
% Obtain names of sub-directories from splitting corresponding paths into cell arrays
% (first and last cell of cell array are empty, if directory is of from "/home/.../test/')
if isIMA
	%dirParts	= regexp(dirString, filesep, 'split');
	dirParts	= strsplit(dirString, filesep);
	nameSpec	= dirParts{end-1};
	%nameSpec	= dirParts{length(dirParts)-1};
end

if isIMA_w
	if isempty(dirString_w)
		error('%s: Error: Name of directory for unsuppressed water signal %s is empty!\n\n', sFunctionName, dirString_w);
	else
		%dirParts_w	= regexp(dirString_w, filesep, 'split');
		dirParts_w	= strsplit(dirString_w, filesep);
		name_w		= dirParts_w{end-1};
		%name_w		= dirParts_w{length(dirParts_w)-1};
	end
end

% Make a new directory for the output report and figures each, if not already existent,
% and if desired
if reportSwitch == 1
	%mkdir([outDirString nameSpec '/report']);
	%mkdir([outDirString nameSpec '/report/figs']);
    reportDirStr		= [nameSpec '_report/'];
    reportFigDirStr		= [nameSpec '_report/figs/'];	
	if ~exist([outDirString reportDirStr], 'dir' )
		mkdir([outDirString reportDirStr]);
	end
	if ~exist( [outDirString reportFigDirStr], 'dir' )
		mkdir([outDirString reportFigDirStr]);
	end
end


% % Create filenames for saving of processed output depending on sequence type
% if( strcmp(seqType, 'MEGA-PRESS') )
% 	% Assuming that MEGA-PRESS editing is usually performed without OVS
% 	% Assuming that no bad average removal is applied for any of the water signals
% 	outFileName				= [nameSpec, sprintf('_%.1f', nSD)];
% 	outFileName_w			= [name_w, '_w', sprintf('_%.1f', nSD)];
% 	outFileName_ref_ECC		= [nameSpec, '_ref_ECC', sprintf('_%.1f', nSD)];
% 	outFileName_ref_Quant	= [nameSpec, '_ref_Quant', sprintf('_%.1f', nSD)];
% else
% 	% Water signals and/or MR spectra were acquired with or without OVS
% 	% Assuming that water reference signals for ECC are always acquired in same way as 
% 	% the MR spectrum,i.e. if OVS is automatically turned on (as it is for dkd_sLaser), 
% 	% then water reference signals for ECC are acquired with OVS ('wOVS'); 
% 	% Assuming that water reference signals for quantification are always acquired without
% 	% OVS ('woutOVS') to avoid MT effects
% 	% Assuming that no bad average removal is applied for any of the water signals
% 	outFileName				= [nameSpec, '_', strOVS, sprintf('_%.1f', nSD)];
% 	outFileName_w			= [name_w, '_w', '_', strOVS];
% 	outFileName_ref_ECC		= [nameSpec, '_ref_ECC', sprintf('%d', 8), '_', strOVS];
% 	outFileName_ref_Quant	= [nameSpec, '_ref_Quant', sprintf('%d', 8), '_', 'woutOVS'];
% 	%outFileName_w			= [name_w, '_w', '_', strOVS, sprintf('_%.1f', nSD)];
% 	%outFileName_ref_ECC		= [nameSpec, '_ref_ECC', '_', strOVS, sprintf('_%.1f', nSD)];
% 	%outFileName_ref_Quant	= [nameSpec, '_ref_Quant', '_', strOVS, sprintf('_%.1f', nSD)];
% end


%% Set parameters for figure display
% h		= figure('position', [left bottom width height]);
fig_left	= 20;
fig_bottom	= 50;
fig_width	= 560;
fig_height	= 420;
fig_dist_b	= 510;
fig_dist_l	= 570;


%% Depending on type of data provided load all corresponding data files
% Init # of shots (averages) for water unsuppressed signals and water reference signals 
noAvg_w			= 0;
noAvg_ref_ECC	= 0;
noAvg_ref_Quant	= 0;

% Display info
fprintf('%s:\n', sFunctionName);
switch dataType 
	case 'mrs'
		% Only MR spectrum
		disp('Data is MR spectrum only ...');
		with_water				= false;
		with_ref				= false;
		out_w					= struct([]);
		out_w_noproc			= struct([]);
		out_ref_ECC				= struct([]);
		out_ref_ECC_noproc		= struct([]);
		out_ref_Quant			= struct([]);
		out_ref_Quant_noproc	= struct([]);
	case 'mrs_w'
		% MR spectrum is provided together with unsuppressed water signal
		% FLAG: Modified
		% Check on filename or directory name depending on MRS data type
        if(isIMA_w)
            if isempty(dirString_w)
                error('%s: Error: Directory name for unsuppresed water signal %s is empty!\n\n', sFunctionName, dirString_w);
            end
        else
            if isempty(filename_w)
                error('%s: Error: Filename for unsuppresed water signal %s is empty!\n\n', sFunctionName, filename_w);
            end
        end
		disp('Data is MR spectrum with additional unsuppressed water signal ...');
		with_water				= true;
		with_ref				= false;
		out_ref_ECC				= struct([]);
		out_ref_ECC_noproc		= struct([]);
		out_ref_Quant			= struct([]);
		out_ref_Quant_noproc	= struct([]);
	case 'mrs_w_ref'
		% MR spectrum is provided together with unsuppressed water signal and reference
		% scans
		% FLAG: Modified
		% Check on filename or directory name depending on MRS data type
        if(isIMA_w)
            if isempty(dirString_w)
                error('%s: Error: Directory name for unsuppresed water signal %s is empty!\n\n', sFunctionName, dirString_w);
            end
        else
            if isempty(filename_w)
                error('%s: Error: Filename for unsuppresed water signal %s is empty!\n\n', sFunctionName, filename_w);
            end
        end
		disp('Data is MR spectrum with additional unsuppressed water signal and with reference scans ...');
		with_water				= true;
		with_ref				= true;
	case 'mrs_ref'
		% MR spectrum is provided together with reference scans	
		disp('Data is MR spectrum with reference scans ...');
		with_water				= false;
		with_ref				= true;
		out_w					= struct([]);
		out_w_noproc			= struct([]);
	case 'water'
		% MR spectrum is unsuppressed water signal ('water') itself
		disp('MR spectrum is unsuppressed water signal itself ...');
		with_water				= false;
		with_ref				= false;
		out_w					= struct([]);
		out_w_noproc			= struct([]);
		out_ref_ECC				= struct([]);
		out_ref_ECC_noproc		= struct([]);
		out_ref_Quant			= struct([]);
		out_ref_Quant_noproc	= struct([]);
	case 'water_ref'
		% MR spectrum is unsuppressed water signal ('water') itself with reference scans
		disp('MR spectrum is unsuppressed water signal itself with reference scans...');
		with_water				= false;
		with_ref				= true;
		out_w					= struct([]);
		out_w_noproc			= struct([]);		
		
	otherwise
		error('%s: Unknown MRS dataType = %s!', sFunctionName, dataType);
end		% End of switch dataType 
				
% if( strcmp(dataType, 'mrs_w') )
% 	% MR spectrum is provided together with unsuppressed water signal
% 	if isempty(filename_w)
% 		error('%s: Error: Filename for unsuppresed water signal %s is empty!\n\n', sFunctionName, filename_w);
% 	end
% 	disp('Data is MR spectrum with additional unsuppressed water signal ...');
% 	with_water		= true;
% else
% 	% Either MR spectrum is unsuppressed water signal ('water') itself or MR spectrum is 
% 	% provided without unsuppressed water signal
% 	if( strcmp(dataType, 'water') )
% 		disp('MR spectrum is unsuppressed water signal itself ...');
% 	else
% 		disp('Data is MR spectrum without additional unsuppressed water signal ...');
% 	end
% 	with_water		= false;
%     out_w			= struct([]);
%     out_w_noproc	= struct([]);
% end
disp(sMsg_newLines);

% Read in the 'main' data together with possibly existing reference scans and, 
% if available, the additional unsuppressed water signal
% FLAG: Modified
% Selects data loading function depending on data type
%out_raw				= io_loadspec_twix([dirString filename]);
if(isIMA)
    [out_raw, out_ref_raw]  = io_loadspec_IMA_s(dirString);
else
    [out_raw, out_ref_raw]	= io_loadspec_twix_s([dirString filename]);
end


% Convert single precision data (default format used my mapVBVD.m for imaging data) 
% into double precision for processing, if not empty
if ~isempty(out_raw)
	out_raw.fids		= double(out_raw.fids);
	out_raw.specs		= double(out_raw.specs);
end
if ~isempty(out_ref_raw)
	out_ref_raw.fids	= double(out_ref_raw.fids);
	out_ref_raw.specs	= double(out_ref_raw.specs);
end
disp(newline);

if with_water
	%disp('***WITH ADDITIONAL WATER UNSUPPRESSED DATA***');
	fprintf('%s: ***WITH ADDITIONAL WATER UNSUPPRESSED DATA***\n', sFunctionName);
    
    % FLAG: Modified
    % Selects data loading function depending on data type
    if isIMA_w
        if isempty(dirString_w)
            error('%s: Error: Name of directory for unsuppressed water signal %s is empty!\n\n', sFunctionName, dirString_w);
		else
			out_w_raw		= io_loadspec_IMA_s(dirString_w);
		end
    else
        if isempty(dirString_w)
            out_w_raw		= io_loadspec_twix_s([dirString filename_w]);
        else
            out_w_raw		= io_loadspec_twix_s([dirString_w filename_w]);
        end
    end
	
	% Convert single precision data (default format used my mapVBVD.m for imaging data) 
	% into double precision for processing, if not empty
	if ~isempty(out_w_raw)
		out_w_raw.fids		= double(out_w_raw.fids);
		out_w_raw.specs		= double(out_w_raw.specs);
	end
	
	% Store # of shots (averages) for water signal
	noAvg_w					= out_w_raw.averages;
else
	%disp('***WITHOUT ADDITIONAL WATER UNSUPPRESSED DATA***');
	fprintf('%s: ***WITHOUT ADDITIONAL WATER UNSUPPRESSED DATA***\n', sFunctionName);
    %out_w			= struct([]);
    %out_w_noproc	= struct([)];
end
disp(sMsg_newLines);


%% Preprocess all data according to sequence type and types of data provided
% Display info
fprintf('%s: Start of preprocessing of MRS data ... \n\n\n', sFunctionName);
switch seqType
	case {'PRESS', 'STEAM', 'sLASER'}
		%% Is this Dinesh K. Deelchand's single voxel (sLaser) sequence from CMRR, U Minnesota
		isSVSdkd_seq = contains(out_raw.seq,'svs_slaser_dkd');
		
		% NOTE:
		% If it is Dinesh's sequence , leftshift all FIDs, i.e. remove specifc # of 
		% points from the beginning of each FID (it is programmed this way to avoid
		% effect of digital filter on top of echo (first point of FID) by staritng data
		% acquisition slightly earlier)
		% Specific # of points in this case should be equal to 3 for oversampled raw 
		% data (.dat on Siemens)
		% Leftshift MR spectrum and, if existent, also the unsuppresed water signal
		if isSVSdkd_seq
			%leftshift	= 3;
			out_raw		= op_leftshift(out_raw, leftshift);
			if with_water
				out_w_raw		= op_leftshift(out_w_raw, leftshift_w);
			end		% End of if with_water			
		end		% End of if isSVSdkd_seq
		
		
		%% Check whether reference scans are included in the MRS data
		if with_ref
			% If svs_slaser_dkd sequence was used, create two data objects for reference
			% scans: one for ECC and one for quantification (Quant)
			if isSVSdkd_seq
				% Create two copies of original data object for reference scans, since
				% several attributes and all flags will remain the same
				out_ref_ECC_raw		= out_ref_raw;
				out_ref_Quant_raw	= out_ref_raw;
				
				% Determine # of all reference scans assuming the information about 
				% averages in data struct for reference scans already includes possibly
				% existing subspectra (see routine io_loadspec_twix_s.m)
				%sz_ref			= size(out_ref_raw.fids);
				noRefScans_all	= out_ref_raw.averages; 
				
				% Determine indices for references scnas for ECC and for Quant
				% Reference scans for ECC for this sequence are stored as scans 
				% 1, 2, ..., (# of all ref scans/4) and as 
				% (# of all ref scans/2+1, ..., (# of all ref scans/2 + # of all ref scans/4) 
				% with the last index = (3/4*# of all ref scans)
				indRefScans_ECC		= [1:(noRefScans_all/4) (noRefScans_all/2+1):(3/4*noRefScans_all)];
				
				% Reference scans for Quant for this sequence are stored as scans 
				% (# of all ref scans/4+1), ...,  (# of all ref scans/2) and as 
				% (3/4*# of all ref scans/2+1, ..., (# of all ref scans) 
				indRefScans_Quant	= [(noRefScans_all/4+1):(noRefScans_all/2) (3/4*noRefScans_all+1):noRefScans_all];

				% Extract data for reference scans for ECC into corresponding data struct
				% Use substruct indexing to extract selected data independent of # of
				% dimensions of data object for all reference scans
				% Define indexing structure for data object (array)
				% S.type is character vector or string scalar containing (), {}, or .,
				% specifying the subscript type; here it is '()' that is used for indexing
				% S.subs is cell array, character vector, or string scalar containing the
				% actual subscripts; here it is cell array of {':'} for each data dimension
				S.type			= '()';
				S.subs			= repmat({':'}, 1, ndims(out_ref_raw.fids));
				
				% Select subscripts (indices) in dimension for averages of data object
				% for all reference scans
				S.subs{out_ref_raw.dims.averages}	= indRefScans_ECC;
				out_ref_ECC_raw.fids				= subsref(out_ref_raw.fids, S);
				out_ref_ECC_raw.specs				= subsref(out_ref_raw.specs, S);
				
				% Adjust information about data object for reference scans for ECC
				% # of averages here should be equal to # of averages for all reference 
				% scans divided by two; same for rawAverages
				% Note: # of subspectra and information about dimensions (dims) should 
				% remain the same
				out_ref_ECC_raw.sz				= size(out_ref_ECC_raw.fids);
				out_ref_ECC_raw.averages		= out_ref_raw.averages/2;
				out_ref_ECC_raw.rawAverages		= out_ref_raw.rawAverages/2;
				noAvg_ref_ECC					= out_ref_ECC_raw.averages;
				
				% Extract data for reference scans for ECC into corresponding data struct
				% Same substruct indexing object can be used, since both types of
				% reference scans, ECC and Quant, are stored in same dimension(s)
				% Select subscripts (indices) in dimension for averages of data object
				% for all reference scans
				S.subs{out_ref_raw.dims.averages}	= indRefScans_Quant;
				out_ref_Quant_raw.fids				= subsref(out_ref_raw.fids, S);
				out_ref_Quant_raw.specs				= subsref(out_ref_raw.specs, S);
				
				% Adjust information about data object for reference scans for Quant
				% # of averages here should be equal to # of averages for all reference 
				% scans divided by two; same for rawAverages
				% Note: # of subspectra and information about dimensions (dims) should 
				% remain the same
				out_ref_Quant_raw.sz			= size(out_ref_Quant_raw.fids);
				out_ref_Quant_raw.averages		= out_ref_raw.averages/2;
				out_ref_Quant_raw.rawAverages	= out_ref_raw.rawAverages/2;
				noAvg_ref_Quant					= out_ref_Quant_raw.averages;
				
				% Leftshift FIDs of all reference scans by specific # of points
				out_ref_ECC_raw		= op_leftshift(out_ref_ECC_raw, leftshift);
				out_ref_Quant_raw	= op_leftshift(out_ref_Quant_raw, leftshift);
				
			else
				warning('%s: Reference scans option for sequence "%s" not yet implemented!\nReference scans will NOT be processed!', sFunctionName, out_raw.seq);
				% Create empty output structs for reference scans and set boolean
				% parameter with_ref to false to avaoid error messages and to enable 
				% processing of the MR spectrum or water signal itself
				out_ref_ECC				= struct([]);
				out_ref_ECC_noproc		= struct([]);
				out_ref_Quant			= struct([]);
				out_ref_Quant_noproc	= struct([]);
				with_ref				= false;
			end		% End of if isSVSdkd_seq 
			
		end		% End of if with_ref
		disp(sMsg_newLines);
		
		
		%% Create filenames for saving of processed output depending on sequence type
		% Information that is used here is only available after loading the data
		% Water signals and/or MR spectra were acquired with or without OVS
		% Assuming that water reference signals for ECC are always acquired in same way as
		% the MR spectrum, i.e. if OVS is automatically turned on (as it is for 
		% svs_dkd_sLaser), then water reference signals for ECC are als acquired with 
		% OVS ('wOVS');
		% Assuming that water reference signals for quantification are always acquired 
		% without OVS ('woutOVS') to avoid MT effects
		% Assuming that no bad average removal is applied for any of the water signals
		outFileName				= [nameSpec, '_', strOVS, sprintf('_%.1f', nSD)];
		outFileName_w			= [name_w, '_w', '_', strOVS_w];
		%outFileName_w			= [name_w, '_w', sprintf('%d', noAvg_w), '_', strOVS_w];
		outFileName_ref_ECC		= [nameSpec, '_ref_ECC', sprintf('%d', noAvg_ref_ECC), '_', strOVS];
		outFileName_ref_Quant	= [nameSpec, '_ref_Quant', sprintf('%d', noAvg_ref_Quant), '_', 'woutOVS'];
		%outFileName_w			= [name_w, '_w', '_', strOVS, sprintf('_%.1f', nSD)];
		%outFileName_ref_ECC		= [nameSpec, '_ref_ECC', '_', strOVS, sprintf('_%.1f', nSD)];
		%outFileName_ref_Quant	= [nameSpec, '_ref_Quant', '_', strOVS, sprintf('_%.1f', nSD)];

		
		%% Combine signals from different coil elements
        % FLAG: Modified
        
        % Zeroth step is to look if MRS data is in .IMA (DICOM) format which is
        % usually already combined on the scanner, and, thus coil combination can be
        % skipped for DICOM data
		% Index for coil dimension in data struct should then also be zero, since coil
		% dimension does not exist; 
		% Thus, perform coil combination, if NRS data is not DICOM and index for coil
		% dimension is non-zero
		% if ~(isIMA && isIMA_w)
        if ~(isIMA && isIMA_w) && (out_raw.dims.coils ~= 0)
		    % First step should be to combine coil channels. For this find the coil phases 
		    % from water unsuppressed data, if available; otherwise from the MR spectra
		    % Arguments referring to the nPos_ccth point of the FID and weighting of channels 
		    % based on maximum signal ('w') are ignored, when coilcombos are provided as input
		    nPos_cc				= 1;
		    nPos_cc_w			= 1;
		    nPos_cc_ref_ECC		= 1;
		    nPos_cc_ref_Quant	= 1;
		    
		    % If separate water scans and reference scans were acquired, get coil phases from
		    % both types of signals, and for coil combination use 
		    % if reference scans exist,
		    %	coil phases from reference scans for reference scans and MR spectra, since
		    %	reference scans were acquired together with MR spectra
		    %	coil phases from water scans for water scans
		    %
		    % if only water scans exist, 
		    %	coil phases from water scans for water scans and for MR spectra
		    %
		    % if neither reference scans nor water scans exist,
		    %	coil phases from MR spectra for MR spectra
		    
		    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		    %
		    % NOTE: If reference scans exist, it is assumed that the dkd sequence from CMRR 
		    % was used, i.e. two types of reference data exist: one for ECC and one for 
		    % quantification (Quant); 
		    % Use reference scans for ECC also for processing of spectra, but combine
		    % reference scans for quantification using their own coil phases
		    %
		    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		    
		    % First, obtain all possible coil phases from existing types of data
		    if with_ref && ~isIMA
			    % Obtain coil phases and amplitudes from reference scans for ECC and for Quant
			    coilcombos_ref_ECC		= op_getcoilcombos(out_ref_ECC_raw,nPos_cc_w,'w');
			    coilcombos_ref_Quant	= op_getcoilcombos(out_ref_Quant_raw,nPos_cc_w,'w');
			    % Combine reference scans using respective coil phases
			    [out_ref_ECC_cc,fid_ref_ECC_pre,spec_ref_ECC_pre,ph_ref_ECC,sig_ref_ECC]			= ...
				    op_addrcvrs(out_ref_ECC_raw,nPos_cc_ref_ECC,'w',coilcombos_ref_ECC);
			    [out_ref_Quant_cc,fid_ref_Quant_pre,spec_ref_Quant_pre,ph_ref_Quant,sig_ref_Quant]	= ...
				    op_addrcvrs(out_ref_Quant_raw,nPos_cc_ref_Quant,'w',coilcombos_ref_Quant);
		    end		% End of if with_ref		
		    if with_water && ~isIMA_w
			    % Obtain coil phases and amplitudes from unsuppressed water signal
			    %coilcombos		= op_getcoilcombos(out_w_raw,1);
			    coilcombos_w	= op_getcoilcombos(out_w_raw,nPos_cc_w,'w');
			    % Combine water scans using respective coil phases
			    [out_w_cc,fid_w_pre,spec_w_pre,ph_w,sig_w]	= op_addrcvrs(out_w_raw,nPos_cc_w,'w',coilcombos_w);
            end			% % End of if with_water
		    
		    % Obtain coil phases and amplitudes from (averaged) MR spectra
		    %coilcombos_mrs	= op_getcoilcombos(op_averaging(out_raw),1);
            if ~isIMA
		        coilcombos_mrs	= op_getcoilcombos(op_averaging(out_raw),nPos_cc,'w');
            end

		    % Select coil phases and amplitudes for coil combination of MR spectra depending
		    % on available signals
		    if with_ref && ~isIMA
			    coilcombos	= coilcombos_ref_ECC;
		    else
			    if with_water
				    coilcombos	= coilcombos_w;
			    else
				    coilcombos	= coilcombos_mrs;
			    end			% % End of if with_water
		    end		% End of if with_ref
			    
		    % Combine coil channels before and after signal averaging for comparison and
		    % plotting
            if ~isIMA
		        [out_cc,fid_pre,spec_pre,ph,sig]	= op_addrcvrs(out_raw,nPos_cc,'w',coilcombos);
		        [out_av_cc,fid_av_pre,spec_av_pre]	= op_addrcvrs(op_averaging(out_raw),nPos_cc,'w',coilcombos);   
            end
            out_raw_av							= op_averaging(out_raw);
		    
		    % Generate unprocessed spectrum or spectra, respectively
            if ~isIMA
		        out_noproc			= op_averaging(out_cc);
            else
                out_cc				= out_raw;
                out_noproc			= op_averaging(out_cc);
                spec_av_pre			= out_noproc.specs;
            end
            
            if with_ref
                if ~isIMA
			        out_ref_ECC_noproc		= op_averaging(out_ref_ECC_cc);
			        out_ref_Quant_noproc	= op_averaging(out_ref_Quant_cc);
                else
                    out_ref_ECC_cc			= out_ref_ECC_raw;
                    out_ref_Quant_cc		= out_ref_Quant_raw;
			        out_ref_ECC_noproc		= op_averaging(out_ref_ECC_cc);
			        out_ref_Quant_noproc	= op_averaging(out_ref_Quant_cc);
                end
            end

		    if with_water
                if ~isIMA_w
			        out_w_noproc	= op_averaging(out_w_cc);
                else
                    out_w_cc		= out_w_raw;
                    out_w_noproc	= op_averaging(out_w_cc);
                    spec_w_pre		= out_w_raw.specs;
                end
		    end
    
		    % Generate plots showing coil channels before and after phase alignment
		    % Only display figure(s), if selected
		    if plotSwitch == 1
			    h1	= figure('position',[fig_left fig_bottom fig_width fig_height]);
		    else
			    h1	= figure('visible','off');
		    end
		    subplot(1,2,1);
		    plot(out_raw_av.ppm,real(out_raw_av.specs(:,:,1)));xlim([1 5]);
		    set(gca,'FontSize',8);
		    set(gca,'XDir','reverse');
		    xlabel('Frequency (ppm)','FontSize',10);
		    ylabel('Amplitude (a.u.)','FontSize',10);
		    title('Before correction','FontSize',12);
		    box off;
		    subplot(1,2,2);
		    plot(out_raw_av.ppm,real(spec_av_pre(:,:,1)));xlim([1 5]);
		    set(gca,'FontSize',12);
		    set(gca,'XDir','reverse');
		    xlabel('Frequency (ppm)','FontSize',10);
		    ylabel('Amplitude(a.u.)','FontSize',10);
		    title('After correction','FontSize',12);
		    box off;
		    set(h1,'PaperUnits','centimeters');
		    set(h1,'PaperPosition',[0 0 20 10]);
		    
		    % Save figures, if report switch is turned ON
		    if reportSwitch == 1
			    %saveas(h1,[outDirString nameSpec '/report/figs/coilReconFig'],'jpg');
			    %saveas(h1,[outDirString nameSpec '/report/figs/coilReconFig'],'fig');
			    saveas(h1,[outDirString reportFigDirStr 'coilReconFig'],'jpg');
			    saveas(h1,[outDirString reportFigDirStr 'coilReconFig'],'fig');
			    %close(h1);
		    end
		    
		    if with_water
			    % Only display figure(s), if selected
			    if plotSwitch == 1
				    h2	= figure('position',[fig_left (fig_bottom+fig_dist_b) fig_width fig_height]);
			    else
				    h2	= figure('visible','off');
			    end
			    subplot(2,1,1);
			    plot(out_w_raw.ppm,real(out_w_raw.specs(:,:,1,1)));xlim([4 5]);
			    xlabel('Frequency (ppm)');
			    ylabel('Amplitude (a.u.)');
			    title('Multi-channel (water unsupp.) data before phase correction');
			    subplot(2,1,2);
			    plot(out_w_raw.ppm,real(spec_w_pre(:,:,1,1)));xlim([4 5]);
			    xlabel('Frequency (ppm)');
			    ylabel('Amplitude (a.u.)');
			    title('Multi-channel (water unsupp.) data after phase correction');
			    box off;
			    set(h2,'PaperUnits','centimeters');
			    set(h2,'PaperPosition',[0 0 20 10]);
			    
			    % Save figures, if report switch is turned ON
			    if reportSwitch == 1
				    %saveas(h2,[outDirString nameSpec '/report/figs/w_coilReconFig'],'jpg');
				    %saveas(h2,[outDirString nameSpec '/report/figs/w_coilReconFig'],'fig');
				    saveas(h2,[outDirString reportFigDirStr 'w_coilReconFig'],'jpg');
				    saveas(h2,[outDirString reportFigDirStr 'w_coilReconFig'],'fig');
			    end
			    close(h2);
		    end
		    close(h1);
        else
            % Move "raw" DICOM (.IMA) data to coil combined "cc" data, as DICOM data are 
			% already coil combined, to execute same preprocssing code as for MRS data
			% that needs to be coil combined, e.g. MRS raw data (.dat)
            out_cc		= out_raw;
            out_noproc	= op_averaging(out_cc);
            if with_water
                out_w_cc				= out_w_raw;
                out_w_noproc			= op_averaging(out_w_cc);
            end
            if with_ref
                out_ref_ECC_cc			= out_ref_ECC_raw;
				out_ref_ECC_noproc		= op_averaging(out_ref_ECC_cc);
                out_ref_Quant_cc		= out_ref_Quant_raw;		    
			    out_ref_Quant_noproc	= op_averaging(out_ref_Quant_cc);
			end
		%end		% End of  if ~(isIMA && isIMA_w)
		end		% End of if ~(isIMA && isIMA_w) && (out_raw.dims.coils == 0)
		
        
		%% Remove bad averages from MRS data
		%%%%%%%% OPTIONAL REMOVAL OF BAD AVERAGES FROM DATASET %%%%%%%%%%%%%%%%%%%%
		close all;
		out_cc2			= out_cc;
		nBadAvgTotal	= 0;
		nbadAverages	= 1;
		rmbadav			= 'y';
		
		% Do not remove bad averages, if either not selected or if dimension of 
		% averages does not exist (index for dimension of averages = 0), 
		% e.g. when data are already averaged
		%if rmbadav=='n' || rmbadav=='N'
		if rmbadav=='n' || rmbadav=='N' || out_cc.dims.averages == 0
			out_rm		= out_cc;
			%nSD			= 'N/A';
		else
			sat		='n';
			while sat=='n' || sat=='N'
				%nSD=4; % Setting the number of standard deviations;
				iter			= 1;
				nbadAverages	= 1;
				nBadAvgTotal	= 0;
				out_cc2			= out_cc;
				while nbadAverages>0
					[out_rm,metric{iter},badAverages]	= op_rmbadaverages(out_cc2,nSD,'t');
					badAverages;
					nbadAverages	= length(badAverages);
					nBadAvgTotal	= nBadAvgTotal+nbadAverages;
					out_cc2			= out_rm;
					iter			= iter+1;
					disp([num2str(nbadAverages) ' bad averages removed on this iteration.']);
					disp([num2str(nBadAvgTotal) ' bad averages removed in total.']);
					close all;
				end
				
				% Create  figure to show pre-post removal of averages
				% Only display figure(s), if selected
				if plotSwitch == 1
					h3	= figure('position',[fig_left fig_bottom fig_width fig_height]);
				else
					h3	= figure('visible','off');
				end
				subplot(1,2,1);
				plot(out_cc.ppm,real(out_cc.specs(:,:)));xlim([1 5]);
				set(gca,'FontSize',8);
				set(gca,'XDir','reverse');
				xlabel('Frequency (ppm)','FontSize',10);
				ylabel('Amplitude(a.u.)','FontSize',10);
				title('Before','FontSize',12);
				box off;
				subplot(1,2,2);
				plot(out_rm.ppm,real(out_rm.specs(:,:)));xlim([1 5]);
				set(gca,'FontSize',8);
				set(gca,'XDir','reverse');
				xlabel('Frequency (ppm)','FontSize',10);
				ylabel('Amplitude(a.u.)','FontSize',10);
				title('After','FontSize',12);
				box off;
				set(h3,'PaperUnits','centimeters');
				set(h3,'PaperPosition',[0 0 20 15]);
				%saveas(h3,[outDirString nameSpec '/report/figs/rmBadAvg_prePostFig'],'jpg');
				%saveas(h3,[outDirString nameSpec '/report/figs/rmBadAvg_prePostFig'],'fig');
				%saveas(h3,[outDirString reportFigDirStr 'rmBadAvg_prePostFig'],'jpg');
				%saveas(h3,[outDirString reportFigDirStr 'rmBadAvg_prePostFig'],'fig');
				%close(h3);
				
				% Only display figure(s), if selected
				if plotSwitch == 1
					h4	= figure('position',[fig_left (fig_bottom+fig_dist_b) fig_width fig_height]);
				else
					h4	= figure('visible','off');
				end
				plot([1:length(metric{1})],metric{1},'.r',[1:length(metric{iter-1})],metric{iter-1},'x','MarkerSize',16);
				set(gca,'FontSize',8);
				xlabel('Scan Number','FontSize',10);
				ylabel('Deviation Metric','FontSize',10);
				legend('Before rmBadAv','After rmBadAv');
				legend boxoff;
				title('Deviation Metric','FontSize',12);
				box off;
				set(h4,'PaperUnits','centimeters');
				set(h4,'PaperPosition',[0 0 20 10]);
				%saveas(h4,[outDirString nameSpec '/report/figs/rmBadAvg_scatterFig'],'jpg');
				%saveas(h4,[outDirString nameSpec '/report/figs/rmBadAvg_scatterFig'],'fig');
					
				% Save figures, if report switch is turned ON
				if reportSwitch == 1
					saveas(h3,[outDirString reportFigDirStr 'rmBadAvg_prePostFig'],'jpg');
					saveas(h3,[outDirString reportFigDirStr 'rmBadAvg_prePostFig'],'fig');					
					saveas(h4,[outDirString reportFigDirStr 'rmBadAvg_scatterFig'],'jpg');
					saveas(h4,[outDirString reportFigDirStr 'rmBadAvg_scatterFig'],'fig');
				end
				close(h4);
				close(h3);
				
				%sat1=input('are you satisfied with the removal of bad averages? ','s');
				sat='y';
				
			end			% End of while sat=='n' || sat=='N'
			
			% Write a readme file to record the number of dropped avgs to output directory
			% for report instead of input data directory, if report switch is turned ON
			if reportSwitch == 1
				%fid1=fopen([dirString '/readme.txt'],'w+');
				%fid1		= fopen([outDirString outFileName '_readme.txt'],'w+');
				%fid1		= fopen([outDirString reportDirStr 'readme.txt'],'w+');
                %disp([outDirString reportDirStr outFileName '_readme.txt']);
				fid1		= fopen([outDirString reportDirStr outFileName '_readme.txt'],'w+');
				fprintf(fid1,'Original number of averages: \t%5.6f',out_raw.sz(out_raw.dims.averages));
				disp(['Original number of averages:  ' num2str(out_raw.sz(out_raw.dims.averages))]);
				fprintf(fid1,'\nNumber of bad Averages removed:  \t%5.6f',nBadAvgTotal);
				disp(['Number of bad averages removed:  ' num2str(nBadAvgTotal)]);
				fprintf(fid1,'\nNumber of remaining averages in processed dataset:  \t%5.6f',out_rm.sz(out_rm.dims.averages));
				disp(['Number of remaining averages in processed dataset:  ' num2str(out_rm.sz(out_rm.dims.averages))]);
				fclose(fid1);
			end
		end		% End of if rmbadav=='n' || rmbadav=='N' || out_cc.dims.averages == 0
				
		%%%%%%%%%%%%%%%%%%%% END OF BAD AVERAGES REMOVAL %%%%%%%%%%%%%%%%%%%%
		
		
		%% NOW ALIGN AVERAGES:  A.K.A. Frequency Drift Correction
		% If minimum user input selected, frequency drift correction is not optional
		% Ask for user input for frequency correction only, if minimum user input is NOT
		% selected
		% Set parameters for drift correction depending on type of data, i.e. whether MRS
		% data is spectrum or water signal
		% NOTE: Check whether aligning of averages in frequency domain works, if the MR
		% spectrum is water signal itself; if not, simply align averages in time domain 
		%if( strcmp(dataType, 'mrs_w') || strcmp(dataType, 'mrs') )
		switch dataType
			case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
				% MR spectrum is provided together without or with unsuppressed water 
				% signal and/or with reference scans
				ppmmin_fix		= 1.6;
				%ppmmaxarray_fix	= [3.5; 4.0; 5.5];
				ppmmaxarray_fix = [2.4,2.85,3.35,4.2,4.4,5.2];
				iamax			= 6;
			case {'water', 'water_ref'}
				% MR spectrum is water signal itself without or with reference scans
				ppmmin_fix		= 4.2;
				ppmmaxarray_fix	= [5.5 5.5 5.2];
				iamax			= 6;
			
			otherwise
				error('%s: Unknown MRS dataType = %s!', sFunctionName, dataType);
		end		% End of switch dataType
		
		% Do not perform drift correction, if either not selected or if dimension of 
		% averages does not exist (index for dimension of averages = 0), 
		% e.g. when data is already averaged
		driftCorr		= 'y';
		%if driftCorr=='n' || driftCorr=='N'
		if driftCorr=='n' || driftCorr=='N' || out_rm.dims.averages == 0
			out_av		= op_averaging(out_rm);
			if with_water
				out_w_av			= op_averaging(out_w_cc);
			end
			if with_ref
				out_ref_ECC_av		= op_averaging(out_ref_ECC_cc);
				out_ref_Quant_av	= op_averaging(out_ref_Quant_cc);
			end
			fs			= 0;
			phs			= 0;
		else
			if with_water
				%out_w_aa		= op_alignAverages(out_w_cc,tmaxin,'n');
				out_w_aa			= op_alignAverages(out_w_cc,0.2,'n');
			end
			if with_ref
				out_ref_ECC_aa		= op_alignAverages(out_ref_ECC_cc,0.2,'n'); 
				out_ref_Quant_aa	= op_alignAverages(out_ref_Quant_cc,0.2,'n'); 
			end
			sat			= 'n';
			out_rm2		= out_rm;
			while sat=='n' || sat=='N'
				fsPoly		= 100;
				phsPoly		= 1000;
				fscum		= zeros(out_rm2.sz(out_rm2.dims.averages),1);
				phscum		= zeros(out_rm2.sz(out_rm2.dims.averages),1);
				iter		= 1;
				while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<iterin
					%iter			= iter+1
					close all
					%tmax			= 0.25+0.03*randn(1);
					%ppmmin			= 1.6+0.1*randn(1);
					%ppmmaxarray	= [3.5+0.1*randn(1,2),4+0.1*randn(1,3),5.5+0.1*randn(1,1)];
					%ppmmax			= ppmmaxarray(randi(6,1));
					tmax			= tmaxin+0.03*randn(1);
					ppmmin			= ppmmin_fix+0.1*randn(1);
					ppmmaxarray		= [ppmmaxarray_fix(1)+0.1*randn(1,2),ppmmaxarray_fix(2)+0.1*randn(1,3),ppmmaxarray_fix(3)+0.1*randn(1,1)];
					ppmmax			= ppmmaxarray(randi(iamax,1));
					switch aaDomain
						case 't'
							[out_aa,fs,phs]		= op_alignAverages(out_rm2,tmax,'y');
							%[out_aa,fs,phs]		= op_alignAverages(out_rm2);
						case 'f'
							[out_aa,fs,phs]		= op_alignAverages_fd(out_rm2,ppmmin,ppmmax,tmax,'y');
							%[out_aa,fs,phs]		= op_alignAverages_fd(out_rm2,ppmmin,ppmmax,tmax,'n');
						otherwise
							error('%s: ERROR: avgAlignDomain %s not recognized!', sFunctionName, aaDomain);
					end
					
					fsPoly		= polyfit([1:out_aa.sz(out_aa.dims.averages)]',fs,1)
					phsPoly		= polyfit([1:out_aa.sz(out_aa.dims.averages)]',phs,1)
					%iter
					%disp( sprintf('Aligning averages iteration %d', iter) );
					fprintf(1, 'Aligning averages iteration %d\n', iter) ;
					
					fscum		= fscum+fs;
					phscum		= phscum+phs;
					
					if driftCorr=='y' || driftCorr=='Y'
						out_rm2		= out_aa;
					end
					iter			= iter+1;
				end		% End of while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<iterin
				
				% Only display figure(s), if selected
				if plotSwitch == 1
					h5	= figure('position',[fig_left fig_bottom fig_width fig_height]);
				else
					h5	= figure('visible','off');
				end
				subplot(1,2,1);
				plot(out_rm.ppm,real(out_rm.specs(:,:)));xlim([1 5]);
				set(gca,'FontSize',8);
				set(gca,'XDir','reverse');
				xlabel('Frequency (ppm)','FontSize',10);
				ylabel('Amplitude(a.u.)','FontSize',10);
				title('Before','FontSize',12);
				box off;
				subplot(1,2,2);
				plot(out_aa.ppm,real(out_aa.specs(:,:)));xlim([1 5]);
				set(gca,'FontSize',8);
				set(gca,'XDir','reverse');
				xlabel('Frequency (ppm)','FontSize',10);
				ylabel('Amplitude(a.u.)','FontSize',10);
				title('After','FontSize',12);
				box off;
				set(h5,'PaperUnits','centimeters');
				set(h5,'PaperPosition',[0 0 20 15]);
				%saveas(h5,[outDirString nameSpec '/report/figs/alignAvgs_prePostFig'],'jpg');
				%saveas(h5,[outDirString nameSpec '/report/figs/alignAvgs_prePostFig'],'fig');
				%saveas(h5,[outDirString reportFigDirStr 'alignAvgs_prePostFig'],'jpg');
				%saveas(h5,[outDirString reportFigDirStr 'alignAvgs_prePostFig'],'fig');
				%close(h5);
				
				% Only display figure(s), if selected
				if plotSwitch == 1
					h6	= figure('position',[fig_left (fig_bottom+fig_dist_b) fig_width fig_height]);
				else
					h6	= figure('visible','off');
				end
				plot([1:out_aa.sz(out_aa.dims.averages)],fscum,'.-','LineWidth',2);
				set(gca,'FontSize',8);
				xlabel('Scan Number','FontSize',10);
				ylabel('Frequency Drift [Hz]','FontSize',10);
				box off;
				legend('Frequency Drift','Location','SouthEast');
				legend boxoff;
				title('Estimated Frequency Drift','FontSize',12);
				set(h6,'PaperUnits','centimeters');
				set(h6,'PaperPosition',[0 0 10 10]);
				%saveas(h6,[outDirString nameSpec '/report/figs/freqDriftFig'],'jpg');
				%saveas(h6,[outDirString nameSpec '/report/figs/freqDriftFig'],'fig');
				%saveas(h6,[outDirString reportFigDirStr 'freqDriftFig'],'jpg');
				%saveas(h6,[outDirString reportFigDirStr 'freqDriftFig'],'fig');
				%close(h6);
			
				% Only display figure(s), if selected
				if plotSwitch == 1
					h7	= figure('position',[(fig_left+fig_dist_l) fig_bottom fig_width fig_height]);
				else
					h7	= figure('visible','off');
				end
				plot([1:out_aa.sz(out_aa.dims.averages)],phscum,'.-','LineWidth',2);
				set(gca,'FontSize',8);
				xlabel('Scan Number','FontSize',10);
				ylabel('Phase Drift [Deg.]','FontSize',10);
				box off;
				legend('Phase Drift','Location','SouthEast');
				legend boxoff;
				title('Estimated Phase Drift','FontSize',12);
				set(h7,'PaperUnits','centimeters');
				set(h7,'PaperPosition',[0 0 10 10]);
				%saveas(h7,[outDirString nameSpec '/report/figs/phaseDriftFig'],'jpg');
				%saveas(h7,[outDirString nameSpec '/report/figs/phaseDriftFig'],'fig');
				
				% Save figures, if report switch is turned ON
				if reportSwitch == 1
					saveas(h5,[outDirString reportFigDirStr 'alignAvgs_prePostFig'],'jpg');
					saveas(h5,[outDirString reportFigDirStr 'alignAvgs_prePostFig'],'fig')
					saveas(h6,[outDirString reportFigDirStr 'freqDriftFig'],'jpg');
					saveas(h6,[outDirString reportFigDirStr 'freqDriftFig'],'fig');
					saveas(h7,[outDirString reportFigDirStr 'phaseDriftFig'],'jpg');
					saveas(h7,[outDirString reportFigDirStr 'phaseDriftFig'],'fig');
				end
				close(h7);
				close(h6);
				close(h5);
				
				sat='y';
				if sat=='n'
					iter		= 0;
					p1			= 100;
					fscum		= zeros(out_rm.sz(2:end));
					phscum		= zeros(out_rm.sz(2:end));
					fs2cum		= zeros(out_cc.sz(2:end));
					phs2cum		= zeros(out_cc.sz(2:end));
					out_rm2		= out_rm;
					out_cc2		= out_cc;
				end
				totalFreqDrift		= mean(max(fscum)-min(fscum));
				totalPhaseDrift		= mean(max(phscum)-min(phscum));
				close all
			end		% End of while sat=='n' || sat=='N'
			
			% Now average the aligned averages
			out_av		= op_averaging(out_aa);
			if with_water
				out_w_av			= op_averaging(out_w_aa);
			end
			if with_ref
				out_ref_ECC_av		= op_averaging(out_ref_ECC_aa);
				out_ref_Quant_av	= op_averaging(out_ref_Quant_aa); 
			end
		end		% End of if driftCorr=='n' || driftCorr=='N' || out_rm.dims.averages == 0
		
		
		%% Combine any subspectra, if existent
		% If at this stage of processing, subspectra are still present, i.e. index to 
		% corresponding dimension in data struct is non-zero, combine subspectra to
		% avoid any error when writing processing results to file
		% Usually, susbspectra should be comnbined before processing of averages, but if
		% different type of input data is unwantingly interpreted as subspectra by input 
		% loading routine, then these can be still present at this point
		if out_av.dims.subSpecs ~= 0
			disp(sMsg_newLines);
			warning('%s: Subspectra still present in spectrum! # of subsprectra = out_av.sz(out_av.dims.subSpecs) = %d', sFunctionName, out_av.sz(out_av.dims.subSpecs));
			disp(sMsg_newLine);
			disp('Averaging subSpecs of spectrum ...');
			disp(sMsg_newLines);
			out_av_tmp				= out_av;
			out_av					= op_average_subSpecs_s(out_av_tmp);
			%out_av		= op_combinesubspecs(out_av_tmp, 'diff');
			clear out_av_tmp;
			
			% Since unprocessed spectra will also be written to file in LCModel format,
			% the unprocessed subSpectra have to be averaged as well to avoid any errors
			% in output routine io_writelcm.m due to isISIS flag being still set to 1
			if out_noproc.dims.subSpecs ~= 0
				out_noproc_tmp			= out_noproc;
				out_noproc				= op_average_subSpecs_s(out_noproc_tmp);
				clear out_noproc_tmp;
			end		% End of if out_noproc.dims.subSpecs ~= 0
		end		% End of if out_av.dims.subSpecs ~= 0
		
		if with_water
			if out_w_av.dims.subSpecs ~= 0
				disp(sMsg_newLines);
				warning('%s: Subspectra still present in unsuppressed water signal! # of subsprectra = out_w_av.sz(out_w_av.dims.subSpecs) = %d', sFunctionName, out_w_av.sz(out_w_av.dims.subSpecs));
				disp(sMsg_newLine);
				disp('Averaging subSpecs of water signal ...');
				disp(sMsg_newLines);
				out_w_av_tmp	= out_w_av;
				out_w_av		= op_average_subSpecs_s(out_w_av_tmp);
				clear out_w_av_tmp;
				
				% Since unprocessed water signals will also be written to file in LCModel 
				% format, the unprocessed subSpectra have to be averaged as well to avoid 
				% any errors in output routine io_writelcm.m due to isISIS flag being 
				% still set to 1
				if out_w_noproc.dims.subSpecs ~= 0
					out_w_noproc_tmp		= out_w_noproc;
					out_w_noproc			= op_average_subSpecs_s(out_w_noproc_tmp);
					clear out_w_noproc_tmp;
				end		% End of if out_w_noproc.dims.subSpecs ~= 0
			end		% End of if out_w_av.dims.subSpecs ~= 0
		end		% End of if with_water
		
		if with_ref
			% Though it is highly unlikely that for only one type of reference scans
			% subspectra are present and not for the other one, both types of reference
			% scans from Dinesh's CMRR sequence are treated independently to be on the
			% safe side
			if out_ref_ECC_av.dims.subSpecs ~= 0
				disp(sMsg_newLines);
				warning('%s: Subspectra still present in reference (water) scans for ECC! # of subsprectra = out_ref_ECC_av.sz(out_ref_ECC_av.dims.subSpecs) = %d', sFunctionName, out_ref_ECC_av.sz(out_ref_ECC_av.dims.subSpecs));
				disp(sMsg_newLine);
				disp('Averaging subSpecs of reference (water) scans for ECC ...');
				disp(sMsg_newLines);
				out_ref_ECC_av_tmp		= out_ref_ECC_av;
				out_ref_ECC_av			= op_average_subSpecs_s(out_ref_ECC_av_tmp);
				clear out_ref_ECC_av_tmp;
				
				% Since unprocessed refrence scans will also be written to file in LCModel
				% format, the unprocessed subSpectra have to be averaged as well to avoid
				% any errors in output routine io_writelcm.m due to isISIS flag being
				% still set to 1
				if out_ref_ECC_noproc.dims.subSpecs ~= 0
					out_ref_ECC_noproc_tmp		= out_ref_ECC_noproc;
					out_ref_ECC_noproc			= op_average_subSpecs_s(out_ref_ECC_noproc_tmp);
					clear out_ref_ECC_noproc_tmp;
				end		% End of if out_ref_ECC_noproc.dims.subSpecs ~= 0
			end		% End of if out_ref_ECC_av.dims.subSpecs ~= 0
			
			if out_ref_Quant_av.dims.subSpecs ~= 0
				disp(sMsg_newLines);
				warning('%s: Subspectra still present in reference (water) scans for Quant! # of subsprectra = out_ref_Quant_av.sz(out_ref_Quant_av.dims.subSpecs) = %d', sFunctionName, out_ref_Quant_av.sz(out_ref_Quant_av.dims.subSpecs));
				disp(sMsg_newLine);
				disp('Averaging subSpecs of reference (water) scans for Quant ...');
				disp(sMsg_newLines);
				out_ref_Quant_av_tmp	= out_ref_Quant_av;
				out_ref_Quant_av		= op_average_subSpecs_s(out_ref_Quant_av_tmp);
				clear out_ref_Quant_av_tmp;
				
				% Since unprocessed refrence scans will also be written to file in LCModel
				% format, the unprocessed subSpectra have to be averaged as well to avoid
				% any errors in output routine io_writelcm.m due to isISIS flag being
				% still set to 1
				if out_ref_Quant_noproc.dims.subSpecs ~= 0
					out_ref_Quant_noproc_tmp	= out_ref_Quant_noproc;
					out_ref_Quant_noproc		= op_average_subSpecs_s(out_ref_Quant_noproc_tmp);
					clear out_ref_Quant_noproc_tmp;
				end		% End of if out_ref_Quant_noproc.dims.subSpecs ~= 0
			end		% End of if out_ref_Quant_av.dims.subSpecs ~= 0
			
		end		% End of if with_ref
		
		
		%% Perform eddy current correction (ECC), if selected
		if bECC		
			% Copy struct of MR spectrum for ECC and init signal used for ECC
			out_av_ECC_In			= out_av;
			out_water_av_ECC_In		= struct([]);
			% Use reference (water) signals for ECC, if acquired
			% If not, then use an unsuppressed water signal, if acquired
			% If no reference and no water signals are acquired, check whether MR spectrum
			% is water signal itself; and if it is, use it for ECC
			% (Same calls of ECC routine are used with same input variables, but different
			% output arguments to facilitate correct assignment of output variables
			% for further processing depending on which if case is actually invoked)
			if with_ref
				out_water_av_ECC_In			= out_ref_ECC_av;
				[out_av, out_ref_ECC_av]	= op_ecc_s(out_av_ECC_In, out_water_av_ECC_In);
				% Perform ECC for water reference signals for quantification using these
				% signals themselves
				% (output signals are then also both the same, since ECC is performed on
				% both input signals;
				% use different output variables to avoid confusion, though  
				% second one is unused)) 
				out_ref_Quant_av_ECC_In			= out_ref_Quant_av;
				out_ref_Quant_water_av_ECC_In	= out_ref_Quant_av;
				[out_ref_Quant_av, out_ref_Quant_av_ECC_w]	= op_ecc_s(out_ref_Quant_av_ECC_In, out_ref_Quant_water_av_ECC_In);
			else
				if with_water
					% If no reference signals, use water signals for ECC, if acquired
					out_water_av_ECC_In		= out_w_av;
					[out_av, out_w_av]		= op_ecc_s(out_av_ECC_In, out_water_av_ECC_In);
				else
					if strcmp(dataType, 'water')
						% (dataType = 'water_ref' is covered by 'with_ref')
						% MR spectrum (= water signal itself) also used for ECC
						% (Input and output arguments are both the same, respectively; 
						% use different output variables to avoid confusion, though  
						% second one is unused)) 
						out_water_av_ECC_In		= out_av;
						[out_av, out_av_ECC_w]	= op_ecc_s(out_av_ECC_In, out_water_av_ECC_In);
					else
						% No reference and no water signals and MR spectrum is not water
						% signal itself => ECC not possible
						error('%s: No reference and no water signals and MR spectrum is not water signal itself (dataType = %s) => ECC not possible!', sFunctionName, dataType);
					end		% End of if strcmp(dataType, 'water')
				end		% End of if with_water
			end		% End of if with_ref
		end		% End of if bECC
		
		
		%% Perform phase correction and frequency shifting, if selected
		% Set parameters for phase correction and frequency shifting inluding reference 
		% frequencies and factors for zeropadding
		freqs_Cr		= [2.9; 3.1; 3.027];
		freqs_w			= [4.0; 5.5; 4.65];
		zp_factor		= 16;
		zp_factor_w		= 16;
		zp_factor_ref	= 16;
		
		if bPhaseCorrFreqShift
			% Create empty figure, so that op_addphase(...) and other routines plot to 
			% this figure and do not overwrite previous plots (e.g. from ECC)
			h_figTmp1	= figure;
			
			% Left shift the MRS signals to eliminate first order phase, perform 
			% zero-order phase correction, and shift data in frequency to obtain 
			% reference peaks at known positions
			% Now left shift, if not already done before
			% MR spectrum
			if out_av.flags.leftshifted == 0
				out_ls					= op_leftshift(out_av,out_av.pointsToLeftshift);
			else
				out_ls					= out_av;
			end
			
			if with_water
				% Water signal
				if out_w_av.flags.leftshifted == 0
					out_w_ls			= op_leftshift(out_w_av,out_w_av.pointsToLeftshift);
				else
					out_w_ls			= out_w_av;
				end
			end		% End of if with_water
			
			if with_ref
				% Reference (water) signals
				if out_ref_ECC_av.flags.leftshifted == 0
					out_ref_ECC_ls		= op_leftshift(out_ref_ECC_av,out_ref_ECC_av.pointsToLeftshift);
				else
					out_ref_ECC_ls		= out_ref_ECC_av;
				end
				if out_ref_Quant_av.flags.leftshifted == 0
					out_ref_Quant_ls	= op_leftshift(out_ref_Quant_av,out_ref_Quant_av.pointsToLeftshift);
				else
					out_ref_Quant_ls	= out_ref_Quant_av;
				end
			end		% End of if with_ref
			
			% Perform automatic zero-order phase correction
			% If data is MR spectrum, use creatine peak at 3.027 ppm:
			%if( strcmp(dataType, 'mrs_w') || strcmp(dataType, 'mrs') )
			switch dataType
				case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
					% MR spectrum is provided together without or with unsuppressed water
					% signal and/or with reference scans
					out_ls_zp		= op_zeropad(out_ls,zp_factor);
					%[out_ph,ph0]	= op_autophase(out_ls,2.9,3.1);
					[out_ph,ph0]	= op_autophase(out_ls,freqs_Cr(1),freqs_Cr(2));
					out_ls_zp		= op_addphase(out_ls_zp,ph0);
					% And now for water unsuppressed data (use water peak):
					if with_water
						out_w_ls_zp	= op_zeropad(out_w_ls,zp_factor_w);
						%indexw		= find(abs(out_w_ls_zp.specs) == max(abs(out_w_ls_zp.specs(out_w_ls_zp.ppm>4 & out_w_ls_zp.ppm<5.5))));
						index_w		= find(abs(out_w_ls_zp.specs) == max(abs(out_w_ls_zp.specs(out_w_ls_zp.ppm>freqs_w(1) & out_w_ls_zp.ppm<freqs_w(2)))));
						ph0_w		= -phase(out_w_ls_zp.specs(index_w))*180/pi;
						out_w_ph	= op_addphase(out_w_ls,ph0_w);
						out_w_ls_zp	= op_addphase(out_w_ls_zp,ph0_w);
					end
				case {'water', 'water_ref'}
					% MR spectrum is water signal itself without or with reference scans
					out_ls_zp	= op_zeropad(out_ls,zp_factor);
					index_w		= find(abs(out_ls_zp.specs) == max(abs(out_ls_zp.specs(out_ls_zp.ppm>freqs_w(1) & out_ls_zp.ppm<freqs_w(2)))));
					ph0			= -phase(out_ls_zp.specs(index_w))*180/pi;
					out_ph		= op_addphase(out_ls,ph0);
					out_ls_zp	= op_addphase(out_ls_zp,ph0);
					
				otherwise
					error('%s: Unknown MRS dataType = %s!', sFunctionName, dataType);
			end		% End of switch dataType
			
			% Reference scans can be processed the same way for all cases of the above 
			% switch statement assuming that these are always water signals
			if with_ref
				% Reference (water) signals
				out_ref_ECC_ls_zp	= op_zeropad(out_ref_ECC_ls,zp_factor_ref);
				index_ref_ECC		= find(abs(out_ref_ECC_ls_zp.specs) == max(abs(out_ref_ECC_ls_zp.specs(out_ref_ECC_ls_zp.ppm>freqs_w(1) & out_ref_ECC_ls_zp.ppm<freqs_w(2)))));
				ph0_ref_ECC			= -phase(out_ref_ECC_ls_zp.specs(index_ref_ECC))*180/pi;
				out_ref_ECC_ph		= op_addphase(out_ref_ECC_ls,ph0_ref_ECC);
				out_ref_ECC_ls_zp	= op_addphase(out_ref_ECC_ls_zp,ph0_ref_ECC);
				
				out_ref_Quant_ls_zp	= op_zeropad(out_ref_Quant_ls,zp_factor_ref);
				index_ref_Quant		= find(abs(out_ref_Quant_ls_zp.specs) == max(abs(out_ref_Quant_ls_zp.specs(out_ref_Quant_ls_zp.ppm>freqs_w(1) & out_ref_Quant_ls_zp.ppm<freqs_w(2)))));
				ph0_ref_Quant		= -phase(out_ref_Quant_ls_zp.specs(index_ref_Quant))*180/pi;
				out_ref_Quant_ph	= op_addphase(out_ref_Quant_ls,ph0_ref_Quant);
				out_ref_Quant_ls_zp	= op_addphase(out_ref_Quant_ls_zp,ph0_ref_Quant);
			end		% End of if with_ref
			
			% Perform same phase corection (left shift and zero-order) on unprocessed data
			% Left shift only, if not already done before
			% MR spectrum
			if out_noproc.flags.leftshifted == 0
				out_noproc					= op_addphase(op_leftshift(out_noproc,out_noproc.pointsToLeftshift),ph0);
			else
				out_noproc					= op_addphase(out_noproc,ph0);
			end
			if with_water
				% Water signal
				if out_w_noproc.flags.leftshifted == 0
					out_w_noproc			= op_addphase(op_leftshift(out_w_noproc,out_w_noproc.pointsToLeftshift),ph0_w);
				else
					out_w_noproc			= op_addphase(out_w_noproc,ph0_w);
				end
			end		% End of if with_water
			
			if with_ref
				% Reference (water) signals
				if out_ref_ECC_noproc.flags.leftshifted == 0
					out_ref_ECC_noproc		= op_addphase(op_leftshift(out_ref_ECC_noproc,out_ref_ECC_noproc.pointsToLeftshift),ph0_ref_ECC);
				else
					out_ref_ECC_noproc		= op_addphase(out_ref_ECC_noproc,ph0_ref_ECC);
				end
				if out_ref_Quant_noproc.flags.leftshifted == 0
					out_ref_Quant_noproc	= op_addphase(op_leftshift(out_ref_Quant_noproc,out_ref_Quant_noproc.pointsToLeftshift),ph0_ref_Quant);
				else
					out_ref_Quant_noproc	= op_addphase(out_ref_Quant_noproc,ph0_ref_Quant);
				end
			end		% End of if with_ref
			
			% Frequency shift all spectra
			% If data is MR spectrum, so that creatine appears at 3.027 ppm
			%if( strcmp(dataType, 'mrs_w') || strcmp(dataType, 'mrs') )
			switch dataType
				case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
					% MR spectrum is provided together without or with unsuppressed water
					% signal and/or with reference scans
					%[~,frqShift]	= op_ppmref(out_ls_zp,2.9,3.1,3.027);
					[~,frqShift]	= op_ppmref(out_ls_zp,freqs_Cr(1),freqs_Cr(2),freqs_Cr(3));
					out				= op_freqshift(out_ph,frqShift);
					out_noproc		= op_freqshift(out_noproc,frqShift);
					% For water unsuppressed data, use water peak and set it to 4.65 ppm
					if with_water
						%[~,frqShift_w]	= op_ppmref(out_w_ls_zp,4,5.5,4.65);
						[~,frqShift_w]	= op_ppmref(out_w_ls_zp,freqs_w(1),freqs_w(2),freqs_w(3));
						out_w			= op_freqshift(out_w_ph,frqShift_w);
						out_w_noproc	= op_freqshift(out_w_noproc,frqShift_w);
					end		% End of if with_water
				case {'water', 'water_ref'}
					% MR spectrum is water signal itself without or with reference scans,
					% so use water peak and set it to 4.65 ppm
					%[~,frqShift]	= op_ppmref(out_ls_zp,4,5.5,4.65);
					[~,frqShift]	= op_ppmref(out_ls_zp,freqs_w(1),freqs_w(2),freqs_w(3));
					out				= op_freqshift(out_ph,frqShift);
					out_noproc		= op_freqshift(out_noproc,frqShift);
					
				otherwise
					error('%s: Unknown MRS dataType = %s!', sFunctionName, dataType);
			end		% End of switch dataType
			
			% Reference scans can be processed the same way for all cases of the above 
			% switch statement assuming that these are always water signals
			if with_ref
				% Reference (water) signals
				[~,frqShift_ref_ECC]	= op_ppmref(out_ref_ECC_ls_zp,freqs_w(1),freqs_w(2),freqs_w(3));
				out_ref_ECC				= op_freqshift(out_ref_ECC_ph,frqShift_ref_ECC);
				out_ref_ECC_noproc		= op_freqshift(out_ref_ECC_noproc,frqShift_ref_ECC);
				
				[~,frqShift_ref_Quant]	= op_ppmref(out_ref_Quant_ls_zp,freqs_w(1),freqs_w(2),freqs_w(3));
				out_ref_Quant			= op_freqshift(out_ref_Quant_ph,frqShift_ref_Quant);
				out_ref_Quant_noproc	= op_freqshift(out_ref_Quant_noproc,frqShift_ref_Quant);
			end		% End of if with_ref
			% Close selected figure(s)
			close(h_figTmp1);
		else	% No phase correction and no frequency shifting performed
			% Assign output struct names to struct names after last processing step to be 
			% able to continue with same code below for further processing
			% (not needed for unprocessed data, since struct names remain the same during
			% phase correction and frequency shifting)
			% MR spectrum
			out					= out_av;
			if with_water
				% Water signal
				out_w			= out_w_av;
			end		% End of if with_water
			
			if with_ref
				% Reference (water) signals
				out_ref_ECC		= out_ref_ECC_av;
				out_ref_Quant	= out_ref_Quant_av;
			end		% End of if with_ref
		end		% End of if bPhaseCorrFreqShift
		
						
		%% Display final MR spectra and save corresponding figures
		% Close all figures
		%close all;
		% Close selected figure(s)
		%close(h_figTmp1);
		
		%h=figure('visible','off');
		%plot(out.ppm,real(out.specs),'linewidth',2);xlim([0.2 5.2]);
		%xlabel('Frequency (ppm)','FontSize',10);
		%ylabel('Amplitude(a.u.)','FontSize',10);
		
		% Set figure parameters
		% Set plotting resolution and properties of axes depending on data type
		resolution		= 600;
		%if( strcmp(dataType, 'mrs_w') || strcmp(dataType, 'mrs') )
		switch dataType
			case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
				% MR spectrum is provided together without or with unsuppressed water
				% signal and/or with reference scans
				xLimValues1		= [0.0 5.5];
				xLimValues2		= [0.2 4.2];
				xTickValues2	= [0.5:0.5:4.0];
				strTitle_mrs	= sprintf('Preprocessed MR Spectrum');
			case {'water', 'water_ref'}
				% MR spectrum is water signal itself without or with reference scans,
				xLimValues1		= [3.3 5.9];
				xLimValues2		= [4.2 5.1];
				xTickValues2	= [4.2:0.2:5.0];
				strTitle_mrs	= sprintf('Preprocessed Water Spectrum');
				
			otherwise
				error('%s: Unknown MRS dataType = %s!', sFunctionName, dataType);
		end		% End of switch dataType
		
		h_mrs			= figure('visible','on');
		plot(out.ppm,real(out.specs),'linewidth',1.5);xlim(xLimValues1);
		set(gca,'FontSize',12, 'FontWeight','bold');
		set(gca,'XDir','reverse');
		set(gca,'XAxisLocation', 'origin');
		xlabel('ppm','FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle');
		ylabel('Amplitude(a.u.)','FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','baseline');
		box off;
		%title('Result: Preprocessed Spectrum','FontSize',12);
		title(strTitle_mrs,'FontSize',12);
		xf1		= xLimValues1(1) + (xLimValues1(2) - xLimValues1(1))/2;
		yf1		= 1.0*min(get(gca, 'ylim'));
		set(get(gca,'XLabel'),'Position', [xf1, yf1], 'VerticalAlignment', 'Top');
		
		% Save figure of spectrum
		% Create name for figure files
		% Extracting the digits before and after decimal point assumes that there is only
		% one digit after the decimal point
		% fix, i.e. rounding towards zero, is used to correctly handle negative numbers
		% Using round for digits 2 and 4 avoids getting zero if difference is 0.1
		%strFigName_add	= '_processed_lcm_0_0_5_5ppm';
		%digits = [fix(xLimValues1(1)) fix(abs(xLimValues1(1)-fix(xLimValues1(1)))*10) fix(xLimValues1(2)) fix(abs(xLimValues1(2)-fix(xLimValues1(2)))*10)];
		digits = [fix(xLimValues1(1)) round(abs(xLimValues1(1)-fix(xLimValues1(1)))*10) fix(xLimValues1(2)) round(abs(xLimValues1(2)-fix(xLimValues1(2)))*10)];
		strFigName_add	= sprintf( '_processed_lcm_%d_%d_%d_%dppm', digits(1), digits(2), digits(3), digits(4) );
		figureName_fig	= [outFileName, strFigName_add, '.fig'];
		figureName_png	= [outFileName, strFigName_add, '.png'];
		saveFigure_s(h_mrs, outDirString, figureName_fig, 'fig', resolution);
		saveFigure_s(h_mrs, outDirString, figureName_png, 'png', resolution);
		
		% Copy and modify figure to show different spectral range and save modified figure
		% as well
		% For figure creation, i.e. if spectrum should be used as figure for display or in 
		% paper, use only ppm range from 0.2 to 4.2 and 'whiten' y-axis
		% For regular use, leave y-axis as is to indicate signal strength
		ax_h_mrs		= gca;
		h_mrs2			= figure('visible','on');
		ax_h_mrs2		= copyobj(ax_h_mrs,h_mrs2);
		%set(gca, 'XLim', xLimValues, 'XTick',[0.5:0.5:4.0], 'XTickLabel',[0.5:0.5:4.0], 'YTickLabel',[], 'YColor',[1 1 1]);
		%set(gca, 'XLim', xLimValues2, 'XTick',[0.5:0.5:4.0], 'XTickLabel',[0.5:0.5:4.0]);
		set(gca, 'XLim', xLimValues2, 'XTick',xTickValues2, 'XTickLabel',xTickValues2);
		xf2		= xLimValues2(1) + (xLimValues2(2) - xLimValues2(1))/2;
		yf2		= 1.0*min(get(gca, 'ylim'));
		set(get(gca,'XLabel'),'Position', [xf2, yf2], 'VerticalAlignment', 'Top');
		% Create filenames for saving of modified figure
		%strFigName_add	= '_processed_lcm_0_2_4_2ppm';
		%digits = [fix(xLimValues2(1)) fix(abs(xLimValues2(1)-fix(xLimValues2(1)))*10) fix(xLimValues2(2)) fix(abs(xLimValues2(2)-fix(xLimValues2(2)))*10)];
		digits = [fix(xLimValues2(1)) round(abs(xLimValues2(1)-fix(xLimValues2(1)))*10) fix(xLimValues2(2)) round(abs(xLimValues2(2)-fix(xLimValues2(2)))*10)];
		strFigName_add	= sprintf( '_processed_lcm_%d_%d_%d_%dppm', digits(1), digits(2), digits(3), digits(4) );
		figureName_fig	= [outFileName, strFigName_add, '.fig'];
		figureName_png	= [outFileName, strFigName_add, '.png'];
		saveFigure_s(h_mrs2, outDirString, figureName_fig, 'fig', resolution);
		saveFigure_s(h_mrs2, outDirString, figureName_png, 'png', resolution);
		
		
		%% Display water signals and save correspoding figures, if available
		xLimValues3		= [3.3 5.9];
		if with_water
			% Show and save water signal
			h_w				= figure('visible','on');
			plot(out_w.ppm,real(out_w.specs),'linewidth',1.5);xlim(xLimValues3);
			axes_h_w		= get(h_w,'CurrentAxes');
			set(axes_h_w,'FontSize',12, 'FontWeight','bold');
			set(axes_h_w,'XDir','reverse');
			set(axes_h_w,'XAxisLocation', 'origin');
			xlabel('ppm','FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle');
			ylabel('Amplitude(a.u.)','FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','baseline');
			box off;
			title('Preprocessed Water Signal','FontSize',12);
			xf3		= xLimValues3(1) + (xLimValues3(2) - xLimValues3(1))/2;
			yf3		= 1.0*min(get(axes_h_w, 'ylim'));
			set(get(axes_h_w,'XLabel'),'Position', [xf3, yf3]);
			% Save figure of water signal
			% Create filenames for saving of processed output
			%strFigName_addW	= '_w_processed_lcm_3_3_5_9ppm';
			digits = [fix(xLimValues2(1)) fix(abs(xLimValues3(1)-fix(xLimValues3(1)))*10) fix(xLimValues3(2)) fix(abs(xLimValues3(2)-fix(xLimValues3(2)))*10)];
			strFigName_add_w	= sprintf( '_processed_lcm_%d_%d_%d_%dppm', digits(1), digits(2), digits(3), digits(4) );
			figureName_w_fig	= [outFileName_w, strFigName_add_w, '.fig'];
			figureName_w_png	= [outFileName_w, strFigName_add_w, '.png'];
			saveFigure_s(h_w, outDirString, figureName_w_fig, 'fig', resolution);
			saveFigure_s(h_w, outDirString, figureName_w_png, 'png', resolution);
		end		% End of if with_water
		
		
		%% Display (water) reference signals and save correspoding figures, if available
		xLimValues4		= [3.3 5.9];
		if with_ref
			% Show and save (water) reference signals
			% Reference signal for ECC
			h_ref_ECC			= figure('visible','on');
			plot(out_ref_ECC.ppm,real(out_ref_ECC.specs),'linewidth',1.5);xlim(xLimValues4);
			axes_h_ref_ECC		= get(h_ref_ECC,'CurrentAxes');
			set(axes_h_ref_ECC,'FontSize',12, 'FontWeight','bold');
			set(axes_h_ref_ECC,'XDir','reverse');
			set(axes_h_ref_ECC,'XAxisLocation', 'origin');
			xlabel('ppm','FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle');
			ylabel('Amplitude(a.u.)','FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','baseline');
			box off;
			title('Preprocessed Reference Signal for ECC','FontSize',12);
			xf3		= xLimValues4(1) + (xLimValues4(2) - xLimValues4(1))/2;
			yf3		= 1.0*min(get(axes_h_ref_ECC, 'ylim'));
			set(get(axes_h_ref_ECC,'XLabel'),'Position', [xf3, yf3]);
			% Save figure of reference signal for ECC
			% Create filenames for saving of processed output
			%strFigName_addW	= '_w_processed_lcm_3_3_5_9ppm';
			digits = [fix(xLimValues2(1)) fix(abs(xLimValues4(1)-fix(xLimValues4(1)))*10) fix(xLimValues4(2)) fix(abs(xLimValues4(2)-fix(xLimValues4(2)))*10)];
			strFigName_add_ref_ECC	= sprintf( '_processed_lcm_%d_%d_%d_%dppm', digits(1), digits(2), digits(3), digits(4) );
			figureName_ref_ECC_fig	= [outFileName_ref_ECC, strFigName_add_ref_ECC, '.fig'];
			figureName_ref_ECC_png	= [outFileName_ref_ECC, strFigName_add_ref_ECC, '.png'];
			saveFigure_s(h_ref_ECC, outDirString, figureName_ref_ECC_fig, 'fig', resolution);
			saveFigure_s(h_ref_ECC, outDirString, figureName_ref_ECC_png, 'png', resolution);
			
			% Reference signal for metabolite quantification (Quant)
			h_ref_Quant			= figure('visible','on');
			plot(out_ref_Quant.ppm,real(out_ref_Quant.specs),'linewidth',1.5);xlim(xLimValues4);
			axes_h_ref_Quant	= get(h_ref_Quant,'CurrentAxes');
			set(axes_h_ref_Quant,'FontSize',12, 'FontWeight','bold');
			set(axes_h_ref_Quant,'XDir','reverse');
			set(axes_h_ref_Quant,'XAxisLocation', 'origin');
			xlabel('ppm','FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle');
			ylabel('Amplitude(a.u.)','FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','baseline');
			box off;
			title('Preprocessed Reference Signal for Quantification','FontSize',12);
			xf3		= xLimValues4(1) + (xLimValues4(2) - xLimValues4(1))/2;
			yf3		= 1.0*min(get(axes_h_ref_Quant, 'ylim'));
			set(get(axes_h_ref_Quant,'XLabel'),'Position', [xf3, yf3]);
			% Save figure of reference signal for metabolite quantification (Quant)
			% Create filenames for saving of processed output
			%strFigName_addW	= '_w_processed_lcm_3_3_5_9ppm';
			digits = [fix(xLimValues2(1)) fix(abs(xLimValues4(1)-fix(xLimValues4(1)))*10) fix(xLimValues4(2)) fix(abs(xLimValues4(2)-fix(xLimValues4(2)))*10)];
			strFigName_add_ref_Quant	= sprintf( '_processed_lcm_%d_%d_%d_%dppm', digits(1), digits(2), digits(3), digits(4) );
			figureName_ref_Quant_fig	= [outFileName_ref_Quant, strFigName_add_ref_Quant, '.fig'];
			figureName_ref_Quant_png	= [outFileName_ref_Quant, strFigName_add_ref_Quant, '.png'];
			saveFigure_s(h_ref_Quant, outDirString, figureName_ref_Quant_fig, 'fig', resolution);
			saveFigure_s(h_ref_Quant, outDirString, figureName_ref_Quant_png, 'png', resolution);	
		end		% End of if with_ref
		
		
		%% Write processed and unprocessed data (spectra water, and reference signals) to file, if user selected
		% Use .RAW format for LCModel
		% Ask for user input for writing results to file only, if minimum user input is 
		% NOT selected
		% If minimum user input selected, writing results to file is not optional
		wrt				= 'y';
		if strcmp(strMinUserIn,'n') || strcmp(strMinUserIn,'N')
			wrt			= input('Write results to file? ','s');
		end
		if wrt=='y' || wrt=='Y'
			disp(sMsg_newLines);
			fprintf('%s: Writing results to file ...\n\n', sFunctionName);
			RF = io_writelcm(out,[outDirString outFileName '_processed_lcm' '.RAW'],out.te);
			RF = io_writelcm(out_noproc,[outDirString outFileName '_unprocessed_lcm' '.RAW'],out_noproc.te);
			if with_water
				% Save water signal
				RF_w = io_writelcm(out_w,[outDirString  outFileName_w '_processed_lcm' '.RAW'],out_w.te);
				RF_w = io_writelcm(out_w_noproc,[outDirString outFileName_w '_unprocessed_lcm' '.RAW'],out_w_noproc.te);
			end		% End of if with_water
			if with_ref
				% Save (water) reference signals
				RF_ref_ECC   = io_writelcm(out_ref_ECC,[outDirString  outFileName_ref_ECC '_processed_lcm' '.RAW'],out_ref_ECC.te);
				RF_ref_ECC	 = io_writelcm(out_ref_ECC_noproc,[outDirString outFileName_ref_ECC '_unprocessed_lcm' '.RAW'],out_ref_ECC_noproc.te);
				RF_ref_Quant = io_writelcm(out_ref_Quant,[outDirString  outFileName_ref_Quant '_processed_lcm' '.RAW'],out_ref_Quant.te);
				RF_ref_Quant = io_writelcm(out_ref_Quant_noproc,[outDirString outFileName_ref_Quant '_unprocessed_lcm' '.RAW'],out_ref_Quant_noproc.te);				
			end		% End of if with_ref
		end
		
		
		%% Make figure(s) of spectrum for html report, if report switch is turned ON
		if reportSwitch == 1
			h_mrs_html			= figure('visible','off');
			plot(out.ppm,real(out.specs),'linewidth',2);xlim(xLimValues1);
			%set(gca,'FontSize',12, 'FontWeight','bold');
			set(gca,'FontSize',8);
			set(gca,'XDir','reverse');
			set(gca,'XAxisLocation', 'origin');
			%xlabel('ppm','FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle');
			%ylabel('Amplitude(a.u.)','FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','baseline');
			xlabel('ppm','FontSize',10, 'HorizontalAlignment','center', 'VerticalAlignment','middle');
			ylabel('Amplitude(a.u.)','FontSize',10, 'HorizontalAlignment','center', 'VerticalAlignment','baseline');
			%legend('Real(MRS Data)');
			%legend boxoff;
			box off;
			%title('Result: Preprocessed Spectrum','FontSize',12);
			title(strTitle_mrs,'FontSize',12);
			xf1		= xLimValues1(1) + (xLimValues1(2) - xLimValues1(1))/2;
			yf1		= 1.0*min(get(gca, 'ylim'));
			set(get(gca,'XLabel'),'Position', [xf1, yf1], 'VerticalAlignment', 'Top');
			set(h_mrs_html,'PaperUnits','centimeters');
			set(h_mrs_html,'PaperPosition',[0 0 20 10]);
			saveas(h_mrs_html, fullfile(outDirString, reportFigDirStr, 'finalSpecFig'), 'jpg');
			saveas(h_mrs_html, fullfile(outDirString, reportFigDirStr, 'finalSpecFig'), 'fig');
			
			% Copy and modify figure to show different spectral range and save modified figure
			% as well
			% For figure creation, i.e. if spectrum should be used as figure for display or in
			% paper, use only ppm range from 0.2 to 4.2 and 'whiten' y-axis
			% For regular use, leave y-axis as is to indicate signal strength
			ax_h_mrs_html	= gca;
			h_mrs_html2		= figure('visible','off');
			ax_h_mrs_html2	= copyobj(ax_h_mrs_html,h_mrs_html2);
			set(gca, 'XLim', xLimValues2, 'XTick',xTickValues2, 'XTickLabel',xTickValues2);
			xf2		= xLimValues2(1) + (xLimValues2(2) - xLimValues2(1))/2;
			yf2		= 1.0*min(get(gca, 'ylim'));
			set(get(gca,'XLabel'),'Position', [xf2, yf2], 'VerticalAlignment', 'Top');
			set(h_mrs_html2,'PaperUnits','centimeters');
			set(h_mrs_html2,'PaperPosition',[0 0 20 10]);
			saveas(h_mrs_html2, fullfile(outDirString, reportFigDirStr, 'finalSpecFig_narrow'), 'jpg');
			saveas(h_mrs_html2, fullfile(outDirString, reportFigDirStr, 'finalSpecFig_narrow'), 'fig');
			
			
			%% Write an html report, if report switch is turned ON
			% Adjust directory information and name report according to each case
			% Use relative paths from directory of html report to figure files, so that html
			% report can be opened on any computer/browser
			%fid2			= fopen(fullfile(outDirString, reportDirStr, 'report.html'),'w+');
			fid2			= fopen(fullfile(outDirString, reportDirStr, [outFileName '_report' '.html']),'w+');
			fprintf(fid2,'<!DOCTYPE html>');
			fprintf(fid2,'\n<html>');
			%logoPath	= which('FID-A_LOGO.jpg');
			%fprintf(fid2,'\n<img src= " %s " width="120" height="120"></body>',logoPath);
			logoFileName	= 'FID-A_LOGO.jpg';
			logoPath		= which(logoFileName);
			% Copy FID-A logo into output directory, if it does not already exist
			if exist(fullfile(reportDirStr, logoFileName), 'file') ~= 2
				%[status,msg] = copyfile(logoPath, fullfile(reportFigsDir, logoFileName));
				%[status,msg] = copyfile(logoPath, fullfile(outDirString, reportFigDirString, logoFileName));
				[status,msg] = copyfile(logoPath, fullfile(outDirString, logoFileName));
				if status ~= 1
					error('%s: Error copying logo file %s!\n\n%s', sFunctionName, logoFileName, msg);
				end
			end
			fprintf(fid2,'\n<img src= " %s " width="120" height="120"></body>', fullfile('./../', logoFileName));
			fprintf(fid2,'\n<h1>FID-A Processing Report</h1>');
			fprintf(fid2,'\n<h2>Processing pipeline applied to %s MRS data using %s </h2>', seqType, sFunctionName);
			fprintf(fid2,'\n<p>FILENAME: %s </p>', fullfile(dirString, filename));
			fprintf(fid2,'\n<p>Data type: %s </p>', dataType);
			fprintf(fid2,'\n<p>DATE: %s </p>',date);
			fprintf(fid2,'\n\n<p> </p>');
			
			% FLAG: Modified
			% Write info about coil combination into report depending on MRS data type
            %if ~(isIMA && isIMA_w)
			if ~(isIMA && isIMA_w) && (out_raw.dims.coils == 0)
			    fprintf(fid2,'\n\n<h2>Results of multi-coil combination:</h2>');
			    %fprintf(fid2,'\n<img src= " %s%scoilReconFig.jpg " width="800" height="400"></body>', outDirString, reportFigDirStr);
			    fprintf(fid2,'\n<img src= " %s " width="800" height="400"></body>', fullfile('./figs/', 'coilReconFig.jpg'));
			    fprintf(fid2,'\n\n<p> </p>');
			else
				fprintf(fid2,'\n\n<h2>Multi-coil combination was not performed, since MRS DICOM data already coil combined.</h2>');
				fprintf(fid2,'\n\n<p> </p>');
			%end		% End of  if ~(isIMA && isIMA_w)
			end		% End of if ~(isIMA && isIMA_w) && (out_raw.dims.coils == 0)
			
			% Write results from removing bad averages into html report, only if step was
			% actually performed, i.e. it was selected and dimension of averages existed
			if (rmbadav=='y' || rmbadav=='Y') && out_cc.dims.averages > 0
				fprintf(fid2,'\n\n<h2>Results of removal of bad averages:</h2>');
				fprintf(fid2,'\n<p>Original number of averages: \t%5.6f </p>', out_raw.sz(out_raw.dims.averages));
				fprintf(fid2,'\n<p>Number of bad Averages removed:  \t%5.6f </p>',nBadAvgTotal);
				fprintf(fid2,'\n<p>Number of remaining averages in processed dataset:  \t%5.6f </p>',out_rm.sz(out_rm.dims.averages));
				fprintf(fid2,'\n<p>Bad Averages Removal Threshold was:  \t%2.2f </p>', nSD);
				%fprintf(fid2,'\n<img src= " %s%srmBadAvg_prePostFig.jpg " width="800" height="600"><img src= " %s%srmBadAvg_scatterFig.jpg " width="800" height="400">', outDirString, reportFigDirStr, outDirString, reportFigDirStr);
				fprintf(fid2,'\n<img src= " %s " width="800" height="600"><img src= " %s " width="800" height="400">', fullfile('./figs/','rmBadAvg_prePostFig.jpg'), fullfile('./figs/','rmBadAvg_scatterFig.jpg'));
			else
				fprintf(fid2,'\n\n<h2>Removal of bad averages was NOT performed!</h2>');
				fprintf(fid2,'\n<p>rmbadav = %s </p>', rmbadav);
				fprintf(fid2,'\n<p>out_cc.dims.averages = %d </p>', out_cc.dims.averages);
			end
			fprintf(fid2,'\n\n<p> </p>');
			fprintf(fid2,'\n\n<p> </p>');
			
			% Write results from drift correction into html report, only if this step was
			% actually performed, i.e. it was selected and dimension of averages existed
			if (driftCorr=='y' || driftCorr=='Y') && out_rm.dims.averages > 0
				fprintf(fid2,'\n\n<h2>Results of spectral registration:</h2>');
				fprintf(fid2,'\n<p>Total frequency drift was: \t%5.6f </p>',max(totalFreqDrift));
				fprintf(fid2,'\n<p>Total phase drift was: \t%5.6f </p>',max(totalPhaseDrift));
				%fprintf(fid2,'\n<img src= " %s%salignAvgs_prePostFig.jpg " width="800" height="600">', outDirString, reportFigDirStr);
				fprintf(fid2,'\n<img src= " %s " width="800" height="600">', fullfile('./figs/','alignAvgs_prePostFig.jpg'));
				fprintf(fid2,'\n\n<p> </p>');
				%fprintf(fid2,'\n<img src= " %s%sfreqDriftFig.jpg " width="400" height="400"><img src="%s/%s/report/figs/phaseDriftFig.jpg " width="400" height="400">', outDirString, reportFigDirStr, outDirString, reportFigDirStr);
				fprintf(fid2,'\n<img src= " %s " width="400" height="400"><img src=" %s " width="400" height="400">', fullfile('./figs/','freqDriftFig.jpg'), fullfile('./figs/','phaseDriftFig.jpg'));
			else
				fprintf(fid2,'\n\n<h2>Drift correction for averages was NOT performed!</h2>');
				fprintf(fid2,'\n<p>driftCorr = %s </p>', driftCorr);
				fprintf(fid2,'\n<p>out_rm.dims.averages = %d </p>', out_rm.dims.averages);
			end
			fprintf(fid2,'\n\n<p> </p>');
			fprintf(fid2,'\n\n<h2>Final Result:</h2>');
			%fprintf(fid2,'\n<img src= " %s%sfinalSpecFig.jpg " width="800" height="400">', outDirString, reportFigDirStr);
			%fprintf(fid2,'\n<img src= " %s " width="800" height="400">', fullfile('./figs/','finalSpecFig.jpg'));
			fprintf(fid2,'\n<img src= " %s " width="800" height="400"><img src= " %s " width="800" height="400">', fullfile('./figs/','finalSpecFig.jpg'), fullfile('./figs/','finalSpecFig_narrow.jpg'));
			fclose(fid2);
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		end
		

	case 'SPECIAL'
		disp('NOT YET COMPLETED ...');
		
		% Create filenames for saving of processed output depending on sequence type
		% Adjust filenames for use of OVS with SPECIAL
		
	case 'MEGA-PRESS'
		disp('Coming soon ...');
		
		% Create filenames for saving of processed output depending on sequence type
		% Assuming that no bad average removal is applied for any of the water signals
		outFileName				= [nameSpec, sprintf('_%.1f', nSD)];
		outFileName_w			= [name_w, '_w'];
		outFileName_ref_ECC		= [nameSpec, '_ref_ECC', sprintf('_%.1f', nSD)];
		outFileName_ref_Quant	= [nameSpec, '_ref_Quant', sprintf('_%.1f', nSD)];
		%outFileName_w			= [name_w, '_w', sprintf('_%.1f', nSD)];
		%outFileName_ref_ECC		= [nameSpec, '_ref_ECC', sprintf('_%.1f', nSD)];
		%outFileName_ref_Quant	= [nameSpec, '_ref_Quant', sprintf('_%.1f', nSD)];
		
		% Figure display and saving
		% Can be realized by calling a separate routine?
		% Use code from run_specialproc_CBF.m
		
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType);
		
end		% End of switch seqType


