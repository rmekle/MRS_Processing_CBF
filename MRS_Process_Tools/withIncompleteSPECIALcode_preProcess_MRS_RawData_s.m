%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% preProcess_MRS_RawData_s.m
%
%% Function to measure the signal-to-noise ratio (SNR) and the linewidth (LW)/FWHM of a 
%  selected resonance/peak in single volume magnetic resonance spectroscopy (MRS) data
%  using functions from the MRS processing toolkit FID-A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USAGE
% [out,out_w,out_noproc,out_w_noproc] = preProcess_MRS_RawData_s(dirString,outDirString,filename,filenamew,seqType,dataType,strOVS,nSD,aaDomain,tmaxin,iterin,plotSwitch,strMinUserIn,reportSwitch);
% 
% DESCRIPTION:
% Processing script for Siemens MRS data in .dat format (twix raw data).  
% Includes combination of receiver channels, removal of bad averages, 
% frequency drift correction, and leftshifting.
% 
% INPUTS:
% dirString    = String variable for the name of the directory containing
%                   the water suppressed .dat file and the water unsuppressed .dat file.
% outDirString = String variable for the name of the output directory, i.e. the directory,
%					where all output files are saved to
% filename     = String variable for the name of the water suppressed .dat file,
%					e.g. 3T_20170510_PetraO_meas_MID00091_FID153257_rm_special_RF_ACC.dat
%					or specialDLPFC.dat
% filenamew    = String variable for the name of the water unsuppressed .dat file,
%					e.g. 3T_20170510_PetraO_meas_MID00091_FID153259_rm_special_RF_ACC.dat
%					or specialDLPFC_w.dat
% seqType	   = String specifying the MRS sequence type used for data acquisition, e.g.
%					'PRESS', 'STEAM', 'SPECIAL', 'sLASER', 'MEGA-PRESS'
% dataType	   = String describing the type of MRS data:
%					'water' = water signal, 'mrs' = MR spectrum without water signal, 
%					'mrs_w'= MR spectrum with respective water signal
% strOVS	   = (optional) String that specifies whether water acquired with OVS ('wOVS')
%					 or withoutOVS ('woutOVS') is used for processing. 
%					 Default is 'woutOVS', which means that OVS was not used. 
% nSD		   = (Optional) # of standard deviations for bad average removal. Default
%					value is 3.2.
% aaDomain     = (Optional) Perform the spectral registration (drift correction) using
%                   the full spectrum ('t'), or only a limited frequency
%                   range ('f').  Default is 'f'.
% tmaxin       = (Optional).  Duration (in sec.) of the time domain signal
%                   used in the spectral registration (drift correction).
%                   Default is 0.2 sec.
% iterin       = (Optional).  Maximum number of allowed iterations for the spectral
%                   registration to converge. Default is 20.
% plotSwitch   = (Optional)	Switch for displaying plots: 1 = ON, 0 = OFF. Default is 0
% strMinUserIn = (Optional) String that specifies whether user input/interaction should be
%					 minimized or not; 'y' or 'Y' lead to minimization, 'n' or 'N' do not
%					Default is 'y'.
% reportSwitch = (Optional) Switch for generating an html report with corresponding  
%					figures and a readme file: 1 = ON, 0 = OFF. Default is 1. 
% 
% OUTPUTS:
% out          = Fully processed, water suppressed output spectrum.
% out_w        = Fully processed, water unsuppressed output spectrum.
% out_noproc   = Water suppressed output spectrum without pre-
%                   processing (No bad-averages removal, no frequency drift
%                   correction).
% out_w_noproc = Water unsuppressed output spectrum without pre-
%                   processing.
%
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2021; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out,out_w,out_noproc,out_w_noproc] = preProcess_MRS_RawData_s(dirString,outDirString,filename,filenamew,seqType,dataType,strOVS,nSD,aaDomain,tmaxin,iterin,plotSwitch,strMinUserIn,reportSwitch)

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'preProcess_MRS_RawData_s.m';
sMsg_newLines		= sprintf('\n\n');
sMsg_newLine		= sprintf('\n');
disp(sMsg_newLines);


%% Check on # of input arguments and assign default values to variables, if required
maxNargin	= 14;
if nargin<maxNargin
	reportSwitch = 1;
	if nargin<(maxNargin-1)
		strMinUserIn = 'y';
		if nargin<(maxNargin-2)
			plotSwitch = 0;
			if nargin<(maxNargin-3)
				iterin = 20;
				if nargin<(maxNargin-4)
					tmaxin = 0.2;
					if nargin<(maxNargin-5)
						aaDomain = 'f';
						if nargin<(maxNargin-6)
							nSD = 3.2;
							if nargin<(maxNargin-7)
								strOVS = 'woutOVS';
							end
						end
					end
				end
			end
		end
	end
end

% Use filenames for the MRS datasets provided as input arguments
% NOTE: It is assumed that water suppressed and water unsuppressed data are in same
% directory!

% Ensure that input and output directory string end with file separator, i.e. '/' or '\'
% (Windows handles the Unix '/' just fine)
if( ~strcmp( dirString(end), filesep ) )
	dirString	= [dirString, filesep];
end

if( ~strcmp( outDirString(end), filesep ) )
	outDirString	= [outDirString, filesep];
end

% Obtain different parts of filenames
[sPathStrSpec,nameSpec,extSpec] 	= fileparts(filename);
[sPathStrW,nameW, extW] 			= fileparts(filenamew);

% Create filenames for saving of processed output depending on sequence type
if( strcmp(seqType, 'MEGA-PRESS') )
	% Assuming that MEGA-PRESS editing is usually performed without OVS
	outFileName		= [nameSpec, sprintf('_%.1f', nSD)];
	outFileNameW	= [nameW, '_w', sprintf('_%.1f', nSD)];
else
	% Water signals and/or MR spectra were acquired with or without OVS
	outFileName		= [nameSpec, '_', strOVS, sprintf('_%.1f', nSD)];
	outFileNameW	= [nameW, '_w', '_', strOVS, sprintf('_%.1f', nSD)];
end


%% Set parameters for figure display
% h		= figure('position', [left bottom width height]);
fig_left	= 20;
fig_bottom	= 50;
fig_width	= 560;
fig_height	= 420;
fig_dist_b	= 510;
fig_dist_l	= 570;


%% Depending on type of data provided load all corresponding data files
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
		if isempty(filenamew)
			error('%s: Error: Filename for unsuppresed water signal %s is empty!\n\n', sFunctionName, filenamew);
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
		if isempty(filenamew)
			error('%s: Error: Filename for unsuppresed water signal %s is empty!\n\n', sFunctionName, filenamew);
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
% 	if isempty(filenamew)
% 		error('%s: Error: Filename for unsuppresed water signal %s is empty!\n\n', sFunctionName, filenamew);
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
%out_raw				= io_loadspec_twix([dirString filename]);
[out_raw, out_ref_raw]	= io_loadspec_twix_s([dirString filename]);

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
	disp('***WITH ADDITIONAL WATER UNSUPPRESSED DATA***');
	out_w_raw		= io_loadspec_twix([dirString filenamew]);
	
	% Convert single precision data (default format used my mapVBVD.m for imaging data) 
	% into double precision for processing, if not empty
	if ~isempty(out_w_raw)
		out_w_raw.fids		= double(out_w_raw.fids);
		out_w_raw.specs		= double(out_w_raw.specs);
	end
else
	disp('***WITHOUT ADDITIONAL WATER UNSUPPRESSED DATA***');    
    %out_w			= struct([]);
    %out_w_noproc	= struct([)];
end

% Make a new directory for the output report and figures each, if not already existent,
% and if desrired
if reportSwitch == 1
	%mkdir([outDirString nameSpec '/report']);
	%mkdir([outDirString nameSpec '/report/figs']);
	reportDirStr		= [nameSpec '_report/'];
	reportFigDirStr		= [nameSpec '_report/figs/'];
	if ~exist( [outDirString reportDirStr], 'dir' )
		mkdir([outDirString reportDirStr]);
	end
	if ~exist( [outDirString reportFigDirStr], 'dir' )
		mkdir([outDirString reportFigDirStr]);
	end
end
disp(sMsg_newLines);


%% Pre-process all data according to sequence type and types of data provided
switch seqType
	case {'PRESS', 'STEAM', 'sLASER'}
		%% Is this Dinesh K. Deelchand's single voxel (sLaser) sequence from CMRR, U Minnesota
		isSVSdkd_seq = contains(out_raw.seq,'svs_slaser_dkd');
		
		
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
				
			else
				error('%s: Reference scans option for sequence %s not yet implemented!', sFunctionName, out.seq);
			end		% End of if isSVSdkd_seq 
			
		end		% End of if with_ref

		
		%% Combine signals from different coil elements
		% First step should be to combine coil channels. For this find the coil phases 
		% from the water unsuppressed data, if available; otherwise from the MR spectra
		% Arguments referring to the nPos_ccth point of the FID and weighting of channels 
		% based on maximum signal ('w') are ignored, when coilcombos are provided as input
		nPos_cc		= 1;
		nPos_cc_w	= 1;
		if with_water
			% Obtain coil phases and amplitudes from unsuppressed water signal
			%coilcombos	= op_getcoilcombos(out_raw_w,1);
			coilcombos	= op_getcoilcombos(out_raw_w,nPos_cc_w,'w');
			[outw_cc,fidw_pre,specw_pre,phw,sigw]	= op_addrcvrs(out_w_raw,nPos_cc_w,'w',coilcombos);
		else
			if with_ref
				% Obtain coil phases and amplitudes from reference scans (water signal)
				if isSVSdkd_seq		% dkd sequence from CMRR
					% Obtain coil phases and amplitudes from reference scans for ECC
					coilcombos	= op_getcoilcombos(out_ref_ECC_raw,nPos_cc_w,'w');
				else
					coilcombos	= op_getcoilcombos(out_ref_raw,nPos_cc_w,'w');
				end
			else
				%coilcombos	= op_getcoilcombos(op_averaging(out_raw),1);
				coilcombos	= op_getcoilcombos(op_averaging(out_raw),nPos_cc,'w');
			end
		end
		% Combine coil channels before and after signal averaging for comparison and
		% plotting
		[out_cc,fid_pre,spec_pre,ph,sig]	= op_addrcvrs(out_raw,nPos_cc,'w',coilcombos);
		[out_av_cc,fid_av_pre,spec_av_pre]	= op_addrcvrs(op_averaging(out_raw),nPos_cc,'w',coilcombos);
		out_raw_av							= op_averaging(out_raw);
		
		% Generate unprocessed spectrum or spectra, respectively
		out_noproc			= op_averaging(out_cc);
		if with_water
			outw_noproc		= op_averaging(outw_cc);
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
			plot(out_w_raw.ppm,out_w_raw.specs(:,:,1,1));xlim([4 5]);
			xlabel('Frequency (ppm)');
			ylabel('Amplitude (a.u.)');
			title('Multi-channel (water unsupp.) data before phase correction');
			subplot(2,1,2);
			plot(out_w_raw.ppm,spec_w_pre(:,:,1,1));xlim([4 5]);
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
		
		
		%% Remove bad averages from data
		%%%%%%%% OPTIONAL REMOVAL OF BAD AVERAGES FROM DATASET %%%%%%%%%%%%%%%%%%%%
		close all;
		out_cc2			= out_cc;
		nBadAvgTotal	= 0;
		nbadAverages	= 1;
		rmbadav			= 'y';
		
		% Do not remove bad averages, if either not selected or if dimension of 
		% averages does not exist (index for dimension of averages = 0), 
		% e.g. when data is already averaged
		%if rmbadav=='n' || rmbadav=='N'
		if rmbadav=='n' || rmbadav=='N' || out_cc.dims.averages == 0
			out_rm		= out_cc;
			%nSD			= 'N/A';
		else
			sat		='n'
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
				
			end
			
			% Write a readme file to record the number of dropped avgs to output directory
			% for report instead of input data directory, if report switch is turned ON
			if reportSwitch == 1
				%fid1=fopen([dirString '/readme.txt'],'w+');
				%fid1		= fopen([outDirString outFileName '_readme.txt'],'w+');
				%fid1		= fopen([outDirString reportDirStr 'readme.txt'],'w+');
				fid1		= fopen([outDirString reportDirStr outFileName '_readme.txt'],'w+');
				fprintf(fid1,'Original number of averages: \t%5.6f',out_raw.sz(out_raw.dims.averages));
				disp(['Original number of averages:  ' num2str(out_raw.sz(out_raw.dims.averages))]);
				fprintf(fid1,'\nNumber of bad Averages removed:  \t%5.6f',nBadAvgTotal);
				disp(['Number of bad averages removed:  ' num2str(nBadAvgTotal)]);
				fprintf(fid1,'\nNumber of remaining averages in processed dataset:  \t%5.6f',out_rm.sz(out_rm.dims.averages));
				disp(['Number of remaining averages in processed dataset:  ' num2str(out_rm.sz(out_rm.dims.averages))]);
				fclose(fid1);
			end
		end
				
		%%%%%%%%%%%%%%%%%%%% END OF BAD AVERAGES REMOVAL %%%%%%%%%%%%%%%%%%%%
		
		
		%% NOW ALIGN AVERAGES:  A.K.A. Frequency Drift Correction
		% If minimum user input selected, frequency drift correction is not optional
		% Ask for user input for frequency correction only, if minimum user input is NOT
		% selected
		% Set parameters for drift correction depending on type of data, i.e. whether MRS
		% data is spectrum or water signal
		% NOTE: Check whether aligning of averages in frequency domain works, if the MR
		% spectrum is water signal itself; if not, simply align averages in time domain 
		if( strcmp(dataType, 'mrs_w') || strcmp(dataType, 'mrs') )
			% MR spectrum is provided together with or without unsuppressed water signal 
			ppmmin_fix		= 1.6;
			ppmmaxarray_fix	= [3.5; 4.0; 5.5];
			iamax			= 6;
		else
			% MR spectrum is water signal itself
			ppmmin_fix		= 4.2;
			ppmmaxarray_fix	= [5.5 5.5 5.2];
			iamax			= 6;
		end
		
		% Do not perform drift correction, if either not selected or if dimension of 
		% averages does not exist (index for dimension of averages = 0), 
		% e.g. when data is already averaged
		driftCorr		= 'y';
		%if driftCorr=='n' || driftCorr=='N'
		if driftCorr=='n' || driftCorr=='N' || out_rm.dims.averages == 0
			out_av		= op_averaging(out_rm);
			if with_water
				outw_av		= op_averaging(outw_cc);
			end
			fs			= 0;
			phs			= 0;
		else
			if with_water
				%outw_aa		= op_alignAverages(outw_cc,tmaxin,'n');
				outw_aa		= op_alignAverages(outw_cc,0.2,'n');
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
						case 'f'
							[out_aa,fs,phs]		= op_alignAverages_fd(out_rm2,ppmmin,ppmmax,tmax,'y');
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
				end
				
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
				title('Estimated Freqeuncy Drift','FontSize',12);
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
			end
			% Now average the aligned averages
			out_av		= op_averaging(out_aa);
			if with_water
				outw_av		= op_averaging(outw_aa);
			end
		end
		
		
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
			if outw_av.dims.subSpecs ~= 0
				disp(sMsg_newLines);
				warning('%s: Subspectra still present in unsuppressed water signal! # of subsprectra = outw_av.sz(outw_av.dims.subSpecs) = %d', sFunctionName, outw_av.sz(outw_av.dims.subSpecs));
				disp(sMsg_newLine);
				disp('Averaging subSpecs of water signal ...');
				disp(sMsg_newLines);
				outw_av_tmp	= outw_av;
				outw_av				= op_average_subSpecs_s(outw_av_tmp);
				clear outw_av_tmp;
				
				% Since unprocessed water signals will also be written to file in LCModel 
				% format, the unprocessed subSpectra have to be averaged as well to avoid 
				% any errors in output routine io_writelcm.m due to isISIS flag being 
				% still set to 1
				if out_w_noproc.dims.subSpecs ~= 0
					out_w_noproc_tmp		= out_w_noproc;
					out_w_noproc			= op_average_subSpecs_s(out_w_noproc_tmp);
					clear out_w_noproc_tmp;
				end		% End of if out_noproc.dims.subSpecs ~= 0
			end		% End of if outw_av.dims.subSpecs ~= 0
		end		% End of if with_water
		
		
		%% Perform phase and frequency corrections
		% Left shift the MRS signals to eliminate first order phase, perform zero-order 
		% phase correction and shift data in frequency to obtain reference peaks at known
		% positions
		% Now left shift
		out_ls		= op_leftshift(out_av,out_av.pointsToLeftshift);
		if with_water
			outw_ls		= op_leftshift(outw_av,outw_av.pointsToLeftshift);
		end
		
		% Perform automatic zero-order phase correction
		% If data is MR spectrum, use creatine peak at 3.027 ppm:
		if( strcmp(dataType, 'mrs_w') || strcmp(dataType, 'mrs') )
			out_ls_zp		= op_zeropad(out_ls,16);
			[out_ph,ph0]	= op_autophase(out_ls,2.9,3.1);
			out_ls_zp		= op_addphase(out_ls_zp,ph0);
			% And now for water unsuppressed data (use water peak):
			if with_water
				outw_ls_zp	= op_zeropad(outw_ls,16);
				indexw		= find(abs(outw_ls_zp.specs) == max(abs(outw_ls_zp.specs(outw_ls_zp.ppm>4 & outw_ls_zp.ppm<5.5))));
				ph0w		= -phase(outw_ls_zp.specs(indexw))*180/pi;
				outw_ph		= op_addphase(outw_ls,ph0w);
				outw_ls_zp	= op_addphase(outw_ls_zp,ph0w);
			end
		else
			% MR spectrum is water signal itself, so use water peak:
			out_ls_zp	= op_zeropad(out_ls,16);
			indexw		= find(abs(out_ls_zp.specs) == max(abs(out_ls_zp.specs(out_ls_zp.ppm>4 & out_ls_zp.ppm<5.5))));
			ph0			= -phase(out_ls_zp.specs(indexw))*180/pi;
			out_ph		= op_addphase(out_ls,ph0);
			out_ls_zp	= op_addphase(out_ls_zp,ph0);
		end
		
		% Perform same phase corection on unprocessed data
		out_noproc		= op_addphase(op_leftshift(out_noproc,out_noproc.pointsToLeftshift),ph0);
		if with_water
			outw_noproc		= op_addphase(op_leftshift(outw_noproc,outw_noproc.pointsToLeftshift),ph0w);
		end

		
		% Frequency shift all spectra
		% If data is MR spectrum, so that creatine appears at 3.027 ppm
		if( strcmp(dataType, 'mrs_w') || strcmp(dataType, 'mrs') )
			[~,frqShift]	= op_ppmref(out_ls_zp,2.9,3.1,3.027);
			out				= op_freqshift(out_ph,frqShift);
			out_noproc		= op_freqshift(out_noproc,frqShift);
			% And now for water unsuppressed data (user water peak and set it to 4.65 ppm)
			if with_water
				[~,frqShiftw]	= op_ppmref(outw_ls_zp,4,5.5,4.65);
				outw			= op_freqshift(outw_ph,frqShiftw);
				outw_noproc		= op_freqshift(outw_noproc,frqShiftw);
			end
		else
			% MR spectrum is water signal itself, so use water peak and set it to 4.65 ppm
			[~,frqShift]	= op_ppmref(out_ls_zp,4,5.5,4.65);
			out				= op_freqshift(out_ph,frqShift);
			out_noproc		= op_freqshift(out_noproc,frqShift);
		end
		close all;		
		
		
		
		%% Make figure(s) to show the final spectra and save results
		
			
		%% Display spectra and save corresponding figures
		%h=figure('visible','off');
		%plot(out.ppm,real(out.specs),'linewidth',2);xlim([0.2 5.2]);
		%xlabel('Frequency (ppm)','FontSize',10);
		%ylabel('Amplitude(a.u.)','FontSize',10);
		
		% Set figure parameters
		% Set plotting resolution and properties of axes depending on data type
		resolution		= 600;
		% If data is MR spectrum mor MR spectrum and a water signal
		if( strcmp(dataType, 'mrs_w') || strcmp(dataType, 'mrs') )
			xLimValues1		= [0.0 5.5];
			xLimValues2		= [0.2 4.2];
			xTickValues2	= [0.5:0.5:4.0];
		else
			% MR spectrum is water signal itself
			xLimValues1		= [3.3 5.9];
			xLimValues2		= [4.2 5.1];
			xTickValues2	= [4.2:0.2:5.0];
		end
		h_mrs			= figure('visible','on');
		plot(out.ppm,real(out.specs),'linewidth',1.5);xlim(xLimValues1);
		set(gca,'FontSize',12, 'FontWeight','bold');
		set(gca,'XDir','reverse');
		set(gca,'XAxisLocation', 'origin');
		xlabel('ppm','FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle');
		ylabel('Amplitude(a.u.)','FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','baseline');
		box off;
		title('Result: Preprocessed Spectrum','FontSize',12);
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
		
		% Modify figure to show different spectral range and save modified figure es well
		% For figure creation, i.e. if spectrum should be used as figure for display or in 
		% paper, use only ppm range from 0.2 to 4.2 and 'whiten' y-axis
		% For regular use, leave y-axis as is to indicate signal strength
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
		saveFigure_s(h_mrs, outDirString, figureName_fig, 'fig', resolution);
		saveFigure_s(h_mrs, outDirString, figureName_png, 'png', resolution);
		
		
		%% Display water signals and save correspoding figures, if available
		xLimValues3		= [3.3 5.9];
		if with_water
			% Show and save water signal
			h_w				= figure('visible','on');
			plot(outw_subSpec2.ppm,real(outw.specs),'linewidth',1.5);xlim(xLimValues3);
			axes_h_w		= get(h_w,'CurrentAxes');
			set(axes_h_w,'FontSize',12, 'FontWeight','bold');
			set(axes_h_w,'XDir','reverse');
			set(axes_h_w,'XAxisLocation', 'origin');
			xlabel('ppm','FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','middle');
			ylabel('Amplitude(a.u.)','FontSize',16, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','baseline');
			box off;
			title('Result: Preprocessed Water Signal','FontSize',12);
			xf3		= xLimValues3(1) + (xLimValues3(2) - xLimValues3(1))/2;
			yf3		= 1.0*min(get(axes_h_w, 'ylim'));
			set(get(axes_h_w,'XLabel'),'Position', [xf3, yf3]);
			% Save figure of water signal
			% Create filenames for saving of processed output
			%strFigName_addW	= '_w_processed_lcm_3_3_5_9ppm';
			digits = [fix(xLimValues2(1)) fix(abs(xLimValues3(1)-fix(xLimValues3(1)))*10) fix(xLimValues3(2)) fix(abs(xLimValues3(2)-fix(xLimValues3(2)))*10)];
			strFigName_addW	= sprintf( '_processed_lcm_%d_%d_%d_%dppm', digits(1), digits(2), digits(3), digits(4) );
			figureNameW_fig	= [outFileNameW, strFigName_addW, '.fig'];
			figureNameW_png	= [outFileNameW, strFigName_addW, '.png'];
			saveFigure_s(h_w, outDirString, figureNameW_fig, 'fig', resolution);
			saveFigure_s(h_w, outDirString, figureNameW_png, 'png', resolution);

		end		% End of if with_water
		
		
		%% If minimum user input selected, writing results to file is not optional
		% Write processed and unprocessed data (spectra and water) to file, if user selected
		% Use .RAW format for LCModel
		% Ask for user input for writing results to file only, if minimum user input is NOT
		% selected
		wrt				= 'y';
		if strcmp(strMinUserIn,'n') || strcmp(strMinUserIn,'N')
			wrt			= input('Write results to file? ','s');
		end
		if wrt=='y' || wrt=='Y'
			disp(sMsg_newLines);
			disp('Writing results to file ...');
			RF=io_writelcm(out,[outDirString outFileName '_processed_lcm' '.RAW'],out.te);
			RF=io_writelcm(out_noproc,[outDirString outFileName '_unprocessed_lcm' '.RAW'],out_noproc.te);
			if with_water
				RF=io_writelcm(out_w,[outDirString  outFileNameW '_processed_lcm' '.RAW'],out_w.te);
				RF=io_writelcm(out_w_noproc,[outDirString outFileNameW '_unprocessed_lcm' '.RAW'],out_w_noproc.te);
			end
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
			title('Result: Preprocessed Spectrum','FontSize',12);
			xf1		= xLimValues1(1) + (xLimValues1(2) - xLimValues1(1))/2;
			yf1		= 1.0*min(get(gca, 'ylim'));
			set(get(gca,'XLabel'),'Position', [xf1, yf1], 'VerticalAlignment', 'Top');
			set(h_mrs_html,'PaperUnits','centimeters');
			set(h_mrs_html,'PaperPosition',[0 0 20 10]);
			saveas(h_mrs_html, fullfile(outDirString, reportFigDirStr, 'finalSpecFig'), 'jpg');
			saveas(h_mrs_html, fullfile(outDirString, reportFigDirStr, 'finalSpecFig'), 'fig');
			
			
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
			fprintf(fid2,'\n<p>DATE: %s </p>',date);
			fprintf(fid2,'\n\n<p> </p>');
			fprintf(fid2,'\n\n<h2>Results of multi-coil combination:</h2>');
			%fprintf(fid2,'\n<img src= " %s%scoilReconFig.jpg " width="800" height="400"></body>', outDirString, reportFigDirStr);
			fprintf(fid2,'\n<img src= " %s " width="800" height="400"></body>', fullfile('./figs/', 'coilReconFig.jpg'));
			fprintf(fid2,'\n\n<p> </p>');
			
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
			fprintf(fid2,'\n<img src= " %s " width="800" height="400">', fullfile('./figs/','finalSpecFig.jpg'));
			fclose(fid2);
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		end
		

	case 'SPECIAL'
		disp('NOT YET COMPLETED ...');
		% Create filenames for saving of processed output including OVS options
		outFileName		= [nameSpec, '_', strOVS, sprintf('_%.1f', nSD)];
		outFileNameW	= [nameW, '_w', '_', strOVS, sprintf('_%.1f', nSD)];
		
		
		%% now combine the coil channels:
		
		% Find the coil phases either from additional unsuppressed watrer signal or from MR
		% spectrum
		% CBF: Note that here the second point of the time domain signal is used
		% 'diff' option for combining spectra means addition (see manual)!
		if with_water
			coilcombos=op_getcoilcombos(op_combinesubspecs(out_w_raw,'diff'),2);
		else
			coilcombos=op_getcoilcombos_specReg(op_combinesubspecs(op_averaging(out_raw),'diff'),0,0.01,2);
		end
		
		% CBF: Arguments referring to the 2nd point of the FID and weighting of channels based
		% on maximum signal ('w') are ignored, when coilcombos are provided as input
		[out_cc,fid_pre,spec_pre,ph,sig]=op_addrcvrs(out_raw,2,'w',coilcombos);
		if with_water
			[out_w_cc,fid_w_pre,spec_w_pre,ph_w,sig_w]=op_addrcvrs(out_w_raw,2,'w',coilcombos);
		end
		
		%make the un-processed spectra, which may be optionally output from this function:
		% CBF: Non-processed SPECIAL data - i.e. ISIS combination to obtain fully localized
		% signals is carried out - are generated for comparison;
		% only (manual) phase correction is applied below, and the non-processed data are saved as
		% well
		out_noproc=op_combinesubspecs(op_averaging(out_cc),'diff');
		if with_water
			out_w_noproc=op_combinesubspecs(op_averaging(out_w_cc),'diff');
		end
		
		%plot the data before and after coil phasing:
		figure('position',[0 50 560 420]);
		subplot(2,1,1);
		plot(out_raw.ppm,out_raw.specs(:,:,1,1));xlim([-1 7]);
		xlabel('Frequency (ppm)');
		ylabel('Amplitude (a.u.)');
		title('Multi-channel (water supp.) data before phase correction');
		subplot(2,1,2);
		plot(out_raw.ppm,spec_pre(:,:,1,1));xlim([-1 7]);
		xlabel('Frequency (ppm)');
		ylabel('Amplitude (a.u.)');
		title('Multi-channel (water supp.) data after phase correction');
		
		if with_water
			figure('position',[0 550 560 420]);
			subplot(2,1,1);
			plot(out_w_raw.ppm,out_w_raw.specs(:,:,1,1));xlim([4 5]);
			xlabel('Frequency (ppm)');
			ylabel('Amplitude (a.u.)');
			title('Multi-channel (water unsupp.) data before phase correction');
			subplot(2,1,2);
			plot(out_w_raw.ppm,spec_w_pre(:,:,1,1));xlim([4 5]);
			xlabel('Frequency (ppm)');
			ylabel('Amplitude (a.u.)');
			title('Multi-channel (water unsupp.) data after phase correction');
		end
		%disp('Press any key to continue...');
		%pause;
		close all;
		
		
		%% CBF: If minimum user input selected, alignment of subspectra is not optional
		
		%%%%%%%%%%%%%%%%%%%%%OPTIONAL ALIGNMENT OF SUBSPECTRA%%%%%%%%%%%%%%%%
		fs_ai			= [];
		phs_ai			= [];
		% CBF: Ask for user input for alignment only, if minimum user input is NOT selected
		alignISIS		= 'y';
		if strcmp(strMinUserIn,'n') || strcmp(strMinUserIn,'N')
			alignISIS		= input('would you like to align subspectra?  ','s');
		end
		if strcmp(alignISIS,'y') || strcmp(alignISIS,'Y')
			%What we're actually doing is aligning the averages, then aligning the
			%subspectra, then aligning the averages again, and then aligning the
			%subspectra again.
			% CBF: Spectral registration in the time domain using the first 0.4 s of each signal
			% and the average of the averages as reference to align with ('y');
			% additional outputs are frequency shifts (in Hz) and phase shifts (in degrees)
			[out_ai,fs_temp,phs_temp]=op_alignAverages(out_cc,0.4,'y');
			fs_ai=fs_temp;
			phs_ai=phs_temp;
			% CBF: Apply spectral registration to align ISIS subspectra prior to subtraction
			% using the first 0.4 s of each signal and add shifts from this alignment to
			% second column of shifts from previous alignment(s)
			[out_ai,fs_temp,phs_temp]=op_alignISIS(out_ai,0.4);
			fs_ai(:,2)=fs_ai(:,2)+fs_temp;
			phs_ai(:,2)=phs_ai(:,2)+phs_temp;
			% CBF: Spectral registration in the time domain using the first 0.4 s of each signal
			% and the average of the averages as reference to align with ('y');
			% add shifts from this alignment to shifts from previous alignment(s)
			[out_ai,fs_temp,phs_temp]=op_alignAverages(out_ai,0.4,'y');
			fs_ai=fs_ai+fs_temp;
			phs_ai=phs_ai+phs_temp;
			% CBF: Apply spectral registration to align ISIS subspectra prior to subtraction
			% using the first 0.4 s of each signal and add shifts from this alignment to
			% second column of shifts from previous alignment(s)
			[out_ai,fs_temp,phs_temp]=op_alignISIS(out_ai,0.4);
			fs_ai(:,2)=fs_ai(:,2)+fs_temp;
			phs_ai(:,2)=phs_ai(:,2)+phs_temp;
			
			%for fs_ai and phs_ai, take the average across both subspecs:
			fs_ai=mean(fs_ai,2);
			phs_ai=mean(phs_ai,2);
			
			if with_water
				%Now repeat above for water unsuppressed data:
				[out_w_ai,fs_w_temp,phs_w_temp]=op_alignAverages(out_w_cc,0.4,'y');
				fs_w_ai=fs_w_temp;
				phs_w_ai=phs_w_temp;
				[out_w_ai,fs_w_temp,phs_w_temp]=op_alignISIS(out_w_ai,0.4);
				fs_w_ai(:,2)=fs_w_ai(:,2)+fs_w_temp;
				phs_w_ai(:,2)=phs_w_ai(:,2)+phs_w_temp;
				[out_w_ai,fs_w_temp,phs_w_temp]=op_alignAverages(out_w_ai,0.4,'y');
				fs_w_ai=fs_w_ai+fs_w_temp;
				phs_w_ai=phs_w_ai+phs_w_temp;
				[out_w_ai,fs_w_temp,phs_w_temp]=op_alignISIS(out_w_ai,0.4);
				fs_w_ai(:,2)=fs_w_ai(:,2)+fs_w_temp;
				phs_w_ai(:,2)=phs_w_ai(:,2)+phs_w_temp;
				
				%for fs_w_ai and phs_w_ai, take the average across both subspecs:
				fs_w_ai=mean(fs_w_ai,2);
				phs_w_ai=mean(phs_w_ai,2);
			end
			
			%Now check the plots to make sure that they look okay:
			%First make combined subspecs plots:
			out_cc_temp=op_combinesubspecs(out_cc,'diff');
			out_ai_cc_temp=op_combinesubspecs(out_ai,'diff');
			if with_water
				out_w_cc_temp=op_combinesubspecs(out_w_cc,'diff');
				out_w_ai_cc_temp=op_combinesubspecs(out_w_ai,'diff');
			end
			
			%Now plot them
			close all
			figure('position',[0 0 560 420]);
			subplot(1,2,1);
			plot(out_cc_temp.ppm,out_cc_temp.specs);
			xlim([0 5]);
			xlabel('Frequency (ppm)');
			ylabel('Amplitude (a.u.)');
			title('Subspecs not aligned: (all averages)');
			subplot(1,2,2);
			plot(out_ai_cc_temp.ppm,out_ai_cc_temp.specs);
			xlim([0 5]);
			xlabel('Frequency (ppm)');
			ylabel('Amplitude (a.u.)');
			title('Subspecs aligned: (all averages)');
			
			if with_water
				figure('position',[0 550 560 420]);
				subplot(1,2,1);
				plot(out_w_cc_temp.ppm,out_w_cc_temp.specs);
				xlim([3.7 5.7]);
				xlabel('Frequency (ppm)');
				ylabel('Amplitude (a.u.)');
				title('Subspecs not aligned: (all averages)');
				subplot(1,2,2);
				plot(out_w_ai_cc_temp.ppm,out_w_ai_cc_temp.specs);
				xlim([3.7 5.7]);
				xlabel('Frequency (ppm)');
				ylabel('Amplitude (a.u.)');
				title('Subspecs aligned: (all averages)');
			end
			
			figure('position',[570 50 560 420]);
			subplot(2,1,1);
			plot([1:length(fs_ai)],fs_ai);
			xlabel('Scan Number');
			ylabel('Frequency Drift (Hz)');
			title('Estimated Frequency Drift');
			subplot(2,1,2);
			plot([1:length(phs_ai)],phs_ai);
			xlabel('Scan Number');
			ylabel('Phase Drift (Degrees)');
			title('Estimated Phase Drift');
			
			sat		= 'y';
			%sat=input('are you satisfied with alignment of subspecs? ','s');
			
			if strcmp(fsPolysat,'n') || strcmp(sat,'N')
				out_ai=out_cc;
				if with_water
					out_w_ai=out_w_cc;
				end
			end
			
		else
			out_ai=out_cc;
			if with_water
				out_w_ai=out_w_cc;
			end
			fs_ai=zeros(size(out_cc.fids,out_cc.dims.averages),1);
			phs_ai=zeros(size(out_cc.fids,out_cc.dims.averages),1);
		end
		
		% Now combine the subspecs to obtain fully localized SPECIAL spectra
		out_cs=op_combinesubspecs(out_ai,'diff');
		if with_water
			out_w_cs=op_combinesubspecs(out_w_ai,'diff');
		end
		
		%%%%%%%%%%%%%%%%%%%%%END OF ALIGNMENT OF SUBSPECTRA%%%%%%%%%%%%%%%%%%
		
		
		%% CBF: If minimum user input selected, removal of bad averages is not optional
		
		%%%%%%%%%%%%%%%%%%%%%OPTIONAL REMOVAL OF BAD AVERAGES%%%%%%%%%%%%%%%%
		close all
		figure('position',[0 50 560 420]);
		plot(out_cs.ppm,out_cs.specs);
		xlabel('Frequency (ppm)');
		ylabel('Amplitude (a.u.)');
		title('Water suppressed spectra (all averages)');
		
		out_cs2=out_cs;
		nBadAvgTotal=0;
		nbadAverages=1;
		allAveragesLeft=[1:out_cs.sz(out_cs.dims.averages)]';
		allBadAverages=[];
		% CBF: Ask for user input for removal of bad averages only, if minimum user input is NOT
		% selected
		rmbadav		= 'y';
		if strcmp(strMinUserIn,'n') || strcmp(strMinUserIn,'N')
			rmbadav=input('would you like to remove bad averages?  ','s');
		end
		
		close all;
		if rmbadav=='n' || rmbadav=='N'
			out_rm=out_cs;
		else
			% CBF: Remove bad averages; choose # of standard deviations by which averages have to
			% differ to be classified as 'bad'
			% (in _auto version of this script, i.e. with essentially no user input, 4 is used)
			sat='n'
			while sat=='n'||sat=='N'
				% CBF: Ask for user input for # of standard deviations only, if minimum user
				% input is NOT selected
				% 'nSD' was initialized at beginning of routine
				if strcmp(strMinUserIn,'n') || strcmp(strMinUserIn,'N')
					nSD		= input('input number of standard deviations.  ');
				end
				iter=1;
				nbadAverages=1;
				nBadAvgTotal=0;
				allAveragesLeft=[1:out_cs.sz(out_cs.dims.averages)]';
				allBadAverages=[];
				out_cs2=out_cs;
				while nbadAverages>0
					[out_rm,metric{iter},badAverages]=op_rmbadaverages(out_cs2,nSD,'t');
					badAverages;
					allBadAverages=[allBadAverages; allAveragesLeft(badAverages)];
					badavMask_temp=zeros(length(allAveragesLeft),1);
					badavMask_temp(badAverages)=1;
					allAveragesLeft=allAveragesLeft(~badavMask_temp);
					nbadAverages=length(badAverages)*out_raw.sz(out_raw.dims.subSpecs);
					nBadAvgTotal=nBadAvgTotal+nbadAverages;
					out_cs2=out_rm;
					iter=iter+1;
					disp([num2str(nbadAverages) ' bad averages removed on this iteration.']);
					disp([num2str(nBadAvgTotal) ' bad averages removed in total.']);
					%disp('Press any key to continue...');
					%pause
					close all;
				end
				figure('position',[0 50 560 420]);
				subplot(1,2,1);
				plot(out_cs.ppm,out_cs.specs);xlim([1 5]);
				xlabel('Frequency (ppm)');
				ylabel('Amplitude (a.u.)');
				title('Before removal of bad averages:');
				subplot(1,2,2);
				plot(out_rm.ppm,out_rm.specs);xlim([1 5]);
				xlabel('Frequency (ppm)');
				ylabel('Amplitude (a.u.)');
				title('After removal of bad averages:');
				figure('position',[0 550 560 420]);
				plot([1:length(metric{1})],metric{1},'.r',[1:length(metric{iter-1})],metric{iter-1},'x');
				xlabel('Scan Number');
				ylabel('Unlikeness metric (a.u.)');
				title('Unlikeness metrics before and after bad averages removal');
				legend('before','after');
				legend boxoff;
				sat		= 'y';
				%sat=input('are you satisfied with bad averages removal? ','s');
			end
		end
		
		
% 		% CBF: Create filenames for saving of processed output below
% 		outFileName		= [nameSpec, '_', strOVS, sprintf('_%.1f', nSD)];
% 		outFileNameW	= [nameW, '_w', '_', strOVS, sprintf('_%.1f', nSD)];
		
		% write a readme file to record the number of dropped avgs
		% CBF: Write readme file to output directory instead of data directory
		%fid=fopen([dirString '/readme.txt'],'w+');
		%fid		= fopen([outDirString outFileName '_readme.txt'],'w+');
		fid		= fopen([outDirString filename '/readme.txt'],'w+');
		fprintf(fid,'Original number of averages: \t%5.6f',out_raw.sz(out_raw.dims.averages)*2);
		disp(['Original number of averages:  ' num2str(out_raw.sz(out_raw.dims.averages)*2)]);
		fprintf(fid,'\nNumber of bad Averages removed:  \t%5.6f',nBadAvgTotal);
		disp(['Number of bad averages removed:  ' num2str(nBadAvgTotal)]);
		fprintf(fid,'\nNumber of remaining averages in processed dataset:  \t%5.6f',out_rm.sz(out_rm.dims.averages)*2);
		disp(['Number of remaining averages in processed dataset:  ' num2str(out_rm.sz(out_rm.dims.averages)*2)]);
		fclose(fid);
		
		close all;
		
		%Now remove the entries from fs_ai and phs_ai that correspond to
		%the bad-averages that were removed.
		BadAvgMask=zeros(length(fs_ai),1);
		BadAvgMask(allBadAverages)=1;
		size(BadAvgMask)
		size(fs_ai)
		fs_ai=fs_ai(~BadAvgMask);
		phs_ai=phs_ai(~BadAvgMask);
		
		%%%%%%%%%%%%%%%%%%%%END OF BAD AVERAGES REMOVAL%%%%%%%%%%%%%%%%%%%%
		
		
		%% CBF: If minimum user input selected, frequency drift correction is not optional
		% now align averages;
		% CBF: Ask for user input for frequency correction only, if minimum user input is NOT
		% selected
		driftCorr		= 'y';
		if strcmp(strMinUserIn,'n') || strcmp(strMinUserIn,'N')
			driftCorr=input('Would you like to perform frequency drift correction?  ','s');
		end
		if driftCorr=='n'|| driftCorr=='N'
			out_aa		= out_rm;
			if with_water
				out_w_aa	= out_w_cs;
			end
		else
			% CBF: Align averages in time or in frequency domain; if in time domain, obtain tmax
			% for drift correction as user input
			% (in _auto version of this routine, i.e. with almost no user input, tmax is also
			% still a user input)
			sat			= 'n'
			while sat=='n' || sat=='N'
				out_rm2			= out_rm;
				fsPoly			= 100;
				phsPoly			= 1000;
				fsCum			= fs_ai;fsPoly
				phsCum			= phs_ai;
				iter			= 1;
				while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<iterin
					close all
					if aaDomain=='t' || aaDomain=='T'
						tmax=input('input tmax for drift correction: ');
						[out_aa,fs,phs]=op_alignAverages(out_rm2,tmax,'n');
					elseif aaDomain=='f' || aaDomain=='F'
						tmax=tmaxin+0.04*randn(1);
						fmin=1.8+0.1*randn(1);
						fmaxarray=[2.4,2.85,3.35,4.2,4.4,5.2];
						fmax=fmaxarray(randi(6,1))
						[out_aa,fs,phs]=op_alignAverages_fd(out_rm2,fmin,fmax,tmax,'n');
					end
					if with_water
						[out_w_aa,fs_w,phs_w]=op_alignAverages(out_w_cs,5*tmax,'n');
					end
					
					fsCum=fsCum+fs;
					phsCum=phsCum+phs;
					fsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',fs,1)
					phsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',phs,1)
					iter
					out_rm2=out_aa;
					if with_water
						out_w_cs=out_w_aa;
					end
					iter=iter+1;
				end
				%Now plot the cumulative frequency drift correction:
				figure('position',[0 550 1125 420]);
				subplot(2,1,1);
				plot(out_rm.ppm,out_rm.specs);xlim([0 6]);
				xlabel('Frequency (ppm)');
				ylabel('Amplitude (a.u.)');
				title('Before drift correction:');
				subplot(2,1,2);
				plot(out_aa.ppm,out_aa.specs);xlim([0 6]);
				xlabel('Frequency (ppm)');
				ylabel('Amplitude (a.u.)');
				title('After drift correction:');
				
				figure('position',[0 50 560 420]);
				plot([1:out_aa.sz(out_aa.dims.averages)],phsCum);
				xlabel('Scan Number');
				ylabel('Phase Drift (deg.)');
				title('Estimated Phase Drift');
				
				figure('position',[570 50 560 420]);
				plot([1:out_aa.sz(out_aa.dims.averages)],fsCum);
				xlabel('Scan Number');
				ylabel('Frequency Drift (Hz)');
				title('Estimated Frequency Drift');
				
				sat		= 'y';
				%sat=input('Are you satisfied with the drift correction? ','s');
			end
			
		end
		
		
		%% now do the averaging and left shift to get rid of first order phase:
		out_av=op_leftshift(op_averaging(out_aa),out_aa.pointsToLeftshift);
		if with_water
			out_w_av=op_leftshift(op_averaging(out_w_aa),out_w_aa.pointsToLeftshift);
		end
		
		
		%% CBF: Do manual phase correction only, if minimum user input is NOT selected
		if strcmp(strMinUserIn,'n') || strcmp(strMinUserIn,'N')
			%Do a manual phase correction:
			SpecTool(out_av,0.05,-2,7);
			ph0		= input('input 0 order phase correction: ');
			ph1		= input('input 1st order phase correction: ');
		else
			% Set zero-order and first-order phases to zero for LCModel processing
			ph0		= 0;
			ph1		= 0;
		end
		% Apply phase correction
		out			= op_addphase(out_av,ph0,ph1);
		out_noproc	= op_addphase(out_noproc,ph0,ph1);
		
		
		% Now do a manual phase correction on the water unsuppressed data, if minimum user input
		% is not selected
		if with_water
			if strcmp(strMinUserIn,'n') || strcmp(strMinUserIn,'N')
				SpecTool(out_w_av,0.05,-2,7);
				ph0_w	= input('input 0 order phase correction: ');
				ph1_w	= input('input 1st order phase correction: ');
			else
				% Set zero-order and first-order phases to zero for LCModel processing
				ph0_w	= 0;
				ph1_w	= 0;
			end
			% Apply phase correction
			out_w			= op_addphase(out_w_av,ph0_w,ph1_w);
			out_w_noproc	= op_addphase(out_w_noproc,ph0_w,ph1_w);
		end
		
		
		%% Auto phase correction and frequency shifting for peak alignment
		
		
		
		%% Make figure(s) to show the final spectra and save results
		
		% Can be realized by calling a separate routine?
		% Use code from run_megapressproc_CBF.m
		
		
		
		
		
	case 'MEGA-PRESS'
		disp('Coming soon ...');
		
		% Figure display and saving
		% Can be realized by calling a separate routine?
		% Use code from run_specialproc_CBF.m
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType);
end



