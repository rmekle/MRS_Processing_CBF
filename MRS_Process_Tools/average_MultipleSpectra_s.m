%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% average_MultipleSpectra_s.m
%
%% Script to average multiple spectra of the same format
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2023;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility
sFunctionName		= 'average_MultipleSpectra_s';


%% Set parameters for loading and averaging multiple spectra
% Input data directory
dirString_In_Base		= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Trauma_MMs_Analysis/';
dirString_In_AddOn1		= 'MMs_PCG_IMA_FID-A_SD_3_2_ECCref_ls1_SR1';
%dirString_In_AddOn1		= 'MMs_PCG_IMA_FID-A_SD_3_2_ECCref_ls1_SR2';
dirString_In_AddOn2		= '_MRS_jmrui';
dirString_In			= [dirString_In_Base, dirString_In_AddOn1, dirString_In_AddOn2, filesep];
filename_In				= '';

% Data format and # of input data files
dataType_MRS			= 'mrs_ref';		% 'mrs_w_ref';		'mrs_w';	% 'mrs_ref';
dataFormat_In			= 'mrui';
structFileListing		= dir([dirString_In, '*.', dataFormat_In]);
noEntriesListing		= length( structFileListing );

% Output data directory and output data format
dirString_Out			= [dirString_In, 'Avgd/'];
dataFormat_Out			= 'mrui';

% If output directory is non-existent, create it
if ~exist( dirString_Out, 'dir' )
	mkdir(dirString_Out);
% else
% 	% Output directory already exists, ask user whether to overwrite or not
% 	% (should help to avoid accidentally overwriting previosuly processed data)
% 	strOverwrite	= input('\n\nDo you want to overwrite the existing output directory (y/n)?  ','s');
% 	if strOverwrite == 'n' || strOverwrite  == 'N'
% 		fprintf('\n%s: Already existing output directory is not overwritten! Processing is aborted!\n\n', sFunctionName)
% 		return;
% 	else
% 		fprintf('\n%s: Already existing output directory is overwritten! Processing is continued!\n\n', sFunctionName)
% 	end		% End of if strOverwrite == 'n' || strOverwrite  == 'N'
end		% End of if ~exist( dirString_Out, 'dir' )


%% Parameters for processing
% Parameters for aligning multiple spectra with each other
bAlignSpectra			= 1;

% Parameters for spectral registration (aligning of averages/frequency and phase drift
% correction) performed in either frequency or time domain
%driftCorr			= 'y';		% 'y';		'n';
iterin					= 20;
aaDomain				= 'f';		% 'f';		't';
tmaxin					= 0.2;		% 0.2;		0.1;
bTmaxset				= 1;
ppmOption				= 6;
medin					= 'y';		% 'y';	'n';	'a';	'ref';
refin					= 'f';		% 'f';	'a';
modein					= 'fp';		% 'fp';		''f';		'p';
strSpecReg				= 'SRs8';	% To distinguish settings for spectral registration

% Set parameters for drift correction depending on type of data, i.e. whether MRS
% data is spectrum or water signal
% NOTE: Check whether aligning of averages in frequency domain works, if the MR
% spectrum is water signal itself; if not, simply align averages in time domain
switch dataType_MRS
	case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
		% MR spectrum is provided together without or with unsuppressed water
		% signal and/or with reference scans
		%ppmmin_fix			= 1.6;		% 1.6;		1.8;
		switch ppmOption
			case 1
				% For MR spectra
				ppmmin_fix			= 1.6;		% 1.6;		1.8;
				ppmmaxarray_fix		= [2.4,2.85,3.35,4.2,4.4,5.2];
			case 2
				% For MR spectra
				ppmmin_fix			= 1.6;
				ppmmaxarray_fix		= [3.5; 4.0; 5.5];
			case 3
				ppmmin_fix			= 4.2;
				ppmmaxarray_fix		= [5.5 5.5 5.2];
			case 4
				% Wide range to always inlcude water resonance
				ppmmin_fix			= 1.6;
				ppmmaxarray_fix		= [5.5 5.5 5.2];
			case 5
				% For MMs signals
				ppmmin_fix			= 0.2;
				ppmmaxarray_fix		= [3.35,4.2,4.4];
			case 6
				% For MMs signals
				ppmmin_fix			= 0.2;
				ppmmaxarray_fix		= [3.35,4.0,4.1];

			otherwise
				error('%s: Unknown ppmOption = %d!', sFunctionName, ppmOption);
		end			% End of switch ppmOption
	case {'water', 'water_ref'}
		% MR spectrum is water signal itself without or with reference scans
		ppmmin_fix			= 4.2;
		ppmmaxarray_fix		= [5.5 5.5 5.2];

	otherwise
		error('%s: Unknown MRS dataType_MRS = %s!', sFunctionName, dataType_MRS);
end		% End of switch dataType_MRS
%alignSS				= 2;		% For aligning subspectra (e.g. in SPECIAL)
plotSwitch				= 1;
reportSwitch			= 1;

% Select indices of spectra to be averaged
indSpectra_All			= [1:noEntriesListing];
%indSpectra_Sel			= [1:noEntriesListing];
indSpectra_Sel			= [3:noEntriesListing];		% [1:2];	[3:noEntriesListing];
noSpectra_All			= length(indSpectra_All);
noSpectra_Sel			= length(indSpectra_Sel);
indexStart				= 1;
indexStep				= 1;


%% Determine output filename
filename_Out_Base		= dirString_In_AddOn1;
if noSpectra_Sel == noSpectra_All
	% All spectra were averaged
	filename_Out_AddOn1		= sprintf('_All%d_avg', noSpectra_All);
else
	% Only selected spectra were averaged
	filename_Out_AddOn1		= sprintf('_Sel%d_avg', noSpectra_Sel);
end		% End of if noSpectra_Sel == noSpectra_All

% Indicate alignment of multiple spectra in output filename
if bAlignSpectra
	filename_Out_AddOn2		= sprintf('_align_%s_%s', aaDomain, strSpecReg);
else
	filename_Out_AddOn2		= '';
end
filename_Out			= [filename_Out_Base, filename_Out_AddOn1, filename_Out_AddOn2];

% Set subdirectory for report (figures)
reportFigDirStr			= [filename_Out, filesep];
reportDir				= [dirString_Out reportFigDirStr];

% Make a new directory for the output report figures, if not already existent,
% and if desired
if reportSwitch == 1
	%     reportDirStr		= [nameSpec '_report/'];
	%     reportFigDirStr		= [nameSpec '_report/figs/'];
	% 	if ~exist([outDirString reportDirStr], 'dir' )
	% 		mkdir([outDirString reportDirStr]);
	% 	end
	if ~exist( reportDir, 'dir' )
		mkdir(reportDir);
	else
		% Output directory already exists, ask user whether to overwrite or not
		% (should help to avoid accidentally overwriting previosuly processed data)
		strOverwrite	= input('\n\nDo you want to overwrite the existing report directory (y/n)?  ','s');
		if strOverwrite == 'n' || strOverwrite  == 'N'
			fprintf('\n%s: Already existing report directory is not overwritten! Processing is aborted!\n\n', sFunctionName)
			return;
		else
			fprintf('\n%s: Already existing report directory is overwritten! Processing is continued!\n\n', sFunctionName)
		end		% End of if strOverwrite == 'n' || strOverwrite  == 'N'
	end		% End of if ~exist( reportDir, 'dir' )
end		% End of if reportSwitch == 1


%% Set parameters for figure display
% For figure(s) of (avaraged) MR spectra
strTitle				= 'Averaged 3T MMs Spectra';
strFigFile_Add			= '';

% For figures of aligned MR spectra
% h		= figure('position', [left bottom width height]);
fig_left	= 20;
fig_bottom	= 50;
fig_width	= 560;
fig_height	= 420;
fig_dist_b	= 510;
fig_dist_l	= 570;


%% Load selected datasets (spectra) according to the selected data format
% Here, datasets are loaded into FID-A toolkit data structures for MR spectra
ceSpectra_Sel	= cell(noSpectra_Sel, 1);
switch dataFormat_In
	case 'dat'
		error('%s: ERROR: Averaging of multiple spectra not yet implemented for dataFormat_In = %s!', sFunctionName, dataFormat_In);
	case 'IMA'
		error('%s: ERROR: Averaging of multiple spectra not yet implemented for dataFormat_In = %s!', sFunctionName, dataFormat_In);
	case 'mrui'
		for ind=indexStart : indexStep : noSpectra_Sel
			filename_In				= structFileListing(indSpectra_Sel(ind)).name;
			ceSpectra_Sel{ind}		= io_loadjmrui( fullfile(dirString_In, filename_In) );
		end		% End of ind=indexStart : indexStep : noSpectra_Selected

	otherwise
		error('%s: ERROR: Unknown data format = %s!', sFunctionName, dataFormat_In);
end		% End of switch dataFormat_In


%% Align multiple spectra in time domain or frequency domain, if selected
% Determine # of initial values for ppmmax
noVals_ppmmax_fix		= length(ppmmaxarray_fix);
if bAlignSpectra
	fprintf('\n\n');
	%out_rm2		= out_rm;
	ceSpectra_Sel2		= ceSpectra_Sel;
	fsPoly				= 100;
	phsPoly				= 1000;
	fscum				= zeros(noSpectra_Sel,1);
	phscum				= zeros(noSpectra_Sel,1);
	iter				= 1;
	while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<iterin
		fprintf(1, 'Aligning multiple spectra (scans): iteration iter = %d\n', iter);
		close all
		%tmax			= tmaxin+0.03*randn(1);	 % From run_pressProc_auto.m
		tmax			= tmaxin+0.04*randn(1);  % From run_specialProc_auto.m
		ppmmin			= ppmmin_fix+0.1*randn(1);
		switch noVals_ppmmax_fix
			case 3
				% Generate array of ppmmax values using random number variations
				ppmmaxarray		= [ppmmaxarray_fix(1)+0.1*randn(1,2),ppmmaxarray_fix(2)+0.1*randn(1,3),ppmmaxarray_fix(3)+0.1*randn(1,1)];
			case 6
				% Generate array of ppmmax values using given values
				ppmmaxarray		= ppmmaxarray_fix;

			otherwise
				error('%s: No option for noVals_ppmmax_fix = %d!', sFunctionName, noVals_ppmmax_fix);
		end		% End of switch noVals_ppmmax_fix
		% Select value for ppmmax used in this iteration
		iamax			= length(ppmmaxarray);
		ppmmax			= ppmmaxarray(randi(iamax,1));
		fprintf('\n');
		fprintf('ppmmaxarray = [%s]\n', join(string(ppmmaxarray), ' '));
		fprintf('ppmmin = %f \t ppmmax = %f\n\n\n', ppmmin, ppmmax);

		switch aaDomain
			case 't'
				% Perform alignment of spectra (scans) in time domain
				% using given values for tmax and other parameters
				[ceSpectra_Sel_as,phs,fs]		= op_alignAllScans(ceSpectra_Sel2,tmax,refin,modein);
			case 'f'
				% Perform alignment of spectra (scans) in frequency domain
				[ceSpectra_Sel_as,phs,fs]		= op_alignAllScans_fd(ceSpectra_Sel2,ppmmin,ppmmax,tmax,refin,modein);

			otherwise
				error('%s: ERROR: alignDomain %s for spectra (scans) not recognized!', sFunctionName, aaDomain);
		end

		fsPoly		= polyfit([1:noSpectra_Sel]',fs,1)
		phsPoly		= polyfit([1:noSpectra_Sel]',phs,1)

		fscum		= fscum+fs;
		phscum		= phscum+phs;

		%if driftCorr=='y' || driftCorr=='Y'
		if bAlignSpectra
			ceSpectra_Sel2		= ceSpectra_Sel_as;
		end
		iter			= iter+1;
	end		% End of while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<iterin

	% For plotting, extract all spectra from a cell array of spectra into one array
	% (only works, if each element of the cell array only contains one spectrum)
	%specsSpectra_Sel	= zeros(length(ceSpectra_Sel{1}.specs), noSpectra_Sel);
	specsSpectra_Sel		= zeros(ceSpectra_Sel{1}.sz(1), noSpectra_Sel);
	specsSpectra_Sel_as		= zeros(ceSpectra_Sel_as{1}.sz(1), noSpectra_Sel);
	for j=1 : 1 : noSpectra_Sel
		specsSpectra_Sel(:, j)		= ceSpectra_Sel{j}.specs;
		specsSpectra_Sel_as(:, j)	= ceSpectra_Sel_as{j}.specs;
	end		% End of for J=1 : 1 noSpectra_Sel

	% Only display figure(s), if selected
	if plotSwitch == 1
		h5	= figure('position',[fig_left fig_bottom fig_width fig_height]);
	else
		h5	= figure('visible','off');
	end
	subplot(1,2,1);
	plot(ceSpectra_Sel{1}.ppm,real(specsSpectra_Sel(:,:)));xlim([0.2 4.5]);
	set(gca,'FontSize',8);
	set(gca,'XDir','reverse');
	xlabel('Frequency (ppm)','FontSize',10);
	ylabel('Amplitude(a.u.)','FontSize',10);
	title('Before','FontSize',12);
	box off;
	subplot(1,2,2);
	plot(ceSpectra_Sel_as{1}.ppm,real(specsSpectra_Sel_as(:,:)));xlim([0.2 4.5]);
	set(gca,'FontSize',8);
	set(gca,'XDir','reverse');
	xlabel('Frequency (ppm)','FontSize',10);
	ylabel('Amplitude(a.u.)','FontSize',10);
	title('After','FontSize',12);
	box off;
	set(h5,'PaperUnits','centimeters');
	set(h5,'PaperPosition',[0 0 20 15]);

	% Only display figure(s), if selected
	if plotSwitch == 1
		h6	= figure('position',[fig_left (fig_bottom+fig_dist_b) fig_width fig_height]);
	else
		h6	= figure('visible','off');
	end
	plot([1:noSpectra_Sel],fscum,'.-','LineWidth',2);
	set(gca,'FontSize',8);
	xlabel('Scan Number','FontSize',10);
	ylabel('Frequency Drift [Hz]','FontSize',10);
	box off;
	legend('Frequency Drift','Location','SouthEast');
	legend boxoff;
	title('Estimated Frequency Drift','FontSize',12);
	set(h6,'PaperUnits','centimeters');
	set(h6,'PaperPosition',[0 0 10 10]);

	% Only display figure(s), if selected
	if plotSwitch == 1
		h7	= figure('position',[(fig_left+fig_dist_l) fig_bottom fig_width fig_height]);
	else
		h7	= figure('visible','off');
	end
	plot([1:noSpectra_Sel],phscum,'.-','LineWidth',2);
	set(gca,'FontSize',8);
	xlabel('Scan Number','FontSize',10);
	ylabel('Phase Drift [Deg.]','FontSize',10);
	box off;
	legend('Phase Drift','Location','SouthEast');
	legend boxoff;
	title('Estimated Phase Drift','FontSize',12);
	set(h7,'PaperUnits','centimeters');
	set(h7,'PaperPosition',[0 0 10 10]);
	
	% Save figures, if report switch is turned ON
	if reportSwitch == 1
		saveas(h5,[reportDir 'alignScans_prePostFig'],'jpg');
		saveas(h5,[reportDir 'alignScans_prePostFig'],'fig')
		saveas(h6,[reportDir 'freqDriftFig'],'jpg');
		saveas(h6,[reportDir 'freqDriftFig'],'fig');
		saveas(h7,[reportDir 'phaseDriftFig'],'jpg');
		saveas(h7,[reportDir 'phaseDriftFig'],'fig');
	end
	%close(h7);
	%close(h6);
	%close(h5);
	% Clear all figure handles from workspace to avoid warning when saving workspace
	vars		= {'h5','h6','h7'};
	clear(vars{:});

	% Compute total frequency and phase drifts
	totalFreqDrift		= mean(max(fscum)-min(fscum));
	totalPhaseDrift		= mean(max(phscum)-min(phscum));
	%close all

	% Pass on aligned cell array of selected spectra (scans) to ouput cell array
	ceSpectra_Sel_out	= ceSpectra_Sel_as;
else
	% Pass on original cell array of selected spectra (scans) to ouput cell array
	ceSpectra_Sel_out	= ceSpectra_Sel;
end		% End of if bAlignSpectra


%% Average multiple spectra
% First, add all spectra using a FID-A routine and then divide results by the # of
% selected spectra
out					= ceSpectra_Sel_out{1};
for j=2 : 1 : noSpectra_Sel
	out		= op_addScans(out, ceSpectra_Sel_out{j}, 0);
end
% Divide FIDs of added spectra by # of selected spectra for averaging
out.fids	= out.fids./noSpectra_Sel;

% Re-calculate spectra (specs) using fft
% Flags of FID_A data struct do not need to be updated for these last steps
out.specs	= fftshift(ifft(out.fids,[],out.dims.t),out.dims.t);


%% Plot averaged spectra
if (plot_MR_Spectrum_s(out, dataType_MRS, reportDir, filename_Out, strTitle, strFigFile_Add, 600, 1)) == -1
	error('%s: ERROR: Plotting of MR spectrum and/or saving of figures was NOT successful!', sFunctionName);
end


%% Save result in selected output data format
wrt				= 'y';		% 'y';		'n';
if wrt=='y' || wrt=='Y'
	fprintf('\n\n');
	fprintf('%s: Writing results to file ...\n\n', sFunctionName);
	switch dataFormat_Out
		case 'dat'
			error('%s: ERROR: Saving of averaged spectra not yet implemented for dataFormat_Out = %s!', sFunctionName, dataFormat_Out);
		case 'IMA'
			error('%s: ERROR: Saving of averaged spectra spectra not yet implemented for dataFormat_Out = %s!', sFunctionName, dataFormat_Out);
		case 'mrui'
			% jMRUI output format
			RF_jmrui = io_writejmrui(out,[dirString_Out, filename_Out, '.mrui']);

		otherwise
			error('%s: ERROR: Unknown dataFormat_Out = %s!', sFunctionName, dataFormat_Out);
	end		% End of switch dataFormat_Out
end		% End of if wrt=='y' || wrt=='Y'


%% Save variables of workspace to file
% Obtain current date and time in specific format
dt		= datestr(now,'yyyymmdd_HH_MM_SS');

% Save workspace into output directory (optional with user input)
% (Extension".mat" in filename explicitly required, so that Matlab can correctly load 
% workspace file with a "." in its filename)
%strSavedWorkspaceFileName		= ['workspace_', sFunctionName, '_', seqType_MRS, '_', dataType_MRS, '_', dt];
%strSavedWorkspaceFileNameFull	= [dirString_Out, strSavedWorkspaceFileName, sprintf('_SD_%.1f.mat', noSD_In)];
strSavedWorkspaceFileName		= [filename_Out, '__workspace_', dt];
strSavedWorkspaceFileNameFull	= [dirString_Out, strSavedWorkspaceFileName, '.mat'];
%strSaveWorkspace	= input('Would you like to save all variables of the workspace to file?  ', 's');
strSaveWorkspace	= 'y';
if strcmp(strSaveWorkspace,'y') || strcmp(strSaveWorkspace,'Y')
	save(strSavedWorkspaceFileNameFull);
end

