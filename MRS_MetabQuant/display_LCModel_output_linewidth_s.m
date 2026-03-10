% display_LCModel_output_linewidth_s.m
%
% Script to display output from LCModel in .coord file including spectrum, fit from 
% LCModel, fit residuals, and metabolite signals
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clc; clear; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bdwidth = 5; topbdwidth = 30; set(0,'Units','pixels'); scnsize = get(0,'ScreenSize');
% pos1	= [bdwidth, 2/3*scnsize(4) + bdwidth, scnsize(3)/2 - 2*bdwidth, scnsize(4)/2.2 - (topbdwidth + bdwidth)];
% pos2	= [pos1(1) + scnsize(3)/2,pos1(2),pos1(3),pos1(4)];
% pos3	= [bdwidth, bdwidth, scnsize(3)/2 - 2*bdwidth, scnsize(4)/2.2 - (topbdwidth + bdwidth)];
% pos4	= [pos1(1) + scnsize(3)/2, bdwidth, scnsize(3)/2 - 2*bdwidth, scnsize(4)/2.2 - (topbdwidth + bdwidth)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Set string for name of routine
sFunctionName		= 'display_LCModel_output_linewidth_s';
fprintf('\n\n');

%% Choose display/plotting options
% Same scale
strBoField							= '3T';			%	'3T';	'7T';
strUseSameScale						= 'YES';		%	'YES';	'NO';
strUseSameScaleFitPlot				= 'NO';
strShowSpectrumAndFit				= 'YES';
strShowSpectrumAndFitFigs			= 'YES';
strShowMetaboliteFits				= 'NO';
strLinewidthFit						= 'NO';
strAllLinewidthFitsWoutResiduals	= 'NO';
strShowLinewidthData				= 'NO';
strCloseLinewidthFits				= 'NO';
strSaveLinewidthFits				= 'YES';
strEstimateSignal					= 'YES';
strShowEstimateSignal				= 'NO';
strSaveEstimateSignal				= 'YES';


%% Set info for file to be opened

% Obtain Bo field strength and name of volunteer/patient via user input
% Bo field strength
%strBoField	= questdlg('Which Bo field strength was used?', ...
%		'Bo Field Strength', '3T', '7T', 'Other', '3T');
	
% Volunteer/Patient
% acellCases	= {'7T_fMRS_RL'; 'V1081_Acc'; 'V1079_Occipital'; 'V1111_2_Occipital';};
% acellCases	= {'AureliaP'; 'DanielleT'; 'GregoryA'; 'GuillaumeS'; 'KatrinP'; ...
%	'NilsR'; 'T2_NilsR';};
% acellCases	= {'AureliaP'; 'DanielleT'; 'GregoryA'; 'GuillaumeS'; 'KatrinP'; ...
% 	'NilsR'; 'MihaelaB';};

% [indCase, ok]	= listdlg('PromptString', 'Select a case:', 'SelectionMode', 'single', ...
%                 'ListString', acellCases);

     
% File selection via UI and user
% Move to selected directory
%cd '/home/allgemein/projects/7T_2Ch_Sf_Coil_fMRS_Right_Left/results/LCModel_Analysis_CasesNewLamp/';
%cd 'LCModel_ver631_150303_Cases_130619-140209_Left_2_140424Re_MRS_coord/';
file_ext				= '.coord'; 
strFilterSpec			= strcat('*', file_ext);
%[file_name, file_path]	= uigetfile(strFilterSpec, 'Select file to be analyzed')
%file_name				= '3T_SBAM_0088_20230413_meas_MID00262_PCG_Con8_Fig.coord';
%file_path				= '/home/mekler/CSB_NeuroRad/mekler/Ralf/Papers/ISMRM_2026_05/Figures_Abstract_MRS_MPM/Fig_1_MRS/Fig_1_MRS_Aux/3T_SBAM_0088_20230413_PCG_LCM_Out_PCG_ref_Quant_Con8_Fig/';
%file_name				= '3T_SBAM_0088_20230413_meas_MID00252_HC_Con8_Fig.coord';
%file_path				= '/home/mekler/CSB_NeuroRad/mekler/Ralf/Papers/ISMRM_2026_05/Figures_Abstract_MRS_MPM/Fig_1_MRS/Fig_1_MRS_Aux/3T_SBAM_0088_20230413_HC_LCM_Out_HC_ref_Quant_Con8_Fig/';
file_name				= '3T_TGA_21_1_20250420_DICOM_SVS_SLASER_DKD_HC_RIGHT_TE23_WS128_0028_wOVS_3_0_SR1_processed_lcm.coord';
file_path				= '/home/mekler/CSB_NeuroRad/mekler/Ralf/Papers/Papers_Co_Author/Charite_202510_Charite_Psych_MRS_Networks/3T_TGA_Data_Display/3T_MRS_TGA_Analysis/HC_IMA_FID-A_ls1_SD3_0_SR1_ECCref/HC_LCM_Out_HC_ref_Quant_Conc8_43206/';
fileName				= [file_path file_name];
caseName				= file_name(1:(strfind(file_name, file_ext) - 1));

% Allocate arrays
acellMetabolites		= cell(30, 1);
%sz_acellMetabolites = size(acellMetabolites)


%% Set parameters specific to acquisition, macromolecules (MMs), LCModel file(s), and plotting
seqType_MRS		= 'sLASER';		% 'SPECIAL';	% 'MEGA-PRESS'; % 'sLASER';
bIs_fMRS		= 0;
bAcquired_MMs	= 1;

% Obtain info about the coord file
% Assume that data values are stored in 10 columns (by LCModel)
noColumns		= 10;

% Set plotting resolution
resolution		= 600;


%% Set parameters for (output) file naming
saveDir         = file_path;
strSaveName		= '';
outDir			= file_path;
outFileName     = 'out_LCModel';
out_ext			= '.txt';

% Select output filename for parameter estimates according to selected options
% Linewidth fits
outFileName_AddOn1	= '';
if( strcmp(strSaveLinewidthFits, 'YES') )
	outFileName_AddOn1	= '_linewidthFits_Cr';
end		% End of if( strcmp(strSaveLinewidthFits, 'YES') )

% Select output filename for parameter estimates according to selected options
% Signal estimates
outFileName_AddOn2	= '';
if( strcmp(strSaveEstimateSignal, 'YES') )
	outFileName_AddOn2	= '_signalEstimates';
end		% End of if( strcmp(strSaveEstimateSignal, 'YES') )

% Complete output filename and full absolute output filename
outFileName			= [outFileName, outFileName_AddOn1, outFileName_AddOn2];
fullOutFileName		= strcat(outDir, outFileName, out_ext);


%% Specific case: 7T_fMRS_RL
if bIs_fMRS
	% Set output file information for saving
	% Add "Left" or "Right" to filename, if directory name contains one of these words
	% (last part of absolute path, i.e. part following second last slash)
	slashInd		= strfind(file_path, '/');
	strDirectory	= file_path( (slashInd(length(slashInd)-1)+1):end );
	if ~isempty(strfind( strDirectory, 'Left' ))
		outFileName		= strcat(outFileName, '_7T_fMRS_RL_Left');
	end
	if ~isempty(strfind( strDirectory, 'Right' ))
		outFileName		= strcat(outFileName, '_7T_fMRS_RL_Right');
	end
	%saveDir         = strcat('/echo2/allgemein/tmp/RalfM/Z_MRS_Processing/7T_fMRS_RL/');
	%strSaveName     = strcat('_7T_fMRS_RL_', caseName);
end		% End of if bIs_fMRS


%% Use different directory, if all data are plotted at same scale
if( strcmp(strUseSameScale, 'YES') )
	parentDir	= saveDir;
	%saveDir		= strcat(saveDir, 'SameScale\');
    saveDir		= strcat(saveDir, 'SameScale/');
end % End of if( strcmp(strUseSameScale, 'YES') )


%% Determine # of macromolecules (MMs)
if bAcquired_MMs
	% Acquired macromolecules = 'mac' in basis set
	noMacromolecules	= 1;
else
	% Macromolecules modelled by LCModel depending on sequence type
	if strcmp(seqType_MRS, 'MEGA-PRESS')
		% For MEGA-PRESS
		% Simulated macromolecules = 'MM09'; 'MM20';'MM12';'MM14'; 'MM17'; 'MM30'
		noMacromolecules	= 6;
	else
		% Simulated macromolecules = 'MM09'; 'MM20';'MM12';'MM14'; 'MM17'
		noMacromolecules	= 5;
	end			% End of if strcmp(seqType_MRS, 'MEGA-PRESS')
end		% End of if bAcquired_MMs


%% Set parameters specific to Bo field strength
switch( strBoField )
    case '3T'
        % 3T
        % 1 ppm in Hz at 3T (according to LCModel)
        ppm_Hz			= 123.26; 		
    case '7T'
        % 7T
        % 1 ppm in Hz at 7T (according to LCModel)
        ppm_Hz			= 297.21;
        
    otherwise
        error('\n%s: No option for chosen field strength strBoField = %s!\n\n', sFunctionName, strBoField);
end		% End of switch( strBoField )
       
% Set parameters specific to sequence type
switch( seqType_MRS )
	case {'PRESS', 'STEAM', 'sLASER', 'SPECIAL'}
		strSaveName		= '';
	case 'MEGA-PRESS'
		strSaveName		= '_M_16_1';		% '_M_01_2';

	otherwise
		error('\n%s: No option for chosen sequence type seqType_MRS = %s!\n\n', sFunctionName, seqType_MRS);
end		% End of switch( seqType_MRS )


%% Read information about metabolite fits from file
% % Read the entire text file into string
% fileText = fileread(fileName);
% 
% % Search for the line of file that includes 'points on ppm-axis'
% % each line is separated by a newline ('\n')
% expr            = '[^\n]*points on ppm-axis[^\n]*';
% fileread_info   = regexp(fileText,expr,'match')

% Open file and search for the line of file that includes 'points on ppm-axis' and count
% lines that were read up to this line to determine start index of data in file
% (data starts in subsequent line)
[fid, errmsg]		= fopen(fileName);
if fid == -1
	error('%s: Invalid file identifier fid = %d with error message = %s!\n\n', sFunctionName, fid, errmsg);
end
tline           = fgetl(fid);
lLineCounter    = 1;
%disp(tline)

while isempty(strfind(tline, 'points on ppm-axis'))
    tline           = fgetl(fid);
    lLineCounter    = lLineCounter + 1;
    %disp(tline)
end
%disp(tline)
%disp(lLineCounter)

% Extract # of data values for each type of data from last line read and determine start
% index into file after which data is stored
%noValues                = 1218
noValues                = sscanf(tline, '%d');
fileIndexStartValues    = lLineCounter;

% Determine end of data values in file (this is the first line after all data) and extract
% all metabolite names that are listed in the coord file
lMetaboliteInd        = 0;
while isempty(strfind(tline, 'following diagnostic table'))
    tline           = fgetl(fid);
    lLineCounter    = lLineCounter + 1;
    if ~isempty( strfind(tline, 'Conc.') ) % Line with metabolite name
        lMetaboliteInd                      = lMetaboliteInd + 1;
        acellMetabolites{lMetaboliteInd}    = strtok(tline);      
    end
end		% End of while isempty(strfind(tline, 'following diagnostic table'))
fileIndexEndValues  = lLineCounter;

% Remove empty cell elements
acellMetabolites = acellMetabolites(~cellfun('isempty',acellMetabolites)); 
%sz_acellMetabolites = size(acellMetabolites)

% Close input file
if( fclose(fid) == -1 )
	error('%s: Error closing input file!', sFunctionName);
end

% Find strings for metabolite fits for Cr and PCr
lInd_Cr		= find(ismember(acellMetabolites,'Cr'));
lInd_PCr	= find(ismember(acellMetabolites,'PCr'));

% Find strings for metabolite fits for GABA and Glu
% (essential for editing using 'MEGA-PRESS'; might be empty, if string is not found)
lInd_GABA	= find(ismember(acellMetabolites,'GABA'));
lInd_Glu	= find(ismember(acellMetabolites,'Glu'));

% Determine step size for indices (# of lines for data of each metabolite plus line of 
% description, such as "NY phased data points follow" or "Cr       Conc. = 2.98E+00") and
% determine indices into the file where specific data columns are stored and how many
% metabolites are included 
% (after start index of data, the first 4 sets of values included are ppm axis, phased
% data, fit to the data, and background; after that the values for metabolites follow)
stepIndices        = ceil(noValues/noColumns) + 1;
fileIndices        = fileIndexStartValues + [0; 1; 2; 3; (3+lInd_Cr); (3+lInd_PCr);] .* stepIndices;
if strcmp(seqType_MRS, 'MEGA-PRESS')
	% For MEGA-PRESS
	% Use indices for metabolite fits for GABA and Glu
	fileIndices        = fileIndexStartValues + [0; 1; 2; 3; (3+lInd_GABA); (3+lInd_Glu);] .* stepIndices
end		% End of if strcmp(seqType_MRS, 'MEGA-PRESS')
noMetabolites      = (fileIndexEndValues - (fileIndexStartValues + 4*stepIndices)) / stepIndices;

% Display some info
fprintf('\n\n');
fprintf('# of data points for each component = noValues = %d\n\n', noValues);
fprintf('Last line of data values = (fileIndexEndValues-1) = %d\n\n', (fileIndexEndValues-1));
fprintf('# of metabolites = noMetabolites = %d\n\n', noMetabolites);
disp('acellMetabolites = '); disp(acellMetabolites);


%% Open input file and read data
[fid, errmsg]		= fopen(fileName);
if fid == -1
	error('%s: Invalid file identifier fid = %d with error message = %s!\n\n', sFunctionName, fid, errmsg);
end
fid_pos		= 0;

%% fileIndices(1)  n points on ppm-axis = NY = 'noValues'
status			= fseek(fid, 0, 'bof');
[C, fid_pos]	= textscan(fid, '%f', noValues, 'headerlines', fileIndices(1));
data = C{1,1}; ppm = data;

%% fileIndices(2)   NY phased data points follow
status			= fseek(fid, 0, 'bof');
[C, fid_pos]	= textscan(fid, '%f', noValues, 'headerlines', fileIndices(2));
data = C{1,1}; raw_spectrum = data;
%sz_data         = size(data)

%% fileIndices(3)  NY points of the fit to the follow
status			= fseek(fid, 0, 'bof');
[C, fid_pos]	= textscan(fid, '%f', noValues, 'headerlines', fileIndices(3));
data = C{1,1}; spectrum_fit = data;

%% fileIndices(4)   NY background values follow
status			= fseek(fid, 0, 'bof');
[C, fid_pos]	= textscan(fid, '%f', noValues, 'headerlines', fileIndices(4));
data = C{1,1}; background = data;

%% Select subsequent metabolite fits depending on sequence type
if strcmp(seqType_MRS, 'MEGA-PRESS')
	% For MEGA-PRESS
	% Use indices for metabolite fits for GABA and Glu
	%% fileIndices(5)	GABA values follow
	status			= fseek(fid, 0, 'bof');
	[C, fid_pos]	= textscan(fid, '%f', noValues, 'headerlines', fileIndices(5));
	data = C{1,1}; metabolite_GABA_withBG = data;

	%% fileIndices(6)	Glutamate (Glu) values follow
	status			= fseek(fid, 0, 'bof');
	[C, fid_pos]	= textscan(fid, '%f', noValues, 'headerlines', fileIndices(6));
	data = C{1,1}; metabolite_Glu_withBG = data;

else
	%% fileIndices(5)	Creatine (Cr) values follow
	status			= fseek(fid, 0, 'bof');
	[C, fid_pos]	= textscan(fid, '%f', noValues, 'headerlines', fileIndices(5));
	data = C{1,1}; metabolite_Cr_withBG = data;

	%% fileIndices(6)	Phosphocreatine (PCr) values follow
	status			= fseek(fid, 0, 'bof');
	[C, fid_pos]	= textscan(fid, '%f', noValues, 'headerlines', fileIndices(6));
	data = C{1,1}; metabolite_PCr_withBG = data;
end		% End of if strcmp(seqType_MRS, 'MEGA-PRESS')


%% Metabolites
%noMetabolites		= length(acellMetabolites)
indMetabolites		= zeros(noMetabolites, 1);
indMetabolites		= fileIndices(4) + [1:noMetabolites].' .* stepIndices;
fitsMetabolites		= zeros(noMetabolites, noValues);
sz_raw_spectrum     = size(raw_spectrum);
sz_background		= size(background);

% Read in data for fits of metabolites and subtract background (baseline) values
for lI=1 : 1 : noMetabolites
	%[lI indMetabolites(lI)]
	status			= fseek(fid, 0, 'bof');
	[C, fid_pos]	= textscan(fid, '%f', noValues, 'headerlines', indMetabolites(lI));
	data			= C{1,1};
    %sz_data_InLoop  = size(data)
    %sz_fitsMetabolites = size(fitsMetabolites(lI, :))
	fitsMetabolites(lI, :) = (data - background).';
end

% Close input file
if( fclose(fid) == -1 )
	error('%s: Error closing input file!', sFunctionName);
end

% Display some info
fprintf('\n\n');
fprintf('Size of raw_spectrum = sz_raw_spectrum = %d\t%d\n', sz_raw_spectrum);
fprintf('Size of background = sz_background = %d\t%d\n', sz_background);


%% Obtain minimum and maximum of ppm values only with the precision of two digits after
%  decimal point; for minimum use floor and for maximum use ceiling
min_ppmFloor	= floor( min(ppm)*100 ) / 100;
max_ppmCeil		= ceil( max(ppm)*100 ) / 100;

%%%%%%%%%%%%%     RESIDUALS    %%%%%%%%%%%%
residuals		= raw_spectrum - spectrum_fit;
min_spectrum	= min(raw_spectrum);
max_spectrum	= max(raw_spectrum);
SIrange			= (max_spectrum - min_spectrum);
ranges_figR		= [min_ppmFloor max_ppmCeil -0.25*SIrange 0.25*SIrange];
%ranges_figR		= [min(ppm) max(ppm) -0.25*SIrange 0.25*SIrange];
%ranges_fig		= [min(ppm) max(ppm) min_spectrum 1.25*max_spectrum];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine minimim and maximum of macromolecules (MMs)
if bAcquired_MMs
	% Acquired macromolecules = 'mac' in basis set
	% (assume that for acquired macromolecules are in first row of array for metabolite
	% fits)
	min_Mac		= min(fitsMetabolites(1, :));
	max_Mac		= max(fitsMetabolites(1, :));
else
	% Sum up all modelled macromolecule contributions, if not acquired
	% (assume that modelled macromolecules are in last rows of array for metabolite fits)
	fitSummedMM	= zeros(1, noValues);
	for lI=0 : 1 : (noMacromolecules-1)
		fitSummedMM		= fitSummedMM + fitsMetabolites((noMetabolites-lI), :);
	end
	min_Mac		= min(fitSummedMM);
	max_Mac		= max(fitSummedMM);
end		% End of if bAcquired_MMs


%% Plot spectrum, fit from LCModel, and fit residuals into one figure

% % Plot into one figure with 3 subplots
% figure;
% subplot(3,1,2)
% plot(ppm, raw_spectrum); set(gca,'Xdir','reverse'); title('RAW SPECTRUM'); axis(ranges_fig)
% 
% subplot(3,1,3)
% plot(ppm, spectrum_fit);set(gca,'Xdir','reverse'); title('SPECTRUM FIT'); axis(ranges_fig)
% 
% subplot(3,1,1)
% plot(ppm, residuals); set(gca,'Xdir','reverse'); title('RESIDUALS '); axis(ranges_figR)


% Init index into array of figure handles
indFigs	= 0;
% Set figure properties for spectrum and fit of spectrum
% ppm_range		= [min(ppm) max(ppm)]
font			= 12;
plotLineWidth	= 1.2;
ppmTextPos      = [0.9, -0.11, 0];
ppmPlotDelta    = 0.10;				% 0.05;		% 0.10;		% 0.15;
ppmPlotRange	= [(min_ppmFloor-ppmPlotDelta) (max_ppmCeil+ppmPlotDelta)];

% Determine a suitable plotting scale, also for the case that all data are plotted at
% same scale
if( min_spectrum < 0 )
	% Water signal or other part of spectrum negative
	spectrumPlotRange		= [(1.05*abs(min_spectrum)*sign(min_spectrum)) 1.05*max_spectrum];
	spectrumPlotRangeSame	= [(min_spectrum-0.05*max_spectrum) 1.05*max_spectrum];
else
	% Water signal and spectrum positive
	spectrumPlotRange		= [0 1.05*max_spectrum];
	spectrumPlotRangeSame	= [-0.05*max_spectrum 1.05*max_spectrum];
end
ranges_fig				= [ppmPlotRange spectrumPlotRange];
ranges_figSameScale		= [ppmPlotRange spectrumPlotRangeSame];

% Set figure properties for plot of residuals
font_figR			= 10;
plotFactor_figR		= 1.5;		% 1.05;		% 1.1;	% 1.5;
max_residuals		= max(abs(residuals));
ranges_figR			= [ppmPlotRange -plotFactor_figR*max_residuals plotFactor_figR*max_residuals];

% Display some info
fprintf('\n\n');
fprintf('min_spectrum = %.2f\t\tmax_spectrum = %.2f\n', min_spectrum, max_spectrum);
disp('ppmPlotRange = '); disp(ppmPlotRange);
disp('spectumPlotRange = '); disp(spectrumPlotRange);
disp('spectumPlotRangeSame = '); disp(spectrumPlotRangeSame);

% Plot into one figure with 2 subplots, if selected
if( strcmp(strShowSpectrumAndFit, 'YES') )
	indFigs				= indFigs + 1;
	h_figures(indFigs)	= figure;

	% Fit residuals
	%subplot(2,1,1);
	subplot('Position',[0.13 0.85 0.775 0.1]);
	hplot				= plot(ppm, residuals, 'k', 'LineWidth', 1.0);
	set(gca,'Xdir','reverse', 'XTick', [],  'FontName', 'Arial', ...
		'FontSize', font_figR, 'LineWidth', 0.7, 'FontWeight', 'bold', 'XColor', [0 0 0], ...
		'YColor', [0 0 0], 'TickDir', 'out', 'Box', 'on', 'Linewidth', 0.7, ...
		'YTickLabelMode', 'Manual', 'YTick', [-2000 0 2000], ...
		'YTickLabel', {'-2*10^3' '0' '2*10^3'}, ...
		'YGrid', 'on', 'GridLineStyle', '--', 'GridColor', [0.0 0.0 0.0], ...
		'GridAlpha', 0.5, 'GridLineWidth', 1.0);		%'YGrid', 'off');
	%if( strcmp(strUseSameScale, 'YES') )
	if( strcmp(strUseSameScaleFitPlot, 'YES') )
		axis(ranges_figSameScale);
	else
		axis(ranges_figR);
	end
	%title('Fit Residuals');
	% xlabel('ppm', 'FontSize', (font+4), 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
	% 	'Units', 'normalized', 'Position', ppmTextPos);

	% Spectrum and fit of spectrum
	%subplot(2,1,2);
	subplot('Position',[0.13 0.10 0.775 0.74]);
	hplot = plot(ppm, raw_spectrum, 'k', ppm, spectrum_fit, 'r', 'LineWidth', plotLineWidth);
	%hplot				= plot(ppm, raw_spectrum, 'k', 'LineWidth', plotLineWidth);
	set(gca,'Xdir','reverse', 'XTick', [0.5:0.5:5], 'XTickLabel', {'', '1', '', '2', '', '3', '', '4', '', '5'}, ...
		'FontSize', font, 'LineWidth', 2, 'FontWeight', 'bold', 'XColor', [0 0 0], ...
		'YColor', [0 0 0], 'TickDir', 'out', 'Box', 'off', 'FontName', 'Arial', ...
		'YTickLabelMode', 'Manual', 'YTick', [0 20000 60000 100000], ...
		'YTickLabel', {'0' '2*10^4' '6*10^4' '10*10^4'}, ...
		'YColor', 'none', ...		% Hide y-axis using set(gca, ...);
		'YGrid', 'off');
	%'YTickLabelMode', 'Manual', 'YTick', []);
	if( strcmp(strUseSameScale, 'YES') )
		axis(ranges_figSameScale);
	else
		axis(ranges_fig);
	end
	%title('Raw Spectrum');
	xlabel('ppm', 'FontSize', (font+2), 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
	% xlabel('ppm', 'FontSize', (font+4), 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
	% 	'Units', 'normalized', 'Position', ppmTextPos);
end		% End of if( strcmp(strShowSpectrumAndFit, 'YES') )


%% Plot spectrum, fit from LCModel, fit residuals, and baseline into separate figures, if selected

% Set figure properties
% ppm_range		= [min(ppm) max(ppm)]
font			= 16;
plotLineWidth	= 1.2;
ppmPlotRange	= [(min_ppmFloor-ppmPlotDelta) (max_ppmCeil+ppmPlotDelta)];
% % HACK for metabolite T2s
% if( indCase == 7 )
% 	ppmPlotRange = [(1.8-ppmPlotDelta) (max(ppm)+ppmPlotDelta)];
% end
% Determine a suitable plotting scale, also for the case that all data are plotted at
% same scale
if( min_spectrum < 0 )
	% Water signal negative
	spectrumPlotRange		= [(1.05*abs(min_spectrum)*sign(min_spectrum)) 1.05*max_spectrum];
	spectrumPlotRangeSame	= [(min_spectrum-0.05*max_spectrum) 1.05*max_spectrum];
else
% Water signal positive
	spectrumPlotRange		= [0 1.05*max_spectrum];
	spectrumPlotRangeSame	= [-0.05*max_spectrum 1.05*max_spectrum];
end
ranges_fig				= [ppmPlotRange spectrumPlotRange];
ranges_figSameScale		= [ppmPlotRange spectrumPlotRangeSame];

% Plot spectrum and fits into separate figures, if selected
if( strcmp(strShowSpectrumAndFitFigs, 'YES') )
	indFigs				= indFigs + 1;
	h_figures(indFigs)	= figure;
	% hplot = plot(ppm, raw_spectrum, 'k', ppm, spectrum_fit, '--r', 'LineWidth', plotLineWidth);
	hplot				= plot(ppm, raw_spectrum, 'k', 'LineWidth', plotLineWidth);
	set(gca,'Xdir','reverse', 'XTick', [0.5:0.5:5], 'XTickLabel', {'', '1', '', '2', '', '3', '', '4', '', '5'}, ...
		'FontSize', font, 'LineWidth', 2, 'FontWeight', 'bold', 'XColor', [0 0 0], ...
		'YColor', [0 0 0], 'TickDir', 'out', 'TickLength', [0.02, 0.025], 'Box', 'off', ...
		'YTickLabelMode', 'Manual', 'YTick', [0], 'YTickLabel', {'0'}, 'YGrid', 'on');
	%'YTickLabelMode', 'Manual', 'YTick', []);
	if( strcmp(strUseSameScale, 'YES') )
		axis(ranges_figSameScale);
	else
		axis (ranges_fig);
	end
	title('Raw Spectrum');
	xlabel('ppm', 'FontSize', (font+4), 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
		'Units', 'normalized', 'Position', ppmTextPos);

	% Plot fit of spectrum from LCModel
	indFigs				= indFigs + 1;
	h_figures(indFigs)	= figure;
	hplot				= plot(ppm, spectrum_fit, 'k', 'LineWidth', plotLineWidth);
	set(gca,'Xdir','reverse', 'XTick', [0.5:0.5:5], 'XTickLabel', {'', '1', '', '2', '', '3', '', '4', '', '5'}, ...
		'FontSize', font, 'LineWidth', 2, 'FontWeight', 'bold', 'XColor', [0 0 0], ...
		'YColor', [0 0 0], 'TickDir', 'out', 'TickLength', [0.02, 0.025], 'Box', 'off', ...
		'YTickLabelMode', 'Manual', 'YTick', [0], 'YTickLabel', {'0'}, 'YGrid', 'on');
	if( strcmp(strUseSameScale, 'YES') )
		axis(ranges_figSameScale);
	else
		axis (ranges_fig);
	end
	title('Fit of Spectrum');
	xlabel('ppm', 'FontSize', (font+4), 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
		'Units', 'normalized', 'Position', ppmTextPos);

	% Plot fit residuals
	ranges_figR			= [ppmPlotRange -0.1*SIrange 0.1*SIrange];
	indFigs				= indFigs + 1;
	h_figures(indFigs)	= figure;
	hplot				= plot(ppm, residuals, 'k', 'LineWidth', plotLineWidth);
	set(gca,'Xdir','reverse', 'XTick', [0.5:0.5:5], 'XTickLabel', {'', '1', '', '2', '', '3', '', '4', '', '5'}, ...
		'FontSize', font, 'LineWidth', 2, 'FontWeight', 'bold', 'XColor', [0 0 0], ...
		'YColor', [0 0 0], 'TickDir', 'out', 'TickLength', [0.02, 0.025], 'Box', 'off', ...
		'YTickLabelMode', 'Manual', 'YTick', [0], 'YTickLabel', {'0'}, 'YGrid', 'off');
	if( strcmp(strUseSameScale, 'YES') )
		axis(ranges_figSameScale);
	else
		axis (ranges_figR);
	end
	title('Fit Residuals');
	xlabel('ppm', 'FontSize', (font+4), 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
		'Units', 'normalized', 'Position', ppmTextPos);

	% Plot background (baseline) signal
	ranges_figBG		= [ppmPlotRange (sign(min(background))*(abs(min(background))+0.1*max(background))) (1.1*max(background))];
	indFigs				= indFigs + 1;
	h_figures(indFigs)	= figure;
	hplot				= plot(ppm, background, 'k', 'LineWidth', plotLineWidth);
	set(gca,'Xdir','reverse', 'XTick', [0.5:0.5:5], 'XTickLabel', {'', '1', '', '2', '', '3', '', '4', '', '5'}, ...
		'FontSize', font, 'LineWidth', 2, 'FontWeight', 'bold', 'XColor', [0 0 0], ...
		'YColor', [0 0 0], 'TickDir', 'out', 'TickLength', [0.02, 0.025], 'Box', 'off', ...
		'YTickLabelMode', 'Manual', 'YTick', [0], 'YTickLabel', {'0'}, 'YGrid', 'off');
	if( strcmp(strUseSameScale, 'YES') )
		axis(ranges_figSameScale);
    else
		axis(ranges_figBG);
	end
	title('Background (Baseline) Signal');
	xlabel('ppm', 'FontSize', (font+4), 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
		'Units', 'normalized', 'Position', ppmTextPos);
end		% End of if( strcmp(strShowSpectrumAndFitFigs, 'YES') )


%% Plot fits for metabolite signals into separate figures, if selected
if( strcmp(strShowMetaboliteFits, 'YES') )
	% Set index into figure handles to specific value
	% (to have a definite starting point for handles to figures for metabolite fits)
	firstMetabFig	= 6;
	indFigs			= firstMetabFig - 1;
	for lI=1 : 1 : noMetabolites
		min_Metabolite	= min(fitsMetabolites(lI, :));
		max_Metabolite	= max(fitsMetabolites(lI, :));
		% Assume that 'max_Metabolite' > 0
		if( min_Metabolite  > 0.1*max_Metabolite )
			ranges_figMetab		= [ppmPlotRange 0 (1.1*max_Metabolite)];
		else
			ranges_figMetab		= [ppmPlotRange ...
				(min_Metabolite-0.1*max_Metabolite) (1.1*max_Metabolite)];
		end
		indFigs				= indFigs + 1;
		h_figures(indFigs)	= figure;
		hplot				= plot(ppm, fitsMetabolites(lI, :), 'k', 'LineWidth', plotLineWidth);
		set(gca,'Xdir','reverse', 'XTick', [0.5:0.5:5], 'XTickLabel', {'', '1', '', '2', '', '3', '', '4', '', '5'}, ...
			'FontSize', font, 'LineWidth', 2, 'FontWeight', 'bold', 'XColor', [0 0 0], ...
			'YColor', [0 0 0], 'TickDir', 'out', 'TickLength', [0.02, 0.025], 'Box', 'off', ...
			'YTickLabelMode', 'Manual', 'YTick', [0], 'YTickLabel', {'0'}, 'YGrid', 'off');
		if( strcmp(strUseSameScale, 'YES') )
			axis(ranges_figSameScale);
		else
			axis(ranges_figMetab);
		end
		title(sprintf('Fit for %s\n', acellMetabolites{lI}));
		xlabel('ppm', 'FontSize', (font+4), 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
			'Units', 'normalized', 'Position', ppmTextPos);
	end

	% Plot summed up macromolecule contributions for 3T data
	if( strcmp(strBoField, '3T') )
		% Assume that 'max_Mac' > 0
		if( min_Mac  > 0.1*max_Mac )
			ranges_figSummedMM	= [ppmPlotRange 0 (1.1*max_Mac)];
		else
			ranges_figSummedMM	= [ppmPlotRange ...
				(min_Mac-0.1*max_Mac) (1.1*max_Mac)];
		end
		indFigs				= indFigs + 1;
		h_figures(indFigs)	= figure;
		indFigSummedMM		= indFigs;
		hplot				= plot(ppm, fitSummedMM, 'k', 'LineWidth', plotLineWidth);
		set(gca,'Xdir','reverse', 'XTick', [0.5:0.5:5], 'XTickLabel', {'', '1', '', '2', '', '3', '', '4', '', '5'}, ...
			'FontSize', font, 'LineWidth', 2, 'FontWeight', 'bold', 'XColor', [0 0 0], ...
			'YColor', [0 0 0], 'TickDir', 'out', 'TickLength', [0.02, 0.025], 'Box', 'off', ...
			'YTickLabelMode', 'Manual', 'YTick', [0], 'YTickLabel', {'0'}, 'YGrid', 'off');
		if( strcmp(strUseSameScale, 'YES') )
			axis(ranges_figSameScale);
		else
			axis(ranges_figSummedMM);
		end
		title('Fit for Summed MM (Macromolecules)');
		xlabel('ppm', 'FontSize', (font+4), 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
			'Units', 'normalized', 'Position', ppmTextPos);
	end
end		% End of if( strcmp(strShowMetaboliteFits, 'YES') )


%% Save plots as Matlab figures and EPS files and PNG files

% Save plots of spectrum, fit from LCModel, fit residuals, and baseline, if selected
% CHECK ON INDICES INTO ARRAY FOR FIGURE HANDLES
if( strcmp(strShowSpectrumAndFit, 'YES') )
	acellNames		= {'spectrum', 'fitOfSpectrum', 'fitResiduals', 'baseline'};
	answer = questdlg('Do you want to save the figures for spectrum, fit, residuals, and baseline?', ...
		'Saving of Spectrum Figures', 'Yes', 'Default', 'No', 'No');
	switch answer
		case 'Yes'
			fprintf('\n\nSaving of spectrum figures ...\n\n\n')
			% If directory for saving does not exist, create it
			if( ~exist(saveDir, 'dir') )
				[success,message,messageID] = mkdir(parentDir, 'SameScale');
				if( success == 0 )
					disp(message);
					disp(messageID);
				end
			end
			for lI=1 : 1 : 4
				figureName	= strcat(acellNames{lI}, strSaveName);
				saveFigure_s(h_figures(lI+1), saveDir, figureName, 'fig', resolution);
				saveFigure_s(h_figures(lI+1), saveDir, figureName, 'eps', resolution);
                saveFigure_s(h_figures(lI+1), saveDir, figureName, 'png', resolution);
				%saveFigure_s(h_figures(lI+1), saveDir, figureName, 'png', 300);
			end
		case 'Default'
			fprintf('\n\nSpectrum figures were not saved!\n\n\n');
		case 'No'
			fprintf('\n\nSpectrum figures were not saved!\n\n\n');
	end		% End of switch answer

	% Save plots of fits for metabolite signals, if selected
	% CHECK ON INDICES INTO ARRAY FOR FIGURE HANDLES
	if( strcmp(strShowMetaboliteFits, 'YES') )
		answer = questdlg('Do you want to save the figures for the metabolite fits?', ...
			'Saving of Metabolite Fit Figures', 'Yes', 'Default', 'No', 'No');
		switch answer
			case 'Yes'
				fprintf('\n\nSaving of metabolite fit figures ...\n\n\n')
				% If directory for saving does not exist, create it
				if( ~exist(saveDir, 'dir') )
					[success,message,messageID] = mkdir(parentDir, 'SameScale');
					if( success == 0 )
						disp(message);
						disp(messageID);
					end
				end
				for lI=1 : 1 : noMetabolites
					figureName	= strcat('fitSig_', acellMetabolites{lI}, strSaveName);
					saveFigure_s(h_figures((lI-1)+firstMetabFig), saveDir, figureName, 'fig', resolution);
					saveFigure_s(h_figures((lI-1)+firstMetabFig), saveDir, figureName, 'eps', resolution);
                    saveFigure_s(h_figures((lI-1)+firstMetabFig), saveDir, figureName, 'png', resolution);
					%saveFigure_s(h_figures((lI-1)+firstMetabFig), saveDir, figureName, 'png', 300);
				end
				% Save figure for summed up macromolecule contributions for 3T data
				if( strcmp(strBoField, '3T') )
					figureName	= strcat('fitSig_', 'SummedMM', strSaveName);
					saveFigure_s(h_figures(indFigSummedMM), saveDir, figureName, 'fig', resolution);
					saveFigure_s(h_figures(indFigSummedMM), saveDir, figureName, 'eps', resolution);
					saveFigure_s(h_figures(indFigSummedMM), saveDir, figureName, 'png', resolution);
				end
			case 'Default'
				fprintf('\n\nMetabolite figures were not saved!\n\n\n');
			case 'No'
				fprintf('\n\nMetabolite figures were not saved!\n\n\n');
		end
	end		% End of if( strcmp(strShowMetaboliteFits, 'YES') )
end		% End of if( strcmp(strShowSpectrumAndFit, 'YES') )


%% Creatine for linewidth estimation, if selected
if( strcmp(strLinewidthFit, 'YES') )
    
    % Generate selected signals for plotting and peak fitting
    % Cr fit after subtracting baseline
    metabolite_Cr		= metabolite_Cr_withBG - background;
    % PCr fit after subtracting baseline
    metabolite_PCr		= metabolite_PCr_withBG - background;
	
    % Cr part of spectrum that is left after subtracting all other fits  
	% from the spectrum, except for Cr
	sum_Metabolites				= (sum(fitsMetabolites, 1)).';
	sum_Metabolites_wout_Cr		= sum_Metabolites - metabolite_Cr;
	spectrum_Cr					= raw_spectrum - background - sum_Metabolites_wout_Cr;
	if( strcmp(strAllLinewidthFitsWoutResiduals, 'YES') )
		spectrum_Cr				= spectrum_Cr - residuals;
	end
	
    % tCr (= Cr+PCr) part of spectrum that is left after subtracting all other fits 
	% from the spectrum, except for tCr
    sum_Metabolites_wout_tCr	= sum_Metabolites - metabolite_Cr - metabolite_PCr;
	spectrum_tCr				= raw_spectrum - background - sum_Metabolites_wout_tCr;
	if( strcmp(strAllLinewidthFitsWoutResiduals, 'YES') )
		spectrum_tCr				= spectrum_tCr - residuals;
	end
    
    % Plot data for linewidth fit, if desired
    if( strcmp(strShowLinewidthData, 'YES') )
        
        % Plot Cr metabolite fit
        delta_plot      = 0.1;
        if (min_spectrum < 0)
            % Part of MRS spectrum negative
            peaksCrPlotRange	= [(2*abs(min(metabolite_Cr))*sign(min(metabolite_Cr))) (1+delta_plot)*max(metabolite_Cr)];
            peaks_tCrPlotRange	= [(2*abs(min(spectrum_tCr))*sign(min(spectrum_tCr))) (1+delta_plot)*max(spectrum_tCr)];
        else
            % MRS spectrum only positive
            peaksCrPlotRange	= [-delta_plot*max(metabolite_Cr) (1+delta_plot)*max(metabolite_Cr)];
            peaks_tCrPlotRange	= [-delta_plot*max(spectrum_tCr) (1+delta_plot)*max(spectrum_tCr)];
        end
        ranges_fig_Cr		= [ppmPlotRange peaksCrPlotRange];
        indFigs				= indFigs + 1;
        h_figures(indFigs)	= figure;
        hplot				= plot(ppm, metabolite_Cr, 'k', 'LineWidth', 1.5);
        set(gca,'Xdir','reverse', 'XTick', [0.5:0.5:5], 'XTickLabel', {'', '1', '', '2', '', '3', '', '4', '', '5'}, ...
            'FontSize', font, 'LineWidth', 2, 'FontWeight', 'bold', 'XColor', [0 0 0], ...
            'YColor', [0 0 0], 'TickDir', 'out', 'TickLength', [0.02, 0.025], 'Box', 'off')%, ...
        %'YTickLabelMode', 'Manual', 'YTick', [0], 'YTickLabel', {'0'}, 'YGrid', 'on');
        %'YTickLabelMode', 'Manual', 'YTick', []);
        axis (ranges_fig_Cr);
        title('Cr');
        xlabel('ppm', 'FontSize', (font+4), 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
            'Units', 'normalized', 'Position', ppmTextPos);
        
        % Plot Cr part of spectrum that is left after subtracting all other fits from the
        % spectrum, except for Cr
        indFigs				= indFigs + 1;
        h_figures(indFigs)	= figure;
        hplot				= plot(ppm, spectrum_Cr, 'k', 'LineWidth', 1.5);
        set(gca,'Xdir','reverse', 'XTick', [0.5:0.5:5], 'XTickLabel', {'', '1', '', '2', '', '3', '', '4', '', '5'}, ...
            'FontSize', font, 'LineWidth', 2, 'FontWeight', 'bold', 'XColor', [0 0 0], ...
            'YColor', [0 0 0], 'TickDir', 'out', 'TickLength', [0.02, 0.025], 'Box', 'off')%, ...
        %'YTickLabelMode', 'Manual', 'YTick', [0], 'YTickLabel', {'0'}, 'YGrid', 'on');
        axis (ranges_fig_Cr);
        title('Cr Part of Spectrum');
        xlabel('ppm', 'FontSize', (font+4), 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
            'Units', 'normalized', 'Position', ppmTextPos);
                
        % Plot tCr (= Cr+PCr) part of spectrum that is left after subtracting all other fits
        % from the spectrum, except for tCr
        ranges_fig_tCr		= [ppmPlotRange peaks_tCrPlotRange];
        indFigs				= indFigs + 1;
        h_figures(indFigs)	= figure;
        hplot				= plot(ppm, spectrum_tCr, 'k', 'LineWidth', 1.5);
        set(gca,'Xdir','reverse', 'XTick', [0.5:0.5:5], 'XTickLabel', {'', '1', '', '2', '', '3', '', '4', '', '5'}, ...
            'FontSize', font, 'LineWidth', 2, 'FontWeight', 'bold', 'XColor', [0 0 0], ...
            'YColor', [0 0 0], 'TickDir', 'out', 'TickLength', [0.02, 0.025], 'Box', 'off')%, ...
        %'YTickLabelMode', 'Manual', 'YTick', [0], 'YTickLabel', {'0'}, 'YGrid', 'on');
        axis (ranges_fig_tCr);
        title('tCr Part of Spectrum');
        xlabel('ppm', 'FontSize', (font+4), 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
            'Units', 'normalized', 'Position', ppmTextPos);
	end		% End of if( strcmp(strShowLinewidthData, 'YES') )	
    
    % Fit peaks of the different versions of the fits for Cr and tCr
	% Methyl peaks around 3.03 ppm are selected
	peakFit_position	= 3.03;
	peakFit_range		= 1.2;
	noPeaks				= 1;
	peakFit_shape		= 2;	% 2 = unconstrained Lorentzian
	peakFit_extra		= 0;	% Specifies the value of 'extra', used in the Pearson, 
								% exponentially-broadened Gaussian, Gaussian/Lorentzian 
								% blend, bifurcated Gaussian, and Breit-Wigner-Fano shapes
								% to fine-tune the peak shape
	NumTrials			= 10;	% Restarts the fitting process "NumTrials" times and 
								% selects the best one (with lowest fitting error). 
								% NumTrials can be any positive integer (default is 1). 
								% In may cases, NumTrials=1 will be sufficient, but if 
								% that does not give consistent results, increase 
								% NumTrials until the result are stable
	
	% Peak fit for Cr metabolite fit from LCModel
	indFigs				= indFigs + 1;
    h_figures(indFigs)	= figure;
    fitMatrix_metabolite_Cr					= [ppm.'; metabolite_Cr.'];
    [FitResults_Cr,FitError_Cr]				= peakfit(fitMatrix_metabolite_Cr, peakFit_position, peakFit_range, noPeaks, peakFit_shape, 0, 10);
	peakWidth_metabolite_Cr_Hz				= FitResults_Cr(4) * ppm_Hz;
    
	% Peak fit for Cr part of spectrum
	indFigs				= indFigs + 1;
    h_figures(indFigs)	= figure;
	fitMatrix_spectrum_Cr					= [ppm.'; spectrum_Cr.'];
    [FitResults_specCr,FitError_specCr]		= peakfit(fitMatrix_spectrum_Cr, peakFit_position, peakFit_range, noPeaks, peakFit_shape, 0, 10);
	peakWidth_spectrum_Cr_Hz				= FitResults_specCr(4) * ppm_Hz;
	
	% Peak fit for tCr part of spectrum
	indFigs				= indFigs + 1;
    h_figures(indFigs)	= figure;
	fitMatrix_spectrum_tCr					= [ppm.'; spectrum_tCr.'];
    [FitResults_spectCr,FitError_spectCr]	= peakfit(fitMatrix_spectrum_tCr, peakFit_position, peakFit_range, noPeaks, peakFit_shape, 0, 10);
	peakWidth_spectrum_tCr_Hz				= FitResults_spectCr(4) * ppm_Hz;
	
	% Close plots of peak fits and adjust index for figures accordingly, if selected
	% (plotting of these fits automatically accurs in external routine peakfit.m)
	if( strcmp(strCloseLinewidthFits, 'YES') )
		close(h_figures(indFigs));
		close(h_figures(indFigs-1));
		close(h_figures(indFigs-2));
		indFigs		= indFigs - 3;
	end		% End of if( strcmp(strCloseLinewidthFits, 'YES') )
    
    % Display results from fits
	fprintf('\n\n\n');
	disp('Peak fitting results are for metabolite_Cr, spectrum_Cr, and spectrum_tCr:');
	disp('      Peak number  Position     Height      Width       Peak area	 Width_Hz  RMS_Err%	 R2');
	disp([FitResults_Cr peakWidth_metabolite_Cr_Hz FitError_Cr]);
	disp([FitResults_specCr peakWidth_spectrum_Cr_Hz FitError_specCr]);
    disp([FitResults_spectCr peakWidth_spectrum_tCr_Hz FitError_spectCr]);
	    
    % Write results from linewidth fits to file, if desired
	% File is saved in same directory, from which the data is read 
    if( strcmp(strSaveLinewidthFits, 'YES') )
		% Create full absolute output filename
		%fullOutFileName		= strcat(outDir, outFileName, out_ext);
		
		% Open file for appending output
		[out_fid, out_message]		= fopen(fullOutFileName, 'a');
		if out_fid == -1
			error('%s: Invalid file identifier out_fid = %d with error message = %s!\n\n', sFunctionName, out_fid, out_message);
		end
		
		% If file is empty, add information about peak fit and output as first line,
		% if not empty, set file pointer to end of file for appending data
		strLinewidthInfo	= sprintf('Fit_Data \t\tPosition_ppm\tHeight\tWidth_ppm\tPeak_area\tWidth_Hz\tRMS_Err%%\tR2');
		if fseek(out_fid, 1, 'bof') == -1
			% Empty file
			nbytes	= fprintf(out_fid, 'Peak_fitting_parameters:\t Position= %.2f Range= %.1f noPeaks= %d PeakShape= %d extra= %d NumTrials= %d\n\n', ...
				peakFit_position, peakFit_range, noPeaks, peakFit_shape, peakFit_extra, NumTrials);
			%fprintf(out_fid, 'Case\t%s\n\n', caseName);
			nbytes	= fprintf(out_fid, 'File\t%s\n\n', file_name);
			%fprintf(out_fid, '%s\t %s\t %s\n\n', strLinewidthInfo, strLinewidthInfo, strLinewidthInfo);
			nbytes	= fprintf(out_fid, '%s\t\n', strLinewidthInfo);
		else
			% File non-empty
			if fseek(out_fid, 0, 'eof') == -1
				error('Error setting file pointer out_fid to end of file!');
			end
		end
		
		% Write output results and some information formatted to textfile
		nbytes	= fprintf(out_fid, 'metabolite_Cr \t\t%.3f\t\t%.2f\t%.4f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.4f\n', ...
			FitResults_Cr(2:end), peakWidth_metabolite_Cr_Hz, FitError_Cr);
		nbytes	= fprintf(out_fid, 'spectrum_Cr \t\t%.3f\t\t%.2f\t%.4f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.4f\n', ...
			FitResults_specCr(2:end), peakWidth_spectrum_Cr_Hz, FitError_specCr);
		nbytes	= fprintf(out_fid, 'spectrum_tCr \t\t%.3f\t\t%.2f\t%.4f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.4f\n', ...
			FitResults_spectCr(2:end), peakWidth_spectrum_tCr_Hz, FitError_spectCr);
			%fprintf(out_fid, 'spectrum_tCr \t\t%.3f %.2f %.4f %.2f %.2f %.2f %.4f\n', ...
			%[FitResults_spectCr(2:end) peakWidth_spectrum_tCr_Hz FitError_spectCr]);
		nbytes	= fprintf(out_fid, '\n\n\n');

		% Close output file
		if( fclose(out_fid) == -1 )
			error('%s: Error closing output file for linewidth fit!', sFunctionName);
		end
	end

% 	% Save metabolite and other signals to file
% 	answer = questdlg('Do you want to save the signals for linewidth fits?', ...
% 		'Saving of Signals for Linewidth Fits', 'Yes','Default','No','No');
% 	switch answer
% 		case 'Yes'
% 			fullOutFileName	= strcat(outDir, outFileName);
% 			save(fullOutFileName, 'ppm', 'metabolite_Cr', 'background', 'raw_spectrum');
% 		case 'Default'
% 			disp('Signals for linewidth fits were not saved!');
% 		case 'No'
% 			disp('Signals for linewidth fits were not saved!');
% 	end
end		% End of if( strcmp(strLinewidthFit, 'YES') )


%% Use fit for NAA to estimate signal for the signal-to-noise ratio (SNR), if selected
if( strcmp(strEstimateSignal, 'YES') )
	% Find index for fit for NAA in cell array of metabolites
	lIndSig		= -1;
	lI			= 1;
	bFound_NAA	= 0;
	while( lIndSig < 0 && lI <= noMetabolites )
		if( strcmp(acellMetabolites{lI}, 'NAA') )
			lIndSig		= lI;
			bFound_NAA	= 1;
		else
			lI	= lI + 1;
		end
	end % End of while( lIndSig < 0 && lI <= noMetabolites )

	% Plot fit for NAA, if selected and NAA found
	if( strcmp(strShowEstimateSignal, 'YES') && bFound_NAA == 1 )
		min_Metabolite	= min(fitsMetabolites(lIndSig, :));
		max_Metabolite	= max(fitsMetabolites(lIndSig, :));
		% Assume that 'max_Metabolite' > 0
		if( min_Metabolite  > 0.1*max_Metabolite )
			ranges_figMetab		= [ppmPlotRange 0 (1.1*max_Metabolite)];
		else
			ranges_figMetab		= [ppmPlotRange ...
				(min_Metabolite-0.1*max_Metabolite) (1.1*max_Metabolite)];
		end
		indFigs				= indFigs + 1;
		h_figures(indFigs)	= figure;
		hplot				= plot(ppm, fitsMetabolites(lIndSig, :), 'k', 'LineWidth', plotLineWidth);
		set(gca,'Xdir','reverse', 'XTick', [0.5:0.5:5], 'XTickLabel', {'', '1', '', '2', '', '3', '', '4', '', '5'}, ...
			'FontSize', font, 'LineWidth', 2, 'FontWeight', 'bold', 'XColor', [0 0 0], ...
			'YColor', [0 0 0], 'TickDir', 'out', 'TickLength', [0.02, 0.025], 'Box', 'off', ...
			'YGrid', 'on', 'XGrid', 'on');
		%'YTickLabelMode', 'Manual', 'YTick', [0], 'YTickLabel', {'0'}, 'YGrid', 'off');
		if( strcmp(strUseSameScale, 'YES') )
			axis(ranges_figSameScale);
		else
			axis(ranges_figMetab);
		end
		title(sprintf('Fit for %s\n', acellMetabolites{lIndSig}));
		xlabel('ppm', 'FontSize', (font+4), 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', ...
			'Units', 'normalized', 'Position', ppmTextPos);
	end		% End of if( strcmp(strShowEstimateSignal, 'YES') && bFound_NAA == 1 )

	% Estimate signal by determining peak height of fit for NAA and maximum of fit for
	% spectrum and maximum of raw süectrum (i.e. spectrum including residuals, which for a
	% perfect fit should only be noise)
	fprintf('\n\n\n');
	if bFound_NAA
		peakHeight_NAA_fit	= max(fitsMetabolites(lIndSig, :));
	else
		warning('\n%s: Fit for NAA not found in data!\n\n', sFunctionName);
		peakHeight_NAA_fit	= 0;
	end
	SI_max_spectrum_fit	= max(spectrum_fit);
	SI_max_raw_spectrum	= max_spectrum;

	% Display info about signal estimates
	fprintf('Signal estimates for \n\n');
	fprintf('File\t%s\n\n', file_name);
	fprintf('peakHeight_NAA_fit \t=\t %.1f\n', peakHeight_NAA_fit);
	fprintf('SI_max_spectrum_fit \t=\t %.1f\n', SI_max_spectrum_fit);
	fprintf('SI_max_raw_spectrum \t=\t %.1f\n', SI_max_raw_spectrum);
	fprintf('\n\n');

	% Save info about signal estimates to output file, if selected
	if( strcmp(strSaveEstimateSignal, 'YES') )
		% Open file for appending output
		[out_fid, out_message]		= fopen(fullOutFileName, 'a');
		if out_fid == -1
			error('%s: Invalid file identifier out_fid = %d with error message = %s!\n\n', sFunctionName, out_fid, out_message);
		end

		% Write signal estimates to output file
		nbytes	= fprintf(out_fid, 'Signal estimates for \n\n');
		nbytes	= fprintf(out_fid, 'File\t%s\n\n', file_name);
		nbytes	= fprintf(out_fid, 'peakHeight_NAA_fit  \t=\t %.1f\n', peakHeight_NAA_fit);
		nbytes	= fprintf(out_fid, 'SI_max_spectrum_fit \t=\t %.1f\n', SI_max_spectrum_fit);
		nbytes	= fprintf(out_fid, 'SI_max_raw_spectrum \t=\t %.1f\n\n', SI_max_raw_spectrum);

		% Close output file
		if( fclose(out_fid) == -1 )
			error('%s: Error closing output file for signal estimates!', sFunctionName);
		end
	end		% End of if( strcmp(strSaveEstimateSignal, 'YES') )
end		% End of if( strcmp(strEstimateSignal, 'YES') )


%% Move all figures to selected position on screen
% Position 'center" is not favourable for double monitor setup
for lI=1 : 1 : indFigs
	%movegui(h_figures(lI), 'center');
	movegui(h_figures(lI), 'west');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Old code
% Selection based on field strength
% % SPECIAL_SVS_Study_3T_7T and other studies
% switch( strBoField )
% 	case '3T'
% 		% 3T
% 		dir				= 'E:\Ralf\MRI\SPECIAL\Results\SPECIAL_SVS_Study_3T_7T\LC_Model_Analysis\3T_Trio_at_CHUV\';
% 		outDir	= 'E:\Ralf\MRI\SPECIAL\Results\SPECIAL_SVS_Study_3T_7T\Linewidth_Measurements\Linewidths_3T_Trio_at_CHUV\';
% 		
% 		% 3T analyzed between 0.2 - 4.2 ppm
% 		fileIndices		= [45; 147; 249; 351; 759];
% 		noValues		= 1010;
% 		stepIndices		= 102;
% 		ppmTextPos		= [0.94, -0.11, 0];
% 		ppmPlotDelta	= 0.05;
% 		
% 		% Metabolites in 3T basis (for range 0.2 - 5.2 ppm, metabolites are
% 		% different, e.g. no 'Ala') For range 0.2 - 4.2 ppm (varies for
% 		% each volunteer) For AureliaP, DanielleT, GuillaumeS, KatrinP,
% 		% NilsR
% 		acellMetabolites = {'Ala'; 'Asp'; 'PCho'; 'Cr'; 'PCr'; 'GABA'; 'Gln'; 'Glu'; ...
% 			'GSH'; 'Gly'; 'Ins'; 'Lac'; 'NAA'; 'Scyllo'; 'Tau'; 'Glc'; 'NAAG'; 'GPC'; ...
% 			'PE'; 'MM09'; 'MM20';'MM12';'MM14'; 'MM17';};
% 		noMacromolecules	= 5;
% 		
% 		switch(indCase)
% 			case 1
% 				% Aurelia Portmann
% 				fileName	= strcat(dir, 'Human_Siemens_3T_Trio_at_CHUV_AureliaP_09_04_08\coord');
% 				outFileName	= 'signal_Cr_3T_spectrum_LCModel_Fit_AureliaP';
% 				saveDir		= strcat('E:\Ralf\Papers\MR_Spectroscopy\SPECIAL_SVS_Study_3T_7T_05_08\Figures\Z_Support_Figures\', ...
% 					'Fits_and_Spectrum_fromLCModel_3T_SPECIAL_AureliaP_Sf_scan16_0_2_4_2ppm\');
% 				strSaveName	= '_3T_Sf_AureliaP_scan16_SPECIAL_TE6_NA128';
% 			case 2
% 				% Danielle Tendall
% 				fileName	= strcat(dir, 'Human_Siemens_3T_Trio_at_CHUV_DanielleT_24_04_08\coord');
% 				outFileName	= 'signal_Cr_3T_spectrum_LCModel_Fit_DanielleT';
% 				saveDir		= strcat('E:\Ralf\Papers\MR_Spectroscopy\SPECIAL_SVS_Study_3T_7T_05_08\Figures\Z_Support_Figures\', ...
% 					'Fits_and_Spectrum_fromLCModel_3T_SPECIAL_DanielleT_Sf_scan21_0_2_4_2ppm\');
% 				strSaveName	= '_3T_Sf_DanielleT_scan21_SPECIAL_TE6_NA128';
% 			case 3
% 				% Gregory Aubert
% 				fileName	= strcat(dir, 'Human_Siemens_3T_Trio_at_CHUV_GregoryA_21_05_08\coord');
% 				outFileName	= 'signal_Cr_3T_spectrum_LCModel_Fit_GregoryA';
% 				saveDir		= strcat('E:\Ralf\Papers\MR_Spectroscopy\SPECIAL_SVS_Study_3T_7T_05_08\Figures\Z_Support_Figures\', ...
% 					'Fits_and_Spectrum_fromLCModel_3T_SPECIAL_GregoryA_Sf_scan19_0_2_4_2ppm\');
% 				strSaveName	= '_3T_Sf_GregoryA_scan19_SPECIAL_TE6_NA128';
% 				
% 				% Metabolites in 3T basis For range 0.2 - 4.2 ppm For
% 				% GregoryA
% 				acellMetabolites = {'Ala'; 'Asp'; 'PCho'; 'Cr'; 'PCr'; 'GABA'; 'Gln'; ...
% 					'Glu'; 'GSH'; 'Gly'; 'Ins'; 'NAA'; 'Scyllo'; 'Tau'; 'Glc'; 'NAAG'; ...
% 					'GPC'; 'PE'; 'MM09'; 'MM20';'MM12';'MM14'; 'MM17';};
% 			case 4
% 				% Guillaume Saurais
% 				fileName	= strcat(dir, 'Human_Siemens_3T_Trio_at_CHUV_GuillaumeS_28_03_08\coord');
% 				outFileName	= 'signal_Cr_3T_spectrum_LCModel_Fit_GuillaumeS';
% 				saveDir		= strcat('E:\Ralf\Papers\MR_Spectroscopy\SPECIAL_SVS_Study_3T_7T_05_08\Figures\Z_Support_Figures\', ...
% 					'Fits_and_Spectrum_fromLCModel_3T_SPECIAL_GuillaumeS_Sf_scan21_0_2_4_2ppm\');
% 				strSaveName	= '_3T_Sf_GuillaumeS_scan21_SPECIAL_TE6_NA128';
% 			case 5
% 				% Katrin Petermann
% 				fileName	= strcat(dir, 'Human_Siemens_3T_Trio_at_CHUV_KatrinP_18_03_08\coord');
% 				outFileName	= 'signal_Cr_3T_spectrum_LCModel_Fit_KatrinP';
% 				saveDir		= strcat('E:\Ralf\Papers\MR_Spectroscopy\SPECIAL_SVS_Study_3T_7T_05_08\Figures\Z_Support_Figures\', ...
% 					'Fits_and_Spectrum_fromLCModel_3T_SPECIAL_KatrinP_Sf_scan25_0_2_4_2ppm\');
% 				strSaveName	= '_3T_Sf_KatrinP_scan25_SPECIAL_TE6_NA128';
% 			case 6
% 				% Nils Rettby
% 				fileName	= strcat(dir, 'Human_Siemens_3T_Trio_at_CHUV_NilsR_03_04_08\coord');
% 				outFileName	= 'signal_Cr_3T_spectrum_LCModel_Fit_NilsR';
% 				saveDir		= strcat('E:\Ralf\Papers\MR_Spectroscopy\SPECIAL_SVS_Study_3T_7T_05_08\Figures\Z_Support_Figures\', ...
% 					'Fits_and_Spectrum_fromLCModel_3T_SPECIAL_NilsR_Sf_scan17_0_2_4_2ppm\');
% 				strSaveName	= '_3T_Sf_NilsR_scan17_SPECIAL_TE6_NA128';
% 			case 7
% 				% MRS_Schizophrenia_Study - Mihaela Badica
% 				fileIndices		= [47; 149; 251; 353; 760];	
% 				
% 				acellMetabolites = {'Ala'; 'Asp'; 'PCho'; 'Cr'; 'PCr'; 'GABA'; 'Gln'; 'Glu'; ...
% 					'GSH'; 'Gly'; 'Ins'; 'Lac'; 'NAA'; 'Scyllo'; 'Tau'; 'Glc'; 'NAAG'; 'GPC'; ...
% 					'PE'; 'MM09'; 'MM20';'MM12';'MM14'; 'MM17';};
% 				
% 				dir			= 'E:\Ralf\MRI\MRS_Schizophrenia_Study\Results\LCModel_Analysis\LNAC_Study_Patients_with_DSI\';
% 				fileName	= strcat(dir, 'MRS_Schizophrenia_3T_at_CHUV_MihaelaB_21_08_09\coord_water_withOVS');
% 				outFileName	= 'signal_Cr_3T_spectrum_LCModel_Fit_MihaelaB';
% 				saveDir		= strcat('E:\Ralf\MRI\MRS_Schizophrenia_Study\Results\LNAC_Study_Patients_with_DSI\MRS_Schizophrenia_3T_at_CHUV_MihaelaB_21_08_09\', ...
% 					'Fits_and_Spectrum_fromLCModel_3T_SPECIAL_MihaelaB_TEM_scan15_0_2_4_2ppm\');
% 				strSaveName	= '_3T_TEM_MihaelaB_scan15_SPECIAL_TE6_NA148';
% 				
% 			otherwise
% 				error('No options for chosen volunteer/patient!');
% 		end
% 		
% 	case '7T'
		% 7T
        % From PTB Hoby box e81174
        %dir				= '/home/allgemein/projects/7T_2Ch_Sf_Coil_fMRS_Right_Left/results/LCModel_Analysis_CasesNewLamp/';
        % From PTB number cruncher e81151
        %dir				= '/echo2/allgemein/projects/7T_2Ch_Sf_Coil_fMRS_Right_Left/results/LCModel_Analysis_CasesNewLamp/';
        %dir				= '/echo2/mekle01/Ralf/PTB_Berlin/7T_MRS_MRI/Results/7T_MRS_24Ch_NovaMedical_Coil/LCModel_Analysis/';
		%dir				= 'E:\Ralf\MRI\SPECIAL\Results\SPECIAL_SVS_Study_3T_7T\LC_Model_Analysis\7T_HeadOnly_at_EPFL\';
		%outDir	= 'E:\Ralf\MRI\SPECIAL\Results\SPECIAL_SVS_Study_3T_7T\Linewidth_Measurements\Linewidths_7T_HeadOnly_at_EPFL\';
        %outDir	= dir';
				
% 		% 7T data analyzed between 0.2 - 4.2 ppm
%         % Some info about the coord file
%         % Assume that data values are stored in 10 columns
%         noColumns               = 10;
% 		fileIndexStartValues	= 41;
% 		noValues                = 1218
% 		stepIndices             = ceil(noValues/noColumns) + 1;
%         fileIndices             = fileIndexStartValues + [0; 1; 2; 3; 7;] .* stepIndices

		% % 7T analyzed between 0.2 - 5.2 ppm
		% fileIndices		= [40; 194; 348; 502; 1272];
		% noValues		= 1522;
		% stepIndices		= 154;
		% ppmTextPos		= [0.94, -0.11, 0];
		% ppmPlotDelta	= 0.05;

		% Metabolites in 7T basis
		% (for range 0.2 - 5.2 ppm, metabolites might be different, e.g. no 'Ala')
		% For range 0.2 - 4.2 ppm 
		% (varies for each volunteer)
		% For GregoryA, GuillaumeS, KatrinP
		% acellMetabolites = {'Mac'; 'Ala'; 'Asp'; 'PCho'; 'Cr'; 'PCr'; 'GABA'; 'Gln'; ...
		% 	'Glu'; 'GSH'; 'Gly'; 'Ins'; 'Lac'; 'NAA'; 'Scyllo'; 'Tau'; 'bHB'; ...
		%	'NAAG'; 'GPC'; 'PE';};	
		% noMacromolecules	= 1;
	
%         % 7T_fMRS_RL
%         % Obtain case name from filename that does NOT include the path
%         [path,caseName,ext] = fileparts(file_name); clear path; clear ext;
%         outFileName     = strcat('signal_Cr_7T_spectrum_LCModel_Fit_', caseName);
%         saveDir         = strcat('/echo2/allgemein/tmp/RalfM/Z_MRS_Processing/7T_fMRS_RL/');
%         strSaveName     = strcat('_7T_fMRS_RL_', caseName);
%         outDir    = saveDir;
        
        % Metabolites in 7T basis
        % For range 0.2 - 4.2 ppm
        % For 7T_fMRS_RL
%         noMacromolecules	= 1;
%         acellMetabolites = {'mac'; 'Asp'; 'PCho'; 'Cr'; 'PCr'; 'GABA'; ...
%             'Gln'; 'Glu'; 'GSH'; 'Ins'; 'Lac'; 'NAA'; 'Scyllo'; 'Tau'; ...
%             'Glc'; 'NAAG'; 'GPC'; 'bHB'; 'PE';};
%         sz_acellMetabolites = size(acellMetabolites)
        
        
        
%         switch(indCase)
%             case 1
%                 % 7T_fMRS_RL
%                 caseName        = '631__MDC-0021,V_1291__20131211_39';
%                 fileName        = strcat(dir, 'LCModel_ver631_150303_Cases_130619-140209_Left_2_140424Re_MRS_coord/', caseName, '.coord');
%                 outFileName     = strcat('signal_Cr_7T_spectrum_LCModel_Fit_', caseName);
%                 saveDir         = strcat('/echo2/allgemein/tmp/RalfM/Z_MRS_Processing/7T_fMRS_RL/');
%                 strSaveName     = strcat('_7T_fMRS_RL_', caseName);
%                 outDir    = saveDir;
%                 
%                 % Metabolites in 7T basis
%                 % For range 0.2 - 4.2 ppm
%                 % For 7T_fMRS_RL
%                 acellMetabolites = {'mac'; 'Asp'; 'PCho'; 'Cr'; 'PCr'; 'GABA'; ...
%                     'Gln'; 'Glu'; 'GSH'; 'Ins'; 'Lac'; 'NAA'; 'Scyllo'; 'Tau'; ...
%                     'Glc'; 'NAAG'; 'GPC'; 'bHB'; 'PE';};
%             case 2
%                 % V1081_Acc
%                 fileName	= strcat(dir, 'LCModel_V1081_Acc_09_06_11/V1081_Acc_coord_water_withOVS');
%                 outFileName	= 'signal_Cr_7T_spectrum_LCModel_Fit_V1081_Acc';
%                 saveDir		= strcat('/echo2/allgemein/tmp/RalfM/7T_MRS_MRI/Display_V1081_Acc_09_06_11/');
%                 strSaveName	= '_7T_24Ch_V1081_Acc_SPECIAL_TE9_NA64';
%                 
%                 % Metabolites in 7T basis
%                 % For range 0.2 - 4.2 ppm
%                 % For V1081_Acc
%                 acellMetabolites = {'Mac'; 'Asp'; 'PCho'; 'Cr'; 'PCr'; 'GABA'; ...
%                     'Gln'; 'Glu'; 'GSH'; 'Gly'; 'Ins'; 'Lac'; 'NAA'; 'Scyllo'; 'Tau'; ...
%                     'Asc'; 'bHB'; 'NAAG'; 'GPC'; 'PE';};
%             case 3
%                 % V1079_Occipital
%                 fileName	= strcat(dir, 'LCModel_V1079_Occipital_09_06_11/V1079_Occipital_coord_water_withOVS');
%                 outFileName	= 'signal_Cr_7T_spectrum_LCModel_Fit_V1079_Occipital';
%                 saveDir		= strcat('/echo2/allgemein/tmp/RalfM/7T_MRS_MRI/Display_V1079_Occipital_09_06_11/');
%                 strSaveName	= '_7T_24Ch_V1079_Occipital_SPECIAL_TE9_NA64';
%                 
%                 % Metabolites in 7T basis
%                 % For range 0.2 - 4.2 ppm
%                 % For V1079_Occipital
%                 acellMetabolites = {'Mac'; 'Ala'; 'Asp'; 'PCho'; 'Cr'; 'PCr'; 'GABA'; ...
%                     'Gln'; 'Glu'; 'GSH'; 'Gly'; 'Ins'; 'Lac'; 'NAA'; 'Scyllo'; 'Tau'; ...
%                     'Asc'; 'bHB'; 'Glc'; 'NAAG'; 'GPC'; 'PE';};
%             case 4
%                 % V1111_2_Occipital
%                 fileName	= strcat(dir, 'LCModel_V1111_2_Occipital_06_07_11/V1111_2_Occipital_coord_water_withOVS');
%                 outFileName	= 'signal_Cr_7T_spectrum_LCModel_Fit_V1111_2_Occipital';
%                 saveDir		= strcat('/echo2/allgemein/tmp/RalfM/7T_MRS_MRI/Display_V1111_2_Occipital_06_07_11/');
%                 strSaveName	= '_7T_24Ch_V1111_2_Occipital_SPECIAL_TE10_NA64';
%                 
%                 % Metabolites in 7T basis
%                 % For range 0.2 - 4.2 ppm
%                 % For V1111_2_Occipital
%                 acellMetabolites = {'Mac'; 'Asp'; 'PCho'; 'Cr'; 'PCr'; 'GABA'; ...
%                     'Gln'; 'Glu'; 'GSH'; 'Gly'; 'Ins'; 'Lac'; 'NAA'; 'Scyllo'; 'Tau'; ...
%                     'Asc'; 'bHB'; 'Glc'; 'NAAG'; 'PE';};
%                 
%             otherwise
%                 error('No options for chosen volunteer/patient!');
%         end
%         
% 	otherwise
% 		error('No options for chosen field strength!');
% end



% Previous cases
% 		switch(indCase)
% 			case 1
% 				% Aurelia Portmann
% 				fileName	= strcat(dir, 'Siemens_7T_HeadOnly_at_EPFL_AureliaP_09_04_08\coord');
% 				outFileName	= 'signal_Cr_7T_spectrum_LCModel_Fit_AureliaP';
% 				saveDir		= strcat('E:\Ralf\Papers\MR_Spectroscopy\SPECIAL_SVS_Study_3T_7T_05_08\Figures\Z_Support_Figures\', ...
% 					'Fits_and_Spectrum_fromLCModel_7T_SPECIAL_AureliaP_Sf_scan27_0_2_4_2ppm\');
% 				strSaveName	= '_7T_Sf_AureliaP_scan27_SPECIAL_TE6_NA64';
% 				
% 				% Metabolites in 7T basis
% 				% For range 0.2 - 4.2 ppm
% 				% For AureliaP
% 				acellMetabolites = {'Mac'; 'Ala'; 'Asp'; 'PCho'; 'Cr'; 'PCr'; 'GABA'; ...
% 					'Gln'; 'Glu'; 'GSH'; 'Gly'; 'Ins'; 'Lac'; 'NAA'; 'Scyllo'; 'Tau'; ...
% 					'Asc'; 'bHB'; 'NAAG'; 'GPC'; 'PE';};
% 			case 2
% 				% Danielle Tendall
% 				fileName	= strcat(dir, 'Siemens_7T_HeadOnly_at_EPFL_DanielleT_24_04_08\coord');
% 				outFileName	= 'signal_Cr_7T_spectrum_LCModel_Fit_DanielleT';
% 				saveDir		= strcat('E:\Ralf\Papers\MR_Spectroscopy\SPECIAL_SVS_Study_3T_7T_05_08\Figures\Z_Support_Figures\', ...
% 					'Fits_and_Spectrum_fromLCModel_7T_SPECIAL_DanielleT_Sf_scan27_0_2_4_2ppm\');
% 				strSaveName	= '_7T_Sf_DanielleT_scan27_SPECIAL_TE6_NA64';
% 				
% 				% Metabolites in 7T basis
% 				% For range 0.2 - 4.2 ppm
% 				% For DanielleT
% 				acellMetabolites = {'Mac'; 'Asp'; 'PCho'; 'Cr'; 'PCr'; 'GABA'; 'Gln'; ...
% 					'Glu'; 'GSH'; 'Ins'; 'Lac'; 'NAA'; 'Scyllo'; 'Tau'; 'Asc'; 'bHB'; ...
% 					'NAAG'; 'GPC'; 'PE';};
% 			case 3
% 				% Gregory Aubert
% 				fileName	= strcat(dir, 'Siemens_7T_HeadOnly_at_EPFL_GregoryA_21_05_08\coord');
% 				outFileName	= 'signal_Cr_7T_spectrum_LCModel_Fit_GregoryA';
% 				saveDir		= strcat('E:\Ralf\Papers\MR_Spectroscopy\SPECIAL_SVS_Study_3T_7T_05_08\Figures\Z_Support_Figures\', ...
% 					'Fits_and_Spectrum_fromLCModel_7T_SPECIAL_GregoryA_Sf_scan24_0_2_4_2ppm\');
% 				strSaveName	= '_7T_Sf_GregoryA_scan24_SPECIAL_TE6_NA64';
% 			case 4
% 				% Guillaume Saurais
% 				fileName	= strcat(dir, 'Siemens_7T_HeadOnly_at_EPFL_GuillaumeS_31_03_08\coord');
% 				outFileName	= 'signal_Cr_7T_spectrum_LCModel_Fit_GuillaumeS';
% 				saveDir		= strcat('E:\Ralf\Papers\MR_Spectroscopy\SPECIAL_SVS_Study_3T_7T_05_08\Figures\Z_Support_Figures\', ...
% 					'Fits_and_Spectrum_fromLCModel_7T_SPECIAL_GuillaumeS_Sf_scan25_0_2_4_2ppm\');
% 				strSaveName	= '_7T_Sf_GuillaumeS_scan25_SPECIAL_TE6_NA64';
% 			case 5
% 				% Katrin Petermann
% 				fileName	= strcat(dir, 'Siemens_7T_HeadOnly_at_EPFL_KatrinP_19_03_08\coord');
% 				outFileName	= 'signal_Cr_7T_spectrum_LCModel_Fit_KatrinP';
% 				saveDir		= strcat('E:\Ralf\Papers\MR_Spectroscopy\SPECIAL_SVS_Study_3T_7T_05_08\Figures\Z_Support_Figures\', ...
% 					'Fits_and_Spectrum_fromLCModel_7T_SPECIAL_KatrinP_Sf_scan25_0_2_4_2ppm\');
% 				strSaveName	= '_7T_Sf_KatrinP_scan25_SPECIAL_TE6_NA64';
% 			case 6
% 				% Nils Rettby
% 				fileName	= strcat(dir, 'Siemens_7T_HeadOnly_at_EPFL_NilsR_02_04_08\coord');
% 				outFileName	= 'signal_Cr_7T_spectrum_LCModel_Fit_NilsR';
% 				saveDir		= strcat('E:\Ralf\Papers\MR_Spectroscopy\SPECIAL_SVS_Study_3T_7T_05_08\Figures\Z_Support_Figures\', ...
% 					'Fits_and_Spectrum_fromLCModel_7T_SPECIAL_NilsR_Sf_scan22_0_2_4_2ppm\');
% 				strSaveName	= '_7T_Sf_NilsR_scan22_SPECIAL_TE6_NA64';
% 				
% 				% Metabolites in 7T basis
% 				% For range 0.2 - 4.2 ppm
% 				% For NilsR
% 				acellMetabolites = {'Mac'; 'Ala'; 'Asp'; 'PCho'; 'Cr'; 'PCr'; 'GABA'; ...
% 					'Gln'; 'Glu'; 'GSH'; 'Gly'; 'Ins'; 'Lac'; 'NAA'; 'Scyllo'; 'Tau'; ...
% 					'Asc'; 'bHB'; 'NAAG'; 'GPC'; 'PE';};
% 			case 7
% 				% Metabolite T2s for Nils Rettby
% 				fileIndices	= [46; 169; 292; 415; 1030];
% 				dir			= 'E:\Ralf\MRI\Metabolite_T2_Relaxation_Times\Results\LCModel_Analysis\Siemens_7T_HeadOnly_at_EPFL_NilsR_02_04_08\';
% 				fileName	= strcat(dir, 's5_isise80_coord');
% 				outFileName	= 'signal_Cr_7T_spectrum_LCModel_Fit_NilsR';
% 				saveDir		= strcat(dir, 'Displayed_LCModel_Output_T2s_NilsR_02_04_08\');
% 				strSaveName	= '_7T_Sf_NilsR_SPECIAL_TE80_NA64';
% 				
% 				% Metabolites in 7T basis
% 				% For range 0.2 - 4.2 ppm
% 				% For NilsR
% 				% For TE = 6 ms
% % 				acellMetabolites = {'Mac'; 'Asp'; 'PCho'; 'Cr_CH2'; 'Cr'; ...
% % 					'PCr_CH'; 'PCr'; 'GABA'; 'Gln'; 'Glu'; 'GSH'; 'Gly'; 'Ins'; 'Lac'; ...
% % 					'NAA'; 'NAA_m'; 'Scyllo'; 'Tau'; 'Asc'; 'bHB'; 'NAAG'; ...
% % 					'GPC'; 'PE';};
% 				
% 				% For TE = 20 ms
% % 				acellMetabolites = {'Mac'; 'Ala'; 'Asp'; 'Cr_CH2'; 'Cr'; ...
% % 					'PCr_CH'; 'PCr'; 'GABA'; 'Gln'; 'Glu'; 'GSH'; 'Gly'; 'Ins'; 'Lac'; ...
% % 					'NAA'; 'NAA_m'; 'Scyllo'; 'Tau'; 'bHB'; 'Glc'; 'NAAG'; ...
% % 					'GPC'; 'PE';};
% 
% 				% For TE = 80 ms
% 				acellMetabolites = {'Mac'; 'Asp'; 'PCho'; 'Cr_CH2'; 'Cr'; ...
% 					'PCr_CH'; 'PCr'; 'GABA'; 'Gln'; 'Glu'; 'GSH'; 'Gly'; 'Ins'; 'Lac'; ...
% 					'NAA'; 'NAA_m'; 'Scyllo'; 'Tau'; 'Asc'; 'bHB'; 'Glc'; 'NAAG'; 'PE';};
% 				
% 				% For TE = 110 ms
% % 				acellMetabolites = {'Mac'; 'Asp'; 'PCho'; 'Cr_CH2'; 'Cr'; ...
% % 					'PCr_CH'; 'PCr'; 'GABA'; 'Glu'; 'GSH'; 'Gly'; 'Ins'; 'Lac'; ...
% % 					'NAA'; 'NAA_m'; 'Scyllo'; 'Tau'; 'bHB'; 'Glc'; 'NAAG'; 'GPC';};
% 			otherwise
% 				error('No options for chosen volunteer/patient!');
% 		end

