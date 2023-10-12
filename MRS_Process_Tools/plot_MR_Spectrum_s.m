function [h_mrs, h_mrs2] = plot_MR_Spectrum_s(mrs, dataType, outDirString, outFileName, resolution, bSaveFigures)
%
% DESCRIPTION:
%  Function to plot an magnetic resonance (MR) spectrum stored in a data structure of the
%  FID-A toolkit for two different plotting ranges


%% Set string for name of routine and display blank lines for enhanced output visibility
sFunctionName		= 'plot_MR_Spectrum_s';

% Set figure parameters
% Set plotting resolution and properties of axes depending on data type
%resolution		= 600;
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
plot(mrs.ppm,real(mrs.specs),'linewidth',1.5);xlim(xLimValues1);
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
if bSaveFigures
	saveFigure_s(h_mrs, outDirString, figureName_fig, 'fig', resolution);
	saveFigure_s(h_mrs, outDirString, figureName_png, 'png', resolution);
end

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
if bSaveFigures
	saveFigure_s(h_mrs2, outDirString, figureName_fig, 'fig', resolution);
	saveFigure_s(h_mrs2, outDirString, figureName_png, 'png', resolution);
end

end