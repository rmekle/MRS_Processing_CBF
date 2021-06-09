function saveFigure_s(h_fig, dir, figureName, format, resolution);
%  
% saveFigure_s --	Save MATLAB figure in specified format at selected resolution.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Usage
%	saveFigure_s(h_figs, dir, filename, format)
%
%  Inputs
%	h_fig			handle of figure to be saved
%	dir				directory for saving figure
%	figureName		filename for figure saving	
%	format			format for figure saving
%	resolution		resolution (number) for figure saving depending on chosen file format
%						(e.g. in dpi for tif or in jpeg quality (1-100), see MATLAB help)
%		
%  Outputs
%
%  Description
%		Save MATLAB figure according to specified format at selected resolution.
%
%  See Also
%    mri_sim_main_s, saveMRIsimFigures_s, signal_analysis_s, updateDisplayTitle_s 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save figure into specified directory and with specified filename according to 
% selected format and resolution
% Generate complete filename for figure saving
filename		= strcat(dir, figureName);
switch(format)
	case('jpeg')
		% As .jpg file
		print(h_fig, '-djpeg', strcat('-r', num2str(resolution)), filename);
	case('fig')
		% As .fig file
		saveas(h_fig, filename, 'fig');
	case('tiffnocompression')
		% As .tif file without compression (figures only)
		%print(h_fig, '-dtiffnocompression', strcat('-r', num2str(resolution)), filename);
		print(h_fig, '-dtiffn', strcat('-r', num2str(resolution)), filename);
	case('eps')
		% As .eps black and white file (including a TIFF preview at 72 dpi resolution)
		print(h_fig, '-deps', '-tiff', strcat('-r', num2str(resolution)), filename);
	case('epsc')
		% As .eps color file (including a TIFF preview at 72 dpi resolution)
		print(h_fig, '-depsc', '-tiff', strcat('-r', num2str(resolution)), filename);
	case('emf')
		% As .emf color file 
		print(h_fig, '-dmeta', strcat('-r', num2str(resolution)), filename);
	case('png')
		% As .png color file 
		print(h_fig, '-dpng', strcat('-r', num2str(resolution)), filename);
		
	otherwise
		error('Specified file format not known (saveFigure_s)!');	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by Ralf Mekle, 2002; modified 2003
% rm197@columbia.edu






