%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% copyMetabResults_IntoExcel_s.m
%
%% Function to copy results from metabolite quantification of MRS data into MS Excel 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USAGE
% [resultsTable_final] = copyMetabResults_IntoExcel_s(fullFileName_In,noLastRowsToDel,bDeleteCols,indColsToDel,bMoveCols,noColsToMove,newColsTable,strOutDir,fileExcel,strSheet,strRange,bWriteVariableNames,bUseTemplate,fullPath_Template);
% 
% DESCRIPTION:
% Function that reads in results from metabolite quantification from a (.csv) file into a
% table, adds additional elements (columns) provided as input to table, reorders and saves
% final table to Excel file.
% If (Excel) template is provided, output is written to formatted Excel file.
% 
% INPUTS:
% fullFileName_In	= String variable for the full filename (path) to the (.csv) file
%						that contains the original results from metabolite quantification
% noLastRowsToDel	= Integer, # of last rows to be deleted from original table of results
% bDeleteCols		= Boolean that specifies whether to delete columns from original table
%						of results
% indColsToDel		= Indices of columns to be deleted from original table of results
% bMoveCols			= Boolean that specifies whether to move columns during reordering of 
%						new table
% noColsToMove		= integer, # of colums to be moved during reordering of new table
% newColsTable		= Integer, # of columns to be added to original table of results
% strOutDir			= String variable for the name of the output directory, i.e. the 
%						directory, where all output files are saved to
% fileExcel			= string variable, filename for output Excel file
% strSheet			= String variable that specifies name of sheet in Excel file, 
%						to which final table of results is saved/copied to
% strRange			= String variable that specifies range of fields or starting field on
%						sheet in Excel file, 
%						to which final table of results is saved/copied to
% bWriteVariableNames = Boolean that specifies whether variable names for data in table
%						are also written to spreadsheet in Excel
% bUseTemplate		= Boolean that specifies whether a (Excel) template file is copied to 
%					 
% fullPath_Template	= (Optional) String variable for full path of (Excel) template file
%						Default is '', an empty string assuming no template is used.
% 
% OUTPUTS:
% resultsTable_final = Final table of results from metabolite qunatification
%						(initial table of results plus modifications)
%
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2022, 2023; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [resultsTable_final] = copyMetabResults_IntoExcel_s(fullFileName_In,noLastRowsToDel,bDeleteCols,indColsToDel,bMoveCols,noColsToMove,newColsTable,strOutDir,fileExcel,strSheet,strRange,bWriteVariableNames,bUseTemplate,fullPath_Template)

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'copyMetabResults_IntoExcel_s.m';
sMsg_newLines		= sprintf('\n\n');
sMsg_newLine		= sprintf('\n');
disp(sMsg_newLines);


%% Check on # of input arguments and assign default values to variables, if required
maxNargin	= 14;
if nargin<maxNargin
	fullPath_Template = '';
end

% Logical check on input arguments
if bUseTemplate && isempty(fullPath_Template)
	error('%s: ERROR: bUseTemplate = %d, but empty fullPath_Template = %s!', sFunctionName, bUseTemplate, fullPath_Template);
end


%% Save/copy results from .csv file also as formatted Excel file
% For that, read contents of .csv file into Matlab table and then write metabolite
% information from table to (formatted) .xlsx file
% Detect import options of input file and keep empty lines of file in table
opts_read						= detectImportOptions(fullFileName_In);
opts_read.EmptyLineRule			= 'read';
opts_read.VariableNamingRule	= 'preserve';
resultsTable					= readtable(fullFileName_In, opts_read);
%resultsTable_orig				= resultsTable;

% Delete selected # of last rows in table of results, e.g rows that contain a blank line, 
% mean and STD computed by Matlab, since this information should be recomputed in Excel
sz_Results					= size(resultsTable);
% resultsTable_withSD		= resultsTable;
resultsTable((sz_Results(1)-(noLastRowsToDel-1)):sz_Results(1), :)		= [];
sz_Results					= size(resultsTable);


% Delete specified columns from table of results, if selected
if bDeleteCols
	resultsTable(:, indColsToDel)	= [];
	sz_Results						= size(resultsTable);
end

% Init final table of results
%resultsTable_final	= resultsTable;

% Insert selected # of columns for additional parameters, e.g. water linewidth (LW_H2O),
% into table of results by concatenating new input table with table of results
% and use # of colums to be moved from original table of results for reordering of columns
%  within resulting table to obtain desired order of columns, if selected
if bMoveCols
	sz_newCols			= size(newColsTable);
	noNewCols			= sz_newCols(2);
	resultsTable_final	= [newColsTable resultsTable];
	sz_final			= size(resultsTable_final);
	resultsTable_final	= resultsTable_final(:, [[noNewCols+(1:noColsToMove)] [1:noNewCols] [(noColsToMove+noNewCols+1):sz_final(2)]]);
else
	resultsTable_final	= resultsTable;
end
sz_final			= size(resultsTable_final);


% Write table of results to (formatted) Excel file
% Generate full path for Excel file; if corresponding input variable is empty,use filename
% from.csv file (without file extension); otherwise use given input Excel filename
if isempty(fileExcel)
	[filepath_In,name_In,ext_In]	= fileparts(fullFileName_In);
	fullFileName_Out				= fullfile(strOutDir, [name_In, '.xlsx']);
else
	fullFileName_Out				= fullfile(strOutDir, fileExcel);
end

% To use desired formatting in Excel, copy formatted Excel template to make it
% formatted Excel file for results, if selected
if bUseTemplate
	[status, msg]	= copyfile(fullPath_Template, fullFileName_Out);
	if ~status
		error('%s: ERROR: Copying of template file failed! msg = %s', sFunctionName, msg);
	end
end

% Write table to selected sheet and range of newly copied/created Excel file; variable
% names are included or not based on user input
% Set 'AutoFitWidth' option to false to preserve column width of the Excel template file
% (Default setting for this option is true, but since only values, i.e. numbers, and no
% variable names are written to the template file, the column width is then adjusted to
% only fit the numbers)
%writetable(resultsTable_final, fullFileName_Out, 'Sheet', strSheet, 'Range', strRange, 'WriteVariableNames', bWriteVariableNames);
writetable(resultsTable_final, fullFileName_Out, 'Sheet', strSheet, 'Range', strRange, 'WriteVariableNames', bWriteVariableNames, 'AutoFitWidth', false);


