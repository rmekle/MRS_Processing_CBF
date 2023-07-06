%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% getParamsFrom_LCM_Output_s.m
%
%% Function to extract values of selected parameter(s) from LCModel output files 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USAGE
% [params] = getParamsFrom_LCM_Output_s(strDirIn, strFileType, listParams, listInds);
% 
% DESCRIPTION:
% Function that reads in a specific type of LCModel output files and extracts
% values of selected parameter(s) into output variable(s).
% 
% INPUTS:
% strDirIn			= Character array or string for the name of the input directory,
%						where all LCModel output files of a specific type are stored
% strFileType		= Type of LCModel output file (character array or string):
%						'coord'	- .coord file 
%						'coraw' - .coraw file
%						'csv'	- .csv file
%						'pdf'	- .pdf file
%						'print'	- .print file
%						'ps'	- .ps file
%						'table'	- .table file
% listParams		= List of selected parameters, for which values are extracted; can be
%						specified as character array for a single parameter; as array of 
%						character arrays or cell array for multiple parameters
% listInds			= List of indices, with which the occurrence of each of the selected
%						parameters can be specifed, from which value should be extracted,
%						is numerical array
% 
% OUTPUTS:
% params			= Values of selected parameter(s) extracted from LCModel output files
%
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2023; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [params] = getParamsFrom_LCM_Output_s(strDirIn, strFileType, listParams, listInds)

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'getParamsFrom_LCM_Output_s.m';
sMsg_newLines		= sprintf('\n\n');
sMsg_newLine		= sprintf('\n');
disp(sMsg_newLines);


%% Check on # of input arguments and assign default values to variables, if required
% maxNargin	= 9;
% if nargin<maxNargin
% 	fullPath_Template = '';
% end

% Logical check on input arguments
if isempty(strDirIn)
	error('%s: ERROR: Empty input directory = %s!', sFunctionName, strDirIn);
end


%% Determine all LCModel output files in input directrory
% (assuming that all data files are included in the same directory)
% (On Linux, file list in Matlab also includes the two directories "." and "..", which
% means that the actual # of files in the directory is (# of entries in list - 2;
% however, if dir is used to list specific files, e.g. using a file extension, these two 
% directories are not included in the resulting list; 
% NOTE that name of directory has to have a file separator, i.e. a slash, at the end for 
% the latter to work)
%cd(strDirIn);
fileExt					= ['.', strFileType];
structFileListing		= dir([strDirIn, ['*', fileExt]]);
noEntriesListing		= length( structFileListing );


%% Prepare output variables to store information about parameters
% List of parameters is string/character array or cell array
% Determine # of selected parameters from list of corresponding indices and allocate array
% for output parameter values
noParams	= length(listInds);
params		= zeros(noEntriesListing, noParams);

% params		= [];
% if ischar(listParams)
% 	% Assume that values for only one parameter is to be extracted
% 	noParams	= 1;
% 	params		= zeros(noEntriesListing, noParams);
% else
% 	% Values for multiple parameters are to be extracted
% 	error('%s: ERROR: Extraction of values for multiple parameters not yet implemented!', sFunctionName);
% end


%% Extravt values for selected parameter(s) from specific LCModel output files
% Read in each LCModel output file in directory, find parameter name(s), and extract 
% values for selected parameter(s)
indexStart		= 1;
indexStep		= 1;
for ind=indexStart : indexStep : noEntriesListing
	fullFile_LCM_Out	= [strDirIn, structFileListing(ind).name];
	% Read in LCModel output file depending on specific file type
	switch strFileType
		case {'coord', 'coraw', 'print', 'table'}
			% All file types are LCModel output textfiles
			encoding	= "";
			filetext	= fileread(fullFile_LCM_Out, Encoding=encoding);
	
		case 'csv'
			error('%s: ERROR: Reading in of LCModel output file type strFileType = %s not yet implemented!', sFunctionName, strFileType);
		case 'pdf'
			error('%s: ERROR: Reading in of LCModel output file type strFileType = %s not yet implemented!', sFunctionName, strFileType);
		case 'ps'
			error('%s: ERROR: Reading in of LCModel output file type strFileType = %s not yet implemented!', sFunctionName, strFileType);
		otherwise
			% Specified LCModel file type not known
			error('%s: ERROR: Unknown LCModel output file type strFileType = %s', sFunctionName, strFileType);
	end		% End of switch(strFileType)

	% Extract parameter values for all parameters
	for indP=1: 1 : noParams
		%if noParams == 1
		switch listParams(indP, :)
			case 'FCALIB'
				% Parameter FCALIB is included in LCModel output .print files
				% At first occurrence, FCALIB= 1.00000000 (always); hence search for a
				% second occurrence in .print files and if present, use this value
				params(ind, indP)		= 1;

				% Define the text to search for
				expr		= '[^\n]*FCALIB[^\n]*';

				%Find and return all lines that contain the text 'FCALIB'.
				matches		= regexp(filetext, expr, 'match');
				no_matches	= length(matches);

				% Display matching lines
				for indMatch=1 : 1 : no_matches
					disp(matches{indMatch});
				end		% End of for indMatch=1 : 1 : no_matches

				% Extract numerical value for 'FCALIB' from occurrence that is specified
				% by list of indices
				% Extract numerical content of string from matching character array and 
				% convert resulting string to numerical value for selected parameter
				strParam			= regexp(char(matches{listInds(indP)}),'[-+]?([0-9]*[.])?[0-9]+([eE][-+]?\d+)?','match');
				params(ind, indP)	= str2double(strParam);


			otherwise
				error('%s: ERROR: Unknown listParams = %s!', sFunctionName, listParams);
		end		% End of switch listParams
	%else
	%	% Values for multiple parameters are to be extracted
	%	error('%s: ERROR: Extraction of values for multiple parameters not yet implemented!', sFunctionName);
	%end		% End of if noParams == 1
	end		% End of for indP=1: 1 : noParams
end		% End of for ind=indexStart : indexStep : noEntriesListing


