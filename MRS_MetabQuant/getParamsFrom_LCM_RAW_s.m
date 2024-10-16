%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% getParamsFrom_LCM_RAW_s.m
%
%% Function to extract values of selected parameter(s) from LCModel .RAW files 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USAGE
% [params] = getParamsFrom_LCM_RAW_s(strDirIn, strDataType, listParams);
% 
% DESCRIPTION:
% Function that reads in MRS data processed for LCModel stored in .RAW files and extracts
% values of selected parameter(s) into output variable(s).
% 
% INPUTS:
% strDirIn			= Character array variable for the name of the input directory, where 
%						all .RAW files are stored
% strDataType		= Type of LCModel raw file (character array):
%						'rda' - .raw file generated from Siemens RDA file
%						'dat' - .raw file generated by FID-A using io_writelcm.
%						'sim' - .raw file generated from FID-A simulated data.
%						'raw' - not sure about this one.
% listParams		= List of selected parameters, for which values are extracted; can be
%						specified as character array for a single parameter; as array of 
%						character arrays or cell array for multiple parameters
% 
% OUTPUTS:
% params			= Values of selected parameter(s) extracted from (processed) LCModel
%						.RAW files (saved using FID-A toolkit)
%
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2022; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [params] = getParamsFrom_LCM_RAW_s(strDirIn, strDataType, listParams)

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'getParamsFrom_LCM_RAW_s.m';
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


%% Determine all files saved in LCModel .RAW files in input directory
% (assuming that all data files are included in the same directory)
% (On Linux, file list in Matlab also includes the two directories "." and "..", which
% means that the actual # of files in the directory is (# of entries in list - 2;
% however, if dir is used to list specific files, e.g. using a file extension, these two 
% directories are not included in the resulting list; 
% NOTE that name of directory has to have a file separator, i.e. a slash, at the end for 
% the latter to work)
%cd(strDirIn);
structFileListing		= dir([strDirIn, '*.RAW']);
noEntriesListing		= length( structFileListing );


%% Prepare output variables to store information about parameters
% List of parameters is string
params		= [];
if ischar(listParams)
	% Assume that values for only one parameter is to be extracted
	noParams	= 1;
	params		= zeros(noEntriesListing, noParams);
	%params		= zeros(noEntriesListing, 1);
else
	% Values for multiple parameters are to be extracted
	error('%s: ERROR: Extraction of values for multiple parameters not yet implemented!', sFunctionName);
end


%% Extravt values for selected parameter(s) from .RAW files
% Load each .RAW file in directory, read in (processed) MRS data, and extract values for
% selected parameter(s)
indexStart		= 1;
indexStep		= 1;
for ind=indexStart : indexStep : noEntriesListing
	fullFile_LCM_RAW	= [strDirIn, structFileListing(ind).name];
	[out]		= io_readlcmraw(fullFile_LCM_RAW, strDataType);
	if noParams == 1
		switch listParams
			case 'hzppm'
				% Parameter hzppm is listed in the $SEQPAR header of the (processed)
				% LCModel .RAW file; the value of this parameter is stored as out.txfrq
				% using txfrq=hzpppm*1e6; out.txfrq=txfrq;
				% (check FID-A routine io_readlcmraw(...) for details)
				params(ind)		= out.txfrq/1e6;
				
			otherwise
				error('%s: ERROR: Unknown listParams = %s!', sFunctionName, listParams);
		end		% End of switch listParams
	else
		% Values for multiple parameters are to be extracted
		error('%s: ERROR: Extraction of values for multiple parameters not yet implemented!', sFunctionName);
	end		% End of if noParams == 1
end		% for ind=indexStart : indexStep : noEntriesListing


