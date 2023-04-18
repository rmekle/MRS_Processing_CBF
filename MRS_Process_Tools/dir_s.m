%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% dir_s.m
%
%% Function to list the files in folderPath ignoring the '.' and '..' paths
%  
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2023;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USAGE
% listing = dir_s(folderPath)
% 
% DESCRIPTION:
% Function to list the files in folderPath ignoring the '.' and '..' paths that are
% included in listing of a directory on Linux/Unix systems
% 
% INPUTS:
% folderPath    = Character vector or scalar string with file or folder name
% 
% OUTPUTS:
% listing		= Structure array with file attributes for file or folder contents
%					(see dir function in MATLAB)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function listing = dir_s(folderPath)
% DIR_S lists the files in folderPath ignoring the '.' and '..' paths

if nargin < 1; 
	folderPath	= '.'; 
elseif nargin == 1
	listing		= dir(folderPath);
	listing		= listing(~ismember({listing.name},{'.','..'}));
end
