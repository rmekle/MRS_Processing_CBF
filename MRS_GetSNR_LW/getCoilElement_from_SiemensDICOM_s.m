%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% getCoilELement_from_SiemensDICOM_s.m
%
%% Function to obtain number of a coil element from the DICOM header of a file generated 
%  on a Siemens MRI/NMR scanner 
%  using functions from the MRS processing toolkit FID-A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USAGE
% [coil_element, coil_string, info] = getCoilELement_from_SiemensDICOM_s(fullFilename_NMR, displaySwitch);
% 
% DESCRIPTION:
% Extract information about the coil elment, with which a DICOM image was acquired on a   
% Siemens MRI scanner. 
% 
% INPUTS:
% fullFilename_NMR	= String variable for the full filename (i.e. full path of filename) 
%					  of the DICOM file generated on a Siemens MRI/NMR scanner
% displaySwitch		= Switch to display outputs; 1 = ON, 0 = OFF
% 
% OUTPUTS:
% coil_element      = Number of coil element, with which DICOM input file was acquired on
%						a Siemens MRI/NMR scanner, if images of individual coil elements
%						were saved; empty, if DICOM file is a combined NMR data file
% coil_string		= String that contains information about RF coil or coil element, with
%					  which NMR data file was acquired:
%					  equal to # of coil element or
%					  equal to 'X', if input DICOM file was combined MR image or
%					  empty (''), if input DICOM file was combined MRS data
% info				= struct that contains the metadata from the compliant DICOM or DICOS 
%					  file specified as MRS input data file obtained using the Matlab
%					  routine dicominfo plus a .csa field that is a struct that contains 
%					  values from the 'SIEMENS CSA HEADER' private group in a DICOM file 
%					  produced by a Siemens MR scanner
%
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2021; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coil_element, coil_string, info] = getCoilElement_from_SiemensDICOM_s(fullFilename_NMR, displaySwitch)

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'getCoilElement_from_SiemensDICOM_s';
sMsg_newLines		= sprintf('\n\n');
disp(sMsg_newLines);


%% Check input arguments and read in DICOM header including Siemens csa field to obtain 
% information about coiul element, with which DICOM input file was acquired

% % Obtain different parts of filename
% [sPathStrNMR,nameNMR,extNMR] 	= fileparts(fullFilename_NMR);

% Check on required input arguments
if( isempty(fullFilename_NMR) )
    error('%s: Error: Empty DICOM input filername  = "%s!"\n\n', sFunctionName, fullFfilename_NMR);
end

% Init coil string and coil element
% (for the case that 'info.csa.ICE_Dims' is empty, e.g. for combined MRS data)
coil_string		= '';
coil_element	= [];

% Read in DICOM header of input file using specific routine form FID-A toolkit that in 
% turn uses the Matlab routine dicominfo to read in DICOM header and additionally creates  
% a .csa field that is a struct that contains values from the 'SIEMENS CSA HEADER' provate
% group in a DICOM file produced by a Siemens MR scanner
info			= SiemensCsaParse(fullFilename_NMR);

% Obtain information about coil element, with which Siemens input DICOM data was acquired,
% if corresponding string is not empty
if ~isempty(info.csa.ICE_Dims)
	indices 		= strfind(info.csa.ICE_Dims, '_');
	coil_string		= info.csa.ICE_Dims(1:indices(1)-1);
	coil_element 	= str2num( coil_string );
end

% Display results on screen, if selected
if displaySwitch == 1	
	fprintf('info.csa.ICE_Dims = %s\t\tcoil_string = %s\t\tcoil_element = %d\n', info.csa.ICE_Dims, coil_string, coil_element);
end



