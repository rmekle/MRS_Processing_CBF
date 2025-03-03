%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% initParams_MRS_s.m
%
%% Function to initialize settings for processing magnetic resonance spectroscopy (MRS) data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USAGE
% [strDir_Out] = completeDirName_MRS_processed_s(strDir_Out_Base, fileExt, strVOI, dataType, signals, leftshift, avgBlockSize, rmbadav, noSD, strSpecReg, driftCorr, bECC, strProcessTool)
% 
% DESCRIPTION:
% Function to complete directory name for processed magnetic resonance spectroscopy (MRS)
% data depending on voxel location, data type, # of SDs, and other options used for 
% pre-processing of MR spectra or acquired macromolecules (MMs).
% 
% INPUTS:
% strDir_Out_Base    = String variable for the name of the base directory containing all
%						processed MR spectra for a specific study
%
% OUTPUTS:
% strDir_Out		= Completed directory name for processed MRS data
%
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2025;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [paramsMRS_struct] = initParams_MRS_s()


end		% End of function