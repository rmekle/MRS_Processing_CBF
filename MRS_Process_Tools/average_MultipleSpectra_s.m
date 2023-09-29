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
dirString_In_AddOn2		= '_MRS_jmrui';
dirString_In			= [dirString_In_Base, dirString_In_AddOn1, dirString_In_AddOn2, filesep];
filename_In				= '';

% Data format and # of input data files
dataType_MRS			= 'mrs_ref';		% 'mrs_w_ref';		'mrs_w';	% 'mrs_ref';	
dataFormat_In			= 'mrui';
structFileListing		= dir([dirString_In, '*.', dataFormat_In]);
noEntriesListing		= length( structFileListing );

% Output data directory and output data format
dirString_Out			= dirString_In;
dataFormat_Out			= 'mrui';


%% Parameters for processing
% Parameters for aligning multiple spectra with each other
bAlignSpectra			= 0;

% Parameters for spectral registration (aligning of averages/frequency and phase drift
% correction) performed in either frequency or time domain
%driftCorr_In			= 'y';		% 'y';		'n';
iterin_In				= 20;
aaDomain_In				= 'f';		% 'f';		't';
tmaxin_In				= 0.2;		% 0.2;		0.1;
bTmaxset_In				= 1;
medin_In				= 'y';		% 'y';	'n';	'a';	'ref';
% Set parameters for drift correction depending on type of data, i.e. whether MRS
% data is spectrum or water signal
% NOTE: Check whether aligning of averages in frequency domain works, if the MR
% spectrum is water signal itself; if not, simply align averages in time domain
switch dataType_MRS
	case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
		% MR spectrum is provided together without or with unsuppressed water
		% signal and/or with reference scans
		ppmmin_fix_In		= 1.6;		% 1.6;		1.8;
		ppmmaxarray_fix_In	= [3.5; 4.0; 5.5];
		%ppmmaxarray_fix_In	= [2.4,2.85,3.35,4.2,4.4,5.2];
	case {'water', 'water_ref'}
		% MR spectrum is water signal itself without or with reference scans
		ppmmin_fix_In		= 4.2;
		ppmmaxarray_fix_In	= [5.5 5.5 5.2];

	otherwise
		error('%s: Unknown MRS dataType_MRS = %s!', sFunctionName, dataType_MRS);
end		% End of switch dataType_MRS
%alignSS_In				= 2;		% For aligning subspectra (e.g. in SPECIAL)
strSpecReg				= 'SRs1';	% To distinguish settings for spectral registration

% Select indices of spectra to be averaged
indSpectra_All			= [1:noEntriesListing];
indSpectra_Select		= [1:noEntriesListing];
noSpectra_All			= length(indSpectra_All);
noSpectra_Select		= length(indSpectra_Select);
indexStart				= 1;
indexStep				= 1;


%% Determine output filename
filename_Out_Base		= dirString_In_AddOn1;
if noSpectra_Select == noSpectra_All
	% All spectra were averaged
	filename_Out_AddOn1		= sprintf('_All%d_avg', noSpectra_All);
else
	% Only selected spectra were averaged
	filename_Out_AddOn1		= sprintf('_Sel%d_avg', noSpectra_Select);
end		% End of if noSpectra_Select == noSpectra_All

% Indicate alignment of multiple spectra in output filename
if bAlignSpectra
	filename_Out_AddOn2		= sprintf('_%s_%s', strSpecReg, aaDomain_In);
else
	filename_Out_AddOn2		= '';
end
filename_Out			= [filename_Out_Base, filename_Out_AddOn1, filename_Out_AddOn2];


%% Load selected datasets (spectra) according to the selected data format
% Here, datasets are loaded into FID-A toolkit data structures for MR spectra
cellSpectra_Select	= cell(noSpectra_Select, 1);
switch dataFormat_In
	case 'dat'
		error('%s: ERROR: Averaging of multiple spectra not yet implemented for dataFormat_In = %s!', sFunctionName, dataFormat_In);
	case 'IMA'
		error('%s: ERROR: Averaging of multiple spectra not yet implemented for dataFormat_In = %s!', sFunctionName, dataFormat_In);
	case 'mrui'
		for ind=indexStart : indexStep : noSpectra_Select
		filename_In					= structFileListing(indSpectra_Select(ind)).name;
		cellSpectra_Select{ind}		= io_loadjmrui( fullfile(dirString_In, filename_In) );
		end		% End of ind=indexStart : indexStep : noSpectra_Selected
	otherwise
		error('%s: ERROR: Unknown data format = %s!', sFunctionName, dataFormat_In);
end		% End of switch dataFormat_In


%% Align multiple spectra in time domain or frequency domain, if selected
if bAlignSpectra
% 	[out,ph,frq]=op_alignAllScans(in,tmax,ref,mode);
% 	[out,ph,frq]=op_alignAllScans_fd(in,fmin,fmax,tmax,ref,mode);
	

end		% End of if bAlignSpectra


%% Average multiple spectra
% First, add all spectra using a FID-A routine and then divide results by the # of
% selected spectra
cellSpectra_Select_as	= cellSpectra_Select;
out						= cellSpectra_Select_as{1};
for j=2 : 1 : noSpectra_Select
	out		= op_addScans(out, cellSpectra_Select_as{j}, 0); 
end
% Divide FIDs of added spectra by # of selected spectra for averaging
out.fids	= out.fids./noSpectra_Select;

% Re-calculate spectra (specs) using fft
% Flags of FID_A data struct do not need to be updated for this last step
out.specs	= fftshift(ifft(out.fids,[],out.dims.t),out.dims.t);


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

