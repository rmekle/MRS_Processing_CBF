% op_alignAverages_spectXcorr_s.m
% Jamie Near, McGill University 2014.
% Edits from Ralf Mekle (RM), Charite, 2025.
% 
% USAGE:
% [out,fs,phs] = op_alignAverages_spectXcorr_s(in,dataFlag,minppmSC,maxppmSC,refSC,filterFlagSC,plotFlagSC,X_nuclOffset);
% 
% DESCRIPTION:
% Perform frequency-domain spectral cross-correlation using a limited range of
% frequencies to correct frequency and phase drifts of MRS data.  As described in 
% Deelchand, D.K. et al.  Simultaneous frequency and phase corrections of single-shot MRS 
% data using cross-correlation.  Magn Reson Med 2025; 93(1):8-17.
% 
% INPUTS:
% in        = Input data structure
% dataFlag	= Flag to indicate whether original data was read as conjugate complex or not
% minppmSC	= Minimum of frequency range (ppm) used for spectral cross-correlation
% maxppmSC	= Maximum of frequency range (ppm) used for spectral cross-correlation
% refSC     = Character array to choose reference signal for spectral cross-correlation
%             - use 1st transient (='f') or mean all spectra ('m') as reference, 
%					default is 'f'
% filterFlagSC	= Flag whether to apply apodization (LB=5 and GF=0.12) to data before SC, 
%					default is off (=0)
% plotFlagSC	= Flag wehether to plot spectra and offsets, default is 0
% XnuclOffsetSC	= Offset in ppm for X-nucleus relative to water. e.g. = 4.65 for 1H
%
% OUTPUTS:
% out       = Output following alignment of averages.  
% fs        = Vector of frequency shifts (in Hz) used for alignment.
% phs       = Vector of phase shifts (in degrees) used for alignment.

function [out,fs,phs] = op_alignAverages_spectXcorr_s(in,dataFlag,minppmSC,maxppmSC,refSC,filterFlagSC,plotFlagSC,XnuclOffsetSC)

%% Set string for name of routine and display blank lines for enhanced output visibility
sFunctionName		= 'op_alignAverages_spectXcorr_s';


%% Check on MRS input data and input arguments
% Check whether MRS input data have already been coil combined
if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

% Check on (missing) input arguments and assign default values
maxNargin	= 8;
if nargin<maxNargin
	% Default value for 1H
	XnuclOffsetSC	= 4.65;
	if nargin<(maxNargin-1)
		plotFlagSC = 0;
		if nargin<(maxNargin-2)
			filterFlagSC = 0;
			if nargin<(maxNargin-3)
				refSC = 'f';
				if nargin<(maxNargin-4)
					maxppmSC = 3.6;
					if nargin<(maxNargin-5)
						minppmSC = 1.8;
						if nargin<(maxNargin-6)
							dataFlag = 'conj';
							if nargin<(maxNargin-7)
								error('%s: MRS input data missing. Aborting!', sFunctionName)
							end
						end
					end
				end
			end
		end
	end
end


%% Perform frequency and phase drift correction using spectral cross-correlation 
% Derive required parameters for spectral cross-correlation from MRS input data and input
% arguments:
% swSC				= Spectral width in Hz
% txfrq_ppmInHzSC	=  1 ppm in Hz from the transmitter frequency (txfrq) of the MRS data
% (minppmSC and maxppmSC are separate input parameters to preserve some  similarity with 
% routine op_alignAverages_fd(...) from FID-A)
swSC				= 1/in.dwelltime;
txfrq_ppmInHzSC		= in.txfrq/1e6;
chemicalRangeSC		= [minppmSC, maxppmSC];

% Determine whether MRS data contains subspectra
if in.dims.subSpecs==0
    B=1;
else
    B=in.sz(in.dims.subSpecs);
end

% Allocate arrays for frequency and phase shifts for all averages and all subspectra and 
% for extracted FIDs in time domain for all time points, all averages, and all subspectra
% and extract FIDs from input data structure
fs		= zeros(in.sz(in.dims.averages),B);
phs		= zeros(in.sz(in.dims.averages),B);
fids	= zeros(in.sz(in.dims.t),in.sz(in.dims.averages),B);
fids	= in.fids(:,:,:);

% Take conjugate complex of FID data, if data flag = 'conj'
% This is done to largely use original code and sign conventions in spectxCorr_(...),
% where to generate spectra from FIDs,
% the forward FFT is used (as is e.g. in Osprey), whereas
% the inverse FFT is used in the FID-A toolkit (probably since Siemens DICOM data (.IMA) 
% is read in using the conjugate Mode in io_loadspec_IMA_s(...) and Siemens raw data
% (twix) is read in using mapVBVD()... in io_loadspec_twix_s(...) that was originally 
% devised to read in k-space data of MR images)
if strcmp(dataFlag, 'conj')
	fids	= conj(fids);
end		% End of if strcmp(dataFlag, 'conj')

% For each subspectrum, extract FIDs in time-domain for all averages and
% perform spectral cross-correlation for extracted FIDs and
% extract frequency and phase shifts from output values
for m=1:1:B
	%fids(:,:,m)				= in.fids(:,:,m);
	[fids(:,:,m), outVal]	= spectXcorr_s(fids(:,:,m), chemicalRangeSC, refSC, filterFlagSC, plotFlagSC, swSC, txfrq_ppmInHzSC, XnuclOffsetSC);
	%fids(:,:,m)				= fidsCor;
	fs(:,m)					= outVal(:,1);
	phs(:,m)				= outVal(:,2);
end		% End of for m=1:1:B

% Take conjugate complex of FID data, if data flag = 'conj'
% to reverse the same operation performed prior to spectral cross-correlation, in order to
% be able to continue using operations from FID-A
if strcmp(dataFlag, 'conj')
	fids	= conj(fids);
end		% End of if strcmp(dataFlag, 'conj')


%% Fill in MRS data and complete all flag settings for output data structure
%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.freqcorrected=1;
   

end		% End of function
