% op_alignAverages_spectXcorr_s.m
% Jamie Near, McGill University 2014.
% Edits from Ralf Mekle (RM), Charite, 2025.
% 
% USAGE:
% [out,fs,phs] = op_alignAverages_spectXcorr_s(in,minppmSC,maxppmSC,refSC,filterFlagSC,plotFlagSC,X_nuclOffset);
% 
% DESCRIPTION:
% Perform frequency-domain spectral cross-correlation using a limited range of
% frequencies to correct frequency and phase drifts.  As described in Near
% et al.  Frequency and phase drift correction of magnetic resonance 
% spectroscopy data by spectral registration in the time domain. Magn Reson 
% Med 2015; 73(1):44-50.
% 
% INPUTS:
% in        = Input data structure.
% minppmSC	= Minimum of frequency range (ppm).
% maxppmSC	= Maximum of frequency range (ppm).
% refSC     = An externally provided reference spectrum that you would like
%             to align everything to (Required only if med = 'r').  
%
% OUTPUTS:
% out       = Output following alignment of averages.  
% fs        = Vector of frequency shifts (in Hz) used for alignment.
% phs       = Vector of phase shifts (in degrees) used for alignment.

function [out,fs,phs] = op_alignAverages_spectXcorr_s(in,minppmSC,maxppmSC,refSC,filterFlagSC,plotFlagSC,X_nuclOffset)

%% Set string for name of routine and display blank lines for enhanced output visibility
sFunctionName		= 'op_alignAverages_spectXcorr_s';


%% Check on MRS input data and input arguments
% Check whether MRS input data have already been coil combined
if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

% Check on (missing) input arguments and assign default values
maxNargin	= 7;
if nargin<maxNargin
	% Default value for 1H
	X_nuclOffset	= 4.65;
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
							error('%s: MRS input data missing. Aborting!', sFunctionName)
						end
					end
				end
			end
		end
	end
end


%% Perform frequency and phase drift correction using spectral cross-correlation 
% Determine whether MRS data contains subspectra
if in.dims.subSpecs==0
    B=1;
else
    B=in.sz(in.dims.subSpecs);
end

% Allocate arrays for frequency and phase shifts and extract FIDs in time domain for one
% subspectrum and for all time points and all averages
fs		= zeros(in.sz(in.dims.averages),B);
phs		= zeros(in.sz(in.dims.averages),B);
fids	= zeros(in.sz(in.dims.t),in.dims.averages,B);
for m=1:1:B
	% For each subspectrum, extract FIDs in time-domain for all averages and
	% perform spectral cross-correlation for extracted FIDs
    fids(:,:,m)		= in.fids(:,:,m);
	

end		% End of for m=1:1:B


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
