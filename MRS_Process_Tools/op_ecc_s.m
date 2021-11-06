% op_ecc_s.m
% Jamie Near, McGill University 2014.
% Edits from Ralf Mekle (RM), Charite, 2021.
% 
% USAGE:
% [out,outw]=op_ecc_s(in,inw);
% 
% DESCRIPTION:
% Perform an eddy current correction by estimating any non-linearity in the
% phase of the water unsuppressed data in the time domain and applying the
% appropriate correction to both the water suppressed and water unsuppressed 
% data.
% 
% INPUTS:
% in     = water suppressed input data in matlab structure format.
% inw    = water unsuppressed input data in matlab structure format.
%
% OUTPUTS:
% out    = Water suppressed output following eddy current correction  
% outw   = Water unsuppressed output following eddy current correction

function [out,outw]=op_ecc_s(in,inw);

% RM:
sFunctionName		= 'op_ecc_s';

% NOTE: Think about removing this condition, if ECC should also be performed on individual
% averages and/or for individual coil elements
if inw.dims.coils~=0 || inw.dims.averages~=0 || inw.dims.subSpecs~=0
    error('ERROR:  Must combine receivers, averages and subspecs prior to running ecc!! ABORTING!!');
end


% %save the phase as a vector of hard numbers
% % RM: Adjust plotting 
% % Show phase of water FID before ECC
% % (Note that plotting might have to be adjusted for multi-dimensional arrays of FIDs)
% inph=double(phase(inw.fids));
% figure;
% plot(inw.t,inph);

% RM:
% Extract FIDs from input structs and reorder dimensions to be able to use Dinesh
% Deelchand's ECC routine that uses FIDs stored in dimensions of the order 
% [time points, averages, coils] = [npX,nt,nbCoils] in Dinesh's code
fidw			= inw.fids;
fidm			= in.fids;

% Obtain info about size and dimensions of extracted FIDs
sz_fidw			= size(fidw);
sz_fidm			= size(fidm);
nDims_fidw		= ndims(fidw);
nDims_fidm		= ndims(fidm);

% Check if # of dimensions of extracted FIDs is equal to 2, i.e. only vectors of FID 
% values are left after coil combination and averaging 
% (CHECK on dimension of subspectra, if present!)
if nDims_fidm == 2 && nDims_fidw == 2
	% Only vectors of FID values left
	% Dinesh's routine for ECC and DC correction can be called without reshaping arrays 
	% of extracted FIDs
	% ECC and DC correction for metabolite FIDs
	fidm_ECC	= eddyCor2(fidw,fidm);
	% ECC and DC correction for accompanying water signal FIDs
	% (use water signal itself for this correction)
	fidw_ECC	= eddyCor2(fidw,fidw);	
else
	% NOT YET!
	error('%s: ECC not yet implemented for dimension of FIDs nDims_fidm = %d!', sFunctionName, nDims_fidm);
	% Reshape array of extracted FIDs to match order that is used in Dinesh's routine
	
	% Perform ECC and DC correction
	
	% Reshape array of corrected FIDs to match order that is used in FID-A structs
	
end


% RM:
% Insert corrected FIDs into corresponding output structs and calculate spectra
% MR spectra
out			= in;
out.fids	= fidm_ECC;
out.specs	= fftshift(ifft(out.fids,[],1),1);

% Water signals used for ECC
outw		= inw;
outw.fids	= fidw_ECC;
outw.specs=fftshift(ifft(outw.fids,[],1),1);

% % Show phase of water FID after ECC
% % (Note that plotting might have to be adjusted for multi-dimensional arrays of FIDs)
% figure;
% plot(outw.t,phase(outw.fids));





% Original code in op_ecc_s(...)
% %choose the part of the phase function that is most linear
% tmin=input('input min t value: ');
% tmax=input('input max t value: ');
% %figure;
% 
% %now fit a straight line to the linear part of the phase function
% p=polyfit(inw.t(inw.t>tmin & inw.t<tmax),inph(inw.t>tmin & inw.t<tmax)',1)
% 
% %now fit a spline to approximate a smooth version of the phase function
% pp=splinefit(inw.t,inph,150);
% 
% %now subtract the line from the spline to get the eddy current related
% %phase offset:
% ecphase=ppval(pp,inw.t)'-polyval(p,inw.t)';
% sz=size(in.fids);
% ecphase_rep=repmat(ecphase,[1 sz(2:end)]);
% size(ecphase_rep);
% figure;
% plot(inw.t,ecphase);
% 
% 
% %Now apply the eddy current correction to both the water suppressed and the
% %water unsuppressed data:
% out=in;
% out.fids=out.fids.*exp(1i*-ecphase_rep);
% out.specs=fftshift(ifft(out.fids,[],1),1);
% size(ecphase);
% % RM: Display first value of eddy current related phase offset
% %ecphase(1);
% fprintf('ecphase(1)\t = %f\n', ecphase(1));
% % RM: Call op_addphase(...) with default values for input variables 3 & 4 to be able to 
% % use input variable 5 to supress plotting of op_addphase(...)
% %out=op_addphase(out,180*ecphase_rep(1)/pi);
% out=op_addphase(out,180*ecphase_rep(1)/pi,0,4.65,1);
% 
% outw=inw;
% outw.fids=outw.fids.*exp(1i*-ecphase);
% outw.specs=fftshift(ifft(outw.fids,[],1),1);
% % RM: Call op_addphase(...) with default values for input variables 3 & 4 to be able to 
% % use input variable 5 to supress plotting of op_addphase(...)
% %outw=op_addphase(outw,180*ecphase(1)/pi);
% outw=op_addphase(outw,180*ecphase(1)/pi,0,4.65,1);
% figure;
% plot(outw.t,phase(outw.fids));