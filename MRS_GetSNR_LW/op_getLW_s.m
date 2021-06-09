% op_getLW_s.m
% Jamie Near, McGill University 2014, Ralf Mekle (RM), Charite, 2021.
% 
% USAGE:
% [FWHM]=op_getLW_s(in,Refppmmin,Refppmmax,zpfactor,displaySwitch);
% 
% DESCRIPTION:
% Estimates the linewidth of a reference peak in the spectrum.  By default, 
% the reference peak is water, between 4.4 and 5.0 ppm.  Two methods are
% used to estimate the linewidth:  1.  FWHM is measured by simply taking the
% full width at half max of the reference peak.  2.  The FWHM is measured by
% fitting the reference peak to a lorentzian lineshape and determine the FWHM of the
% best fit.  The output FWHM is given by the average of these two measures.
% 
% INPUTS:
% in         = input spectrum in structure format.
% Refppmmin  = Min of frequency range (ppm) in which to search for reference peak.
%                  (Optional.  Default = 4.4 ppm);
% Refppmmax  = Max of frequency range to (ppm) in which search for reference peak
%                  (Optional.  Default = 5/3 ppm per Tesla B0);
% zpfactor   = zero-padding factor (used for method 1.)
%                  (Optional.  Default = 8);
% Change RM
% displaySwitch = Switch to turn on/off any plots or user interaction
%                  (Optional.  Default = 1);
%
% OUTPUTS:
% FWHM       = Estimated linewidth of the input spectrum (in Hz).


function [FWHM]=op_getLW_s(in,Refppmmin,Refppmmax,zpfactor,displaySwitch);

% Change RM
if nargin<5
    displaySwitch = 1;
    if nargin<4
        zpfactor=8;
        if nargin<3
            Refppmmax=5.0;
            if nargin<2
                Refppmmin=4.4;
            end
        end
    end
end

%in=io_readlcmraw(filestring,'dat');
in=op_zeropad(in,zpfactor);

%FIRST FIND FWHM USING TWO METHODS:

%METHOD 1:  ACTUALLY MEAUSURE FWHM OF WATER PEAK
Refwindow=in.specs(in.ppm>Refppmmin & in.ppm<Refppmmax);
ppmwindow=in.ppm(in.ppm>Refppmmin & in.ppm<Refppmmax);

maxRef_index=find(abs(real(Refwindow))==max(abs(real((Refwindow)))));
maxRef=real(Refwindow(maxRef_index));

% Change RM
% Only plot figure(s), if displaySwitch is on (= 1)
if displaySwitch == 1
% Add figure command to create figure for plotting
%plot(ppmwindow,abs(real(Refwindow)),'.');
figure, plot(ppmwindow,abs(real(Refwindow)),'.');
end

gtHalfMax=find(abs(real(Refwindow)) >= 0.5*abs(maxRef));

FWHM1=ppmwindow(gtHalfMax(1)) - ppmwindow(gtHalfMax(end));
FWHM1=FWHM1*(42.577*in.Bo);  %Assumes proton.


%METHOD 2:  FIT WATER PEAK TO DETERMINE FWHM PARAM
sat='n';
waterFreq=ppmwindow(maxRef_index);
while sat=='n'
    parsGuess=zeros(1,5);
    parsGuess(1)=maxRef; %AMPLITUDE
    parsGuess(2)=(5*in.Bo/3)/(42.577*in.Bo); %FWHM.  Assumes Proton.  LW = 5/3 Hz/T.
    parsGuess(3)=waterFreq; %FREQUENCY
    parsGuess(4)=0; %Baseline Offset
    parsGuess(5)=0; %Phase
    
    yGuess=op_lorentz(parsGuess,ppmwindow);
    parsFit=nlinfit(ppmwindow,real(Refwindow'),@op_lorentz,parsGuess);
    yFit=op_lorentz(parsFit,ppmwindow);
    
    % Change RM
    % Only plot figure(s), if displaySwitch is on (= 1)
    if displaySwitch == 1
        figure;
        plot(ppmwindow,real(Refwindow),'.',ppmwindow,real(yGuess),':',ppmwindow,yFit);
        legend('data','guess','fit');
        
        % Change RM
        % Only display diakog box for user input, if displaySwitch is on (= 1)
        sat=input('are you satisfied with fit? y/n [y] ','s');
    else
        % Set sat to 'y', so that no user interaction is required
        sat     = 'y';
    end
    if isempty(sat)
        sat='y';
    end
    if sat=='n';
        waterFreq=input('input new water frequency guess: ');
    end

end


FWHM2=abs(parsFit(2));
FWHM2=FWHM2*(42.577*in.Bo);  %Assumes Proton.

FWHM=mean([FWHM1 FWHM2]);  

disp(['The calculated linewidth is:  ' num2str(FWHM) ' Hz.' ]);
