%% with phase correction method
function fid = eddyCor2(fidw,fidm)
%
% function fid=eddyCor2(waterfid,metabfid,flag)
% DC, ECC and zero-order phase corrections
%
% Dinesh Deelchand 21 March 2008
%

% duplicate fidw to have same nt as fidm
repFlag = 0;
myFlagRun = 1;
if (size(fidw,2)==1)
    [~,nt] = size(fidm);
    fidw = repmat(fidw,[1 nt]);
    repFlag=1;
end


% DC correction
% set glitch to zero - 15 July 2014
glitchPts = 200;
fidw (end-glitchPts:end,:,:) = 0;

[npX,nt,nbCoils] =size(fidm);
myPointsUse = round(npX*20/100);
for jcal = 1:nt
    for ical = 1:nbCoils
        tempM = mean(fidm(length(fidm)-glitchPts-myPointsUse:length(fidm)-glitchPts,jcal,ical));
        fidm(:,jcal,ical) = fidm(:,jcal,ical) - tempM;
    end
end
nbCoils = size(fidw,2);
for ical = 1:nbCoils
    tempW = mean(fidw(length(fidw)-glitchPts-myPointsUse:length(fidw)-glitchPts,ical));
    fidw(:,ical) = fidw(:,ical) - tempW;
end

% ECC
[~,nt,nbCoils]=size(fidm);
fid = complex(zeros(size(fidm)));
for ical=1:nbCoils
    for jcal=1:nt
        if nbCoils==1
            % 2D matrix
            tempM =  fidm(:,jcal,ical);
            tempW = fidw(:,jcal);
            phiM = angle(tempM) - angle(tempM(1));
            % extrapolate phase of wref1
            if (myFlagRun==1)
                myTempEP = phzEP(tempW);
                phiW = myTempEP - myTempEP(1);
                %run only once if repmat used above
                if repFlag==1
                    myFlagRun=0;
                end
            end
            fid(:,jcal,ical) = abs(tempM).*(exp(1i*(phiM-phiW)));
        else
            % 3D matrix
            tempM =  fidm(:,jcal,ical);
            tempW = fidw(:,ical);
            phiM = angle(tempM) - angle(tempM(1));
            myTempEP = phzEP(tempW);
            phiW = myTempEP - myTempEP(1);
            fid(:,jcal,ical) = abs(tempM).*(exp(1i*(phiM-phiW)));
        end
    end
end

% Redo DC correction
fid = fid - mean(fid(length(fid)-myPointsUse:length(fid)));

return




function phzCor = phzEP(fidw)
% function to extrapolate phase
% based on where the phase deviated
% Output: phase info in radians
% Dinesh Deelchand, 17 July 2014

% unwrap phase info
myphz = unwrap(angle(fidw));
NP=length(myphz);

% determine if slope is straight or curved
% This is the simplest way to do this
Slope = diff(myphz);
myTol = max(abs(Slope(2:200)));  %ignore 1st point
indx = find(abs(Slope)>myTol);

% index
if (length(indx)>length(myphz)/2) %noisy data
    myindx = indx(2);
else
    myindx = indx(2)-100;   % ignore 1st point
    % ignore extra 100points to make sure it's not on the curvature part
end

% Use part of phase which is linear
myphzUse = myphz(1:myindx);
xtempUse = (1:length(myphzUse))';
xReq = (xtempUse(end)+1:length(myphz))';


% extrapolate by fitting linear part
p = polyfit(xtempUse,myphzUse,1);
phzNew = polyval(p,xReq);

% merge original linear part with extrapolated term
final_phz = [myphzUse; phzNew];

phzCor = final_phz;
return