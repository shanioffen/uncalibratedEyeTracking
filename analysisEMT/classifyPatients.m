function classifyPatients


toggle = 3;

[datafolder filefolder matfolder] = setvariables(toggle);
cd(datafolder)


vals2measure = {'A1x','B1x','C1x','R2x','A1y','B1y','C1y','R2y'};
nmeasures = length(vals2measure);


patientDM = mkDataMatrixPatients;
controlDM = mkDataMatrixControls;

patientRvalsX = sqrt(patientDM(:,4));
patientRvalsY = sqrt(patientDM(:,8));

controlRvalsX = sqrt(controlDM(:,4));
controlRvalsY = sqrt(controlDM(:,8));

patientRvals = [patientRvalsX;patientRvalsY];
controlRvals = [controlRvalsX;controlRvalsY];

figure,
subplot(221), hist(patientRvals), title('Patient R vals');
subplot(222), hist(controlRvals), title('Control R vals');

controlRvals(length(controlRvals)+1:length(patientRvals)) = NaN;

[p,table] = kruskalwallis([patientRvals controlRvals]);

patmean = nanmean(patientRvals)
patSTD = nanstd(patientRvals)
conmean = nanmean(controlRvals)
conSTD = nanstd(controlRvals)

patientFisherZ = (.5)*log((1+patientRvals)./(1-patientRvals));
controlFisherZ = (.5)*log((1+controlRvals)./(1-controlRvals));

figure, subplot(221),
hist(patientFisherZ), title('Patient Fisher Z vals')
subplot(222), hist(controlFisherZ), title('Control Fisher Z vals');

[p,h] = ttest2(patientFisherZ,controlFisherZ)



ZpatMean = nanmean(patientFisherZ)
ZpatSTD = nanstd(patientFisherZ)
ZconMean = nanmean(controlFisherZ)
ZconSTD = nanstd(controlFisherZ)

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction
%**************************************setvariables******************8***************
function [datafolder filefolder matfolder] = setvariables(toggle)
if nargin == 0
    toggle = 3;
end
switch toggle
    case 1
        fileprefix = '~/';
        
        
        
    case 2
        fileprefix = '~/Dropbox/eyetracking/';
        
        
    case 3
        fileprefix = '/Volumes/ShaniWSQBackupHD/Dropbox/eyetracking/';
        
end

datafolder = [fileprefix 'data/'];
figurefolder = [fileprefix 'figures/'];
filefolder = [fileprefix 'measures/'];
matfolder = [fileprefix 'matFiles/'];

function discard


LL = patMean-2*patSTD;
UL = patMean+2*patSTD;

correctRejectAlt = (.5)*length(find(controlFisherZ>=zTH)); % control labeled healthy
FAsAlt = (.5)*length(find(controlFisherZ<zTH));            % control labeled impaired

missesAlt = (.5)*length(find(patientFisherZ>=zTH));        % patient labeled healthy
hitsAlt = (.5)*length(find(patientFisherZ<zTH));           % patient labeled impaired

testSensitivityAlt = hitsAlt/(hitsAlt+missesAlt)
testSpecificityAlt = correctRejectAlt/(FAsAlt+correctRejectAlt)



zTH = 2;

correctReject = (.5)*length(find(controlFisherZ>=zTH)); % control labeled healthy
FAs = (.5)*length(find(controlFisherZ<zTH));            % control labeled impaired

misses = (.5)*length(find(patientFisherZ>=zTH));        % patient labeled healthy
hits = (.5)*length(find(patientFisherZ<zTH));           % patient labeled impaired

testSensitivity = hits/(hits+misses)
testSpecificity = correctReject/(FAs+correctReject)


keyboard

patMean = nanmean(patientFisherZ)
patSTD = nanstd(patientFisherZ)

LL = 0; %patMean-2*patSTD
UL = 2; %patMean+2*patSTD

1-length(find(and(patientFisherZ>=LL,patientFisherZ<=UL)))/length(patientFisherZ)

1-length(find(and(controlFisherZ>=LL,controlFisherZ<=UL)))/length(controlFisherZ)

