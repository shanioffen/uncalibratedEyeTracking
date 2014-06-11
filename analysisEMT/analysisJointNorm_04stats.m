function [zMat, hMat] = analysisForMethodsPaper_04stats

% ************************************************************************
onKirkwood = 1;

% ************************************************************************
if(onKirkwood)
    filePrefixKeep = '/Volumes/ShaniWSQBackupHD/Dropbox/';
else
    filePrefixKeep = '~/Dropbox/';
end
% ************************************************************************

% For each test time-course, calculate z-score for 4 measures based
% on control time-courses

%%%%%%%%%%%% load control distributions %%%%%%%%%%%%%%%%%%%%%
load([filePrefixKeep 'figures/figsForMayPaper/data/Distributions.mat']);
clear filePrefix
% this loads d, d{1} and d{2}, {1} is test, {2} is control
% the fields of d are
%         xVar: [1x5000 double]
%         yVar: [1x5000 double]
%         xR2: [1x5000 double]
%         yR2: [1x5000 double]
%         xR2TTA: [1x5000 double]
%         yR2TTA: [1x5000 double]

controlVar = d{2}.xyMeanVar;
controlInt = d{2}.xyIntegrity;


clear d


%%%%%%%%%%%% load test measure values %%%%%%%%%%%%%%%%%%%%%
load([filePrefixKeep 'figures/figsForMayPaper/data/TESTmeasures.mat']);
clear filePrefix
% this loads list, which is the list of test observers
% this loads m, which are measures for test observers
%
%        percentData: [1x77 double]
%           xmeanVar: [1x77 double]
%           ymeanVar: [1x77 double]
%          xymeanVar: [1x77 double]
%            maxxVar: [1x77 double]
%            maxyVar: [1x77 double]
%           maxxyVar: [1x77 double]
%      xTTAintegrity: [1x77 double]
%      yTTAintegrity: [1x77 double]
%     xyTTAintegrity: [1x77 double]
%             xR2TTA: [1x77 double]
%             yR2TTA: [1x77 double]
%            xyR2TTA: [1x77 double]
%                xR2: [1x77 double]
%                yR2: [1x77 double]
%               xyR2: [1x77 double]

nSubj = length(list);

muVec = [nanmean(controlVar) nanmean(controlInt)];
muMat = repmat(muVec,nSubj,1);

sigmaVec = [nanstd(controlVar) nanstd(controlInt)];
sigmaMat = repmat(sigmaVec,nSubj,1);

testMat = [m.xymeanVar' m.xyTTAintegrity'];

th = 10;
zMat = (testMat-muMat)./sigmaMat;
hMat = or(zMat>th,zMat<-th);


% notes:
% method 1: p val based on count of higher/lower
%         varXcontrol = repmat(controlVarX,nSubj,1);
%         varXtest = repmat(m.xVar',1,length(controlVarX));
%         difVarX = varXtest-varXcontrol;
%         binVarX = (difVarX>0);
%         countVarX = sum(binVarX,2)/size(binVarX,2);
%         hVarXm1 = or(countVarX>.95,countVarX<.05);
%
%
%
% method 2: z-score using control distribution mu and sigma
%         muVarX = nanmean(controlVarX);
%         sigmaVarX = nanstd(controlVarX);
%
%         pValsTestVarX = (m.xVar - muVarX)/sigmaVarX;
%         hVarXm2 = or(pValsTestVarX>2,pValsTestVarX<-2);
%
% they are the same, so I will use method 2
