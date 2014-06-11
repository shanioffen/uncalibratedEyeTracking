function [allTCcontrols,allTCcontrolsNeg] = calculateNeglect(toggle)
% function trackerbatchanalysis_so(toggle,listnum)
% toggle: 1 for eyetracker, 2 for most places, 3 for kirkwood

if nargin==0
    toggle = 3;
end

list = loadControlMatFiles;

[datafolder filefolder matfolder] = setvariables(toggle);
cd(datafolder)
nsubj = length(list);

vals2measure = {'A1x','B1x','C1x','R2x','A1y','B1y','C1y','R2y'};
nmeasures = length(vals2measure);
datamatrix = nan*ones(nsubj,nmeasures);
negdatamatrix = datamatrix;

argOutCount = 0;
argOutCountCal = 0;
goodTCcontrols = [];
goodTCcontrolsCal = [];
allTCcontrols = [];
allTCcontrolsCal = [];

redoDataMatrix=1;
if(redoDataMatrix)
    h=waitbar(0,'making data matrix');
    for isubj = 1:nsubj
        waitbar(isubj/nsubj,h);
        subj = list{isubj};
        matName = subj{1};
        s(isubj).matName = matName;
        s(isubj).subjid = stripext(stripext(matName));
        s(isubj).subjcode = s(isubj).subjid;
        s(isubj).edfname = [s(isubj).subjid '.edf'];
        s(isubj).subjeye = 'R';
        s(isubj).subjdelay = 0;
        
        matfilename = [matfolder matName];
        load(matfilename); % this loads t
        
        [newx newy negx negy] = raw2clean(t,s,isubj);
        
        x = (newx - nanmean(newx))/nanstd(newx);
        y = -(newy - nanmean(newy))/nanstd(newy);        
        xind = find(~isnan(x));
        yind = find(~isnan(y));        
        datamatrix = getStats(x,xind,y,yind,isubj,datamatrix);
       
        nnegx = (negx - nanmean(negx))/nanstd(negx);
        nnegy = -(negy - nanmean(negy))/nanstd(negy);        
        negxind = find(~isnan(nnegx));
        negyind = find(~isnan(nnegy));        
        negdatamatrix = getStats(nnegx,negxind,nnegy,negyind,isubj,negdatamatrix);
        
        allTCcontrols{isubj} = [x' y'];
        allTCcontrolsNeg{isubj} = [nnegx' nnegy'];
      
        clear t newx newy negx negy  xind yind negxind negyind   
        clear nnegx nnegy 
    end
    
    save allTCcontrols.mat allTCcontrols
    save allTCcontrolsNeg.mat allTCcontrolsNeg
 

    
    measname = ['dataMatrixControlsUncalibrated'];
    xlswrite([filefolder '/' measname],datamatrix);
    save([filefolder '/' measname '.mat'],'datamatrix')
    
    negmeasname = ['dataMatrixModelNeglect'];
    xlswrite([filefolder '/' negmeasname],negdatamatrix);
    save([filefolder '/' negmeasname '.mat'],'negdatamatrix')
  
    
else
    measname = ['sinAnalysisControlsNormalized4figs'];
    load([filefolder '/' measname '.mat'])
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions
%**************************************setvariables******************8***************
        function datamatrix = getStats(x,xind,y,yind,isubj,datamatrix);
        [fox1, gx1] = fit(xind',x(xind)','sin1');
        [foy1, gy1] = fit(yind',y(yind)','sin1');
        
        x1vals = coeffvalues(fox1);
        y1vals = coeffvalues(foy1);
        
        datamatrix(isubj,1) = x1vals(1);
        datamatrix(isubj,2) = x1vals(2);
        datamatrix(isubj,3) = x1vals(3);
        datamatrix(isubj,4) = gx1.adjrsquare;
        
        datamatrix(isubj,5) = y1vals(1);
        datamatrix(isubj,6) = y1vals(2);
        datamatrix(isubj,7) = y1vals(3);
        datamatrix(isubj,8) = gy1.adjrsquare;


%**************************************setvariables******************8***************
function [datafolder filefolder matfolder] = setvariables(toggle)

if nargin == 0
    toggle = 2;
end
switch toggle
    case 1
        fileprefix = '~/';
        
        datafolder = [fileprefix 'data/'];
        figurefolder = [fileprefix 'figures/'];
        filefolder = [fileprefix 'measures/'];
        matfolder = [fileprefix 'matFiles/'];
        
    case 2
        fileprefix = '~/Dropbox/sharedEyetrackingFolder/';
        datafolder = [fileprefix 'DATA/'];
        figurefolder = [fileprefix 'FIGURES/'];
        filefolder = [fileprefix 'MEASURES/'];
        matfolder = [fileprefix 'matFiles/'];
        
    case 3
        fileprefix = '/Volumes/ShaniWSQBackupHD/Dropbox/sharedEyetrackingFolder/';
        datafolder = [fileprefix 'DATA/'];
        figurefolder = [fileprefix 'FIGURES/'];
        filefolder = [fileprefix 'MEASURES/'];
        matfolder = [fileprefix 'matFiles/'];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction
%********************** makejulymeasures *********************
function [newx newy neglectx neglecty] = raw2clean(t,bigs,isubj)

f = t.fixations;
c = t.saccades;
b = t.blinks;
g = t.gaze;

xcal = g.x;
ycal = g.y;
x = g.rawx;
y = g.rawy;
p = g.pupil;
gt = g.time;

s=bigs(isubj);
if(isempty(s.subjdelay))
    s.subjdelay = 0;
end
droplength = 11+s.subjdelay;
mpersec = 500;
cycledirections = {'cw','cw','cw','cw','cw'};
cyclelengthsec =   [40;   40;  40;   40;   40];
cyclestarttimess =  [0;   40;  80;  120;  160];
cycleendtimess =   [40;   80; 120;  160;  200];

cyclestart = mpersec*cyclestarttimess + 1;
cycleend =    mpersec*cycleendtimess;
numcycles = length(cycledirections);
numcycletypes = 1; % only cw

nummeasurements = sum(cyclelengthsec)*mpersec;
% start after the 10s countdown
firstel = (mpersec)*droplength+1;
% account \ possibility that recording ended early or late:
lastel = min(size(x,2),nummeasurements+firstel-1);
keeprange = firstel:lastel;

if(length(keeprange)<nummeasurements)
    xgaze(length(keeprange)+1:nummeasurements) = nan;
end

gazetimes = gt(1,keeprange);
xgaze = x(1,keeprange);
ygaze = y(1,keeprange);
pgaze = p(1,keeprange);
xcalgaze = xcal(1,keeprange);
ycalgaze = ycal(1,keeprange);



% pad with nans if didnt watch whole thing
if(length(keeprange)<nummeasurements)
    xgaze(1,length(keeprange)+1:nummeasurements) = nan;
    ygaze(1,length(keeprange)+1:nummeasurements) = nan;
    pgaze(1,length(keeprange)+1:nummeasurements) = nan;
    xcalgaze(1,length(keeprange)+1:nummeasurements) = nan;
    ycalgaze(1,length(keeprange)+1:nummeasurements) = nan;
    gazetimes(1,length(keeprange)+1:nummeasurements)= nan;
end


newx = removeblinks(xgaze,b,gazetimes); newx = newx';
newy = removeblinks(ygaze,b,gazetimes); newy = newy';
newcalx = removeblinks(xcalgaze,b,gazetimes); newcalx = newcalx';
newcaly = removeblinks(ycalgaze,b,gazetimes); newcaly = newcaly';

neglectrange = [1:2*mpersec 29*mpersec+1:42*mpersec 69*mpersec+1:82*mpersec 109*mpersec:122*mpersec 149*mpersec+1:162*mpersec 189*mpersec+1:200*mpersec];

neglectx = newx-nanmean(newx);
neglectx(neglectrange) = NaN;
neglectx(neglectrange) = 2*nanmedian(neglectx)*randn(size(neglectrange));
neglecty = newy-nanmean(newy);
neglecty(neglectrange) = 2*nanmedian(neglecty)*randn(size(neglectrange));



% subfunction
%********************** makejulymeasures *********************
function newvec = removeblinks(invec,blinks,gazetimes);

bwin = 200;
ntpnts = length(invec);
t1 = gazetimes(1);
bstart = ((blinks.startTime - t1)/2)+1;
bend = ((blinks.endTime-t1)/2) + 1;
btimes = ((gazetimes-t1)/2)+1;

tmax = nanmax(btimes);
tmaxind = find(btimes == tmax);

nobmatrix = nan*ones(length(invec),length(bstart));
keepst = nan*ones(1,length(bstart));
keepet = nan*ones(1,length(bstart));

if(length(bstart)>150)
    bwin = 50;
end

for iblink = 1:length(bstart)
    startt = nanmax(1,bstart(iblink)-bwin);
    keepst(iblink) = startt;
    endt = nanmin(tmax,bend(iblink)+bwin);
    keepet(iblink) = endt;
    nobmatrix(:,iblink) = invec;
    nobmatrix(startt:endt,iblink) = nan;
end %  iblink

newvec = mean(nobmatrix,2); % take out all the blinks

