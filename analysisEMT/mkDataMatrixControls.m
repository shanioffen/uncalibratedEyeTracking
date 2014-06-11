function datamatrix = mkDataMatrixControls
% function trackerbatchanalysis_so(toggle,listnum)
% toggle: 1 for eyetracker, 2 for most places, 3 for kirkwood

toggle = 3;


list = loadControlMatFiles;

[datafolder filefolder matfolder] = setvariables(toggle);
cd(datafolder)
nsubj = length(list);

vals2measure = {'A1x','B1x','C1x','R2x','A1y','B1y','C1y','R2y'};
nmeasures = length(vals2measure);
datamatrix = nan*ones(nsubj,nmeasures);

redoDataMatrix=0;
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
        
        [newx newy] = raw2clean(t,s,isubj);
        
        x = (newx - nanmean(newx))/nanstd(newx);
        y = -(newy - nanmean(newy))/nanstd(newy);
        
 
        xind = find(~isnan(x));
        yind = find(~isnan(y));
        
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
        
        
        clear t newx newy xind yind fftx ffty fox1 fox2 foy1 foy2 gx1 gx2 gy1 gy2 x y        
    end
   
    
    measname = ['dataMatrixControls'];
    xlswrite([filefolder '/' measname],datamatrix);
    save([matfolder '/' measname '.mat'],'datamatrix')
else
    measname = ['dataMatrixControls'];
    load([matfolder '/' measname '.mat'])
end


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
        fileprefix = '~/Dropbox/NYU/eyetracking/';
  
        
    case 3
        fileprefix = '/Volumes/ShaniWSQBackupHD/Dropbox/eyetracking/';
     
end

        datafolder = [fileprefix 'data/'];
        figurefolder = [fileprefix 'figures/'];
        filefolder = [fileprefix 'measures/'];
        matfolder = [fileprefix 'matFiles/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction
%********************** makejulymeasures *********************
function [newx newy] = raw2clean(t,bigs,isubj)

f = t.fixations;
c = t.saccades;
b = t.blinks;
g = t.gaze;

x = g.rawx;
y = g.rawy;
p = g.pupil;
gt = g.time;

s=bigs(isubj);
if(isempty(s.subjdelay))
    s.subjdelay = 0;
end
droplength = 10+s.subjdelay;
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

% pad with nans if didnt watch whole thing
if(length(keeprange)<nummeasurements)
    xgaze(1,length(keeprange)+1:nummeasurements) = nan;
    ygaze(1,length(keeprange)+1:nummeasurements) = nan;
    pgaze(1,length(keeprange)+1:nummeasurements) = nan;
    gazetimes(1,length(keeprange)+1:nummeasurements)= nan;
end


newx = removeblinks(xgaze,b,gazetimes); newx = newx';
newy = removeblinks(ygaze,b,gazetimes); newy = newy';
newp = removeblinks(pgaze,b,gazetimes); newp = newp';
newpdif = newp(2:end) - newp(1:end-1);
newpdif(end+1) = nan;


% ntpnts = length(newx);
%
% meanmatx = repmat(nanmean(newx),1,ntpnts);
% meanmaty = repmat(nanmean(newy),1,ntpnts);
%
% newx = newx-meanmatx; % don't want to change the variance
% newy = newy-meanmaty;

% normx = (newx./meanmatx)-1;
% normy = (newy./meanmaty)-1;
%
% varmat1x = repmat(nanstd(newx),1,ntpnts);
% varmat1y = repmat(nanstd(newy),1,ntpnts);
%
% x1 = ((newx-medmatx)./varmat1x);
% y1 = ((newy-medmaty)./varmat1y);
%
% varmat2 = repmat(nanstd([newx newy]),1,ntpnts);
% x2 = ((newx-medmatx)./varmat2);
% y2 = ((newy-medmaty)./varmat2);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction
%********************** measures *********************
function d = getvariousmeasures(x,y,p,pdif);
%vals2measure =
%{'edfname','xvar','yvar','xRangeRaw','yRangeRaw','xRangeJoint','yRangeJoint','pRangeRaw','pDifRangeRaw','meanP','meanPdif'};
d = [];

% clip to leave out outliers
newX = x - nanmean(x);
tempX = newX(~isnan(newX));
sortX = sort(tempX);
xTH = nanmedian(sortX(ceil(9*end/10):end));
clipX = tempX(abs(tempX)<xTH);

newY = y - nanmean(y);
tempY = newY(~isnan(newY));
sortY = sort(tempY);
yTH = nanmedian(sortY(ceil(9*end/10):end));
clipY = tempY(abs(tempY)<yTH);

newP = p - nanmean(p);
tempP = newP(~isnan(newP));
sortP = sort(tempP);
pTH = nanmedian(sortP(ceil(9*end/10):end));
clipP = tempP(abs(tempP)<pTH);

newPD = pdif - nanmean(pdif);
tempPD = newPD(~isnan(newPD));
sortPD = sort(tempPD);
pdTH = nanmedian(sortPD(ceil(9*end/10):end));
clipPD = tempPD(abs(tempPD)<pdTH);


d.xvar = nanstd(clipX);
d.yvar = nanstd(clipY);
d.pvar = nanstd(clipP);
d.pdvar = nanstd(clipPD);

d.meanP = nanmean(clipP);
d.meanPdif = nanmean(clipPD);

d.xRangeRaw = nanmax(clipX)-nanmin(clipX);
d.yRangeRaw = nanmax(clipY)-nanmin(clipY);
d.pRangeRaw = nanmax(clipP)-nanmin(clipP);
d.pDifRangeRaw = nanmax(clipPD)-nanmin(clipPD);

% variance across cycles, and TTA:

nummperp = 500*10; % this is (mpersec*secperpart)
numcycles = 5;
cycleTime = length(x)/numcycles;
sideTime = cycleTime/4;

xByTrials = reshape(newX,cycleTime,numcycles);
yByTrials = reshape(newY,cycleTime,numcycles);

d.xvarAcrossCycles= nanmedian(nanvar(xByTrials'));
d.yvarAcrossCycles = nanmedian(nanvar(yByTrials'));

xtta = nanmean(xByTrials');
ytta = nanmean(yByTrials');

% vertical stripiness, boxiness:

topind = 1:sideTime;
rightind = 1*sideTime+1:2*sideTime;
botind =   2*sideTime+1:3*sideTime;
leftind =  3*sideTime+1:4*sideTime;

% figure, plot(1:sideTime,-ytta(topind),'g.'), hold on, plot(1:sideTime,ytta(botind),'r.');
% figure, plot(1:sideTime,-ytta(topind),'g.'), ...
%     hold on, plot(1:sideTime,-ytta(botind),'r.'), ...
%     hold on, plot(xtta(rightind),1:sideTime,'b.'), ...
%     hold on, plot(xtta(leftind),1:sideTime,'m.');

yTopMatrix = yByTrials(topind,:);
yTop = -yTopMatrix(:);
yBotMatrix = yByTrials(botind,:);
yBot = -yBotMatrix(:);

yDif = yBot-yTop;

yDif = sort(yBot)-sort(yTop);

minYtop = min(yTop);
maxYbot = max(yBot);

d.stripiness = sum(and(-clipY<maxYbot,-clipY>minYtop))/length(clipY);

d = gettrackerboxiness(d,xtta,ytta,topind,rightind,botind,leftind);


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = gettrackerboxiness(d,xtta,ytta,topind,rightind,botind,leftind);

x = xtta/nanstd(xtta); % normalize for the model comparison
y = ytta/nanstd(ytta);

xtop = x(topind);
xright = x(rightind);
xbottom = x(botind);
xleft = x(leftind);

ytop = y(topind);
yright = y(rightind);
ybottom = y(botind);
yleft = y(leftind);

xmodeltop = linspace(-2,2,length(topind)); % x value on top of box
xmodelright = repmat(2,1,length(topind)); % x value on right of box
xmodelbottom = linspace(2,-2,length(topind)); % x value on bottom of box
xmodelleft = repmat(-2,1,length(topind)); % x value on left of box
xmodel = [xmodeltop xmodelright xmodelbottom xmodelleft]; % one cycle of x values

ymodeltop = linspace(-2,2,length(topind)); % x value on top of box
ymodelright = repmat(2,1,length(topind)); % x value on right of box
ymodelbottom = linspace(2,-2,length(topind)); % x value on bottom of box
ymodelleft = repmat(-2,1,length(topind)); % x value on left of box
ymodel = [ymodeltop ymodelright ymodelbottom ymodelleft]; % one cycle of x values


xerrtop = xmodeltop - xtop;
yerrtop = ymodeltop - ytop;

xerrrit = xmodelright - xright;
yerrrit = ymodelright - yright;

xerrbot = xmodelbottom - xbottom;
yerrbot = ymodelbottom - ybottom;

xerrlef = xmodelleft - xleft;
yerrlef = ymodelleft - yleft;



xerroverall = xmodel - x;
yerroverall = ymodel - y;


d.xr2 = 1-(nanvar(xerroverall)/nanvar(x));
d.yr2 = 1-(nanvar(yerroverall)/nanvar(y));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function datamatrix = makedatamatrix(d,s,isubj,datamatrix);
edfname = stripext(s(isubj).edfname);
%vals2measure =
%{'edfname','xvar','yvar','xvarAcrossCycles','yvarAcrossCycles','xr2','yr2','stripiness','aspectRatio','meanP','meanPdif'};

datamatrix(isubj,1) = str2double(edfname(3:end));
datamatrix(isubj,2) = d.xvar;
datamatrix(isubj,3) = d.yvar;
datamatrix(isubj,4) = d.xvarAcrossCycles;
datamatrix(isubj,5) = d.yvarAcrossCycles;
datamatrix(isubj,6) = d.xr2;
datamatrix(isubj,7) = d.yr2;
datamatrix(isubj,8) = d.stripiness;
datamatrix(isubj,9) = d.yRangeRaw/d.xRangeRaw;
datamatrix(isubj,10) = d.meanP;
datamatrix(isubj,11) = d.meanPdif;


%____________________________________________________%


function plotsubplots(plotvals, ifile, x1, y1, x2, y2, newp, newpdif)
%____________________________________________________%
whichplot = toplot.whichplot;
xtta1 = reshape(x1,5,length(x1)/5);
ytta1 = reshape(y1,5,length(y1)/5);
xtta2 = reshape(x2,5,length(x2)/5);
ytta2 = reshape(y2,5,length(y2)/5);

switch whichplot
    case 1
        numfiles=toplot.numfiles;
        numfilesperfig=toplot.numfilesperfig;
        numfigs=toplot.numfigs;
        numrows = toplot.numrowsperfig;
        numcols = toplot.numcolsperfig;
        plotlistnums = 1:2:10;
        currsubplot = plotlistnums(mod(ifile-1,numfilesperfig)+1);
        currfig = ceil(mod(ifile/numfilesperfig,numfigs+1));
        
        figure(currfig)
        subplot(numrows,numcols,currsubplot)
        plot(x1,y1)
        title('norm individually')
        ylabel(sprintf('file num = %s',num2str(ifile)));
        xlim([-3,3])
        ylim([-3,3])
        axis equal
        
        
        subplot(numrows,numcols,currsubplot+1);
        plot(x2,y2)
        title('norm together')
        xlim([-3,3])
        ylim([-3,3])
        axis equal
        
    case 2
        
        
        figure(ifile)
        subplot(221)
        % full time course
        plot([x2,y2]);
        title('full time course')
        
        subplot(222)
        % tta time course
        plot([nanmean(xtta),nanmean(ytta)])
        title('tta time course')
        
        subplot(223)
        % pupil diameter
        plot(newp)
        title('pupil diameter')
        
        subplot(224)
        % pupil change
        plot(newpdif)
        title('change in pupil diameter')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = getaspectratios(d,x,y);

topind = 1:length(x1)/4;
rightind = (length(x1)/4)+1:length(x1)/2;
botind = (length(x1)/2)+1:3*length(x1)/4;
leftind = (3*length(x1)/4)+1:length(x1);

d.height1 = nanmax([y1(topind) y1(botind)]) - nanmin([y1(topind) y1(botind)]);
d.width1 = nanmax([x1(rightind) x1(leftind)]) - nanmin([x1(rightind) x1(leftind)]);
d.ar1 = d.height1/d.width1;

d.height2 = nanmax([y2(topind) y2(botind)]) - nanmin([y2(topind) y2(botind)]);
d.width2 = nanmax([x2(rightind) x2(leftind)]) - nanmin([x2(rightind) x2(leftind)]);
d.ar2 = d.height2/d.width2;



%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction
%********************** measures *********************
function d = OLDgetvariousmeasures(x,y, newp, newpdif);

d = [];
d = getreplicability(d,x1,y1,x2,y2); % gives d.xreplicability d.yreplicability

xtta1 = nanmean(reshape(x1,5,length(x1)/5));
ytta1 = nanmean(reshape(y1,5,length(y1)/5));
xtta2 = nanmean(reshape(x2,5,length(x2)/5));
ytta2 = nanmean(reshape(y2,5,length(y2)/5));

d = gettrackerboxiness(d,xtta1,ytta1); % gives d.xr2, d.yr2
d.xvar1 = nanvar(xtta1);      % gives d.xvar1, d.xvar2, d.yvar1, d.yvar2
d.xvar2 = nanvar(xtta2);
d.yvar1 = nanvar(ytta1);
d.yvar2 = nanvar(ytta2);


d.pvar = nanvar(newp);
d.pdifvar = nanvar(newpdif);
d.pdifmean = nanmean(newpdif);
d.pdifmedian = nanmedian(newpdif);

d = getaspectratios(d,xtta1, ytta1, xtta2, ytta2); % gets d.height, d.width, d.ar x2
