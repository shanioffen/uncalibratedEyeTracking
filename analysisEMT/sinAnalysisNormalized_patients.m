function [goodTCpatients, badPhaseTC, badFitTC] = sinAnalysisNormalized_patients(toggle,makeFigs)
% function trackerbatchanalysis_so(toggle,listnum)
% toggle: 1 for eyetracker, 2 for most places, 3 for kirkwood

toggle = 3;
makeFigs = 0;

[datafolder filefolder matfolder] = setvariables(toggle);
cd(datafolder)

%s.listlabels = {'EDF filename','which eye','De-ID code','ID code','delay','time','sitting','chinrest','inpatient'};
redomatwrite = 0;
if(redomatwrite)
    
    listmakemat = loadVAlist(3);
    nmatsubj = length(listmakemat);
    
    for isubj = 1:nmatsubj
        subj = listmakemat{isubj};
        s(isubj).edfname = subj{1};
        if(length(subj)>5)
            s(isubj).subjeye = subj{2};
            s(isubj).subjid = subj{3};
            s(isubj).subjcode = subj{4};
            s(isubj).subjdelay = subj{5};
            s(isubj).subjtime = subj{6};
            s(isubj).subjsit = subj{7};
            s(isubj).subjchinrest = subj{8};
            s(isubj).subjtype = subj{9};
        else
            s(isubj).subjdelay = subj{2};
            s(isubj).subjeye = subj{3};
            s(isubj).subjid = subj{4};
            s(isubj).subjcode = subj{5};
            s(isubj).subjsit = 1;
            s(isubj).subjchinrest=1;
        end
        
        
        if(isempty(s(isubj).subjdelay)), s(isubj).subjdelay = 0; end
        t = mgleyelinkedfread1(s(isubj).edfname);
        if(nanmin(t.gaze.rawx) == nanmax(t.gaze.rawx))
            t = mgleyelinkedfread(s(isubj).edfname);
        end
        if(nanmin(t.gaze.rawx) == nanmax(t.gaze.rawx))
            (sprintf('%s is not raw data; isubj = %i',s(isubj).edfname,isubj))
            return
        end
        save([matfolder s(isubj).subjid '_' s(isubj).edfname '.mat'], 't')
        clear t
    end
    clear s
end


list = loadVAlist(4);
nsubj = length(list);
vals2measure = {'A1x','B1x','C1x','R2x','A1y','B1y','C1y','R2y'};
nmeasures = length(vals2measure);
datamatrix = nan*ones(nsubj,nmeasures);

goodOutCount = 0;
goodTCpatients = [];

badPhaseCount = 0;
badPhaseTC = [];

badFitCount = 0;
badFitTC = [];

redoDataMatrix=0;
if(redoDataMatrix)
    h=waitbar(0,'making data matrix');
    for isubj = 1:nsubj
        subj = list{isubj};
        s(isubj).edfname = subj{1};
        
        if(length(subj)>5)
            s(isubj).subjeye = subj{2};
            s(isubj).subjid = subj{3};
            s(isubj).subjcode = subj{4};
            s(isubj).subjdelay = subj{5};
            s(isubj).subjtime = subj{6};
            s(isubj).subjsit = subj{7};
            s(isubj).subjchinrest = subj{8};
            s(isubj).subjtype = subj{9};
        else
            s(isubj).subjdelay = subj{2};
            s(isubj).subjeye = subj{3};
            s(isubj).subjid = subj{4};
            s(isubj).subjcode = subj{5};
            s(isubj).subjsit = 'yes';
            s(isubj).subjchinrest='yes';
        end
        
        if(isempty(s(isubj).subjdelay)), s(isubj).subjdelay = 0; end
        waitbar(isubj/nsubj,h);
        subjid = s(isubj).subjid;
        edfname = s(isubj).edfname;
        issit = strcmp(s(isubj).subjsit,'yes');
        ischin = strcmp(s(isubj).subjchinrest,'yes');
        
        matfilename = [matfolder subjid '_' edfname '.mat'];
        
        load(matfilename); % this loads t
        [newx newy] = raw2clean(t,s,isubj);
        
        % check to make sure at least 60% data:
        if sum(~isnan(newx))/length(newx)<.6
            datamatrix(isubj,1) = NaN;
            datamatrix(isubj,2) = NaN;
            datamatrix(isubj,3) = NaN;
            datamatrix(isubj,4) = NaN;
            
            datamatrix(isubj,5) = NaN;
            datamatrix(isubj,6) = NaN;
            datamatrix(isubj,7) = NaN;
            datamatrix(isubj,8) = NaN;
            
        elseif ~(or(issit,ischin))
            datamatrix(isubj,1) = NaN;
            datamatrix(isubj,2) = NaN;
            datamatrix(isubj,3) = NaN;
            datamatrix(isubj,4) = NaN;
            
            datamatrix(isubj,5) = NaN;
            datamatrix(isubj,6) = NaN;
            datamatrix(isubj,7) = NaN;
            datamatrix(isubj,8) = NaN;
        else
            
            
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
            phaseDif = y1vals(3)-x1vals(3);
            
            if(makeFigs)
                if(and(gx1.adjrsquare>.95,gy1.adjrsquare>.95))
                    goodOutCount = goodOutCount+1;
                    goodTCpatients{goodOutCount} = [x' y'];
                    
                elseif(or(gx1.adjrsquare < .5, gy1.adjrsquare < .5))
                    badFitCount = badFitCount +1;
                    badFitTC{badFitCount} = [x' y'];
                end
                
                if(abs(phaseDif-45)>2.5)
                    badPhaseCount = badPhaseCount + 1;
                    badPhaseTC{badPhaseCount} = [x' y'];
                    
                end
            end
            
        end
        

        
        clear t newx newy x y fox1 gx1 foy1 gy1 xind yind
        
    end
    
        save goodTCpatients.mat goodTCpatients
        save badFitTC.mat badFitTC
        save badPhaseTC.mat badPhaseTC
        clear goodTCpatients badFitTC badPhaseTC
    
    measname = ['sinAnalysisPatientsNormalized'];
    xlswrite([filefolder '/' measname],datamatrix);
    save([filefolder '/' measname '.mat'],'datamatrix')
else
    measname = ['sinAnalysisPatientsNormalized'];
    load([filefolder '/' measname '.mat'])
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction
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
