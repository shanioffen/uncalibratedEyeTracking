function analysisForMethodsPaper_06plotControlsWstats

onKirkwood = 1;

if(onKirkwood)
    filePrefixKeep = '/Volumes/ShaniWSQBackupHD/Dropbox/';
else
    filePrefixKeep = '~/Dropbox/';
end
summaryFolder = [filePrefixKeep 'figures/figsForMayPaper/summary/controlsWstats/'];

listType = 'control';
load([filePrefixKeep 'figures/figsForMayPaper/data/CONTROLmeasures.mat'])

numSubj = length(list);

%--------------
% MEASURES
% need to bootstrap the summary stats
% but for each control, need to leave out the observer being analyzed

for iSubj = 1:numSubj
    
    if iSubj == 1
        subIndex = 2:numSubj;
    elseif iSubj == numSubj
        subIndex = 1:numSubj-1;
    else
        subIndex = [1:(numSubj-1) (numSubj+1):numSubj];
    end
    
    nboots = 5000;
    nsamples = 40;
    
    controlVarSub = bootstrapDist(m.xymeanVar(subIndex),nboots,nsamples);
    controlIntSub = bootstrapDist(m.xyTTAintegrity(subIndex),nboots,nsamples);
    
    % For each control time-course, calculate z-score using leave-one-out
    muVec = [nanmean(controlVarSub) nanmean(controlIntSub)];
    sigmaVec = [nanstd(controlVarSub) nanstd(controlIntSub)];
    testVals = [m.xymeanVar(iSubj) m.xyTTAintegrity(iSubj)];
    
    th = 10;
    zMat = (testVals-muVec)./sigmaVec;
    hMatControl(iSubj,:) = or(zMat>th,zMat<-th);
    zMatControl(iSubj,:) = zMat;
    
    clear subIndex controlVarSub controlIntSub muVec sigmaVec testVals zmat
    
end

% save out as .xls file
xlswrite([filePrefixKeep 'figures/figsForMayPaper/data/controlZTable.xls'],zMatControl);
xlswrite([filePrefixKeep 'figures/figsForMayPaper/data/controlHTable.xls'],hMatControl);


% clip at +/- 10 so can see range in bar form:
zMatControl(zMatControl>=25) = 25;
zMatControl(zMatControl<=-25) = -25;




% make plots with stats
load([filePrefixKeep 'figures/figsForMayPaper/data/CONTROLgraphingdata.mat']);

numFilesPerFig = 6;
numFigs = ceil(numSubj/numFilesPerFig);

count = 0;
for iFig = 1:numFigs-1
    figure(iFig)
    for iPlot = 1:3:numFilesPerFig*3
        
        count = count + 1;
        TC = TCs{count};
        title1 = list{count}{6};
        title2 = sprintf('Var = %0.2g, Int = %0.2g',...
            m.xymeanVar(count), m.xyTTAintegrity(count));
        title3 = sprintf('hVar = %0.2g, hInt = %0.2g',...
            hMatControl(count,1), hMatControl(count,2));
        
        subplot(3,6,iPlot),
        h = plot(TC.stackedXcw(:,1), -TC.stackedYcw(:,1),'r.');
        hold on
        plot(TC.stackedXcw(:,2), -TC.stackedYcw(:,2),'g.');
        plot(TC.stackedXcw(:,3), -TC.stackedYcw(:,3),'c.');
        plot(TC.stackedXcw(:,4), -TC.stackedYcw(:,4),'m.');
        plot(TC.stackedXcw(:,5), -TC.stackedYcw(:,5),'b.');
        xlim([-2.5 2.5]);
        ylim([-2.5 2.5]);
        axis equal
        title(title1)
        axis off
        
        
        subplot(3,6,iPlot+1)
        g = plot([TC.TTAxCW,-TC.TTAyCW]);
        title(title2)
        
        subplot(3,6,iPlot+2)
        bar(hMatControl(count,:)'.*zMatControl(count,:)')
        title(title3)
        ylim([-10 10])
        
    end
    
    sumFigName = ['controlSummaryPlots_p' num2str(iFig)]
    print('-djpeg',[summaryFolder sumFigName  '.jpg'])
    close all
    
end
% may be fewer than 6 left, check:
numLeft = numSubj-count;

figure
for iPlot= 1:3:numLeft*3
    count = count + 1;
    TC = TCs{count};
    title1 = list{count}{6};
    title2 = sprintf('Var = %0.2g, Int = %0.2g',...
        m.xymeanVar(count), m.xyTTAintegrity(count));
    title3 = sprintf('hVar = %0.2g, hInt = %0.2g',...
        hMatControl(count,1), hMatControl(count,2));
    
    
    subplot(3,6,iPlot),
    h = plot(TC.stackedXcw(:,1), -TC.stackedYcw(:,1),'r.');
    hold on
    plot(TC.stackedXcw(:,2), -TC.stackedYcw(:,2),'g.');
    plot(TC.stackedXcw(:,3), -TC.stackedYcw(:,3),'c.');
    plot(TC.stackedXcw(:,4), -TC.stackedYcw(:,4),'m.');
    plot(TC.stackedXcw(:,5), -TC.stackedYcw(:,5),'b.');
    xlim([-2.5 2.5]);
    ylim([-2.5 2.5]);
    axis square
    title(title1)
    axis off
    
    subplot(3,6,iPlot+1)
    g = plot([TC.TTAxCW,-TC.TTAyCW]);
    title(title2)
    
    
    subplot(3,6,iPlot+2)
    bar(hMatControl(count,:)'.*zMatControl(count,:)')
    ylim([-10 10])
end

sumFigName = ['controlSummaryPlots_p' num2str(numFigs)]
print('-djpeg',[summaryFolder sumFigName  '.jpg'])
close all


function newdist = bootstrapDist(cVar1,nboots,nsamples);

numOrigMeas = length(cVar1);
toPickFrom = cVar1;

indexMatrix = rand(nboots,nsamples);
indexMatrix = ceil(numOrigMeas*indexMatrix);

newpicks = toPickFrom(indexMatrix);
newdist = nanmean(newpicks');

