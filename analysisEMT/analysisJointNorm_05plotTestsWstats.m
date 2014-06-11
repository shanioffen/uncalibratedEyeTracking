afunction analysisForMethodsPaper_05plotTestsWstats

onKirkwood = 1;

if(onKirkwood)
    filePrefixKeep = '/Volumes/ShaniWSQBackupHD/Dropbox/';
else
    filePrefixKeep = '~/Dropbox/';
end






[zMat hMat] = analysisForMethodsPaper_04stats;
load([filePrefixKeep 'figures/figsForMayPaper/data/TESTgraphingdata.mat'])
clear filePrefix

% save out as .xls file
xlswrite([filePrefixKeep 'figures/figsForMayPaper/data/testZTable.xls'],zMat);
xlswrite([filePrefixKeep 'figures/figsForMayPaper/data/testHTable.xls'],hMat);


% trim so can see both in bar form:
zMat(zMat>=25) = 25;
zMat(zMat<=-25) = -25; 

summaryFolder = [filePrefixKeep 'figures/figsForMayPaper/summary/testsWstats/'];

numFilesPerFig = 6;
numFiles = length(TESTlist);
numFigs = ceil(numFiles/numFilesPerFig);
t = [6, 7, 5];

count = 0;
for iFig = 1:numFigs-1
    figure(iFig)
    for iSub = 1:3:numFilesPerFig*3
        
        count = count + 1;
        TC = TCs{count};
        d = TESTlist{count};
        d{10} = '';
        title1 = [d{t(1)} ' ' d{t(2)}];
        title2 = d{t(3)};
        
        
        subplot(3,6,iSub),
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
        
        
        subplot(3,6,iSub+1)
        g = plot([TC.TTAxCW,-TC.TTAyCW]);
        title(title2)
        
        subplot(3,6,iSub+2)
        bar(hMat(count,:)'.*zMat(count,:)')
        ylim([-10 10])
        
    end
    
    sumFigName = ['testSummaryPlots_p' num2str(iFig)]
    print('-djpeg',[summaryFolder sumFigName  '.jpg'])
    close all
    
end
% may be fewer than 6 left, check:
numLeft = numFiles-count;

figure
for iSub = 1:3:numLeft*3
    count = count + 1;
    TC = TCs{count};
    d = TESTlist{count};
    d{10} = '';
    title1 = [d{t(1)} ' ' d{t(2)}];
    title2 = d{t(3)};
    
    
    subplot(3,6,iSub),
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
    
    subplot(3,6,iSub+1)
    g = plot([TC.TTAxCW,-TC.TTAyCW]);
    title(title2)

    
    subplot(3,6,iSub+2)
    bar(hMat(count,:)'.*zMat(count,:)')
    ylim([-10 10])
end

sumFigName = ['testSummaryPlots_p' num2str(numFigs)]
print('-djpeg',[summaryFolder sumFigName  '.jpg'])
close all

