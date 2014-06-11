function analysisJointNorm_01preprocess
% ************************************************************************
onKirkwood = 1;
redoTestNormalization = 1;
redoControlNormalization = 0;
% ************************************************************************
if(onKirkwood)
    filePrefix = '/Volumes/ShaniWSQBackupHD/Dropbox/';
else
    filePrefix = '~/Dropbox/';
end
% ************************************************************************

% TEST
% ************************************************************************
doPlotTestSummary = 0;
doPlotBoxesTestIndividual = 0;

figureFolderTest = [filePrefix 'figures/figsForMayPaper/individuals/tests/'];
dataFolder = [filePrefix 'VA/ETdata/VA/'];
cd(dataFolder)

% load list of test population
TESTlist = loadTestPop(onKirkwood);
% structure is {}{(1) EDF filename, (2) test#, (3) movie, (4)subj code, (5) date, (6) initials, (7) which eye, (8) diagnosis, (9), notes}

if(redoTestNormalization)
    % Get normalized time courses for each:
    for iTest = 1:length(TESTlist)
        testData = TESTlist{iTest};
        fileName = testData{1};
        TCs{iTest} = preProcessWaka(fileName,onKirkwood); %4 fields: TTAxCW; TTAyCW; stackedXcw; stackedYcw; %columns contain TC for each repeat
        
        if(doPlotBoxesTestIndividual)
            plotIndividual(TCs{iTest})
            figureName = [testData{6} '_' testData{4} '_' testData{5} '_' testData{7}];
            print('-djpeg',[figureFolderTest figureName  '.jpg'])
            close all
        end
    end
    
    clear iTest testData fileName
    
    save([filePrefix 'figures/figsForMayPaper/data/TESTgraphingdata.mat'],'TESTlist','TCs')
else
    load([filePrefix 'figures/figsForMayPaper/data/TESTgraphingdata.mat'])
end


if(doPlotTestSummary)
   
    plotSummaries(TESTlist,TCs,filePrefix,'tests')
end

clear TESTlist TCs 
close all
% ************************************************************************
%
% CONTROLS
% ************************************************************************
doPlotControlSummary = 0;
doPlotBoxesControlIndividual = 0;

figureFolderControl = [filePrefix 'figures/figsForMayPaper/individuals/controls/'];
dataFolder = [filePrefix 'VA/ETdata/WSQ/waka/'];
cd(dataFolder)
% load list of control population
CONTROLlist = loadControlPop(onKirkwood);
% structure is {}{(1) EDF filename, (2) control#, (3) movie, (4)subj type, (5) year, (6) file#, (7) which eye, (8) light or dark, (9), delay}

if(redoControlNormalization)
    % Get normalized time courses for each:
    for iControl = 1:length(CONTROLlist)
        controlData = CONTROLlist{iControl};
        fileName = controlData{1};
        TCs{iControl} = preProcessWaka(fileName,onKirkwood); %4 fields: TTAxCW; TTAyCW; stackedXcw; stackedYcw; %columns contain TC for each repeat
        
        if(doPlotBoxesTestIndividual)
            plotIndividual(TCs{iControl})
            figureName = [testData{2} '_' testData{8} '_' testData{6}];
            print('-djpeg',[figureFolderTest figureName  '.jpg'])
            close all
        end
    end
    
    clear iControl controlData fileName
    
    save([filePrefix 'figures/figsForMayPaper/data/CONTROLgraphingdata.mat'],'CONTROLlist','TCs')
else
    load([filePrefix 'figures/figsForMayPaper/data/CONTROLgraphingdata.mat'])
end

clear CONTROLlist TCs 

close all






%*******************************************%
%**        **SUBFUNCTIONS **              **%
%*******************************************%

%___________________________________________%
function plotIndividual(TCs)
%___________________________________________%

figure
subplot(221),
plot(TCs.TTAxCW, -TCs.TTAyCW,'k.');
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
axis square
title(['TTA'])

subplot(222),
plot(TCs.stackedXcw(:,1), -TCs.stackedYcw(:,1),'r.');
hold on
plot(TCs.stackedXcw(:,2), -TCs.stackedYcw(:,2),'g.');
plot(TCs.stackedXcw(:,3), -TCs.stackedYcw(:,3),'c.');
plot(TCs.stackedXcw(:,4), -TCs.stackedYcw(:,4),'m.');
plot(TCs.stackedXcw(:,5), -TCs.stackedYcw(:,5),'b.');
xlim([-2.5 2.5]);
ylim([-2.5 2.5]);
axis square
title(['Five repeats'])

subplot(223)
plot([TCs.TTAxCW,-TCs.TTAyCW])
title('TTA time course')
ylabel('% signal change')
xlabel('Time in 2ms')

subplot(224)
plot([TCs.stackedXcw(:),-TCs.stackedYcw(:)])
title('Full time course')
ylabel('% signal change')
xlabel('Time in 2ms')






%____________________________________________________%

function plotSummaries(list,TCs,filePrefix,listType)
%____________________________________________________%

summaryFolder = [filePrefix 'figures/figsForMayPaper/summary/' listType '/'];

numFilesPerFig = 6;
numFiles = length(list);
numFigs = ceil(numFiles/numFilesPerFig);

% structure is {}{(1) EDF filename, (2) control#, (3) movie, (4)subj type, (5) year, (6) file#, (7) which eye, (8) light or dark, (9), delay}
switch listType
    case 'tests'
        t = [6, 7, 5];
    case 'controls'
        t = [6,10,10]
end

count = 0;
for iFig = 1:numFigs-1
    figure(iFig)
    for iSub = 1:2:numFilesPerFig*2
        
        count = count + 1;
        TC = TCs{count};
        d = list{count};
        d{10} = '';
        title1 = [d{t(1)} ' ' d{t(2)}];
        title2 = d{t(3)};
        
        
        subplot(3,4,iSub),
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
        
        
        subplot(3,4,iSub+1)
        g = plot([TC.TTAxCW,-TC.TTAyCW]);
        title(title2)
        axis off
    end
    
    sumFigName = [listType 'SummaryPlots_p' num2str(iFig)]
    print('-djpeg',[summaryFolder sumFigName  '.jpg'])
    close all
    
end
% may be fewer than 6 left, check:
numLeft = numFiles-count;

figure
for iSub = 1:2:numLeft*2
    count = count + 1;
    TC = TCs{count};
    d = list{count};
    d{10} = '';
    title1 = [d{t(1)} ' ' d{t(2)}];
    title2 = d{t(3)};
    
    
    subplot(3,4,iSub),
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
    
    subplot(3,4,iSub+1)
    g = plot([TC.TTAxCW,-TC.TTAyCW]);
    title(title2)
    axis off
    
end

    sumFigName = [listType 'SummaryPlots_p' num2str(numFigs)]
    print('-djpeg',[summaryFolder sumFigName  '.jpg'])
    close all














