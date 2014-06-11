
function analysisForMethodsPaper_alt03dist
% ************************************************************************
onKirkwood = 0;
filePrefix = setPrefix(onKirkwood);
doPlot = 1;
jF = figure;

for iList = 1:2 % load m, 15 fields: x,y and xy: mean var, max var, integrity, TTA-R2 and R2
    switch iList
        case 1
            listType = 'test';
            load([filePrefix 'figures/figsForMayPaper/data/TESTmeasures.mat'])
        case 2
            listType = 'control';
            load([filePrefix 'figures/figsForMayPaper/data/CONTROLmeasures.mat'])
    end
    
    numSubj = length(list);
    
    %--------------
    % MEASURES
    % need to bootstrap the summary stats
    nboots = 5000;
    nsamples = 40;    
    
    d{iList}.xyMeanVar = bootstrapDist(m.xVar(:),nboots,nsamples);
    d{iList}.xyIntegrity = bootstrapDist(m.yR2TTA(:),nboots,nsamples);    
    
    
    if(doPlot)
        % just for visualization:
        % leave out extreme outliers so we can see the details
        sortMeanVxy = sort(d{iList}.xyMeanVar,1,'descend');
        sortMeanVxy = sortMeanVxy(~isnan(sortMeanVxy));
        sortMeanVxy = sortMeanVxy(1000:end);
        
        sortIntegrityxy = sort(d{iList}.xyIntegrity,1,'descend');
        sortIntegrityxy = sortIntegrityxy(~isnan(sortIntegrityxy));
        sortIntegrityxy = sortIntegrityxy(1000:end);
        

        
        
        switch iList
            case 1
                figure(jF)
                
                subplot(2,2,1)
                hist(sortMeanVxy);
                title([listType ' Distribution: Mean Variance']);
                ylim([0 4000])
                legend(sprintf('Median = %0.1g', (nanmedian(d{iList}.xyMeanVar))))
                
                subplot(2,2,2)
                hist(sortIntegrityxy);
                title([listType ' Distribution: Box Integrity']);
                ylim([0 4000])
                legend(sprintf('Median = %0.1g', (nanmedian(d{iList}.xyIntegrity))))
                

                
            case 2
                figure(jF)
                
                subplot(2,2,3)
                hist(sortMeanVxy);
                title([listType ' Distribution: Mean Variance']);
                ylim([0 4000])
                legend(sprintf('Median = %0.1g', (nanmedian(d{iList}.xyMeanVar))))
                
                subplot(2,2,4)
                hist(sortIntegrityxy);
                title([listType ' Distribution: Box Integrity']);
                ylim([0 4000])
                legend(sprintf('Median = %0.1g', (nanmedian(d{iList}.xyIntegrity))))
                
                
        end
        
        
    end
    figure(jF)
    print('-djpeg',[filePrefix 'figures/figsForMayPaper/summary/compareDists.jpg']);
    
    
    
    
    
    
    
    % ************************************************************************
end %**** for list = 1:2 ********* %

[s.hV,s.pV] = ttest2(d{1}.xyMeanVar,d{2}.xyMeanVar,.05);
[s.hI,s.pI] = ttest2(d{1}.xyIntegrity,d{2}.xyIntegrity,.05);


meanTest = [mean(d{1}.xyMeanVar) mean(d{1}.xyIntegrity)];
meanControl = [mean(d{2}.xyMeanVar) mean(d{2}.xyIntegrity)];


medianTest = [median(d{1}.xyMeanVar) median(d{1}.xyIntegrity)];
medianControl = [median(d{2}.xyMeanVar) median(d{2}.xyIntegrity)];


% ************************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([filePrefix 'figures/figsForMayPaper/data/Distributions.mat'], 'd','s');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataTable(1,1:2) = [s.hV, s.hI];
dataTable(2,1:2) = [s.pV, s.pI];
dataTable(3,1:2) = meanTest;
dataTable(4,1:2) = meanControl;
dataTable(5,1:2) = medianTest;
dataTable(6,1:2) = medianControl;

xlswrite([filePrefix 'figures/figsForMayPaper/data/summaryTable.xls'],dataTable)

%___________________________________________%
function filePrefix = setPrefix(onKirkwood)
if(onKirkwood)
    filePrefix = '/Volumes/ShaniWSQBackupHD/Dropbox/';
else
    filePrefix = '~/Dropbox/';
end



function newdist = bootstrapDist(cVar1,nboots,nsamples);

numOrigMeas = length(cVar1);
toPickFrom = cVar1;

indexMatrix = rand(nboots,nsamples);
indexMatrix = ceil(numOrigMeas*indexMatrix);

newpicks = toPickFrom(indexMatrix);
newdist = nanmean(newpicks');

function save deleted bits

[s.hxV,s.pxV] = ttest2(d{1}.xVar,d{2}.xVar,.05);
[s.hyV,s.pyV] = ttest2(d{1}.yVar,d{2}.yVar,.05);
[s.hxR,s.pxR] = ttest2(d{1}.xR2,d{2}.xR2,.05);
[s.hyR,s.pyR] = ttest2(d{1}.yR2,d{2}.yR2,.05);
[s.hxRT,s.pxRT] = ttest2(d{1}.xR2TTA,d{2}.xR2TTA,.05);
[s.hyRT,s.pyRT] = ttest2(d{1}.yR2TTA,d{2}.yR2TTA,.05)

meanTest = [mean(d{1}.xVar) mean(d{1}.yVar) mean(d{1}.xR2) mean(d{1}.yR2) mean(d{1}.xR2TTA) mean(d{1}.yR2TTA)];
meanControl = [mean(d{2}.xVar) mean(d{2}.yVar) mean(d{2}.xR2) mean(d{2}.yR2) mean(d{2}.xR2TTA) mean(d{2}.yR2TTA)];
medianTest = [median(d{1}.xVar) median(d{1}.yVar) median(d{1}.xR2) median(d{1}.yR2) median(d{1}.xR2TTA) median(d{1}.yR2TTA)];
medianControl = [median(d{2}.xVar) median(d{2}.yVar) median(d{2}.xR2) median(d{2}.yR2) median(d{2}.xR2TTA) median(d{2}.yR2TTA)];

