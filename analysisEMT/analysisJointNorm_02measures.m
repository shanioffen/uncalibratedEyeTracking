function analysisForMethodsPaper_02measures
% ************************************************************************
onKirkwood = 1;
filePrefix = setPrefix(onKirkwood);
% have to first run analysisForMethodsPaper_01preprocess & load variables
for listNum = 1:2;
    switch listNum
        case 1
            listType = 'test';
            disp('calculating and saving measures for test data\n')
            load([filePrefix 'figures/figsForMayPaper/data/TESTgraphingdata.mat']) %'TESTlist','TCs','filePrefix'
            list = TESTlist;
            clear TESTlist filePrefix
        case 2
            listType = 'control';
            disp('calculating and saving measures for control data\n')
            load([filePrefix 'figures/figsForMayPaper/data/CONTROLgraphingdata.mat']) %'CONTROLlist','TCs','filePrefix'
            list = CONTROLlist;
            clear CONTROLlist filePrefix
    end
    filePrefix = setPrefix(onKirkwood);
    
    
    % ************************************************************************
    nSubj = length(list);
    for iSubj = 1:nSubj;
        d = TCs{iSubj}; % has 4 fields: TTAxCW, TTAyCW, stackedXcw, stackedYcw
        
        xTS = d.stackedXcw;
        yTS = d.stackedYcw;
        xTTA = d.TTAxCW;
        yTTA = d.TTAyCW;
        clear d
        
        xyTS = xTS + yTS;
        pcData = (length(find(~isnan(xTS(:)))))/length(xTS(:));
        m.percentData(iSubj) = round(pcData*100);
        clear pcData
        
        numMeas = length(xTS(:));
        numCycles = size(xTS,2);
        numMperC = numMeas/numCycles;
        numPartsPerCycle = 4; % UL->UR, UR->LR, LR->LL, LL->UL
        numMperP = numMperC/numPartsPerCycle;
        
        %----------------------------------------%
        %  MEASURE 1: VARIANCE ACROSS 5 REPEATES %
%_______________________________________%
        
        maxwindow = 500;
        
        xVarVec = nanvar(xTS')';
        yVarVec = nanvar(yTS')';
        xyVarVec = nanvar(xTS'+yTS');
        
        m.xmeanVar(iSubj) = nanmean(xVarVec);
        m.ymeanVar(iSubj) = nanmean(yVarVec);
        m.xymeanVar(iSubj) =(((m.xmeanVar(iSubj)).^2) + ((m.ymeanVar(iSubj)).^2)).^(1/2);
        
%_______________________________________%
 
        xmax = nanmax(xVarVec);
        xmaxind = find(xVarVec == xmax);
        if(length(xmaxind)>1)
            xmaxind = xmaxind(end/2);
        end
            
        firstEl = max(1,xmaxind-maxwindow);
        lastEl = min(length(xVarVec),xmaxind+maxwindow);
        xmaxrange = firstEl:lastEl;
        xmaxvals = xVarVec(xmaxrange);
        m.maxxVar(iSubj) = nanmean(xmaxvals);
        
        %_______________________________________%
        
        ymax = nanmax(yVarVec);
        ymaxind = find(yVarVec == ymax);
        if(length(ymaxind)>1)
            ymaxind = ymaxind(floor(end/2));
        end
            
        firstEl = max(1,ymaxind-maxwindow);
        lastEl = min(length(yVarVec),ymaxind+maxwindow);
        ymaxrange = firstEl:lastEl;
        ymaxvals = yVarVec(ymaxrange);
        m.maxyVar(iSubj) = nanmean(ymaxvals);
        %_______________________________________%
         xymax = nanmax(xyVarVec);
        xymaxind = find(xyVarVec == xymax);
        if(length(xymaxind)>1)
            xymaxind = xymaxind(floor(end/2));
        end
            
        firstEl = max(1,xymaxind-maxwindow);
        lastEl = min(length(xyVarVec),xymaxind+maxwindow);
        maxrange = firstEl:lastEl;
        xymaxvals = xyVarVec(maxrange);
        m.maxxyVar(iSubj) = nanmean(xymaxvals);
        
%_______________________________________%
        
        %----------------------------------------%
        %  MEASURE 2: BOXY-NESS %
        
        xModelpart1 = linspace(-2,2,numMperP);
        xModelpart2 = repmat(2,1,numMperP);
        xModelpart3 = linspace(2,-2,numMperP);
        xModelpart4 = repmat(-2,1,numMperP);
        xModelTTA = [xModelpart1 xModelpart2 xModelpart3 xModelpart4];
        
        yModelpart1 = repmat(-2,1,numMperP);
        yModelpart2 = linspace(-2,2,numMperP);
        yModelpart3 = repmat(2,1,numMperP);
        yModelpart4 = linspace(2,-2,numMperP);
        yModelTTA = [yModelpart1 yModelpart2 yModelpart3 yModelpart4];
        
        xModel = repmat(xModelTTA,1,numCycles);
        yModel = repmat(yModelTTA,1,numCycles);

        
        xErr = xModel - xTS(:)';
        yErr = yModel - yTS(:)';
        xyErr = (xModel+yModel) - (xTS(:)'+yTS(:)');
        
        xErrTTA = xModelTTA' - xTTA;
        yErrTTA = yModelTTA' - yTTA;
        xyErrTTA = (xModelTTA' + yModelTTA') - (xTTA + yTTA);
        
        xyTTA = xTTA+yTTA;
        
        m.xTTAintegrity(iSubj) = 1 - sqrt(nansum(xErrTTA.^2))/sqrt(nansum(xTTA.^2));
        m.yTTAintegrity(iSubj) = 1 - sqrt(nansum(yErrTTA.^2))/sqrt(nansum(yTTA.^2));
        m.xyTTAintegrity(iSubj) = 1-sqrt(nansum(xyErrTTA.^2))/sqrt(nansum(xyTTA.^2));
        
        m.xR2TTA(iSubj) = 1-(nanvar(xErrTTA)/nanvar(xTTA));
        m.yR2TTA(iSubj) = 1-(nanvar(yErrTTA)/nanvar(yTTA));
        m.xyR2TTA(iSubj) = 1 - (nanvar(xyErrTTA)/nanvar(xyTTA));
        
        m.xR2(iSubj) = 1-(nanvar(xErr)/nanvar(xTS(:)));
        m.yR2(iSubj) = 1-(nanvar(yErr)/nanvar(yTS(:)));
        m.xyR2(iSubj) = 1-(nanvar(xyErr)/nanvar(xTS(:)+yTS(:)));
    end
    switch listNum
        case 1
            save([filePrefix 'figures/figsForMayPaper/data/TESTmeasures.mat'],'list','m')
        case 2
            save([filePrefix 'figures/figsForMayPaper/data/CONTROLmeasures.mat'],'list','m')
    end
 
    dataMatrix(:,1) = m.xmeanVar';
    dataMatrix(:,2) = m.ymeanVar';
    dataMatrix(:,3) = m.xymeanVar';
    
    dataMatrix(:,4) = m.maxxVar';
    dataMatrix(:,5) = m.maxyVar';
    dataMatrix(:,6) = m.maxxyVar';
    
    dataMatrix(:,7) = m.xTTAintegrity';
    dataMatrix(:,8) = m.yTTAintegrity';
    dataMatrix(:,9) = m.xyTTAintegrity';
    
    dataMatrix(:,10) = m.xR2TTA';
    dataMatrix(:,11) = m.yR2TTA';
    dataMatrix(:,12) = m.xyR2TTA';
    
    dataMatrix(:,13) = m.xR2';
    dataMatrix(:,14) = m.yR2';
    dataMatrix(:,15) = m.xyR2';
            
    xlswrite([filePrefix 'figures/figsForMayPaper/data/' listType 'measureMatrix.xls'],dataMatrix)
    
    clear TCs m dataMatrix
end

    
    % ************************************************************************
 %**** for list = 1:2 ********* %
% ************************************************************************


%*******************************************%
%**        **SUBFUNCTIONS **              **%
%*******************************************%

%___________________________________________%
function filePrefix = setPrefix(onKirkwood)
if(onKirkwood)
    filePrefix = '/Volumes/ShaniWSQBackupHD/Dropbox/';
else
    filePrefix = '~/Dropbox/';
end




