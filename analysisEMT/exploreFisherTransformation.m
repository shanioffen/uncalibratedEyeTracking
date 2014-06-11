function exploreFisherTransformation

% I want to make sure I understand the details of the transformation.

% step 1: create a null distribution to mimic our EMT timecourses and
% correlation to a sinusoid

countInd = 0;
for vectorLength = [100 500 1000];
    vectorIndex(:,1) = 1:vectorLength;
    countInd = countInd+1;

    sampleSize = 100;
    
    data = rand(vectorLength,sampleSize+1);
    rvals = NaN*ones(sampleSize,1);
    
    cormat = corr(data);
    corvals = abs(cormat(1,2:end));
    
    % do analysis just like for the paper, to double check
    
    for iVec = 1:sampleSize
        
        
        [f, g] = fit(vectorIndex,data(:,iVec),'sin1');
        
        
        rsquarevals(iVec) = g.rsquare;
        rvals(iVec) = sqrt(abs(g.rsquare));
        
        clear f g
        
        
    end
    
%     figure(vectorLength)
%     subplot(221)
%     hist(rvals)
%     title('r vals')
%     subplot(222)
%     hist(rsquarevals)
%     title('rsquare vals')
%     subplot(223)
%     hist(corvals)
%     title('correlation vals')
%     subplot(224)
%     hist(sqrt(corvals.^2))
%     title('sqrt of squares of corr vals')
    
      
    % do fisher transformation on r values:
    
    fisherValsR = (.5)*log((1+rvals)./(1-rvals));
    fisherValsCorr = (.5)*log((1+corvals)./(1-corvals));
    
    
    
    % now compare means and std for the 4 distributions (r, corr, fisherR,
    % fisherCorr) and check on %age inside the std's
    
    rMeans(countInd) = mean(rvals);
    corMeans(countInd) = mean(corvals);
    fisherRmeans(countInd) = mean(fisherValsR);
    fisherCormeans(countInd) = mean(fisherValsCorr);
    
    rSTD(countInd) = std(rvals);
    corSTD(countInd) = std(corvals);
    fishRstd(countInd) = std(fisherValsR);
    fishCorstd(countInd) = std(fisherValsCorr);
    
    
    % figure out %-iles
    % within 1 STD
    temp = length(find(and(corvals<=(corMeans(countInd)+corSTD(countInd)),...
                           corvals>=(corMeans(countInd)-corSTD(countInd)))));
    z1PercentCorvals(countInd) = temp/length(corvals);
    clear temp
    
    temp = length(find(and(rvals<=(rMeans(countInd)+rSTD(countInd)),...
                           rvals>=(rMeans(countInd)-rSTD(countInd)))));
    z1PercentRvals(countInd) = temp/length(rvals);
    clear temp
    
    temp = length(find(and(fisherValsCorr<=(fisherCormeans(countInd)+fishCorstd(countInd)),...
                           fisherValsCorr>=(fisherCormeans(countInd)-fishCorstd(countInd)))));
    z1PercentFishCorvals(countInd) = temp/length(fisherValsCorr);
    clear temp
    
    temp = length(find(and(fisherValsR<=(fisherRmeans(countInd)+fishRstd(countInd)),...
                           fisherValsR>=(fisherRmeans(countInd)-fishRstd(countInd)))));
    z1PercentFishRvals(countInd) = temp/length(fisherValsR);
    clear temp
    
    
    % within 2 std
    temp = length(find(and(corvals<=(corMeans(countInd)+2*corSTD(countInd)),...
                           corvals>=(corMeans(countInd)-2*corSTD(countInd)))));
    z2PercentCorvals(countInd) = temp/length(corvals);
    clear temp
    
    temp = length(find(and(rvals<=(rMeans(countInd)+2*rSTD(countInd)),...
                           rvals>=(rMeans(countInd)-2*rSTD(countInd)))));
    z2PercentRvals(countInd) = temp/length(rvals);
    clear temp
    
    temp = length(find(and(fisherValsCorr<=(fisherCormeans(countInd)+2*fishCorstd(countInd)),...
                           fisherValsCorr>=(fisherCormeans(countInd)-2*fishCorstd(countInd)))));
    z2PercentFishCorvals(countInd) = temp/length(fisherValsCorr);
    clear temp
    
    temp = length(find(and(fisherValsR<=(fisherRmeans(countInd)+2*fishRstd(countInd)),...
                           fisherValsR>=(fisherRmeans(countInd)-2*fishRstd(countInd)))));
    z2PercentFishRvals(countInd) = temp/length(fisherValsR);
    clear temp
    
    
    
    
    
    figure(vectorLength+1)
    subplot(221)
    hist(fisherValsCorr)
    title(['Fisher Corr z1= ' num2str(z1PercentFishCorvals(countInd))...
                      ' z2= ' num2str(z2PercentFishCorvals(countInd))])
    subplot(222)
    hist(fisherValsR)
    title(['Fisher R  z1= ' num2str(z1PercentFishRvals(countInd))...
                   ' z2= ' num2str(z2PercentFishRvals(countInd))])
    subplot(223)
    hist(corvals)
    title(['original corr  z1= ' num2str(z1PercentCorvals(countInd))...
                        ' z2= ' num2str(z2PercentCorvals(countInd))])
    subplot(224)
    hist(rvals)
    title(['original R   z1= ' num2str(z1PercentRvals(countInd))...
                      ' z2= ' num2str(z2PercentRvals(countInd))])
    
    
  
    
    
    clear data vectorIndex 

end

meansMat = [rMeans' corMeans' fisherRmeans' fisherCormeans']
stdMat = [rSTD' corSTD' fishRstd' fishCorstd']

keyboard









