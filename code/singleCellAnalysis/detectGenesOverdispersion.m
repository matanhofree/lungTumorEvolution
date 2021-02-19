function [countCV,dataMean,dataStd,piMom,dispNorm] = detectGenesOverdispersion(xData,inCV,inOpts)

    defaultOpts.resampleN = 20;
    defaultOpts.resampleFrac = 0.75;
    
    defaultOpts.doZIP = 1;      
    defaultOpts.nBins = 20; 
    defaultOpts.dispNormThr = 0.9;
    
    defaultOpts.minMeanExp = 0;
    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
        
    [D,N] = size(xData);        
    runN = opts.resampleN;            
    piMom = 0;
    dispNormNnz = 0;
    smpN = floor(opts.resampleFrac(1)*N);
    if isempty(inCV)
        dM = zeros(D,runN);
        dStd  = zeros(D,runN);

        zT = tic;        
        progressbar();
        for i = 1:runN
            if (mod(i,10) == 0)
                fprintf('Subsample %d (time %d)\n',i,toc(zT));zT = tic;
            end

            smpCols = randsample(N,smpN);
            xDataSmp = full(xData(:,smpCols)'); 

            % zTot =  mean(xDataSmp);
            
            dM(:,i) = mean(xDataSmp);
            % zvarE = sqrt(sum((xDataSmp - zmean).^2)/(smpN-1));
            % zvarE = std(xDataSmp);
            try 
                dStd(:,i) = std(xDataSmp);
            catch 
                dStd(:,i) = (sum(xDataSmp.^2) - dM(:,i))/(smpN-1);
            end
                    
            progressbar(i/runN);
        end
        %%
        if opts.doZIP == 1       
            lambdaMom = (dStd.^2 + dM.^2 - dM)./dM;
            piMom = (dStd.^2 - dM)./(dStd.^2 + dM.^2 - dM); 

            dataMean = (1-piMom).*lambdaMom;
            dataStd = sqrt(dataMean.*(1+lambdaMom.*piMom));
        else
            dataMean = dM;
            dataStd = dStd;
        end
    
        meanCV = (dataStd)./dataMean;        
        gMean = nanmean(dataMean,2);
        
        dropMinVal = isnan(gMean) | gMean <= opts.minMeanExp;
        
        meanCV(dropMinVal,:) = [];
        gMean(dropMinVal) = [];
        
        
%         
%        binStep = 1/opts.nBins;
%        binEdges = unique([ quantile(nonzeros(gMean),binStep:binStep:(1-binStep)) inf]);    
%        binEdges = setdiff(binEdges,0);
%
%         [zCnt,~,smpBins] = histcounts(gMean,binEdges);
% %%
%         dBinDispMedian = cell2mat(arrayfun(@(x)nanmedian(meanCV(smpBins == x,:)),full(unique(smpBins)),'uniformoutput',0));
%         
%         % calculate normalized dispersion
%         dBinDispDev = meanCV - dBinDispMedian(smpBins+1,:);
%         dBinDispAbsDev = abs(dBinDispDev);
% 
%         dBinDispMAD = cell2mat(arrayfun(@(x)nanmedian(dBinDispAbsDev(smpBins == x)),full(unique(smpBins)),'uniformoutput',0));
% 
% 
%         dispNorm = dBinDispDev ./ dBinDispMAD(smpBins+1,:);
%%
       binStep = 1/opts.nBins;
       binEdges = [0 quantile(nonzeros(gMean),binStep:binStep:1) ];
       % binEdges = setdiff(binEdges,0);

       [zCnt,~,smpBins] = histcounts(gMean,binEdges);             
       
       [smpBinsList,pIdx,pRef,zCntBin,binPos] = fastUnique(smpBins);
       
%%
        dBinDispMedian = cell2mat(cellfun(@(x)nanmedian(meanCV(x,:),1),binPos,'uniformoutput',0)');
        
%%      % calculate normalized dispersion
        dBinDispDev = meanCV - dBinDispMedian(pRef,:);
        dBinDispAbsDev = abs(dBinDispDev);
%%       
        dBinDispMAD = cell2mat(cellfun(@(x)nanmedian(dBinDispAbsDev(x,:),1),binPos,'uniformoutput',0)');

        dispNormNnz = dBinDispDev ./ dBinDispMAD(pRef,:);


    else
        
        error('Not implemented');
        dataMean = inCV.dataMean;
        dataStd = inCV.dataStd;
        piMom = inCV.piMom;
        runN = size(dataMean,2);
        meanCV = inCV.dataStd./inCV.dataMean;
        dispNormNnz = inCV.dispNorm;
    end
    
    dispNormThr = quantile(dispNormNnz,opts.dispNormThr);    
    countCVnnz = sum(bsxfun(@gt,dispNormNnz,dispNormThr),2);
    
    countCV = nan(D,1);
    countCV(~dropMinVal) = countCVnnz;
        
    dispNorm = nan(D,runN);
    dispNorm(~dropMinVal,:) = dispNormNnz;
    % Plot: 
    % zSel = countCVnnz>10; figure; hold on; scatter(gMean,nanmean(meanCV,2));scatter(gMean(zSel),nanmean(meanCV(zSel,:),2));
    % set(gca,'xscale','log','yscale','log')
end
    


