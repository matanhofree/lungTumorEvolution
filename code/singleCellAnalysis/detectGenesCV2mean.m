function [countCV,dataMean,dataStd,piMom] = detectGenesCV2mean(xData,inCV,inOpts)

    defaultOpts.doPlot = 1;
    defaultOpts.power = 1; 
    defaultOpts.resampleN = 20;
    defaultOpts.resampleFrac = 0.75;
    defaultOpts.doZIP = 1;
    % defaultOpts.windowSize = 200;
    defaultOpts.windowSize = 100;
    defaultOpts.windowStep = 50;
    defaultOpts.outlierMethod = 2; 
    % 1 - zscore (x - median)/iqr > factor
    % 2 - outlier ( x > 75 + factor * iqr 
    % 3 - zscore simple (x - mean)/std > factor
    % 
    defaultOpts.rfactor = 2;    
    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
        
    [D,N] = size(xData);        
    runN = opts.resampleN;            
    
    if isempty(inCV)
        dM = zeros(D,runN);
        dStd  = zeros(D,runN);

        zT = tic;        
        progressbar();
        for i = 1:runN
            if (mod(i,10) == 0)
                fprintf('Subsample %d (time %d)\n',i,toc(zT));zT = tic;
            end

            smpCols = randsample(N,floor(opts.resampleFrac(1)*N));
            xDataSmp = xData(:,smpCols); 

            dM(:,i) = nanmean(xDataSmp,2);
            dStd(:,i) = nanstd(xDataSmp,0,2);        
            progressbar(i/runN);
        end

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
    else
        dataMean = inCV.dataMean;
        dataStd = inCV.dataStd;
        piMom = inCV.piMom;
        runN = size(dataMean,2);
        meanCV = inCV.dataStd./inCV.dataMean;   
    end
    countCV = zeros(D,runN);
    
    progressbar();
    %% 
    for i = 1:runN
       [meanV,meanVIdx] = sort(dataMean(:,i));
       meanVIdx(isnan(meanV)) = [];
       meanV(isnan(meanV)) = [];       
       curCV = zeros(D,1);
       countC = zeros(D,1);
       
       
       zj = 1;
       
       %%
       while (zj < length(meanVIdx))
           %
           zpos = max(zj - opts.windowSize,1):min(zj+opts.windowSize-1,length(meanVIdx)); 
           zidx = meanVIdx(zpos);
           
           switch opts.outlierMethod
               case 1
                    statMedian = quantile(meanCV(zidx,i),[ 0.25 0.5 0.75 ]);
                    ratioR = (meanCV(zidx,i) - statMedian(2))./(statMedian(3) - statMedian(1));
                    curCV(zidx) = nansum([ curCV(zidx) ratioR ],2);                    
                    
                    countC(zidx) = countC(zidx) + ~isnan(ratioR);    
               case 2 
                    statMedian = quantile(meanCV(zidx,i),[ 0.25 0.5 0.75 ]);                    
                    curCV(zidx) = curCV(zidx) + (meanCV(zidx,i) > (statMedian(3) + opts.rfactor*(statMedian(3) - statMedian(1))));
                    countC(zidx) = countC(zidx) + 1; 
               case 3 
                    statMean = nanmean(meanCV(zidx,i));
                    statStd = nanstd(meanCV(zidx,i));
                    curCV(zidx) = curCV(zidx) + (meanCV(zidx,i) - statMean)./statStd;
                    countC(zidx) = countC(zidx) + 1;                    
               otherwise
                   error('Outlier method does not exist');
                    
           end
           zj = zj + opts.windowStep;
           
       end
       %%
       
       switch opts.outlierMethod
               case 1
                    countCV(:,i) = curCV./countC > opts.rfactor;
               case 2 
                    countCV(:,i) = curCV > 0;
               case 3
                    countCV(:,i) = curCV./countC > opts.rfactor;
               otherwise
                   error('Outlier method does not exist');

       end

       progressbar(i/runN);
    end
    
    %%
    estMean = nanmean(dataMean,2);
    mCV = mean(meanCV,2);
    zz = sum(countCV,2)>1;
    
    if opts.doPlot
        figure()
        hold on
        scatter(estMean(~zz),mCV(~zz)); 
        scatter(estMean(zz),mCV(zz),'k'); 
        set(gca,'yscale','log','xscale','log')
        
        text(1,1,sprintf('Num diff = %d',sum(zz)));
    else
        fprintf('Num diff = %d',sum(zz));
    end

%%
% Use smooth 
%     %%
%     for i = 1:runN
%         zT = tic;
%         smoothCV(:,i) = smooth(dataMean(:,i),meanCV(:,i),opts.windowSize, 'rlowess');        
%         fprintf('Subsample %d: ',i);toc(zT);
%     end
%     %%
%     %%
%     tic
%     smoothAlt = lowess_matrix(dataMean(:,1),meanCV(:,1), opts.windowSize, 'lowess',1, 5);
%     toc
%     %%
%     
%     estCV = nanmean(smoothCV,2);
%     estMean = nanmean(dataMean,2);
%     mCV = mean(meanCV,2);
%     
%     %%
%     if opts.doPlot
%         figure()
%         [zz,zidx] = sort(estMean);
%         zidx(isnan(zz)) = [];
%         zz(isnan(zz)) = [];
%         
%         hold on
%     
%         scatter(estMean(zidx(1:end-1)),mCV(zidx(1:end-1)));
%         plot(estMean(zidx(1:end-1)),estCV(zidx(1:end-1)),'-k');
%         set(gca,'yscale','log','xscale','log')
%     end
   

end
    
    
%% 

