function [gMean,gStd,meanCV,dM,dStd] = zip_moment(xData,dim,inOpts)
    
    defaultOpts.resampleN = 20;
    defaultOpts.resampleFrac = 0.75;
    defaultOpts.doZIP = 1;

        
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
        
    if nargin < 2 
        dim = 2;
    elseif dim == 1        
        xData = xData';
    end
    
    [D,N] = size(xData);        
    
    runN = opts.resampleN;            
    piMom = 0;    
    
    smpN = floor(opts.resampleFrac(1)*N);

    dM = zeros(D,runN);
    % dStd  = zeros(D,runN);
    dVar  = zeros(D,runN);

    zT = tic;        
    progressbar();
    for i = 1:runN
        if (mod(i,10) == 0)
            fprintf('Subsample %d (time %d)\n',i,toc(zT));zT = tic;
        end

        smpCols = randsample(N,smpN);
        xDataSmp = full(xData(:,smpCols)'); 

        dM(:,i) = mean(xDataSmp);

        % dStd(:,i) = std(xDataSmp);
        dVar(:,i) = var(xDataSmp,1);
        progressbar(i/runN);
    end
    %%
    if opts.doZIP == 1       
        lambdaMom = (dVar + dM.^2 - dM)./dM;
        % lambdaMom = (dVar+ dM.^2)./dM - 1;
        piMom = (dVar - dM)./(dVar + dM.^2 - dM); 

        dataMean = (1-piMom).*lambdaMom;
        dataStd = sqrt(dataMean.*(1+lambdaMom.*piMom));
    else
        dataMean = dM;
        dataStd = sqrt(dVar);
    end
    
    dStd = sqrt(dVar);
    meanCV = (dataStd)./dataMean;        
    gMean = nanmean(dataMean,2);
    gStd = nanmean(dataStd,2);


end