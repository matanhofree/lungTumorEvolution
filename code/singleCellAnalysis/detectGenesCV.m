function outCV = detectGenesCV(xData,batchID,CVin,inOpts)

    defaultOpts.doPlot = 0;
    defaultOpts.power = 1; 
    defaultOpts.resampleN = 100;
    defaultOpts.resampleFrac = 0.75;
    defaultOpts.doZIP = 1;
    defaultOpts.useNormalizedOverDisp = 1;
    
    % Normalized over disp params
    defaultOpts.nBins = 30; 
    defaultOpts.dispNormThr = 0.9;
    
    % CV 2 mean deviation (older version)   
    % defaultOpts.windowSize = 200;
    defaultOpts.windowSize = 100;
    defaultOpts.windowStep = 50;
    
    defaultOpts.outlierMethod = 2; 
    % 1 - zscore (x - median)/iqr > factor
    % 2 - outlier ( x > 75 + factor * iqr 
    % 3 - zscore simple (x - mean)/std > factor
    % 
    defaultOpts.rfactor = 2;
    
    
    defaultOpts.doGlobal = 0;
    defaultOpts.globalPropSample = 0;
    
    defaultOpts.minDiffFreq = 0.1;
    defaultOpts.diffThr = 1;
    defaultOpts.minGroup = 50;
    
    defaultOpts.keepSubSetData = 0;
    defaultOpts.minExp = 0;
    
    defaultOpts.ncores = 1;

    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    disp(opts);
    clear defaultOpts;
    
    if opts.ncores > 1 && isempty(gcp('nocreate')) && ~ismac()
        % distPool(opts.ncores);
        parpool(opts.ncores);
    end
    
    
    [D,N] = size(xData);    
    if ~exist('batchID','var') || isempty(batchID)
        batchID = cell(N,1);
        batchID(:) = { 'global' };
        
        if N < opts.minGroup            
            fprintf(1,'Minimum group size is set to 0');
            opts.minGroup = 0;            
        end
        
        
    end
        
    [batchNames,~,~,bCnt,bPos] = fastUnique(batchID);
    
    dropI = bCnt < opts.minGroup;
    batchNames(dropI) = [];
    bPos(dropI) = [];
    bCnt(dropI) = [];    
    
    
%     if (opts.doGlobal)
%         batchNamesList =  [ setdiff(batchNames,'global'); 'global'];
%     else
%         batchNamesList = unique(batchNames);
%     end
    
    if ~exist('CVin','var') 
        CVin = [];
    end
   
    batchN = length(batchNames);    
    
    zbt = tic();
    if opts.ncores > 1 && isempty(CVin)
        
        xDataSlice = cell(batchN,1);
        for i = 1:batchN
            cRunData{i} = xData(:,bPos{i});             
        end
        
        
        parfor i = 1:batchN
            
            cBatch = batchNames{i};
            
            fprintf('Running %d/%d -- %s\n',i,batchN,cBatch)
            zbt = tic();
            
%             if (strcmp(cBatch,'global'))
%                 selectV = true(1,N);
%             else
%                 zi = strcmp(batchNames,cBatch);
%                 selectV = trueV(bPos{zi},N);            
%             end         
% 
%             cvdata = [];
%             if ~isempty(CVin) && isfield(CVin,cBatch)
%                 fprintf('Using exising cvdata for %s\n',cBatch);
%                 cvdata = CVin.(cBatch);            
%             end
%             cRunData = xData(:,selectV);   
            
            subSetOutCV{i} = calcDisp(cRunData{i},[],opts);
            % typeID{i} = cBatch;
            toc(zbt);

        end    
        clear cRunData;
        
    else
        for i = 1:batchN
            cBatch = batchNames{i};
            
            fprintf('Running %d/%d -- %s\n',i,batchN,cBatch)
            zbt = tic();
            
            if (strcmp(cBatch,'global'))
                selectV = true(1,N);
            else
                zi = strcmp(batchNames,cBatch);
                selectV = trueV(bPos{zi},N);            
            end         

            cvdata = [];
            if ~isempty(CVin) && isfield(CVin,cBatch)
                fprintf('Using exising cvdata for %s\n',cBatch);
                cvdata = CVin.(cBatch);            
            end
            cRunData = xData(:,selectV);   
            
            subSetOutCV{i} = calcDisp(cRunData,cvdata,opts);
            % typeID{i} = cBatch;
            toc(zbt);
        end
    end
    
    %% Fix parallel reformat    
    outCV.typeID = batchNames;
    outCV.countCV = nan(D,batchN);
    outCV.geneNnzFreq = nan(D,batchN);
    outCV.medianDispNorm = nan(D,batchN);
    
    for i = 1:batchN    
        cBatch = batchNames{i};
        
        outCV.countCV(:,i) = subSetOutCV{i}.countCV;
        outCV.geneNnzFreq(:,i) = subSetOutCV{i}.geneNnzFreq;
        outCV.medianDispNorm(:,i) = subSetOutCV{i}.medianDispNorm;
        
        if opts.keepSubSetData == 1
            outCV.(cBatch) = subSetOutCV{i}.subSetData;
        end
        
    end
    %%
    diffThr = max(opts.minDiffFreq*opts.resampleN,1);
    outCV.isDiff = any(outCV.countCV>=diffThr,2);
    
    
    
end
    

function outCV = calcDisp(cRunData,cvdata,opts)

    if opts.useNormalizedOverDisp            
        [zCV.countCV,zCV.dataMean,zCV.dataStd,zCV.piMom,zCV.dispNorm] = detectGenesOverdispersion(cRunData,cvdata,opts);
        outCV.countCV = zCV.countCV;
        outCV.medianDispNorm = nanmedian(zCV.dispNorm,2);            
    else
        [zCV.countCV,zCV.dataMean,zCV.dataStd,zCV.piMom] = detectGenesCV2mean(cRunData,cvdata,opts);
        outCV.countCV = sum(zCV.countCV>=opts.diffThr,2);
    end

    if opts.keepSubSetData 
        outCV.subSetData = zCV;
    end
    
    outCV.geneNnzFreq = sum(cRunData>opts.minExp,2)/size(cRunData,2);

end
