function [freqO,expO,isNanRow] = quantMatchExpByBinMultiNoref(cDataIn,batchID,inOpts)
 
    defaultOpts.geneVar = 'geneID';
    defaultOpts.dataMat = 'normTPM';
        
    defaultOpts.outExp = 1;
    defaultOpts.outFreq = 0;
    defaultOpts.mergeExpVect = 1;
    defaultOpts.mergeExpGrid = 100;
    defaultOpts.sparseOut = 0; 
    defaultOpts.sparseOutput = 1; 
    
    defaultOpts.meanFunc = 0; % 0 Mean
                              % 1 Median
                              % 2 Quantile 0.75      
                              % 3 trimean
                              
    defaultOpts.qDim = 0.75;
    defaultOpts.normalizeZeroRate =  0;
    defaultOpts.centerFreqDist = 0;
    
    defaultOpts.minNumForFreq = 20;
    defaultOpts.removeNan = 1;
        
    defaultOpts.minRefVal = 4;
    defaultOpts.pb = 1;
        
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;  
    disp(opts)    
    
    switch opts.meanFunc
        case 0
            mFunc = @(x)nanmean(x,1);
        case 1
            mFunc = @(x)nanmedian(x,1);
        case 2 
            mFunc = @(x)quantile(x,opts.qDim,1);
        case 3 
            mFunc = @(x)trimean(x,1);
        otherwise
            mFunc = @(x)nanmean(x,1);
    end
          
    % geneA = cDataIn.(opts.geneVar); 
    if isstruct(cDataIn)        
        dataM = (cDataIn.(opts.dataMat));
    else
        dataM = cDataIn;
    end
 
    [zD,nA] = size(dataM);
    nnzA = nnz(dataM);

 
    if nA ~= length(batchID)
        error('Batchid does not match');
    end
    
    if opts.outExp 
        if opts.sparseOut            
            expO = sparse([],[],[],zD,nA,nnzA);
        else       
            expO = zeros(zD,nA);
        end  
    else
        expO = [];
    end
    
    if opts.outFreq
        if opts.sparseOut
            freqO = sparse([],[],[],zD,nA,nnzA);
        else
            freqO = zeros(zD,nA);        
        end  
    else
        freqO = [];
    end
    
    [batchList,~,nIdx,bCnt,bPos] = fastUnique(batchID);            
    bVect = cellfun(@(x)trueV(x,nA),bPos,'uniformoutput',0);
       
    nBatch = length(bVect);
    mergeGrid = linspace(0, 1, opts.mergeExpGrid);

    doExp = 0;
    doFreq  = opts.outFreq;
    if opts.outExp && nargout > 1
        doExp = 1;
        if opts.mergeExpVect
            mergeGridMat = zeros(nBatch,length(mergeGrid));
        end
    end
    doExpMerge = opts.mergeExpVect;       
    
    % progressbar();
    zt = tic();    
    if opts.pb == 1 && ismac()
        progressbar();
    elseif opts.pb == 2
        [~,parForMonFile] = parfor_progress(zD);
    end
    for zi = 1:zD        
        %%
        cFreq = zeros(1,nA);
        batchZeroFreq = zeros(nBatch,1);
        batchMeanFreq = zeros(nBatch,1);
        mergeGridMat = zeros(nBatch,length(mergeGrid));
        
        for zb = 1:nBatch
            bFull = full(dataM(zi,bVect{zb}));                         
            
            rankB = tiedrank(bFull);
            fB = rankB - min(rankB);
            fB = fB/max(max(fB),1);
            
            cFreq(1,bVect{zb}) = fB;
            % fBsplit{zb} = fB;
            
            if doExp && doExpMerge
                [eB,zidx] = unique(bFull);
                
                
                if length(zidx) < opts.minRefVal
                    % fprintf('Skiping %d) %d, meanE %f\n',zi,length(zidx),full(mean(dataM(zi,:))));
                    continue;
                end
                
               
                if length(eB) > opts.minNumForFreq
                    xfB = fB(zidx);
                    xB = bFull(zidx);  
                    
                    nnzFreq = nonzeros(xfB);
                    batchZeroFreq(zb) = min(nnzFreq);
                    batchMeanFreq(zb) = mean(nonzeros(fB));
                    mergeGridMat(zb,:) = max(interp1(xfB,xB,mergeGrid,'linear','extrap'),0); 
                else
                    if (eB ~= 0)
                        warning('Unique nnz expB (%d)',eB);
                    end

                    mergeGridMat(zb,:) = nan;
                end
            end
        end
        
        %%                        
        if doExp            
            if doExpMerge == 0
                error('Not implemented');
                % expO(zi,:) = interp1(xfA,xA,fA,'linear','extrap');                         
                %expO(zi,:) = interp1(xfA,xA,full(freqB(zi,:)),'linear','extrap');
            else
                if opts.normalizeZeroRate 
                    zLim = mergeGrid<batchZeroFreq;
                    mergeGridMat(zLim) = 0;
                end                
                mExp = mFunc(mergeGridMat);
                
                if opts.centerFreqDist  
                    %%
                    
                    centerCorrection = batchMeanFreq - median(batchMeanFreq);
                    cFreqCorr = cFreq - centerCorrection(nIdx)';
                    
                    cFreqCorr = max(min(cFreqCorr,1),0);
                    cFreqCorr(cFreq == 0) = 0;
                    
                    expO(zi,:) = interp1(mergeGrid,mExp,cFreqCorr,'linear','extrap');
                    
                    
                    % cExp = interp1(mergeGrid,mExp,cFreq,'linear','extrap');
                    
                    % cExpMean = cellfun(@(x)mean(cExp(x)),bVect);
                    % cExpMeanC =  cExpMean - mean(cExpMean) ;
                    
                    % expO(zi,:) = cExp - cExpMeanC(nIdx)';
                else
                    expO(zi,:) = interp1(mergeGrid,mExp,cFreq,'linear','extrap');
                end
                
                % expO(zi,:) = interp1(mergeGrid,mExp,fA,'linear','extrap');                         
                % expO(zi,:) = interp1(mergeGrid,mExp,cFreq,'linear','extrap');
            end           
        end
        if doFreq
            freqO(zi,:) = cFreq;
        end
        % progressbar(zi/zD);
        
        if opts.pb == 1 && ismac()
            progressbar(zi/zD);
        elseif opts.pb == 2
            parfor_progress(-1,parForMonFile);             
        elseif mod(zi,500) == 0
            fprintf('%d) %2f\n',zi,zi/zD*100);
            toc(zt);
            zt = tic(); 
        end
        
    end
    
    isNanRow = [];
    if doFreq
        if opts.sparseOutput == 1
            freqO = sparse(freqO);
        end
    end
    
    if doExp
        if opts.removeNan == 1
            isNanRow = false(zD,1);
            for i = 1:zD
                cnanV = isnan(expO(i,:));
                cN = sum(cnanV);                
                if cN > 0                    
                    fprintf('Note: nans found in gene %d. In-total %d nans (%f).\n',i,cN,cN/nA);  
                    isNanRow(i) = 1;
                    expO(i,cnanV) = 0;
                end
            end
        end
            
        if opts.sparseOutput == 1
            expO = sparse(expO);
        end
        
        if doFreq == 0 && nargout == 1
            fprint('Returning expression as output\n');
            freqO = expO;
        end
        nV = numel(dataM);
        fprintf('Done gene quantile matching. Normalized %d batches.\nNNZ before: %d (%f), after: %d (%f).\n',nBatch,nnzA,nnzA/nV,nnz(expO),nnz(expO)/nV)
    end
    
end

