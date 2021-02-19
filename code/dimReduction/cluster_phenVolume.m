function [phenV,clustList,qData,isZero] = cluster_phenVolume(cData,clustV,inOpts)
    
     
    defaultOpts.useShrink = 1;
    defaultOpts.mahalOnly = 1;
    defaultOpts.sampleNum = 200;
    defaultOpts.sampleReplace = 0;
    defaultOpts.repNum = 20;
    defaultOpts.subSmpFrac = 0.8;
    
    defaultOpts.minSampleEdges = [];
    
    defaultOpts.minCl = 30;
    defaultOpts.dType = 1;
    defaultOpts.dist = 'euclid';
    defaultOpts.mFunc = @(x)nanmean(x(:));
    
    defaultOpts.nPC = round(defaultOpts.sampleNum/2);
    
    defaultOpts.discBins = 10;
    defaultOpts.useNMI = 1;
    
    defaultOpts.useCdfByHist = 1;
    defaultOpts.addNnz = 1;
    
    defaultOpts.pb = 1;
               
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    
    if opts.pb == 1 && ~ismac()
        opts.pb = 2;        
    end
          
    if ~isempty(opts.minCl) 
       
        [bTypes,~,~,bCnt,bPos] = fastUnique(clustV);
        
        cDrop = bCnt < opts.minCl;
        if sum(cDrop) > 0
            
            fprintf('Filtering the following batches dues to inssufficient cells in channel:'); 
            disp(bTypes(cDrop));
            
            selectFilter = trueV([ bPos{~cDrop} ],length(clustV));            
            tabFilter(clustV,selectFilter);
            
            cData = cData(:,selectFilter);
            clustV = clustV(selectFilter);
        end        
    end
    
   if ~isempty(opts.minSampleEdges)
        mSmp = round(nchoosek(round(opts.minCl*opts.subSmpFrac),2)*opts.minSampleEdges);
        fprintf('Sampling the %d edges from each rep sample\n',mSmp);
    end
    
    dType = opts.dType;
    
    if dType == 4
        if opts.useCdfByHist ==  1
           fprintf('CDF discritizing data into %d bins by histogram method\n',opts.discBins);
           optsT.bkgNumBins = opts.discBins
           qData = expToQuant(cData,optsT);
        else
           fprintf('CDF discritizing data into %d bins by equal quantiles\n',opts.discBins);
           optsT.numBins = opts.discBins;
           optsT.addNnz = opts.addNnz;
           [qData,~,isZero] = mapToQuantileFix(cData,[],optsT);
           
           if any(isZero)
               fprintf('Dropping all zero values %d\n',sum(isZero))
               qData = qData(~isZero,:);
               cData = cData(~isZero,:);
           end
           assert(sum(all(qData==0,2)) == 0);
        end
        fprintf('Done disc nnz %d - %f\n',nnz(qData),nnz(qData)/numel(qData));
    end
      
    
 	[D,N] = size(cData);
    [clustList,~,~,~,clPos] = fastUnique(clustV);
    nG = length(clustList);    

    nRep = opts.repNum;
    phenV = nan(nRep,nG);
    

    nT = nG*nRep;
    zp = 1;
    
    if opts.pb == 1
        progressbar();
    elseif opts.pb == 2
        [~,parForMonFile] = parfor_progress_mod(nT,[],10);
    end
    
    for i = 1:nG
        
        ciFull = clPos{i};               
        
        for j = 1:nRep
        
            if opts.sampleReplace 
                ci = randsample(ciFull,opts.sampleNum,1);
            else 
                np = round(length(ciFull)*opts.subSmpFrac);
                if np > opts.sampleNum
                    ci = randsample(ciFull,opts.sampleNum);
                else
                    if j == 1
                        fprintf('Note: insufficient samples for %s -- sampling %d instead of %d\n',clustList{i},np,opts.sampleNum);
                    end
                    
                    ci = randsample(ciFull,np);
                end
                
            end
            
            switch dType                                                        
                case 1
                    iCov = est_cov(cData(:,ci)',opts.useShrink);
                    covS = svd(iCov);
                    %%
                    logCovS = log(covS);
                    logCovS(logCovS<0) = [];
                    %

                    phenV(j,i) = sum(logCovS);
                case 2

                    zD = pdist(cData(:,ci)',opts.dist);           

                    phenV(j,i) = opts.mFunc(zD);
                case 3

                    [~, zPC] = pca(full(cData(:,ci)'),'NumComponents',opts.nPC);
                    zD = pdist(zPC,opts.dist);           

                    phenV(j,i) = opts.mFunc(zD);
                    
                case 4 
                    
                    cX = full(qData(:,ci));                       
                    zD = miAllPairs(cX);                                             
                    
                    if opts.useNMI 
                        %% 
                        zH = sqrt(arrayfun(@(x)h(cX(:,x)),1:size(cX,2)));
                        % zH = (zH);
                        %%

                        zD = zD./zH;
                        zD = zD./(zH');
                                        
                    end
                    
                                            
                    zD(isnan(zD)) = [];
                    
                    if ~isempty(opts.minSampleEdges)
                        zD = zD(randsample(length(zD),mSmp));
                    end
                
                    phenV(j,i) = opts.mFunc(zD);                    

                otherwise 
                    error('Distance type not found');

            
            end
            
            
            if opts.pb == 1                
                progressbar(zp/nT); zp = zp + 1;
            elseif opts.pb == 2
                parfor_progress_mod(-1,parForMonFile);
            end
           
        end
        
    end


end

function iCov = est_cov(cData,covType)
    
    
    switch covType
        case 1
            iCov = shrinkage_cov(full(cData));                        
        case 2 
            iCov = shrinkage_cov(full(cData),'rblw');                        
        case 3 
            [n,p] = size(cData);
            if n < 2*p
                iCov = shrinkage_cov(full(cData));                        
            else
                iCov = robustcov(full(cData));
            end
        otherwise
            iCov = full(cov(cData));            
    end
        
end