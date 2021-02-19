function [outConsensusNMF,zfig] = extractConsensusNMF_graphClust(baseVect,cDataM,batchID,inOpts)

    defaultOpts.knnFilterRank = 20;
    
    defaultOpts.maxRankQ = 0.25;
    
    defaultOpts.doPlot = 1;
%    defaultOpts.testKRange = [ 3 30 ];
%    defaultOpts.testK = [];
%    defaultOpts.nrepSubClust = 5;    
    defaultOpts.chooseBestHueristic = 1; % 1 ReconstKneePt
                                         % 2 Mean Sil

                                         
    defaultOpts.maxGraphClustK = 50;
    defaultOpts.filterOutlierBatch = 1;
    
    defaultOpts.mergedDistType = 'cosine';
%     defaultOpts.wmiCluster = 0;
    
    defaultOpts.keepInt = 1;
    defaultOpts.fixW = 1;
    defaultOpts.nmfSolver = 'bp';
    
    defaultOpts.tRange = 10.^(-2:0.2:1);
    defaultOpts.Laplacian = 'Combinatorial';
    defaultOpts.maxClust = 500;
    defaultOpts.minClust = 2;
    defaultOpts.fullOptimizeNMF = 1;
    defaultOpts.bp_nmf = struct();
    
    defaultOpts.doAlphaKernelKnn = 0;
    defaultOpts.knnSimK = 30;
    defaultOpts.knnSimAlpha = 40;
    defaultOpts.knnDistFun = 'cosine';
    
    defaultOpts.knnExcludeHardRank = [];
    
    
    defaultOpts.maxResampleF = 0.8;
    defaultOpts.minSamplePerSet = 50;
    
        
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    disp(opts)
    %%
    
    if opts.fixW == 1
        [dH,dW] = size(baseVect{1});
        if dH < dW
            baseVect = cellfun(@(x)x',baseVect,'uniformoutput',0);
        end  
    end    
    
    mergeFull = cell2mat(baseVect);
    disp(size(mergeFull));
    
    mergeFullWeight = [];
    switch opts.mergedDistType
        case 'euclid'
            fprintf('Using euclidean distance\b');
            cDist = squareform(pdist(mergeFull','euclid'));
            [cDistSorted] = sort(cDist,'ascend');
        case 'cosine'            
            fprintf('Using cosine distance\b');
            cDist = squareform(pdist(mergeFull','cosine'));
            [cDistSorted] = sort(cDist,'ascend');
%            mergeFullWeight = mergeFull;
        case 'wmiCorr' 
            fprintf('Weighted nearest MI\n');            
            mergeFullWeight = cell2mat(cellfun(@(x)weightedNNcontrast(x),baseVect,'uniformoutput',0));
            
            cDist = squareform(pdist(mergeFullWeight','corr'));
            [cDistSorted] = sort(cDist,'ascend');            
        otherwise
            error('Unkown distance function');
                        
    end
    
    knnRankV = sort(cDistSorted(opts.knnFilterRank,:),'descend');        
    % kneeIdx = knee_pt(knnRankV);
    % kneeThr = knnRankV(kneeIdx);
    
    if isempty(opts.maxRankQ)
        [kneeThr,kneeIdx] = findElbow(knnRankV);
    else
        zt = ceil(length(knnRankV)*opts.maxRankQ);        
        [kneeThr,kneeIdx] = findElbow(knnRankV(1:zt));
    end

    if opts.doPlot
        zfig{1} = figure();
        plot(knnRankV,'linewidth',3)
        addLine(kneeIdx);           
        
        if ~isempty(opts.knnExcludeHardRank)
            addLine(opts.knnExcludeHardRank);
        end
        
    end

    if isempty(opts.knnExcludeHardRank)
        dropW = cDistSorted(opts.knnFilterRank,:) > kneeThr;
        fprintf('Dropping %d outlier gene programs\n',sum(dropW));
    else
        dropW = cDistSorted(opts.knnFilterRank,:) > cDistSorted(opts.knnFilterRank,opts.knnExcludeHardRank);
        fprintf('Dropping %d outlier gene programs\n',sum(dropW));
    end

    cMergeClean = mergeFull(:,~dropW);         

    if opts.doAlphaKernelKnn == 1       
        fprintf('Using phate kernel function');
        tic()
        simKnn = compute_alpha_kernel_sparse(cMergeClean, 'k', opts.knnSimK, 'a', opts.knnSimAlpha, 'distfun', opts,knnDistFun);
        toc()
    else
        fprintf('Using simple knn sim kernel');
        sKnnOpts.knnNum = opts.knnSimK;
        sknnOpts.distType =  opts.knnDistFun;
        
        simKnn = knndist2simMatrix(cMergeClean,sKnnOpts);
        simKnn = (simKnn + simKnn')./2;        
    end
    
    %%
    [subStability.S, subStability.N, subStability.VI, subStability.C] = stability(simKnn,opts.tRange,'v','p','Laplacian',opts.Laplacian);
    subStability = structSubSelectMat(subStability,subStability.N > opts.minClust & subStability.N < opts.maxClust);
    testK = subStability.N;
    skipIdx = false(size(subStability.N));
    if ~isempty(opts.maxGraphClustK) && max(subStability.N) > opts.maxGraphClustK
        fprintf('Ignoring clustering solutions larger than %d\n',opts.maxGraphClustK);
        skipIdx = subStability.N > opts.maxGraphClustK;
        testK(skipIdx) = [];
    end    
    
    cInitNorm = norm(cDataM,'fro');
    zi = 1;

    cMeanSil = nan(length(testK),1);
    cReconstError = nan(length(testK),1);     
    % cBestH = cell(length(testK),1);

    cReconstErrorPostOpt = nan(length(testK),1);

    for zk = 1:length(subStability.N)
        zt = tic;
        if skipIdx(zk)
            continue;
        end
        
        cSubCl = subStability.C(:,zk)+1;

        [zD] = silhouette(cMergeClean',cSubCl,'euclidean');
        cMeanSil(zi) = mean(zD);       
        
        cMergedW = cell2mat(arrayfun(@(x)median(cMergeClean(:,cSubCl == x),2),unique(cSubCl),'uniformoutput',0)');   
        % disp(cMergedW);        
               
        fprintf('Completing NMF by NNLS K=%d\n',subStability.N(zk));
        Hinf = nnlsm_nmf(cMergedW, cDataM, [],opts.nmfSolver);        
        cReconstError(zi) = norm((cMergedW*Hinf)-cDataM,'fro')/cInitNorm;
        
        if opts.fullOptimizeNMF
            fprintf('Reoptimizing by ccNMF\n');
%            if exist('batchID','var') && ~isempty(batchID)         
%                 [cMergedW,Hinf] = reoptimizeNMF(cDataM,cMergedW,Hinf,batchID,opts);
%            else
            [cMergedW,Hinf] = bp_nmf_opt(cDataM,size(cMergedW,2),cMergedW,Hinf,opts.bp_nmf);
            cReconstErrorPostOpt(zi) = norm((cMergedW*Hinf)-cDataM,'fro')/cInitNorm;
%           end
        end  
        
        mergedW{zi} = cMergedW;
        extrapH{zi} = Hinf;
       
        
        fprintf('(%d) Time %f., K=%d, Sil: %f, Err: %f, PostOpt: %f\n',zi,toc(zt),size(Hinf,1),cMeanSil(zi),cReconstError(zi),cReconstErrorPostOpt(zi));
        zi = zi+1;
    end
%%
    if opts.chooseBestHueristic == 1
        if opts.fullOptimizeNMF
            [~,bestIdx] = knee_pt(cReconstErrorPostOpt,testK);
        else
            [~,bestIdx] = knee_pt(cReconstError,testK);
        end
        % bestIdx = find(bestIdx == testK,1,'first');
        fprintf('Reconstruction knee point at %d\n',testK(bestIdx))
    elseif opts.chooseBestHueristic == 2
        [~,bestIdx] = knee_pt(cMeanSil,testK);
        % bestIdx = find(bestIdx == testK,1,'first');
        fprintf('Reconstruction knee point at %d\n',testK(bestIdx))
    end
    %%
    %%
    
    cMergedW = mergedW{bestIdx};
    Hinf = extrapH{bestIdx};
    
    %%
    if opts.doPlot
        figure('Position',[0 0 1800 900]);
        subplot(1,2,1);
        yyaxis left 
        plot(testK(:),cMeanSil);
        xlabel('GP cluster K');
        ylabel('Mean Sil');

        yyaxis right
        plot(testK(:),cReconstError);
        ylabel('Error');
        
        subplot(1,2,2);        

        cSubCl = subStability.C(:,bestIdx)+1;
        [~,zfig{2}] = silhouette(cMergeClean',cSubCl,'euclidean');       
       
    end
    

    cHinfNorm = Hinf./sum(Hinf);
    
    outConsensusNMF.bestConsH = cHinfNorm;
    outConsensusNMF.consW = cMergedW;  
    
    
    outConsensusNMF.testK = testK;
    
    outConsensusNMF.reconstError = cReconstError;
    if opts.fullOptimizeNMF
        outConsensusNMF.reconstErrorPostOpt = cReconstErrorPostOpt;
    end
    
    outConsensusNMF.cMeanSil = cMeanSil;
    
    
    outConsensusNMF.subStability = subStability;
    outConsensusNMF.baseMergeClean = cMergeClean;
    
    outConsensusNMF.mergedW = mergedW;
    outConsensusNMF.extrapH = extrapH;
    

end