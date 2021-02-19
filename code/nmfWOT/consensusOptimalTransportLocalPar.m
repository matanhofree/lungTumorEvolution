function [otMat,otVar,outStat] = consensusOptimalTransportLocalPar(Xorig,weightMat,stageV,testSet,inOpts)

    defaultOpts.doNorm = 0;
    defaultOpts.epsV = 1e-6;
    
    defaultOpts.pcaMethod = 'NeNMF'; % Alternative: fb_pca, wpca, standard
    defaultOpts.seed = 42;
    defaultOpts.resampleNMF = 5;
 
    defaultOpts.niter = 25;
    defaultOpts.normFun = @(X)nnzcenter(X,nnzmean(X,1));

    defaultOpts.nstarts = 1;
    defaultOpts.fb_iter = 5;
    defaultOpts.verbose = 1;
    
    defaultOpts.tsneMethod = 'resample'; % resampling - Resampling based approach to pca distance
                                           % fast       - bh tsne
                                           % cosine     - cosine similarity
                                           % dist       - Use externa distance function 
    
    defaultOpts.resampleN = 20;
    defaultOpts.resampleFrac = [ 0.8 0.8 ];
%     defaultOpts.numNeighbors = [ 10 30 ]; 
    defaultOpts.dimPCA = [ 7 30 ];

    defaultOpts.distType = 'cosine';
    defaultOpts.distMPower = 1;
    
    defaultOpts.otLambda1 = 1;
    defaultOpts.otLambda2 = 25;
    defaultOpts.otEpsilon = 0.1;
    defaultOpts.otNIter = 1000;
    defaultOpts.otGrowthRate = 1;
    

    defaultOpts.graphDist = 0;

    defaultOpts.pb = 1;
    defaultOpts.verbose = 0;
   
    defaultOpts.normCost = 1;
    defaultOpts.globalDist = 0;
    defaultOpts.pot_eps = 1e-6;
    
    defaultOpts.ncores = 1;
    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;

    outStat = [];
    disp(opts);
    %%
    [D,N] = size(Xorig);
    
    % Normalize     
    X = normalize_expMatrix(Xorig,[],opts);
    
    % Prepare data structs Compute pairwise distance matrix
    for zj = 1:size(testSet,1)
        subsetA{zj} = (ismember(stageV,testSet{zj,1}));
        subsetB{zj} = (ismember(stageV,testSet{zj,2}));
        
        xA{zj} = X(:,subsetA{zj});
        xB{zj} = X(:,subsetB{zj});
    end

    zT = tic;
    if opts.pb
        progressbar();
    end

    zc = 1;
    zTotal = length(testSet)*opts.resampleN;
    
    %% Start pool
    if opts.ncores > 1
        
        if isempty(gcp('nocreate'))
            maxCore = feature('numCores');
            parpool(min(maxCore,opts.ncores));
        end
        
        parfor zj = 1:length(testSet)
            fprintf('Testing: %s to %s\n',testSet{zj,1},testSet{zj,2});     
            [countMat{zj},countMatSq{zj},countDenom{zj}] = extractOT(xA{zj},xB{zj},opts);
        end          
    else       
        for zj = 1:length(testSet)
            fprintf('Testing: %s to %s\n',testSet{zj,1},testSet{zj,2});     
            [countMat{zj},countMatSq{zj},countDenom{zj}] = extractOT(xA{zj},xB{zj},opts);
        end  
    end
    outStat.countMat = countMat;
    outStat.countMatSq = countMatSq;
    outStat.countDenom = countDenom;
    outStat.opts = opts;        

    if opts.globalDist == 1        
        error('Not Implemented');
    end
    %%
    otMat = cellfun(@(x,d)x./d,countMat,countDenom,'uniformoutput',0);
    otVar = cellfun(@(x,x2,d)(x2 - (x.^2)./d)./(max(d-1,1)),countMat,countMatSq,countDenom,'uniformoutput',0);

end


function [countMat,countMatSq,countDenominator] = extractOT(xA,xB,opts)

    nA = size(xA,2);
    nB = size(xB,2);

    vA = ones(nA,1);
    vB = ones(nB,1);
    
    countMat = nan(nA,nB);
    countMatSq = nan(nA,nB);
    countDenominator = zeros(nA,nB);

        
    for i = 1:opts.resampleN

        costMat = estimateCostMatrix(xA,xB,opts);
        
        if ~isempty(opts.distMPower) & opts.distMPower ~= 0
            costMat = costMat.^opts.distMPower;
        end

        if (opts.normCost)
            costMat = costMat/median(costMat(:));
        end 
       
        outTransport = computeTransportStable(vA,vB,costMat,opts.otLambda1,opts.otLambda2,opts.otEpsilon,opts.otNIter,opts.otGrowthRate,opts.verbose);

        countMat = nansum(cat(3,countMat,outTransport),3);
        countMatSq = nansum(cat(3,countMatSq,outTransport.^2),3);
        countDenominator = countDenominator + 1;

%         if opts.pb
%             progressbar(zc/zTotal);
%         end
%         zc = zc + 1;
    end
end

function costMat = estimateCostMatrix(dA,dB,opts)

    opts.doNorm = 0;
    opts.resampleN = opts.resampleNMF;
    
    nA = size(dA,2);
    nB = size(dB,2);
      
    xData = [dA dB];
    selA = false(nA+nB,1);
    selA(1:nA) = 1;    
    
    [~,XscoreOut] = consensusGraphDistNMFConcat(xData,opts);


    switch opts.distType
        case 'cosine'
            xScoreMerged = cell2mat(cellfun(@(x)x./sum(x,2),XscoreOut,'uniformoutput',0));        
            costMat = pdist2(xScoreMerged(selA,:),xScoreMerged(~selA,:),opts.distType);            
        case 'logEucledian'
            xScoreMerged = cell2mat(cellfun(@(x)(x+opts.pot_eps)./sum(x+opts.pot_eps,2),XscoreOut,'uniformoutput',0));    
            xScoreMerged = -log(xScoreMerged);
            costMat = pdist2(xScoreMerged(selA,:),xScoreMerged(~selA,:),'euclidean');            
        case 'euclidean'
            costMat = pdist2(xScoreMerged(selA,:),xScoreMerged(~selA,:),opts.distType);            
        otherwise 
            error('Metric not implemented');
    end
end
