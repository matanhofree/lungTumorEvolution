function [outScore,outSignalScore,outBkgScore,outScaledDiff,rnkPos,clustE] = genSignatureScoreSimple(countData,signatureList,groupID,clusterSet,inOpts)

    defaultOpts.bkgScoreSmpNum = 1000;
    defaultOpts.bkgNumBins = 20;
    
    defaultOpts.aggrFunc = @(x)nanmean(x,1);   
    defaultOpts.varFunc = @(x)nanstd(x,0,1);
    % defaultOpts.varFunc = @(x)mad(x,1,1);
    
    defaultOpts.minGenes = 3;
    defaultOpts.maxGenes = [];
    
    defaultOpts.inputCounts = 'normTPM';
    defaultOpts.normType = 0; % Ansecombie normalization
                              % log normalization
                              
    defaultOpts.patientMaxExp = [];
    defaultOpts.calcRank = 0;
    defaultOpts.zscore = 1;
    defaultOpts.softMax = 1;
    
    defaultOpts.stdMin = 0.1;
        
    defaultOpts.scoreType = 1; % Simple mean
    
    defaultOpts.patientMaxExp = [];
    defaultOpts.filterByFreq = 0.01;
    
    defaultOpts.dropNegativeCorr = 0.5;
    defaultOpts.useIterative = 1;
    defaultOpts.geneID = 'geneID';
    
    defaultOpts.minGeneUseMean = 1;

    
    defaultOpts.clusterSetScore = 1;               
    defaultOpts.ncores = 1;
    defaultOpts.binsPreSample = 0;
    defaultOpts.pb = 0;
    defaultOpts.stdMeanThr = [];
    
    defaultOpts.minCorrPos = 0;
    defaultOpts.stdOnNNZ = 0;
           
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    disp(opts);            
    
    if opts.pb > 0 & ~ismac()
        opts.pb = 2;
    end
    
    if nargin < 4 || isempty(clusterSet)
        opts.clusterSetScore = 0;
    end
    
    if iscell(clusterSet) || isnumeric(clusterSet)
        clSet.cl = clusterSet;
        clusterSet = clSet;
        clear clSet;
    end
    
    if isempty(groupID)
        groupID = ones(length(countData.sampleID),1);
    else        
        error('Not implemented -- Implementation is likely wrong');
        groupID = grp2idx(groupID(:));
    end        
    
    expGenes.expData = countData.(opts.inputCounts);
    expGenes.geneID = countData.(opts.geneID);         
    
    % Normalize data
    if opts.normType == 1
        fprintf(1,'Using anscombie\n');
        expGenes.expData = 2*sqrt(expGenes.expData+3/8);
    elseif opts.normType == 2
        fprintf(1,'Using log+1\n');
        expGenes.expData = log1p(expGenes.expData);
    elseif opts.normType == 0
        fprintf(1,'Using unnormalized data!\n');
        % expGenes.expData = full(expGenes.expData+3/8);
    else 
        error('Normtype not implemented');
    end
    
    % Clean outlier patients 
    if ~isempty(opts.patientMaxExp)
        fprintf(1,'Applying max patient Q%f threshold cap\n',opts.patientMaxExp);
        patientQThr = quantile(expGenes.expData,opts.patientMaxExp);

        for zi = 1:size(expGenes.expData,2)
            expGenes.expData(expGenes.expData(:,zi) > patientQThr(zi),zi) = patientQThr(zi);
        end
    end
        
    % Prepare data 
    cMat = expGenes.expData';
    [N,D] = size(cMat);
    expGenes.meanExp = full(mean(cMat))';
  
    if opts.stdOnNNZ == -1
        fprintf('Skipping std calculation');
        expGenes.stdExp = expGenes.meanExp;
        expGenes.stdExp(:) = 1;
    elseif opts.stdOnNNZ == 1
        fprintf('Std calculation on nonzeros');
        expGenes.stdExp = arrayfun(@(x)std(full(nonzeros(cMat(:,x)))),1:D)';        
    else
        expGenes.stdExp = arrayfun(@(x)std(full(cMat(:,x))),1:D)';
    end
    
    if ~isempty(opts.stdMeanThr)         
         expGenes.stdExp = max(expGenes.stdExp,expGenes.meanExp*opts.stdMeanThr);
    end
    
    expGenes.freqExp = full(sum(cMat>0))'/N;
    clear cMat;       
    
    fprintf('Done with mean/std');
    if ~isempty(opts.filterByFreq)
        zSel = expGenes.freqExp > opts.filterByFreq;
        fprintf(1,'Keeping only genes with exp freq > %f (%d/%d)\n',opts.filterByFreq,sum(zSel),length(zSel));
        expGenes = structSubSelectMat(expGenes,zSel);
    end
    
    if opts.zscore == 1 
        fprintf(1,'Z normalizing the inputs (data will be made dense)\n');        
        expGenes.expData = full(expGenes.expData);
        expGenes.expData = bsxfun(@minus,expGenes.expData,expGenes.meanExp);
        expGenes.expData = bsxfun(@rdivide,expGenes.expData,max(expGenes.stdExp,opts.stdMin));    
    elseif opts.zscore == 2 
        fprintf(1,'Scaling expression by std\n');        
        expGenes.expData = full(expGenes.expData);
        expGenes.expData = bsxfun(@rdivide,expGenes.expData,max(expGenes.stdExp,opts.stdMin))
    end
    if opts.softMax == 1
        fprintf(1,'Sigmoid transform normalizing');
        expGenes.expData = (1+exp(-expGenes.expData)).^-1;
    end
      
    %% Generate background bin set    
    binStep = 1/opts.bkgNumBins;
    binEdges = unique([ -inf quantile(nonzeros(expGenes.meanExp),binStep:binStep:(1-binStep)) inf]);    
    binEdges = setdiff(binEdges,0);
    
    [zCnt,~,expGenes.smpBins] = histcounts(expGenes.meanExp,binEdges);    
    opts.binNum = 1:length(zCnt);
    
    expGenes.smpBinsIdx = arrayfun(@(xi)find(expGenes.smpBins == xi),opts.binNum,'uniformoutput',0);   
    %%
    fieldNames = [];
    if isstruct(signatureList)
        fieldNames = fieldnames(signatureList);
        signatureList = struct2cell(signatureList);
    end
    
    [idxSigKey,sigListMat] = geneSet2IdxMat(signatureList);   

    overLapGeneID = ismember(idxSigKey,expGenes.geneID);
    cntOverlap = sum(overLapGeneID);
    fprintf('Overlap berween gene set ID and gene IDs: %d (%f)\n',cntOverlap,cntOverlap/length(idxSigKey));
    if cntOverlap < 1
        error('Gene name space of signature set and data are not matched.');
    end
%%
    sigListMat = sigListMat(:,overLapGeneID);
    sigListSize = full(sum(sigListMat,2));
%%
    cDrop = false(size(sigListMat,1),1);
    if ~isempty(opts.minGenes)
        cDrop =  sigListSize <= opts.minGenes & cDrop;
        fprintf('Dropping %d sets due to being too small\n',sum(cDrop));
    end 

    if ~isempty(opts.maxGenes)        
        cDrop = sigListSize >= opts.maxGenes & cDrop;
        fprintf('Dropping %d sets due to being too large or too small\n',sum(cDrop));
    end 
    
    if sum(cDrop) > 0
        fieldNames(cDrop) = [];
        signatureList(cDrop) = [];
    end
  %%
        
    groupList = unique(groupID);

    nG = length(groupList);
    % nN = length(groupID);
    nS = length(signatureList);
    if opts.clusterSetScore > 0
        nCl = length(fieldnames(clusterSet));        
    end

    if opts.ncores == 1 
        outScore = nan(N,nS);  
        outSignalScore = nan(N,nS);
        outBkgScore = nan(N,nS);
        outScaledDiff = nan(N,nS);
        clustE = nan(1,nS);
        if opts.calcRank
            rnkPos = nan(N,nS);
        else 
            rnkPos = [];
        end
        if opts.clusterSetScore > 0        
            clustE = nan(nCl,nS);
        end
                
        if nG > 1
            for gi = 1:nG
                selID = groupID == groupList(gi);
                if all(selID)
                    subExpGenes = expGenes;
                else
                    subExpGenes = structSubSelectMat(expGenes,selID);
                end
                for zi = 1:nS
                    sigGenes = ismember(expGenes.geneID,signatureList{zi}); 
                    if sum(sigGenes) < opts.minGenes
                        fprintf('Skipping %d - found insufficient genes (%d)\n',zi,opts.minGenes);                                        
                        continue;
                    end                    
                    [ outScore(selID,zi),outSignalScore(selID,zi),outBkgScore(selID,zi),outScaledDiff(selID,zi),rnkPos(selID,zi) ] = calcSignatureList(subExpGenes,sigGenes,clusterSet,opts);                     
                end
            end 
        else
            for zi = 1:nS
                %%
                fprintf('Testing %d\n',zi);
                sigGenes = ismember(expGenes.geneID,signatureList{zi}); 
                if sum(sigGenes) < opts.minGenes
                    if opts.minGeneUseMean == 1
                        fprintf('Found (%d) insufficient genes (%d) -- reporting mean\n',zi,sum(sigGenes));
                        outScore(:,zi) = mean(expGenes.expData(sigGenes,:),1);
                        outSignalScore(:,zi) = outScore(:,zi);
                        outBkgScore(:,zi) = 0;

                    else
                        fprintf('Skipping %d - found insufficient genes (%d)\n',zi,sum(sigGenes));
                    end
                    continue
                end 
                if opts.calcRank
                    [ outScore(:,zi),outSignalScore(:,zi),outBkgScore(:,zi),outScaledDiff(:,zi),rnkPos(:,zi),clustE(:,zi)] = calcSignatureList(expGenes,sigGenes,clusterSet,opts);    
                else
                    [ outScore(:,zi),outSignalScore(:,zi),outBkgScore(:,zi),outScaledDiff(:,zi),~,clustE(:,zi)] = calcSignatureList(expGenes,sigGenes,clusterSet,opts);    
                end
            end
        end
        
        outScore = num2cell(outScore,1);  
        outSignalScore = num2cell(outSignalScore,1);
        outBkgScore = num2cell(outBkgScore,1);
        outScaledDiff = num2cell(outScaledDiff,1);
        if opts.calcRank
            rnkPos = num2cell(rnkPos,1);
        end
        
        clustE = num2cell(clustE,1); 
        
        if ~isempty(fieldNames)        
            outScore = cell2struct(outScore,fieldNames,2);
            outSignalScore = cell2struct(outSignalScore,fieldNames,2);
            outBkgScore = cell2struct(outBkgScore,fieldNames,2);
            outScaledDiff = cell2struct(outScaledDiff,fieldNames,2); 
            if opts.calcRank
                rnkPos = cell2struct(rnkPos,fieldNames,2);
            end
            clustE = cell2struct(clustE,fieldNames,2);
        end
        
    else       
        parfor zi = 1:nS 
            fprintf('Testing %d\n',zi);
            sigGenes = ismember(expGenes.geneID,signatureList{zi}); 
            [outScore{zi},outSignalScore{zi},outBkgScore{zi},outScaledDiff{zi},rnkPos{zi},clustE{zi}] = calcSignatureList(expGenes,sigGenes,clusterSet,opts);    
        end 
        
        %%
        outScore = cellfun(@(x)x',outScore,'uniformoutput',0);  
        outSignalScore = cellfun(@(x)x',outSignalScore,'uniformoutput',0); 
        outBkgScore = cellfun(@(x)x',outBkgScore,'uniformoutput',0); 
        outScaledDiff = cellfun(@(x)x',outScaledDiff,'uniformoutput',0); 
        rnkPos = cellfun(@(x)x',rnkPos,'uniformoutput',0); 
        clustE = cellfun(@(x)x',clustE,'uniformoutput',0); 
        %%
        if ~isempty(fieldNames)        
            outScore = cell2struct(outScore,fieldNames,2);
            outSignalScore = cell2struct(outSignalScore,fieldNames,2);
            outBkgScore = cell2struct(outBkgScore,fieldNames,2);
            outScaledDiff = cell2struct(outScaledDiff,fieldNames,2);    
            rnkPos = cell2struct(rnkPos,fieldNames,2);
            clustE = cell2struct(clustE,fieldNames,2);
        end
    end
        
    
    if nargout == 1
        sampleID = [];
        if isfield(countData,'sampleID')
            sampleID = countData.sampleID;
        end
        outScoreOut = structpack([],sampleID,outScore,outSignalScore,outBkgScore,outScaledDiff,rnkPos,clustE)
        outScore = outScoreOut;
    end  
end

function [ outScore,outSignalScore,outBkgScore,outScoreDiff,rnkPos,clustE ] = calcSignatureList(expGenes,sigGenes,clusterSet,opts)    
    
    outScore = [];
    outSignalScore = [];
    outBkgScore = [];
    outScoreDiff = [];
    rnkPos = [];
    clustE = [];

    if sum(sigGenes) < opts.minGenes
        fprintf('Found insufficient number of genes from signature %d\n',sum(sigGenes));
        return 
    end
    %%  
    sigGenesList = find(sigGenes);
    nG = length(sigGenesList);    
    sigMat = expGenes.expData(sigGenes,:); 
    
    if opts.dropNegativeCorr > 0
        zD = corr(sigMat');
        corrFreq = sum(zD > opts.minCorrPos)./nG;
        
        remGenes = full(corrFreq < opts.dropNegativeCorr);
        if sum(remGenes) > 0
            fprintf('Dropping %d/%d genes (cutoff %f)\n',sum(remGenes),nG,opts.dropNegativeCorr);               
            if (nG - sum(remGenes)) > opts.minGenes
                sigGenes(sigGenesList(remGenes)) = 0;           
                sigGenesList(remGenes) = [];
                sigMat(remGenes,:) = [];
                nG = length(sigGenesList);
            else
                fprintf('Note: unable to drop genes by corr filter - too few will be left\n');
            end
        end
    end
   
      
    outSignalScore = opts.aggrFunc(sigMat);            
    runBin = expGenes.smpBins(sigGenes);
    N = size(expGenes.expData,2);
        
    binC = histc(runBin,opts.binNum);    
    
        
    if opts.clusterSetScore        
        [outE,outP] = calcClusterEnrichment(outSignalScore,clusterSet);
        outBkgClustE = nan(opts.bkgScoreSmpNum,length(outE));
    end
    
    %%
    if opts.binsPreSample
       bkgSmpIdx = sampleBkgGenes(expGenes,nG,runBin,binC,opts.bkgScoreSmpNum)
    end 
    

    %%
    if opts.useIterative || opts.binPreSample == 0
        %%
        outBkgScore = zeros(1,N);
        outBkgScoreM2 = zeros(1,N);
                
        % outBkgScoreStd = zeros(1,nN);
        if opts.calcRank
            nM = luniq(clusterSet.cl);
            rnkPos = nan(1,nM);
        else 
            rnkPos = [];
        end
                
        nRep = opts.bkgScoreSmpNum;    
        if opts.pb == 1
            progressbar();
        elseif opts.pb == 2
            [~,parForMonFile] = parfor_progress(nRep);
        end
        for zi = 1:nRep            
            
            if  opts.binsPreSample
                mV = opts.aggrFunc(expGenes.expData(bkgSmpIdx(:,zi),:));
            else
                bkgSmpIdx = sampleBkgGenes(expGenes,nG,runBin,binC,1);
                mV = opts.aggrFunc(expGenes.expData(bkgSmpIdx,:));
            end
            outBkgScore = outBkgScore + mV;
            outBkgScoreM2 = outBkgScoreM2 + mV.^2;                        
            if opts.calcRank
                error('Not implemented');
            end
            if opts.clusterSetScore 
                outBkgClustE(zi,:) = calcClusterEnrichment(mV,clusterSet);
            end

            if opts.pb == 1
                progressbar(zi/nRep);
            elseif opts.pb == 2
                parfor_progress(-1,parForMonFile);
            end
        end        
        
        outBkgScoreStd = (outBkgScoreM2 - (outBkgScore.^2)./nRep)./(nRep-1);
        outBkgScoreStd = sqrt(outBkgScoreStd);
        outBkgScore = outBkgScore./nRep;             
    else
        smpData = squeeze(opts.aggrFunc(reshape(full(expGenes.expData(bkgSmpIdx,:)),nG,opts.bkgScoreSmpNum,N)));
        outBkgScore = opts.aggrFunc(smpData);
        outBkgScoreStd = opts.varFunc(smpData);
        
        if opts.calcRank
            rnkPos = min(sum(bsxfun(@ge,outSignalScore,smpData)),sum(bsxfun(@le,outSignalScore,smpData)));
        else 
            rnkPos = nan;
        end
    end

    zStd = min(nonzeros(outBkgScoreStd));
    if isempty(zStd)
        outBkgScoreStd(outBkgScoreStd == 0) = 0.01;
    else
        outBkgScoreStd(outBkgScoreStd == 0) = min(nonzeros(outBkgScoreStd));
    end

    outScore = outSignalScore - outBkgScore;
    outScoreDiff = outScore./outBkgScoreStd;
    
    if opts.clusterSetScore
        clustE = sum(outE < outBkgClustE);
    else 
        clustE = nan;
    end
    
  
end

function bkgSmpIdx = sampleBkgGenes(expGenes,nGenes,runBin,binC,nSmp)

    
    bkgSmpIdx = nan(nGenes,nSmp);

    for i = 1:length(binC)
       if (binC(i) == 0)
            continue;
       end
       
       zidx = runBin == i;
       if nSmp > 1
           bkgSmpIdx(zidx,:) = reshape(randsample(expGenes.smpBinsIdx{i},binC(i)*nSmp,1),binC(i),nSmp);      
       else
           bkgSmpIdx(zidx,:) = randsample(expGenes.smpBinsIdx{i},binC(i),1);      
       end
    end

end


function [outE,outP] = calcClusterEnrichment(outSignalScore,clusterSet)

    fClNames = fieldnames(clusterSet);

    for i = 1:length(fClNames)

        clName = fClNames{i};
        cl = clusterSet.(clName);
        
        [outP(i),zTbl] = anova1(outSignalScore,cl,'off');
 
        outE(i) = zTbl{2,5};
    end

end