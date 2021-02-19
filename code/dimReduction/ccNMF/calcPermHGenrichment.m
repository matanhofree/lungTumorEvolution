function outE = calcPermHGenrichment(geneSet,scoreIdx,geneID,expMat,inOpts)

    defaultOpts.permN = 1000;            
    defaultOpts.doPerm = 1;
    
    defaultOpts.bkgNumBins = 20;
    defaultOpts.minTestSize = 3;
    defaultOpts.minSetSize = 5;
    defaultOpts.maxSetSize = 2000;
    defaultOpts.pb = 1;
    defaultOpts.permTestByFDR = 0.1;
    
    defaultOpts.effectiveN = -1;
    defaultOpts.testGenesetOnlyGenes = 0;
    
                
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;    
    zfig = [];
    
    if opts.pb > 0 && ~ismac()
        opts.pb = 2;
    end
    
    nG = length(geneID);
    if size(scoreIdx,1) ~= nG
        error('Number of genes does not match scoreIdx');
    end
    
    if ~islogical(scoreIdx)
        scoreIdx = scoreIdx > 0;        
    end
    nT = sum(scoreIdx);
    
    if opts.effectiveN == -1
        opts.effectiveN = nG;
        fprintf('Effective N set to genome size (N=%d)\n',opts.effectiveN)
    end
    
    
    if ~exist('expMat','var')
        expMat = [];
    end
    
    if ~isempty(expMat) && size(expMat,1) ~= nG
        error('Number of genes does not match expMat');
    end
    
    % Prepare geneSet
    zix = ismember(geneSet.idxGeneKey,geneID);  
    geneSetSub = structSubSelectMat(geneSet,zix);   
    if ~isempty(opts.minSetSize)
        gsTotal = sum(geneSetSub.idxMat,2);        
        geneSetSub = structSubSelectMat(geneSetSub,gsTotal>opts.minSetSize);
    end
    if ~isempty(opts.maxSetSize)
        gsTotal = sum(geneSetSub.idxMat,2);        
        geneSetSub = structSubSelectMat(geneSetSub,gsTotal<opts.maxSetSize);    
    end
    clear geneSet;
    
    gsAnnot = sum(geneSetSub.idxMat,1)==0;
    if any(gsAnnot)
        geneSetSub.idxMat(:,gsAnnot) = [];
        geneSetSub.idxGeneKey(gsAnnot) = [];
    end
    
    geneSetSub.gsTotal = sum(geneSetSub.idxMat,2);
    
    effGeneSize = min(nG,length(geneSetSub.idxGeneKey));    
    if opts.testGenesetOnlyGenes == 1
        geneSetSub.effGeneSize = effGeneSize;
        fprintf('Note: testing only for HG enrichment of genes within geneSets - effective N %d\n',geneSetSub.effGeneSize);
    else        
        if isempty(opts.effectiveN) || effGeneSize > opts.effectiveN            
            geneSetSub.effGeneSize = effGeneSize;
            fprintf('Note: testing only for HG enrichment of genes in the genome or geneSet reference - effective N %d\n',geneSetSub.effGeneSize);
        else
            geneSetSub.effGeneSize = opts.effectiveN;
            fprintf('Note: testing only for HG enrichment within a heuristic genome size - effective N %d\n',geneSetSub.effGeneSize);
        end        
    end
    
    disp(geneSetSub)

    nCol = size(scoreIdx,2);
    nSets = length(geneSetSub.setGeneList);
      
    pval = nan(nSets,nCol);
    fdr = nan(nSets,nCol);
    outOverlap = nan(nSets,nCol);
    
    if ~isempty(expMat) && opts.permN > 0
        permPval = nan(nSets,nCol);
        permMean = nan(nSets,nCol);
        permStd = nan(nSets,nCol);
        opts.doPerm = 1;
    else
        opts.doPerm = 0;
    end
    
    opts.preBin = 0;
    if opts.doPerm && size(expMat,2) == 1
        opts.preBin = 1;
        [smpBins,smpBinsIdx,binCnt] = assignGeneToExpBin(expMat,opts);        
    end
    
    
    if opts.pb == 1
        progressbar();
    elseif opts.pb == 2
        [~,parForMonFile] = parfor_progress(nCol);
    end
    
    for zi = 1:nCol        
        if nT(zi) < opts.minTestSize
            fprintf('Skipping column %d -- %d too few genes for HG enrichment\n',zi,nT(zi));
            continue;
        end
        
        fprintf('Test column %d -- %d genes tested for HG enrichment\n',zi,nT(zi));
        
        
        if opts.doPerm       
            if ~opts.preBin 
                cMeanExp = expMat(:,zi);
                cSel = isnan(cMeanExp);
                cMeanExp(cSel) = [];
                
                geneC = geneID(~cSel);
                
                [smpBins,smpBinsIdx,binCnt] = assignGeneToExpBin(cMeanExp,opts);     
                
                cScoreIdx = scoreIdx(~cSel,zi);
            else 
                geneC = geneID;
                cScoreIdx = scoreIdx(:,zi);
            end
            [pval(:,zi),fdr(:,zi),outOverlap(:,zi),permPval(:,zi),permMean(:,zi),permStd(:,zi)] = testHGerichment(geneSetSub,geneC,cScoreIdx,smpBins,smpBinsIdx,opts);
        else
            [pval(:,zi),fdr(:,zi),outOverlap(:,zi)] = testHGerichment(geneSetSub,geneID,scoreIdx(:,zi),[],opts);
        end
        
        if opts.pb == 1
            progressbar(zi/nCol);
        elseif opts.pb == 2
            parfor_progress(-1,parForMonFile);
        end
    end

    outE.geneSetName = geneSetSub.setNames;
    outE.pval = pval;    
    outE.FDR = fdr;
    outE.overlap = outOverlap;
    outE.setSize = geneSetSub.gsTotal;
    
    if opts.doPerm 
        outE.permP = permPval;
        outE.permMean = permMean;
        outE.permStd = permStd;
    end
   

end

function [smpBins,smpBinsIdx,cntB] = assignGeneToExpBin(expMat,opts)

    binStep = 1/opts.bkgNumBins;
    binEdges = unique([ -inf quantile(expMat,binStep:binStep:(1-binStep)) inf]);    
    binEdges = setdiff(binEdges,0);

    [cntB,~,smpBins] = histcounts(expMat,binEdges); 
    
    smpBinsIdx = arrayfun(@(xi)find(smpBins == xi),unique(smpBins),'uniformoutput',0);   
end

function bkgSmpIdx = sampleBkgGenes(smpBinsIdx,nGenes,runBin,binC,nSmp)
    
    bkgSmpIdx = nan(nGenes,nSmp);

    for i = 1:length(binC)
       if (binC(i) == 0)
            continue;
       end
       
       zidx = runBin == i;
       if nSmp > 1
           bkgSmpIdx(zidx,:) = reshape(randsample(smpBinsIdx{i},binC(i)*nSmp,1),binC(i),nSmp);      
       else
           bkgSmpIdx(zidx,:) = randsample(smpBinsIdx{i},binC(i),1);      
       end
    end

end


function [outPval,outFdr,outOverlap,outPermP,outPermMean,outPermStd] = testHGerichment(geneSet,geneID,scoreV,smpBins,smpBinsIdx,opts)
        
    outP = [];
    outPermP = [];
    outPermMean = [];
    outPermStd = [];
    
    gsTotal = geneSet.gsTotal;    
    Ntotal = geneSet.effGeneSize;
    
    topGene = geneID(scoreV);        
    zidx = ismember(geneSet.idxGeneKey,topGene);
    
    outOverlap = sum(geneSet.idxMat(:,zidx)>0,2);                    
    if opts.testGenesetOnlyGenes == 1
        nG = sum(zidx);
    else 
        nG = sum(scoreV);
    end   

    outPval = hygecdf(outOverlap-1,Ntotal,gsTotal,nG,'upper');        
    outFdr = mafdr(outPval,'BHfdr',1);


    if opts.doPerm        
        nTest = length(geneSet.setNames);

        if isempty(opts.permTestByFDR)
            
            negLogP = -log(outPval);
            
           [outPermP,outPermMean,outPermStd] = permTest(geneSet,negLogP,geneID,scoreV,smpBins,smpBinsIdx,opts);
        else            
            outPermP = nan(nTest,1);
            outPermMean = nan(nTest,1);
            outPermStd = nan(nTest,1); 

            if opts.permTestByFDR < 1
                zSubSet = outFdr < opts.permTestByFDR;
            else                
                [cP,cPvalIdx] = sort(outPval);
                
                cPvalIdx(isnan(cP) | cP == 1) = [];                
                
                if ~isempy(cPvalIdx)
                    zSubSet = trueV(cPvalIdx,nTest);
                end                                
            end

            if sum(zSubSet) > 0 
                negLogP = -log(outPval(zSubSet));
                geneSetSub = structSubSelectMat(geneSet,zSubSet);

                [outPermP(zSubSet),outPermMean(zSubSet),outPermStd(zSubSet)] = permTest(geneSetSub,negLogP,geneID,scoreV,smpBins,smpBinsIdx,opts);
            end

        end
    end
end

function [outPermP,outPermMean,outPermStd] = permTest(geneSetSub,negLogP,geneID,scoreV,smpBins,smpBinsIdx,opts)
           
    nTest = length(geneSetSub.setNames);
    outBkgScore = zeros(nTest,1);
    outBkgScoreM2 = zeros(nTest,1);
    outBkgRank = zeros(nTest,1);       
    
    Ntotal = geneSetSub.effGeneSize;
    
    topGene = geneID(scoreV);        
    zidx = ismember(geneSetSub.idxGeneKey,topGene);
      
    testGSonly = opts.testGenesetOnlyGenes;

    if testGSonly ~= 1
        nGlobal = sum(zidx);
    else 
        nGlobal = sum(scoreV);
    end   
    
    gsTotal = geneSetSub.gsTotal;
    
    fprintf('Testing by permutation -- %d\n',nTest);
    runBin = smpBins(scoreV);
    binC = histc(runBin,1:max(runBin));            

    nRep = opts.permN;
    for i = 1:nRep

        geneSmpIdx = sampleBkgGenes(smpBinsIdx,nGlobal,runBin,binC,1);
        geneSmp = geneID(geneSmpIdx);
        
        zidx = ismember(geneSetSub.idxGeneKey,geneSmp);           
        permOverlap = nansum(geneSetSub.idxMat(:,zidx)>0,2);    
        
        permNegLog = -log(hygecdf(permOverlap-1,Ntotal,gsTotal,nGlobal,'upper'));  
        permNegLog(isnan(permNegLog)) = 0;


        outBkgScore = outBkgScore + permNegLog;
        outBkgScoreM2 = outBkgScoreM2 + permNegLog.^2;

        outBkgRank = outBkgRank + (permNegLog > negLogP);
    end

    outBkgScoreStd = (outBkgScoreM2 - (outBkgScore.^2)./nRep)./(nRep-1);

    % outP.permP = (outBkgRank+1)./(nrep+1);
    outPermMean = outBkgScore/nRep;
    outPermStd= sqrt(outBkgScoreStd);

    outPermP = (outBkgRank+1)./(nRep+1);
    fprintf('Done testing by permutation\n');

end