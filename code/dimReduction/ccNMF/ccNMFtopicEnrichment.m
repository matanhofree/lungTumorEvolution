function ccEnrichment = ccNMFtopicEnrichment(ccNMF,outPathStub,countData,geneSet,ydata,outCl,timeP,inOpts)

    defaultOpts.basisExpansion = 1;
    defaultOpts.plotOverview = 1;
    defaultOpts.plotTopicEnrichment = 1;
    defaultOpts.plotTopN = 50;

    defaultOpts.outGeneTable = 1;
    defaultOpts.outGeneSetTable = 1;
    
    defaultOpts.geneTableTopN = 50;
    defaultOpts.geneSetTableTopN = 2000;
    defaultOpts.useDynamic = 200;
    defaultOpts.dynamicMaxSetSize = 100;
    

    defaultOpts.geneID = 'geneID';
    defaultOpts.printType = '-dpng';

    defaultOpts.trimWeightByQfrac = 0.1;
    defaultOpts.trimWeightByQmean = 0.1;
    defaultOpts.filterRiboMt = 1;
    defaultOpts.filterGeneSetOnly = 1;
    
    % defaultOpts.normH = 0;

    defaultOpts.enrichmentThr = 200;
    defaultOpts.fdrThr = 0.1;
    defaultOpts.maxTestFisher = 100;
    
    defaultOpts.doPerm = 1;

    defaultOpts.sortBy = 'Pval'; % or Jaccard

    defaultOpts.plotSetSize = 1;

    defaultOpts.minSetSize = 10;
    defaultOpts.maxSetSize = [];
    defaultOpts.plotFigViz = 'off';
    
    defaultOpts.corrType = 'spe';    
    defaultOpts.includeTopTF = 5;
    defaultOpts.markSig = 1;
    defaultOpts.pvalThr = 0.05;

    defaultOpts.useMacc = 0;
    defaultOpts.doTF = 0;
    
    defaultOpts.useJaccardSetClustering = 1;
    defaultOpts.jaccardSetThreshold = 0;
    defaultOpts.jaccardSetThresholdFactor = 1.33;
    defaultOpts.jaccardSetBinSort = 1;
    defaultOpts.permFDR_thr = 0.1;
    
    defaultOpts.saveIntermediate = 1;
    defaultOpts.forceRecalc = 0;
    
    defaultOpts.applyZscore = 1;
    
    defaultOpts.normDataVar = 'normTPM';
    
    
    defaultOpts.subAxisOpts = { 'SpacingVertical',0.02,'SpacingHorizontal',0.08,'MarginLeft',.08,'MarginRight',.08};
    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts

    if nargin < 2
        error('Must input ccNMF and outPath')
    end

    if nargin < 3
        countData = [];
    end

    if nargin < 4
        geneSet = [];
    end

    if nargin < 5
        ydata = [];
    end

    if nargin < 6
        outCl = [];
    end

    if nargin < 7
        timeP = [];
    end


    % Expand basis as needed
    if opts.basisExpansion && ~isfield(ccNMF,'expandedW') && ~isempty(countData)
        error('Not yet implemented');
    end
    
    if length(countData.sampleID) ~= length(ccNMF.sampleID)
        bothID = intersect(countData.sampleID,ccNMF.sampleID);
        zSel = ismember(countData.sampleID,bothID);
        if ~all(zSel)
            countData = structSubSelectMat(countData,zSel);
        end
        
        zSel = ismember(ccNMF.sampleID,bothID);
        if ~all(zSel)
            ccNMF = structSubSelectMat(ccNMF,zSel);
        end
        
    end

    % Select output to plot
    if isfield(ccNMF,'bestConsH')
        outH = ccNMF.bestConsH;
        % zi = find(ccNMF.testK == size(outH,1));
        if isfield(ccNMF,'Kselect')
            zi = ccNMF.Kselect;        
            outW = ccNMF.expandedW{zi};

        elseif isfield(ccNMF,'extrapH')
            zi = find(cellfun(@(x)isequal(outH,x),ccNMF.extrapH));
            outW = ccNMF.expandedW{zi};
        else
            error('This looks like a problem');
            % outW = ccNMF.expandedW{1};
        end
    else
        outH = ccNMF.H;
        outW = ccNMF.W;        
    end
        
        
    geneID = ccNMF.(opts.geneID);
    %%
    % Todo make sure normTPM exists    
    dataM = countData.(opts.normDataVar);
    
    if opts.Hnorm == 1
        fprintf('H matrix is scaled to sum to one\n');
        outH = outH./sum(outH);
    elseif opts.Hnorm == 2
        fprintf('H matrix is scaled to unit norm (normalized by root mean square)\n');
        outH = outH./sqrt(sum(outH.^2));
    else
        fprintf('H matrix is not normalized');        
    end
          
    % Plot cc enrichment overview
    if opts.plotOverview
        opts.titleText = 'CC nmf %d ';
        zfig = plot_tsne_scatter_multi(ydata,outH,[],[],opts);
        if ~isempty(opts.printType)
            for i = 1:length(zfig)
                outP = sprintf('%s_ccNMF_overview_%02d',outPathStub,i);
                print(zfig{i},outP,opts.printType);
            end
        end
    end

    %% Write table of top genes
    outWeight = weightedNNcontrast(outW);
    cHsto = outH./sum(outH,2);
    clMean = full(dataM*cHsto');    

    weightFilter = outWeight;    
    

    if ~isempty(opts.trimWeightByQfrac)
        fprintf('Filtering by quantile exp fraction\n');
        clFrac = (dataM>0)*cHsto';
        zFilter = clFrac > quantile(clFrac,opts.trimWeightByQfrac);
        weightFilter(~zFilter) = nan;
        
        clear clFrac;
    end

    if ~isempty(opts.trimWeightByQmean)
        fprintf('Filtering by mean expression\n');        
        zFilter = clMean > quantile(clMean,opts.trimWeightByQmean);
        weightFilter(~zFilter) = nan;
    end
    clear cHsto    
    
    if opts.filterRiboMt == 1
        fprintf('Filtering out ribosomal and mitochondrial\n');
        zFilter = strgrepi(geneID,'RPS.*') | strgrepi(geneID,'RPL.*') | strgrepi(geneID,'MT-.*');
        weightFilter(zFilter,:) = nan;
    end

    if opts.filterGeneSetOnly == 1
        zFilter = ismember(geneID,geneSet.idxGeneKey);
        fprintf('Filtering to include only geneset genes (%d)\n',sum(zFilter));
        weightFilter(~zFilter,:) = nan;
    end
 
    [D,K] = size(weightFilter);
    [cVal,cidxWord] = sort(weightFilter,'descend','MissingPlacement','last');
    
    if ~isempty(opts.useDynamic) && opts.useDynamic > 0 
        fprintf('Using dynamic elbow based threshold rule.\n');
        cValThr = findElbow(cVal(1:opts.useDynamic,:)) ;
        cValThr = max(cValThr,cVal(opts.dynamicMaxSetSize,:));                
    else 
        fprintf('Using hard rank threshold rule.\n');
        
        cValThr = cVal(opts.geneTableTopN,:);
    end
    
    % Filter by nnzmax 
    if any(cValThr <= 0)
        fprintf('The following NMFs result in 0 threshold\n');
        cFzero = find(cValThr == 0);
        cFzero = cFzero(:)';
        disp(cFzero);
        %%
        cPos = @(x)x(x>0);
        cValThrNnz = arrayfun(@(x)min(cPos(cVal(:,x))),1:K,'unif',0);
        cValThrNnz(isemptycell(cValThrNnz)) = { nan };
        cValThrNnz = cell2mat(cValThrNnz);
       %%
        cValThr = max(cValThr,cValThrNnz);
        cValThr(cValThr == 0) = nan;
        
    end
    %%
    cidxGeneSelect = weightFilter >= cValThr;
    cidxGeneSelect(isnan(weightFilter)) = 0;    
    
    cidxGeneSelectTotal = sum(cidxGeneSelect);
    disp(cidxGeneSelectTotal)       

    outTopGenes = cell(max(cidxGeneSelectTotal),K);
    outTopGenes(:) = { '' };
    %%
    for i = 1:K
        % outTopGenes(1:cidxGeneSelectTotal(i),i) = geneID(cidxGeneSelect(:,i));
        cSelTop = 1:cidxGeneSelectTotal(i);
        outTopGenes(1:cidxGeneSelectTotal(i),i) = geneID(cidxWord(cSelTop,i));
    end
    
    if opts.outGeneTable
        outP = sprintf('%s_ccNMF_topGene.xlsx',outPathStub);
        
        outPtable = cell2table(outTopGenes);
        
        writetable(outPtable,outP);

        ccEnrichment.topGenes = outTopGenes;
        ccEnrichment.idxGeneSelect = cidxGeneSelect;
    end
    %%
    
    % clMean(isnan(weightFilter)) = nan;
    [D,K] = size(weightFilter);
    
    % Perform GS enrichment
    if isempty(geneSet)
        return;
    end
        
    %% GS enrichment
    outGSfile = sprintf('%s_ccNMF_GStest.mat',outPathStub);   
    
    if exist(outGSfile,'file') && opts.forceRecalc == 0
        fprintf('Loading intermediate GStest results\n');
        load(outGSfile);
        
        geneSetNames = ccEnrichment.outE.geneSetName;

        geneSetSub = ismember(geneSet.setNames,geneSetNames);

        if any(~geneSetSub)
            geneSet = structSubSelectMat(geneSet,geneSetSub);
        end

    else
        if opts.doPerm == 1                
            % TODO: add per cluster mean
            ccEnrichment.outE = calcPermHGenrichment(geneSet,cidxGeneSelect,geneID,clMean,opts);
            
            geneSetNames = ccEnrichment.outE.geneSetName;
            geneSetSub = ismember(geneSet.setNames,geneSetNames);

            if any(~geneSetSub)
                geneSet = structSubSelectMat(geneSet,geneSetSub);
            end
        else
            
            error('Untested -- probably not working')
%             warning('Null model needs some fixing -- the set of geneIDs is filtered');
%             [ccEnrichment.outE,geneSet] = calcHGenrichment(geneSet,outWeight,cidxWord,cVal,geneID,opts);
        end
        %%

        if opts.doTF
            % tfRefSet = geneSet.idxGeneKey;
            ccEnrichment.outE.tfCorr = nan(length(geneSet.setNames),K);    
            [~,zia,zib] = intersect(geneID,geneSet.setNames);

            if ~isempty(zia)
                zDataSub = dataM(zia,:)';

                if opts.useMacc
                    ccEnrichment.outE.tfCorr(zib,:) =  madccMat(outH',zDataSub,opts.useMacc)'; 
                else
                    ccEnrichment.outE.tfCorr(zib,:) =  corr(outH',zDataSub,'type',opts.corrType)'; 
                end
            else
                warning('TF names and geneID have no overlap corellation is not useful');
            end
        end

        if opts.saveIntermediate                                 
            fprintf('Saving intermediate GStest results\n');
            save(outGSfile,'ccEnrichment');                    
        end
    end     
       
    %
    % Write summary table of GS enrichment
    %
    if opts.outGeneSetTable
        %%
        zT = opts.enrichmentThr;
        nG = length(geneID);
        nD = size(outWeight,2);

        outP = sprintf('%s_ccNMF_GSE.xlsx',outPathStub);
        cFDRsig = ccEnrichment.outE.FDR < opts.fdrThr;

        fprintf('Significant sets after FDR correction:\n');
        disp(sum(cFDRsig));


        for zi = 1:nD
            
            cListSig = cFDRsig(:,zi);
            
            if opts.doTF && opts.includeTopTF               
                fprintf('Include top %d TFs\n',opts.includeTopTF);
                [zz,zidxJ] = sort(ccEnrichment.outE.tfCorr(:,zi),'descend','MissingPlacement','last');

                cAddIdx = zidxJ(1:opts.includeTopTF);
                cAddIdx(isnan(1:opts.includeTopTF)) = 0;

                cListSig(cAddIdx) = 1;                    
            end

            if sum(cListSig) == 0
                fprintf('No pathways are signifincant');
                [zz,zMin] = sort(ccEnrichment.outE.pval(:,zi),'ascend','MissingPlacement','last');
                cListSig(zMin(1:10)) = 1;
                isPathwaySig = 0;
            else
                isPathwaySig = 1;
            end          
            
            geneSetSub = structSubSelectMat(geneSet,cListSig);

%             cSelTop = cidxWord(1:zT,zi);
%             cSelTop(cVal(1:zT,zi) <= 0) = [];
%             cSelTop = trueV(cSelTop,nG);

            cSelTop = ccEnrichment.idxGeneSelect(:,zi);
            
            cSelSubE = trueV(zi,nD);
            enSub = structSubSelectMat(ccEnrichment.outE,cSelSubE);
            % enSub = structSubSelectMat(enSub,cListSig);
 %%           
            outPlot = sprintf('%s_ccNMF_GSplot_%02d',outPathStub,zi);
            outSheet = sprintf('cc%02d',zi);
                        
            if sum(cSelTop) == 0
                continue;
            end
            enSub = writePlotGeneSet(enSub,geneSet,geneID,cSelTop,opts);
            
            ccEnrichment.enSub{zi} = enSub;
            writetable(struct2table(enSub),outP,'Sheet',outSheet);
        
            % Plot per topic top GS and top Genes
            if opts.plotTopicEnrichment
                %%
                Hstat = outH(zi,:);
                % Plot Other 
                cSelGene = cidxWord(1:opts.plotTopN,zi);
                outTopG = geneID(cSelGene);
                % outTopV = cVal(1:opts.plotTopN,zi);
                outTopWeight = weightFilter(cSelGene,zi);
                outTopW = outW(cSelGene,zi);
                %
                if opts.doTF
                    zfig = plot_tf(enSub,outPlot,zi,Hstat,ydata,outCl,timeP,outTopG,outTopWeight,opts);
                else
                    zfig = plot_geneSet(enSub,outPlot,zi,Hstat,ydata,outCl,timeP,outTopG,outTopWeight,outTopW,opts);
                end
            end

        end            
    end

end

function [enSub,geneSetSub] = writePlotGeneSet(enSub,geneSet,geneID,cSelTop,opts)

        nTestSig = length(enSub.geneSetName);
        assert(isequal(enSub.geneSetName,geneSet.setNames))

        enSub.fisherP = nan(nTestSig,1);
        enSub.overlapGenes = cell(nTestSig,1);   
        enSub.coverage = enSub.overlap(:)./enSub.setSize(:);

        % Pick the correct sorting procedure
        [enSub,geneSetSub] = sortGStable(enSub,geneSet,opts);

        [~,zia,zib] = intersect(geneID,geneSetSub.idxGeneKey);
        cSetGene = false(length(geneID),length(geneSetSub.setNames));
        cSetGene(zia,:) = (geneSetSub.idxMat(:,zib)>0)';

        nTestSig = length(enSub.geneSetName);

        for zg = 1:nTestSig
            cIdxGene = cSetGene(:,zg);

            if zg<opts.maxTestFisher
                [~,enSub.fisherP(zg)] = fishertest(crosstab(cSelTop,cIdxGene));
            end

            subGene = geneID(cIdxGene & cSelTop);
            enSub.overlapGenes{zg} = strjoin(subGene,',');
        end              

        if opts.useJaccardSetClustering && nTestSig > 10
            jThr = opts.jaccardSetThreshold;                    
            if opts.useJaccardSetClustering == 1
                fprintf('Collapse gene sets by overlap jaccard (set genes only), threshold=%f\n',opts.jaccardSetThreshold);                       
                cSetGeneOverlap = cSetGene(cSelTop,:)';                        
            else
                fprintf('Collapse gene sets by overlap jaccard (Global!), threshold=%f\n',opts.jaccardSetThreshold);                        
                cSetGeneOverlap = cSetGene';                        
            end                  

            [enSub,zidxJ] = ccNMFtopicEnr_collapseByGraphSim(enSub,cSetGeneOverlap,opts);                    
            if opts.jaccardSetBinSort  
                zord = argsort(zidxJ);
                enSub = structSortMat(enSub,zord);  
                geneSetSub = structSortMat(geneSetSub,zord);                        
            end                    
        end                               

end

function [enSub,geneSetSub] = sortGStable(enSub,geneSetSub,opts)

    switch opts.sortBy
        case { 'Corr', 'corr' } 
            fprintf('Sort by TF topic corr \n');
            [~,zidxJ] = sort(enSub.tfCorr,'descend','MissingPlacement','last');
            zidxJ(opts.geneSetTableTopN+1:end) = nan;

            enSub = structSortMat(enSub,zidxJ);
            geneSetSub = structSortMat(geneSetSub,zidxJ);


        case { 'Pval', 'pval' }                        
            fprintf('Sort by Pval \n');
            [zzS,zidxJ] = sort(enSub.pval,'ascend','MissingPlacement','last');
            zidxJ(opts.geneSetTableTopN+1:end) = nan;

            enSub = structSortMat(enSub,zidxJ);
            geneSetSub = structSortMat(geneSetSub,zidxJ);

        case { 'PermP','permp' }

            fprintf('Sort by Pval \n');
            [zzS,zidxJ] = sort(enSub.pval,'ascend','MissingPlacement','last');
            zidxJ(opts.geneSetTableTopN+1:end) = nan;

            enSub = structSortMat(enSub,zidxJ);
            geneSetSub = structSortMat(geneSetSub,zidxJ);
            if opts.doPerm
                fprintf('Sort by Pval \n');
                [zzS,zidxJ] = sort(enSub.permP,'ascend','MissingPlacement','last');

                enSub = structSortMat(enSub,zidxJ);
                geneSetSub = structSortMat(geneSetSub,zidxJ);                             
            end
        case { 'Coverage', 'coverage' }                        
            fprintf('Sort by coverage \n');
            [zzS,zidxJ] = sort(enSub.coverage,'descend','MissingPlacement','last');
            zidxJ(opts.geneSetTableTopN+1:end) = nan;

            enSub = structSortMat(enSub,zidxJ);
            geneSetSub = structSortMat(geneSetSub,zidxJ);                        

        case { 'Jaccard', 'jaccard' }
            error('Not implemented');

        otherwise
            fprintf('Keeping unsorted\n');
    end
    

end

function zfig = plot_tf(enSub,outPlot,zi,Hstat,ydata,outCl,timeP,outTopG,outTopV,opts)

    nSetFound = length(enSub.geneSetName);
    if nSetFound > opts.plotTopN
        enSub = structSubSelectMat(enSub,trueV(1:opts.plotTopN,nSetFound));
        nSetFound = length(enSub.geneSetName);
    end

    zfig = figure('Position',[0 0 2000 1980'],'visible',opts.plotFigViz);                    

    % Plot tSNE                    
    % Hstat = outH(zi,:);
    zopts = opts;
    zopts.cmapLimits = quantile(Hstat,[0.01 0.99]);
    zopts.maxRowPlot = luniq(outCl) + 1;
    %%
    if isempty(timeP)

        zopts.axisObj{1} = subplot(2,4,1:2);
        zopts.axisObj{2} = subplot(2,4,3:4);

        zopts.inLabel = sprintf('GeneProgram - %d',zi);
        zopts.yAxis_unit = '';

        plot_sigvalue_overlay([],[],ydata,outCl,Hstat,zopts);

    else
        %%
        zopts.newPlot = 0;
        subplot(2,4,3:4);
        plot_tsne_scatter(ydata,Hstat,[],[],zopts);

        subplot(2,4,1:2);
        plot_summary_time(Hstat,outCl,timeP,[],zopts);

    end

    zMx = min(length(enSub.geneSetName),opts.plotTopN);
    outTopG = enSub.geneSetName(1:zMx);
    outTopV = enSub.tfCorr(1:zMx);
    outTopV(isnan(outTopV)) = 0;

    cDrop = isnan(outTopV);
    outTopG(cDrop) = [];
    outTopV(cDrop) = [];

    % Plot top genes 
    subplot(2,4,5);
    zg = barh(outTopV);
    set(gca,'yticklabel',outTopG,'ytick',1:length(outTopG));
    xlabel('TF topic corr');

    % Plot top Gene Sets 
    subplot(2,4,8);

    zTopGS = enSub.geneSetName;                    
    zTopGS = regexprep(zTopGS,'_','-');

    if opts.markSig 
        zIsSig = enSub.fisherP*length(enSub.fisherP)<opts.pvalThr;
        if sum(zIsSig) > 0                    
            zTopGS(zIsSig) = regexprep(zTopGS(zIsSig),'(.*)','$1^*');
        end
    end

    zCval = full([ enSub.overlap enSub.setSize ]);

    if opts.plotSetSize 
        barh(zCval,'stacked')
        set(gca,'YTick',1:nSetFound,'YTickLabel',zTopGS);
        xlabel(sprintf('Overlap genes among top %d',opts.enrichmentThr));
        set(gca,'xscale','log');
    else

        barh(zCval(:,1))
        set(gca,'YTick',1:nSetFound,'YTickLabel',zTopGS);
        xlabel(sprintf('Overlap genes among top %d',opts.enrichmentThr));
    end
   
    print(outPlot,'-dpng');      
    close(gcf)
 
end