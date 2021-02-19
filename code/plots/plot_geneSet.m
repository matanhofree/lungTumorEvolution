function zfig = plot_geneSet(enSubIn,outPlot,zi,Hstat,ydata,outCl,timeP,outTopG,outTopWeight,outTopW,opts)

    
    nSetFound = length(enSubIn.geneSetName);        
    
    if ~isfield(enSubIn,'permFDR') && isfield(enSubIn,'permP')
        enSubIn.permFDR = mafdr(enSubIn.permP,'BHFDR',1);
    end
    
    if nSetFound > opts.plotTopN
        enSub = structSubSelectMat(enSubIn,trueV(1:opts.plotTopN,nSetFound));
        nSetFound = length(enSub.geneSetName);
    else 
        enSub = enSubIn;
    end
%%
    zfig = figure('Position',[0 0 2000 1980'],'visible',opts.plotFigViz);

    % Plot tSNE                    
    % Hstat = outH(zi,:);
    zopts = opts;
    zopts.cmapLimits = quantile(Hstat,[0.01 0.99]);
    zopts.maxRowPlot = luniq(outCl) + 1;
    
    if isempty(timeP)
%         zopts.axisObj{1} = subplot(4,4,1:2);
%         zopts.axisObj{2} = subplot(4,4,3:4);
%
        zopts.axisObj{1} = subaxis(4,4,1,1,3,1,opts.subAxisOpts{:}); 
        zopts.axisObj{2} = subaxis(4,4,4,1,1,1,opts.subAxisOpts{:});

%
        zopts.inLabel = sprintf('GeneProgram - %d',zi);
        zopts.yAxis_unit = '';
        
        plot_sigvalue_overlay([],[],ydata,outCl,Hstat,zopts);
    else
        zopts.newPlot = 0;
        % subplot(4,4,3:4);
        subaxis(4,4,1,1,3,1,opts.subAxisOpts{:})
        plot_tsne_scatter(ydata,Hstat,[],[],zopts);

        % subplot(4,4,1:2);
        subaxis(4,4,4,1,1,1,opts.subAxisOpts{:})
        plot_summary_time(Hstat,outCl,timeP,[],zopts);
    end

    cDrop = isnan(outTopWeight);
    outTopG(cDrop) = [];
    outTopWeight(cDrop) = [];
    outTopW(cDrop) = [];

    %% Plot top genes
    % subplot(2,4,5);
    subaxis(4,4,1,2,1,3,opts.subAxisOpts{:})
    zg = barh(outTopWeight);
    set(gca,'yticklabel',outTopG,'ytick',1:length(outTopG));
    xlabel('Scaled weight');
    
    subaxis(4,4,2,2,1,3,opts.subAxisOpts{:})
    zg = barh(outTopW);
    set(gca,'ytick',[]);
    xlabel('Raw NMF weight');

    % Plot top Gene Sets
    % subplot(2,4,8);
    subaxis(4,4,4,2,1,3,opts.subAxisOpts{:})

    zTopGS = enSub.geneSetName;
    if isfield(enSub,'sortBin')        
        zTopGS = mergeStringPair(enSub.sortBin,zTopGS);        
    end
    
    if isfield(enSub,'permFDR')
        selSig = enSub.permFDR < opts.permFDR_thr;
        zTopGS(selSig) = mergeStringPair('*%s%s','',zTopGS(selSig));
    end
    
        
    zTopGS = regexprep(zTopGS,'_','-');

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
%%
    print(zfig,outPlot,'-dpng');
    close(gcf)
    
end