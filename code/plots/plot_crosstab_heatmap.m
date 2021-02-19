function [figOut,outMat,zOrdTX,zOrdTY] = plot_crosstab_heatmap(clustMain,clustStack,cmap,inOpts)

    defaultOpts.axisObj =[];
    defaultOpts.clustTxt = 1;
    defaultOpts.XTickLabelRotation = 60;
    defaultOpts.useFrac = 1;
    defaultOpts.doSortX = 1;
    defaultOpts.doSortGroup = 1;
    defaultOpts.doLeafOptimalOrder = 0;
    defaultOpts.doLeafOptimalOrderX = 0;
    defaultOpts.externXOrder = 0;
    defaultOpts.colormap = flipud(cptcmap('ylwhbl_sym'));
    defaultOpts.maxColorLim = 0.95;
    defaultOpts.fontSize = 18;
    defaultOpts.axisFontSize = 18;
    defaultOpts.leafDistTrans = 'linear';
    defaultOpts.addDendrogram = 0;
    defaultOpts.linkageF = 'average';
    defaultOpts.pdist = 'Spearman';
    defaultOpts.pseudoCount = 1;
    defaultOpts.colorType = 'pearson';
    
    defaultOpts.colorHard = [];
    
    defaultOpts.xRef = [];
    defaultOpts.doPval = 1;
            
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    % disp(opts);
    zfig = [];
   
    if isempty(opts.axisObj)     
        zfig = figure('color','w','Position', [500 0 1700 1300]);        

        opts.axisObj = gca;    
    end           
  
    if exist('cmap','var') && ~isempty(cmap)
        colormap(cmap)
    end
    
    [crossTable,~,pStat,headerC] = crosstab(clustMain,clustStack);
%     if opts.doPval == 2
%         pStat = fishertest(crossTable)
%     end
%     
%     if opts.useFrac
%         crossTable = bsxfun(@rdivide,crossTable,sum(crossTable,2));
%     end

    xName = headerC(:,1);
    xName = xName(~isemptycell(xName));
     
    if opts.doSortX
        if isnumeric(clustMain)
            [~,xOrder] = sort(str2double(xName));
        else
            [~,xOrder] = sort(xName);
        end
        xName = xName(xOrder);
        crossTable = crossTable(xOrder,:);
    elseif opts.externXOrder 
        fprintf('Using external X order\n');
        xOrder = opts.externXOrder;
        xName = xName(xOrder);
        crossTable = crossTable(xOrder,:);
    end
    
    zRefName = headerC(:,2);
    zRefName = zRefName(~isemptycell(zRefName));

    if opts.doSortGroup
        if isnumeric(clustStack)
            [zRefName,refOrder] = sort(str2double(zRefName));
        else
            [zRefName,refOrder] = sort(zRefName);
        end
        
        crossTable = crossTable(:,refOrder);
    end
    %%    
    pTable = crossTable + opts.pseudoCount;
    cTotal = sum(pTable(:));
    pTableN = pTable./cTotal;
    
    cExp = (sum(pTableN,2)*sum(pTableN))*cTotal;
    
    switch opts.colorType
        case 'pearson'
            cPR = (pTable-cExp)./sqrt(cExp);
        case 'logOdds'
            cPR = log2(pTable./cExp);
        otherwise 
            error('Pick coloring by - pearson or logOdds');                
    end
    

    %%
    % Alternative sort method 
    if opts.doLeafOptimalOrder
        if opts.doLeafOptimalOrder == 1 
            zDY = pdist(crossTable./sum(crossTable,2),opts.pdist);
            zOrdTY = linkage(zDY,opts.linkageF);
        elseif opts.doLeafOptimalOrder == 2 
            zDY = pdist(cPR,opts.pdist);
            zOrdTY = linkage(zDY,opts.linkageF);            
        end
        
        leafOrderY = optimalleaforder(zOrdTY,zDY,'Transformation',opts.leafDistTrans);
        cPR = cPR(leafOrderY,:)
        crossTable = crossTable(leafOrderY,:);
        xNameO = xName(leafOrderY);
    else
        xNameO = xName;
    end
        
    if opts.doLeafOptimalOrderX
        if opts.doLeafOptimalOrderX == 1 
            zDX = pdist((crossTable./sum(crossTable,1))',opts.pdist);
            zOrdTX = linkage(zDX,opts.linkageF);
        elseif opts.doLeafOptimalOrderX == 2 
            zDX = pdist(cPR',opts.pdist);
            zOrdTX = linkage(zDX,opts.linkageF);            
        end
        
        leafOrderX = optimalleaforder(zOrdTX,zDX,'Transformation',opts.leafDistTrans);
        cPR = cPR(:,leafOrderX);
        crossTable = crossTable(:,leafOrderX);
        zRefNameO = zRefName(leafOrderX);
    else
        zRefNameO = zRefName;
    end
        
   
%%
    if opts.maxColorLim 
        if opts.maxColorLim > 1
            clim = opts.maxColorLim;
        else
            clim = quantile(abs(cPR(:)),opts.maxColorLim)            
        end
    else
        clim = max(abs(cPR(:)));
    end
    if ~isempty(opts.colorHard)
        clim = opts.colorHard;
    end
   
    %%
    if ~isempty(opts.xRef)
        [~,zia,zib] = intersect(opts.xRef,zRefNameO);
        
        cPRtmp = zeros(size(cPR,1),length(opts.xRef));
        crossTableTmp = zeros(size(cPR,1),length(opts.xRef));
        
        cPRtmp(:,zia) = cPR(:,zib)
        crossTableTmp(:,zia) = crossTable(:,zib);
        
        zRefNameO = opts.xRef;
        
        cPR = cPRtmp;
        crossTable = crossTableTmp;
    end
        %%
    
    [hImage, hText, hXText] = heatmapTXT(cPR,zRefNameO, xNameO, crossTable, 'TickAngle', opts.XTickLabelRotation,...
        'ShowAllTicks', true, 'TickFontSize', opts.axisFontSize,'FontSize',opts.fontSize,'Colormap',opts.colormap,'Colorbar',1,'mincolorvalue',-clim,'maxcolorvalue',clim);
%%    
    if opts.addDendrogram  > 0
        %%
        
        zHMpos = zfig.Position;
        zHMpos(1) = zHMpos(1)-300;
        zHMpos(3) = 300;
        zf = figure('Position',zHMpos);
        zfd = dendrogram(zOrdTY,0,'Orientation','left','Reorder',flipud(leafOrderY(:)),'Labels',(regexprep(xName(:),'_','-')));
        axis tight        
    end    
    if opts.addDendrogram > 2
        
        zHMpos = zfig.Position;
        zHMpos(2) = zHMpos(4)+zHMpos(1);
        zHMpos(4) = 400
        
        %
        zfX = figure('Position',zHMpos);
        zfdX = dendrogram(zOrdTX,0,'Orientation','Top','Reorder',(leafOrderX(:)),'Labels',(regexprep(zRefName(:),'_','-')));
        axis tight        
        
        zfX.Children.XTickLabelRotation = 90;%opts.XTickLabelRotation;
    end
 
    if opts.doPval == 1
        
        title(sprintf('\\chi^2 P-val = %.3g',pStat));
    elseif opts.doPval == 2
        [~,pStatF] = fishertest(crossTable);
        title(sprintf('Fisher P-val = %.3g',pStatF));
    elseif opts.doPval == 3
        [~,pStatF] = fishertest(crossTable);
        title(sprintf('-Log(Fisher P-val) = %g',-log(pStatF)));
        
    end
        
    
    
    
    if nargout > 1
        if ~iscell(zRefNameO)
            zRefNameO = matlab.lang.makeValidName(num2cellstr(zRefNameO));
        end
        %%
        if ~iscell(xNameO)
            xNameO = matlab.lang.makeValidName(num2cellstr(xNameO));
        end


        outMat.cTable = array2table(crossTable,'VariableNames',zRefNameO','rownames',xNameO);

        %%
        outMat.cPR = array2table(cPR,'variablenames',zRefNameO, 'rownames',xNameO);

        outMat.refName = zRefName;
        outMat.xName = xName;
    end
    %%
    plot_dump([mfilename() '_main'],zfig);
    
    figOut{1} = zfig;
    zi = 2;
    if exist('zf','var')     
        figOut{zi} = zf; zi = zi + 1;
        plot_dump([mfilename() '_dY'],zf);
    end
    if exist('zfX','var')     
        figOut{zi} = zfX; zi = zi + 1;
        plot_dump([mfilename() '_dX'],zf);
    end
    
    
end

%     bar(crossTable,'stacked')
%     zz = gca;
%     
% 
%     set(zz,'xtick',1:length(xName));
%     set(zz,'xticklabels', regexprep(xName,'_','-'));
%     zz.XTickLabelRotation = opts.XTickLabelRotation;
% 
%     % zRefName = zHnames % { 'Inflammatory' 'Goblet'  'Enterocytes' 'Transit-Amplifying' 'STEM' };
% 
%     
%     legend(zRefName,'Location','northeastoutside','box','off','interpreter','none')
%     if opts.useFrac
%         ylim([ 0 1 ]);
%         ylabel('Cell fraction');
%     else
%         ylabel('Cell count');
%     end
%     box off
%     
%     xlabel('Subtype');    

%     %%
%     h = heatmap(zRefName,xName,cPR,'Colormap',opts.colormap);
%     
%     clim = max(abs(h.ColorLimits));
%     h.ColorLimits = [ -clim clim ];
