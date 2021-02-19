function [zfig,outMat] = plot_crosstab_bar(clustMain,clustStack,cmap,inOpts)

    defaultOpts.axisObj =[];
    defaultOpts.clustTxt = 1;
    defaultOpts.XTickLabelRotation = 60;
    defaultOpts.useFrac = 1;
    defaultOpts.doSortX = 1;
    defaultOpts.doSortGroup = 1;
    defaultOpts.doLeafOptimalOrder = 0;
    defaultOpts.doLeafOptimalOrderGroup = [];
    defaultOpts.externXOrder = [];
    defaultOpts.normTotal = 0;
            
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    % disp(opts);
    zfig = [];
   
    outMat = [];
    if isempty(opts.axisObj)     
        zfig = figure('color','w','Position', [100, 100, 2000, 800]);        

        opts.axisObj = gca;    
    end           
  
    if exist('cmap','var') && ~isempty(cmap)
        colormap(cmap)
    end
    
    [zTab,~,~,headerC] = crosstab(clustMain,clustStack);
    
    if opts.normTotal 
        zTab = bsxfun(@rdivide,zTab,sum(zTab,1));
    end
    
    if opts.useFrac
        zTab = bsxfun(@rdivide,zTab,sum(zTab,2));
    end

    xName = headerC(:,1);
    xName = xName(~isemptycell(xName));
    
    if opts.doSortX
        [xName,xOrder] = sort(xName);
        zTab = zTab(xOrder,:);
    elseif ~isempty(opts.externXOrder) 
        fprintf('Using external X order\n');
        xOrder = opts.externXOrder;
        xName = xName(xOrder);
        zTab = zTab(xOrder,:);
    end
    
    zRefName = headerC(:,2);
    zRefName = zRefName(~isemptycell(zRefName));

    if opts.doSortGroup
        [zRefName,refOrder] = sort(zRefName);
        zTab = zTab(:,refOrder);
    end
    %%
    if opts.doLeafOptimalOrder
        
        if isempty(opts.doLeafOptimalOrderGroup)
        
            zD = pdist(zTab);
            zOrdT = linkage(zD,'average');

            leafOrder = optimalleaforder(zOrdT,zD,'Transformation','inverse');

            zTab = zTab(leafOrder,:);
            xName = xName(leafOrder,:);
        else
            zGroup = opts.doLeafOptimalOrderGroup(xName);
            [zGroupList,~,~,zGcnt,zGroupPos] = fastUnique(zGroup);                        
            
            for i = 1:length(zGroupList)
                cidx = zGroupPos{i};
                if zGcnt(i) > 1                
                    zD = pdist(zTab(cidx,:));
                    zOrdT = linkage(zD,'average');

                    cSubOrd = optimalleaforder(zOrdT,zD,'Transformation','inverse');
                    leafOrder{i} = cidx(cSubOrd);
                else
                    leafOrder{i} = cidx;
                end
            end            
            leafOrder = cell2mat(leafOrder);
            
            zTab = zTab(leafOrder,:);
            xName = xName(leafOrder,:);
        end
    end
    %%
    bar(zTab,'stacked')
    zz = gca;
    
    outMat.cTab = zTab;
    outMat.xName = xName;
    
    set(zz,'xtick',1:length(xName));
    set(zz,'xticklabels', regexprep(xName,'_','-'));
    zz.XTickLabelRotation = opts.XTickLabelRotation;

    % zRefName = zHnames % { 'Inflammatory' 'Goblet'  'Enterocytes' 'Transit-Amplifying' 'STEM' };

    
    legend(zRefName,'Location','northeastoutside','box','off','interpreter','none')
    if opts.useFrac
        ylim([ 0 1 ]);
        ylabel('Cell fraction');
    else
        ylabel('Cell count');
    end
    box off
    
    xlabel('Subtype');
    
    if exist('cmap','var') && ~isempty(cmap) 
       
        zz.ColorOrder = cmap;
    end
  
end
