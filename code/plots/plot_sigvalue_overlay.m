function zfig = plot_sigvalue_overlay(zData,zGene,ydata,zClustOrig,overlayValue,inOpts)

    defaultOpts.yAxis_unit = ' log2(TPX)';
    defaultOpts.cmap = 'jet';
    defaultOpts.verbose = 1;
    defaultOpts.inLabel = [];
    defaultOpts.axisObj = [];
    defaultOpts.doScatter = 1;
    defaultOpts.doBox = 1;
    defaultOpts.dropNan = 1;
    defaultOpts.selectCluster = [];
    defaultOpts.maxScatter = [];
    defaultOpts.minScatter = [];  
    defaultOpts.xTickLabelRotation = 60;
    defaultOpts.collapseFunc = @(x)nansum(x,1);
    
    defaultOpts.minCnt = []; 
    defaultOpts.maxCnt = [];
    defaultOpts.doText = 0;
    defaultOpts.txtPos = 1.05;
    defaultOpts.clustTxt = 0;
    defaultOpts.markOutlier = 0;
    defaultOpts.pSize = 30;
    defaultOpts.cmapLimits = [];
    defaultOpts.plotSize = [1, 1, 2400, 800];
    defaultOpts.showPlot = 1;
    defaultOpts.plotViolin = 2;
    defaultOpts.selectV = [];
    
            
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    % disp(opts);
    zfig = [];

    clustNames = [];
    if iscellstr(zClustOrig)
        [zClust,clustNames] = grp2idx(zClustOrig);
    else
        zClust = zClustOrig;
        clustNames = num2cellstr(unique(zClustOrig));
    end


    zDrop = [];
    if opts.dropNan
        zDrop = isnan(zClust) | zClust < 0;

        zData(:,zDrop) = [];
        ydata(zDrop,:) = [];
        zClust(zDrop) = [];
        zClustOrig(zDrop) = [];
    end
    

    if ischar(overlayValue) 
        zGname = overlayValue;
        % zi = find(strgrepi(zGene,overlayValue));
        zi = find(strcmp(zGene,overlayValue));        
        if (isempty(zi));            
            fprintf('Gene name not found. Try:\n');
            disp(zGene(strgrepi(zGene,overlayValue)));
            return;
        end
           
        if numel(zi) > 1
            fprintf('Found multiple:\n');
            disp(zGene(zi));
            
            zVal = opts.collapseFunc(zData(zi,:));
        else                
            zi = zi(1);
            disp(zGene(zi));                    
            zVal = zData(zi,:);
        end
    else         
        zGname = inputname(5);
        zVal = overlayValue;
        zVal(zDrop) = [];
    end

    if ~isempty(opts.selectCluster)
        if iscell(opts.selectCluster)
            opts.selectCluster = find(ismember(clustNames,opts.selectCluster));
        end
        zDropClust = ~ismember(zClust,opts.selectCluster);
        
        zData(:,zDropClust) = [];
        ydata(zDropClust,:) = [];
        zClust(zDropClust) = [];
        zClustOrig(zDropClust) = [];
        zVal(zDropClust) = [];        
    end

    if ~isempty(opts.inLabel)
        zGname = opts.inLabel;
    end
    
    if isempty(opts.axisObj)
        if (opts.doBox && opts.doScatter)
            if opts.showPlot == 0
                zfig = figure('Position', opts.plotSize,'visible','off');        
            else
                zfig = figure('color','w','Position', opts.plotSize);        
            end
                
            opts.axisObj{1} = subplot(1,2,1);
            opts.axisObj{2} = subplot(1,2,2);
        elseif opts.doBox || opts.doScatter
            zfig = figure('Position',[ 10 10 900 800] ); 
        else
            error('Need to plot either a box or a scatter');
        end    
    end           
  
    
    %%
    if opts.doBox == 1
%%
        if ~isempty(opts.axisObj{1})
            if isempty(zfig)
                zfig = gcf;
            end
            zfig.CurrentAxes = opts.axisObj{1};
            % axes(opts.axisObj{1});
        end
        
        opts.newPlot = 0;
        
        if ~isempty(opts.selectV)
            cSel = opts.selectV;
            plot_summary_notboxplot(full(zVal(cSel)),zClustOrig(cSel),[],opts);
        else
            plot_summary_notboxplot(full(zVal),zClustOrig,[],opts);   
        end
        

        % zClustOrig = regexprep(zClustOrig,'[0-9]_','');
        zPos = get(gca,'OuterPosition');
        zPos(2) = zPos(2) + 0.05;
        zPos(4) = zPos(4)*0.95;        
        set(gca,'OuterPosition',zPos);

        zz = gca();
        zz.XTickLabel = regexprep(zz.XTickLabel,'^[0-9]-','');
        zz.XTickLabelRotation = opts.xTickLabelRotation;
    
%     if ~isempty(clustNames)
%         set(gca,'xticklabel',clustNames,'TickLabelInterpreter','none');
%     end
        ylabel(sprintf('%s %s',zGname,opts.yAxis_unit),'interpreter','none')
        
    else 
        fprintf('No box?');
    end
    
    if opts.doScatter
        
        if ~isempty(opts.axisObj)
            if isempty(zfig)
                zfig = gcf;
            end
            zfig.CurrentAxes = opts.axisObj{2};
            
            % axes(opts.axisObj{2});
            if (~isempty(opts.axisObj{1}))
                opts.axisObj{1}.Position(3) = opts.axisObj{1}.Position(3) + 0.05;
                opts.axisObj{2}.Position(3) = opts.axisObj{2}.Position(3) + 0.05;
            end
        end
        
        if ~isempty(opts.maxScatter)
            zVal = min(zVal,opts.maxScatter);
        end
        
        if ~isempty(opts.minScatter)
            zVal = max(zVal,opts.minScatter);
        end
        cSel = [];
        if ~isempty(opts.selectV)
            cSel = opts.selectV;
        end
        
        opts.newPlot = [];
        opts.isCategory = 0;
        % opts.pSize = ;
        opts.yAxis_unit = sprintf('%s %s',zGname,opts.yAxis_unit);
        zz = plot_tsne_scatter(ydata,zVal,[],cSel,opts);
        
        
        if (opts.clustTxt) && length(unique(zClust)) < 50
            
            meanCt_x = arrayfun(@(x)median(ydata(zClust == x,1)),unique(zClust));
            meanCt_y = arrayfun(@(x)median(ydata(zClust == x,2)),unique(zClust));
            % text(meanCt_x,meanCt_y,num2cellstr(unique(zClust)),'fontsize',30,'HorizontalAlignment','center');
            
            disp(clustNames);
            zz = text(meanCt_x,meanCt_y,clustNames,'fontsize',30,'HorizontalAlignment','center','interpreter','none');
        end
        
%        [~,zord] = sort(zVal);       
%         zf = scatter(ydata(zord,1),ydata(zord,2),60,zVal(zord),'filled');
%         
%          colormap(opts.cmap);
%          c = colorbar();
% %         ax = opts.axisObj{2};
% %         c.Position(3) = c.Position(3)*0.5;
% %         ax.Parent.Children(2).Position(3) = ax.Parent.Children(2).Position(3)*1.001;
% sss
%         %%
%         xlabel('tSNE 1');
%         ylabel('tSNE 2');
%         ylabel(c,sprintf('%s %s',zGname,opts.yAxis_unit),'interpreter','none');
    end
end
