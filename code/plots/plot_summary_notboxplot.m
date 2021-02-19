function zfig = plot_summary_notboxplot(dataQ,groupV,isOutlier,inOpts)
    

    defaultOpts.newPlot = 1;
    defaultOpts.minCnt = []; 
    defaultOpts.maxCnt = [];
    defaultOpts.doText = 0;
    defaultOpts.sortGroup = 1;
    defaultOpts.txtPos = 1.1;
    defaultOpts.markOutlier = 0;
    defaultOpts.rfactor = 1.5;
    defaultOpts.doLog = 1;
   
    defaultOpts.violinParam = [];
    
    defaultOpts.pullOutGroup = [];
    
    defaultOpts.ylabel = '';
    defaultOpts.fpos = [100, 100, 1200, 600];
    defaultOpts.plotViolin = 0;
    defaultOpts.doBar = 0;
    defaultOpts.doBox = 1;
    defaultOpts.jitterC = 0.2;
    defaultOpts.title = '';
    defaultOpts.xTickLabelRotation = 60;
    defaultOpts.figVisible = 'on';
    defaultOpts.maxRowPlot = 50;
    defaultOpts.colorMap = [];
       
    if nargin <= 2 
        isOutlier = [];
    elseif nargin == 3 && isstruct(isOutlier)
        inOpts = isOutlier;
        clear('isOutlier');             
        isOutlier = [];
    end
    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end        
    clear defaultOpts;
            
    [zGrp,GN]=grp2idx(groupV);
    
    if opts.sortGroup
        if isnumeric(groupV)
            [~,zidx] = sort(str2double(GN));
            outGroupNames = GN(zidx);
        else
            [outGroupNames,zidx] = sort(GN);
        end

        zMap = containers.Map(zidx,1:length(zidx));
        zGrpSort = nanvalues(zMap,zGrp);
    else
        zidx = 1:length(GN);
        outGroupNames = GN;
        zGrpSort = zGrp;
    end
    
    if ~isempty(opts.pullOutGroup)
        isSideGrp = strgrep(outGroupNames,opts.pullOutGroup);
        zidxNew = [ find(~isSideGrp); find(isSideGrp)];
        zMap = containers.Map(zidxNew,1:length(zidx));
        %%
        zGrpSort = nanvalues(zMap,zGrpSort);
        outGroupNames = outGroupNames(zidxNew);
    end 
    
    if length(GN) < opts.maxRowPlot
        
        if opts.newPlot
            zfig = figure('color','w','Position', opts.fpos,'visible',opts.figVisible);        
        end
    
        zfig = plot_summary_subplot(dataQ,zGrpSort,outGroupNames,isOutlier,opts);
    else
        
        totalGrp = length(GN);
        
        nRow = ceil(totalGrp/opts.maxRowPlot);
        nPerRow = ceil(totalGrp/nRow);
        
        if opts.newPlot
            opts.newPlot = 0;
            opts.fpos(4) = (opts.fpos(4)-opts.fpos(2))*nRow + opts.fpos(2);
            zfig = figure('color','w','Position', opts.fpos,'visible',opts.figVisible);        
        end
        
        for i = 1:nRow            
            subplot(nRow,1,i);                
            zSel = ((i-1)*nPerRow+1):min(i*nPerRow,totalGrp);
            zIsSub = ismember(zGrpSort,zSel);
            if isempty(isOutlier)
                plot_summary_subplot(dataQ(zIsSub),grp2idx(zGrpSort(zIsSub)),outGroupNames(zSel),[],opts)
            else
                plot_summary_subplot(dataQ(zIsSub),grp2idx(zGrpSort(zIsSub)),outGroupNames(zSel),isOutlier(zIsSub),opts)
            end
        end
    end
    
  
    
end 

function zfig = plot_summary_subplot(dataQ,zGrpSort,outGroupNames,isOutlier,opts)
    zfig = [];

    if min(size(dataQ)) > 1
        dataQ = sum(dataQ>0);               
    end
        
    outG = [];
    
    if opts.plotViolin == 1
        
       distributionPlot(dataQ(:),'groups',zGrpSort,'color',[0.8314    0.7255    0.8549],'showMM',0);
       % distributionPlot(dataQ(:),'groups',zGrpSort,'color',[0.6,0.6,1],'showMM',0,'histOpt',1.1);
       boxData = notBoxPlot(dataQ,zGrpSort,opts.jitterC/2);  
    elseif opts.plotViolin == 2
        
       % distributionPlot(dataQ(:),'groups',zGrpSort,'color',[0.8314    0.7255    0.8549],'showMM',0);
       % distributionPlot(dataQ(:),'groups',zGrpSort,'color',[0.6,0.6,1],'showMM',0,'histOpt',1.1);       
       % boxData = notBoxPlot(dataQ,zGrpSort,opts.jitterC/2); 
       %%
        if ~isempty(opts.violinParam)
            vParams = struct2param(opts.violinParam);
            zfig = violinplot(full(dataQ(:)),zGrpSort,opts.colorMap,vParams{1},vParams{2});
        else 
            zfig = violinplot(full(dataQ(:)),zGrpSort,opts.colorMap);
        end
    elseif opts.doBar == 1
        %% 
        if isnumeric(zGrpSort)
            grpSortTxt = mergeStringPair('Cl%s%d','',zGrpSort);
        else
            grpSortTxt = zGrpSort
        end
        ropts = opts;
        ropts.newPlot = 0;
        ropts.reorderX = 0;
        [zfig,outG] = plot_bar_simple(grpSortTxt,dataQ,[],opts);
       
    elseif opts.doBox == 1
        boxData = notBoxPlot(dataQ,zGrpSort,opts.jitterC); 

        if opts.markOutlier
            boxMarkers = [ boxData.data ];

            hold on;
            for zi = 1:length(boxMarkers)
                xData = boxMarkers(zi).XData;
                yData = boxMarkers(zi).YData;

                if opts.markOutlier == 1
                    [isOutlier] = findOutlierMedianIQR(yData,[],opts.rfactor);
                    % disp(sum(isOutlier))
                elseif opts.markOutlier == 2
                    isOutlier = abs(zscore(yData)) > opts.rfactor;
                end

                plot(xData(isOutlier),yData(isOutlier),'xb');
            end
        end        
    end   
    
    if isempty(outG)
        set(gca,'XTick',1:length(outGroupNames),'xTickLabel',regexprep(outGroupNames,'_','-'));
        ylabel(opts.ylabel);
        
        if exist('isOutlier','var') && ~isempty(isOutlier) % && opts.plotViolin == 1      
            hold on;
            boxMarkers = [ boxData.data ];

            for zi = 1:length(boxMarkers)
                cOutlier = isOutlier(zGrpSort == zi);

                xData = boxMarkers(zi).XData;
                yData = boxMarkers(zi).YData;

                plot(xData(cOutlier),yData(cOutlier),'xb');
            end
            zKeep = isOutlier;        
        else
            zKeep = false(size(dataQ));

            zxl = xlim();
            zyl = ylim();

            if ~isempty(opts.minCnt)
                zKeep = dataQ > opts.minCnt;
                line(zxl,[opts.minCnt opts.minCnt],'lineStyle','--','color','k')
            end

            if ~isempty(opts.maxCnt)
                zKeep = zKeep & dataQ < opts.maxCnt;
                line(zxl,[opts.maxCnt opts.maxCnt],'lineStyle','--','color','k')
            else
                if isempty(opts.minCnt)
                    zKeep(:) = 1;
                end
            end
            opts.maxCnt = zyl(2)*0.9;
        end   
    end

    %% keepMatOrig = fastUnique(groupV,zKeep);    
    if opts.doText
        [groupVlist,~,~,groupCnt,groupCntPos] = fastUnique(outGroupNames(zGrpSort));
        keepCnt = cellfun(@(x)sum(full(zKeep(x))),groupCntPos);    

        % keepMat = keepMatOrig(zidx,:);
        % disp(keepMat);
        if isnumeric(groupVlist)
            groupVlist = num2cellstr(double(groupVlist));
        end
        groupCnt = double(full(groupCnt));
        keepCnt = double(full(keepCnt));
        
        groupCnt = groupCnt(:);
        keepCnt = keepCnt(:);

        [~,zia,zib] = intersect(outGroupNames,groupVlist);
   
        groupCnt(zia) = groupCnt(zib);
        keepCnt(zia) = keepCnt(zib);
    
        if opts.doText == 1          
            txtVect = arrayfun(@(x,y)sprintf('%d/%d',x,y),keepCnt,groupCnt,'uniformoutput',0);
        elseif opts.doText == 2            
            txtVect = arrayfun(@(x)sprintf('%.2f',x),keepCnt./groupCnt,'uniformoutput',0);
        elseif opts.doText == 3
            txtVect = arrayfun(@(y)sprintf('%d',y),keepCnt,'uniformoutput',0);            
        else 
            txtVect = arrayfun(@(x)sprintf('%d',x),keepMat(:,2)','uniformoutput',0);
        end
        
        
        zyl = ylim();

        maxCnt = zyl(2)*0.9;
        text(1:length(outGroupNames),ones(length(outGroupNames),1)*maxCnt*opts.txtPos,txtVect,'fontSize',14,'HorizontalAlignment','center');
    end
    
    if isempty(zfig)
        zfig = gcf;
    end
    if ~isempty(opts.title)
        title(opts.title,'interpreter','none');        
    end
%     if ~isempty(opts.xTickLabelRotation)        
%         zfig.Children(1).XTickLabelRotation = opts.xTickLabelRotation;
%     end
    if ~isempty(opts.xTickLabelRotation)        
        % zfig.Children(1).XTickLabelRotation = opts.xTickLabelRotation;
        zzGCA = gca;
        zzGCA.XTickLabelRotation = opts.xTickLabelRotation;
    end
end