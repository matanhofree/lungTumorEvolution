function zfig = plot_summary_time(dataQ,groupV,timepoint,isOutlier,inOpts)
    

    defaultOpts.newPlot = 1;
    defaultOpts.minCnt = []; 
    defaultOpts.maxCnt = [];
    defaultOpts.doText = 0;
    defaultOpts.sortGroup = 1;
    defaultOpts.txtPos = 1.1;
    defaultOpts.markOutlier = 0;
    defaultOpts.rfactor = 1.5;
    defaultOpts.doLog = 1;
    
    defaultOpts.pullOutGroup = [];
    
    defaultOpts.ylabel = '';
    defaultOpts.plotSize = [100, 100, 1200, 600];
    defaultOpts.plotViolin = 0;
    defaultOpts.doBox = 1;
    defaultOpts.jitterC = 0.2;
    defaultOpts.title = '';
    defaultOpts.xTickLabelRotation = [];
    defaultOpts.figVisible = 'on';
    
    defaultOpts.lineWidth = 3;
       
    if nargin <= 2 
        isOutlier = [];
    elseif nargin == 4 && isstruct(isOutlier)
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
    
    %%
    [cTime,timeNames]=grp2idx(timepoint);    
    [timeNames,zidx] = sort(timeNames);
    [~,revIdx] = sort(zidx);
    cTime = revIdx(cTime);
    %%
    [cGrp,groupNames]=grp2idx(groupV);
    
%     if opts.sortGroup
%         if isnumeric(groupV)
%             [~,zidx] = sort(str2double(groupNames));
%             outGroupNames = groupNames(zidx);
%         else
%             [outGroupNames,zidx] = sort(groupNames);
%         end
% 
%         zMap = containers.Map(zidx,1:length(zidx));
%         zGrpSort = nanvalues(zMap,cGrp);
%     else
%         zidx = 1:length(groupNames);
%         outGroupNames = groupNames;
%         zGrpSort = cGrp;
%     end
%     
%     if ~isempty(opts.pullOutGroup)
%         isSideGrp = strgrep(outGroupNames,opts.pullOutGroup);
%         zidxNew = [ find(~isSideGrp); find(isSideGrp)];
%         zMap = containers.Map(zidxNew,1:length(zidx));
%         %%
%         zGrpSort = nanvalues(zMap,zGrpSort);
%         outGroupNames = outGroupNames(zidxNew);
%     end 

    if opts.newPlot 
        zfig = figure('Position',opts.plotSize,'visible',opts.figVisible,'color','w'); 
    else
        zfig = gcf;
    end
    % zaxes = gca();
    
    hold on;
    
    for gi = 1:length(groupNames)
        cSel = cGrp == gi;
        
        cStat = full(dataQ(cSel));
        cStat = cStat(:);
        cTimeSub = cTime(cSel);
        
        [cTimeSubList,~,~,~,cTPpos] = fastUnique(cTimeSub);
        
        %%
        cStatS = cellfun(@(x)median(cStat(x)),cTPpos);
        cStatS = cStatS(:);
        cStatSvar = cell2mat(cellfun(@(x)quantile(cStat(x),[ 0.25 0.75 ]),cTPpos,'uniformoutput',0)');
        cStatSvar = cStatSvar-cStatS(:);
        %%
        
        [cTimeSubList,zidx] = sort(cTimeSubList);
        %%
        zf(gi) = plot(cTimeSubList,cStatS(zidx),'-o','Linewidth',opts.lineWidth);
        errorbar(cTimeSubList,cStatS(zidx),cStatSvar(:,1),cStatSvar(:,2))
    end
    
    %%
    legend(zf,groupNames,'interpreter','none')
    
    set(gca,'XTick',1:length(timeNames),'xTickLabel',regexprep(timeNames,'_','-'));
    if opts.xTickLabelRotation 
        set(gca,'xTickLabelRotation',opts.xTickLabelRotation);
    end
    
end 

