function zx = plot_corr_scatter_pairs(cdataA,clabelA,cdataB,clabelB,inOpts)


    defaultOpts.newFigure = 1;
    defaultOpts.doSub = 1;
    defaultOpts.hexres = 30;
    defaultOpts.hexPlot = 0;
    defaultOpts.corrType = 'pearson';
    defaultOpts.xlim = [];
    defaultOpts.ylim = [];
    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    disp(opts);

    if opts.hexPlot == 1 && ~isdeployed()
        addpath('~/regevdata/local/external-matlab/bean-matlab-toolkit');
    end
    
    nRow = size(cdataA,2);
    nCol = size(cdataB,2);
    
    if opts.newFigure == 1 
        zfig = figure('Position',[100 100 2000 2000 ]);
        set(gcf,'color','w');
        set(0,'DefaultAxesFontSize',12);   
    end    
    
   
    zAll = [];
    subplot1(nRow,nCol,'Gap',[0.025 0.025]);
    rc = 1;
    
    for i = 1:nRow
        for j = 1:nCol
        
        

            subplot1(rc);            
            rc = rc + 1;
            
            if opts.hexPlot == 1
                zx{i,j} = hexscatter(cdataA(:,i),cdataB(:,j),'res',opts.hexres);
            else
                zx{i,j} = scatter(cdataA(:,i),cdataB(:,j));
            end
            axis square
            
            if ~isempty(clabelA)
                xlabel(clabelA{i},'interpreter','none');
                ylabel(clabelB{j},'interpreter','none');
            end
            
            if ~isempty(opts.xlim)
                xlim(opts.xlim);
            end                
            if ~isempty(opts.ylim)
                lim(opts.ylim);
            end

            if ~isempty(opts.corrType)
                [cR] = corr(cdataA(:,i),cdataB(:,j),'rows','pairwise','type',opts.corrType);
                yl = ylim();
                xl = xlim();

                text(diff(xl)*0.1+xl(1),diff(yl)*0.9+yl(1),sprintf('R^2=%.2f',cR.^2));
            end
            
            addLine();
        end
    end    
  
    
end