function [figOut,outMat,zOrdTX,zOrdTY] = plot_pairwise_corr(ccMat,cAnnot,inOpts)

    defaultOpts.corrShow = 'pears';
    defaultOpts.useCorr = 1;
    defaultOpts.corrNan = 'pairwise'
    defaultOpts.pdistX = [];
    defaultOpts.pdistY = [];
        
    defaultOpts.position = [ 1 1 1800 1800 ];
    defaultOpts.colorSym = 1;
                
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    % disp(opts);
    zfig = [];
    
    if isempty(opts.pdistX)
        opts.pdistX = opts.corrShow;
    end
    
    if isempty(opts.pdistY)
        opts.pdistY = opts.corrShow;
    end
        
    cSelZero = sum(ccMat ~= 0) == 0;
    if any(cSelZero)
        ccMat(:,cSelZero) = [];
        cAnnot(cSelZero) = [];
    end
    
    [D,N] = size(ccMat);
    fprintf('Running matrix of size (%d,%d)\n',D,N);       
    
    if ~exist('cAnnot','var') || isempty(cAnnot)
        cAnnot = mergeStringPair('Cl',1:N);
    end    
    
    isCorr = 1;
    if strgrep(opts.corrShow,'pears|spe|kend')
        fprintf('Using %s corr\n',opts.corrShow);
        cCorr = corr(ccMat,'type',opts.corrShow,'rows','pairwise');   
    else
        fprintf('Using %s distance\n',opts.corrShow);
        isCorr = 0;
        if strgrep(opts.corrShow,'mahal')
            cDist = pdist(ccMat',opts.corrShow,shrinkage_cov(full(ccMat')));      
            cCorr = (cDist - median(cDist))./iqr(cDist);
            cCorr = squareform(cCorr);
        else
            cDist = pdist(ccMat',opts.corrShow);
            cCorr = squareform(cDist);
        end       
    end
    
%     if any(isnan(cCorr))
%         [i,j] = find(isnan(cCorr))
%         
%         cCorr(i,:) = [];
%         cCorr(:,j) = [];
%         
%         if exist('cDist','var')
%             cDist(i,:) = [];
%             cDist(:,j) = [];
%         end
%         
%         cAnnot(i) = [];
% 
%     end
     
    if opts.useCorr
        if isCorr           
            cDist = squareform(1 - cCorr);
            opts.doLeafOptimalOrderX_data.pdist = cDist;
            opts.doLeafOptimalOrderY_data.pdist = cDist;
        else            
            opts.doLeafOptimalOrderX_data.pdist = cDist;
            opts.doLeafOptimalOrderY_data.pdist = cDist;
        end                
    else
        opts.doLeafOptimalOrderX_data.pdist = pdist(ccMat',opts.pdistX);
    
        if strcmp(opts.pdistX,opts.pdistY)
            opts.doLeafOptimalOrderY_data.pdist = opts.doLeafOptimalOrderX_data.pdist;
        end
    end
    %%
    [figOut,outMat,zOrdTX,zOrdTY] = plot_heatmap_annot(cCorr,cAnnot,cAnnot,[],[],opts)


end