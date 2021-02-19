    function outKnnMat = knndist2simMatrix(X,inOpts)


    defaultOpts.doNorm = 0;
    defaultOpts.distType = 'cosine';
    defaultOpts.knnNum = 30;
    defaultOpts.distMPower = 1;
    
    defaultOpts.scaleMaxGreaterOne = 1;
    

    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    disp(opts);


    [D,N] = size(X);
    X = normalize_expMatrix(X,[],opts);

    fprintf(1,'Begin Knn\n');
    zt = tic;
    [idx,outKernal] = knnsearch(X',X','Distance',opts.distType,'K',opts.knnNum+1);    
    %toc(zt);
    fprintf(1,'Done Knn\n');
    
    
    if opts.distMPower > 1
        fprintf('Raising to %d power',opts.distMPower);
        if strcmp(opts.distType,'cosine') || strcmp(opts.distType,'corr')
            outKernal = 1 - (1-outKernal).^opts.distMPower;
        else
            outKernal = outKernal.^opts.distMPower;
        end
    end

    %%
    zt = tic;
    idx(:,1) = [];
    outKernal(:,1) = [];
    idxRow = repmat((1:N)',1,size(idx,2));


    outKernal = outKernal(:);
    maxK = max(outKernal);

    if maxK > 1
        warning('Distance measure results in values greater than 1. Distances will be scaled to the [0,1] range.');
        % outKernal = outKernal - min(outKernal);
        if opts.scaleMaxGreaterOne == 1
            outKernal = outKernal/max(outKernal);
        % else opts.scaleMaxGreaterOne == 1
            % outKernal = ecdf(outKernal);
        else
            error('Not implmented');
        end
    end


    outKnnMat = sparse(idxRow(:),idx(:),1-outKernal,N,N);
    toc(zt)

end