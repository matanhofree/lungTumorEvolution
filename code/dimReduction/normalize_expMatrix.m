function X = normalize_expMatrix(X,W,inOpts)

    defaultOpts.doNorm = 2;    
    defaultOpts.verbose = 0;    
    defaultOpts.zipMean.robStdFact = 2;
    defaultOpts.zipMean.doZIP = 0;
    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;

    if ~exist('W','var') || isempty(W)
        W = [];
        % Normalize input data
        switch opts.doNorm
            case 1
                fprintf(1,'Doing rescaling only.\n!!WARNING: no centering!!).\n');
                X = X - min(X(:));
                X = X / max(X(:));
            case 2
                fprintf(1,'Doing rescaling and centering.\n');
                X = X - min(X(:));
                X = X / max(X(:));
                xVmean = nanmean(X,2);
                X = bsxfun(@minus,X,xVmean);
            case 3
                fprintf(1,'Centering only .\n');
                xVmean = nanmean(X,2);
                X = bsxfun(@minus,X,xVmean);
            case 4
                fprintf(1,'Externally defined normalization function on X.\n');
                disp(opts.normFun);
                X = opts.normFun(X,nanmean(X,2));
            case 5
                fprintf(1,'Centering and robust Z using ZIP w/ moment method\n');
                X = X - min(X(:));
                X = X / max(X(:));
                
                [gMean,gStd] = zip_moment(X,2,opts.zipMean);
                
                zFact = opts.zipMean.robStdFact;
                robStd = min(max(gStd,gMean/zFact),gMean*zFact);
                
                X = bsxfun(@minus,X,gMean);
                X = bsxfun(@rdivide,X,robStd);
                
            case 6
                fprintf(1,'Scaling samples\n');
                X = X./max(X);
                X = X./sum(X);
            case 7
                %%
                fprintf(1,'Quantile norm\n');
                Xn = full(X);
                Xn(X==0) = nan;
                Xn = quantilenorm(Xn);
                Xn(X==0) = 0;
                X = sparse(Xn);
                clear Xn;
            case 8
                fprintf(1,'Quantile norm\n');
                Xn = full(X);
                Xn(X==0) = nan;
                Xn = quantilenorm(Xn);
                Xn(X==0) = 0;
                X = sparse(Xn);
                clear Xn;
            case 9
                fprintf(1,'Quantile norm columns right aligned\n');                
                X = quantnormColSimple(X);
            otherwise
                fprintf('No normalization is applied\n');
        end
    else
       % Normalize input data
        switch opts.doNorm 
        case 1
            fprintf(1,'Doing rescaling only.\n!!WARNING: no centering!!).\n');
            X = X - min(X(:));
            X = X / max(X(:));
        case 2
            fprintf(1,'Doing rescaling and weighted centering.\n');
            X = X - min(X(:));
            X = X / max(X(:));   
            W = bsxfun(@rdivide,W,nansum(W,2));
            xVmean = nansum(X.*W,2); 
            X = bsxfun(@minus,X,xVmean); 
        case 3
            fprintf(1,'Weighted centering only .\n');
            W = bsxfun(@rdivide,W,nansum(W,2));
            xVmean = nansum(X.*W,2); 
            X = bsxfun(@minus,X,xVmean); 
        case 4
            fprintf(1,'Externally defined normalization function on X.\n');
            disp(opts.normFun);
            X = opts.normFun(X,nanmean(X,2));             
        otherwise
        end        
    end

end
