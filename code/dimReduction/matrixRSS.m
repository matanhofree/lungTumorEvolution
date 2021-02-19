function [outRSS,fNorm] = matrixRSS(dataM,W,H,fNorm,inOpts)

    defaultOpts.blockMatrixCalc = 0;
    defaultOpts.nnzOnly = 1;
    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end

    if opts.blockMatrixCalc > 0
        if opts.nnzOnly == 1
            error('Not implemented');
        end
        RSS = 0
        bl = 1;
        N = size(H,2);
        while bl < N
            sel = bl:min(bl+opts.blockMatrixCalc-1,N);
            RSS = RSS + norm(dataM(:,sel) - W*H(:,sel),'fro').^2;    
            bl = bl + opts.blockMatrixCalc;
        end
    else        
        if opts.nnzOnly == 1
            cSel = dataM > 0;
            RSSmat = W*H;
            RSSmat = dataM(cSel) - RSSmat(cSel);
            RSS = sum(RSSmat(:).^2);                        
        else
            RSS = norm(dataM - W*H,'fro').^2;
        end         
    end

    if isempty(fNorm)         
        if opts.nnzOnly == 1
            fNorm = sum(dataM(cSel).^2);
        else
            fNorm = norm(dataM,'fro').^2;
        end
    end
    
    outRSS = full(1 - (RSS/fNorm));
end