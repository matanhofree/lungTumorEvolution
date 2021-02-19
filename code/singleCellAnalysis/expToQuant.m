function [eData,outS] = expToQuant(cData,inOpts)
    
    defaultOpts.bkgNumBins = 10;
    defaultOpts.pb = 1;
    defaultOpts.doSparse = 0;
           
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    disp(opts);            
    
    [D,N] = size(cData);
    if opts.doSparse & isparse(cData)        
        eData = sparse([],[],D,N,nnz(cData));
    else
        eData = zeros(D,N);
    end
    
    
    nBin = opts.bkgNumBins;
    
    nM = zeros(D,opts.bkgNumBins);
    nE = zeros(D,opts.bkgNumBins+1);
    
    if opts.pb == 1
        progressbar()
    end
    
    for i = 1:D
        
        rowX = cData(i,:);
        [zN,zE,zC] = histcounts(rowX,nBin);
        
        nM(i,1:length(zN)) = zN;
        nE(i,1:length(zE)) = zE;
        eData(i,:) = zC-1;
        if opts.pb == 1
            progressbar(i/D);
        end
    end
    
    outS.countN = nM;
    outS.bEdge = nE;
    
end