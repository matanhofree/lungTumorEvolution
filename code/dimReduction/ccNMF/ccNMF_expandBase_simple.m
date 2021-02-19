function [WexpBase,HlistOut,cErr,outListVarOut] = ccNMF_expandBase_simple(ccNMF,Wlist,Hlist,dataFull,listVar,listVarExclude,inOpts)

    defaultOpts.randW = 0;
    defaultOpts.solver = 'bp';
    defaultOpts.runFullNMF = 0;  
    
    defaultOpts.Hnorm = 0;
    defaultOpts.Wnorm = 0;
    
    defaultOpts.reWeightExpand = 0;
   
    defaultOpts.reWeightExpandDyn = 500;
    defaultOpts.reWeightExpandHard = 100;
    defaultOpts.reWeightExpandAddListVar = 1;
    defaultOpts.pb = 1;
    defaultOpts.calcRss = 1;
    defaultOpts.forceEmbedIncludeAll = 0;
    defaultOpts.mergedIsForced = [];
    
     
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts % inOpts;

    [nD,nN] = size(dataFull);
    
    dataM = dataFull(listVar,:);
    meanData = mean(dataFull,2);
    
    if ~isempty(ccNMF)
        if isfield(ccNMF,'extrapH')
            Hlist = ccNMF.extrapH;
        else
            Hlist = { ccNMF.bestConsH };
        end

        if isfield(ccNMF,'mergedW')
            Wlist = ccNMF.mergedW;
        else
            Wlist = { ccNMF.consW };
        end           
    end
    
    if ~exist('Hlist','var')
        Hlist = [];
    end
        
    if opts.Hnorm == 1
        fprintf('Rescaling H to unit norm before expanding\n');
    elseif opts.Hnorm == 2
        fprintf('Rescaling H to sum to one before expanding\n');
    else
        fprintf('No rescaling of H\n');
    end
    
    if opts.Wnorm == 1
        fprintf('Rescaling W to unit norm before expanding\n');
    else
        fprintf('No rescaling of W\n');
    end
    
    if opts.reWeightExpand == 1
        fprintf('Reweight H and W after expanding best on top WMI genes');   
    end
    
    cN = max(length(Wlist),length(Hlist));
    fNorm = [];
    cErr = nan(cN,1);
    % outListVar = listVar;
    outListVarOut = [];
    HlistOut = [];
    WexpBase = [];
    
    isW = 0;
    if ~isempty(Wlist)
        isW = 1;
    end
    
    isH = 0;
    if ~isempty(Hlist)
        isH = 1;
    end
    
    if isH == 0 && isW == 0
        error('Unable to expand base');
    end
    
    cMergedForced = [];
    if ~isempty(opts.mergedIsForced) && opts.forceEmbedIncludeAll == 1
        cMergedForced = opts.mergedIsForced;
    end
    
    if opts.pb && ismac()
        progressbar();
    end
    for zi = 1:cN
        zt = tic;
        if isW
            W = Wlist{zi};        
            nK = size(W,2);
        end        
        
        if ~isempty(Hlist)
            H = Hlist{zi};
            nK = size(H,1);
        end
        
        
        if ~isempty(Hlist)
            fprintf('Using input H\n');
            
            if opts.Hnorm == 1                
                H = H./sqrt(sum(H.^2));
            elseif opts.Hnorm == 2
                H = H./sum(H);
            end
        else
            if opts.Wnorm == 1
                wNorm = sqrt(sum(W.^2));
                W = W./wNorm;

                if ~isempty(H)
                    H = H.*(wNorm');
                end
            end
            
            fprintf('Reconstructing H\n');
            H = nnlsm_nmf(W, dataM, [],opts.solver);  
        end
        
        if opts.randW == 1
            Wext = rand(nD,nK).*meanData;
            
            Wext(listVar,:) = W;
            
            if opts.Wnorm == 1
                wNorm = sqrt(sum(Wext.^2));
                Wext = Wext./wNorm;
                        
                H = H.*(wNorm');
            end
            
        elseif opts.randW == 2 
            zHsto = H./sum(H,2);
            Wext = dataFull*zHsto'; 
            
            Wext(listVar,:) = W;
            
            if opts.Wnorm == 1
                wNorm = sqrt(sum(Wext.^2));
                Wext = Wext./wNorm;
                        
                H = H.*(wNorm');
            end
        else
            Wext = [];            
        end
        
        % outExpBase{zi} = nnlsm_nmf(zH', dataFull', Wext',opts.solver)';                        
        %        
        Wexp = nnlsm_nmf(H', dataFull', Wext',opts.solver)';
        
        if opts.Wnorm == 1
            wNorm = sqrt(sum(Wexp.^2));
            Wexp = Wexp./wNorm;
                        
            H = H.*(wNorm');
        end
        
        if opts.reWeightExpand == 1
            
            [~,~,outListVar] = ccNMFextractTopGenes(Wexp,[],opts.reWeightExpandDyn,opts.reWeightExpandHard);
           
            if opts.reWeightExpandAddListVar
                outListVar = outListVar | listVar;                
            end
            fprintf('Using %d gene in reWeighting\n',sum(outListVar));
            
            H = nnlsm_nmf(Wexp(outListVar,:), dataFull(outListVar,:), H,opts.solver);                    
            if opts.Hnorm == 1
                H = H./sqrt(mean(H.^2));
            elseif opts.Hnorm == 2
                H = H./sum(H);
            end
            
            outListVarOut{zi} = outListVar;
            
            WexpP = nnlsm_nmf(H', dataFull', Wexp',opts.solver)';
            fprintf('Diff %f\n',norm(Wexp-WexpP,'fro'));
            Wexp = WexpP;            
        end
        
        %%
        if opts.runFullNMF == 1
            %%
            fprintf('Running complete NMF based on CC start');
            [Wexp,H] = bp_nmf(dataFull,nK,'W_init',Wexp,'H_init',H);                       
        
            if opts.forceEmbedIncludeAll == 1
                fprintf('Re-running complete NMF after fixing W matrix with forced input');
                
                fixDim = cMergedForced{zi};
                assert(length(fixDim) == nK,'Missmatch betweek W dim and fixed position indicator variable');
                
                Wexp(listVar,fixDim) = Wlist{zi}(:,fixDim);
                if opts.Wnorm == 1
                    wNorm = sqrt(sum(Wexp.^2));
                    Wexp = Wexp./wNorm;
                        
                    H = H.*(wNorm');
                end
                
                copts = opts.bp_nmf;
                copts.fixDim = fixDim;
                [Wexp,H] = bp_nmf_opt_fsubset(cDataM,nK,Wexp,H,copts);
                
                Wexp(listVar,fixDim) = Wlist{zi}(:,fixDim);
                if opts.Wnorm == 1
                    wNorm = sqrt(sum(Wexp.^2));
                    Wexp = Wexp./wNorm;
                        
                    H = H.*(wNorm');
                end
                H = nnlsm_nmf(Wexp, cDataM,H,opts.solver);  
            end        
        end
        
        WexpBase{zi} = Wexp;
        HlistOut{zi} = H;
        
        if opts.pb && ismac()
            progressbar(zi/cN);
        end        
        if opts.calcRss
            [cErr(zi),fNorm] = matrixRSS(dataFull,Wexp,H,fNorm);
            fprintf('%d) Completed time %g, var exp %f\n',zi,toc(zt),cErr(zi));
        else
            fprintf('%d) Completed time %g\n',zi,toc(zt));
        end
        
    end
    
    if isempty(outListVarOut)
        outListVarOut = listVar;
    end
end