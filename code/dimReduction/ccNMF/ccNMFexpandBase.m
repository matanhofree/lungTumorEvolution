function outExpBase = ccNMFexpandBase(ccNMF,dataFull,listVar,inOpts)

    defaultOpts.randW = 0;
    defaultOpts.solver = 'bp';
    defaultOpts.runFullNMF = 0;    
    defaultOpts.Hnorm = 0;
    defaultOpts.reWeightExpand = 0;
   
    defaultOpts.reWeightExpandDyn = 500;
    defaultOpts.reWeightExpandHard = 100;
    defaultOpts.reWeightExpandAddListVar = 1;
    
     
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts % inOpts;

    [nD,nN] = size(dataFull);
    
    dataM = dataFull(listVar,:);
    meanData = mean(dataFull,2);

    if isfield(ccNMF,'extrapH')
        HList = ccNMF.extrapH;
    else
        HList = { ccNMF.bestConsH };
    end
    
    if isfield(ccNMF,'mergedW')
        WList = ccNMF.mergedW;
    else
        WList = { ccNMF.consW };
    end
    
    if opts.Hnorm == 1
        fprintf('Rescaling H to uni norm before expanding');
    elseif opts.Hnorm == 2
        fprintf('Rescaling H to sum to one before expanding');
    else
        fprintf('No rescaling of H');
    end
    
    if opts.reWeightExpand == 1
        fprintf('Reweight H and W after expanding best on top WMI genes');   
    end
    
    
    for zi = 1:length(ccNMF.testK)
        H = HList{zi};       
        W = WList{zi};        
        nK = size(H,1);
        
        if opts.Hnorm == 1
            H = H./sqrt(mean(H.^2));
        elseif opts.Hnorm == 2
            H = H./sum(H);
        end
        
        if opts.randW == 1
            Wext = rand(nD,nK).*meanData;
            
            Wext(listVar,:) = W;
        elseif opts.randW == 2 
            zHsto = H./sum(H,2);
            Wext = dataFull*zHsto'; 
            
            Wext(listVar,:) = W;
        else
            Wext = [];            
        end
        
        % outExpBase{zi} = nnlsm_nmf(zH', dataFull', Wext',opts.solver)';                        
        %
        Wexp = nnlsm_nmf(H', dataFull', Wext',opts.solver)';
        if opts.reWeightExpand == 1
            
            [~,~,outTopMerged] = ccNMFextractTopGenes(Wexp,[],opts.reWeightExpandDyn,opts.reWeightExpandHard)
            if opts.reWeightExpandAddListVar
                outTopMerged = outTopMerged | listVar;                
            end
            fprintf('Using %d gene in reWeighting\n',sum(outTopMerged));
            
            H = nnlsm_nmf(Wexp(outTopMerged,:), dataFull(outTopMerged,:), H,opts.solver);                    
            if opts.Hnorm == 1
                H = H./sqrt(mean(H.^2));
            elseif opts.Hnorm == 2
                H = H./sum(H);
            end
            
            WexpP = nnlsm_nmf(H', dataFull', Wexp',opts.solver)';
            fprintf('Diff %f\n',norm(Wexp-WexpP,'fro'));
            Wexp = WexpP;            
        end
        
        %%
        if opts.runFullNMF == 1
            %%
            fprintf('Running complete NMF based on CC start');
            [Wfull,Hfull] = bp_nmf(dataFull,nK,'W_init',Wexp,'H_init',H);
            outExpBase{zi} = Wfull;
        else
            
            outExpBase{zi} = Wexp;
        end
        
    end
    
    
end