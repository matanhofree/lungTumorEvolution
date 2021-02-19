function [quantData,quantBin,isZero] = mapToQuantileFix(zData,quantLims,inOpts)
    
    defaultOpts.numBins = 5;
    defaultOpts.dropZeroQuant = 1;
    defaultOpts.addNnz = 0;
    
    defaultOpts.useHistC = 0;
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    
    %%
    if exist('quantLims','var') == 0 || isempty(quantLims)
        quantLims = 1/opts.numBins;        
        quantLims = quantLims:quantLims:(1-quantLims);
    end
    
%     if opts.dropZeroQuant
%         isZero = zData == 0;
%         zData(isZero) = nan;
%     end
%%       
    
    [nD,nN] = size(zData);
    if issparse(zData)        
        nV = nnz(zData);
        quantData = sparse([],[],[],nN,nD,nV);
    else
        quantData = zeros(nN,nD);        
    end
    
    quantBin = zeros(nD,length(quantLims));
    isZero = full(all(zData==0,2));
    
    if opts.useHistC
        for i = 1:nD
            cD = zData(i,:);
            if opts.dropZeroQuant
                isZero = cD == 0;
                cD(isZero) = nan;
            end
            [~,quantBin(i,:),quantData(:,i)] = histcounts(cD(i,:),opts.numBins);           
        end
    else
        %%
        if opts.dropZeroQuant
            progressbar()
            for i = 1:nD
                gData = full(nonzeros(zData(i,:)));

                cD = (gData);
                cD = cD(:)';
                quantBin(i,:) = quantile(cD,quantLims,2);
                
                progressbar(i/nD)

            end
        else
            quantBin = quantile(zData,quantLims,2);
        end
        
        progressbar()
        quantBin = quantBin';
        for i = 1:nD  
            gData = full(zData(i,:));

            qB = quantBin(:,i);
            if any(isnan(qB))
                fprintf('Skipping %d\n',i);                
                disp(qB);
                isZero(i) = 1;
                continue;
            end
            
            
            % [~,~,quantData(i,:)] = histcounts(gData,quantBin(i,:));
            quantData(:,i) = sum(gData>qB);
            
            %%
            % [~,quantData(i,:)] = histc(full(zData(i,:)),quantBin(i,:));
            progressbar(i/nD)
        end
        
%         if opts.addNnz == 1
%             quantData(zData > 0) = 1;
%         end

    end
    
    quantData = quantData';
end