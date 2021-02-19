function [expO,freqO] = quantnormColSimple(expMat,batchID,inOpts)

        
    defaultOpts.outExp = 1;
    defaultOpts.mergeExpVect = 1;
    defaultOpts.mergeExpGrid = 1000;
    defaultOpts.sparseOut = 0; 
    
    defaultOpts.meanFunc = 3; % 0 Mean
                              % 1 Median
                              % 3 Quantile 0.75      
    defaultOpts.qDim = 0.75;
    defaultOpts.normalizeZeroRate =  0;
    defaultOpts.centerFreqDist = 0;
    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;  
    disp(opts)    
    
    switch opts.meanFunc
        case 0
            mFunc = @(x)mean(x,1);
        case 1
            mFunc = @(x)median(x,1);
        case 2 
            mFunc = @(x)quantile(x,opts.qDim,1);
        otherwise
            mFunc = @(x)mean(x,1);
    end
          
    [zD,nA] = size(expMat);
    nnzA = nnz(expMat);

%  
%     if nA ~= length(batchID)
%         error('Batchid does not match');
%     end
    
    if opts.sparseOut
        freqO = sparse(zD,nA,nnzA);
        expO = sparse(zD,nA,nnzA);
    else
        freqO = zeros(zD,nA);
        expO = zeros(zD,nA);
    end  
    
%     
%     [batchList,~,nIdx,bCnt,bPos] = fastUnique(batchID);            
%     bVect = cellfun(@(x)trueV(x,nA),bPos,'uniformoutput',0);
%        
%     nBatch = length(bVect);
%     %%
%     mergeGrid = linspace(0, 1, opts.mergeExpGrid);
%     %%
% 
%     doExp = 0;
%     if opts.outExp && nargout > 1
%         doExp = 1;
%         if opts.mergeExpVect
%             mergeGridMat = zeros(nBatch,length(mergeGrid));
%         end
%     end
%     doExpMerge = opts.mergeExpVect;
%     
%     

   mergeGrid = linspace(0, 1, opts.mergeExpGrid);
   progressbar();
   zt = tic();    
   mergeGridMat = zeros(nA,length(mergeGrid));
       
   for zi = 1:nA
        
        bFull = full(expMat(:,zi));                         

        rankB = tiedrank(bFull);
        fB = rankB - min(rankB);
        fB = fB/max(max(fB),1);
        
        [eB,zidx] = unique(bFull);
        xfB = fB(zidx);
        xB = bFull(zidx);  
        mergeGridMat(zi,:) = max(interp1(xfB,xB,mergeGrid,'linear','extrap'),0); 
       
%         % cFreq(1,bVect{zb}) = fB;
%         % fBsplit{zb} = fB;
% 
%         if doExp && doExpMerge
%             [eB,zidx] = unique(bFull);
% 
% 
%             if length(eB) > 1
%                 xfB = fB(zidx);
%                 xB = bFull(zidx);  
% 
%                 nnzFreq = nonzeros(xfB);
%                 batchZeroFreq(zb) = min(nnzFreq);
%                 batchMeanFreq(zb) = mean(nonzeros(fB));
%                 mergeGridMat(zb,:) = max(interp1(xfB,xB,mergeGrid,'linear','extrap'),0); 
%             else
%                 if (eB ~= 0)
%                     warning('Unique nnz expB (%d)',eB);
%                 end
% 
%                 mergeGridMat(zb,:) = eB;
%             end
%         end
            
        freqO(:,zi) = fB;
        progressbar(zi/nA);
   end
   
   %%
   mergeGridMatSum = mFunc(mergeGridMat);   
   selNnz = freqO > 0;
   
   expO(selNnz) = max(interp1(mergeGrid,mergeGridMatSum,freqO(selNnz),'linear','extrap'),0); 
   
   
%    for zi = 1:nA
%        
%    end                              
%         if doExp            
%             if doExpMerge == 0
%                 error('Not implemented');
%                 % expO(zi,:) = interp1(xfA,xA,fA,'linear','extrap');                         
%                 %expO(zi,:) = interp1(xfA,xA,full(freqB(zi,:)),'linear','extrap');
%             else
%                 if opts.normalizeZeroRate 
%                     zLim = mergeGrid<batchZeroFreq;
%                     mergeGridMat(zLim) = 0;
%                 end                
%                 mExp = mFunc(mergeGridMat);
%                 
%                 if opts.centerFreqDist  
%                     %%
%                     
%                     centerCorrection = batchMeanFreq - median(batchMeanFreq);
%                     cFreqCorr = cFreq - centerCorrection(nIdx)';
%                     
%                     cFreqCorr = max(min(cFreqCorr,1),0);
%                     cFreqCorr(cFreq == 0) = 0;
%                     
%                     expO(zi,:) = interp1(mergeGrid,mExp,cFreqCorr,'linear','extrap');
%                     
%                     
%                     % cExp = interp1(mergeGrid,mExp,cFreq,'linear','extrap');
%                     
%                     % cExpMean = cellfun(@(x)mean(cExp(x)),bVect);
%                     % cExpMeanC =  cExpMean - mean(cExpMean) ;
%                     
%                     % expO(zi,:) = cExp - cExpMeanC(nIdx)';
%                 else
%                     expO(zi,:) = interp1(mergeGrid,mExp,cFreq,'linear','extrap');
%                 end
%                 
%                 % expO(zi,:) = interp1(mergeGrid,mExp,fA,'linear','extrap');                         
%                 %expO(zi,:) = interp1(mergeGrid,mExp,cFreq,'linear','extrap');
%             end           
%         end
%         
% %         freqO(zi,:) = cFreq;
% %         % progressbar(zi/zD);
%         
%         if mod(zi,500) == 0
%             fprintf('%d) %2f\n',zi,zi/nA*100);
%             toc(zt);
%         end
%         
%     end
    
    freqO = sparse(freqO);
    expO = sparse(expO);
end

