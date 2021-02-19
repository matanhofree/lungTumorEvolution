function [outStat,XscoreOut,baseVectOut] = consensusGraphDistNMFConcat(Xorig,inOpts)


    defaultOpts.doNorm = 0;
    defaultOpts.epsV = 1e-6;
    
    defaultOpts.pcaMethod = 'NeNMF'; % Alternative: fb_pca, wpca, standard
    defaultOpts.seed = 42;
    defaultOpts.max_epoch = 200;
    defaultOpts.verbose = 1;
 
    defaultOpts.niter = 25;
    defaultOpts.normFun = @(X)nnzcenter(X,nnzmean(X,1));

    defaultOpts.nstarts = 1;
    defaultOpts.fb_iter = 5;
    defaultOpts.verbose = 1;
    
    defaultOpts.tsneMethod = 'resample'; % resampling - Resampling based approach to pca distance
                                           % fast       - bh tsne
                                           % cosine     - cosine similarity
                                           % dist       - Use externa distance function 
                                                                                          
    defaultOpts.resampleN = 100;
    defaultOpts.resampleFrac = [ 0.8 0.8 ];
    defaultOpts.dimR = [ 15 40 ];

    defaultOpts.distType = 'cosine';
    defaultOpts.distMPower = 0;
    defaultOpts.graphDist = 1;

    defaultOpts.pb = 0;
    defaultOpts.verbose = 0;
    defaultOpts.normCost = 1;
    
    defaultOpts.zipMean.robStdFact = 2;
    defaultOpts.zipMean.doZIP = 0;
    
    defaultOpts.doDistance = 0;
    defaultOpts.inferFull = 1;
    defaultOpts.nmf_max_iter = 5000;
    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;

%     if (opts.minCnt<1)
%         opts.minCnt = ceil(opts.resampleN*opts.minCnt);
%     end
    disp(opts);    
    
    if ~exist('W','var') || isempty(W)
        W = [];
    end
    
    outStat = [];
    XscoreOut = [];
    baseVectOut = [];
    %%
    [D,N] = size(Xorig);
    
    % Normalize     
    X = normalize_expMatrix(Xorig,W,opts);
        
    zT = tic;
    if opts.pb
        progressbar();
    end
    
    if opts.doDistance 
        outKernalSum = nan(N,N);
        outKernalSumSq = nan(N,N);
        outKernalDenom = zeros(N,N);
    end
    
    for i = 1:opts.resampleN
        if opts.resampleN > 1  && ~isempty(opts.resampleFrac)
            fprintf('Subsample %d (time %.2f)\n',i,toc(zT));zT = tic;
            if opts.resampleFrac(1) < 1
                smpCols = randsample(N,floor(opts.resampleFrac(1)*N)); 
           
            else
                smpCols = 1:N; 
            end
            if opts.resampleFrac(2) < 1
                smpRows = randsample(D,floor(opts.resampleFrac(2)*D));
            else
                smpRows = 1:D;
            end

            if numel(opts.dimR) == 1
                nDim = opts.dimR;
            elseif numel(opts.dimR) == 2
                nDim = randsample(opts.dimR(1):opts.dimR(2),1);            
            else
                error('Number of elements must be fixed or a nonegative range');
            end

            if isempty(W)
                [Xscore,baseVect] = extract_embedding(X(smpRows,smpCols),[],nDim,opts);                        
            else
                error('Not implemented');
                % Xscore = extract_pca_embedding(X(smpRows,smpCols),W(smpRows,smpCols),nDim,opts);                        
            end
        else
            fprintf(1,'Note: no sampling');
            
            nDim = round(mean(opts.dimR)); 
            if isempty(W)
                [Xscore,baseVect] = extract_embedding(X,[],nDim,opts);                        
            else
                [Xscore,baseVect] = extract_embedding(X,W,nDim,opts);                        
            end
            
            smpCols = true(N,1);
        end

        if opts.doDistance 
            outKernal = pdist(Xscore,opts.distType); 
            if opts.distMPower > 1
                if strcmp(opts.distType,'cosine')
                    outKernal = 1 - (1-outKernal).^opts.distMPower;
                else
                    outKernal = outKernal.^opts.distMPower;
                end            
            end
            outKernal = squareform(outKernal);

            outKernalSum(smpCols,smpCols) = nansum(cat(3,outKernalSum(smpCols,smpCols),outKernal),3);
            outKernalSumSq(smpCols,smpCols) = nansum(cat(3,outKernalSumSq(smpCols,smpCols),outKernal.^2),3);

            nSmp = length(smpCols);
            outKernalDenom(smpCols,smpCols) =  nansum(cat(3,outKernalDenom(smpCols,smpCols),ones(nSmp)),3);        

            if opts.pb == 1
                progressbar(i/opts.resampleN);
            %elseif opts.pb == 2
            %    parfor_progress(-1,parForMonFile);
            end
        end
        
        if opts.inferFull    
            
            if strgrepi(opts.pcaMethod,'PCA')
                
                fprintf('Completeing matrix by matrix complition\n');

                % zH = X(smpRows,:);
                
                zH = (baseVect'*baseVect)\(baseVect'*X(smpRows,:));
                %%
                zW = (zH*zH')\zH*X';
                %% 
                XscoreOut{i} = zH';
                baseVectOut{i} = zW;

            else
                fprintf('Completeing matrix by NNLSM\n');
                % H = nnlsm_activeset(W'*W, W'*V, 0, 1, H);
                % W = nnlsm_activeset(H*H', H*V', 0, 1, W');
                % W = W';                        
                % xScoreInf = fcnnls(X(smpRows,:), baseVect);
                % zH = nnlsm_activeset(baseVect'*baseVect, baseVect'*X(smpRows,:), 0, 0);
                zH = nnlsm_activeset(baseVect'*baseVect, baseVect'*X(smpRows,:), 0, 1);
                % zH = nnlsm_blockpivot(baseVect'*baseVect, baseVect'*X(smpRows,:), 0, 1);

                %%
                zW = nnlsm_activeset(zH*zH', zH*X', 0, 1);
                % zW = nnlsm_blockpivot(zH*zH', zH*X', 0, 1);

                XscoreOut{i} = zH';
                baseVectOut{i} = zW;
            end
        else
            baseVectOut{i} = baseVect;
            XscoreOut{i} = Xscore;            
        end
        
    end
    
    if opts.doDistance
        outStat.kernal = outKernalSum./outKernalDenom;
        outStat.var = (outKernalSumSq - (outKernalSum.^2)./outKernalDenom)./(max(outKernalDenom-1,1));
    end
    
end


function [Xscore,baseVect] = extract_embedding(X,W,nDim,opts)
        switch opts.pcaMethod
            case 'bwpca'                
                %%
                fprintf('Preprocessing data using Bailey weighted PCA (d=%d).\n',nDim);

                [zV,zU,zScoreWeights,lambda,totVar] = bwpca(X',W',opts.initial_dims,opts.nstarts,opts.doSmooth,opts.epsV,opts.niter,0,0);

                tV = reshape(zV,size(X,1),nDim);
                tU = reshape(zU,size(X,2),nDim)*diag(sqrt(lambda).^-1);               
                Xscore = X'*tV;               
                
            case 'fb_pca'  
                %%
                fprintf('Preprocessing data using fb PCA (d=%d).\n',nDim);
                
                [tU,tS,tV] = fb_pca(X',nDim,1);
                
                % Xscore = X*tU;
                Xscore = X'*tV;
                % clear tU tS tV;
                baseVect = tV;
                
            case 'SS_ALS'
                fprintf('Preprocessing data using subsampling ALS-PCA (d=%d).\n',nDim);
                %%
                W = bsxfun(@rdivide,W,sum(W,1));
                %%
                [D,N] = size(W);
                Xnan = X;
                
                tt = tic;
                zIsZero = X == 0;
                for i = 1:N
                   cZero = find(zIsZero(:,i));
                   nSmpZeros = floor(length(cZero)*opts.dropZeros);
                   zIdx = randsample(cZero,nSmpZeros,1,W(cZero,i));
                   Xnan(zIdx,i) = nan;
                end
                fprintf(1,'SS took %d\n',toc(tt));
                
                zopts = struct( 'maxiters', 250,...
                               'algorithm', 'vb',...
                               'earlystop', 1,...
                               'rmsstop',[ 25 1e-4 1e-3 ],...
                               'verbose', 1,...
                               'display', 1);
           
                %
                Xnorm = bsxfun(@minus,Xnan,nanmean(Xnan,2));
                %
                
                tic;
                [ A, S, Mu, V, CV, HP, LC ] = pca_diag( Xnorm, nDim,zopts );
                toc;
                %%
                XnormT = Xnorm;
                XnormT(isnan(XnormT)) = 0;
                %%
                % zXscore = XnormT*S';
                Xscore = XnormT'*A;
                %% 
                % figure; scatter(zXscore(:,1),Xscore(:,1))
            case 'NMF_bpp'
                
                
                zoptions.alg = 'anls_bpp';
                zoptions.max_epoch = opts.max_epoch;
                zoptions.verbose = opts.verbose;
                [nmfOut, nmfIno] = nmf_anls(X, nDim,zoptions);
                
                Xscore = nmfOut.H';
            case 'NeNMF'
                               
                [baseVect,H]=NeNMF(X,nDim,'MAX_ITER',opts.nmf_max_iter);
                
                Xscore = H';

            otherwise
                %%
                fprintf('Preprocessing data using simple  PCA (d=%d).\n',nDim);
                
                [zCoeff,zScore,zLatent] = pca(X','NumComponents',nDim,'Algorithm',opts.pcaAlg,'centered',false);
                %%
                zSigma = sqrt(zLatent*size(X,2));                
                Xscore = X*(zScore*diag(zSigma(1:nDim).^-1))';
                %%
                
                % Xscore = X*(zCoeff*diag(zSigma(1:nDim).^-1))';
                % Xscore = X'*tV;
                
                
        end

end 
