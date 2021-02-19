function [ydata,Xscore] = run_embed_mcTSNE(Xorig, inOpts )
%TSNE Performs symmetric t-SNE on dataset X
%
    defaultOpts.doNorm = 0;
    
    defaultOpts.pcaMethodType = 'skip';
    defaultOpts.initial_dims = 30;
    defaultOpts.doNorm = 0;
    defaultOpts.covDist = 1;
    
    defaultOpts.no_dims = 2;
    defaultOpts.perplexity = 30;
    defaultOpts.theta = 0.5;
    defaultOpts.multiCore = -1;
    defaultOpts.max_iter = 1000;
    defaultOpts.random_state = 42;
    defaultOpts.verbose = 1;
    defaultOpts.early_aggeration = 12;
    defaultOpts.learning_rate = 200;
    defaultOpts.distance = 1; % 1 - Sq euclid
                             % 0 - Euclid
    
    defaultOpts.n_iter_exag = 250;
    
    defaultOpts.tsneMethod = 'mcTSNE_new';

    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;

    end
    clear defaultOpts;
    disp(opts);
    
    if opts.multiCore == -1
        opts.multiCore = feature('numcores')
        fprintf('Attempting to use %d cores\n',opts.multiCore);
    end

    [D,N] = size(Xorig);

    fprintf('Using MC-tSNE on data (%d,%d).\n',D,N);
    
    Xscore = normalize_expMatrix(Xorig,[],opts);
    if ~strcmp(opts.pcaMethodType,'skip')
        Xscore = extract_pca_embedding(Xscore,W,opts);
    end

    if opts.covDist 
        disp('Dividing rows by norm value (using cosine distances)\n');

        Xnorm = sqrt(nansum(Xscore.^2,2));
        Xscore = bsxfun(@rdivide,Xscore,Xnorm);
    end
    
    zT = tic;
    switch opts.tsneMethod
        case { 'mcTSNE', 'fastMC_tSNE' }
            fprintf('Using original mcTSNE procedure\n');
            ydata = fastMC_tSNE(Xscore',opts.no_dims,opts.perplexity,opts.theta,opts.multiCore,opts.max_iter);
        case { 'mcTSNE_new', 'fastMC_tSNE_new' }
            fprintf('Using updated mcTSNE procedure\n');
            ydata = fastMC_tSNE_new(Xscore,opts.no_dims,opts.perplexity,opts.theta,opts.multiCore,opts.max_iter,opts.random_state,opts.verbose,opts.early_aggeration,opts.learning_rate,opts.distance,opts.n_iter_exag);
        case { 'mcTSNE_sparse', 'fastMC_tSNE_sparse' }
            fprintf('Using sparse distabce based mcTSNE procedure\n');
            warnning('Needs more testing');
            ydata = fastMC_tSNE_sparse(X_orig,opts.no_dims,opts.perplexity,opts.theta,opts.multiCore,opts.max_iter,opts.random_state,opts.verbose,opts.early_aggeration,opts.learning_rate,opts.distance);
        otherwise
            error('tSNE method not implemented');
    end


end


function Xscore = extract_pca_embedding(X,W,opts)

    % Taken from /Users/mhofree/projects/cancer_SC/code_sync/dimReduction/tSNE/tsne_weighted_subsampling.m
    fprintf('Doing PCA ?!');

    switch opts.methodType
        case 'bwpca'
            %%
            disp('Preprocessing data using Bailey weighted PCA!');

            [zV,zU,zScoreWeights,lambda,totVar] = bwpca(X',W',opts.initial_dims,opts.nstarts,opts.doSmooth,opts.epsV,opts.niter,0,0);


            tV = reshape(zV,size(X,1),opts.initial_dims);
            tU = reshape(zU,size(X,2),opts.initial_dims)*diag(sqrt(lambda).^-1);

            Xscore = X'*tV;


        case 'fb_pca'
            %%
            disp('Preprocessing data using fb PCA.');

            [tU,tS,tV] = fb_pca(X',opts.initial_dims,1);

            Xscore = X'*tV;

        case 'SS_ALS'
            disp('Preprocessing data using subsampling ALS-PCA.');
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

            tic;
            [ A, S, Mu, V, CV, HP, LC ] = pca_diag( Xnorm, opts.initial_dims,zopts );
            toc;
            %%
            XnormT = Xnorm;
            XnormT(isnan(XnormT)) = 0;
            %%
            % zXscore = XnormT*S';
            Xscore = XnormT'*A;
            %%
            % figure; scatter(zXscore(:,1),Xscore(:,1))
        case 'ICA'
            %%
            disp('Preprocessing data using fast ICA');

            [Xscore, A, W] = fastica(X,'numOfIC',opts.initial_dims,'stabilization','on');

            Xscore = Xscore';

        case 'robust_cov'
            disp('Preprocessing data using robust cov method');

            opts.dropMinFrac = 0.1;
            opts.dropMinFracSample = 0.05;

            dataM = full(X)';

            zFrac = sum(dataM>0)/size(dataM,1);
            zDrop = zFrac < opts.dropMinFrac;
            fprintf('Dropping %d genes (frac<%f)\n',sum(zDrop),opts.dropMinFrac);
            dataM(:,zDrop) = [];

            fprintf('Using robus corr pca\n');
            zFrac = sum(dataM>0,2)/size(dataM,2);
            zSel = zFrac > opts.dropMinFracSample;
            try
                [sig,mu,~,covOutliers , rCovS] = robustcov(dataM(zSel,:),'Method','olivehawkins');
            catch e
                fprintf('Trying OGK method\n');
                [sig,mu,~,covOutliers , rCovS] = robustcov(dataM(zSel,:),'Method','ogk');
            end

            [coeff,latent] = pcacov(sig);
            dataM = bsxfun(@minus,dataM,mu);
            %%
            Xscore = dataM*coeff(:,1:opts.initial_dims);

        case 'skip'
            fprintf('No PCA preprocessing.');
            Xscore = X';

        otherwise
            %%
            disp('Preprocessing data using simple  PCA (wig...');

            [zCoeff,zScore,zLatent] = pca(X','NumComponents',opts.initial_dims,'Algorithm',opts.pcaAlg,'centered',false);
            %%
            zSigma = sqrt(zLatent*size(X,2));
            Xscore = X*(zScore*diag(zSigma(1:opts.initial_dims).^-1))';
            %%

            % Xscore = X*(zCoeff*diag(zSigma(1:opts.initial_dims).^-1))';
            % Xscore = X'*tV;


    end

end
