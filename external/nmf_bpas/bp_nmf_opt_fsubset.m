% Nonnegative Matrix Factorization by Alternating Nonnegativity Constrained Least Squares
%      using Block Principal Pivoting/Active Set method
%
% This software solves one the following problems: given A and k, find W and H such that
%     (1) minimize 1/2 * || A-WH ||_F^2 
%     (2) minimize 1/2 * ( || A-WH ||_F^2 + alpha * || W ||_F^2 + beta * || H ||_F^2 )
%     (3) minimize 1/2 * ( || A-WH ||_F^2 + alpha * || W ||_F^2 + beta * (sum_(i=1)^n || H(:,i) ||_1^2 ) )
%     where W>=0 and H>=0 elementwise.
%
% Reference:
%  [1] For using this software, please cite:
%          Jingu Kim and Haesun Park, Toward Faster Nonnegative Matrix Factorization: A New Algorithm and Comparisons,
%                 In Proceedings of the 2008 Eighth IEEE International Conference on Data Mining (ICDM'08), 353-362, 2008
%  [2] If you use 'nnls_solver'='as' (see below), please cite:
%          Hyunsoo Kim and Haesun Park, Nonnegative Matrix Factorization Based on Alternating Nonnegativity Constrained Least Squares and Active Set Method,
%                 SIAM Journal on Matrix Analysis and Applications, 2008, 30, 713-730
%
% Written by Jingu Kim (jingu@cc.gatech.edu)
% Copyright 2008-2009 by Jingu Kim and Haesun Park,
%                        School of Computational Science and Engineering,
%                        Georgia Institute of Technology
%
% Check updated code at http://www.cc.gatech.edu/~jingu
% Please send bug reports, comments, or questions to Jingu Kim.
% This code comes with no guarantee or warranty of any kind.
%
% Last modified Feb-20-2010
%
% <Inputs>
%        A : Input data matrix (m x n)
%        k : Target low-rank
%
%        (Below are optional arguments: can be set by providing name-value pairs)
%        TYPE : 'plain' to use formulation (1)
%               'regularized' to use formulation (2)
%               'sparse' to use formulation (3)
%               Default is 'regularized', which is recommended for quick application testing unless 'sparse' or 'plain' is explicitly needed.
%               If sparsity is needed for 'W' factor, then apply this function for the transpose of 'A' with formulation (3).
%                      Then, exchange 'W' and 'H' and obtain the transpose of them.
%               Imposing sparsity for both factors is not recommended and thus not included in this software.
%        NNLS_SOLVER : 'bp' to use the algorithm in [1]
%                      'as' to use the algorithm in [2]
%                      Default is 'bp', which is in general faster.
%        ALPHA : Parameter alpha in the formulation (2) or (3). 
%                Default is the average of all elements in A. No good justfication for this default value, and you might want to try other values.
%        BETA : Parameter beta in the formulation (2) or (3).
%               Default is the average of all elements in A. No good justfication for this default value, and you might want to try other values.
%        MAX_ITER : Maximum number of iterations. Default is 100.
%        MIN_ITER : Minimum number of iterations. Default is 20.
%        MAX_TIME : Maximum amount of time in seconds. Default is 100,000.
%        W_INIT : (m x k) initial value for W.
%        H_INIT : (k x n) initial value for H.
%        TOL : Stopping tolerance. Default is 1e-3. If you want to obtain a more accurate solution, decrease TOL and increase MAX_ITER at the same time.
%        VERBOSE : 0 (default) - No debugging information is collected.
%                  1 (debugging purpose) - History of computation is returned by 'HIS' variable.
%                  2 (debugging purpose) - History of computation is additionally printed on screen.
% <Outputs>
%        W : Obtained basis matrix (m x k)
%        H : Obtained coefficients matrix (k x n)
%        iter : Number of iterations
%        HIS : (debugging purpose) History of computation
% <Usage Examples>
%        nmf(A,10)
%        nmf(A,20,'verbose',2)
%        nmf(A,30,'verbose',2,'nnls_solver','as')
%        nmf(A,5,'verbose',2,'type','sparse')
%        nmf(A,60,'verbose',1,'type','plain','w_init',rand(m,k))
%        nmf(A,70,'verbose',2,'type','sparse','nnls_solver','bp','alpha',1.1,'beta',1.3)

function [W,H,iter,HIS,final]=bp_nmf_opt_fsubset(A,k,W_init,H_init,inOpts)
    
    [m,n] = size(A);
    

    % Default configuration
    defaultOpts.m = m;
    defaultOpts.n = n;
    defaultOpts.type = 'regularized';
    defaultOpts.nnls_solver = 'bp';
    defaultOpts.alpha = [];
    defaultOpts.beta = [];
    defaultOpts.max_iter = 100;
    defaultOpts.min_iter = 20;
    defaultOpts.max_time = 1e6;
    defaultOpts.tol = 1e-4;
    defaultOpts.verbose = -1;
    defaultOpts.ST_RULE = 1;
    defaultOpts.fixDim = [];
               
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    disp(opts);
    clear defaultOpts;
    
    ST_RULE = opts.ST_RULE;
    if ~exist('W_init','var') || isempty(W_init)        
        W = rand(m,k);
    else
        W = W_init;        
    end
    
    if ~exist('H_init','var') || isempty(H_init)                
        H = rand(k,n);
    else
        H = H_init;        
    end

    if size(W,2) < k
        if isempty(opts.fixDim) && opts.fixDim == 1
            opts.fixDim = trueV(1:size(W,2),k);
        end
        
        kdiff = k - size(W,2);
        W = [ W rand(m,kdiff) ];
    end
    
    if size(H,1) < k
        kdiff = k - size(H,1);
        H = [ H; rand(kdiff,n) ];
    end
    fixDim = opts.fixDim;
%     % Read optional optsameters
%     if (rem(length(varargin),2)==1)
%         error('Optional optsameters should always go by pairs');
%     else
%         for i=1:2:(length(varargin)-1)
%             switch upper(varargin{i})
%                 case 'TYPE',                opts.type = varargin{i+1};
%                 case 'NNLS_SOLVER',         opts.nnls_solver = varargin{i+1};
%                 case 'ALPHA',               argAlpha = varargin{i+1};,opts.alpha = argAlpha;
%                 case 'BETA',                argBeta = varargin{i+1};,opts.beta = argBeta;
%                 case 'MAX_ITER',            opts.max_iter = varargin{i+1};
%                 case 'MIN_ITER',            opts.min_iter = varargin{i+1};
%                 case 'MAX_TIME',            opts.max_time = varargin{i+1};
%                 case 'W_INIT',              W = varargin{i+1};
%                 case 'H_INIT',              H = varargin{i+1};
%                 case 'TOL',                 opts.tol = varargin{i+1};
%                 case 'VERBOSE',             opts.verbose = varargin{i+1};
%                 otherwise
%                     error(['Unrecognized option: ',varargin{i}]);
%             end
%         end
%     end

    if isempty(opts.alpha) 
        opts.alpha = mean(A(:));
    end
    if isempty(opts.beta) 
        opts.beta = mean(A(:));
    end

    % for regularized/soptsse case
    if strcmp(opts.type,'regularized')
        
        salphaI = sqrt(opts.alpha)*eye(k);
        sbetaI = sqrt(opts.beta)*eye(k);
        zerokn = zeros(k,n);
        zerokm = zeros(k,m);
    elseif strcmp(opts.type,'soptsse')
        salphaI = sqrt(opts.alpha)*eye(k);
        sbetaE = sqrt(opts.beta)*ones(1,k);
        betaI = opts.beta*ones(k,k);
        zero1n = zeros(1,n);
        zerokm = zeros(k,m);
    elseif ~strcmp(opts.type,'plain')
        error(['Unrecognized type: use ''plain'', ''regularized'', or ''soptsse''.']);
    end
    
    if ~strcmp(opts.nnls_solver,'bp') && ~strcmp(opts.nnls_solver,'as')
        error(['Unrecognized nnls_solver: use ''bp'' or ''as''.']);
    end
    
    if opts.verbose 
        display(opts); 
    end
    
    HIS = [];
    if opts.verbose > 0         % collect information for analysis/debugging
        if any(fixDim)
            [gradW,gradH] = getGradient(A,W(:,~fixDim),H(~fixDim,:),opts.type,opts.alpha,opts.beta);
            initGrNormW = norm(gradW,'fro');
            initGrNormH = norm(gradH,'fro');
            initNorm = norm(A,'fro'); 
            numSC = 3;
            initSCs = zeros(numSC,1);

            for j=1:numSC
                initSCs(j) = getInitCriterion(j,A,W(:,~fixDim),H(~fixDim,:),opts.type,opts.alpha,opts.beta,gradW,gradH);
            end
        else
            [gradW,gradH] = getGradient(A,W,H,opts.type,opts.alpha,opts.beta);
            initGrNormW = norm(gradW,'fro');
            initGrNormH = norm(gradH,'fro');
            initNorm = norm(A,'fro'); 
            numSC = 3;
            initSCs = zeros(numSC,1);

            for j=1:numSC
                initSCs(j) = getInitCriterion(j,A,W,H,opts.type,opts.alpha,opts.beta,gradW,gradH);
            end
        end
        
%---(1)------(2)--------(3)--------(4)--------(5)---------(6)----------(7)------(8)-----(9)-------(10)--------------(11)-------
% iter # | elapsed | totalTime | subIterW | subIterH | rel. obj.(%) | NM_GRAD | GRAD | DELTA | W density (%) | H density (%)
%------------------------------------------------------------------------------------------------------------------------------
        HIS = zeros(1,11);
        HIS(1,[1:5])=0;
        ver.initGrNormW = initGrNormW;
        ver.initGrNormH = initGrNormH;
        ver.initNorm = initNorm;                            HIS(1,6)=ver.initNorm;
        ver.SC1 = initSCs(1);                               HIS(1,7)=ver.SC1;
        ver.SC2 = initSCs(2);                               HIS(1,8)=ver.SC2;
        ver.SC3 = initSCs(3);                               HIS(1,9)=ver.SC3;
        ver.W_density = length(find(W>0))/(m*k);            HIS(1,10)=ver.W_density;
        ver.H_density = length(find(H>0))/(n*k);            HIS(1,11)=ver.H_density;
        if opts.verbose == 2, display(ver);, end
        tPrev = cputime;
    end
    
    tStart = cputime;, tTotal = 0;
    
    if any(fixDim)
        initSC = getInitCriterion(ST_RULE,A,W(:,~fixDim),H(~fixDim,:),opts.type,opts.alpha,opts.beta);
    else
        initSC = getInitCriterion(ST_RULE,A,W,H,opts.type,opts.alpha,opts.beta);
    end
    
    SCconv = 0; SC_COUNT = 3;

    for iter=1:opts.max_iter
        switch opts.type
            case 'plain'
                % [H,gradHX,subIterH] = nnlsm(W,A,H,opts.nnls_solver);
                [H,~,subIterH] = nnlsm(W,A,H,opts.nnls_solver);
                [W,gradW,subIterW] = nnlsm(H',A',W',opts.nnls_solver);
                W=W';
                gradW=gradW';
                gradH = (W'*W)*H - W'*A;
            case 'regularized'
                % [H,gradHX,subIterH] = nnlsm([W;sbetaI],[A;zerokn],H,opts.nnls_solver);
                [H,~,subIterH] = nnlsm([W;sbetaI],[A;zerokn],H,opts.nnls_solver);
                [W,gradW,subIterW] = nnlsm([H';salphaI],[A';zerokm],W',opts.nnls_solver);
                W=W';
                gradW=gradW';
                gradH = (W'*W)*H - W'*A + opts.beta*H;
            case 'soptsse'
                % [H,gradHX,subIterH] = nnlsm([W;sbetaE],[A;zero1n],H,opts.nnls_solver);
                [H,~,subIterH] = nnlsm([W;sbetaE],[A;zero1n],H,opts.nnls_solver);
                [W,gradW,subIterW] = nnlsm([H';salphaI],[A';zerokm],W',opts.nnls_solver);
                W=W';
                gradW=gradW';
                gradH = (W'*W)*H - W'*A + betaI*H;
        end
        
        if any(fixDim) 
            if iter == 1
                fprintf('Fixing %d W dim\n',sum(fixDim));
            end
            W(:,fixDim) = W_init(:,fixDim);     
            if opts.verbose ~= 0
                fprintf('*');
            end
        end
        
        if opts.verbose > 0          % collect information for analysis/debugging
            elapsed = cputime-tPrev;
            tTotal = tTotal + elapsed;
            % ver = 0;
            idx = iter+1;
%---(1)------(2)--------(3)--------(4)--------(5)---------(6)----------(7)------(8)-----(9)-------(10)--------------(11)-------
% iter # | elapsed | totalTime | subIterW | subIterH | rel. obj.(%) | NM_GRAD | GRAD | DELTA | W density (%) | H density (%)
%------------------------------------------------------------------------------------------------------------------------------
            ver.iter = iter;                                    HIS(idx,1)=iter;
            ver.elapsed = elapsed;                              HIS(idx,2)=elapsed;
            ver.tTotal = tTotal;                                HIS(idx,3)=tTotal;
            ver.subIterW = subIterW;                            HIS(idx,4)=subIterW;
            ver.subIterH = subIterH;                            HIS(idx,5)=subIterH;
            ver.relError = norm(A-W*H,'fro')/initNorm;          HIS(idx,6)=ver.relError;
            
            if any(fixDim)
                ver.SC1 = getStopCriterion(1,A,W(:,~fixDim),H(~fixDim,:),opts.type,opts.alpha,opts.beta,gradW(:,~fixDim),gradH(~fixDim,:))/initSCs(1);     HIS(idx,7)=ver.SC1;
                ver.SC2 = getStopCriterion(2,A,W(:,~fixDim),H(~fixDim,:),opts.type,opts.alpha,opts.beta,gradW(:,~fixDim),gradH(~fixDim,:))/initSCs(2);     HIS(idx,8)=ver.SC2;
                ver.SC3 = getStopCriterion(3,A,W(:,~fixDim),H(~fixDim,:),opts.type,opts.alpha,opts.beta,gradW(:,~fixDim),gradH(~fixDim,:))/initSCs(3);     HIS(idx,9)=ver.SC3;
            else
            
                ver.SC1 = getStopCriterion(1,A,W,H,opts.type,opts.alpha,opts.beta,gradW,gradH)/initSCs(1);     HIS(idx,7)=ver.SC1;
                ver.SC2 = getStopCriterion(2,A,W,H,opts.type,opts.alpha,opts.beta,gradW,gradH)/initSCs(2);     HIS(idx,8)=ver.SC2;
                ver.SC3 = getStopCriterion(3,A,W,H,opts.type,opts.alpha,opts.beta,gradW,gradH)/initSCs(3);     HIS(idx,9)=ver.SC3;
            end
            ver.W_density = length(find(W>0))/(m*k);            HIS(idx,10)=ver.W_density;
            ver.H_density = length(find(H>0))/(n*k);            HIS(idx,11)=ver.H_density;
            if opts.verbose == 2, display(ver);, end
            tPrev = cputime;
        elseif opts.verbose == -1
            fprintf('.');
        end
        
        if (iter > opts.min_iter)
            if any(fixDim)
                SC = getStopCriterion(ST_RULE,A,W(:,~fixDim),H(~fixDim,:),opts.type,opts.alpha,opts.beta,gradW(:,~fixDim),gradH(~fixDim,:));                
            else
                SC = getStopCriterion(ST_RULE,A,W,H,opts.type,opts.alpha,opts.beta,gradW,gradH);
            end
            if (opts.verbose>0 && (tTotal > opts.max_time)) || (~(opts.verbose>0) && ((cputime-tStart)>opts.max_time))
                break;
            end
            
            if (SC/initSC <= opts.tol)
                SCconv = SCconv + 1;
                if (SCconv >= SC_COUNT)
                    break;
                end
            else
                SCconv = 0;
            end
        end
    end
    [m,n]=size(A);
    norm2=sqrt(sum(W.^2,1));
    toNormalize = norm2>0;
    W(:,toNormalize) = W(:,toNormalize)./repmat(norm2(toNormalize),m,1);
    H(toNormalize,:) = H(toNormalize,:).*repmat(norm2(toNormalize)',1,n);
    
    final.iterations = iter;
    if any(fixDim)
        fprintf('BP-NMF (%s, %s, %d fixed dims) converged after %d iterations.\n',opts.nnls_solver,opts.type,sum(fixDim),iter);
    else
        fprintf('BP-NMF (%s, %s) converged after %d iterations.\n',opts.nnls_solver,opts.type,iter);
    end
    if opts.verbose > 0
        final.elapsed_total = tTotal;
    else
        final.elapsed_total = cputime-tStart;
    end
    final.relative_error = norm(A-W*H,'fro')/norm(A,'fro');
    final.W_density = length(find(W>0))/(m*k);
    final.H_density = length(find(H>0))/(n*k);
    final.opts = opts;
    
    if ~isempty(HIS)
        HIS = array2table(HIS,'VariableNames',{'iter' 'elapsed' 'totalTime' 'subIterW' 'subIterH' 'rel_obj_percent' 'NM_GRAD' 'GRAD' 'DELTA' 'W density' 'H density' });
        final.history = HIS;
    end
    if opts.verbose > 0
        display(HIS);
    end
end

%------------------------------------------------------------------------------------------------------------------------
%                                    Utility Functions 
%------------------------------------------------------------------------------------------------------------------------
function [X,grad,iter] = nnlsm(A,B,init,solver)
    switch solver
        case 'bp'
            [X,grad,iter] = nnlsm_blockpivot(A,B,0,init);
        case 'as'
            [X,grad,iter] = nnlsm_activeset(A,B,1,0,init);
    end
end    
%-------------------------------------------------------------------------------
function retVal = getInitCriterion(stopRule,A,W,H,type,alpha,beta,gradW,gradH)
% STOPPING_RULE : 1 - Normalized proj. gradient
%                 2 - Proj. gradient
%                 3 - Delta by H. Kim
%                 0 - None (want to stop by MAX_ITER or MAX_TIME)
    if nargin~=9
        [gradW,gradH] = getGradient(A,W,H,type,alpha,beta);
    end
    [m,k]=size(W);, [k,n]=size(H);, numAll=(m*k)+(k*n);
    switch stopRule
        case 1
            retVal = norm([gradW; gradH'],'fro')/numAll;
        case 2
            retVal = norm([gradW; gradH'],'fro');
        case 3
            retVal = getStopCriterion(3,A,W,H,type,alpha,beta,gradW,gradH);
        case 0
            retVal = 1;
    end
end
%-------------------------------------------------------------------------------
function retVal = getStopCriterion(stopRule,A,W,H,type,alpha,beta,gradW,gradH)
% STOPPING_RULE : 1 - Normalized proj. gradient
%                 2 - Proj. gradient
%                 3 - Delta by H. Kim
%                 0 - None (want to stop by MAX_ITER or MAX_TIME)
    if nargin~=9
        [gradW,gradH] = getGradient(A,W,H,type,alpha,beta);
    end

    switch stopRule
        case 1
            % pGradW = gradW(gradW<0|W>0);
            % pGradH = gradH(gradH<0|H>0);
            pGrad = [gradW(gradW<0|W>0); gradH(gradH<0|H>0)];
            pGradNorm = norm(pGrad);
            retVal = pGradNorm/length(pGrad);
        case 2
            % pGradW = gradW(gradW<0|W>0);
            % pGradH = gradH(gradH<0|H>0);
            pGrad = [gradW(gradW<0|W>0); gradH(gradH<0|H>0)];
            retVal = norm(pGrad);
        case 3
            resmat=min(H,gradH); resvec=resmat(:);
            resmat=min(W,gradW); resvec=[resvec; resmat(:)]; 
            deltao=norm(resvec,1); %L1-norm
            num_notconv=length(find(abs(resvec)>0));
            retVal=deltao/num_notconv;
        case 0
            retVal = 1e100;
    end
end
%-------------------------------------------------------------------------------
function [gradW,gradH] = getGradient(A,W,H,type,alpha,beta)
    switch type
        case 'plain'
            gradW = W*(H*H') - A*H';
            gradH = (W'*W)*H - W'*A;
        case 'regularized'
            gradW = W*(H*H') - A*H' + alpha*W;
            gradH = (W'*W)*H - W'*A + beta*H;
        case 'soptsse'
            k=size(W,2);
            betaI = beta*ones(k,k);
            gradW = W*(H*H') - A*H' + alpha*W;
            gradH = (W'*W)*H - W'*A + betaI*H;
    end
end
