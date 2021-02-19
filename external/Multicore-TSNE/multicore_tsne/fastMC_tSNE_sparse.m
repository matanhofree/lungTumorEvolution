function ydata = fastMC_tSNE_sparse(distMat,ndim,preplexity,theta,openMPcores,iterations,random_state,verbose,early_aggeration,learning_rate,distance)
%
% dataMat -- Sample by feature matrix (N x D). Transposed for compatibility
% with C code. 
%
% %
%     
%     random_state = -1;
%     verbose = 1;
%     early_exaggeration = 12;
%     learning_rate = 200;
%     distance = 1;
    fprintf('Please check input is a dist matrix');
    
    %%
    N = length(distMat);
    %%
    K = preplexity*3;
    
    if ~issymmetric(distMat)
        error('This should be a symmetric distance matrix');
    end    
    %%
    [P, beta] = d2p_matlab(distMat, preplexity);
    P = 0.5 * (P + P');                                 % symmetrize P-values
    P = max(P ./ sum(P(:)), realmin);
    
    %% Turn into sim sparse 
    [pOrd,pIdx] = sort(P,'descend');
    
    %%
    i = 0:K:(N*K);
    j = pIdx(1:K,:)-1;    
    v = pOrd(1:K,:);
    
    i = (i);
    j = j(:);
    v = v(:);
    
    %%
    
    ydata = mexTSNE_new_step2(i,j,v,size(P,1),ndim,preplexity,theta,openMPcores,iterations,...
        random_state,verbose,early_exaggeration,learning_rate,distance);
    
    % row-major to colMajor 
    ydata = reshape(ydata(:),fliplr(size(ydata)))';

end