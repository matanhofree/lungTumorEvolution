function ydata = fastMC_tSNE_new(dataMat,ndim,preplexity,theta,openMPcores,iterations,random_state,verbose,early_exaggeration,learning_rate,distance,n_iter_exag)
%
% dataMat -- Sample by feature matrix (N x D). Transposed for compatibility
% with C code. 
%
%
%     
%     random_state = -1;
%     verbose = 1;
%     early_exaggeration = 12;
%     learning_rate = 200;
%     distance = 1;
    

    ydata = mexTSNE_new(dataMat,ndim,preplexity,theta,openMPcores,iterations,...
        random_state,verbose,early_exaggeration,learning_rate,distance,n_iter_exag);
    
    % row-major to colMajor 
    ydata = reshape(ydata(:),fliplr(size(ydata)))';

end