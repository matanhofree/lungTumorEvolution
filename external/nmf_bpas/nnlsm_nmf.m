function [X,grad,iter,pinvCounter] = nnlsm_nmf(A,B,init,solver)
    
    if nargin < 3
        init = [];        
        solver = 'as';
    elseif nargin < 4
        solver = 'as';
    end
    pinvCounter = 0;

    if isempty(init)
        switch solver
            case 'bp'
                [X,grad,iter,pinvCounter] = nnlsm_blockpivot(A,B,0);
            case 'as'
                [X,grad,iter] = nnlsm_activeset(A,B,1,0);
        end
    else
        switch solver
            case 'bp'
                [X,grad,iter,pinvCounter] = nnlsm_blockpivot(A,B,0,init);
            case 'as'
                [X,grad,iter] = nnlsm_activeset(A,B,0,0,init);
        end
    end
end 