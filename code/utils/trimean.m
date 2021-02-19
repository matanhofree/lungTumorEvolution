function y = trimean(x,dim,p)

    if ~exist('p','var') || isempty(p)
        p = [ 0.25 0.5 0.75 ];
    end
    
    if ~exist('dim','var') || isempty(dim)
        dim = 1;
    end
    
    cQ = quantile(x,p,dim);
    cProp = [ 0.25 0.5 0.25 ]';
    if dim == 1           
        
        y = sum(cQ.*cProp);
    else
        cProp = cProp';
        y = sum(cQ.*cProp,dim);
        
    end
    

end