function outS = data2simpleStruct(data,colnames)
    
    if ismatrix(data)
        if length(colnames) == size(data,1)
            for i = 1:length(colnames)
                outS.(colnames{i}) = data(i,:)';
            end
        elseif length(colnames) == size(data,2)
            for i = 1:length(colnames)
                 outS.(colnames{i}) = data(:,i);
            end
        else
            error('Colnames and data dim missmatch');
        end
    else iscell(data)
        error('Not implemented');
    end
        
end