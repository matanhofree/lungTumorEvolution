function optionFinal = removeString(option)
% Merge two struct options into one struct
% Usage:
% optionFinal = mergeOption(option,optionDefault)
    optionFinal=option;
    if isempty(option)
        return;
    end

    names=fieldnames(option);
    for i=1:numel(names)
        cfield = option.(names{i});
        
        if (isstring(cfield))            
            if numel(cfield) == 1
                optionFinal.(names{i}) = char(cfield);        
            else
                optionFinal.(names{i}) = arrayfun(@(x)char(x),cfield,'uniformoutput',0);
                
            end
        end
    end

end