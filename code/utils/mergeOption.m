function optionFinal = mergeOption(option,optionDefault)
% Merge two struct options into one struct
% Usage:
% optionFinal = mergeOption(option,optionDefault)
    optionFinal=optionDefault;
    if isempty(option)
        return;
    end

    names=fieldnames(option);
    for i=1:numel(names)
        cfield = option.(names{i});
        if isstring(cfield)
            cfield = cellstr(cfield);
        end
        if (isstruct(cfield) && isfield(optionDefault,names{i}))
            cfield = mergeOption(cfield,optionDefault.(names{i}));
            if isstring(cfield)
                cfield = cellstr(cfield);
            end
            
            optionFinal.(names{i}) = cfield;
        else
            optionFinal.(names{i}) = cfield;
        end
    end

end