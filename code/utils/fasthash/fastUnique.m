function [b, m, n, cnt,dup_list,dupVect] = fastUnique(A)
% A unique function that doesn't break the ordering things are given in
% Outputs:
% b unique list
% b = A(m)
% A = b(n)
% cnt = count(b(i)))
% dup_list = list of occurances of each item in original array!
    if iscategorical(A)
        A = cellstr(A);
    end

    if iscell(A) && iscellstr(A)
        [b, m, n, cnt,dup_list] = fastUniqueStr(A);
    elseif isnumeric(A) | islogical(A)
        switch class(A)
            case 'double'
                if all(mod(A,1) == 0) && min(A) >= 0
                    A = uint64(A);
                    
                    % [b, m, n, cnt,dupVect] = fastUniqueIntMultiPass(A);
                    [b, m, n, cnt,dupVect] = fastUniqueIntHash(A,uint64(length(unique(A))));
                    
                    if nargout > 4
                        dupListIdx = [ 1; cumsum(cnt)+1 ];
                        dup_list = arrayfun(@(zi)dupVect(dupListIdx(zi):dupListIdx(zi+1)-1)',1:length(dupListIdx)-1,'uniformoutput',0);
                    end

                else
                    error('Array contains nans or non-integer values');
                end                
            case 'uint64'
                % [b, m, n, cnt,dup_list] = fastUniqueInt(A);
                %% 
%                 [b, m, n, cnt,dupListIdx,dupVect] = fastUniqueIntFlat(A)
%                 dup_list = cell(size(b));
%                 dup_list(1:end-1) = arrayfun(@(zi)dupVect(dupListIdx(zi):dupListIdx(zi+1)-1)',1:length(dupListIdx)-1,'uniformoutput',0);
%                 dup_list{end+1} = dupVect(dupListIdx(end):end)';                    
                [b, m, n, cnt,dupVect] = fastUniqueIntHash(A,uint64(length(unique(A))));

                if nargout > 4
                    dupListIdx = [ 1; cumsum(cnt)+1 ];
                    dup_list = arrayfun(@(zi)dupVect(dupListIdx(zi):dupListIdx(zi+1)-1)',1:length(dupListIdx)-1,'uniformoutput',0);
                end

            case { 'uint32', 'logical' }
                [b, m, n, cnt,dupVect] = fastUniqueIntHash(uint64(A),uint64(length(unique(A))));

                if nargout > 4
                    dupListIdx = [ 1; cumsum(cnt)+1 ];
                    dup_list = arrayfun(@(zi)dupVect(dupListIdx(zi):dupListIdx(zi+1)-1)',1:length(dupListIdx)-1,'uniformoutput',0);
                end

            otherwise
                error('Unimplemented');
        end
    else
        error('Only cell or numeric inputs are supported');
    end
    
%     b = double(b);
%     m = double(m);
%     n = double(n);
end
       
        