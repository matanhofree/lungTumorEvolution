function ifprintf(ifprint, varargin)
    if ifprint
        fprintf(varargin{:});
    end
end