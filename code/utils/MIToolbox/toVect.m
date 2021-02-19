function zm = toVect(zmat)
    %n = length(zmat);
    %zm = zmat(triu(true(n),1));
    %zm = zm(:)';

    if (issparse(zmat))
        zm = fastToVect(2,zmat);
    elseif isa(zmat,'double')
        zm = fastToVect(1,zmat);
    else
        n = length(zmat);
        zm = zmat(triu(true(n),1));
        zm = zm(:);
    end

end