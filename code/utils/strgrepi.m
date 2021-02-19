function [tselect,gout] = strgrep(strIn,regE)
    
    if ~iscell(strIn)
        strIn = {strIn};
    end
    tt = regexpi(strIn,regE);
        
    tselect = ~isemptycell(tt);
    gout = strIn(tselect);
    
end
    