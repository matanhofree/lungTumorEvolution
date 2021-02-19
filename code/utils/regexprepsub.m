function outF = regexprepsub(inF,selectV,exp,rep)
        
    if ischar(selectV) 
        selI = strgrep(inF,selectV);
    elseif islogical(selectV) && length(inF) == length(selectV)
        selI = selectV;
    end
    
    outF = inF;
    
    outF(selI) = regexprep(outF(selI),exp,rep);
end