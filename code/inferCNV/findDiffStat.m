function [outDiffStat,AUC] = findDiffStat(xIn,ctrlRef);

%%
    caseIdx = ~ctrlRef(:);   
    xIn = xIn(:);
    
    [X,Y,T,AUC,OPTROCPT] = perfcurve(caseIdx,xIn,'true');    
    zPos = find(X == OPTROCPT(1) & Y == OPTROCPT(2));    
    outDiffStat = T(zPos);
    
    fprintf('AUC = %f\n',AUC);
  %%  
    % caseDrop = (caseIdx & xIn > outDiffStat) | (~caseIdx);
    
    
end