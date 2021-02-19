function listVar = extractOverDispersed(outCV,minFreq,minDup,dupTest)

    if nargin < 4 
        dupTest = [];
    end

    zDiffStat = outCV.medianDispNorm;
    zDiffStat(outCV.geneNnzFreq<minFreq) = 0;
    zDiffStat(zDiffStat<0) = 0;

    zVal_medianDispNorm = findElbow(zDiffStat,0);
    zVarMat = zDiffStat > zVal_medianDispNorm;
    
    nD = size(outCV.countCV,1);
    
    if isempty(dupTest)
        listVar = sum(zVarMat,2) >= min(size(zVarMat,2),minDup);
    else
        listVar = false(nD,1);        
        for i = 1:length(dupTest)        
            cTest = dupTest{i};

            selSub = strgrep(outCV.typeID,dupTest{i});
            
            if sum(selSub) == 0
                fprintf('Skipping %s as it is not found\n',dupTest{i}');
                continue;
            end

            zAdd = sum(zVarMat(:,selSub),2) >= min(sum(selSub),minDup);
            listVar = listVar | zAdd;
            
            fprintf('%s) Found %d over-dispersed variables (adding %d)\n',cTest,sum(listVar),sum(zAdd));
        end
        
        if sum(listVar) == 0
            fprintf('No over-dispersed variables. Trying a global test\n');
            listVar = sum(zVarMat,2) >= min(size(zVarMat,2),minDup);
        end                        
    end 
    fprintf('Found a total of %d over-dispersed genes for next step\n',sum(listVar));
    
end
