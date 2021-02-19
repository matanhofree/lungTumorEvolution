function [listVar,listVarMat] = extractOverDispersedAndFilter(outCV,minFreq,minDup,dupTest,countData,inOpts)

    defaultOpts.filterExpQ = 0.995;
    defaultOpts.batchID_var = 'batchID';
    defaultOpts.expVar = 'normTPM';
    
    defaultOpts.removeMT = 1;
    defaultOpts.removeHLA = 1;
    
    defaultOpts.removeExtList = [];
    
    defaultOpts.globalDup = [];
    
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    disp(opts);

    if nargin < 4 
        dupTest = [];
    end

    cDiffStat = outCV.medianDispNorm;
    cDiffStat(outCV.geneNnzFreq<minFreq) = 0;
    cDiffStat(cDiffStat<0) = 0;
    
    if ~isempty(opts.filterExpQ)    
        fprintf('Removing top exp quantile (%f)\n',opts.filterExpQ);
        [cMeanExp,~,cBatch] = summarize_subset_value(countData.(opts.expVar),countData.(opts.batchID_var));        
        
        zSel = ~ismember(cBatch,outCV.typeID);
        
        cMeanExp(:,zSel) = [];
        cBatch(zSel) = [];
        %%
        [~,zia,zib] = intersect(outCV.typeID,cBatch);
        
        zord = zib(argsort(zia));
        %
        cMeanExp = cMeanExp(:,zord);
        cBatch = cBatch(zord);
        %%
        if ~isequal(cBatch,outCV.typeID)
            error('Expression batch do not match the overdisperssion batch. Fix!');            
        end
        %%
        cThr = quantile(cMeanExp,opts.filterExpQ);
                
        cDiffStat(cMeanExp > cThr) = 0;       
    end   
    
    if opts.removeMT
        cSel = strgrepi(countData.geneID,'MT-');
        fprintf('Filtering the following mito genes:\n');        
        disp(countData.geneID(cSel));
        cDiffStat(cSel,:) = 0;       
    end
    
    if opts.removeHLA 
        cSel = strgrepi(countData.geneID,'hla-');
        fprintf('Filtering the following HLA genes:\n');        
        disp(countData.geneID(cSel));        
        cDiffStat(cSel,:) = 0;       
    end
    
    zVal_medianDispNorm = findElbow(cDiffStat,0);
    zVarMat = cDiffStat > zVal_medianDispNorm;    
    nD = size(outCV.countCV,1);
    
    if isempty(dupTest)
        listVar = sum(zVarMat,2) >= min(size(zVarMat,2),minDup);
        listVarMat = zVarMat;
    else
        
        listVarMat = false(nD,length(dupTest));
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
            listVarMat(:,i) = zAdd;
            fprintf('%s) Found %d over-dispersed variables (adding %d)\n',cTest,sum(listVar),sum(zAdd));
        end
        
        if sum(listVar) == 0
            fprintf('No over-dispersed variables. Trying a global test\n');
            listVar = sum(zVarMat,2) >= min(size(zVarMat,2),minDup);
        end                                        
    end 
    fprintf('Found a total of %d over-dispersed genes for next step\n',sum(listVar));
    
    if ~isempty(opts.globalDup) && ~isempty(dupTest)
        listVar = sum(listVarMat,2) > opts.globalDup;
    end
end
