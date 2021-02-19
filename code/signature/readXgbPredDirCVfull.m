function [outPredS,outScore] = readXgbPredDirCVfull(inPath,assignThr,ambigThr,doScale,inOpts)
    
    if ~exist('assignThr','var') || isempty(assignThr)
        assignThr = 0.9;
    end

    if ~exist('ambigThr','var') 
        ambigThr = 0.5;
    end
    
    if ~exist('doScale','var') 
        doScale = 0;
    end
               
    defaultOpts.refList = [];
    defaultOpts.fullID = 'full';
    defaultOpts.useFullPred = 0;
            
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
    
    zfile = fileList([ inPath '/*_trainLabel.tsv']);
    
    zMatHead = fastTxtRead(zfile{1},'\t');
    zMatRaw = fastTxtRead(zfile{1},'\t',[],[],1);
    trainLabels =data2simpleStruct(zMatRaw,[ 'sampleID' zMatHead(1,:)]);
    %%
    trainLabelsV = trainLabels.(zMatHead{1});
    batchID = trainLabels.(zMatHead{2});
    
    %%
    inPathSummary = sprintf( [ inPath '/all_predSummary.mat']);
    
    if exist(inPathSummary,'file')
        fprintf('Found existing intermediate: %s\n',inPathSummary);
        
        inScore = load(inPathSummary);
    else
        fprintf('Loading raw data...\n');
        [~,inScore] = readXgbPredDir(inPath,assignThr,ambigThr);
        
        save(inPathSummary,'-v7.3','-struct','inScore');
    end
    
    outScore = rmfield(inScore,intersect(fieldnames(inScore),{'predMat','predMatOrig','predAmbig','pred','predAmbigNA'}));        
    
    N = length(outScore.sampleID)';
    outScore.batchID = batchID;
    assert(length(batchID) == N);
    
    if ~isempty(opts.fullID)        
        selFull = strgrep(outScore.fName,opts.fullID);        
        
        if opts.useFullPred         
            fprintf('Using full (non-cv) predictions\n');
            outScore.predMat = inScore.predMat(:,selFull);                        
        else
            fprintf('Using only CV predictions\n');
            outScore.predMat = nan(N,sum(selFull));
        end
        outScore.typeH = inScore.typeH(selFull);
        if ~isempty(setdiff(setdiff(unique(trainLabelsV),'NA'),outScore.typeH))
            warning('Missmatch in label list and typeList');
        end
    else
        fprintf('Full CV not found - using average of CV\n');
        [zB,~,~,~,zCn ] = fastUnique(inScore.typeH);
        outScore.typeH = zB;
        outScore.predMat = cell2mat(cellfun(@(x)median(inScore.predMat(:,x),2),zCn,'unif',0)');
    end        
    
    refList = unique(batchID);   

    for fi = 1:length(refList)
        cf = refList{fi};
        cb = strcmp(batchID,cf);
        selCol = strgrep(inScore.fName,[ '_drop' cf '_' ]);
        
        fprintf('Replacing values for %s - %d (numCol %d)\n',cf,sum(cb),sum(selCol))
        
        inSub = structSubSelectMat(inScore,cb);
        inSub = structSubSelectMat(inSub,selCol);
        
        [~,zia,zib] = intersect(inSub.typeH,outScore.typeH);
        [~,zpa,zpb] = intersect(inSub.sampleID,outScore.sampleID);
        
        outScore.predMat(zpb,zib) = inSub.predMat(zpa,zia);
    end    
    
    [outPredS,outScore] = readXgbPredDir(outScore,assignThr,ambigThr,doScale);    
    cTypes = outScore.typeH;
    
    cSel = ~strcmp(trainLabelsV,'NA');
    trainLabelsV(~cSel) = [];
  
    tabFilter(trainLabelsV)
  
    for zi = 1:length(cTypes)
        cPred = outScore.predMat(cSel,zi);
        cLabel = ~strcmp(trainLabelsV,cTypes{zi});
        [cThr(zi),cAuc(zi)] = findDiffStat(cPred,cLabel);
    end
    outScore.cThr = cThr;
    outScore.cAuc = cAuc;

    cvPred = outScore.predMat > cThr;
    sMult = sum(cvPred,2);
    tabFilter(sMult)

    outScore.predCV = outScore.pred;
    for  zi = 1:length(cTypes)
        zSel = cvPred(:,zi) == 1 & sMult == 1;
        outScore.predCV(zSel) = outScore.typeH(zi);
    end    
    outScore.predCV_NA = outScore.predCV;
    outScore.predCV_NA(sMult >= 2) = { 'NA' };  
    
    %%
    if isfield(trainLabels,'sampleID')
        outScore = addSampleSlice(outScore,trainLabels);
    end
end