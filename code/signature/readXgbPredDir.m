function [outPredS,outScore] = readXgbPredDir(inPath,assignThr,ambigThr,doScale)

    if ~exist('assignThr','var') || isempty(assignThr)
        assignThr = 0.9;
    end

    if ~exist('ambigThr','var') 
        ambigThr = 0.5;
    end
    
    if ~exist('doScale','var') 
        doScale = 0;
    end
    
    if ~isstruct(inPath)
        if exist(inPath,'dir')
            inPath = regexprep(inPath,'/*$','');
            [zFile,~,zName] = fileList(inPath,'/*_pred.csv');

            %%

            zName = regexprep(zName,'([^_]*)_xgb_pred.csv','$1');

            %%

            [~,zi] = sort(zName);
            zName = zName(zi);
            zFile = zFile(zi);

            %%

            [zData,~,zH,zC] = fastMatRead(zFile{1},',');

            %%

            zN = length(zName);

            outScore.sampleID = zH;
            outScore.predMat = nan(length(zH),zN);
            outScore.predMat(:,1) = zData;
            outScore.typeH(1) = zC;
            outScore.fName = zName;

            progressbar()
            for i = 2:zN
                [zData,~,zH,zC] = fastMatRead(zFile{i},',');

                assert(isequal(zH,outScore.sampleID));

                outScore.predMat(:,i) = zData;        
                outScore.typeH(i) = zC;
                progressbar((i-1)/(zN-1));
            end
        elseif exist(inPath,'file')
            % [outScore.predMat,~,outScore.sampleID,outScore.typeH] = fastMatRead('/ahg/regevdata/projects/sawyersProstateSC/results/2019_11_08_mouse_post_rev/outReClust/mmProst10x_lowLevel_xgb_normTPM_DEG.tsv')
            % [outScore.predMat,~,outScore.sampleID,outScore.typeH] = fastMatRead('/ahg/regevdata/projects/sawyersProstateSC/results/2019_11_08_mouse_post_rev/outReClust/mmProst10x_lowLevel_xgb_normTPM_DEG.tsv')
            error('Needs fixuing');
            outScore.fName = outScore.typeH;
        else
            error('Inptu path:%s not found',inPath);
        end
    else
        fprintf('Using precalculated struct\n');
        outScore = inPath;
    end
    
    if doScale
        fprintf('Rescaling');
        if ~isfield(outScore,'predMatOrig')         
            outScore.predMatOrig = outScore.predMat;
        end
        
        outScore.predMat = outScore.predMatOrig./sum(outScore.predMatOrig,2);
    end
    
    [pScore,pOrd] = sort(outScore.predMat,2,'descend','MissingPlacement','last');
    
    outPred = cell(length(outScore.sampleID),1);
    outPred(:) = {'NA'};
    
    cSel = pScore(:,1) > assignThr;
    outPred(cSel) = outScore.typeH(pOrd(cSel));
    
    outScore.pred = outPred;
    outPredS.sampleID = outScore.sampleID;
    outPredS.pred = outPred;
    
    if ~isempty(ambigThr)
        cSel = cSel & pScore(:,2) > ambigThr;
        outPred = regexprepsub(outPred,cSel,'(.*)','ambig_$1');
        outScore.predAmbig = outPred;
        outPredS.predAmbig = outPred;
        
        outPred(cSel) = {'NA'};
        
        outScore.predAmbigNA = outPred;
        outPredS.predAmbigNA = outPred;

    end

end