function outCluster = wrapper_ccNMF_enrichment(dataFile,subClustDataFile,outPath,paramJson,outStubName,inCluster)
    
    totalTime = tic();

    defaultOpts.ncores = 1;
    defaultOpts.saveOut = 1;
    defaultOpts.loadInt = 1;
    
    defaultOpts.dataMat = 'rawCount';
    defaultOpts.loadVariables = { 'geneID' 'ensgID' 'sampleID' 'normTPM' 'batchID' }; % 
    defaultOpts.normTPM_val = 'normTPM';

    defaultOpts.loadSubClustDataVars = { 'sampleID' 'outConsensusNMF' 'listVar' 'ydata_NMF_hBase' };
    defaultOpts.doDiagPlots = 1;
    
    defaultOpts.geneSymbolVar = 'geneSymbol';
    defaultOpts.geneIDVar = 'geneID';
    
    defaultOpts.batchID = 'subClust.batchID';
    
    defaultOpts.batchAwareCC = 1;
    
    defaultOpts.optimalValSelect = [ 6 ]; % Use default;
    
    defaultOpts.forceAll = 0;
    
    defaultOpts.minK = 3;
    defaultOpts.maxK = 50;
    
    
    % defaultOpts.gsEnrichmentFile = '/ahg/regevdata/projects/colonSC/analysis/2018_08_04_c163_realign/geneSetCollection/geneSetMerged/geneSetAllMsigDB_hs.mat';     
    defaultOpts.gsEnrichmentFile = '';   
    defaultOpts.tfEnrichmentFile = '';
    
    defaultOpts.minSetSize = 3;
    defaultOpts.maxSetSize = [];
    
    defaultOpts.reExpandW = 1;
    defaultOpts.Hnorm = 2; % 1 - Unit norm 
                           % 2 - Sum to 1
                           % 0 - no change
       
    defaultOpts.yData = 'subClust.ydata_NMF_hBase';
    
    % defaultOpts.clusterFile = '';
    defaultOpts.setFigureOpts = 1;
    defaultOpts.figureFont = 16;
    
    defaultOpts.filterMinClust = 20;
    defaultOpts.filterMinClustVar = [];
    

    
    defaultOpts.clustVarName = [];
    defaultOpts.timeVarName = [];
    defaultOpts.recalcCCnmf = 0;
    
    defaultOpts.sortByTF = 'corr';
    defaultOpts.ccNMF = [];
    
    defaultOpts.versionCode = '0.0.5';
    
    defaultOpts.writeNMFraw = 1;
    defaultOpts.writeNMFcorr = 1;
    defaultOpts.forceRedo = 0;
    defaultOpts.forceRedoSub = [];
    defaultOpts.maxMahalH = 10000;
    
    outCluster = [];

    if exist('paramJson','var')
        if isstruct(paramJson)
            inOpts = paramJson;
        elseif exist(paramJson,'file')
            fprintf('Loading config file: %s\n',paramJson);
            inOpts = loadjson(paramJson);
        % else
            % error('Configuration file not found %s\n',paramsJson);
        end
    end
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end
    clear defaultOpts;
          
    if exist('varGenesFile','var')                    
        opts.detectGenesCV.useExternalList = varGenesFile;
    end
            
    if exist('smpFilterListFile','var')       
        opts.sampleFilter = smpFilterListFile;
    end       
    
    disp(opts);
    structdisp(opts);

    if opts.setFigureOpts
        set(0, 'DefaultFigureColor', 'White', ...
                'DefaultFigurePaperType', 'a4letter', ...
                'DefaultAxesColor', 'white', ...
                'DefaultAxesFontUnits', 'points', ...
                'DefaultAxesFontSize', opts.figureFont, ...
                'DefaultAxesFontName', 'Ariel', ...
                'DefaultAxesGridLineStyle', ':', ...
                'DefaultAxesInterruptible', 'on', ...
                'DefaultAxesLayer', 'Bottom', ...
                'DefaultAxesNextPlot', 'replace', ...
                'DefaultAxesUnits', 'normalized', ...
                'DefaultAxesXcolor', [0, 0, 0], ...
                'DefaultAxesYcolor', [0, 0, 0], ...
                'DefaultAxesZcolor', [0, 0, 0], ...
                'DefaultAxesVisible', 'on', ...
                'DefaultLineColor', 'Red', ...
                'DefaultLineLineStyle', '-', ...
                'DefaultLineLineWidth', 1, ...
                'DefaultLineMarker', 'none', ...
                'DefaultLineMarkerSize', 8, ...
                'DefaultTextColor', [0, 0, 0], ...
                'DefaultTextFontUnits', 'Points', ...
                'DefaultTextFontSize', opts.figureFont, ...
                'DefaultTextFontName', 'Ariel', ...
                'DefaultTextVerticalAlignment', 'middle', ...
                'DefaultTextHorizontalAlignment', 'left');
    end
    
    if ~exist('outStubName','var')
        outStubName = [];
    else
        fprintf('Using stub name:\n%s\n',outStubName);
    end

    %% Start pool
    % if opts.ncores > 1 && isempty(gcp('nocreate')) && ~ismac()
    %     maxCore = feature('numCores');
    %     parpool(min(maxCore,opts.ncores));
    % end

    %% Load data
    if ~isstruct(dataFile)        
        fprintf('Attempting to load data from file: %s\n',dataFile);
        
        if exist(dataFile,'file')
            loadVar = union(opts.loadVariables,opts.dataMat);
            % loadVar = union(loadVar,opts.batchVar);
            
            countData = loadSafe(dataFile,loadVar);

            if isempty(outStubName)
                [~,outStubName] = fileparts(dataFile);
            end
        else
            error('File not found');
        end
    else
        countData = dataFile;
        if isempty(outStubName)
            outStubName = inputname(1);
        end
    end
    
    % NormTPM
    if ~isfield(countData,'normTPM')
        fprintf('Calculating normTPM\n');
        if isequal(opts.normTPM_val,'normTPM')
            fprintf('Calculating normTPM from %s\n',opts.dataMat);
            countData.normTPM = calcNormTPM(countData.(opts.dataMat));
        else
            fprintf('Alternative normalized data - %s\n',opts.normTPM_val);
            countData.normTPM = countData.(opts.normTPM_val);
        end
    end

    %% Load ccNMF precalculated variables
    timeofrun = datetime;
    
    if ~isempty(subClustDataFile)
        if isstruct(subClustDataFile)
            subClust = subClustDataFile;
        else
            fprintf('Loading ccNMF data:\n',subClustDataFile)
            subClust = loadSafe(subClustDataFile,opts.loadSubClustDataVars);
        end
    else
        inPath = [ outPath filesep 'outClust_' outStubName '.mat' ];
        subClust = loadSafe(inPath,opts.loadSubClustDataVars);
    end
    
    if ~isequal(countData.sampleID,subClust.sampleID)
        zidx = ismember(countData.sampleID,subClust.sampleID);
        if ~all(zidx)
            countData = structSubSelectMat(countData,zidx);
        end
        
        if ~isequal(countData.sampleID,subClust.sampleID)
            error('Data does not match after subsetting unable to continue'); 
        end        
    end

    
    %% Start 
    fprintf('Starting\n');
    disp(countData);
    disp(subClust);
    
    %% Output path
    % if ~exist(outPath,'dir')
    outPath = [ outPath filesep outStubName '_ccNMF' ];        
    if ~exist(outPath,'dir')
        mkdir(outPath);
    end
    % end
    % outPath = [ outPath filesep 'ccNMF_' outStubName ];   
    
    N = length(countData.sampleID);
    D = length(countData.geneID);

    outFilePmat = [ outPath filesep 'ccNMF_' outStubName '_outCons.mat' ]; 

    if opts.batchAwareCC
        eval(sprintf('batchID = %s;',opts.batchID));
        fprintf('Batch info:\n');
        tabFilter(batchID);
    end
    
    %%
    reCalcCC = 0;
    if isfield(subClust,'outConsensusNMF') && opts.recalcCCnmf == 0 
        fprintf('Using consensus NMF from subClust object\n');
        outConsensusNMF = subClust.outConsensusNMF;
        outConsensusNMF.sampleID = subClust.sampleID;
    else      
        if exist(outFilePmat,'file') == 2 && opts.recalcCCnmf == 0               
            fprintf('Using out consensus from output directory:\n\t%s\n',outFilePmat);
            outConsensusNMF = load(outFilePmat);
        else
            reCalcCC = 1;
            
            if opts.batchAwareCC
                fprintf('Caclulating consensus NMF -- batch aware\n');            
                                
                [outConsensusNMF,cfig] = extractConsensusNMF_graphClust(subClust.wBase,countData.normTPM(subClust.listVar,:),batchID,opts.ccNMF);
                outConsensusNMF.sampleID = subClust.sampleID;    
            else
                fprintf('Caclulating consensus NMF\n');            
                [outConsensusNMF,cfig] = extractConsensusNMF_graphClust(subClust.wBase,countData.normTPM(subClust.listVar,:),[],opts.ccNMF);
                outConsensusNMF.sampleID = subClust.sampleID;   

            end 
            
            if ~isempty(cfig) 
                if ~isempty(cfig{1}) && isa(cfig{1},'matlab.ui.Figure')
                    outFileP = sprintf('%s/ccNMF_%s_kTopicDrop',outPath,outStubName);
                    print(cfig{1},outFileP,'-dpng');                    
                end
                if ~isempty(cfig{2}) && isa(cfig{2},'matlab.ui.Figure')
                    outFileP = sprintf('%s/ccNMF_%s_bestKSilPlot',outPath,outStubName);
                    print(cfig{2},outFileP,'-dpng');                    
                end
            end
        end
    end
    
    
    %%
    if opts.doDiagPlots
        
        outPlotFile = [ outPath filesep 'ccNMF_' outStubName '_summary' ];
        
        zfig = figure('Position',[0 0 1800 900]);
        
        yyaxis left 
        plot(outConsensusNMF.testK,outConsensusNMF.cMeanSil);
        xlabel('GP cluster K');
        ylabel('Mean Sil');

        yyaxis right
        plot(outConsensusNMF.testK,outConsensusNMF.reconstError);
        ylabel('Error');
        
        print(zfig,outPlotFile,'-dpng');
        
        close(zfig);
    end
    
    %%
    cDataM = countData.normTPM;
    
    if ~isfield(outConsensusNMF,'expandedW') || opts.reExpandW 
        outConsensusNMF.expandedW = ccNMFexpandBase(outConsensusNMF,cDataM,subClust.listVar,opts);           
        outConsensusNMF.geneID = countData.(opts.geneIDVar);
        if isfield(countData,opts.geneSymbolVar)
            outConsensusNMF.geneSymbol = countData.(opts.geneSymbolVar);
        end
        outConsensusNMF.listVar = subClust.listVar;
    end        
       
    %% Todo save intermediate file 
    if reCalcCC 
        save(outFilePmat,'-struct','outConsensusNMF');        
    end
    
    %%
    clustV = [];    
    timeP = [];
        
    if exist('inCluster','var')        
        if iscell(inCluster) || (ismatrix(inCluster) && isnumeric(inCluster))
            clustV = inCluster(:);
        else            
            if ~isstruct(inCluster) && exist(inCluster,'file')
                fprintf('Loading cluster file: %s\n',inCluster);
                inCluster = load(inCluster);
            end
            
            if isstruct(inCluster)
                [~,zia,zib] = intersect(outConsensusNMF.sampleID,inCluster.sampleID);
                
                if isempty(opts.clustVarName)
                    zFnames = setdiff(fieldnames(inCluster),'sampleID');            
                    fprintf('Using incluster struct -- %s\n',zFnames{1});
                    clustV = inCluster.(zFnames{1});                            
                else
                    fprintf('Using incluster struct -- %s\n',opts.clustVarName);
                    clustV = inCluster.(opts.clustVarName);
                end

                if ~isequal(zia,zib)
                    clustAll = clustV;
                    if iscell(clustV)
                        clustV = cell(N,1);
                        clustV(:) = { 'NA' };
                        clustV(zia) = clustAll(zib);
                    else
                        clustV = zeros(N,1);
                        clustV(zia) = clustAll(zib);
                    end
                end
                
                if ~isempty(opts.timeVarName)
                    timeP = inCluster.(opts.timeVarName);
                    
                    if ~isequal(zia,zib)
                        timePall = timeP;
                        if iscell(timeP)
                            timeP = cell(N,1);
                            timeP(:) = { 'NA' };
                            timeP(zia) = timePall(zib);
                        else
                            timeP = zeros(N,1);
                            timeP(zia) = timePall(zib);
                        end
                    end
                end                    
            else
                error('This should be a struct');
            end                        
        end
    end
    
    if isempty(clustV)
        fprintf('Using batch as cluster\n');
        clustV = subCluster.batchID;      
    end          
    %%
    %     if isempty(opts.gsEnrichmentFile) && isempty(opts.tfEnrichmentFile)
    %         fprintf('Completed. Nothing to test for enrichment\n');
    %         return;
    %     end
    
    gsTestSet = [];
    if ~isempty(opts.tfEnrichmentFile) 
        if ~isstruct(opts.gsEnrichmentFile)
            gsTestSet = load(opts.gsEnrichmentFile);
        else
            gsTestSet = opts.gsEnrichmentFile;
        end
    end
        
    tfTestSet = [];
    if ~isempty(opts.tfEnrichmentFile)
        tfTestSet = opts.tfEnrichmentFile;        
        if ~isstruct(tfTestSet)
            tfTestSet = load(opts.tfEnrichmentFile);
        else 
            tfTestSet = opts.tfEnrichmentFile;
            
        end
    end  
    
    if isnumeric(opts.yData)
        ydata = opts.yData;
    else
        eval(sprintf('ydata = %s;',opts.yData));
    end
        %%
  
    if opts.filterMinClust
        fprintf('Filtering minimal cluster assignment %d\n',opts.filterMinClust);
        
        if isempty(opts.filterMinClustVar)
            mCl = clustV;
        else
            eval(sprintf('mCl = %s;',opts.filterMinClustVar));
            % mCl = regexprep(clustV,'_T[0-9][0-9]_.*','');
        end
        %%
        [clList,~,~,clCnt,clPos] = fastUnique(mCl);

        clPosSub = concatCell(clPos(clCnt < opts.filterMinClust));

        clSub = ~trueV(clPosSub,length(mCl));
        
        outConsensusNMF = structSubSelectMat(outConsensusNMF,clSub,1);
        clustV = clustV(clSub);
        if size(ydata,1) == length(clSub)
            ydata = ydata(clSub,:);
        end
        countData = structSubSelectMat(countData,clSub);
        
        if ~isempty(timeP)
            timeP = timeP(clSub);
        end
    end

    if opts.forceAll
        cListOut = outConsensusNMF.testK;
        opts.optimalValSelect = -1*ones(length(cListOut));
    else
        cListOut = opts.optimalValSelect;
    end

    cSubList = [];
    for i = 1:length(cListOut)
      
        cVal = opts.optimalValSelect(i);
        [outConsensusNMF,bestIdx,kSelect,outSt] = selectOptType(outConsensusNMF,cVal,i,opts);

        if any(cSubList == kSelect)
            fprintf('Skipping %s -- already processed %d\n',outSt,kSelect);
            continue;
        end
        
        if kSelect > opts.maxK
            fprintf('Skipping %s -- too large %d\n',outSt,kSelect);
            continue;
        end
        %%        
        cSubList(end+1) = bestIdx;       
        
        if ~isempty(gsTestSet)            
            outFileP = [ outPath filesep 'ccNMF_' outStubName '_' outSt ];              
            if luniq(fileList(outFileP,'*'))>0 & ~ismember(kSelect,opts.forceRedoSub)
                fprintf('Found results from:\n%s\nskipping...\n',outFileP);
                continue;
            end            
            
            ccEnrichment.(outSt) = ccNMFtopicEnrichment(outConsensusNMF,outFileP,countData,gsTestSet,ydata,clustV,timeP,opts);                        
        end                
       
        if ~isempty(tfTestSet)            
            ropts = opts;
            ropts.doTF = 1;
            ropts.sortBy = opts.sortByTF;
            ropts.plotOverview = 0;
            outFileP = [ outPath filesep 'ccNMF_' outStubName '_TFcorr_' outSt ]; 
            ccTF.(outSt) = ccNMFtopicEnrichment(outConsensusNMF,outFileP,countData,tfTestSet,ydata,clustV,timeP,ropts)
        end 
        
        if opts.writeNMFcorr
            %%
            cBaseName = mergeStringPair('H%s%d','',1:kSelect);
            ropts = opts;
            ropts.corrShow = 'spe';                              
            [figOut,outMat,zOrdTX,zOrdTY] = plot_pairwise_corr(outConsensusNMF.bestConsH',cBaseName,ropts);            
            %%
            outFileP = sprintf('%s/ccNMF_%s_%s_H_speCorr',outPath,outStubName,outSt);
            print(figOut,outFileP,'-dpng');               

            %%
            cBaseName = mergeStringPair('W%s%d','',1:kSelect);
            ropts = opts;
            ropts.corrShow = 'spe';  
            [figOut,outMat,zOrdTX,zOrdTY] = plot_pairwise_corr(outConsensusNMF.consW,cBaseName,ropts);            
            %%
            outFileP = sprintf('%s/ccNMF_%s_%s_W_speCorr',outPath,outStubName,outSt);
            print(figOut,outFileP,'-dpng');               
            
            %%
            try      
                cBaseName = mergeStringPair('W%s%d','',1:kSelect);
                ropts = opts;
                ropts.corrShow = 'mahal';
                % [figOut,outMat,zOrdTX,zOrdTY] = plot_pairwise_corr(zOutH.Hinf(:,zSel)',zTypeNames,zopts);            
                % [zfig,outMat] = plot_heatmap_annot(outConsensusNMF.bestConsH',cBaseName,cBaseName,[],[],ropts);
                [figOut,outMat,zOrdTX,zOrdTY] = plot_pairwise_corr(outConsensusNMF.consW,cBaseName,ropts);            
                %%
                outFileP = sprintf('%s/ccNMF_%s_%s_W_mahalScov',outPath,outStubName,outSt);
                print(figOut,outFileP,'-dpng');  
            catch e 
                warning('Failed on mahal');
            end
               
            try 
                cBaseName = mergeStringPair('H%s%d','',1:kSelect);
                ropts = opts;
                ropts.corrShow = 'mahal';    
                    
                if opts.maxMahalH < N
                    fprintf('Subsampling cells for cov calculations\n');
                    
                    [figOut,outMat,zOrdTX,zOrdTY] = plot_pairwise_corr(outConsensusNMF.bestConsH(:,randsample(N,opts.maxMahalH))',cBaseName,ropts);                            
                else
                    [figOut,outMat,zOrdTX,zOrdTY] = plot_pairwise_corr(outConsensusNMF.bestConsH',cBaseName,ropts);            
                end
                
                outFileP = sprintf('%s/ccNMF_%s_%s_H_mahalScov',outPath,outStubName,outSt);
                print(figOut,outFileP,'-dpng');   
            catch e
                warning('Failed on mahal');
            end
            
        end
        
        if opts.writeNMFraw 
            fprintf('Writing raw NMF components');
  
            zW = outConsensusNMF.expandedW{bestIdx};            
            outFileP = sprintf('%s/ccNMF_%s_%s_W_top.csv',outPath,outStubName,outSt);             
            
            %% Todo addd the ccNMF indices. 
            %% if exists(ccEnrichment)
            %% end
            cBaseName = mergeStringPair('W%s%d','',1:kSelect);
            writeTableFile(zW(outConsensusNMF.listVar,:),outConsensusNMF.geneSymbol(outConsensusNMF.listVar),cBaseName,outFileP);            
            outFileP = sprintf('%s/ccNMF_%s_%s_W_full.csv',outPath,outStubName,outSt); 
            writeTableFile(zW,outConsensusNMF.geneSymbol,cBaseName,outFileP);
            %%
            zW = weightedNNcontrast(zW);   
            %%
            outFileP = sprintf('%s/ccNMF_%s_%s_W_weighted_top.csv',outPath,outStubName,outSt); 
            writeTableFile(zW(outConsensusNMF.listVar,:),outConsensusNMF.geneSymbol(outConsensusNMF.listVar),cBaseName,outFileP);                                                
            outFileP = sprintf('%s/ccNMF_%s_%s_W_weighted_full.csv',outPath,outStubName,outSt); 
            writeTableFile(zW,outConsensusNMF.geneSymbol,cBaseName,outFileP);
            %%
            cBaseName = mergeStringPair('H%s%d','',1:kSelect);
            zH = outConsensusNMF.bestConsH;
            zH = zH./sum(zH);
            outFileP = sprintf('%s/ccNMF_%s_%s_H_norm.csv',outPath,outStubName,outSt);
            writeTableFile(zH',outConsensusNMF.sampleID,cBaseName,outFileP);
            
        end
        
        try 
            close all 
        catch e            
            warning('No figures were found\n');
        end        
    end
    
%%
% outFileP = [ outPath filesep 'ccNMF_' outStubName ]; 
    if exist('ccEnrichment','var') && ~isempty(ccEnrichment) && exist('ccTF','var') && ~isempty(ccTF)
        save(outFilePmat,'-append','ccEnrichment','ccTF');
    elseif exist('ccEnrichment','var') && ~isempty(ccEnrichment)
        save(outFilePmat,'-append','ccEnrichment');
    elseif exist('ccTF','var') && ~isempty(ccTF)
        save(outFilePmat,'-append','ccTF');
    end
%     
    whosQ();
    toc(totalTime);
    fprintf('About to error exit due to matlab figure writing bug\n');
    if ~ismac()
        error('Exiting at complition (with error).');
    end
end

function [outConsensusNMF,bestIdx,kSelect,outSt] = selectOptType(outConsensusNMF,cVal,i,opts)

        switch cVal
            case 1
                fprintf('Using default selection\n');                
                outSt = 'default';    
            case 2              
                [~,bestIdx] = knee_pt(outConsensusNMF.reconstError,outConsensusNMF.testK);          
                fprintf('%% Using reconst error selection -- %d\n',outConsensusNMF.testK(bestIdx));
                outConsensusNMF.bestConsH = outConsensusNMF.extrapH{bestIdx};
                outConsensusNMF.consW = outConsensusNMF.mergedW{bestIdx};
                outSt = 'reconstError';    
            case 3     
                cSub = find(outConsensusNMF.testK > opts.minK & outConsensusNMF.testK < opts.maxK);
                
                [~,bestIdx] = max(outConsensusNMF.cMeanSil(cSub));
                bestIdx = cSub(bestIdx);
                fprintf('%% Using max mean sil -- %d\n',outConsensusNMF.testK(bestIdx));                
                
                outConsensusNMF.bestConsH = outConsensusNMF.extrapH{bestIdx};
                outConsensusNMF.consW = outConsensusNMF.mergedW{bestIdx};
                outSt = 'maxSil';                    
            case 4                
                cSub = find(outConsensusNMF.testK > opts.minK & outConsensusNMF.testK < opts.maxK);

                cRecError = cellfun(@(X,Y)norm(X*Y-cDataM,'fro'),outConsensusNMF.expandedW(cSub),outConsensusNMF.extrapH(cSub));                                
                
                [~,bestIdx] = knee_pt(cRecError,outConsensusNMF.testK(cSub));
                bestIdx = cSub(bestIdx);                
                outConsensusNMF.bestConsH = outConsensusNMF.extrapH{bestIdx};
                outConsensusNMF.consW = outConsensusNMF.mergedW{bestIdx};
                outSt = 'fullReconstError';                    
                
                fprintf('%% Using fullReoncstruction error -- %d\n',outConsensusNMF.testK(bestIdx));                
            case 5 
                %%
                [~,bestIdx] = min(abs(outConsensusNMF.testK - luniq(clustV)));
                
                outConsensusNMF.bestConsH = outConsensusNMF.extrapH{bestIdx};
                outConsensusNMF.consW = outConsensusNMF.mergedW{bestIdx};
                outSt = 'closestToClust'; 
                fprintf('%% UsingmclosestToClust -- %d\n',outConsensusNMF.testK(bestIdx));
            case 6
                
                if isfield(outConsensusNMF,'reconstErrorPostOpt')
                    [~,bestIdx] = knee_pt(outConsensusNMF.reconstErrorPostOpt,outConsensusNMF.testK);
                    outSt = 'reconstErrorPostOpt';  
                else
                    [~,bestIdx] = knee_pt(outConsensusNMF.reconstError,outConsensusNMF.testK);
                    outSt = 'reconstError';  
                end
                fprintf('%% Using reconst error selection -- %d\n',outConsensusNMF.testK(bestIdx));
                outConsensusNMF.bestConsH = outConsensusNMF.extrapH{bestIdx};
                outConsensusNMF.consW = outConsensusNMF.mergedW{bestIdx};
            case -1
                outSt = sprintf('forceTestAll%02d',i);
                
                fprintf('%% Force testing all cc clust outcomes -- %d\n',outConsensusNMF.testK(i));
                bestIdx = i;
                outConsensusNMF.bestConsH = outConsensusNMF.extrapH{i};
                outConsensusNMF.consW = outConsensusNMF.mergedW{i};
                
                
            otherwise           
                error('Not implemented');
    
        end
        
        kSelect = size(outConsensusNMF.bestConsH,1);
        outSt = sprintf('%s_K%02d',outSt,kSelect);
        
end