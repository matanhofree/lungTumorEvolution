function outWot = wrapper_calcRunWOT(dataFile,outPath,paramJson,externalTime,clId,outStubName)

    totalTime = tic();

    defaultOpts.ncores = 1;
    defaultOpts.saveOut = 1;
    defaultOpts.loadInt = 1;

    defaultOpts.dataMat = 'rawCount';
    defaultOpts.dataMatFilter = {};
    defaultOpts.sampleFilter = '';
    defaultOpts.sampleFilterReverseSelection = 0;
    defaultOpts.externalTime = '';
    defaultOpts.doNormTPM = 1;

    defaultOpts.keyVarList = '';

    defaultOpts.loadVariables = { 'geneID' 'ensgID' 'sampleID' 'batchID' 'cluster' 'normTPM'};

    defaultOpts.clusterSelect = {};
    
    defaultOpts.testSet = '';
% 
%     defaultOpts.minClustSize = 2;
%     defaultOpts.calcClusterStats = 1;
%     defaultOpts.calcClusterPval = 0;
% 
%     defaultOpts.saveSigScore = 1;
%     defaultOpts.saveSigDiffScore = 1;
% 
%     defaultOpts.saveClusterSummaryScore = 1;    
%     defaultOpts.saveClusterSummaryDiffScore = 1;    
% 
    defaultOpts.doSortOutScore = 1;
    defaultOpts.closePool = 1;
    
    defaultOpts.time_id_string = 'countData.batchID';

     defaultOpts.pcaMethod = 'NeNMF'; % Alternative: fb_pca, wpca, standard
    defaultOpts.seed = 42;
    defaultOpts.resampleNMF = 5;
 
    defaultOpts.niter = 25;
    defaultOpts.normFun = @(X)nnzcenter(X,nnzmean(X,1));

    defaultOpts.nstarts = 1;
    defaultOpts.fb_iter = 5;
    defaultOpts.verbose = 1;
    
%     defaultOpts.tsneMethod = 'resample'; % resampling - Resampling based approach to pca distance
%                                            % fast       - bh tsne
%                                            % cosine     - cosine similarity
%                                            % dist       - Use externa distance function 
%     
    defaultOpts.resampleN = 20;
    defaultOpts.resampleFrac = [ 0.8 0.8 ];

    defaultOpts.distType = 'logEucledian';
    defaultOpts.distMPower = 2;
    
    defaultOpts.otLambda1 = 1;
    defaultOpts.otLambda2 = 25;
    defaultOpts.otEpsilon = 0.1;
    defaultOpts.otNIter = 1000;
    defaultOpts.otGrowthRate = 1;
    

    defaultOpts.graphDist = 0;

    defaultOpts.pb = 0;
    defaultOpts.verbose = 0;
   
    defaultOpts.normCost = 1;
    defaultOpts.globalDist = 0;
    defaultOpts.pot_eps = 1e-6;
    
    defaultOpts.ncores = 1;
    
    %%

    defaultOpts.versionCode = '0.0.2';
    % savejson('',rmfield(defaultOpts,{'normFun'}),'/ahg/regevdata/projects/Lung_ca_het/analysis/2018_11_07_cohen2018/qsub_WOT/calcRunWOT_basic.json')
    if exist('paramJson','var')
        if isstruct(paramJson)
            inOpts = paramJson;
        elseif exist(paramJson,'file')
            fprintf('Loading config file: %s\n',paramJson);
            inOpts = loadjson(paramJson,'FastArrayParser',0,'SimplifyCell',1);
        end
    end
    if (exist('inOpts','var') == 1)
        opts = mergeOption(inOpts,defaultOpts);
    else
        opts = defaultOpts;
    end

    opts = removeString(opts);
    clear defaultOpts;
    disp(opts);
    % structdisp(opts);

    if nargin < 5
        clId = [];   
    end
    
    if ~isempty(clId) && ischar(clId)
        clId = str2double(clId);
        fprintf('Runnsing subset: %d\n',clId);
    end
        
    
    if ~exist('outStubName','var')
        outStubName = [];
    else
        fprintf('Using stub name:\n%s\n',outStubName);
    end

    %% Start pool
    if opts.ncores > 1 && isempty(gcp('nocreate'))
        maxCore = feature('numCores');
        parpool(min(maxCore,opts.ncores));
    end

    %% Load data
    if ~isstruct(dataFile)       
        fprintf('Attempting to load data from file: %s\n',dataFile)
        if exist(dataFile,'file')

            loadVar = union(opts.loadVariables,opts.dataMat);
            countData = loadSafe(dataFile,loadVar);

            if isempty(outStubName)
                [~,outStubName] = fileparts(dataFile);
            end
        else
            error('File not found');
        end
    else
        fprintf('Using workspace variable %s\n',inputname(1));
        loadVar = union(opts.loadVariables,opts.dataMat);
        countData = structSelectField(dataFile,loadVar);

        if isempty(outStubName)
            outStubName = inputname(1);
        end    
    end
    
    varGeneIdx = extractVarGeneList(countData,opts);
    if isempty(varGeneIdx)
        error('Variable gene list appears to be empty');
    end
    
    %% Filter batches
    if ~isempty(opts.dataMatFilter)
        zBatchVector = countData.(opts.batchVar);

        zDrop = false(size(zBatchVector));
        for i = 1:length(opts.dataMatFilter)
            disp(opts.dataMatFilter{i});
            zDrop = zDrop | strgrep(zBatchVector,opts.dataMatFilter{i});
        end

        fprintf(1,'Dropping %d samples due to channel filters\n',sum(zDrop));
        countData = structSubSelectMat(countData,~zDrop);
    end

    %% Filtering by id
    if ~isempty(opts.sampleFilter) && exist(opts.sampleFilter,'file')

        smpFilter = fastTxtRead(opts.sampleFilter);
        zDrop = ismember(countData.sampleID,smpFilter);

        if opts.sampleFilterReverseSelection == 1
            fprintf(1,'Selecting %d samples due sampleID filter\n',sum(zDrop));
            countData = structSubSelectMat(countData,zDrop);        
        else
            fprintf(1,'Dropping %d samples due sampleID filter\n',sum(~zDrop));
            countData = structSubSelectMat(countData,~zDrop);
        end
    end

    %% Output path
    if ~exist('outPath','var') || isempty(outPath)
        outPath = './';
    end
    if ~exist(outPath,'dir')
        mkdir(outPath);
    end

    outPathOrig = [ outPath filesep 'out_' outStubName ];

    %% Start
    fprintf('Starting\n');
    disp(countData);
    clusterSet = [];

    if ~exist('externalTime','var')  
         if ~isempty(opts.externalTime) 
             externalTime = opts.externalTime;
         else
             externalTime = load(externalTime);
         end
    end
    
    timeVar = [];
    if ~isempty(externalTime)
        timeVar = externalTime;
    else
        eval(sprintf('timeVar = %s;',opts.time_id_string));
    end
    
    % Begin processing
    [nD,N] = size(countData.(opts.dataMat));
    sampleID = countData.sampleID;

    if opts.doNormTPM
        fprintf(1,'Normalizing input matrix (%s) to log2 to TPX',opts.dataMat);

        countData.normTPM = calcNormTPM(countData.(opts.dataMat),[],10000);
        countData = rmfield(countData,opts.dataMat); 

        opts.inputCounts = 'normTPM';
    else    
        countData.(opts.dataMat) = countData.(opts.dataMat);
        % countData = rmfield(countData,opts.dataMat);    
        opts.inputCounts = opts.dataMat;
    end        

    if ~isempty(varGeneIdx)
        fprintf('Filtering by key variable genes.\n');               
        countData = structSubSelectMat(countData,varGeneIdx)
    else
        error('Must provide list of genes to filter');       
    end
    
    outPathId = sprintf('%s_%02d',outPath,clId);
    
    testSet = extractTestSet(timeVar,opts);
    if isempty(testSet)
        error('Ill defined test time points');
    end
    
    if isempty(clId)
        fprintf('Note: clId empty running everything.');
        testSetSub = testSet;
    else
        testSetSub = testSet(clId,:);
    end
    
    fprintf('Running the folloing time points in WOT:');
    disp(testSetSub)
    
    outWot = extractWOTbyTP(countData,outPathId,timeVar,testSetSub,clId,opts);    


    ppool = gcp('nocreate');
    if ~isempty(ppool) && opts.closePool 
        fprintf('Trying to delete parallel pool\n');
        delete(ppool);
    end

    fprintf('Total run time:%d\n',toc(totalTime));

end

function testSet = extractTestSet(timeVar,opts)

    if isempty(opts.testSet)
        zStageList = unique(timeVar)        

        clear zTestList
        for zi = 1:length(zStageList)-1
            testSet{zi}{1} = zStageList{zi};
            testSet{zi}{2} = zStageList{zi+1};    
        end
    else
        zStageList = unique(timeVar)        
        testSet = opts.testSet;
        
        checkMissing = cellfun(@(x)all(ismember(x,zStageList)),testSet);
        
        testSet(~checkMissing) = [];        
    end

end


function varGeneIdx = extractVarGeneList(countData,opts)

    [~,cName,ext] = fileparts(opts.keyVarList);

    if strcmp(ext,'.mat')            
        varGenes = load(opts.keyVarList)
                
        if isfield(varGenes,'listVar')
            varGeneIdx = varGenes.listVar;
            
            assert(length(varGeneIdx) == length(countData.geneID),'Error size missmatch');
            
        elseif isfield(varGenes,'geneList')
            varGeneIdx = ismember(countData.geneID,varGenes.geneID);
        else
            error('Variable gene param not recognized');            
        end
        
    elseif strcmp(ext,'.tsv') || strcmp(ext,'.csv') || strcmp(ext,'.txt')
        if strcmp(ext,'.tsv')
            varGenes = fastTxtRead(opts.keyVarList,'\t');            
        else            
            varGenes = fastTxtRead(opts.keyVarList,',');
        end
        
        varGeneIdx = ismember(countData.geneID,varGenes);
        
    else
        error('Var gene file not found');
    end 

end

function outWot = extractWOTbyTP(countData,outPath,timeVar,testSet,clId,opts)       

    %% Filter data 
    subSliceList = concatCell_to_cell(testSet);
    selIdxSet = ismember(timeVar,subSliceList);
    
    timeVar = timeVar(selIdxSet);
    countData = structSubSelectMat(countData,selIdxSet);
    
    %%
    if ~isfield(countData,'normTPM')
        countData = calcNormTPM(countData.rawCount);
    end
    
    
    [outWot.otMat,outWot.otVar,outWot.outStat] = consensusOptimalTransportLocalPar(countData.normTPM,[],timeVar,testSet,opts);
    outWot.timeVar = timeVar;
    outWot.testSet = testSet;
    outWot.clId = clId;
    outWot.opts = opts;
    outWot.sampleID = countData.sampleID;        
    
    %%      
                            
    if ~isempty(outPath)                
        save(outPath,'-v7.3','-struct','outWot');        
    end    
end

