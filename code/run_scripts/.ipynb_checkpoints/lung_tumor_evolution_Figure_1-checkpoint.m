% ---
% jupyter:
%   jupytext:
%     formats: ipynb,mfiles//m:percent
%     text_representation:
%       extension: .m
%       format_name: percent
%       format_version: '1.3'
%       jupytext_version: 1.4.0
%   kernelspec:
%     display_name: Matlab
%     language: matlab
%     name: matlab
% ---

% %%
%%

codeRoot = []
if isempty(codeRoot)
    codeRoot = matlab.desktop.editor.getActiveFilename();
    codeRoot = regexprep(codeRoot,'/[^/]*$','')    
end

cLoadPath = sprintf('%s/../utils/loadProjectPath.m',codeRoot);
run(cLoadPath)

cd(envVar.outDir);


% %%
mfilename('')

% %%
%% Data

envVar.reRunAnalysis = 0;

% %%
%% Data

mmLungPlate_fQC = load([ envVar.dataRoot 'mmLungKP_timecourse.mat' ]);

%%

timeID = mmLungPlate_fQC.timesimple;
zSel = strcmp(timeID,'02_KorKP_early_ND');
timeID(zSel) = mmLungPlate_fQC.typeID(zSel);
tabFilter(timeID);

% %%
%% Ref

cSubClustNmf = load([ envVar.dataRoot 'mmLungKP_timecourse_subClustNmf.mat' ])

% %%
%% Selecting over-dispersed genes 

cOpts.minGroup = 50;
cOpts.resampleN = 250;
cOpts.minFreq = 0.0500;
cOpts.minDup = 2;
cOpts.minExpGenes = 10;
cOpts.targetGene = 1000;
  

runSection = 0; 
if envVar.reRunAnalysis || runSection 

    outCV = detectGenesCV(mmLungPlate_fQC.rawCount,mmLungPlate_fQC.timesimple,[],cOpts);
    
    listVar = extractOverDispersed(outCV,cOpts.minFreq,cOpts.minDup);
else
    
    listVar = extractOverDispersed(cSubClustNmf.outCV,cOpts.minFreq,cOpts.minDup);
end 

% %%
%% CV vs Mean

expMean = mean(mmLungPlate_fQC.rawCount,2);
expCV = std(mmLungPlate_fQC.rawCount,[],2)./expMean;

figure;
gscatter(expMean,expCV,listVar,[],[],10);
set(gca,'xscale','log')

ylabel('Read-count Coefficient of variation (CV)');
xlabel('Read-count Mean');

% %%
%% Filter by freq and sex-chromosomes

geneChrMap = load([ envVar.dataRoot 'GRCm38_geneChrMap.mat' ])

%%
cEnsgID = regexprep(mmLungPlate_fQC.ensgID,'\..*$','');
listVarFilter = ismember(cEnsgID,geneChrMap.ensgID.chrX);
listVarFilter = listVarFilter | ismember(cEnsgID,geneChrMap.ensgID.chrY);
listVarFilter = listVarFilter | ismember(cEnsgID,geneChrMap.ensgID.chrM);
sum(listVarFilter)

listVarFilter = listVarFilter | ismember(mmLungPlate_fQC.geneSymbol,geneChrMap.symbol.chrX);
listVarFilter = listVarFilter | ismember(mmLungPlate_fQC.geneSymbol,geneChrMap.symbol.chrY);
listVarFilter = listVarFilter | ismember(mmLungPlate_fQC.geneSymbol,geneChrMap.symbol.chrM);
sum(listVarFilter)

%%

listVarAfterFilter = listVar;
sum(listVarAfterFilter)
listVarAfterFilter(listVarFilter) = 0;
sum(listVarAfterFilter)
sum(listVarAfterFilter&cSubClustNmf.listVar)/sum(listVarAfterFilter|cSubClustNmf.listVar)

% %%
%% Run NMF procedure on data 

runSection = 0;
if envVar.reRunAnalysis || runSection

    cOpts.type = 'concatTPM';
    cOpts.pcaMethod = 'NeNMF';
    cOpts.distMPower = 1;
    cOpts.doNorm = 7;
    cOpts.resampleN = 100;
    cOpts.batchID = 'typeID_num';
    cOpts.useODonly = 1;
    cOpts.saveExpMat = 0;
    cOpts.useGlobalFreq = 1;

    [~,hScore,wBase] = consensusGraphDistNMFConcat(mmLungPlate_fQC.normTPM(listVarAfterFilter,:),cOpts)

    %% Aggregate NMF laodings and generate a tSNE plot

    zH = cell2mat(hScore)';
else
    
    hScore = cSubClustNmf.hScore;
    wBase = cSubClustNmf.wBase;
    
    zH = cell2mat(hScore)';
    
end

% %%
%% Plot tSNE

runSection = 0;
if envVar.reRunAnalysis || runSection
    zT = tic;
    [ydata_tsne] = run_embed_mcTSNE(zH)
    toc(zT);
else
    ydata_tsne = mmLungPlate_fQC.annot.tSNE;
end

%%

plot_tsne_scatter(ydata_tsne,mmLungPlate_fQC.timesimple,colSet.cmapFull)

plot_tsne_scatter(ydata_tsne,mmLungPlate_fQC.annot.clusterK12,colSet.cmapClust)

% %%
%% Generate a phate map 

runSection = 0;
if envVar.reRunAnalysis || runSection

    zDataT = mmLungPlate_fQC.normTPM(listVarAfterFilter,:);

    zt = tic
    y_phate_2D = phate(zDataT','k',15,'npca',30,'pot_method','sqrt');
    toc
else
    y_phate_2D = mmLungPlate_fQC.annot.phate;
end

%%

plot_tsne_scatter(y_phate_2D,mmLungPlate_fQC.timesimple,colSet.cmapFull)

plot_tsne_scatter(y_phate_2D,mmLungPlate_fQC.annot.clusterK12)

% %%
%% Cluster

runSection = 0;
if envVar.reRunAnalysis || runSection

    zopts = cSubClustNmf.opts

    zopts.knn2dist
    %%
    zSimKnn = knndist2simMatrix(hScore {1}',zopts.knn2dist);
    for i = 2:length(hScore )
        zSimKnn = zSimKnn + knndist2simMatrix(hScore {i}',zopts.knn2dist);
    end

    zSimKnn = (zSimKnn + zSimKnn')./2;

    %%

    zT = 10.^[-1:0.2:1]

    [zStabClMerge.S, zStabClMerge.N, zStabClMerge.VI, zStabClMerge.C] = stability(zSimKnn,zT,'v','Laplacian','combinatorial','plot');

    %%

    zn = knee_pt(zStabClMerge.S)
    zStabClMerge.N(zn)


    %%

    zNewCl = zStabClMerge.C(:,zn);
    tabFilter(zNewCl)

    %%

    plot_crosstab_heatmap(zNewCl,mmLungPlate_fQC.annot.clusterK12)
end

% %%
%% Cluster reordering
if envVar.reRunAnalysis || runSection
    zTypeNum = str2double(regexprep(mmLungPlate_fQC.typeID,'^([0-9]+)_.*','$1'));
    tabFilter(zTypeNum);
    %%

    zopts = []
    zopts.equalFreqWeight = 1;
    zopts.weightedMedian = 1;
    zNewClReord = clusterReorderByVal(zNewCl,zTypeNum,[],zopts);

    %% Compare new clustering to reference (pub version)

    plot_crosstab_heatmap(zNewClReord,mmLungPlate_fQC.annot.clusterK12)
    
    zCl = zNewClReord;
else
    zCl = mmLungPlate_fQC.annot.clusterK12;
end

% %%
%% Figure 1C

zTimeSimple = mmLungPlate_fQC.timesimple;
zListType = unique(zTimeSimple)
zopts = [];
zopts.pSize = 20;
zopts.doPosTxt = 0;
zopts.newPlot = 0;

figure('Position',[1 1 3000 1800]);

for zi = 1:length(zListType)
    zSel = strgrep(zTimeSimple,zListType{zi});

    subplot(2,4,zi)
    plot_tsne_scatter(y_phate_2D,zTimeSimple,colSet.cmapFull,zSel,zopts);

    xlabel('Phate X'); ylabel('Phate Y');
end

% %%
%% Figure 1D 

zMouseID = mergeStringPair(regexprep(mmLungPlate_fQC.typeID,'([^_]+)_.*','$1'),mmLungPlate_fQC.mouseID);
zClTxt = mergeStringPair('%s%02d','Cl',zCl);

plot_crosstab_bar(zClTxt,timeID,colSet.cmapFull)

% %%
%% Figure 1E -- Normalized mutual infromation
% Selecting informative genes 
runSection = 0;
if envVar.reRunAnalysis || runSection
    % Todo: Add DEG with tweeDE model
    zDEG_all_ND = load([envVar.dataRoot 'mmLungKP_timecourse_plate_DEG_tweeDE.mat' ])

    %% 
    % Todo: Add wilcoxon based DEG


    zOutE = load([ envVar.dataRoot 'mmLungKP_timecourse_plate_DEG_ranksum.mat' ])

    %% 

    zDEG_pos = zDEG_all_ND.AUC > 0.65 & zDEG_all_ND.pval_adjust < 0.1 & zDEG_all_ND.log2fc > log2(1.25) & zDEG_all_ND.TRUE_freq > 0.1 & zOutE.maxLogR > log2(1);
    sum(zDEG_pos)
    sum(any(zDEG_pos,2))

    %%
    zopts.summaryGenes = 200;
    zopts.filterSubset = 10;   
    zopts.doWrite = 0;    
    [~,zSummaryTable] = extractWriteDEGtableByFlatIn(zOutE,mmLungPlate_fQC.geneSymbol,zDEG_pos,zDEG_all_ND.AUC,[],zopts)

    %% 

    zGeneList = table2cell(zSummaryTable);
    zGeneList = zGeneList(:);
    zGeneList(isemptycell(zGeneList)) = [];
    zGeneListAll = zGeneList;

    %%

    zSelGenes = ismember(mmLungPlate_fQC.geneSymbol,zGeneListAll);
    sum(zSelGenes)

    %%

    zTP = mmLungPlate_fQC.timesimple;
    tabFilter(zTP)

    %%



    zopts = [];
    zopts.sampleNum = 100;
    zopts.repNum = 50;
    zopts.dType = 4;
    zopts.nPC = 50;

    [zDistTP.phenoV,zDistTP.clustNames] = cluster_phenVolume(mmLungPlate_fQC.normTPM(zSelGenes,:),zTP,zopts);
else 
    
    zDistTP = load([ envVar.dataRoot 'mmLungKP_timecourse_NMI_TP.mat' ])
end


% %%
%% Figure 1E

[~,zi] = sort(zDistTP.clustNames)
zDistTPord = structSortMat(zDistTP,zi)
zClmat = repmat(zDistTPord.clustNames',size(zDistTPord.phenoV,1),1);


figure('Position',[1 1 1007 635])
g = gramm('x',zClmat(:),'y',zDistTPord.phenoV(:));

g.stat_boxplot();
g.set_text_options('base_size',18);
g.set_names('x','Time-point','y','NMI')
zg = g.draw()


zg.facet_axes_handles.YLim(1) = 0;

zg.facet_axes_handles.XTickLabelRotation = 60;
% %%
%% Figure 1F

zCl = mergeStringPair('Cl%s%02d','',mmLungPlate_fQC.annot.clusterK12);
zMouseID = mergeStringPair(regexprep(mmLungPlate_fQC.typeID,'([^_]+)_.*','$1'),mmLungPlate_fQC.mouseID);

%%

plot_crosstab_bar(zMouseID,zCl,colSet.cmapClust)

%% Figure 1G 

% Todo: add inferCNV 



