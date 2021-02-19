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

codeRoot = './code'
% if isempty(codeRoot)
%     codeRoot = matlab.desktop.editor.getActiveFilename();
%     codeRoot = regexprep(codeRoot,'/[^/]*$','')    
% end

cLoadPath = sprintf('%s/utils/loadProjectPath.m',codeRoot);
run(cLoadPath)

cd(envVar.outDir);

%% Settings

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
%% Load gene sets 

zCellHallmark = load([ envVar.dataRoot 'refData/geneSetCollectionMouse/msigDB_v62_symbol_homMM_h_all.mat'])
zCellMarker = load([ envVar.dataRoot 'refData/geneSetCollectionMouse/cellMarker_mouse_set_mm.mat'])
zMCAset = load([ envVar.dataRoot 'refData/geneSetCollectionMouse/mca_STable4_top.mat'])

%%
clear zCellSelected 

zCellSelected.AV2 = zMCAset.Alveolar_type_II_cell;
zCellSelected.AV1 = zCellMarker.s1147_pm24739965_Lung_Normal_cell_Type_I_pneumocyte_24739965
zCellSelected.Intestinal_epithelial_cell = zMCAset.Intestinal_epithelial_cell;
zCellSelected.Intestinal_TA = zCellMarker.s972_pm26287467_Intestinal_crypt_Normal_cell_Transit_amplifying;
zCellSelected.Hepatocyte = zMCAset.Hepatocyte;
zCellSelected.Mouse_embryonic_fibroblast = zMCAset.MEF_Cultured;
zCellSelected.Embryonic_stem_cell = zMCAset.Embryonic_stem_cell
zCellSelected.EMT = zCellHallmark.s30_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION

%%

runSection = 0;
if envVar.reRunAnalysis || runSection

    zopts = [];
    zopts.geneID = 'geneSymbol';
    zopts.bkgScoreSmpNum = 1000;
    zGStestCellSelected = genSignatureScoreSimple(mmLungPlate_fQC,zCellSelected,[],mmLungPlate_fQC.annot.clusterK12,zopts);

else
    zGStestCellSelected = load([ envVar.dataRoot 'mmLungKP_timecourse_GStestCellType.mat' ])
end
%%

zP = (struct2array(zGStestCellSelected.clustE) + 1)/(1001)

%

zSel = true(size(zP));
sum(zSel)
zSelSig = mafdr(zP,'BH',1) < 0.1

%

zCl = mergeStringPair('Cl%s%02d','',mmLungPlate_fQC.annot.clusterK12);

zGSscore = cell2mat(struct2cell(zGStestCellSelected.outScaledDiff)');
[zGSmean,~,zBid] = summarize_subset_value((zGSscore)',zCl)

%

zFlist = fieldnames(zGStestCellSelected.outScaledDiff);
zFlist = zFlist(zSel);

zGSmean = zGSmean(zSel,:);
zGSnorm = zscore(zGSmean,[],2)

%
fixNames = @(x)regexprep(x,'_','-');

% %%
%plot -s 1350 1050

zopts =[];

zopts.doSortX = 0;
zopts.doSortY = 0;
zopts.addDendrogramX =0
zopts.addDendrogramY = 0
zopts.doLeafOptimalOrderY = 0
zopts.doLeafOptimalOrderX = 0

zf = plot_heatmap_annot(zGSnorm,fixNames(zFlist),fixNames(zBid),[],[],zopts)

zf.Position = [81 81 1341 1052]

% %%
%% Generate Embryo gene signatures - Cao2019
% Based on supplementary tables of differentially expressed genes

runSection = 0;
if envVar.reRunAnalysis || runSection

    zCaoEmbryo = readtable([ envVar.dataRoot 'refData/Cao2019/41586_2019_969_MOESM3_ESM_S2.csv' ]);
    zCaoEmbryoTypeTable = readtable([ envVar.dataRoot 'refData/Cao2019/41586_2019_969_MOESM3_ESM_S3.csv' ]);

    %%

    zCaoEmbryoTypeTable = zCaoEmbryoTypeTable(1:38,:)
    zCaoEmbryoTypeTableMarker = arrayfun(@(x)strsplit(zCaoEmbryoTypeTable.MarkersUsedForCellTypeIdentification{x},', '),1:38,'unif',0);

    %% Global

    zSel = strgrep(zCaoEmbryo.Group,'Main');
    zCaoEmbryoSub = zCaoEmbryo(zSel,:);
    zCaoEmbryoSub = zCaoEmbryoSub(zCaoEmbryoSub.qval < 0.1,:);

    %%

    topN = @(x,n)x(1:min(length(x),n));

    %%

    zNames = mergeStringPair('Cao2019Embryo_main',matlab.lang.makeValidName(zCaoEmbryoTypeTable.CellType))

    %%

    tabFilter(zCaoEmbryoSub.max_cluster)
    [zClList,~,~,zCnt,zPos] = fastUnique(zCaoEmbryoSub.max_cluster);
    [zClList,zidx] = sort(zClList);
    zPos = zPos(zidx);
    zCnv = zCnt(zidx);

    %%
    zCaoGeneSetMainKnee = [];

    for zi = 1:40

         zCsel = zPos{zi};
         zSubGene = zCaoEmbryoSub.gene_short_name(zCsel);
         zSubStat = log2(zCaoEmbryoSub.max_expr(zCsel)./zCaoEmbryoSub.second_expr(zCsel));
         [zz,zidx] = sort(zSubStat,'descend','missingplacement','last');
         zSubGene = zSubGene(zidx);

        if length(zz)>5
            [zKneePt,zidx] = knee_pt(zz);
        else
            zKneePt=5;
            zidx = nan;
        end

        zTopList = topN(zSubGene,max(zidx,20));       
        zCaoGeneSetMainKnee{zi} = zTopList;
    end

    %%

    zSel = setdiff(1:40,[22 27])
    zCaoGeneSetMainKnee      = zCaoGeneSetMainKnee(zSel);

    for zi = 1:38
        zl = [ zCaoGeneSetMainKnee{zi}(:); zCaoEmbryoTypeTableMarker{zi}(:) ];
        length(zl)
        zl = fastUnique(zl);
        length(setdiff(zl,union(zCaoGeneSetMainKnee{zi},zCaoEmbryoTypeTableMarker{zi}(:)))) 
        zCaoGeneSetMainKneeAug{zi} = setdiff(zl,'NA');
    end

    zCaoGeneSetMainKnee = cell2struct(zCaoGeneSetMainKneeAug(:),zNames(:))
    
    %% Epithelial cells only
    
    zCaoS5 = readtable([ envVar.dataRoot 'refData/Cao2019/41586_2019_969_MOESM3_ESM_S5.csv' ]);

    zCaoS5_sub = zCaoS5(zCaoS5.Main_cluster_number == 6,:);
    zCaoS5_sub = sortrows(zCaoS5_sub,5);

    zCaoEpiMarker = arrayfun(@(x)strsplit(zCaoS5_sub.Minimum_defining_markers{x},','),1:29,'unif',0)

    zSel = strgrep(zCaoEmbryo.Group,'6');
    zCaoEmbryoSub = zCaoEmbryo(zSel,:);
    zCaoEmbryoSub = zCaoEmbryoSub(zCaoEmbryoSub.qval < 0.1,:);

    tabFilter(zCaoEmbryoSub.max_cluster)
    [zClList,~,~,zCnt,zPos] = fastUnique(zCaoEmbryoSub.max_cluster);

    [zClList,zidx] = sort(zClList);
    zPos = zPos(zidx);
    zCnv = zCnt(zidx);

    zCaoGeneSetEpiKnee = [];
    %

    for zi = 1:length(zClList)

         zCsel = zPos{zi};
         zSubGene = zCaoEmbryoSub.gene_short_name(zCsel);
         zSubStat = log2(zCaoEmbryoSub.max_expr(zCsel)./zCaoEmbryoSub.second_expr(zCsel));
         [zz,zidx] = sort(zSubStat,'descend','missingplacement','last');
         zSubGene = zSubGene(zidx);

        zCaoGeneSetEpiTopAll{zi} = zSubGene;
        if length(zz)>5
            [zKneePt,zidx] = knee_pt(zz);
        else
            zKneePt=5;
            zidx = nan;
        end

        zTopList = topN(zSubGene,max(zidx,20));       
        zCaoGeneSetEpiKnee{zi} = zTopList;

        fprintf('Cl %d) %f\t%f\t%f\t Knee=%f KneeIdx=%d/%d\n',zi,min(topN(zz,50)),min(topN(zz,100)),min(topN(zz,200)),zKneePt,zidx,length(zz));       
    end

    %%     
    
    zNames = mergeStringPair('Cao2019Embryo_epi',matlab.lang.makeValidName(zCaoS5_sub.Sub_cluster_name))
    zNames{end} =   'Cao2019Embryo_epi_NA2'

    for zi = 1:29
        zl = [ zCaoGeneSetEpiKnee{zi}(:); zCaoEpiMarker{zi}(:) ];
        length(zl)
        zl = fastUnique(zl);
        length(setdiff(zl,union(zCaoGeneSetEpiKnee{zi},zCaoEpiMarker{zi}(:)))) 
        zCaoGeneSetEpiKneeAug{zi} = setdiff(zl,'NA');
    end
    %% 
    zCaoGeneSetEpiKnee = cell2struct(zCaoGeneSetEpiKneeAug',zNames)
    
else
    zCaoGeneSetMainKnee = load([ envVar.dataRoot 'refData/embryoGeneSetsCao' ]);
    zCaoGeneSetEpiKnee = load([ envVar.dataRoot 'refData/embryoGeneSetsCaoEpi' ]);
end
% %%

%% Nowotschin2019, DEG were recalculated based on published data using a a Wilcoxon Ranksum test.

zDEG_embryoNowotschin2019 = load([ envVar.dataRoot 'refData/Nowotschin2019/DEG/mmEmbryo10x_DEG_cellTypeTP_by_slice.mat' ])

%%

zTPlist = setdiff(fieldnames(zDEG_embryoNowotschin2019),'geneID');
clear zEmbryoTopM

for j = 1:length(zTPlist)

    zTp = zTPlist{j}

    zDEG = zDEG_embryoNowotschin2019.(zTp);
    
    %%
    if isfield(zDEG,'maxLogR') && ~isempty(zDEG.maxLogR)
        zDEG_pos = zDEG.auc > 0.60 & zDEG.fdr < 0.1 & zDEG.logR > log2(1.25) & zDEG.expT > 0.1 & zDEG.maxLogR > log2(1);
        sum(zDEG_pos)
    else
        zDEG_pos = zDEG.auc > 0.60 & zDEG.fdr < 0.1 & zDEG.logR > log2(1) & zDEG.expT > 0.1;
        sum(zDEG_pos)
    end

    %%
    zopts.summaryGenes = 100;
    %zopts.filterSubset = 200;   
    zopts.doWrite = 0;    

    [~,zSummaryTable] = extractWriteDEGtableByFlatIn(zDEG,zDEG_embryoNowotschin2019.geneID,zDEG_pos,zDEG.auc,[],zopts);

    %% 
    zBatchID = mergeStringPair(zTp,matlab.lang.makeValidName(zDEG.outClustNames));

    for i = 1:length(zBatchID)
        zSub = table2cell(zSummaryTable(:,i));
        zSub(isemptycell(zSub)) = [];
        if ~isempty(zSub)
            zEmbryoTopM.(zBatchID{i}) = zSub;
        end
    end
    
end

%%

zEmbryoTopExtCao = zEmbryoTopM;
zEmbryoTopExtCao.Cao2019Embryo_main_EpithelialCells = zCaoGeneSetMainKnee.Cao2019Embryo_main_EpithelialCells;
zEmbryoTopExtCao.Cao2019Embryo_main_Hepatocytes = zCaoGeneSetMainKnee.Cao2019Embryo_main_Hepatocytes;

zEmbryoTopExtCao.Cao2019Embryo_epi_LungEpithelialTrajectory_1_of_1 = zCaoGeneSetEpiKnee.Cao2019Embryo_epi_LungEpithelialTrajectory_1_of_1;
zEmbryoTopExtCao.Cao2019Embryo_epi_Midgut_HindgutEpithelialTrajectory_1_of_3 = zCaoGeneSetEpiKnee.Cao2019Embryo_epi_Midgut_HindgutEpithelialTrajectory_1_of_3;
zEmbryoTopExtCao.Cao2019Embryo_epi_Midgut_HindgutEpithelialTrajectory_2_of_3 = zCaoGeneSetEpiKnee.Cao2019Embryo_epi_Midgut_HindgutEpithelialTrajectory_2_of_3;
zEmbryoTopExtCao.Cao2019Embryo_epi_Midgut_HindgutEpithelialTrajectory_3_of_3 = zCaoGeneSetEpiKnee.Cao2019Embryo_epi_Midgut_HindgutEpithelialTrajectory_3_of_3;
%%
zEmbryoTopExtCao

% %%
%%
runSection = 0;
if envVar.reRunAnalysis || runSection
    
    zopts = [];
    zopts.geneID = 'geneSymbol';
    zopts.bkgScoreSmpNum = 10000;

    zGStestEmbryoPeer = genSignatureScoreSimple(mmLungPlate_fQC,zEmbryoTopExtCao,[],mmLungPlate_fQC.annot.clusterK12,zopts);

    %%
else
    zGStestEmbryoPeer = load([ envVar.dataRoot 'mmLungKP_timecourse_GStest_embryo.mat' ])
    zGStestEmbryoPeer = zGStestEmbryoPeer.zGStestEmbryoPeer;
end

% %%
%plot -s 1600,1600
%% Plot heatmap for significant embryo signatures

zP = (struct2array(zGStestEmbryoPeer.clustE) + 1)/(10001)

zSel = mafdr(zP,'BH',1) < 0.1
sum(zSel)

zCl = mergeStringPair('Cl%s%02d','',mmLungPlate_fQC.annot.clusterK12);

zGSscore = cell2mat(struct2cell(zGStestEmbryoPeer.outScaledDiff)');
[zGSmean,~,zBid] = summarize_subset_value((zGSscore)',zCl)

%%

zFlist = fieldnames(zGStestEmbryoPeer.outScaledDiff);
zFlist = zFlist(zSel)
zGSmean = zGSmean(zSel,:);

fixNames = @(x)regexprep(x,'_','-');


[zFlistSort,zIdx] = sort(zFlist);
zIdx = [ zIdx(5:end); zIdx(3); zIdx(4); zIdx(2); zIdx(1)];
zIdx = flipud(zIdx)

zFixedList = zFlist(zIdx)
zGSnorm = zscore(zGSmean(zIdx,:),[],2)
[zz,zi] = sort(zGSnorm,'descend')


zz = unique([ setdiff(unique(zi(1,:)),6) 4 3 13 7 8])

zSel = trueV(zz,length(zFixedList))
zFixedList = zFixedList(zSel)
zGSnorm = zGSnorm(zSel,:);


%% 
zopts =[];

zopts.doSortX = 0;
zopts.doSortY = 0;
zopts.addDendrogramX =0
zopts.addDendrogramY = 0
zopts.doLeafOptimalOrderY = 0
zopts.doLeafOptimalOrderX = 0

zf = plot_heatmap_annot(zGSnorm,zFixedList,fixNames(zBid),[],[],zopts)
zf.Position = [87 1 1664 1704];

% %%
%plot -s 2000,450
%% Figure 2C

zGeneList = { 'Nkx2-1' 'Hnf4a' 'Hmga2' };

[~,zia,zib] = intersect(mmLungPlate_fQC.geneSymbol,zGeneList);
zord = zia(argsort(zib))

mmLungPlate_fQC.geneSymbol(zord)

%%

zExpM = mmLungPlate_fQC.normTPM(zord,:)';
zopts = [];
zopts.nRow = 1;
zopts.nCol = 3;
zopts.plotSize = [1 1 1999 443]
zopts.titleText = mmLungPlate_fQC.geneSymbol(zord);

zf = plot_tsne_scatter_multi(mmLungPlate_fQC.annot.phate,zExpM,[],[],zopts)

% %%
%plot -s 1600,1600
%% Figure 2D -- Consenus NMF programs 


runSection = 0;
if envVar.reRunAnalysis || runSection

    zOutP = [ envVar.outDir 'ccNMF_out/all_plate_ND' ];
    mkdir(fileparts(zOutP));

    zopts = [];

    % zopts.gsEnrichmentFile = '/ahg/regevdata/projects/Lung_ca_het/analysis/2018_11_14_plate_recluster/subCluster/geneSet_selected_mm.mat';
    % zopts.tfEnrichmentFile = '/ahg/regevdata/projects/Lung_ca_het/analysis/2018_11_14_plate_recluster/subCluster/TF_testSet_regnetworkweb.mat';

    zopts.optimalValSelect = -1;
    zopts.forceAll = 1
    zopts.enrichmentThr = 100;
    zopts.batchAwareCC = 0;

    zStub = 'mmLungPlate_ccNMF_top100';

    zopts.yData = mmLungPlate_fQC.annot.phate;
    zopts.Hnorm = 0;
    zopts.reWeightExpand = 1;

    zOutCCnmf = wrapper_ccNMF_enrichment(mmLungPlate_fQC,zSubClustNmf,zOutP,zopts,zStub,zClustIsT)
   
    zExpM = zOutCCnmf.bestConsH';
    
    zHeadType = mergeStringPair('NMF',size(zExpM,2));
    
    
else
    
    zOutCCnmf = load('/Users/mhofree/scratch/lungKP/2020_05_22_finalize_data/dataPrep/ccNMF_mmLungPlate_ccNMF_top100_byAlt_elbow_timesimple_outCons.mat')

    zHeadType = { 'AT1_AT2'
        'Cycling'
        'Low_conf_1'
        'Biosynthetic_mixed'
        'Stressed'
        'Highly_mixed'
        'EMT'
        'GI_epithelium_like'
        'Hepatocyte_like'
        'Low_conf_2'
        'AT2_like'}
    
    zExpM = zOutCCnmf.extrapH{4}';
    zopts = [];


end

zopts.titleText = mergeStringPair('NMF ',zHeadType);

zf = plot_tsne_scatter_multi(mmLungPlate_fQC.annot.phate,zExpM,[],[],zopts)

% %%
%plot -s 1600,1600
%% Figure 2E -- 

zGeneList = { 'Lyz2' 'Hopx' 'Cldn4' 'Cldn2' 'Zeb2' };

[~,zia,zib] = intersect(mmLungPlate_fQC.geneSymbol,zGeneList);
zord = zia(argsort(zib))

mmLungPlate_fQC.geneSymbol(zord)

%

zExpM = mmLungPlate_fQC.normTPM(zord,:)';
zopts = [];
zopts.nRow = 2;
zopts.nCol = 3;
zopts.plotSize = [1 1 1999 843]
zopts.titleText = mmLungPlate_fQC.geneSymbol(zord);

zf = plot_tsne_scatter_multi(mmLungPlate_fQC.annot.phate,zExpM,[],[],zopts)

% %%
%plot -s 1600,1600
%% Figure 2G

zopts = [];
zopts.maxRow = 12;
zopts.reorderX = 1;
zopts.widthV = 0.5;
zopts.widthBox = 0.1;
zopts.doJitter = 0.2;
zopts.dodge = 0;
zopts.cmap = cmapClust;
zopts.npoints = 30;
zCl = mergeStringPair('Cl%s%02d','',mmLungPlate_fQC.annot.clusterK12);

% plot_violin_simple(zCl,zOutCCnmf.extrapH{4}(6,:)',zCl,zopts)

plot_violin_simple(zCl,zOutCCnmf.bestConsH(6,:)',zCl,zopts)

% %%
%plot -s 1600,1600
%% Figure 2H

runSection = 0;
if envVar.reRunAnalysis || runSection
    % Todo: Add DEG with tweeDE model
    zDEG_all_ND = load([envVar.dataRoot 'mmLungKP_timecourse_plate_DEG_tweeDE.mat' ])

    %% 
    % Todo: Add wilcoxon based DEG

    zOutE = load([ envVar.dataRoot 'mmLungKP_timecourse_plate_DEG_ranksum.mat' ])    

    %%

    zDEG_pos = zDEG_all_ND.AUC > 0.6 & zDEG_all_ND.pval_adjust < 0.1 & zDEG_all_ND.log2fc > log2(1.25) & zDEG_all_ND.TRUE_freq > 0.05 & zOutE.maxLogR > log2(1);
    sum(zDEG_pos)
    sum(any(zDEG_pos,2))

    %%
    zopts = [];
    zopts.summaryGenes = 50;
    zopts.filterSubset = 10;   
    zopts.doWrite = 0;    
    [~,zSummaryTable] = extractWriteDEGtableByFlatIn(zOutE,mmLungPlate_fQC.geneSymbol,zDEG_pos,zDEG_all_ND.AUC,[],zopts)

    %% 

    zGeneListSt = table2cell(zSummaryTable);
    zGeneListSt = arrayfun(@(x)zGeneListSt(:,x),1:12,'unif',0);
    zGeneListSt = cell2struct(zGeneListSt',mergeStringPair('Cl%s%02d','',1:12))


else

    zGeneListSt = load([ envVar.dataRoot 'mmLungKP_timecourse_plate_DEG_top50_AUC.mat' ])
end

% %%
%plot -s 500,500
%% Gene signatures enrichment 
runSection = 0;
if envVar.reRunAnalysis || runSection
    
    zopts = [];
    zopts.geneID = 'geneSymbol';
    zopts.bkgScoreSmpNum = 1000;
    zGStestCellSelected = genSignatureScoreSimple(mmLungPlate_fQC,zGeneListSt,[],mmLungPlate_fQC.annot.clusterK12,zopts);
else
    zGStestCellSelected = load([ envVar.dataRoot 'mmLungKP_timecourse_GStest_clustSig.mat' ]);
end

%% 

figure;
hexscatter(zscore(zGStestCellSelected.outScore.Cl05),zscore(zOutCCnmf.extrapH{4}(6,:)'))

%%
% %%
%plot -s 2000,2000
%% Gene signatures enrichment 

zDataA = zscore(struct2array(zGStestCellSelected.outScore));
zDataB = zscore(zOutCCnmf.extrapH{4}');
zCa = mergeStringPair('Cl%s%02d','',1:12) 
zCb = zHeadType

%%

zopts = [];
zopts.hexPlot = 1;
zx = plot_corr_scatter_pairs(zDataA,zCa,zDataB,zCb,zopts)

%%

