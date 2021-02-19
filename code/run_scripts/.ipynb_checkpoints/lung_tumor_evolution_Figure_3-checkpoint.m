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

zTimeSimple = mmLungPlate_fQC.timesimple;
zStageList = unique(zTimeSimple)

%%

zTestList = [];
zi = 1;

zTestList{zi}{1} = '01_T_early_ND'; %T - K2
zTestList{zi}{2} = '02_KorKP_early_ND';
zi = zi + 1;

% K only
zTestList{zi}{1} = '02_KorKP_early_ND'; %K2 - K12
zTestList{zi}{2} = '04_K_12w_ND'; 
zi = zi + 1;

 
zTestList{zi}{1} = '04_K_12w_ND'; %K12 - K30
zTestList{zi}{2} = '05_K_30w_ND';
zi = zi + 1;

 % KP Only
zTestList{zi}{1} = '02_KorKP_early_ND'; %K3 - K12
zTestList{zi}{2} = '06_KP_12w_ND';
zi = zi + 1;

zTestList{zi}{1} = '06_KP_12w_ND'; %K12 - K20
zTestList{zi}{2} = '07_KP_20w_ND';
zi = zi + 1;

zTestList{zi}{1} = '07_KP_20w_ND'; %K20 - K30w
zTestList{zi}{2} = '08_KP_30w_ND';
zi = zi + 1;

% %%
%% Note this takes a while if re-running

runSection = 0;
if envVar.reRunAnalysis || runSection
    
    zParam = loadjson([ envVar.dataRoot 'calcRunWOT_basic_plate_withTime.json' ]);

    
    zOutPath = './WOT_redo/'
    mkdir(zOutPath)

    otModelNew = wrapper_calcRunWOT(mmLungPlate_fQC,zOutPath,zParam,[])

else 

    [zFiles,~,zNames] = fileList([ envVar.dataRoot 'WOT/*' ])
    zNames = regexprep(zNames,'\.mat','')

    clear otModelNew
    for i = 1:length(zFiles)
        zOtSingle = load(zFiles{i})
        zf = setdiff(fieldnames(zOtSingle),{'outStat','clId','testSet','opts','sampleID'});
        %%
        for zj = 1:length(zf)
            zVar = zOtSingle.(zf{zj});

            if iscell(zVar)
                otModelNew.(zf{zj}){i} = zVar{1};
            else
                otModelNew.(zf{zj}){i} = zVar;
            end
        end

        otModelNew.clId(i) = zOtSingle.clId;

        otModelNew.outStat.coutMat{i} = zOtSingle.outStat.countMat{1};
        otModelNew.outStat.countMatSq{i} = zOtSingle.outStat.countMatSq{1};
        otModelNew.outStat.countDenom{i} = zOtSingle.outStat.countDenom{1};
        otModelNew.testSet{i} = zOtSingle.testSet{1};
    end

end

% %%
%%

zTimeSimple = mmLungPlate_fQC.timesimple;
tabFilter(zTimeSimple)

zSel = strcmp(mmLungPlate_fQC.timesimple,'01_T_early_ND');
zCluster = mmLungPlate_fQC.annot.clusterK12;
zCluster(zSel) = 0;

zopts = [];
zopts.minFrac = 0.1;
zopts.minClusterSize = 5;
zopts.normOut = 0;
zopts.normOutDir = 1;
zopts.connectMST = 0;

edgeList = ot_to_ancestorGraph_splitSummary(otModelNew.otMat,otModelNew.outStat,zCluster,otModelNew.testSet,zTimeSimple,zopts);

% %%
%plot -s 1500,1500
%% Generate a global graph

logit = @(x)log(x) - log(1-x);

edgeList_KP = edgeList;
zTab = struct2table(edgeList_KP);

% Focus on cell states
zSource = regexprep(edgeList_KP.sourceList,'.*_(cl.*)','$1');
zTarget = regexprep(edgeList_KP.targetList,'.*_(cl.*)','$1');

%
zSourceType = regexprep(edgeList_KP.sourceList,'(.*)_(cl.*)','$1');
% [zST,zN] = grp2idx(zSourceType);

zTargetType = regexprep(edgeList_KP.targetList,'(.*)_(cl.*)','$1');
% [zTT,zNT] = grp2idx(zTargetType);
% [zST,zN] = grp2idx(zTarget);
zTargetSortType = cellfun(@(x,y)sprintf('%s_%s',x,y),zSourceType,zTargetType,'uniformoutput',0);

%%

zTab = table(zSource,zTarget,zSourceType,zTargetType,zTargetSortType,edgeList_KP.valList);

%%

zTabSort = sortrows(zTab);
% [zTT,zNT] = grp2idx(zTabSort.zTargetType);
zNT = unique(zTabSort.zTargetType);
zMap = containers.Map(zNT,1:length(zNT));

zTT = nanvalues(zMap,zTabSort.zTargetType);
tabFilter(zTT)

%% Betweeness

zG = digraph(zTabSort.zSource,zTabSort.zTarget,zTabSort.Var6);

%%
clear zEtry;

zEtry.EndNodes = [ zTabSort.zSource zTabSort.zTarget];
zEtry.Weight = zTabSort.Var6;
zEtry.Names = zTabSort.zTargetSortType;
zEtry.TT = zTT;

zG = digraph(struct2table(zEtry))

%
figure('Position',[20 20 1220 1220]);
p = plot(zG,'layout','layered');

% Alteranative weighting scheme
% p.LineWidth = zG.Edges.Weight*20;
% zS = logit(zscore(zG.Edges.Weight));
% zS = zS - min(zS) + 1;

zS = zG.Edges.Weight;
[zX,zFx,zSout] = ecdfQuantile(zS);

p.LineWidth = zSout*10;

p.EdgeColor = cmapFull_TN(zG.Edges.TT,:);

zImp = zG.Edges.Weight;

wbc = centrality(zG,'pagerank','importance',zImp)

zV = round(wbc*150)+4;
p.MarkerSize = zV ;

p.NodeFontSize = 24
p.ArrowSize = 14;

% %%
% plot -s 2000,1200
%% 

for zi = 1:6

    zG = digraph(struct2table(zEtry))

    %
    figure('Position',[20 20 1220 1220]);
    p = plot(zG,'layout','layered');


    zS = zG.Edges.Weight;
    [zX,zFx,zSout] = ecdfQuantile(zS);

    p.LineWidth = zSout*10;

    zTT = zG.Edges.TT;
    
    cmapFull_mono = cmapFull_TN;
    cmapFull_mono(:) = 1;
    cmapFull_mono(zi,:) = cmapFull_TN(zi,:);

    p.EdgeColor = cmapFull_mono(zTT,:)

    % 
    zImp = zG.Edges.Weight;

    wbc = centrality(zG,'pagerank','importance',zImp)


    zV = round(wbc*150)+4;
    p.MarkerSize = zV ;

    p.NodeFontSize = 24
    p.ArrowSize = 14;

    title(sprintf('%s --> %s',zTestList{zi}{1},zTestList{zi}{2}),'Interpreter','none')

end



% %%
%%
zClList = unique(zTarget)

for zi = 1:12
%%
    digraph(struct2table(zEtry));
    %
    figure('Position',[20 20 1220 1220]);
    p = plot(zG,'layout','layered');

    zS = zG.Edges.Weight;
    [zX,zFx,zSout] = ecdfQuantile(zS);

    p.LineWidth = zSout*10;

    
    %%
    zEsel = strcmp(zG.Edges.EndNodes(:,1),zClList{zi});
    zSsel = strcmp(zG.Edges.EndNodes(:,2),zClList{zi});
    zRemE = zEsel | zSsel;
    %%
    p.EdgeColor = zeros(length(zRemE),3);
    p.EdgeColor(~zRemE,:) = 1
    %%
    p.EdgeColor(zSsel,:) = ones(sum(zSsel),3).*colorSet({'#bcbddc'});

    %
    zImp = zG.Edges.Weight;

    wbc = centrality(zG,'pagerank','importance',zImp)


    zV = round(wbc*150)+4;
    p.MarkerSize = zV ;

    p.NodeFontSize = 24
    p.ArrowSize = 14;

end


