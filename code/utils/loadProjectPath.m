
fprintf('Loading Lung Tumor Evolution pkg enviornment\n');

defEnvVar.pkgRoot = regexprep(mfilename('fullpath'),'.code.utils.loadProjectPath','');
defEnvVar.codeDir = [ defEnvVar.pkgRoot filesep 'code/' ];
defEnvVar.externalCodeDir = [ defEnvVar.pkgRoot filesep 'external/' ];
defEnvVar.dataRoot = [ defEnvVar.pkgRoot filesep 'data/' ];
defEnvVar.intRefDir = [ defEnvVar.dataRoot 'intRefDir/' ];
defEnvVar.outDir = [ defEnvVar.pkgRoot filesep 'output/' ];

if exist('envVar','var') == 1
    if exist('mergeOption')   
        envVar = mergeOption(envVar,defEnvVar);
    end        
else
    envVar = defEnvVar;
end
clear defEnvVar;  
disp(envVar)    

addpath([ envVar.codeDir '/utils']);
addPathExcl(envVar.codeDir);
addPathExcl(envVar.externalCodeDir);

%%

set_figure_defaults()

%%

colorSet = @(zCol)cell2mat(cellfun(@(x)hexColor(x(2:end)),zCol,'uniformoutput',0))';
fixNames = @(x)regexprep(x,'_','-');

clear colSet

colSet.rainbow14={'#882E72', '#B178A6', '#D6C1DE', '#1965B0', '#5289C7', '#7BAFDE', '#4EB265', '#90C987', '#CAE0AB', '#F7EE55', '#F6C141', '#F1932D', '#E8601C', '#DC050C'};
colSet.rainbow15={'#114477', '#4477AA', '#77AADD', '#117755', '#44AA88', '#99CCBB', '#777711', '#AAAA44', '#DDDD77', '#771111', '#AA4444', '#DD7777', '#771144', '#AA4477', '#DD77AA'};
colSet.rainbow18={'#771155', '#AA4488', '#CC99BB', '#114477', '#4477AA', '#77AADD', '#117777', '#44AAAA', '#77CCCC', '#777711', '#AAAA44', '#DDDD77', '#774411', '#AA7744', '#DDAA77', '#771122', '#AA4455', '#DD7788'};
colSet.rainbow21={'#771155', '#AA4488', '#CC99BB', '#114477', '#4477AA', '#77AADD', '#117777', '#44AAAA', '#77CCCC', '#117744', '#44AA77', '#88CCAA', '#777711', '#AAAA44', '#DDDD77', '#774411', '#AA7744', '#DDAA77', '#771122', '#AA4455', '#DD7788'};

colSet = structfun(@(x)colorSet(x),colSet,'unif',0);
%%

zGreens = brewermap(3,'Greens');
zBlue = brewermap(7,'Blues');
zRed = brewermap(9,'Reds');
% zCol = brewermap(8,'Dark2');
% zOranges = brewermap(3,'YlOrRd');

cmapFull = [];
zi = 1;
cmapFull(zi,:) = zGreens(3,:);
zi = zi+1;
cmapFull(zi,:) = zBlue(2,:);
zi = zi+1;
cmapFull(zi,:) = zRed(3,:);

zi = zi+1;
cmapFull(zi,:) = zBlue(4,:);
zi = zi+1;
cmapFull(zi,:) = zBlue(7,:);

zi = zi+1;
cmapFull(zi,:) = zRed(5,:);
zi = zi+1;
cmapFull(zi,:) = zRed(7,:);
zi = zi+1;
cmapFull(zi,:) = zRed(8,:);


colSet.cmapFull = cmapFull;


colSet.cmapClust = colSet.rainbow14(1:12,:);

%%

zGreens = brewermap(3,'Greens');
zBlue = brewermap(6,'Blues');
zRed = brewermap(6,'Reds');
zOranges = brewermap(3,'YlOrRd');

cmapFull_TN = [];
zi = 1;
cmapFull_TN(zi,:) = zGreens(2,:);

zi = zi+1;
cmapFull_TN(zi,:) = zBlue(2,:);
zi = zi+1;
cmapFull_TN(zi,:) = zBlue(4,:);

zi = zi+1;
cmapFull_TN(zi,:) = zOranges(1,:);
zi = zi+1;
cmapFull_TN(zi,:) = zOranges(2,:);
zi = zi+1;
cmapFull_TN(zi,:) = zOranges(3,:);

colSet.cmapClust_TN = cmapFull_TN;

clear cmapFull zi zGreens zBlue zRed cmapClust zOranges cmapClust_TN
