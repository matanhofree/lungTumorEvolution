%%

cd ~/projects/cancer_SC/external/tSNE/multicore_tSNE_update/Multicore-TSNE/multicore_tsne

%%

mex -v -largeArrayDims mexTSNE_new.cpp splittree.cpp

%% 

open /ahg/regevdata/users/mhofree/projects/cancer_SC/external/tSNE/multicore_tSNE_update/Multicore-TSNE/multicore_tsne/fastMC_tSNE_new.m

%%

mex -v -largeArrayDims mexTSNE_new_step2.cpp splittree.cpp


% ############################# Old 
%% -O3 -fPIC -ffast-math -funroll-loops
%% GCC compile 

!/usr/local/Cellar/gcc/6.2.0/bin/g++-6 -c -DMATLAB_MEX_FILE  -I"/Applications/MATLAB_R2016a.app/extern/include" -I"/Applications/MATLAB_R2016a.app/simulink/include" -fno-common -arch x86_64 -std=c++11 -O3 -fPIC -ffast-math -funroll-loops -fwrapv -fopenmp -DNDEBUG /ahg/regevdata/users/mhofree/projects/cancer_SC/external/tSNE/multicore_tSNE/Multicore-TSNE/multicore_tsne/mexTSNE.cpp -o mexTSNE.o
%%
!/usr/local/Cellar/gcc/6.2.0/bin/g++-6 -c -DMATLAB_MEX_FILE  -I"/Applications/MATLAB_R2016a.app/extern/include" -I"/Applications/MATLAB_R2016a.app/simulink/include" -fno-common -arch x86_64 -std=c++11 -O3 -fPIC -ffast-math -funroll-loops -fwrapv -fopenmp -DNDEBUG /ahg/regevdata/users/mhofree/projects/cancer_SC/external/tSNE/multicore_tSNE/Multicore-TSNE/multicore_tsne/quadtree.cpp -o quadtree.o
%
!/usr/local/Cellar/gcc/6.2.0/bin/g++-6 -Wl,-twolevel_namespace -undefined error -arch x86_64 -bundle -Wl,-exported_symbols_list,"/Applications/MATLAB_R2016a.app/extern/lib/maci64/mexFunction.map" -O -Wl,-exported_symbols_list,"/Applications/MATLAB_R2016a.app/extern/lib/maci64/mexFunction.map" mexTSNE.o quadtree.o -L"/Applications/MATLAB_R2016a.app/bin/maci64" -L/usr/local/Cellar/gcc/6.2.0/lib/gcc/6/ -lmx -lmex -lmat -lgomp -o mexTSNE.mexmaci64


%% 

zMtag = reshape(zMtag(:),fliplr(size(zM)))';

%%

zz = zM(1:20,1:5);

%%
tic;
ydataTag = mexTSNE(zM',2,30,0.5,4,1000);
toc

%%

tic;
ydataTag = fastMC_tSNE(zM,2,30,0.5,4,3000);
toc


%%


zz = figure('position',[ 20 20 1020 920]);
set(gcf,'color','w'); 

gscatter(ydataTag(:,1),ydataTag(:,2),zBatch,[],[],15);
zz.Children(1).Interpreter = 'none';

ylabel('tSNE 2');
xlabel('tSNE 1');

legend('boxoff');
box off 

set(gca,'YTick',[-50 0 50]);

%%
ydataT = reshape(ydataTag(:),fliplr(size(ydataTag)))';

%%


zz = figure('position',[ 20 20 1020 920]);
set(gcf,'color','w'); 

gscatter(ydataT(:,1),ydataT(:,2),zBatch,[],[],15);
zz.Children(1).Interpreter = 'none';

ylabel('tSNE 2');
xlabel('tSNE 1');

legend('boxoff');
box off 

set(gca,'YTick',[-50 0 50]);

%%
zz = figure('position',[ 20 20 1020 920]);
set(gcf,'color','w'); 

gscatter(ydata(:,1),ydata(:,2),zBatch,[],[],15);
zz.Children(1).Interpreter = 'none';

ylabel('tSNE 2');
xlabel('tSNE 1');

legend('boxoff');
box off 

set(gca,'YTick',[-50 0 50]);

%%



clear zopts
zopts.methodType = 'fb_pca';
zopts.weightedPCA = 0;
zopts.initial_dims = 15;
zopts.perplexity = 20;
zopts.doNorm = 2;
zopts.doSmooth = 0;
% zopts.normFun = @(X,W)nnzcenter(X,nnzmean(X,1));
zopts.initial_solution = zM;
zopts.niter = 25;
zopts.useBHtsne = 1;
zopts.doWrite = 1;
%

% zidxCD = true(size(mm_colon_10x_fQC.batchID));
% zidxCD = ~strcmp(mm_colon_10x_fQC.batchID,'T1D_unsorted');
zT = tic;
% [ydata,zM,zW] = tsne_weighted(full(mm_colon_10x_fQC.normTPM(zSelect,zidxCD)'), [],[],mm_colon_10x_fQC.weightMat(zSelect,zidxCD)',zopts);
[ydata] = tsne_weighted([], [],[],-1,zopts);

toc(zT);

%% Linux version

!g++ -c -DMATLAB_MEX_FILE  -I"/Applications/MATLAB_R2016a.app/extern/include" -I"/Applications/MATLAB_R2016a.app/simulink/include" -fno-common -arch x86_64 -std=c++11 -O3 -fPIC -ffast-math -funroll-loops -fwrapv -fopenmp -DNDEBUG /ahg/regevdata/users/mhofree/projects/cancer_SC/external/tSNE/multicore_tSNE/Multicore-TSNE/multicore_tsne/mexTSNE.cpp -o mexTSNE.o
!g++ -c -DMATLAB_MEX_FILE  -I"/Applications/MATLAB_R2016a.app/extern/include" -I"/Applications/MATLAB_R2016a.app/simulink/include" -fno-common -arch x86_64 -std=c++11 -O3 -fPIC -ffast-math -funroll-loops -fwrapv -fopenmp -DNDEBUG /ahg/regevdata/users/mhofree/projects/cancer_SC/external/tSNE/multicore_tSNE/Multicore-TSNE/multicore_tsne/quadtree.cpp -o quadtree.o
%
!g++ -Wl,-twolevel_namespace -undefined error -arch x86_64 -bundle -Wl,-exported_symbols_list,"/Applications/MATLAB_R2016a.app/extern/lib/maci64/mexFunction.map" -O -Wl,-exported_symbols_list,"/Applications/MATLAB_R2016a.app/extern/lib/maci64/mexFunction.map" mexTSNE.o quadtree.o -L"/Applications/MATLAB_R2016a.app/bin/maci64" -L/usr/local/Cellar/gcc/6.2.0/lib/gcc/6/ -lmx -lmex -lmat -lgomp -o mexTSNE.mexmaci64

%% 

mex -v -largeArrayDims mexTSNE.cpp quadtree.cpp

%%

!/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -I"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2016b/extern/include" -I"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2016b/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11  -O3 -fPIC -ffast-math -funroll-loops -fwrapv -fopenmp -DNDEBUG /ahg/regevdata/users/mhofree/projects/cancer_SC/external/tSNE/multicore_tSNE/Multicore-TSNE/multicore_tsne/mexTSNE.cpp -o mexTSNE.o
!/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -I"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2016b/extern/include" -I"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2016b/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11  -O3 -fPIC -ffast-math -funroll-loops -fwrapv -fopenmp -DNDEBUG /ahg/regevdata/users/mhofree/projects/cancer_SC/external/tSNE/multicore_tSNE/Multicore-TSNE/multicore_tsne/quadtree.cpp -o quadtree.o
!/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -I"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2016b/extern/include" -I"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2016b/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11  -O3 -fPIC -ffast-math -funroll-loops -fwrapv -fopenmp -DNDEBUG /broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2016b/extern/version/cpp_mexapi_version.cpp -o cpp_mexapi_version.o
!/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/bin/g++ -pthread -Wl,--no-undefined  -shared -O -Wl,--version-script,"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2016b/extern/lib/glnxa64/c_exportsmexfileversion.map" mexTSNE.o quadtree.o scpp_mexapi_version.o   -Wl,-rpath-link,/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2016b/bin/glnxa64 -L"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2016b/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -lgomp -o mexTSNE.mexa64

%%

%% GCC compile 

!/usr/local/Cellar/gcc/7.1.0/bin/g++-7 -c -DMATLAB_MEX_FILE  -I"/Applications/MATLAB_R2016a.app/extern/include" -I"/Applications/MATLAB_R2016a.app/simulink/include" -fno-common -arch x86_64 -std=c++11 -O3 -fPIC -ffast-math -funroll-loops -fwrapv -fopenmp -DNDEBUG /ahg/regevdata/users/mhofree/projects/cancer_SC/external/tSNE/multicore_tSNE/Multicore-TSNE/multicore_tsne/mexTSNE.cpp -o mexTSNE.o
%%
!/usr/local/Cellar/gcc/7.1.0/bin/g++-7 -c -DMATLAB_MEX_FILE  -I"/Applications/MATLAB_R2016a.app/extern/include" -I"/Applications/MATLAB_R2016a.app/simulink/include" -fno-common -arch x86_64 -std=c++11 -O3 -fPIC -ffast-math -funroll-loops -fwrapv -fopenmp -DNDEBUG /ahg/regevdata/users/mhofree/projects/cancer_SC/external/tSNE/multicore_tSNE/Multicore-TSNE/multicore_tsne/quadtree.cpp -o quadtree.o
%%
!/usr/local/Cellar/gcc/7.1.0/bin/g++-7 -Wl,-twolevel_namespace -undefined error -arch x86_64 -bundle -Wl,-exported_symbols_list,"/Applications/MATLAB_R2016a.app/extern/lib/maci64/mexFunction.map" -O -Wl,-exported_symbols_list,"/Applications/MATLAB_R2016a.app/extern/lib/maci64/mexFunction.map" mexTSNE.o quadtree.o -L"/Applications/MATLAB_R2016a.app/bin/maci64" -L/usr/local/Cellar/gcc/7.1.0/lib/gcc/7/ -lmx -lmex -lmat -lgomp -o mexTSNE.mexmaci64


%% 
zopts.verbose = 0;
zy = run_embed_mcTSNE(rand(50,100),zopts);


%% Ubuntu 18.04 (on GCE)

 mex -v -largeArrayDims -lgomp CXXFLAGS="$CXXFLAGS -fopenmp" mexTSNE_new.cpp splittree.cpp