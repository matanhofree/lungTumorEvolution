%% Make for fastTable Read function 
cd /ahg/regevdata/users/mhofree/local/toolbox-matlab/fastTableRead
%%
% mex -v ~/regev/local/toolbox-matlab/fastTableRead/fastTableRead.cpp  -I/cellar/users/mhofree/projects/cancer_strat/external/clixo/clixo_release_0.1 -largeArrayDims -outdir /cellar/users/mhofree/projects/collab_cilia/code/util/fastTableRead/ CXXFLAGS="-fexceptions -fPIC -fno-omit-frame-pointer -pthread -O3 -lang-c-c++-comments -std=c++0x" COPTIMFLAGS="-O3 -DNDEBUG"

%% Mac 

% mex -v /Users/mhofree/regev/local/toolbox-matlab/fastTableRead/fastTableRead.cpp -largeArrayDims -I/usr/local/include/ -outdir /Users/mhofree/regev/local/toolbox-matlab/fastTableRead -L/usr/local/Cellar/boost/1.58.0/lib/ -lboost_system -lboost_thread
% mex -v /Users/mhofree/regev/local/toolbox-matlab/fastTableRead/fastTableRead.cpp -largeArrayDims -I/usr/local/opt/boost/include/ -I/usr/local/opt/boost/lib/ -outdir /Users/mhofree/regev/local/toolbox-matlab/fastTableRead -L/usr/local/opt/boost/lib/ CXXFLAGS="$CXXFLAGS -fexceptions -arch x86_64 -DMATLAB_BGL_LARGE_ARRAYS -fPIC -fno-omit-frame-pointer -pthread -O3 -lang-c-c++-comments -std=c++11" COPTIMFLAGS="-O3 -DNDEBUG"

% mex -v /Users/mhofree/regev/local/toolbox-matlab/fastTableRead/fastTableRead.cpp -largeArrayDims -outdir /Users/mhofree/regev/local/toolbox-matlab/fastTableRead -L/usr/local/opt/boost/lib/ CXXFLAGS="$CXXFLAGS -fexceptions -DMATLAB_BGL_LARGE_ARRAYS -fPIC -fno-omit-frame-pointer -pthread -O3 -lang-c-c++-comments -std=c++11"


%mex -v /Users/mhofree/regev/local/toolbox-matlab/fastTableRead/fastTableRead.cpp -largeArrayDims -outdir /Users/mhofree/regev/local/toolbox-matlab/fastTableRead -I/Users/mhofree/tmp/boost_1_49_0


mex -v /Users/mhofree/regev/local/toolbox-matlab/fastTableRead/fastTableRead.cpp -largeArrayDims -I/usr/local/include/ -outdir /Users/mhofree/regev/local/toolbox-matlab/fastTableRead 

% mex -v /Users/mhofree/regev/local/toolbox-matlab/fastTableRead/fastTableReadMac.cpp -largeArrayDims -outdir /Users/mhofree/regev/local/toolbox-matlab/fastTableRead

%% Linux 

mex -v /ahg/regev/users/mhofree/local/toolbox-matlab/fastTableRead/fastTableRead.cpp -I/home/unix/mhofree/local/boost_1_59/include -largeArrayDims  -outdir /ahg/regev/users/mhofree/local/toolbox-matlab/fastTableRead/
% CXXFLAGS="-fexceptions -fPIC -fno-omit-frame-p'ointer -pthread -O3 -lang-c-c++-comments -std=c++0x" COPTIMFLAGS="-O3 -DNDEBUG"

%%

mex -v /Users/mhofree/regev/local/toolbox-matlab/fastTableRead/fastTableReadV2.cpp -largeArrayDims -I/usr/local/include/ -outdir /Users/mhofree/regev/local/toolbox-matlab/fastTableRead -lboost_iostreams


%%


mex -v /ahg/regevdata/users/mhofree/local/toolbox-matlab/fastTableRead/fastTableReadV2.cpp -largeArrayDims -I/usr/local/include/ -outdir /ahg/regevdata/users/mhofree/local/toolbox-matlab/fastTableRead -I/ahg/regevdata/users/mhofree/local_unix/boost_install/boost_1_60_0/boost  -lboost_iostreams

%%


mex -v fastTableReadV2.cpp -largeArrayDims -I/usr/local/include/ -lboost_iostreams


%%

fname = '/ahg/regevdata/users/mhofree/local/toolbox-matlab/fastTableRead/test_data.small.txt';
[zdata,zerr] = fastTableReadV2(2,fname,[],char(9));
%%

fname = '/ahg/regevdata/users/mhofree/local/toolbox-matlab/fastTableRead/test_data.small.txt';
[zdata,zerr,zRowH,zColH] = fastTableReadV2(2,fname,[],char(9),0,1,1);

%%

fname = '/ahg/regevdata/users/mhofree/local/toolbox-matlab/fastTableRead/test_data.150by100.txt';
[zdata,zerr,zRowH,zColH] = fastTableReadV2(2,fname,[],char(9),0,1,1);

%%

fname = '/broad/hptmp/mhofree/active/kras_tumor_redo/analysis_b1to6_deNovo/summary/T1.mouse2.kallisto.genes.TPM.matrix';
[zdata,zerr,zRowH,zColH] = fastTableReadV2(2,fname,[],char(9),0,1,1);


%%
fname = '/ahg/regevdata/users/mhofree/local/toolbox-matlab/fastTableRead/test_data.small.txt';
[zz,zerr] = fastTxtRead(fname,'\t');

%% In linux UGER - build thus 
% The crux is to tstatically link libstdc++

/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE  -I/usr/local/include/ -I/ahg/regevdata/users/mhofree/local_unix/boost_1_59/include  -I"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2016a/extern/include" -I"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2016a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -std=c++11 -O3 -DNDEBUG /ahg/regevdata/users/mhofree/local/toolbox-matlab/fastTableRead/fastTableReadV2.cpp -o fastTableReadV2.o
/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/bin/g++ -pthread -Wl,--no-undefined  -shared -O3 -Wl,--version-script,"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2016a/extern/lib/glnxa64/mexFunction.map" -L/ahg/regevdata/users/mhofree/local_unix/boost_1_59/lib   -Wl,-rpath-link,/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2016a/bin/glnxa64 -L"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2016a/bin/glnxa64" -lmx -lmex -lmat -lm fastTableReadV2.o /broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.2.0/lib64/libstdc++.a -o /ahg/regevdata/users/mhofree/local/toolbox-matlab/fastTableRead/fastTableReadV2.mexa64
