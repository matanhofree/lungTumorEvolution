%% Make for fastTable Read function 

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

mex -v /Users/mhofree/regev/local/toolbox-matlab/MIToolbox/fastToVect.c -largeArrayDims -I/usr/local/include/ -outdir /Users/mhofree/regev/local/toolbox-matlab/MIToolbox/

%%

mex -v miAllPairsMex.c MutualInformation.c Entropy.c CalculateProbability.c ArrayOperations.c

%%

output = miAllPairs(randi(5,5,10));

%% 

mex -v miAllPairsMexSub.c MutualInformation.c Entropy.c CalculateProbability.c ArrayOperations.c

%%

output = miCondAllPairs(randi(5,5,10),randi(5,5,10));

