%% 

cd /ahg/regevdata/users/mhofree/projects/cancer_SC/results/bwpca

%%
cd ~/projects/cancer_SC/code/bwpca
mex -v /ahg/regevdata/users/mhofree/projects/cancer_SC/code/bwpca/bwpca.cpp ... 
    -I${HOME}/local/armadillo-6.400.3/include ... 
    -L${HOME}/local/armadillo-6.400.3/lib ...
    -L/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/lib64 ...
    -largeArrayDims -lopenblas -llapack -larmadillo -lstdc++... 
    CXXFLAGS="-std=c++11 -O2 -fPIC -fno-omit-frame-pointer -pthread -m64"

%%

mex -v /ahg/regevdata/users/mhofree/projects/cancer_SC/code/bwpca/bwpca.cpp ... 
    -I${HOME}/local/armadillo-6.400.3/include ...
    -L${HOME}/local/armadillo-6.400.3/lib .../usr/lib64/atlas/
    -L/usr/lib64/atlas/
    -largeArrayDims -larmadillo -lgfortran -lopenblas -llapac
    % CXXFLAGS="-std=c++11 -O2 -fPIC -fno-omit-frame-pointer -pthread -m64"
%%

mex -v /ahg/regevdatas/users/mhofree/projects/cancer_SC/code/bwpca/bwpca.cpp ... 
    -I${HOME}/local/armadillo-6.400.3/include ...
     -largeArrayDims -larmadillo -lgfortran -lopenblas -llapac
    % CXXFLAGS="-std=c++11 -O2 -fPIC -fno-omit-frame-pointer -pthread -m64"
    
%% Mac 

mex -v /ahg/regev/users/mhofree/projects/cancer_SC/code/bwpca/bwpca.cpp ... 
    -I/Users/mhofree/tmp/arma/armadillo-6.500.4/include ... 
    -L/Users/mhofree/tmp/arma/armadillo-6.500.4/ ... 
    -L/usr/local/Cellar/openblas/0.2.15/lib ...
    -largeArrayDims -lgfortran  -lopenblas -larpack -larmadillo -lstdc++... 
    CXXFLAGS="-std=c++11 -O2 -fPIC -fno-omit-frame-pointer -pthread -m64"

%%

mex -v /ahg/regevdata/users/mhofree/projects/cancer_SC/code/bwpca/bwpca.cpp ... 
    -I/Users/mhofree/tmp/arma/armadillo-6.500.4/include ... 
    -L/Users/mhofree/tmp/arma/armadillo-6.500.4/ ... 
    -L/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Versions/Current/ ...
    -L/usr/local/Cellar/arpack/3.2.0/lib/libarpack.dylib ...
    -largeArrayDims -lgfortran -larpack -lBLAS -lLAPACK -larmadillo 
      
%% g++ -c -I/opt/matlab/extern/include -I/opt/matlab/simulink/include -DMATLAB_MEX_FILE -ansi -fPIC -fno-omit-frame-pointer -pthread -DMX_COMPAT_32 -DNDEBUG -DARMA_NO_DEBUG "toy_example.cpp"

!xcrun  -sdk macosx10.11  clang++ -c  -I/Users/mhofree/tmp/arma/armadillo-6.500.4/include -I/Applications/MATLAB_R2014b.app/extern/include -I/Applications/MATLAB_R2014b.app/simulink/include -DMATLAB_MEX_FILE -std=c++11 -O2 -fPIC -fno-omit-frame-pointer -pthread -m64  -O2 -DNDEBUG  "/ahg/regev/users/mhofree/projects/cancer_SC/code/bwpca/bwpca.cpp"
!xcrun -sdk macosx10.11 clang -O -arch x86_64 -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk -mmacosx-version-min=10.11 -bundle -Wl,-exported_symbols_list,/Applications/MATLAB_R2014b.app/extern/lib/maci64/mexFunction.map -o  "bwpca.mexmaci64"  bwpca.o /usr/local/lib/libgfortran.a -L/Users/mhofree/tmp/arma/armadillo-6.500.4/ -L/usr/local/Cellar/openblas/0.2.15/lib -lopenblas -larpack -larmadillo -lstdc++ -L/Applications/MATLAB_R2014b.app/bin/maci64 -lmx -lmex -lmat -lstdc++


%%

!xcrun  -sdk macosx10.11  clang++ -c  -I/Users/mhofree/tmp/arma/armadillo-6.500.4/include -I/Applications/MATLAB_R2014b.app/extern/include -I/Applications/MATLAB_R2014b.app/simulink/include -DMATLAB_MEX_FILE -fno-common -fexceptions -arch x86_64 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk -mmacosx-version-min=10.11  -O2 -DNDEBUG  "/ahg/regevdata/users/mhofree/projects/cancer_SC/code/bwpca/bwpca.cpp"

!xcrun -sdk macosx10.11 clang -O -arch x86_64 -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk -mmacosx-version-min=10.11 -bundle -Wl,-exported_symbols_list,/Applications/MATLAB_R2014b.app/extern/lib/maci64/mexFunction.map -o  "bwpca.mexmaci64"  bwpca.o /usr/local/lib/libarpack.a  -L/Users/mhofree/tmp/arma/armadillo-6.500.4/ -L/System/Library/Frameworks/Accelerate.framework//Frameworks/vecLib.framework/Versions/Current/ -L/usr/local/lib/ -lgfortran -larpack -lBLAS -lLAPACK -larmadillo -L/Applications/MATLAB_R2014b.app/bin/maci64 -lmx -lmex -lmat -lstdc++

%%

mat = [ randn(10,5) randn(10,5)+5 ]; % matrix( c(rnorm(5*10,mean=0,sd=1), rnorm(5*10,mean=5,sd=1)), 10, 10)  # random matrix
matw = abs([ randn(10,5) randn(10,5)+5 ]);  % matrix( c(rnorm(5*10,mean=0,sd=1), rnorm(5*10,mean=5,sd=1)), 10, 10)  # random weight matrix
% ' matw <- abs(matw)/max(matw)
% ' base.pca.weighted <- bwpca(mat, matw)  # weighted pca

matw = ones(size(mat));
%%

mat = bsxfun(@minus,mat,sum(mat));

%% Mex free compilation

!/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/bin/g++ -c -DARMA_BLAS_LONG_LONG -D_GNU_SOURCE -DMATLAB_MEX_FILE  -I${HOME}/local/armadillo-6.400.3/include  -I"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2014b/extern/include" -I"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2014b/simulink/include" -std=c++11 -O2 -fPIC -fno-omit-frame-pointer -pthread -m64 -O -DNDEBUG /ahg/regevdata/users/mhofree/projects/cancer_SC/code/bwpca/bwpca.cpp -o bwpca.o
%%
!/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/bin/g++  -DARMA_BLAS_LONG_LONG -pthread -Wl,--no-undefined  -shared -O -Wl,--version-script,"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2014b/extern/lib/glnxa64/mexFunction.map" bwpca.o /broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/lib64/libstdc++.a   -larmadillo  -lopenblas  -llapack -L${HOME}/local/armadillo-6.400.3/lib  -Wl,-rpath,/broad/software/free/Linux/redhat_6_x86_64/pkgs/gcc_5.1.0/lib64/ -L"/broad/software/nonfree/Linux/redhat_6_x86_64/pkgs/matlab_2014b/bin/glnxa64" -lmx -lmex -lmat -lm -o bwpca.mexa64

%%

[zEigenv,zScores,zScoreWeights,zVar,totVar] = bwpca(mat,matw,5,1,0,1e-6,25,0,0)

%%
addpath('~/projects/cancer_SC/external/fpca/matlab/')
%%
[tU,tS,tV] = fb_pca(mat,5,1);


%%  Mac try 

!xcrun  -sdk macosx10.11  clang++ -c  -I/Users/mhofree/tmp/arma/armadillo-6.500.4/include -I/Applications/MATLAB_R2014b.app/extern/include -I/Applications/MATLAB_R2014b.app/simulink/include -DMATLAB_MEX_FILE -fno-common -fexceptions -arch x86_64 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk -mmacosx-version-min=10.11  -std=c++11 -O2 -fPIC -fno-omit-frame-pointer -pthread -m64 -O2 -DNDEBUG  "/ahg/regevdata/users/mhofree/projects/cancer_SC/code/bwpca/bwpca.cpp"
%%
        
!xcrun -v -sdk macosx10.11 clang -O -arch x86_64 -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk -mmacosx-version-min=10.11 -bundle -Wl,-exported_symbols_list,/Applications/MATLAB_R2014b.app/extern/lib/maci64/mexFunction.map -o  "bwpca.mexmaci64"  bwpca.o -L/Users/mhofree/tmp/arma/armadillo-6.500.4/ -L/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Versions/Current/  -lBLAS -lLAPACK -larmadillo -L/Applications/MATLAB_R2014b.app/bin/maci64 -lmx -lmex -lmat -lstdc++

%%

mex -v bwpca.cpp -lgfortran -largeArrayDims -larmadillo -L/usr/local/Cellar/armadillo/6.700.4/lib/ -I/usr/local/Cellar/armadillo/6.700.4/include

%%

! xcrun  -sdk macosx10.11  clang++ -c  -I/usr/local/Cellar/armadillo/6.700.4/include -I/Applications/MATLAB_R2014b.app/extern/include -I/Applications/MATLAB_R2014b.app/simulink/include -DMATLAB_MEX_FILE  -DARMA_BLAS_LONG_LONG -fno-common -fexceptions -arch x86_64 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk -mmacosx-version-min=10.11 -std=c++11 -O2 -fPIC -fno-omit-frame-pointer -pthread -m64  -O2 -DNDEBUG  "bwpca.cpp"
%%
! xcrun -sdk macosx10.11 clang++ -O -arch x86_64 -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk -mmacosx-version-min=10.11 -bundle -Wl,-exported_symbols_list,/Applications/MATLAB_R2014b.app/extern/lib/maci64/mexFunction.map -o  "bwpca.mexmaci64"  bwpca.o /usr/local/Cellar/arpack/3.2.0/lib/libarpack.a -larmadillo -L/usr/local/Cellar/armadillo/6.700.4/lib/ -L/Applications/MATLAB_R2014b.app/bin/maci64 -lmx -lmex -lmat -lstdc++ -L/Applications/MATLAB_R2014b.app/sys/os/maci64/ -lgfortran -framework Accelerate

%%

%%

mex -v -L/Users/mhofree/local/arma/armadillo-7.200.1_mex/usr/local/lib -I/Users/mhofree/local/arma/armadillo-7.200.1_mex/usr/local/include -L/usr/local/lib/ -larmadillo -lgfortran -largeArrayDims -DARMA_BLAS_LONG_LONG bwpca.cpp
%% 
setenv('DYLD_LIBRARY_PATH','/usr/local/lib/:/Applications/MATLAB_R2014b.app/sys/os/maci64:/Applications/MATLAB_R2014b.app/bin/maci64/../../Contents/MacOS:/Applications/MATLAB_R2014b.app/bin/maci64:/Applications/MATLAB_R2014b.app/extern/lib/maci64:/Applications/MATLAB_R2014b.app/runtime/maci64:/Applications/MATLAB_R2014b.app/sys/java/jre/maci64/jre/lib/./native_threads:/Applications/MATLAB_R2014b.app/sys/java/jre/maci64/jre/lib/./server:/Applications/MATLAB_R2014b.app/sys/java/jre/maci64/jre/lib/./lib/jli:/Applications/MATLAB_R2014b.app/sys/os/maci64:/Applications/MATLAB_R2014b.app/bin/maci64/../../Contents/MacOS:/Applications/MATLAB_R2014b.app/bin/maci64:/Applications/MATLAB_R2014b.app/extern/lib/maci64:/Applications/MATLAB_R2014b.app/runtime/maci64:/Applications/MATLAB_R2014b.app/sys/java/jre/maci64/jre/lib/./native_threads:/Applications/MATLAB_R2014b.app/sys/java/jre/maci64/jre/lib/./server:/Applications/MATLAB_R2014b.app/sys/java/jre/maci64/jre/lib/./lib/jli');
%%

!install_name_tool -change libarmadillo.7.dylib /Users/mhofree/local/arma/armadillo-7.200.1_mex/usr/local/lib/libarmadillo.7.dylib  bwpca.mexmaci64

%%
!install_name_tool -change /usr/local/lib/libgfortran.3.dylib /Users/mhofree/local/arma/armadillo-7.200.1_mex/gfortran/libgfortran.3.dylib  bwpca.mexmaci64

%%

!xcrun  -sdk macosx10.11  clang++ -c  -I/Users/mhofree/local/arma/armadillo-7.200.1_mex/usr/local/include -I/Applications/MATLAB_R2014b.app/extern/include -I/Applications/MATLAB_R2014b.app/simulink/include -DMATLAB_MEX_FILE -fno-common -fexceptions -arch x86_64 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk -mmacosx-version-min=10.11  -DARMA_BLAS_LONG_LONG -O2 -DNDEBUG  "bwpca.cpp"

%%
!xcrun -sdk macosx10.11 clang -O -arch x86_64 -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk -mmacosx-version-min=10.11 -bundle -Wl,-exported_symbols_list,/Applications/MATLAB_R2014b.app/extern/lib/maci64/mexFunction.map -o  "bwpca.mexmaci64"  bwpca.o  -L/Users/mhofree/local/arma/armadillo-7.200.1_mex/usr/local/lib -larmadillo -lgfortran -L/Applications/MATLAB_R2014b.app/bin/maci64 -lmx -lmex -lmat -lstdc++ -framework Accelerate -headerpad_max_install_names

% Note: To make this work on macosx I had to relink the malab gfortran
% libraries to the ones in /usr/local/lib/libgfortran.dylib

%% Updated 
xcrun  -sdk macosx10.11  clang++ -c  -I/ahg/regevdata/users/mhofree/local/external-matlab/mexplus/include -I/Users/mhofree/local/arma/armadillo-7.200.1_mex/usr/local/include -I/Applications/MATLAB_R2016a.app/extern/include -I/Applications/MATLAB_R2016a.app/simulink/include -DMATLAB_MEX_FILE -fno-common -fexceptions -arch x86_64 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk -mmacosx-version-min=10.11  -DARMA_BLAS_LONG_LONG -O2 -DNDEBUG -std=c++11  "bwpca.cpp"

xcrun -sdk macosx10.11 clang -O -arch x86_64 -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk -mmacosx-version-min=10.11 -bundle -Wl,-exported_symbols_list,/Applications/MATLAB_R2016a.app/extern/lib/maci64/mexFunction.map -o  "bwpca.mexmaci64"  bwpca.o  -L/Users/mhofree/local/arma/armadillo-7.200.1_mex/usr/local/lib -larmadillo -lgfortran -L/Applications/MATLAB_R2016a.app/bin/maci64 -lmx -lmex -lmat -lstdc++ -framework Accelerate -headerpad_max_install_names
