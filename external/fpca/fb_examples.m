% This script provides example usage of the codes for "An implementation
% of a randomized algorithm for principal component analysis."


%
% Set the dimensions of the m x n matrix to be approximated.
%
m = 100;
n = m;

%
% Set the rank of the approximation.
%
k = 10;

%
% Set the spectral-norm accuracy of the best rank-approximation.
%
best = 1e-12

%
% Construct the eigenvalues.
%
for j = 1:k
  d(j) = exp(log(best)*(j-1)/(k-1));
end
for j = (k+1):n
  d(j) = best*(n-j)/(n-k-1);
end

%
% Construct the eigenvectors.
%
V = randn(n);
[V,R] = qr(V,0);

%
% Form the matrix A to be approximated.
%
A = V*diag(d)*V';

%
% Construct rank-k approximations to A.
%
[V1,D1] = eigenn(A,k);
[V2,D2] = eigens(A,k);
[U,D3,V3] = pca(A,k,true);

%
% Check the accuracies of the approximations.
%
its = 10000;
errn = diffsnorms(A,V1,D1,its)
errs = diffsnorms(A,V2,D2,its)
erra = diffsnorm(A,V3,D3,V3,its)

if(errn > 100*best)
  error('errn is too big.');
end
if(errs > 100*best)
  error('errs is too big.');
end
if(erra > 100*best)
  error('erra is too big.');
end
