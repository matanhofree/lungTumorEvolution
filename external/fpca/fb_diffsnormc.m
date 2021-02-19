function snorm = diffsnormc(A,U,S,V,its)
%DIFFSNORMC  2-norm approx. error to a matrix upon centering.
%
%
%   snorm = diffsnormc(A,U,S,V)  computes an estimate snorm of the
%           spectral norm (the operator norm induced by the Euclidean
%           vector norm) of C(A) - USV', using 20 iterations of the
%           power method started with a random vector, where C(A) refers
%           to the matrix A from the input, after centering its columns.
%
%   snorm = diffsnormc(A,U,S,V,its)  finds an estimate snorm of the
%           spectral norm (the operator norm induced by the Euclidean
%           vector norm) of C(A) - USV', using its iterations of the
%           power method started with a random vector, where C(A) refers
%           to the matrix A from the input, after centering its columns;
%           its must be a positive integer.
%
%
%   Increasing its improves the accuracy of the estimate snorm of the
%   spectral norm of C(A) - USV', where C(A) refers to A after centering
%   its columns.
%
%
%   Note: DIFFSNORMC invokes RANDN. To obtain repeatable results, reset
%         the seed for the pseudorandom number generator.
%
%
%   inputs (the first four are required):
%   A -- first matrix in the column-centered A - USV' whose spectral
%        norm is being estimated
%   U -- second matrix in the column-centered A - USV' whose spectral
%        norm is being estimated
%   S -- third matrix in the column-centered A - USV' whose spectral
%        norm is being estimated
%   V -- fourth matrix in the column-centered A - USV' whose spectral
%        norm is being estimated
%   its -- number of iterations of the power method to conduct;
%          its must be a positive integer, and defaults to 20
%
%   output:
%   snorm -- an estimate of the spectral norm of the column-centered A
%            - USV' (the estimate fails to be accurate with
%            exponentially low probability as its increases; see
%            references 1 and 2 below)
%
%
%   Example:
%     A = rand(1000,2)*rand(2,1000);
%     A = A/normest(A);
%     [U,S,V] = pca(A,2,false);
%     diffsnormc(A,U,S,V)
%
%     This example produces a rank-2 approximation USV' to the column-
%     centered A such that the columns of U are orthonormal, as are the
%     columns of V, and the entries of S are all nonnegative and are
%     zero off the diagonal.
%     diffsnormc(A,U,S,V) outputs an estimate of the spectral norm
%     of the column-centered A - USV', which should be close to the
%     machine precision.
%
%
%   References:
%   [1] Jacek Kuczynski and Henryk Wozniakowski, Estimating the largest
%       eigenvalues by the power and Lanczos methods with a random
%       start, SIAM Journal on Matrix Analysis and Applications, 13 (4):
%       1094-1122, 1992.
%   [2] Edo Liberty, Franco Woolfe, Per-Gunnar Martinsson, Vladimir
%       Rokhlin, and Mark Tygert, Randomized algorithms for the low-rank
%       approximation of matrices, Proceedings of the National Academy
%       of Sciences (USA), 104 (51): 20167-20172, 2007. (See the
%       appendix.)
%   [3] Franco Woolfe, Edo Liberty, Vladimir Rokhlin, and Mark Tygert,
%       A fast randomized algorithm for the approximation of matrices,
%       Applied and Computational Harmonic Analysis, 25 (3): 335-366,
%       2008. (See Section 3.4.)
%
%
%   See also DIFFSNORM, PCA, NORMEST, NORM.
%

%   Copyright 2014 Mark Tygert.


%
% Check the number of inputs.
%
if(nargin < 4)
  error('MATLAB:diffsnormc:TooFewIn',...
        'There must be at least 4 inputs.')
end

if(nargin > 5)
  error('MATLAB:diffsnormc:TooManyIn',...
        'There must be at most 5 inputs.')
end

%
% Check the number of outputs.
%
if(nargout > 1)
  error('MATLAB:diffsnormc:TooManyOut',...
        'There must be at most 1 output.')
end

%
% Set the input its to its default value 20, if necessary.
%
if(nargin == 4)
  its = 20;
end

%
% Check the input arguments.
%
if(~isfloat(A))
  error('MATLAB:diffsnormc:In1NotFloat',...
        'Input 1 must be a floating-point matrix.')
end

if(isempty(A))
  error('MATLAB:diffsnormc:In1Empty',...
        'Input 1 must not be empty.')
end

if(~isfloat(U))
  error('MATLAB:diffsnormc:In2NotFloat',...
        'Input 2 must be a floating-point matrix.')
end

if(~isfloat(S))
  error('MATLAB:diffsnormc:In3NotFloat',...
        'Input 3 must be a floating-point matrix.')
end

if(isempty(S))
  error('MATLAB:diffsnormc:In3Empty',...
        'Input 3 must not be empty.')
end

if(~isfloat(V))
  error('MATLAB:diffsnormc:In4NotFloat',...
        'Input 4 must be a floating-point matrix.')
end

if(size(its,1) ~= 1 || size(its,2) ~= 1)
  error('MATLAB:diffsnormc:In5Not1x1',...
        'Input 5 must be a scalar.')
end

if(~(its > 0))
  error('MATLAB:diffsnormc:In5NonPos',...
        'Input 5 must be > 0.')
end

%
% Retrieve the dimensions of A, U, S, and V.
%
[m n] = size(A);
[m2 k] = size(U);
[k2 l] = size(S);
[n2 l2] = size(V);

%
% Make sure that the dimensions of A, U, S, and V are commensurate.
%
if(m ~= m2)
  error('MATLAB:diffsnormc:In1In2BadDim',...
        'The 1st dims. of Inputs 1 and 2 must be equal.')
end

if(k ~= k2)
  error('MATLAB:diffsnormc:In2In3BadDim',...
        'The 2nd dim. of Input 2 must equal the 1st dim. of Input 3.')
end

if(l ~= l2)
  error('MATLAB:diffsnormc:In3In4BadDim',...
        'The 2nd dims. of Inputs 3 and 4 must be equal.')
end

if(n ~= n2)
  error('MATLAB:diffsnormc:In1In4BadDim',...
        'The 2nd dim. of Input 1 must equal the 1st dim. of Input 4.')
end

%
% Calculate the average of the entries in every column.
%
c = sum(A)/m;


if(m >= n)

%
% Generate a random vector x.
%
  if(isreal(A) && isreal(U) && isreal(S) && isreal(V))
    x = randn(n,1);
  end

  if(~isreal(A) || ~isreal(U) || ~isreal(S) || ~isreal(V))
    x = randn(n,1) + i*randn(n,1);
  end

  x = x/norm(x);

%
% Run its iterations of the power method, starting with the random x.
%
  for it = 1:its
%
%   Set y = (A-ones(m,1)*c-USV')x.
%
    y = (x'*V)';
    y = S*y;
    y = U*y;
    y = A*x-ones(m,1)*(c*x)-y;
%
%   Set x = (A'-c'*ones(1,m)-VS'U')y.
%
    x = (y'*U)';
    x = (x'*S)';
    x = V*x;
    x = (y'*A-(y'*ones(m,1))*c)'-x;
%
%   Normalize x, memorizing its Euclidean norm.
%
    snorm = norm(x);
    if(snorm == 0)
      break;
    end
    x = x/snorm;
  end

  snorm = sqrt(snorm);

  clear x y;

end


if(m < n)

%
% Generate a random vector y.
%
  if(isreal(A) && isreal(U) && isreal(S) && isreal(V))
    y = randn(1,m);
  end

  if(~isreal(A) || ~isreal(U) || ~isreal(S) || ~isreal(V))
    y = randn(1,m) + i*randn(1,m);
  end

  y = y/norm(y);

%
% Run its iterations of the power method, starting with the random y.
%
  for it = 1:its
%
%   Set x = y(A-ones(m,1)*c-USV').
%
    x = y*U;
    x = x*S;
    x = (V*x')';
    x = y*A-(y*ones(m,1))*c-x;
%
%   Set y = x(A'-c'*ones(1,m)-VS'U').
%
    y = x*V;
    y = (S*y')';
    y = (U*y')';
    y = (A*x'-ones(m,1)*(c*x'))'-y;
%
%   Normalize y, memorizing its Euclidean norm.
%
    snorm = norm(y);
    if(snorm == 0)
      break;
    end
    y = y/snorm;
  end

  snorm = sqrt(snorm);

  clear x y;

end
