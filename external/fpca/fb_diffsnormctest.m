%
% This script tests diffsnormc on dense matrices.
%

pass = [];

for m = [100 200]
  for n = [100 200]
    for isreal = [true false]

      if(isreal)
        A = randn(m,n);
      end
      if(~isreal)
        A = randn(m,n) + 1i*randn(m,n);
      end

      c = sum(A)/m;
      Ac = A - ones(m,1)*c;

      [U,S,V] = svd(Ac,'econ');
      snorm = diffsnormc(A,U,S,V);
      pass = [pass snorm<.1d-10*S(1,1)];

    end
  end
end


if(all(pass))
  disp('All tests succeeded.');
end

if(~all(pass))
  error('A test failed.');
end
