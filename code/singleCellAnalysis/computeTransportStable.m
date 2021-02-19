function outTransport = computeTransportStable(p,q,C,lambda1,lambda2,epsilon,maxIter,g,verbose)
% simple unbalanced transport solver
% this solver has stabilized numerics
% inputs:
%         p,q:      uniform distributions on input and output cells
%         g:        growth value for input cells
%         C:        cost matrix to transport cell i to cell j (distance)
%     epsilon:      entropy parameter
%     maxIter:      number of scaling iterations
%     lambda1:       regularization parameter for marginal constraint for p.
%     lambda2:       regularization parameter for marginal constraint for q.

  p = p/sum(p);
  q = q/sum(q);

  p = p.*g;
  q = q*(sum(g)/length(g));

  u = zeros(length(p),1);
  v = zeros(length(q),1);
  a = ones(length(p),1);
  b = ones(length(q),1);
  
  K0 = exp(-C/epsilon);
  K = K0;
  alpha1 = lambda1/(lambda1+epsilon);
  alpha2 = lambda2/(lambda2+epsilon);
  
  a_old = a;
  b_old = b;
  i = 0;
  
  while (i < maxIter)
    i = i+1;
    % scaling iteration
    a = (p./(K*b)).^alpha1 .* (exp(-u/(lambda1+epsilon)));
    b = (q./(K'*a)).^alpha2 .* (exp(-v/(lambda2+epsilon)));
    
    % print(paste0("iter: ",i," a=",max(a)," b=",max(b)))
    ifprintf(verbose,'Iter %d: max(a)=%f, max(b)=%f\n',i,max(a),max(b));
    ifprintf(verbose,'\t\tdiff(a)=%f, diff(b)=%f\n',norm(a-a_old),norm(b-b_old));
    a_old = a;
    b_old = b;
    
    % stabilization
    if max(max(abs(a(:))),max(abs(b(:)))) > 1e6
      u = u + epsilon*log(a);
      v = v + epsilon*log(b); % absorb
      ifprintf(verbose,"Stabilization!\n")
      K = diag(exp(u/epsilon))*K0*diag(exp(v/epsilon));
      b = ones(length(q),1);
      a = ones(length(q),1);
    end
  end
  A = diag(a); 
  B = diag(b);
  outTransport = (A*K*B);
end
