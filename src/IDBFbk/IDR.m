function [sk,T] = IDR(fun,xx,kk,rank_or_tol,srand)
  
  % set default parameters
  if nargin < 3 || isempty(srand)
    srand = 1;
  end

  % check inputs
  assert(rank_or_tol >= 0,'FLAM:id:negativeRankOrTol', ...
         'Rank or tolerance must be nonnegative.')

  % initialize
  m = size(xx,1);
  n = size(kk,1);

  % return if matrix is empty
  if isempty(xx) || isempty(kk)
    sk = [];
    rd = 1:n;
    T = zeros(0,n);
    return
  end

  % sample against Gaussian matrix if too rectangular
  if srand && m > 2*n
    A = randn(n+16,m)*fun(xx,kk);
  else
    A = fun(xx,kk);
  end

  % compute ID
  [~,R,E] = qr(A,0);
  if rank_or_tol < 1
    k = sum(abs(diag(R)) > abs(R(1))*rank_or_tol);
  else
    k = min(rank_or_tol,n);
  end
  sk = E(1:k);
  rd = E(k+1:end);
  T1 = R(1:k,1:k)\R(1:k,k+1:end);
  
  T = zeros(length(sk),n);
  
  for i = 1:length(rd)
      T(:,rd(i)) = T1(:,i);
  end
  
  for j = 1:length(sk)
      T(j,sk(j)) = 1;
  end
  
  sk = kk(sk,:);
end