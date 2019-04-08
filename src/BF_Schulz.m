function X = BF_Schulz(Afun,Xfun,f,maxIter)
% use Schulz iteration to approximate matrix inverse
%
% Copyright 2018 by Qiyuan Pang, Ken Ho, and Haizhao Yang

if nargin < 4, maxIter = 5; end
if maxIter < 2
    X = 2*Xfun(f)-Xfun(Afun(Xfun(f)));
else
    Xfun2 = @(f) BF_Schulz(Afun,Xfun,f,maxIter-1);
    X = 2*Xfun2(f)-Xfun2(Afun(Xfun2(f)));
end
end