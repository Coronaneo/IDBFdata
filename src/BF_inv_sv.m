% Copyright 2018 by Qiyuan Pang, Ken Ho, and Haizhao Yang

% Solve with LU factors of sparse embedding.
function y = BF_inv_sv(L,U,P,Q,x)
    n = size(x,1);
    x = [x; zeros(size(L,1)-n,size(x,2))];
    y = Q*(U\(L\(P*x)));
    y = y(1:n,:);
end