function [tr1,tr2] = IDBF_multitr2(levels)
if  levels > 1
    m = 2^levels;
    k = 2^(levels+1);
    n = m*k;
    tr1 = zeros(n,1);
    tr2 = zeros(n,1);
    [tr,ts] = IDBF_multitr2(levels-1);
    tr2(1:n/4) = ts;
    tr1(1:n/4) = tr;
    tr2(n/4+1:n/2) = m/2 + ts;
    tr1(n/4+1:n/2) = tr;
    tr2(n/2+1:3*n/4) = ts;
    tr1(n/2+1:3*n/4) = k/2 + tr;
    tr2(3*n/4+1:n) = m/2 + ts;
    tr1(3*n/4+1:n) = k/2 + tr;
end
if  levels == 1
    m = 2^levels;
    k = 2^(levels+1);
    n = m*k;
    tr1 = zeros(n,1);
    tr2 = zeros(n,1);
    tr2 = [1 1 2 2 1 1 2 2]';
    tr1 = [1 2 1 2 3 4 3 4]';
end
if  levels == 0
    tr2 = [1 1]';
    tr1 = [1 2]';
end