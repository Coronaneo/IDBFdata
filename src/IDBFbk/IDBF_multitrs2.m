function [tr1,tr2] = IDBF_multitrs2(n,m,levels)
N = n*2^levels;
M = m*2^(levels+1);
tr1 = zeros(N*M,1);
tr2 = zeros(N*M,1);
[tr,ts] = IDBF_multitr2(levels);
for i = 1:m
    for j = 1:n
        tr2((i-1)*N*2^(levels+1)+(j-1)*2^(2*levels+1)+1:(i-1)*N*2^(levels+1)+j*2^(2*levels+1))=...
            (j-1)*2^levels+ts;
        tr1((i-1)*N*2^(levels+1)+(j-1)*2^(2*levels+1)+1:(i-1)*N*2^(levels+1)+j*2^(2*levels+1))=...
            (i-1)*2^(levels+1)+tr;
    end
end
end