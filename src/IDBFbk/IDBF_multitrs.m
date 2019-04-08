function [tr1,tr2] = IDBF_multitrs(n,m,levels)
N = n*2^levels;
M = m*2^levels;
tr1 = zeros(N*M,1);
tr2 = zeros(N*M,1);
[tr,ts] = IDBF_multitr(levels);
for i = 1:n
    for j = 1:m
        tr1((i-1)*M*2^levels+(j-1)*2^(levels*2)+1:(i-1)*M*2^levels+j*2^(levels*2))=...
            (i-1)*2^levels+tr;
        tr2((i-1)*M*2^levels+(j-1)*2^(levels*2)+1:(i-1)*M*2^levels+j*2^(levels*2))=...
            (j-1)*2^levels+ts;
    end
end
end