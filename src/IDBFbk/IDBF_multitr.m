function [tr1,tr2] = IDBF_multitr(levels)
if  levels > 1
    m = 2^levels;
    n = m^2;
    tr1 = zeros(n,1);
    tr2 = zeros(n,1);
    [tr11,tr21] = IDBF_multitr(levels-1);
    tr1(1:n/4) = tr11;
    tr2(1:n/4) = tr21;
    
    tr1(n/4+1:n/2) = tr11;
    tr2(n/4+1:n/2) = tr21 + m/2;
    
    tr1(n/2+1:3*n/4) = tr11 + m/2;
    tr2(n/2+1:3*n/4) = tr21;
    
    tr1(3*n/4+1:n) = tr11 + m/2;
    tr2(3*n/4+1:n) = tr21 + m/2;
end

if  levels == 1
    m = 2^levels;
    n = m^2;
    tr1 = zeros(n,1);
    tr2 = zeros(n,1);
    tr1 = [1 1 2 2]';
    tr2 = [1 2 1 2]';
end

if  levels == 0
    tr1 = 1;
    tr2 = 1;
end

end