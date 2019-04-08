% Copyright 2018 by Qiyuan Pang, Ken Ho, and Haizhao Yang

% Construct sparse embedding
%
%       [   U   ]
%   A = [V    -I]
%       [  -I  S]
%
% where S itself is expanded in the same form.
function S = BF_IDBF_spembed(F)
    L = length(F.U);
    assert(length(F.V) == L)  % check size

    % total number of nonzeros
    nz = nnz(F.S);
    for l = 1:L
        nz = nz + nnz(F.U{l}) + nnz(F.V{l}) + size(F.U{l},2) + size(F.V{l},1);
    end

    % allocate sparse storage and fill
    I = zeros(nz,1);
    J = zeros(nz,1);
    V = zeros(nz,1);
    M  = 0;
    N  = 0;
    nz = 0;
    for l = 1:L
        [Im,Jm,Vm] = find(F.U{l});
        [In,Jn,Vn] = find(F.V{l});
        [um,un] = size(F.U{l});
        [vm,vn] = size(F.V{l});
        % U block (1,2)
        nv = length(Vm);
        I(nz+(1:nv)) = M + Im;
        J(nz+(1:nv)) = N + vn + Jm;
        V(nz+(1:nv)) = Vm;
        nz = nz + nv;
        % V block (2,1)
        nv = length(Vn);
        I(nz+(1:nv)) = M + um + In;
        J(nz+(1:nv)) = N + Jn;
        V(nz+(1:nv)) = Vn;
        nz = nz + nv;
        % -I block (2,3)
        I(nz+(1:vm)) = M + um + (1:vm);
        J(nz+(1:vm)) = N + vn + un + (1:vm);
        V(nz+(1:vm)) = -ones(vm,1);
        nz = nz + vm;
        % -I block (3,2)
        I(nz+(1:un)) = M + um + vm + (1:un);
        J(nz+(1:un)) = N + vn + (1:un);
        V(nz+(1:un)) = -ones(un,1);
        nz = nz + un;
        M  = M  + um + vm;
        N  = N  + vn + un;
    end
    % S block (3,3)
    [Is,Js,Vs] = find(F.S);
    nv = length(Vs);
    I(nz+(1:nv)) = M + Is;
    J(nz+(1:nv)) = N + Js;
    V(nz+(1:nv)) = Vs;
    M = M + size(F.S,1);
    N = N + size(F.S,2);
    S = sparse(I,J,V,M,N);
end


