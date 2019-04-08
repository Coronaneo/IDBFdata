function [U] = IDBF_iterU1(Ucell,rv)
% called in IDBBF_arrg.m
%
% Copyright 2018 by Qiyuan Pang

prodnpx = size(Ucell,1);
prodnpk = size(Ucell,2);

if prodnpx*prodnpk > 1
   k = prodnpx*prodnpk/4;
   rv1 = rv(1:k,:);
   rv2 = rv(k+1:2*k,:);
   rv3 = rv(2*k+1:3*k,:);
   rv4 = rv(3*k+1:4*k,:);
   U1 = IDBF_iterU1(Ucell(1:prodnpk/2,1:prodnpx/2,:),rv1);
   U2 = IDBF_iterU1(Ucell(prodnpk/2+1:prodnpk,1:prodnpx/2,:),rv2);
   U3 = IDBF_iterU1(Ucell(1:prodnpk/2,prodnpx/2+1:prodnpx,:),rv3);
   U4 = IDBF_iterU1(Ucell(prodnpk/2+1:prodnpk,prodnpx/2+1:prodnpx,:),rv4);
   row = [size(U1,1) size(U2,1) size(U3,1) size(U4,1)];
   col = [size(U1,2) size(U2,2) size(U3,2) size(U4,2)];
   U = sparse(max(row(1),row(2))+max(row(3),row(4)),sum(col));
   U(1:row(1),1:col(1)) = U1;
   U(1:row(2),col(1)+1:col(1)+col(2)) = U2;
   U(max(row(1),row(2))+1:row(3)+max(row(1),row(2)),col(1)+col(2)+1:col(1)+col(2)+col(3)) = U3;
   U(max(row(1),row(2))+1:row(4)+max(row(1),row(2)),col(1)+col(2)+col(3)+1:col(1)+col(2)+col(3)+col(4)) = U4;
else
   lvl = size(Ucell,3);
   row = [];
   col = [];
   for i = 1:lvl
       row = [row size(Ucell{1,1,i},1)];
       col = [col size(Ucell{1,1,i},2)];
   end
   U = sparse(sum(row),sum(col));
   U(1:row(1),1:col(1)) = Ucell{1,1,1};
   for i = 2:lvl
       U(sum(row(1:i-1))+1:sum(row(1:i)),sum(col(1:i-1))+1:sum(col(1:i))) = Ucell{1,1,i};
   end
end
end