function [V] = IDBF_iterV1(Vcell,cv)
% called in IDBBF_arrg.m
%
% Copyright 2018 by Qiyuan Pang

prodnpx = size(Vcell,1);
prodnpk = size(Vcell,2);

if prodnpx*prodnpk > 1
   k = prodnpx*prodnpk/4;
   cv1 = cv(1:k,:);
   cv3 = cv(k+1:2*k,:);
   cv2 = cv(2*k+1:3*k,:);
   cv4 = cv(3*k+1:4*k,:);
   V1 = IDBF_iterV1(Vcell(1:prodnpk/2,1:prodnpx/2,:),cv1);
   V2 = IDBF_iterV1(Vcell(prodnpk/2+1:prodnpk,1:prodnpx/2,:),cv2);
   V3 = IDBF_iterV1(Vcell(1:prodnpk/2,prodnpx/2+1:prodnpx,:),cv3);
   V4 = IDBF_iterV1(Vcell(prodnpk/2+1:prodnpk,prodnpx/2+1:prodnpx,:),cv4);
   row = [size(V1,1) size(V2,1) size(V3,1) size(V4,1)];
   col = [size(V1,2) size(V2,2) size(V3,2) size(V4,2)];
   V = sparse(sum(row),max(col(1),col(2))+max(col(3),col(4)));
   V(1:row(1),1:col(1)) = V1;
   V(row(1)+1:row(1)+row(2),1:col(2)) = V2;
   V(row(1)+row(2)+1:row(1)+row(2)+row(3),max(col(1),col(2))+1:col(3)+max(col(1),col(2))) = V3;
   V(row(1)+row(2)+row(3)+1:row(1)+row(2)+row(3)+row(4),max(col(1),col(2))+1:col(4)+max(col(1),col(2))) = V4;
else
   lvl = size(Vcell,3);
   row = [];
   col = [];
   for i = 1:lvl
       row = [row size(Vcell{1,1,i},1)];
       col = [col size(Vcell{1,1,i},2)];
   end
   V(1:row(1),1:col(1)) = Vcell{1,1,1};
   for i = 2:lvl
       V(sum(row(1:i-1))+1:sum(row(1:i)),sum(col(1:i-1))+1:sum(col(1:i))) = Vcell{1,1,i};
   end
end
end