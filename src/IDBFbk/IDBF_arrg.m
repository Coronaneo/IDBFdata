function [Factor] = IDBF_arrg(Faccell,Ucell,Vcell,flag)
% called in IDBBF.m
%
% Copyright 2018 by Qiyuan Pang
if flag > 0
   prodnpx = size(Faccell,1);
   prodnpk = size(Faccell,2);
   levels = size(Faccell{1,1}.U,2)+1;
   Factor = struct('U',[],'S',[],'V',[]);
   Factor.U = cell(1,levels);
   Factor.V = cell(1,levels);
   [Factor.S,rv,cv] = IDBF_iterS(Faccell,flag);
   Factor.U{1} = IDBF_iterU1(Ucell,fftshift(rv));
   Factor.V{1} = IDBF_iterV1(Vcell,cv);
      
   for k = 2:levels
       row = [];
       col = [];
       for i = 1:prodnpx*prodnpk
           row = [row size(Faccell{rv(i,1),rv(i,2)}.U{k-1},1)];
           col = [col size(Faccell{rv(i,1),rv(i,2)}.U{k-1},2)];
       end
       Factor.U{k} = sparse(sum(row),sum(col));
       Factor.U{k}(1:row(1),1:col(1)) = Faccell{rv(1,1),rv(1,2)}.U{k-1};
       for i = 2:prodnpx*prodnpk
           Factor.U{k}(sum(row(1:i-1))+1:sum(row(1:i)),sum(col(1:i-1))+1:sum(col(1:i))) = Faccell{rv(i,1),rv(i,2)}.U{k-1};
           
       end
       
       row = [];
       col = [];
       for i = 1:prodnpx*prodnpk
           row = [row size(Faccell{cv(i,1),cv(i,2)}.V{k-1},1)];
           col = [col size(Faccell{cv(i,1),cv(i,2)}.V{k-1},2)];
       end
       Factor.V{k} = sparse(sum(row),sum(col));
       Factor.V{k}(1:row(1),1:col(1)) = Faccell{cv(1,1),cv(1,2)}.V{k-1};
       for i = 2:prodnpx*prodnpk
           Factor.V{k}(sum(row(1:i-1))+1:sum(row(1:i)),sum(col(1:i-1))+1:sum(col(1:i))) = Faccell{cv(i,1),cv(i,2)}.V{k-1};
       end
   end
else
   Factor = struct('U',[],'S',[],'V',[]);
   Factor.U = cell(1);
   Factor.V = cell(1);
   [Factor.S,rv,cv] = IDBF_iterS(Faccell,flag);
   Factor.U{1} = IDBF_iterU1(Ucell,fftshift(rv));
   Factor.V{1} = IDBF_iterV1(Vcell,cv);
end