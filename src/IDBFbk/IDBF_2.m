function [Factor] = IDBF_2(fun,xx,kk,r,levels,tol)
% called in IDBBF.m
%
% Copyright 2018 by Qiyuan Pang

if levels > 1

   Factor = struct('U',[],'S',[],'V',[]);
   Factor.U = cell(1,levels);
   Factor.V = cell(1,levels);
   %Factor.S = cell(1);
   NN = max(size(xx));
   MM = max(size(kk));
   N = 0;
   for i = 1:NN
       N = N + size(xx{i},1);
   end
   M = 0;
   for i = 1:MM
       M = M + size(kk{i},1);
   end
   
   sk = cell(2,MM);
   Factor.V{1} = sparse(2*N,M);
   fvr = 0;
   fvc = 0;
   fvcj = 0;
   for j = 1:2
       fvc = fvc + fvcj;
       for i = 1:2
           fvcj = 0;
           for k = 1:MM/2               
               x = [];
               for ii = 1:NN/2
                   x = [x;xx{(i-1)*NN/2+ii}];
               end
               [UR,sk1] = CURR(fun,x,kk{(j-1)*MM/2+k},r,tol);
               Factor.V{1}(fvr+1:fvr+size(UR,1),fvc+fvcj+1:fvc+fvcj+size(UR,2)) = UR;
               fvr = fvr + size(UR,1);
               fvcj = fvcj + size(UR,2);
               sk{i,(j-1)*MM/2+k} = sk1;
           end
       end
   end
   Factor.V{1} = Factor.V{1}(1:fvr,:);
   
   

   sx = cell(2,NN);
   Factor.U{1} = sparse(N,2*M);
   fur = 0;
   fuc = 0;
   furi = 0;
   for i = 1:2
       fur = fur + furi;
       for j = 1:2
           furi = 0;
           for k = 1:NN/2
               idk = [];
               for ii = 1:MM/2
                   idk = [idk;sk{i,(j-1)*MM/2+ii}];
               end
               [CU,sx1] = CURL(fun,xx{(i-1)*NN/2+k},idk,r,tol);
               Factor.U{1}(fur+furi+1:fur+furi+size(CU,1),fuc+1:fuc+size(CU,2)) = CU;
               furi = furi + size(CU,1);
               fuc = fuc + size(CU,2);
               sx{j,(i-1)*NN/2+k} = sx1;
           end
       end
   end
   Factor.U{1} = Factor.U{1}(:,1:fuc);
   
   idx = cell(1,NN/2);
   idk = cell(1,MM/2);
   for i = 1:NN/2
       idx{i} = sx{1,i};
   end
   for j = 1:MM/2
       idk{j} = sk{1,j};
   end
   [Factor1] = IDBF_2(fun,idx,idk,r,levels-1,tol);
      
   idx = cell(1,NN/2);
   idk = cell(1,MM/2);
   for i = 1:NN/2
       idx{i} = sx{2,i};
   end
   for j = 1:MM/2
       idk{j} = sk{1,MM/2+j};
   end
   [Factor2] = IDBF_2(fun,idx,idk,r,levels-1,tol);
   
   idx = cell(1,NN/2);
   idk = cell(1,MM/2);
   for i = 1:NN/2
       idx{i} = sx{1,NN/2+i};
   end
   for j = 1:MM/2
       idk{j} = sk{2,j};
   end
   [Factor3] = IDBF_2(fun,idx,idk,r,levels-1,tol);
   
   
   
   idx = cell(1,NN/2);
   idk = cell(1,MM/2);
   for i = 1:NN/2
       idx{i} = sx{2,NN/2+i};
   end
   for j = 1:MM/2
       idk{j} = sk{2,MM/2+j};
   end
   [Factor4] = IDBF_2(fun,idx,idk,r,levels-1,tol);
   
   for i = 2:levels
       row = [size(Factor1.U{i-1},1) size(Factor2.U{i-1},1) size(Factor3.U{i-1},1) size(Factor4.U{i-1},1)];
       col = [size(Factor1.U{i-1},2) size(Factor2.U{i-1},2) size(Factor3.U{i-1},2) size(Factor4.U{i-1},2)];
       Factor.U{i} = sparse(sum(row),sum(col));
       Factor.U{i}(1:row(1),1:col(1)) = Factor1.U{i-1};
       Factor.U{i}(row(1)+1:sum(row(1:2)),col(1)+1:sum(col(1:2))) = Factor2.U{i-1};
       Factor.U{i}(sum(row(1:2))+1:sum(row(1:3)),sum(col(1:2))+1:sum(col(1:3))) = Factor3.U{i-1};
       Factor.U{i}(sum(row(1:3))+1:sum(row),sum(col(1:3))+1:sum(col)) = Factor4.U{i-1};
       
       row = [size(Factor1.V{i-1},1) size(Factor3.V{i-1},1) size(Factor2.V{i-1},1) size(Factor4.V{i-1},1)];
       col = [size(Factor1.V{i-1},2) size(Factor3.V{i-1},2) size(Factor2.V{i-1},2) size(Factor4.V{i-1},2)];
       Factor.V{i} = sparse(sum(row),sum(col));
       Factor.V{i}(1:row(1),1:col(1)) = Factor1.V{i-1};
       Factor.V{i}(row(1)+1:sum(row(1:2)),col(1)+1:sum(col(1:2))) = Factor3.V{i-1};
       Factor.V{i}(sum(row(1:2))+1:sum(row(1:3)),sum(col(1:2))+1:sum(col(1:3))) = Factor2.V{i-1};
       Factor.V{i}(sum(row(1:3))+1:sum(row),sum(col(1:3))+1:sum(col)) = Factor4.V{i-1};
   end
   
   row = [size(Factor1.S,1) size(Factor2.S,1) size(Factor3.S,1) size(Factor4.S,1)];
   col = [size(Factor1.S,2) size(Factor2.S,2) size(Factor3.S,2) size(Factor4.S,2)];
   Factor.S = sparse(sum(row),sum(col));
   Factor.S(1:row(1),1:col(1)) = Factor1.S;
   Factor.S(row(1)+1:row(1)+row(2),col(1)+col(3)+1:sum(col(1:3))) = Factor2.S;
   Factor.S(row(1)+row(2)+1:sum(row(1:3)),col(1)+1:col(1)+col(3)) = Factor3.S;
   Factor.S(sum(row(1:3))+1:end,sum(col(1:3))+1:end) = Factor4.S;
    
end

if levels == 1
      
   Factor = struct('U',[],'S',[],'V',[]);
   Factor.U = cell(1);
   Factor.V = cell(1);
   %Factor.S = cell(1);
   NN = max(size(xx));
   MM = max(size(kk));
   N = 0;
   for i = 1:NN
       N = N + size(xx{i},1);
   end
   M = 0;
   for i = 1:MM
       M = M + size(kk{i},1);
   end
   
   sk = cell(2,MM);
   Factor.V{1} = sparse(2*N,M);
   fvr = 0;
   fvc = 0;
   fvcj = 0;
   for j = 1:2
       fvc = fvc + fvcj;
       for i = 1:2
           fvcj = 0;
           for k = 1:MM/2               
               x = [];
               for ii = 1:NN/2
                   x = [x;xx{(i-1)*NN/2+ii}];
               end
               [UR,sk1] = CURR(fun,x,kk{(j-1)*MM/2+k},r,tol);
               Factor.V{1}(fvr+1:fvr+size(UR,1),fvc+fvcj+1:fvc+fvcj+size(UR,2)) = UR;
               fvr = fvr + size(UR,1);
               fvcj = fvcj + size(UR,2);
               sk{i}{(j-1)*MM/2+k} = sk1;
           end
       end
   end
   Factor.V{1} = Factor.V{1}(1:fvr,:);
   
   

   sx = cell(2,NN);
   Factor.U{1} = sparse(N,2*M);
   fur = 0;
   fuc = 0;
   furi = 0;
   for i = 1:2
       fur = fur + furi;
       for j = 1:2
           furi = 0;
           for k = 1:NN/2
               idk = [];
               for ii = 1:MM/2
                   idk = [idk;sk{i}{(j-1)*MM/2+ii}];
               end
               [CU,sx1] = CURL(fun,xx{(i-1)*NN/2+k},idk,r,tol);
               Factor.U{1}(fur+furi+1:fur+furi+size(CU,1),fuc+1:fuc+size(CU,2)) = CU;
               furi = furi + size(CU,1);
               fuc = fuc + size(CU,2);
               sx{j}{(i-1)*NN/2+k} = sx1;
           end
       end
   end
   Factor.U{1} = Factor.U{1}(:,1:fuc);
   
   
   R = zeros(4,1);
   C = zeros(4,1);
   for i = 1:2
       for j = 1:NN/2
           R(i) = R(i) + size(sx{i}{j},1);
           R(2+i) = R(i+2) + size(sx{i}{NN/2+j},1);
       end
   end
   for i = 1:2
       for j = 1:MM/2
           C(i) = C(i) + size(sk{i}{j},1);
           C(2+i) = C(i+2) + size(sk{i}{MM/2+j},1);
       end
   end
   Factor.S = zeros(sum(R),sum(C));
   
   sr = 0;
   sv = 0;
   for i = 1:NN/2
       svj = 0;
       for j = 1:MM/2
           Factor.S(sr+1:sr+size(sx{1}{i},1),sv+svj+1:sv+svj+size(sk{1}{j},1)) = fun(sx{1}{i},sk{1}{j});
           svj = svj + size(sk{1}{j},1);
       end
       sr = sr + size(sx{1}{i},1);
   end

   sr = R(1);
   sv = C(1)+C(2);
   for i = 1:NN/2
       svj = 0;
       for j = 1:MM/2
           Factor.S(sr+1:sr+size(sx{2}{i},1),sv+svj+1:sv+svj+size(sk{1}{MM/2+j},1)) = fun(sx{2}{i},sk{1}{MM/2+j});
           svj = svj + size(sk{1}{MM/2+j},1);
       end
       sr = sr + size(sx{2}{i},1);
   end

   sr = R(1)+R(2);
   sv = C(1);
   for i = 1:NN/2
       svj = 0;
       for j = 1:MM/2
           Factor.S(sr+1:sr+size(sx{1}{NN/2+i},1),sv+svj+1:sv+svj+size(sk{2}{j},1)) = fun(sx{1}{NN/2+i},sk{2}{j});
           svj = svj + size(sk{2}{j},1);
       end
       sr = sr + size(sx{1}{NN/2+i},1);
   end

   sr = R(1)+R(2)+R(3);
   sv = C(1)+C(2)+C(3);
   for i = 1:NN/2
       svj = 0;
       for j = 1:MM/2
           Factor.S(sr+1:sr+size(sx{2}{NN/2+i},1),sv+svj+1:sv+svj+size(sk{2}{MM/2+j},1)) = fun(sx{2}{NN/2+i},sk{2}{MM/2+j});
           svj = svj + size(sk{2}{MM/2+j},1);
       end
       sr = sr + size(sx{2}{NN/2+i},1);
   end
end

end