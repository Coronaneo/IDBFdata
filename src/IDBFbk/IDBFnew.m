function [Factor] = IDBFnew(fun_org,xx,kk,r,tol,tag)
% This is the main file to compute the ID butterfly factorization. It
% compress the butterfly matrix from both leaf levels towards the middle
% level.
%
% Copyright 2018 by Qiyuan Pang

if size(xx,1)<size(xx,2), error('xx should be a tall skinny matrix'); end
if size(kk,1)<size(kk,2), error('kk should be a tall skinny matrix'); end
if nargin == 5
    tag = 'Regular';
end

[fun,xx,kk,xbox,kbox,npx,npk] = BF_prepbox(fun_org,xx,kk,tag);
% make sure that the smallest dimension of the low-rank submatrices is
% larger than max(8,r)
levels = max(0,min(floor(log2([npx npk])-max(3,ceil(log2(r))))))+1;
Dim = size(xx,2);

xxset = cell(levels+1,1);

xxboxidx = zeros(size(xx));
npxx = npx*2^levels;
xxset{levels+1} = cell(prod(npxx),1);
for i = 1:Dim
    edges = linspace(xbox(1,i),xbox(2,i),npxx(i)+1);
    [~,xxboxidx(:,i)] = histc(xx(:,i),edges);
end
[xxboxidx,xxidx] = sortrows(xxboxidx,Dim:-1:1);
[xxC,xxIA,~] = unique(xxboxidx,'rows','stable');
xxC = [xxC;zeros(1,size(xxC,2))];
xxIA = [xxIA;size(xx,1)+1];
for itx = 1:prod(npxx)
    x = BF_idx2vec(npxx,itx);
    if any(x ~= xxC(itx,:))
        xxC = [xxC(1:itx-1,:);x;xxC(itx:end,:)];
        xxIA = [xxIA(1:itx-1);xxIA(itx);xxIA(itx:end)];
    end
end
nx = size(xx,1);
iidx = zeros(nx,1);
ax = 0;
for itx = 1:prod(npxx)
    xxset{levels+1}{itx} = xx(xxidx(xxIA(itx):xxIA(itx+1)-1),:);
    IIdx = xxidx(xxIA(itx):xxIA(itx+1)-1);
    for i = 1:length(IIdx)
        iidx(IIdx(i)) = ax + i;
    end
    ax = ax + length(IIdx);
end
for lvl = levels-1:-1:0
    npxx = npx*2^lvl;
    xxset{lvl+1} = cell(prod(npxx),1);
    for itx = 1:prod(npxx)
        xxset{lvl+1}{itx} = [];
        for itc = BF_childidx(npxx,itx)
            xxset{lvl+1}{itx} = [ xxset{lvl+1}{itx}; xxset{lvl+2}{itc}];
        end
    end
end

kkset = cell(levels+1,1);

kkboxidx = zeros(size(kk));
npkk = npk*2^levels;
kkset{levels+1} = cell(prod(npkk),1);
for i = 1:Dim
    edges = linspace(kbox(1,i),kbox(2,i),npkk(i)+1);
    [~,kkboxidx(:,i)] = histc(kk(:,i),edges);
end
[kkboxidx,kkidx] = sortrows(kkboxidx,Dim:-1:1);
[kkC,kkIA,~] = unique(kkboxidx,'rows','stable');
kkC = [kkC;zeros(1,size(kkC,2))];
kkIA = [kkIA;size(kk,1)+1];
for itk = 1:prod(npkk)
    k = BF_idx2vec(npkk,itk);
    if any(k ~= kkC(itk,:))
        kkC = [kkC(1:itk-1,:);k;kkC(itk:end,:)];
        kkIA = [kkIA(1:itk-1);kkIA(itk);kkIA(itk:end)];
    end
end
nk = size(kk,1);
iidk = zeros(nk,1);
ak = 0;
for itk = 1:prod(npkk)
    kkset{levels+1}{itk} = kk(kkidx(kkIA(itk):kkIA(itk+1)-1),:);
    IIdk = kkidx(kkIA(itk):kkIA(itk+1)-1);
    for i = 1:length(IIdk)
        iidk(IIdk(i)) = ak + i;
    end
    ak = ak + length(IIdk);
end

for lvl = levels-1:-1:0
    npkk = npk*2^lvl;
    kkset{lvl+1} = cell(prod(npkk),1);
    for itk = 1:prod(npkk)
        kkset{lvl+1}{itk} = [];
        for itc = BF_childidx(npkk,itk)
            kkset{lvl+1}{itk} = [ kkset{lvl+1}{itk}; kkset{lvl+2}{itc}];
        end
    end
end

npxx = npx*2^levels;
npkk = npk*2^levels;

flag = prod(npx)/prod(npk);

   if flag >= 1
      Faccell = cell(flag);
      for i = 1:flag
           Faccell{i} = IDBF_2(fun,xxset{levels+1}((i-1)*prod(npxx)/flag+1:i*prod(npxx)/flag),...
               kkset{levels+1}(:),r,levels,tol);
      end
   else 
      flag = prod(npk)/prod(npx);
      Faccell = cell(flag);
      for i = 1:flag
           Faccell{i} = IDBF_2(fun,xxset{levels+1}(:),...
               kkset{levels+1}((i-1)*prod(npkk)/flag+1:i*prod(npkk)/flag),r,levels,tol);
      end
   end
flag = prod(npx)/prod(npk);

%flag = prod(npx)/prod(npk);
if flag == 1
   Factor = Faccell{1};   
end
if flag > 1
   Factor = struct('U',[],'S',[],'V',[]);
   Factor.U = cell(1,levels);
   Factor.V = cell(1,levels);
   for i = 1:levels
       row = [];
       col = [];
       for j = 1:flag
           row = [row size(Faccell{j}.U{i},1)];
           col = [col size(Faccell{j}.U{i},2)];
       end
       Factor.U{i} = sparse(sum(row),sum(col));
       Factor.U{i}(1:row(1),1:col(1)) = Faccell{1}.U{i};
       for j = 2:flag
           Factor.U{i}(sum(row(1:j-1))+1:sum(row(1:j)),sum(col(1:j-1))+1:sum(col(1:j))) = Faccell{j}.U{i};
       end
   end
   row = [];
   col = [];
   for j = 1:flag
       row = [row size(Faccell{j}.S,1)];
       col = [col size(Faccell{j}.S,2)];
   end
   Factor.S = sparse(sum(row),sum(col));
   Factor.S(1:row(1),1:col(1)) = Faccell{1}.S;
   for j = 2:flag
       Factor.S(sum(row(1:j-1))+1:sum(row(1:j)),sum(col(1:j-1))+1:sum(col(1:j))) = Faccell{j}.S;
   end
   for i = 2:levels
       row = [];
       col = [];
       for j = 1:flag
           row = [row size(Faccell{j}.V{i},1)];
           col = [col size(Faccell{j}.V{i},2)];
       end
       Factor.V{i} = sparse(sum(row),sum(col));
       Factor.V{i}(1:row(1),1:col(1)) = Faccell{1}.V{i};
       for j = 2:flag
           Factor.V{i}(sum(row(1:j-1))+1:sum(row(1:j)),sum(col(1:j-1))+1:sum(col(1:j))) = Faccell{j}.V{i};
       end
   end
   row = [];
   col = [];
   for j = 1:flag
       row = [row size(Faccell{j}.V{1},1)];
       col = [col size(Faccell{j}.V{1},2)];
   end
   Factor.V{1} = sparse(sum(row),max(col));
   Factor.V{1}(1:row(1),1:col(1)) = Faccell{1}.V{1};
   for j = 2:flag
       Factor.V{1}(sum(row(1:j-1))+1:sum(row(1:j)),1:col(j)) = Faccell{j}.V{1};
   end
end

if flag < 1
   flag = prod(npk)/prod(npx);
   Factor = struct('U',[],'S',[],'V',[]);
   Factor.U = cell(1,levels);
   Factor.V = cell(1,levels);
   for i = 1:levels
       row = [];
       col = [];
       for j = 1:flag
           row = [row size(Faccell{j}.V{i},1)];
           col = [col size(Faccell{j}.V{i},2)];
       end
       Factor.V{i} = sparse(sum(row),sum(col));
       Factor.V{i}(1:row(1),1:col(1)) = FF{1}.V{i};
       for j = 2:flag
           Factor.V{i}(sum(row(1:j-1))+1:sum(row(1:j)),sum(col(1:j-1))+1:sum(col(1:j))) = Faccell{j}.V{i};
       end
   end
   row = [];
   col = [];
   for j = 1:flag
       row = [row size(Faccell{j}.S,1)];
       col = [col size(Faccell{j}.S,2)];
   end
   Factor.S = sparse(sum(row),sum(col));
   Factor.S(1:row(1),1:col(1)) = FF{1}.S;
   for j = 2:flag
       Factor.S(sum(row(1:j-1))+1:sum(row(1:j)),sum(col(1:j-1))+1:sum(col(1:j))) = Faccell{j}.S;
   end
   for i = 2:levels
       row = [];
       col = [];
       for j = 1:flag
           row = [row size(Faccell{j}.U{i},1)];
           col = [col size(Faccell{j}.U{i},2)];
       end
       Factor.U{i} = sparse(sum(row),sum(col));
       Factor.U{i}(1:row(1),1:col(1)) = FF{1}.U{i};
       for j = 2:flag
           Factor.U{i}(sum(row(1:j-1))+1:sum(row(1:j)),sum(col(1:j-1))+1:sum(col(1:j))) = Faccell{j}.U{i};
       end
   end
   row = [];
   col = [];
   for j = 1:flag
       row = [row size(Faccell{j}.U{1},1)];
       col = [col size(Faccell{j}.U{1},2)];
   end
   Factor.U{1} = sparse(max(row),sum(col));
   Factor.U{1}(1:row(1),1:col(1)) = FF{1}.U{1};
   for j = 2:flag
       Factor.U{1}(1:row(j),sum(col(1:j-1))+1:sum(col(1:j))) = Faccell{j}.U{1};
   end
end

Factor.U{1} = Factor.U{1}(iidx,:);
Factor.V{1} = Factor.V{1}(:,iidk);
end