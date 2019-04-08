function [Factor] = IDBF(fun_org,xx,kk,r,tol,tag,srand)

if size(xx,1)<size(xx,2), error('xx should be a tall skinny matrix'); end
if size(kk,1)<size(kk,2), error('kk should be a tall skinny matrix'); end
if nargin == 5
    tag = 'Regular';
    srand = 1;
end
if nargin == 6
    srand = 1;
end

[fun,xx,kk,xbox,kbox,npx,npk] = BF_prepbox(fun_org,xx,kk,tag);
% make sure that the smallest dimension of the low-rank submatrices is
% larger than max(8,r)
levels = max(0,min(floor(log2([npx npk])-max(3,ceil(log2(r))))))+1;
Dim = size(xx,2);

xxboxidx = zeros(size(xx));
npxx = npx*2^levels;
xxset = cell(1,prod(npxx));
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
    xxset{itx} = xx(xxidx(xxIA(itx):xxIA(itx+1)-1),:);
    IIdx = xxidx(xxIA(itx):xxIA(itx+1)-1);
    for i = 1:length(IIdx)
        iidx(ax + i) = IIdx(i);
    end
    ax = ax + length(IIdx);
end

kkboxidx = zeros(size(kk));
npkk = npk*2^levels;
kkset = cell(prod(npkk),1);
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
    kkset{itk} = kk(kkidx(kkIA(itk):kkIA(itk+1)-1),:);
    IIdk = kkidx(kkIA(itk):kkIA(itk+1)-1);
    for i = 1:length(IIdk)
        iidk(ak + i) = IIdk(i);
    end
    ak = ak + length(IIdk);
end

k1 = [];
k2 = [];
for i = 1:prod(npkk)/2
    k1 = [k1;kkset{i}];
    k2 = [k2;kkset{prod(npkk)/2+i}];
end

Factor = struct('U',[],'S',[],'V',[]);
Factor.U = cell(1,levels);
Factor.V = cell(1,levels);

sx = cell(2*prod(npxx),1);
sc = 0;
fsr = 0;
XT = [];
YT = [];
ST = [];
for i = 1:2
    for j = 1:2
        sr = 0;
        for k = 1:prod(npxx)/2
            x = xxset{(i-1)*prod(npxx)/2+k};
            if j == 1
               kj = k1;
            else
               kj = k2;
            end
            [CU,sx1] = CURL(fun,x,kj,r,tol);
            dk = randsample(size(kj,1),min(r,size(kj,1)));
            kk = kj(dk,:);
            err = norm(fun(x,kk)-CU*fun(sx1,kk))/norm(fun(x,kk));
            if err > 10*tol
               [sx1,CU] = IDL(fun,x,kj,tol,srand); 
            end
            r1 = size(CU,1);
            c1 = size(CU,2);
            S = CU(:);
            ST = [ST;S];
            X = iidx(fsr+sr+1:fsr+sr+r1);
            Y = ones(r1,1);
            for jj = 1:c1
                XT = [XT;X];
                YT = [YT;(sc+jj)*Y];
            end
            sx{(i-1)*prod(npxx)+(j-1)*prod(npxx)/2+k} = sx1;
            sc = sc + c1;
            sr = sr + r1;
        end
    end
    fsr = fsr + sr;
end
Factor.U{1} = sparse(XT,YT,ST,fsr,sc);

x1 = [];
x2 = [];
x3 = [];
x4 = [];
for i = 1:prod(npxx)/2
    x1 = [x1;sx{i}];
    x2 = [x2;sx{prod(npxx)/2+i}];
    x3 = [x3;sx{prod(npxx)+i}];
    x4 = [x4;sx{3*prod(npxx)/2+i}];
end

sk = cell(2*prod(npkk),1);
sr = 0;
fsc = 0;
XT = [];
YT = [];
ST = [];
for j = 1:2
    for i = 1:2
        sc = 0;
        for ii = 1:prod(npkk)/2
            k = kkset{(j-1)*prod(npkk)/2+ii};
            if j == 1 && i == 1
               xij = x1;
            end
            if j == 1 && i == 2
                xij = x3;
            end
            if j == 2 && i == 1
                xij = x2;
            end
            if j == 2 && i == 2
                xij = x4;
            end
            [UR,sk1] = CURR(fun,xij,k,r,tol);
            dx = randsample(size(xij,1),min(r,size(xij,1)));
            x = xij(dx,:);
            err = norm(fun(x,k)-fun(x,sk1)*UR)/norm(fun(x,k));
            if err > 10*tol
               [sk1,UR] = IDR(fun,xij,k,tol,srand); 
            end
            r1 = size(UR,1);
            c1 = size(UR,2);
            S = UR(:);
            ST = [ST;S];
            X = [sr+1:sr+r1].';
            Y = ones(r1,1);
            for jj = 1:c1
                XT = [XT;X];
                YT = [YT;iidk(fsc+sc+jj)*Y];
            end
            sk{(j-1)*prod(npkk)+(i-1)*prod(npkk)/2+ii} = sk1;
            sr = sr + r1;
            sc = sc + c1;
        end
    end
    fsc = fsc + sc;
end
Factor.V{1} = sparse(XT,YT,ST,sr,fsc);

%levels = min(log2(prod(npxx))-1,log2(prod(npkk))-1);

if levels > 1
Factor1 = IDBFpre(fun,sx(1:prod(npxx)/2),sk(1:prod(npkk)/2),r,tol,levels-1,srand);
Factor2 = IDBFpre(fun,sx(prod(npxx)/2+1:prod(npxx)),sk(prod(npkk)+1:3*prod(npkk)/2),r,tol,levels-1,srand);
Factor3 = IDBFpre(fun,sx(prod(npxx)+1:3*prod(npxx)/2),sk(prod(npkk)/2+1:prod(npkk)),r,tol,levels-1,srand);
Factor4 = IDBFpre(fun,sx(3*prod(npxx)/2+1:2*prod(npxx)),sk(3*prod(npkk)/2+1:2*prod(npkk)),r,tol,levels-1,srand);
end

for lvl = 2:levels
    totalR = [Factor1.U{lvl-1}.totalR Factor2.U{lvl-1}.totalR Factor3.U{lvl-1}.totalR Factor4.U{lvl-1}.totalR];
    totalC = [Factor1.U{lvl-1}.totalC Factor2.U{lvl-1}.totalC Factor3.U{lvl-1}.totalC Factor4.U{lvl-1}.totalC];
    XT = [Factor1.U{lvl-1}.XT;totalR(1)+Factor2.U{lvl-1}.XT;...
        totalR(1)+totalR(2)+Factor3.U{lvl-1}.XT;sum(totalR(1:3))+Factor4.U{lvl-1}.XT];
    YT = [Factor1.U{lvl-1}.YT;totalC(1)+Factor2.U{lvl-1}.YT;...
        totalC(1)+totalC(2)+Factor3.U{lvl-1}.YT;sum(totalC(1:3))+Factor4.U{lvl-1}.YT];
    ST = [Factor1.U{lvl-1}.ST;Factor2.U{lvl-1}.ST;Factor3.U{lvl-1}.ST;Factor4.U{lvl-1}.ST];
    Factor.U{lvl} = sparse(XT,YT,ST,sum(totalR),sum(totalC));
    
    totalR = [Factor1.V{lvl-1}.totalR Factor3.V{lvl-1}.totalR Factor2.V{lvl-1}.totalR Factor4.V{lvl-1}.totalR];
    totalC = [Factor1.V{lvl-1}.totalC Factor3.V{lvl-1}.totalC Factor2.V{lvl-1}.totalC Factor4.V{lvl-1}.totalC];
    XT = [Factor1.V{lvl-1}.XT;totalR(1)+Factor3.V{lvl-1}.XT;...
        totalR(1)+totalR(2)+Factor2.V{lvl-1}.XT;sum(totalR(1:3))+Factor4.V{lvl-1}.XT];
    YT = [Factor1.V{lvl-1}.YT;totalC(1)+Factor3.V{lvl-1}.YT;...
        totalC(1)+totalC(2)+Factor2.V{lvl-1}.YT;sum(totalC(1:3))+Factor4.V{lvl-1}.YT];
    ST = [Factor1.V{lvl-1}.ST;Factor3.V{lvl-1}.ST;Factor2.V{lvl-1}.ST;Factor4.V{lvl-1}.ST];
    Factor.V{lvl} = sparse(XT,YT,ST,sum(totalR),sum(totalC));
end

if levels > 1
totalR = [Factor1.S.totalR Factor2.S.totalR Factor3.S.totalR Factor4.S.totalR];
totalC = [Factor1.S.totalC Factor3.S.totalC Factor2.S.totalC Factor4.S.totalC];
XT = [Factor1.S.XT;totalR(1)+Factor2.S.XT;...
        totalR(1)+totalR(2)+Factor3.S.XT;sum(totalR(1:3))+Factor4.S.XT];
YT = [Factor1.S.YT;totalC(1)+totalC(2)+Factor2.S.YT;...
        totalC(1)+Factor3.S.YT;sum(totalC(1:3))+Factor4.S.YT];
ST = [Factor1.S.ST;Factor2.S.ST;Factor3.S.ST;Factor4.S.ST];
Factor.S = sparse(XT,YT,ST,sum(totalR),sum(totalC));
end

if levels == 1
XT = [];
YT = [];
ST = [];
sr = 0;
npxx = prod(npxx);
npkk = prod(npkk);
for i = 1:npxx/2
    sc = 0;
    for j = 1:npkk/2
        C = fun(sx{i},sk{j});
        r1 = size(C,1);
        c1 = size(C,2);
        S = C(:);
        ST = [ST;S];
        X = [sr+1:sr+r1]';
        Y = ones(r1,1);
        for k = 1:c1
            XT = [XT;X];
            YT = [YT;(sc+k)*Y];
        end
        sc = sc + c1;
    end
    sr = sr + size(sx{i},1);
end
totalR1 = sr;
totalC1 = sc;

YT1 = [];
for i = 1:npxx/2
    sc = 0;
    for j = 1:npkk/2
        C = fun(sx{npxx/2+i},sk{npkk+j});
        r1 = size(C,1);
        c1 = size(C,2);
        S = C(:);
        ST = [ST;S];
        X = [sr+1:sr+r1]';
        Y = ones(r1,1);
        for k = 1:c1
            XT = [XT;X];
            YT1 = [YT1;(sc+k)*Y];
        end
        sc = sc + c1;
    end
    sr = sr + size(sx{npxx/2+i},1);
end
totalR2 = sr - totalR1;
totalC2 = sc;

for i = 1:npxx/2
    sc = 0;
    for j = 1:npkk/2
        C = fun(sx{npxx+i},sk{npkk/2+j});
        r1 = size(C,1);
        c1 = size(C,2);
        S = C(:);
        ST = [ST;S];
        X = [sr+1:sr+r1]';
        Y = ones(r1,1);
        for k = 1:c1
            XT = [XT;X];
            YT1 = [YT1;(sc+k)*Y+totalC1];
        end
        sc = sc + c1;
    end
    sr = sr + size(sx{npxx+i},1);
end
totalR3 = sr - totalR1 - totalR2;
totalC3 = sc;

YT1(1:totalR2*totalC2) = YT1(1:totalR2*totalC2) + totalC1 + totalC3;
YT = [YT;YT1];

for i = 1:npxx/2
    sc = 0;
    for j = 1:npkk/2
        C = fun(sx{3*npxx/2+i},sk{3*npkk/2+j});
        r1 = size(C,1);
        c1 = size(C,2);
        S = C(:);
        ST = [ST;S];
        X = [sr+1:sr+r1]';
        Y = ones(r1,1);
        for k = 1:c1
            XT = [XT;X];
            YT = [YT;(sc+k)*Y+totalC1+totalC2+totalC3];
        end
        sc = sc + c1;
    end
    sr = sr + size(sx{3*npxx/2+i},1);
end
totalC4 = sc;
Factor.S = sparse(XT,YT,ST,sr,totalC1 + totalC2 + totalC3 + totalC4);
end
end