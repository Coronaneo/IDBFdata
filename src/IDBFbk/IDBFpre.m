function [Factor] = IDBFpre(fun,xx,kk,r,tol,levels,srand)
if  levels > 1
Factor = struct('U',[],'S',[],'V',[]);
Factor.U = cell(1,levels);
Factor.V = cell(1,levels);

npxx = size(xx,1);
npkk = size(kk,1);

k1 = [];
k2 = [];
for i = 1:npkk/2
    k1 = [k1;kk{i}];
    k2 = [k2;kk{npkk/2+i}];
end

sx = cell(2*npxx,1);
sc = 0;
fsr = 0;
XT = [];
YT = [];
ST = [];
for i = 1:2
    for j = 1:2
        sr = 0;
        for k = 1:npxx/2
            x = xx{(i-1)*npxx/2+k};
            if j == 1
               kj = k1;
            else
               kj = k2;
            end
            [CU,sx1] = CURL(fun,x,kj,r,tol);
            dk = randsample(size(kj,1),min(r,size(kj,1)));
            kkk = kj(dk,:);
            err = norm(fun(x,kkk)-CU*fun(sx1,kkk))/norm(fun(x,kkk));
            if err > 10*tol
               [sx1,CU] = IDL(fun,x,kj,tol,srand); 
            end
            r1 = size(CU,1);
            c1 = size(CU,2);
            S = CU(:);
            ST = [ST;S];
            X = [fsr+sr+1:fsr+sr+r1]';
            Y = ones(r1,1);
            for jj = 1:c1
                XT = [XT;X];
                YT = [YT;(sc+jj)*Y];
            end
            sx{(i-1)*npxx+(j-1)*npxx/2+k} = sx1;
            sc = sc + c1;
            sr = sr + r1;
        end
    end
    fsr = fsr + sr;
end
Factor.U{1} = struct('XT',[],'YT',[],'ST',[],'totalR',[],'totalC',[]);
Factor.U{1}.XT = XT;
Factor.U{1}.YT = YT;
Factor.U{1}.ST = ST;
Factor.U{1}.totalR = fsr;
Factor.U{1}.totalC = sc;

x1 = [];
x2 = [];
x3 = [];
x4 = [];
for i = 1:npxx/2
    x1 = [x1;sx{i}];
    x2 = [x2;sx{npxx/2+i}];
    x3 = [x3;sx{npxx+i}];
    x4 = [x4;sx{3*npxx/2+i}];
end

sk = cell(2*npkk,1);
sr = 0;
fsc = 0;
XT = [];
YT = [];
ST = [];
for j = 1:2
    for i = 1:2
        sc = 0;
        for ii = 1:npkk/2
            k = kk{(j-1)*npkk/2+ii};
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
                YT = [YT;(fsc+sc+jj)*Y];
            end
            sk{(j-1)*npkk+(i-1)*npkk/2+ii} = sk1;
            sr = sr + r1;
            sc = sc + c1;
        end
    end
    fsc = fsc + sc;
end
Factor.V{1} = struct('XT',[],'YT',[],'ST',[],'totalR',[],'totalC',[]);
Factor.V{1}.XT = XT;
Factor.V{1}.YT = YT;
Factor.V{1}.ST = ST;
Factor.V{1}.totalR = sr;
Factor.V{1}.totalC = fsc;

Factor1 = IDBFpre(fun,sx(1:npxx/2),sk(1:npkk/2),r,tol,levels-1,srand);
Factor2 = IDBFpre(fun,sx(npxx/2+1:npxx),sk(npkk+1:3*npkk/2),r,tol,levels-1,srand);
Factor3 = IDBFpre(fun,sx(npxx+1:3*npxx/2),sk(npkk/2+1:npkk),r,tol,levels-1,srand);
Factor4 = IDBFpre(fun,sx(3*npxx/2+1:2*npxx),sk(3*npkk/2+1:2*npkk),r,tol,levels-1,srand);

for lvl = 2:levels
    totalR = [Factor1.U{lvl-1}.totalR Factor2.U{lvl-1}.totalR Factor3.U{lvl-1}.totalR Factor4.U{lvl-1}.totalR];
    totalC = [Factor1.U{lvl-1}.totalC Factor2.U{lvl-1}.totalC Factor3.U{lvl-1}.totalC Factor4.U{lvl-1}.totalC];
    XT = [Factor1.U{lvl-1}.XT;totalR(1)+Factor2.U{lvl-1}.XT;...
        totalR(1)+totalR(2)+Factor3.U{lvl-1}.XT;sum(totalR(1:3))+Factor4.U{lvl-1}.XT];
    YT = [Factor1.U{lvl-1}.YT;totalC(1)+Factor2.U{lvl-1}.YT;...
        totalC(1)+totalC(2)+Factor3.U{lvl-1}.YT;sum(totalC(1:3))+Factor4.U{lvl-1}.YT];
    ST = [Factor1.U{lvl-1}.ST;Factor2.U{lvl-1}.ST;Factor3.U{lvl-1}.ST;Factor4.U{lvl-1}.ST];
    Factor.U{lvl} = struct('XT',[],'YT',[],'ST',[],'totalR',[],'totalC',[]);
    Factor.U{lvl}.XT = XT;
    Factor.U{lvl}.YT = YT;
    Factor.U{lvl}.ST = ST;
    Factor.U{lvl}.totalR = sum(totalR);
    Factor.U{lvl}.totalC = sum(totalC);
    
    totalR = [Factor1.V{lvl-1}.totalR Factor3.V{lvl-1}.totalR Factor2.V{lvl-1}.totalR Factor4.V{lvl-1}.totalR];
    totalC = [Factor1.V{lvl-1}.totalC Factor3.V{lvl-1}.totalC Factor2.V{lvl-1}.totalC Factor4.V{lvl-1}.totalC];
    XT = [Factor1.V{lvl-1}.XT;totalR(1)+Factor3.V{lvl-1}.XT;...
        totalR(1)+totalR(2)+Factor2.V{lvl-1}.XT;sum(totalR(1:3))+Factor4.V{lvl-1}.XT];
    YT = [Factor1.V{lvl-1}.YT;totalC(1)+Factor3.V{lvl-1}.YT;...
        totalC(1)+totalC(2)+Factor2.V{lvl-1}.YT;sum(totalC(1:3))+Factor4.V{lvl-1}.YT];
    ST = [Factor1.V{lvl-1}.ST;Factor3.V{lvl-1}.ST;Factor2.V{lvl-1}.ST;Factor4.V{lvl-1}.ST];
    Factor.V{lvl} = struct('XT',[],'YT',[],'ST',[],'totalR',[],'totalC',[]);
    Factor.V{lvl}.XT = XT;
    Factor.V{lvl}.YT = YT;
    Factor.V{lvl}.ST = ST;
    Factor.V{lvl}.totalR = sum(totalR);
    Factor.V{lvl}.totalC = sum(totalC);
end

totalR = [Factor1.S.totalR Factor2.S.totalR Factor3.S.totalR Factor4.S.totalR];
totalC = [Factor1.S.totalC Factor3.S.totalC Factor2.S.totalC Factor4.S.totalC];
XT = [Factor1.S.XT;totalR(1)+Factor2.S.XT;...
        totalR(1)+totalR(2)+Factor3.S.XT;sum(totalR(1:3))+Factor4.S.XT];
YT = [Factor1.S.YT;totalC(1)+totalC(2)+Factor2.S.YT;...
        totalC(1)+Factor3.S.YT;sum(totalC(1:3))+Factor4.S.YT];
ST = [Factor1.S.ST;Factor2.S.ST;Factor3.S.ST;Factor4.S.ST];
Factor.S = struct('XT',[],'YT',[],'ST',[],'totalR',[],'totalC',[]);
Factor.S.XT = XT;
Factor.S.YT = YT;
Factor.S.ST = ST;
Factor.S.totalR = sum(totalR);
Factor.S.totalC = sum(totalC);

end

if  levels == 1
Factor = struct('U',[],'S',[],'V',[]);
Factor.U = cell(1,levels);
Factor.V = cell(1,levels);

npxx = size(xx,1);
npkk = size(kk,1);

k1 = [];
k2 = [];
for i = 1:npkk/2
    k1 = [k1;kk{i}];
    k2 = [k2;kk{npkk/2+i}];
end

sx = cell(2*npxx,1);
sc = 0;
fsr = 0;
XT = [];
YT = [];
ST = [];
for i = 1:2
    for j = 1:2
        sr = 0;
        for k = 1:npxx/2
            x = xx{(i-1)*npxx/2+k};
            if j == 1
               kj = k1;
            else
               kj = k2;
            end
            [CU,sx1] = CURL(fun,x,kj,r,tol);
            dk = randsample(size(kj,1),min(r,size(kj,1)));
            kkk = kj(dk,:);
            err = norm(fun(x,kkk)-CU*fun(sx1,kkk))/norm(fun(x,kkk));
            if err > 10*tol
               [sx1,CU] = IDL(fun,x,kj,tol,srand); 
            end
            r1 = size(CU,1);
            c1 = size(CU,2);
            S = CU(:);
            ST = [ST;S];
            X = [fsr+sr+1:fsr+sr+r1]';
            Y = ones(r1,1);
            for jj = 1:c1
                XT = [XT;X];
                YT = [YT;(sc+jj)*Y];
            end
            sx{(i-1)*npxx+(j-1)*npxx/2+k} = sx1;
            sc = sc + c1;
            sr = sr + r1;
        end
    end
    fsr = fsr + sr;
end
Factor.U{1} = struct('XT',[],'YT',[],'ST',[],'totalR',[],'totalC',[]);
Factor.U{1}.XT = XT;
Factor.U{1}.YT = YT;
Factor.U{1}.ST = ST;
Factor.U{1}.totalR = fsr;
Factor.U{1}.totalC = sc;

x1 = [];
x2 = [];
x3 = [];
x4 = [];
for i = 1:npxx/2
    x1 = [x1;sx{i}];
    x2 = [x2;sx{npxx/2+i}];
    x3 = [x3;sx{npxx+i}];
    x4 = [x4;sx{3*npxx/2+i}];
end

sk = cell(2*npkk,1);
sr = 0;
fsc = 0;
XT = [];
YT = [];
ST = [];
for j = 1:2
    for i = 1:2
        sc = 0;
        for ii = 1:npkk/2
            k = kk{(j-1)*npkk/2+ii};
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
                YT = [YT;(fsc+sc+jj)*Y];
            end
            sk{(j-1)*npkk+(i-1)*npkk/2+ii} = sk1;
            sr = sr + r1;
            sc = sc + c1;
        end
    end
    fsc = fsc + sc;
end
Factor.V{1} = struct('XT',[],'YT',[],'ST',[],'totalR',[],'totalC',[]);
Factor.V{1}.XT = XT;
Factor.V{1}.YT = YT;
Factor.V{1}.ST = ST;
Factor.V{1}.totalR = sr;
Factor.V{1}.totalC = fsc;

XT = [];
YT = [];
ST = [];
sr = 0;
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

Factor.S = struct('XT',[],'YT',[],'ST',[],'totalR',[],'totalC',[]);
Factor.S.XT = XT;
Factor.S.YT = YT;
Factor.S.ST = ST;
Factor.S.totalR = sr;
Factor.S.totalC = totalC1+ totalC2+ totalC3+ totalC4;

end