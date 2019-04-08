function x = BF_inv_spemd(Blu,szZ,b,B,opt)
% if nargin > 3
% [Blu{1},Blu{2},Blu{3},Blu{4},Blu{5}] = lu(B([N+1:end,1:N],[N+1:end,1:N]));
% else
% [Blu{1},Blu{2},Blu{3},Blu{4},Blu{5}] = lu(B);
% end
if nargin < 5, opt = 0; end
if nargin < 4
    if iscell(Blu)
        switch length(Blu)
            case 5
                % [Blu{1},Blu{2},Blu{3},Blu{4},Blu{5}] = lu(Z);
                c = [b;zeros(size(Blu{2},2)-szZ,size(b,2))];
                z = Blu{2}\(Blu{1}\(Blu{3}*((Blu{5}\c))));
                y = Blu{4}*z;
            case 4
                % [Blu{1},Blu{2},Blu{3},Blu{4}] = lu(Z);
                c = [b;zeros(size(Blu{2},2)-szZ,size(b,2))];
                z = Blu{2}\(Blu{1}\(Blu{3}*(c)));
                y = Blu{4}*z;
            case 3
                % [Blu{1},Blu{2},Blu{3}] = lu(Z);
                c = [b;zeros(size(Blu{2},2)-szZ,size(b,2))];
                y = Blu{2}\(Blu{1}\(Blu{3}*(c)));
        end
    else
        c = [b;zeros(size(Blu,2)-szZ,size(b,2))];
        y=Blu\c;
    end
    x = y(1:szZ,:);
    nitp = 1;
else
    % solve B*y=c
    c = [zeros(size(Blu{2},2),size(b,2));b];
    U = B(1:szZ,szZ+1:end);
    V = B(szZ+1:end,1:szZ);
    Mfun = @(f) appM(Blu,U,V,f,szZ);
    restart = 50; tol = 1e-8; maxit = 100;
    [y,~,~,iter,~] = gmres(Mfun,c,restart,tol,maxit);
    nitp = (iter(1)-1)*restart+iter(2);
    x = y(end-szZ+1:end,:);
end
if opt == 1
    fprintf('num iterations: %d\n',nitp);
end
end

function y = appM(B,U,V,f,szZ)
f1 = f(1:end-szZ,:);
f2 = f(end-szZ+1:end,:);
switch length(B)
    case 5
        y1 = f1+ B{1}\(B{3}*(B{5}\(V*f2)));
        % method 1
        y2 = B{5}\(V*f2);
        y2 = B{1}\(B{3}*y2);
        y2 = -U*(B{4}*(B{2}\y2));
        % method 2
        %y2 = U*(B{4}*(B{2}\f1));
    case 4
        y1 = f1+ B{1}\(B{3}*(V*f2));
        % method 1
        y2 = V*f2;
        y2 = B{1}\(B{3}*y2);
        y2 = -U*(B{4}*(B{2}\y2));
        % method 2
        %y2 = U*(B{4}*(B{2}\f1));
    case 3
        y1 = f1+ B{1}\(B{3}*(V*f2));
        % method 1
        y2 = V*f2;
        y2 = B{1}\(B{3}*y2);
        y2 = -U*(B{2}\y2);
        % method 2
        %y2 = U*(B{2}\f1);
end
y = [y1;y2];
end