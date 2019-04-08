function [U,V] = BF_Phase_Correction_Mat_1(U,V,posu,posv,disPosu,disPosv,disPosSubu,disPosSubv,tau,dim)
% O(N log N) operation and memory complexity.
% Reference: H. Yang, A Unified Framework for Oscillatory Integral 
% Transform: When to use NUFFT or Butterfly factorization? preprint, 2018.


% U for rows and V for columns, all talk skinny matrices
[mu,nu] = size(U);
[mv,nv] = size(V);
switch dim
    case 1
        % correct first row in each continuous sector of rows
        for cntd = 1:numel(disPosu)
            posc = disPosSubu(cntd);
            V(:,posc) = BF_Phase_Correction_Vec_1(V(:,posc),disPosv,tau);
        end

        % correct first three columns in each continuous sector of columns
        for cntd = 1:numel(disPosv)
            for cnt1 = 1:3
                pos = cnt1+disPosSubv(cntd)-1;
                val = V(disPosv(cntd)+cnt1-1,disPosSubu);
                U(:,pos) = BF_Phase_Correction_Vec_2(U(:,pos),tau,disPosu,val);
            end
            % for loop for each continuous sector of rows
            for cntdd = 1:numel(disPosu)
                app = U(disPosu(cntdd),disPosSubv(cntd)+1) + U(disPosu(cntdd)+1,disPosSubv(cntd)) - U(disPosu(cntdd),disPosSubv(cntd));
                slope = round((app-U(disPosu(cntdd)+1,disPosSubv(cntd)+1))/tau);
                if cntdd == numel(disPosu)
                    len = mv-disPosu(cntdd)+1;
                    ed = mv;
                else
                    len = disPosu(cntdd+1)-disPosu(cntdd)+1;
                    ed = disPosu(cntdd+1)-1;
                end
                st = disPosu(cntdd);
                vecAdd = tau*(0:(len-1))*slope;
                U(st:ed,disPosSubv(cntd)+1) = U(st:ed,disPosSubv(cntd)+1) + vecAdd';
                app = U(disPosu(cntdd),disPosSubv(cntd)+2) + U(disPosu(cntdd)+1,disPosSubv(cntd)+1) - U(disPosu(cntdd),disPosSubv(cntd)+1);
                slope = round((app-U(disPosu(cntdd)+1,disPosSubv(cntd)+2))/tau);
                vecAdd = tau*(0:(len-1))*slope;
                U(:,disPosSubv(cntd)+2) = U(:,disPosSubv(cntd)+2) + vecAdd';
            end
        end
        
        %correct all other rows
        for cntd = 1:numel(disPosu)
            if cntd == numel(disPosu)
                ed = nv;
            else
                ed = disPosSubu(cntd+1)-1;
            end
            for cnt1 = disPosSubu(cntd)+1:ed
                pos1 = disPosv; pos2 = posv(cnt1);
                val1 = U(pos2,disPosSubv); val2 = U(pos2,disPosSubv+1); val3 = U(pos2,disPosSubv+2);
                scnDer = (val1+val3-2*val2);
                V(:,cnt1) = BF_Phase_Correction_Vec_3(V(:,cnt1),tau,pos1,val1,pos1+1,val2,pos1+2,val3,scnDer);
            end
        end
        
        % correct all other columns
        for cntd = 1:numel(disPosv)
            if cntd == numel(disPosv)
                ed = nu;
            else
                ed = disPosSubv(cntd+1)-1;
            end
            for cnt1 = disPosSubv(cntd)+3:ed
                pos1 = disPosu; pos2 = posu(cnt1);
                val1 = V(pos2,disPosSubu); val2 = V(pos2,disPosSubu+1); val3 = V(pos2,disPosSubu+2);
                scnDer = (val1+val3-2*val2);
                U(:,cnt1) = BF_Phase_Correction_Vec_3(U(:,cnt1),tau,pos1,val1,pos1+1,val2,pos1+2,val3,scnDer);
            end
        end
        
    case 2
        % correct first row
        cnt1 = 1;
        for cnt2 = 1:mv^(1/dim)
            if cnt2>1
                V((cnt2-1)*n+1,cnt1) = round((V((cnt2-2)*n+1,cnt1)-V((cnt2-1)*n+1,cnt1)+tau/2)/tau)*tau + V((cnt2-1)*n+1,cnt1)-tau;
            end
          %  V(((cnt2-1)*n+1):(cnt2*n),cnt1) = BF_Phase_Correction_Vec_1(V(((cnt2-1)*n+1):(cnt2*n),cnt1),tau);
        end
        % correct all columns
        for cnt1 = 1:nu
            for cnt2 = 1:mu^(1/dim)
                if cnt2>1
                    U((cnt2-1)*n+1,cnt1) = round((U((cnt2-2)*n+1,cnt1)-U((cnt2-1)*n+1,cnt1)+tau/2)/tau)*tau + U((cnt2-1)*n+1,cnt1)-tau;
                end
          %      U(((cnt2-1)*n+1):(cnt2*n),cnt1) = BF_Phase_Correction_Vec_1(U(((cnt2-1)*n+1):(cnt2*n),cnt1),tau);
            end
        end
        %correct all other rows
        for cnt1 = 1:nv
            for cnt2 = 1:mv^(1/dim)
                if cnt2>1
                    V((cnt2-1)*n+1,cnt1) = round((V((cnt2-2)*n+1,cnt1)-V((cnt2-1)*n+1,cnt1)+tau/2)/tau)*tau + V((cnt2-1)*n+1,cnt1)-tau;
                end
            %    V(((cnt2-1)*n+1):(cnt2*n),cnt1) = BF_Phase_Correction_Vec_1(V(((cnt2-1)*n+1):(cnt2*n),cnt1),tau);
            end
        end
        
        
    case 3
        
end

end


