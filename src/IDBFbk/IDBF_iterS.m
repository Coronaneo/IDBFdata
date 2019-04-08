function [FactorS,rv,cv] = IDBF_iterS(Faccell,flag)
% called in IDBBF_arrg.m
%
% Copyright 2018 by Qiyuan Pang
   prodnpx = size(Faccell,1);
   prodnpk = size(Faccell,2);
   if prodnpx*prodnpk > 1
     [FactorS1,rv1,cv1] = IDBF_iterS(Faccell(1:prodnpx/2,1:prodnpk/2),flag);
     [FactorS2,rv2,cv2] = IDBF_iterS(Faccell(1:prodnpx/2,prodnpk/2+1:prodnpk),flag);
     [FactorS3,rv3,cv3] = IDBF_iterS(Faccell(prodnpx/2+1:prodnpx,1:prodnpk/2),flag);
     [FactorS4,rv4,cv4] = IDBF_iterS(Faccell(prodnpx/2+1:prodnpx,prodnpk/2+1:prodnpk),flag);
     row = [size(FactorS1,1) size(FactorS2,1) size(FactorS3,1) size(FactorS4,1)];
     col = [size(FactorS1,2) size(FactorS3,2) size(FactorS2,2) size(FactorS4,2)];
     FactorS = sparse(sum(row),sum(col));
     FactorS(1:row(1),1:col(1)) = FactorS1;
     FactorS(row(1)+1:row(1)+row(2),col(1)+col(2)+1:col(1)+col(2)+col(3)) = FactorS2;
     FactorS(row(1)+row(2)+1:row(1)+row(2)+row(3),col(1)+1:col(1)+col(2)) = FactorS3;
     FactorS(row(1)+row(2)+row(3)+1:sum(row),col(1)+col(2)+col(3)+1:sum(col)) = FactorS4;
   
     rv2(:,2) = rv2(:,2) + prodnpk/2;
     rv3(:,1) = rv3(:,1) + prodnpx/2;
     rv4(:,1) = rv4(:,1) + prodnpx/2;
     rv4(:,2) = rv4(:,2) + prodnpk/2;
     cv2(:,2) = cv2(:,2) + prodnpk/2;
     cv3(:,1) = cv3(:,1) + prodnpx/2;
     cv4(:,1) = cv4(:,1) + prodnpx/2;
     cv4(:,2) = cv4(:,2) + prodnpk/2;
     rv = [rv1;rv2;rv3;rv4];
     cv = [cv1;cv3;cv2;cv4];
   else
     if flag > 0
        FactorS = Faccell{1,1}.S;
     else
        FactorS = Faccell{1,1};
     end
     rv = [1 1];
     cv = [1 1];
  end

    
end