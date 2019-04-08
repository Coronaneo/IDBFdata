m = 32;
M = 10;
r = 40;
n = 10;
N = 2^n;
ww = (-N/2:N/2-1)';
ns = rand(N,1)/2;
ww = ww + ns;
rr = ((0:N-1)/N)';
ns = rand(N,1)/N/2;
rr = rr + ns;
Z = @(x,k)funFT(rr(x),ww(k));
A = Z([1:N]',[1:N]');
B = @(x,k)funFT(rr,ww);
F = BF_IDBF(Z,[1:N]',[1:N]',m,r,1E-6,'rand',20,1);
Factor = IBF_Cheby(B,ww,rr,10,1E-8);
x = randn(N,1);
err = norm(BF_apply(Factor,x)-B(rr,ww)*x)/norm(B(rr,ww)*x)
accuracy = norm(BF_apply(F,x)-A*x)/norm(A*x)

train_data = zeros(N,M);
target_data = zeros(N,M);
for i = 1:M
    x = randn(N,1);
    train_data(:,i) = x;
    target_data(:,i) = A*x;
end

save train_data.mat train_data 
save target_data.mat target_data 


