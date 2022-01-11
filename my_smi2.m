function [x,Y] = my_smi2(b,U,V,A0inv,epsilon)
Y = [];
M = size(V,2);
%(i)
x = A0inv*b;
%(ii)
if M <= 0
    return
end
Y = A0inv*U;
%(iii)
for l=1:M-1
%     f = 1+epsilon+V(:,l)'*Y(:,l);
    x = x - Y(:,l) * ((V(:,l)'*x) / (1+V(:,l)'*Y(:,l)+epsilon));
    Y(:,l+1:M) = Y(:,l+1:M)- Y(:,l) * ((V(:,l)'*Y(:,l+1:M)) / (1+V(:,l)'*Y(:,l)+epsilon));
end
%(iv)
x = x - Y(:,M) * ((V(:,M)'*x) / (1+V(:,M)'*Y(:,M)+epsilon));