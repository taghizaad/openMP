function [x,Y] = my_smi(b,U,V,A0)
Y = [];
M = size(V,2);
%(i)
x = inv(A0)*b;
%(ii)
if M <= 0
    return
end
% cols = size(V,2);
% P = randperm(cols);
% U = U(:,P);
% V= V(:,P);

Y = inv(A0)*U;
%(iii)
count = 0;
for l=1:M-1
    den = 1+V(:,l)'*Y(:,l);
    l
    Y
     V(:,l)'*Y(:,l)
%     V'*Y
    if(den<1e-5)
        count = count+1;
    end
    x = x - Y(:,l) * ((V(:,l)'*x) / den);
    Y(:,l+1:M) = Y(:,l+1:M)- Y(:,l) * ((V(:,l)'*Y(:,l+1:M)) / den);
end
%(iv)
l=M
den = 1+V(:,l)'*Y(:,l);
 V(:,l)'*Y(:,l)
if(den<1e-5)
    count = count+1;
end
x = x - Y(:,l) * ((V(:,l)'*x) / den);
count