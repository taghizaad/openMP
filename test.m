clear
clc
non = 21;
nos = 18;
len = 2^nos;
b = ones(non,1);
dd = zeros(1,len);
Aall=zeros(non*len,non);
Uall=zeros(non*len,non);
V = eye(non);
for i=1:len
    A = rand(non);
    A0 = diag(diag(A));
    U = A-A0;
    [xSM,Y] = my_smi(b,U,V,A0);
    x = inv(A)*b;
    dd(i) = norm(xSM-x);
    Aall((i-1)*non+1:i*non,:) = A;
    Uall((i-1)*non+1:i*non,:) = U;
end
[val,ind]= max(dd)
% semilogy(dd,'*')
% xlabel('matrix number')
% ylabel('|| xSM - x ||')



