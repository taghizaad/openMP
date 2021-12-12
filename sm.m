function Ainv = sm(baseTT,cTT,AbaseInv,switchNodeMat,gon,numOfNodes)

sw = find(cTT-baseTT ~= 0);
%admittance = ttd(sw,i);
S = length(sw); %number of switch mismatch
% what if S==0
A0Inv = AbaseInv;
for j=1:S
    %obtain u and v
    [u,v]=generateUandV(baseTT,cTT,switchNodeMat,gon,numOfNodes,sw(j));
    uv = u*v;
    %calculate 2*2 matrices D and invD
    vA0Invu=v*A0Inv*u;
    D = eye(size(u,2))+vA0Invu;
    dett = det(D);
    Dinv = inv(D);
    %update Ainv
    q = Dinv * v* A0Inv;
    p = A0Inv * u * q;
    A0Inv = A0Inv - p;
    
end
Ainv = A0Inv;
end

