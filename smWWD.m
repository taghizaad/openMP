function result = smWWD(baseTT,cTT,AbaseInv,switchNodeMat,gon,numOfNodes)

sw = find(cTT-baseTT ~= 0);
u = zeros(numOfNodes,numOfNodes);
v = zeros(numOfNodes,numOfNodes);
% what if S==0
for j=1:length(sw)
    nodes = find(switchNodeMat(sw(j),:)==1);
    if(length(nodes) == 1)
        u(nodes) = u(nodes) + (cTT(sw(j)) - baseTT(sw(j))) * gon;
        v(nodes) = 1;
    else
        u(nodes(1),nodes(1)) = u(nodes(1),nodes(1)) + (cTT(sw(j)) - baseTT(sw(j))) * gon;
        u(nodes(1),nodes(2))= u(nodes(1),nodes(2)) + (cTT(sw(j)) - baseTT(sw(j))) * (-gon);
        u(nodes(2),nodes(1))= u(nodes(2),nodes(1)) + (cTT(sw(j)) - baseTT(sw(j))) * (-gon);
        u(nodes(2),nodes(2))= u(nodes(2),nodes(2)) + (cTT(sw(j)) - baseTT(sw(j))) * gon;
        v(nodes(1),nodes(1))=1;
        v(nodes(2),nodes(2))=1;
    end
end

    %prune u and v
    u( :, ~any(u,1) ) = [];
    v( :, ~any(v,1) ) = [];
    uv = u*v.';
    %calculate 2*2 matrices D and invD
    vAbaseInvu=v.'*AbaseInv*u;
    D = eye(size(u,2),size(u,2))+vAbaseInvu;
    dett = det(D);
    Dinv = inv(D);
    result = AbaseInv - AbaseInv * u * Dinv * v.'* AbaseInv;
    
end

