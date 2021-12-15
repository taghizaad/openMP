clc
first=202437;
last=202437;
diff=zeros(last,3);
baseTT = base1;
AbaseInv = Abase1Inv;
tic
for i=first:last
    curTT = TT(:,i);
    %selectedBase = chooseBase(base1,base2,curTT);
%     baseTT = base1;
%     if(all(baseTT == base1))
%         AbaseInv = Abase1Inv;
%     else
%         AbaseInv = Abase2Inv;
%     end
    
    %-----------------smwwd impl---------------------------
    sw = find(curTT-baseTT ~= 0);
    u = zeros(numOfNodes,numOfNodes);
    v = zeros(numOfNodes,numOfNodes);
    % what if S==0
    for j=1:length(sw)
        nodes = find(switchNodeMat(sw(j),:)==1);
        if(length(nodes) == 1)
            u(nodes,nodes) = u(nodes,nodes) + (curTT(sw(j)) - baseTT(sw(j))) * gon;
            v(nodes,nodes) = 1;
        else
            u(nodes(1),nodes(1)) = u(nodes(1),nodes(1)) + (curTT(sw(j)) - baseTT(sw(j))) * gon;
            u(nodes(1),nodes(2))= u(nodes(1),nodes(2)) + (curTT(sw(j)) - baseTT(sw(j))) * (-gon);
            u(nodes(2),nodes(1))= u(nodes(2),nodes(1)) + (curTT(sw(j)) - baseTT(sw(j))) * (-gon);
            u(nodes(2),nodes(2))= u(nodes(2),nodes(2)) + (curTT(sw(j)) - baseTT(sw(j))) * gon;
            v(nodes(1),nodes(1))=1;
            v(nodes(2),nodes(2))=1;
        end
    end

    %prune u and v
    u( :, ~any(u,1) ) = [];
    v( :, ~any(v,1) ) = [];
    uv = u*v.';
    diff(i,1) = max(max(abs(uv - (HH(21*i-20:21*i,:)-HH(21*1-20:21*1,:)))));
    %calculate 2*2 matrices D and invD
    vAbaseInvu=v.'*AbaseInv*u;
    D = eye(size(u,2),size(u,2))+vAbaseInvu;
    diff(i,2) = cond(D);
    Dinv = inv(D);
    Ainv = AbaseInv - AbaseInv * u * Dinv * v.'* AbaseInv;
    %----------------------------------------------------------------
    diff(i,3) = max(max(abs(IH(21*i-20:21*i,:)-Ainv)));
    %fprintf('***%d***\n', i);
end
toc
max(diff(:,1))
max(diff(:,2))
max(diff(:,3))


%bar = IH(21*first-20:21*last,:);
%max(max(abs(bar-Ainv)))
%fprintf('*****************%d*************************\n', i);

aa = zeros(length(TT),3);
for i=1:length(TT)
    bb = HH(21*i-20:21*i,:)-HH(21*1-20:21*1,:);
    aa(i,1) = rank(bb);
    bb( :, ~any(bb,1) ) = [];
    aa(i,2) = size(bb,2);

end
