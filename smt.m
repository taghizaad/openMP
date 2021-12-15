clc
dif=zeros(length(TT),1);
totInv = zeros(size(HH,1),size(HH,2));
totInv(21*1-20:21*1,:) = IH(21*1-20:21*1,:);
tic
for i=100:100
    curTT = TT(:,i);
    if(sum(curTT-base1 ==0)>= sum(curTT-base2 ==0))
        baseTT = base1;
        AbaseInv = Abase1Inv;
    else
%         baseTT = base1;
%         AbaseInv = Abase1Inv;
        baseTT = base2;
        AbaseInv = Abase2Inv;
    end
    uv=zeros(numOfNodes,numOfNodes);
    %-----------------smwwd impl---------------------------
    sw = find(curTT-baseTT ~= 0);
    A0Inv = AbaseInv;
    for j=1:length(sw)
        u = zeros(numOfNodes,1);
        v = zeros(numOfNodes,1);
        nodes = find(switchNodeMat(sw(j),:)==1);
        if(length(nodes) == 1)
            u(nodes) = gon;
            v(nodes) = curTT(sw(j)) - baseTT(sw(j));
        else
            u(nodes(1)) = gon;
            u(nodes(2)) = -gon;
            v(nodes(1)) = curTT(sw(j)) - baseTT(sw(j));
            v(nodes(2)) = baseTT(sw(j)) - curTT(sw(j));
        end
        uv = uv+ u*v.';
        vA0Invu=v.'*A0Inv*u;
        D = 1+vA0Invu;
        Dinv = 1/D;
        %update Ainv
        A0Inv = A0Inv - A0Inv * u *  Dinv * v.' * A0Inv;
    end

    dif(i,1) = max(max(abs(IH(21*i-20:21*i,:)-A0Inv)));
    totInv(21*i-20:21*i,:) = A0Inv;
end
toc
% max(max(abs(totInv-IH)))
%calculate 2*2 matrices D and invD
% vAbaseInvu = v.'*AbaseInv*u;
% D = 1 + vAbaseInvu;
% diff(i,2) = cond(D);
% Dinv = 1/D;
% Ainv = AbaseInv - AbaseInv * u * Dinv * v.'* AbaseInv;
% %----------------------------------------------------------------
% diff(i,3) = max(max(abs(IH(21*i-20:21*i,:)-Ainv)));
%fprintf('***%d***\n', i);




%bar = IH(21*first-20:21*last,:);
%max(max(abs(bar-Ainv)))
%fprintf('*****************%d*************************\n', i);

% aa = zeros(length(TT),3);
% for i=1:length(TT)
%     bb = HH(21*i-20:21*i,:)-HH(21*1-20:21*1,:);
%     aa(i,1) = rank(bb);
%     bb( :, ~any(bb,1) ) = [];
%     aa(i,2) = size(bb,2);
%
% end
