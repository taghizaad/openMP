clc
dif=zeros(length(TT),1);
totInv = zeros(size(HH,1),size(HH,2));
totInv(21*1-20:21*1,:) = IH(21*1-20:21*1,:);
tic
for i=230646:230646
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
    %-----------------smwwd impl---------------------------
    sw = find(curTT-baseTT ~= 0);
    U = zeros(numOfNodes,length(sw));
    V = zeros(numOfNodes,length(sw));
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
        U(:,j)=u;
        V(:,j)=v;
    end
    p = 1:length(sw);
    A0Inv = AbaseInv;
    for k=1:length(sw)
        t = 1;
        while (1 + V(:,p(k)).'*A0Inv * U(:,p(k)) < 1e-10 && k+t<=length(p))
            temp = p(k);
            p(k) = p(k+t);
            p(k+t) = temp;
            t = t+1;
        end
        vA0Invu= V(:,p(k)).'*A0Inv * U(:,p(k));
        D = 1+vA0Invu;
        Dinv = 1/D;
        A0Inv = A0Inv - A0Inv * U(:,p(k)) *  Dinv * V(:,p(k)).' * A0Inv;

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
