
%  loadData;
clc
base = ones(numOfSwitches,1);
for t=1:numOfSwitches-1
    base(t+1) = 2^(t)+1;
end

b = ones(numOfNodes,1);

start = 1;
stop = len;

difSM = zeros(1,stop-start+1);
% difInv = zeros(1,stop-start+1);
% difReg = zeros(1,stop-start+1);

tic
for i=start:stop
    curttd = ttd(:,i);
    if base(1)<=i && i< base(2)
        selectedBase = base(1);
    elseif base(2)<=i && i< base(3)
        selectedBase = base(2);
    elseif base(3)<=i && i< base(4)
        selectedBase = base(3);
    elseif base(4)<=i && i< base(5)
        selectedBase = base(4);
    elseif base(5)<=i && i< base(6)
        selectedBase = base(5);
    elseif base(6)<=i && i< base(7)
        selectedBase = base(6);
    elseif base(7)<=i && i< base(8)
        selectedBase = base(7);
    elseif base(8)<=i && i< base(9)
        selectedBase = base(8);
    elseif base(9)<=i && i< base(10)
        selectedBase = base(9);
    elseif base(10)<=i && i< base(11)
        selectedBase = base(10);
    elseif base(11)<=i && i< base(12)
        selectedBase = base(11);
    elseif base(12)<=i && i< base(13)
        selectedBase = base(12);
    elseif base(13)<=i && i< base(14)
        selectedBase = base(13);
    elseif base(14)<=i && i< base(15)
        selectedBase = base(14);
    elseif base(15)<=i && i< base(16)
        selectedBase = base(15);
    elseif base(16)<=i && i< base(17)
        selectedBase = base(16);
    elseif base(17)<=i && i< base(18)
        selectedBase = base(17);
    else 
        selectedBase = base(18);
    end

    basettd = TT(:,selectedBase);
    A0 = HH(selectedBase*numOfNodes-(numOfNodes-1):selectedBase*numOfNodes,:);

    %-----------------smwwd impl---------------------------
    sw = find(TT(:,i) ~= TT(:,selectedBase));
    U = zeros(numOfNodes,length(sw));
    V = zeros(numOfNodes,length(sw));
    for j=1:length(sw)
        nodes = find(switchNodeMat(sw(j),:)==1);
        if(length(nodes) == 1)
            U(nodes,j) = curttd(sw(j)) - basettd(sw(j));
            V(nodes,j) = 1;
        else
            U(nodes(1),j) = curttd(sw(j)) - basettd(sw(j));
            V(nodes(1),j) = 1;
            U(nodes(2),j) = basettd(sw(j)) - curttd(sw(j));
            V(nodes(2),j) = -1;
        end
    end
    [x_ls,Y] = my_smi(b,U,V,A0);
    x = IH((i-1)*numOfNodes+1:i*numOfNodes,:)*b;
    difSM(i)=norm(x_ls-x);

    %     inverse = eye(numOfNodes);
    %     for j=size(V,2):-1:1
    %         inverse = inverse * (eye(numOfNodes) - (Y(:,j)*V(:,j)' / (1+V(:,j)'*Y(:,j))));
    %     end
    %     inverse = inverse * A0inv;
    %     difInv(i) = norm(inverse- IH((i-1)*numOfNodes+1:i*numOfNodes,:));

    %%% Regularization %%%
    %     A = HH((i-1)*numOfNodes+1:i*numOfNodes,:);
    %     g = A'*b;
    %     [x_reg, info] = shermanMorrisonIteration(g,A,Z0);
    %     difReg(i)=norm(x_reg-x);

    %%%%%%%%%%%%%%%%%%%%%%
end
toc
[valX,indX]=max(difSM)
% [valReg,indReg]=max(difReg);
% [valInv,indInv]=max(difInv);

% close all
% subplot(2,1,1)
% axis = 1:len;
% plot(axis,difSM)
% title('norm(x - xSM)')
%
% subplot(2,1,2)
% plot(axis,difReg)
% title('norm(x - xReg)')

