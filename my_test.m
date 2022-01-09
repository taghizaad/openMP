
% loadData;
clc
base1Ind = 1;
base2Ind = 1;
base1 = ttd(:,base1Ind);
A01 = HH(base1Ind*numOfNodes-(numOfNodes-1):base1Ind*numOfNodes,:);
A01inv = IH(base1Ind*numOfNodes-(numOfNodes-1):base1Ind*numOfNodes,:);
base2 = ttd(:,base2Ind);
A02 = HH(base2Ind*numOfNodes-(numOfNodes-1):base2Ind*numOfNodes,:);
A02inv = IH(base2Ind*numOfNodes-(numOfNodes-1):base2Ind*numOfNodes,:);

b = ones(numOfNodes,1);

difSM = zeros(1,len);
difInv = zeros(1,len);
difReg = zeros(1,len);

start = 244765;
stop = 244765;

tic
for i=start:stop
    curttd = ttd(:,i);
    if(sum(abs(curttd-base1) == 0) >= sum(abs(curttd-base2) == 0))
        baseInd = base1Ind;
        basettd = base1;
        A0 = A01;
        A0inv = A01inv;
    else
        baseInd = base2Ind;
        basettd = base2;
        A0 = A02;
        A0inv = A02inv;
    end
    %-----------------smwwd impl---------------------------
    sw = find(TT(:,i) ~= TT(:,baseInd));
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
    A = HH((i-1)*numOfNodes+1:i*numOfNodes,:);
    g = A'*b;
    [x_reg, info] = shermanMorrisonIteration(g,A,Z0);
    difReg(i)=norm(x_reg-x);

    %%%%%%%%%%%%%%%%%%%%%%
end
toc
[valX,indX]=max(difSM)
[valReg,indReg]=max(difReg);
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

