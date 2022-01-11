
%  loadData;
clc
b = ones(numOfNodes,1);

start = 1;
stop = 1;
% tic
for i=start:stop
    curttd = ttd(:,i);
    selectedBase = 2;
    basettd = ttd(:,selectedBase);
    A0 = HH(selectedBase*numOfNodes-(numOfNodes-1):selectedBase*numOfNodes,:);
    A = HH(i*numOfNodes-(numOfNodes-1):i*numOfNodes,:);
    %     A0inv = IH(selectedBase*numOfNodes-(numOfNodes-1):selectedBase*numOfNodes,:);
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
    % b <-- A'*b
    % U <-- A'
    % V <-- A'
    % A0 <-- lambda^2*L'*L
    lambda = 1;
%     L = eye(numOfNodes);
    [x_reg, info] = my_smi(A'*b,A',A',lambda^2*(L'*L));
    difReg(i)=norm(x_reg-x);
    %%%%%%%%%%%%%%%%%%%%%%
end
% toc
valX=max(difSM)
valReg=max(difReg)
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

