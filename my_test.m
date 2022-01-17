
loadData
clc
b = ones(numOfNodes,1);
start = 1;
stop = len;
difSM=[];
% allDen = zeros(18,len);
selectedBase = len; %bin2dec('111111111000000000')+1;
basettd = ttd(:,selectedBase);
% tic
for i=start:stop
    curttd = ttd(:,i);
    A0 = HH(selectedBase*numOfNodes-(numOfNodes-1):selectedBase*numOfNodes,:);
    A = HH((i-1)*numOfNodes+1:i*numOfNodes,:);
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
    [xSM,Y] = my_smi(b,U,V,A0);
    x = IH((i-1)*numOfNodes+1:i*numOfNodes,:)*b;
    difSM(i)=norm(xSM-x);
%     difSMrel(i)=norm(x_ls-x)/norm(x);
%     AllCount(:,i)= counters;
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
%     lambda = 1;
%     L = eye(numOfNodes);
%     [x_reg, info] = my_smi(A'*b,A',A',lambda^2*(L'*L));
%     difReg(i)=norm(x_reg-x);
    %%%%%%%%%%%%%%%%%%%%%%
end
% toc
valSM=max(difSM)
% minDen=zeros(1,len);
% for t=start:stop
%     if(t==selectedBase)
%         continue;
%     end
%     foo = allDen(:,t);
%     minDen(t) = min(foo(foo>0));
% end
figure
semilogy(0:len-1,difSM,'*')
xlim([0 len-1])
ylim([1e-8 1])
% % title(dec2bin(selectedBase-1))
% xlabel('matrix number')
% ylabel('|| xSM - x ||')
% bin2dec('111111111000000000')

% detVal = zeros(1,len);
% for t=1:len
%     detVal(t) = det(HH((t-1)*numOfNodes+1:t*numOfNodes,:));
% end
