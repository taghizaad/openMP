
% loadData;
clc
base1Ind = 1;
base2Ind = len;
base1 = ttd(:,base1Ind);
A01 = HH(base1Ind*numOfNodes-(numOfNodes-1):base1Ind*numOfNodes,:);
A01inv = IH(base1Ind*numOfNodes-(numOfNodes-1):base1Ind*numOfNodes,:);
base2 = ttd(:,base2Ind);
A02 = HH(base2Ind*numOfNodes-(numOfNodes-1):base2Ind*numOfNodes,:);
A02inv = IH(base2Ind*numOfNodes-(numOfNodes-1):base2Ind*numOfNodes,:);

b = ones(numOfNodes,1);

difX = zeros(1,len);
difInv = zeros(1,len);
difReg = zeros(1,len);

start = 1;
stop = len;

lambda = 1e+5;
L = eye(numOfNodes);
Z0 = 0*L;
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
    Y = my_smi3(b,U,V,A0inv,0);
%     x = IH((i-1)*numOfNodes+1:i*numOfNodes,:)*b;
    difX(i)=norm(Y-IH((i-1)*numOfNodes+1:i*numOfNodes,:));
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
      Y2 = my_smi3(b,U,V,A0inv,1e-17);
    difReg(i)=norm(Y2-IH((i-1)*numOfNodes+1:i*numOfNodes,:));

    %%%%%%%%%%%%%%%%%%%%%%
end
toc
[valX,indX]=max(difX);
[valInv,indInv]=max(difInv);

close all
subplot(3,1,1)
axis = 1:len;
plot(axis,difX)
title('norm(x - xSM)')


subplot(3,1,2)
plot(axis,difReg)
title('norm(x - xReg)')

subplot(3,1,3)
plot(axis,difX-difReg)
title('difX-difReg')

mean(difReg<difX)
figure
plot(difX,difReg-difX,'.')
xlabel('difX')
ylabel('difReg-difX')

foo2=zeros(1,len);
foo = zeros(1,len);
for t = 1:len
    A=HH((t-1)*numOfNodes+1:t*numOfNodes,:);
    foo(t) = max(max(abs( inv(inv(A))-A)));
    foo2(t) = max(max(abs(A)));
end
figure
plot(log10(foo))
figure
plot(log10(foo./foo2))

%%
warning off
lambda = 1e-11;

foo = zeros(1,len);
foo2=zeros(1,len);
foo3=zeros(1,len);
for t = 1:10
    A=HH((t-1)*numOfNodes+1:t*numOfNodes,:);
%     if mod(t,1000)==0
%         display(num2str(t))
%     end
    foo(t) = max(max(abs(inv(inv(A'*A+lambda*eye(numOfNodes))*A')-A)));
    foo2(t) = max(max(abs(A)));
    invMat = double(vpa(sym(double(vpa(inv(sym(A)),100))) - inv(sym(A)),100));
    invInvMat = double(vpa(sym(double(vpa(inv(sym(invMat)),100))) - inv(sym(invMat)),100));
    foo3(t) = max(max(abs(invInvMat)));
end
close all
figure
plot(log10(foo))
ylim([-6 -5])
figure
plot(log10(foo./foo2))
figure
plot(log10(foo2))
warning on



