
loadData;
clc
base1Ind = 1;
base2Ind = len;
base1 = ttd(:,base1Ind);
Z01inv = HH(base1Ind*numOfNodes-(numOfNodes-1):base1Ind*numOfNodes,:);
Z01 = IH(base1Ind*numOfNodes-(numOfNodes-1):base1Ind*numOfNodes,:);
base2 = ttd(:,base2Ind);
Z02inv = HH(base2Ind*numOfNodes-(numOfNodes-1):base2Ind*numOfNodes,:);
Z02 = IH(base2Ind*numOfNodes-(numOfNodes-1):base2Ind*numOfNodes,:);

X_ls = zeros(numOfNodes,len);
X = zeros(numOfNodes,len);
dif_ls = zeros(1,len);
dif = zeros(1,len);


start = 1;
stop = len;
tic
for i=start:stop
    curttd = ttd(:,i);
    b = B(:,i);
    if(sum(abs(curttd-base1) == 0) >= sum(abs(curttd-base2) == 0))
        baseInd = base1Ind;
        basettd = base1;
        Z0 = Z01;
    else
        baseInd = base2Ind;
        basettd = base2;
        Z0 = Z02;
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
    [x_ls,info] = smi(b,U,V,Z0);
    X_ls(:,i) = x_ls;
    X(:,i) = IH((i-1)*numOfNodes+1:i*numOfNodes,:)*b;
    dif_ls(i) = norm(x_ls-x);
    dif(i)= norm(X(:,i) -x);
    
end
toc
[val_ls,ind_ls]=max(dif_ls);
[val,ind]=max(dif);
Z = info.Z;
I = eye(numOfNodes);
den = (1+V(:,1)'*Z(:,1));
num = Z(:,1)*V(:,1)';
invi = I - num/den;
for i=2:size(V,2)
    den = 1+V(:,i)'*Z(:,i);
    num = Z(:,i)*V(:,i)';
    invi = invi * (I - num/den);
end

foo = invi * Z0;
bar = IH((start-1)*numOfNodes+1:start*numOfNodes,:);

norm(foo-bar)
