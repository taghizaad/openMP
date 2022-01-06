
difRandom=zeros(1,to);
seg =zeros(nbdiode,to);
difAllRandom = zeros(len*21,21);

tic
for i=from:to
    curTT = TT(:,i);
    AcurInv = IH(21*(i)-20:21*(i),:);

    if(sum(curTT-base1 ==0)>= sum(curTT-base2 ==0))
        base = base1;
        AbaseInv = Abase1Inv;
        Abase = Abase1;
        baseNumber = base1Len;
    else
        base = base2;
        AbaseInv = Abase2Inv;
        Abase = Abase2;
        baseNumber = base2Len;
    end

    %-----------------smwwd impl---------------------------
    sw = find(curTT-base ~= 0);
    numOfSw = length(sw);
    A0Inv = AbaseInv;
    U = zeros(numOfNodes,numOfSw);
    V = zeros(numOfNodes,numOfSw);
    for j=1:numOfSw
        u = zeros(numOfNodes,1);
        v = zeros(numOfNodes,1);
        nodes = find(switchNodeMat(sw(j),:)==1);
        if(length(nodes) == 1)
            u(nodes) = ttd(sw(j),i)-ttd(sw(j),baseNumber);
            v(nodes) = 1;
        else
            u(nodes(1)) = ttd(sw(j),i)-ttd(sw(j),baseNumber);
            u(nodes(2)) = -1 * (ttd(sw(j),i)-ttd(sw(j),baseNumber));
            v(nodes(1)) = 1;
            v(nodes(2)) =-1;
        end
        U(:,j)=u;
        V(:,j)=v;
    end
    UU = U;
    VV = V;
    series = 1:numOfSw;
    random = series(randperm(length(series)));
    p=0;
    while p<numOfSw
        p=p+1;
        denom = 1 + VV(:,random(p)).' * A0Inv * UU(:,random(p));
        denomInv = 1/denom;
        A0Inv = A0Inv - A0Inv *  UU(:,random(p))*  denomInv * VV(:,random(p)).' * A0Inv;
        seg(p,i) = max(max(abs(A0Inv-IH(21*i-20:21*i,:))));
    end
    difAllRandom(21*i-20:21*i,:) = A0Inv;
    if(p ~=0)
        difRandom(1,i) = seg(p,i);
    end

end
toc
% save('difAllRandom.mat','difAllRandom')
% save('difRandom.mat','difRandom')


