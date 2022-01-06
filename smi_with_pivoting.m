clc
len = length(TT);
base1Len = 1;
base2Len = len-base1Len+1;
%  base1Len = base2Len;

base1 = TT(:,base1Len);
Abase1 = HH(21*base1Len-20:21*base1Len,:);
Abase1Inv = IH(21*base1Len-20:21*base1Len,:);
base2 = TT(:,base2Len);
Abase2 = HH(21*base2Len-20:21*base2Len,:);
Abase2Inv = IH(21*base2Len-20:21*base2Len,:);

from = 1;
to = len;
difPivit=zeros(1,to);
segPivit =zeros(nbdiode,to);
difAllPivit = zeros(len*21,21);

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
    p=0;
    while p<numOfSw
        p=p+1;
        di = 1 + diag(VV(:,series(p:end)).' * A0Inv * UU(:,series(p:end)));
        [~,index] = max(abs(di));
        index= index + p - 1;
        temp = series(p);
        series(p) = index;
        series(index) = temp;
        denom = 1 + VV(:,series(p)).' * A0Inv * UU(:,series(p));
        denomInv = 1/denom;
        A0Inv = A0Inv - A0Inv *  UU(:,series(p))*  denomInv * VV(:,series(p)).' * A0Inv;
        segPivit(p,i) = max(max(abs(A0Inv-IH(21*i-20:21*i,:))));
    end
    difAllPivit(21*i-20:21*i,:) = A0Inv;
    if(p ~=0)
        difPivit(1,i) = segPivit(p,i);
    end

end
toc
save('difAllPivit.mat','difAllPivit')
save('difPivit.mat','difPivit')


