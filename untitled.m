clc
Ainv=Abase1Inv;
first=5537;
last=5537;
diff=zeros(last,1);
tic
for i=first:last
    curTT = TT(:,i);
    %selectedBase = chooseBase(base1,base2,curTT);
    selectedBase = base1;
    if(all(selectedBase == base1))
        AbaseInv = Abase1Inv;
    else
        AbaseInv = Abase2Inv;
    end
    AbaseInvResh = reshape(AbaseInv.',1,[]);
    %fprintf('*****************%d*************************\n', i);
    AinvResh = smwWithWholeD(numOfNodes,numOfSwitches,gon,selectedBase,curTT,switchNodeMatResh,AbaseInvResh);
    %ss(i,:) = [i invArr selectedBase(1)];
    tempAinv = transpose(reshape(AinvResh,size(AbaseInv,2),[]));
    diff(i) = max(max(abs(IH(21*i-20:21*i,:)-tempAinv)));
    Ainv = [Ainv;tempAinv];
    %fprintf('***%d***\n', i);
end
toc
max(diff)
%bar = IH(21*first-20:21*last,:);
%max(max(abs(bar-Ainv)))