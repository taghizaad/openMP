clear
clc
All3levelmatrixML
load("switchNodeMat.mat");
switchNodeMatResh = reshape(switchNodeMat.',1,[]);
numOfSwitches = size(TT,1);
numOfNodes = size(HH,2);
base2 = TT(:,end);
base1 = TT(:,1);
Abase1Inv = IH(1:21,:);
Abase2Inv = IH(end-20:end,:);
Abase1 = HH(1:21,:);
Abase2 = HH(end-20:end,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ss = zeros(length(TT)-2,4);
ss(1,:)=[];
Ainv=[];
tic
%for i=2:length(TT)-1
for i=2:100
    curTT = TT(:,i);
    selectedBase = chooseBase(base1,base2,curTT);
    %selectedBase = base1;
    if(all(selectedBase == base1))
        AbaseInv = Abase1Inv;
    end
    if(all(selectedBase == base2))
        AbaseInv = Abase2Inv;
    end
    AbaseInvResh = reshape(AbaseInv.',1,[]);
    %fprintf('*****************%d*************************\n', i);
    [AinvResh,invArr] = smwWithWholeD(numOfNodes,numOfSwitches,gon,selectedBase,curTT,switchNodeMatResh,AbaseInvResh);
    %ss(i,:) = [i invArr selectedBase(1)];
    Ainv = [Ainv;transpose(reshape(AinvResh,size(AbaseInv,2),[]))];
end

toc

 dim1Prob = find(ss(:,2));
 dim2Prob = find(ss(:,3));
%--------------------------------------------------------------------------------------------------------------%

clc
cTT = TT(:,dim2Prob(2));
baseTT = base1;
AbaseInv = Abase1Inv;
Ainv =sm(baseTT,cTT,AbaseInv,switchNodeMat,gon,numOfNodes);





