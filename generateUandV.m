function [u,v] = generateUandV(baseTT,cTT,switchNodeMat,gon,numOfNodes,switchNum)


nodes = find(switchNodeMat(switchNum,:)==1);

u = zeros(numOfNodes,length(nodes));
v = zeros(length(nodes),numOfNodes);
if(length(nodes) == 1)
    u(nodes) = (cTT(switchNum) - baseTT(switchNum)) * gon;
    v(nodes) = 1;
else
    u(nodes(1),1)=(cTT(switchNum) - baseTT(switchNum)) * gon;
    u(nodes(1),2)=(cTT(switchNum) - baseTT(switchNum)) * (-gon);
    u(nodes(2),1)=(cTT(switchNum) - baseTT(switchNum)) * (-gon);
    u(nodes(2),2)=(cTT(switchNum) - baseTT(switchNum)) * gon;


    v(1,nodes(1))=1;
    v(2,nodes(2))=1;
end
end



