clear
clc
All3levelmatrixML
load("switchNodeMat.mat");
switchNodeMatResh = reshape(switchNodeMat.',1,[]);
numOfSwitches = size(TT,1);
numOfNodes = size(HH,2);
len = length(ttd);
