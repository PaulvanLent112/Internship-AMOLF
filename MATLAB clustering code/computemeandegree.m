% Compute the mean degree of the network
% Input:
% mydistmat: a distancematrix as obtained from pdist
% mythresholds: thresholds for connecting nodes, typically, mylinkmat(i,3)

function averagedegree=computemeandegree(mydistmat,mythresholds)

N=(1+sqrt(1+8*length(mydistmat)))/2;
mts=length(mythresholds);
if N~=mts+1, disp('Thresholds not equal to merging points!'), end

averagedegree=nan(mts,1);
for i=1:mts
    tmp=(mydistmat<mythresholds(i));
    averagedegree(i)=sum(tmp)/N*2;
end