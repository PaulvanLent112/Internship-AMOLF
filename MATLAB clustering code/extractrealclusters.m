%Function to extract all instances where a branch is above the thresholds
%needed for extractallclusters6

function realclusters=extractrealclusters(position,linkmat,gn,...
    clustersize,meansizethresh,falseposcontrol,abovethresh)

%compute the fraction 'falseposcontrol' of the clustersize
sizethreshpos=round(falseposcontrol*clustersize(position,3))-1;
sizethreshpos(sizethreshpos<1)=1;
if linkmat(position,3)<meansizethresh(sizethreshpos) && abovethresh(position)>0
    realclusters=zeros(1,gn-1);
    realclusters(position)=1;
%not below the mean size threshold and both children are single points
elseif linkmat(position,1)<gn+1 && linkmat(position,2)<gn+1
    realclusters=zeros(1,gn-1);
elseif linkmat(position,1)<gn+1 %only one child is single point
    realclusters=zeros(1,gn-1);
    realclusters=realclusters+extractrealclusters(linkmat(position,2)-gn,linkmat,gn,...
        clustersize,meansizethresh,falseposcontrol,abovethresh);
elseif linkmat(position,2)<gn+1 %only one child is single point
    realclusters=zeros(1,gn-1);
    realclusters=realclusters+extractrealclusters(linkmat(position,1)-gn,linkmat,gn,...
        clustersize,meansizethresh,falseposcontrol,abovethresh);
else, realclusters=extractrealclusters(linkmat(position,2)-gn,linkmat,gn,...
        clustersize,meansizethresh,falseposcontrol,abovethresh)+...
        extractrealclusters(linkmat(position,1)-gn,linkmat,gn,...
        clustersize,meansizethresh,falseposcontrol,abovethresh);
end

realclusters=logical(realclusters);

end