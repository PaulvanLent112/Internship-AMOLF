% Function to track largest cluster sizes
% NOTE: there was a bug in all previous versions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [allclustersizes,maxelements]=trackclustersizes(clustersize,gn,Nmaxsizes)

allclustersizes=zeros(gn,gn); allclustersizes(1,1)=gn;
maxelements=zeros(gn,Nmaxsizes); maxelements(1,:)=ones(1,Nmaxsizes);
for i=1:gn-1
    allclustersizes(i+1,:)=allclustersizes(i,:);
    allclustersizes(i+1,clustersize(i,1))=allclustersizes(i+1,clustersize(i,1))-1;
    allclustersizes(i+1,clustersize(i,2))=allclustersizes(i+1,clustersize(i,2))-1;
    allclustersizes(i+1,clustersize(i,3))=allclustersizes(i+1,clustersize(i,3))+1;
    myallclustersizes=allclustersizes(i+1,1:i+1); allpossiblesizes=(1:(i+1));
    allpossiblesizes=allpossiblesizes(myallclustersizes>0);
    myallclustersizes=myallclustersizes(myallclustersizes>0);
    temp=[]; macssize=size(myallclustersizes,2);
    for j=0:(min(macssize,Nmaxsizes)-1)
        temp=[temp(:)', allpossiblesizes((macssize-j)*...
            ones(1,myallclustersizes(macssize-j)))];
    end
    if size(temp,2)>Nmaxsizes, temp=temp(1:Nmaxsizes); end
    maxelements(i+1,1:size(temp,2))=temp;
end

end

