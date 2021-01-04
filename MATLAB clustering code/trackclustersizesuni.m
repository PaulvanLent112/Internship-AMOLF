% Function to track largest unique cluster sizes
% NOTE: this is a little more robust to determine the percolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [allclustersizes,maxelements]=trackclustersizesuni(clustersize,gn,Nmaxsizes)

allclustersizes=zeros(gn,gn); allclustersizes(1,1)=gn;
maxelements=zeros(gn,Nmaxsizes); maxelements(1,Nmaxsizes)=1;

currentlylargest=1;
for i=1:gn-1
    allclustersizes(i+1,:)=allclustersizes(i,:);
    allclustersizes(i+1,clustersize(i,1))=allclustersizes(i+1,clustersize(i,1))-1;
    allclustersizes(i+1,clustersize(i,2))=allclustersizes(i+1,clustersize(i,2))-1;
    allclustersizes(i+1,clustersize(i,3))=allclustersizes(i+1,clustersize(i,3))+1;
    
    if currentlylargest<clustersize(i,3), currentlylargest=clustersize(i,3); end
    [~,temp]=find(allclustersizes(i+1,1:currentlylargest)>0,Nmaxsizes,'last');
%     [~,temp]=find(allclustersizes(i+1,1:i+1)>0,Nmaxsizes,'last');

    maxelements(i+1,Nmaxsizes-size(temp,2)+1:Nmaxsizes)=temp;
end

maxelements=maxelements(:,end:-1:1);

