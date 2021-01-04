%Function to automatically extract all clusters using a global percolation
%threshold
%
%Output:
%-trackclustersize: cluster ID of largest cluster at each merging
%-percolind: index of last merging before the percolation transition
%-allclusters: logical marking the last merging point that forms each
%cluster (when cluster is completed)
%
%Input:
%-gn: number of clustered elements (genes)
%-clustersize: output of computeclustersize, tracking all cluster sizes
%-mylinkmat: output of linkage, a hierarchical linkage matrix
%-percolthresh: a percolationthreshold (output of computesizepercolation2)
%-meansizethresh: expected largest cluster size for the reference
%-stdsizethresh: standard deviation about the expected largerst cluster
%size for the reference
%-falseposcontrol: trade-off between false positive and false negative
%rate, threshold for the fraction between the cluster size and the largest
%background cluster
%-output: 0-no output, >0-number of found clusters, >1-plot of
%cluster growth relative to reference
%
%based on extractallclusters
%
%%%%%

function [trackclustersize,percolind,allclusters]...
    = extractallclusters6(gn,clustersize,mylinkmat,percolthresh,...
    meansizethresh,stdsizethresh,falseposcontrol,output)

% last merging before percolation
percolind=find(mylinkmat(:,3)>percolthresh,1,'first')-1;
%number of pairs originated till percolation:
maxClnum=sum(clustersize(1:percolind,3)==2);

%% Track cluster origin, keep ID of largest cluster at merging
trackclustersize=nan(gn-1,1);
clustersatpercolation=zeros(1,maxClnum);
icount=0;
for i=1:percolind
    if clustersize(i,3)==2 %new cluster
        icount=icount+1;
        trackclustersize(i)=icount;
        clustersatpercolation(icount)=1;
    elseif clustersize(i,1)<clustersize(i,2) %keep ID of larger cluster
        trackclustersize(i)=trackclustersize(mylinkmat(i,2)-gn);
        if mylinkmat(i,1)>gn %if smaller one was a cluster, delete from clusteratpercolation list
            clustersatpercolation(trackclustersize(mylinkmat(i,1)-gn))=0;
        end
    elseif clustersize(i,1)>clustersize(i,2) %keep ID of larger cluster
        trackclustersize(i)=trackclustersize(mylinkmat(i,1)-gn);
        if mylinkmat(i,2)>gn %if smaller one was a cluster, delete from clusteratpercolation list
            clustersatpercolation(trackclustersize(mylinkmat(i,2)-gn))=0;
        end
    %if both merging clusters have the same size, take older cluster [new]
    elseif trackclustersize(mylinkmat(i,1)-gn)>trackclustersize(mylinkmat(i,2)-gn)
        trackclustersize(i)=trackclustersize(mylinkmat(i,2)-gn);
        clustersatpercolation(trackclustersize(mylinkmat(i,1)-gn))=0;
    else, trackclustersize(i)=trackclustersize(mylinkmat(i,1)-gn);
        clustersatpercolation(trackclustersize(mylinkmat(i,2)-gn))=0;
    end
end
capnumber=sum(clustersatpercolation);
% disp([num2str(capnumber),' clusters at percolation.'])

%Check for all the clusters that were independently present at some point
%whether they were above the size threshold before each merging
abovethresh=zeros(1,gn-1);
sizethresh=meansizethresh-3*stdsizethresh;
sizethresh(1)=0; %don't use clusters of size 2
for i=1:icount
    temp=(trackclustersize==i);
    abovethresh(temp)=(cumsum(mylinkmat(temp,3)<...
            sizethresh(clustersize(temp,3)-1))>0);
end

%Consider each branch that exists at the percolation independently
capind=find(clustersatpercolation);
realclusterstart=zeros(1,gn-1);
for i=1:capnumber
    branchstart=find(trackclustersize==capind(i),1,'last');
    realclusterstart=realclusterstart+...
        extractrealclusters(branchstart,mylinkmat,gn,...
                    clustersize,meansizethresh,falseposcontrol,abovethresh);
end
allclusters=find(realclusterstart);

if output>0
    disp(['We have found ', num2str(sum(realclusterstart)), ' clusters.'])
end

if output>1
    figure(output), clf
    sizepercol=find(meansizethresh<percolthresh,1,'last');
    plot(sizethresh(1:sizepercol),(2:sizepercol+1)'), hold on
    plot(meansizethresh(1:sizepercol),(2:sizepercol+1)')
    for i=1:sum(realclusterstart)
        temp=(trackclustersize==trackclustersize(allclusters(i)));
        plot(mylinkmat(temp,3),clustersize(temp,3))
    end
end

end



