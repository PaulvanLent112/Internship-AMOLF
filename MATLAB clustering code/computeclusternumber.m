%Function computing the number of clusters at each merging
%'minsize' specifies the minimum size of a cluster

function clusternumber=computeclusternumber(clustersize,gn,minsize)

clusternumber=nan(gn,1); clusternumber(1)=0;
for i=1:gn-1
    if clustersize(i,1)>minsize-1 && clustersize(i,2)>minsize-1 %two clusters merging
        clusternumber(i+1)=clusternumber(i)-1; %cluster number decreases
    elseif clustersize(i,1)>minsize-1 && clustersize(i,2)<minsize %existing cluster grows (first child)
        clusternumber(i+1)=clusternumber(i); %nothing changes
    elseif clustersize(i,1)<minsize && clustersize(i,2)>minsize-1 %cluster grows (second child)
        clusternumber(i+1)=clusternumber(i); %nothing change
    elseif clustersize(i,3)>minsize-1 %two children merge to form new cluster
        clusternumber(i+1)=clusternumber(i)+1; %cluster number increases
    else clusternumber(i+1)=clusternumber(i); %mergings not generating cluster
    end
end
clusternumber = clusternumber(2:end);

end