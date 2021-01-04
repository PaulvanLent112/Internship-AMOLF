 %function checking whether the merging clusters are already part of a real
%cluster above or below this merge point: 
%0-both are not part of a real cluster (check for clustering)
%1-both are part of the same cluster (check for subclustering)
%2-one branch leads to a real cluster below (do nothing)
%4-both branches lead to a real cluster below (do nothing)
%
%helper function for extractallclusters8 to extractallclusters11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function clusterflag=computingclusterflag(mergepoint,myclusterIDreal)

    clusterflag=zeros(1,2);
    if ~isnan(myclusterIDreal(1))
        if myclusterIDreal(1)<mergepoint, clusterflag(1)=2;
        else, clusterflag(1)=1;
        end
    end
    if ~isnan(myclusterIDreal(2))
        if myclusterIDreal(2)<mergepoint, clusterflag(2)=2;
        else, clusterflag(2)=1;
        end
    end
    
    clusterflag=sum(clusterflag);
end