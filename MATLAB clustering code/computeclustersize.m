%Function to compute the sizes of parent and child clusters at each merging
%Output: [size child 1, size child 2, size parent]
%Input: output from the linkage function and gn=size(mylinkmat,1)+1

function clustersize=computeclustersize(mylinkmat,gn)

clustersize=nan(gn-1,3);
for i=1:gn-1
    mytest=(mylinkmat(i,1)>gn)+2*(mylinkmat(i,2)>gn);
    if mytest==0
        clustersize(i,:)=[1,1,2];
    elseif mytest==1
        clustersize(i,:)=[clustersize(mylinkmat(i,1)-gn,3),1,...
            clustersize(mylinkmat(i,1)-gn,3)+1];
    elseif mytest==2
        clustersize(i,:)=[1,clustersize(mylinkmat(i,2)-gn,3),...
            clustersize(mylinkmat(i,2)-gn,3)+1];
    elseif mytest==3
        clustersize(i,:)=[clustersize(mylinkmat(i,1)-gn,3),...
            clustersize(mylinkmat(i,2)-gn,3),...
            clustersize(mylinkmat(i,1)-gn,3)+clustersize(mylinkmat(i,2)-gn,3)];
    end
end

end