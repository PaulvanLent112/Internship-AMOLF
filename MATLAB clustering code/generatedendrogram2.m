%Function that generates a linkage matrix and tracks cluster memberships
%for clusters of a specified minimum size of members
%Input:
%- clustersize: cluster sizes at each merging (from computeclustersize)
%- gn: number of clustered elements
%- mylinkmat: hierarchical linkage matrix
%- inisize: minimum size of initial cluster
%Output:
%- generateID: classification of merging types
%- gnCl: number of initial clusters of at least "inisize"s members
%- linkmatcluster: new reduced linkage matrix
%- mergcluster1/mergcluster2: memberships of the two merging clusters


function [generateID,gnCl,linkmatcluster,mergcluster1,mergcluster2]=...
    generatedendrogram2(clustersize,gn,mylinkmat,inisize)

inisize=inisize-1;
%generateID (vector, for each merging): 0 - parent no cluster, 10 - parent cluster
%but non child (new cluster), 11 - parent cluster and only first child,
%12 - parent cluster and only second child, 13 - parent cluster and both children
generateID=(clustersize(:,3)>inisize)*10 ...
    + (clustersize(:,1)>inisize) + (clustersize(:,2)>inisize)*2;

%newclusterID (vector, for each merging): 0 - generates no cluster, number
%up to gnCl - ID of initial cluster, number larger gnCl - ID of derived cluster
newclusterID=zeros(gn-1,1);
gnCl=sum(generateID==10);
newclusterID(generateID==10)=1:gnCl; %only parent is cluster -> new cluster
newclusterID(generateID==13)=(1:sum(generateID==13))+gnCl;
addingID=find(logical((generateID==11)+(generateID==12)));
for i=1:size(addingID,1)
    if generateID(addingID(i))==11
        newclusterID(addingID(i))=newclusterID(mylinkmat(addingID(i),1)-gn);
    else
        newclusterID(addingID(i))=newclusterID(mylinkmat(addingID(i),2)-gn);
    end
end
%linkmatcluster: basically linkmat, just for the derived clusters
linkmatcluster=[newclusterID(mylinkmat(generateID==13,1)-gn),...
    newclusterID(mylinkmat(generateID==13,2)-gn),...
    mylinkmat(generateID==13,3)];


%mergcluster1,2: logical, which genes belong to both clusters at merging
tempclusters=zeros(gn,1);
mergcluster1=zeros(gn,gnCl-1); mergcluster2=zeros(gn,gnCl-1); myclustcount=0;
for i=1:gn-1
    if mylinkmat(i,1)>gn
        mytempclusters=(tempclusters==mylinkmat(i,1));
        if generateID(i)==13 %two clusters merge
            myclustcount=myclustcount+1;
            mergcluster1(:,myclustcount)=mytempclusters;
        end
        tempclusters(mytempclusters)=gn+i;
    else, tempclusters(mylinkmat(i,1))=gn+i;
    end
    if mylinkmat(i,2)>gn
        mytempclusters=(tempclusters==mylinkmat(i,2));
        if generateID(i)==13 %two clusters merge
            mergcluster2(:,myclustcount)=mytempclusters;
        end
        tempclusters(mytempclusters)=gn+i;
    else, tempclusters(mylinkmat(i,2))=gn+i;
    end
end
mergcluster1=logical(mergcluster1); mergcluster2=logical(mergcluster2);

end

