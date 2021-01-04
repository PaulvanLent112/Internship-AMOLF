%Function to perform all the clustering steps at once
%Input: data set to be analyzed
%Output: logical matrix of clusters cluster

function clusteroutput=completeclustering4(inputdata,Nref,falsecluster,...
    falsemembers)

mydistmat=acos(1-pdist(inputdata,'correlation'))/pi;

mylinkmat=linkage(mydistmat,'single');
gn=size(mylinkmat,1)+1; %number of clustered elements
%Matrix with size of parent and child clusters at each merging:
%[size child 1, size child 2, size parent]
clustersize=computeclustersize(mylinkmat,gn);
%Vector with cluster number at each merging (for the epsilon plot)
clusternumber3=computeclusternumber(clustersize,gn,3);
%Number of clusters as a function of cluster distance (former 'epsilon scan plot')
%additional peaks can come from zeros, effectively lower dimensional genes!
figure(201), clf, hold on, box on
plot(mylinkmat(:,3),clusternumber3,'k')
set(gca,'Fontsize',16,'Linewidth',2)
ylabel('Number of clusters','Fontsize',16)
xlabel('Distance','Fontsize',16)

%% All cluster sizes and the largest clusters at each epsilon
%allclustersizes: occurance of all cluster sizes for each epsilon
%[5,3,1,4,0,0,0,....] means 5 clusters of size 1, 3 clusters of size 2 and so on 
%maxelements: sizes of the four largest clusters as a function of epsilon
Nmaxsizes=3;
[~,maxelements]=trackclustersizes(clustersize,gn,Nmaxsizes);
%Plot Maximum cluster size with epsilon to see percolation transition
%below critical point: <log(gn), at critical point: ~n^2/3
%above: ~n, second largest <log(gn)
figure(202), clf, box on
semilogy(mylinkmat(:,3),maxelements(2:end,1),'r'), hold on %largest cluster
plot(mylinkmat(:,3),maxelements(2:end,2),'g') %second largest cluster
plot(mylinkmat(:,3),maxelements(2:end,3),'b')
xlim([mylinkmat(1,3),mylinkmat(end,3)])
set(gca,'Fontsize',16,'Linewidth',2)
ylabel('Size of cluster','Fontsize',16)
xlabel('Distance','Fontsize',16)

%% Extract all clusters
[~,allclusters]...
    = extractallclusters12(gn,clustersize,mylinkmat,size(inputdata,2),...
    Nref,falsecluster,falsemembers,1);
%% Extract cluster
clusteroutput=nan(gn,size(allclusters,2));
for i=1:size(allclusters,2)
    clusteroutput(:,i)=extractcluster(allclusters(i),mylinkmat,gn)';
end

end