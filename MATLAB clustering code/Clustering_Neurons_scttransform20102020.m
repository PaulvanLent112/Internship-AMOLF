% Analysis clustering paper Fig. 3
%
% like AnalyseClusteringPaperFig301 & 302, just more polished
% 
% SW, 24-3-20
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear variables;
%% ---- Load data
path2data='C:\Users\Paul\Documents\Internship AMOLF\Matlab';
set(0,'DefaultAxesColorOrder',brewermap(NaN,'Set1'))
BS1=brewermap(NaN,'Set1'); BS1(6,:)=[];
BS1ng=BS1; BS1ng(3,:)=[]; %no green
BG=brewermap(NaN,'Greens');
% set(0, 'DefaultAxesColorOrder', 'factory')
% load('GeneExpressionData/DataSpeed3.mat');
%CT_expr= csvread('Input_data/Multiple_exprs.csv',1,1);
%normdataS3=zscore(CT_expr,1,2);
%% Pick data
mydistmat= csvread('Input_Data\Neuron_clustering_scttransform20102020.csv',1,1);

%%
mydistmat(logical(eye(size(mydistmat)))) = 0;
mydistmat=squareform(mydistmat);

%% ---- constructing hierarchical cluster tree
%Matrix of hierarchical tree: [index child 1, index child 2, epsilon]
mylinkmat=linkage(mydistmat,'single');
gn=size(mylinkmat,1)+1; %number of clustered elements
%Matrix with size of parent and child clusters at each merging:
%[size child 1, size child 2, size parent]
clustersize=computeclustersize(mylinkmat,gn);
%Vector with cluster number at each merging (for the epsilon plot)
clusternumber3=computeclusternumber(clustersize,gn,3);
%Number of clusters as a function of cluster distance (former 'epsilon scan plot')
%additional peaks can come from zeros, effectively lower dimensional genes!
figure(1), clf, hold on, box on
plot(mylinkmat(:,3),clusternumber3,'k')
set(gca,'Fontsize',16,'Linewidth',2)
ylabel('Number of clusters','Fontsize',16)
xlabel('Distance','Fontsize',16)

%% All cluster sizes and the largest clusters at each epsilon
%allclustersizes: occurance of all cluster sizes for each epsilon
%[5,3,1,4,0,0,0,....] means 5 clusters of size 1, 3 clusters of size 2 and so on 
%maxelements: sizes of the largest clusters as a function of epsilon
Nmaxsizes=18;
[allclustersizes,maxelements]=trackclustersizes(clustersize,gn,Nmaxsizes);
figure(1), clf, box on
semilogy(mylinkmat(:,3),maxelements(2:end,1)/gn), hold on %largest cluster
plot(mylinkmat(:,3),maxelements(2:end,2)/gn) %second largest cluster
plot(mylinkmat(:,3),maxelements(2:end,3)/gn)
% plot(mylinkmat(:,3),maxelements(2:end,4),'m')
plot(mylinkmat(:,3),maxelements(2:end,10)/gn)
xlim([mylinkmat(1,3),mylinkmat(end,3)])
set(gca,'Fontsize',16,'Linewidth',2)
ylabel('Size of cluster','Fontsize',16)
xlabel('Distance','Fontsize',16)

%% Extract all clusters: 
[trackclustersize,allclusters]=extractallclusters13(gn,...
    clustersize,mylinkmat,52412,100,2,1,1);

%% Check for cluster, all at once
inisize=3;
[generateID,gnCl,linkmatcluster,mergcluster1,mergcluster2]=...
    generatedendrogram2(clustersize,gn,mylinkmat,inisize);
[mycf,dendrogramhandles,~,mergcluster1corr,mergcluster2corr,figureclick]...
    = plotdendrogram3(61,mylinkmat,linkmatcluster,gnCl,0,generateID,...
    mergcluster1,mergcluster2);
set(dendrogramhandles,'color',BS1ng(7,:));
symbollist={'or';'og';'ob';'om';'oc';'^r';'^g';'^b';'^m';'^c';...
    'sr';'sg';'sb';'sm';'sc';'*r';'*g';'*b';'*m';'*c';'dr';'dg';'db';'dm';'dc'};
collist={'-r';'-g';'-b';'-m';'-c';'--r';'--g';'--b';'--m';'--c';':r';':g';':b';':m';':c'};
linestylelist={'-','--',':','-.'}; colorlist=BS1ng(1:6,:); clsize=size(colorlist,1);
sizelog=nan(1,size(allclusters,2));
for i=1:size(allclusters,2)
    mycluster=extractcluster(allclusters(i),mylinkmat,gn);
    sizelog(i)=sum(mycluster);
%     mytextinputnew=mytextinput(mycluster,:);
%     plotoverlapdendrogram2(get(mycf,'Number'),dendrogramhandles,gnCl,mergcluster1corr,...
%         mergcluster2corr,mytextinput,mytextinputnew(:,1),8,symbollist{i},1)
%     hold on, plot([0, gnCl+1],[1,1]*mylinkmat(allclusters(i),3),collist{i})
    plotoverlapdendrogram4(get(mycf,'Number'),dendrogramhandles,gnCl,...
        mergcluster1corr,mergcluster2corr,mycluster',mylinkmat(allclusters(i),3),...
        inisize,colorlist(mod(i-1,clsize)+1,:),linestylelist{ceil(i/clsize)})
end
disp(sizelog)

ylim([0,0.3])
set(gca,'Fontsize',24,'TickLength',[.03 .03],'Layer','top')
ylabel('Distance d','Fontsize',24)




%% Sorting results
% mysortvec=[1,6,10,9,2,3,7,8,4,5,11,12,13];
mysortvec=[1,6,9,10,2,3,7,8,12,5,4,11,13];
allclusterssort=allclusters(mysortvec);
%% Save results
warning('off','MATLAB:xlswrite:AddSheet')
myclustersizes=nan(1,size(allclusters,2));
highcountsS3=[1:17189]'
for i=1:size(allclusters,2)
  res=extractcluster(allclusters(i),mylinkmat,gn);
  res=highcountsS3(res);
  writematrix(res,'Results/Cluster_Neurons_sctransform_2sd.xlsx','Sheet',i);
end





