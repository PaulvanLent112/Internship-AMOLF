%Function to determine the percolation transition based on the peaks in
%cluster sizes
%
%Output:
%- allclustersizes: occurance of all cluster sizes for each epsilon
%[5,3,1,4,0,0,0,....] means 5 clusters of size 1, 3 clusters of size 2 ...
%- maxelements: sizes of all 'Nmaxsizes' largest clusters as a function of
%distance (control parameter)
%- percolthresh: distance of the percolation transition, based on the peak
%of nth-largest cluster, which approaches a constant value for large n
%
%Input:
%- clustersize: sizes of all clusters (from computeclustersize)
%- gn: number of clustered elements
%- Nmaxsizes: number of largest clusters to determine percolation transition
%(needs to be <max(clusternumber); if visually no saturation, reduce it)
%- linkmat: hierarchical linkage matrix (from linkage)
%- plotnumber: figure number, if ==0 no control plot
%
%based on computesizepercolation
%%%%

function [allclustersizes, maxelements, percolthresh]=...
    computesizepercolation2(clustersize,gn,Nmaxsizes,linkmat,plotnumber)

%Finds all cluster sizes
%Unique cluster sizes might be slightly more robust, yet for simplicity, we
%use all cluster sizes [can be changed here]
% [allclustersizes,maxelements]=trackclustersizesuni(clustersize,gn,Nmaxsizes);
[allclustersizes,maxelements]=trackclustersizes(clustersize,gn,Nmaxsizes);

%Find (last) maximum for all cluster sizes
[~,maxindex]=max(maxelements(end:-1:2,:),[],1); maxindex=gn-maxindex;

%Check saturation and fit exponential
stepdiff=abs(diff(maxindex));
% stepdiff=abs(diff(linkmat(maxindex,3)));
meanprevstepdiff=cumsum(stepdiff(1:end-1))./(1:Nmaxsizes-2);
if sum(stepdiff(3:end)>2*meanprevstepdiff(2:end))>0
    endfit=find(stepdiff(3:end)>2*meanprevstepdiff(2:end),1)+2;
else
    endfit=Nmaxsizes;
end
f = fittype('a*exp(-b*x)+c');
fitres = fit((1:endfit)',linkmat(maxindex(1:endfit),3),f,...
    'StartPoint',[linkmat(maxindex(1),3),1,linkmat(maxindex(endfit),3)]);

%Control plot
if plotnumber>0
    figure(plotnumber), clf, box on,
    plot(1:endfit,linkmat(maxindex(1:endfit),3)','r'), hold on
    plot(endfit:Nmaxsizes,linkmat(maxindex(endfit:end),3)','*r')
    plot(1:Nmaxsizes,fitres.a*exp(-fitres.b*(1:Nmaxsizes))+fitres.c,'--k')
%     plot(2:endfit,linkmat(maxindex(2:endfit),3)','.r'), hold on
%     plot(2:Nmaxsizes,fitres.a*exp(-fitres.b*(2:Nmaxsizes))+fitres.c,'--k')
    plot([1 Nmaxsizes],[1,1]*fitres.c,'k')
    ylim([0 inf])
    set(gca,'Fontsize',16,'Linewidth',2)
    ylabel('Peak of size of largest clusters','Fontsize',16)
    xlabel('Ranking of cluster size','Fontsize',16)
end
    
percolthresh=fitres.c;

end
