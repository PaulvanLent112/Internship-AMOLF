%Function to plot dendrogram for clusters with a given minimum size,
%needs output from generatedendrogram2 as input.
%
%Input:
%- plotnum: figure number
%- mylinkmat: original (full) hierarchical linkage matrix
%- linkmatcluster: new reduced linkage matrix
%- gnCl: (reduced) number of clustered elements
%- percolthresh: percolation threshol (0 if none)
%- generateID: classification of types of merging events
%- mergcluster1/mergcluster2: memberships of the two merging cluster
%
%Output:
%- myfc: figure handle
%- dendrogramhandles: handles for lines of dedrogram
%- iniclustthres: distance at which cluster comes into existence
%- mergcluster1corr/mergcluster2corr: corrected mergcluster1/mergcluster2
%such that mergcluster1corr is always the left branch in dendrogram
%- figureclick: handle of data cursor, for reading out clicked point  


function [mycf,dendrogramhandles,iniclustthres,mergcluster1corr,...
    mergcluster2corr,figureclick] = plotdendrogram3(plotnum,mylinkmat,...
    linkmatcluster,gnCl,percolthresh,generateID,mergcluster1,mergcluster2)

%handles are sorted by distance value
mycf = figure(plotnum); clf,
[dendrogramhandles,~,clusterdefaultpos]=dendrogram(linkmatcluster,gnCl);
if percolthresh(1)>0
    for i=1:length(percolthresh)
        hold on, plot([0, gnCl+1],[1,1]*percolthresh(i),'r--')
    end
end
set(gca,'Xtick',[]), box on
set(gca,'Fontsize',16,'Linewidth',2)
ylabel('Distance','Fontsize',16)
ylim([0 inf])

%cut tree, such that branches start at initial cluster size
iniclustthres=mylinkmat(generateID==10,3);
for i=1:gnCl-1
    cutnumber=sum(linkmatcluster(i,1:2)<(gnCl+1)); %one or two cuts
    if cutnumber==1
        myydata=get(dendrogramhandles(i),'YData');
        cutpos=(myydata==0);
        myxdata=get(dendrogramhandles(i),'XData');
        myydata(cutpos)=iniclustthres(clusterdefaultpos(myxdata(cutpos)));
        set(dendrogramhandles(i),'YData',myydata)
    elseif cutnumber==2
        myydata=get(dendrogramhandles(i),'YData');
        cutpos=logical([1,0,0,1]);
        myxdata=get(dendrogramhandles(i),'XData');
        myydata(cutpos)=iniclustthres(clusterdefaultpos(myxdata(cutpos)));
        set(dendrogramhandles(i),'YData',myydata)
    end
end

%make sure mergcluster1corr is left in dendrogram of mergecluster2corr
sortvector=nan(gnCl,1);
sortvector(clusterdefaultpos)=(1:gnCl);
linkmatclustertrace=linkmatcluster(:,1:2);
for i=1:gnCl-1
    if linkmatclustertrace(i,1)>gnCl
        linkmatclustertrace(i,1)=linkmatclustertrace(linkmatclustertrace(i,1)-gnCl,1);
    end
    if linkmatclustertrace(i,2)>gnCl
        linkmatclustertrace(i,2)=linkmatclustertrace(linkmatclustertrace(i,2)-gnCl,1);
    end
end
postrace=sortvector(linkmatclustertrace); posswitch=postrace(:,1)>postrace(:,2);
mergcluster1corr=mergcluster1; mergcluster2corr=mergcluster2;
mergcluster1corr(:,posswitch)=mergcluster2(:,posswitch);
mergcluster2corr(:,posswitch)=mergcluster1(:,posswitch);
figureclick = datacursormode(mycf);

end