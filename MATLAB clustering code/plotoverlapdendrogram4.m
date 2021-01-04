%Function to color the dendrogram based on cluster overlap
%based on plotoverlapdendrogram3
%Needs output from plotdendrogram3 and generatedendrogam2 as input
%
%Input:
%-figurenumber: number of figure
%-dendrogramhandles: handles of the dendrogram
%-gnCl: (reduced) number of clustered elements
%-mergcluster1corr/mergcluster2corr: corrected mergcluster1/mergcluster2
%-logicalclusterind: clustermembership as logical
%-clusterstart: last merging event of the cluster
%-linecolor: color of the line
%-linestyle: line style (e.g. solid, dashed)
%%%%%%%%%%%%%%%%%%%%%%%


function plotoverlapdendrogram4(figurenumber,dendrogramhandles,gnCl,...
    mergcluster1corr,mergcluster2corr,logicalclusterind,clusterend,...
    minclsize,linecolor,linestyle)

mti2size=sum(logicalclusterind);
figure(figurenumber), hold on
i=1; %go through all merging points
while i<gnCl
    
    temp=0;
    geneoverlap=sum(mergcluster1corr(:,i).*logicalclusterind);
    if (geneoverlap>0)&&(geneoverlap<mti2size)
        set(dendrogramhandles(i),'color',linecolor)
        set(dendrogramhandles(i),'LineStyle',linestyle)
    elseif geneoverlap==mti2size
        myxdata=get(dendrogramhandles(i),'XData');
        myydata=get(dendrogramhandles(i),'YData');
        if  myxdata(1)<myxdata(4)
            if myydata(1)<clusterend
                plot([1,1]*myxdata(1),[myydata(1),clusterend],...
                    'color',linecolor,'LineStyle',linestyle)
                myydata(1)=clusterend;
                set(dendrogramhandles(i),'YData',myydata)
            elseif (myydata(1)==clusterend) && ...
                    (sum(mergcluster1corr(:,i))<2*minclsize)
                plot(myxdata(1),myydata(1),'o','color',linecolor)
            end
        else
            if myydata(4)<clusterend
                plot([1,1]*myxdata(4),[myydata(4),clusterend],...
                    'color',linecolor,'LineStyle',linestyle)
                myydata(4)=clusterend;
                set(dendrogramhandles(i),'YData',myydata)
            elseif (myydata(4)==clusterend) &&...
                    (sum(mergcluster1corr(:,i))<2*minclsize)
                plot(myxdata(4),myydata(4),'o','color',linecolor)
            end
        end
        temp=1;
    end
    
    geneoverlap=sum(mergcluster2corr(:,i).*logicalclusterind);
    if (geneoverlap>0)&&(geneoverlap<mti2size)
        set(dendrogramhandles(i),'color',linecolor)
        set(dendrogramhandles(i),'LineStyle',linestyle)
    elseif geneoverlap==mti2size
        myxdata=get(dendrogramhandles(i),'XData');
        myydata=get(dendrogramhandles(i),'YData');
        if  myxdata(4)<myxdata(1)
            if myydata(1)<clusterend
                plot([1,1]*myxdata(1),[myydata(1),clusterend],...
                    'color',linecolor,'LineStyle',linestyle)
                myydata(1)=clusterend;
                set(dendrogramhandles(i),'YData',myydata)
            elseif (myydata(1)==clusterend) &&...
                    (sum(mergcluster2corr(:,i))<2*minclsize)
                plot(myxdata(1),myydata(1),'o','color',linecolor)
            end
        else
            if myydata(4)<clusterend
                plot([1,1]*myxdata(4),[myydata(4),clusterend],...
                    'color',linecolor,'LineStyle',linestyle)
                myydata(4)=clusterend;
                set(dendrogramhandles(i),'YData',myydata)
            elseif (myydata(4)==clusterend) &&...
                    (sum(mergcluster2corr(:,i))<2*minclsize)
                plot(myxdata(4),myydata(4),'o','color',linecolor)
            end
        end
        temp=1;
    end
    
    i=i+1;
    if temp==1, i=gnCl; end     
end