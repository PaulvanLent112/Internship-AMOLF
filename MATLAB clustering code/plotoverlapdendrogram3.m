%Function to show overlap between a given cluster and the cluster at the
%branch of dendrogram, same as plotoverlapdendrogram2, just using logical
%input instead of text comparison

function plotoverlapdendrogram3(figurenumber,dendrogramhandles,gnCl,...
    mergcluster1corr,mergcluster2corr,logicalclusterind,markersize,marker,earlystop)

mti2size=sum(logicalclusterind);
figure(figurenumber), hold on
i=1;
while i<gnCl
    myxdata=get(dendrogramhandles(i),'XData'); myxdata=myxdata(2:3);
    if myxdata(1)>myxdata(2), myxdata=myxdata([2,1]); end
    myydata=get(dendrogramhandles(i),'YData'); myydata=myydata(2);
    
    temp=0;
    geneoverlap=sum(mergcluster1corr(:,i).*logicalclusterind);
    if geneoverlap>0
        plot(myxdata(1),myydata,marker,'MarkerSize',markersize)
    end
    if (earlystop==1) && (geneoverlap==mti2size), temp=1; end
    
    geneoverlap=sum(mergcluster2corr(:,i).*logicalclusterind);
    if geneoverlap>0
        plot(myxdata(2),myydata,marker,'MarkerSize',markersize)
    end
    if (earlystop==1) && (geneoverlap==mti2size), temp=1; end
    
    i=i+1;
    if temp==1, i=gnCl; end     
end