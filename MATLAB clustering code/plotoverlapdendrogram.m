%Function to show overlap between a given cluster and the cluster at the
%branch of dendrogram, same as plotoverlapdendrogram2, just with a relative
%marker size

function plotoverlapdendrogram(figurenumber,dendrogramhandles,gnCl,...
    mergcluster1corr,mergcluster2corr,mytextinput,mytextinput2,markersizefract,marker,earlystop)

mti2size=size(mytextinput2,1);
figure(figurenumber), hold on
i=1;
while i<gnCl
    myxdata=get(dendrogramhandles(i),'XData'); myxdata=myxdata(2:3);
    if myxdata(1)>myxdata(2), myxdata=myxdata([2,1]); end
    myydata=get(dendrogramhandles(i),'YData'); myydata=myydata(2);
    
    temp=0;
    [geneoverlap,~]=IdentifyGenes(mytextinput(mergcluster1corr(:,i)),mytextinput2);
    myval=sum(~isnan(geneoverlap))/size(geneoverlap,1);
    if sum(~isnan(geneoverlap))>0
        plot(myxdata(1),myydata,marker,'MarkerSize',markersizefract*myval)
    end
    if (earlystop==1) && (sum(~isnan(geneoverlap))==mti2size), temp=1; end
    
    [geneoverlap,~]=IdentifyGenes(mytextinput(mergcluster2corr(:,i)),mytextinput2);
    myval=sum(~isnan(geneoverlap))/size(geneoverlap,1);
    if sum(~isnan(geneoverlap))>0
        plot(myxdata(2),myydata,marker,'MarkerSize',markersizefract*myval)
    end
    if (earlystop==1) && (sum(~isnan(geneoverlap))==mti2size), temp=1; end
    
    i=i+1;
    if temp==1, i=gnCl; end     
end