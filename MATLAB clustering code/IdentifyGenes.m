function [matchlookup12,matchlookup21]=IdentifyGenes(data1,data2)
%this function matches the gene names in two data sets
%before: data1=JN

ds1=size(data1,1);
ds2=size(data2,1);
matchlookup12=nan(ds1,1); matchlookup21=nan(ds2,1);
textinputtemp=data1; numtemp=1:ds1;
for i=1:ds2
    j=0; ds1temp=size(numtemp,2);
    while j<ds1temp
        j=j+1;
        if strcmp(textinputtemp{j},data2{i})
            matchlookup12(numtemp(j))=i;
            matchlookup21(i)=numtemp(j);
            textinputtemp(j)=[]; numtemp(j)=[];
            j=ds1temp;
        end
    end
end
disp([num2str(ds1),'-',num2str(ds2),'->',num2str(sum(~isnan(matchlookup21)))])