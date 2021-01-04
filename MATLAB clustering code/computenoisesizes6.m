%Function to compute the behavior of the average size + std of the second
%largest cluster in noise data as a threshold for the clustering
%
%like computenoisesizes5 but allowing for restriction to cap
%

function [sizethresh,sizethreshmean,sizethreshstd]=computenoisesizes6(...
    datasize1,datasize2,repeats,capfrac)

sizethresh=nan(datasize1-1,repeats);
for j=1:repeats
    %% - Generate random variables
    if capfrac==1
        myrand=randn(datasize1,datasize2);
%         myrandnorm=sqrt(sum(myrand.^2,2));
%         myrand=myrand./myrandnorm(:,ones(1,datasize2));
        myrand=zscore(myrand,1,2);
    elseif capfrac>1
        myrand=randn(datasize1*capfrac,datasize2);
        myrand=zscore(myrand,1,2);
        dist2ref=acos(1-pdist2(myrand,myrand(1,:),'correlation'))/pi;
        [~,res]=sort(dist2ref); myrand=myrand(res,:);
        myrand=myrand(1:datasize1,:);
    end

    mydistmat=acos(1-pdist(myrand,'correlation'))/pi;
    %Matrix of hierarchical tree: [index child 1, index child 2, epsilon]
    mylinkmat=linkage(mydistmat,'single');
    gn=size(mylinkmat,1)+1; %number of clustered elements
    %Matrix with size of parent and child clusters at each merging:
    %[size child 1, size child 2, size parent]
    clustersize=computeclustersize(mylinkmat,gn);
    
    %extracting maximum cluster size
%     maxelements=nan(gn,2); maxelements(1,:)=[1,0]; %for keeping all
%     for i=1:gn-1
%         if clustersize(i,3)>maxelements(i,1)
%             maxelements(i+1,:)=[clustersize(i,3),mylinkmat(i,3)];
%         else maxelements(i+1,:)=maxelements(i,:);
%         end
%     end
    maxelements=nan(gn-1,2); maxelement=1; %for only when it changes
    for i=1:gn-1
        if clustersize(i,3)>maxelement
            maxelement=clustersize(i,3);
            maxelements(i,:)=[clustersize(i,3),mylinkmat(i,3)];
        end
    end
    maxelements=maxelements(~isnan(maxelements(:,1)),:);

    sizethresh(:,j)=interp1(maxelements(:,1),maxelements(:,2),...
        (2:datasize1)','linear');
%     sizethresh(:,j)=interp1(maxelements(:,1),maxelements(:,2),...
%         (2:datasize1)','next');
    if mod(j,10)==0, disp(j); end
end

sizethreshmean=mean(sizethresh,2);
sizethreshstd=std(sizethresh,0,2);

end