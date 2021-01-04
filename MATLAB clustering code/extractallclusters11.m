%Function to extract all clusters using a local percolation criterion
%
%Output:
%-trackclustersize: ID of larger cluster at each merging point
%-allclusters: logical marking the last merging point that forms each
%cluster
%-sizethreshlog: expected largest sizes and std as a function of distance
%for different dimensions [additionally saved in supportfolder]
%
%Input:
%-gn: number of clustered elements (genes)
%-clustersize: output of computeclustersize, tracking all cluster sizes
%-mylinkmat: output of linkage, a hierarchical linkage matrix
%-mydim: number of dimensions (samples) of the data
%-Nref: number of simulated instances of the reference null model
%-falseposcontrol: trade-off between false positive and false negative
%rate, threshold for the fraction between the cluster size and the largest
%background cluster
%-output: 0-no output, >0-number of found clusters, >1-plot of
%cluster growth relative to reference
%
%
%based on extractallclusters10
%like extractallclusters7 for inhomogeneous background but without taking a
%global percolation transition into account.
%%%%%



function [trackclustersize,allclusters,sizethreshlog]...
    = extractallclusters11(gn,clustersize,mylinkmat,mydim,Nref,falseposcontrol,output)

%Make support folder if not exists 
supportfolder='PercolationReference';
if ~isfolder(supportfolder)
    mkdir(supportfolder)
end

%% Load reference for percolation transition
if isfile([supportfolder,'\PercolationDimensionReference.mat'])
    load([supportfolder,'\PercolationDimensionReference.mat'],...
        'dimpercolsmooth','mydimlist');
else
    disp('Warning: Referenze file missing.')
    disp('Assuming percolation at degree 1, not valid for D<20.')
    dimpercolsmooth=[1,1]; mydimlist=[4,100];
end
[~,inflectionpoint]=min(-diff(dimpercolsmooth)./diff(mydimlist));
dimpercolsmooth=interp1([mydimlist(1:inflectionpoint),10^12],...
    [dimpercolsmooth(1:inflectionpoint),1],4:mydim);

%% Load reference for largest cluster sizes if exist
if isfile([supportfolder,'\sizethreshlog',num2str(gn),'-',num2str(Nref),'.mat'])
    load([supportfolder,'\sizethreshlog',num2str(gn),'-',num2str(Nref),'.mat']...
        ,'sizethreshlog')
    stlsize=size(sizethreshlog,1);
    if stlsize<mydim-3
        sizethreshlog(stlsize+1:mydim-3,:)=cell(mydim-3-stlsize,2);
    end
    sizethreshlist=zeros(1,size(sizethreshlog,1)); 
    for i=1:size(sizethreshlog,1)
        sizethreshlist(i)=~isempty(sizethreshlog{i,1});
    end
else
    sizethreshlog=cell(mydim-3,2); sizethreshlist=zeros(1,mydim-3);
end

%% Track cluster origin, keep ID of largest cluster at merging
trackclustersize=nan(gn-1,1); %track cluster ID by size
icount=0;
for i=1:gn-1
    if clustersize(i,3)==2 %new cluster
        icount=icount+1;
        trackclustersize(i)=icount;
    elseif clustersize(i,1)<clustersize(i,2) %keep ID of larger cluster
        trackclustersize(i)=trackclustersize(mylinkmat(i,2)-gn);
    elseif clustersize(i,1)>clustersize(i,2) %keep ID of larger cluster
        trackclustersize(i)=trackclustersize(mylinkmat(i,1)-gn);
    %if both merging clusters have the same size, take older cluster
    elseif trackclustersize(mylinkmat(i,1)-gn)>trackclustersize(mylinkmat(i,2)-gn)
        trackclustersize(i)=trackclustersize(mylinkmat(i,2)-gn);
    else
        trackclustersize(i)=trackclustersize(mylinkmat(i,1)-gn);
    end
end


%% Actual clustering
%what is the maximum percolation treshold given the dimensionality of the data
maxpercolthreshold=fzero(@(x) gn*betainc((sin(x*pi)).^2,(mydim-2)/2,1/2)/2-...
    dimpercolsmooth(mydim-3),mylinkmat(end,3));
[~,maxpercolind]=min(abs(maxpercolthreshold-mylinkmat(:,3)));
clusterIDreal=nan(1,icount); %distance at which the cluster of this ID was found to be real
for i=maxpercolind:-1:3 %scan from last merging to first merging
    mergepoints=mylinkmat(i,1:2);
    if sum(mergepoints>gn)==2 %only if clusters merge this can be a percolation point
        %For the two branches, since when are they real clusters
        myclusterIDreal=clusterIDreal(trackclustersize(mergepoints-gn));
        %checks whether ID of both branches is associated with real cluster
        clusterflag=computingclusterflag(i,myclusterIDreal);
        if clusterflag<2 %either both are not part of any real cluster or
            %both are part of the same real cluster
            
            %check which one of the clusters survives
            if trackclustersize(mergepoints(1)-gn)==trackclustersize(i)
                smind=2; largind=1;
            else
                smind=1; largind=2;
            end
            
            if clusterflag==1
                clusterIDrealref=clusterIDreal(trackclustersize(mergepoints(largind)-gn));
            end

            %compute effective dimension
            [~,mydim]=min(abs(gn*betainc((sin(mylinkmat(i,3)*pi)).^2,...
                ((4:mydim)-2)/2,1/2)/2-dimpercolsmooth(1:mydim-3)));
            mydim=mydim+3;
            
            %simulate reference profile
            if sizethreshlist(mydim-3)==0
                [~,sizethreshmean,sizethreshstd]=computenoisesizes5(gn,mydim,Nref);
                sizethreshlog{mydim-3,1}=sizethreshmean;
                sizethreshlog{mydim-3,2}=sizethreshstd;
                sizethreshlist(mydim-3)=1;
            else
                sizethreshmean=sizethreshlog{mydim-3,1};
                sizethreshstd=sizethreshlog{mydim-3,2};
            end
            sizethresh=sizethreshmean-3*sizethreshstd; %!!!!!!!!!!!!!!!!!!!!!!!!!!
            sizethresh(1)=0; %don't use clusters of size 2
            
            %check whether smaller cluster includes one or more real clusters
            checkIDlist=(mergepoints(smind)-gn); checkflag=0;
            while ~isempty(checkIDlist)
                %Check whether smaller cluster was above size threshold before merging
                abovethresh1=zeros(1,checkIDlist(1));
                %find all instaces with the same ID
                temp=(trackclustersize(1:checkIDlist(1))==trackclustersize(checkIDlist(1)));
                mymergings=mylinkmat(1:checkIDlist(1),3);
                myclustersize=clustersize(1:checkIDlist(1),3);
                abovethresh1(temp)=(cumsum(mymergings(temp)<...
                        sizethresh(myclustersize(temp)-1))>0);
%                 if i==1976%1836%1990%1828%1858%1828%2001
%                     mydim
%                     figure(1010),  hold on
%                     semilogy(sizethreshmean-3*sizethreshstd,1:gn-1,'b'), hold on
%                     plot(sizethreshmean,1:gn-1,'b')
%                     plot(mymergings(temp),myclustersize(temp),'r')
%                     error('b')
%                 end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if sum(abovethresh1)>0 %only if the first criterion is met
                    sizethreshpos=round(falseposcontrol*myclustersize(temp))-1;
                    sizethreshpos(sizethreshpos<1)=1;
                    abovethresh2=zeros(1,checkIDlist(1));
                    abovethresh2(temp)=(mymergings(temp)<sizethreshmean(sizethreshpos));
                    realclustpos=find((abovethresh1.*abovethresh2)==1,1,'last');
                    if ~isempty(realclustpos) %there was a real cluster found
                        clusterIDreal(trackclustersize(checkIDlist(1)))=realclustpos;
                        checkflag=1;
                    %if criterion is not met, set realclusterpos to beginning
                    %(this allows to identify small real clusters on long
                    %branches (very tight) that however merge into a larger
                    %a non-real cluster before the percolation threshold)
                    else
                        realclustpos=find(trackclustersize==...
                            trackclustersize(checkIDlist(1)),1,'first');
                    end
                %if criterion is not met, set realclusterpos to beginning of cluster
                else
                    realclustpos=find(trackclustersize==...
                        trackclustersize(checkIDlist(1)),1,'first');
                end
                
                if realclustpos<checkIDlist(1) %there are branchings above
                    %find all branchings above
                    tempnew=temp.*[zeros(realclustpos,1);...
                        ones(checkIDlist(1)-realclustpos,1)];
                    if sum(tempnew)==0, error('tempnew'), end %%
                    %Find all clusters (!) branching in above
                    mymergingspos=mylinkmat(1:checkIDlist(1),1:2);
                    mymergingspos=mymergingspos(tempnew>0,:);
                    mymergingspos=mymergingspos(sum(mymergingspos>gn,2)==2,:);
                    if ~isempty(mymergingspos) %clusters branching in
                        checkIDlisttemp=reshape((mymergingspos-gn),[],1);
                        %remove branches that have been just considered
                        checkIDlisttemp(trackclustersize(checkIDlisttemp)...
                            ==trackclustersize(checkIDlist(1)))=[];
                        checkIDlist=[checkIDlist(:);checkIDlisttemp];
                        
                    end
                end
                
                checkIDlist(1)=[];
            end
            

            if checkflag==1 %small cluster was real
                %check whether larger cluster includes one or more real clusters
                checkIDlist=(mergepoints(largind)-gn);
                while ~isempty(checkIDlist)
                    %Check whether larger cluster was above size threshold before merging
                    abovethresh1=zeros(1,checkIDlist(1));
                    temp=(trackclustersize(1:checkIDlist(1))==...
                        trackclustersize(checkIDlist(1)));
                    mymergings=mylinkmat(1:checkIDlist(1),3);
                    myclustersize=clustersize(1:checkIDlist(1),3);
                    abovethresh1(temp)=(cumsum(mymergings(temp)<...
                            sizethresh(myclustersize(temp)-1))>0);
                    if sum(abovethresh1)>0 %only if the first criterion is met
                        sizethreshpos=round(falseposcontrol*myclustersize(temp))-1;
                        sizethreshpos(sizethreshpos<1)=1;
                        abovethresh2=zeros(1,checkIDlist(1));
                        abovethresh2(temp)=(mymergings(temp)<sizethreshmean(sizethreshpos));
                        realclustpos=find((abovethresh1.*abovethresh2)==1,1,'last');
                        if ~isempty(realclustpos) %there was a real cluster found
                            clusterIDreal(trackclustersize(checkIDlist(1)))=realclustpos;
                        else
                            realclustpos=find(trackclustersize==...
                                trackclustersize(checkIDlist(1)),1,'first');    
                        end
                    else
                        realclustpos=find(trackclustersize==...
                            trackclustersize(checkIDlist(1)),1,'first');
                    end
                    
                    if realclustpos<checkIDlist(1) %there are branchings above
                        %find all branchings above
                        tempnew=temp.*[zeros(realclustpos,1);...
                            ones(checkIDlist(1)-realclustpos,1)];
                        if sum(tempnew)==0, error('tempnew'), end %%
                        %Find all clusters (!) branching in above
                        mymergingspos=mylinkmat(1:checkIDlist(1),1:2);
                        mymergingspos=mymergingspos(tempnew>0,:);
                        mymergingspos=mymergingspos(sum(mymergingspos>gn,2)==2,:);
                        if ~isempty(mymergingspos) %clusters branching in
                            checkIDlisttemp=reshape((mymergingspos-gn),[],1);
                            %remove branches that have been just considered
                            checkIDlisttemp(trackclustersize(checkIDlisttemp)...
                                ==trackclustersize(checkIDlist(1)))=[];
                            checkIDlist=[checkIDlist(:);checkIDlisttemp];
                        end
                    end
                    
                    checkIDlist(1)=[];
                end  
            end
            
            %if there was a subclustering but the large cluster did not
            %change, remove the large cluster event above
            if checkflag>0 && clusterflag==1
                if clusterIDrealref==clusterIDreal(trackclustersize(mergepoints(largind)-gn))
                    clusterIDreal(trackclustersize(mergepoints(largind)-gn))=nan;
                end
            end
            
        end
    end
    
        
end
allclusters=clusterIDreal(~isnan(clusterIDreal));

%% Remove later
myclusterIDs=trackclustersize(allclusters);
myclusterIDsUnique=unique(myclusterIDs);
myclusterIDsUniquesize=size(myclusterIDsUnique,1);
if myclusterIDsUniquesize<sum(~isnan(clusterIDreal))
    error('Clusters not distinct!')
end

%Save reference for largest cluster sizes
save([supportfolder,'\sizethreshlog',num2str(gn),'-',num2str(Nref),'.mat'],'sizethreshlog')

%%
if output>0
    disp(['We have found ', num2str(sum(~isnan(clusterIDreal))),' clusters.'])
end

end