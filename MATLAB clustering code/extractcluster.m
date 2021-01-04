%Function to extract all cluster members at a specific merging position

function clustermembers=extractcluster(position,linkmat,gn)

if linkmat(position,1)<gn+1 && linkmat(position,2)<gn+1
    clustermembers=zeros(1,gn);
    clustermembers(linkmat(position,1))=1;
    clustermembers(linkmat(position,2))=1;
elseif linkmat(position,1)<gn+1
    clustermembers=zeros(1,gn);
    clustermembers(linkmat(position,1))=1;
    clustermembers=clustermembers+...
        extractcluster(linkmat(position,2)-gn,linkmat,gn);
elseif linkmat(position,2)<gn+1
    clustermembers=zeros(1,gn);
    clustermembers(linkmat(position,2))=1;
    clustermembers=clustermembers+...
        extractcluster(linkmat(position,1)-gn,linkmat,gn);
else, clustermembers=extractcluster(linkmat(position,1)-gn,linkmat,gn)+...
        extractcluster(linkmat(position,2)-gn,linkmat,gn);
end

clustermembers=logical(clustermembers);

end