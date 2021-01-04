%Function to extract branches of the dendrogram by hand

function myclusterpick=pickbranchesdendrogram(figureclick,...
    linkmatcluster,dendrogramhandles,mergcluster1corr,mergcluster2corr)

flag=0;
myclusterpick=[];
while flag==0
    
    prompt = 'Extract branch? y/n [y]: ';
    str = input(prompt,'s');
    if ~isempty(str) && str ~= 'y'
        flag = 1;
    else
        figureinfo = getCursorInfo(figureclick);
        ind1=figureinfo.Position(1); ind2=figureinfo.Position(2);
        [~,mergind]=min(abs(linkmatcluster(:,3)-ind2));
        xcheck=get(dendrogramhandles(mergind),'XData'); xcheck=xcheck(2:3);
        if xcheck(1)>xcheck(2), xcheck=xcheck(2:-1:1); end
        if xcheck(1)==ind1, leftright=0; else, leftright=1; end
        if leftright==0
            mergclustercorr=mergcluster1corr(:,mergind);
        else
            mergclustercorr=mergcluster2corr(:,mergind);
        end
        myclusterpick=[myclusterpick(:,:),mergclustercorr];
    end
      
end

myclusterpick=logical(myclusterpick);

end
