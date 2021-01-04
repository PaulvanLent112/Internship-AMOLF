function distance=myclusterdistance2(XI,XJ)
%distance=1-joint members/size of smaller group

m2=size(XJ,1);
sum1=sum(XI,2);
sum2=sum(XJ,2);
distance=1-sum(XI(ones(1,m2),:).*XJ,2)./min(sum1,sum2);

end