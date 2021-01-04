function distance=myclusterdistance1(XI,XJ)
%distance=1-joint members/mean members

m2=size(XJ,1);
sum1=sum(XI,2);
sum2=sum(XJ,2);
distance=1-2*sum(XI(ones(1,m2),:).*XJ,2)./(sum1+sum2);

end