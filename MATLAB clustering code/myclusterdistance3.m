function distance=myclusterdistance3(XI,XJ)
% We define a distance based on the Fisher test, which uses the hypergeometric
% distribution to determine the probability of a random overlap.
% The probability of an overlap of t genes  between two sets of m genes and
% n genes out of a total of N genes can be mapped to the probibability to
% have x white balls in a draw of n balls, when the urn contains m white
% balls and N-m black balls. The probability density is p=hygepdf(x,N,m,n)
%(dhyper in R has a different definition: dhyper(x,m,N-m,n))
% The probability to get at least an overlap of x if n<=m is:
% sum(hygepdf(x:n,N,m,n))

N=size(XI,2);
ClN=size(XJ,1);
overlapx=sum(XI(ones(1,ClN),:).*XJ,2);
mn=sort([sum(XI,2)*ones(ClN,1),sum(XJ,2)],2);
distance=ones(ClN,1);
for i=1:ClN
    if overlapx(i)>0
        distance(i)=sum(hygepdf(overlapx(i):mn(i,1),N,mn(i,2),mn(i,1)),2);
    end
end

end