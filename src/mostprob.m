function [x,cts]=mostprob(X,nbins)
X=sort(X(:));
w=(max(X)-min(X))/nbins;
bins=min(X)+(0:nbins)*w;
bins(nbins+1)=max(X)+1;
cts=zeros(size(bins));
for k=1:numel(X)
    for n=1:nbins
        if X(k)>=bins(n) && X(k)<bins(n+1)
            cts(n)=cts(n)+1;
            break
        end
    end
end
[~,x]=max(cts);
bins(nbins+1)=max(X);
x=0.5*(bins(x)+bins(x+1));
return