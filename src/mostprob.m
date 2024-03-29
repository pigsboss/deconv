function [x_ll,x_ul,x_mean,x_median,cts]=mostprob(X,nbins)
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
x_ll=bins(x);
x_ul=bins(x+1);
x_mean=mean(X(sum(cts(1:x-1))+1:sum(cts(1:x))));
x_median=median(X(sum(cts(1:x-1))+1:sum(cts(1:x))));
return