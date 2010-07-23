function M=immoment(p,q,I)
[X,Y]=meshgrid(1:size(I,2),1:size(I,1));
X=X-(1+size(I,2))/2;
Y=Y-(1+size(I,1))/2;
M=sum(sum((X.^p).*(Y.^q).*I));
return