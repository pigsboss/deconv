function M=immomentc(p,q,I)
[X,Y]=meshgrid(1:size(I,2),1:size(I,1));
X=X-(1+size(I,2))/2;
Y=Y-(1+size(I,1))/2;
C10=immoment(1,0,I)/immoment(0,0,I);
C01=immoment(0,1,I)/immoment(0,0,I);
M=sum(sum(((X-C10).^p).*((Y-C01).^q).*I));
return