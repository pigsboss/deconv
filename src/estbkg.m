function [mu,sigma,F]=estbkg(X,NUMIT)
mu=zeros(NUMIT,1);
sigma=zeros(NUMIT,1);
X=X(:);
F=ones(size(X));
for k=1:NUMIT
    mu(k)=sum(F.*X)/sum(F);
    sigma(k)=sqrt(sum(((X-mu(k)).*F).^2)/sum(F));
    F=double(X<=(mu(k)+3*sigma(k)));
end
return