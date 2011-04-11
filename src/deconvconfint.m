function [ul,ll]=deconvconfint(x,X,level)
%DECONVCONFINT Confidence interval of deconvolution estimate
% ul    - upper limit, scalar or row matrix
% ll    - lower limit, scalar or row matrix
% x     - point estimate, scalar or column vector
% X     - all samples, row vector or matrix
% level - confident level
[~,d2]=size(x);
if d2>1
    error('x must be either scalar or column vector.');
end
[d1,~]=size(X);
[d2,~]=size(x);
if d1~=d2
    error('x and X must have the same number of rows.');
end
X=sort(X,2);
ul=zeros(size(x));
ll=ul;
[d1,d2]=size(X);
for k=1:d1
    m=sum(double(X(k,:)<=x(k)));
    ll(k)=X(k,max(1,m-int32(0.5*level*d2)));
    ul(k)=X(k,min(d2,m+int32(0.5*level*d2)));
end
return