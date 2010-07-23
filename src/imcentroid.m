function [cx,cy]=imcentroid(I)
[X,Y]=meshgrid(1:size(I,2),1:size(I,1));
X=X-(1+size(I,2))/2;
Y=Y-(1+size(I,1))/2;
M00=sum(sum(I));
cx=sum(sum(X.*I))/M00;
cy=sum(sum(Y.*I))/M00;
if nargout<2
    cx=[cx,cy];
end
return