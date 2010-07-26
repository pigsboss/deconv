function I=multframe(A,H,sigma,NUMFRM)
if numel(A)>=numel(H)
    siz=size(A);
else
    siz=size(H);
end
I=zeros(siz(1),siz(2),NUMFRM);
for k=1:NUMFRM
    I(:,:,k)=abs(normrnd(0,sigma,siz)+ifftn(fftn(A,siz).*fftn(ifftshift(H),siz)));
end
return