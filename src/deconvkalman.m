function [J,sigma]=deconvkalman(I,PSF,R,NUMIT)
%DECONVKALMAN Deconvolution based on Kalman filter.

delta=1e-6;
if numel(I)>=numel(PSF)
    siz=size(I);
else
    siz=size(PSF);
end
Z=fftn(I,siz);
H=fftn(ifftshift(PSF),siz);
sigma=zeros(NUMIT,1);
P=ones(siz);
X=zeros(siz);
J=zeros(siz);
for k=1:NUMIT
    S=H.*P.*conj(H);
    K=P.*conj(H)./(S.*double(abs(S)>delta)+double(abs(S)<=delta)+R);
    X=X+K.*(Z-H.*X);
    J_new=real(ifftn(X));
    J=J_new.*double(J_new>=0)+J.*double(J_new<0);
    X=fftn(J);
    P=(1-K.*H).*P;
    sigma(k)=sqrt(mean(mean(abs(ifftn(X.*H)-I))));
end
return