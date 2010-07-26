function [J,sigma]=deconvkalmannoc(I,PSF,R,NUMIT)
%DECONVKALMAN Deconvolution based on Kalman filter.

if numel(I)>=numel(PSF)
    siz=size(I);
else
    siz=size(PSF);
end

% Prepare Kalman filtering iteration:
Z=fftn(I,siz);
H=fftn(ifftshift(PSF),siz);
sigma=zeros(NUMIT,1);
P=ones(siz);
X=fftn(I);
R_k=R/var(I(:));
% =========================================================================

% Configure physical constraints parameters:
% [XI,YI]=meshgrid(1:siz(2),1:siz(1));
J=I;
% [I_cx,I_cy]=imcentroid(I);
% [PSF_cx,PSF_cy]=imcentroid(PSF);
% J_cx=I_cx+PSF_cx;
% J_cy=I_cy+PSF_cy;
% =========================================================================
for k=1:NUMIT
    % Update gain factor:
    S=H.*P.*conj(H);
    K=P.*conj(H)./(S+R_k);

    % Update estimates:
    X=X+K.*(Z-H.*X);
    J=real(ifftn(X));
    X=fftn(J);

    % Update covariance:
    R_k=R/var(J(:));
    P=abs((1-K.*H).*P);

    % Record std:
    sigma(k)=sqrt(mean(mean((real(ifftn(X.*H))-I).^2)));
end
return