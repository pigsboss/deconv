function [J,sigma]=deconvkalmanmf(I,PSF,R)
%DECONVKALMAN Deconvolution based on Kalman filter.

if numel(I(:,:,1))>=numel(PSF)
    siz=size(I(:,:,1));
else
    siz=size(PSF);
end
[~,~,NUMIT]=size(I);
% Prepare Kalman filtering iteration:
H=fftn(ifftshift(PSF),siz);
sigma=zeros(NUMIT,1);
P=ones(siz);
X=fftn(I(:,:,1),siz);
R_k=R/var(I(:));
% =========================================================================

% Configure physical constraints parameters:
% [XI,YI]=meshgrid(1:siz(2),1:siz(1));
J=I(:,:,1);
m_0=sum(sum(I(:,:,1)));
g=0.0;
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
    Z=fftn(I(:,:,k),siz);
    X=X+K.*(Z-H.*X);

    % Physical constraints list:
        J_new=real(ifftn(X));
        % Non-negative:
        J_new=J_new.*double(J_new>=0)+J.*double(J_new<0);
        J=g*J + (1-g)*J_new;
        % First-order moment:
%         [cx,cy]=imcentroid(J);
%         dx=cx-J_cx;
%         dy=cy-J_cy;
%         J=interp2(XI,YI,J,XI+dx,YI+dy,'linear',mean(mean(J)));
        % Zero-order moment:
        J=J*m_0/sum(sum(J));
    X=fftn(J);
        
    J=real(ifftn(X));
    X=fftn(J);

    % Update covariance:
    R_k=R/var(J(:));
    P=abs((1-K.*H).*P);

    % Record std:
    sigma(k)=sqrt(mean(mean((real(ifftn(X.*H))-I(:,:,k)).^2)));
end
return