function [J,sigma]=deconvkalman(I,PSF,R,NUMIT)
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
m_0=sum(sum(I));
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

    % Update covariance:
    R_k=R/var(J(:));
    P=abs((1-K.*H).*P);

    % Record std:
    sigma(k)=sqrt(mean(mean((real(ifftn(X.*H))-I).^2)));
end
return