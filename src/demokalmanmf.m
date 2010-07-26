function J=demokalmanmf(A,PSF,sigma,NUMIT)
%DECONVKALMAN Deconvolution based on Kalman filter.

if numel(A)>=numel(PSF)
    siz=size(A);
else
    siz=size(PSF);
end

% Prepare Kalman filtering iteration:
H=fftn(ifftshift(PSF),siz);
I=abs(normrnd(0,sigma,siz)+ifftn(fftn(A,siz).*H));
J=I;
X=fftn(J);
P=ones(siz);
% =========================================================================

% Configure physical constraints parameters:
g=0;
% [XI,YI]=meshgrid(1:siz(2),1:siz(1));
% [I_cx,I_cy]=imcentroid(I);
% [PSF_cx,PSF_cy]=imcentroid(PSF);
% J_cx=I_cx+PSF_cx;
% J_cy=I_cy+PSF_cy;
% =========================================================================
for k=1:NUMIT
    % Measure:
    I=abs(normrnd(0,sigma,siz)+ifftn(fftn(A,siz).*H));
    Z=fftn(I);
    m00=sum(sum(I));

    % Update gain factor:
    R=sigma/var(I(:));
    S=H.*P.*conj(H);
    K=P.*conj(H)./(S+R);

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
        J=J*m00/sum(sum(J));
    X=fftn(J);

    % Update covariance:
    P=abs((1-K.*H).*P);
end
return