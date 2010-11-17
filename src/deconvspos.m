function [J,sigma]=deconvspos(I,PSF,NUMIT,W)
%DECONVSPOS Simple positive Iterative deconvolution algorithm.

bg = 10;
G = fftn(I);
H = fftn(ifftshift(PSF));
if nargin < 4
    W = conj(H);
end
R = G;
C = G.*W;
F = G.*W;
sigma = sum(abs(R(:)).^2);
for k=1:NUMIT
    disp(k)
    disp(sigma)
    Rnew = R-C.*H;
    Cnew = W.*Rnew.*(2-W.*H);
    Fnew = F+Cnew;
    J = real(ifftn(Fnew));
    Jp = J.*double(J >= bg)+bg*double(J < bg);
    Jn = (J-bg).*double(J < bg);
    Fp = fftn(Jp);
    Fn = fftn(Jn);
    rho = -1.0*sum(Jn(:))/sum(Jp(:));
    sigmaNew = sum(abs(Rnew(:)).^2);
    if sigmaNew >= sigma
        break
    end
    R = Rnew;
    C = Cnew-Fn-rho*Fp;
    F = Fnew-Fn-rho*Fp;
    sigma = sigmaNew;
end
J = real(ifftn(F));
return