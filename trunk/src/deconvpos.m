function [J,sigma]=deconvpos(I,PSF,W)
%DECONVPOS Positive Iterative deconvolution algorithm (orginally by Pruksch
%and Fleischmann).

bg = 10;
delta = 1e-6;
G = fftn(I);
H = fftn(ifftshift(PSF));
if nargin < 3
    W = conj(H);
end
%  outer algorithm
R = G;
C = G.*W;
F = G.*W;
sigma = sum(abs(R(:)-R(:).*H(:).*W(:)).^2);
n = 2^20;
inc = n/2;
Hp = double(abs(H) > delta).*H.*W + double(abs(H) <= delta);
while 1
    disp(n)
    disp(sigma)
    Rnew = R-C.*H;
    Cnew = double(abs(H) > delta).*W.*Rnew.*(1-(1-Hp).^(n+1))./Hp + double(abs(H) <= delta).*W.*Rnew*(n+1);
    Fnew = F+Cnew;
    J = real(ifftn(Fnew));
    Jp = J.*double(J >= bg)+bg*double(J < bg);
    Jn = (J-bg).*double(J < bg);
    Fp = fftn(Jp);
    Fn = fftn(Jn);
    rho = -1.0*sum(Jn(:))/sum(Jp(:));
    sigmaNew = sum(abs(Rnew(:)-(Cnew(:)-Fn(:)-rho*Fp(:)).*H(:)).^2);
    if sigmaNew < sigma
        n = n+inc;
        R = Rnew;
        C = Cnew-Fn-rho*Fp;
        F = Fnew-Fn-rho*Fp;
    else
        if inc > 1
            inc = floor(inc/2);
        end
        if n > 0
            n = floor(n/2);
        else
            break
        end
        continue
    end
    if abs(sigmaNew-sigma)/sigma < delta 
        break
    elseif sigmaNew < sigma
        sigma = sigmaNew;
    end
end
J = real(ifftn(F));
return