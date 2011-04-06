function J=deconvintbkg(I,PSF,R,NUMIT,NSIG)
% DECONVINTBKG Deconvolution regularized by interval background
% constraints (use interval estimation instead of point estimation).
J=I;
IPSF=rot90(PSF,2);
H=fftn(ifftshift(PSF));
IH=fftn(ifftshift(IPSF));
BG=mean(R(:));
sig=std(R(:));
for k=1:NUMIT
    J=J.*exp(real(ifftn(fftn((I)./real(ifftn(H.*fftn(J)))-1).*IH)));
    J=J.*double(J>=(BG+NSIG*sig))+BG.*double(J<(BG+NSIG*sig));
end
return