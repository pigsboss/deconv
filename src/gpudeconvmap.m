function J=gpudeconvmap(I,PSF,NUMIT,BG)
J=I;
U=ones(size(I));
IPSF=rot90(PSF,2);
H=fftn(ifftshift(PSF));
IH=fftn(ifftshift(IPSF));
gJ=gpuArray(complex(J));
gI=gpuArray(complex(I));
gU=gpuArray(complex(U));
gH=gpuArray(H);
gIH=gpuArray(IH);
gBG=gpuArray(BG);
tic
for k=1:NUMIT
    gJ=gJ.*exp(real(ifft2(fft2(gI./real(ifft2(gH.*fft2(gJ)))-gU).*gIH)));
    gJ=gJ.*ge(gJ,gBG)+gBG.*lt(gJ,gBG);
end
toc
J=real(gather(gJ));
return