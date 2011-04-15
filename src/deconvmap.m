function J=deconvmap(I,PSF,NUMIT,BG)
%DECONVMAP Deconvolution by maximum a posteriori estimate.
J=I;
IPSF=rot90(PSF,2);
H=fftn(ifftshift(PSF));
IH=fftn(ifftshift(IPSF));
tic
for k=1:NUMIT
    J=J.*exp(real(ifft2(fft2(I./real(ifft2(H.*fft2(J)))-1).*IH)));
    J=J.*double(J>=BG)+BG.*double(J<BG);
end
toc
return