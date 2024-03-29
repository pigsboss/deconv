function J=deconvmap(I,PSF,BG,NUMIT)
%DECONVMAP Deconvolution by maximum a posteriori estimate.
J=I;
IPSF=rot90(PSF,2);
H=fftn(ifftshift(PSF));
IH=fftn(ifftshift(IPSF));
%sigma=zeros(NUMIT,1);
%tic
for k=1:NUMIT
    J=J.*exp(real(ifft2(fft2(I./real(ifft2(H.*fft2(J)))-1).*IH)));
    J=J.*double(J>=BG)+BG.*double(J<BG);
%    sigma(k)=std(I(:)-reshape(imconv(J,PSF),numel(I),[]));
end
%toc
return