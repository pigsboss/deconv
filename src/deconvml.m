function J=deconvml(I,PSF,NUMIT)
IPSF=rot90(PSF,2);
OTF=fftn(ifftshift(PSF));
IOTF=fftn(ifftshift(IPSF));
J=I;
for k=1:NUMIT
    J=J.*real(ifftn(fftn(I./real(ifftn(fftn(J).*OTF))).*IOTF));
end
return