function J=deconvml(I,PSF,NUMIT)
IPSF=rot90(PSF,2);
OTF=fftn(ifftshift(PSF));
IOTF=fftn(ifftshift(IPSF));
J=I;
for k=1:NUMIT
    G=fftn(J);
    J=J.*real(ifftn(fftn(I./real(ifftn(G.*OTF))).*IOTF));
    J=J.*double(J>1)+double(J<=1);
    J=J/sum(sum(J))*sum(sum(I));
end
return