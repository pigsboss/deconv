function J=deconvmap(I,PSF,NUMIT,BG)
J=I;
IPSF=rot90(PSF,2);
H=fftn(ifftshift(PSF));
IH=fftn(ifftshift(IPSF));
for k=1:NUMIT
    J=J.*exp(real(ifftn(fftn(I./real(ifftn(H.*fftn(J)))-1).*IH)));
    J=J.*double(J>=BG)+BG.*double(J<BG);
%   J=J*(sum(I(:))/sum(J(:)));
end
return