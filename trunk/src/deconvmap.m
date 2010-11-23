function J=deconvmap(I,PSF,NUMIT,W,BG)
% J=normrnd(I,sqrt(I));
J=I;
IPSF=rot90(PSF,2);
H=fftn(ifftshift(PSF));
IH=fftn(ifftshift(IPSF));
for k=1:NUMIT
    I_new=real(ifftn(J).*IH);
    R=I-I_new;
    R=dwt(R,W);
    J=J.*exp(real(ifftn(fftn((I_new+R)./real(ifftn(H.*fftn(J)))-1).*IH)));
    J=J.*double(J>=BG)+BG.*double(J<BG);
end
return