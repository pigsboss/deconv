function [O_est,J]=gaussianml(I,PSF,gamma,NUMIT)
IPSF=rot90(PSF,2);
OTF=fftn(ifftshift(PSF));
IOTF=fftn(ifftshift(IPSF));
O_est=I;
J=zeros(1,NUMIT);
for k=1:NUMIT
    G=fftn(O_est);
    J(k)=sum(sum(abs(I-real(ifftn(OTF.*G))).^2));
    O_est=O_est+gamma*real(ifftn(IOTF.*fftn(I-real(ifftn(OTF.*G)))));
    O_est=O_est.*double(O_est>0);
    O_est=O_est/sum(sum(O_est))*sum(sum(I));
end
return