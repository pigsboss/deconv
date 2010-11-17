function J=deconvrl(I,PSF,NUMIT,BKG)
IPSF=rot90(PSF,2);
OTF=fft2(ifftshift(PSF));
IOTF=fft2(ifftshift(IPSF));
% mu=3*sqrt(I);
J=poissrnd(I);
% T=0.5;
% Gamma=0.9;
% J=I;
for k=1:NUMIT
    G=real(ifft2(fft2(J).*OTF));
    J=J.*real(ifft2(fft2(I./(G+double(abs(G)<=eps))).*IOTF));
    J=J.*double(abs(G)>eps);
%    J=mu.*double(J<=(0.5+max(mu,BKG)-T))+J.*double(J>(0.5+max(mu,BKG)-T));
%    J=mu.*double(J<=(0.5+max(mu,BKG)-T))+(mu+1).*double(J>=(0.5+mu+T))+J.*double(J<(0.5+mu+T)).*double(J>(0.5+max(mu,BKG)-T));
    J=J.*double(J>=BKG)+(BKG*double(J<BKG));
%    T=T*Gamma;
end
return