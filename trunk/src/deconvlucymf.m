function J=deconvlucymf(I,PSF,NUMIT)
[~,~,NUMFRM]=size(I);
siz=size(I(:,:,1));
J=zeros(siz);
for k=1:NUMFRM
    J=J+deconvlucy(I(:,:,k),PSF,NUMIT);
end
J=J/NUMFRM;
return