function J=demolucymf(A,PSF,sigma,NUMIT,NUMFRM)
if numel(A)>=numel(PSF)
    siz=size(A);
else
    siz=size(PSF);
end
J=zeros(siz);
for k=1:NUMFRM
    I=imfilter(A,PSF,'circular')+normrnd(0,sigma,siz);
    J=J+deconvlucy(I,PSF,NUMIT);
end
J=J/NUMFRM;
return