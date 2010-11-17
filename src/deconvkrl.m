function [J,R,P]=deconvkrl(I,PSF,NUMIT,BKG,NUMK)
NUMR=10;
z=zeros(size(I));
for k=1:NUMR
    disp(k)
    z=z+deconvrl(I,PSF,NUMIT,BKG);
end
z=z/NUMR;
R=zeros(size(I));
for k=1:NUMR
    disp(k)
    R=R+(deconvrl(I,PSF,NUMIT,BKG)-z).^2;
end
R=BKG+R/NUMR;
J=I;
P=I;
for k=1:NUMK
    disp(k)
    z=deconvrl(I,PSF,NUMIT,BKG);
    y=z-J;
    S=P+R;
    K=P./S;
    J=J+K.*y;
    P=(ones(size(I))-K).*P;
end
return