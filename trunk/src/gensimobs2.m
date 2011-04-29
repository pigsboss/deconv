function [I,O,A,x,y]=gensimobs2(num_src,flux_bkg,flux_src_mean,mag_std,loc_std,exposure,siz)
O=flux_bkg*ones(siz);
x=mod(round(normrnd(0.5*(siz(2)+1),loc_std,[num_src,1])),siz(2)-1)+1;
y=mod(round(normrnd(0.5*(siz(1)+1),loc_std,[num_src,1])),siz(1)-1)+1;
mag_src_mean=log10(flux_src_mean)*-2.5;
A=(10*ones([num_src,1])).^(normrnd(mag_src_mean,mag_std,[num_src,1])/(-2.5));
for k=1:num_src
    O(y(k),x(k))=A(k);
end
I=normrnd(exposure*O,sqrt(exposure*O));
return