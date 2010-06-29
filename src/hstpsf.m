function I=hstpsf(r,N)
% r: arcsec per pixel
% N: size(I)=N x N
% author: HUO

r_0=0.0686;
alpha=3;
beta=3.69;
[x,y]=meshgrid(1:N);
x=x-0.5*(1+N);
y=y-0.5*(1+N);
x=x*r;
y=y*r;
R=(x.^2+y.^2).^0.5;
I=(1+(R/r_0).^alpha).^(-beta/alpha);
I=I/sum(sum(I));
return