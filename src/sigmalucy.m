function sigma=sigmalucy(I,H,NUMIT)
sigma=zeros(NUMIT,1);
for k=1:NUMIT
    J=deconvlucy(I,H,k);
    sigma(k)=sqrt(mean(mean((imfilter(J,H,'circular')-I).^2)));
end
return