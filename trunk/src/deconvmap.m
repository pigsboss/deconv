function J=deconvmap(I,PSF,NUMIT,BG)
J=I;
%Js=zeros([numel(I),NUMIT]);
%Jregs=zeros([numel(I),NUMIT]);
IPSF=rot90(PSF,2);
H=fftn(ifftshift(PSF));
IH=fftn(ifftshift(IPSF));
for k=1:NUMIT
    J=J.*exp(real(ifftn(fftn(I./real(ifftn(H.*fftn(J)))-1).*IH)));
%    Js(:,k)=J(:);
%    subplot(4,6,(k-1)*2+1),plot(J)
%    drawnow
    J=J.*double(J>=BG)+BG.*double(J<BG);
%    Jregs(:,k)=J(:);
%    subplot(4,6,(k-1)*2+2),plot(J)
%    drawnow
%   J=J*(sum(I(:))/sum(J(:)));
end
%Js=reshape(Js,[size(I),NUMIT]);
%Jregs=reshape(Js,[size(I),NUMIT]);
return