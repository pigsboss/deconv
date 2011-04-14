function [Js,Rs]=deconvmapana(I,REF,PSF,BG,NUMIT)
%DECONVMAPANA MAP Deconvolution Analysis Routine.

J=I;
R=REF;
Js=zeros([numel(I),2*NUMIT]);
Rs=zeros([numel(I),2*NUMIT]);
IPSF=rot90(PSF,2);
H=fftn(ifftshift(PSF));
IH=fftn(ifftshift(IPSF));

for k=1:NUMIT
    J=J.*exp(real(ifftn(fftn(I./real(ifftn(H.*fftn(J)))-1).*IH)));
    R=R.*exp(real(ifftn(fftn(REF./real(ifftn(H.*fftn(R)))-1).*IH)));
    Js(:,2*(k-1)+1)=J(:);
    Rs(:,2*(k-1)+1)=R(:);
%   subplot(4,6,(k-1)*2+1),plot(J-R)
    subplot(4,5,k),plot(J-R),title(num2str(k))
    drawnow
    J=J.*double(J>=BG)+BG.*double(J<BG);
    R=R.*double(R>=BG)+BG.*double(R<BG);
    Js(:,2*(k-1)+2)=J(:);
    Rs(:,2*(k-1)+2)=R(:);
%   subplot(4,6,(k-1)*2+2),plot(J-R)
    drawnow
%   J=J*(sum(I(:))/sum(J(:)));
end
%Js=reshape(Js,[size(I),NUMIT]);
%Jregs=reshape(Js,[size(I),NUMIT]);
return