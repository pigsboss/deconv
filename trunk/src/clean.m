function [cmap,rmap,sigma,mu]=clean(dmap,psf,gain,N)
%CLEAN Implementation of CLEAN algorithm.
rmap=dmap;
cmap=zeros(size(dmap));
sigma=zeros(N,1);
mu=zeros(N,1);
for k=1:N
    beam=zeros(size(dmap));
    [~,x]=max(max(rmap,[],1));
    [~,y]=max(max(rmap,[],2));
    beam(y,x)=rmap(y,x);
    if k>1
        if rmap(y,x)<3*sigma(k-1)+mu(k-1)
            break
        end
    end
    cmap=cmap+gain*beam;
    beam=imconv(beam,psf);
    rmap=rmap-gain*beam;
    sigma(k)=std(rmap(:));
    mu(k)=mean(rmap(:));
end
return