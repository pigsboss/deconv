function [cmap,rmap,sigma]=clean(dmap,psf,gain,N)
rmap=dmap;
cmap=zeros(size(dmap));
sigma=zeros(N,1);
for k=1:N
    beam=zeros(size(dmap));
    [~,x]=max(max(rmap,[],1));
    [~,y]=max(max(rmap,[],2));
    beam(y,x)=rmap(y,x);
    cmap=cmap+gain*beam;
    beam=imfilter(beam,psf,'circular');
    rmap=rmap-gain*beam;
    sigma(k)=std(rmap);
end
return