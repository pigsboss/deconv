function c=dwt(x,J)
c=x;
if J==0
    return
end
N=length(x);
if J > 0
    for k=0:J-1
        h=zeros(N,1);
        h(1)=3/8;
        h(mod(2^k,N)+1)=1/4;
        h(mod(2^(k+1),N)+1)=1/16;
        h(mod(N-2^k,N)+1)=1/4;
        h(mod(N-2^(k+1),N)+1)=1/16;
        h=h*h';
        c=real(ifftn(fftn(c).*fftn(h)));
    end
end
return