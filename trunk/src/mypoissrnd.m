function k=mypoissrnd(lambda,m,n)
k=ones(m,n);
for x=1:m
    for y=1:n
        p=1;
        l=exp(-lambda);
        while p>l
            p=p*rand(1);
            k(x,y)=k(x,y)+1;
        end
    end
end
k=k-1;
return