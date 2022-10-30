
bern(p,n,r)=binomial(n,r)*p^r*(1-p)^(n-r)
berncum(p,n,a,b)=sum([bern(p,n,i) for i=a:b])

function stableVec(p)
    n = size(p)[1]
    b = [1;zeros(n-1)]
    a = [ones(1,n);I-transpose(p)]
    a = a[1:n,:]
    transpose(a^-1*b)
end
