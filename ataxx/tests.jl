using Test

include("rectUtils.jl")

res0 = [1  2  3  4  4  4;
        3  3  3  3  4  5];

res1 = neighbor2Locs(2,5,5);
res2 = neighborNLocs(2,5,5,2);

function nbrsEq(n1,n2)
    if length(n1)!=length(n2)
        return false
    end
    res = true;
    (_,w) = size(n1);
    i = 1;
    while res && i<=w
        res = res && any(all(n1.==n2[:,i],dims=1));
        i += 1;
    end
    return res
end


@test nbrsEq(res0,res1)
@test nbrsEq(res0,res2)
