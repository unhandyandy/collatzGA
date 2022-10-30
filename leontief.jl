
function regCol( vec )
    s = sum( vec )
    r = rand()
    r/s * vec
end

function randDiag( n )
    res = zeros(n,n)
    for i = 1:n
        r = rand(Float64)
        if rand()<0.5; r = 1/r end
        res[i,i] = r
    end
    res
end

function randLeontiev( n )
    res = rand(n,n)
    for i = 1:n
        res[:,i] = regCol( res[:,i] )
    end
    m = randDiag(n)
    m*res*m^-1
end
    
