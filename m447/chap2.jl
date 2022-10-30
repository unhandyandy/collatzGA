
## 2.5 #8

popmat(s) = [ 0 0 10 25; s 0 0 0; 0 .2 .2 0; 0 0 .4 .2 ]

function growthadult( mat )
    eigs = eig( mat )
    #println( eigs )
    lambda, ind = findmax( map( real, eigs[1] ) )
    eigvec = eigs[2][:,ind]
    #@hiya eigvec
    tot = sum( eigvec )
    adult = real( eigvec[4,1] / tot )
    lambda, adult
end

#fplot( x -> growthadult( popmat( x ))[1], [.05 .2] )
#fplot( x -> growthadult( popmat( x ))[2], [.05 .2] )


## 2.5 #9

function surv( popvec )
    pop = sum( popvec )
    2 / ( 2 + pop )
end

function currentPopMat( popvec )
    s = surv( popvec )
    [ 0 2 5; s 0 0; 0 s 0 ]
end

function nextT( popvec )
    currentPopMat( popvec ) * popvec
end

function popNormalize( vec )
    tot = sum( vec )
    vec / tot
end

function step1( vec )
    vec |> nextT |> popNormalize
end

function iterate{T}( f::Function, n::Int, x::T )
    n==0 ? x : iterate( f, n - 1, f( x ) )
end
