
include("../mydefs.jl")

roll() = rand(1:6)

function rollN( n )
    rollSeq = [ roll() for i = 1:n ]
    #[ count( x -> x==j, rollSeq ) for j = 1:6 ]
end

## function countSums( rolls::Array{Int,1} )
##     prods = [ 0 for i=1:12 ]
##     for i=1:6
##         prods[2*i] += binomial( rolls[i], 2 )
##         for j=i+1:6
##             prods[i+j] += rolls[i]*rolls[j]
##         end  end
##     prods
## end

function pairs( rolls::Array{Int,1} )
    total = sum( rolls )
    res = Set{Array{Int,1}}()
    for i = 2:4
        x = rolls[1] + rolls[i]
        y = total - x
        push!( res, [x,y] )
    end
    res
end

function filterPairs( prs::Set{Array{Int,1}}, active::Array{Int,1} )
    f1( pr::Array{Int,1} ) = filter( x -> ( x in active ), pr )
    preres = mapset( f1, prs )
    filter( x -> (length( x ) > 0 ), preres)
end

function track( x::Int )
    13 - 2 * abs( x - 7 )
end

function progress( x::Int )
    1/track(x)
end
function progress( x::Int, y::Int )
    progress(x) + progress(y)
end

function progress( prs::Set{Array{Int,1}} )
    if length( prs ) == 0
        return 0
    end
    maximum( [ progress( x... ) for x in prs ] )
end
    
    function makeData( active )
        [ progress( filterPairs( pairs( rollN(4) ), active ) ) for i = 1:1000000 ]
    end

function getStats( active )
    data = makeData( active )
    prob = .001 * count( x -> ( x > 0 ), data )
    ave = mean( data ) / prob
    println( "prob: ", prob )
    println( " ave: ", ave )
    prob, ave
end

function stopSign( active )
    prob, ave = getStats( active )
    ave * prob / ( 1 - prob )
end
