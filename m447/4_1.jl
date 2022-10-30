
# Ex 4.1

function randRate( )
    rnd = rand()
    if rnd < .3
        -.01
    elseif rnd < .4
        0
    elseif rnd < .8
        .01
    else
        .02
    end
end

function randPop()
    pop = 100
    for i = 1:20
        pop *= 1 + randRate()
    end
    pop
end

function generateData( n )
    [ randPop() for i = 1:n ]
end

function countRange( data, a, b )
    count( x -> a <= x < b, data )
end

function makeTable( data, start=85, last=140, incr=5 )
    n = length( data )
    res = {}
    for a = start:incr:last
        b = a + incr
        cnt = countRange( data, a, b )
        push!( res, { string( "[", a, ",", b, ")" ), cnt, cnt/n } )
    end
    res
end

# function printTable( tab )
#     n = length( tab )
#     for i = 1:n
#         rng, val, prb = tab[i]
#         @printf( "%9s      %7d      %f\n", rng, val, prb )
#     end
# end


## Ex 4.2

type ftState
    last::Bool
    rate::Float64
end

function ftState()
    ftState( false, 0.5 )
end

function update!( fts::ftState, hit::Bool )
    fts.last = hit
    fts.rate = hit ? 0.5 : 0.7
    fts
end

function ftOne( fts::ftState )
    rand() < fts.rate
end

function ftIterate!( fts::ftState )
    hit = ftOne( fts )
    update!( fts, hit )
    hit, fts
end

function ftIterN( n::Int )
    fts = ftState()
    cnt = 0
    for i = 1:n
        hit, _ = ftIterate!( fts )
        cnt += hit ? 1 : 0
    end
    cnt
end

function ftRunsMN( m::Int, n::Int )
    [ ftIterN( n ) for i = 1:m ]
end


#Ex 4.3

type catsCol
    c::Bool
    a::Bool
    t::Bool
    s::Bool
end

function catsCol()
    catsCol( false, false, false, false )
end

function randLetter()
    r = rand()
    r < .3 ? "c": r < .4 ? "a" : r < .7 ? "t" : "s"
end

function update!( cc::catsCol, l::String )
    if l == "c"
        cc.c = true
    elseif l == "a"
        cc.a = true
    elseif l == "t"
        cc.t = true
    else
        cc.s = true
    end
    cc
end

function checkCats( cc::catsCol )
    cc.c && cc.a && cc.t && cc.s
end

function catsOne (  )
    cc = catsCol()
    cnt = 0
    while !checkCats( cc )
        cnt += 1
        update!( cc, randLetter() )
    end
    cnt
end

function catsRunN( n::Int )
    [ catsOne() for i = 1:n ]
end

#Ex 4.7

function randDemand()
    r = rand() 
    r < .2 ? 100 : r < .4 ? 125 : r < .7 ? 150 : r < .9 ? 175 : r < .95 ? 200 : 225 
end

function kitProfit( ord::Int )
    dem = randDemand()
    10 * ord - 20 * ( dem < ord ? ord - dem : 0 )
end

function meanSD( fun, n::Int )
    data = Int[ fun() for i = 1:n ]
    #@hiya data
    mean( data ), std( data )
end

function kitStats( ord:: Int )
    m, sd = meanSD( function() kitProfit( ord ) end, 2000 )
    ord, m, sd
end

kitData = [ [kitStats(i)...] for i = 100:25:225 ]

#Ex 4.8

function randWeather()
    r = rand()
    r < .5 ? "s" : r < .8 ? "c" : "r"
end

function randSales( w )
    r = rand()
    cut = w == "s" ? .5 : w == "c" ? .8 : .3
    r < cut ? 1 : 0
end

function sell( n )
    [ randSales( randWeather() ) for i = 1:n ]
end

function sellRunMN( m, n )
    data = Int[ sell( m ) |> sum for i = 1:n ]
    makeTable( data, 0, 30, 1 ) |> printTable
    mean( data ), std( data )
end

    
