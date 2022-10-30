
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

function printTable( tab )
    n = length( tab )
    for i = 1:n
        rng, val, prb = tab[i]
        @printf( "%9s      %7d      %f\n", rng, val, prb )
    end
end

function putRowsAsCols( tab )
    rn = length(tab)
    cn = length(tab[1])
    [ tab[j][i] for i = 1:cn, j = 1:rn ]
end



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


#Ex 4.9

function dieRoll()
    rand( 1:6 )
end

type JJgame
    jess::Int
    joe::Int
end

function JJgame()
    JJgame( 2, 5 )
end

function JJnewRoll!( jjg::JJgame )
    r = dieRoll()
    if r <= 2
        jjg.jess -= 1
        jjg.joe += 1
    else
        jjg.jess += 1
        jjg.joe -= 1
    end
    jjg
end

function JJcheck( jjg::JJgame )
    jjg.jess == 0 ? ( true, "Joe" ) : jjg.joe == 0 ? ( true, "Jessica" ) : ( false, "" )
end

function JJrun1()
    jjg = JJgame()
    cnt = 0
    while !JJcheck( jjg )[1]
        cnt += 1
        JJnewRoll!( jjg )
    end
    cnt, jjg.jess == 7 ? 1 : 0
end

function JJrunN( n::Int )
    data = [ JJrun1() for i = 1:n ]
    times = map( first, data )
    # jessWins = map( x->x[2], data ) |> sum
    # mean( times ), std( times ), jessWins / n
    tab = makeTableJJ( data, 0, 80, 5 )
    printTable( tab )
end

function makeTableJJ( data, start, last, incr=5 )
    n = length( data )
    res = {}
    jessWins = data[ Bool[ x[2] == 1 for x in data ] ]
    for a = start:incr:last
        b = a + incr
        cnt = count( x -> a <= x[1] < b, data )
        cntJess = count( x -> a <= x[1] < b, jessWins )
        push!( res, { string( "[", a, ",", b, ")" ), cnt, cntJess } )
    end
    res
end

#Ex 4.10

function attemptShot( prev::Bool )
    r = rand()
    prev ? r < .6 : r < .3
end

function basketballRun()
    res = Bool[]
    prev = rand() < .5
    push!( res, prev )
    for i = 2:10
        curr = attemptShot( prev )
        push!( res, curr )
        prev = curr
    end
    res
end

function basketballRunN( n )
    data = Int[ count( x->x, basketballRun() ) for i = 1:n ]
    dh = hist( data )
    tab = { [dh[1]][2:end], dh[2] }
    println( "mean: ", mean( data ), "  SD: ", std( data ) )
    putRowsAsCols( tab )
end

## Ex 4.11

pmat = [ .6  .3  0  .1; .5  0. .5 0; 0  .2  .2  .6; .2  0  .8  0 ]

function randSt( pvec )
    len = length( pvec )
    cums = [ sum( pvec[1:i] ) for i = 1:len ]
    r = rand()
    len + 1 - count( x -> r < x, cums )
end

function MCstep( mat, st )
    pvec = mat[st,:]
    randSt( pvec )
end

function MCstepN( mat, st, n )
    curst = st
    for i = 1:n
        curst = MCstep( mat, curst )
    end
    curst
end

function MCrunsN( mat, m )
    data = [ MCstepN( mat, i, 5 ) for i = 1:4, j = 1:m ]
    newmat = [ count( x -> x == j, data[i,:] ) for i = 1:4, j = 1:4 ]
    1/m * newmat
end

## Ex 4.12

x0 = [ .1 .6 .1 .2 ]

function MCfromXN( mat, xinit, m )
    data = [ MCstepN( mat, randSt( xinit ), 5 ) for i = 1:m ]
    newmat = [ count( x -> x == i, data ) for i = 1:4 ]
    1/m * newmat
end

## Ex 4.13

absMat = [ 1 0 0 0; .5 0 .5 0; 0 .2 .2 .6; .2 0 .8 0 ]

function timeAbs( mat, st0, absset )
    cnt = 0
    curst = st0
    while !( curst in absset )
        cnt += 1
        curst = randSt( mat[curst,:] )
    end
    cnt
end

function absRunN( mat, st0, absset, n )
    data = [ timeAbs( mat, st0, absset ) for i = 1:n ]
    printTable( makeTable( data, 0, 120, 10 ) )
    mean(data), std(data)
end

## Ex 4.14

function randUnif( a, b )
    r = rand()
    a + r*(b - a)
end

function randRand( funs, weights )
    ind = randSt( weights )
    funs[ind]()
end

randRand( [ () -> randUnif(0,1), () -> randUnif(1,3) ], [.5,.5] )

using Roots

function randForF( Fun )
    r = rand()
    fzero( x -> Fun(x) - r, 1.0 )
end

## Ex 4.19

function coffeeDistr( x )
    x < 0  ? 0 :
    x < 15 ? (.25/450)*x^2 :
    x < 45 ? .125 + (1/60)*(x-15) :
    x < 90 ? 1 - (.25/1350)*(90-x)^2 : 1
end

## 4.5: Customer Flow

using Distributions

function calcNP( mu, sigma )
    p = 1 - sigma^2 / mu
    n = round( mu / p )
    p = mu / n
    n, p
end

distr = [ Binomial( calcNP( 10, 2.3 )... )
          Binomial( calcNP( 15, 2.8 )... )
          Binomial( calcNP(  5, 1.6 )... )
          Binomial( calcNP(  5, 1.6 )... )
          Binomial( calcNP(  5, 1.6 )... )
          Binomial( calcNP(  5, 1.6 )... )
          Binomial( calcNP( 15, 2.8 )... )
          Binomial( calcNP( 15, 2.8 )... )
          Binomial( calcNP( 15, 2.8 )... )
          Binomial( calcNP( 10, 2.3 )... ) ]

mall = [ 100 150 50 50 50 50 150 150 150 100 ]

function calcPcasual( n )
    .002 * n + .01
end

function calcPcasual( )
    .01
end

function simHH( k, prev )
    hr = floor( (k + 1)/2 )
    pcas = k == 1 ? calcPcasual() : calcPcasual( prev )
    core = rand( distr[ hr ] )
    distrCasual = Binomial( mall[ hr ], pcas )
    casual = rand( distrCasual )
    core, casual, core + casual
end

function simDay()
    prev = 0
    data = Any[]
    for k = 1:20
        new = simHH( k, prev )
        push!( data, [ k new... ] )
        prev = new[3]
    end
    #map( println, data )
    data
end

function simDaysN( n )
    data = [ simDay() for t = 1:n ]
    avedata = sum(data)/n
    map( println, avedata );
    sum( [ avedata[i][4] for i = 1:20 ] )
end


## Exercise bikes

dist1 = Exponential( 1 )
dist2 = Exponential( 2/3 )
dist3 = Exponential( 3/2 )

function bikeDist( t )
    t < 30 ? dist1 : t < 60 ? dist2 : dist3
end

usageT = 8
numBikes = 5


type bikeQueue
    length::Int
    bikes::Array{Float64,1}
    time::Float64
end

function newArrival( t::Number, bq::bikeQueue )
    newbq = bikeQueue( 0, [], t )
    delt = t - bq.time
    newbikes = filter( x -> x > delt, bq.bikes )
    oldbikes = filter( x -> x <= delt, bq.bikes )
    newbikes -= delt
    newopen = length( oldbikes )
    dequeue = min( bq.length, newopen )
    newbq.length = bq.length - dequeue
    for i = 1:dequeue
        push!( newbikes, usageT - delt + oldbikes[i] )
    end
    occ = length( newbikes )
    bikeless = occ == numBikes
    if !bikeless
        push!( newbikes, usageT )
    end
    newbq.bikes = newbikes
    disappointed = bikeless && bq.length == 4
    if bikeless && !disappointed
        newbq.length += 1
    end
    if disappointed
        wait = -1
    elseif bikeless
        wait = newbikes[ newbq.length ]
    else
        wait = 0
    end
    newbq, disappointed, wait
end

function bikeSim()
    bq = bikeQueue( 0, [], 0 )
    t = 0
    disappointed = 0
    customers = 0
    waiters = 0
    totalwait = 0
    while t < 90
        ## glitch at border!
        dist = bikeDist( t )
        delt = rand( dist )
        t += delt
        customers += 1
        bq, dQ, wait = newArrival( t, bq )
        if wait > 0
            waiters += 1
            totalwait += wait
        end
        if dQ; disappointed += 1 end
    end
    [ disappointed customers waiters ( totalwait / waiters ) ]
end

    
function bikeSimN( n )
    data = [ bikeSim() for i = 1:n ]
    sumdata = sum( data )
    [ sumdata[2] sumdata[1] sumdata[4] ] / n
end

function bikeEcon( prob, numCust )
    dist = Binomial( 10, prob )
    pQuit = 0.5 * ccdf( dist, 3.5 )
    numQuit = pQuit * numCust
    cost = numQuit * 30
    [ pQuit cost ]
end

function bikeEconM( m )
    global numBikes = m
    numCust, percent, wait = bikeSimN( 1000 )
    pQuit, cost = bikeEcon( percent / 100, numCust )
    bikecost = ( m - 5 ) * 100
    revenue = numCust * 30 / 10 * 30 - cost - bikecost
    [ m percent pQuit cost bikecost revenue ]
end

function bikeCompare()
    data = [ bikeEconM( m ) for m = 5:12 ]
#    datapr = Float64[ data[i][j] for i = 1:8, j = 1:6 ]
    map( println, data )
#    print( datapr )
end

## Epidemiology

distStateT = Dict{ Int, Distribution }()
distStateT[1] = Binomial(  5, 0.5 )
distStateT[2] = Binomial(  5, 0.5 )
distStateT[3] = Binomial( 10, 0.3 )
distContacts  = Binomial(  5, 0.2 )
probInfect    = 0.2
population   = 1000
quarantineRate = 0.0

type infectionState
    stage::Int
    tRemaining::Number
    quarantine::Bool
end

function infectionState( stg::Int )
    if stg == 0
        res = infectionState( stg, Inf, false )
    elseif 0 < stg < 4
        res = infectionState( stg, rand( distStateT[ stg ] ), false )
    else
        res = infectionState( stg, Inf, false )
    end
    if stg == 3
        res.quarantine = rand() < quarantineRate
    end
    if res.tRemaining > 0 
        res
    else
        infectionState( stg + 1 )
    end
end

type popEpidemic
    infectedN::Int
    immuneN::Int
    #quarantineN::Int
    infectedPop::Set{ infectionState }
end

function popEpidemic()
    popEpidemic( 0, 0, Set{ infectionState }() )
end

function nextDay( st::infectionState, pop::popEpidemic=popEpidemic() )
    ## Contact
    newInfect = 0
    if !st.quarantine && 2 <= st.stage <= 3
        numContacts = rand( distContacts ) 
        probSusceptible = ( population - pop.infectedN - pop.immuneN ) / population
        numSusceptible = rand( Binomial( numContacts, probSusceptible ) )
        newInfect = rand( Binomial( numSusceptible, probInfect ) )
    end
    ## Stage
    newimmune = 0
    newst = infectionState( st.stage, st.tRemaining, false )
    if st.tRemaining > 1
        newst.tRemaining -= 1
    else
        newst = infectionState( st.stage + 1 )
        newimmune += ( newst.stage == 4 ) ? 1 : 0
    end
    newst, newInfect, newimmune
end

function genInitInf()
    day = rand( -5:-1 )
    newinf = infectionState( 1 )
    #pop = popEpidemic()
    for d = day:-1
        newinf, _ = nextDay( newinf )
    end
    newinf
end

function genInitInfN( n )
    res = Set{ infectionState }()
    for i = 1:n
        push!( res, genInitInf() )
    end
    res
end

function countStage( popset::Set{ infectionState }, stg::Int )
    count( x -> x.stage == stg, popset )
end

function initPop( n )
    infected = genInitInfN( n )
    immune = countStage( infected, 4 )
    popEpidemic( n, immune, infected )
end

function nextDay( pop::popEpidemic )
    newinf = pop.infectedN
    newinfected = Set{ infectionState }()
    newimmune = pop.immuneN
    for sick in pop.infectedPop
        newsick, delnewinf, delimmune = nextDay( sick, pop )
        newinf += delnewinf
        newimmune += delimmune
        if newsick.stage != 4
            push!( newinfected, newsick )
        end
    end
    for i = pop.infectedN + 1 : newinf 
        newsick = infectionState( 1 )
        push!( newinfected, newsick )
        if newsick.stage == 4; newimmune += 1 end
    end
    popEpidemic( newinf, newimmune, newinfected )
end

function epidemicSim( m, qrate )
    global quarantineRate = qrate
    pop = initPop( m )
    duration = 1
    while pop.infectedN > pop.immuneN
        duration += 1
        pop = nextDay( pop )
        # if pop.infectedN != pop.immuneN + length( pop.infectedPop ) && duration > 1
        #     @hiya duration
        #     @hiya pop.infectedN
        #     @hiya pop.immuneN
        #     @hiya length( pop.infectedPop )
        #     throw( pop )
        # end            
    end
    [ pop.immuneN, duration ]
end

function epidemicSimN( m, qrate, n )
    data = [ epidemicSim( m, qrate ) for i = 1:n ]
    mean( [ data[i][1] for i = 1:n ]  ), mean( [ data[i][2] for i = 1:n ] )
end

function printEpidemicTab( n )
    @printf( "     %2d                 %2d                %2d\n", 5, 10, 20 )
    for i = 0:10:100
        @printf( "%3d   %5.1f   %5.1f     %5.1f   %5.1f     %5.1f   %5.1f\n", i, epidemicSimN( 5, i/100, n )..., epidemicSimN( 10, i/100, n )..., epidemicSimN( 20, i/100, n )...)
    end
end
