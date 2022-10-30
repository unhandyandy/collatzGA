
include("../mydefs.jl")

import Base.print

roll() = rand(1:6)

#rollN( n ) = [ roll() for i=1:n ]

function rollN( n )
    rollSeq = [ roll() for i = 1:n ]
    [ count( x -> x==j, rollSeq ) for j = 1:6 ]
end

function resolveOne( offense, defense )
    sum = offense + defense
    sum == 7 ? 0 : sum > 7 ? 1 : -1
end

type battleState
    score::Array{Int,1}
    scoreUlt::Array{Int,1}
    army::Array{Int,2}
    init::Int
end

function clone( bs::battleState )
    battleState( [ bs.score[1], bs.score[2] ],
                [ bs.scoreUlt[1], bs.scoreUlt[2] ],
                [ bs.army[:,1] bs.army[:,2] ],
                bs.init )
end

function battleState( army1, army2 )
    score = [0,0]
    battleState( score, copy(score), [ army1 army2 ], 1 )
end

type battleMove
    offense::Bool
    dice::Array{Int,1}
end

function printScore( bs::battleState, p::Int )
    @printf( "%d::%d\n", bs.scoreUlt[p], bs.score[p] )
    #@printf( "::%d\n", bs.score[p] )
end

function print( bs::battleState )
    printScore( bs, 1)
    println( bs.army[:,1] )
    println()
    println( bs.army[:,2] )
    printScore( bs, 2)
end
function print( bs::battleState, mv::battleMove )
    printScore( bs, 1)
    println( bs.army[:,1] )
    bs.init == 1 ? print( mv ) : println()
    bs.init == 2 ? print( mv ) : println()
    println( bs.army[:,2] )
    printScore( bs, 2)
end
function print( bs::battleState, mv1::battleMove, mv2::battleMove )
    printScore( bs, 1)
    println( bs.army[:,1] )
    bs.init == 1 ? print( mv1 ) : print( mv2, false )
    bs.init == 1 ? print( mv2, false ) : print( mv1 )
    println( bs.army[:,2] )
    printScore( bs, 2)
end

function letterQ( c::Char )
    'A' <= c <= 'z'
end

function print( mv::battleMove, off::Bool=true )
    local str
    if off
        str = mv.offense ? "O: " : "D: "
    else
        str = "   "
    end
    println( str, mv.dice )
end

    
function battleMove( dice::Array{Int,1} )
    battleMove( false, dice )
end

function (==)( m1::battleMove, m2::battleMove )
    #m1.offense == m2.offense && length( m1.dice ) == length( m2.dice ) && m1.dice[1] == m2.dice[1] && sum( m1.dice - m2.dice ) == 0
    m1.offense == m2.offense && m1.dice == m2.dice
end


function parseMoveStr( str::String )
    letQ = letterQ( str[1] ) 
    off = letQ ? str[1] == 'o' : true
    dice = letQ ? str[ 2 : end - 1 ] : str[ 1 : end - 1 ]
    dicelist = [ Int(d) - 48 for d in dice ]
    battleMove( off, dicelist )
end

function generateMoves( bs::battleState, player::Int )
    res = battleMove[]
    dice = bs.army[ :, player ]
    ## if ( dice[3] > 1 )
    ##     push!( res, battleMove( [ 3, 3 ] ) )
    ## end
    stop = bs.init == player ? 1 : 0
    for i=0:stop
        b = i==0
        if ( dice[3] > 0 && dice[4] > 0 )
            push!( res, battleMove( b, [ 3, 4 ] ) )
        end
        if ( dice[2] > 0 && dice[5] > 0 )
            push!( res, battleMove( b, [ 2, 5 ] ) )
        end
        if ( dice[1] > 0 && dice[6] > 0 )
            push!( res, battleMove( b, [ 1, 6 ] ) )
        end
        for j=1:6
            if ( dice[j] > 0 )
                push!( res, battleMove( b, [ j ] ) )
            end
            if ( dice[j] > 1 )
                push!( res, battleMove( b, [ j, j ] ) )
            end
        end
    end
    res
end

function parseTwo( off::Bool, dice::Array{Int,1} )
    if length( dice ) != 2
        return dice[1]
    else
        d1, d2 = dice
        return off ? ( d1 == d2 ? d1 + 1 : 7 ) : ( d1 == d2 ? d1 - 1 : 0 )
    end
end
                                               
function resolveMoves( m1::battleMove, m2::battleMove, bs::battleState )
    res = clone( bs )
    p1 = bs.init
    p2 = 3 - p1
    offplayer = m1.offense ? p1 : p2
    defplayer = 3 - offplayer
    ## d1 = ( length( m1.dice ) == 1 || p1 == offplayer ) ? sum( m1.dice ) : 0
    ## d2 = ( length( m2.dice ) == 1 || p2 == offplayer ) ? sum( m2.dice ) : 0
    d1 = parseTwo( p1 == offplayer, m1.dice )
    d2 = parseTwo( p2 == offplayer, m2.dice )
    total = d1 + d2
    value = length( m1.dice ) + length( m2.dice )

    for d in m1.dice
        res.army[ d, p1 ] -= 1
    end
    for d in m2.dice
        res.army[ d, p2 ] -= 1
    end
    endQ, left = checkEnd( res )

    if ( total < 7 )
        res.init = defplayer
        res.score[defplayer] += value
        if endQ
            res.scoreUlt[defplayer] += 1
        end
    elseif ( total > 7 )
        res.init = offplayer
        res.score[offplayer] += value
        if endQ
            res.scoreUlt[offplayer] += 1
        end
    end
    #res.army[ m1.dice[1], p1 ] -= length( m1.dice )
    #res.army[ m2.dice[1], p2 ] -= length( m2.dice )
 
    if endQ
        if left > 0
            res.score[1] +=  left
            if total == 7 ; res.scoreUlt[1] += 1 end
        elseif left < 0
            res.score[2] += -left
            if total == 7 ; res.scoreUlt[2] += 1 end
        end
    end
    return res
end

function checkEnd( bs::battleState )
    s1, s2 = map( sum,  [ bs.army[ :, i ] for i = 1:2 ] )
    s1 * s2 == 0, s1 - s2
end

function bestMove( bs::battleState )
    p1 = bs.init
    p2 = 3 - p1
    mvs = generateMoves( bs, p1 )
    if length( mvs ) == 0
        return false, bs.score[1] - bs.score[2]
    end
    data = [ bestMove( bs, m )[2] for m in mvs ]
    opt = p1 == 1 ? maximum( data ) : minimum( data )
    pos = findin( data, [opt] )[1]
    return mvs[ pos ], opt
end

function bestMove( bs::battleState, mv::battleMove )
    p1 = bs.init
    p2 = 3 - p1
    mvs = generateMoves( bs, p2 )
    if length( mvs ) == 0
        return false, bs.score[1] - bs.score[2]
    end
    resolves = [ resolveMoves( mv, m, bs ) for m in mvs ]
    data = [ bestMove( r )[2] for r in resolves ]
    opt = p2 == 1 ? maximum( data ) : minimum( data )
    pos = findin( data, [opt] )[1]
    return mvs[ pos ], opt
end

function BDloop( armySize::Int=6, init::Int=1, score::Array{Int,1}=[0,0], scoreUlt::Array{Int,1}=[0,0] )
    a1 = rollN( armySize )
    a2 = rollN( armySize )
    bs = battleState( a1, a2 )
    bs.init = init
    bs.score = score
    bs.scoreUlt = scoreUlt
    print( bs )
    newscore, newscoreUlt = gameloop( bs )
    a1, a2 = a2, a1
    bs = battleState( a1, a2 )    
    bs.init = 3 - init
    bs.score = newscore
    bs.scoreUlt = newscoreUlt
    print( bs )
    newscore, newscoreUlt = gameloop( bs )
    BDloop( armySize, init, newscore, newscoreUlt )
end

function getUserMove( bs::battleState, p::Int=2 )
    println( "Your move:" )
    mv = parseMoveStr( readline(STDIN) )
    mvs = generateMoves( bs, p )
    if mv in mvs
        return mv
    else
        return getUserMove( bs, p )
    end
end

function gameloop( bs::battleState )
    local m1, m2
    if bs.init == 1
        m1, val = bestMove( bs )
    else
        m1 = getUserMove( bs )
    end
    println(); println()
    print( bs, m1 )
    if bs.init == 1
        m2 = getUserMove( bs )
    else
        m2, val = bestMove( bs, m1 )
    end
    println(); println()
    print( bs, m1, m2 )
    newbs = resolveMoves( m1, m2, bs )
    print( newbs )
    if checkEnd( newbs )[1]
        println( "Round Over" )
        println()
        return newbs.score, newbs.scoreUlt
    else
        gameloop( newbs )
    end
end

type diceDuelState
    bs::battleState
    pos::Array{Int,1}
end

function clone( dds::diceDuelState )
    diceDuelState( clone( dds.bs ), copy( pos ) )
end

type diceDuelMove
    slideQ::Bool
    player::Int
    move::Union(Int,battleMove)
end

floorLength = 8

function generateMoves( dds::diceDuelState, p::Int, response::Bool )
    res = diceDuelMove[]
    s1, s2 = dds.pos
    if ( s1 - s2 == 1 )
        bmvs = generateMoves( dds.bs, p )
        for m in bmvs
            push!( res, diceDuelMove( false, p, m ) )
        end
    end
    army = dds.bs.army[:,p]
    a, b = ( p == 1 ) ? ( s2 - s1 + 1, floorLength - s1 ) : ( - s2, s1 - s2 - 1 )
    if !response
        for m = 1:6
            if ( army[m] > 0 || army[7 - m] > 0 )
                if m <= b
                    push!( res, diceDuelMove( true, p, m ) )
                end
                if -m >= a 
                    push!( res, diceDuelMove( true, p, -m ) )
                end
            end
        end
    end
    return res
end

function resolveMoves( m::diceDuelMove, dds::diceDuelState )
    res = clone( dds )
    if m.slideQ
        res.pos[m.player] += m.move
    end
    res.bs.init = m.player == 1 ? 2 : 1
    return res
end

function resolveMoves( m1::diceDuelMove, m2::diceDuelMove, dds::diceDuelState )
    res = clone( dds )
    newbs = resolveMoves( m1.move, m2.move, dds.bs )
    res.bs = newbs
    return res
end

function bestMove( dds::diceDuelState )
    p1 = dds.bs.init
    p2 = 3 - p1
    mvs = generateMoves( dds, p1, false )
    if length( mvs ) == 0
        return false, dds.bs.score[1] - dds.bs.score[2]
    end
    data = [ bestMove( dds, m )[2] for m in mvs ]
    opt = p1 == 1 ? maximum( data ) : minimum( data )
    pos = findin( data, [opt] )[1]
    return mvs[ pos ], opt
end

function bestMove( dds::diceDuelState, mv::diceDuelMove )
    p1 = dds.bs.init
    p2 = 3 - p1
    mvs = generateMoves( dds, p2, true )
    if length( mvs ) == 0
        return false, dds.bs.score[1] - dds.bs.score[2]
    end
    resolves = [ resolveMoves( mv, m, dds ) for m in mvs ]
    data = [ bestMove( r )[2] for r in resolves ]
    opt = p2 == 1 ? maximum( data ) : minimum( data )
    pos = findin( data, [opt] )[1]
    return mvs[ pos ], opt
end

    
## army1 = rollN(7)
## army2 = rollN(7)
## bs = battleState( army1, army2 )
army1 = rollN(6)
army2 = rollN(6)
bs = battleState( army1, army2 )
