#-*- julia; eval: (auto-fill-mode); eval: (auto-complete-mode) -*-

include("/home/dabrowsa/lang/julia/mine/mydefs.jl")
#reload("/media/office/lang/julia/mine/mydefs.jl")

function roll()
    rand(1:6)
end

function rollN( n )
    rollSeq = [ roll() for i = 1:n ]
    [ count( x -> x==j, rollSeq ) for j = 1:6 ]
end

type targSum
    domain::Range
    val::Int
end

function targSum( v::Int )
    targSum( 1:3, v )
end

type targCom
    val::Array{Int,1}
end

target = Union{ targSum, targCom }

function numDice( targ::targCom )
    sum( targ.val )
end

function captureP( targ::targSum, rollNums::Array{Int,1} )
    weights = [ i in targ.domain ? i : 0 for i = 1:6 ]
    res = dot( weights, rollNums ) >= targ.val
    num = 0
    rnums = [ i in targ.domain ? rollNums[i] : 0 for i = 1:6 ]
    cnt = 0
    if res
        sum = 0
        ind = 6
        while sum < targ.val
            if rnums[ind] > 0
                rnums[ind] -= 1
                sum += ind
                cnt += 1
            else
                ind -= 1
            end
        end
    end
    res, cnt
end

function captureP( targ::targCom, rollNums::Array{Int,1} )
    nums = zip( targ.val, rollNums )
    res = reduce( true, nums ) do x, pr
        x && pr[1] <= pr[2]
    end
    res, numDice( targ )
end
 
function estProb( targ::target, n::Int, m::Int )
    data = [ captureP( targ, rollN( n ) )[1] for i = 1:m ]
    succ = count( x -> x, data )
    succ / m 
end

type castle
    walls::Array{target,1}
end

function clone( cst::castle )
    castle( copy( cst.walls ) )
end

function attack( cst::castle, rnums::Array{Int,1} )
    test = false
    ind = 0
    len = length( cst.walls )
    numD = sum( rnums )
    dnum = 0
    while !test && ind < len
        ind += 1
        test, dnum = captureP( cst.walls[ind], rnums )
    end
    if test
        #newwalls = [ cst.walls[ 1 : ind - 1 ]..., cst.walls[ ind + 1 : end ] ]
        newwalls = copy( cst.walls )
        deleteat!( newwalls, ind )
        newcst = castle( newwalls )
        numD -= dnum
    else
        newcst = cst
        numD -= 1
    end
    test, newcst, numD
end

function besiege( cst::castle, ndice::Int )
    if length( cst.walls ) == 0
        return true
    elseif ndice == 0
        return false
    end
    rnums = rollN( ndice )
    #@hiya rnums
    oldcst = clone( cst )
    sortWalls!( oldcst, ndice )
    test, newcst, nd = attack( oldcst, rnums )
    besiege( newcst, nd )
end
    
function estProb( cst::castle, n::Int, m::Int )
    data = [ besiege( cst, n )[1] for i = 1:m ]
    succ = count( x -> x, data )
    succ / m 
end

function harder( wall1::target, wall2::target, n::Int )
    estProb( wall1, n, 100 ) < estProb( wall2, n, 100 )
end

function sortWalls!( cst::castle, n::Int )
    function sortFun( w1::target, w2::target )
        harder( w1, w2, n )
    end
    sort!( cst.walls, lt=sortFun )
end

soloCastle01 = castle( target[ 
       targCom([0,0,0,0,0,2]),
       targCom([0,0,0,1,0,0]),
       targCom([0,0,0,0,1,0]),
       targSum(4)
       ] )

soloCastle02 = castle( target[ 
       targCom([0,0,0,0,0,2]),
       targSum(4),
       targCom([0,0,0,1,0,0]),
       targCom([0,0,0,0,1,0])
       ] )

soloCastle03 = castle( target[ 
       targCom([0,0,0,0,0,2]),
       targCom([0,0,0,1,0,0]),
       targSum(4),
       targCom([0,0,0,0,1,0])
       ] )

blueCastle1 = castle( target[ 
       targCom([0,0,0,0,2,0]),
       targCom([0,0,0,0,0,1])
       ] )

blueCastle2 = castle( target[ 
       targCom([0,0,0,0,0,1]),
       targCom([0,0,0,0,1,0]),
       targCom([0,0,0,1,0,0]),
       targSum( 3 )
       ] )
sortWalls!( blueCastle2, 7 )

blueCastle3 = castle( target[ 
       targCom([0,0,0,1,1,0]),
       targCom([0,0,0,1,1,0]),
       targSum( 3 )
       ] )
sortWalls!( blueCastle3, 7 )

