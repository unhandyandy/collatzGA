
reload("/home/dabrowsa/lang/julia/mine/mydefs.jl")
#reload("/mnt/office/lang/julia/mine/mydefs.jl")

function dist(i,j)
    len = size( votemat )[ 1 ]
    ( len - dot( votemat[:,i], votemat[:,j] ) )/2
end
    
function disk( i, r )
    res = Set{Int}()
    for j=1:7
        if dist(i,j) <= r
            push!( res, j )
        end
    end
    res
end

function clique( i, r )
    dr = disk( i, r )
    for x in dr
        dr = intersect( dr, disk( x, r ) )
    end
    dr
end

function checkClique( clq, r )
    testclq = copy( clq )
    for x in clq
        testclq = intersect( testclq, disk( x, r ) )
    end
    testclq == clq
end

function order( clq )
    ord = 0
    while !checkClique( clq, ord )
        ord += 1
    end
    ord
end

function makeTable( fun )
    ss = subsets( 7 )
    if length( doubles ) == 2
        ss = filter( ss ) do t doubleCnt( t ) == 1 end
    end
    tab = [ ( s, fun( s ) ) for s in ss ]
    sort( tab, by=second )
end

function makeSuspectsTable( n )
    merlins = length( doubles ) == 2 ? doubles : Set([ 1:7 ])
    ss = subsets( Set([1:7]), n )
    if length( doubles ) == 2
        ss = filter( ss ) do t doubleCnt( t ) == 1 end
    end
    tab = [ ( s, m, order(s), scoreSuspects( s, m ) ) for s in ss, m in merlins ]
    tabflat =  vec( tab )
    filter( sort( tabflat, by=@fn(x[4]) ) ) do x x[4] > -Inf end
end


function maxAtOrd( sss, ord )
    level = filter( x -> x[2]==ord, sss )
    level = map( first, level )
    max = maximum( map( length, level ) )
    res = filter( x -> length( x )==max, level )
    [ ( x, order( comp( x, 7 ) ) ) for x in res ]
end

function comp( s, n )
    res = Set{Int}()
    for i = 1:n
        if !( i in s )
            res = union( res, Set{Int}(i) )
        end
    end
    res
end

function charVec( pset )
    [ ( i in pset ) ? 1 : 0 for i = 1:7 ]
end
    
function scoreTeam( vote, team, suspects )
    len = length( suspects )
    dirty = length( intersect( team, suspects ) ) > 0 ? 1 : -1
    cvec = charVec( suspects ) .* transpose( vote )
    score = ( len + sum( cvec ) * dirty ) / 2
end

function scoreTeam( vote, team, suspects, merlin )
    len = length( suspects )
    dirtCnt = length( intersect( team, suspects ) )
    dirty = dirtCnt > 0 ? 1 : -1
    cvec = charVec( suspects ) .* transpose( vote )
    mvote = vote[merlin]
    merlscr = - dirty * mvote + 1
    score = ( len + sum( cvec ) * dirty + merlscr ) / 2 
end

function scoreQuest( quest, suspects )
    expvote = length( intersect( quest[1], suspects ) )
    actvote = quest[2]
    actvote > expvote ? -Inf : 2 * actvote - expvote
end

function scoreSuspects( suspects )
    len = length( voteTeams )
    numquests = length( quests )
    teamscore = sum( [ scoreTeam( votemat[i,:], voteTeams[i], suspects ) for i = 1:len ] )
    questscore = sum( [ scoreQuest( quests[i], suspects ) for i = 1 : numquests ] )
    teamscore + questscore
end

function scoreSuspects( suspects, merlin )
    if merlin in suspects return -Inf end
    len = length( voteTeams )
    numquests = length( quests )
    teamscore = sum( [ scoreTeam( votemat[i,:], voteTeams[i], suspects, merlin )
                      for i = 1:len ] )
    questscore = sum( [ scoreQuest( quests[i], suspects )
                       for i = 1 : numquests ] )
    propscore = sum( [ scoreProposal( voteTeams[i], proposer(i), suspects, merlin )
                      for i = 1:len ] )
    teamscore + questscore + propscore
end

function scoreProposal( team, author, suspects, merlin )
    dirtCnt = length( intersect( team, suspects ) )
    if author == merlin
        if dirtCnt == 0
            1
        else
            -1
        end
    elseif author in team
        1
    elseif author in suspects
        if dirtCnt > 0
            1
        else
            -1
        end
    else
        0
    end
end

function countUnam( group )
    cnt = 0
    grpsize = length( group )
    len, num = size( votemat )
    for i = 1 : len
        test = true
        vote = votemat[i,:]
        grpvote = mapset( x -> vote[x], group )
        test = length( grpvote ) == 1
        cnt += test ? 1 : 0
    end
    cnt
end

function doubleCnt( group )
    length( intersect( doubles, group ) )
end

function filterTab( tab, rng1, rng2, suspects=false )
    if suspects
        filter( tab ) do x
            doubleCnt( s )==1 && length( x[1] ) in rng1 && x[2] in rng2
        end
    else
        filter( tab ) do x
            length( x[1] ) in rng1 && x[2] in rng2
        end
    end
end

function maxScoreFor( n )
    teams = filter( subsets( 7 ) ) do x
        n in x &&
        length( x ) == 3
    end
    scores = map( scoreSuspects, teams )
    maximum( scores )
end

function teamSuspicion( team )
    teamarray = [ team... ]
    scores = map( maxScoreFor, teamarray )
    sum( scores )
end

function bestTeam( m )
    ## teams = filter( subsets( 7 ) ) do x
    ##     length( x ) == m
    ## end
    teams = subsets( Set([1:7]), m )
    sort( [ teams... ], by=teamSuspicion )
end

function scoreFun( n )
    function( t )
        scoreSuspects( t, n )
    end
end

function proposer( i )
    res = mod( i, 7 )
    res + ( res==0 ? 7 : 0 )
end

y = 1; n = -1
votemat =
    [ y y n n y n y;
      n y n n n n n;
      n n y n n n y;
      n n n y n n n;
      n y n n y n n;
      0 0 0 0 0 0 0;
      n n n n n n y;
      n n n n n n y;
      n y n n n n n;
      0 0 0 0 0 0 0 ]

voteTeams = [ Set{Int32}([5,1])
              Set{Int32}([2,4,7])
              Set{Int32}([3,4,6])
              Set{Int32}([4,5,6])
              Set{Int32}([3,5,7])
              Set{Int32}([1,3,6])
              Set{Int32}([3,4,6])
              Set{Int32}([1,2,4])
              Set{Int32}([1,2,7])
              Set{Int32}([1,3,4]) ]

quests = [ ( Set{Int32}([5,1]), 1 )
           ( Set{Int32}([1,3,6]), 1 )
           ( Set{Int32}([1,3,4]), 0 ) ]

doubles = Set{Int}([2,3])

    
