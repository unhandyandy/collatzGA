
include( "mydefs.jl" )

type disjunct
    vars::Array{Bool,1}
    signs::Array{Bool,1}
    constant::(Bool,Bool)
    comment::String
end

truedj = disjunct( [], [], ( true, true ), "True" )
falsedj = disjunct( [], [], ( true, false ), "False" )


function disjunct( vars::Array{Int,1}, sgns::Array{Bool,1}, len::Int, com::String="" )
    newdisj = disjunct( [ false for i = 1 : len ],
                        [ false for i = 1 : len ],
                        ( false, false ),
                        com )
    function assign( i::Int, x::Bool )
        newdisj.vars[ i ] = true
        newdisj.signs[ i ] = x
    end
    map( assign, vars, sgns )
    newdisj
end

function disjunct( vars::Array{(Int,Bool),1}, len::Int, com::String="" )
    vs = [ x[1] for x in vars ]
    ss = [ x[2] for x in vars ]
    disjunct( vs, ss, len, com )
end

function (+)( p::disjunct, q::disjunct )
    if p.constant[1]
        if p.constant[2]
            return truedj
        else
            return q
        end
    elseif q.constant[1]
        if q.constant[2]
            return truedj
        else
            return p
        end
    end
    pl, ql, len = level( p, q )
    newdisjunct = disjunct( [ false for i = 1 : len ],
                           [ false for i = 1 : len ],
                           ( false, false ),
                           p.comment + " + " + q.comment )
    trivial = true
    for i = 1 : len
         if p.vars[ i ] && q.vars[ i ] && ( p.signs[ i ] == !q.signs[ i ] )
            return truedj
        end
        newvar = p.vars[ i ] || q.vars[ i ]
        trivial = trivial && !newvar
        if newvar
            newdisjunct.signs[ i ] = p.signs[ i ] || q.signs[ i ]
            newdisjunct.vars[ i ] = true
        end
    end
    if trivial
        falsedj
    else 
        newdisjunct
    end
end
        
    

function level( p::disjunct, q::disjunct )    
    lenp, lenq = length( p.vars ), length( q.vars )
    if lenp == lenq ; return p, q, lenp end

    lenmax, lenmin = max( lenp, lenq ), min( lenp, lenq )
    diff = lenmax - lenmin
    short, long = ( lenp == lenmin ) ? (p,q) : (q,p)
    newshortvars = [ short.vars..., [ false for i = 1 : diff ]... ]
    newshortsigns = [ short.signs..., [ false for i = 1 : diff ]... ]
    newshort = disjunct( newshortvars, newshortsigns, (false,false) )
    return newshort, long, lenmax
end

function toVarList( p::disjunct )
    res = (Int,Bool)[]
    len = length( p.vars )
    for i = 1 : len
        if p.vars[ i ]
            push!( res, ( i, p.signs[ i ] ) )
        end
    end
    res
end

function varPairToStr( pr::(Int,Bool) )
    b26 = reverse( digits( pr[ 1 ], 26 ) )
    charlst = map( x -> string(char( x + 97 )), b26 )
    str = foldl( *, "", charlst )
    pr[ 2 ] ? str : "-" * str
end
    

function toString( p::disjunct )
    if p.constant[ 1 ]
        return p.constant[ 2 ]
    end
    varlst = toVarList( p )
    strlst = map( varPairToStr, varlst )
    res = foldl( (x,y) -> (x * "|" * y), "", strlst )
    res[2:end]
end

function clone( p::disjunct )
    disjunct( copy( p.vars ), copy( p.signs ), copy( p.constant ), p.comment )
end

function resolve( p::disjunct, q::disjunct, ind::Int )
    pc, qc = clone( p ), clone( q )
    ( pc.vars[ ind ], pc.signs[ ind ] ) = ( false, false )
    ( qc.vars[ ind ], qc.signs[ ind ] ) = ( false, false )
    pc + qc
end

function (*)( p::disjunct, q::disjunct )
    if p.constant[1]
        if p.constant[2]
            return q
        else
            return falsedj
        end
    elseif q.constant[1]
        if q.constant[2]
            return p
        else
            return falsedj
        end
    end
    pl, ql, len = level( p, q )
    for i = 1 : len
        if pl.vars[ i ] && ql.vars[ i ] && ( pl.signs[ i ] == ! ql.signs[ i ] )
            return resolve( pl, ql, i )
        end
    end
    return truedj 
end

function (==)( p::disjunct, q::disjunct )
    if p.constant[1]
        if q.constant[1]
            return p.constant[2] == q.constant[2]
        else
            return false
        end
    elseif q.constant[1]
        return false
    end
    pl, ql, len = level( p, q )
    return all( map( ==, pl.vars, ql.vars ) ) && all( map( ==, pl.signs, ql.signs ) )
end

function strToTerm( str::String )
    sign = true
    numstr = trimWhite( str )
    if numstr[1] == '-'
        sign = false
        numstr = str[2:end]
    end
    digstr = map( x -> ( x < 107 ?  x - 49 : x - 10 ), numstr )
    #@hiya digstr
    ind = parseint( digstr, 26 ) 
    return ind, sign
end

function strToDisj( str::String )
    strlst = split( str, "|" )
    termlst = map( strToTerm, strlst )
    len = maximum( map( first, termlst ) )
    disjunct( termlst, len )
end

        
function eliminate!( dj::disjunct, ind::Int, b::Bool )
    len = length( dj.vars )
    if !dj.constant[1] && ind <= len && dj.vars[ ind ]
        if dj.signs[ ind ] != b
            dj.vars[ ind ], dj.signs[ ind ] = false, false
            checkEmpty!( dj )
        else 
            dj.constant = ( true, true )
        end
        dj.comment *= " - " * string( ind )
    end
end

function eliminate( dj::disjunct, ind::Int, b::Bool )
    djcopy = clone( dj )
    eliminate!( djcopy, ind, b )
    djcopy
end

function eliminateAll( djset::Set{disjunct}, ind::Int, b::Bool )
    function mapper( dj::disjunct )
        eliminate( dj, ind, b )
    end
    mapset( mapper, djset )
end
  
function checkEmpty!( dj::disjunct )
    if !dj.constant[1]
        if !any( dj.vars )
            dj.constant = ( true, false )
        end
    end
end

function checkSingle( dj::disjunct )
    if dj.constant[1]
        return -1, true
    end
    cnt = 0
    ind = -1
    len = length( dj.vars )
    for i = 1 : len 
        if dj.vars[ i ]
            cnt += 1
            if cnt > 1
                return -1, false
            end
            ind = i
        end
    end
    if cnt == 1
        return ind, dj.signs[ ind ]
    else
        return -1, true
    end
end
 
function countVars( dj::disjunct )
    if dj.constant[1]
        0
    end
    count( first, dj.vars )
end

function checkInvariance( djset::Set{disjunct}, ind::Int, b::Bool )
    function getsign( dj::disjunct )
        len = length( dj.vars )
        ind > len || dj.constant[1] || !dj.vars[ ind ] || dj.signs[ ind ] == b
    end
    all( map( getsign, djset ) )
end

function removeInvariant( djset::Set{disjunct}, ind::Int )
    function replacer( dj::disjunct )
        len = length( dj.vars )
        !dj.constant[1] && ind <= len && dj.vars[ ind ] ? truedj : dj
    end
    mapset( replacer, djset )
end

function removeInvariantsAll( djset::Set{disjunct} )
    len = maximum( map( x -> length( x.vars ), djset ) )
    res = djset
    for i = 1 : len
        for b in [ true, false ]
            if checkInvariance( djset, i, b )
                res = removeInvariant( res, i )
            end
        end
    end
    res
end

function eliminateSinglesAll!( djset::Set{disjunct} )
    for dj in djset
        ind, b = checkSingle( dj )
        if ind > 0 
            for x in djset; eliminate!( x, ind, b ) end
        end
    end
end

function checkTrivial( djset::Set{disjunct} )
    for dj in djset
        if !dj.constant[1] 
            return false
        end
    end
    return true
end

function findLiteral( djset::Set{disjunct} )
    for dj in djset
        if !dj.constant[1]
            return findfirst( x->x, dj.vars )
        end
    end
    return -1
end

function clone( djset::Set{disjunct} )
    mapset( clone, djset )
end

function falseP( dj::disjunct )
    dj.constant == ( true, false )
end

function checkSetFalse( djset::Set{disjunct} )
    #println( "checking for false..." )
    for dj in djset
        if falseP( dj )
            #println( "...false!" )
            return true, dj.comment
        end
    end
    #println( "... no." )
    return false, ""
end

function dpll( djset::Set{disjunct} )
    checkFalse, com = checkSetFalse( djset )
    if checkFalse 
        println( "Violated: ", com )
        return false
    elseif checkTrivial( djset )
        #println( "Trivial." )
        return true
    end
    djsetcopy = clone( djset )
    eliminateSinglesAll!( djsetcopy )
    djsetcopy = removeInvariantsAll( djsetcopy )
    #println( "Looking for fresh literals." )
    ind = findLiteral( djsetcopy )
    if ind < 0 
        #println( "No literals left." )
        return dpll( djsetcopy )
    end
    #println( "int = ", ind )
    nextgen1 = eliminateAll( djsetcopy, ind, true )
    nextgen2 = eliminateAll( djsetcopy, ind, false )
    return dpll( nextgen1 ) || dpll( nextgen2 )
end

function simplify( djset::Set{disjunct} )
    djsetcopy = clone( djset )
    #eliminateSinglesAll!( djsetcopy )
    if checkSetFalse( djsetcopy )[1]
        return Set{disjunct}({ falsedj })
    elseif checkTrivial( djsetcopy )
        return Set{disjunct}({ truedj })
    end
    return djsetcopy
end

