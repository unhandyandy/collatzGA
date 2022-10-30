#-*- mode: julia; eval: (auto-fill-mode) -*-
-
# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

include("mydefs.jl")

type cSet
    old :: Set{Int}
    new :: Set{Int}
end

# <codecell>

# testCSet = cSet(Set(2,3,7),Set(9))

# <codecell>

cSort(s::cSet) = begin
    s = union(s.old,s.new)
    res = [ x for x in s ]
    sort( res )
end     

# <codecell>

cAll(s::cSet) = union(s.old,s.new)

# <codecell>

# cAll(testCSet)

# <codecell>

backOne(n::BigInt)=begin
    r = n % 6   
    if r == 4 
        res = Set{BigInt}(div((n-1),3)) 
    else
        res = Set{BigInt}()
    end
    union!(res,Set{BigInt}(2*n))
    res
end

# collD = Dict{BigInt,Int}()
# collD[big(1)]=1
# collD[big(0)]=1
# collNew = Set(big(1))
collD = Dict{BigInt,Int}()
collD[big(1)]=0
flowD = Dict{BigInt,BigInt}()
flowD[big(1)]=4


backOne()=begin
    global collNew
    new = Set{BigInt}()
    #old = union(cs.old,cs.new)
    for n in collNew
        union!( new, backOne(n) )
    end
    filter!( n -> !haskey(collD,n), new)
    t = collD[big(0)] + 1
    for n in new
        collD[n] = t
    end
    collD[big(0)] = t
    collNew = new
    nothing
end


function forwardOne{T<:Number}(n::T) 
    ( n % 2 == 0 ) ? div(n,2) : 3*n + 1
end

function findOrder{T<:Number}(n::T) 
    nxt = get!( flowD, n ) do
        forward(n)
    end 
    ord = get!( collD, n ) do
        1 + findOrder( nxt )[1]
    end
    ord, nxt
end

connected(x::Int,y::Int) = forwardOne(x) == y || forwardOne(y) == x
    


iterate(f::Function,n::Int,x...) = (n==0) ? x : iterate(f,n-1,apply(f,x...))    
iterate(f::Function,n::Int) = (n==0) ? nothing : ( f(); iterate(f,n-1) )

# <codecell>

#cs1 = cSet(Set{BigInt}(),Set{BigInt}(big(1)))
cs1 = Set(big(1))

BtoN(x::Bool) = x ? 1 : 0

# <codecell>

#csHist100 = [ BtoN(in(i,cAll(iterate(backOne,j,cs1)))) for i=1:3, j=0:2 ]

# <codecell>

# iterMemory = [cs1]

iterMemo(n) = begin
    len = collD[big(0)]
    if n > len
        iterate( backOne, n - len )
    end
end

# <codecell>

#iterMemo(20)



# csHist(n) = begin
#     iterMemo(n)
#     [ BtoN( haskey(collD,i) && collD[i] <= j ) 
#        for i=big(2*n):-1:1, j=1:n ]
# end

collGraph(r::Range1; showflow=false, xtransf=(z->z) ) = begin
    hold(false)
    x = [r]
    y = map( x->findOrder(x)[1], x )
    x = map( xtransf, x )
    plot( x, y, ".")
    if showflow
        for i=r
            ord, nxt = findOrder(i)
            if nxt <= r[end]
                nxtord, = findOrder(nxt)
                oplot( [xtransf(i), xtransf(nxt)], [ord, nxtord], "m-" )
            end
        end
    end
    oplot()
end

collGraph(n::Int; showflow=false, xtransf=(z->z) ) = collGraph( 1:n; showflow=showflow, xtransf=xtransf )

trailGraph(n::Int, color="k") = begin
    ord, nxt = findOrder(n)
    x = [ord:-1:0]
    y = nestList( n->findOrder(n)[2], ord, n )
    plot( x, y, color * "-" )
end

trailGraph( nset::Set{Int} ) = begin
    colors = [ "k", "m", "c", "r", "g", "b", "y" ]
    len = length( colors )
    hold(false)
    plot()
    hold( true )
    colind = 0
    for n in nset
        col = colors[ ( colind % len ) + 1 ]
        trailGraph( n, col )
        colind += 1
    end
    plot()
end


#(==)(a::cSet,b::cSet) = a.old == b.old && a.new == b.new

using Winston

# <codecell>

#imagesc(csHist(5))

#makeBar(m,n) = 


randForward( n::Int ) = if ( rand() < 0.5 )
    del = ( n%2 == 0 ) ? 0 : ( rand()<0.5 ? 0 : 1 )
    div(n,2) + del
else
    3*n + 1
end

toggleRand() = begin
    global forward, collD, flowD
    forward = ( forward == forwardOne ) ? randForward : forwardOne
    collD = Dict{Int,Int}()
    collD[1]=0
    flowD = Dict{Int,Int}()
    flowD[1]=4
end;


type affineDiaph
    m :: Rational
    b :: Rational
end

outAD( f::affineDiaph, n::Int ) = begin
    val = f.m * n + f.b
    ( val == int( val ) ) ? Set(int(val)) : Set{Int}()
end

type genColl
    funs :: Array{affineDiaph,1}
    component :: cSet
    orderDB :: Dict{Int,Int}
    order :: Int
end

initGenColl( f... ) = begin
    funs = [f...]
    cmp = cSet(Set{Int}(),Set(0,1))
    dct =  Dict{Int,Int}()
    new = genColl( funs, cmp, dct, 0 )
    updateDB( new )
end

updateDB( gc::genColl ) = begin
    t = gc.order
    for n in gc.component.new
        gc.orderDB[ n ] = t
    end
    gc
end


iterColl( gc::genColl ) = begin
    new = Set{Int}()
    old = gc.component.old
    for n in gc.component.new 
        for f in gc.funs
            union!( new, outAD( f, n ) )
        end
    end
    union!( old, gc.component.new )
    setdiff!( new, old )
    gc.component = cSet( old, new )
    gc.order += 1
    updateDB( gc )
    cSort( gc.component )
end

function findOrder{T<:Number}( n::T, gc::genColl )
    get!( gc.orderDB, n ) do
    iterColl( gc )
    findOrder( n, gc )
    end
end


inverse( f::affineDiaph ) = affineDiaph( 1//f.m, -f.b//f.m )

inverse( gc::genColl ) = begin
    funs = [ inverse(f) for f in gc.funs ]
    apply( initGenColl, funs )
end


function collOF{T<:Number}(n::oddFactor{T}; fun=(d->(3*d+1)),
                           trim1=true )
    d = n.odd
    p = n.power - ( trim1 ? 1 : 0 )
    newd = fun(d)
    oddFactor( T(1), p ) * oddFactor( newd )
end


    
function tostring{T<:Number}( m::oddFactor{T}, width=400 ) 
    str =  bin( m.odd )
    pad = width - m.power - length( str )
    if ( pad < 0 ) pad = 0 end
    #fun = x->x
    res = repeat( " ", pad ) * str * repeat(" ", m.power )
    res = replace( res, "0", " " )
    res = replace( res, "1", "X" )
    #res * " : " * string( reduceOddLeft( m.odd ) ) * ", " * string(
    #reduceOddRight( m.odd ) )
    #res *  " \n   " * string( metric( m ) )
end

function tovector{T<:Number}( m::oddFactor{T} ) 
    str = tostring( m )
    [ ( str[i] == 'X' ? 1 : 0 ) for i=1:length(str) ]'
end

function binFlow{T<:Number}(n::T, metric = x -> genHalfFT( bin( x.odd ), "01", 'r' )' ) 
    ord, = findOrder( n )
    map( x -> tostring( x ), nestList( collOF, div(ord,2), oddFactor(n) ) )
end

function showFlowRand( wid )
    len = int( wid / 0.4 ) + 1
    cols =  wid + int( 0.6 * len )
    start = randLargeOdd( wid )
    showFlow( start, len, cols )
end

function showFlow{T<:Number}( start::oddFactor{T}, len::Int, cols::Int )
    #len = int( log( 2, start.odd ) ) + 2
    data = map( x -> tostring( x, cols), nestList( collOF, len, start ))
    mat = [ ( data[i][j] == 'X' ? 0 : 1 ) for i=len+1:-1:1, j=1:cols ]
    imagesc( mat )
end

function showFlowColor{T<:Number}( start::oddFactor{T}, len::Int,
                                  cols::Int )
    function tolist( x::oddFactor ) 
        str = bin( x.odd )
        colors = cwallColor( str )
        wid = length( str )
        left = cols - wid - x.power
        if left<0 left=0 end
        #@hiya wid
        #@hiya left
        vcat( zeros( left ), colors, zeros(x.power ) )
    end
    data = map( tolist, nestList( collOF, len, start ))
    #@hiya data
    mat::Array{Float64}
    mat = [ float(data[i][j]) for i=len+1:-1:1, j=1:cols ]
    # @hiya mat
    imagesc( mat, (0.0,1.0) )
end

function showFlow( str::String, showFun=showFlow )
    start = parseint( BigInt, str, 2 )
    wid = int( log( 2, start ) ) + 2
    len = int( wid / 0.4 ) + 1
    cols =  wid + int( 0.6 * len )
    showFun( oddFactor( start ), len, cols )
end

forward = forwardOne

makeBinAltRight( len::Int ) = begin
    res = 1
    pow = 2
    while pow < len
        res += 2^pow
        pow += 2
    end
    res
end
makeBinAltLeft( len::Int ) = begin
    res = makeBinAltRight( len )
    if len % 2 == 1
        res
    else
        res * 2
    end
end

intLog2( n::Int ) =  int( floor( ( log(2,n) ) ) )

reduceOddOneLeft(n::Int) = begin
    len = intLog2( n ) + 1
    n $ makeBinAltLeft( len )
end

reduceOddOneRight(n::Int) = begin
    len = int( floor( ( log(2,n) ) ) ) + 1
    removePowers( n $ makeBinAltRight( len ) )
end

removePowers( n::Int )  = begin
    if n == 0
        return 0
    end
    res = n
    while res % 2 == 0
        res = div( res, 2 )
    end
    res
end


reduceOddLeft( n::Int ) = begin
    if n == 0
        0
    else
        1 + reduceOddLeft( reduceOddOneLeft( n ) )
    end
end

reduceOddRight( n::Int ) = begin
    #@hiya n
    if n == 0
        0
    else
        1 + reduceOddRight( reduceOddOneRight( n ) )
    end
end

reduceRightEdge( ns::String ) = begin
    res, l = reduceLeftEdge( reverse( ns ) )
    reverse( res ), l
end

reduceLeftEdge( ns::String ) = begin
    # @hiya ns
    len = 0
    lenns =  length( ns )
    msk = makeAltStr( ns[1], lenns )
    while len < lenns && ns[len+1] == msk[len+1] 
        # @hiya len
        len += 1
    end
    ns[len+1:end], len
end


makeAltStr( s::Char, len::Int ) = begin
    l = div( len, 2 )
    rep = s == '1' ? "10" : "01"
    term = len % 2 == 0 ? "" : string( s )
    repeat( rep, l ) * term
end

reduceSymEdges( ns::String, side::Char ) = begin
    if length( ns ) == 0
        return Int[]
    end
    if side == 'r'
        cur, l = reduceRightEdge( ns )
    else
        cur, l = reduceLeftEdge( ns )
    end
    #newside = side == 'r' ? 'l' : 'r'
    newside = 'r'
    res = reduceSymEdges( cur, newside )
    push!( res, l )
    res
end

sideFun( s::Char ) = s == 'l' ? x->x : reverse
oppSide( s::Char ) = s == 'l' ? 'r' : 'l'

removeZeros( str::String, side::Char ) = begin
    #@hiya str
    if side == 'b'
        return removeZeros( removeZeros( str, 'l' )[1], 'r' )[1]
    end
    if length( str ) == 0
        return "", 0
    end
    fun = sideFun( side )
    cur = fun( str )
    cur, len = ( cur[1] == '1' ) ? (cur,-1) : removeZeros( cur[2:end], 'l' )
    # len = 0
    # while cur[1] == '0'
    #     cur = cur[2:end]
    #     len += 1
    # end
    fun( cur ), len + 1
end

xorStr( s1::String, s2::String ) = begin
    #@hiya (s1,s2)
    len = min( length( s1 ), length( s2 ) )
    fun( i::Int ) = 
    if i <= len
        ( s1[i]==s2[i] ? "0" : "1" ) * fun( i + 1 )
    else
        ""
    end
    fun(1)
end

eraseProg( s1::String, s2::String ) = begin
    if length( s1 ) == 0
        return "", 0
    end
    cur, len = ( s1[1] == s2[1] ) ? eraseProg( s1[2:end], s2[2:end] ) : (s1,-1)
    cur, len + 1
end

negStr( s1::String ) = begin
    replace( s1, r"[01]", x -> ( x=="0" ? '1' : '0' ) )
end

halfFT( ns::String, side::Char ) = begin
    fun = sideFun( side )
    #@hiya ns
    len = length( ns )
    if len == 0
        return Int[]
    end
    msk = fun( makeAltStr( '1', len ) )
    #@hiya msk
    cur = xorStr( ns, msk )
    cur, scr = removeZeros( cur, side )
    # shave other side?
    # cur, = removeZeros( cur, oppSide( side ) )
    res = halfFT( cur, side )
    push!( res, scr )
    fun( res )
end

fillWPat( pat::String, len::Int ) = begin
    patlen = length( pat )
    l = div( len, patlen )
    rem = len % patlen
    repeat( pat, l ) * pat[1:rem]
end

iterPatDbl( pat::String ) = begin
    dblarr = [ repeat( string(pat[i]), 2 ) for i=1:length(pat) ]
    reduce( (*), dblarr )
end

iterPatOne( pat::String ) = begin
    spt = pat[1] == '1' ? "10" : "01"
    rpt = iterPatDbl( spt )
    replace( pat, spt, rpt )
end

iterPatInv( pat::String ) = begin
    spt = pat[1] == '1' ? "10" : "01"
    replace( pat, spt, "" )
end


genHalfFT( ns::String, pat::String, side::Char ) = begin
    fun = sideFun( side )
    cur = fun( ns )
    pator = fun( pat )
    #@hiya pator
    #@hiya ns
    len = length( ns )
    if len == 0
        return Int[]
    end
    msk = fillWPat( pator, len )
    #@hiya msk
    curmsk, len = eraseProg( cur, msk )
    # if len > 0
    #     #newpat = iterPatOne( pator ) 
    #     #newpat = length( pator ) == 2 ? "1100" : "10"
    #     newpat = length( pator ) == 2 ? "1" : "10"
    #     #if evenQ( len )  newpat = reverse( newpat )  end
    # else
    #     newpat = negStr( pator )
    # end
    newpat = len==0 ? negStr(pator) : pator
    res = genHalfFT( curmsk, newpat, 'l' )
    if ( len > 0 ) push!( res, len ) end
    res        
end

function waveMetric( str::String )
    genHalfFT( str, "01", 'r' )
end

function carries( wm::Array{Int} )
    ch = 1
    i = 1
    res = Int[]
    while i <= length( wm )
        push!( res, ch )
        if evenQ( wm[i] )  ch = 1-ch end
        i += 1
    end
    reverse( res )
end

function carriesAll( str::String )
    strrv = reverse( str )
    ch = 1
    i = 1
    res = Int[]
    look = '0'
    len =  length( strrv )
    while i <= len
        push!( res, ch )
        if i < len && strrv[i] == look && strrv[i+1] == look
            ch = 1-ch 
            look = look=='0' ? '1' : '0'
        end
        i += 1
    end
    reverse( res )
end
   

function cwMetric( str::String )
    wm = waveMetric( str )
    cr = carries( reverse(wm) )
    sgns = 2 * cr .- 1
    #[ reverse( wm )', reverse( sgns )' ]
    ( wm .* sgns )'
end

function cwallColor( str::String )
    ca = carriesAll( str )
    con = 0.1
    val( pos::Int ) = int( string( str[pos] ) )
    [ ca[i] - (2*ca[i]-1)*0.1 * val(i)  for i=1:length(str) ]
end

trimEdges( ns::String ) = begin
    len = length( ns )
    msk = fillWPat( "10", len )
    res, left = eraseProg( ns, msk )
    res, right = eraseProg( reverse( res ), msk )
    [ left, length( res ), right ]
end
    
sumOnes( str::String ) = reduce( (+), [ int(string(c)) for c in str ]
                                )
countXO( str::String ) = begin
    len = length( str ) 
    nXs = sumOnes( str ) 
    nOs = len - nXs
    nOs, nXs
end


function randLargeOdd( n::Int )
    rndarr = [ ( rand() < 0.5 ? "0" : "1" ) for i=1:n ] 
    rndstr = reduce( *, rndarr ) * "1"
    oddFactor( parseint( BigInt, rndstr, 2 ) )
end

         
