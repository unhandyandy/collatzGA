#-*- julia; eval: (auto-complete-mode) -*-

include( "mydefs.jl" )

#using Winston

type cState
    bit::Bool
    carry::Bool
    stop::Bool
    channel::Bool
end

stop = cState( false, true, true, false )

function (==)( cs1::cState, cs2::cState )
    cs1.bit==cs2.bit && cs1.carry==cs2.carry &&
    cs1.stop==cs2.stop && cs1.channel==cs2.channel
end

function cState( b::Bool, c::Bool )
    cState( b, c, false, false )
end

function cState( b::Bool )
    cState( b, false )
end

function tostring( cs::cState )
    ( cs.bit ? "1" : "0" ) * ( cs.channel ? ( cs.carry ? "+" : "-" ) : " " )
end

function tocolor( cs::cState )
    if cs.stop
        -1
    else
        ( cs.channel ? 1 : 0 ) + ( cs.channel && cs.carry ? 1 : 0 ) + ( cs.channel ? -1 : 1 ) * 0.2 * ( cs.bit ? 1 : 0 )
    end
end

function nextBit( csLeft::cState, cs::cState, csRight::cState )
    if cs.stop
        false
    else
        csLeft.bit $ cs.bit $ csRight.carry
    end
end

function curCarry( csLeft::cState, cs::cState, csRight::cState )
    ( csLeft.bit && cs.bit ) || (  csLeft.bit && csRight.carry ) ||
    ( cs.bit && csRight.carry )
end

function curChannel( csLeft::cState, cs::cState, csRight::cState )
    csLeft.bit != cs.bit || cs.bit != csRight.bit
end

function updateRowBits( row::Array{cState,1}, grow=false )
    currow = vcat( [ cState(false), cState(false) ], copy( row ), [ stop ] )
    len = length( row )
    newrow = [ nextBit( currow[i-1],currow[i],currow[i+1] ) for i=2:len+2 ]
    if !grow || !newrow[1]
        newrow = newrow[2:end]
    end
    map( cState, newrow )
end

function addRowCarries( row::Array{cState,1} )
    currow = vcat( [cState(false)], copy( row ), [ stop ] )
    len = length( row )
    for i = len + 1 : -1 : 2
        if currow[i+1].stop && !currow[i].bit
            currow[i] = stop
        else
            currow[i].carry = curCarry( currow[i-1],currow[i],currow[i+1] )
            currow[i].channel = curChannel( currow[i-1],currow[i],currow[i+1] )
        end
    end
    currow[2:len+1]
end

function nextRow( row::Array{cState,1}, grow=false )
    addRowCarries( updateRowBits( row, grow ) )
end

function nextRowExact( row::Array{cState,1} )
    newrow = nextRow( row, true )
    prepRow( newrow )
end
   

function rowBoolToCstate( row::Array{Bool,1} )
    addRowCarries( map( cState, row ) )
end

function rowStrToCstate( str::String )
    brow = [ c=='1' for c in str ]
    rowBoolToCstate( brow )
end

function rowCstateToStr( row::Array{cState,1} )
    strarr = [ ( cs.bit ? "1" : "0" ) for cs in row ]
    reduce( *, strarr )
end

function flow( row::Array{cState,1}, tm::Int )
    res = Array{cState,1}[]
    push!( res, row )
    for i=1:tm
        #@hiya rowCstateToStr( res[i] )
        newrow = nextRow( res[i] )
        #@hiya rowCstateToStr( newrow )
        push!( res, newrow )
    end
    [ res[i][j] for i=tm+1:-1:1, j=1:length(row) ]
end

function arrToImage1( arr::Array{cState,2} )
    convert( Array{Float64,2}, mapMat( tocolor, arr ) )
end

function showFlowCA( str::String, time = -1 )
    wid = length( str )
    len = ( time > 0 ) ? time : int( wid / 0.4 ) + 1
    #cols =  wid + int( 0.6 * len )
    #right = int( 0.6 * wid )
    left = int( 0.6 * len )
    startstr = repeat( "0", left ) * str
    flowarr = flow( rowStrToCstate( startstr ), len )
    #printStrArr( map( parseToStr, parseFlow( flowarr ) ) )
    flowmat = arrToImage1( flowarr )
    display( imagesc( flowmat, (-1,2) ) )
end

function eraseProg( s1::Array{cState,1}, pat::Array{cState,1} )
    #@hiya pat
    unit = length( pat )
    if length( s1 ) < unit
        return s1, 0
    end
    cur, len = ( s1[1:unit] == pat ) ? eraseProg( s1[ unit + 1 : end ], pat ) : (s1,-1)
    cur, len + 1
end

function neg( cs::cState )
    cState( !cs.bit, cs.carry, cs.stop, cs.channel )
end

function stripStoppers( ns::Array{cState,1} )
    if length( ns ) == 0 ; return ns end
    ns[1].stop ? stripStoppers( ns[2:end] ) : ns
end

function stripZeros( ns::Array{cState,1} )
    if length( ns ) == 0 ; return ns end
    !ns[1].bit ? stripZeros( ns[2:end] ) : ns
end

function stripZeros( str::String )
    if length( str ) == 0 ; return str end
    str[1]=='0' ? stripZeros( str[2:end] ) : str
end


function prepRow(  ns::Array{cState,1} )
    res = reverse( stripZeros( ns ) )
    reverse( stripZeros( stripStoppers( res ) ) )
end

function prepRow( str::String  )
    res = reverse( stripZeros( str ) )
    reverse( stripZeros( res ) )
end


function parseL( ns::Array{cState,1},
               pat::Array{cState,1}=[ cState( true, true ),
                                     cState( false, true ) ],
               flipped=false )
    len = length( ns )
    if len == 0
        return String[]
    end
    carry = pat[1].carry
    unit = length( pat )
    newns, len = eraseProg( ns, pat )
    if len==0
        if flipped
            newflipped = false
            if unit == 2
                newpat = [ cState( true, carry ) ]
            else
                newpat = [ cState( true ), cState( false ) ]
                for cs in newpat
                    cs.carry = !carry
                end
            end
        else
            newflipped = true
            newpat = map( neg, pat )
        end
    else
        newpat = [ cState( true, true ), cState( false, true ) ]
        newflipped = false
    end
    #@hiya length( newns )
    res = parseL( newns, newpat, newflipped )
    bitstrarr = reverse( map( x -> ( x.bit ? "1" : "0" ), pat ) )
    bitstr = reduce( *, bitstrarr ) * ( carry ? "+" : "-" )
    newstr = bitstr * @sprintf( "%2s", len * unit )
    if ( len > 0 ) push!( res, newstr ) end
    res
end

function parseR( ns::Array{cState,1} )
    parseL( prepRow( ns ) )
end

#function unparse( prs::Array{String,1} )

function parseFlow( flow::Array{cState,2} )
    r,c = size( flow )
    [ parseR( [ flow[i,j] for j=1:c ] ) for i=r : -1 : 1 ]
end

function parseToStr( prs::Array{String,1} )
    res = ""
    for p in prs
        res *= @sprintf( "%7s", p )
    end
    #reduce( (x,y) -> ( x * ", " * y ), prs )
    res
end

function printStrArr( arr::Array{Array{String,1},1} )
    # len = maximum( map( length, arr ) )
    # frm = @sprintf( "%%%ds", 7 * len )
    for s in arr
        @printf( "%80s \n", printToStr( s ) )
    end
end

function mutate( str::String, p=(1/length(str)), q=p )
    #println( "mutating..." )
    res = replace( str, r"(.)", x -> ( rand() > p ? x : ( x=="1" ? "0" : "1" ) ) )
    pos = rand( 1:length( str ) )
    if rand() < 2 * q
        res = res[1:pos-1] * res[pos+1:end]
    end
    left = ( rand() < q ) ? "1" : ""
    right = ( rand() < q ) ? "1" : ""
    left * res * right
end

function findOrder( str::String, trunc::Bool=false )
    findOrder( rowStrToCstate( str ); trunc=trunc )
end

function findOrder( row::Array{cState,1}, cnt=0, 
                   hist=Set{Array{cState,1}}(); trunc::Bool=false )
    #trimmed = prepRow( row )
    if length( row ) < ( trunc ? 63 : 2 )
        cnt
    else
        if row in hist
            println( "Cycle!!!!" )
            Inf
        else
            findOrder( nextRowExact( row ), cnt + 1, trunc=trunc )
        end
    end
end



function procreate( str::String )
    beststr = str
    bestord = findOrder( beststr )
    newstr = str
    neword = bestord
    while 0 < neword <= bestord
        cand = beststr
        while cand == beststr
            cand = mutate( beststr )
        end
        newstr = cand
        neword = findOrder( newstr )
        display( oplot() )
    end
    showFlowCA( newstr, neword )
    #sleep(0.5)
    println( "good offspring: " * newstr )
    @printf( "order: %d \n", neword )
    newstr
end


function progress( str::String )
    newstr = str
    while true
        newstr = procreate( newstr )
    end
end


function haveSex( str1::String, str2::String )
    len1, len2 = length( str1 ), length( str2 )
    len = max( len1, len2 )
    function elm( i::Int )
        if i <= len1 && i <= len2
            ( rand() < 0.5 ) ? str1[i] : str2[i]
        elseif i <= len1
            str1[i]
        else
            str2[i]
        end
    end
    babyarr = [ string(elm( i )) for i=1:len ]
    baby = reduce( *, babyarr )
    prepRow( mutate( baby, 0.01 ) )
end

function haveSex1( row1::Array{cState,1}, row2::Array{cState,1} )
    gen1, gen2 = splitR( row1 ), splitR( row2 )
    len1, len2 = length( gen1 ), length( gen2 )
    len = max( len1, len2 )
    function elm( i::Int )
        res = ""
        if i <= len1 && i <= len2
            res = ( rand() < 0.5 ) ? geneToCSrow( gen1, i, row1 ) : geneToCSrow( gen2, i, row2 )
        elseif i <= len1
            res = geneToCSrow( gen1, i, row1 )
        else
            res = geneToCSrow( gen2, i, row2 )
        end
        reduce( *, map( x -> x.bit ? "1" : "0", res ) )
    end
    babyarr = [ elm( i ) for i=len:-1:1 ] 
    baby  = reduce( *, babyarr )
    prepRow( mutate( baby, 0.01 ) )
end



function orgy!{T<:String}( fitlst::Dict{ T, Float64 } )
    pmonster = 0.1
    maxfit = maximum( values( fitlst ) )
    strs = collect( keys( fitlst ) )
    ppl = length( strs )
    newfit = 0
    minstr, minfit = findMin( fitlst )
    newstr, neword = "", 0
    while 0 <= newfit <= minfit || newstr in strs
        mate1 = strs[ rand( 1:ppl ) ]
        mate2 = strs[ rand( 1:ppl ) ]
        matefun = ( rand() > pmonster ) ? haveSex1 : haveSex2
        newstr = matefun( rowStrToCstate(mate1), rowStrToCstate(mate2) )

        #neword = findOrder( newstr, true )
        newfit, neword = fitnessTruncOrder( newstr )
        display( oplot() )
    end
    fitlst[ newstr ] = newfit
    delete!( fitlst, minstr )
    showFlowCA( newstr, neword )
    if newfit > maxfit ; maxfit = newfit end
    #print( "\n good offspring: " )
    showGenome( newstr )
    @printf( "max: %.4f                              ", maxfit )
    @printf( "fitness: %.4f > %.4f \n", newfit, minfit )
    #println( fitlst )
end

function findMin{T,N<:Real}( dct::Dict{ T, N } )
    min = Inf
    local key::T
    for s in collect( keys( dct ) )
        val = dct[ s ]
        if val < min
            min = val
            key = s
        end
    end
    key, min
end

function evolution( adam::String )
    bestHist = [ adam ]
    evolution( bestHist )
end

function evolution{T<:String}( eden::Array{T,1} )
    # bestHist = copy( eden )
    # newstr = eden[end]
    # neword = findOrder( newstr )
    fitD = Dict{ T, Float64 }()
    for str in eden
        fit = findOrder( str ) / length( str )
        fitD[str] = fit
    end
    while true
        orgy!( fitD )
    end
end

function evolution!{T<:String}( fitD::Dict{ T, Float64 }, minpopsize=30 )
    #minpopsize = 30
    strs = keys( fitD )
    len = maximum( map( length, strs ) )
    while length( keys( fitD ) ) < minpopsize
        newstr = "1" * randBinStr( len - 2 ) * "1"
        fitD[ newstr ] = findOrder( newstr ) / len
    end
    while true
        orgy!( fitD )
    end
end

function randBinStr( len::Int )
    res = repeat( "0", len )
    replace( res, r"(.)", x -> ( rand() < 0.5 ? "0" : "1" ) )
end

function evolution( len::Int, eden::Array{String,1}=String[] )
    popsize = 30
    while length( eden ) < popsize
        push!( eden, "1" * randBinStr( len - 2 ) * "1" )
    end
    evolution( eden )
end

function monster( str1::String, str2::String )
    pos = rand( 1:length( str1 ) )
    str3 = str1[1:pos] * str2
    haveSex( str1, str3 )
end

function fitnessTruncOrder( str::ASCIIString )
    len = length( str )
    ord = findOrder( str, true )
    len > 62 ? ord / ( len ) : 0 , ord
end   

function enterDict{T<:String}( str::T, dict::Dict{T,Float64} )
    dict[ str ], _ = fitnessTruncOrder( str )
end

function enterDict{T,N<:Real}( genlst::T, dict::Dict{T,N} )
    dict[ genlst ] = fitnessOrth( genlst )[1]
end

function saveEden( sym::Symbol )
    name, val = string( sym ), eval( sym )
    open( x -> write( x, @sprintf( "%s = %s", name, val ) ),
          @sprintf("%s.jl", name ), "w" )
end



type gene
    kind::Char
    length::Int
end

function splitL( row::Array{cState,1}, kinder::Function, 
                res::Array{gene,1}=gene[], prev::Char=' ' )
    if length( row ) == 0 
        return res
    end
    cur = row[1]
    curkind = kinder( cur )
    if curkind == prev
        curgene = res[end]
        curgene.length += 1
        newprev = prev
    else
        newgene = gene( curkind, 1 )
        push!( res, newgene )
        newprev = curkind
    end
    #@hiya ( cur )
    splitL( row[2:end], kinder, res, newprev )
end

function channelKind( cs::cState ) 
    cs.channel ? 'a' : ( cs.bit ? 'b' : 'c' )
end

function splitR( row::Array{cState,1}, kinder::Function=channelKind )
    splitL( prepRow( row ), kinder )
end

function geneToCSrow( prs::Array{gene,1}, gnum::Int, row::Array{cState,1} )
    revgen, revrow = reverse( prs ), reverse( row )
    st = 1 + reduce( +, map( x->x.length, revgen[1:gnum-1] ) )
    fn = st + revgen[ gnum ].length - 1
    reverse( revrow[st:fn] )
end

function showGenome( genome::Array{cState,1} )
    prs = splitR( genome )
    println()
    for i = length( prs ) : -1 : 1
        strgen = geneToCSrow( prs, i, genome ) |> rowCstateToStr 
        print( strgen * " " )
    end
    println()
end

function showGenome( str::ASCIIString )
    str |> rowStrToCstate |> showGenome
end

function findLastStop( row::Array{cState,1}, res=0 )
    #@hiya row
    if length( row ) != 0 && row[end].stop 
        findLastStop( row[ 1 : end - 1 ], res + 1 )
    else
        res
    end
end

function calcStopSeq( flow::Array{cState,2} )
    r,c = size( flow )
    [ findLastStop( vec( flow[ i, 1:end ] ) ) for i = 1 : r ] |> reverse
end

function stopSeqToGenome( seq::Array{Int,1}, gen::Array{Int,1}=Int[], prev::Int=-1 )
    if length( seq ) == 0
        gen[ 2 : end ]
    else
        cur = seq[1]
        if cur == prev
            gen[end] += 1
        else
            push!( gen, cur - prev - 1 )
            push!( gen, 0 )
        end
        stopSeqToGenome( seq[2:end], gen, cur )
    end
end

function strToOrthGenome( str::ASCIIString )
    rowCstateToOrthGenome( rowStrToCstate( str ) )
end    

function rowCstateToOrthGenome( row::Array{cState,1} )
    len = length( row )
    flow( row, len ) |> calcStopSeq |> stopSeqToGenome
end    

function calcLastRowFromGenome( gen::Array{Int,1} )
    len = length( gen )
    d = 0
    for i = 2 : 2 : len
        d += 1 + gen[i]
    end
    str = "1" * repeat( "0", d )
    rowStrToCstate( str )
end

function nextVert( row::Array{cState,1}, v::Int )
    [ row for i = 1 : v ]
end

function nextHoriz( row::Array{cState,1}, h::Int )
    newrow = vcat( row, cState[ stop for i = 1 : h + 1 ] )
    [ newrow for i = 1:1 ]
end

initRow = [ [ cState(true,true,false,true) ] for i=1:1 ]

function startRowsFromGenome( 
             gen::Array{Int,1}, 
             rows::Array{Array{cState,1},1}=initRow, 
             nextFun::Function=nextVert )
    if length( gen ) == 0 
        rows
    else
        newnext = nextFun==nextVert ? nextHoriz : nextVert
        funapp = nextFun( rows[end], gen[1] )
        typeof( funapp )
        newrows = [ rows, funapp ]
        newgen = gen[ 2 : end ]
        startRowsFromGenome( newgen, newrows, newnext )
    end
end

function completePredecessor!( row::Array{cState,1}, pred::Array{cState,1} )
    #@hiya pred
    lenp, lenr = length( pred ), length( row )
    if lenr >= lenp
        pos = lenr - lenp + 1
        cs = cState( false )
        cs.bit = row[pos].bit $ pred[1].carry $ pred[1].bit
        cs.carry = ( cs.bit && pred[1].bit ) ||
                   ( cs.bit && pred[1].carry ) ||
                   ( pred[1].bit && pred[1].carry )
        unshift!( pred, cs )
        completePredecessor!( row, pred )
    end
end

function completeAllPredecessors( rows::Array{Array{cState,1},1} )
    res = copy( rows )
    for i = length( rows ) : -1 : 2
        completePredecessor!( res[ i ], res[ i - 1 ] )
    end
    res
end

function startPred( currow::Array{cState,1}, gen::Array{Int,1}, horiz::Bool )
    len = horiz ? pop!( gen ) + 1 : 0
    change = true
    if horiz 
        if gen[end] == 0
            pop!( gen )
            change = false
        end
    else 
        gen[end] -= 1
        if gen[end] == 0
            pop!( gen )
        else
            change = false
        end
    end
    newcur = [ currow..., [ stop for i = 1 : len ]... ]
    #@hiya newcur
    newcur, copy( initRow[1] ), horiz $ change
end

function quickGrow( gen::Array{Int,1} )
    genc = copy( gen )
    hor = evenQ( length( gen ) )
    pred = copy( initRow[1] )
    local currow::Array{cState,1}
    while length( genc ) > 0
        currow, pred, hor = startPred( pred, genc, hor )
        completePredecessor!( currow, pred )
        #currow = pred
    end
    prepRow( currow )
end

function quickFitness( gen::Array{Int,1} )
    _, time = orthGenSlope( gen )
    row1 = addRowCarries( quickGrow( gen ) )
    row2 = iterate( nextRowExact, time, row1 )
    len = max( length( row1 ), length( row2 ) )
    matchLength( reverse(row1), reverse(row2) ) / len, row1, time
end



function growFromOrthGenome( gen::Array{Int,1} )
    prepRow( ( gen |> startRowsFromGenome |> completeAllPredecessors )[1] )
end

function mate( p1::Array{Int,1}, p2::Array{Int,1} )
    len1, len2 = length( p1 ), length( p2 )
    lenmin, lenmax = min( len1, len2 ), max( len1, len2 )
    big = ( len1 == lenmax ) ? p1 : p2
    function elm( i::Int )
        if i <= lenmin
            ( rand() < 0.5 ) ? p1[i] : p2[i]
        else
            big[i]
        end
    end
    Int[ elm( i ) for i=1:lenmax ]
end

function mutate!( genlst::Array{Int,1}, p=0.5 )
    if rand() >= p
        len = length( genlst )
        pos = rand( 0 : len + 1 )
        #@hiya pos
        #res = copy( genlst )
        if pos == 0
            genlst = genlst[3:end]
        elseif pos <= len
            genlst[ pos ] += ( rand() < 0.5 ) ? -1 : 1
            genlst[ pos ] = max( 0, genlst[pos] )
        else
            push!( genlst, 0 )
        end
    end
end

function haveSex2( row1::Array{cState,1}, row2::Array{cState,1} )
    p1, p2 = map( rowCstateToOrthGenome, ( row1, row2 ) )
    offgen = mate( p1, p2 )
    mutate!( offgen )
    off = growFromOrthGenome( offgen )
    offstr = rowCstateToStr( off )
end

function fitnessOrth( genome::Array{Int,1} )
    #rowMetric = (x,y) -> - distRow(x,y)
    rowMetric = matchLength
    fllst = ( genome |> startRowsFromGenome |> completeAllPredecessors )
    len = length( fllst )
    row1 = addRowCarries( prepRow( fllst[1] ) )
    #@hiya rowCstateToStr( row1 )
    row2 = iterate( nextRowExact, len - 1, row1 )
    #fullflow = flow( row1, len - 1 )
    #row2 = prepRow( vec( fullflow[1,1:end] ) )
    #@hiya rowCstateToStr( row2 )
    rowMetric( reverse(row1), reverse(row2) ), row1, len
end

function distCS( cs1::cState, cs2::cState )
    ( cs1.bit==cs2.bit ? 0 : 1 ) +
    ( cs1.carry==cs2.carry ? 0 : 1 ) +
    ( cs1.channel==cs2.channel ? 0 : 1 )
end
    
function distRow( row1::Array{cState,1}, row2::Array{cState,1} )
    difcon = 4
    len1, len2 = length( row1 ), length( row2 )
    lenmin, lenmax = sort( [ len1, len2 ] )
    dist = difcon * ( lenmax - lenmin )
    #revrow1, revrow2 = reverse( row1 ), reverse( row2 )
    for i = 1 : lenmin
        dist += distCS( row1[i], row2[i] )
    end
    dist
end

function matchLength( row1::Array{cState,1}, row2::Array{cState,1}, cnt::Int=0 )
    length( row2 )==0 || length( row2 )==0 || row1[1]!=row2[1] ? cnt : matchLength( row1[2:end], row2[2:end], cnt + 1 )
end

function demographics{T<:Real,C}( fitD::Dict{ Array{C,1}, T } )
    maxfit = maximum( values( fitD ) )
    ppl = length( fitD )
    mean = calcMean( fitD )
    ppl, mean, maxfit
end 


function orgy!{T<:Real,C}( fitD::Dict{ Array{C,1}, T } )
    pmonster = 0.0
    #fitfun = fitnessOrth
    #fitfun = coreFitness
    fitfun = quickFitness
    #maxfit = maximum( values( fitD ) )
    keycol = collect( keys( fitD ) )
    #ppl = length( keycol )
    newfit = -Inf
    #minfit = calcMean( fitD )
    ppl, minfit, maxfit = demographics( fitD )
    local newkey::Array{C,1}, newrow::Array{cState,1}, len::Int, newrowlen::Int=0, newkeycand::Array{Int,1}
    while newfit <= minfit || newkey in keycol || newrowlen < 63
        mate1 = keycol[ rand( 1:ppl ) ]
        mate2 = keycol[ rand( 1:ppl ) ]
        #minfit = ( fitD[ mate1 ] + fitD[ mate2 ] ) / 2
        #matefun = ( rand() > pmonster ) ? mate : mate
        #matefun = coreMate
        matefun = mate
        newkey = matefun( mate1, mate2 )
        #@hiya typeof(newkey)
        #newkey = correctOrthGen( newkeycand )
        #mutateFun = coreMutate!
        mutateFun = mutate!
        mutateFun( newkey )
        newfit, newrow, len = fitfun( newkey )
        newrowlen = length( newrow )
        #display( oplot() )
    end
    #@hiya newkeycand'
    #@hiya newkey'
    fitD[ newkey ] = newfit
    #delete!( fitD, minkey )
    #showFlowCA( rowCstateToStr(newrow), len )
    if newfit > maxfit ; maxfit = newfit end
    #@hiya typeof( newkey )
    #printList( newkey )
    @printf( "max: %.4f                 ", maxfit )
    @printf( "row length: %d;    fitness: %.4f > %.4f \n", newrowlen, newfit, minfit )
end





function pOrgy!{T<:Real,C}( fitD::Dict{ Array{C,1}, T } )
    fitfun = quickFitness
    keycol = collect( keys( fitD ) )
    newfit = -Inf
    ppl, minfit, maxfit = demographics( fitD )
    local newkey::Array{C,1}, newrow::Array{cState,1}, newrowlen::Int=0

    nproc = nprocs()
    offlst = Array{C,1}[]
    for i = 1 : nproc
        mate1 = keycol[ rand( 1:ppl ) ]
        mate2 = keycol[ rand( 1:ppl ) ]
        matefun = mate
        newkey = matefun( mate1, mate2 )
        mutateFun = mutate!
        mutateFun( newkey )
        push!( offlst, newkey )
    end

    fitlst = pmap( fitfun, offlst )

    for i = 1 : nproc
        newfit, newrow, len = fitlst[ i ]
        newrowlen = length( newrow )
        if newfit > minfit
            fitD[ newkey ] = newfit
            if newfit > maxfit ; maxfit = newfit end
            @printf( "max: %.4f                 ", maxfit )
            @printf( "row length: %d;    fitness: %.4f > %.4f \n", newrowlen, newfit, minfit )
        end
    end
end





function evolution!{T<:Real}( fitD::Dict{ Array{Int,1}, T }, minpopsize=32 )
    fitfun = fitnessOrth
    randkey = ( x -> Int[ rand( 0:9 ) for i = 1 : x ] )
    keytuple = keys( fitD )
    #len = maximum( map( length, keytuple ) )
    while length( keys( fitD ) ) < minpopsize
        newkey = correctOrthGen( randkey( rand( 40:80 ) ) )
        fitD[ newkey ] = fitfun( newkey )[1]
    end
    while true
        orgy!( fitD )
    end
end

function printList{T<:Number}( seq::Array{T,1} )
    strlst = map( x -> string(x) * " ", seq )
    println( "\n" * reduce( *, strlst ) )
end
function printList( seq::Array{Bool,1} )
    strlst = map( x -> x ? "1" : "0", seq )
    println( "\n" * reduce( *, strlst ) )
end


function makeRandDict( len::Int, pop::Int )
    newdict = Dict{ ASCIIString, Float64 }()
    for i = 1 : pop
        enterDict( "1" * randBinStr( len ) * "1", newdict )
    end
    newdict
end

targetSlope = log( 2, 3 ) - 1

function orthGenSlope( genlst::Array{Int,1} )
    vrt, hrz = 0, 0
    len = length( genlst )
    for i = 1 : 2 : len
        vrt += genlst[ i ] + 1
    end
    for i = 2 : 2 : len
        hrz += genlst[ i ] + 1
    end
    hrz, vrt
end

function correctOrthGen( genlst::Array{Int,1} )
    res = copy( genlst )
    len = length( genlst )
    hrz, vrt = orthGenSlope( genlst )
    hrztarg = floor( targetSlope * vrt )
    hrztarg += hrz==hrztarg ? 0 : ( hrz < hrztarg ? 0 : 1 )
    vrttarg = floor( hrz / targetSlope )
    vrttarg += vrt==vrttarg ? 0 : ( vrt < vrttarg ? 0 : 1 ) 
    #@hiya hrztarg, vrttarg
    if hrz != hrztarg
        #@hiya hrz, vrt
        if evenQ( len ) 
            if res[end] + hrztarg - hrz >= 0
                res[end] += hrztarg - hrz
            else
                push!( res, vrttarg - vrt - 1 )
            end
        else
            if res[end] + vrttarg - vrt >= 0
                res[end] += vrttarg - vrt
            else
                push!( res, hrztarg - hrz - 1 )
            end
        end
    end
    res
end


function discreteTarget( len::Int )
    vrt, hrz = 0, 0
    res = Int[]
    for i = 1 : len
        if evenQ( i )
            newhrz = ( hrz < targetSlope * vrt ) ? 1 : 0
            push!( res, newhrz )
            hrz += newhrz + 1
        else
            newvrt = ( hrz < targetSlope * vrt ) ? 0 : 1
            push!( res, newvrt )
            vrt += newvrt + 1
        end
    end
    res
end

function randOrthGen( len::Int, lim::Int, correct::Bool=true )
    res = Int[]
    lowlim = round( targetSlope * lim )
    for i = 1 : len
        push!( res, rand( 0:lim ) )
        push!( res, rand( 0:lowlim ) )
    end
    correct ? correctOrthGen( res ) : res
end

function randOrthEden( len::Int, lim::Int, pop::Int )
    fitnessFun = quickFitness
    res = Dict{ Array{Int,1}, Int }()
    dt = discreteTarget( len )
    res[ dt ] = fitnessFun( dt ) [1]
    for i = 1 : pop - 1
        newgen = randOrthGen( len, lim )
        res[ newgen ] = fitnessFun( newgen )[1]
    end
    res
end



function calcMean{T,N<:Real}( dct::Dict{T,N} )
    sum( values( dct ) ) / length( keys( dct ) )
end

function cull!{T,N<:Real}( dct::Dict{T,N} )
    mn = calcMean( dct )
    for (k, v) in dct
        if v < mn
            delete!( dct, k )
        end
    end
end


function showOrthList( genlst::Array{Int,1} )
    fit, row0, len = fitnessOrth( genlst )
    println( "fitness: ", fit )
    showFlowCA( rowCstateToStr( row0 ), len - 1 )
end

### core functions

function coreStrToCstate( str::String )
    len = length( str )
    res = [ cState( str[i] == '1' ) for i = 1 : len ]
    #res[1].carry = res[1].bit $ res[2].bit
    res
end

function coreAddCarries!( row::Array{cState,1} )
    len = length( row )
    for i = len - 1 : -1 : 1
        row[i].carry = ( row[ i ].bit && row[ i + 1 ].bit ) || 
                       ( row[ i ].bit && row[ i + 1 ].carry ) || 
                       ( row[ i + 1 ].bit && row[ i + 1 ].carry )
    end
end

function coreNextBits( row::Array{cState,1}, right::cState )
    res = cState[ right ]
    revrow = reverse( row )
    for i = 2 : length( row ) - 1
        newsc = cState( revrow[i].bit $ revrow[i].carry $ revrow[ i + 1 ].bit )
        push!( res, newsc )
    end
    push!( res, cState( revrow[end].bit $ revrow[end].carry ) )
    if revrow[end].bit && revrow[end].carry
        push!( res, cState( true ) )
    end
    reverse( res )
end 

function coreNextRCarry!( row::Array{cState,1}, nextRight::cState )
    left = length( row ) > 1 ? row[ end - 1 ].bit : false
    row[end].carry = left $ row[end].bit $ nextRight.bit
end

function coreGrow( core::Array{cState,1} )
    len = length( core )
    res = Array{cState,1}[ cState[ core[1] ] ]
    #@hiya res
    for i = 1 : len - 1
        newright = core[ i + 1 ]
        coreNextRCarry!( res[i], newright )
        coreAddCarries!( res[i] )
        newrow = coreNextBits( res[ i ], newright )
        push!( res, newrow )
    end
    res
end

function coreGenToRow0Str( core::Array{cState,1} )
    coreGrow( core )[end] #|> rowCstateToStr
end

function coreFitness( coreB::Array{Bool,1} )
    core = coreFromBoolList( coreB )
    row0 = coreGenToRow0Str( core ) |> addRowCarries
    ord = findOrder( row0; trunc=true )
    len = length( row0 )
    len > 62 ? ord / (  len - 62 ) : 0 , row0, ord
end

    
        
function pEvolve!{T<:Real,C}( fitD::Dict{ Array{C,1}, T }, time::Int=Inf )
    for i = 1 : time
        pOrgy!( fitD )
    end
end

function coreFromBoolList( bl::Array{Bool,1} )
    map( cState, bl )
end

function coreMate( m1::Array{Bool,1}, m2::Array{Bool,1} )
    len1, len2 = length( m1 ), length( m2 )
    lenmin, lenmax = sort( [ len1, len2 ] )
    small, large = len1==lenmin ? (m1,m2) : (m2,m1)
    res = Bool[ true for i = 1 : lenmax ]
    for i = 1 : lenmax
        upper = lenmax - lenmin 
        res[i] = ( rand() < 0.5 || lenmin < i <= upper ) ? large[i] : ( i <= lenmin ? small[i] : small[ end + i - lenmax ] )
    end
    #@hiya typeof( res )
    res
end


function coreRand( len::Int )
    [ true, [ rand() < 0.5 for i = 1 : len - 2 ]..., true ]
end


function coreMutate!( cg::Array{Bool,1}, p::Float64=0.1 )
    if rand() < p
        len = length( cg )
        pos = rand( 2 : len - 1 )
        if rand() < 0.5
            cg[pos] = !cg[pos]
        else
            if rand() < 0.5
                cg = [ cg[ 1 : pos - 1 ]..., cg[ pos + 1, end ]... ]
            else
                cg = [ cg[ 1 : pos ]..., rand() < 0.5, cg[ pos + 1, end ]... ]
            end
        end
    end
end

function coreRandDict( len::Int, pop::Int )
    dct = Dict{ Array{Bool,1}, Float64 }()
    for i = 1 : pop
        newcore = coreRand( len )
        dct[ newcore ] = coreFitness( newcore )[1]
    end
    dct
end

function getQGenFit( gen::Array{Int,1}  )
    res = quickFitness( gen )[1]
    println( "fitness calculated..." )
    res
end

function pAddQGen( n::Int, dct::Dict{ Array{Int,1}, Float64 }, len::Int=4000, lim::Int=5 )
    lst = [ randOrthGen( len, lim, false ) for i = 1 : n ]
    function enter( gen:: Array{Int,1}, fit::Float64 )
        println( "adding entry..." )
        dct[ gen ] = fit
        println( "...done." )
    end
    fits = pmap( getQGenFit, lst )
    println( "...fitnesses all calculated." )
    map( enter, lst, fits )
    println( "All done." )
end

