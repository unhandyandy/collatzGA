#-*- julia; eval: (auto-fill-mode); eval: (auto-complete-mode) -*-

macro hiya(ex)
    quote
        begin
            local val = $ex
            print( $(string(ex)), ": ", $ex, "\n" )
            val
        end
    end
end

macro fn( ex )
    vars = [ "x0" "x1" "x2" "x3" "x4" "x5" "x6" "x7" "x8" "x9" ]
    used = filter( vars ) do v
        typeof( match( Regex(v), string(ex) )) == RegexMatch
    end
    if length( used ) > 0
        ## args = "("
        ## for v in used
        ##     args *= v * ","
        ## end
        ## funstr = "( " * args[1:end-1] * ") -> " * string(ex) * " )"
        ## return parse( funstr )
        vars = map( parse, used )
        args = Expr( :tuple, vars... )
        funex = Expr( :->, args, ex )
        return funex
    else
        return :( x -> $ex )
    end
end

function nestList( f::Function, n::Int, x0 ) 
    t = typeof( x0 )
    res = t[]
    cur = x0
    for i=0:n
        push!( res, cur )
        cur = f( cur )
    end
    return res
end

function evenQ( n::Int ) 
    n % 2 == 0
end

function mapAt( f::Function, x, lev::Int )
    if lev == 0
        f(x)
    else
        map( y -> mapAt( f, y, lev - 1 ), x )
    end
end

function mapMat( f::Function, arr::Array{Any,2} )
    r, c = size( arr )
    [ f( arr[i,j] ) for i=1:r, j=1:c ]
end

# function iterate( f::Function, n::Int, x )
#     n==0 ? x : iterate( f, n - 1, f( x ) )
# end

function iterateList( f::Function, n::Int, x::T, acc::Array{T,1}=T[] ) where {T}
    newacc = acc
    push!( newacc, x )
    n==0 ? newacc : iterateList( f, n - 1, f( x ), newacc )
end

function trimWhiteL( str::AbstractString )
    if str[1] == ' '
        return trimWhiteL( str[2:end] )
    else
        return str
    end
end

function trimWhiteR( str::AbstractString )
    reverse( trimWhiteL( reverse( str ) ) )
end

function trimWhite( str::AbstractString )
    trimWhiteL( str ) |> trimWhiteR
end

function mapset( f::Function, s::Set )
    if length( s ) == 0; return Set() end
    el = first( s )
    t = typeof( f( el ) )
    newset = Set{t}()
    for x in s
        push!( newset, f( x ) )
    end
    newset
end

function haseq( s::Set{T}, el::T ) where {T}
    for x in s
        if x == el
            return true
        end
    end
    return false
end


function subsets( s::Set{T}, n::Int ) where {T}
    card = length( s )
    if card < n 
        return Set()
    end        
    if n == 0 
        return Set(Any[Set{T}()])
    end
    if card == 0
        return Set()
    end
    el = first( s )
    #@hiya el
    rest = copy( s )
    delete!( rest, el )
    #@hiya rest
    ans = sizehint( subsets( rest, n ), binomial( card, n ) )
    #ans = subsets( rest, n )
    smalls = subsets( rest, n - 1 )
    #@hiya ans, smalls
    for ss in smalls
        push!( ss, el )
    end
    union!( ans, smalls )
    ans
end


function subsets( n )
    if n==0
        Set{Set{Int}}([Set{Int}()])
    else
        smaller = subsets( n - 1 )
        bigger = mapset( s->union(s,Set(n)), smaller )
        union( smaller, bigger )
    end
end

function second( lst )
    lst[2]
end
function third( lst )
    lst[3]
end

function rand_element(list::Array{T,1},probs::Array{N,1})::T where {T,N <: Number}
    len = length(list)
    if len > length(probs)
        println("probs are too short")
        return(rand(list))
    end
    tot = sum(probs[1:len])
    if tot == 0
        return rand(list)
    end
    cnt = 0
    ind = 0
    r = rand() * tot
    while cnt <= r
        ind += 1
        if (ind>len) println(probs,tot,r,cnt) end
        cnt += probs[ind]
    end
    return list[ind]
end

## approx. Dirichlet noise in dim n with alpha = 2.6?
function dir_noise1(n::Int64)::Array{Float64,1}
    noise = rand(n)
    tot = sum(noise)
    return noise / tot
end

## alpha = 2?
function dir_noise2(n::Int64)::Array{Float64,1}
    brks = sort(rand(n-1))
    b1,b2 = vcat([0],brks),vcat(brks,[1])
    return b2 - b1
end

randBool(p::Float64=0.5) = rand()<p

## read dict of Expr,val pairs into current scope
macro read_dict(d)
    quote
        for (a,v) in $(esc(d))
            ex = Expr(:(=),a,v)
            eval(ex)
        end
    end
end

macro print_dict(d)
    quote
        println($(esc(d)))
        for (a,v) in $(esc(d))
            println(a," => ",v)
        end
    end
end


function eval_with_dict(dict,ex::Expr)
    @read_dict(dict)
    return eval(ex)
end

# function histogram(data::Array{Dict,1},ex::Expr,rng)
#     header = string("cut "," # "," Prob")

function df_get_entry(f,v,c)
    return f[f[:variable].==v,c][1]
end

function find_min_local(fun::Function,left::Float64,right::Float64,d::Float64)::Float64
    l,r = left,right
    while r - l > d
        ml = 2/3 * l + 1/3 * r
        mr = 1/3 * l + 2/3 * r
        if fun(ml) > fun(mr)
            l,r = ml,r
        else
            l,r = l,mr
        end
    end
    return 0.5 * (l + r)
end

function fixed_point(f::Function,init,dist::Function,tol::Float64)
    old = new = init
    d = Inf
    while d > tol
        new = f(old)
        d = dist(old,new)
        old = new
    end
    return new
end

function interpolation_lin(pairlist::Array{Array{Float64,1},1})::Function
    function(t::Float64)
        if t <= pairlist[1][1] || t >= pairlist[end][1]
            0
        else
            ind = findfirst(p -> p[1]>=t, pairlist) - 1
            #println(ind)
            a,b = pairlist[ind],pairlist[ind+1]
            a[2] + (b[2]-a[2])*(t-a[1])/(b[1]-a[1])
        end
    end
end

function interpolation_lin_int(pairlist::Array{Array{Float64,1},1})#::Function
    intlist = Array{Float64,1}[]
    cum = 0
    oldi = 0
    oldp = pairlist[1]
    for i = 1:length(pairlist)
        #println(pairlist[i])
        newp = pairlist[i]
        add = 0.5*(newp[2] + oldp[2])*(newp[1] - oldp[1])
        newi = oldi + add
        oldp,oldi = newp,newi
        #println(add)
        push!(intlist,[newp[1],newi])
    end
    intlin = interpolation_lin(pairlist)
    function(t::Float64)
        if t <= pairlist[1][1]
            0
        elseif t >= pairlist[end][1]
            intlist[end][2]
        else
            ind = findfirst(p -> p[1]>=t, pairlist) - 1
            #println(ind)
            p1,b = pairlist[ind],intlin(t)
            intlist[ind][2] + 0.5*(p1[2]+b)*(t-p1[1])
        end
    end
end

function array_transpose(tupes::Array{Array{T,1},1})::Array{Array{T,1},1} where T
    r,c = length(tupes),length(tupes[1])
    [ [ tupes[i][j] for i=1:r ] for j=1:c ]
end

function invert_dict(d)
    Dict(values(d) .=> keys(d))
end

;
