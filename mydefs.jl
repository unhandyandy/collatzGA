#-*- julia; eval: (auto-fill-mode) -*-

macro hiya(ex)
    quote
        begin
            local val = $ex
            print( $(string(ex)), ": ", $ex, "\n" )
            val
        end
    end
end

function nestList{T}( f::Function, n::Int, x0::T ) 
    t = typeof( x0 )
    res = t[]
    cur = x0
    for i=0:n
        push!( res, cur )
        cur = f( cur )
    end
    return res
end

function evenQ{T<:Number}( n::T ) 
    n % 2 == 0
end

function mapAt( f::Function, x, lev::Int )
    if lev == 0
        f(x)
    else
        map( y -> mapAt( f, y, lev - 1 ), x )
    end
end

function mapMat{T<:Any}( f::Function, arr::Array{T,2} )
    r, c = size( arr )
    [ f( arr[i,j] ) for i=1:r, j=1:c ]
end

function iterate{T}( f::Function, n::Int, x::T )
    n==0 ? x : iterate( f, n - 1, f( x ) )
end
