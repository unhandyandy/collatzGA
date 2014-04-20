# make parametric
type oddFactor{T<:Number}
    odd::T
    power::Int
end

# function oddFactor{T<:Number}( odd::T, pow::Int )
#     oddFactor{T}( odd, pow )
# end


# oddFactor( m::Int, p::Int ) = begin
#     conQ = floor(log(2,m)) > floor(log(2,m-1))
#     oddFactor( m, p, conQ )
# end

function (*){T<:Number}( m::oddFactor{T}, n::oddFactor{T} ) 
    oddFactor{T}( m.odd*n.odd, m.power + n.power )
end


# oddFactor{T<:Number}(n::T) = 
# if n%2 == 1
#     oddFactor{T}(n,0)
# else
#     oddFactor{T}(convert(T,1),1) * oddFactor{T}( div(n,2) )
# end

function oddFactor{T<:Number}(n::T) 
    if n%2 == 1
        oddFactor{T}(n,0)
    else
        #fe = int(factor(n)[T(2)])
        fo, fe = n, 0
        while evenQ( fo )
            fo = div( fo, 2 )
            fe += 1
        end
        oddFactor{T}( fo, fe )
    end
end


