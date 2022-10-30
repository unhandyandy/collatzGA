using Dates,ImageView,Test

struct cn2
    bits::Vector{Bool}
    fitness::Union{Missing,Float16}
    len::UInt16
end

function cn2(n::Int)::cn2
    m = n
    while m & 1 == 0
        m = m >> 1
    end
    bits = Bool[]
    len = 0
    while m>0
        len += 1
        push!(bits,m&1)
        m = m >> 1
    end
    cn2(reverse(bits),missing,len)
end
            
# function copy(cn::cn2)::cn2
#     cn2(cn.bits,cn.fitness,cn.len)
# end

function nextGen(cn::cn2)::cn2
    bits = cn.bits
    len = cn.len + 1
    newbits = zeros(Bool, len)
    prev = true
    carry = true
    for i in 1:len - 2
        current = bits[end - i]
        odd = current ⊻ prev ⊻ carry
        carry = (current && (prev || carry)) || (prev && carry)
        newbits[end - i + 1] = odd
        prev = current
    end
    odd = prev ⊻ carry
    carry = prev && carry
    newbits[2] = odd
    newbits[1] = carry
    if !newbits[1]
        newbits = newbits[2:end]
        len -= 1
    end
    while !newbits[end]
        newbits = newbits[1:end-1]
        len -= 1
    end
    cn2(newbits,missing,len)
end

function findOrder(cn::cn2)::Int64
    nn = cn
    o = 0
    while nn.len > 1
        o += 1
        nn = nextGen(nn)
    end
    o
end

function setFitness(cn::cn2)::cn2
    o = findOrder(cn)
    cn2(cn.bits,o,cn.len)
end

function randomCN2(len::Int)::cn2
    bits = rand(Bool,len)
    bits[1] = true
    bits[end] = true
    cn = cn2(bits,missing,len)
    setFitness(cn)
end

function mutate1(cn::cn2)::cn2
    bits = deepcopy(cn.bits)
    idx = rand(2:cn.len-1)
    bits[idx] = !bits[idx]
    nn = cn2(bits,missing,cn.len)
    setFitness(nn)
end

function choice(w1::Number,w2::Number)::Bool
    rand() < w1/(w1+w2)
end

function MH1(start::cn2,t::TimePeriod)::Tuple{cn2,cn2}
    nn = deepcopy(start)
    best = deepcopy(start)
    endtime = now() + t
    while now() < endtime
        nnn = mutate1(nn)
        if nnn.fitness > best.fitness
            best = nnn
            println("Best: ",best.fitness)
        end
        if choice(nnn.fitness,nn.fitness)
            nn = nnn
        end
    end
    best, nn
end

function evolve(cn::cn2)::Vector{cn2}
    evo = cn2[cn]
    nn = cn
    while nn.len>1
        nn = nextGen(nn)
        push!(evo,nn)
    end
    evo
end

function grevo(evo::Vector{cn2})::Matrix{Bool}
    width = max(map(c->c.len,evo)...)
    height = length(evo)
    mat = zeros(height,width)
    for r = 1:height
        mat[r,1:evo[r].len] = evo[r].bits
    end
    mat
end




###  Tests

@test findOrder(cn2(37))==6
@test findOrder(cn2(63))==39
