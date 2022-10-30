
struct position
    pos::Matrix{Int8}
    mov::Int8
    val::Int8
    n::Int8
end

position0 = position(zeros(3,3),1,0,1)

function getChildren(p::position)::Array{position,1}
    mat = p.pos
    locs =  findall(x -> x==0,p.pos)
    n = length(locs)
    children = Array{position}(undef,n)
    newmov = p.n==9 ? 0 : 3 - p.mov
    newn = p.n + 1
    for i in 1:n
        tmp = copy(mat)
        tmp[locs[i]] = p.mov
        tmp0 = position(tmp,newmov,0,newn)
        children[i] = position(tmp,newmov,staticVal(tmp0),newn)
    end
    children
end

function getLines(p::position)::Matrix{Int8}
    mat = p.pos
    lines = zeros(Int8,8,3)
    lines[1,:] = mat[1,:]
    lines[2,:] = mat[2,:]
    lines[3,:] = mat[3,:]
    lines[4,:] = mat[:,1]
    lines[5,:] = mat[:,2]
    lines[6,:] = mat[:,3]
    lines[7,:] = [mat[1,1] mat[2,2] mat[3,3]]
    lines[8,:] = [mat[1,3] mat[2,2] mat[3,1]]
    lines
end

function staticVal(p::position)::Int8
    lines = getLines(p)
    if any(all(x -> x==1,lines,dims=2))
        val = 1
    elseif any(all(x -> x==2,lines,dims=2))
        val = -1
    else
        val = 0
    end
    val
end

function setStaticVal(p::position)::position
    position(p.pos,p.mov,staticVal(p),p.n)
end

function findBest(p::position)::Tuple{Int8,position}
    pc = getChildren(p)
    sgn = -2*p.mov + 3
    if p.n==9
        val = pc[1].val
        best = pc[1]
    else
        vs = map(x -> x.val,pc)
        m,idx = findmax(sgn*vs)
        if m==1
            val = sgn
            best = pc[idx]
        else
            bv = -2
            numc = 10 - p.n
            for i in 1:numc
                c = pc[i]
                v,_ = findBest(c)
                if sgn*v > bv
                    bv = sgn*v
                    val = v
                    best = c
                end
                if v==sgn
                    break
                end
            end
        end
    end
    return val, best
end
