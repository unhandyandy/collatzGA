
using LinearAlgebra
using Debugger,Dates,JLD2

maxDepth = 10

struct Position
    board::Matrix{Int8}
    mover::Int8
    staticVal::Float16
    id::String
    status::Int8
    term::Bool
end

mutable struct Node
    position::Position
    n::Int8
    hist::Array{Position,1}
    dynVal::Union{Float16,Missing}
    children::Union{Array{Node,1},Missing}
    depth::Union{Int8,Missing}
    ## desired depth
    dd::Union{Int8,Missing}
    parent::Union{Node,Missing}
    best::Union{String,Missing}
end

function isequal(p1::Position,p2::Position)::Bool
    p1.id==p2.id 
end

function findHeight(bd::Matrix{Int8},c::Integer)::Int8
    count(x->x==0,bd[1:6,c])
end

coefficients = ones(Float16,138)./16

function scoreLine(l::Vector{Int8},p::Integer,q::Integer)::Tuple{Float16,Int8}
    stat::Int8 = 0
    if all(x->x==q,l)
        stat = -1
    end
    if all(x->x!=p,l)
        return (count(x->x==q,l)/4,stat)
    elseif all(x->x!=q,l)
        return (count(x->x==p,l)/4,stat)
    else
        return (0,0)
    end
end

function staticValue(bd::Matrix{Int8},player::Integer)::Tuple{Float16,Int8,
                                                           Vector{Float16}}
    p,q = player,toggleMover(player)
    vec = zeros(Float16,138)
    numzs = count(x->x==0,bd)
    idx = mod(numzs,2)==0 ? 0 : 69
    stat::Int8 = 0
    ## check verticals
    for j = 1:7
        for i = 1:3
            line = bd[i:i+3,j]
            idx += 1
            vec[idx],st = scoreLine(line,p,q)
            stat = min(stat,st)
        end
    end
    ## check horizontals
    for i = 1:6
        for j = 1:4
            line = bd[i,j:j+3]
            idx += 1
            vec[idx],st = scoreLine(line,p,q)
            stat = min(stat,st)
        end
    end
    ## check down diags
    for j = 1:4
        for i = 1:3
            line = diag(bd[i:i+3,j:j+3])
            idx += 1
            vec[idx],st = scoreLine(line,p,q)
            stat = min(stat,st)
        end
    end
    ## check up diags
    for j = 1:4
        for i = 4:6
            line = diag(reverse(bd[i-3:i,j:j+3],dims=1))
            idx += 1
            vec[idx],st = scoreLine(line,p,q)
            stat = min(stat,st)
         end
    end
    global coefficients
    return (vec⋅coefficients,stat,vec)
end

function toggleMover(m::Integer)::Int8
    Int8(3) - m
end

function Position(mat::Matrix{Int8},m::Integer)::Position
    val,stat = staticValue(mat,m)
    term = stat!=0 || count(x->x==0,mat[1,:])==0
    id = string(m)
    for x in mat
        id *= string(x)
    end
    Position(mat,m,-val,id,stat,term)
end

position0 = Position(zeros(Int8,6,7),1)

function getChildren(p::Position)::Array{Position,1}
    mat = p.board
    newmover = toggleMover(p.mover)
    n = count(x->x==0,mat[1,1:7])
    children = Array{Position}(undef,n)
    cnt = 1
    for j=1:7
        if mat[1,j]==0
            tmp = copy(mat)
            i = findHeight(mat,j)
            tmp[i,j] = p.mover
            children[cnt] = Position(tmp,newmover)
            cnt += 1
        end
    end
    children
end

function display(p::Position)
    b = p.board
    dmat = Array{String,2}(undef,7,7)
    for j=1:7
        dmat[1,j] = string(j)
    end
    for i=1:6
        for j=1:7
            dmat[i+1,j] = b[i,j]==0 ? " " : b[i,j]==1 ? "X" : "O"
        end
    end
    pmat = Array{String,1}(undef,7)
    for i=1:7
        pmat[i] = join(dmat[i,1:7]," ")
    end
    Base.display(pmat)
    println(p.mover," to move")
    println("status: ",p.status)
end

###### memory dict

mutable struct memCell
    depth::Int8
    val::Union{Missing,Float16}
    best::Union{Missing,String}
end

memCell0 = memCell(Int8(0),missing,missing)

memory = Dict{String,memCell}()

#### tree

function Node(p::Position,n::Integer,hist::Array{Position,1})::Node
    global memory
    dynval = p.staticVal
    mem = get(memory,p.id,memCell0)
    if mem.depth==0 && ismissing(mem.val)
        if p.term
            val = p.status
        else
            val = p.staticVal
        end
    else
        val = mem.val
    end
    ########      dynVal,childr,depth,dd,    best
    Node(p,n,hist,val,missing,Int8(0),mem.depth,missing,mem.best)
end

function copyNodeToMem(node::Node)
    global memory
    #println("old cell: ",get(memory,node.position.id,memCell0).depth)
    mem = memCell(node.depth,node.dynVal,node.best)
    #println("new cell: ",mem.depth)
    memory[node.position.id] = mem
end

function makeStartNode()::Node
    Node(position0,Int8(1),Position[])
end

node0 = makeStartNode()

function display(n::Node)
    display(n.position)
    # println(n.position.status)
    println("dynVal: ",n.dynVal)
    println(n.depth)
    # println(n.dd)
    println(n.position.id)
    # par = n.parent
    # s = ismissing(par) ? par : par.position.id
    # println("parent: ",s)
    # if ismissing(n.children)
    #     println("children: missing")
    # else
    #     for c in n.children
    #         println(c.dynVal)
    #     end
    # end
end

function bearChildren!(node::Node)
    n::Int8 = node.n + 1
    hist = [node.hist; node.position]
    childPoss = getChildren(node.position)
    children = map(p->Node(p,n,hist),childPoss)
    map(n->n.parent=node,children)
    node.children = children
end

function valAtDepth(node::Node,dep::Integer)::Tuple{Union{Missing,Float16},
                                                 Union{Missing,String}}
    ###     using memory                       
    if (node.depth >= dep) || node.position.term
        ### using search           terminal node
        return (node.dynVal,node.best)
    else
        return (missing,missing)
    end
end
    
function deepen1!(node::Node,dep::Integer,alpha::Number,beta::Number)::Float16
    val,_ = valAtDepth(node,dep)
    if !ismissing(val)
        return val
    end
    if ismissing(node.children)
        bearChildren!(node)
    end
    bs = missing
    if size(node.children)[1]==0
        bv = node.position.status
    else
        b,a = -beta,-alpha
        sort!(node.children,by=x->x.dynVal)
        bs = node.children[1].position.id
        for c in node.children
            v = deepen1!(c,Int8(dep - 1),b,a)
            if v < a
                a = v
                bs = c.position.id
            end
            if a<=b
                a = -1
                bs = c.position.id
                break
            end
        end
    end
    node.dynVal = -a
    if alpha==-1 && beta==1
        node.depth = dep
        # node.dd = node.depth
        node.best = bs
    end
    return node.dynVal
end

## diagnostic
function countNodes(node)
    tot = 1
    if !ismissing(node.children)
        for c in node.children
            tot += countNodes(c)
        end
    end
    return tot
end
#### 

function findBest(node::Node)::Node
    idx = findfirst(c->c.dynVal==-node.dynVal,node.children)
    node.children[idx]
end

function strToChild(node::Node,str::String)::Node
    if ismissing(node.children)
        bearChildren!(node)
    end
    idx = findfirst(x->x.position.id==str,node.children)
    newnode = node.children[idx]
    newnode.parent = node
    return newnode
end

## Miles on StackOverflow
function randItem(items, weights)
    cs = cumsum(weights)
    cs /= cs[end]
    return items[findfirst(cs .> rand())]
end

function chooseMove(node::Node,dep::Integer,rnd::Bool=false)::Node
    v,b = valAtDepth(node,dep)
    if !ismissing(b)
        # if !ismissing(node.parent)
        #     node.parent.dynVal = -v
        #     copyNodeToMem(node.parent)
        # end
        # copyNodeToMem(node)
        return strToChild(node,b)
    end
    depth = ismissing(node.dd) ? dep : node.dd>dep ? node.dd : dep
    # println(node.depth+1)
    # println(depth)
    for d = node.depth+1:depth
        deepen1!(node,d,-1,1)
    end
    #@bp
    ## parent&memory update
    # if !ismissing(node.parent)
    #     node.parent.dynVal = -node.dynVal
    #     copyNodeToMem(node.parent)
    # end
    # copyNodeToMem(node)
    if rnd
        chids = map(x->x.position.id,node.children)
        wts = map(x->exp(x).dynVal,node.children)
        # println("wts: ",wts)
        chc = randItem(chids,wts)
    else
        chc = node.best
    end
    newnode = strToChild(node,chc)
    return newnode
end

function selfplay1(d::Integer)::Node
    n = makeStartNode()
    while n.position.status>-1 && n.n<43
        n = chooseMove(n,d,true)
        # display(n)
    end
    # display(strToChild(n,n.best))
    return n
end

## diagnostic
function strToPosition(str::String)::Position
    mat = zeros(Int8,6,7)
    idx = 2
    for c = 1:7
        for r = 1:6
            s = str[idx:idx]
            mat[r,c] = s=="1" ? Int8(1) : s=="2" ? Int8(2) : Int8(0)
            idx += 1
        end
    end
    mover = str[1:1]=="1" ? Int8(1) : Int8(2)
    Position(mat,mover)
end

### diagnostic
function backTraverse(n::Node)
    while !ismissing(n.parent)
        println(n.position.status)
        n = n.parent
    end
end
    
function postMortem(en::Node,comp)
    local node = en.position.mover==comp ? en : en.parent
    # diff::Float16 = -Inf
    # worst = node
    remdep = Inf
    sign = -1
    while !ismissing(node.parent) 
        curdiff = node.parent.dynVal + node.dynVal
        #println("diff: ",curdiff*sign)
        # if curdiff >= diff
        #     diff = curdiff
        #     worst = node.parent
        # end
        #println("node depth: ",node.depth)
        #println("   node dd: ",node.dd)
        if remdep <= node.depth
            node.depth = 0
        elseif curdiff*sign > 0 
            # node.parent.dynVal = -node.dynVal
            remdep = 0
        end
        #println("remdep: ",remdep)
        copyNodeToMem(node)
        node = node.parent
        sign *= -1
        remdep += 1
    end
    if remdep <= node.depth
        node.depth = 0
    end
    copyNodeToMem(node)
    # if dep<maxDepth
    #     dep = max(worst.dd,worst.depth)
    #     worst.dd = dep + 2
    #     worst.best = missing
    #     worst.dynVal = missing
    #     copyNodeToMem(worst)
    # end
end

function parent(node::Node,ind::Integer)::Node
    ind==0 ? node : parent(node,ind - 1).parent
end

function selfImprove(t::TimePeriod,dep)
    endtime = now() + t
    while now() < endtime
        len = length(memory)
        en = selfplay1(dep)
        display(en)
        if length(memory) == len
            postMortem(en,1)
            postMortem(en,2)
        end
        updateCoefficients()
    end
end

function updateCoefficients()
    if length(memory) < 300
        return
    end
    global coefficients
    idx = 0
    tab = zeros(Float16,length(memory),138)
    vals = zeros(Float16,length(memory))
    for (s,m) in memory
        if ismissing(m.val)
            return
        end
        p = strToPosition(s)
        _,_,v = staticValue(p.board,p.mover)
        idx += 1
        tab[idx,:] = v
        vals[idx] = m.val
    end
    cov = tab'*tab
    # covinv = inv(cov)
    # rnk = rank(cov)
    # println("rank: ",rnk)
    println("Recalculating coefficients")
    coefficients = solveLinEq(cov,tab'*vals)
end

function solveLinEq(mat,vec)
    nmat = copy(UpperTriangular(mat))
    pmat = copy(LowerTriangular(mat))
    for i=1:138
       pmat[i,i]=0
    end
    c = zeros(138,1)
    for i=1:100
        c = nmat\(vec-pmat*c)
    end
    return c
end

function makeMove(p::Position,j::Integer)
    if p.board[1,j]==0
        tmp = copy(p.board)
        i = findHeight(tmp,j)
        tmp[i,j] = p.mover
        return Position(tmp,toggleMover(p.mover))
    else
        return missing
    end
end
function playGame(human::Integer)
    playGame(human,makeStartNode())
end
function playGame(human::Integer,en::Node,dep)::Node
    #@bp
    if en.position.status == -1
        winner = en.position.mover==human ? "I" : "You"
        println(winner," win!")
        # postMortem(n)
        # updateCoefficients()
        return en
    elseif human==en.position.mover
        println("Your move: ")
        move = readline()
        j = parse(Int8,move)
        newpos = makeMove(en.position,j)
        if ismissing(newpos)
            println("Illegal move.")
            return playGame(human,en,dep)
        else
            newen = Node(newpos,en.n+1,en.hist)
        end
        println(en.position.id)
    else
        newen = chooseMove(en,Int8(dep))
    end
    newen.parent = en
    # display(newen)
    playGame(human,newen,dep)
end        

function displayChildren(node)
    map(c->c.dynVal,node.children)
end

function Node(p::Position)
    Node(p,1,Position[])
end

## utility
function saveParams()
    @JLD2.save "memory.jld2" memory
    @JLD2.save "coefficients.jld2" coefficients
end

function loadParams()
    global memory, coefficients
    @JLD2.load "memory.jld2"
    @JLD2.load "coefficients.jld2"
end
