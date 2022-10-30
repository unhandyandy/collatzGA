using Flux

include("../mydefs.jl")


function policyIter!(numiters::Int64,numeps::Int64,numwus::Int64,policy)
    examples = []
    #newpolicy = solve(solver, policy.mdp)
    local newpolicy
    for i = 1:numiters
        println("iter: ",i)
        #Flux.testmode!(predict,true)
        # println("params 1: ",params(newpolicy.mdp.net)[1][1:5,1])
        # println(" end - 1: ",params(newpolicy.mdp.net)[end - 1][1:5])
        #println("     end: ",params(newpolicy.mdp.net)[end])
        println("tp6: ",predict(testpos6))
        println("tp8: ",predict(testpos8))
        println("tp9: ",predict(testpos9))
        println(" st: ",predict(startpos))
        #Flux.testmode!(predict,false)
        for j = 1:numeps
        if mod(j,5)==0 println("ep. #: ",j) end
            newpolicy = solve(solver, policy.mdp)
            newex =  execEp(newpolicy, numwus,ceil(Int64,i/2))
            #if j > numeps / 2
            examples = vcat(examples,newex)
            #end
        end
        println(" ex: ",map(e -> e[2][end],filter(e -> POMDPs.isequal(e[1],testpos8),examples)))
        trainNet(examples,newpolicy)
    end
    return newpolicy
end

function param_dist(p1,p2)
    ps = collect(zip(p1,p2))
    ds = map(p -> Flux.mse(p[1],p[2]),ps)
    return sum(ds)
end

function policyIter!(numiters::Int,numeps::Int,numwus::Int64)
    newnet = makeNN()
    global mdp =  TTTGame(3,newnet)
    global policy = solve(solver,mdp)
    global predict = mapleaves(Flux.Tracker.data, policy.mdp.net)
    Flux.testmode!(predict, true)
    return policyIter!(numiters,numeps,numwus,policy)
end

function prep(state::TTTState)::Tuple{Array{Float64,1},Float64}
    return  vec(state.board),state.player
end

# function make_predictors(policy)
#     predpol = mapleaves(Flux.Tracker.data, predict[1])
#     predval = mapleaves(Flux.Tracker.data, predict[2])
#     Flux.testmode!(predpol, true)
#     Flux.testmode!(predval, true)
#     return s -> predpol(vec(s.board),s.player), s -> predval(vec(s.board),s.player)
# end



# struct TTTexample
#     state::TTTState
#     policy::Array{Float64,1}
#     value::Float64
# end

function get_pis(policy::MCTSPlanner,state::TTTState)::Array{Float64,1}
    wts = get_wts(policy,state)
    tot = sum(wts)
    if tot == 0
        tot = 9
        wts = ones(9)
    end
    wts /= tot
    return wts
end

function get_wts(policy::MCTSPlanner,state::TTTState)::Array{Float64,1}
    node = policy.tree[state]
    sanodes = node.sanodes
    mvs = [ n.action for n in sanodes ]
    function getN(m::TTTMove)::Int64
        if m in mvs
            ind = findfirst(n -> n.action==m,sanodes)
            return sanodes[ind].N
        else
            return 0
        end
    end
    return map( m -> getN(m), action_list)
end

function rand_move(policy::MCTSPlanner,state::TTTState)::TTTMove
    #mvs = actions(mdp,state)
    wts = get_wts(policy,state)
    # if (wts[1]===NaN) println(state) end
    return rand_element(action_list,wts)
end

# function symm_examples(exs::Array{Any,1})::Array{Any,1}
#     res = [ (apply_sym(e[1],s), for s in symms for m in

function execEp(policy::MCTSPlanner,numwus::Int,depth::Int64=100,start::TTTState=startpos)::Array{Any,1}
    examples = []
    s = start
    while true
        #println(numwus)
        for w = 1:numwus
            MCTS.simulate(policy,s,9)
        end
        pis = get_pis(policy,s)
        # println(pis)
        # if (sum(pis)==0) println(s) end
        examples = vcat(examples,(s,pis,0))
        m = rand_move(policy,s)
        s,r = next_state(s,m)
        if s.done
            examples = assign_rewards(examples,r)
            d = min(depth,length(examples)) - 1
            examples = examples[end - d:end]
            examples = [ apply_sym(s,e) for e in examples for s in symms ]
            #println(examples,"\n\n")
            return examples
        end
    end
end

function assign_rewards(examples::AbstractArray{Any,1},r::Int64)::AbstractArray{Any,1}
    len = length(examples) + 1
    #termplyr = last(examples).state.player
    endreward = r > 0 ? 1.0 : 0.0
    return [ (examples[i][1],
              examples[i][2],
              endreward * (-1)^(len - i - 1))
             for i = 1:(len - 1) ]
end

function get_net_policy(net,state::TTTState)
    return net(state)[1].tracker.data
end

# function get_net_best_move(net,state::TTTState)
#     pis = get_net_policy(net[1:2],state)
#     bv = max(pis...)
#     bi = findfirst(bv,pis)
#     return actions(mdp,state)[bi]
# end

function MCTS.best_sanode_UCB(snode,c::Float64,state::TTTState,policy)
    # ::MCTS.StateActionNode{TTTMove}
    ## state field is not implemented!!
    #state = snode.state
    best_UCB = -Inf
    best_sanode = snode.sanodes[1]
    sN = snode.N
    ind = 0
    pis = 0.8 * predict(state)[1] + 0.2 * dir_noise2(9)
    for sanode in snode.sanodes
        # if (sN == 1 && sanode.N == 0) || c == 0.0
        #     UCB = sanode.Q
        # else
            #println(c)
            #UCB = sanode.Q + c*sqrt(log(sN)/sanode.N)
            ## predict calls policy.mdp.net
            #UCB = sanode.Q + c * policy.mdp.net(state)[1][ind] * sqrt(sanode.N)/(sN + 1)
            ind = POMDPs.action_index(policy.mdp,sanode.action)
            UCB = sanode.Q + c * pis[ind] * sqrt(sanode.N)/(sN + 1)
            # if POMDPs.isequal(state,TTTState([1 0 0; 0 2 0; 0 0 0], 1, false))
            #println(UCB)
            # end
        #end
        @assert !isnan(UCB)
        @assert !isequal(UCB, -Inf)
        #if UCB == best_UCB println("nonrandom: ",state) end
        if UCB > best_UCB #|| (UCB == best_UCB && rand() < 0.5)
            best_UCB = UCB
            best_sanode = sanode
        end
    end
    return best_sanode
end

function denseL(innum::Int64,outnm::Int64)::Chain
    Chain(Dense(innum,outnm,relu,initb=Flux.glorot_normal), BatchNorm(outnm))
    #Chain(Dense(innum,outnm), BatchNorm(outnm))
    #Chain(Dense(innum,outnm), BatchNorm(outnm,relu))
    #Chain(Dense(innum,outnm,relu,initb=Flux.glorot_normal))
end

# opt = ADAM(params(scanner, encoder))
# evalcb = () -> @show testloss()

# Flux.train!(loss, train, opt, cb = throttle(evalcb, 10))

# struct Trunk
#     mixerfunc
#     remlayers
# end

# Flux.treelike(Trunk)

# (t::Trunk)(x, y) = t.remlayers(t.mixerfunc(x, y))

#mixer(m::TTTmodel,x,y) = vcat(x,y)

struct TTTmodel
    in_bd
    in_pl
    mixer
    mtrunk
    outpol
    outval
end

Flux.treelike(TTTmodel)

function (m::TTTmodel)(s::TTTState)
    return (m.outpol(m.mtrunk(m.mixer(m.in_bd(vec(s.board)),m.in_pl(s.player)))),
            m.outval(m.mtrunk(m.mixer(m.in_bd(vec(s.board)),m.in_pl(s.player)))))
end

function makeNN()::TTTmodel
    return TTTmodel(Dense(9,40, relu,initb=Flux.glorot_normal),
                    Dense(1,10, relu,initb=Flux.glorot_normal),
                    vcat,
                    Chain(denseL(50,80),denseL(80,160),denseL(160,80)),
                    Chain(denseL(80,40),denseL(40,9),softmax),
                    Chain(denseL(80,40),Dense(40,1,initb=Flux.glorot_normal),x -> tanh.(x)))
end
                    

# function makeNN()
#     in_bd = Dense(9,40, relu)
#     in_pl = Dense(1,10, relu)
#     mixer(x,y) = vcat(in_bd(x),in_pl(y))
#     #trunk = mixer(x,y) |> denseL(50,80) |> denseL(80,80) |> denseL(80,40)
#     maintr = Trunk(mixer,Chain(denseL(50,80),denseL(80,80),denseL(80,40)))
#     #outmv(x,y) = trunk(x,y) |> denseL(40,20) |> denseL(20,9) |> softmax
#     outmv = Trunk(maintr,Chain(denseL(40,20),denseL(20,9),softmax))
#     #outvl(x,y) = trunk(x,y) |> denseL(40,20) |> denseL(20,1) |> σ
#     outvl = Trunk(maintr,Chain(denseL(40,20),Dense(20,1),x -> tanh.(x)))
#     #println(params(in_bd,in_pl,maintr,outmv,outvl))
#     return (outmv,outvl,params(in_bd,in_pl,maintr,outmv,outvl))
# end

l2nm(prms) =  sqrt(sum(map(p -> Flux.mse(p,0),prms)))

function trainNet(exs::AbstractArray{Any,1},policy)
    prms = params(policy.mdp.net)
    old_params = [ Flux.data.(p) for p in prms ]
    #outs = policy.mdp.net[1:2]
    #net = call_net(outs...)
    net = policy.mdp.net
    ### ADAM SGD
    opt = ADAM(prms)
    function loss(state::TTTState,actpol::Array{Float64,1},actval::Float64)
        netpol,netval = net(state)
        #actpol,actval = act
        ## include norm of parameters?
        # println()
        # println(state)
        # println(netpol)
        # println(actpol)
        # println(Flux.crossentropy(netpol,actpol))
        l = 0.9 * Flux.crossentropy(netpol,actpol) + 0.1 * Flux.mse(netval,actval) + 0.01 * l2nm(prms)
        # println("netpol = ",netpol)
        # println("actpol = ",actpol)
        # println("l = ",l)
        return l
    end
    #println(exs)
    #data = map( e -> (e.state,(e.policy,e.value)),exs)
    Flux.@epochs 10 Flux.train!(loss,exs,opt)
    new_params = params(policy.mdp.net)
    new_params = [ Flux.data.(p) for p in new_params ]
    println("prms d: ",param_dist(old_params,new_params))
    #println("old prms: ",old_params)
    return 0
end

    
# function call_net(polnet,valnet)::Function
#     function fun(state::TTTState)
#         inbd,inpl = vec(state.board),state.player
#         return polnet(inbd,inpl),valnet(inbd,inpl)
#     end
# end

#import MCTS.simulate

function MCTS.simulate(policy::AbstractMCTSPlanner, state, depth::Int64)
    # model parameters
    mdp = policy.mdp
    discount_factor = discount(mdp)
    rng = policy.rng

    # once depth is zero return
    if depth == 0 || isterminal(policy.mdp, state)
        return 0.0
    end

    # if unexplored state add to the tree and run rollout
    if !MCTS.hasnode(policy, state)
        newnode = MCTS.insert_node!(policy, state)
        return estimate_value(policy.solved_estimate, policy.mdp, state, depth)
    end
    # if previously visited node
    snode = MCTS.getnode(policy, state)

    # pick action using UCT
    snode.N += 1 # increase number of node visits by one
    sanode = MCTS.best_sanode_UCB(snode, policy.solver.exploration_constant,state,policy)

    # transition to a new state
    sp, r = generate_sr(mdp, state, sanode.action, rng)
    
    if policy.solver.enable_tree_vis
        MCTS.record_visit(policy, sanode, sp)
    end

    q = r + discount_factor * MCTS.simulate(policy, sp, depth - 1)
    sanode.N += 1
    sanode.Q += ((q - sanode.Q) / (sanode.N)) # moving average of Q value
    return q
end

function play_against_net(policy,sp)
    #outs = policy.mdp.net[1:2]
    pos = copy(startpos)
    #rep = pos
    p = sp
    inv = sym0
    while !pos.done
        if p==1
            println("Your move: ")
            mstr = readline()
            m = TTTMove(parse(Int,mstr[1]),parse(Int,mstr[3]))
            pos,+ = next_state(pos,m)
        else
            m = get_net_best(policy,pos)
            pos,_ = next_state(pos,m)
            #pos = apply_sym(inv,nst)
        end
        println(m)
        #rep,inv = choose_rep(pos)
        #println(pos)
        p = 3 - p
        printpos(pos)
    end
end

function get_net_best(p,s::TTTState)
    pol = p.mdp.net(s)[1]
    mvs = actions(p.mdp,s)
    z = collect(zip(pol[1:length(mvs)],mvs))
    sort!(z,by=first)
    return z[end][2]
end


#### init

mdp = TTTGame(3,makeNN())
action_list = actions(mdp)
policy = solve(solver, mdp)

predict = mapleaves(Flux.Tracker.data, mdp.net)
Flux.testmode!(predict, true)

;

#### debugging

# TTTexample[TTTexample(TTTState([0 0 0; 0 0 0; 0 0 0], 1, false), [0.999882, 1.47048e-5, 1.47048e-5, 1.47048e-5, 1.47048e-5, 1.47048e-5, 1.47048e-5, 1.47048e-5, 1.47048e-5], 1.0), TTTexample(TTTState([1 0 0; 0 0 0; 0 0 0], 2, false), [0.000132607, 0.000169443, 0.000221012, 6.63037e-5, 5.89366e-5, 0.999322, 1.47342e-5, 1.47342e-5, 0.0], -1.0), TTTexample(TTTState([1 0 2; 0 0 0; 0 0 0], 1, false), [0.999911, 1.48085e-5, 1.48085e-5, 1.48085e-5, 1.48085e-5, 1.48085e-5, 1.48085e-5, 0.0, 0.0], 1.0), TTTexample(TTTState([1 0 2; 1 0 0; 0 0 0], 2, false), [0.961997, 1.48216e-5, 1.48216e-5, 1.48216e-5, 1.48216e-5, 0.0379434, 0.0, 0.0, 0.0], -1.0), TTTexample(TTTState([1 0 2; 1 0 0; 2 0 0], 1, false), [6.98553e-6, 0.999972, 6.98553e-6, 6.98553e-6, 6.98553e-6, 0.0, 0.0, 0.0, 0.0], 1.0), TTTexample(TTTState([1 0 2; 1 1 0; 2 0 0], 2, false), [4.89929e-6, 0.999985, 4.89929e-6, 4.89929e-6, 0.0, 0.0, 0.0, 0.0, 0.0], -1.0), TTTexample(TTTState([1 0 2; 1 1 0; 2 2 0], 1, false), [3.77306e-6, 0.999996, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 1.0)]

#testBN = Chain(Dense(4,8,relu), BatchNorm(8))
