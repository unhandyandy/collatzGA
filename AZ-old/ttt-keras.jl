using Keras

import Keras.Layers: Dense, Activation, Concatenate, BatchNormalization

# def policyIterSP(game):
#     nnet = initNNet()                                       # initialise random neural network
#     examples = []    
#     for i in range(numIters):
#         for e in range(numEps):
#             examples += executeEpisode(game, nnet)          # collect examples from this game
#         new_nnet = trainNNet(examples)                  
#         frac_win = pit(new_nnet, nnet)                      # compare new net with previous net
#         if frac_win > threshold: 
#             nnet = new_nnet                                 # replace with new net            

function policyIter(numiters::Int,numeps::Int,net::Model)::Model
    newnet = net
    examples = []
    for i = 1:numiters
        for j = 1:numeps
            examples = vcat(examples, execEp(newnet::Model))
        end
        return  trainNet(examples)
    end
end

function policyIter(numiters::Int,numeps::Int)::Model
    newnet = makeNN()
    policyIter(numiters,numeps,newnet)
end


                           
# def executeEpisode(game, nnet):
#     examples = []
#     s = game.startState()
#     mcts = MCTS()                                           # initialise search tree
        
#     while True:
#         for _ in range(numMCTSSims):
#             mcts.search(s, game, nnet)
#         examples.append([s, mcts.pi(s), None])              # rewards can not be determined yet 
#         a = random.choice(len(mcts.pi(s)), p=mcts.pi(s))    # sample action from improved policy
#         s = game.nextState(s,a)
#         if game.gameEnded(s):
#             examples = assignRewards(examples, game.gameReward(s)) 
# return examples

struct TTTexample
    state::TTTState
    policy::Array{Float64,1}
    value::Float64
end

function get_pis(policy::MCTSPlanner,state::TTTState)::Array{Float64,1}
    node = policy.tree[state]
    totnum = node.N
    sanodes = node.sanodes
    res = map( n -> n.N, sanodes)
    return res / totnum
end

function rand_move(state::TTTState,wts::AbstractArray{Float64,1})::TTTMove
    mvs = actions(mdp,state)
    inds = map(m -> action_index(mdp,m),mvs)
    possmvs = map(i -> action_index_inv(mdp,i),inds)
    posswts = map(i -> wts[i],inds)
    return rand_element(possmvs,posswts)
end

function execEp(policy::MCTSPlanner,numeps::Int)::Array{TTTexample,1}
    examples = []
    s = startpos
    while true
        for _ = 1:numeps
            MCTS.simulate(policy,s,9)
        end
        pis = get_pis(policy,s)
        examples = vcat(examples,TTTexample(s,pis,0))
        m = rand_move(s,pis)
        s = next_state(s,m)
        if s.done
            examples = assign_rewards(examples)
            return examples
        end
    end
end

function assign_rewards(examples::AbstractArray{TTTexample,1})::AbstractArray{TTTexample,1}
    len = length(examples) + 1
    termstate = last(examples).state
    endreward = reward(termstate)
    return [ TTTexample(examples[len - i].state,
                        examples[len - i].policy,
                        termreward * (-1)^(len - i -1))
             for i = 1:len ]
end

function get_net_policy(net::Model,state::TTTState)::AbstractArray{Float64,1}
    return net.predict(state.board,state.player)[1]
end

function best_sanode_UCB(snode::StateNode, state::TTTState,c::Float64,net::Model)::StateActionNode
    best_UCB = -Inf
    best_sanode = snode.sanodes[1]
    sN = snode.N
    for sanode in snode.sanodes
        if (sN == 1 && sanode.N == 0) || c == 0.0
            UCB = sanode.Q
        else
            #UCB = sanode.Q + c*sqrt(log(sN)/sanode.N)
            UCB = sanode.Q + c * get_net_policy(net,state) * sqrt(sanode.N)/(sN + 1)
        end
        @assert !isnan(UCB)
        @assert !isequal(UCB, -Inf)
        if UCB > best_UCB
            best_UCB = UCB
            best_sanode = sanode
        end
    end
    return best_sanode
end


function makeNN()::Model
    
