

using POMDPs,POMDPToolbox,GraphPlot,MCTS,GraphPlot,D3Trees,Distributions,Winston


# struct ChainMDP <: MDP{Int, Symbol}
#     len::Int
#     p_success::Float64
#     discount::Float64
# end

# # the transition function returns the distribution of sp given that
# # action a is taken in state s
# function transition(mdp::ChainMDP, s::Int, a::Symbol)
#     if a == :right
#         success = min(s+1, mdp.len)
#         failure = max(s-1, 1)
#     else # a == :left
#         success = max(s-1, 1)
#         failure = min(s+1, mdp.len)
#     end
#     return SparseCat([success, failure], [mdp.p_success, 1.0-mdp.p_success])
# end

# function reward(mdp::ChainMDP, s::Int, a::Symbol)
#     if s == mdp.len
#         return 1.0
#     else 
#         return 0.0
#     end
# end


# discount(mdp::ChainMDP) = mdp.discount
# actions(::ChainMDP) = [:left, :right]

# history = sim(mdp, initial_state=1, max_steps=10) do s
#     if s == 2
#         return :left
#     end
#     return :right
# end

# solver = MCTSSolver()
# mdp = ChainMDP(4, 0.8, 0.95)



# solver = MCTSSolver(enable_tree_vis=true)
# planner = solve(solver, mdp);

### GridWorld

struct GridWorldState 
    x::Int64 # x position
    y::Int64 # y position
    done::Bool # are we in a terminal state?
end

# initial state constructor
GridWorldState(x::Int64, y::Int64) = GridWorldState(x,y,false)
# checks if the position of two states are the same
posequal(s1::GridWorldState, s2::GridWorldState) = s1.x == s2.x && s1.y == s2.y

# the grid world mdp type
type GridWorld <: MDP{GridWorldState, Symbol} # Note that our MDP is parametarized by the state and the action
    size_x::Int64 # x size of the grid
    size_y::Int64 # y size of the grid
    reward_states::Vector{GridWorldState} # the states in which agent recieves reward
    reward_values::Vector{Float64} # reward values for those states
    tprob::Float64 # probability of transitioning to the desired state
    discount_factor::Float64 # discount factor
end

# we use key worded arguments so we can change any of the values we pass in 
function GridWorld(;sx::Int64=10, # size_x
                    sy::Int64=10, # size_y
                    rs::Vector{GridWorldState}=[GridWorldState(4,3), GridWorldState(4,6), GridWorldState(9,3), GridWorldState(8,8)], # reward states
                    rv::Vector{Float64}=rv = [-10.,-5,10,3], # reward values
                    tp::Float64=0.7, # tprob
                    discount_factor::Float64=0.9)
    return GridWorld(sx, sy, rs, rv, tp, discount_factor)
end

# we can now create a GridWorld mdp instance like this:
mdp = GridWorld()
mdp.reward_states # mdp contains all the defualt values from the constructor

function POMDPs.states(mdp::GridWorld)
    s = GridWorldState[] # initialize an array of GridWorldStates
    # loop over all our states, remeber there are two binary variables:
    # done (d)
    for d = 0:1, y = 1:mdp.size_y, x = 1:mdp.size_x
        push!(s, GridWorldState(x,y,d))
    end
    return s
end;

mdp = GridWorld()
state_space = states(mdp);

POMDPs.actions(mdp::GridWorld) = [:up, :down, :left, :right];

# transition helpers
function inbounds(mdp::GridWorld,x::Int64,y::Int64)
    if 1 <= x <= mdp.size_x && 1 <= y <= mdp.size_y
        return true
    else
        return false
    end
end

inbounds(mdp::GridWorld, state::GridWorldState) = inbounds(mdp, state.x, state.y);

function POMDPs.transition(mdp::GridWorld, state::GridWorldState, action::Symbol)
    a = action
    x = state.x
    y = state.y
    
    if state.done
        return SparseCat([GridWorldState(x, y, true)], [1.0])
    elseif state in mdp.reward_states
        return SparseCat([GridWorldState(x, y, true)], [1.0])
    end

    neighbors = [
        GridWorldState(x+1, y, false), # right
        GridWorldState(x-1, y, false), # left
        GridWorldState(x, y-1, false), # down
        GridWorldState(x, y+1, false), # up
        ] # See Performance Note below
    
    targets = Dict(:right=>1, :left=>2, :down=>3, :up=>4) # See Performance Note below
    target = targets[a]
    
    probability = fill(0.0, 4)

    if !inbounds(mdp, neighbors[target])
        # If would transition out of bounds, stay in
        # same cell with probability 1
        return SparseCat([GridWorldState(x, y)], [1.0])
    else
        probability[target] = mdp.tprob

        oob_count = sum(!inbounds(mdp, n) for n in neighbors) # number of out of bounds neighbors

        new_probability = (1.0 - mdp.tprob)/(3-oob_count)

        for i = 1:4 # do not include neighbor 5
            if inbounds(mdp, neighbors[i]) && i != target
                probability[i] = new_probability
            end
        end
    end

    return SparseCat(neighbors, probability)
end;

function POMDPs.reward(mdp::GridWorld, state::GridWorldState, action::Symbol, statep::GridWorldState) #deleted action
    if state.done
        return 0.0
    end
    r = 0.0
    n = length(mdp.reward_states)
    for i = 1:n
        if posequal(state, mdp.reward_states[i])
            r += mdp.reward_values[i]
        end
    end
    return r
end;



POMDPs.n_states(mdp::GridWorld) = 2*mdp.size_x*mdp.size_y
POMDPs.n_actions(mdp::GridWorld) = 4

POMDPs.discount(mdp::GridWorld) = mdp.discount_factor;

function POMDPs.state_index(mdp::GridWorld, state::GridWorldState)
    sd = Int(state.done + 1)
    return sub2ind((mdp.size_x, mdp.size_y, 2), state.x, state.y, sd)
end
function POMDPs.action_index(mdp::GridWorld, act::Symbol)
    if act==:up
        return 1
    elseif act==:down
        return 2
    elseif act==:left
        return 3
    elseif act==:right
        return 4
    end
    error("Invalid GridWorld action: $act")
end;

POMDPs.isterminal(mdp::GridWorld, s::GridWorldState) = s.done

function MCTS.node_tag(s::GridWorldState)
    if s.done
        return "done"
    else
        return "[$(s.x),$(s.y)]"
    end
end

#MCTS.node_tag(a::GridWorldAction) = "[$a]"

n_iter = 100000
depth = 10
ec = 10.0

solver = MCTSSolver(n_iterations=n_iter,
    depth=depth,
    exploration_constant=ec,
    enable_tree_vis=true
)
mdp = GridWorld()

policy = solve(solver, mdp)
#state = initial_state(mdp, MersenneTwister(4))
state = GridWorldState(5, 5, false)

a = action(policy, state)

