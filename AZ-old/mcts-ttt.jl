using POMDPs,MCTS,D3Trees,POMDPToolbox
import Base.rand,Base.copy

struct TTTState 
    board::Array{Int,2}  # board
    player::Int # on move, 1 or 2
    done::Bool
end

function copy(s::TTTState)::TTTState
    TTTState(copy(s.board),s.player,s.done)
end

struct TTTMove
    x::Int
    y::Int
end

###### Symmetry

#len = mdp.size + 1
len = 4
function sym0(i::Int,j::Int)::Tuple{Int,Int}  return i,j              end
function sym1(i::Int,j::Int)::Tuple{Int,Int}  return j,len - i        end
function sym2(i::Int,j::Int)::Tuple{Int,Int}  return len - i,len - j  end
function sym3(i::Int,j::Int)::Tuple{Int,Int}  return len - j,i        end
function sym4(i::Int,j::Int)::Tuple{Int,Int}  return j,i              end
function sym5(i::Int,j::Int)::Tuple{Int,Int}  return len - i,j        end
function sym6(i::Int,j::Int)::Tuple{Int,Int}  return len - j,len - i  end
function sym7(i::Int,j::Int)::Tuple{Int,Int}  return i,len - j        end

symms    = [sym0, sym1, sym2, sym3, sym4, sym5, sym6, sym7]
invsymms = [sym0, sym3, sym2, sym1, sym4, sym5, sym6, sym7]

function get_inv(sym::Function)::Function
    ind = findfirst(symms,sym)
    return invsymms[ind]
end

function apply_sym(sym::Function,state::TTTState)::TTTState
    bd = copy(state.board)
    newbd = [ bd[sym(i,j)...] for i=1:3,j=1:3 ]
    return TTTState(newbd,state.player,state.done)
end
function apply_sym(sym::Function,m::TTTMove)::TTTMove
    x,y = sym(m.x,m.y)
    return TTTMove(x,y)
end

function apply_sym(sym::Function,state::TTTState,pol::Array{Float64,1})::Array{Float64,1}
    len = 9
    prs = collect(zip(action_list,pol))
    prs = map(p -> (apply_sym(sym,p[1]),p[2]),prs)
    sort!(prs,by=(p -> POMDPs.action_index(mdp,p[1])))
    newpol = map(second,prs)
    return newpol
end

function apply_sym(sym,e::Tuple{TTTState,Array{Float64,1},Float64})::Tuple{TTTState,Array{Float64,1},Float64}
    syminv = get_inv(sym)
    return (apply_sym(syminv,e[1]),apply_sym(sym,e[1],e[2]),e[3])
end

### ennd symms


                                           
# the grid world mdp struct
struct TTTGame <: MDP{TTTState, TTTMove} 
    size::Int
    net##::AbstractTTTmodel
end


POMDPs.actions(::TTTGame) = vec([TTTMove(i,j) for i=1:3, j=1:3 ]);

function POMDPs.actions(::TTTGame,s::TTTState)::Array{TTTMove,1}
    mvs = filter(m -> inbounds(s,m), action_list)
    # sts = map(m -> next_state(s,m)[1],mvs)
    # prs = collect(zip(sts,mvs))
    # prs = unique(p -> choose_rep(p[1])[1],prs)
    # return map(second,prs)
    return mvs
end

# transition helpers
function inbounds(state::TTTState,x::Int,y::Int)::Bool
    return 1 <= x <= 3 && 1 <= y <= 3 && state.board[x,y]==0
end

inbounds(state::TTTState, move::TTTMove) = inbounds(state, move.x, move.y);

function next_state(state::TTTState, move::TTTMove)::Tuple{TTTState,Int64}
    if !inbounds(state, move)
        return state, 0
    else
        player = state.player
        newboard = copy(state.board)
        newboard[move.x,move.y] = player
        newstate = TTTState(newboard,3 - player,false)
        checkwin = check_lines(newboard,player)
        #println(checkwin)
        moves = filter(m -> inbounds(newstate,m), action_list)
        #println(moves)
        newterm = checkwin != 0 || length(moves) == 0
        newstate = TTTState(newboard,3 - player,newterm)
        return newstate, checkwin
    end
end;

function POMDPs.transition(mdp::TTTGame, state::TTTState, move::TTTMove)::POMDPToolbox.SparseCat{Tuple{TTTState},Tuple{Float64}}
    nst,_ = next_state(state,move)
    #newstate,_ = choose_rep(nst)
    return SparseCat((nst,), (1.0,))
end;

lines = [ [ [1, 1], [2, 2], [3, 3] ],
          [ [3, 1], [2, 2], [1, 3] ],
          [ [3, 1], [3, 2], [3, 3] ],
          [ [2, 1], [2, 2], [2, 3] ],
          [ [1, 1], [1, 2], [1, 3] ],
          [ [3, 3], [2, 3], [1, 3] ],
          [ [3, 2], [2, 2], [1, 2] ],
          [ [3, 1], [2, 1], [1, 1] ] ];

function check_constant( l::Array{Int,1})::Bool
    h, t = l[1], l[2:3]
    test = map( x -> x==h, t)
    return reduce(&,test)
end

function check_one_line(b::Array{Int,2},l::Array{Array{Int,1},1})::Int
    pces = map(c -> b[c...], l)
    check_constant(pces) ? (return pces[1]) : (return 0)
end

function check_lines(b::Array{Int,2},p::Int)::Int
    q = 3 - p
    tests = map( l -> check_one_line(b,l), lines)
    return (q in tests) ? q : (p in tests) ? p : 0
end  

function choose_rep(state::TTTState)::Tuple{TTTState,Function}
    ss = map( s -> apply_sym(s,state), symms)
    rep = first(sort(ss,by=hash))
    ind = findfirst(ss,rep)
    inv = invsymms[ind]
    return rep,inv
end


function POMDPs.reward(mdp::TTTGame, state::TTTState, move::TTTMove)::Int
    p = state.player
    _,winner = next_state(state,move)
    return (winner==p) ? 1 : (winner == 3 - p) ? -1 : 0
end;
function POMDPs.reward(mdp::TTTGame, state::TTTState, move::TTTMove, statep::TTTState)::Int
    return -reward(statep)
end;
function reward(state::TTTState)::Int64
    p = state.player
    winner = check_lines(state.board,p)
    return (winner==p) ? 1 : (winner == 3 - p) ? -1 : 0
end


POMDPs.n_states(::TTTGame) = 3^9
POMDPs.n_states() = 3^9
POMDPs.n_actions(mdp::TTTGame) = 9
POMDPs.discount(mdp::TTTGame) = -1;

function POMDPs.state_index(state::TTTState)
    sd = Int(state.done + 1)
    pces = [ 1 + state.board[ind2sub((3,3),i)...] for i=1:9 ]
    return sub2ind((3,3,3,3,3,3,3,3,3), pces...)
end

POMDPs.state_index(::TTTGame, state::TTTState) =  state_index(state)

function POMDPs.action_index(::TTTGame, move::TTTMove)
    return sub2ind((3,3), move.x, move.y)
end;

function action_index_inv(mdp::TTTGame, ind::Int64)::TTTMove
    x,y = ind2sub((3,3), ind)
    return TTTMove(x,y)
end;

# checks if the position of two states are the same
POMDPs.isequal(s1::TTTState, s2::TTTState) = s1.board == s2.board && s1.player == s2.player

POMDPs.hash(state::TTTState,n::Int64) = mod(state_index(state) + n,
                                            n_states() )
POMDPs.hash(state::TTTState) = POMDPs.hash(state,0)


    
POMDPs.isterminal(::TTTGame, s::TTTState)::Bool = s.done 

function MCTS.node_tag(s::TTTState)
    return "$s.board"
end

# function rand(rng::AbstractRNG,state::TTTState)::TTTMove
#     moves = filter(m -> inbounds(state,m), action_list)
#     return rand(rng,moves)
# end

MCTS.node_tag(m::TTTMove) = "[$m]"

n_iter = 100
depth = 9
ec = 4.0

# function get_value(outs,state::TTTState)
#     fun = outs[2]
#     return fun(vec(state.board),state.player)[1]
# end

function init_value(mdp::TTTGame,state::TTTState,move::TTTMove)::Float64
    newstate,_ = next_state(state,move)
    # println(state)
    # println(move)
    predict(newstate)[2][1]
end

solver = MCTSSolver(n_iterations=n_iter,
                    depth=depth,
                    exploration_constant=ec,
                    enable_tree_vis=true,
                    init_Q=init_value
                    )
    ## must define makeNN!!

#action(policy, startpos)

printpos(s::TTTState) = map(println,[ [s.board[i,j] for j=1:3] for i=1:3 ])
printpol(p::Array{Float64,1}) = map(println,[ [p[3(i-1)+j] for j=1:3] for i=1:3 ])

function play_against(policy,sp)
    pos = copy(startpos)
    rep = pos
    p = sp
    inv = sym0
    while !pos.done
        if p==1
            println("Your move: ")
            mstr = readline()
            m = TTTMove(parse(Int,mstr[1]),parse(Int,mstr[3]))
            pos,+ = next_state(pos,m)
        else
            m = action(policy,rep)
            nst,_ = next_state(rep,m)
            pos = apply_sym(inv,nst)
        end
        println(m)
        rep,inv = choose_rep(pos)
        #println(pos)
        p = 3 - p
        printpos(pos)
    end
end


# hashred(a::Int,b::Int)::Int = 3*a + b

# state_hash(s::TTTState)::Int = reduce(hashred,vec(s.board))


    #### debugging


startpos = TTTState([0 0 0; 0 0 0; 0 0 0],1,false)
testpos1 = choose_rep(TTTState([1 1 1; 2 2 0; 0 0 0],2,true))[1]
testpos2 = choose_rep(TTTState([1 1 0; 2 2 2; 1 0 0],1,true))[1]
testpos3 = choose_rep(TTTState([0 0 0; 0 1 0; 0 0 0],2,false))[1]
testpos4 = choose_rep(TTTState([0 2 0; 0 1 0; 0 0 0],1,false))[1]
testpos5 = choose_rep(TTTState([2 0 0; 0 1 0; 0 0 1],2,false))[1]
testpos6 = choose_rep(TTTState([2 1 0; 2 1 2; 1 2 1],1,false))[1]
testpos7 = choose_rep(TTTState([2 1 0; 2 1 0; 1 2 1],2,false))[1]
testpos8 = choose_rep(TTTState([1 0 0; 1 2 0; 2 0 0],1,false))[1]
testpos9 = choose_rep(TTTState([1 2 1; 0 1 0; 0 2 0],2,false))[1]
;
