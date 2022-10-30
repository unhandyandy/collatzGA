
include("dpll.jl")

#create deck
#ptot = # people playing
function createDeck( ptot::Int )
    res = Set{(Int,Int)}()
    for p = 0 : ptot - 1, j = 0 : 3
        push!( res, ( p, j ) )
    end
    global numProps = ptot * ptot * 4
    global numPlayers = ptot
    res
end

function createHands( deck::Set{(Int,Int)} )
    dict = Dict{(Int,Int),Int}()
    for c in deck
        dict[ c ] = c[1]
    end
    dict
end

function getHand( dict::Dict{(Int,Int),Int}, p::Int )
    hand = Set{(Int,Int)}()
    for (k,v) in dict
        if v == p
            push!( hand, k )
        end
    end
    hand
end

function codeAtomicProp( card::(Int,Int), suit::Int )
    p, r = card
    4 * numPlayers * p + numPlayers * r + suit + 1
end

function decodeAtomicProp( code::Int )
    q, s = divrem( code - 1, numPlayers )
    p, r = divrem( q, 4 )
    ( ( p, r ), s )
end

function makeCardAxioms( card::(Int,Int) )
    p, r = card
    lst4 = [ ( card, i ) for i = 0 : numPlayers - 1 ]
    codelst4 = map( x -> codeAtomicProp( x... ), lst4 )
    com4 = "Suited"
    dsj4 = disjunct( codelst4, [ true for i = 0 : numPlayers - 1 ], numProps, com4 )
    axset = Set{disjunct}( {dsj4} )
    for i = 0 : numPlayers - 2, j = i + 1 : numPlayers - 1
        com = string(card) * "UniqueSuit" * string(i) * string(j)
        dsj = disjunct( [ codeAtomicProp( card, i ), 
                         codeAtomicProp( card, j ) ],
                       [ false, false ],
                       numProps,
                       com )
        push!( axset, dsj )
    end
    axset
end

function makeCardAxiomsAll( deck::Set{(Int,Int)} )
    res = Set{disjunct}()
    for c in deck
        union!( res, makeCardAxioms( c ) )
    end
    res
end

# p = player #
function makeDealAxiom( p::Int, suit::Int )
    cardlist = [ ( p, i ) for i = 0 : 3 ]
    proplist = Int[ codeAtomicProp( c, suit ) for c in cardlist ]
    signlist = [ false for i = 1 : 4 ]
    #@hiya proplist, signlist
    com = "Misdeal"
    disjunct( proplist, signlist, numProps, com )
end

function makeDealAxiomsAll( )
    res = Set{disjunct}()
    for p = 0 : numPlayers - 1, s = 0 : numPlayers - 1
        push!( res, makeDealAxiom( p, s ) )
    end
    res
end

function makeSuitAxioms( deck )
    axioms = Set{disjunct}()
    for i = 0 : numPlayers - 1
        union!( axioms, djsetSuitMinSize( deck, i ) )
    end
    axioms
end

function makeAxiomsInit( deck )
    union( makeDealAxiomsAll( ),
           makeCardAxiomsAll( deck ),
           makeSuitAxioms( deck ) )
end

function ask( hands::Dict{(Int,Int),Int}, p::Int, suit::Int )
    hand = getHand( hands, p )
    com = string(p) * "Ask" * string(suit)
    disjunct( Int[ codeAtomicProp( c, suit ) for c in hand ],
              [ true for c in hand ],
             numProps,
             com )
end

function ask( axioms::Set{disjunct}, hands::Dict{(Int,Int),Int}, p::Int, suit::Int )
    newaxioms = clone( axioms )
    push!( newaxioms, ask( hands, p, suit ) )
    if !dpll( newaxioms )
        printLoss( p )
        return false, simplify( newaxioms )
    else
        return true, simplify( newaxioms )
    end
end

function propSameSuit( cset::Set{(Int,Int)}, suit::Int )
    proplist = [ codeAtomicProp( c, suit ) for c in cset ]
    newaxioms = Set{disjunct}()
    com = "SameSuit" * string(suit)
    for pc in proplist
        push!( newaxioms, disjunct( [ pc ], [ true ], numProps, com ) )
    end
    newaxioms
end

function propHasSuit( cset::Set{(Int,Int)}, suit::Int, com::String=("HasSuit" * string(suit)) )
    proplist = Int[ codeAtomicProp( c, suit ) for c in cset ]
    #@hiya proplist
    disjunct( proplist, [ true for c in cset ], numProps, com )
end

function djsetSuitNone( cset::Set{(Int,Int)}, suit::Int )
    res = Set{disjunct}()
    com = "OmitsSuit" * string(suit)
    for  c in cset
        push!( res, disjunct( [ codeAtomicProp( c, suit ) ], [ false ], numProps, com ) )
    end
    res
end

function djsetSuitMinSize( cset::Set{(Int,Int)}, suit::Int, size::Int=4 )
    propsize = length( cset ) - size + 1
    propsets = subsets( cset, propsize )
    res = Set{disjunct}()
    com = "FullSuit" * string(suit)
    for ps in propsets
        push!( res, propHasSuit( ps, suit, com ) )
    end
    res
end

function layoff( hands::Dict{(Int,Int),Int}, p::Int, suit::Int )
    hand = getHand( hands, p )
    djsetSuitMinSize( hand, suit )
end

function layoff( axioms::Set{disjunct}, hands::Dict{(Int,Int),Int}, p::Int, suit::Int )
    hand = getHand( hands, p )
    newhands = copy( hands )
    sets = subsets( hand, 4 )
    for s in sets
        newaxioms = union( axioms, propSameSuit( s, suit ) )
        if dpll( newaxioms )
            for c in s; newhands[ c ] = -1; end
            return true, simplify( newaxioms ), newhands
        end
    end
    printLoss( p )
    return false, simplify( axioms ), newhands
end

function printLoss( p::Int )
    @printf( "\n Player %d loses! \n", p )
end
function printWin( p::Int )
    @printf( "\n Player %d wins! \n", p )
end

function findSuit( axioms::Set{disjunct}, cset::Set{(Int,Int)}, suit::Int )
    for c in cset
        notc = disjunct( [ codeAtomicProp( c, suit ) ], [ false ], numProps )
        newaxioms = simplify( union( axioms, Set{disjunct}({notc}) ) )
        if !dpll( newaxioms )
            cax = disjunct( [ codeAtomicProp( c, suit ) ], [ true ], numProps )
            newaxioms = union( axioms, Set{disjunct}({cax}) )
            return c, newaxioms
        end
    end
    for c in cset
        cax = disjunct( [ codeAtomicProp( c, suit ) ], [ true ], numProps )
        newaxioms = simplify( union( axioms, Set{disjunct}({cax}) ) )
        if dpll( newaxioms )
            return c, newaxioms
        end
    end
end
        

# p = asker, q = answerer
function answer( askaxioms::Set{disjunct}, hands::Dict{(Int,Int),Int}, p::Int, suit::Int, q::Int, response::Bool )
    #askaxioms = ask( axioms, hands, p, suit )
    newhands = copy( hands )
    if response
        rsp, respaxioms = ask( askaxioms, hands, q, suit )
        if !rsp
            printLoss( q )
            return false, respaxioms, newhands
        else
            handq = getHand( hands, q )
            card, newaxioms = findSuit( respaxioms, handq, suit )
            #@hiya card
            newhands[ card ] = p
            return true, newaxioms, newhands
        end
    else
        handq = getHand( hands, q )
        respaxioms = simplify( union( askaxioms, djsetSuitNone( handq, suit ) ) )
        if !dpll( respaxioms )
            printLoss( q )
            return false, respaxioms, newhands
        else
            return true, respaxioms, newhands
        end
    end
end

type uiModel
    cardNumbers::Array{Int,1}
    suit::Int
    player::Int
    phase::Int
    respondent::Int
    answer::Bool
    over::Bool
    winner::Int
end

type gameState
    axioms::Set{disjunct}
    hands::Dict{(Int,Int),Int}}
    numberPlayers::Int
    activePlayer::Int
end

