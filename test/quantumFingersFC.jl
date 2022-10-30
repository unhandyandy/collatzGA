
module factCheckQF

using FactCheck

facts("2 Players") do
    include( "/home/dabrowsa/lang/julia/mine/quantumFingers.jl" )
    deck = createDeck( 2 )
    axioms0 = makeAxiomsInit( deck )
    @fact numProps => 16
    @fact numPlayers => 2
    @fact dpll( axioms0 ) => true
    hands0 = createHands( deck )
    rsp, axioms1 = ask( axioms0, hands0, 0, 0 )
    @fact rsp => true
    @fact dpll( axioms1 ) => true
    rsp, axioms2, hands2 = answer( axioms1, hands0, 0, 0, 1, true )
    @fact rsp => true
    @fact dpll( axioms2 ) => true
    @fact length( getHand( hands2, 0 ) ) => 5
    @fact length( getHand( hands2, 1 ) ) => 3
    rsp, axioms3 = ask( axioms2, hands2, 0, 1 )
    @fact rsp => true
    @fact dpll( axioms3 ) => true
    rsp, axioms4, hands4 = answer( axioms3, hands2, 0, 1, 1, false )
    @fact rsp => false
    @fact dpll( axioms4 ) => false
    rsp, axioms4, hands4 = answer( axioms3, hands2, 0, 1, 1, true )
    @fact rsp => true
    @fact dpll( axioms4 ) => true
    @fact length( getHand( hands4, 0 ) ) => 6
    @fact length( getHand( hands4, 1 ) ) => 2
    rsp, axioms5, hands5 = layoff( axioms4, hands4, 0, 0 )
    @fact rsp => true
    @fact length( getHand( hands5, 0 ) ) => 2
    rsp, axioms6 = ask( axioms5, hands5, 0, 1 )
    @fact rsp => true
    rsp, axioms7, hands7 = answer( axioms6, hands5, 0, 1, 1, true )
    @fact rsp => true
    @fact length( getHand( hands7, 0 ) ) => 3
    @fact length( getHand( hands7, 1 ) ) => 1
    rsp, axioms8 = ask( axioms7, hands7, 0, 1 )
    @fact rsp => true
    rsp, axioms9, hands9 = answer( axioms8, hands7, 0, 1, 1, true )
    @fact rsp => true
    @fact length( getHand( hands9, 0 ) ) => 4
    @fact length( getHand( hands9, 1 ) ) => 0
    @fact dpll( axioms9 ) => true
    rsp, axioms10, hands10 = layoff( axioms9, hands9, 0, 0 )
    @fact rsp => false
    rsp, axioms10, hands10 = layoff( axioms9, hands9, 0, 1 )
    @fact rsp => true
    @fact length( getHand( hands10, 0 ) ) => 0
    @fact length( getHand( hands10, 1 ) ) => 0
end

facts("3 Players") do
    include( "/home/dabrowsa/lang/julia/mine/quantumFingers.jl" )
    deck = createDeck( 3 )
    axioms0 = makeAxiomsInit( deck )
    @fact numProps => 36
    @fact numPlayers => 3
    @fact dpll( axioms0 ) => true
    hands0 = createHands( deck )
    rsp, axioms1 = ask( axioms0, hands0, 0, 0 )
    #@fact rsp => true
    @fact dpll( axioms1 ) => true
    rsp, axioms2, hands2 = answer( axioms1, hands0, 0, 0, 1, false )
    @fact rsp => true
    #@fact dpll( axioms2 ) => true
    @fact length( getHand( hands2, 0 ) ) => 4
    @fact length( getHand( hands2, 1 ) ) => 4    
    @fact length( getHand( hands2, 2 ) ) => 4    
    rsp, axioms2, hands2 = answer( axioms1, hands0, 0, 0, 1, true )
    @fact rsp => true
    #@fact dpll( axioms2 ) => true
    @fact length( getHand( hands2, 0 ) ) => 5
    @fact length( getHand( hands2, 1 ) ) => 3
    @fact length( getHand( hands2, 2 ) ) => 4    
    rsp, axioms3 = ask( axioms2, hands2, 0, 1 )
    @fact rsp => true
    #@fact dpll( axioms3 ) => true
    rsp, axioms4, hands4 = answer( axioms3, hands2, 0, 1, 1, false )
    @fact rsp => true
    #@fact dpll( axioms4 ) => true
    rsp, axioms4, hands4 = answer( axioms3, hands2, 0, 1, 1, true )
    @fact rsp => true
    #@fact dpll( axioms4 ) => true
    @fact length( getHand( hands4, 0 ) ) => 6
    @fact length( getHand( hands4, 1 ) ) => 2
    rsp, axioms5, hands5 = layoff( axioms4, hands4, 0, 0 )
    @fact rsp => true
    @fact length( getHand( hands5, 0 ) ) => 2
    @fact length( getHand( hands5, 1 ) ) => 2
    @fact length( getHand( hands5, 2 ) ) => 4
    rsp, axioms6, hands6 = layoff( axioms5, hands5, 0, 1 )
    @fact rsp => false
    rsp, axioms6 = ask( axioms5, hands5, 0, 1 )
    @fact rsp => true
    rsp, axioms7, hands7 = answer( axioms6, hands5, 0, 1, 1, false )
    @fact rsp => true
    @fact length( getHand( hands7, 0 ) ) => 2
    @fact length( getHand( hands7, 1 ) ) => 2
    @fact length( getHand( hands7, 2 ) ) => 4
    rsp, axioms8 = ask( axioms7, hands7, 1, 0 )
    @fact rsp => false
    rsp, axioms8 = ask( axioms7, hands7, 1, 1 )
    @fact rsp => false
    rsp, axioms8 = ask( axioms7, hands7, 1, 2 )
    @fact rsp => true
    rsp, axioms9, hands9 = answer( axioms8, hands7, 1, 2, 2, false )
    @fact rsp => false
    rsp, axioms9, hands9 = answer( axioms8, hands7, 1, 2, 2, true )
    @fact rsp => true
    @fact length( getHand( hands9, 0 ) ) => 2
    @fact length( getHand( hands9, 1 ) ) => 3
    @fact length( getHand( hands9, 2 ) ) => 3
    
end

end
