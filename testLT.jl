1+1

using Winston

display( imagesc( rand( 2,2 )))

display( imagesc( rand( 4,4 )))

imagesc( rand( 8,8 ))

showFlowCA( repeat( "1100", 20 ), 100 )

include( "eden.jl")

evolution!( eden )
