using Test

include("engine.jl")

positionA = Position(Int8[0 0 0 0 0 0 0;
       0 1 1 1 1 0 0;
       0 2 2 2 2 0 0;
       0 1 2 1 2 0 0;
       0 2 1 2 1 1 2;
       0 1 2 1 2 1 1],Int8(1))
positionB = Position(Int8[0 0 0 0 0 0 0;
       0 0 0 0 0 0 0;
       0 2 2 2 1 0 0;
       0 1 2 1 2 0 0;
       0 2 1 2 1 1 2;
       0 1 2 1 2 1 1],Int8(1))
positionC= Position(Int8[0 0 0 0 0 0 0;
       0 0 0 0 0 0 0;
       0 0 2 2 1 0 0;
       0 1 2 1 2 0 0;
       0 2 1 2 1 1 2;
       0 1 2 1 2 1 1],Int8(1))

@test staticValue(position0.board,Int8(1))==0
@test staticValue(positionA.board,Int8(1))==1
@test staticValue(positionA.board,Int8(2))==1
@test staticValue(positionB.board,Int8(1))==1
@test staticValue(positionB.board,Int8(2))==1
@test staticValue(positionC.board,Int8(1))==1
@test staticValue(positionC.board,Int8(2))==Float16(1/18)

bearChildren!(node0)
@test node0.children[1].position.id == "2000001000000000000000000000000000000000000"
