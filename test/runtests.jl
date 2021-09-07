using LiftedMaps
using Test

import LiftedMaps.LinearMaps

@testset "LiftedMaps.jl" begin
    A = rand(10,10)
    B = rand(3,3)

    U = 1:13
    V = 1:13
    
    Alifted = LiftedMap(A, 1:10, 1:10, U, V)
    Blifted = LiftedMap(B, 11:13, 11:13, U, V)

    @test size(Alifted) == (13,13)
    @test size(Blifted) == (13,13)

    C = Alifted + Blifted

    @test size(C) == (13,13)
    @test C isa LinearMaps.LinearMap
    @test C isa LinearMaps.LinearCombination
    @test isblockdiagonal(C)

    D = Alifted + 2 * Blifted
    @test D isa LinearMaps.LinearMap
    @test D isa LinearMaps.LinearCombination
    @test isblockdiagonal(D)

end
