using LiftedMaps
using BlockArrays
using Test

import LiftedMaps.LinearMaps

@testset "LiftedMaps.jl" begin
    A = rand(10,10)
    B = rand(3,3)

    U = Base.OneTo(13)
    V = Base.OneTo(13)

    I = 1:10
    J = 11:13
    
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

    x = rand(13)
    x1 = x[I]
    x2 = x[J]

    y1 = A * x1
    y2 = 2 * B * x2
    y = [y1; y2]

    z = D * x
    @test y ≈ z

    Amatrix = Matrix(Alifted)
    @test Amatrix isa Matrix
    @test Amatrix[I,I] == A

    ax1 = blockedrange([10,3])
    ax2 = blockedrange([10,3])
    Alifted = @inferred LiftedMap(A,Block(1),Block(1),ax1,ax2)
    Blifted = @inferred LiftedMap(B,Block(2),Block(2),ax1,ax2)
    @test axes(Alifted) === (ax1,ax2)
    @test axes(Blifted) === (ax1,ax2)

    D = Alifted + 2 * Blifted
    z = @inferred D * x
end