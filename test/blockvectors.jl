@testitem "blockvectors" begin
    using BlockArrays
    using LinearAlgebra

    N1, N2 = 20_000, 30_000

    A = rand(N1, N1)
    ax = blockedrange([N1,N2])
    L = LiftedMap(A, Block(1), Block(1), ax, ax)

    a = rand(N1+N2)
    b = BlockVector(a, [N1,N2])
    p = BlockedVector(a, [N1,N2])

    axb = axes(b,1)
    axp = axes(p,1)

    I = Block(1)
    B = axb[I]
    P = axp[I]
    @test B == P

    La = L * a
    Lb = L * b
    Lp = L * p

    @time L * a
    @time L * b
    @time L * p

    @test norm(La - Lb) < 1e-10
    @test norm(La - Lp) < 1e-10

    aI = view(a, B)
    pI = view(p, P)
    bI = view(b, B)

    @show typeof(bI)
    @show typeof(pI)
end