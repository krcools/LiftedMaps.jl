using BlockArrays
using LiftedMaps
using LinearAlgebra

ax1 = blockedrange((2,2,3))
ax2 = blockedrange((1,1,5))

block = rand(2,5)

A = LiftedMap(block, Block(1), Block(3), ax1, ax2)
B = zeros((ax1,ax2))

mul!(B, 1, A, 1, 1)

# BIJ = view(B, A.I, A.J)
# mul!(BIJ, 1, A.A, 1, 1)

@assert Matrix(A) ≈ Matrix(B)