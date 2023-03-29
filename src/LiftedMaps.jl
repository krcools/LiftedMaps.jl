module LiftedMaps

export LiftedMap, isondiagonal, isblockdiagonal

using LinearMaps
using LinearAlgebra
# using BlockArrays

struct LiftedMap{T,TI,TJ,TM,TN} <: LinearMap{T}
    A::LinearMap{T}
    I::TI
    J::TJ
    M::TM
    N::TN
end

LiftedMap(A::AbstractMatrix, I, J, M, N) = LiftedMap(LinearMap(A), I, J, M, N)
LinearMaps.MulStyle(A::LiftedMap) = LinearMaps.FiveArg()

Base.size(A::LiftedMap) = (length(A.M), length(A.N))
Base.axes(A::LiftedMap) = (A.M, A.N)
Base.axes(A::LiftedMap, i::Int) = axes(A)[i]

LinearAlgebra.adjoint(A::LiftedMap) = LiftedMap(adjoint(A.A), A.J, A.I, A.N, A.M)
LinearAlgebra.transpose(A::LiftedMap) = LiftedMap(transpose(A.A), A.J, A.I, A.N, A.M)

function LinearMaps._unsafe_mul!(y::AbstractVector, L::LiftedMap,
    x::AbstractVector, α::Number=true, β::Number=false)

    y .*= β
    temp = y

    I = L.I
    J = L.J

    P = getindex(axes(L,1), I)
    Q = getindex(axes(L,2), J)

    # yI = view(y, I)
    # xJ = view(x, J)
    yI = view(y, P)
    xJ = view(x, Q)
    AIJ = L.A


    LinearMaps._unsafe_mul!(yI, AIJ, xJ, α, true)
    return temp
end


function LinearMaps._unsafe_mul!(Y::AbstractMatrix, X::LiftedMap, c::Number, a::Number=true, b::Number=false)

    LinearAlgebra.rmul!(Y, b)
    temp = Y

    I = X.I
    J = X.J

    P = getindex(axes(X,1), I)
    Q = getindex(axes(X,2), J)

    YIJ = view(Y, P, Q)
    XIJ = X.A
    LinearMaps._unsafe_mul!(YIJ, X.A, c, a, true)
    return temp
end

function isondiagonal(A::LiftedMaps.LiftedMap)
    return A.I == A.J
end

function isondiagonal(A::LinearMaps.ScaledMap)
    return isondiagonal(A.lmap)
end

function isondiagonal(A::LinearMaps.LinearMap)
    return false
end

function isblockdiagonal(A::LinearMaps.LinearCombination)
    for map in A.maps
        isondiagonal(map) || return false
    end
    return true
end

function isblockdiagonal(A::LinearMaps.LinearMap)
    return false
end

end