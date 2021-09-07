module LiftedMaps

export LiftedMap, isondiagonal, isblockdiagonal

using LinearMaps
using LinearAlgebra

struct LiftedMap{T,TI,TJ,TM,TN} <: LinearMap{T}
    A::LinearMap{T}
    I::TI
    J::TJ
    M::TM
    N::TN
end

LiftedMap(A::AbstractMatrix, I, J, M, N) = LiftedMap(LinearMap(A), I, J, M, N)

Base.size(A::LiftedMap) = (length(A.M), length(A.N))

function LinearAlgebra.mul!(y::AbstractVector, L::LiftedMap,
    x::AbstractVector, α::Number, β::Number)

    yI = view(y, L.I)
    xJ = view(x, L.J)
    AIJ = L.A
    LinearAlgebra.mul!(yI, AIJ, xJ, α, β)
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