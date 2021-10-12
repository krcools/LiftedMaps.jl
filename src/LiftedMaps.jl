module LiftedMaps

export LiftedMap, isondiagonal, isblockdiagonal

using LinearMaps
using LinearAlgebra
using BlockArrays

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

function LinearAlgebra.mul!(y::AbstractVector, L::LiftedMap,
    x::AbstractVector, α::Number, β::Number)

    bvy = PseudoBlockVector(y, blocksizes(L.M)...)
    bvx = PseudoBlockVector(x, blocksizes(L.N)...)

    yI = view(bvy, L.I)
    xJ = view(bvx, L.J)
    AIJ = L.A

    y .*= β
    LinearAlgebra.mul!(yI, AIJ, xJ, α, 1)
    return y
end



function LinearAlgebra.mul!(y::AbstractVector, L::LiftedMap, x::AbstractVector)

    bvy = PseudoBlockVector(y, blocksizes(L.M)...)
    bvx = PseudoBlockVector(x, blocksizes(L.N)...)

    yI = view(bvy, L.I)
    xJ = view(bvx, L.J)
    AIJ = L.A

    fill!(y,0)
    LinearAlgebra.mul!(yI, AIJ, xJ)
    return y
end

function Base.:(*)(A::LiftedMap, x::AbstractVector)
    axes(A,2) == axes(x,1) ||
        throw(DimensionMismatch("second dimension of left factor, $(size(A, 2)), " *
                "does not match first dimension of right factor, $(size(x, 1))"))

    T = promote_type(eltype(A), eltype(x))
    y = PseudoBlockVector{T}(undef, BlockArrays.blocksizes(A,1))
    fill!(y, zero(T))
    LinearAlgebra.mul!(y, A, x)
end


function Base.Matrix{T}(A::LiftedMap) where {T}
    # M = zeros(T, size(A))
    # M = PseudoBlockMatrix{T}(undef, BlockArrays.blocksizes(A)...)
    M = similar(Array{T}, axes(A))
    fill!(M,0)
    m = Matrix(A.A)
    M[A.I,A.J] .= m
    return Matrix{T}(M)
end

function Base.Matrix{T}(A::LinearMaps.ScaledMap{<:Any, <:Any, <:LiftedMap}) where {T}
    M = Matrix{T}(A.lmap)
    M .*= A.λ
    return M
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