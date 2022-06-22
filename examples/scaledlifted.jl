using LinearMaps

struct CustomMap{T} <: LinearMap{T}
    A::Matrix{T}
end

LinearMaps.MulStyle(C::CustomMap) = LinearMaps.FiveArg()
LinearMaps.size(C::CustomMap) = size(C.A)
LinearMaps.mul!(y::AbstractVector, C::CustomMap, x::AbstractVector, a::Number, b::Number) = LinearMaps.mul!(y, C.A, x, a, b)
function LinearMaps.mul!(Y::AbstractMatrix, C::CustomMap, X::AbstractMatrix, a::Number, b::Number)
    println("User supplied BLAS-3")
    LinearMaps.mul!(Y, C.A, X, a, b)
end

A = rand(3,3)
C = CustomMap(A)
S = 2C

X = rand(3,4)
Y = zeros(3,4)

LinearMaps.mul!(Y,C,X,1,1);
LinearMaps.mul!(Y,S,X,1,1);
