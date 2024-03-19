module TestHelpers

using ..SimplicialSets
using StructEqualHash: @struct_equal_hash
using LinearCombinations
using LinearCombinations: diff
using SimplicialSets: d, s

export BasicSimplex, undo_basic

#
# BasicSimplex
#

# used to test that functions only use the basic operations dim, d, s

struct BasicSimplex{T<:AbstractSimplex} <: AbstractSimplex
    x::T
end

# Base.:(==)(y::BasicSimplex, z::BasicSimplex) = y.x == z.x
# Base.hash(y::BasicSimplex, h::UInt) = hash(y.x, h)

@struct_equal_hash BasicSimplex{T} where T
Base.copy(y::BasicSimplex) = BasicSimplex(copy(y.x))

SimplicialSets.dim(y::BasicSimplex) = dim(y.x)
SimplicialSets.d(y::BasicSimplex, k::Integer) = BasicSimplex(d(y.x, k))
SimplicialSets.s(y::BasicSimplex, k::Integer) = BasicSimplex(s(y.x, k))

SimplicialSets.:⋅(ys::BasicSimplex...) = BasicSimplex(⋅(map(y -> y.x, ys)...))
Base.:/(y::BasicSimplex, z::BasicSimplex) = BasicSimplex(y.x/z.x)
Base.inv(y::BasicSimplex) = BasicSimplex(inv(y.x))
Base.one(::Type{BasicSimplex{T}}, n...) where T = BasicSimplex(one(T, n...))
Base.one(y::BasicSimplex, n...) = BasicSimplex(one(y.x, n...))

undo_basic(x::AbstractSimplex) = x
undo_basic(y::BasicSimplex) = y.x
undo_basic(x::ProductSimplex) = ProductSimplex((undo_basic(y) for y in x)...)
undo_basic(x::AbstractTensor) = Tensor((undo_basic(y) for y in x)...)

# undo_basic(a::Linear) = Linear(undo_basic(x) => c for (x, c) in a)
@linear undo_basic

#
# twisting cochain
#

export check_twc

"""
    check_twc(f, a)

Check that `f` satisfies the twisting cochain identity for the argument `a`,
```
diff(f(a)) + f(diff(a)) == coprod(a) |> Tensor(f, f) |> TensorSplat(*)
```
The returned value is the left-hand side minus the right-hand side.
"""
function check_twc(f, a)
    b1 = diff(f(a)) + f(diff(a))
    b2 = coprod(a) |> Tensor(f, f) |> TensorSplat(*)
    b1-b2
end

end # module TestHelpers
