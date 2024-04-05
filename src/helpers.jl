module TestHelpers

using Test

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

SimplicialSets.:⋄(ys::BasicSimplex...) = BasicSimplex(⋄(map(y -> y.x, ys)...))
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
# twisting function
#

export test_twf

"""
    test_twf(f, x::AbstractSimplex)

Test that `f` satisfies the defining identites for a twisting functions for the argument `x`.
"""
function test_twf(f, x::AbstractSimplex)
    n = dim(x)
    y = f(s(x, 0))
    @test dim(y) == n && isone(y)
    if n == 0
        @test_throws Exception f(x)
    else
        y = f(x)
        @test dim(y) == n-1
        for k in 0:n-1
            @test s(y, k) == f(s(x, k+1))
            n == 1 && continue
            if k == 0
                @test d(y, k) == inv(f(d(x, 0))) ⋄ f(d(x, 1))
            else
                @test d(y, k) == f(d(x, k+1))
            end
        end
    end
end

#
# twisting cochain
#

export test_twc

"""
    test_twc(f, a)

Test that `f` satisfies the twisting cochain identity for the argument `a`,
```
diff(f(a)) + f(diff(a)) == coprod(a) |> Tensor(f, f) |> TensorSplat(*)
```
"""
function test_twc(f, a)
    b1 = diff(f(a)) + f(diff(a))
    b2 = coprod(a) |> Tensor(f, f) |> TensorSplat(*)
    @test b1 == b2
end

end # module TestHelpers
