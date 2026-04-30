module TestHelpers

using Test

using ..SimplicialSets
using StructEqualHash: @struct_equal_hash
using LinearCombinations
using LinearCombinations: diff
using SimplicialSets: d, s

export BasicSimplex, undo_basic, test_simplex, test_groupsimplex, test_twf, test_twc

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
# simplex
#

function test_simplex(x::AbstractSimplex, n)
    @test dim(x) isa Integer
    @test n == dim(x) >= 0

    @test hash(x) isa UInt
    xc = @inferred copy(x)
    @test x == xc
    @test hash(x) == hash(xc)

    @test_throws Exception d(x, -1)
    @test_throws Exception d(x, n+1)
    n == 0 && @test_throws Exception d(x, 0)
    @test_throws Exception s(x, -1)
    @test_throws Exception s(x, n+1)

    for i in 0:n
        if n >= 1
            y = @inferred d(x, i)
            @test dim(y) == n-1 && typeof(y) == typeof(x)
        end
        y = @inferred s(x, i)
        @test dim(y) == n+1 && typeof(y) == typeof(x)
        @test @inferred isdegenerate(y)
        @test @inferred isdegenerate(y, i)
        !isdegenerate(x) && @test all(0:n) do j
            @inferred(isdegenerate(y, j)) == (j == i)
        end
    end
    if n >= 2
        for j in 0:n, i in 0:j-1
           @test d(d(x, j), i) == d(d(x, i), j-1)
        end
    end
    for j in 0:n, i in 0:j
        @test s(s(x, j), i) == s(s(x, i), j+1)
    end
    for j in 0:n, i in 0:j-1
        @test d(s(x, Int8(j)), BigInt(i)) == s(d(x, Int16(i)), Int32(j-1))
    end
    for j in 0:n
        @test d(s(x, j), j) == d(s(x, j), j+1) == x
    end
    for j in 0:n, i in j+2:n+1
        @test d(s(x, Int16(j)), Int8(i)) == s(d(x, Int32(i-1)), BigInt(j))
    end

    # test d(x, kk)
    for R in (Int8, Int, BigInt)
        kk = R[k for k in 0:n if rand(Bool)]
        length(kk) == n+1 && popfirst!(kk)
        @test d(x, kk) == undo_basic(d(BasicSimplex(x), kk))
    end

    # test s(x, kk)
    for R in (Int8, Int, BigInt)
        kk = R[]
        l = 0
        for i in 0:div(n, 3)
            k = rand(l:n+i)
            push!(kk, k)
            l = k+1
        end
        @test s(x, kk) == undo_basic(s(BasicSimplex(x), kk))
    end

    # test r(x, kk)
    @test_throws Exception r(x, [])
    @test_throws Exception r(x, [(1,2)])
    for R in (Int8, Int, BigInt)
         kk = UnitRange{R}[]
         k2 = 0
         while k2 <= (n <= 2 ? n-1 : n-2)
            k1 = rand(k2:n)
            k2 = rand(k1:n)
            push!(kk, R(k1):R(k2))
        end
        if !isempty(kk)
            y = SimplicialSets.r(x, kk)
            @test dim(y) == sum(map(length, kk))-1
            @test y == undo_basic(SimplicialSets.r(BasicSimplex(x), kk))
        end
    end
end

#
# group simplex
#

function test_groupsimplex(x::T, is_commutative) where T <: AbstractSimplex
    n = dim(x)

    onex = @inferred one(x)
    test_simplex(onex, n)
    @test isone(onex)
    @test one(x, Int8(n+1)) == one(T, BigInt(n+1))
    @test one(T) == one(T, 0)
    @test_throws Exception one(x, -1)
    for i in 0:n
        n > 0 && @test d(onex, i) == one(x, n-1)
        @test s(onex, i) == one(x, n+1)
    end

    @test ⋄(x) == x
    @test x ⋄ onex == x == onex ⋄ x
    @test isone(x ⋄ inv(x)) && isone(inv(x) ⋄ x)
    @test_throws Exception x⋄one(T, n+1)

    #=
    @test_broken @inferred(x^0) == onex
    @test_broken @inferred(x^1) == x
    @test_broken @inferred(x^3) == x⋄x⋄x
    @test_broken @inferred(x^(-1)) == inv(x)
    @test_broken @inferred(x^(-2)) == inv(x)^2
    =#

    invx = @inferred inv(x)
    @test dim(invx) == n
    for i in 0:n
        n > 0 && @test d(invx, i) == inv(d(x, i))
        @test s(invx, i) == inv(s(x, i))
    end

    a = Linear(x => Int8(1), inv(x) => Int8(2))
    @test @inferred(one(a)) == Linear(one(T) => one(Int8))
    @test isone(one(a))
    @test a * one(a) == a == one(a) * a
    if iseven(n)
        @test_broken a*a*a == @inferred a^3
    else
        is_commutative && @test iszero(a*a)
    end
end

function test_groupsimplex(x::T, y::T, is_commutative) where T <: AbstractSimplex
    k, l = dim(x), dim(y)

    if k == l
        xy = @inferred x⋄y
        @test dim(xy) == k
        for i in 0:k
            k > 0 && @test d(xy, i) == d(x, i)⋄d(y, i)
            @test s(xy, i) == s(x, i)⋄s(y, i)
        end
        @test xy⋄x == x⋄y⋄x == x⋄(y⋄x)
        if is_commutative
	    @test y⋄x == xy
            @test x/y == x⋄inv(y)
	else
            @test_throws Exception x/y
        end
    else
        @test_throws Exception x⋄y
    end

    a = Linear{T,BigInt}(x => 2)
    b = Linear{T,Float32}(y => 3)
    ab = @inferred a*b
    @test coefftype(ab) == promote_type(BigInt, Float32) && termtype(ab) == T
    !iszero(ab) && @test deg(ab) == k+l

    is_commutative && @test b*a == (-1)^(k*l) * ab

    @test diff(ab) == diff(a)*b + (-1)^k*a*diff(b)
end

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
