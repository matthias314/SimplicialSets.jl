#
# Interval
#

const Interval = Union{Tuple{Integer,Integer}, Pair{<:Integer,<:Integer}, UnitRange{<:Integer}}

interval_length(k::Interval) = last(k)-first(k)+1

#
# AbstractSimplex datatype
#

export AbstractSimplex, isdegenerate

abstract type AbstractSimplex end

@linear_broadcastable AbstractSimplex

# for a new concrete subtype NewSimplex of AbstractSimplex, at least
# the following methods must be defined:
#     dim(::NewSimplex)
#     d(:NewSimplex, ::Int)
#     s(:NewSimplex, ::Int)

"""
    dim(x::AbstractSimplex) -> Int

Return the dimension of the simplex `x`.
"""
function dim end

function d!(x::AbstractSimplex, kk::AbstractVector{<:Integer})
    for k in reverse(kk)
        d!(x, k)
    end
    x
end

"""
    d(x::T, k) where T <: AbstractSimplex -> T
    d(x::T, ks::Integer...) where T <: AbstractSimplex -> T

In the first form, return the `k`-th facet of `x`. In the second form,
apply `d` repeatedly to `x`, starting with the last element of `ks`.

TODO CHECK
```
"""
@generated function d(x::T, kk) where T <: AbstractSimplex
    if hasmethod(d!, (T, Int))
        :(d!(copy(x), kk))
    else
        :(foldl(d, reverse(kk), init = x))
    end
end

function s!(x::AbstractSimplex, kk::AbstractVector{<:Integer})
    for k in kk
        s!(x, k)
    end
    x
end

@generated function s(x::T, kk) where T <: AbstractSimplex
    if hasmethod(s!, (T, Int))
        :(s!(copy(x), kk))
    else
        quote
            for k in kk
                x = s(x, k)
            end
            x
        end
    end
end

#=
function r!(x::AbstractSimplex, kk)
    l = dim(x)
    for k in reverse(kk)
        if !(0 <= k <= l)
            error("index outside the allowed range 0:$l")
        end
        d!(x, k+1:l)
        l = k-1
    end
    d!(x, 0:l)
end
=#

function r(x::T, kk) where T <: AbstractSimplex
# kk is a tuple/vector/iterator of Tuple{Integer,Integer}
    isempty(kk) && error("at least one interval must be given")

    kk[end] isa Interval || error("second argument must contain integer tuples, pairs or unit ranges")
    y = x
    for k in dim(x):-1:last(kk[end])+1
        y = d(y, k)
    end
    for i in length(kk)-1:-1:1
        kk[i] isa Interval || error("second argument must contain integer tuples, pairs or unit ranges")
        kki1 = first(kk[i+1])
        kki = last(kk[i])
        if kki == kki1
            y = s(y, kki)
        else
            for k in kki1-1:-1:kki+1
                y = d(y, k)
            end
        end
    end
    for k in first(kk[1])-1:-1:0
        y = d(y, k)
    end
    y
end

"""
    isdegenerate(x::AbstractSimplex, k::Integer) -> Bool

Return `true` if `x` is degenerate at position `k` and `false` otherwise.
The index `k` must be between `0` and `dim(x)-1`.

A simplex `x` of positive dimension is degenerate at position `k` if `x == s(d(x, k), k)`.
"""
function isdegenerate(x::AbstractSimplex, k::Integer)
    @boundscheck if k < 0 || k >= dim(x)
        error("index outside the allowed range 0:$(dim(x)-1)")
    end
    @inbounds x == s(d(x, k), k)
end

"""
    isdegenerate(x::AbstractSimplex, k0::Integer, k1::Integer) -> Bool

Return `true` if `x` is degenerate on the interval `k0:k1` and `false` otherwise.
The indices must satisfy `0 <= k0 <= k1 <= dim(x)`.

A simplex is degenerate on the interval `k0:k1` if it degenerate at some position `k`
between `k0` and `k1-1`.
"""
@propagate_inbounds isdegenerate(x::AbstractSimplex, k0::Integer, k1::Integer) = any(k -> isdegenerate(x, k), k0:k1-1)

"""
    isdegenerate(x::AbstractSimplex) -> Bool

Return `true` if `x` is degenerate and `false` otherwise.
"""
isdegenerate(x::AbstractSimplex) = @inbounds isdegenerate(x, 0, dim(x))

"""
    LinearCombinations.linear_filter(x::AbstractSimplex) -> Bool

Return `true` if `x` is non-degenerate and `false` otherwise.
The effect of this is that linear combinations of simplices represent
elements of the *normalized* chain complex of the corresponding simplicial set.

See also `LinearCombinations.linear_filter`.

# Examples
```jldoctest
julia> using LinearCombinations; using SimplicialSets: s

julia> x = SymbolicSimplex(:x, 2)
x[0,1,2]

julia> Linear(s(x, 1) => 1)
0
```
"""
linear_filter(x::AbstractSimplex) = !isdegenerate(x)

"""
    deg(x::AbstractSimplex) -> Int

Return the degree of `x`, which is defined as its dimension.

See also [`dim`](@ref).
"""
deg(x::AbstractSimplex) = dim(x)

"""
    diff(x::T) where T <: AbstractSimplex -> Linear{T}

Return the differential or boundary of `x` as a linear combination.
By default, the coefficients are of type `Int`.

This functions supports the keyword arguments `coefftype`, `addto`,
`coeff` and `is_filtered` as described for `@linear`.

TODO

See also `LinearCombinations.@linear`.

# Examples
```jldoctest
julia> using LinearCombinations; using LinearCombinations: diff

julia> x = SymbolicSimplex(:x, 2)
x[0,1,2]

julia> a = diff(x)
x[1,2]-x[0,2]+x[0,1]

julia> b = zero(a); diff(x; addto = b, coeff = 2)
2*x[1,2]-2*x[0,2]+2*x[0,1]

julia> b
2*x[1,2]-2*x[0,2]+2*x[0,1]
```
"""
@linear_kw function diff(x::T;
        coefftype = Int,
        addto = zero(Linear{T,unval(coefftype)}),
        coeff = ONE,
        is_filtered = false) where T <: AbstractSimplex
    if iszero(coeff)
        return addto
    end
    n = dim(x)
    if n != 0
        for k in 0:n
            addmul!(addto, d(x, k), signed(k, coeff))
        end
    end
    addto
end

function ^(g::AbstractSimplex, n::Integer)
    if n >= 32
        # square-and-multiply
        s = g
        m = one(g)
        while n != 0
            if isodd(n)
                m *= s
            end
            s *= s
            n >>= 1
        end
        m
    elseif n > 0
        *(ntuple(Returns(g), n)...)
    elseif n == 0
        one(g)
    else
        inv(g)^(-n)
    end
end
