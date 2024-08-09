#
# simplicial bar construction
#

export BarSimplex

import Base: *, /, ^, one, isone, inv,
    eltype, length, iterate

"""
    BarSimplex{T} <: AbstractSimplex

A type representing simplices in a simplicial bar construction. If `T`
is a subtype of `AbstractSimplex`, then it is assumed to be a simplicial group;
otherwise `T` is assumed to be a discrete group.

Iterating over a `BarSimplex` means iterating over its components.
"""
struct BarSimplex{T} <: AbstractSimplex
    g::Vector{T}
end

function BarSimplex(g::Vector{T}) where T <: AbstractSimplex
    @boundscheck all(enumerate(g)) do (i, x)
        dim(x) == i-1
    end || error("simplices have incorrect dimensions", g)
    BarSimplex{T}(g)
end

@propagate_inbounds function BarSimplex(iter; op::Union{typeof(*),typeof(+)} = *)
    if op isa typeof(+)
        BarSimplex(AddToMul.(iter))
        # todo: this leads to additional allocations
    else
        BarSimplex(collect(iter))
    end
end

show(io::IO, x::BarSimplex) = print(io, '[', join(x.g, ','), ']')

@struct_equal_hash BarSimplex{T} where T

copy(x::BarSimplex{T}) where T = BarSimplex{T}(copy(x.g))

eltype(::Type{BarSimplex{T}}) where T = T

"""
    length(x::BarSimplex) -> Int

The length of a `BarSimplex` is the number of its components.
"""
length(x::BarSimplex) = length(x.g)

iterate(x::BarSimplex, state...) = iterate(x.g, state...)

dim(x::BarSimplex) = length(x)

# note: multiplication of bar simplices only makes sense for commutative groups

"""
    one(x::BarSimplex{T}, n::Integer = dim(x)) where T -> BarSimplex{T}

Return the identity element in the group of `n`-simplices
in the simplicial bar construction containing `x`.
Here `T` is assumed to be a commutative (simplicial) group.
"""
one(x::BarSimplex, n::Integer = dim(x)) = one(typeof(x), n)

isone(x::BarSimplex) = all(isone, x.g)

"""
    inv(x::BarSimplex{T}) where T -> BarSimplex{T}

Return the inverse element of `x` in the group of `n`-simplices in the simplicial
bar construction containing that simplex.
Here `T` is assumed to be a commutative (simplicial) group.
"""
inv(x::BarSimplex{T}) where T = BarSimplex(inv.(x.g))

"""
    *(x::BarSimplex{T}...) where T -> BarSimplex{T}

Return the product of the given simplices, which must all have the same dimension.
Here `T` is assumed to be a commutative (simplicial) group.
"""
function *(x::BarSimplex{T}, xs::BarSimplex{T}...) where T
    xs = (x, xs...)
    @boundscheck allequal(map(dim, xs)) || error("illegal arguments")
    BarSimplex(map(*, map(x -> x.g, xs)...))
end

"""
    ^(x::BarSimplex{T}, n::Integer) where T -> BarSimplex{T}

Return the `n`-th power of the simplex `x`.
Here `T` is assumed to be a commutative (simplicial) group.
"""
function ^(x::BarSimplex, n::Integer)
    BarSimplex(x.g .^ n)
end

"""
    /(x::BarSimplex{T}, y::BarSimplex{T}) where T -> BarSimplex{T}

Return the quotient of the `x` by `y` in the commutative group of `n`-simplices
in the simplicial bar construction, where `n` is the common dimension of `x` and `y`.
Here `T` is assumed to be a commutative (simplicial) group.
"""
function /(x::BarSimplex{T}, y::BarSimplex{T}) where T
    @boundscheck dim(x) == dim(y) || error("illegal arguments")
    BarSimplex(x.g ./ y.g)
    # BarSimplex(map(splat(/), zip(x.g, y.g)))
end

#
# bar construction for discrete groups
#

function d(x::BarSimplex{T}, k::Integer) where T
    @boundscheck if k < 0 || k > (n = dim(x)) || n == 0
        error("illegal arguments")
    end
    g = x.g
    @inbounds gg = T[if i < k
            g[i]
        elseif i == k
            g[i]*g[i+1]
        else
            g[i+1]
        end
        for i in 1:length(g)-1]
    BarSimplex(gg)
end

function s(x::BarSimplex{T}, k::Integer) where T
    @boundscheck if k < 0 || k > dim(x)
        error("illegal arguments")
    end
    g = x.g
    k += 1
    @inbounds gg = T[if i < k
            g[i]
        elseif i == k
            one(T)
        else
            g[i-1]
        end
        for i in 1:length(g)+1]
    BarSimplex(gg)
end

@inline function isdegenerate(x::BarSimplex, k::Integer)
    @boundscheck if k < 0 || k >= dim(x)
        error("illegal arguments")
    end
    @inbounds isone(x.g[k+1])
end

function one(::Type{BarSimplex{T}}, n::Integer = 0) where T
    n >= 0 || error("dimension must be non-negative")
    BarSimplex(T[one(T) for k in 0:n-1])
end

#
# bar construction for simplicial groups
#

export twf_bar

function d(x::BarSimplex{T}, k::Integer) where T <: AbstractSimplex
    @boundscheck if k < 0 || k > (n = dim(x)) || n == 0
        error("illegal arguments")
    end
    g = x.g
    @inbounds gg = T[if i < k
            g[i]
        elseif i == k
            g[i]*d(g[i+1], 0)
        else
            d(g[i+1], i-k)
        end
        for i in 1:length(g)-1]
    @inbounds BarSimplex(gg)
end

function s(x::BarSimplex{T}, k::Integer) where T <: AbstractSimplex
    @boundscheck if k < 0 || k > dim(x)
        error("illegal arguments")
    end
    g = x.g
    k += 1
    @inbounds gg = [if i < k
            g[i]
        elseif i == k
            one(T, k-1)
        else
            s(g[i-1], i-k-1)
        end
        for i in 1:length(g)+1]
    @inbounds BarSimplex(gg)
end

@inline function isdegenerate(x::BarSimplex{T}, k::Integer) where T <: AbstractSimplex
    @boundscheck if k < 0 || k >= dim(x)
        error("illegal arguments")
    end
    g = x.g
    @inbounds isone(g[k+1]) && all(i -> isdegenerate(g[i], i-k-2), k+2:dim(x))
end

function one(::Type{BarSimplex{T}}, n::Integer = 0) where T <: AbstractSimplex
    n >= 0 || error("dimension must be non-negative")
    BarSimplex(T[one(T, k) for k in 0:n-1])
end

# twisting function
function twf_bar(x::BarSimplex{T}) where T <: AbstractSimplex
    if deg(x) == 0
        error("twisting function not defined for 0-simplices")
    else
        return x.g[end]
    end
end
