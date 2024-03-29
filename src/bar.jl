#
# simplicial bar construction
#

export BarSimplex

import Base: *, /, ^, one, isone, inv

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

length(x::BarSimplex) = length(x.g)

iterate(x::BarSimplex, state...) = iterate(x.g, state...)

dim(x::BarSimplex) = length(x)

# note: multiplication of bar simplices only makes sense for commutative groups

one(x::BarSimplex, n::Integer = dim(x)) = one(typeof(x), n)

isone(x::BarSimplex) = all(isone, x.g)

inv(x::BarSimplex{T}) where T = BarSimplex(inv.(x.g))

function *(x::BarSimplex{T}, ys::BarSimplex{T}...) where T
    @boundscheck all(==(dim(x)) ∘ dim, ys) || error("illegal arguments")
    BarSimplex(.*(x.g, map(y -> y.g, ys)...))
end

function ^(x::BarSimplex, n::Integer)
    BarSimplex(x.g .^ n)
end

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
