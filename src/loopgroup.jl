#
# loop group
#

import Base: one, isone, inv, *, /, ^

# LoopGroupGenerator

struct LoopGroupGenerator{T <: AbstractSimplex}
    gen::T
    inv::Bool
end

function show(io::IO, u::LoopGroupGenerator)
    show(io, u.gen)
    if u.inv
        # print(io, "^{-1}")
        print(io, "⁻¹")
    end
end

@struct_equal_hash LoopGroupGenerator{T} where T <: AbstractSimplex

dim(u::LoopGroupGenerator) = dim(u.gen)-1

@propagate_inbounds d(u::LoopGroupGenerator, k) = LoopGroupGenerator(d(u.gen, k), u.inv)

@propagate_inbounds s(u::LoopGroupGenerator, k) = LoopGroupGenerator(s(u.gen, k), u.inv)

@propagate_inbounds isdegenerate(u::LoopGroupGenerator, k) = isdegenerate(u.gen, k)

inv(u::LoopGroupGenerator) = LoopGroupGenerator(u.gen, !u.inv)

areinverse(u::LoopGroupGenerator{T}, v::LoopGroupGenerator{T}) where T = u.gen == v.gen && u.inv != v.inv

# LoopGroupSimplex

export LoopGroupSimplex, twf_loop

"""
    LoopGroupSimplex{T<:AbstractSimplex} <: AbstractSimplex

    LoopGroupSimplex(x::T) where T <: AbstractSimplex

This type represents simplices in Kan loop groups.

The constructor returns the simplex in the loop group determined by
the simplex `x`, which must be of strictly positive dimension.

# Examples
```jldoctest
julia> using SimplicialSets: s, d

julia> x = SymbolicSimplex(:x, 2); y = LoopGroupSimplex(x)
⟨x[0,1,2]⟩

julia> d(y, 1), d(y, 0)
(⟨x[0,1]⟩, ⟨x[1,2]⁻¹,x[0,2]⟩)

julia> LoopGroupSimplex(s(x, 0))
⟨⟩
```
"""
struct LoopGroupSimplex{T<:AbstractSimplex} <: AbstractSimplex
    gens::Vector{LoopGroupGenerator{T}}
    dim::Int
end

function LoopGroupSimplex(x::T) where T <: AbstractSimplex
    n = dim(x)
    if n == 0
         error("illegal argument")
    end
    if isdegenerate(x, 0)
        LoopGroupSimplex(LoopGroupGenerator{T}[], n-1)
    else
        LoopGroupSimplex([LoopGroupGenerator(x, false)], n-1)
    end
end

function show(io::IO, g::LoopGroupSimplex)
    print(io, '⟨')
    join(io, g.gens, ',')
    print(io, '⟩')
end

dim(g::LoopGroupSimplex) = g.dim

copy(g::LoopGroupSimplex) = LoopGroupSimplex(copy(g.gens), g.dim)

length(g::LoopGroupSimplex) = length(g.gens)

iterate(g::LoopGroupSimplex, state...) = iterate(g.gens, state...)

function one(::Type{LoopGroupSimplex{T}}, n::Integer = 0) where T <: AbstractSimplex
    n >= 0 || error("dimension must be non-negative")
    LoopGroupSimplex{T}(LoopGroupGenerator{T}[], n)
    # we need the parameter in "LoopGroupSimplex{T}" to get the automatic conversion of n to Int
end

one(g::T, n = dim(g)) where T <: LoopGroupSimplex = one(T, n)

isone(g::LoopGroupSimplex) = length(g) == 0

@struct_equal_hash LoopGroupSimplex{T} where T <: AbstractSimplex

function append!(g::LoopGroupSimplex{T}, u::LoopGroupGenerator{T}, ::Val{nondeg0}) where {T, nondeg0}
    gens = g.gens
    @boundscheck if dim(g) != dim(u)
        error("illegal arguments")
    end
    @inbounds if nondeg0 || !isdegenerate(u, 0)
        if isempty(gens) || !areinverse(gens[end], u)
            push!(gens, u)
        else
            pop!(gens)
        end
    end
    g
end

function d(g::T, k::Integer) where T <: LoopGroupSimplex
    n = dim(g)
    m = n == 0 ? -1 : n
    0 <= k <= m || error("index outside the allowed range 0:$m")
    h = one(T, dim(g)-1)
    sizehint!(h.gens, (k == 0 ? 2 : 1) * length(g))
    @inbounds for u in g.gens
        v = d(u, k+1)
        if k != 0
            append!(h, v, Val(false))
        else
            w = inv(d(u, 0))
            append!(h, u.inv ? v : w, Val(false))
            append!(h, u.inv ? w : v, Val(false))
        end
    end
    h
end

function s(g::LoopGroupSimplex, k::Integer)
    n = dim(g)
    0 <= k <= n || error("index outside the allowed range 0:$n")
    LoopGroupSimplex(map(u -> @inbounds(s(u, k+1)), g.gens), n+1)
end

"""
    Base.inv(g::L) where L <: LoopGroupSimplex -> L

Return the inverse of the simplex `g` in the loop group.
"""
function inv(g::LoopGroupSimplex)
    gens = similar(g.gens)
    for (k, u) in enumerate(g.gens)
        @inbounds gens[end+1-k] = inv(u)
    end
    LoopGroupSimplex(gens, g.dim)
end

"""
    mul!(g::L, hs::L...) where L <: LoopGroupSimplex

Multiply the simplex `g` from the right by the simplices given as other arguments in-place
and return `g`.
"""
function mul!(g::LoopGroupSimplex{T}, hs::LoopGroupSimplex{T}...) where T <: AbstractSimplex
    all(==(dim(g)) ∘ dim, hs) || error("illegal arguments")
    sizehint!(g.gens, length(g)+sum(length, hs; init = 0))
    for h in hs, v in h.gens
        append!(g, v, Val(true))
    end
    g
end

"""
    *(g::L...) where L <: LoopGroupSimplex -> L

Multiply the given simplices in the loop group. At least one simplex must be given, and they
must all have the same dimension.
"""
*(g::LoopGroupSimplex{T}, hs::LoopGroupSimplex{T}...) where T <: AbstractSimplex = mul!(copy(g), hs...)

# twisting function
twf_loop(x::AbstractSimplex) = LoopGroupSimplex(x)
