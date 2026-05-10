#
# AbstractProductSimplex
#

abstract type AbstractProductSimplex{T<:Tuple} <: AbstractSimplex end

"""
    Tuple(x::AbstractProductSimplex{T}) where T -> T <: Tuple{Vararg{AbstractSimplex}}
    components(x::AbstractProductSimplex{T}) where T -> T <: Tuple{Vararg{AbstractSimplex}}

Return the tuple of component simplices of `x`.

!!! note
    The function `components` is deprecated. Use `Tuple` instead.
"""
Tuple(x::AbstractProductSimplex), components

"""
    fieldtypes(::Type{P}) where P <: AbstractProductSimplex -> Tuple

Return the types of the components of `P` as a tuple.

# Example
```jldoctest
julia> fieldtypes(ProductSimplex{Tuple{SymbolicSimplex{Symbol},}})
(Char, String)
```
"""
Base.fieldtypes(::Type{<:AbstractProductSimplex{T}}) where T <: Tuple = fieldtypes(T)

"""
    length(x::AbstractProductSimplex) -> Int

Return the number of components (or factors) of `x`.
"""
Base.length(x::AbstractProductSimplex) = length(Tuple(x))

Base.firstindex(x::AbstractProductSimplex) = 1
Base.lastindex(x::AbstractProductSimplex) = length(x)

Base.iterate(x::AbstractProductSimplex, state...) = iterate(Tuple(x), state...)

@propagate_inbounds Base.getindex(x::AbstractProductSimplex, k) = Tuple(x)[k]

#
# ProductSimplex datatype
#

export ProductSimplex

using Base: @__MODULE__ as @MODULE

# struct ProductSimplex{T<:Tuple{Vararg{AbstractSimplex}}} <: AbstractSimplex
"""
    ProductSimplex{T<:Tuple{Vararg{AbstractSimplex}}} <: AbstractSimplex

    ProductSimplex{T}{t::Tuple{Vararg{AbstractSimplex}} [; dim::Integer]}

    ProductSimplex(t::Tuple{Vararg{AbstractSimplex}} [; dim::Integer])
    ProductSimplex(xs::AbstractSimplex... [; dim::Integer])

A type representing an element in the product of simplicial sets. Empty products are allowed.
The component simplices must all be of the same dimension.
They may be given as a tuple or as individual arguments.

In the case of the empty product, the keyword argument `dim` is required to determine the
dimension of the resulting simplex. Otherwise `dim` is optional, but if present, it must be correct.

See also [`Tuple(x::ProductSimplex)`](@ref).

# Examples
```jldoctest
julia> x, y = SymbolicSimplex(:x, 2), SymbolicSimplex(:y, 2)
(x[0,1,2], y[0,1,2])

julia> z = ProductSimplex(x, y)
(x[0,1,2],y[0,1,2])

julia> w = ProductSimplex(dim = 2)
()

julia> dim(w)
2

julia> ProductSimplex(x, y; dim = 1)
ERROR: dimensions of simplices do not match
[...]
```
"""
struct ProductSimplex{T<:Tuple} <: AbstractProductSimplex{T}
    xl::T
    dim::Int

    # we need `@propagate_inbounds` instead of `@inline`, see julia/#30411
    @propagate_inbounds ProductSimplex{T}(xl::T; dim::Union{Integer,Missing} = missing) where T<:Tuple{Vararg{AbstractSimplex}} = begin
        if dim === missing
            if isempty(xl)
                error("use 'dim' to specify the dimension of an empty product simplex")
            else
                dim = (@MODULE).dim(xl[1])
            end
        end
        @boundscheck begin
            dim >= 0 || error("dimension must be non-negative")
            all(==(dim) ∘ (@MODULE).dim, xl) || error("dimensions of simplices do not match")
        end
        new{T}(xl, dim)
    end
end

@propagate_inbounds ProductSimplex(xl::T; dim::Union{Integer,Missing} = missing) where T <: Tuple  =
    ProductSimplex{T}(xl; dim)

@propagate_inbounds ProductSimplex(x::AbstractSimplex...; kw...) = ProductSimplex(x; kw...)

function show(io::IO, x::ProductSimplex)
    print(io, '(', join(map(repr, Tuple(x)), ','), ')')
end

Base.Tuple(x::AbstractProductSimplex) = x.xl
@deprecate components(x::AbstractProductSimplex) Tuple(x)

# copy(x::ProductSimplex) = ProductSimplex(copy(x.xl))
copy(x::ProductSimplex) = x

Base.convert(::Type{P}, x::ProductSimplex) where P <: ProductSimplex = @inbounds P(Tuple(x); dim = dim(x))

@struct_equal_hash ProductSimplex{T} where T
# @struct_equal_hash ProductSimplex
# TODO: should we take the tuple type T into account?

dim(x::ProductSimplex) = x.dim

@propagate_inbounds function d(x::ProductSimplex, k::Integer)
    n = dim(x)
    @boundscheck begin
        m = n == 0 ? -1 : n
        0 <= k <= m || error("index outside the allowed range 0:$m")
    end
    ProductSimplex(map(y -> @inbounds(d(y, k)), x.xl); dim = n-1)
end

@propagate_inbounds function s(x::ProductSimplex, k::Integer)
    n = dim(x)
    @boundscheck begin
        0 <= k <= n || error("index outside the allowed range 0:$n")
    end
    ProductSimplex(map(y -> @inbounds(s(y, k)), x.xl); dim = n+1)
end

@propagate_inbounds function r(x::ProductSimplex, kk)
    ProductSimplex(map(y -> r(y, kk), x.xl))
end

@propagate_inbounds function r(x::ProductSimplex{Tuple{}}, kk)
    isempty(kk) && error("at least one interval must be given")
    ProductSimplex((); dim = sum(map(interval_length, kk))-1)
end

@inline function isdegenerate(x::ProductSimplex, k::Integer)
    @boundscheck if k < 0 || k >= dim(x)
        error("index outside the allowed range 0:$(dim(x))")
    end
    @inbounds all(y -> isdegenerate(y, k), Tuple(x))
end

# group operations

function Base.one(::Type{ProductSimplex{T}}, n::Integer = 0) where T
    ProductSimplex(map(Fix2(one, n), fieldtypes(T)); dim = n)
end

Base.one(g::T, n::Integer = dim(g)) where T <: ProductSimplex = one(T, n)

function ⋅(gs::Vararg{ProductSimplex{<:NTuple{N,AbstractSimplex}},M}) where {N,M}
    ProductSimplex(map(⋅, map(Tuple, gs)...); dim = dim(gs[1]))
end

Base.inv(g::ProductSimplex) = ProductSimplex(map(inv, Tuple(g)); dim = dim(g))

# concatenating and flattening ProductSimplex

using LinearCombinations: tuple_cat
import LinearCombinations: cat, flatten, tuple_flatten

"""
    SimplicialSets.cat(x::AbstractProductSimplex...) -> ProductSimplex

Return the product simplex that is the concatenation of the simplices given as arguments.

This function is linear. Also note that it is overloaded from the package `LinearCombinations`,
not from `Base`. If one wants to use the short form `cat`, then one needs to import the function
via `using` or `import`.

See also [`flatten`](@ref).

# Example
```jldoctest
julia> using SimplicialSets: cat   # or: using LinearCombinations: cat

julia> u = ProductSimplex(SymbolicSimplex(:x, 2), SymbolicSimplex(:y, 2))
(x[0,1,2],y[0,1,2])

julia> v = ProductSimplex(SymbolicSimplex(:z, 2), SymbolicSimplex(:w, 2))
(z[0,1,2],w[0,1,2])

julia> cat(u, v)
(x[0,1,2],y[0,1,2],z[0,1,2],w[0,1,2])
```
"""
cat(x::AbstractProductSimplex...) = ProductSimplex(tuple_cat(x...); dim = dim(x[1]))

tuple_flatten(x::AbstractProductSimplex) = tuple_cat(map(tuple_flatten, Tuple(x))...)

"""
    SimplicialSets.flatten(x::AbstractProductSimplex) -> ProductSimplex

Return the product simplex that is obtained by recursively flattening all product simplices
appearing within `x`.

This function is linear. Also note that it is overloaded from the package `LinearCombinations`.

See also [`LinearCombinations.Regroup`](@ref), [`SimplicialSets.cat`](@ref).

# Examples
```jldoctest
julia> using SimplicialSets: flatten   # or: using LinearCombinations: flatten

julia> u = ProductSimplex(SymbolicSimplex(:x, 2), SymbolicSimplex(:y, 2))
(x[0,1,2],y[0,1,2])

julia> v = ProductSimplex(SymbolicSimplex(:z, 2), SymbolicSimplex(:w, 2))
(z[0,1,2],w[0,1,2])

julia> flatten(ProductSimplex(u, v))
(x[0,1,2],y[0,1,2],z[0,1,2],w[0,1,2])

julia> flatten(ProductSimplex(ProductSimplex(u, v), u))
(x[0,1,2],y[0,1,2],z[0,1,2],w[0,1,2],x[0,1,2],y[0,1,2])
```
"""
flatten(x::AbstractProductSimplex) = ProductSimplex(tuple_flatten(x); dim = dim(x))

#
# regrouping
#

using LinearCombinations: regroup_check_arg, regroup_eval_expr, regroup_getindex

"""
    swap(z::AbstractProductSimplex{Tuple{S,T}}) where {S <: AbstractSimplex, T <: AbstractSimplex} -> ProductSimplex{Tuple{T,S}}

Swap the two components of the `AbstractProductSimplex` `z` and return the resulting `ProductSimplex`.

This function is linear. Also note that it is overloaded from the package `LinearCombinations`.

See also [`LinearCombinations.Regroup`](@ref).
"""
swap(::AbstractProductSimplex{Tuple{S,T}}) where {S <: AbstractSimplex, T <: AbstractSimplex}

"""
    (rg::LinearCombinations.Regroup)(z::AbstractProductSimplex) -> ProductSimplex

Apply the `Regroup` object `rg` to `z` and return the result. This allows to permute and restructure
the components of a product simplex in an arbitrary way (without dropping any component).

This functions is linear and supports the keyword arguments `coefftype`, `addto`,
`coeff` and `is_filtered` as described for `@linear`.

See `LinearCombinations.@linear`, `LinearCombinations.@regroup_str`, [`swap`](@ref), [`flatten`](@ref).

# Example

```@jldoctest
julia> using LinearCombinations

julia> rg = regroup"((1, 2), 3) -> (2, (3, 1))"
Regroup{((1, 2), 3), (2, (3, 1))}

julia> x, y, z = SymbolicSimplex(:x, 2), SymbolicSimplex(:y, 2), SymbolicSimplex(:z, 2)
(x[0,1,2], y[0,1,2], z[0,1,2])

julia> w = ProductSimplex(ProductSimplex(x, y), z)
((x[0,1,2],y[0,1,2]),z[0,1,2])

julia> rg(w)
(y[0,1,2],(z[0,1,2],x[0,1,2]))
```
"""
function (rg::Regroup{A})(x::T) where {A,T<:AbstractProductSimplex}
    regroup_check_arg(AbstractProductSimplex, typeof(A), T) ||
        error("argument type $(typeof(x)) does not match first Regroup parameter $A")
    @inbounds regroup_eval_expr(rg, regroup_getindex, ProductSimplex, x)
end
