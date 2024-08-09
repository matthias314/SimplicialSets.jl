#
# ProductSimplex datatype
#

export ProductSimplex

using Base: @__MODULE__ as @MODULE

import Base: length, iterate, convert

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
struct ProductSimplex{T<:Tuple} <: AbstractSimplex
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

"""
    Tuple(x::ProductSimplex{T}) where T <: Tuple{Vararg{AbstractSimplex}} -> T

Return the tuple of component simplices of `x`.
"""
Base.Tuple(x::ProductSimplex) = x.xl

"""
    length(x::ProductSimplex) -> Int

Return the number of components (or factors) of `x`.
"""
length(x::ProductSimplex) = length(Tuple(x))

firstindex(x::ProductSimplex) = 1
lastindex(x::ProductSimplex) = length(x)

iterate(x::ProductSimplex, state...) = iterate(Tuple(x), state...)

@propagate_inbounds getindex(x::ProductSimplex, k) = Tuple(x)[k]

# copy(x::ProductSimplex) = ProductSimplex(copy(x.xl))
copy(x::ProductSimplex) = x

convert(::Type{P}, x::ProductSimplex) where P <: ProductSimplex = @inbounds P(Tuple(x); dim = dim(x))

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

# concatenating and flattening ProductSimplex

using LinearCombinations: _cat
import LinearCombinations: cat, flatten, _flatten

"""
    SimplicialSets.cat(x::ProductSimplex...) -> ProductSimplex

Return the product simplex that is the concatenation of the simplices given as arguments.

This function is linear. Also note that it is overloaded from the package `LinearCombinations`,
not from `Base`. If one wants to use the short form `cat`, then one needs to import the function
via `using` or `import`.

See also [`flatten`](@ref).

# Example
```jldoctest
julia> using SimplicialSets: cat   # or: using LinearCombinations: cat

julia> u = ProductSimplex(SymbolicSimplex(:x, 2), SymbolicSimplex(:y, 2))
(x[0,1],y[0,1])

julia> v = ProductSimplex(SymbolicSimplex(:z, 2), SymbolicSimplex(:w, 2))
(z[0,1],w[0,1])

julia> cat(u, v)
(x[0,1],y[0,1],z[0,1],w[0,1])
```
"""
cat(x::ProductSimplex...) = ProductSimplex(_cat(x...); dim = dim(x[1]))

_flatten(x::ProductSimplex) = _cat(map(_flatten, Tuple(x))...)

"""
    SimplicialSets.flatten(x::ProductSimplex) -> ProductSimplex

Return the product simplex that is obtained by recursively flattening all product simplices
appearing within `x`.

This function is linear. Also note that it is overloaded from the package `LinearCombinations`.

See also [`LinearCombinations.Regroup`](@ref), [`SimplicialSets.cat`](@ref).

# Examples
```jldoctest
julia> using SimplicialSets: flatten   # or: using LinearCombinations: flatten

julia> u = ProductSimplex(SymbolicSimplex(:x, 2), SymbolicSimplex(:y, 2))
(x[0,1],y[0,1])

julia> v = ProductSimplex(SymbolicSimplex(:z, 2), SymbolicSimplex(:w, 2))
(z[0,1],w[0,1])

julia> flatten(ProductSimplex(u, v))
(x[0,1],y[0,1],z[0,1],w[0,1])

julia> flatten(ProductSimplex(ProductSimplex(u, v), u))
(x[0,1],y[0,1],z[0,1],w[0,1],x[0,1],y[0,1])
```
"""
flatten(x::ProductSimplex) = ProductSimplex(_flatten(x); dim = dim(x))

#
# regrouping
#

using LinearCombinations: regroup_check_arg, regroup_eval_expr
import LinearCombinations: _length, _getindex

_length(::Type{<:ProductSimplex{T}}) where T <: Tuple = _length(T)

@propagate_inbounds _getindex(::Type{T}, i) where T <: ProductSimplex = _getindex(T.parameters[1], i)

"""
    swap(z::ProductSimplex{Tuple{S,T}}) where {S <: AbstractSimplex, T <: AbstractSimplex} -> ProductSimplex{Tuple{T,S}}

Swap the two components of the `ProductSimplex` `z` and return the resulting `ProductSimplex`.

This function is linear. Also note that it is overloaded from the package `LinearCombinations`.

See also [`LinearCombinations.Regroup`](@ref).
"""
swap(::ProductSimplex{Tuple{S,T}}) where {S <: AbstractSimplex, T <: AbstractSimplex}

"""
    (rg::LinearCombinations.Regroup)(z::ProductSimplex) -> ProductSimplex

Apply the `Regroup` element `rg` to `z` and return the result. This allows to permute and restructure
the components of a product simplex in an arbitrary way (without dropping any component).

This functions is linear and supports the keyword arguments `coefftype`, `addto`,
`coeff` and `is_filtered` as described for `@linear`.

See `LinearCombinations.@linear`, `LinearCombinations.regroup`, [`swap`](@ref), [`flatten`](@ref).

# Example

```@jldoctest
julia> using LinearCombinations

julia> rg = regroup(:( ((1, 2), 3) ), :( (2, (3, 1))  ))
Regroup{((1, 2), 3),(2, (3, 1))}

julia> x, y, z = SymbolicSimplex(:x, 2), SymbolicSimplex(:y, 2), SymbolicSimplex(:z, 2)
(x[0,1,2], y[0,1,2], z[0,1,2])

julia> w = ProductSimplex(ProductSimplex(x, y), z)
((x[0,1,2],y[0,1,2]),z[0,1,2])

julia> rg(w)
(y[0,1,2],(z[0,1,2],x[0,1,2]))
```
"""
function (rg::Regroup{A})(x::T) where {A,T<:ProductSimplex}
    regroup_check_arg(ProductSimplex, typeof(A), T) ||
        error("argument type $(typeof(x)) does not match first Regroup parameter $A")
    @inbounds regroup_eval_expr(rg, _getindex, ProductSimplex, x)
end
