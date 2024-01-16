#
# ProductSimplex datatype
#

export ProductSimplex, components

using Base: @__MODULE__ as @MODULE

import Base: length, iterate, convert

# struct ProductSimplex{T<:Tuple{Vararg{AbstractSimplex}}}
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
    print(io, '(', join(map(repr, components(x)), ','), ')')
end

components(x::ProductSimplex) = x.xl

length(x::ProductSimplex) = length(components(x))

firstindex(x::ProductSimplex) = 1
lastindex(x::ProductSimplex) = length(x)

iterate(x::ProductSimplex, state...) = iterate(components(x), state...)

@propagate_inbounds getindex(x::ProductSimplex, k) = components(x)[k]

# copy(x::ProductSimplex) = ProductSimplex(copy(x.xl))
copy(x::ProductSimplex) = x

convert(::Type{P}, x::ProductSimplex) where P <: ProductSimplex = @inbounds P(components(x); dim = dim(x))

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
    @inbounds all(y -> isdegenerate(y, k), components(x))
end

# concatenating and flattening ProductSimplex

using LinearCombinations: _cat
import LinearCombinations: cat, flatten, _flatten

cat(x::ProductSimplex...) = ProductSimplex(_cat(x...); dim = dim(x[1]))

_flatten(x::ProductSimplex) = _cat(map(_flatten, components(x))...)

flatten(x::ProductSimplex) = ProductSimplex(_flatten(x); dim = dim(x))

#
# regrouping
#

using LinearCombinations: regroup_check_arg, regroup_eval_expr
import LinearCombinations: _length, _getindex

_length(::Type{<:ProductSimplex{T}}) where T <: Tuple = _length(T)

@propagate_inbounds _getindex(::Type{T}, i) where T <: ProductSimplex = _getindex(T.parameters[1], i)

function (rg::Regroup{A})(x::T) where {A,T<:ProductSimplex}
    regroup_check_arg(ProductSimplex, typeof(A), T) ||
        error("argument type $(typeof(x)) does not match first Regroup parameter $A")
    @inbounds regroup_eval_expr(rg, _getindex, ProductSimplex, x)
end
