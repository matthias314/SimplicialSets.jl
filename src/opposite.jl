#
# OppositeSimplex datatype
#

export OppositeSimplex, opposite

"""
    OppositeSimplex{T<:AbstractSimplex} <: AbstractSimplex

`OppositeSimplex(x)` is the 'opposite' simplex to `x` in the sense
that `i`-th face operator corresponds to the `(n-i)`-th face operator
for `x`, and likewise for degeneracy operators.

The functions `*`, `\\`, `inv` and `one` with argument(s) of type
`OppositeSimplex` work on the underlying simplices.

Note that the linear extension of `OppositeSimplex` would not be a chain map
from the chain complex of the simplicial set containing the terms `x`
to the opposite simplicial set. For this purpose there is the function `opposite`.

See also [`opposite`](@ref).

# Examples
```jldoctest
julia> using SimplicialSets: d, s

julia> x = SymbolicSimplex(:x, 3)
x[0,1,2,3]

julia> y = OppositeSimplex(x)
OppositeSimplex(x[0,1,2,3])

julia> d(y, 1)
OppositeSimplex(x[0,1,3])

julia> s(y, 1)
OppositeSimplex(x[0,1,2,2,3])

julia> z = OppositeSimplex(LoopGroupSimplex(x))
OppositeSimplex(⟨x[0,1,2,3]⟩)

julia> inv(z)
OppositeSimplex(⟨x[0,1,2,3]⁻¹⟩)
```
"""
struct OppositeSimplex{T<:AbstractSimplex} <: AbstractSimplex
    x::T
end

function show(io::IO, y::OppositeSimplex)
    print(io, "OppositeSimplex(", repr(y.x), ')')
end

@struct_equal_hash OppositeSimplex{T} where T

copy(y::OppositeSimplex) = OppositeSimplex(copy(y.x))

dim(y::OppositeSimplex) = dim(y.x)
d(y::OppositeSimplex, k::Integer) = OppositeSimplex(d(y.x, dim(y)-k))
s(y::OppositeSimplex, k::Integer) = OppositeSimplex(s(y.x, dim(y)-k))

*(ys::OppositeSimplex...) = OppositeSimplex(prod(map(y -> y.x, ys)))
/(y::OppositeSimplex, z::OppositeSimplex) = OppositeSimplex(y.x/z.x)
inv(y::OppositeSimplex) = OppositeSimplex(inv(y.x))
one(::Type{OppositeSimplex{T}}, n...) where T = OppositeSimplex(one(T, n...))
one(y::OppositeSimplex, n...) = OppositeSimplex(one(y.x, n...))

"""
    opposite(x::T) where T <: AbstractSimplex -> OppositeSimplex
    opposite(x::T) where {S, T <: OppositeSimplex{S}} -> S
    opposite(x::T) where T <: ProductSimplex -> ProductSimplex

Return an "interpreted version" of the opposite simplex of `x`:

* If `x` is already an `OppositeSimplex`, return the underlying simplex.
* If `x` is a `ProductSimplex`, apply `opposite` to the components and
  return the corresponding `ProductSimplex`.
* In all other cases, return `OppositeSimplex(x)`.

See also
[`OppositeSimplex`](@ref),
[`opposite(::AbstractTensor)`](@ref),
[`opposite(::AbstractLinear)`](@ref).
"""
opposite(x::AbstractSimplex) = OppositeSimplex(x)
opposite(y::OppositeSimplex) = y.x
opposite(x::ProductSimplex) = ProductSimplex(map(opposite, Tuple(x)))

"""
    opposite(t::AbstractTensor) -> Tensor

Apply `opposite` to the components of `t` and return their tensor product.

See also
[`opposite(::AbstractSimplex)`](@ref),
[`opposite(::AbstractLinear)`](@ref).
"""
opposite(x::AbstractTensor) = Tensor(map(opposite, Tuple(x)))

q2(x::AbstractSimplex) = ifelse((dim(x)+1) & 2 == 0, 0, 1)
# this is 0 if dim(x) == 0 or 3 mod 4, and 1 if dim(x) == 1 or 2 mod 4

q2(t::AbstractTensor) = sum0(map(q2, Tuple(t)))

"""
    opposite(a::AbstractLinear{T,R}) -> Linear

Apply `opposite` to the terms in `a` and return a linear combination of them.
If the degree of a term `x` is congruent to `1` or `2` mod `4`, then it is
transformed to `-opposite(x)` and otherwise to `opposite(x)`.

Note that this function is **not** the linear extension of `opposite(x::AbstractSimplex)` or `OppositeSimplex`.
The signs are chosen such that this function is a chain map from the chain complex
of the simplicial set containing the terms `x` to the opposite simplicial set.

See also
[`opposite(::AbstractSimplex)`](@ref),
[`opposite(::AbstractTensor)`](@ref).

# Examples
```jldoctest
julia> using LinearCombinations; using LinearCombinations: diff

julia> a = Linear(SymbolicSimplex(:x, 1) => 2)
Linear{SymbolicSimplex{Symbol}, Int64} with 1 term:
2*x[0,1]

julia> b = opposite(a)
Linear{OppositeSimplex{SymbolicSimplex{Symbol}}, Int64} with 1 term:
-2*OppositeSimplex(x[0,1])

julia> diff(b) == opposite(diff(a))
true

julia> c = Linear(OppositeSimplex(y) => c for (y, c) in a)
Linear{OppositeSimplex{SymbolicSimplex{Symbol}}, Int64} with 1 term:
2*OppositeSimplex(x[0,1])

julia> diff(c) == opposite(diff(a))
false
```
"""
opposite(::AbstractLinear)

@linear_kw function opposite(a::AbstractLinear{T,R};
        coefftype = R,
        addto = zero(Linear{return_type(opposite, T),unval(coefftype)}),
        coeff = ONE) where {T,R}
    iszero(coeff) && return addto
    for (x, c) in a
        addmul!(addto, opposite(x), coeff*signed(q2(x), c); is_filtered = true)
    end
    addto
end
