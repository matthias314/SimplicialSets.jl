#
# OppositeSimplex datatype
#

export OppositeSimplex, opposite, qsign

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

*(ys::OppositeSimplex...) = OppositeSimplex(*((y.x for y in ys)...))
/(y::OppositeSimplex, z::OppositeSimplex) = OppositeSimplex(y.x/z.x)
inv(y::OppositeSimplex) = OppositeSimplex(inv(y.x))
one(::Type{OppositeSimplex{T}}, n...) where T = OppositeSimplex(one(T, n...))
one(y::OppositeSimplex, n...) = OppositeSimplex(one(y.x, n...))

opposite(x::AbstractSimplex) = OppositeSimplex(x)
opposite(y::OppositeSimplex) = y.x
opposite(x::ProductSimplex) = ProductSimplex((opposite(y) for y in x)...)
opposite(x::Tensor) = Tensor((opposite(y) for y in x)...)

q2(x::AbstractSimplex) = ifelse((dim(x)+1) & 2 == 0, 0, 1)
# this is 0 if dim(x) == 0 or 3 mod 4, and 1 if dim(x) == 1 or 2 mod 4

q2(t::Tensor) = sum0(map(q2, factors(t)))

@linear_kw function opposite(a::Linear{T,R};
        coefftype = R,
        addto = zero(Linear{return_type(opposite, T),coefftype}),
        coeff = ONE) where {T,R}
    iszero(coeff) && return addto
    for (x, c) in a
        addmul!(addto, opposite(x), coeff*signed(q2(x), c); is_filtered = true)
    end
    addto
end
