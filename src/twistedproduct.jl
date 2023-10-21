#
# twisted Cartesian products
#

export TwistedProductSimplex

struct TwistedProductSimplex{TWF, B<:AbstractSimplex, F<:AbstractSimplex} <: AbstractSimplex
    b::B
    f::F
    twf::TWF
    @inline function TwistedProductSimplex(b::B, f::F, twf::TWF) where {B <: AbstractSimplex, F <: AbstractSimplex, TWF}
        @boundscheck if dim(b) != dim(f)
            error("simplices must be of the same dimension")
        end
        new{TWF, B, F}(b, f, twf)
    end
end

show(io::IO, x::TwistedProductSimplex) = print(io, "($(x.b),$(x.f))")

==(x::T, y::T) where T <: TwistedProductSimplex{TWF,B,F} where {TWF,B,F} = x.b == y.b && x.f == y.f

hash(x::TwistedProductSimplex, h::UInt) = hash(x.f, hash(x.b, h))

dim(x::TwistedProductSimplex) = dim(x.b)

function d(x::TwistedProductSimplex, k)
    if k == 0
        f = x.twf(x.b) * d(x.f, 0)
    else
        f = d(x.f, k)
    end
    TwistedProductSimplex(d(x.b, k), f, x.twf)
end

function s(x::TwistedProductSimplex, k)
    TwistedProductSimplex(s(x.b, k), s(x.f, k), x.twf)
end

@propagate_inbounds isdegenerate(x::TwistedProductSimplex, k::Int) = isdegenerate(x.b, k) && isdegenerate(x.f, k)
