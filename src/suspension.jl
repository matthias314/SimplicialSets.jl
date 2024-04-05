#
# simplicial suspension
#

export SuspensionSimplex

struct SuspensionSimplex{T <: AbstractSimplex} <: AbstractSimplex
    i::IntervalSimplex
    x::T
    function SuspensionSimplex{T}(i::IntervalSimplex) where T <: AbstractSimplex
        isabovevertex(i) ? new{T}(i) : error("illegal arguments")
    end
    function SuspensionSimplex(x::T, i::IntervalSimplex) where T <: AbstractSimplex
        if dim(x) == dim(i)
            isabovevertex(i) ? new{T}(i) : new{T}(i, x)
        else
            error("illegal arguments")
        end
    end
end

SuspensionSimplex(x::T, p::Int, q::Int = dim(x)+1-p) where T <: AbstractSimplex = SuspensionSimplex(x, IntervalSimplex(p, q))

SuspensionSimplex(x::ProductSimplex{Tuple{T, IntervalSimplex}}) where T <: AbstractSimplex = SuspensionSimplex(Tuple(x)...)

function show(io::IO, x::SuspensionSimplex)
    print(io, isabovevertex(x) ? "Σ($(x.i))" : "Σ($(x.x),$(x.i))")
end

dim(x::SuspensionSimplex) = dim(x.i)

isabovevertex(x::SuspensionSimplex) = isabovevertex(x.i)

function ==(x::SuspensionSimplex{T}, y::SuspensionSimplex{T}) where T <: AbstractSimplex
    x.i == y.i && (isabovevertex(x) || x.x == y.x)
end

function hash(x::SuspensionSimplex, h::UInt)
    h = hash(x.i, h)
    if isabovevertex(x)
        h
    else
        hash(x.x, h)
    end
end

function d(x::SuspensionSimplex{T}, k::Integer) where T <: AbstractSimplex
    if isabovevertex(x)
        SuspensionSimplex{T}(d(x.i, k))
    else
        SuspensionSimplex(d(x.x, k), d(x.i, k))
    end
end

function s(x::SuspensionSimplex{T}, k::Integer) where T <: AbstractSimplex
    if isabovevertex(x)
        SuspensionSimplex{T}(s(x.i, k))
    else
        SuspensionSimplex(s(x.x, k), s(x.i, k))
    end
end

function isdegenerate(x::SuspensionSimplex, k::Integer)
    if isabovevertex(x)
        dim(x) > 0
    else
        isdegenerate(x.i, k) && isdegenerate(x.x, k)
    end
end

function isdegenerate(x::SuspensionSimplex)
    if isabovevertex(x)
        dim(x) > 0
    else
        any(k -> isdegenerate(x.i, k) && isdegenerate(x.x, k), 0:dim(x)-1)
    end
end
