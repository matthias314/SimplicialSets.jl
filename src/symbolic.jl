#
# SymbolicSimplex datatype
#

export SymbolicSimplex, dim, vertices

# implementation using UInt128

const Label = Union{Symbol,Char}

struct SymbolicSimplex{L<:Label} <: AbstractSimplex
    label::L
    dim::Int
    v::UInt128
    function SymbolicSimplex(label::L, dim, v) where L <: Label
        @boundscheck if dim > 24
            error("SymbolicSimplex is limited to dimension at most 24")
        end
        new{L}(label, dim, v)
    end
end

function SymbolicSimplex(label::Label, w::AbstractVector{<:Integer})
    v = UInt128(0)
    for k in reverse(w)
        0 <= k < 32 || error("vertex numbers must be between 0 and 31")
        v = 32*v+k
    end
    SymbolicSimplex(label, length(w)-1, v)
end

SymbolicSimplex(label::Label, n::Integer) = SymbolicSimplex(label, 0:n)

dim(x::SymbolicSimplex) = x.dim

# copy(x::SymbolicSimplex) = SymbolicSimplex(x.label, x.dim, x.v)
copy(x::SymbolicSimplex) = x

to_uint(x::Symbol) = objectid(x)
to_uint(x::Char) = UInt(1073741831)*UInt(x)

# hash(x::SymbolicSimplex, h::UInt) = hash(x.v, hash(256*x.dim+Int(x.label), h))
function Base.hash(x::SymbolicSimplex, h::UInt)
    m = UInt(1073741827)*(x.dim % UInt) + to_uint(x.label)
    # the constants are primes with product 0x1000000280000015
    # which is greater than the 60 = 12*5 bits needed for 12 vertices
    if x.dim <= 12
        hash(x.v % UInt + m, h)
    else
        hash(x.v + m, h)
    end
end

function Base.:(==)(x::SymbolicSimplex, y::SymbolicSimplex)
    # x.label == y.label &&
    x.dim == y.dim && x.v == y.v
end

function vertices(x::SymbolicSimplex)
    d = x.v
    m = dim(x)+1
    s = Vector{Int}(undef, m)
    for k in 1:m
        d, r = divrem(d, UInt128(1) << 5)
        s[k] = r
    end
    s
end

function show(io::IO, x::SymbolicSimplex)
    print(io, x.label, '[')
    join(io, vertices(x), ',')
    print(io, ']')
end

const bitmask = [[UInt128(0)]; [UInt128(1) << (5*k) - UInt128(1) for k in 1:25]]

function d(x::SymbolicSimplex, k::Integer)
    n = dim(x)
    @boundscheck if k < 0 || k > n || n == 0
        error("index outside the allowed range 0:$(n == 0 ? -1 : n)")
    end
    @inbounds w = (x.v & bitmask[k+1]) | (x.v & ~bitmask[k+2]) >> 5
    @inbounds SymbolicSimplex(x.label, n-1, w)
end

# missing: d for kk a collections

function r(x::SymbolicSimplex, kk)
    isempty(kk) && error("at least one interval must be given")
    n = zero(first(kk[1]))
    w = x.v & bitmask[last(kk[end])+2]
    for i in length(kk)-1:-1:1
        ka = last(kk[i])
        kb = first(kk[i+1])
        w = (w & bitmask[ka+2]) | ((w & ~bitmask[kb+1]) >> (5*(kb-ka-1)))
        n += last(kk[i+1])-kb+1
    end
    w >>= 5*first(kk[1])
    n += last(kk[1])-first(kk[1])
    # println(n, typeof(n))
    SymbolicSimplex(x.label, n, w)
end

function s(x::SymbolicSimplex, k::Integer)
    n = dim(x)
    if n == 24
        error("SymbolicSimplex is limited to dimension at most 24")
    end
    @boundscheck if k < 0 || k > n
        error("index outside the allowed range 0:$n")
    end
    @inbounds w = (x.v & bitmask[k+2]) | (x.v & ~bitmask[k+1]) << 5
    @inbounds SymbolicSimplex(x.label, n+1, w)
end

function s(x::SymbolicSimplex, kk::AbstractVector{<:Integer})
    w = UInt128(0)
    v = x.v
    l = length(kk)
    n = dim(x)
    if n+l > 24
        error("SymbolicSimplex is limited to dimension at most 24")
    end
    for i in 1:l
        @inbounds k = kk[i]+1
        @boundscheck if k > n+i
            error("indices outside the allowed range")
        end
        @inbounds w |= v & bitmask[k+1]
        @inbounds v = (v & ~bitmask[k]) << 5
    end
    @inbounds SymbolicSimplex(x.label, n+l, w | v)
end

@inline function isdegenerate(x::SymbolicSimplex, k::Integer)
    @boundscheck if k < 0 || k >= dim(x)
        error("illegal arguments")
    end
    @inbounds xor(x.v, x.v >> 5) & bitmask[k+2] & ~bitmask[k+1] == 0
end
