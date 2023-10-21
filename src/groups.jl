#
# group stuff
#

export AddToMul, Lattice

import Base: *, /, ^, +, -, zero, iszero, one, isone

#
# AddToMul
#

struct AddToMul{T}
    x::T
end

show(io::IO, a::AddToMul) = print(io, a.x)

one(::Type{AddToMul{T}}) where T = AddToMul(zero(T))
one(a::AddToMul) = AddToMul(zero(a.x))

isone(a::AddToMul) = iszero(a.x)

*(a::AddToMul{T}...) where T = AddToMul(mapreduce(b -> b.x, +, a))

/(a::AddToMul{T}, b::AddToMul{T}) where T = AddToMul(a.x-b.x)

inv(a::AddToMul) = AddToMul(-a.x)

^(a::AddToMul, n::Integer) = AddToMul(n * a.x)

#
# Lattice
#

struct Lattice{N}
    v::NTuple{N, Int}
end

Lattice(ii::Int...) = Lattice(ii)

show(io::IO, g::Lattice) = print(io, g.v)

length(::Lattice{N}) where N = N

iterate(g::Lattice, state...) = iterate(g.v, state...)

zero(::Type{Lattice{N}}) where N = Lattice(ntuple(Returns(0), N))
zero(::T) where T <: Lattice = zero(T)

# this doesn't seem to be faster than the default g.v == zero(g.v)
iszero(g::Lattice) = all(iszero, g.v)

+(g::Lattice{N}...) where N = Lattice(map(+, map(h -> h.v, g)...))

-(g::Lattice) = Lattice(.- g.v)

-(g::Lattice{N}, h::Lattice{N}) where N = Lattice(g.v .- h.v)

*(n::Integer, g::Lattice) = Lattice(n .* g.v)
