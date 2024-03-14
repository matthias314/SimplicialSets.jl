#
# group stuff
#

export AddToMul, Lattice

import Base: *, /, ^, +, -, zero, iszero, one, isone

#
# AddToMul
#

"""
    AddToMul{T}

A wrapper to turn an additive group structure defined for the type `T`
into a multiplicative structure.

# Examples
```jldoctest
julia> x, y = AddToMul(2), AddToMul(3)
(2, 3)

julia> x*y
5

julia> x^2
4

julia> inv(x)
-2

julia> x/y
-1

julia> one(x)
0
```
"""
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

"""
    Lattice{N}

    Lattice(t::NTuple{N,Integer}) -> Lattice{N}
    Lattice(x::Integer...) -> Lattice

A type representing elements in a lattice (free abelian group) of rank `N`.

# Examples
```jldoctest
julia> x, y = Lattice(1, 2, 3), Lattice(0, -1, 5)
((1, 2, 3), (0, -1, 5))

julia> x+y
(1, 1, 8)

julia> x+y, x-y
((1, 1, 8), (1, 3, -2))

julia> 2*x
(2, 4, 6)

julia> zero(x)
(0, 0, 0)

julia> length(x)
3
```
"""
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
