"""
    SimplicialSets

This packages provides functions to work with simplicial sets.
Various kinds of simplicial sets are supported, including symbolic simplices,
products, bar constructions and Kan loop groups.
The Eilenberg-Zilber maps and interval cut operations are also implemented.

The package uses `LinearCombinations.jl` to represent formal linear combinations
of simplices and other objects.
"""
module SimplicialSets

using StructEqualHash, LinearCombinations

using LinearCombinations: Sign, ONE, signed, Zero, sum0, return_type, @Function, unval
using LinearCombinations: coefftype as coeff_type

import LinearCombinations: linear_filter, deg, diff, coprod, hastrait

using Base: Fix1, Fix2, @propagate_inbounds

import Base: show, ==, hash, copy, one, isone, *, /, ^, inv,
    length, firstindex, lastindex, getindex, setindex!

if VERSION < v"1.11.0-"
    const Memory = Vector
end

include("abstract.jl")
include("product.jl")
include("opposite.jl")
include("interval.jl")
include("suspension.jl")
include("symbolic.jl")
include("loopgroup.jl")
include("bar.jl")
include("groups.jl")
include("twistedproduct.jl")
# include("szczarba.jl")
include("ez.jl")
include("surj.jl")
include("helpers.jl")

end
