module SimplicialSets

using StructEqualHash, LinearCombinations

using LinearCombinations: Sign, ONE, signed, Zero, sum0, return_type, @Function, unval
using LinearCombinations: coefftype as coeff_type

import LinearCombinations: linear_filter, deg, diff, coprod, hastrait

using Base: Fix1, Fix2, @propagate_inbounds

import Base: show, ==, hash, copy, one, isone, *, /, ^, inv,
    length, firstindex, lastindex, getindex, setindex!

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
