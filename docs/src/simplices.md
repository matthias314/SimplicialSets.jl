# Simplices

## `AbstractSimplex`

```@docs
AbstractSimplex
```

## `SymbolicSimplex`

```@docs
SymbolicSimplex
vertices
```

## `ProductSimplex`

```@docs
ProductSimplex
Base.Tuple(::ProductSimplex)
length(x::ProductSimplex)
SimplicialSets.cat
SimplicialSets.flatten
swap
LinearCombinations.Regroup
```

## `BarSimplex`

```@docs
BarSimplex
length(::BarSimplex)
one(::BarSimplex)
isone(::BarSimplex)
*(::BarSimplex{T}, ::BarSimplex{T}...) where T
inv(::BarSimplex)
/(::BarSimplex{T}, ::BarSimplex{T}) where T
^(::BarSimplex, ::Integer)
```

## `LoopGroupSimplex`

```@docs
SimplicialSets.LoopGroupGenerator
LoopGroupSimplex
length(::LoopGroupSimplex)
isone(::LoopGroupSimplex)
inv(::LoopGroupSimplex)
SimplicialSets.mul!
*(::LoopGroupSimplex{T}, ::LoopGroupSimplex{T}...) where T <: AbstractSimplex
```

## Other simplices

```@docs
IntervalSimplex
OppositeSimplex
opposite
```
