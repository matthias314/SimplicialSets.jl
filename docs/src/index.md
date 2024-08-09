```@meta
DocTestSetup = doctestsetup
```

# Overview

```@docs
SimplicialSets
```

## Simplices

### `AbstractSimplex`

```@docs
AbstractSimplex
```

### `SymbolicSimplex`

```@docs
SymbolicSimplex
vertices
```

### `ProductSimplex`

```@docs
ProductSimplex
Base.Tuple(::ProductSimplex)
length(x::ProductSimplex)
SimplicialSets.cat
SimplicialSets.flatten
swap
LinearCombinations.Regroup
```

### `BarSimplex`

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

### `LoopGroupSimplex`

```@docs
SimplicialSets.LoopGroupGenerator
LoopGroupSimplex
length(::LoopGroupSimplex)
isone(::LoopGroupSimplex)
inv(::LoopGroupSimplex)
SimplicialSets.mul!
*(::LoopGroupSimplex{T}, ::LoopGroupSimplex{T}...) where T <: AbstractSimplex
```

### Other simplices

```@docs
IntervalSimplex
OppositeSimplex
opposite
```

## Basic functions for simplices

```@docs
dim
deg(::AbstractSimplex)
SimplicialSets.d
SimplicialSets.s
isdegenerate(::AbstractSimplex)
isdegenerate(::AbstractSimplex, ::Integer)
isdegenerate(::AbstractSimplex, ::Integer, ::Integer)
LinearCombinations.linear_filter(::AbstractSimplex)
LinearCombinations.diff(::AbstractSimplex)
diag
LinearCombinations.coprod
```

## Eilenberg–Zilber maps

```@docs
ez
aw
shih
shih_opp
```

## Surjection operad and interval cuts

```@docs
Surjection
arity
deg(::Surjection)
SimplicialSets.is_surjection
isdegenerate(::Surjection)
LinearCombinations.linear_filter(::Surjection)
LinearCombinations.diff(::Surjection)
```

## Helper types for groups

```@docs
AddToMul
Lattice
Tuple(::Lattice)
```
