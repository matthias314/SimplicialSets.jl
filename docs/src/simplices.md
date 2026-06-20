# Simplices

## `AbstractSimplex`

```@docs
AbstractSimplex
⋅
*(::AbstractSimplex, ::AbstractSimplex...)
```

## `SymbolicSimplex`

```@docs
SymbolicSimplex
vertices
```

## Product simplices

```@docs
AbstractProductSimplex
Base.Tuple(::AbstractProductSimplex)
length(::AbstractProductSimplex)
fieldtypes
```

### `ProductSimplex`

```@docs
ProductSimplex
SimplicialSets.cat
SimplicialSets.flatten
swap
LinearCombinations.Regroup
```

### Twisted product simplices

```@docs
LeftTwistedProductSimplex
RightTwistedProductSimplex
lefttwistedproductsimplex
righttwistedproductsimplex
```

## `BarSimplex`

```@docs
BarSimplex
length(::BarSimplex)
one(::BarSimplex)
isone(::BarSimplex)
⋅(::BarSimplex{T}, ::BarSimplex{T}...) where T
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
⋅(::LoopGroupSimplex{T}, ::LoopGroupSimplex{T}...) where T <: AbstractSimplex
```

## Other simplices

```@docs
IntervalSimplex
OppositeSimplex
opposite
```
