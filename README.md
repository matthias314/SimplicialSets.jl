# SimplicialSets.jl

This packages provides types and functions to work with simplicial sets. Various kinds
of simplicial sets are supported, including symbolic simplices, products,
bar constructions and Kan loop groups.
The Eilenberg–Zilber maps and interval cut operations are also implemented.

The package uses [LinearCombinations.jl](https://github.com/matthias314/LinearCombinations.jl)
to represent formal linear combinations of simplices. By default, coefficients are of type `Int`.

The package comes with an extensive [documentation](https://matthias314.github.io/SimplicialSets.jl/dev/).

## Examples

### Basic operations and symbolic simplices
```julia
julia> using LinearCombinations, SimplicialSets

julia> x = SymbolicSimplex(:x, 4)
x[0,1,2,3,4]

julia> dim(x)
4

julia> using SimplicialSets: d, s

julia> d(x, 2)
x[0,1,3,4]

julia> s(x, 2)
x[0,1,2,2,3,4]

julia> isdegenerate(x), isdegenerate(s(x, 2))
(false, true)

julia> using LinearCombinations: diff

julia> diff(x)
x[0,1,2,3]-x[0,2,3,4]+x[1,2,3,4]+x[0,1,3,4]-x[0,1,2,4]
```

### Loop groups
```julia
julia> x, y = SymbolicSimplex(:x, 3), SymbolicSimplex(:y, 3)
(x[0,1,2,3], y[0,1,2,3])

julia> u, v = LoopGroupSimplex(x), LoopGroupSimplex(y)
(⟨x[0,1,2,3]⟩, ⟨y[0,1,2,3]⟩)

julia> w = u*v
⟨x[0,1,2,3],y[0,1,2,3]⟩

julia> inv(w)
⟨y[0,1,2,3]⁻¹,x[0,1,2,3]⁻¹⟩

julia> diff(w)
-⟨x[0,1,3],y[0,1,3]⟩+⟨x[0,1,2],y[0,1,2]⟩
+⟨x[1,2,3]⁻¹,x[0,2,3],y[1,2,3]⁻¹,y[0,2,3]⟩
```

### Eilenberg–Zilber maps
```julia
julia> x, y = SymbolicSimplex(:x, 2), SymbolicSimplex(:y, 2)
(x[0,1,2], y[0,1,2])

julia> z = ProductSimplex(x, y)
(x[0,1,2],y[0,1,2])

julia> ez(x, y)   # shuffle map
-(x[0,0,1,1,2],y[0,1,1,2,2])+(x[0,1,2,2,2],y[0,0,0,1,2])
+(x[0,0,1,2,2],y[0,1,1,1,2])+(x[0,0,0,1,2],y[0,1,2,2,2])
+(x[0,1,1,1,2],y[0,0,1,2,2])-(x[0,1,1,2,2],y[0,0,1,1,2])

julia> aw(z)   # Alexander–Whitney map
x[0,1]⊗y[1,2]+x[0]⊗y[0,1,2]+x[0,1,2]⊗y[2]

julia> shih(z)    # Eilenberg–MacLane homotopy
-(x[0,1,1,2],y[0,1,2,2])+(x[0,0,1,1],y[0,1,1,2])
-(x[0,0,0,1],y[0,1,2,2])+(x[0,0,1,2],y[0,2,2,2])
```
Let's check that `shih` is indeed a homotopy from the identity to `ez∘aw`:
```julia
julia> diff(shih(z)) + shih(diff(z)) == ez(aw(z)) - z
true
```
Let's verify the "side conditions" for the Eilenberg–Zilber maps:
```julia
julia> shih(ez(x, y)), aw(shih(z)), shih(shih(z))
(0, 0, 0)
```
Let's check that the shuffle map is commutative:
```julia
julia> x, y = SymbolicSimplex(:x, 1), SymbolicSimplex(:y, 3)
(x[0,1], y[0,1,2,3])

julia> t = tensor(x, y)
x[0,1]⊗y[0,1,2,3]

julia> ez(t)
-(x[0,0,1,1,1],y[0,1,1,2,3])+(x[0,0,0,1,1],y[0,1,2,2,3])
+(x[0,1,1,1,1],y[0,0,1,2,3])-(x[0,0,0,0,1],y[0,1,2,3,3])

julia> a = swap(ez(t))
-(y[0,1,2,3,3],x[0,0,0,0,1])-(y[0,1,1,2,3],x[0,0,1,1,1])
+(y[0,0,1,2,3],x[0,1,1,1,1])+(y[0,1,2,2,3],x[0,0,0,1,1])

julia> swap(t)
-y[0,1,2,3]⊗x[0,1]

julia> b = ez(swap(t))
-(y[0,1,2,3,3],x[0,0,0,0,1])-(y[0,1,1,2,3],x[0,0,1,1,1])
+(y[0,0,1,2,3],x[0,1,1,1,1])+(y[0,1,2,2,3],x[0,0,0,1,1])

julia> a == b
true
```

### Interval cut operations
```julia
julia> sj = Surjection([1,3,2,1])
Surjection{3}([1, 3, 2, 1])

julia> x = SymbolicSimplex(:x, 2)
x[0,1,2]

julia> sj(x)
-x[0,1,2]⊗x[1]⊗x[0,1]-x[0,2]⊗x[1,2]⊗x[0,1]-x[0,1,2]⊗x[0,1]⊗x[0]
-x[0,1,2]⊗x[2]⊗x[1,2]+x[0,2]⊗x[0,1,2]⊗x[0]+x[0,2]⊗x[2]⊗x[0,1,2]
-x[0,1,2]⊗x[1,2]⊗x[1]

julia> diff(sj)
Surjection{3}([2, 3, 1])-Surjection{3}([1, 2, 3])

julia> diff(sj(x)) == diff(sj)(x) + (-1)^deg(sj) * sj(diff(x))
true
```
