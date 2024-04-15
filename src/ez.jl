#
# Eilenberg-Zilber maps
#

export ez, aw, shih_opp, shih_eml, shih

#
# shuffle map
#

multinomial() = 1
multinomial(k1) = 1
multinomial(k1, ks...) = binomial(k1+sum(ks)::Int, k1)*multinomial(ks...)
# multinomial(k::Int...) = length(k) <= 1 ? 1 : binomial(sum(k)::Int, k[1])*multinomial(k[2:end]...)

function foreach_shuffle_simplex(f, p::Int, q::Int, x0::S, y0::T, m::Int = 0;
        xv = Memory{S}(undef, q+1),
        yv = Memory{T}(undef, q+1),
        kv = Memory{Int}(undef, q+1),
        sv = Memory{Int}(undef, q+1))  where {S <: AbstractSimplex, T <: AbstractSimplex}
    i = 1
    xv[1] = x0
    yv[1] = y0
    kv[1] = 0
    sv[1] = p*q
    @inbounds while i != 0
        if i <= q && kv[i] < p+i
            kv[i+1] = kv[i]+1
            xv[i+1] = s(xv[i] , kv[i]+m)
            yv[i+1] = yv[i]
            sv[i+1] = sv[i]
            i += 1
        else
            if i > q
            # this means i == q+1
                yy = yv[q+1]
                for kk in kv[q+1]:p+q-1
                    yy = s(yy, kk+m)
                end
                f(xv[q+1], yy, sv[q+1])
            end
            i -= 1
            if i != 0
                yv[i] = s(yv[i], kv[i]+m)  # TODO: this fails for SymbolicSimplex in dim p+q close to 24
                kv[i] += 1
                sv[i] -= q+1-i
            end
        end
    end
    nothing
end

_ez(f, addto, coeff, x) = addmul!(addto, x, coeff; is_filtered = true)

function _ez(f::F, addto, coeff, x::S, y::T, z...) where {F,S,T}
    p = dim(x)
    q = dim(y)

    # stack for foreach_shuffle_simplex
    # this is not necessary, but avoids allocations
    xv = Memory{S}(undef, q+1)
    yv = Memory{T}(undef, q+1)
    kv = Memory{Int}(undef, q+1)
    sv = Memory{Int}(undef, q+1)

    # foreach_shuffle_simplex(p, q, x, y) do xx, yy, ss
    foreach_shuffle_simplex(p, q, x, y; xv, yv, kv, sv) do xx, yy, ss
        _ez(f, addto, signed(ss, coeff), f(xx, yy), z...)
    end
end

ez_prod() = ProductSimplex(; dim = 0)
ez_prod(x) = ProductSimplex(x)
ez_prod(x, y) = ProductSimplex(x..., y)

_ez_term_type(f, P) = P
_ez_term_type(f, P, T, U...) = _ez_term_type(f, return_type(f, P, T), U...)

ez_term_type(f) = return_type(f)
ez_term_type(f, T) = return_type(f, T)
ez_term_type(f, T, U...) = _ez_term_type(f, ez_term_type(f, T), U...)

"""
    ez(x::AbstractSimplex...) -> ProductSimplex

    ez(AbstractTensor{T}) where T <: Tuple{Vararg{AbstractSimplex}} -> ProductSimplex{T}

Return the image of the given simplices under the Eilenberg-Zilber or shuffle map.
The first version is multilinear in the simplices and the second one linear in the
tensor argument. Any number of simplices is allowed, including zero.

This function supports the keyword arguments `coefftype`, `addto`, `coeff`,
`sizehint` and `is_filtered` as described for the macro `@linear`.

See also [`aw`](@ref), [`shih`](@ref), `LinearCombinations.@linear`.

# Examples
```jldoctest
julia> x, y = SymbolicSimplex(:x, 1), SymbolicSimplex(:y, 2)
(x[0,1], y[0,1,2])

julia> a = ez(x, y)
(x[0,1,1,1],y[0,0,1,2])-(x[0,0,1,1],y[0,1,1,2])+(x[0,0,0,1],y[0,1,2,2])

julia> ez(Tensor(x, y); addto = a, coeff = -1)
0

julia> z = SymbolicSimplex(:z, 0)
z[0]

julia> ez(x, y, z)
-(x[0,0,1,1],y[0,1,1,2],z[0,0,0,0])+(x[0,0,0,1],y[0,1,2,2],z[0,0,0,0])+(x[0,1,1,1],y[0,0,1,2],z[0,0,0,0])

julia> ez(x, y, z) == ez(ez(x, y), z)
false

julia> ez(x)
(x[0,1])

julia> ez(Tensor())
()
```
"""
@linear_kw function ez(xs::AbstractSimplex...;
        f = ez_prod,
        coefftype = Int,
        addto = zero(Linear{ez_term_type(f, typeof.(xs)...),unval(coefftype)}),
        coeff = one(coeff_type(addto)),
        sizehint = true,
        is_filtered = false)

    isempty(xs) && return addmul!(addto, f(), coeff; is_filtered = true)

    is_filtered || all(!isdegenerate, xs) || return addto

    sizehint && sizehint!(addto, length(addto) + multinomial(map(dim, xs)...))

    x1, xr... = xs
    _ez(f, addto, coeff, f(x1), xr...)
    addto
end

ez(t::AbstractTensor; kw...) = ez(t...; kw...)

hastrait(::typeof(ez), prop::Val, ::Type{<:AbstractTensor{T}}) where T <: Tuple = hastrait(ez, prop, fieldtypes(T)...)

@multilinear ez

deg(::typeof(ez)) = Zero()

keeps_filtered(::typeof(ez), ::Type) = true

#
# Alexander-Whitney map
#

_aw(addto, coeff, ::Tuple{}) = addmul!(addto, Tensor(), coeff)

function _aw(addto, coeff, x::Tuple{AbstractSimplex}, z...)
    isdegenerate(x[1]) || addmul!(addto, Tensor(x[1], z...), coeff; is_filtered = true)
end

@inline function _aw(addto, coeff, x::Tuple, z...)
    n = dim(x[end])
    for k in 0:n
        y = r(x[end], ((k, n),))
        isdegenerate(y) || _aw(addto, coeff, map(Fix2(r, ((0, k),)), x[1:end-1]), y, z...)
    end
end

"""
    aw(x::ProductSimplex{T}) where T <: Tuple{Vararg{AbstractSimplex}} -> Linear{Tensor{T}}

Return the image of the given product simplex under the Alexander-Whitney map.
The product simplex may have any number of components, including zero. The number
of components in the resulting tensor product is the same as the number of
components of the product simplex.

This function is linear and supports the keyword arguments `coefftype`, `addto`, `coeff`
and `is_filtered` as described for the macro `@linear`.

See also [`coprod`](@ref), [`ez`](@ref), [`shih`](@ref), `LinearCombinations.@linear`.

# Examples
```jldoctest
julia> x, y = SymbolicSimplex(:x, 2), SymbolicSimplex(:y, 2)
(x[0,1,2], y[0,1,2])

julia> aw(ProductSimplex(x, y))
x[0]⊗y[0,1,2]+x[0,1,2]⊗y[2]+x[0,1]⊗y[1,2]

julia> z = SymbolicSimplex(:z, 2); aw(ProductSimplex(x, y, z))
x[0,1]⊗y[1]⊗z[1,2]+x[0,1,2]⊗y[2]⊗z[2]+x[0,1]⊗y[1,2]⊗z[2]+x[0]⊗y[0]⊗z[0,1,2]+x[0]⊗y[0,1]⊗z[1,2]+x[0]⊗y[0,1,2]⊗z[2]

julia> aw(ProductSimplex(x))
x[0,1,2]

julia> aw(ProductSimplex(; dim = 0))
()
```
"""
@linear_kw function aw(x::ProductSimplex{T};
        coefftype = Int,
        addto = zero(Linear{Tensor{T},unval(coefftype)}),
        coeff = ONE,
        is_filtered = false) where T <: Tuple{Vararg{AbstractSimplex}}
    if !iszero(coeff) && (is_filtered || !isdegenerate(x))
        _aw(addto, coeff, Tuple(x))
    end
    addto
end

@linear aw

deg(::typeof(aw)) = Zero()

#
# homotopy
#

"""
    shih_opp(z::ProductSimplex{T}) where T <: Tuple{AbstractSimplex,AbstractSimplex} -> Linear{T}

Return the image of the product simplex `z` under the *opposite* Eilenberg-MacLane homotopy
(which is sometimes called the opposite Shih map).

This function is linear and supports the keyword arguments `coefftype`, `addto`, `coeff`,
`sizehint` and `is_filtered` as described for the macro `@linear`.

See also [`aw`](@ref), [`ez`](@ref), [`shih_eml`](@ref), [`opposite`](@ref), [`swap`](@ref),
`LinearCombinations.@linear`.

# Example
```jldoctest
julia> x, y = SymbolicSimplex(:x, 2), SymbolicSimplex(:y, 2)
(x[0,1,2], y[0,1,2])

julia> a = Linear(ProductSimplex(x, y) => 1)   # we use `Linear` to get the correct sign from `opposite` below
(x[0,1,2],y[0,1,2])

julia> shih_opp(a)
-(x[0,1,1,2],y[1,1,2,2])-(x[0,0,0,2],y[0,1,2,2])+(x[0,0,1,2],y[0,1,1,2])+(x[0,0,1,2],y[1,2,2,2])

julia> shih(a)
-(x[0,0,0,1],y[0,1,2,2])+(x[0,0,1,2],y[0,2,2,2])-(x[0,1,1,2],y[0,1,2,2])+(x[0,0,1,1],y[0,1,1,2])

julia> shih_opp(opposite(swap(a))) == opposite(swap(shih(a)))
true
```
"""
@linear_kw function shih_opp(z::ProductSimplex{Tuple{S, T}};
        coefftype = Int,
        addto = zero(Linear{ProductSimplex{Tuple{S,T}},unval(coefftype)}),
        coeff = ONE,
        sizehint = true,
        is_filtered = false) where {S <: AbstractSimplex, T <: AbstractSimplex}

    iszero(coeff) && return addto

    n = dim(z)
    sizehint && sizehint!(addto, length(addto)+(1<<(n+1))-n-2)
    x, y = Tuple(z)

    # stack for foreach_shuffle_simplex
    xv = Memory{S}(undef, n+1)
    yv = Memory{T}(undef, n+1)
    kv = Memory{Int}(undef, n+1)
    sv = Memory{Int}(undef, n+1)

    for p in 0:n-1, q in 0:n-1-p
        xx = r(x, ((0, p), (p+q+1, n)))
        yy = r(y, ((p, n),))

        # if nondeg || !(isdegenerate(xx, 0, p+1) || isdegenerate(yy, 0, q+1))
        if !(isdegenerate(xx, 0, p+1) || isdegenerate(yy, 0, q+1))
            foreach_shuffle_simplex(p, q+1, xx, yy; xv, yv, kv, sv) do xxx, yyy, sss
                @inbounds addmul!(addto, ProductSimplex(xxx, s(yyy, p+q+1); dim = n+1), signed(p+q+sss, coeff); is_filtered = true)
            end
        end
    end

    addto
end

@linear shih_opp

deg(::typeof(shih_opp)) = 1

"""
    shih(z::ProductSimplex{T}) where T <: Tuple{AbstractSimplex,AbstractSimplex} -> Linear{T}
    shih_eml(z::ProductSimplex{T}) where T <: Tuple{AbstractSimplex,AbstractSimplex} -> Linear{T}

Return the image of the product simplex `z` under the Eilenberg-MacLane homotopy
(which is sometimes called the Shih map).

This function is linear and supports the keyword arguments `coefftype`, `addto`, `coeff`,
`sizehint` and `is_filtered` as described for the macro `@linear`.

See also [`aw`](@ref), [`ez`](@ref), [`shih_opp`](@ref), `LinearCombinations.@linear`.

# Examples
```jldoctest
julia> x, y = SymbolicSimplex(:x, 2), SymbolicSimplex(:y, 2); z = ProductSimplex(x, y)
(x[0,1,2],y[0,1,2])

julia> shih(z)
-(x[0,1,1,2],y[0,1,2,2])+(x[0,0,1,1],y[0,1,1,2])-(x[0,0,0,1],y[0,1,2,2])+(x[0,0,1,2],y[0,2,2,2])

julia> shih(ez(x, y))
0

julia> shih(ez(x, y)), aw(shih(z)), shih(shih(z))
(0, 0, 0)
```
"""
@linear_kw function shih_eml(z::ProductSimplex{Tuple{S, T}};
        coefftype = Int,
        addto = zero(Linear{ProductSimplex{Tuple{S,T}},unval(coefftype)}),
        coeff = ONE,
        sizehint = true,
        is_filtered = false) where {S <: AbstractSimplex, T <: AbstractSimplex}

    iszero(coeff) && return addto

    n = dim(z)
    sizehint && sizehint!(addto, length(addto)+(1<<(n+1))-n-2)
    x, y = Tuple(z)

    # stack for foreach_shuffle_simplex
    xv = Memory{S}(undef, n)
    yv = Memory{T}(undef, n)
    kv = Memory{Int}(undef, n)
    sv = Memory{Int}(undef, n)

    for p in 0:n-1, q in 0:n-1-p
        m = n-1-p-q
        xx = r(x, ((0, n-q),))
        yy = r(y, ((0, m), (n-q, n)))

        # if nondeg || !(isdegenerate(xx, m, m+p+1) || isdegenerate(yy, m, m+q+1))
        if !(isdegenerate(xx, m, m+p+1) || isdegenerate(yy, m, m+q+1))
            foreach_shuffle_simplex(p+1, q, s(xx, m), yy, m+1; xv, yv, kv, sv) do xxx, yyy, sss
                @inbounds addmul!(addto, ProductSimplex(xxx, yyy; dim = n+1), signed(m+sss, coeff); is_filtered = true)
            end
        end
    end

    addto
end

@linear shih_eml

deg(::typeof(shih_eml)) = 1

# setting the default version
"""
    const shih = shih_eml

`shih` is a short form for the usual Eilenberg-MacLane homotopy [`shih_eml`](@ref).
"""
const shih = shih_eml

#
# group operations on chains
#

import Base: *, inv

@linear inv

keeps_filtered(::typeof(inv), ::Type) = true

_mul(addto, coeff, t) = ez(Tensor(t); f = *, addto, coeff)

function _mul(addto, coeff, t, a, b...)
    for (x, c) in a
        _mul(addto, coeff*c, (t..., x), b...)
    end
end

function *(a::AbstractLinear{<:AbstractSimplex}...)
# TODO: add addto & coeff
    R = promote_type(map(coefftype, a)...)
    T = return_type(*, map(termtype, a)...)
    addto = zero(Linear{T,R})
    l = prod(map(length, a))
    l == 0 && return addto
    _mul(addto, ONE, (), a...)
    addto
end

#
# diagonal map and coproduct
#

export diag

"""
    diag(x::T) where T <: AbstractSimplex -> ProductSimplex{T,T}

Return the image of `x` under the diagonal map from the simplicial set containing `x`
to the Cartesian product with itself.

See also [`coprod`](@ref).
"""
diag(x::AbstractSimplex) = @inbounds ProductSimplex(x, x)

@linear diag

keeps_filtered(::typeof(diag), ::Type) = true

"""
    coprod(x::T) where T <: AbstractSimplex -> Linear{Tensor{Tuple{T,T}}}

Return the image of the simplex `x` under coproduct (or diagonal map) of the normalized chain complex
containing `x`.

This function is linear and supports the keyword arguments `coefftype`, `addto`, `coeff`
and `is_filtered` as described for the macro `@linear`.

See also [`aw`](@ref), [`diag`](@ref), `LinearCombinations.@linear`.

# Examples
```jldoctest
julia> x = SymbolicSimplex(:x, 2)
x[0,1,2]

julia> coprod(x)
x[0]⊗x[0,1,2]+x[0,1,2]⊗x[2]+x[0,1]⊗x[1,2]

julia> coprod(x) == aw(diag(x))
true
```
"""
coprod(x::AbstractSimplex; kw...) = aw(diag(x); kw...)
