export LeftTwistedProductSimplex, RightTwistedProductSimplex

abstract type AbstractTwistedProductSimplex{X,Y,TWF} <: AbstractProductSimplex{Tuple{X,Y}} end

Base.Tuple(x::AbstractTwistedProductSimplex) = (x.x, x.y)

dim(z::AbstractTwistedProductSimplex) = dim(z.x)

s(z::TWP, k) where TWP <: AbstractTwistedProductSimplex = TWP(z.twf, s(z.x, k), s(z.y, k))

"""
    LeftTwistedProductSimplex{X,Y,TWF}

See also [`RightTwistedProductSimplex`](@ref),  [`lefttwistedproductsimplex`](@ref).
"""
struct LeftTwistedProductSimplex{X,Y,TWF} <: AbstractTwistedProductSimplex{X,Y,TWF}
    twf::TWF
    x::X
    y::Y
end

"""
    RightTwistedProductSimplex{X,Y,TWF}

See also [`LeftTwistedProductSimplex`](@ref),  [`righttwistedproductsimplex`](@ref).
"""
struct RightTwistedProductSimplex{X,Y,TWF} <: AbstractTwistedProductSimplex{X,Y,TWF}
    twf::TWF
    x::X
    y::Y
end

function show(io::IO, z::LeftTwistedProductSimplex)
    print(io, '(', z.x, ',', z.y, ")˱")
end

function show(io::IO, z::RightTwistedProductSimplex)
    print(io, '(', z.x, ',', z.y, ")˲")
end

@struct_equal_hash LeftTwistedProductSimplex
@struct_equal_hash RightTwistedProductSimplex

export lefttwistedproductsimplex, righttwistedproductsimplex

"""
    lefttwistedproductsimplex(twf)

Return a callable object that converts an `AbstractProductSimplex` to a `LeftTwistedProductSimplex`
with twisting function `twf`. The callable object is linear.

See also [`LeftTwistedProductSimplex`](@ref), [`righttwistedproductsimplex`](@ref).
"""
function lefttwistedproductsimplex(twf)
    @linear f
    f(z::AbstractProductSimplex{<:NTuple{2,AbstractSimplex}}) = LeftTwistedProductSimplex(twf, z[1], z[2])
end

"""
    righttwistedproductsimplex(twf)

Return a callable object that converts an `AbstractProductSimplex` to a `RightTwistedProductSimplex`
with twisting function `twf`. The callable object is linear.

See also [`RightTwistedProductSimplex`](@ref), [`lefttwistedproductsimplex`](@ref).
"""
function righttwistedproductsimplex(twf)
    @linear f
    f(z::AbstractProductSimplex{<:NTuple{2,AbstractSimplex}}) = RightTwistedProductSimplex(twf, z[1], z[2])
end

function d(z::LeftTwistedProductSimplex, k)
    if k == 0
        LeftTwistedProductSimplex(z.twf, d(z.x, k)⋅inv(z.twf(z.y)), d(z.y, k))
    else
        LeftTwistedProductSimplex(z.twf, d(z.x, k), d(z.y, k))
    end
end

function d(z::RightTwistedProductSimplex, k)
    if k == 0
        RightTwistedProductSimplex(z.twf, d(z.x, k), z.twf(z.x)⋅d(z.y, k))
    else
        RightTwistedProductSimplex(z.twf, d(z.x, k), d(z.y, k))
    end
end

#
# Shih
#

export Contraction, LeftTwistedProductPerturbation, RightTwistedProductPerturbation,
    ShihContractionF, ShihContractionG, ShihContractionH, ShihTwc,
    shih_contraction, new_diff

struct Contraction{F,G,H}
    f::F
    g::G
    h::H
end

@struct_equal_hash Contraction

function test_contraction(ct::Contraction, a, b)
    ct.g(ct.f(a)) == (a isa Linear ? a : Linear(a => 1)) &&
    diff(ct.h(b)) + ct.h(diff(b)) == ct.f(ct.g(b)) - b &&
    iszero(ct.h(ct.f(a))) &&
    iszero(ct.g(ct.h(b))) &&
    iszero(ct.h(ct.h(b)))
end

struct LeftTwistedProductPerturbation{TWF}
    twf::TWF
end

@struct_equal_hash LeftTwistedProductPerturbation

@linear p::LeftTwistedProductPerturbation

@linear_kw function (p::LeftTwistedProductPerturbation)(z::P;
        coefftype = Int,
        addto = zero(Linear{P,unval(coefftype)}),
        coeff = 1) where P <: ProductSimplex{<:Tuple{AbstractSimplex,AbstractSimplex}}
    dim(z) == 0 && return addto
    x, y = z
    addmul!(addto, ProductSimplex(d(x, 0)⋅inv(p.twf(y)), d(y, 0)), coeff)
    addmul!(addto, d(z, 0), -coeff)
end

struct RightTwistedProductPerturbation{TWF}
    twf::TWF
end

@struct_equal_hash RightTwistedProductPerturbation

@linear p::RightTwistedProductPerturbation

@linear_kw function (p::RightTwistedProductPerturbation)(z::P;
        coefftype = Int,
        addto = zero(Linear{P,unval(coefftype)}),
        coeff = 1) where P <: ProductSimplex{<:Tuple{AbstractSimplex,AbstractSimplex}}
    dim(z) == 0 && return addto
    x, y = z
    addmul!(addto, ProductSimplex(d(x, 0), p.twf(x)⋅d(y, 0)), coeff)
    addmul!(addto, d(z, 0), -coeff)
end

struct ShihContractionF{CT<:Contraction,P}
    ct::CT
    p::P
end

@struct_equal_hash ShihContractionF

# @linear f::ShihContractionF

function (f::ShihContractionF)(x)
    a = f.ct.f(x)
    b = a
    while begin b = f.ct.h(f.p(b)); !iszero(b) end
        add!(a, b)
    end
    a
end

struct ShihContractionG{CT<:Contraction,P}
    ct::CT
    p::P
end

@struct_equal_hash ShihContractionG

# @linear G::ShihContractionG

function (g::ShihContractionG)(y)
    a = y isa Linear ? copy(y) : Linear(y => 1)
    b = a
    while begin b = g.p(g.ct.h(b)); !iszero(b) end
        add!(a, b)
    end
    g.ct.g(a)
end

struct ShihContractionH{CT<:Contraction,P}
    ct::CT
    p::P
end

@struct_equal_hash ShihContractionH

# @linear h::ShihContractionH

function (h::ShihContractionH)(y)
    a = h.ct.h(y)
    b = a
    while begin b = h.ct.h(h.p(b)); !iszero(b) end
        add!(a, b)
    end
    a
end

shih_contraction(ct::Contraction, p) =
    Contraction(ShihContractionF(ct, p), ShihContractionG(ct, p), ShihContractionH(ct, p))

function new_diff(ct::Contraction, p, x)
    diff(x) + ShihContractionG(ct, p)(p(ct.f(x)))
end

struct ShihTwc{CT,TWF}
    ct::CT
    twf::TWF
end

@linear st::ShihTwc

deg(::ShihTwc) = -1

function (st::ShihTwc)(x::T) where T <: AbstractSimplex
    G = return_type(st.twf, T)
    # dim(x) == 0 && return zero(Linear{G,Int})
    # ct = Contraction(ez, aw, shih_opp)
    p = LeftTwistedProductPerturbation(st.twf)
    a = ProductSimplex(one(G, dim(x)), x) |> p |> ShihContractionG(st.ct, p)
    Linear(g => dim(x) == 0 ? c : 0 for ((g, x), c) in a)
end

#
# Szczarba operators
#

export SzczarbaTwc, SzczarbaShuffle

function szczarba(x::AbstractSimplex, ii::AbstractVector, k::Int, m::Int = length(ii), l::Int = 0)
# m is the length of the valid part of ii
# l keeps track of how often the simplicial operator has been derived
    m == 0 && return x
    i1 = ii[m]
    if k < i1
        szczarba(s(d(x, i1-k+l), l), ii, k, m-1, l+1)
    elseif k == i1
        szczarba(x, ii, k, m-1, l+1)
    else
        szczarba(s(x, l), ii, k-1, m-1, l+1)
    end
end

function szczarba_twc(twf, ii::AbstractVector, x::AbstractSimplex)
    n = dim(x)
    @assert n != 0
    @assert length(ii) == n-1
    # m = length(ii)
    g = szczarba(twf(x), ii, 0)
    for k in 1:n-1
        x = d(x, 0)
        g = szczarba(twf(x), ii, k) ⋅ g
    end
    inv(g)
end

function foreach_szczarba(f, n::Int)
# Note: ii is reversed: ii[1] == 0 and 0 <= ii[n] <= n-1
    ii = zeros(Int, n)
    n == 0 && begin f(ii); return nothing end
    while true
        f(ii)
        k = n
        while ii[k] == k-1
            ii[k] = 0
            k -= 1
            k == 0 && return
        end
        ii[k] += 1
    end
end

"""
    SzczarbaTwc(twf)

Construct Szczarba's twisting cochain from the twisting function `twf`.
The resulting twisting cochain can be called with a simplex or a linear combination
of simplices as argument.

# Example
```jldoctest
julia> SymbolicSimplex(:x, 3);

julia> sz = SzczarbaTwc(twf_loop);

julia> sz(x)
Linear{LoopGroupSimplex{SymbolicSimplex{Symbol}}, Int64} with 2 terms:
⟨x[0,1,2,3]⁻¹,x[1,2,2,3]⁻¹,x[2,3,3,3]⁻¹⟩-⟨x[0,1,1,3]⁻¹,x[1,2,3,3]⁻¹,x[2,3,3,3]⁻¹⟩

julia> using LinearCombinations: diff

julia> a = diff(sz(x)) + sz(diff(x))
Linear{LoopGroupSimplex{SymbolicSimplex{Symbol}}, Int64} with 4 terms:
-⟨x[0,1,2]⁻¹,x[1,2,2]⁻¹⟩+⟨x[0,1,2]⁻¹,x[1,2,2]⁻¹,x[2,3,3]⁻¹⟩-⟨x[0,1,1]⁻¹,x[1,2,3]⁻¹,x[2,3,3]⁻¹⟩+⟨x[1,2,3]⁻¹,x[2,3,3]⁻¹⟩

julia> b = x |> coprod |> Tensor(sz, sz) |> TensorSplat(*)
Linear{LoopGroupSimplex{SymbolicSimplex{Symbol}}, Int64} with 4 terms:
-⟨x[0,1,2]⁻¹,x[1,2,2]⁻¹⟩+⟨x[0,1,2]⁻¹,x[1,2,2]⁻¹,x[2,3,3]⁻¹⟩-⟨x[0,1,1]⁻¹,x[1,2,3]⁻¹,x[2,3,3]⁻¹⟩+⟨x[1,2,3]⁻¹,x[2,3,3]⁻¹⟩

julia> a == b
true
```
"""
struct SzczarbaTwc{F}
    twf::F
end

@linear sz::SzczarbaTwc

@linear_kw function (sz::SzczarbaTwc)(x::T;
        coefftype = Int,
        addto = begin
            G = return_type(sz.twf, T)
            zero(Linear{G,unval(coefftype)})
        end,
        coeff = 1) where T <: AbstractSimplex
    coefftype = unval(coefftype)
    n = dim(x)
    if n > 1
        foreach_szczarba(n-1) do ii
            m = has_char2(coefftype) ? Zero() : sum(ii)
            addmul!(addto, szczarba_twc(sz.twf, ii, x), withsign(m, coeff))
        end
    elseif n == 1
        g = inv(sz.twf(x))
        addmul!(addto, g, coeff)
        addmul!(addto, one(g), -coeff)
    end
    addto
end

deg(::SzczarbaTwc) = 1

"""
    SzczarbaShuffle(twf)

Construct Szczarba's twisting shuffle map from the twisting function `twf`.
The resulting map can be called with a `Tensor` or a linear combination
of tensors as argument. Note that the base space corresponds to the first
tensor factor and the fiber to the second.

# Example
```jldoctest
julia> x = SymbolicSimplex(:x, 2); y = one(typeof(twf_loop(x)));

julia> Tensor(x, y) |> SzczarbaShuffle(twf_loop)
Linear{RightTwistedProductSimplex{typeof(twf_loop), SymbolicSimplex{Symbol}, LoopGroupSimplex{SymbolicSimplex{Symbol}}}, Int64} with 2 terms:
-(x[0,0,2],⟨x[0,1,2,2]⁻¹,x[1,2,2,2]⁻¹⟩)˲+(x[0,1,2],⟨x[0,1,1,2]⁻¹,x[1,2,2,2]⁻¹⟩)˲
```
"""
struct SzczarbaShuffle{F}
    twf::F
end

function szczarba_shuffle(twf, ii::AbstractVector, x::AbstractSimplex, y::AbstractSimplex)
    n = dim(x)
    n == 0 && return ProductSimplex(x, y)
    # m = length(ii)
    xx = szczarba(x, ii, 0)
    g = szczarba(twf(x), ii, 1)
    for k in 2:n
        x = d(x, 0)
        g = szczarba(twf(x), ii, k) ⋅ g
    end
    yy = inv(g) ⋅ y
    ProductSimplex(xx, yy)
end

szczarba_mul(twf, xg) = xg

szczarba_mul(twf, (x, g), y) = RightTwistedProductSimplex(twf, x, g⋅y)

@linear szsh::SzczarbaShuffle

function (szsh::SzczarbaShuffle{TWF})(t::RightTwistedTensor{X,Y};
        coefftype = Int,
        addto = zero(Linear{RightTwistedProductSimplex{X,Y,TWF},unval(coefftype)}),
        coeff = 1) where {X <: AbstractSimplex, Y <: AbstractSimplex, TWF}
    G = return_type(szsh.twf, X)
    x, y = t
    n = dim(x)
    g = one(G, n)
    foreach_szczarba(n) do ii
        m = has_char2(coefftype) ? Zero() : sum(ii)
        ez(szczarba_shuffle(szsh.twf, ii, x, g), y;
            addto, coeff = withsign(m, coeff), f = Fix1(szczarba_mul, szsh.twf))
    end
    addto
end

#
# coproduct
#

Eop_right(k) = Surjection([isodd(i) ? k+1 : i÷2 for i in 1:2*k+1])

function FFF_right(twc, x::X, y::Y) where {X <: AbstractSimplex, Y <: AbstractSimplex}
    addto = zero(Linear{Tensor{Tuple{Y,X}},Int})
    for k in 0:dim(x)
        for (xx, c) in Eop_right(k)(x)
            yy = foldr(*, map(twc, xx[1:end-1]); init = y)
            m = deg(xx[end])*deg(y) +
                # (k*(k-1))÷2 + sum(deg(xx[i]) for i in k:-2:1; init = 0) +
                # sum(deg, xx[1:k]; init = 0) - k
                (k*(k+1))÷2 + sum(deg(xx[i]) for i in k-1:-2:1; init = 0)
            addmul!(addto, tensor(yy, xx[end]), withsign(m, c))
        end
    end
    addto
end

@linear_kw function coprod(t::T;
        coefftype = Int,
        addto = zero(Linear{Tensor{Tuple{T, T}}, unval(coefftype)}),
        coeff = ONE) where T <: RightTwistedTensor{<:AbstractSimplex, <:AbstractSimplex, <:SzczarbaTwc}
    for ((x1, x2), c1) in coprod(t.x), ((y1, y2), c2) in coprod(t.y)
        for ((yy, xx), c3) in FFF_right(t.twc, x2, y1)
            t1 = RightTwistedTensor(t.twc, x1, yy)
            t2 = RightTwistedTensor(t.twc, xx, y2)
            addmul!(addto, Tensor(t1, t2), coeff*c1*c2*c3)
        end
    end
    addto
end
