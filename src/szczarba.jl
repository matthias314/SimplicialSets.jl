#
# _Fix1
#

struct _Fix1{F,T}
    f::F
    x::T
end

(f::_Fix1)(xs...; kw...) = f.f(f.x, xs...; kw...)

#
# twisted Cartesian product
#

export LeftTwistedProductSimplex, RightTwistedProductSimplex

abstract type AbstractTwistedProductSimplex <: AbstractSimplex end

struct LeftTwistedProductSimplex{TWF,X,Y} <: AbstractTwistedProductSimplex
    twf::TWF
    x::X
    y::Y
end

struct RightTwistedProductSimplex{TWF,X,Y} <: AbstractTwistedProductSimplex
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

dim(z::AbstractTwistedProductSimplex) = dim(z.x)

s(z::TWP, k) where TWP <: AbstractTwistedProductSimplex = TWP(z.twf, s(z.x, k), s(z.y, k))

function d(z::LeftTwistedProductSimplex, k)
    if k == 0
        LeftTwistedProductSimplex(z.twf, d(z.x, k)⋄inv(z.twf(z.y)), d(z.y, k))
    else
        LeftTwistedProductSimplex(z.twf, d(z.x, k), d(z.y, k))
    end
end

function d(z::RightTwistedProductSimplex, k)
    if k == 0
        RightTwistedProductSimplex(z.twf, d(z.x, k), z.twf(z.x)⋄d(z.y, k))
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
    addmul!(addto, ProductSimplex(d(x, 0)⋄inv(p.twf(y)), d(y, 0)), coeff)
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
    addmul!(addto, ProductSimplex(d(x, 0), p.twf(x)⋄d(y, 0)), coeff)
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
        g = szczarba(twf(x), ii, k) ⋄ g
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
        g = szczarba(twf(x), ii, k) ⋄ g
    end
    yy = inv(g) ⋄ y
    ProductSimplex(xx, yy)
end

szczarba_mul(twf, xg) = xg

szczarba_mul(twf, (x, g), y) = RightTwistedProductSimplex(twf, x, g⋄y)

@linear szsh::SzczarbaShuffle

function (szsh::SzczarbaShuffle{F})(t::Tensor{Tuple{X,Y}};
        coefftype = Int,
        addto = zero(Linear{RightTwistedProductSimplex{F,X,Y},coefftype}),
        coeff = 1) where {F, X <: AbstractSimplex, Y <: AbstractSimplex}
    G = return_type(szsh.twf, X)
    x, y = t
    n = dim(x)
    g = one(G, n)
    foreach_szczarba(n) do ii
        m = has_char2(coefftype) ? Zero() : sum(ii)
        ez(szczarba_shuffle(szsh.twf, ii, x, g), y;
            addto, coeff = withsign(m, coeff), f = _Fix1(szczarba_mul, szsh.twf))
    end
    addto
end
