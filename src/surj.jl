#
# surjections and interval cut operations
#

export Surjection, arity

"""
    Surjection{K}

    Surjection{K}(u [; check = true]) where K
    Surjection(u)

A type representing an element of arity `K` in the surjection operad.

The first constructor return the surjection given by `u` TODO. If the optional
keyword argument is set to `false`, then it is not checked that `u` is indeed a surjection.
The second constructor is equivalent to the first with `K` set to `maximum(u)`.

See also [`arity`](@ref).

# Examples
```jldoctest
julia> Surjection([1, 2, 3, 1])
Surjection{3}([1, 2, 3, 1])

julia> Surjection(Int[])
Surjection{0}(Int64[])
```
"""
struct Surjection{K}          # K is the number of labels
    u::Vector{Int}            # the surjection proper
    # u1::Vector{Int}           # previous interval with same label (0 for first occurrence)
    f::Vector{Bool}           # final intervals
    v::NTuple{K,Vector{Int}}  # intervals for each i in 1:K
end

show(io::IO, surj::Surjection{K}) where K = print(io, "Surjection{$K}($(repr(surj.u)))")

==(surj1::Surjection, surj2::Surjection) = surj1.u == surj2.u

hash(surj::Surjection, h::UInt) = hash(surj.u, h)

deg(surj::Surjection{K}) where K = length(surj.u)-K

"""
    arity(surj::Surjection{K}) where K -> K

Return the arity of the surjection `surj`.
"""
arity(surj::Surjection{K}) where K = K

is_surjection(k::Integer, u::AbstractVector{<:Integer}) = extrema(u; init = (1, 0)) == (1, k) && all(in(u), 2:k-1)

isdegenerate_surjection(u) = any(i -> @inbounds(u[i-1] == u[i]), 2:length(u))

isdegenerate(surj::Surjection) = isdegenerate_surjection(surj.u)

linear_filter(surj::Surjection) = !isdegenerate(surj)

# TODO: do we need to allow empty u? yes
function Surjection{k}(u; check = true) where k
    !check || is_surjection(k, u) || error("argument is not a surjection onto 1:$k")

    l = length(u)

    v = ntuple(_ -> Int[], k)   # we cannot use Returns
    for i in 1:l
        push!(v[u[i]], i)
    end

    # determine final intervals
    f = fill(false, l)
    for i in 1:k
        f[v[i][end]] = true
    end

    Surjection{k}(u, f, v)
end

Surjection(u) = Surjection{maximum(u; init = 0)}(u)

"""
    diff(surj::Surjection{K}) where K -> Linear{Surjection{K}}

Return the differential (or boundary) of `u` in the surjection operad.

This function is linear and supports the keyword arguments `coefftype`, `addto`,
`coeff` and `is_filtered` as described for the macro `@linear`.

See also `LinearCombinations.@linear`.
"""
@linear_kw function diff(surj::Surjection{k};
        coefftype = Int,
        addto = zero(Linear{Surjection{k},unval(coefftype)}),
        coeff = ONE,
        is_filtered = false) where k
    (; u, f, v) = surj
    s = 0
    l = length(u)
    us = Vector{Int}(undef, l)
    @inbounds for i in 1:l
        if !f[i]
            si = s
            s += 1
        else
            vi = v[u[i]]
            length(vi) == 1 && continue
            si = us[vi[end-1]] + 1
        end
        us[i] = si
        if i == 1 || i == l || u[i-1] != u[i+1]
            # this means that w below is non-degenerate
            w = deleteat!(copy(u), i)
            addmul!(addto, Surjection{k}(w; check = false), signed(si, coeff); is_filtered = true)
        end
    end
    addto
end

# interval cuts

struct IC
    n::Int
    k::Int
end

length(a::IC) = binomial(a.n+a.k-1, a.k-1)

function iterate(a::IC)
    ic = zeros(Int, a.k+1)
    @inbounds ic[a.k+1] = a.n
    ic, ic
end

@inline function iterate(a::IC, ic)
    i = a.k
    while i > 1 && @inbounds ic[i] == a.n
        i -= 1
    end
    if i == 1
        nothing
    else
        @inbounds m = ic[i]+1
        for j in i:a.k
            @inbounds ic[j] = m
        end
        ic, ic
    end
end

@propagate_inbounds function isdegenerate_ic(v, ic)
    for pp in v
        i = -1
        for p in pp
            if p != 1 && ic[p] == i
                return true
            end
            i = ic[p+1]
        end
    end
    false
end

@linear_kw function (surj::Surjection{k})(x::T;
        coefftype = Int,
        addto = zero(Linear{Tensor{NTuple{k,T}},unval(coefftype)}),
        coeff = one(coeff_type(addto)),
        is_filtered = false) where {k, T <: AbstractSimplex}

    iszero(coeff) && return addto

    if k == 0
    # in this case the interval cut operation is the augmentation map
        deg(x) == 0 && addmul!(addto, Tensor(), coeff)
        return addto
    end

    u = surj.u
    f = surj.f
    v = surj.v
    l = length(u)

    a = IC(dim(x), l)
    sizehint!(addto, length(addto)+length(a))

    rr = ntuple(i -> Vector{Tuple{Int,Int}}(undef, length(v[i])), k)
    L = Vector{Int}(undef, l)

    @inbounds for ic in a
        isdegenerate_ic(v, ic) && continue

        for i in 1:k, j in 1:length(v[i])
            rr[i][j] = (ic[v[i][j]], ic[v[i][j]+1])
        end

        xl = map(Fix1(r, x), rr)
        any(isdegenerate, xl) && continue
        # TODO: incorporate non-deg testing into iterate

        # compute permutation sign
        s = sum(@inbounds ifelse(f[i], 0, ic[i+1]) for i in 1:l)
        for i in 1:l
            L[i] = Li = ic[i+1] - ic[i] + ifelse(f[i], 0, 1)
            ui = u[i]
            for j in 1:i-1
                if ui < u[j]
                    s += Li*L[j]
                end
            end
        end

        addmul!(addto, Tensor(xl), signed(s, coeff); is_filtered = true)
    end

    addto
end

@linear surj::Surjection
