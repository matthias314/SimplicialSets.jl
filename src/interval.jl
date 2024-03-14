#
# simplicial interval
#

export IntervalSimplex

"""
    IntervalSimplex <: AbstractSimplex

    IntervalSimplex(p::Int, q::Int)
    IntervalSimplex()

A type representing simplices in a simplicial interval (a 1-simplex).
An `n`-simplex in a simplicial interval is determined by two non-negative
integers `p` and `q` such that `p+q-1` equals `n`.

The first constructor above returns the `p+q-1`-simplex with `p` vertices
equal to `0` and `q` vertices equal to `1`. The second constructor is a
short form for `SimplicialInterval(1, 1)` and returns the unique non-degenerate
`1`-simplex.
"""
struct IntervalSimplex <: AbstractSimplex
    p::Int
    q::Int
    function IntervalSimplex(p::Int, q::Int)
        if p >= 0 && q >= 0 && p+q != 0
            new(p, q)
        else
            error("illegal arguments")
        end
    end
end

IntervalSimplex() = IntervalSimplex(1,1)

@struct_equal_hash IntervalSimplex

dim(x::IntervalSimplex) = x.p + x.q - 1

function d(x::IntervalSimplex, k::Int)
    p, q = x.p, x.q
    if 0 <= k < p
        IntervalSimplex(p-1, q)
    elseif 0 <= k < p+q
        IntervalSimplex(p, q-1)
    else
        error("illegal arguments")
    end
end

function s(x::IntervalSimplex, k::Int)
    p, q = x.p, x.q
    if 0 <= k < p
        IntervalSimplex(p+1, q)
    elseif 0 <= k < p+q
        IntervalSimplex(p, q+1)
    else
        error("illegal arguments")
    end
end

isdegenerate(x::IntervalSimplex, k::Integer) = k != x.p-1

isdegenerate(x::IntervalSimplex) = x.p > 1 || x.q > 1

# not exported at the moment
isabovevertex(x::IntervalSimplex) = x.p == 0 || x.q == 0
