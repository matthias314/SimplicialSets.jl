#
# simplicial interval
#

export IntervalSimplex

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
