#
# Szczarba operators
#

export szczarba

szczarba(x::AbstractSimplex, ::Tuple{}, ::Int, ::Int) = x

function szczarba(x::AbstractSimplex, ii::Tuple, k::Int, l::Int = 0)
# l keeps track of how often the simplicial operator has been derived
    # println("x = $x  ii = $ii  k = $k  l = $l")
    i1 = first(ii)
    i2 = Base.tail(ii)
    if k < i1
        szczarba(s(d(x, i1-k+l), l), i2, k, l+1)
    elseif k == i1
        szczarba(x, i2, k, l+1)
    else
        szczarba(s(x, l), i2, k-1, l+1)
    end
end

function szczarba(x::AbstractSimplex, twf, ii::Tuple)
    n = dim(x)
    if x == 0
        error("illegal arguments")
    end
    g = szczarba(inv(twf(x)), ii, 0)
    for k in 1:n-1
        x = d(x, 0)
        g *= szczarba(inv(twf(x)), ii, k)
    end
    g
end

# TODO: sign
function foreach_szczarba(f, n::Int)
    ii = zeros(Int, n)
    while true
        f(ii)
        k = n
        while ii[k] == k-1
            ii[k] = 0
            k -= 1
            if k == 0
                return
            end
        end
        ii[k] += 1
    end
end

function szczarba(x::AbstractSimplex, twf)
    n = dim(x)
    if n > 1
        a = TODO
        foreach_szczarba(n) do ii, ss
            addcoeff!(a, szczarba(x, twf, ii), ss)
        end
        a
    elseif n == 1
        g = inv(twf(x))
        Linear(g => 1, one(g) => -1)
    else
        error("illegal arguments")
    end
end
