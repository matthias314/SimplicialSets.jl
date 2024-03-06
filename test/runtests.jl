using Test, StructEqualHash, LinearCombinations, SimplicialSets

using SimplicialSets: d, s, Interval, interval_length
using SimplicialSets.TestHelpers

using LinearCombinations: diff, signed
const TENSORMAP = Tensor

#
# Interval
#

@testset "Interval" begin
    for a in ( (2, 5), 2 => 5, 2:5 )
        @test a isa Interval
	@test interval_length(a) == 4
    end
end

function test_simplex(x::AbstractSimplex, n)
    @test dim(x) isa Integer
    @test n == dim(x) >= 0
    
    @test hash(x) isa UInt
    xc = @inferred copy(x)
    @test x == xc
    @test hash(x) == hash(xc)
    
    @test_throws Exception d(x, -1)
    @test_throws Exception d(x, n+1)
    n == 0 && @test_throws Exception d(x, 0)
    @test_throws Exception s(x, -1)
    @test_throws Exception s(x, n+1)    
    
    for i in 0:n
        if n >= 1
            y = @inferred d(x, i)
            @test dim(y) == n-1 && typeof(y) == typeof(x)
        end
        y = @inferred s(x, i)
        @test dim(y) == n+1 && typeof(y) == typeof(x)
        @test @inferred isdegenerate(y)
        @test @inferred isdegenerate(y, i)
        !isdegenerate(x) && @test all(0:n) do j
            @inferred(isdegenerate(y, j)) == (j == i)
        end
    end
    if n >= 2
        for j in 0:n, i in 0:j-1
           @test d(d(x, j), i) == d(d(x, i), j-1)
        end
    end
    for j in 0:n, i in 0:j
        @test s(s(x, j), i) == s(s(x, i), j+1)
    end
    for j in 0:n, i in 0:j-1
        @test d(s(x, Int8(j)), BigInt(i)) == s(d(x, Int16(i)), Int32(j-1))
    end
    for j in 0:n
        @test d(s(x, j), j) == d(s(x, j), j+1) == x
    end
    for j in 0:n, i in j+2:n+1
        @test d(s(x, Int16(j)), Int8(i)) == s(d(x, Int32(i-1)), BigInt(j))
    end

    # test d(x, kk)
    for R in (Int8, Int, BigInt)
        kk = R[k for k in 0:n if rand(Bool)]
        length(kk) == n+1 && popfirst!(kk)
        @test d(x, kk) == undo_basic(d(BasicSimplex(x), kk))
    end

    # test s(x, kk)
    for R in (Int8, Int, BigInt)
        kk = R[]
        l = 0
        for i in 0:div(n, 3)
            k = rand(l:n+i)
            push!(kk, k)
            l = k+1
        end
        @test s(x, kk) == undo_basic(s(BasicSimplex(x), kk))
    end

    # test r(x, kk)
    @test_throws Exception r(x, [])
    @test_throws Exception r(x, [(1,2)])
    for R in (Int8, Int, BigInt)
         kk = UnitRange{R}[]
         k2 = 0
         while k2 <= (n <= 2 ? n-1 : n-2)
            k1 = rand(k2:n)
            k2 = rand(k1:n)
            push!(kk, R(k1):R(k2))
        end
	if !isempty(kk)
            y = SimplicialSets.r(x, kk)
            @test dim(y) == sum(map(length, kk))-1
            @test y == undo_basic(SimplicialSets.r(BasicSimplex(x), kk))
	end
    end
end

@testset failfast=true "BasicSimplex" begin
    for n in 0:3
        x = SymbolicSimplex('x', n)
	y = BasicSimplex(x)
	@test dim(y) == dim(x)
        n > 0 && @test all(0:n) do k ; d(y, k) == BasicSimplex(d(y.x, k)) end
        @test all(0:n) do k ; s(y, k) == BasicSimplex(s(y.x, k)) end
        test_simplex(y, n)
    end

    x = BarSimplex([Lattice(1,0)]; op = +)
    y = BarSimplex([Lattice(0,1)]; op = +)
    p = LoopGroupSimplex(SymbolicSimplex('p', 3))
    q = LoopGroupSimplex(SymbolicSimplex('q', 3))
    for (x1, x2) in ( (x, y), (p, q) )
        y1 = BasicSimplex(x1)
        y2 = BasicSimplex(x2)
        @test *(y1) == y1
        @test y1 * y2 * y1 == BasicSimplex(x1 * x2 * x1)
        @test inv(y1) == BasicSimplex(inv(x1))
        if x1 isa BarSimplex{AddToMul{Lattice{2}}}
            # / is defined for commutative groups only
            @test y1 / y2 == BasicSimplex(x1 / x2)
	else
	    @test_throws Exception y1 / y2
	end
        @test one(y1) == BasicSimplex(one(x1))
        @test one(y1, 5) == BasicSimplex(one(x1, 5))
        @test one(typeof(y1)) == BasicSimplex(one(typeof(x1)))
        @test one(typeof(y1), 4) == BasicSimplex(one(typeof(x1), 4))
    end
end

function test_group(x::T, is_commutative) where T <: AbstractSimplex
    n = dim(x)
    
    onex = @inferred one(x)
    test_simplex(onex, n)
    @test isone(onex)
    @test one(x, Int8(n+1)) == one(T, BigInt(n+1))
    @test one(T) == one(T, 0)
    @test_throws Exception one(x, -1)
    for i in 0:n
        n > 0 && @test d(onex, i) == one(x, n-1)
        @test s(onex, i) == one(x, n+1)
    end

    @test *(x) == x
    @test x * onex == x == onex * x
    @test isone(x * inv(x)) && isone(inv(x) *x)
    @test_throws Exception x*one(T, n+1)
    
    @test @inferred(x^0) == onex
    @test @inferred(x^1) == x
    @test @inferred(x^3) == x*x*x
    @test @inferred(x^(-1)) == inv(x)
    @test @inferred(x^(-2)) == inv(x)^2
    
    invx = @inferred inv(x)
    @test dim(invx) == n
    for i in 0:n
        n > 0 && @test d(invx, i) == inv(d(x, i))
        @test s(invx, i) == inv(s(x, i))
    end
    
    a = Linear(x => Int8(1), inv(x) => Int8(2))
    @test @inferred(one(a)) == Linear(one(T) => one(Int8))
    @test isone(one(a))
    @test a * one(a) == a == one(a) * a
    if iseven(n)
        @test a*a*a == @inferred a^3
    else
        is_commutative && @test iszero(a*a)
    end
end

function test_group(x::T, y::T, is_commutative) where T <: AbstractSimplex
    k, l = dim(x), dim(y)
    
    if k == l
        xy = @inferred x*y
        @test dim(xy) == k
        for i in 0:k
            k > 0 && @test d(xy, i) == d(x, i)*d(y, i)
            @test s(xy, i) == s(x, i)*s(y, i)
        end
        @test xy*x == x*y*x == x*(y*x)
        if is_commutative
	    @test y*x == xy
            @test x/y == x*inv(y)
	else
            @test_throws Exception x/y
	end
    else
        @test_throws Exception x*y
    end	

    a = Linear{T,BigInt}(x => 2)
    b = Linear{T,Float32}(y => 3)
    ab = @inferred a*b
    @test coefftype(ab) == promote_type(BigInt, Float32) && termtype(ab) == T
    !iszero(ab) && @test deg(ab) == k+l
    
    is_commutative && @test b*a == (-1)^(k*l) * ab
    
    @test diff(ab) == diff(a)*b + (-1)^k*a*diff(b)
end

@testset failfast=true "SymbolicSimplex" begin
    for n in (0, 1, 2, 14)
        x = SymbolicSimplex(:x, BigInt(n))
        test_simplex(x, n)
        y = SymbolicSimplex('y', 1:2:2*n+1)
        test_simplex(y, n)
	v = sort!(rand(UInt8(0):UInt8(31), n+1))
	z = SymbolicSimplex(:z, v)
	test_simplex(z, n)
	@test isdegenerate(z) == !allunique(v)
    end
end

@testset failfast=true "ProductSimplex" begin
    @test_throws Exception ProductSimplex()
    @test_throws Exception ProductSimplex(())
    
    x = SymbolicSimplex('x', 2)
    y = SymbolicSimplex('y', 3)
    @test_throws Exception ProductSimplex(x, y)
    @test_throws Exception ProductSimplex(x; dim = 3)
    @test_throws Exception ProductSimplex(x, y; dim = 3)
    
    for n in 0:3
        xv = ntuple(k -> BasicSimplex(SymbolicSimplex('a'+k, n)), 4)
	for k in 0:4
            @test_throws Exception ProductSimplex(xv[1:k]; dim = -1)
	    w = @inferred ProductSimplex(xv[1:k]; dim = BigInt(n))
	    @test w == @inferred ProductSimplex(xv[1:k]...; dim = n)
	    if k > 0
                v = @inferred ProductSimplex(xv[1:k])
		@test v == w
	    end
            test_simplex(w, n)
	end
    end
end

@testset failfast=true "LoopGroupSimplex" begin
    x = SymbolicSimplex('x', 3)
    xx = LoopGroupSimplex(x)
    xc = LoopGroupSimplex(copy(x))
    @test xc == copy(xx) !== xx
    m = 3
    for n in 1:4
        # xv = ntuple(k -> LoopGroupSimplex(SymbolicSimplex('a'+k, n)), m)
        xv = ntuple(k -> LoopGroupSimplex(BasicSimplex(SymbolicSimplex('a'+k, n))), m)
	u = xv[1]
	for k in 0:m, l in 0:m
	    v = prod(xv[1:k]; init = one(u))
	    w = prod(xv[m-l+1:m]; init = one(u))
	    test_simplex(v, n-1)
	    test_group(v, false)
	    test_group(v, w, false)
        end	
    end
end

function test_barsimplex(n, m, groupsimplex, is_commutative)
    for k in 0:n
        T = typeof(groupsimplex(0, m))
	v = T[groupsimplex(i-1, m) for i in 1:k]
        x = @inferred BarSimplex(v)
	xc = BarSimplex(copy(v))
        @test xc == copy(x) !== x

	test_simplex(x, k)
	
	is_commutative || continue
	
        test_group(x, true)
	for l in 0:n
            y = BarSimplex(T[groupsimplex(i-1, m) for i in 1:l])
	    test_group(x, y, true)
	end
    end
end

# random_lattice_element(_, m) = AddToMul(Lattice(Tuple(rand(-3:3, m))))
# we need an inferrable return type
# random_lattice_element(_, _)::AddToMul{Lattice{3}} = AddToMul(Lattice(ntuple(_ -> rand(-3:3), 3)))
random_lattice_element(_, _) = AddToMul(Lattice(ntuple(_ -> rand(-3:3), 3)))

struct M
    a::Matrix{Rational{Int}}
end

@struct_equal_hash M

Base.one(::Type{M}) = M([1 0; 0 1])
Base.one(::M) = one(M)
Base.isone(x::M) = isone(x.a)
Base.inv(x::M) = M(inv(x.a))
Base.:*(x::M, ys::M...) = M(*(x.a, map(y -> y.a, ys)...))
# Base.:/(x::M, y::M) = x*inv(y)

function random_matrix_2x2(_, _)
    a11 = rand(1:8)
    a22 = rand(1:8)
    b = a11*a22+1
    a12 = findfirst(i -> rem(b, i) == 0, 2:b) + 1
    a21 = div(b, a12)
    M([a11 a12; a21 a22])
end

@testset failfast=true "BarSimplex discrete" begin
    test_barsimplex(4, 3, random_lattice_element, true)
    test_barsimplex(4, 3, random_matrix_2x2, false)
end

function random_loopgroupsimplex(n, m)
    prod(LoopGroupSimplex(SymbolicSimplex(rand('a':'z'), n+1)) for _ in 1:m)
end

function random_barsimplex(n, m)
    BarSimplex([random_lattice_element(0, m) for i in 1:n])
end

function random_barbarsimplex(n, m)
    BarSimplex([random_barsimplex(i-1, m) for i in 1:n])
end

@testset failfast=true "BarSimplex simplicial" begin
    test_barsimplex(4, 3, random_barsimplex, true)      # test double bar construction of Z^3
    test_barsimplex(3, 3, random_barbarsimplex, true)   # test triple bar construction of Z^3
    test_barsimplex(4, 3, BasicSimplex ∘ random_loopgroupsimplex, false)
    # it seems that the multiplication of chains on the bar construction of a loop group is graded commutative!
end

#=
q2(n) = div(n*(n+1), 2)

qsign(a::Linear) = Linear(x => signed(q2(deg(x)), c) for (x, c) in a)

opp_linear(a::Linear) = qsign(opposite(a))

opp_swap(a::Linear{<:ProductSimplex}) =
    # deg(x) == degree of ProductSimplex(x, y)
    Linear(opposite(ProductSimplex(y, x)) => signed(q2(deg(x)), c) for ((x, y), c) in a)

opp_swap(a::Linear{<:Tensor}) =
    Linear(opposite(Tensor(y, x)) => signed(q2(deg(x)+deg(y)), c) for ((x, y), c) in a)
=#

const opposite_swap = opposite ∘ swap

@testset failfast=true "OppositeSimplex" begin
    for n in (0, 1, 4)
        x = SymbolicSimplex(:x, n)
        test_simplex(x, n)
	test_simplex(BasicSimplex(x), n)
	
	x = random_loopgroupsimplex(n, 2)
	y = random_loopgroupsimplex(n, 2)
        test_simplex(x, n)
	test_group(x, false)
	test_group(x, y, false)
	
        x = random_barsimplex(n, 2)
        y = random_barsimplex(n, 2)
        test_simplex(x, n)
	test_group(x, true)
	test_group(x, y, true)
    end
    
    for n in 0:8
        x = SymbolicSimplex(:x, n)
        a = Linear(x => 1)
	@test opposite(diff(a)) == diff(opposite(a))
    end
end

@testset failfast=true "ez" begin
    @test @inferred(ez(Tensor())) == Linear(ProductSimplex(; dim = 0) => 1)

    for k in 0:8
        t = Tensor(ntuple(Returns(SymbolicSimplex('i', 1)), k))
	@inferred ez(t)
	@inferred ez(t; coefftype = Val(Float16))
    end

    for m in 0:4, n in 0:4
        x = SymbolicSimplex(:x, m)
        y = SymbolicSimplex(:y, n)
        a = tensor(x, y)
        c1 = a |> ez |> opposite_swap
        c2 = a |> opposite_swap |> ez
        @test c1 == c2
    end
    
    n = 3
    
    x = random_loopgroupsimplex(n, 2)
    y = random_loopgroupsimplex(n, 2)
    z = random_loopgroupsimplex(n, 2)
    
    t = Tensor(BasicSimplex.((x, y, z)))
    @test ez(undo_basic(t)) == undo_basic(ez(t))
    
    rgr = regroup( :( (1,(2,3)) ), :( (1,2,3) ) )
    rgl = regroup( :( ((1,2),3) ), :( (1,2,3) ) )
    @test rgr(ez(x, ez(y, z))) == ez(x, y, z) == rgl(ez(ez(x, y), z))
    
    w = SymbolicSimplex('w', 0)
    @test @inferred(ez(Tensor(w, w))) == Linear(ProductSimplex(w, w) => 1)
    @test ez(x, w) == Linear(ProductSimplex(x, s(w, 0:n-1)) => 1)
    @test ez(w, x) == Linear(ProductSimplex(s(w, 0:n-1), x) => 1)
    
    a = Linear(Tensor(x,y,z) => 1)
    @test diff(ez(a)) == ez(diff(a))
end

@testset failfast=true "aw" begin
    @test @inferred(aw(ProductSimplex(; dim = 0))) == Linear(Tensor() => 1)
    @test iszero(aw(ProductSimplex(; dim = 1)))
    
    for k in 0:8
        w = ProductSimplex(ntuple(Returns(SymbolicSimplex('i', 1)), k); dim = 1)
	@inferred aw(w)
	@inferred aw(w; coefftype = Val(Float32))
    end

    for n in 0:8
        x = SymbolicSimplex('x', n)
        y = SymbolicSimplex('y', n)
        w = ProductSimplex(x, y)
        b = Linear(w => 1)
    
        c1 = b |> aw |> opposite_swap
        c2 = b |> opposite_swap |> aw
        @test c1 == c2
    end

    n = 3

    x = random_loopgroupsimplex(n, 2)
    y = random_loopgroupsimplex(n, 2)
    z = random_loopgroupsimplex(n, 2)

    w = ProductSimplex(x, y, z)
    w = ProductSimplex(BasicSimplex.((x, y, z)))
    @test aw(undo_basic(w)) == undo_basic(aw(w))

    rgr, rgri = regroup_inv( :( (1,(2,3)) ), :( (1,2,3) ) )
    rgl, rgli = regroup_inv( :( ((1,2),3) ), :( (1,2,3) ) )
    a = w |> aw
    ar = w |> rgri |> aw |> TENSORMAP(identity, aw) |> rgr
    al = w |> rgli |> aw |> TENSORMAP(aw, identity) |> rgl
    @test ar == a == al
    
    w = SymbolicSimplex('w', 0)
    @test @inferred(aw(ProductSimplex(w, w))) == Linear(Tensor(w, w) => 1)
    @test aw(ProductSimplex(x, s(w, 0:n-1))) == Linear(Tensor(x, w) => 1)
    @test aw(ProductSimplex(s(w, 0:n-1), x)) == Linear(Tensor(w, x) => 1)

    a = Linear(ProductSimplex(x,y,z) => 1)
    @test diff(aw(a)) == aw(diff(a))
end

@testset failfast=true "shih" begin

    for n in 0:8
        x = SymbolicSimplex('x', n)
        y = SymbolicSimplex('y', n)
        w = ProductSimplex(x, y)
        @inferred shih_eml(w; coefftype = Val(Int16))
        @inferred shih_opp(w; coefftype = Val(Int32))

        b = Linear(w => 1)
    
        c1 = b |> shih_eml |> opposite_swap
        c2 = b |> opposite_swap |> shih_opp
        @test c1 == c2
    end

    n = 3

    x = random_loopgroupsimplex(n, 2)
    y = random_loopgroupsimplex(n, 2)
    z = random_loopgroupsimplex(n, 2)

    w = ProductSimplex(x, y)
    b = Linear(w => 1)

    c = Linear(ProductSimplex(BasicSimplex(x), BasicSimplex(y)) => 1)
    u = SymbolicSimplex('u', 0)
    v = ProductSimplex(x, y, z)
    a = Linear(v => 1)

    rgr, rgri = regroup_inv( :( (1,(2,3)) ), :( (1,2,3) ) )
    rgl, rgli = regroup_inv( :( ((1,2),3) ), :( (1,2,3) ) )

    for shih in (shih_eml, shih_opp)
        @test shih(undo_basic(c)) == undo_basic(shih(b))

        @test iszero(shih(shih(b)))
	
	@test iszero(shih(ProductSimplex(u, u)))

        c1 = a |> rgri |> shih |> rgr |> rgli |> shih |> rgl
        c2 = a |> rgli |> shih |> rgl |> rgri |> shih |> rgr
	@test c1 == -c2
    end
end

@testset failfast=true "EZ relations" begin
    for n in (0, 1, 4)
        x = SymbolicSimplex('x', n)
        y = SymbolicSimplex('y', n)
        a = tensor(x, y)
        w = ProductSimplex(x, y)
        b = Linear(w => 1)

        @test aw(ez(a)) == a
        for shih in (shih_eml, shih_opp)
            @test ez(aw(b)) - b == diff(shih(b)) + shih(diff(b))
        end
    end
end

@testset failfast=true "Hirsch formula" begin
    f11, g21 = regroup_inv( :( (1,2,3) ), :( (1,(2,3)) ) )
    f31, g11 = regroup_inv( :( (1,2,3) ), :( ((1,2),3) ) )
    f21 = regroup( :( (1,2,3) ), :( (2,(1,3)) ) )
    g31 = regroup( :( ((2,1),3) ), :( (2,3,1) ) )

    for n in (0, 1, 2, 4)
        x = SymbolicSimplex('x', n)
        y = SymbolicSimplex('y', n)
        z = SymbolicSimplex('z', n)
        a = Linear( ProductSimplex(x, y, z) => 1 )

        u = a |> f11 |> shih |> swap |> g11 |> aw
        v = a |> f21 |> aw |> TENSORMAP(identity,  aw ∘ swap ∘ shih) |> g21
        w = a |> f31 |> aw |> TENSORMAP(aw ∘ swap ∘ shih, identity) |> g31

        @test u == v + w
    end
end

@testset failfast=true "AAFR formula" begin
# we check a formula from
# Alvarez, V.; Armario, J. A.; Frau, M. D.; Real, P.
# Algebra structures on the twisted Eilenberg-Zilber theorem

# TODO: what is the formula?

    n = 3

    x = SymbolicSimplex('x', n)
    y = SymbolicSimplex('y', n)
    z = SymbolicSimplex('z', n)
    w = SymbolicSimplex('w', n)
    xy = ProductSimplex(x, y)
    zw = ProductSimplex(z, w)

    a = Linear(xy => 1)
    b = Linear(zw => 1)

    f = regroup( :( ((1,2),(3,4)) ), :( ((1,3),(2,4)) ) )

    t1 = ez(x, y)
    t2 = shih(b)
    t3 = ez(tensor(t1, t2))
    t4 = f(t3)
    t5 = shih(t4)
    @test iszero(t5)

    t1 = shih(a)
    t2 = ez(z, w)
    t3 = ez(tensor(t1, t2))
    t4 = f(t3)
    t5 = shih(t4)
    @test iszero(t5)

    t1 = shih(a)
    t2 = shih(b)
    t3 = ez(tensor(t1, t2))
    t4 = f(t3)
    t5 = shih(t4)
    @test iszero(t5)
end

@testset failfast=true "surjections" begin
    surj = @inferred Surjection{0}(Int[])
    @test @inferred arity(surj) == 0
    @test @inferred deg(surj) == 0
    
    surj = @inferred Surjection{4}(Int[1,2,3,2,2,1,4])
    @test @inferred arity(surj) == 4
    @test @inferred deg(surj) == 3
    @test @inferred isdegenerate(surj)
    @test iszero(Linear(surj => 1))

    surj = Surjection([1,3,2,1,4,2,1])
    a = Linear(surj => 1)
    b = @inferred diff(a)
    @test typeof(b) == typeof(a)
    @test iszero(diff(b))
    
    c = @inferred diff(a; coeff = 2)
    @test typeof(c) == typeof(a)
    @test c == 2*b
    
    c2 = @inferred diff(a; addto = c, coeff = -2)
    @test c2 === c
    @test iszero(c2)
end

@testset failfast=true "interval cuts" begin
    surj = Surjection(Int[])
    y = SymbolicSimplex('y', 0)
    @test surj(y; coeff = -1) == Linear(Tensor() => -1)    

    s12 = Surjection(1:2)
    s21 = Surjection([2,1])
    f4 = TENSORMAP(coprod, coprod) ∘ coprod
    rg = regroup(:(((1,2),(3,4))), :((1,2,3,4)))

    y = SymbolicSimplex('y', 4)
    z = random_loopgroupsimplex(4, 3)
    w = random_barsimplex(4, 3)
    for x in (BasicSimplex(y), y, ProductSimplex(y, y), z, w)
        surj = Surjection(Int[])	
        @test iszero(surj(x))

        cx = coprod(x)
        @test s12(x) == cx
        @test s21(x) == swap(cx)
    
        s1234 = Surjection(1:4)
        @test rg(f4(x)) == s1234(x)
    
        surj = Surjection([1,3,3,2])
        @test iszero(surj(x))
    
        for n in 0:8
            surj = Surjection(1:n)
	    @inferred surj(x)
        end
    end
end
