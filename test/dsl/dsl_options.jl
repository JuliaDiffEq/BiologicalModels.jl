#! format: off

### Fetch Packages and Set Global Variables ###
using Catalyst, ModelingToolkit, OrdinaryDiffEq, Plots

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)

# Sets globally used variable
@variables t

### Run Tests ###

# Test creating networks with/without options.
let
    @reaction_network begin (k1, k2), A <--> B end
    @reaction_network begin
        @parameters k1 k2
        (k1, k2), A <--> B
    end
    @reaction_network begin
        @parameters k1 k2
        @species A(t) B(t)
        (k1, k2), A <--> B
    end
    @reaction_network begin
        @species A(t) B(t)
        (k1, k2), A <--> B
    end

    @reaction_network begin
        @parameters begin
            k1
            k2
        end
        (k1, k2), A <--> B
    end
    @reaction_network begin
        @species begin
            A(t)
            B(t)
        end
        (k1, k2), A <--> B
    end
    @reaction_network begin
        @parameters begin
            k1
            k2
        end
        @species begin
            A(t)
            B(t)
        end
        (k1, k2), A <--> B
    end

    n1 = @reaction_network name begin (k1, k2), A <--> B end
    n2 = @reaction_network name begin
        @parameters k1 k2
        (k1, k2), A <--> B
    end
    n3 = @reaction_network name begin
        @species A(t) B(t)
        (k1, k2), A <--> B
    end
    n4 = @reaction_network name begin
        @parameters k1 k2
        @species A(t) B(t)
        (k1, k2), A <--> B
    end
    n5 = @reaction_network name begin
        (k1, k2), A <--> B
        @parameters k1 k2
    end
    n6 = @reaction_network name begin
        (k1, k2), A <--> B
        @species A(t) B(t)
    end
    n7 = @reaction_network name begin
        (k1, k2), A <--> B
        @parameters k1 k2
        @species A(t) B(t)
    end
    n8 = @reaction_network name begin
        @parameters begin
            k1
            k2
        end
        (k1, k2), A <--> B
    end
    n9 = @reaction_network name begin
        @species begin
            A(t)
            B(t)
        end
        (k1, k2), A <--> B
    end
    n10 = @reaction_network name begin
        @parameters begin
            k1
            k2
        end
        @species begin
            A(t)
            B(t)
        end
        (k1, k2), A <--> B
    end
    @test all(==(n1), (n2, n3, n4, n5, n6, n7, n8, n9, n10))
end

# Tests that when either @species or @parameters is given, the other is infered properly. 
let
    rn1 = @reaction_network begin
        k*X, A + B --> 0
    end
    @test issetequal(species(rn1), @species A(t) B(t))
    @test issetequal(parameters(rn1), @parameters k X)

    rn2 = @reaction_network begin
        @species A(t) B(t) X(t)
        k*X, A + B --> 0
    end
    @test issetequal(species(rn2), @species A(t) B(t) X(t))
    @test issetequal(parameters(rn2), @parameters k)

    rn3 = @reaction_network begin
        @parameters k
        k*X, A + B --> 0
    end
    @test issetequal(species(rn3), @species A(t) B(t))
    @test issetequal(parameters(rn3), @parameters k X)

    rn4 = @reaction_network begin
        @species A(t) B(t) X(t)
        @parameters k
        k*X, A + B --> 0
    end
    @test issetequal(species(rn4), @species A(t) B(t) X(t))
    @test issetequal(parameters(rn4), @parameters k)

    rn5 = @reaction_network begin
        @parameters k B [isconstantspecies=true]
        k*X, A + B --> 0
    end
    @test issetequal(species(rn5), @species A(t))
    @test issetequal(parameters(rn5), @parameters k B X)
end

# Test inferring with stoichiometry symbols and interpolation.
let
    @parameters k g h gg X y [isconstantspecies = true]
    t = Catalyst.DEFAULT_IV
    @species A(t) B(t) BB(t) C(t)

    rni = @reaction_network inferred begin
        $k*X, $y + g*A + h*($gg)*B + $BB * C --> k*C
    end
    @test issetequal(species(rni), [A, B, BB, C])
    @test issetequal(parameters(rni), [k, g, h, gg, X, y])

    rnii = @reaction_network inferred begin
        @species BB(t)
        @parameters y [isconstantspecies = true]
        k*X, y + g*A + h*($gg)*B + BB * C --> k*C
    end
    @test rnii == rni
end

# Tests that when some species or parameters are left out, the others are set properly.
let
    rn6 = @reaction_network begin
        @species A(t)
        k*X, A + B --> 0
    end
    @test issetequal(species(rn6), @species A(t) B(t))
    @test issetequal(parameters(rn6), @parameters k X)

    rn7 = @reaction_network begin
        @species A(t) X(t)
        k*X, A + B --> 0
    end
    @test issetequal(species(rn7), @species A(t) X(t) B(t))
    @test issetequal(parameters(rn7), @parameters k)

    rn7 = @reaction_network begin
        @parameters B [isconstantspecies=true]
        k*X, A + B --> 0
    end
    @test issetequal(species(rn7), @species A(t))
    @test issetequal(parameters(rn7), @parameters B k X)

    rn8 = @reaction_network begin
        @parameters B [isconstantspecies=true] k
        k*X, A + B --> 0
    end
    @test issetequal(species(rn8), @species A(t))
    @test issetequal(parameters(rn8), @parameters B k X)

    rn9 = @reaction_network begin
        @parameters k1 X1
        @species A1(t) B1(t)
        k1*X1, A1 + B1 --> 0
        k2*X2, A2 + B2 --> 0
    end
    @test issetequal(species(rn9), @species A1(t) B1(t) A2(t) B2(t))
    @test issetequal(parameters(rn9), @parameters k1 X1 k2 X2)

    rn10 = @reaction_network begin
        @parameters k1 X2 B2 [isconstantspecies=true]
        @species A1(t) X1(t)
        k1*X1, A1 + B1 --> 0
        k2*X2, A2 + B2 --> 0
    end
    @test issetequal(species(rn10), @species A1(t) X1(t) B1(t) A2(t))
    @test issetequal(parameters(rn10), @parameters k1 X2 B2 k2)

    rn11 = @reaction_network begin
        @parameters k1 k2
        @species X1(t)
        k1*X1, A1 + B1 --> 0
        k2*X2, A2 + B2 --> 0
    end
    @test issetequal(species(rn11), @species X1(t) A1(t) A2(t) B1(t) B2(t))
    @test issetequal(parameters(rn11), @parameters k1 k2 X2)
end

##Checks that some created networks are identical.
let
    rn12 = @reaction_network name begin (k1, k2), A <--> B end
    rn13 = @reaction_network name begin
        @parameters k1 k2
        (k1, k2), A <--> B
    end
    rn14 = @reaction_network name begin
        @species A(t) B(t)
        (k1, k2), A <--> B
    end
    rn15 = @reaction_network name begin
        @parameters k1 k2
        @species A(t) B(t)
        (k1, k2), A <--> B
    end
    @test all(==(rn12), (rn13, rn14, rn15))
end

# Checks that the rights things are put in vectors. 
let
    rn18 = @reaction_network name begin
        @parameters p d1 d2
        @species A(t) B(t)
        p, 0 --> A
        1, A --> B
        (d1, d2), (A, B) --> 0
    end
    rn19 = @reaction_network name begin
        p, 0 --> A
        1, A --> B
        (d1, d2), (A, B) --> 0
    end
    @test rn18 == rn19

    @parameters p d1 d2
    @species A(t) B(t)
    @test isequal(parameters(rn18)[1], p)
    @test isequal(parameters(rn18)[2], d1)
    @test isequal(parameters(rn18)[3], d2)
    @test isequal(species(rn18)[1], A)
    @test isequal(species(rn18)[2], B)

    rn20 = @reaction_network name begin
        @species X(t)
        @parameters S
        mm(X,v,K), 0 --> Y
        (k1,k2), 2Y <--> Y2
        d*Y, S*(Y2+Y) --> 0
    end
    rn21 = @reaction_network name begin
        @species X(t) Y(t) Y2(t)
        @parameters v K k1 k2 d S
        mm(X,v,K), 0 --> Y
        (k1,k2), 2Y <--> Y2
        d*Y, S*(Y2+Y) --> 0
    end
    rn22 = @reaction_network name begin
        @species X(t) Y2(t)
        @parameters d k1
        mm(X,v,K), 0 --> Y
        (k1,k2), 2Y <--> Y2
        d*Y, S*(Y2+Y) --> 0
    end
    @test all(==(rn20), (rn21, rn22))
    @parameters v K k1 k2 d S
    @species X(t) Y(t) Y2(t)
    @test issetequal(parameters(rn22),[v K k1 k2 d S])
    @test issetequal(species(rn22), [X Y Y2])
end

# Tests that defaults work. 
let
    rn26 = @reaction_network name begin
        @parameters p=1.0 d1 d2=5
        @species A(t) B(t)=4
        p, 0 --> A
        1, A --> B
        (d1, d2), (A, B) --> 0
    end

    rn27 = @reaction_network name begin
    @parameters p1=1.0 p2=2.0 k1=4.0 k2=5.0 v=8.0 K=9.0 n=3 d=10.0
    @species X(t)=4.0 Y(t)=3.0 X2Y(t)=2.0 Z(t)=1.0
        (p1,p2), 0 --> (X,Y)
        (k1,k2), 2X + Y --> X2Y
        hill(X2Y,v,K,n), 0 --> Z
        d, (X,Y,X2Y,Z) --> 0
    end
    u0_27 = []
    p_27 = []

    rn28 = @reaction_network name begin
    @parameters p1=1.0 p2 k1=4.0 k2 v=8.0 K n=3 d
    @species X(t)=4.0 Y(t) X2Y(t) Z(t)=1.0
        (p1,p2), 0 --> (X,Y)
        (k1,k2), 2X + Y --> X2Y
        hill(X2Y,v,K,n), 0 --> Z
        d, (X,Y,X2Y,Z) --> 0
    end
    u0_28 = symmap_to_varmap(rn28, [:p2=>2.0, :k2=>5.0, :K=>9.0, :d=>10.0])
    p_28 = symmap_to_varmap(rn28, [:Y=>3.0, :X2Y=>2.0])
    defs28 = Dict(Iterators.flatten((u0_28, p_28)))

    rn29 = @reaction_network name begin
    @parameters p1 p2 k1 k2 v K n d
    @species X(t) Y(t) X2Y(t) Z(t)
        (p1,p2), 0 --> (X,Y)
        (k1,k2), 2X + Y --> X2Y
        hill(X2Y,v,K,n), 0 --> Z
        d, (X,Y,X2Y,Z) --> 0
    end
    u0_29 = symmap_to_varmap(rn29, [:p1=>1.0, :p2=>2.0, :k1=>4.0, :k2=>5.0, :v=>8.0, :K=>9.0, :n=>3, :d=>10.0])
    p_29 = symmap_to_varmap(rn29, [:X=>4.0, :Y=>3.0, :X2Y=>2.0, :Z=>1.0])
    defs29 = Dict(Iterators.flatten((u0_29, p_29)))

    @test ModelingToolkit.defaults(rn27) == defs29
    @test merge(ModelingToolkit.defaults(rn28), defs28) == ModelingToolkit.defaults(rn27)
end

### Observables ###

# Test basic functionality.
# Tests various types of indexing.
let 
    rn = @reaction_network begin
        @observables begin
            X ~ Xi + Xa
            Y ~ Y1 + Y2
        end
        (p,d), 0 <--> Xi
        (k1,k2), Xi <--> Xa
        (k3,k4), Y1 <--> Y2
    end
    @unpack X, Xi, Xa, Y, Y1, Y2, p, d, k1, k2, k3, k4 = rn

    # Test that ReactionSystem have the correct properties.
    @test length(species(rn)) == 4
    @test length(states(rn)) == 4
    @test length(observed(rn)) == 2
    @test length(equations(rn)) == 6

    @test isequal(observed(rn)[1], X ~ Xi + Xa)
    @test isequal(observed(rn)[2], Y ~ Y1 + Y2)

    # Tests correct indexing of solution.
    u0 = [Xi => 0.0, Xa => 0.0, Y1 => 1.0, Y2 => 2.0]
    ps = [p => 1.0, d => 0.2, k1 => 1.5, k2 => 1.5, k3 => 5.0, k4 => 5.0]

    oprob = ODEProblem(rn, u0, (0.0, 1000.0), ps)
    sol = solve(oprob, Tsit5())
    @test sol[X][end] ≈ 10.0
    @test sol[Y][end] ≈ 3.0
    @test sol[rn.X][end] ≈ 10.0
    @test sol[rn.Y][end] ≈ 3.0
    @test sol[:X][end] ≈ 10.0
    @test sol[:Y][end] ≈ 3.0

    # Tests that observables can be used for plot indexing.
    @test plot(sol; idxs=X).series_list[1].plotattributes[:y][end] ≈ 10.0
    @test plot(sol; idxs=rn.X).series_list[1].plotattributes[:y][end] ≈ 10.0
    @test plot(sol; idxs=:X).series_list[1].plotattributes[:y][end] ≈ 10.0
    @test plot(sol; idxs=[X, Y]).series_list[2].plotattributes[:y][end] ≈ 3.0
    @test plot(sol; idxs=[rn.X, rn.Y]).series_list[2].plotattributes[:y][end] ≈ 3.0
    @test plot(sol; idxs=[:X, :Y]).series_list[2].plotattributes[:y][end] ≈ 3.0
end

# Compares programmatic and DSL system with observables.
let
    # Model declarations.
    rn_dsl = @reaction_network begin
        @observables begin
            X ~ x + 2x2y
            Y ~ y + x2y
        end
        k, 0 --> (x, y)
        (kB, kD), 2x + y <--> x2y
        d, (x,y,x2y) --> 0
    end

    @variables t X(t) Y(t)
    @species x(t), y(t), x2y(t)
    @parameters k kB kD d
    r1 = Reaction(k, nothing, [x], nothing, [1])
    r2 = Reaction(k, nothing, [y], nothing, [1])
    r3 = Reaction(kB, [x, y], [x2y], [2, 1], [1])
    r4 = Reaction(kD, [x2y], [x, y], [1], [2, 1])
    r5 = Reaction(d, [x], nothing, [1], nothing)
    r6 = Reaction(d, [y], nothing, [1], nothing)
    r7 = Reaction(d, [x2y], nothing, [1], nothing)
    obs_eqs = [X ~ x + 2x2y, Y ~ y + x2y]
    @named rn_prog = ReactionSystem([r1, r2, r3, r4, r5, r6, r7], t, [x, y, x2y], [k, kB, kD, d]; observed = obs_eqs)

    # Make simulations.
    u0 = [x => 1.0, y => 0.5, x2y => 0.0]
    tspan = (0.0, 15.0)
    ps = [k => 1.0, kD => 0.1, kB => 0.5, d => 5.0]
    oprob_dsl = ODEProblem(rn_dsl, u0, tspan, ps)
    oprob_prog = ODEProblem(rn_prog, u0, tspan, ps)

    sol_dsl = solve(oprob_dsl, Tsit5(); saveat=0.1)
    sol_prog = solve(oprob_prog, Tsit5(); saveat=0.1)

    # Tests observables equal in both cases.
    @test oprob_dsl[:X] == oprob_prog[:X]
    @test oprob_dsl[:Y] == oprob_prog[:Y]
    @test sol_dsl[:X] == sol_prog[:X]
    @test sol_dsl[:Y] == sol_prog[:Y]
end

# Tests for complicated observable formula.
# Tests using a single observable (without begin/end statement).
# Tests using observable component not part of reaction.
# Tests using parameters in observables formula.
let 
    rn = @reaction_network begin
        @parameters op_1 op_2
        @species X4(t)
        @observables X ~ X1^2 + op_1*(X2 + 2X3) + X1*X4/op_2 + p        
        (p,d), 0 <--> X1
        (k1,k2), X1 <--> X2
        (k3,k4), X2 <--> X3
    end
    
    u0 = Dict([:X1 => 1.0, :X2 => 2.0, :X3 => 3.0, :X4 => 4.0])
    ps = Dict([:p => 1.0, :d => 0.2, :k1 => 1.5, :k2 => 1.5, :k3 => 5.0, :k4 => 5.0, :op_1 => 1.5, :op_2 => 1.5])

    oprob = ODEProblem(rn, u0, (0.0, 1000.0), ps)
    sol = solve(oprob, Tsit5())

    @test sol[:X][1] == u0[:X1]^2 + ps[:op_1]*(u0[:X2] + 2*u0[:X3]) + u0[:X1]*u0[:X4]/ps[:op_2] + ps[:p]  
end

# Checks that ivs are correctly found.
let
    rn = @reaction_network begin
        @ivs t x y
        @species V1(t) V2(t,x) V3(t, y) W1(t) W2(t, y)
        @observables begin
            V ~ V1 + 2V2 + 3V3
            W ~ W1 + W2
        end
    end
    V,W = getfield.(observed(rn), :lhs)
    @test isequal(arguments(ModelingToolkit.unwrap(V)), Any[rn.iv, rn.sivs[1], rn.sivs[2]])
    @test isequal(arguments(ModelingToolkit.unwrap(W)), Any[rn.iv, rn.sivs[2]])
end

# Checks that metadata is written properly.
let
    rn = @reaction_network rn_observed begin
        @observables (X, [description="my_description"]) ~ X1 + X2
        k, 0 --> X1 + X2
    end
    @test getdescription(observed(rn)[1].lhs) == "my_description"
end

# Declares observables implicitly/explicitly.
# Cannot test `isequal(rn1, rn2)` because the two sets of observables have some obscure Symbolics
# substructure that is different.
let 
    # Basic case.
    rn1 = @reaction_network rn_observed begin
        @observables X ~ X1 + X2
        k, 0 --> X1 + X2
    end
    rn2 = @reaction_network rn_observed begin
        @variables X(t)
        @observables X ~ X1 + X2
        k, 0 --> X1 + X2
    end
    @test isequal(observed(rn1)[1].rhs, observed(rn2)[1].rhs)
    @test isequal(observed(rn1)[1].lhs.metadata, observed(rn2)[1].lhs.metadata)
    @test isequal(states(rn1), states(rn2))

    # Case with metadata.
    rn3 = @reaction_network rn_observed begin
        @observables (X,  [description="description"]) ~ X1 + X2
        k, 0 --> X1 + X2
    end
    rn4 = @reaction_network rn_observed begin
        @variables X(t) [description="description"]
        @observables X ~ X1 + X2
        k, 0 --> X1 + X2
    end
    @test isequal(observed(rn3)[1].rhs, observed(rn4)[1].rhs)
    @test isequal(observed(rn3)[1].lhs.metadata, observed(rn4)[1].lhs.metadata)
    @test isequal(states(rn3), states(rn4))
end

# Tests various erroneous declarations throw errors.
let 
    # Independent variable in @compounds.
    @test_throws Exception @eval @reaction_network begin
        @observables X(t) ~ X1 + X2
        k, 0 --> X1 + X2
    end

    # System with observable in observable formula.
    @test_throws Exception @eval @reaction_network begin
        @observables begin
            X ~ X1 + X2
            X2 ~ 2X
        end
        (p,d), 0 <--> X1 + X2
    end

    # Multiple @compounds options
    @test_throws Exception @eval @reaction_network begin
        @observables X ~ X1 + X2
        @observables Y ~ Y1 + Y2
        k, 0 --> X1 + X2
        k, 0 --> Y1 + Y2
    end
    @test_throws Exception @eval @reaction_network begin
        @observables begin
            X ~ X1 + X2
        end
        @observables begin
            X ~ 2(X1 + X2)
        end
        (p,d), 0 <--> X1 + X2
    end

    # Default value for compound.
    @test_throws Exception @eval @reaction_network begin
        @observables (X = 1.0) ~ X1 + X2
        k, 0 --> X1 + X2
    end

    # Forbidden symbols as observable names.
    @test_throws Exception @eval @reaction_network begin
        @observables t ~ t1 + t2
        k, 0 --> t1 + t2
    end
    @test_throws Exception @eval @reaction_network begin
        @observables im ~ i + m
        k, 0 --> i + m
    end

    # Non-trivial observables expression.
    @test_throws Exception @eval @reaction_network begin
        @observables X - X1 ~ X2
        k, 0 --> X1 + X2
    end

    # Occurrence of undeclared dependants.
    @test_throws Exception @eval @reaction_network begin
        @observables X ~ X1 + X2
        k, 0 --> X1
    end

    # A forbidden symbol used as observable name.
    @test_throws Exception @eval @reaction_network begin
        @observables begin
            t ~ X1 + X2
        end
        (p,d), 0 <--> X1 + X2
    end
end

### Differential Equations ###

# Basic checks on simple case with additional differential equations.
# Checks indexing works.
# Checks that non-block form for single equation works.
let
    # Creates model.
    rn = @reaction_network rn begin
        @parameters k
        @equations begin
            D(V) ~ X - k*V
        end 
        (p,d), 0 <--> X
    end
    
    @unpack k, p, d, X, V = rn
    @variables t
    D = Differential(t)
    
    # Checks that the internal structures have the correct lengths.
    @test length(species(rn)) == 1
    @test length(states(rn)) == 2
    @test length(reactions(rn)) == 2
    @test length(equations(rn)) == 3
    @test has_diff_equations(rn)
    @test isequal(diff_equations(rn), [D(V) ~ X - k*V])
    @test !has_alg_equations(rn)
    @test isequal(alg_equations(rn), [])
    
    # Checks that the internal structures contain the correct stuff, and are correctly sorted.
    @test isspecies(states(rn)[1])
    @test !isspecies(states(rn)[2])
    @test equations(rn)[1] isa Reaction
    @test equations(rn)[2] isa Reaction
    @test equations(rn)[3] isa Equation
    @test isequal(equations(rn)[3], D(V) ~ X - k*V)
    
    # Checks that simulations has the correct output
    u0 = Dict([X => 1 + rand(rng), V => 1 + rand(rng)])
    ps = Dict([p => 1 + rand(rng), d => 1 + rand(rng), k => 1 + rand(rng)])
    oprob = ODEProblem(rn, u0, (0.0, 10000.0), ps)
    sol = solve(oprob, Tsit5(); abstol=1e-9, reltol=1e-9)
    @test sol[X][end] ≈ ps[p]/ps[d]
    @test sol[V][end] ≈ ps[p]/(ps[d]*ps[k])
    
    # Checks that set and get index works for variables.
    @test oprob[V] == u0[V]
    oprob[V] = 2.0
    @test oprob[V] == 2.0
    integrator = init(oprob, Tsit5())
    @test integrator[V] == 2.0
    integrator[V] = 5.0
    @test integrator[V] == 5.0

    # Checks that block form is not required when only a single equation is used.
    rn2 = @reaction_network rn begin
        @parameters k
        @equations D(V) ~ X - k*V
        (p,d), 0 <--> X
    end
    @test rn == rn2
end

# Tries complicated set of equations.
# Tries with pre-declaring some variables (and not others).
# Tries using default values.
# Tries using mixing of parameters, variables, and species in different places.
let 
    rn = @reaction_network begin
        @species S(t)=2.0
        @variables Y(t)=0.5 Z(t)
        @parameters p1 p2 p3
        @equations begin
            D(X) ~ p1*S - X
            D(Y) ~ p2 + X - Y 
            D(Z) ~ p + Y - Z 
        end
        (p*T,d), 0 <--> S
        (p*Z,d), 0 <--> T
    end

    u0 = [:X => 1.0, :Z => 10.0, :S => 1.0, :T => 1.0]
    ps = [:p1 => 0.5, :p2 => 1.0, :p3 => 1.0, :p => 1.0, :d => 1.0]
    oprob = ODEProblem(rn, u0, (0.0, 1000.0), ps)
    sol = solve(oprob, Rosenbrock23(); abstol=1e-9, reltol=1e-9)
    @test sol[:S][end] ≈ 4
    @test sol[:T][end] ≈ 4
    @test sol[:X][end] ≈ 2
    @test sol[:Y][end] ≈ 3
    @test sol[:Z][end] ≈ 4
end

# Tries for reaction system without any reactions (just an equation).
# Tries with interpolating a value into an equation.
# Tries using rn.X notation for designating variables.
# Tries for empty parameter vector.
let 
    c = 4.0
    rn = complete(@reaction_network begin
        @equations D(X) ~ $c - X
    end)

    u0 = [rn.X => 0.0]
    ps = []
    oprob = ODEProblem(rn, u0, (0.0, 100.0), ps)
    sol = solve(oprob, Tsit5(); abstol=1e-9, reltol=1e-9)
    @test sol[rn.X][end] ≈ 4.0
end

# Checks hierarchical model.
let 
    base_rn = @reaction_network begin
        @equations begin
            D(V1) ~ X - 2V1
        end 
        (p,d), 0 <--> X
    end
    @unpack X, V1, p, d = base_rn
    
    internal_rn = @reaction_network begin
        @equations begin
            D(V2) ~ X - 3V2
        end 
        (p,d), 0 <--> X
    end
    
    rn = compose(base_rn, [internal_rn])
    
    u0 = [V1 => 1.0, X => 3.0, internal_rn.V2 => 2.0, internal_rn.X => 4.0]
    ps = [p => 1.0, d => 0.2, internal_rn.p => 2.0, internal_rn.d => 0.5]
    oprob = ODEProblem(rn, u0, (0.0, 1000.0), ps)
    sol = solve(oprob, Tsit5(); abstol=1e-9, reltol=1e-9)
    
    @test sol[X][end] ≈ 5.0
    @test sol[V1][end] ≈ 2.5
    @test sol[internal_rn.X][end] ≈ 4.0
    @test sol[internal_rn.V2][end] ≈ 4/3
end

# Tests that various erroneous declarations throw errors.
let 
    # Using = instead of ~ (for equation).
    @test_throws Exception @eval @reaction_network begin
        @equations D(X) = 1 - X
        (p,d), 0 <--> S
    end

    # Using ~ instead of = (for differential).
    @test_throws Exception @eval @reaction_network begin
        @differentials D ~ Differential(t)
        (p,d), 0 <--> S
    end

    # Equation with component undeclared elsewhere.
    @test_throws Exception @eval @reaction_network begin
        @equations D(X) ~ p - X
        (P,D), 0 <--> S
    end

    # Using default differential D and a symbol D.
    @test_throws Exception @eval @reaction_network begin
        @equations D(X) ~ -X
        (P,D), 0 <--> S
    end

    # Declaring a symbol as a differential when it is used elsewhere.
    @test_throws Exception @eval @reaction_network begin
        @differentials d = Differential(t)
        (p,d), 0 <--> S
    end

    # Declaring differential equation using a forbidden variable.
    @test_throws Exception @eval @reaction_network begin
        @equations D(pi) ~ 1 - pi
        (p,d), 0 <--> S
    end

    # Declaring forbidden symbol as differential.
    @test_throws Exception @eval @reaction_network begin
        @differentials pi = Differential(t)
        (p,d), 0 <--> S
    end

    # Differential with respect to a species.
    @test_throws Exception @eval @reaction_network begin
        @equations D(S) ~ 1 - S
        (p,d), 0 <--> S
    end

    # System with derivatives with respect to several independent variables.
    @test_throws Exception @eval @reaction_network begin
        @ivs s t
        @variables X(s) Y(t)
        @differentials begin
            Ds = Differential(s)
            Dt = Differential(t)
        end
        @equations begin
            Ds(X) ~ 1 - X
            Dt(Y) ~ 1 - Y
        end
    end
end

### Algebraic Equations ###

# Checks creation of basic network.
# Check indexing of output solution.
# Check that DAE is solved correctly.
let
    rn = @reaction_network rn begin
        @parameters k
        @variables X(t) Y(t)
        @equations begin
            X + 5 ~ k*S
            3Y + X  ~ S + X*d
        end 
        (p,d), 0 <--> S
    end

    @unpack X, Y, S, k, p, d = rn

    # Checks that the internal structures have the correct lengths.
    @test length(species(rn)) == 1
    @test length(states(rn)) == 3
    @test length(reactions(rn)) == 2
    @test length(equations(rn)) == 4
    @test !has_diff_equations(rn)
    @test isequal(diff_equations(rn), [])
    @test has_alg_equations(rn)
    @test isequal(alg_equations(rn), [X + 5 ~ k*S, 3Y + X  ~ S + X*d])

    # Checks that the internal structures contain the correct stuff, and are correctly sorted.
    @test isspecies(states(rn)[1])
    @test !isspecies(states(rn)[2])
    @test !isspecies(states(rn)[3])
    @test equations(rn)[1] isa Reaction
    @test equations(rn)[2] isa Reaction
    @test equations(rn)[3] isa Equation
    @test equations(rn)[3] isa Equation
    @test isequal(equations(rn)[3], X + 5 ~ k*S)
    @test isequal(equations(rn)[4], 3Y + X  ~ S + X*d)

    # Checks that simulations has the correct output
    u0 = Dict([S => 1 + rand(rng), X => 1 + rand(rng), Y => 1 + rand(rng)])
    ps = Dict([p => 1 + rand(rng), d => 1 + rand(rng), k => 1 + rand(rng)])
    oprob = ODEProblem(rn, u0, (0.0, 10000.0), ps; structural_simplify=true)
    sol = solve(oprob, Tsit5(); abstol=1e-9, reltol=1e-9)
    @test sol[S][end] ≈ ps[p]/ps[d]
    @test sol[X] .+ 5 ≈ sol[k] .*sol[S]
    @test 3*sol[Y] .+ sol[X] ≈ sol[S] .+ sol[X].*sol[d]
end

# Checks that block form is not required when only a single equation is used.
let
    rn1 = @reaction_network rn begin
        @parameters k
        @variables X(t)
        @equations X + 2 ~ k*S
        (p,d), 0 <--> S
    end
    rn2 = @reaction_network rn begin
        @parameters k
        @variables X(t)
        @equations begin 
            X + 2 ~ k*S
        end
        (p,d), 0 <--> S
    end
    @test rn1 == rn2
end

# Tries for reaction system without any reactions (just an equation).
# Tries with interpolating a value into an equation.
# Tries using rn.X notation for designating variables.
# Tries for empty parameter vector.
let 
    c = 6.0
    rn = complete(@reaction_network begin
        @variables X(t)
        @equations 2X ~ $c - X
    end)

    u0 = [rn.X => 0.0]
    ps = []
    oprob = ODEProblem(rn, u0, (0.0, 100.0), ps; structural_simplify=true)
    sol = solve(oprob, Tsit5(); abstol=1e-9, reltol=1e-9)
    @test sol[rn.X][end] ≈ 2.0
end

# Checks hierarchical model.
let 
    base_rn = @reaction_network begin
        @variables V1(t)
        @equations begin
            X*3V1 ~ X - 2
        end 
        (p,d), 0 <--> X
    end
    @unpack X, V1, p, d = base_rn
    
    internal_rn = @reaction_network begin
        @variables V2(t)
        @equations begin
            X*4V2 ~ X - 3
        end 
        (p,d), 0 <--> X
    end
    
    rn = compose(base_rn, [internal_rn])
    
    u0 = [V1 => 1.0, X => 3.0, internal_rn.V2 => 2.0, internal_rn.X => 4.0]
    ps = [p => 1.0, d => 0.2, internal_rn.p => 2.0, internal_rn.d => 0.5]
    oprob = ODEProblem(rn, u0, (0.0, 1000.0), ps; structural_simplify=true)
    sol = solve(oprob, Rosenbrock23(); abstol=1e-9, reltol=1e-9)
    
    @test sol[X][end] ≈ 5.0
    @test sol[X][end]*3*sol[V1][end] ≈ sol[X][end] - 2
    @test sol[internal_rn.X][end] ≈ 4.0
end

# Check for combined differential and algebraic equation.
# Check indexing of output solution using Symbols.
let
    rn = @reaction_network rn begin
        @parameters k
        @variables X(t) Y(t)
        @equations begin
            X + 5 ~ k*S
            D(Y) ~ X + S - 5*Y
        end 
        (p,d), 0 <--> S
    end

    # Checks that the internal structures have the correct lengths.
    @test length(species(rn)) == 1
    @test length(states(rn)) == 3
    @test length(reactions(rn)) == 2
    @test length(equations(rn)) == 4
    @test has_diff_equations(rn)
    @test length(diff_equations(rn)) == 1
    @test has_alg_equations(rn)
    @test length(alg_equations(rn)) == 1

    # Checks that the internal structures contain the correct stuff, and are correctly sorted.
    @test isspecies(states(rn)[1])
    @test !isspecies(states(rn)[2])
    @test !isspecies(states(rn)[3])
    @test equations(rn)[1] isa Reaction
    @test equations(rn)[2] isa Reaction
    @test equations(rn)[3] isa Equation
    @test equations(rn)[3] isa Equation

    # Checks that simulations has the correct output
    u0 = Dict([S => 1 + rand(rng), X => 1 + rand(rng), Y => 1 + rand(rng)])
    ps = Dict([p => 1 + rand(rng), d => 1 + rand(rng), k => 1 + rand(rng)])
    oprob = ODEProblem(rn, u0, (0.0, 10000.0), ps; structural_simplify=true)
    sol = solve(oprob, Tsit5(); abstol=1e-9, reltol=1e-9)
    @test sol[:S][end] ≈ sol[:p]/sol[:d]
    @test sol[:X] .+ 5 ≈ sol[:k] .*sol[:S]
    @test 5*sol[:Y][end] ≈ sol[:S][end] + sol[:X][end]
end

# Tests that various erroneous declarations throw errors.
let 
    # Using = instead of ~ (for equation).
    @test_throws Exception @eval @reaction_network begin
        @variables X(t)
        @equations X = 1 - S
        (p,d), 0 <--> S
    end

    # Equation with component undeclared elsewhere.
    @test_throws Exception @eval @reaction_network begin
        @equations X ~ p - S
        (P,D), 0 <--> S
    end
end

### Events ###

# Compares models with complicated events that are created programmatically/with the DSL.
# Checks that simulations are correct.
# Checks that various simulation inputs works.
# Checks continuous, discrete, preset time, and periodic events.
# Tests event affecting non-species components.

let
    # Creates model via DSL.
    rn_dsl = @reaction_network rn begin
        @parameters thres=1.0 dY_up
        @variables Z(t)
        @continuous_events begin
            t - 2.5 => p ~ p + 0.2
            [X - thres, Y - X] => [X ~ X - 0.5, Z ~ Z + 0.1]
        end
        @discrete_events begin
            2.0 => [dX ~ dX + 0.1, dY ~ dY + dY_up]
            [1.0, 5.0] => [p ~ p - 0.1]
            [Z > Y, Z > X] => [Z ~ Z - 0.1]
        end

        (p, dX), 0 <--> X
        (p, dY), 0 <--> Y
    end

    # Creates model programmatically.
    @variables t Z(t)
    @species X(t) Y(t)
    @parameters p dX dY thres=1.0 dY_up
    rxs = [
        Reaction(p, nothing, [X], nothing, [1])
        Reaction(dX, [X], nothing, [1], nothing)
        Reaction(p, nothing, [Y], nothing, [1])
        Reaction(dY, [Y], nothing, [1], nothing)
    ]
    continuous_events = [
        t - 2.5 => p ~ p + 0.2
        [X - thres, Y - X] => [X ~ X - 0.5, Z ~ Z + 0.1]
    ]
    discrete_events = [
        2.0 => [dX ~ dX + 0.1, dY ~ dY + dY_up]
        [1.0, 5.0] => [p ~ p - 0.1]
        [Z > Y, Z > X] => [Z ~ Z - 0.1]
    ]
    rn_prog = ReactionSystem([rx1, rx2, eq], t; continuous_events, discrete_events, name=:rn)

    # Tests that approaches yield identical results.
    @test isequal(rn_dsl, rn_prog)

    u0 = [X => 1.0, Y => 0.5, Z => 0.25]
    tspan = (0.0, 20.0)
    ps = [p => 1.0, dX => 0.5, dY => 0.5, dY_up => 0.1]

    sol_dsl = solve(ODEProblem(rn_dsl, u0, tspan, ps), Tsit5())
    sol_prog = solve(ODEProblem(rn_prog, u0, tspan, ps), Tsit5())
    @test sol_dsl == sol_prog
end

# Compares DLS events to those given as callbacks.
# Checks that events works when given to SDEs.
let
    # Creates models.
    rn = @reaction_network begin
        (p, d), 0 <--> X
    end
    rn_events = @reaction_network begin
        @discrete_events begin
            [5.0, 10.0] => [X ~ X + 100.0]
        end
        @continuous_events begin
            [X - 90.0] => [X ~ X + 10.0]
        end
        (p, d), 0 <--> X
    end
    cb_disc = PresetTimeCallback([5.0, 10.0], int -> (int[:X] += 100.0))
    cb_cont = ContinuousCallback((u, t, int) -> (u[1] - 90.0), int -> (int[:X] += 10.0))

    # Simulates models,.
    u0 = [:X => 100.0]
    tspan = (0.0, 50.0)
    ps = [:p => 100.0, :d => 1.0]
    sol = solve(SDEProblem(rn, u0, tspan, ps), ImplicitEM();  seed = 1234, callback = CallbackSet(cb_disc, cb_cont))
    sol_events = solve(SDEProblem(rn_events, u0, tspan, ps), ImplicitEM(); seed = 1234)

    @test sol == sol_events
end