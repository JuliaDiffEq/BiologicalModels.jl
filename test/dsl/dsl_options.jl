#! format: off

### Prepares Tests ###

# Fetch packages.
using Catalyst, ModelingToolkit, OrdinaryDiffEqTsit5, OrdinaryDiffEqRosenbrock, StochasticDiffEq, Plots, Test
using Symbolics: unwrap

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)

# Sets the default `t` to use.
t = default_t()
D = default_time_deriv()

### Tests `@parameters`, `@species`, and `@variables` Options ###

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

    n1 = @reaction_network rnname begin (k1, k2), A <--> B end
    n2 = @reaction_network rnname begin
        @parameters k1 k2
        (k1, k2), A <--> B
    end
    n3 = @reaction_network rnname begin
        @species A(t) B(t)
        (k1, k2), A <--> B
    end
    n4 = @reaction_network rnname begin
        @parameters k1 k2
        @species A(t) B(t)
        (k1, k2), A <--> B
    end
    n5 = @reaction_network rnname begin
        (k1, k2), A <--> B
        @parameters k1 k2
    end
    n6 = @reaction_network rnname begin
        (k1, k2), A <--> B
        @species A(t) B(t)
    end
    n7 = @reaction_network rnname begin
        (k1, k2), A <--> B
        @parameters k1 k2
        @species A(t) B(t)
    end
    n8 = @reaction_network rnname begin
        @parameters begin
            k1
            k2
        end
        (k1, k2), A <--> B
    end
    n9 = @reaction_network rnname begin
        @species begin
            A(t)
            B(t)
        end
        (k1, k2), A <--> B
    end
    n10 = @reaction_network rnname begin
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

# Tests that when either @species or @parameters is given, the other is inferred properly.
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

# Checks that some created networks are identical.
let
    rn12 = @reaction_network rnname begin (k1, k2), A <--> B end
    rn13 = @reaction_network rnname begin
        @parameters k1 k2
        (k1, k2), A <--> B
    end
    rn14 = @reaction_network rnname begin
        @species A(t) B(t)
        (k1, k2), A <--> B
    end
    rn15 = @reaction_network rnname begin
        @parameters k1 k2
        @species A(t) B(t)
        (k1, k2), A <--> B
    end
    @test all(==(rn12), (rn13, rn14, rn15))
end

# Checks that the rights things are put in vectors.
let
    rn18 = @reaction_network rnname begin
        @parameters p d1 d2
        @species A(t) B(t)
        p, 0 --> A
        1, A --> B
        (d1, d2), (A, B) --> 0
    end
    rn19 = @reaction_network rnname begin
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

    rn20 = @reaction_network rnname begin
        @species X(t)
        @parameters S
        mm(X,v,K), 0 --> Y
        (k1,k2), 2Y <--> Y2
        d*Y, S*(Y2+Y) --> 0
    end
    rn21 = @reaction_network rnname begin
        @species X(t) Y(t) Y2(t)
        @parameters v K k1 k2 d S
        mm(X,v,K), 0 --> Y
        (k1,k2), 2Y <--> Y2
        d*Y, S*(Y2+Y) --> 0
    end
    rn22 = @reaction_network rnname begin
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
    rn26 = @reaction_network rnname begin
        @parameters p=1.0 d1 d2=5
        @species A(t) B(t)=4
        p, 0 --> A
        1, A --> B
        (d1, d2), (A, B) --> 0
    end

    rn27 = @reaction_network rnname begin
    @parameters p1=1.0 p2=2.0 k1=4.0 k2=5.0 v=8.0 K=9.0 n=3 d=10.0
    @species X(t)=4.0 Y(t)=3.0 X2Y(t)=2.0 Z(t)=1.0
        (p1,p2), 0 --> (X,Y)
        (k1,k2), 2X + Y --> X2Y
        hill(X2Y,v,K,n), 0 --> Z
        d, (X,Y,X2Y,Z) --> 0
    end
    u0_27 = []
    p_27 = []

    rn28 = @reaction_network rnname begin
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

    rn29 = @reaction_network rnname begin
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

# Tests that parameter type designation works.
let
    # Creates a model.
    rn = @reaction_network begin
        @parameters begin
            k1
            l1
            k2::Float64 = 2.0
            l2::Float64
            k3::Int64 = 2, [description="A parameter"]
            l3::Int64
            k4::Float32, [description="Another parameter"]
            l4::Float32
            k5::Rational{Int64}
            l5::Rational{Int64}
        end
        (k1,l1), X1 <--> Y1
        (k2,l2), X2 <--> Y2
        (k3,l3), X3 <--> Y3
        (k4,l4), X4 <--> Y4
        (k5,l5), X5 <--> Y5
    end

    # Checks parameter types.
    @test unwrap(rn.k1) isa SymbolicUtils.BasicSymbolic{Real}
    @test unwrap(rn.l1) isa SymbolicUtils.BasicSymbolic{Real}
    @test unwrap(rn.k2) isa SymbolicUtils.BasicSymbolic{Float64}
    @test unwrap(rn.l2) isa SymbolicUtils.BasicSymbolic{Float64}
    @test unwrap(rn.k3) isa SymbolicUtils.BasicSymbolic{Int64}
    @test unwrap(rn.l3) isa SymbolicUtils.BasicSymbolic{Int64}
    @test unwrap(rn.k4) isa SymbolicUtils.BasicSymbolic{Float32}
    @test unwrap(rn.l4) isa SymbolicUtils.BasicSymbolic{Float32}
    @test unwrap(rn.k5) isa SymbolicUtils.BasicSymbolic{Rational{Int64}}
    @test unwrap(rn.l5) isa SymbolicUtils.BasicSymbolic{Rational{Int64}}

    # Checks that other parameter properties are assigned properly.
    @test !ModelingToolkit.hasdefault(unwrap(rn.k1))
    @test ModelingToolkit.getdefault(unwrap(rn.k2)) == 2.0
    @test ModelingToolkit.getdefault(unwrap(rn.k3)) == 2
    @test ModelingToolkit.getdescription(unwrap(rn.k3)) == "A parameter"
    @test ModelingToolkit.getdescription(unwrap(rn.k4)) == "Another parameter"
    @test !ModelingToolkit.hasdescription(unwrap(rn.k5))
end

# Test @variables in DSL.
let
    rn = @reaction_network tester begin
        @parameters k1
        @variables V1(t) V2(t) V3(t)
        @species B1(t) B2(t)
        (k1*k2 + V3), V1*A + 2*B1 --> V2*C + B2
    end

    @parameters k1 k2
    @variables V1(t) V2(t) V3(t)
    @species A(t) B1(t) B2(t) C(t)
    rx = Reaction(k1*k2 + V3, [A, B1], [C, B2], [V1, 2], [V2, 1])
    @named tester = ReactionSystem([rx], t)
    @test tester == rn

    sts = (A, B1, B2, C, V1, V2, V3)
    spcs = (A, B1, B2, C)
    @test issetequal(unknowns(rn), sts)
    @test issetequal(species(rn), spcs)
end

# Tests error when disallowed name is used for variable.
let
    @test_throws Exception @eval @reaction_network begin
        @variables π(t)
    end
    @test_throws Exception @eval @reaction_network begin
        @variables Γ(t)
    end
    @test_throws Exception @eval @reaction_network begin
        @variables ∅(t)
    end
    @test_throws Exception @eval @reaction_network begin
        @ivs s
        @variables t(s)
    end
end

# Tests that explicitly declaring a single symbol as several things does not work.
# Several of these are broken, but note sure how to test broken-ness on `@test_throws false Exception @eval`.
let
    # Species + parameter.
    @test_throws Exception @eval @reaction_network begin
        @species X(t)
        @parameters X
    end

    # Species + variable.
    @test_throws Exception @eval @reaction_network begin
        @species X(t)
        @variables X(t)
    end

    # Variable + parameter.
    @test_throws Exception @eval @reaction_network begin
        @variables X(t)
        @parameters X
    end

    # Species + differential.
    @test_throws Exception @eval @reaction_network begin
        @species X(t)
        @differentials X = Differential(t)
    end

    # Parameter + differential.
    @test_throws Exception @eval @reaction_network begin
        @parameters X
        @differentials X = Differential(t)
    end

    # Variable + differential.
    @test_throws Exception @eval @reaction_network begin
        @variables X(t)
        @differentials X = Differential(t)
    end

    # Parameter + observable (species/variable + observable is OK, as this e.g. provide additional observables information).
    @test_throws Exception @eval @reaction_network begin
        @species Y(t)
        @parameters X
        @observables X ~ Y
    end

    # Species + compound.
    @test_throws Exception @eval @reaction_network begin
        @species X(t) O(t)
        @compounds begin X(t) ~ 2O end
    end

    # Parameter + compound.
   @test_throws Exception @eval @reaction_network begin
        @species O(t)
        @parameters X
        @compounds begin X(t) ~ 2O end
    end

    # Variable + compound.
    @test_throws Exception @eval @reaction_network begin
        @species O(t)
        @variables X(t)
        @compounds begin X(t) ~ 2O end
    end
end

### Test Independent Variable Designations ###

# Test ivs in DSL.
let
    rn = @reaction_network ivstest begin
        @ivs s x
        @parameters k2
        @variables D(x) E(s) F(s,x)
        @species A(s,x) B(s) C(x)
        k*k2*D, E*A +B --> F*C + C2
    end

    @parameters k k2
    @parameters s x
    @variables D(x) E(s) F(s,x)
    @species A(s,x) B(s) C(x) C2(s,x)
    rx = Reaction(k*k2*D, [A, B], [C, C2], [E, 1], [F, 1])
    @named ivstest = ReactionSystem([rx], s; spatial_ivs = [x])

    @test ivstest == rn
    @test issetequal(unknowns(rn), [D, E, F, A, B, C, C2])
    @test issetequal(species(rn), [A, B, C, C2])
    @test isequal(ModelingToolkit.get_iv(rn), s)
    @test issetequal(Catalyst.get_sivs(rn), [x])
end

### Test Symbolic Variable Inference ###

# Basic checks that that symbolic variables not explicitly declared are correctly inferred.
let
    # Case 1 (a reaction only).
    rn1 = @reaction_network begin
        (p1/(S1+p2) + S2), S1 --> S2
    end
    @test issetequal(species(rn1), [rn1.S1, rn1.S2])
    @test issetequal(parameters(rn1), [rn1.p1, rn1.p2])

    # Case 2 (reactions and equations).
    rn2 = @reaction_network begin
        @equations V1 + log(V2 + S1) ~ V2^2
        (p1/V1 + S1 + log(S2 + V2 + p2)), S1 --> S2
    end
    @test issetequal(species(rn2), [rn2.S1, rn2.S2])
    @test issetequal(nonspecies(rn2), [rn2.V1, rn2.V2])
    @test issetequal(parameters(rn2), [rn2.p1, rn2.p2])

    # Case 3 (reaction and equations with a differential).
    rn3 = @reaction_network begin
        @equations begin
            D(V1) ~ S1 + V1
            V2 + S2 ~ V1^2 + V2^2
        end
        (p1/V1 + S1 + log(S2 + V2 + p2)), S1 --> S2
    end
    @test issetequal(species(rn3), [rn3.S1, rn3.S2])
    @test issetequal(nonspecies(rn3), [rn3.V1, rn3.V2])
    @test issetequal(parameters(rn3), [rn3.p1, rn3.p2])

    # Case 4 (reactions and equations with a pre-declared parameter).
    rn4 = @reaction_network begin
        @parameters p1
        @equations V1 + sin(p1 + S1) ~ S2*V2
        (p1+p2+V1+V2+S1+S2), S1 --> S2
    end
    @test issetequal(species(rn4), [rn2.S1, rn2.S2])
    @test issetequal(nonspecies(rn4), [rn2.V1, rn2.V2])
    @test issetequal(parameters(rn4), [rn2.p1, rn2.p2])

    # Case 5 (algebraic equation containing D, which is pre-declared as a species).
    rn5 = @reaction_network begin
        @species D(t)
        @equations D * (S1 + V1 + V2) ~ S2
        (p1 + p2*(D + V1 + V2 + S2 + S2)), S1 --> S2 + D
    end
    @test issetequal(species(rn5), [rn5.S1, rn5.S2, rn5.D])
    @test issetequal(nonspecies(rn5), [rn5.V1, rn5.V2])
    @test issetequal(parameters(rn5), [rn5.p1, rn5.p2])

    # Case 6 (algebraic equation containing D, which is pre-declared as a parameter).
    rn6 = @reaction_network begin
        @parameters D
        @equations D * (S1 + V1 + V2) ~ S2
        (p1 + p2*(D + V1 + V2 + S2 + S2)), S1 --> S2
    end
    @test issetequal(species(rn6), [rn6.S1, rn6.S2])
    @test issetequal(nonspecies(rn6), [rn6.V1, rn6.V2])
    @test issetequal(parameters(rn6), [rn6.p1, rn6.p2, rn6.D])

    # Case 7 (algebraic equation containing D, which is pre-declared as a variable).
    rn7 = @reaction_network begin
        @variables D(t)
        @equations D * (S1 + V1 + V2) ~ S2
        (p1 + p2*(D + V1 + V2 + S2 + S2)), S1 --> S2
    end
    @test issetequal(species(rn7), [rn7.S1, rn7.S2])
    @test issetequal(nonspecies(rn7), [rn7.V1, rn7.V2, rn7.D])
    @test issetequal(parameters(rn7), [rn7.p1, rn7.p2])

    # Case 8 (reactions, equations, and a custom differential).
    rn8 = @reaction_network begin
        @differentials Δ = Differential(t)
        @equations Δ(V1) + Δ(V2) + log(V2 + S1) ~  S2
        (p1/V1 + S1 + log(S2 + V2 + p2)), S1 --> S2
    end
    @test issetequal(species(rn8), [rn8.S1, rn8.S2])
    @test issetequal(nonspecies(rn8), [rn8.V1, rn8.V2])
    @test issetequal(parameters(rn8), [rn8.p1, rn8.p2])
end

# Checks that various cases where symbolic variables cannot (or shouldn't) be inferred generate errors.
let
    # Species/variables/parameter named after default differential used as function call.
    # In the future, species/variables should be usable this way (designating a time delay).
    @test_throws Exception @eval @reaction_network begin
        @equations D(V) ~ 1 - V
        d, D --> 0
    end
    @test_throws Exception @eval @reaction_network begin
        @variables D(t)
        @equations D(V) ~ 1 - V
        d, X --> 0
    end
    @test_throws Exception @eval @reaction_network begin
        @parameters D
        @equations D(V) ~ 1 - V
        d, X --> 0
    end

    # Symbol only occurring in events.
    @test_throws Exception @eval @reaction_network begin
        @discrete_event (X > 1.0) => [V => V/2]
        d, X --> 0
    end
    @test_throws Exception @eval @reaction_network begin
        @continuous_event [X > 1.0] => [V => V/2]
        d, X --> 0
    end
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
    @test length(unknowns(rn)) == 4
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
    plot(sol; idxs=X).series_list[1].plotattributes[:y][end] ≈ 10.0
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

    @variables X(t) Y(t)
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
    rn_prog = complete(rn_prog)

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

# Checks that models created w/o specifying `@variables` for observables are identical.
# Compares both to model with explicit declaration, and programmatically created model.
let
    # With default ivs.
    rn1 = @reaction_network rn begin
        @variables X(t) X1(t) X2(t)
        @observables X ~ X1 + X2
    end
    rn2 = @reaction_network rn begin
        @variables X1(t) X2(t)
        @observables X ~ X1 + X2
    end
    @variables X(t) X1(t) X2(t)
    rn3 = complete(ReactionSystem([], t, [X1, X2], []; name = :rn, observed = [X ~ X1 + X2]))
    @test isequal(rn1, rn2)
    @test isequal(rn1, rn3)
    @test isequal(rn1.X, rn2.X)
    @test isequal(rn1.X, rn3.X)

    # With non-default ivs.
    rn4 = @reaction_network rn begin
        @ivs τ x
        @variables X(τ,x) X1(τ,x) X2(τ,x)
        @observables X ~ X1 + X2
    end
    rn5 = @reaction_network rn begin
        @ivs τ x
        @variables X1(τ,x) X2(τ,x)
        @observables X ~ X1 + X2
    end
    @parameters τ x
    @variables X(τ,x) X1(τ,x) X2(τ,x)
    rn6 = complete(ReactionSystem([], τ, [X1, X2], []; name = :rn, observed = [X ~ X1 + X2], spatial_ivs = [x]))
    @test isequal(rn4, rn5)
    @test isequal(rn4, rn6)
    @test isequal(rn4.X, rn5.X)
    @test isequal(rn4.X, rn6.X)
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
    @test isequal(Symbolics.sorted_arguments(ModelingToolkit.unwrap(V)), Any[Catalyst.get_iv(rn), Catalyst.get_sivs(rn)[1], Catalyst.get_sivs(rn)[2]])
    @test isequal(Symbolics.sorted_arguments(ModelingToolkit.unwrap(W)), Any[Catalyst.get_iv(rn), Catalyst.get_sivs(rn)[2]])
end

# Checks that metadata is written properly.
let
    rn = @reaction_network rn_observed begin
        @observables (X, [description="my_description"]) ~ X1 + X2
        k, 0 --> X1 + X2
    end
    @test ModelingToolkit.getdescription(observed(rn)[1].lhs) == "my_description"
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
    @test isequal(unknowns(rn1), unknowns(rn2))

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
    @test isequal(unknowns(rn3), unknowns(rn4))
end

# Tests for interpolation into the observables option.
let
    # Interpolation into lhs.
    @species X [description="An observable"]
    rn1 = @reaction_network begin
        @observables $X ~ X1 + X2
        (k1, k2), X1 <--> X2
    end
    @test isequal(observed(rn1)[1].lhs, X)
    @test ModelingToolkit.getdescription(rn1.X) == "An observable"
    @test isspecies(rn1.X)
    @test length(unknowns(rn1)) == 2

    # Interpolation into rhs.
    @parameters n [description="A parameter"]
    @species S(t)
    rn2 = @reaction_network begin
        @observables Stot ~ $S + $n*Sn
        (kB, kD), $n*S <--> Sn
    end
    @unpack Stot, Sn, kD, kB = rn2

    u0 = Dict([S => 5.0, Sn => 1.0])
    ps = Dict([n => 2, kB => 1.0, kD => 1.0])
    oprob = ODEProblem(rn2, u0, (0.0, 1.0), ps)

    @test issetequal(Symbolics.get_variables(observed(rn2)[1].rhs), [S, n, Sn])
    @test oprob[Stot] == u0[S] + ps[n]*u0[Sn]
    @test length(unknowns(rn2)) == 2
end

# Tests specific declaration of observables as species/variables.
let
    rn = @reaction_network begin
        @species X(t)
        @variables Y(t)
        @observables begin
            X ~ X + 2X2
            Y ~ Y1 + Y2
            Z ~ X + Y
        end
        (kB,kD), 2X <--> X2
        (k1,k2), Y1 <--> Y2
    end

    @test isspecies(rn.X)
    @test !isspecies(rn.Y)
    @test !isspecies(rn.Z)
end

# Tests various erroneous declarations throw errors.
let
    # Independent variable in @observables.
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

    # Multiple @observables options
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

    # Interpolation and explicit declaration of an observable.
    @variables X(t)
    @test_throws Exception @eval @reaction_network begin
        @variables X(t)
        @observables $X ~ X1 + X2
        (k1,k2), X1 <--> X2
    end

    # Observable metadata provided twice.
    @test_throws Exception @eval @reaction_network begin
        @species X2 [description="Twice the amount of X"]
        @observables (X2, [description="X times two."]) ~ 2X
        d, X --> 0
    end

end


### Test `@equations` Option for Coupled CRN/Equations Models ###

# Checks creation of basic network.
# Check indexing of output solution.
# Check that DAE is solved correctly.
let
    rn = @reaction_network rn begin
        @parameters k d
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
    @test length(unknowns(rn)) == 3
    @test length(reactions(rn)) == 2
    @test length(equations(rn)) == 4
    @test !has_diff_equations(rn)
    @test isequal(diff_equations(rn), [])
    @test has_alg_equations(rn)
    @test isequal(alg_equations(rn), [X + 5 ~ k*S, 3Y + X  ~ S + X*d])

    # Checks that the internal structures contain the correct stuff, and are correctly sorted.
    @test isspecies(unknowns(rn)[1])
    @test !isspecies(unknowns(rn)[2])
    @test !isspecies(unknowns(rn)[3])
    @test equations(rn)[1] isa Reaction
    @test equations(rn)[2] isa Reaction
    @test equations(rn)[3] isa Equation
    @test equations(rn)[3] isa Equation
    @test isequal(equations(rn)[3], X + 5 ~ k*S)
    @test isequal(equations(rn)[4], 3Y + X  ~ S + X*d)

    # Checks that simulations has the correct output
    u0 = Dict([S => 1 + rand(rng)])
    ps = Dict([p => 1 + rand(rng), d => 1 + rand(rng), k => 1 + rand(rng)])
    oprob = ODEProblem(rn, u0, (0.0, 10000.0), ps; structural_simplify=true)
    sol = solve(oprob, Tsit5(); abstol=1e-9, reltol=1e-9)
    @test sol[S][end] ≈ sol.ps[p]/sol.ps[d]
    @test sol[X] .+ 5 ≈ sol.ps[k] .*sol[S]
    @test 3*sol[Y] .+ sol[X] ≈ sol[S] .+ sol[X].*sol.ps[d]
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
    oprob = ODEProblem(rn, [], (0.0, 100.0); structural_simplify=true)
    sol = solve(oprob, Tsit5(); abstol=1e-9, reltol=1e-9)
    @test sol[rn.X][end] ≈ 2.0
end

# Checks hierarchical model.
let
    base_rn = @network_component begin
        @variables V1(t)
        @equations begin
            X*3V1 ~ X - 2
        end
        (p,d), 0 <--> X
    end
    @unpack X, V1, p, d = base_rn

    internal_rn = @network_component begin
        @variables V2(t)
        @equations begin
            X*4V2 ~ X - 3
        end
        (p,d), 0 <--> X
    end

    rn = complete(compose(base_rn, [internal_rn]))

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
    @unpack X, Y, S, p, d, k = rn

    # Checks that the internal structures have the correct lengths.
    @test length(species(rn)) == 1
    @test length(unknowns(rn)) == 3
    @test length(reactions(rn)) == 2
    @test length(equations(rn)) == 4
    @test has_diff_equations(rn)
    @test length(diff_equations(rn)) == 1
    @test has_alg_equations(rn)
    @test length(alg_equations(rn)) == 1

    # Checks that the internal structures contain the correct stuff, and are correctly sorted.
    @test isspecies(unknowns(rn)[1])
    @test !isspecies(unknowns(rn)[2])
    @test !isspecies(unknowns(rn)[3])
    @test equations(rn)[1] isa Reaction
    @test equations(rn)[2] isa Reaction
    @test equations(rn)[3] isa Equation
    @test equations(rn)[3] isa Equation

    # Checks that simulations has the correct output
    u0 = Dict([S => 1 + rand(rng), Y => 1 + rand(rng)])
    ps = Dict([p => 1 + rand(rng), d => 1 + rand(rng), k => 1 + rand(rng)])
    oprob = ODEProblem(rn, u0, (0.0, 10000.0), ps; structural_simplify=true)
    sol = solve(oprob, Tsit5(); abstol=1e-9, reltol=1e-9)
    @test sol[:S][end] ≈ sol.ps[:p]/sol.ps[:d]
    @test sol[:X] .+ 5 ≈ sol.ps[:k] .*sol[:S]
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

    # Differential equation using a forbidden variable (in the DSL).
    @test_throws Exception @eval @reaction_network begin
        @equations D(π) ~ -1
    end

    # Algebraic equation using a forbidden variable (in the DSL).
    @test_throws Exception @eval @reaction_network begin
        @equations Γ ~ 1 + 3(Γ^2 + Γ)
    end
end

### Other DSL Option Tests ###

# Test that various options can be provided in block and single line form.
# Also checks that the single line form takes maximally one argument.
let
    # The `@equations` option.
    rn11 = @reaction_network rn1 begin
        @equations D(V) ~ 1 - V
    end
    rn12 = @reaction_network rn1 begin
        @equations begin
            D(V) ~ 1 - V
        end
    end
    @test isequal(rn11, rn12)
    @test_throws Exception @eval @reaction_network begin
        @equations D(V) ~ 1 - V D(W) ~ 1 - W
    end

    # The `@observables` option.
    rn21 = @reaction_network rn1 begin
        @species X(t)
        @observables X2 ~ 2X
    end
    rn22 = @reaction_network rn1 begin
        @species X(t)
        @observables begin
            X2 ~ 2X
        end
    end
    @test isequal(rn21, rn22)
    @test_throws Exception @eval @reaction_network begin
        @species X(t)
        @observables X2 ~ 2X X3 ~ 3X
    end

    # The `@compounds` option.
    rn31 = @reaction_network rn1 begin
        @species X(t)
        @compounds X2 ~ 2X
    end
    rn32 = @reaction_network rn1 begin
        @species X(t)
        @compounds begin
            X2 ~ 2X
        end
    end
    @test isequal(rn31, rn32)
    @test_throws Exception @eval @reaction_network begin
        @species X(t)
        @compounds X2 ~ 2X X3 ~ 3X
    end

    # The `@differentials` option.
    rn41 = @reaction_network rn1 begin
        @differentials D = Differential(t)
    end
    rn42 = @reaction_network rn1 begin
        @differentials begin
            D = Differential(t)
        end
    end
    @test isequal(rn41, rn42)
    @test_throws Exception @eval @reaction_network begin
        @differentials D = Differential(t) Δ = Differential(t)
    end

    # The `@continuous_events` option.
    rn51 = @reaction_network rn1 begin
        @species X(t)
        @continuous_events [X ~ 3.0] => [X ~ X - 1]
    end
    rn52 = @reaction_network rn1 begin
        @species X(t)
        @continuous_events begin
            [X ~ 3.0] => [X ~ X - 1]
        end
    end
    @test isequal(rn51, rn52)
    @test_throws Exception @eval @reaction_network begin
        @species X(t)
        @continuous_events [X ~ 3.0] => [X ~ X - 1] [X ~ 1.0] => [X ~ X + 1]
    end

    # The `@discrete_events` option.
    rn61 = @reaction_network rn1 begin
        @species X(t)
        @discrete_events [X > 3.0] => [X ~ X - 1]
    end
    rn62 = @reaction_network rn1 begin
        @species X(t)
        @discrete_events begin
            [X > 3.0] => [X ~ X - 1]
        end
    end
    @test isequal(rn61, rn62)
    @test_throws Exception @eval @reaction_network begin
        @species X(t)
        @discrete_events [X > 3.0] => [X ~ X - 1] [X < 1.0] => [X ~ X + 1]
    end
end

# test combinatoric_ratelaws DSL option
let
    # Test for `@combinatoric_ratelaws false`.
    rn = @reaction_network begin
        @combinatoric_ratelaws false
        (k1,k2), 2A <--> B
    end
    combinatoric_ratelaw = Catalyst.get_combinatoric_ratelaws(rn)
    @test combinatoric_ratelaw == false
    rl = oderatelaw(reactions(rn)[1]; combinatoric_ratelaw)
    @unpack k1, A = rn
    @test isequal(rl, k1*A^2)

    # Test for `@combinatoric_ratelaws true`.
    rn2 = @reaction_network begin
        @combinatoric_ratelaws true
        (k1,k2), 2A <--> B
    end
    combinatoric_ratelaw = Catalyst.get_combinatoric_ratelaws(rn2)
    @test combinatoric_ratelaw == true
    rl = oderatelaw(reactions(rn2)[1]; combinatoric_ratelaw)
    @unpack k1, A = rn2
    @test isequal(rl, k1*A^2/2)

    # Test for interpolation into `@combinatoric_ratelaws`.
    crl = false
    rn3 = @reaction_network begin
        @combinatoric_ratelaws $crl
        (k1,k2), 2A <--> B
    end
    combinatoric_ratelaw = Catalyst.get_combinatoric_ratelaws(rn3)
    @test combinatoric_ratelaw == crl
    rl = oderatelaw(reactions(rn3)[1]; combinatoric_ratelaw)
    @unpack k1, A = rn3
    @test isequal(rl, k1*A^2)

    # Test for erroneous inputs (to few, to many, wrong type).
    @test_throws Exception @eval @reaction_network begin
        @combinatoric_ratelaws
        d, 3X --> 0
    end
    @test_throws Exception @eval @reaction_network begin
        @combinatoric_ratelaws false false
        d, 3X --> 0
    end
    @test_throws Exception @eval @reaction_network begin
        @combinatoric_ratelaws "false"
        d, 3X --> 0
    end
end

# Test whether user-defined functions are properly expanded in equations.
let
    f(A, t) = 2*A*t
    g(A) = 2*A + 2

    # Test user-defined function (on both lhs and rhs).
    rn = @reaction_network begin
        @equations D(A) + g(A) ~ f(A, t)
    end
    @test length(equations(rn)) == 1
    @test equations(rn)[1] isa Equation
    @species A(t)
    @test isequal(equations(rn)[1], D(A) + 2*A + 2 ~ 2*A*t)

    # Test whether expansion happens properly for unregistered/registered functions.
    hill_unregistered(A, v, K, n) = v*(A^n) / (A^n + K^n)
    rn2 = @reaction_network begin
        @parameters v K n
        @equations D(A) ~ hill_unregistered(A, v, K, n)
    end
    @test length(equations(rn2)) == 1
    @test equations(rn2)[1] isa Equation
    @parameters v K n
    @test isequal(equations(rn2)[1], D(A) ~ v*(A^n) / (A^n + K^n))

    hill2(A, v, K, n) = v*(A^n) / (A^n + K^n)
    @register_symbolic hill2(A, v, K, n)
    # Registered symbolic function should not expand.
    rn2r = @reaction_network begin
        @parameters v K n
        @equations D(A) ~ hill2(A, v, K, n)
    end
    @test length(equations(rn2r)) == 1
    @test equations(rn2r)[1] isa Equation
    @parameters v K n
    @test isequal(equations(rn2r)[1], D(A) ~ hill2(A, v, K, n))

    rn3 = @reaction_network begin
        @species Iapp(t)
        @equations begin
            D(A) ~ Iapp
            Iapp ~ f(A,t)
        end
    end
    @test length(equations(rn3)) == 2
    @test equations(rn3)[1] isa Equation
    @test equations(rn3)[2] isa Equation
    @variables Iapp(t)
    @test isequal(equations(rn3)[1], D(A) ~ Iapp)
    @test isequal(equations(rn3)[2], Iapp ~ 2*A*t)

    # Test whether the DSL and symbolic ways of creating the network generate the same system
    @species Iapp(t) A(t)
    eq = [D(A) ~ Iapp, Iapp ~ f(A, t)]
    @named rn3_sym = ReactionSystem(eq, t)
    rn3_sym = complete(rn3_sym)
    @test isequivalent(rn3, rn3_sym)

    # Test more complicated expression involving both registered function and a user-defined function.
    g(A, K, n) = A^n + K^n
    rn4 = @reaction_network begin
        @parameters v K n
        @equations D(A) ~ hill(A, v, K, n)*g(A, K, n)
    end
    @test length(equations(rn4)) == 1
    @test equations(rn4)[1] isa Equation
    @parameters v n
    @test isequal(Catalyst.expand_registered_functions(equations(rn4)[1]), D(A) ~ v*(A^n))
end

# Erroneous `@default_noise_scaling` declaration (other noise scaling tests are mostly in the SDE file).
let
    # Default noise scaling with multiple entries.
    @test_throws Exception @eval @reaction_network begin
        @default_noise_scaling η1 η2
    end
end

# test that @require_declaration properly throws errors when undeclared variables are written.
let
    import Catalyst: UndeclaredSymbolicError

    # Test error when species are inferred
    @test_throws UndeclaredSymbolicError @macroexpand @reaction_network begin
        @require_declaration
        @parameters k
        k, A --> B
    end
    @test_nowarn @macroexpand @reaction_network begin
        @require_declaration
        @species A(t) B(t)
        @parameters k
        k, A --> B
    end

    # Test error when a parameter in rate is inferred
    @test_throws UndeclaredSymbolicError @macroexpand @reaction_network begin
        @require_declaration
        @species A(t) B(t)
        @parameters k
        k*n, A --> B
    end
    @test_nowarn @macroexpand @reaction_network begin
        @require_declaration
        @parameters n k
        @species A(t) B(t)
        k*n, A --> B
    end

    # Test error when a parameter in stoichiometry is inferred
    @test_throws UndeclaredSymbolicError @macroexpand @reaction_network begin
        @require_declaration
        @parameters k
        @species A(t) B(t)
        k, n*A --> B
    end
    @test_nowarn @macroexpand @reaction_network begin
        @require_declaration
        @parameters k n
        @species A(t) B(t)
        k, n*A --> B
    end

    # Test error when a variable in an equation is inferred
    @test_throws UndeclaredSymbolicError @macroexpand @reaction_network begin
        @require_declaration
        @equations V ~ V^2 + 2
    end
    @test_nowarn @macroexpand @reaction_network begin
        @require_declaration
        @variables V(t)
        @equations V ~ V^2 + 2
    end

    # Test error when a variable in an observable is inferred
    @test_throws UndeclaredSymbolicError @macroexpand @reaction_network begin
        @require_declaration
        @variables X1(t)
        @observables X2 ~ X1
    end
    @test_nowarn @macroexpand @reaction_network begin
        @require_declaration
        @variables X1(t) X2(t)
        @observables X2 ~ X1
    end

    # Test when the default differential D is inferred
    @test_throws UndeclaredSymbolicError @macroexpand @reaction_network begin
        @require_declaration
        @variables V(t)
        @equations D(V) ~ 1 - V
    end
    @test_nowarn @macroexpand @reaction_network begin
        @differentials D = Differential(t)
        @variables X1(t) X2(t)
        @observables X2 ~ X1
    end
end
