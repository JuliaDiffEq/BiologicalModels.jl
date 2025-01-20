### Prepares Tests ###

# Fetch packages.
using DiffEqBase, Catalyst, Random, Test
using ModelingToolkit: operation, get_unknowns, get_ps, get_eqs, get_systems,
                       get_iv, nameof
using Symbolics: iscall

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)

# Sets the default `t` to use.
t = default_t()

# Fetch test networks and functions.
include("../test_networks.jl")
include("../test_functions.jl")

### Declares Testing Functions ###

function unpacksys(sys)
    get_eqs(sys), get_iv(sys), get_unknowns(sys), get_ps(sys), nameof(sys), get_systems(sys)
end

opname(x) = iscall(x) ? nameof(operation(x)) : nameof(x)
alleq(xs, ys) = all(isequal(x, y) for (x, y) in zip(xs, ys))

# Gets all the reactants in a set of equations.
function all_reactants(eqs)
    all_reactants = []
    for eq in eqs
        append!(all_reactants, opname.(eq.substrates))
        append!(all_reactants, opname.(eq.products))
    end
    return Set{Symbol}(unique(all_reactants))
end

# Gets all parameters (where every reaction rate is constant).
function all_parameters(eqs)
    return Set(unique(map(eq -> opname(eq.rate), eqs)))
end

# Perform basic tests.
function basic_test(rn, N, unknowns_syms, p_syms)
    eqs, iv, unknowns, ps, name, systems = unpacksys(rn)
    @test length(eqs) == N
    @test opname(iv) == :t
    @test length(unknowns) == length(unknowns_syms)
    @test issetequal(map(opname, unknowns), unknowns_syms)
    @test all_reactants(eqs) == Set(unknowns_syms)
    @test length(ps) == length(p_syms)
    @test issetequal(map(opname, ps), p_syms)
end

## Run Tests ###

# Tests the various network constructors. Test for `@network_component` and `@network_component`.
# Tests for combinations of reactions/no reactions, no name/name/interpolated name.
let
    # Declare comparison networks programmatically.
    @parameters d
    @species X(t)
    rx = Reaction(d, [X], [])

    rs_empty = ReactionSystem([], t; name = :name)
    rs = ReactionSystem([rx], t; name = :name)
    rs_empty_comp = complete(rs_empty)
    rs_comp = complete(rs)

    # Declare empty networks.
    name_sym = :name
    rs_empty_1 = @network_component
    rs_empty_2 = @network_component name
    rs_empty_3 = @network_component $name_sym
    rs_empty_comp_1 = @reaction_network
    rs_empty_comp_2 = @reaction_network name
    rs_empty_comp_3 = @reaction_network $name_sym

    # Check that empty networks are correct.
    isequivalent(rs_empty_1, rs_empty)
    rs_empty_2 == rs_empty
    rs_empty_3 == rs_empty
    isequivalent(rs_empty_comp_1, rs_empty_comp)
    rs_empty_comp_2 == rs_empty_comp
    rs_empty_comp_3 == rs_empty_comp

    # Declare non-empty networks.
    rs_1 = @network_component begin
        d, X --> 0
    end
    rs_2 = @network_component name begin
        d, X --> 0
    end
    rs_3 = @network_component $name_sym begin
        d, X --> 0
    end
    rs_comp_1 = @reaction_network begin
        d, X --> 0
    end
    rs_comp_2 = @reaction_network name begin
        d, X --> 0
    end
    rs_comp_3 = @reaction_network $name_sym begin
        d, X --> 0
    end

    # Check that non-empty networks are correct.
    isequivalent(rs_1, rs)
    rs_2 == rs
    rs_3 == rs
    isequivalent(rs_empty_1, rs_empty)
    rs_empty_2 == rs_empty
    rs_empty_3 == rs_empty
end

# Test basic properties of networks.
let
    basic_test(reaction_networks_standard[1], 10, [:X1, :X2, :X3],
               [:p1, :p2, :p3, :k1, :k2, :k3, :k4, :d1, :d2, :d3])
    @test all_parameters(get_eqs(reaction_networks_standard[1])) ==
          Set([:p1, :p2, :p3, :k1, :k2, :k3, :k4, :d1, :d2, :d3])
    basic_test(reaction_networks_standard[2], 3, [:X1, :X2], [:v1, :K1, :v2, :K2, :d])
    basic_test(reaction_networks_standard[3], 10, [:X1, :X2, :X3, :X4],
               [:v1, :K1, :v2, :K2, :k1, :k2, :k3, :k4, :d])
    basic_test(reaction_networks_standard[4], 8, [:X1, :X2, :X3, :X4],
               [:v1, :K1, :v2, :K2, :v3, :K3, :v4, :K4, :d1, :d2, :d3, :d4])
    basic_test(reaction_networks_standard[5], 8, [:X1, :X2, :X3, :X4],
               [:p, :k1, :k2, :k3, :k4, :k5, :k6, :d])
    @test all_parameters(get_eqs(reaction_networks_standard[5])) ==
          Set([:p, :k1, :k2, :k3, :k4, :k5, :k6, :d])
    basic_test(reaction_networks_hill[1], 4, [:X1, :X2],
               [:v1, :v2, :K1, :K2, :n1, :n2, :d1, :d2])
    basic_test(reaction_networks_conserved[1], 6, [:X1, :X2, :X3],
               [:k1, :k2, :k3, :k4, :k5, :k6])
    basic_test(reaction_networks_real[1], 4, [:X, :Y], [:A, :B])
    basic_test(reaction_networks_weird[1], 2, [:X], [:p, :d])
    basic_test(reaction_networks_weird[2], 4, [:X, :Y, :Z], [:k1, :k2, :k3, :k4])
end

# Compares networks to networks created using different arrow types.
let
    identical_networks_1 = Vector{Pair}()

    different_arrow_1 = @reaction_network begin
        (p1, p2, p3), ∅ > (X1, X2, X3)
        (k1, k2), X2 ↔ X1 + 2X3
        (k3, k4), X1 ⟷ X3
        (d1, d2, d3), (X1, X2, X3) → ∅
    end
    push!(identical_networks_1, reaction_networks_standard[1] => different_arrow_1)

    different_arrow_2 = @reaction_network begin
        mmr(X2, v1, K1), ∅ → X1
        mm(X1, v2, K2), ∅ ↣ X2
        d, X1 + X2 ↦ ∅
    end
    push!(identical_networks_1, reaction_networks_standard[2] => different_arrow_2)

    different_arrow_3 = @reaction_network begin
        mm(X2, v1, K1), ∅ ⇾ X1
        mm(X3, v2, K2), ∅ ⟶ X2
        (k1, k2), X1 ⇄ X3
        (k3, k4), X3 + X2 ⇆ X4 + X1
        d, (X1, X2, X3, X4) ⟼ ∅
    end
    push!(identical_networks_1, reaction_networks_standard[3] => different_arrow_3)

    different_arrow_4 = @reaction_network begin
        mmr(X4, v1, K1), ∅ ⥟ X1
        mmr(X1, v2, K2), ∅ ⥟ X2
        mmr(X2, v3, K3), ∅ ⇀ X3
        mmr(X3, v4, K4), ∅ ⇁ X4
        (d1, d2, d3, d4), (X1, X2, X3, X4) --> ∅
    end
    push!(identical_networks_1, reaction_networks_standard[4] => different_arrow_4)

    # Yes the name is different, I wanted one with several single direction arrows.
    different_arrow_8 = @reaction_network begin
        p, 2X1 < ∅
        k1, X2 ← X1
        (k2, k3), X3 ⟻ X2
        d, ∅ ↼ X3
    end
    push!(identical_networks_1, reaction_networks_standard[8] => different_arrow_8)

    for networks in identical_networks_1
        for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
            u0 = rnd_u0(networks[1], rng; factor)
            p = rnd_ps(networks[1], rng; factor)
            t = rand(rng)

            @test f_eval(networks[1], u0, p, t) ≈ f_eval(networks[2], u0, p, t)
            @test jac_eval(networks[1], u0, p, t) ≈ jac_eval(networks[2], u0, p, t)
            @test g_eval(networks[1], u0, p, t) ≈ g_eval(networks[2], u0, p, t)
        end
    end
end

# Compares simulations for network with different species and parameter names
let
    # Fetches the original network, and also declares it using alternative notation.
    network = reaction_networks_standard[5]
    differently_written_5 = @reaction_network begin
        q, ∅ → Y1
        (l1, l2), Y1 ⟷ Y2
        (l3, l4), Y2 ⟷ Y3
        (l5, l6), Y3 ⟷ Y4
        c, Y4 → ∅
    end

    # Checks that the networks' functions evaluates equally for various randomised inputs.
    @unpack X1, X2, X3, X4, p, d, k1, k2, k3, k4, k5, k6 = network
    for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
        u0_1 = Dict(rnd_u0(network, rng; factor))
        p_1 = Dict(rnd_ps(network, rng; factor))
        u0_2 = [:Y1 => u0_1[X1], :Y2 => u0_1[X2], :Y3 => u0_1[X3], :Y4 => u0_1[X4]]
        p_2 = [:q => p_1[p], :c => p_1[d], :l1 => p_1[k1], :l2 => p_1[k2], :l3 => p_1[k3],
               :l4 => p_1[k4], :l5 => p_1[k5], :l6 => p_1[k6]]
        t = rand(rng)

        @test f_eval(network, u0_1, p_1, t) ≈ f_eval(differently_written_5, u0_2, p_2, t)
        @test jac_eval(network, u0_1, p_1, t) ≈ jac_eval(differently_written_5, u0_2, p_2, t)
        @test g_eval(network, u0_1, p_1, t) ≈ g_eval(differently_written_5, u0_2, p_2, t)
    end
end

# Compares networks to networks written in different ways.
let
    identical_networks_2 = Vector{Pair}()

    # Unfold reactions.
    differently_written_6 = @reaction_network begin
        p1, ∅ → X1
        p2, ∅ → X2
        k1, 2X1 → X3
        k2, X3 → 2X1
        k3, X2 + X3 → 4X4
        k4, 4X4 → X2 + X3
        k5, X4 + X1 → 2X3
        k6, 2X3 → X4 + X1
        d, X1 → ∅
        d, X2 → ∅
        d, X3 → ∅
        d, X4 → ∅
    end
    push!(identical_networks_2, reaction_networks_standard[6] => differently_written_6)

    # Ignore mass action.
    differently_written_7 = @reaction_network begin
        @parameters p1 p2 p3 k1 k2 k3 v1 K1 d1 d2 d3 d4 d5
        (p1, p2, p3), ∅ ⇒ (X1, X2, X3)
        (k1 * X1 * X2^2 / 2, k2 * X4), X1 + 2X2 ⟺ X4
        (mm(X3, v1, K1) * X4, k3 * X5), X4 ⇔ X5
        (d1 * X1, d2 * X2, d3 * X3, d4 * X4, d5 * X5), ∅ ⟽ (X1, X2, X3, X4, X5)
    end
    push!(identical_networks_2, reaction_networks_standard[7] => differently_written_7)

    # Ignore mass action new arrows.
    differently_written_8 = @reaction_network begin
        @parameters p1 p2 p3 k1 k2 k3 v1 K1 d1 d2 d3 d4 d5
        (p1, p2, p3), ∅ => (X1, X2, X3)
        (k1 * X1 * X2^2 / 2, k2 * X4), X1 + 2X2 ⟺ X4
        (mm(X3, v1, K1) * X4, k3 * X5), X4 ⇔ X5
        (d1 * X1, d2 * X2, d3 * X3, d4 * X4, d5 * X5), ∅ <= (X1, X2, X3, X4, X5)
    end
    push!(identical_networks_2, reaction_networks_standard[7] => differently_written_8)

    for networks in identical_networks_2
        for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
            u0 = rnd_u0(networks[1], rng; factor)
            p = rnd_ps(networks[1], rng; factor)
            t = rand(rng)

            @test f_eval(networks[1], u0, p, t) ≈ f_eval(networks[2], u0, p, t)
            @test jac_eval(networks[1], u0, p, t) ≈ jac_eval(networks[2], u0, p, t)
            @test g_eval(networks[1], u0, p, t) ≈ g_eval(networks[2], u0, p, t)
        end
    end
end

# Compares networks to networks written without parameters.
let
    identical_networks_3 = Vector{Pair}()
    parameter_sets = []

    # Different parameter and variable names.
    no_parameters_9 = @reaction_network begin
        (1.5, 1, 2), ∅ ⟶ (X1, X2, X3)
        (0.01, 2.3, 1001), (X1, X2, X3) ⟶ ∅
        (π, 42), X1 + X2 ⟷ X3
        (19.9, 999.99), X3 ⟷ X4
        (sqrt(3.7), exp(1.9)), X4 ⟷ X1 + X2
    end
    push!(identical_networks_3, reaction_networks_standard[9] => no_parameters_9)
    push!(parameter_sets, [:p1 => 1.5, :p2 => 1, :p3 => 2, :d1 => 0.01, :d2 => 2.3, :d3 => 1001,
                           :k1 => π, :k2 => 42, :k3 => 19.9, :k4 => 999.99, :k5 => sqrt(3.7), :k6 => exp(1.9)])

    no_parameters_10 = @reaction_network begin
        0.01, ∅ ⟶ X1
        (3.1, 3.2), X1 → X2
        (0.0, 2.1), X2 → X3
        (901.0, 63.5), X3 → X4
        (7, 8), X4 → X5
        1.0, X5 ⟶ ∅
    end
    push!(identical_networks_3, reaction_networks_standard[10] => no_parameters_10)
    push!(parameter_sets, [:p => 0.01, :k1 => 3.1, :k2 => 3.2, :k3 => 0.0, :k4 => 2.1, :k5 => 901.0,
                           :k6 => 63.5, :k7 => 7, :k8 => 8, :d => 1.0])

    for (networks, p_1) in zip(identical_networks_3, parameter_sets)
        for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
            u0 = rnd_u0(networks[1], rng; factor)
            t = rand(rng)

            @test f_eval(networks[1], u0, p_1, t) ≈ f_eval(networks[2], u0, [], t)
            @test jac_eval(networks[1], u0, p_1, t) ≈ jac_eval(networks[2], u0, [], t)
            @test g_eval(networks[1], u0, p_1, t) ≈ g_eval(networks[2], u0, [], t)
        end
    end
end

# Tests that Reaction System created manually and through macro are identical.
let
    identical_networks_4 = Vector{Pair}()
    @parameters v1 K1 v2 K2 k1 k2 k3 k4 k5 p d t
    @species X1(t) X2(t) X3(t) X4(t) X5(t)

    rxs_1 = [Reaction(p, nothing, [X1], nothing, [2]),
        Reaction(k1, [X1], [X2], [1], [1]),
        Reaction(k2, [X2], [X3], [1], [1]),
        Reaction(k3, [X2], [X3], [1], [1]),
        Reaction(d, [X3], nothing, [1], nothing)]
    @named rs_1 = ReactionSystem(rxs_1, t, [X1, X2, X3], [p, k1, k2, k3, d])
    push!(identical_networks_4, reaction_networks_standard[8] => rs_1)

    rxs_2 = [Reaction(k1, [X1], [X2], [1], [1]),
        Reaction(k2 * X5, [X2], [X1], [1], [1]),
        Reaction(k3 * X5, [X3], [X4], [1], [1]),
        Reaction(k4, [X4], [X3], [1], [1]),
        Reaction(p + k5 * X2 * X3, nothing, [X5], nothing, [1]),
        Reaction(d, [X5], nothing, [1], nothing)]
    @named rs_2 = ReactionSystem(rxs_2, t, [X1, X2, X3, X4, X5], [k1, k2, k3, k4, p, k5, d])
    push!(identical_networks_4, reaction_networks_conserved[3] => rs_2)

    rxs_3 = [Reaction(k1, [X1], [X2], [1], [1]),
        Reaction(0, [X2], [X3], [1], [1]),
        Reaction(k2, [X3], [X4], [1], [1]),
        Reaction(k3, [X4], [X5], [1], [1])]
    @named rs_3 = ReactionSystem(rxs_3, t, [X1, X2, X3, X4, X5], [k1, k2, k3])
    push!(identical_networks_4, reaction_networks_weird[7] => rs_3)

    for networks in identical_networks_4
        @test isequal(get_iv(networks[1]), get_iv(networks[2]))
        @test alleq(get_unknowns(networks[1]), get_unknowns(networks[2]))
        @test alleq(get_ps(networks[1]), get_ps(networks[2]))
        @test ModelingToolkit.get_systems(networks[1]) ==
              ModelingToolkit.get_systems(networks[2])
        @test length(get_eqs(networks[1])) == length(get_eqs(networks[2]))
        for (e1, e2) in zip(get_eqs(networks[1]), get_eqs(networks[2]))
            @test isequal(e1.rate, e2.rate)
            @test isequal(e1.substrates, e2.substrates)
            @test isequal(e1.products, e2.products)
            @test isequal(e1.substoich, e2.substoich)
            @test isequal(e1.prodstoich, e2.prodstoich)
            @test isequal(e1.netstoich, e2.netstoich)
            @test isequal(e1.only_use_rate, e2.only_use_rate)
        end
    end
end

### Tests Usage of Various Symbols ###

# Tests that time is handled properly.
let
    time_network = @reaction_network begin
        (t, k2), X1 ↔ X2
        (k3, t), X2 ↔ X3
        (t, k6), X3 ↔ X1
    end

    for factor in [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]
        τ = rand(rng)
        u = rnd_u0(reaction_networks_conserved[1], rng; factor)
        p_2 = rnd_ps(time_network, rng; factor)
        p_1 = [p_2; reaction_networks_conserved[1].k1 => τ;
               reaction_networks_conserved[1].k4 => τ; reaction_networks_conserved[1].k5 => τ]

        @test f_eval(reaction_networks_conserved[1], u, p_1, τ) ≈ f_eval(time_network, u, p_2, τ)
        @test jac_eval(reaction_networks_conserved[1], u, p_1, τ) ≈ jac_eval(time_network, u, p_2, τ)
        @test g_eval(reaction_networks_conserved[1], u, p_1, τ) ≈ g_eval(time_network, u, p_2, τ)
    end
end

# Check that various symbols can be used as species/parameter names.
let
    @reaction_network begin
        (a, A), n ⟷ N
        (b, B), o ⟷ O
        (c, C), p ⟷ P
        (d, D), q ⟷ Q
        (e, E), r ⟷ R
        (f, F), s ⟷ S
        (g, G), u ⟷ U
        (h, H), v ⟷ V
        (j, J), w ⟷ W
        (k, K), x ⟷ X
        (l, L), y ⟷ Y
        (m, M), z ⟷ Z
    end

    @reaction_network begin (1.0, 1.0), i ⟷ T end

    @reaction_network begin
        (å, Å), ü ⟷ Ü
        (ä, Ä), ñ ⟷ Ñ
        (ö, Ö), æ ⟷ Æ
    end

    @reaction_network begin
        (α, Α), ν ⟷ Ν
        (β, Β), ξ ⟷ Ξ
        (γ, γ), ο ⟷ Ο
        (δ, Δ), Π ⟷ Π
        (ϵ, Ε), ρ ⟷ Ρ
        (ζ, Ζ), σ ⟷ Σ
        (η, Η), τ ⟷ Τ
        (θ, Θ), υ ⟷ Υ
        (ι, Ι), ϕ ⟷ Φ
        (κ, κ), χ ⟷ Χ
        (λ, Λ), ψ ↔ Ψ
        (μ, Μ), ω ⟷ Ω
    end
end

# Test that the `I` symbol works as a quantity name.
let
    rn = @reaction_network begin
        k1, S + I --> 2I
        k2, I --> R
    end
    @species I(t)
    @test any(isequal(I), species(rn))
    @test any(isequal(I), unknowns(rn))
end

# Tests backwards and bi-directional arrows.
let
    rn1 = @reaction_network arrowtest begin
        (a1, a2), C <--> 0
        (k1, k2), A + B <--> C
        b1, 0 <-- B
    end

    rn2 = @reaction_network arrowtest begin
        a1, C --> 0
        a2, 0 --> C
        k1, A + B --> C
        k2, C --> A + B
        b1, B --> 0
    end

    @test rn1 == rn2
end

# Tests arrow variants in `@reaction` macro.
let
    @test isequal((@reaction k, 0 --> X), (@reaction k, X <-- 0))
    @test isequal((@reaction k, 0 --> X), (@reaction k, X ⟻ 0))
    @test isequal((@reaction k, 0 --> X), (@reaction k, 0 → X))
    @test isequal((@reaction k, 0 --> X), (@reaction k, 0 ⥟ X))
end

# Test that symbols with special meanings are handled properly.
let 
    test_network = @reaction_network begin t * k, X --> ∅ end
    @test length(species(test_network)) == 1
    @test length(parameters(test_network)) == 1

    test_network = @reaction_network begin π, X --> ∅ end
    @test length(species(test_network)) == 1
    @test length(parameters(test_network)) == 0
    @test reactions(test_network)[1].rate == π

    test_network = @reaction_network begin pi, X --> ∅ end
    @test length(species(test_network)) == 1
    @test length(parameters(test_network)) == 0
    @test reactions(test_network)[1].rate == pi

    test_network = @reaction_network begin ℯ, X --> ∅ end
    @test length(species(test_network)) == 1
    @test length(parameters(test_network)) == 0
    @test reactions(test_network)[1].rate == ℯ
end

### Error Test ###

# Erroneous `@reaction` usage.
let
    # Bi-directional reaction using the `@reaction` macro.
    @test_throws Exception @eval @reaction (k1,k2), X1 <--> X2

    # Bundles reactions.
    @test_throws Exception @eval @reaction k, (X1,X2) --> 0
end

# Tests that malformed reactions yields errors.
let
    # Checks that malformed combinations of entries yields errors.
    @test_throws Exception @eval @reaction_network begin
        d, X --> 0, Y --> 0
    end
    @test_throws Exception @eval @reaction_network begin
        d, X --> 0, [misc="Ok metadata"], [description="Metadata in (erroneously) extra []."]
    end

    # Checks that incorrect bundling yields error.
    @test_throws Exception @eval @reaction_network begin
        (k1,k2,k3), (X1,X2) --> 0
    end

    # Checks that incorrect stoichiometric expression yields an error.
    @test_throws Exception @eval @reaction_network begin
        k, X^Y --> XY
    end
end

# Check that forbidden symbols correctly generates errors.
let
    # @reaction macro, symbols that cannot be in the rate.
    @test_throws Exception @eval @reaction im, 0 --> X
    @test_throws Exception @eval @reaction nothing, 0 --> X
    @test_throws Exception @eval @reaction Γ, 0 --> X
    @test_throws Exception @eval @reaction ∅, 0 --> X
    
    # @reaction macro, symbols that cannot be a reactant.
    @test_throws Exception @eval @reaction 1, 0 --> im
    @test_throws Exception @eval @reaction 1, 0 --> nothing
    @test_throws Exception @eval @reaction 1, 0 --> Γ
    @test_throws Exception @eval @reaction 1, 0 --> ℯ
    @test_throws Exception @eval @reaction 1, 0 --> pi
    @test_throws Exception @eval @reaction 1, 0 --> π
    @test_throws Exception @eval @reaction 1, 0 --> t
    
    # @reaction_network macro, symbols that cannot be in the rate.
    @test_throws Exception @eval @reaction_network begin im, 0 --> X end
    @test_throws Exception @eval @reaction_network begin nothing, 0 --> X end
    @test_throws Exception @eval @reaction_network begin Γ, 0 --> X end
    @test_throws Exception @eval @reaction_network begin ∅, 0 --> X end
    
    # @reaction_network macro, symbols that cannot be a reactant.
    @test_throws Exception @eval @reaction_network begin 1, 0 --> im end
    @test_throws Exception @eval @reaction_network begin 1, 0 --> nothing end
    @test_throws Exception @eval @reaction_network begin 1, 0 --> Γ end
    @test_throws Exception @eval @reaction_network begin 1, 0 --> ℯ end
    @test_throws Exception @eval @reaction_network begin 1, 0 --> pi end
    @test_throws Exception @eval @reaction_network begin 1, 0 --> π end
    @test_throws Exception @eval @reaction_network begin 1, 0 --> t end

    # Checks that non-supported arrow type usage yields error.
    @test_throws Exception @eval @reaction_network begin
        d, X ⇻ 0
    end
end