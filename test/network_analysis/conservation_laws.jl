### Prepares Tests ###

# Fetch packages.
using Catalyst, JumpProcesses, LinearAlgebra, NonlinearSolve, OrdinaryDiffEqTsit5, OrdinaryDiffEqDefault, OrdinaryDiffEqVerner, SteadyStateDiffEq, StochasticDiffEq, Test

# Sets stable rng number.
using StableRNGs
rng = StableRNG(123456)
seed = rand(rng, 1:100)

# Fetch test networks.
include("../test_networks.jl")

# Except where we test the warnings, we do not want to print this warning.
remove_conserved_warn = false

### Basic Tests ###

# Tests basic functionality on system with known conservation laws.
let
    rn = @reaction_network begin
        (1, 2), A + B <--> C
        (3, 2), D <--> E
        (0.1, 0.2), E <--> F
        6, F --> G
        7, H --> G
        5, 0 --> K
        3, A + B --> Z + C
    end

    S = netstoichmat(rn)
    C = conservationlaws(S)
    @test size(C, 1) == 3
    b = [0, 0, 0, 1, 1, 1, 1, 1, 0, 0]
    @test any(b == C[i, :] for i in 1:size(C, 1))

    # For the A + B <--> C subsystem one of these must occur
    # as a conservation law.
    D = [1 -1 0 0 0 0 0 0 0 0;
         -1 1 0 0 0 0 0 0 0 0
         0 1 1 0 0 0 0 0 0 0]
    @test any(D[j, :] == C[i, :] for i in 1:size(C, 1), j in 1:size(D, 1))

    C = conservationlaws(rn)
    @test size(C, 1) == 3
    @test Catalyst.get_networkproperties(rn).nullity == 3
    @test any(b == C[i, :] for i in 1:size(C, 1))
    @test any(D[j, :] == C[i, :] for i in 1:size(C, 1), j in 1:size(D, 1))
end

# Tests conservation law computation on large number of networks where we know which have conservation laws.
let
    # networks for which we know there is no conservation laws.
    Cs_standard = map(conservationlaws, reaction_networks_standard)
    Cs_hill = map(conservationlaws, reaction_networks_hill)
    @test all(size(C, 1) == 0 for C in Cs_standard)
    @test all(size(C, 1) == 0 for C in Cs_hill)

    # Networks for which there are known conservation laws (stored in `reaction_network_conslaws`).
    function consequiv(A, B)
        rank([A; B]) == rank(A) == rank(B)
    end
    Cs_constraint = map(conservationlaws, reaction_networks_conserved)
    @test all(consequiv.(Matrix{Int}.(Cs_constraint), reaction_network_conslaws))
end

# Tests additional conservation law-related functions.
let
    rn = @reaction_network begin
        (k1, k2), X1 <--> X2
        (k3, k4), X3 <--> X4
    end
    cons_laws = conservationlaws(rn)
    cons_eqs = conservedequations(rn)
    cons_laws_constants = conservationlaw_constants(rn)
    conserved_quantity = conservedquantities(cons_laws[1, :], unknowns(rn)[1])

    @test sum(cons_laws) == 4
    @test size(cons_laws) == (2, 4)
    @test length(cons_eqs) == 2
    @test length(conserved_quantity) == 4
    @test length(cons_laws_constants) == 2
    @test count(isequal.(conserved_quantity, Num(0))) == 2
end

# Tests that `conservationlaws`'s caches something.
let
    # Creates network with/without cached conservation laws.
    rn = @reaction_network rn begin
        (k1,k2), X1 <--> X2
    end
    rn_cached = deepcopy(rn)
    conservationlaws(rn_cached)

    # Checks that equality is correct (currently equality does not consider network property caching).
    @test rn_cached == rn
    @test Catalyst.get_networkproperties(rn_cached) != Catalyst.get_networkproperties(rn)
end

### Simulation & Solving Tests ###

# Test conservation law elimination.
let
    # Declares the model
    rn = @reaction_network begin
        (k1, k2), A + B <--> C
        (m1, m2), D <--> E
        b12, F1 --> F2
        b23, F2 --> F3
        b31, F3 --> F1
    end
    @unpack A, B, C, D, E, F1, F2, F3, k1, k2, m1, m2, b12, b23, b31 = rn
    sps = species(rn)
    u0 = [A => 10.0, B => 10.0, C => 0.0, D => 10.0, E => 0.0, F1 => 8.0, F2 => 0.0,
        F3 => 0.0]
    p = [k1 => 1.0, k2 => 0.1, m1 => 1.0, m2 => 2.0, b12 => 1.0, b23 => 2.0, b31 => 0.1]
    tspan = (0.0, 20.0)

    # Simulates model using ODEs and checks that simulations are identical.
    osys = complete(convert(ODESystem, rn; remove_conserved = true, remove_conserved_warn))
    oprob1 = ODEProblem(osys, u0, tspan, p)
    oprob2 = ODEProblem(rn, u0, tspan, p)
    oprob3 = ODEProblem(rn, u0, tspan, p; remove_conserved = true, remove_conserved_warn)
    osol1 = solve(oprob1, Tsit5(); abstol = 1e-8, reltol = 1e-8, saveat= 0.2)
    osol2 = solve(oprob2, Tsit5(); abstol = 1e-8, reltol = 1e-8, saveat= 0.2)
    osol3 = solve(oprob3, Tsit5(); abstol = 1e-8, reltol = 1e-8, saveat= 0.2)
    @test osol1[sps] ≈ osol2[sps] ≈ osol3[sps]

    # Checks that steady states found using nonlinear solving and steady state simulations are identical.
    nsys = complete(convert(NonlinearSystem, rn; remove_conserved = true, remove_conserved_warn))
    nprob1 = NonlinearProblem{true}(nsys, u0, p)
    nprob2 = NonlinearProblem(rn, u0, p)
    nprob3 = NonlinearProblem(rn, u0, p; remove_conserved = true, remove_conserved_warn)
    ssprob1 = SteadyStateProblem{true}(osys, u0, p)
    ssprob2 = SteadyStateProblem(rn, u0, p)
    ssprob3 = SteadyStateProblem(rn, u0, p; remove_conserved = true, remove_conserved_warn)
    nsol1 = solve(nprob1, NewtonRaphson(); abstol = 1e-8)
    # Nonlinear problems cannot find steady states properly without removing conserved species.
    nsol3 = solve(nprob3, NewtonRaphson(); abstol = 1e-8)
    sssol1 = solve(ssprob1, DynamicSS(Tsit5()); abstol = 1e-8, reltol = 1e-8)
    sssol2 = solve(ssprob2, DynamicSS(Tsit5()); abstol = 1e-8, reltol = 1e-8)
    sssol3 = solve(ssprob3, DynamicSS(Tsit5()); abstol = 1e-8, reltol = 1e-8)
    @test nsol1[sps] ≈ nsol3[sps] ≈ sssol1[sps] ≈ sssol2[sps] ≈ sssol3[sps]

    # Creates SDEProblems using various approaches.
    u0_sde = [A => 100.0, B => 20.0, C => 5.0, D => 10.0, E => 3.0, F1 => 8.0, F2 => 2.0,
        F3 => 20.0]
    ssys = complete(convert(SDESystem, rn; remove_conserved = true, remove_conserved_warn))
    sprob1 = SDEProblem(ssys, u0_sde, tspan, p)
    sprob2 = SDEProblem(rn, u0_sde, tspan, p)
    sprob3 = SDEProblem(rn, u0_sde, tspan, p; remove_conserved = true, remove_conserved_warn)

    # Checks that the SDEs f and g function evaluates to the same thing.
    ind_us = ModelingToolkit.get_unknowns(ssys)
    us = ModelingToolkit.get_unknowns(rn)
    ind_uidxs = findall(in(ind_us), us)
    u1 = copy(sprob1.u0)
    u2 = sprob2.u0
    u3 = copy(sprob3.u0)
    du1 = similar(u1)
    du2 = similar(u2)
    du3 = similar(u3)
    g1 = zeros(length(u1), numreactions(rn))
    g2 = zeros(length(u2), numreactions(rn))
    g3 = zeros(length(u3), numreactions(rn))
    sprob1.f(du1, sprob1.u0, sprob1.p, 1.0)
    sprob2.f(du2, sprob2.u0, sprob2.p, 1.0)
    sprob3.f(du3, sprob3.u0, sprob3.p, 1.0)
    @test du1 ≈ du2[ind_uidxs] ≈ du3
    sprob1.g(g1, sprob1.u0, sprob1.p, 1.0)
    sprob2.g(g2, sprob2.u0, sprob2.p, 1.0)
    sprob3.g(g3, sprob3.u0, sprob3.p, 1.0)
    @test g1 ≈ g2[ind_uidxs, :] ≈ g3
end

# Tests simulations for various input types (using X, rn.X, and :X forms).
# Tests that conservation laws can be generated for system with non-default parameter types.
let
    # Prepares the model.
    rn = @reaction_network rn begin
        @parameters kB::Int64
        (kB,kD), X + Y <--> XY
    end
    sps = species(rn)
    @unpack kB, kD, X, Y, XY = rn

    # Creates `ODEProblem`s using three types of inputs. Checks that solutions are identical.
    u0_1 = [X => 2.0, Y => 3.0, XY => 4.0]
    u0_2 = [rn.X => 2.0, rn.Y => 3.0, rn.XY => 4.0]
    u0_3 = [:X => 2.0, :Y => 3.0, :XY => 4.0]
    ps = (kB => 2, kD => 1.5)
    oprob1 = ODEProblem(rn, u0_1, 10.0, ps; remove_conserved = true, remove_conserved_warn)
    oprob2 = ODEProblem(rn, u0_2, 10.0, ps; remove_conserved = true, remove_conserved_warn)
    oprob3 = ODEProblem(rn, u0_3, 10.0, ps; remove_conserved = true, remove_conserved_warn)
    @test solve(oprob1)[sps] ≈ solve(oprob2)[sps] ≈ solve(oprob3)[sps]
end

# Tests conservation laws in SDE simulation.
let
    # Creates `SDEProblem`s.
    rn = @reaction_network begin
        (k1,k2), X1 <--> X2
    end
    u0 = Dict([:X1 => 100.0, :X2 => 120.0])
    ps = [:k1 => 0.2, :k2 => 0.15]
    sprob = SDEProblem(rn, u0, 10.0, ps; remove_conserved = true, remove_conserved_warn)

    # Checks that conservation laws hold in all simulations.
    sol = solve(sprob, ImplicitEM(); seed)
    @test sol[:X1] + sol[:X2] ≈ sol[rn.X1 + rn.X2] ≈ fill(u0[:X1] + u0[:X2], length(sol.t))
end

# Checks that the conservation law parameter's value can be changed in simulations.
let
    # Prepares `ODEProblem`s.
    rn = @reaction_network begin
        (k1,k2), X1 <--> X2
    end
    osys = complete(convert(ODESystem, rn; remove_conserved = true, remove_conserved_warn))
    u0 = [osys.X1 => 1.0, osys.X2 => 1.0]
    ps_1 = [osys.k1 => 2.0, osys.k2 => 3.0]
    ps_2 = [osys.k1 => 2.0, osys.k2 => 3.0, osys.Γ[1] => 4.0]
    oprob1 = ODEProblem(osys, u0, 10.0, ps_1)
    oprob2 = ODEProblem(osys, u0, 10.0, ps_2)

    # Checks that the solutions generates the correct conserved quantities.
    sol1 = solve(oprob1; saveat = 1.0)
    sol2 = solve(oprob2; saveat = 1.0)
    @test all(sol1[osys.X1 + osys.X2] .== 2.0)
    @test all(sol2[osys.X1 + osys.X2] .== 4.0)
end

# Tests system problem updating when conservation laws are eliminated.
# Checks that the correct values are used after the conservation law species are updated.
# Here is an issue related to the broken tests: https://github.com/SciML/Catalyst.jl/issues/952
let
    # Create model and fetch the conservation parameter (Γ).
    t = default_t()
    @parameters k1 k2
    @species X1(t) X2(t)
    rxs = [
        Reaction(k1, [X1], [X2]),
        Reaction(k2, [X2], [X1])
    ]
    @named rs = ReactionSystem(rxs, t)
    osys = convert(ODESystem, complete(rs); remove_conserved = true, remove_conserved_warn = false)
    osys = complete(osys)
    @unpack Γ = osys

    # Creates an `ODEProblem`.
    u0 = [X1 => 1.0, X2 => 2.0]
    ps = [k1 => 0.1, k2 => 0.2]
    oprob = ODEProblem(osys, u0, (0.0, 1.0), ps)

    # Check `ODEProblem` content.
    @test oprob[X1] == 1.0
    @test oprob[X2] == 2.0
    @test oprob.ps[k1] == 0.1
    @test oprob.ps[k2] == 0.2
    @test oprob.ps[Γ[1]] == 3.0

    # Currently, any kind of updating of species or the conservation parameter(s) is not possible.

    # Update problem parameters using `remake`.
    oprob_new = remake(oprob; p = [k1 => 0.3, k2 => 0.4])
    @test oprob_new.ps[k1] == 0.3
    @test oprob_new.ps[k2] == 0.4
    integrator = init(oprob_new, Tsit5())
    @test integrator.ps[k1] == 0.3
    @test integrator.ps[k2] == 0.4

    # Update problem parameters using direct indexing.
    oprob.ps[k1] = 0.5
    oprob.ps[k2] = 0.6
    @test oprob.ps[k1] == 0.5
    @test oprob.ps[k2] == 0.6
    integrator = init(oprob, Tsit5())
    @test integrator.ps[k1] == 0.5
    @test integrator.ps[k2] == 0.6
end

### Other Tests ###

# Checks that `JumpSystem`s with conservation laws cannot be generated.
let
    rn = @reaction_network begin
        (k1,k2), X1 <--> X2
    end
    @test_throws ArgumentError convert(JumpSystem, rn; remove_conserved = true, remove_conserved_warn)
end

# Checks that `conserved` metadata is added correctly to parameters.
# Checks that the `isconserved` getter function works correctly.
let
    # Creates ODESystem with conserved quantities.
    rs = @reaction_network begin
        (k1,k2), X1 <--> X2
        (k1,k2), Y1 <--> Y2
    end
    osys = convert(ODESystem, rs; remove_conserved = true, remove_conserved_warn)

    # Checks that the correct parameters have the `conserved` metadata.
    @test Catalyst.isconserved(osys.Γ[1])
    @test Catalyst.isconserved(osys.Γ[2])
    @test !Catalyst.isconserved(osys.k1)
    @test !Catalyst.isconserved(osys.k2)
end

# Checks that conservation law elimination warnings are generated in the correct cases.
let
    # Prepare model.
    rn = @reaction_network begin
        (k1,k2), X1 <--> X2
    end
    u0 = [:X1 => 1.0, :X2 => 2.0]
    tspan = (0.0, 1.0)
    ps = [:k1 => 3.0, :k2 => 4.0]

    # Check warnings in system conversion.
    for XSystem in [ODESystem, SDESystem, NonlinearSystem]
        @test_nowarn convert(XSystem, rn)
        @test_logs (:warn, r"You are creating a system or problem while eliminating conserved quantities. Please *") convert(XSystem, rn; remove_conserved = true)
        @test_nowarn convert(XSystem, rn; remove_conserved_warn = false)
        @test_nowarn convert(XSystem, rn; remove_conserved = true, remove_conserved_warn = false)
    end

    # Checks during problem creation (separate depending on whether they have a time span or not).
    @test_broken false # SDEProblems generate anacounted for warnings here. These should be figured out and `warn_initialize_determined = false` removed (https://github.com/SciML/Catalyst.jl/pull/1159).
    for XProblem in [ODEProblem, SDEProblem]
        @test_nowarn XProblem(rn, u0, tspan, ps, warn_initialize_determined = false)
        @test_logs (:warn, r"You are creating a system or problem while eliminating conserved quantities. Please *") XProblem(rn, u0, tspan, ps; remove_conserved = true, warn_initialize_determined = false)
        @test_nowarn XProblem(rn, u0, tspan, ps; remove_conserved_warn = false, warn_initialize_determined = false)
        @test_nowarn XProblem(rn, u0, tspan, ps; remove_conserved = true, remove_conserved_warn = false, warn_initialize_determined = false)
    end
    for XProblem in [NonlinearProblem, SteadyStateProblem]
        @test_nowarn XProblem(rn, u0, ps)
        @test_logs (:warn, r"You are creating a system or problem while eliminating conserved quantities. Please *") XProblem(rn, u0, ps; remove_conserved = true)
        @test_nowarn XProblem(rn, u0, ps; remove_conserved_warn = false)
        @test_nowarn XProblem(rn, u0, ps; remove_conserved = true, remove_conserved_warn = false)
    end
end

# Conservation law simulations for vectorised species.
let
    # Prepares the model.
    t = default_t()
    @species (X(t))[1:2]
    @parameters k[1:2]
    rxs = [
        Reaction(k[1], [X[1]], [X[2]]),
        Reaction(k[2], [X[2]], [X[1]])
    ]
    @named rs = ReactionSystem(rxs, t)
    rs = complete(rs)

    # Checks that simulation reaches known equilibrium.
    @test_broken false # Currently broken on MTK .
    # u0 = [:X => [3.0, 9.0]]
    # ps = [:k => [1.0, 2.0]]
    # oprob = ODEProblem(rs, u0, (0.0, 1000.0), ps; remove_conserved = true)
    # sol = solve(oprob, Vern7())
    # @test sol[X[1]][end] ≈ 8.0
    # @test sol[X[2]][end] ≈ 4.0
end
