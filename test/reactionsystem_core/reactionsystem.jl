### Fetch Packages and Set Global Variables ###

# Fetch packages.
using Catalyst, LinearAlgebra, JumpProcesses, OrdinaryDiffEqTsit5, OrdinaryDiffEqVerner, StochasticDiffEq, Test
const MT = ModelingToolkit

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)

# Sets the default `t` to use.
t = default_t()

# Fetch test functions.
include("../test_functions.jl")

### Creates Basic Test Network ###

# Create the network.
@parameters k[1:20]
@species A(t) B(t) C(t) D(t)
rxs = [Reaction(k[1], nothing, [A]),            # 0 -> A
    Reaction(k[2], [B], nothing),            # B -> 0
    Reaction(k[3], [A], [C]),                  # A -> C
    Reaction(k[4], [C], [A, B]),              # C -> A + B
    Reaction(k[5], [C], [A], [1], [2]),      # C -> A + A
    Reaction(k[6], [A, B], [C]),              # A + B -> C
    Reaction(k[7], [B], [A], [2], [1]),      # 2B -> A
    Reaction(k[8], [A, B], [A, C]),            # A + B -> A + C
    Reaction(k[9], [A, B], [C, D]),            # A + B -> C + D
    Reaction(k[10], [A], [C, D], [2], [1, 1]), # 2A -> C + D
    Reaction(k[11], [A], [A, B], [2], [1, 1]), # 2A -> A + B
    Reaction(k[12], [A, B, C], [C, D], [1, 3, 4], [2, 3]),          # A+3B+4C -> 2C + 3D
    Reaction(k[13], [A, B], nothing, [3, 1], nothing),           # 3A+B -> 0
    Reaction(k[14], nothing, [A], nothing, [2]),               # 0 -> 2A
    Reaction(k[15] * A / (2 + A), [A], nothing; only_use_rate = true), # A -> 0 with custom rate
    Reaction(k[16], [A], [B]; only_use_rate = true),             # A -> B with custom rate.
    Reaction(k[17] * A * exp(B), [C], [D], [2], [1]),              # 2C -> D with non constant rate.
    Reaction(k[18] * B, nothing, [B], nothing, [2]),             # 0 -> 2B with non constant rate.
    Reaction(k[19] * t, [A], [B]),                                # A -> B with non constant rate.
    Reaction(k[20] * t * A, [B, C], [D], [2, 1], [2]),                  # 2A +B -> 2C with non constant rate.
]
@named rs = ReactionSystem(rxs, t, [A, B, C, D], [k])
rs = complete(rs)
odesys = complete(convert(ODESystem, rs))
sdesys = complete(convert(SDESystem, rs))

# Hard coded ODE rhs.
function oderhs(u, kv, t)
    A = u[1]
    B = u[2]
    C = u[3]
    D = u[4]
    k = kv[1]
    du = zeros(eltype(u), 4)
    du[1] = k[1] - k[3] * A + k[4] * C + 2 * k[5] * C - k[6] * A * B + k[7] * B^2 / 2 -
            k[9] * A * B - k[10] * A^2 - k[11] * A^2 / 2 - k[12] * A * B^3 * C^4 / 144 -
            3 * k[13] * A^3 * B / 6 + 2 * k[14] - k[15] * A / (2 + A) - k[16] -
            k[19] * t * A
    du[2] = -k[2] * B + k[4] * C - k[6] * A * B - k[7] * B^2 - k[8] * A * B - k[9] * A * B +
            k[11] * A^2 / 2 - 3 * k[12] * A * B^3 * C^4 / 144 - k[13] * A^3 * B / 6 +
            k[16] + 2 * k[18] * B + k[19] * t * A - 2 * k[20] * t * A * B^2 * C
    du[3] = k[3] * A - k[4] * C - k[5] * C + k[6] * A * B + k[8] * A * B + k[9] * A * B +
            k[10] * A^2 / 2 - 2 * k[12] * A * B^3 * C^4 / 144 -
            2 * k[17] * A * exp(B) * C^2 / 2 - k[20] * t * A * B^2 * C
    du[4] = k[9] * A * B + k[10] * A^2 / 2 + 3 * k[12] * A * B^3 * C^4 / 144 +
            k[17] * A * exp(B) * C^2 / 2 + 2 * k[20] * t * A * B^2 * C
    du
end

# SDE noise coefs.
function sdenoise(u, kv, t)
    A = u[1]
    B = u[2]
    C = u[3]
    D = u[4]
    k = kv[1]
    G = zeros(eltype(u), length(k), length(u))
    z = zero(eltype(u))

    G = [sqrt(k[1]) z z z;
         z -sqrt(k[2] * B) z z;
         -sqrt(k[3] * A) z sqrt(k[3] * A) z;
         sqrt(k[4] * C) sqrt(k[4] * C) -sqrt(k[4] * C) z;
         2*sqrt(k[5] * C) z -sqrt(k[5] * C) z;
         -sqrt(k[6] * A * B) -sqrt(k[6] * A * B) sqrt(k[6] * A * B) z;
         sqrt(k[7] * B^2 / 2) -2*sqrt(k[7] * B^2 / 2) z z;
         z -sqrt(k[8] * A * B) sqrt(k[8] * A * B) z;
         -sqrt(k[9] * A * B) -sqrt(k[9] * A * B) sqrt(k[9] * A * B) sqrt(k[9] * A * B);
         -2*sqrt(k[10] * A^2 / 2) z sqrt(k[10] * A^2 / 2) sqrt(k[10] * A^2 / 2);
         -sqrt(k[11] * A^2 / 2) sqrt(k[11] * A^2 / 2) z z;
         -sqrt(k[12] * A * B^3 * C^4 / 144) -3*sqrt(k[12] * A * B^3 * C^4 / 144) -2*sqrt(k[12] * A * B^3 * C^4 / 144) 3*sqrt(k[12] * A * B^3 * C^4 / 144);
         -3*sqrt(k[13] * A^3 * B / 6) -sqrt(k[13] * A^3 * B / 6) z z;
         2*sqrt(k[14]) z z z;
         -sqrt(k[15] * A / (2 + A)) z z z;
         -sqrt(k[16]) sqrt(k[16]) z z;
         z z -2*sqrt(k[17] * A * exp(B) * C^2 / 2) sqrt(k[17] * A * exp(B) * C^2 / 2);
         z 2*sqrt(k[18] * B) z z;
         -sqrt(k[19] * t * A) sqrt(k[19] * t * A) z z;
         z -2*sqrt(k[20] * t * A * B^2 * C) -sqrt(k[20] * t * A * B^2 * C) +2*sqrt(k[20] * t * A * B^2 * C)]'
    return G
end

### Basic Tests ###

# Test equation only constructor.
let
    @named rs2 = ReactionSystem(rxs, t)
    @test Catalyst.isequivalent(rs, rs2)

    # Test with a type mismatch
    @test Catalyst.isequivalent(rs, "Not a ReactionSystem") == false
end

# Defaults test.
let
    kvals = Float64.(1:length(k))
    def_p = [k => kvals]
    def_u0 = [A => 0.5, B => 1.0, C => 1.5, D => 2.0]
    defs = merge(Dict(def_p), Dict(def_u0))

    @named rs = ReactionSystem(rxs, t, [A, B, C, D], [k]; defaults = defs)
    rs = complete(rs)
    odesys = complete(convert(ODESystem, rs))
    sdesys = complete(convert(SDESystem, rs))
    js = complete(convert(JumpSystem, rs))

    @test ModelingToolkit.get_defaults(rs) ==
          ModelingToolkit.get_defaults(odesys) ==
          ModelingToolkit.get_defaults(sdesys) ==
          ModelingToolkit.get_defaults(js) ==
          defs

    u0map = [A => 5.0]
    kvals[1] = 5.0
    pmap = [k => kvals]
    prob = ODEProblem(rs, u0map, (0, 10.0), pmap)
    @test prob.ps[k[1]] == 5.0
    @test prob.u0[1] == 5.0
end

### Check ODE, SDE, and Jump Functions ###

# Test by evaluating drift and diffusion terms.

let
    u = rnd_u0(rs, rng)
    p = rnd_ps(rs, rng)
    du = oderhs(last.(u), last.(p), 0.0)
    G = sdenoise(last.(u), last.(p), 0.0)
    sdesys = complete(convert(SDESystem, rs))
    sf = SDEFunction{false}(sdesys, unknowns(rs), parameters(rs))
    sprob = SDEProblem(rs, u, (0.0, 0.0), p)
    du2 = sf.f(sprob.u0, sprob.p, 0.0)

    du2 = sf.f(sprob.u0, sprob.p, 0.0)
    @test norm(du - du2) < 100 * eps()
    G2 = sf.g(sprob.u0, sprob.p, 0.0)
    @test norm(G - G2) < 100 * eps()
end

# Test with JumpSystem.
let
    @species A(t) B(t) C(t) D(t) E(t) F(t)
    rxs = [Reaction(k[1], nothing, [A]),            # 0 -> A
        Reaction(k[2], [B], nothing),            # B -> 0
        Reaction(k[3], [A], [C]),                  # A -> C
        Reaction(k[4], [C], [A, B]),              # C -> A + B
        Reaction(k[5], [C], [A], [1], [2]),      # C -> A + A
        Reaction(k[6], [A, B], [C]),              # A + B -> C
        Reaction(k[7], [B], [A], [2], [1]),      # 2B -> A
        Reaction(k[8], [A, B], [A, C]),            # A + B -> A + C
        Reaction(k[9], [A, B], [C, D]),            # A + B -> C + D
        Reaction(k[10], [A], [C, D], [2], [1, 1]), # 2A -> C + D
        Reaction(k[11], [A], [A, B], [2], [1, 1]), # 2A -> A + B
        Reaction(k[12], [A, B, C], [C, D], [1, 3, 4], [2, 3]),          # A+3B+4C -> 2C + 3D
        Reaction(k[13], [A, B], nothing, [3, 1], nothing),           # 3A+B -> 0
        Reaction(k[14], nothing, [A], nothing, [2]),               # 0 -> 2A
        Reaction(k[15] * A / (2 + A), [A], nothing; only_use_rate = true), # A -> 0 with custom rate
        Reaction(k[16], [A], [B]; only_use_rate = true),             # A -> B with custom rate.
        Reaction(k[17] * A * exp(B), [C], [D], [2], [1]),              # 2C -> D with non constant rate.
        Reaction(k[18] * B, nothing, [B], nothing, [2]),             # 0 -> 2B with non constant rate.
        Reaction(k[19] * t, [D], [E]),                                # D -> E with non constant rate.
        Reaction(k[20] * t * A, [D, E], [F], [2, 1], [2]),                  # 2D + E -> 2F with non constant rate.
    ]
    @named rs = ReactionSystem(rxs, t, [A, B, C, D, E, F], [k])
    rs = complete(rs)
    js = complete(convert(JumpSystem, rs))

    midxs = 1:14
    cidxs = 15:18
    vidxs = 19:20
    @test all(map(i -> typeof(equations(js)[i]) <: JumpProcesses.MassActionJump, midxs))
    @test all(map(i -> typeof(equations(js)[i]) <: JumpProcesses.ConstantRateJump, cidxs))
    @test all(map(i -> typeof(equations(js)[i]) <: JumpProcesses.VariableRateJump, vidxs))

    p = rand(rng, length(k))
    pmap = [k => p]
    u0 = rand(rng, 2:10, 6)
    u0map = unknowns(js) .=> u0
    ttt = rand(rng)
    jumps = Vector{Union{ConstantRateJump, MassActionJump, VariableRateJump}}(undef,
                                                                              length(rxs))

    jumps[1] = MassActionJump(p[1], Vector{Pair{Int, Int}}(), [1 => 1])
    jumps[2] = MassActionJump(p[2], [2 => 1], [2 => -1])
    jumps[3] = MassActionJump(p[3], [1 => 1], [1 => -1, 3 => 1])
    jumps[4] = MassActionJump(p[4], [3 => 1], [1 => 1, 2 => 1, 3 => -1])
    jumps[5] = MassActionJump(p[5], [3 => 1], [1 => 2, 3 => -1])
    jumps[6] = MassActionJump(p[6], [1 => 1, 2 => 1], [1 => -1, 2 => -1, 3 => 1])
    jumps[7] = MassActionJump(p[7], [2 => 2], [1 => 1, 2 => -2])
    jumps[8] = MassActionJump(p[8], [1 => 1, 2 => 1], [2 => -1, 3 => 1])
    jumps[9] = MassActionJump(p[9], [1 => 1, 2 => 1], [1 => -1, 2 => -1, 3 => 1, 4 => 1])
    jumps[10] = MassActionJump(p[10], [1 => 2], [1 => -2, 3 => 1, 4 => 1])
    jumps[11] = MassActionJump(p[11], [1 => 2], [1 => -1, 2 => 1])
    jumps[12] = MassActionJump(p[12], [1 => 1, 2 => 3, 3 => 4],
                               [1 => -1, 2 => -3, 3 => -2, 4 => 3])
    jumps[13] = MassActionJump(p[13], [1 => 3, 2 => 1], [1 => -3, 2 => -1])
    jumps[14] = MassActionJump(p[14], Vector{Pair{Int, Int}}(), [1 => 2])

    jumps[15] = ConstantRateJump((u, p, t) -> p[15] * u[1] / (2 + u[1]),
                                 integrator -> (integrator.u[1] -= 1))
    jumps[16] = ConstantRateJump((u, p, t) -> p[16],
                                 integrator -> (integrator.u[1] -= 1; integrator.u[2] += 1))
    jumps[17] = ConstantRateJump((u, p, t) -> p[17] * u[1] * exp(u[2]) * binomial(u[3], 2),
                                 integrator -> (integrator.u[3] -= 2; integrator.u[4] += 1))
    jumps[18] = ConstantRateJump((u, p, t) -> p[18] * u[2],
                                 integrator -> (integrator.u[2] += 2))

    jumps[19] = VariableRateJump((u, p, t) -> p[19] * u[4] * t,
                                 integrator -> (integrator.u[4] -= 1; integrator.u[5] += 1))
    jumps[20] = VariableRateJump((u, p, t) -> p[20] * t * u[1] * binomial(u[4], 2) * u[5],
                                 integrator -> (integrator.u[4] -= 2; integrator.u[5] -= 1; integrator.u[6] += 2))

    unknownoid = Dict(unknown => i for (i, unknown) in enumerate(unknowns(js)))
    dprob = DiscreteProblem(js, u0map, (0.0, 10.0), pmap)
    mtkpars = dprob.p
    jspmapper = ModelingToolkit.JumpSysMajParamMapper(js, mtkpars)
    symmaj = ModelingToolkit.assemble_maj(equations(js).x[1], unknownoid, jspmapper)
    maj = MassActionJump(symmaj.param_mapper(mtkpars), symmaj.reactant_stoch, symmaj.net_stoch,
                         symmaj.param_mapper, scale_rates = false)
    for i in midxs
        @test abs(jumps[i].scaled_rates - maj.scaled_rates[i]) < 100 * eps()
        @test jumps[i].reactant_stoch == maj.reactant_stoch[i]
        @test jumps[i].net_stoch == maj.net_stoch[i]
    end
    for i in cidxs
        crj = ModelingToolkit.assemble_crj(js, equations(js)[i], unknownoid)
        @test isapprox(crj.rate(u0, mtkpars, ttt), jumps[i].rate(u0, p, ttt))
        fake_integrator1 = (u = zeros(6), p = mtkpars, t = 0.0)
        fake_integrator2 = (u = zeros(6), p, t = 0.0)
        crj.affect!(fake_integrator1)
        jumps[i].affect!(fake_integrator2)
        @test fake_integrator1.u == fake_integrator2.u
    end
    for i in vidxs
        crj = ModelingToolkit.assemble_vrj(js, equations(js)[i], unknownoid)
        @test isapprox(crj.rate(u0, mtkpars, ttt), jumps[i].rate(u0, p, ttt))
        fake_integrator1 = (u = zeros(6), p = mtkpars, t = 0.0)
        fake_integrator2 = (u = zeros(6), p, t = 0.0)
        crj.affect!(fake_integrator1)
        jumps[i].affect!(fake_integrator2)
        @test fake_integrator1.u == fake_integrator2.u
    end
end

### Niche Model Declarations ###

# Checks model with vector species and parameters.
# Checks that it works for programmatic/dsl-based modelling.
# Checks that all forms of model input (parameter/initial condition and vector/non-vector) are
# handled properly.
let
    # Declares programmatic model.
    @parameters p[1:2] k d1 d2
    @species (X(t))[1:2] Y1(t) Y2(t)
    rxs = [
        Reaction(p[1], [], [X[1]]),
        Reaction(p[2], [], [X[2]]),
        Reaction(k, [X[1]], [Y1]),
        Reaction(k, [X[2]], [Y2]),
        Reaction(d1, [Y1], []),
        Reaction(d2, [Y2], []),
    ]
    rs_prog = complete(ReactionSystem(rxs, t; name = :rs))

    # Declares DSL-based model.
    rs_dsl = @reaction_network rs begin
        @parameters p[1:2] k d1 d2
        @species (X(t))[1:2] Y1(t) Y2(t)
        (p[1],p[2]), 0 --> (X[1],X[2])
        k, (X[1],X[2]) --> (Y1,Y2)
        (d1,d2), (Y1,Y2) --> 0
    end

    # Checks equivalence.
    rs_dsl == rs_prog

    # Creates all possible initial conditions and parameter values.
    u0_alts = [
        [X => [2.0, 5.0], Y1 => 0.2, Y2 => 0.5],
        [X[1] => 2.0, X[2] => 5.0, Y1 => 0.2, Y2 => 0.5],
        [rs_dsl.X => [2.0, 5.0], rs_dsl.Y1 => 0.2, rs_dsl.Y2 => 0.5],
        [rs_dsl.X[1] => 2.0, X[2] => 5.0, rs_dsl.Y1 => 0.2, rs_dsl.Y2 => 0.5],
        [:X => [2.0, 5.0], :Y1 => 0.2, :Y2 => 0.5]
    ]
    ps_alts = [
        [p => [1.0, 10.0], d1 => 5.0, d2 => 4.0, k => 2.0],
        [p[1] => 1.0, p[2] => 10.0, d1 => 5.0, d2 => 4.0, k => 2.0],
        [rs_dsl.p => [1.0, 10.0], rs_dsl.d1 => 5.0, rs_dsl.d2 => 4.0, rs_dsl.k => 2.0],
        [rs_dsl.p[1] => 1.0, p[2] => 10.0, rs_dsl.d1 => 5.0, rs_dsl.d2 => 4.0, rs_dsl.k => 2.0],
        [:p => [1.0, 10.0], :d1 => 5.0, :d2 => 4.0, :k => 2.0]
    ]

    # Loops through all inputs and check that the correct steady state is reached
    # Target steady state: (X1, X2, Y1, Y2) = (p1/k, p2/k, p1/d1, p2/d2).
    # Technically only one model needs to be check. However, "equivalent" models in MTK can still
    # have slight differences, so checking for both here to be certain.
    for rs in [rs_prog, rs_dsl]
        oprob = ODEProblem(rs, u0_alts[1], (0.0, 10000.), ps_alts[1])
        @test_broken false # Cannot currently `remake` this problem/
        # for rs in [rs_prog, rs_dsl], u0 in u0_alts, p in ps_alts
        #     oprob_remade = remake(oprob; u0, p)
        #     sol = solve(oprob_remade, Vern7(); abstol = 1e-8, reltol = 1e-8)
        #     @test sol[end] ≈ [0.5, 5.0, 0.2, 2.5]
        # end
    end
end

### Miscellaneous Getters ###

# Tests spatial independent variables `get_sivs` and `has_sivs`.
let
    # Creates a `ReactionSystem`.
    @parameters x
    @parameters p d
    @species S(t,x)
    rxs = [
        Reaction(p, [], [S]),
        Reaction(d, [S], []),
    ]
    @named rs = ReactionSystem(rxs, t; spatial_ivs = [x])

    # Test system independent variables (spatial and non-spatial).
    @test Catalyst.has_sivs(rs)
    @test isequal(Catalyst.get_sivs(rs), [x])
    @test isequal(Catalyst.independent_variables(rs), [t])
end

# Tests `numparams` and `numreactions` for reaction system with subsystem.
# Tests where a subssystem is a non-`ReactionSystem`.
let
    # Prepares content
    @parameters k1 k2 k3 p1 p2
    @species X1(t) X2(t) X3(t)
    @variables V1(t) V2(t)
    D = default_time_deriv()

    # Creates a reaction system with a subsystem.
    sub_rxs = [
        Reaction(k1, [X1], []),
        Reaction(k2, [X2], [])
    ]
    @named sub_rs = ReactionSystem(sub_rxs, t)
    sub_eqs = [
        D(V1) ~ p1 - V1,
        D(V2) ~ p2 - V2,
    ]
    @named sub_osys = ODESystem(sub_eqs, t)
    rxs = [
        Reaction(k2, [X2], []),
        Reaction(k3, [X3], [])
    ]
    @named rs = ReactionSystem(rxs, t; systems = [sub_rs, sub_osys])

    # Tests content.
    @test numparams(rs) == 6
    @test numreactions(rs) == 4
end

# Tests `has_sivs` getter.

### Test Show ###

# Basic show test.
let
    io = IOBuffer()
    show(io, rs)
    str = String(take!(io))
    @test count(isequal('\n'), str) < 30
end

# Test printing with arrays is working ok.
# Needs fix for https://github.com/JuliaSymbolics/Symbolics.jl/issues/842.
let
    @parameters a
    @species A(t) B(t) (C(t))[1:2]
    rx1 = Reaction(a, [A, C[1]], [C[2], B], [1, 2], [2, 3])
    io = IOBuffer()
    show(io, rx1)
    str = String(take!(io))
    @test str == "a, A + 2*(C(t))[1] --> 2*(C(t))[2] + 3*B"
end

### Boundary Condition Species Tests ###

# Test for constant and boundary condition species.
function f!(du, u, p, t)
    A = p[1]
    k1 = p[2]
    k2 = p[3]
    B = u[1]
    D = u[2]
    E = u[3]
    C = u[4]
    du[1] = k1 * A - k2 * B
    du[2] = -k1 * C * D + k2 * C * E
    du[3] = k1 * C * D - k2 * C * E
    du[4] = -C
    nothing
end
function fs!(du, u, p, t)
    A = p[1]
    k1 = p[2]
    k2 = p[3]
    B = u[1]
    D = u[2]
    E = u[3]
    C = u[4]
    du[1] = k1 * A - k2 * B
    du[2] = -k1 * C * D + k2 * C * E
    du[3] = k1 * C * D - k2 * C * E
    nothing
end
function gs!(dg, u, p, t)
    A = p[1]
    k1 = p[2]
    k2 = p[3]
    B = u[1]
    D = u[2]
    E = u[3]
    C = u[4]
    dg .= 0.0
    dg[1, 1] = sqrt(k1 * A)
    dg[1, 2] = -sqrt(k2 * B)
    dg[2, 3] = -sqrt(k1 * C * D)
    dg[2, 4] = sqrt(k2 * C * E)
    dg[3, 3] = -dg[2, 3]
    dg[3, 4] = -dg[2, 4]
    nothing
end

# Tests for BC and constant species.
let
    @parameters k1 k2 A [isconstantspecies = true]
    @species B(t) C(t) [isbcspecies = true] D(t) E(t)
    Dt = default_time_deriv()
    eqs = [(@reaction k1, $A --> B),
        (@reaction k2, B --> $A),
        (@reaction k1, $C + D --> E + $C),
        Dt(C) ~ -C,
        (@reaction k2, E + $C --> $C + D)]
    @named rs = ReactionSystem(eqs, t)
    rs = complete(rs)
    @test all(eq -> eq isa Reaction, ModelingToolkit.get_eqs(rs)[1:4])
    osys = complete(convert(ODESystem, rs))
    @test issetequal(MT.get_unknowns(osys), [B, C, D, E])
    @test issetequal(MT.get_ps(osys), [k1, k2, A])

    # test nonlinear systems
    u0 = [1.0, 2.0, 3.0, 4.0]
    p = [2.0, 2.5, 3.5]
    u0map = [B, D, E, C] .=> u0
    pmap = [A, k1, k2] .=> p
    tspan = (0.0, 5.0)
    oprob1 = ODEProblem(osys, u0map, tspan, pmap)
    sts = [B, D, E, C]
    syms = [:B, :D, :E, :C]
    ofun = ODEFunction(f!; sys = ModelingToolkit.SymbolCache(syms))
    oprob2 = ODEProblem(ofun, u0, tspan, p)
    saveat = tspan[2] / 50
    abstol = 1e-10
    reltol = 1e-10
    sol1 = solve(oprob1, Tsit5(); saveat, abstol, reltol)
    sol2 = solve(oprob2, Tsit5(); saveat, abstol, reltol)
    for i in eachindex(sts)
        @test isapprox(sol1[sts[i]], sol2[syms[i]])
    end

    # Test sde systems.
    rxs = [(@reaction k1, $A --> B),
        (@reaction k2, B --> $A),
        (@reaction k1, $C + D --> E + $C),
        (@reaction k2, E + $C --> $C + D)]
    @named rs = ReactionSystem(rxs, t)   # add constraint csys when supported!
    rs = complete(rs)
    ssys = complete(convert(SDESystem, rs))
    @test issetequal(MT.get_unknowns(ssys), [B, C, D, E])
    @test issetequal(MT.get_ps(ssys), [A, k1, k2])
    du1 = zeros(4)
    du2 = zeros(4)
    sprob = SDEProblem(ssys, u0map, tspan, pmap; check_length = false)
    sprob.f(du1, sprob.u0, sprob.p, 1.0)
    fs!(du2, u0, p, 1.0)
    @test isapprox(du1, du2)
    dg1 = zeros(4, 4)
    dg2 = zeros(4, 4)
    sprob.g(dg1, sprob.u0, sprob.p, 1.0)
    gs!(dg2, u0, p, t)
    @test isapprox(dg1, dg2)

    # Test jump systems.
    rxs = [(@reaction k1, $A --> B),
        (@reaction k2, B --> $A),
        (@reaction k1, $C + D --> E + $C),
        (@reaction k2, $C + E --> $C + D),
        (@reaction k1 * t, $A + $C --> B + $C),
        (@reaction k1 * B, 2 * $A + $C --> $C + B)]
    @named rs = ReactionSystem(rxs, t)
    rs = complete(rs)
    jsys = complete(convert(JumpSystem, rs))
    @test issetequal(unknowns(jsys), [B, C, D, E])
    @test issetequal(parameters(jsys), [k1, k2, A])
    majrates = [k1 * A, k1, k2]
    majrs = [[], [C => 1, D => 1], [C => 1, E => 1]]
    majns = [[B => 1], [D => -1, E => 1], [D => 1, E => -1]]
    for (i, maj) in enumerate(equations(jsys).x[1])
        @test isequal(maj.scaled_rates, majrates[i])
        @test issetequal(maj.reactant_stoch, majrs[i])
        @test issetequal(maj.net_stoch, majns[i])
    end
    @test isempty(equations(jsys).x[2])
    vrj1 = equations(jsys).x[3][1]
    @test isequal(vrj1.rate, k2 * B)
    @test issetequal(vrj1.affect!, [B ~ B - 1])
    vrj2 = equations(jsys).x[3][2]
    @test isequal(vrj2.rate, k1 * t * A * C)
    @test issetequal(vrj2.affect!, [B ~ B + 1])
    vrj3 = equations(jsys).x[3][3]
    @test isequal(vrj3.rate, k1 * B * A * (A - 1) / 2 * C)
    @test issetequal(vrj3.affect!, [B ~ B + 1])
end

# Test that jump solutions actually run correctly for constants and BCs.
let
    @parameters k1 A [isconstantspecies = true]
    @species C(t) [isbcspecies = true] B1(t) B2(t) B3(t)
    @named rn = ReactionSystem([(@reaction k1, $C --> B1 + $C),
                                   (@reaction k1, $A --> B2),
                                   (@reaction 10 * k1, ∅ --> B3)], t)
    rn = complete(rn)
    dprob = DiscreteProblem(rn, [A => 10, C => 10, B1 => 0, B2 => 0, B3 => 0], (0.0, 10.0),
                            [k1 => 1.0])
    jprob = JumpProblem(rn, dprob, Direct(); rng, save_positions = (false, false))
    umean = zeros(4)
    Nsims = 40000
    for i in 1:Nsims
        sol = solve(jprob, SSAStepper(), saveat = 10.0)
        umean += sol(10.0, idxs = [B1, B2, B3, C])
    end
    umean /= Nsims
    @test isapprox(umean[1], umean[2]; rtol = 1e-2)
    @test isapprox(umean[1], umean[3]; rtol = 1e-2)
    @test umean[4] == 10
end

### Species Tests ###

# Test various species related checker functions.
let
    # Creates species and parameters.
    @species X(t) Y(t) [isbcspecies=true]
    @parameters x(t) y(t) [isconstantspecies=true]

    # Tests properties.
    @test !isspecies(x)
    @test !isspecies(y)
    @test isspecies(X)
    @test isspecies(Y)
    @test !Catalyst.isbc(x)
    @test !Catalyst.isbc(y)
    @test !Catalyst.isbc(X)
    @test Catalyst.isbc(Y)
    @test !Catalyst.isconstant(x)
    @test Catalyst.isconstant(y)
    @test !Catalyst.isconstant(X)
    @test !Catalyst.isconstant(Y)
end

### Error Tests ###

# Tests various erroneous `ReactionSystem` creations.
let
    # Prepare model inputs.
    @parameters k1 k2 x [isconstantspecies=true] Γ
    @species X1(t) X2(t)
    @variables V(t)

    # System using a forbidden symbol.
    @test_throws Exception rs = ReactionSystem([Reaction(Γ, [X1], [])], t; name = :rs)

    # Species among parameters.
    @test_throws Exception rs = ReactionSystem([Reaction(k1, [X1], [])], t, [X1], [k1, X2])

    # Variable among parameters.
    @test_throws Exception rs = ReactionSystem([Reaction(k1, [X1], [])], t, [X1], [k1, V])

    # Parameter among unknowns.
    @test_throws Exception rs = ReactionSystem([Reaction(k1, [X1], [])], t, [X1, k2], [k1])

    # Constant species-parameter among unknowns.
    @test_throws Exception rs = ReactionSystem([Reaction(k1, [X1], [])], t, [X1, x], [k1])
end

# Tests various erroneous `convert` calls.
let
    # Conversion of non-autonomous `ReactionSystem` to `NonlinearSystem`.
    rs = @reaction_network begin
        (p/(1+t),d), 0 <--> X
    end
    @test_throws Exception convert(NonlinearSystem, rs)

    # Conversion of non-complete system to various system types.
    nc = @network_component begin
        (p,d), 0 <--> X
    end
    @test_throws Exception convert(ODESystem, nc)
    @test_throws Exception convert(SDESystem, nc)
    @test_throws Exception convert(JumpSystem, nc)
    @test_throws Exception convert(NonlinearSystem, nc)
end

# Checks that the same name cannot be used for two different of parameters/species/variables.
let
    # Stores a parameter, a species, and a variable (with identical names) in different variables.
    x_p = let
        only(@parameters x)
    end
    x_sp = let
        only(@species x(t))
    end
    x_v = let
        only(@variables x(t))
    end

    # Checks that creating systems with different in combination produces errors.
    # Currently broken on MTK, potentially fix in Catalyst once sorted out there (https://github.com/SciML/ModelingToolkit.jl/issues/2883).
    @parameters d
    @species X(t)
    rx = Reaction(d, [X], [])
    @test_broken false # (not sure how to mark a `@test_throws` as broken)
    # @test_throws rs1 = ReactionSystem([rx], t, [X, x_sp,], [d, x_p]; name = :rs)
    # @test_throws rs2 = ReactionSystem([rx], t, [X, X, x_v], [d, x_p]; name = :rs)
    # @test_throws rs2 = ReactionSystem([rx], t, [X, x_sp, x_v], [d]; name = :rs)
end

### Other Tests ###

# Test for https://github.com/SciML/ModelingToolkit.jl/issues/436.
let
    @parameters t
    @species S(t) I(t)
    rxs = [Reaction(1, [S], [I]), Reaction(1.1, [S], [I])]
    @named rs = ReactionSystem(rxs, t, [S, I], [])
    rs = complete(rs)
    js = complete(convert(JumpSystem, rs))
    dprob = DiscreteProblem(js, [S => 1, I => 1], (0.0, 10.0))
    jprob = JumpProblem(js, dprob, Direct(); rng)
    sol = solve(jprob, SSAStepper())

    # Test for https://github.com/SciML/ModelingToolkit.jl/issues/1042.
    jprob = JumpProblem(rs, dprob, Direct(); rng, save_positions = (false, false))

    @parameters k1 k2
    @species R(t)
    rxs = [Reaction(k1 * S, [S, I], [I], [2, 3], [2]),
        Reaction(k2 * R, [I], [R])]
    @named rs = ReactionSystem(rxs, t, [S, I, R], [k1, k2])
    rs = complete(rs)
    @test isequal(oderatelaw(equations(rs)[1]),
                  k1 * S * S^2 * I^3 / (factorial(2) * factorial(3)))
    @test_skip isequal(jumpratelaw(equations(eqs)[1]),
                       k1 * S * binomial(S, 2) * binomial(I, 3))
    dep = Set()
    ModelingToolkit.get_variables!(dep, rxs[2], Set(unknowns(rs)))
    dep2 = Set([R, I])
    @test dep == dep2
    dep = Set()
    ModelingToolkit.modified_unknowns!(dep, rxs[2], Set(unknowns(rs)))
    @test dep == Set([R, I])

    isequal2(a, b) = isequal(simplify(a), simplify(b))

    @test isequal2(jumpratelaw(rxs[1]), k1 * S * S * (S - 1) * I * (I - 1) * (I - 2) / 12)
    @test isequal2(jumpratelaw(rxs[1]; combinatoric_ratelaw = false),
                   k1 * S * S * (S - 1) * I * (I - 1) * (I - 2))
    @test isequal2(oderatelaw(rxs[1]), k1 * S * S^2 * I^3 / 12)
    @test isequal2(oderatelaw(rxs[1]; combinatoric_ratelaw = false), k1 * S * S^2 * I^3)

    @named rs2 = ReactionSystem(rxs, t, [S, I, R], [k1, k2]; combinatoric_ratelaws = false)
    rs2 = complete(rs2)

    # Test ODE scaling:
    os = complete(convert(ODESystem, rs))
    @test isequal2(equations(os)[1].rhs, -2 * k1 * S * S^2 * I^3 / 12)
    os = convert(ODESystem, rs; combinatoric_ratelaws = false)
    @test isequal2(equations(os)[1].rhs, -2 * k1 * S * S^2 * I^3)
    os2 = complete(convert(ODESystem, rs2))
    @test isequal2(equations(os2)[1].rhs, -2 * k1 * S * S^2 * I^3)
    os3 = complete(convert(ODESystem, rs2; combinatoric_ratelaws = true))
    @test isequal2(equations(os3)[1].rhs, -2 * k1 * S * S^2 * I^3 / 12)

    # Test ConstantRateJump rate scaling.
    js = complete(convert(JumpSystem, rs))
    @test isequal2(equations(js)[1].rate,
                   k1 * S * S * (S - 1) * I * (I - 1) * (I - 2) / 12)
    js = complete(convert(JumpSystem, rs; combinatoric_ratelaws = false))
    @test isequal2(equations(js)[1].rate, k1 * S * S * (S - 1) * I * (I - 1) * (I - 2))
    js2 = complete(convert(JumpSystem, rs2))
    @test isequal2(equations(js2)[1].rate, k1 * S * S * (S - 1) * I * (I - 1) * (I - 2))
    js3 = complete(convert(JumpSystem, rs2; combinatoric_ratelaws = true))
    @test isequal2(equations(js3)[1].rate,
                   k1 * S * S * (S - 1) * I * (I - 1) * (I - 2) / 12)

    # Test MassActionJump rate scaling.
    rxs = [Reaction(k1, [S, I], [I], [2, 3], [2]),
        Reaction(k2, [I], [R])]
    @named rs = ReactionSystem(rxs, t, [S, I, R], [k1, k2])
    rs = complete(rs)
    js = complete(convert(JumpSystem, rs))
    @test isequal2(equations(js)[1].scaled_rates, k1 / 12)
    js = complete(convert(JumpSystem, rs; combinatoric_ratelaws = false))
    @test isequal2(equations(js)[1].scaled_rates, k1)

    # test building directly from rxs
    @parameters x, y
    rxs = [Reaction(x * t * A * B + y, [A], nothing)]
    @named rs1 = ReactionSystem(rxs, t, [A, B], [x, y])
    @named rs2 = ReactionSystem(rxs, t)
    @test Catalyst.isequivalent(rs1, rs2)

    @species L(t), H(t)
    obs = [Equation(L, 2 * x + y)]
    @named rs3 = ReactionSystem(rxs, t; observed = obs)
    L2 = L
    @unpack L = rs3
    @test isequal(L, L2)
end

# Test that non-integer stoichiometry goes through.
let
    @parameters k b
    @species A(t) B(t) C(t) D(t)
    rx1 = Reaction(k, [B, C], [B, D], [2.5, 1], [3.5, 2.5])
    rx2 = Reaction(2 * k, [B], [D], [1], [2.5])
    rx3 = Reaction(2 * k, [B], [D], [2.5], [2])
    @named mixedsys = ReactionSystem([rx1, rx2, rx3], t, [A, B, C, D], [k, b])
    mixedsys = complete(mixedsys)
    osys = convert(ODESystem, mixedsys; combinatoric_ratelaws = false)
end

# Test balanced_bc_check.
let
    @species A(t) [isbcspecies = true]
    rx = @reaction k, 2 * $A + B --> C + $A
    @test_throws ErrorException ReactionSystem([rx], t; name = :rs)
    @named rs = ReactionSystem([rx], t; balanced_bc_check = false)
end

# Fix for SBML test 305.
let
    @parameters k1 k2 S2 [isconstantspecies = true]
    @species S1(t) S3(t)
    rx = Reaction(k2, [S1], nothing)
    ∂ₜ = default_time_deriv()
    eq = ∂ₜ(S3) ~ k1 * S2
    @mtkbuild osys = ODESystem([eq], t)
    @named rs = ReactionSystem([rx, eq], t)
    rs = complete(rs)
    @test issetequal(unknowns(rs), [S1, S3])
    @test issetequal(parameters(rs), [S2, k1, k2])
    osys = convert(ODESystem, rs)
    @test issetequal(unknowns(osys), [S1, S3])
    @test issetequal(parameters(osys), [S2, k1, k2])
end
let
    @parameters k1 k2 S2 [isconstantspecies = true]
    @species S1(t) S3(t) [isbcspecies = true]
    rx = Reaction(k2, [S1], nothing)
    ∂ₜ = default_time_deriv()
    eq = S3 ~ k1 * S2
    @named rs = ReactionSystem([rx, eq], t)
    rs = complete(rs)
    @test issetequal(unknowns(rs), [S1, S3])
    @test issetequal(parameters(rs), [S2, k1, k2])
    osys = convert(ODESystem, rs)
    @test issetequal(unknowns(osys), [S1, S3])
    @test issetequal(parameters(osys), [S2, k1, k2])
    osys2 = structural_simplify(osys)
    @test length(equations(osys2)) == 1
    @test issetequal(unknowns(osys2), [S1])
    @test issetequal(parameters(osys2), [S2, k1, k2])
end
let
    @parameters k1 k2 S2 [isconstantspecies = true]
    @variables S3(t)
    @species S1(t)
    rx = Reaction(k2, [S1], nothing)
    ∂ₜ = default_time_deriv()
    eq = S3 ~ k1 * S2
    @named rs = ReactionSystem([rx, eq], t)
    rs = complete(rs)
    @test issetequal(unknowns(rs), [S1, S3])
    @test issetequal(parameters(rs), [S2, k1, k2])
    osys = convert(ODESystem, rs)
    @test issetequal(unknowns(osys), [S1, S3])
    @test issetequal(parameters(osys), [S2, k1, k2])
    osys2 = structural_simplify(osys)
    @test length(equations(osys2)) == 1
    @test issetequal(unknowns(osys2), [S1])
    @test issetequal(parameters(osys2), [S2, k1, k2])
end

# Constant species = parameters basic tests.
let
    @parameters k b [isconstantspecies = true] c
    @test_throws ArgumentError @species A(t) B(t) a [isconstantspecies = true]
    @test_throws ArgumentError Reaction(k, [A, c], [B])
    @test_throws ArgumentError Reaction(k, [A], [B, c])
    rx = Reaction(k, [A, b], [B, b], [1, 1], [1, 2])
    @named rs = ReactionSystem([rx], t)
    @test issetequal(unknowns(rs), [A, B])
    @test issetequal(parameters(rs), [k, b])
end

# Test parameteric initial conditions.
let
    @parameters d X0
    @species X(t)=X0
    rx = Reaction(d, [X], nothing, [1], nothing)
    @named rs = ReactionSystem([rx], t)
    rs = complete(rs)
    prob = ODEProblem(rs, [], (0.0, 1.0), [d => 1.0, X0 => 7.6])
    @test prob[X] == 7.6
end

# Test for classification of jump types.
let
    rn = @reaction_network begin
        t, A --> B          # vrj
        1.0, B --> D        # vrj
        k * D, H --> I + H  # vrj
        k2, I --> L         # vrj
        k, E --> F          # maj
        k * E, E --> G      # crj
        k2, G --> H         # maj
        k2, G --> A + B     # maj
    end
    jsys = convert(JumpSystem, rn)
    jumps = Catalyst.assemble_jumps(rn)
    @test count(j -> j isa VariableRateJump, jumps) == 4
    @test count(j -> j isa ConstantRateJump, jumps) == 1
    @test count(j -> j isa MassActionJump, jumps) == 3
    dg = [[1, 2], [2, 3], [4], [4], [5, 6], [5, 6, 7, 8], [3, 7, 8], [1, 2, 7, 8]]
    dgact = Catalyst.get_depgraph(rn)
    @test dg == dgact
end

# Test array metadata for species works.
let
    @species (A(t))[1:20]
    using ModelingToolkit: value
    Av = value(A)
    @test isspecies(Av)
    @test all(i -> isspecies(Av[i]), 1:length(Av))
end

# Test mixed models are formulated correctly.
let
    @parameters k1 k2
    @variables V(t)
    @species A(t) B(t)
    rx = Reaction(k1, [A], [B], [k2], [2])
    D = default_time_deriv()
    eq = D(V) ~ -k1 * k2 * V + A
    @named rs = ReactionSystem([eq, rx], t)
    rs = complete(rs)
    @test length(unknowns(rs)) == 3
    @test issetequal(unknowns(rs), [A, B, V])
    @test length(parameters(rs)) == 2
    @test issetequal(parameters(rs), [k1, k2])
    @test length(species(rs)) == 2
    @test issetequal(species(rs), [A, B])
    @test all(typeof.(ModelingToolkit.get_eqs(rs)) .<: (Reaction, Equation))
    @test length(Catalyst.get_rxs(rs)) == 1
    @test reactions(rs)[1] == rx
    osys = convert(ODESystem, rs)
    @test issetequal(unknowns(osys), [A, B, V])
    @test issetequal(parameters(osys), [k1, k2])
    @test length(equations(osys)) == 3
end

# Test errors for repeated substrates or products
let
    @species A(t) B(t)
    @test_throws ArgumentError Reaction(1.0, [A, A, B], [B])
    @test_throws ArgumentError Reaction(1.0, [B], [A, A])
    @test_throws ArgumentError Reaction(1.0, [A, A], [B, B])
end

# Test order of species and products doesn't matter for equality or hashing
let
    @species A(t) α(t)
    rx = Reaction(1.0, [α, A], [α, A], [2, 3], [4, 5])
    rx2 = Reaction(1.0, [A, α], [A, α], [3, 2], [5, 4])
    @test rx == rx2
    @test hash(rx) == hash(rx2)

    rx = Reaction(1.0, [α, A], [α, A], [2, 3], [4, 5]; netstoich = [α => 2, A => 2])
    rx2 = Reaction(1.0, [A, α], [A, α], [3, 2], [5, 4]; netstoich = [A => 2, α => 2])
    @test rx == rx2
    @test hash(rx) == hash(rx2)
end

# Additional unsorted tests.
let
    rn = @reaction_network begin k, X --> 0 end
    isspecies(species(rn)[1])
    @test Catalyst.has_species(rn)
    @test Catalyst.has_rxs(rn)

    @species X
    @variables Y
    @test isspecies(X)
    @test !isspecies(Y)
    @test isspecies(Catalyst.tospecies(Y))
end

# Tests system metadata.
let
    @test isnothing(ModelingToolkit.get_metadata(rs))
end

# Tests construction of empty reaction networks.
let
    # Using DSL.
    empty_network = @reaction_network
    @test length(ModelingToolkit.get_eqs(empty_network)) == 0
    @test nameof(ModelingToolkit.get_iv(empty_network)) == :t
    @test length(ModelingToolkit.get_unknowns(empty_network)) == 0
    @test length(ModelingToolkit.get_ps(empty_network)) == 0

    # Using `make_empty_network`.
    empty_network = make_empty_network()
    @test length(ModelingToolkit.get_eqs(empty_network)) == 0
    @test nameof(ModelingToolkit.get_iv(empty_network)) == :t
    @test length(ModelingToolkit.get_unknowns(empty_network)) == 0
    @test length(ModelingToolkit.get_ps(empty_network)) == 0
end

# Checks that the `reactionsystem_uptodate` function work. If it does not, the ReactionSystem
# structure's fields have been updated, without updating the `reactionsystem_fields` constant. If so,
# there are several places in the code where the `reactionsystem_uptodate` function is called, here
# the code might need adaptation to take the updated reaction system into account.
let
    @test_nowarn Catalyst.reactionsystem_uptodate_check()
end

# Test that functions using the incidence matrix properly cache it
let
    rn = @reaction_network begin
        k1, A --> B
        k2, B --> C
        k3, C --> A
    end

    nps = Catalyst.get_networkproperties(rn)
    @test isempty(nps.incidencemat) == true

    img = incidencematgraph(rn)
    @test size(nps.incidencemat) == (3,3)

    Catalyst.reset!(nps)
    lcs = linkageclasses(rn)
    @test size(nps.incidencemat) == (3,3)

    Catalyst.reset!(nps)
    sns = subnetworks(rn)
    @test size(nps.incidencemat) == (3,3)

    Catalyst.reset!(nps)
    δ = deficiency(rn)
    @test size(nps.incidencemat) == (3,3)

    Catalyst.reset!(nps)
    δ_l = linkagedeficiencies(rn)
    @test size(nps.incidencemat) == (3,3)

    Catalyst.reset!(nps)
    rev = isreversible(rn)
    @test size(nps.incidencemat) == (3,3)

    Catalyst.reset!(nps)
    weakrev = isweaklyreversible(rn, sns)
    @test size(nps.incidencemat) == (3,3)
end

########## tests related to hybrid systems ##########
let
    t = default_t()
    D = default_time_deriv()
    @parameters λ k
    @variables V(t)
    @species A(t) B(t) C(t)
    rxs = [Reaction(k*V, [], [A]), Reaction(λ*A, [B], nothing),
        Reaction(k, [A, B], nothing), Reaction(λ, [C], [A])]
    eqs = [D(V) ~ λ*V*C]
    cevents = [[V ~ 2.0] => [V ~ V/2, A ~ A/2]]
    @named rs = ReactionSystem(vcat(rxs, eqs), t; continuous_events = cevents)
    rs = complete(rs)
    jinput = JumpInputs(rn, [:A => 0, :V => 1.0], (0.0, 10.0), [:k => 1.0, :λ => .4])
    @test jinput.prob isa ODEProblem
    sys = jinput.sys
    @test sys isa JumpSystem
    @test MT.has_equations(sys)
    eqs = MT.equations(sys)
    @test isempty(eqs.x[1]) 
    @test isempty(eqs.x[2]) 
    @test length(eqs.x[3]) == 1
    @test length(eqs.x[4]) == 2  # dA/dt and dV/dt
    @test length(continuous_events(sys)) == 1
end
