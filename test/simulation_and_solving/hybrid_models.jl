# Fetch packages.
using Catalyst, JumpProcesses, OrdinaryDiffEqTsit5, Statistics, Test

# Sets stable rng number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rng, 1:100)

# test JumpInputs function auto problem selection
let
    rn = @reaction_network begin
        k*(1 + sin(t)), 0 --> A
    end
    jinput = JumpInputs(rn, [:A => 0], (0.0, 10.0), [:k => .5])
    @test jinput.prob isa ODEProblem
    jprob = JumpProblem(jinput; rng)
    sol = solve(jprob, Tsit5())
    @test sol(10.0; idxs = :A) > 0

    rn = @reaction_network begin
        k, 0 --> A
    end
    jinput = JumpInputs(rn, [:A => 0], (0.0, 10.0), [:k => .5])
    @test jinput.prob isa DiscreteProblem
    jprob = JumpProblem(jinput; rng)
    sol = solve(jprob)
    @test sol(10.0; idxs = :A) > 0

    rn = @reaction_network begin
        @parameters λ
        k*V, 0 --> A
        @equations D(V) ~ λ*V
        @continuous_events begin
            [V ~ 2.0] => [V ~ V/2, A ~ A/2]
        end
    end        
    jinput = JumpInputs(rn, [:A => 0, :V => 1.0], (0.0, 10.0), [:k => 1.0, :λ => .4])
    @test jinput.prob isa ODEProblem
    jprob = JumpProblem(jinput; rng)
    sol = solve(jprob, Tsit5()) 
end

    