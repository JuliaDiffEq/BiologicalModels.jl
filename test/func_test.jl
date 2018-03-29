using DiffEqBiological, StochasticDiffEq

@reaction_func new_hill(x, v, k, n) = v*x^n/(k^n+x^n)
@reaction_func new_poly(x) = 3x^2+1
@reaction_func new_exp(x) = exp(x)

network1 = @reaction_network rn begin
    hill(X,4,3,2), X + Y --> Z1
    3X^2+1, X + Y --> Z2
    exp(Y), X + Y --> Z3
end
network2 = @reaction_network rn begin
    new_hill(X,4,3,2), X + Y --> Z1
    new_poly(X), X + Y --> Z2
    new_exp(Y), X + Y --> Z3
end

for i = 1:100
    u = 5*rand(5)
    du1 = 3*rand(5); du2 = du1;
    du1g = 2.5*rand(5,3); du2g = du1g;
    t = 9*rand(1)[1]
    p = []

    @test network1.f(du1,u,p,t) == network2.f(du2,u,p,t)
    @test network1.g(du1g,u,p,t) == network2.g(du2g,u,p,t)
    @test network1.jumps[1].rate(u,p,t) == network2.jumps[1].rate(u,p,t)
    @test network1.jumps[2].rate(u,p,t) == network2.jumps[2].rate(u,p,t)
    @test network1.jumps[3].rate(u,p,t) == network2.jumps[3].rate(u,p,t)
end

model1 = @reaction_network rn begin
    (5,5000), X ↔ 0
end
model2 = @reaction_network rn η begin
    (5,5000), X ↔ 0
end η
model3 = @reaction_network rn η begin
    (1,5000), X ↔ 0
end η
function tmp_std(sol)
    vect = Vector{Float64}(length(sol.u))
    for i = 1:length(sol.u)
        vect[i] = sol.u[i][1]
    end
    return std(vect)
end
prob1 = SDEProblem(model1,[1000.],(0.,200.))
prob2 = SDEProblem(model2,[1000.],(0.,200.), 0.1)
prob3 = SDEProblem(model3,[1000.],(0.,200.),10)
sol1 = solve(prob1, EM(), dt = 0.0001)
sol2 = solve(prob2, EM(), dt = 0.0001)
sol3 = solve(prob3, EM(), dt = 0.0001)
@test tmp_std(sol2) < tmp_std(sol1) < tmp_std(sol3)
