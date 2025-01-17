# [Surface Plasmon Resonance Model Example](@id surface_plasmon_resonance_model_Example)
This tutorial shows how to programmatically construct a ['ReactionSystem'](@ref) corresponding to the chemistry underlying the [Surface Plasmon Resonance Model](https://en.wikipedia.org/wiki/Surface_plasmon_resonance) using [ModelingToolkit](http://docs.sciml.ai/ModelingToolkit/stable/)/[Catalyst](http://docs.sciml.ai/Catalyst/stable/). You can find a simpler constructruction of this model [here](https://docs.sciml.ai/Catalyst/stable/catalyst_applications/parameter_estimation/) in the example of parameter estimation. 

Using the symbolic interface, we will create our states/species, define our reactions, and add a discrete event. This model is simulated and used to generate sample data points with added noise and then used to fit parameters in iterations. This event will correspond to the time at which the antibody stop binding to antigen and switch to dissociating. 

We begin by importing some necessary packages. 
```julia
using ModelingToolkit, Catalyst
using Optimization, OptimizatizationOptimJL
using Plots
```
In this example, the concentration of antigen,`\beta` is varied to determine the constant of proportionality, `\alpha`, `k_{on}`(association rate constant), and `k_{off}` (dissociation rate constant) which characterized the binding interaction between the antigen and the immobilized receptor molecules on the sensor slide. We start by defining our reaction equations, parameters, variables, event, and states. 
```julia
osys = @reaction_network begin
    @variables t 
    @parameters k_on=100.0 k_off α
    @species A(t)B(t)
    
    @discrete_events begin
        t == switch_time => [k_on ~ 0.0] 
    end

    α*k_on, A --> B
    k_off, B --> A
end

switch_time = 2.0

tspan = (0.0, 4.0)
alpha_list = [0.1, 0.2, 0.3, 0.4] #list of concentrations 
```
Iterating over values of `\alpha`, now we create a list of ODE solutions and set the initial conditions of our states and parameters. 
```julia 
results_list = []

u0 = [:A => 10.0, :B => 0.0]
p_real = [k_on => 100.0, k_off => 10.0, α => 1.0]
oprob = ODEProblem(osys, u0, tspan, p_real)
sample_times = range(tspan[1]; stop = tspan[2], length = 1001) 

for alpha in alpha_list
    p_real = [k_on => 100.0, k_off => 10.0, α => alpha]
    oprobr = remake(oprob, p=p_real)
    sol_real = solve(oprobr, Tsit5(); tstops = sample_times)

    push!(results_list, sol_real(sample_times))
end
```
```@example ceq3
default(; lw = 3, framestyle = :box, size = (800, 400))
p = plot()
plot(p, sample_times, results_list[1][2,:])
plot!(p, sample_times, results_list[2][2,:])
plot!(p, sample_times, results_list[3][2,:])
plot!(p, sample_times, results_list[4][2,:])

sample_vals = []
for result in results_list
    sample_val = Array(result)
    sample_val .*= (1 .+ .1 * rand(Float64, size(sample_val)) .- .01)
    push!(sample_vals, sample_val)
end

for val in sample_vals
    scatter!(p, sample_times, val[2,:]; color = [:blue :red], legend = nothing)
end
plot(p)
```
Next, we create a function to fit the parameters. Here we are using `NelderMead()`. In order to achieve a better fit, we are incorporating all of our solutions into the loss function. We will fit seperately the association and dissociation signals so for the first estimate, `tend < switch_time`. 
```julia 
function optimise_p(pinit, tend)
    function loss(p, _)
        newtimes = filter(<=(tend), sample_times)
        solutions = []
        for alpha in alpha_list
            newprob = remake(oprob; tspan = (0.0, tend), p = [k_on => p[1],k_off => p[2],α => alpha])
            sol = Array(solve(newprob, Tsit5(); saveat = newtimes, tstops = switch_time))
            push!(solutions,sol[2,:])
        end
        loss = 0
        for (idx, solution) in enumerate(solutions)
            loss += sum(abs2, p[3]*solution .- sample_vals[idx][2, 1:size(sol,2)])
        end

        return loss
    end

    optf = OptimizationFunction(loss)
    
    optprob = OptimizationProblem(optf, pinit)
    sol = solve(optprob, Optim.NelderMead())

    return sol.u
end

p_estimate = optimise_p([100.0, 10.0, 1.0], 1.5)
```
Finally, we remake the solution using the estimate, fit the entire time span for some `alpha`, and plot the solution. 
```julia
alpha = 0.2
newprob = remake(oprob; tspan = (0.0, 4.0), p = [k_on => p_estimate[1], k_off => p_estimate[2], α => alpha])
newsol = solve(newprob, Tsit5(); tstops = switch_time)

#plotting simulated sample date with new solution
default(; lw = 3, framestyle = :box, size = (800, 400))

scatter!(sample_times, sample_vals; color = [:darkblue :darkred], legend = nothing)

plot!(newsol[2,:]; legend = nothing, color = [:blue :red], linestyle = :dash, tspan= (0.0, 4.0))
```
