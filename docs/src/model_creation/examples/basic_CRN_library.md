# [Library of Basic Chemical Reaction Network Models](@id basic_CRN_library)
```@raw html
<details><summary><strong>Environment setup and package installation</strong></summary>
```
The following code sets up an environment for running the code on this page.
```julia
using Pkg
Pkg.activate(".")
Pkg.add("Catalyst")
Pkg.add("OrdinaryDiffEqDefault")
Pkg.add("StochasticDiffEq")
Pkg.add("JumpProcesses")
Pkg.add("Plots")
```
```@raw html
</details>
```
  \
  
Below we will present various simple and established chemical reaction network (CRN) models. Each model is given some brief background, implemented using the `@reaction_network` DSL, and basic simulations are performed.

## [Birth-death process](@id basic_CRN_library_bd)
The birth-death process is one of the simplest possible CRN models. It consists of a single component ($X$) which is both produced and degraded at linear rates:
```@example crn_library_birth_death
using Catalyst
bd_process = @reaction_network begin
    (p,d), ∅ <--> X
end
```
Next, we define simulation conditions. Note that the initial condition is integer-valued (more natural than decimal numbers for jump simulations).
```@example crn_library_birth_death
u0 = [:X => 1]
tspan = (0.0, 10.0)
ps = [:p => 1.0, :d => 0.2]
nothing # hide
```
We can now simulate our model using all three interpretations. First, we perform a reaction rate equation-based ODE simulation:
```@example crn_library_birth_death
using OrdinaryDiffEqDefault
oprob = ODEProblem(bd_process, u0, tspan, ps)
osol = solve(oprob)
nothing # hide
```
Next, a chemical Langevin equation-based SDE simulation:
```@example crn_library_birth_death
using StochasticDiffEq
sprob = SDEProblem(bd_process, u0, tspan, ps)
ssol = solve(sprob, STrapezoid())
ssol = solve(sprob, STrapezoid(); seed = 12) # hide
nothing # hide
```
Next, a stochastic chemical kinetics-based jump simulation:
```@example crn_library_birth_death
using JumpProcesses
jinput = JumpInputs(bd_process, u0, tspan, ps)
jprob = JumpProblem(jinput)
jsol = solve(jprob)
jsol = solve(jprob; seed = 12) # hide
nothing # hide
```
Finally, we plot the results:
```@example crn_library_birth_death
using Plots
oplt = plot(osol; title = "Reaction rate equation (ODE)")
splt = plot(ssol; title = "Chemical Langevin equation (SDE)")
jplt = plot(jsol; title = "Stochastic chemical kinetics (Jump)")
plot(oplt, splt, jplt; lw = 3, size=(800,700), layout = (3,1))
```

## [Two-state model](@id basic_CRN_library_two_states)
The two-state model describes a component (here called $X$) which can exist in two different forms (here called $X₁$ and $X₂$). It switches between these forms at linear rates. First, we simulate the model using both ODEs and SDEs:
```@example crn_library_two_states
using Catalyst, OrdinaryDiffEqDefault, StochasticDiffEq, Plots
two_state_model = @reaction_network begin
    (k₁,k₂), X₁ <--> X₂
end

u0 = [:X₁ => 50.0, :X₂ => 50.0]
tspan = (0.0, 1.0)
ps = [:k₁ => 2.0, :k₂ => 3.0]

oprob = ODEProblem(two_state_model, u0, tspan, ps)
osol = solve(oprob)
oplt = plot(osol; title = "Reaction rate equation (ODE)")

sprob = SDEProblem(two_state_model, u0, tspan, ps)
ssol = solve(sprob, STrapezoid())
ssol = solve(sprob, STrapezoid(); seed = 12) # hide
splt = plot(ssol; title = "Chemical Langevin equation (SDE)")

plot(oplt, splt; lw = 3, size = (800,550), layout = (2,1))
```
What is interesting about this model is that it has a *conserved quantity*, where $X₁ + X₂$ remains constant throughout the simulation (both in deterministic and stochastic cases). We can show this by instead plotting this conserved quantity.
```@example crn_library_two_states
@unpack X₁, X₂ = two_state_model
oplt = plot(osol; idxs = X₁ + X₂, title = "Reaction rate equation (ODE)")
splt = plot(ssol; idxs = X₁ + X₂, title = "Chemical Langevin equation (SDE)")
plot(oplt, splt; lw = 3, ylimit = (99,101), size = (800,450), layout = (2,1))
```

## [Michaelis-Menten enzyme kinetics](@id basic_CRN_library_mm)
[Michaelis-Menten enzyme kinetics](https://en.wikipedia.org/wiki/Michaelis%E2%80%93Menten_kinetics) is a simple description of an enzyme ($E$) transforming a substrate ($S$) into a product ($P$). Under certain assumptions, it can be simplified to a single function (a Michaelis-Menten function) and used as a reaction rate. Here we instead present the full system model:
```@example crn_library_michaelis_menten
using Catalyst
mm_system = @reaction_network begin
    kB, S + E --> SE
    kD, SE --> S + E
    kP, SE --> P + E
end
```
Next, we perform ODE, SDE, and jump simulations of the model:
```@example crn_library_michaelis_menten
u0 = [:S => 301, :E => 100, :SE => 0, :P => 0]
tspan = (0., 100.)
ps = [:kB => 0.00166, :kD => 0.0001, :kP => 0.1]

using OrdinaryDiffEqDefault
oprob = ODEProblem(mm_system, u0, tspan, ps)
osol  = solve(oprob)

using StochasticDiffEq
sprob = SDEProblem(mm_system, u0, tspan, ps)
ssol = solve(sprob, STrapezoid())
ssol = solve(sprob, STrapezoid(); seed = 12) # hide

using JumpProcesses
jinput = JumpInputs(mm_system, u0, tspan, ps)
jprob = JumpProblem(jinput)
jsol = solve(jprob)
jsol = solve(jprob; seed = 12) # hide

using Plots
oplt = plot(osol; title = "Reaction rate equation (ODE)")
splt = plot(ssol; title = "Chemical Langevin equation (SDE)")
jplt = plot(jsol; title = "Stochastic chemical kinetics (Jump)")
plot(oplt, splt, jplt; lw = 2, size=(800,800), layout = (3,1))
fullplt = plot(oplt, splt, jplt; lw = 2, size=(800,800), layout = (3,1), fmt = :png,
    dpi = 200, bottom_margin = 3Plots.Measures.mm) # hide
Catalyst.PNG(fullplt) # hide
```

## [SIR infection model](@id basic_CRN_library_sir)
The [SIR model](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model) is the simplest model of the spread of an infectious disease. While the real system is very different from the chemical and cellular processes typically modelled with CRNs, it (and several other epidemiological systems) can be modelled using the same CRN formalism. The SIR model consists of three species: susceptible ($S$), infected ($I$), and removed ($R$) individuals, and two reaction events: infection and recovery.
```@example crn_library_sir
using Catalyst
sir_model = @reaction_network begin
    α, S + I --> 2I
    β, I --> R
end
```
First, we perform a deterministic ODE simulation:
```@example crn_library_sir
using OrdinaryDiffEqDefault, Plots
u0 = [:S => 99, :I => 1, :R => 0]
tspan = (0.0, 500.0)
ps = [:α => 0.001, :β => 0.01]

# Solve ODEs.
oprob = ODEProblem(sir_model, u0, tspan, ps)
osol = solve(oprob)
plot(osol; title = "Reaction rate equation (ODE)", size=(800,350))
```
Next, we perform 3 different Jump simulations. Note that for the stochastic model, the occurrence of an outbreak is not certain. Rather, there is a possibility that it fizzles out without a noteworthy peak.
```@example crn_library_sir
using JumpProcesses
jinput = JumpInputs(sir_model, u0, tspan, ps)
jprob = JumpProblem(jinput)

jsol1 = solve(jprob)
jsol2 = solve(jprob)
jsol3 = solve(jprob)
jsol1 = solve(jprob, seed = 1) # hide
jsol2 = solve(jprob, seed = 2) # hide
jsol3 = solve(jprob, seed = 3) # hide

jplt1 = plot(jsol1; title = "Outbreak")
jplt2 = plot(jsol2; title = "Outbreak")
jplt3 = plot(jsol3; title = "No outbreak")
plot(jplt1, jplt2, jplt3; lw = 3, size=(800,700), layout = (3,1))
fullplt = plot(jplt1, jplt2, jplt3; lw = 3, size=(800,700), layout = (3,1), fmt = :png, dpi = 200) # hide
Catalyst.PNG(fullplt) # hide
```

## [Chemical cross-coupling](@id basic_CRN_library_cc)
In chemistry, [cross-coupling](https://en.wikipedia.org/wiki/Cross-coupling_reaction) is when a catalyst combines two substrates to form a product. In this example, the catalyst ($C$) first binds one substrate ($S₁$) to form an intermediary complex ($S₁C$). Next, the complex binds the second substrate ($S₂$) to form another complex ($CP$). Finally, the catalyst releases the now-formed product ($P$). This system is an extended version of the [Michaelis-Menten system presented earlier](@ref basic_CRN_library_mm).
```@example crn_library_cc
using Catalyst
cc_system = @reaction_network begin
    k₁, S₁ + C --> S₁C
    k₂, S₁C + S₂ --> CP
    k₃, CP --> C + P
end
```
Below, we perform a simple deterministic ODE simulation of the system. Next, we plot both:
- The concentration of the substrates and the product.
- The concentration of the catalyst and the intermediaries.

In two separate plots.
```@example crn_library_cc
using OrdinaryDiffEqDefault, Plots
u0 = [:S₁ => 1.0, :C => 0.05, :S₂ => 1.2, :S₁C => 0.0, :CP => 0.0, :P => 0.0]
tspan = (0., 15.)
ps = [:k₁ => 5.0, :k₂ => 5.0, :k₃ => 100.0]

# solve ODEs
oprob = ODEProblem(cc_system, u0, tspan, ps)
osol  = solve(oprob)

plt1 = plot(osol; idxs = [:S₁, :S₂, :P], title = "Substrate and product dynamics")
plt2 = plot(osol; idxs = [:C, :S₁C, :CP], title = "Catalyst and intermediaries dynamics")
plot(plt1, plt2; lw = 3, size = (800,600), layout = (2,1))
```

## [The Wilhelm model](@id basic_CRN_library_wilhelm)
The Wilhelm model was introduced in [*Wilhelm (2009)*](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-3-90) as the smallest CRN model (with constant rates) that exhibits bistability.
```@example crn_library_wilhelm
using Catalyst
wilhelm_model = @reaction_network begin
    k1, Y --> 2X
    k2, 2X --> X + Y
    k3, X + Y --> Y
    k4, X --> 0
end
```
We can simulate the model for two different initial conditions, demonstrating the existence of two different stable steady states.
```@example crn_library_wilhelm
using OrdinaryDiffEqDefault, Plots
u0_1 = [:X => 1.5, :Y => 0.5]
u0_2 = [:X => 2.5, :Y => 0.5]
tspan = (0., 10.)
ps = [:k1 => 8.0, :k2 => 2.0, :k3 => 1.0, :k4 => 1.5]

oprob1 = ODEProblem(wilhelm_model, u0_1, tspan, ps)
oprob2 = ODEProblem(wilhelm_model, u0_2, tspan, ps)
osol1 = solve(oprob1)
osol2 = solve(oprob2)
plot(osol1; lw = 4, idxs = :X, label = "X(0) = 1.5")
plot!(osol2; lw = 4, idxs = :X, label = "X(0) = 2.5", yguide = "X", size = (800,350))
plot!(bottom_margin = 3Plots.Measures.mm) # hide
```

## [Simple self-activation loop](@id basic_CRN_library_self_activation)
The simplest self-activation loop consists of a single species (here called $X$) which activates its own production. If its production rate is modelled as a [Hill function](https://en.wikipedia.org/wiki/Hill_equation_(biochemistry)) with $n>1$, the system may exhibit bistability.
```@example crn_library_self_activation
using Catalyst
sa_loop = @reaction_network begin
    v₀ + hill(X,v,K,n), ∅ --> X
    d, X --> ∅
end
```
A simple example of such a loop is a transcription factor which activates its own gene. Here, $v₀$ represents a basic transcription rate (leakage) in the absence of the transcription factor.

We simulate the self-activation loop from a single initial condition using both deterministic (ODE) and stochastic (jump) simulations. We note that while the deterministic simulation reaches a single steady state, the stochastic one switches between two different states.
```@example crn_library_self_activation
using JumpProcesses, OrdinaryDiffEqDefault, Plots
u0 = [:X => 4]
tspan = (0.0, 1000.0)
ps = [:v₀ => 0.1, :v => 2.0, :K => 10.0, :n => 2, :d => 0.1]

oprob = ODEProblem(sa_loop, u0, tspan, ps)
osol = solve(oprob)

jinput = JumpInputs(sa_loop, u0, tspan, ps)
jprob = JumpProblem(jinput)
jsol = solve(jprob)
jsol = solve(jprob, seed = 12) # hide

fplt = plot(osol; lw = 3, label = "Reaction rate equation (ODE)")
plot!(fplt, jsol; lw = 3, label = "Stochastic chemical kinetics (Jump)", yguide = "X", size = (800,350))
plot!(fplt, bottom_margin = 3Plots.Measures.mm, left_margin=3Plots.Measures.mm) # hide
plot!(fplt; fmt = :png, dpi = 200) # hide
Catalyst.PNG(fplt) # hide
```

## [The Brusselator](@id basic_CRN_library_brusselator)
The [Brusselator](https://en.wikipedia.org/wiki/Brusselator) is a well-known (theoretical) CRN model able to produce oscillations (its name is a portmanteau of "Brussels" and "oscillator").
```@example crn_library_brusselator
using Catalyst
brusselator = @reaction_network begin
    A, ∅ --> X
    1, 2X + Y --> 3X
    B, X --> Y
    1, X --> ∅
end
```
It is generally known to (for reaction rate equation-based ODE simulations) produce oscillations when $B > 1 + A^2$. However, this result is based on models generated when *combinatorial adjustment of rates is not performed*. Since Catalyst [automatically perform these adjustments](@ref introduction_to_catalyst_ratelaws), and one reaction contains a stoichiometric constant $>1$, the threshold will be different. Here, we trial two different values of $B$. In both cases, $B < 1 + A^2$, however, in the second case the system can generate oscillations.
```@example crn_library_brusselator
using OrdinaryDiffEqDefault, Plots
u0 = [:X => 1.0, :Y => 1.0]
tspan = (0., 50.)
ps1 = [:A => 1.0, :B => 1.0]
ps2 = [:A => 1.0, :B => 1.8]

oprob1 = ODEProblem(brusselator, u0, tspan, ps1)
oprob2 = ODEProblem(brusselator, u0, tspan, ps2)
osol1  = solve(oprob1)
osol2  = solve(oprob2)
oplt1 = plot(osol1; title = "No Oscillation")
oplt2 = plot(osol2; title = "Oscillation")

plot(oplt1, oplt2; lw = 3, size = (800,600), layout = (2,1))
```

## [The Repressilator](@id basic_CRN_library_repressilator)
The Repressilator was introduced in [*Elowitz & Leibler (2000)*](https://www.nature.com/articles/35002125) as a simple system that can generate oscillations (most notably, they demonstrated this both in a model and in a synthetic in vivo implementation in *Escherichia col*). It consists of three genes, repressing each other in a cycle. Here, we will implement it using three species ($X$, $Y$, and $Z$) whose production rates are (repressing) [Hill functions](https://en.wikipedia.org/wiki/Hill_equation_(biochemistry)).
```@example crn_library_brusselator
using Catalyst
repressilator = @reaction_network begin
    hillr(Z,v,K,n), ∅ --> X
    hillr(X,v,K,n), ∅ --> Y
    hillr(Y,v,K,n), ∅ --> Z
    d, (X, Y, Z) --> ∅
end
```
Whether the Repressilator oscillates or not depends on its parameter values. Here, we will perform deterministic (ODE) simulations for two different values of $K$, showing that it oscillates for one value and not the other one. Next, we will perform stochastic (SDE) simulations for both $K$ values, showing that the stochastic model can sustain oscillations in both cases. This is an example of the phenomena of *noise-induced oscillation*.
```@example crn_library_brusselator
using OrdinaryDiffEqDefault, StochasticDiffEq, Plots
u0 = [:X => 50.0, :Y => 15.0, :Z => 15.0]
tspan = (0., 200.)
ps1 = [:v => 10.0, :K => 20.0, :n => 3, :d => 0.1]
ps2 = [:v => 10.0, :K => 50.0, :n => 3, :d => 0.1]

oprob1 = ODEProblem(repressilator, u0, tspan, ps1)
oprob2 = ODEProblem(repressilator, u0, tspan, ps2)
osol1  = solve(oprob1)
osol2  = solve(oprob2)
oplt1 = plot(osol1; title = "Oscillation (ODE, K = 20)")
oplt2 = plot(osol2; title = "No oscillation (ODE, K = 50)")

sprob1 = SDEProblem(repressilator, u0, tspan, ps1)
sprob2 = SDEProblem(repressilator, u0, tspan, ps2)
ssol1  = solve(sprob1, STrapezoid())
ssol2  = solve(sprob2, STrapezoid())
ssol1  = solve(sprob1, STrapezoid(); seed = 1) # hide
ssol2  = solve(sprob2, STrapezoid(); seed = 100) # hide
splt1 = plot(ssol1; title = "Oscillation (SDE, K = 20)")
splt2 = plot(ssol2; title = "Oscillation (SDE, K = 50)")

plot(oplt1, oplt2, splt1, splt2; lw = 2, layout = (2,2), size = (800,600))
```

## [The Willamowski–Rössler model](@id basic_CRN_library_wr)
The Willamowski–Rössler model was introduced in [*Willamowski & Rössler (1979)*](https://www.degruyter.com/document/doi/10.1515/zna-1980-0308/html?lang=en) as an example of a simple CRN model which exhibits [*chaotic behaviours*](https://en.wikipedia.org/wiki/Chaos_theory). This means that small changes in initial conditions can produce relatively large changes in the system's trajectory.
```@example crn_library_chaos
using Catalyst
wr_model = @reaction_network begin
    k1, 2X --> 3X
    k2, X --> 2X
    k3, Z + 2X --> 2Z
    k4, Y + X --> 2Y
    k5, Y --> ∅
    k6, 2Z --> ∅
    k7, Z --> ∅
end
```
Here we simulate the model for a single initial condition, showing both time-state space and phase space how it reaches a [*strange attractor*](https://www.dynamicmath.xyz/strange-attractors/).
```@example crn_library_chaos
using OrdinaryDiffEqDefault, Plots
u0 = [:X => 1.5, :Y => 1.5, :Z => 1.5]
tspan = (0.0, 50.0)
p = [:k1 => 2.1, :k2 => 0.7, :k3 => 2.9, :k4 => 1.1, :k5 => 1.0, :k6 => 0.5, :k7 => 2.7]
oprob = ODEProblem(wr_model, u0, tspan, p)
sol = solve(oprob)

plt1 = plot(sol; title = "Time-state space")
plt2 = plot(sol; idxs = (:X, :Y, :Z), title = "Phase space")
plot(plt1, plt2; layout = (1,2), size = (800,400))
plot!(bottom_margin = 3Plots.Measures.mm) # hide
```
