# [The Catalyst DSL - Advanced Features and Options](@id dsl_advanced_options)
Within the Catalyst DSL, each line can represent either *a reaction* or *an option*. The [previous DSL tutorial](@ref dsl_description) described how to create reactions. This one will focus on options. These are typically used to supply a model with additional information. Examples include the declaration of initial condition/parameter default values, or the creation of observables. 

All option designations begin with a declaration starting with `@`, followed by its input. E.g. the `@observables` option allows for the generation of observables. Each option can only be used once within each use of `@reaction_network`. This tutorial will also describe some additional advanced DSL features that do not involve using an option. 

As a first step, we import Catalyst (which is required to run the tutorial):
```@example dsl_advanced_explicit_definitions
using Catalyst
```

## [Explicit specification of network species and parameters](@id dsl_advanced_options_declaring_species_and_parameters)
Previously, we mentioned that the DSL automatically determines which symbols correspond to species and which to parameters. This is done by designating everything that appears as either a substrate or a product as a species, and all remaining quantities as parameters (i.e. those only appearing within rates or [stoichiometric constants](@ref dsl_description_stoichiometries_parameters)). Sometimes, one might want to manually override this default behaviour for a given symbol. I.e. consider the following model, where the conversion of a protein `P` from its inactive form (`Pᵢ`) to its active form (`Pₐ`) is catalysed by an enzyme (`E`). Using the most natural description:
```@example dsl_advanced_explicit_definitions
catalysis_sys = @reaction_network begin
    k*E, Pᵢ --> Pₐ
end
```
`E` (as well as `k`) will be considered a parameter, something we can confirm directly:
```@example dsl_advanced_explicit_definitions
parameters(catalysis_sys)
```
If we want `E` to be considered a species, we can designate this using the `@species` option:
```@example dsl_advanced_explicit_definitions
catalysis_sys = @reaction_network begin
    @species E(t)
    k*E, Pᵢ --> Pₐ
end
parameters(catalysis_sys)
```
!!! note
    When declaring species using the `@species` option, the species symbol must be followed by `(t)`. The reason is that species are time-dependent variables, and this time-dependency must be explicitly specified ([designation of non-`t` dependant species is also possible](@ref dsl_advanced_options_ivs)).

Similarly, the `@parameters` option can be used to explicitly designate something as a parameter:
```@example dsl_advanced_explicit_definitions
catalysis_sys = @reaction_network begin
    @parameters k
    k*E, Pᵢ --> Pₐ
end
```
Here, while `k` is explicitly defined as a parameter, no information is provided about `E`. Hence, the default case will be used (setting `E` to a parameter). The `@species` and `@parameter` options can be used simultaneously (although a quantity cannot be declared *both* as a species and a parameter). They may be followed by a full list of all species/parameters, or just a subset.

While designating something which would default to a parameter as a species is straightforward, the reverse (creating a parameter which occurs as a substrate or product) is more involved. This is, however, possible, and described [here](@ref dsl_advanced_options_constant_species).

Rather than listing all species/parameters on a single line after the options, a `begin ... end` block can be used (listing one species/parameter on each line). E.g. in the following example we use this notation to explicitly designate all species and parameters of the system:
```@example dsl_advanced_explicit_definitions
catalysis_sys = @reaction_network begin
    @species begin 
        E(t)
        Pᵢ(t)
        Pₐ(t)
    end
    @parameters begin
        k
    end
    k*E, Pᵢ --> Pₐ
end
```

A side-effect of using the `@species` and `@parameter` options is that they specify *the order in which the species and parameters are stored*. I.e. lets check the order of the parameters in the parameters in the following dimerisation model:
```@example dsl_advanced_explicit_definitions
dimerisation = @reaction_network begin
    (p,d), 0 <--> X
    (kB,kD), 2X <--> X2
end
parameters(dimerisation)
```
The default order is typically equal to the order with which the parameters (or species) are encountered in the DSL (this is, however, not guaranteed). If we specify the parameters using `@parameters`, the order used within the option is used instead:
```@example dsl_advanced_explicit_definitions
dimerisation = @reaction_network begin
    @parameters kB kD p d
    (p,d), 0 <--> X
    (kB,kD), 2X <--> X2
end
parameters(dimerisation)
```
!!! danger
    Generally, Catalyst and the SciML ecosystem *do not* guarantee that parameter and species order are preserved throughout various operations on a model. Writing programs that depend on these orders is *strongly discouraged*. There are, however, some legacy packages which still depend on order (one example can be found [here](@ref optimization_parameter_fitting_basics)). In these situations, this might be useful. However, in these cases, it is recommended that the user is extra wary, and also checks the order manually. 

!!! note
    The syntax of the `@species` and `@parameters` options is identical to that used by the `@species` and `@parameters` macros [used in programmatic modelling in Catalyst](@ref programmatic_CRN_construction) (for e.g. designating metadata or initial conditions). Hence, if one has learnt how to specify species/parameters using either approach, that knowledge can be transferred to the other one.

Generally, there are four main reasons for specifying species/parameters using the `@species` and `@parameters` options:
1. To designate a quantity, that would otherwise have defaulted to a parameter, as a species.
2. To designate default values for parameters/species initial conditions (described [here](@ref dsl_advanced_options_default_vals)).
3. To designate metadata for species/parameters (described [here](@ref dsl_advanced_options_species_and_parameters_metadata)).
4. To designate a species or parameters that do not occur in reactions, but are still part of the model (e.g a [parametric initial condition](@ref dsl_advanced_options_parametric_initial_conditions))

!!! warning
    Catalyst's DSL automatically infer species and parameters from the input. However, it only does so for *quantities that appear in reactions*. Until now this has not been relevant. However, this tutorial will demonstrate cases where species/parameters that are not part of reactions are used. These *must* be designated using either the `@species` or `@parameters` options (or the `@variables` option, which is described [later](@ref constraint_equations)).

### [Setting default values for species and parameters](@id dsl_advanced_options_default_vals)
When declaring species/parameters using the `@species` and `@parameters` options, one can also assign them default values (by appending them with `=` followed by the desired default value). E.g here we set `X`'s default initial condition value to $1.0$, and `p` and `d`'s default values to $1.0$ and $0.2$, respectively:
```@example dsl_advanced_defaults
using Catalyst # hide
rn = @reaction_network begin
    @species X(t)=1.0
    @parameters p=1.0 d=0.1
    (p,d), 0 <--> X
end
```
Next, if we simulate the model, we do not need to provide values for species or parameters that have default values. In this case all have default values, so both `u0` and `ps` can be empty vectors:
```@example dsl_advanced_defaults
using OrdinaryDiffEqDefault, Plots
u0 = []
tspan = (0.0, 10.0)
p = []
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```
It is still possible to provide values for some (or all) initial conditions/parameters in `u0`/`ps` (in which case these overrides the default values):
```@example dsl_advanced_defaults
u0 = [:X => 4.0]
p = [:d => 0.5]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```
It is also possible to declare a model with default values for only some initial conditions/parameters:
```@example dsl_advanced_defaults
using Catalyst # hide
rn = @reaction_network begin
    @species X(t)=1.0
    (p,d), 0 <--> X
end

tspan = (0.0, 10.0)
p = [:p => 1.0, :d => 0.2]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob)
plot(sol)
```

### [Setting parametric initial conditions](@id dsl_advanced_options_parametric_initial_conditions)
In the previous section, we designated default values for initial conditions and parameters. However, the right-hand side of the designation accepts any valid expression (not only numeric values). While this can be used to set up some advanced default values, the most common use case is to designate a species's initial condition as a parameter. E.g. in the following example we represent the initial condition of `X` using the parameter `X₀`. 
```@example dsl_advanced_defaults
rn = @reaction_network begin
    @species X(t)=X₀
    @parameters X₀
    (p,d), 0 <--> X
end
```
Please note that as the parameter `X₀` does not occur as part of any reactions, Catalyst's DSL cannot infer whether it is a species or a parameter. This must hence be explicitly declared. We can now simulate our model while providing `X`'s value through the `X₀` parameter:
```@example dsl_advanced_defaults
using OrdinaryDiffEqTsit5
u0 = []
p = [:X₀ => 1.0, :p => 1.0, :d => 0.5]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob, Tsit5())
plot(sol)
```
It is still possible to designate $X$'s value in `u0`, in which case this overrides the default value.
```@example dsl_advanced_defaults
u0 = [:X => 0.5]
p = [:X₀ => 1.0, :p => 1.0, :d => 0.5]
oprob = ODEProblem(rn, u0, tspan, p)
sol = solve(oprob, Tsit5())
plot(sol)
```
Please note that `X₀` is still a parameter of the system, and as such its value must still be designated to simulate the model (even if it is not actually used).

### [Designating metadata for species and parameters](@id dsl_advanced_options_species_and_parameters_metadata)
Catalyst permits the user to define *metadata* for species and parameters. This permits the user to assign additional information to these, which can be used for a variety of purposes. Some Catalyst features depend on using metadata (with each such case describing specifically how this is done). Here we will introduce how to set metadata, and describe some common metadata types. 

Whenever a species/parameter is declared using the `@species`/`@parameters` options, it can be followed by a `[]` within which the metadata is given. Each metadata entry consists of the metadata's name, followed by a `=`, followed by its value. E.g. the `description` metadata allows you to attach a [`String`](https://docs.julialang.org/en/v1/base/strings/) to a species/parameter. Here we create a simple model where we add descriptions to all species and parameters.
```@example dsl_advanced_metadata
using Catalyst # hide
two_state_system = @reaction_network begin
    @species Xᵢ(t) [description="X's inactive form"] Xₐ(t) [description=" X's active form"]
    @parameters kA [description="Activation rate"] kD [description="Deactivation rate"]
    (ka,kD), Xᵢ <--> Xₐ
end
```
A metadata can be given to only a subset of a system's species/parameters, and a quantity can be given several metadata entries. To give several metadata, separate each by a `,`. Here we only provide a description for `kA`, for which we also provide a [`bounds` metadata](https://docs.sciml.ai/ModelingToolkit/dev/basics/Variable_metadata/#Bounds),
```@example dsl_advanced_metadata
two_state_system = @reaction_network begin
    @parameters kA [description="Activation rate", bounds=(0.01,10.0)]
    (ka,kD), Xᵢ <--> Xₐ
end
```

It is possible to add both default values and metadata to a parameter/species. In this case, first provide the default value, and next the metadata. I.e. to in the above example set $kA$'s default value to $1.0$ we use
```@example dsl_advanced_metadata
two_state_system = @reaction_network begin
    @parameters kA=1.0 [description="Activation rate", bounds=(0.01,10.0)]
    (ka,kD), Xᵢ <--> Xₐ
end
```

When designating metadata for species/parameters in `begin ... end` blocks the syntax changes slightly. Here, a `,` must be inserted before the metadata (but after any potential default value). I.e. a version of the previous example can be written as
```@example dsl_advanced_metadata
two_state_system = @reaction_network begin
    @parameters begin
        kA, [description="Activation rate", bounds=(0.01,10.0)]
        kD = 1.0, [description="Deactivation rate"]
    end
    (kA,kD), Xᵢ <--> Xₐ
end
```

Each metadata has its own getter functions. E.g. we can get the description of the parameter `kA` using `ModelingToolkit.getdescription`:
```@example dsl_advanced_metadata
ModelingToolkit.getdescription(two_state_system.kA)
```

### [Designating constant-valued/fixed species parameters](@id dsl_advanced_options_constant_species)
Catalyst enables the designation of parameters as `constantspecies`. These parameters can be used as species in reactions, however, their values are not changed by the reaction and remain constant throughout the simulation (unless changed by e.g. the [occurrence of an event](@ref constraint_equations_events). Practically, this is done by setting the parameter's `isconstantspecies` metadata to `true`. Here, we create a simple reaction where the species `X` is converted to `Xᴾ` at rate `k`. By designating `X` as a constant species parameter, we ensure that its quantity is unchanged by the occurrence of the reaction.
```@example dsl_advanced_constant_species
using Catalyst # hide
rn = @reaction_network begin
    @parameters X [isconstantspecies=true]
    k, X --> Xᴾ
end
```
We can confirm that $Xᴾ$ is the only species of the system:
```@example dsl_advanced_constant_species
species(rn)
```
Here, the produced model is actually identical to if $X$ had simply been a parameter in the reaction's rate:
```@example dsl_advanced_constant_species
rn = @reaction_network begin
    k*X, 0 --> Xᴾ
end
```

A common use-case for constant species is when modelling systems where some species are present in such surplus that their amounts the reactions' effect on it is negligible. A system which is commonly modelled this way is the [Brusselator](https://en.wikipedia.org/wiki/Brusselator).

### [Designating parameter types](@id dsl_advanced_options_parameter_types)
Sometimes it is desired to designate that a parameter should have a specific [type](https://docs.julialang.org/en/v1/manual/types/). When supplying this parameter's value to e.g. an `ODEProblem`, that parameter will then be restricted to that specific type. Designating a type is done by appending the parameter with `::` followed by its type. E.g. in the following example we specify that the parameter `n` (the number of `X` molecules in the `Xn` polymer) must be an integer (`Int64`)
```@example dsl_advanced_parameter_types
using Catalyst # hide
polymerisation_network = @reaction_network begin
    @parameters n::Int64
    (kB,kD), n*X <--> Xn
end
nothing # hide
```
Generally, when simulating models with mixed parameter types, it is recommended to [declare parameter values as tuples, rather than vectors](@ref simulation_intro_ODEs_input_forms), e.g.:
```@example dsl_advanced_parameter_types
ps = (:kB => 0.2, :kD => 1.0, :n => 2)
nothing # hide
```

If a parameter has a type, metadata, and a default value, they are designated in the following order:
```@example dsl_advanced_parameter_types
polymerisation_network = @reaction_network begin
    @parameters n::Int64 = 2 [description="Parameter n, an integer with defaults value 2."]
    (kB,kD), n*X <--> Xn
end
nothing # hide
```

### [Vector-valued species or parameters](@id dsl_advanced_options_vector_variables)
Sometimes, one wishes to declare a large number of similar parameters or species. This can be done by *creating them as vectors*. E.g. below we create a [two-state system](@ref basic_CRN_library_two_states). However, instead of declaring `X1` and `X2` (and `k1` and `k2`) as separate entities, we declare them as vectors:
```@example dsl_advanced_vector_variables
using Catalyst # hide
two_state_model = @reaction_network begin
    @parameters k[1:2]
    @species (X(t))[1:2]
    (k[1],k[2]), X[1] <--> X[2]
end
```
Now, we can also declare our initial conditions and parameter values as vectors as well:
```@example dsl_advanced_vector_variables
using OrdinaryDiffEqDefault, Plots # hide
u0 = [:X => [0.0, 2.0]]
tspan = (0.0, 1.0)
ps = [:k => [1.0, 2.0]]
oprob = ODEProblem(two_state_model, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```

### [Turning off species, parameter, and variable inference](@id dsl_advanced_options_require_dec)
In some cases it may be desirable for Catalyst to not infer species and parameters from the DSL, as in the case of reaction networks with very many variables, or as a sanity check that variable names are written correctly. To turn off inference, simply use the `@require_declaration` option when using the `@reaction_network` DSL. This will require every single variable, species, or parameter used within the DSL to be explicitly declared using the `@variable`, `@species`, or `@parameter` options. In the case that the DSL parser encounters an undeclared symbolic, it will error with an `UndeclaredSymbolicError` and print the reaction or equation that the undeclared symbolic was found in. 

```julia
using Catalyst
rn = @reaction_network begin
    @require_declaration
    (k1, k2), A <--> B
end
```
Running the code above will yield the following error: 
```
LoadError: UndeclaredSymbolicError: Unrecognized variables A detected in reaction expression: "((k1, k2), A <--> B)". Since the flag @require_declaration is declared, all species must be explicitly declared with the @species macro.
```
In order to avoid the error in this case all the relevant species and parameters will have to be declared.
```@example dsl_advanced_require_dec
# The following case will not error. 
using Catalyst
t = default_t()
rn = @reaction_network begin
    @require_declaration
    @species A(t) B(t)
    @parameters k1 k2
    (k1, k2), A <--> B
end
```

The following cases in which the DSL would normally infer variables will all throw errors if `@require_declaration` is set and the variables are not explicitly declared.
- Occurrence of an undeclared species in a reaction, as in the example above.
- Occurrence of an undeclared parameter in a reaction rate expression, as in the reaction line `k*n, A --> B`.
- Occurrence of an undeclared parameter in the stoichiometry of a species, as in the reaction line `k, n*A --> B`.
- Occurrence of an undeclared differential variable on the LHS of a coupled differential equation, as in `A` in `@equations D(A) ~ A^2`.
- Occurrence of an undeclared [observable](@ref dsl_advanced_options_observables) in an `@observables` expression, such as `@observables X1 ~ A + B`.

## [Naming reaction networks](@id dsl_advanced_options_naming)
Each reaction network model has a name. It can be accessed using the `nameof` function. By default, some generic name is used:
```@example dsl_advanced_names
using Catalyst # hide
rn = @reaction_network begin
    (p,d), 0 <--> X
end
nameof(rn)
```
A specific name can be given as an argument between the `@reaction_network` and the `begin`. E.g. to name a network `my_network` we can use:
```@example dsl_advanced_names
rn = @reaction_network my_network begin
    (p,d), 0 <--> X
end
nameof(rn)
```

A consequence of generic names being used by default is that networks, even if seemingly identical, by default are not. E.g.
```@example dsl_advanced_names
rn1 = @reaction_network begin
    (p,d), 0 <--> X
end
rn2 = @reaction_network begin
    (p,d), 0 <--> X
end
rn1 == rn2
```
The reason can be confirmed by checking that their respective (randomly generated) names are different:
```@example dsl_advanced_names
nameof(rn1) == nameof(rn2)
```
By designating the networks to have the same name, however, identity is achieved.
```@example dsl_advanced_names
rn1 = @reaction_network my_network begin
    (p,d), 0 <--> X
end
rn2 = @reaction_network my_network begin
    (p,d), 0 <--> X
end
rn1 == rn2
```
If you wish to check for identity, and wish that models that have different names but are otherwise identical, should be considered equal, you can use the [`isequivalent`](@ref) function.

Setting model names is primarily useful for [hierarchical modelling](@ref compositional_modeling), where network names are appended to the display names of subnetworks' species and parameters.

## [Creating observables](@id dsl_advanced_options_observables)
Sometimes one might want to use observable variables. These are variables with values that can be computed directly from a system's state (rather than having their values implicitly given by reactions or equations). Observables can be designated using the `@observables` option. Here, the `@observables` option is followed by a `begin ... end` block with one line for each observable. Each line first gives the observable, followed by a `~` (*not* a `=`!), followed by an expression describing how to compute it.

Let us consider a model where two species (`X` and `Y`) can bind to form a complex (`XY`, which also can dissociate back into `X` and `Y`). If we wish to create a representation for the total amount of `X` and `Y` in the system, we can do this by creating observables `Xtot` and `Ytot`:
```@example dsl_advanced_observables
using Catalyst # hide
rn = @reaction_network begin
    @observables begin
        Xtot ~ X + XY
        Ytot ~ Y + XY
    end
    (kB,kD), X + Y <--> XY
end
```
We can now simulate our model using normal syntax (initial condition values for observables should not, and can not, be provided):
```@example dsl_advanced_observables
using OrdinaryDiffEqTsit5
u0 = [:X => 1.0, :Y => 2.0, :XY => 0.0]
tspan = (0.0, 10.0)
ps = [:kB => 1.0, :kD => 1.5]
oprob = ODEProblem(rn, u0, tspan, ps)
sol = solve(oprob, Tsit5())
nothing # hide
```

Next, we can use [symbolic indexing](@ref simulation_structure_interfacing) of our solution object, but with the observable as input. E.g. we can use 
```@example dsl_advanced_observables
sol[:Xtot]
nothing # hide
```
to get a vector with `Xtot`'s value throughout the simulation. We can also use
```@example dsl_advanced_observables
using Plots
plot(sol; idxs = :Xtot)
plot!(ylimit = (minimum(sol[:Xtot])*0.95, maximum(sol[:Xtot])*1.05)) # hide
```
to plot the observables (rather than the species).

Observables can be defined using complicated expressions containing species, parameters, and [variables](@ref constraint_equations) (but not other observables). In the following example (which uses a [parametric stoichiometry](@ref dsl_description_stoichiometries_parameters)) `X` polymerises to form a complex `Xn` containing `n` copies of `X`. Here, we create an observable describing the total number of `X` molecules in the system:
```@example dsl_advanced_observables
rn = @reaction_network begin
    @observables Xtot ~ X + n*Xn
    (kB,kD), n*X <--> Xn
end
nothing # hide
```
!!! note
    If only a single observable is declared, the `begin ... end` block is not required and the observable can be declared directly after the `@observables` option.

[Metadata](@ref dsl_advanced_options_species_and_parameters_metadata) can be supplied to an observable directly after its declaration (but before its formula). If so, the metadata must be separated from the observable with a `,`, and the observable plus the metadata encapsulated by `()`. E.g. to add a [description metadata](@ref dsl_advanced_options_species_and_parameters_metadata) to our observable we can use
```@example dsl_advanced_observables
rn = @reaction_network begin
    @observables (Xtot, [description="The total amount of X in the system."]) ~ X + n*Xn
    (kB,kD), n*X <--> Xn
end
nothing # hide
```

Observables are by default considered [variables](@ref constraint_equations) (not species). To designate them as a species, they can be pre-declared using the `@species` option. I.e. Here `Xtot` becomes a species:
```@example dsl_advanced_observables
rn = @reaction_network begin
    @species Xtot(t)
    @observables Xtot ~ X + n*Xn  
    (kB,kD), n*X <--> Xn
end
nothing # hide
```

Some final notes regarding observables:
- The left-hand side of the observable declaration must contain a single symbol only (with the exception of metadata, which can also be supplied).
- All quantities appearing on the right-hand side must be declared elsewhere within the `@reaction_network` call (either by being part of a reaction, or through the `@species`, `@parameters`, or `@variables` options).
- Observables may not depend on other observables.
- Observables have their dependent variable(s) automatically assigned as the union of the dependent variables of the species and variables on which it depends.

## [Specifying non-time independent variables](@id dsl_advanced_options_ivs)

Catalyst's `ReactionSystem` models depend on a *time independent variable*, and potentially one or more *spatial independent variables*. By default, the independent variable `t` is used. We can declare another independent variable (which is automatically used as the default one) using the `@ivs` option. E.g. to use `τ` instead of `t` we can use
```@example dsl_advanced_ivs
using Catalyst # hide
rn = @reaction_network begin
    @ivs τ
    (ka,kD), Xᵢ <--> Xₐ
end
nothing # hide
```
We can confirm that `Xᵢ` and `Xₐ` depend on `τ` (and not `t`):
```@example dsl_advanced_ivs
species(rn)
```

It is possible to designate several independent variables using `@ivs`. If so, the first one is considered the default (time) independent variable, while the following one(s) are considered spatial independent variable(s). If we want some species to depend on a non-default independent variable, this has to be explicitly declared:
```@example dsl_advanced_ivs
rn = @reaction_network begin
    @ivs τ x
    @species X(τ) Y(x)
    (p1,d1), 0 <--> X
    (p2,d2), 0 <--> Y
end
species(rn)
```
It is also possible to have species which depends on several independent variables:
```@example dsl_advanced_ivs
rn = @reaction_network begin
    @ivs t x
    @species Xᵢ(t,x) Xₐ(t,x)
    (ka,kD), Xᵢ <--> Xₐ
end
species(rn)
```

!!! note
    Setting spatial independent variables is primarily intended for the modelling of spatial systems on continuous domains. Catalyst's support for this is currently under development. Hence, the utility of specifying spatial independent variables is limited.


## [Setting reaction metadata](@id dsl_advanced_options_reaction_metadata)
It is possible to supply reactions with *metadata*, containing some additional information of the reaction. A reaction's metadata follows after its declaration (first using the metadata's name, then a `=`, then its value) and is encapsulated by `[]` (where individual entries are separated by `,`). Here, we add a `description` metadata to the reactions of a birth-death process:
```@example dsl_advanced_reaction_metadata
using Catalyst # hide
bd_model = @reaction_network begin
    p, 0 --> X, [description="Production reaction"]
    d, X --> 0, [description="Degradation reaction"]
end
nothing # hide
```

When [bundling reactions](@ref dsl_description_reaction_bundling), reaction metadata can be bundled using the same rules as rates. Bellow we re-declare our birth-death process, but on a single line:
```@example dsl_advanced_reaction_metadata
bd_model = @reaction_network begin
    (p,d), 0 <--> X, ([description="Production reaction"], [description="Degradation reaction"])
end
nothing # hide
```

Here we declare a model where we also provide a `misc` metadata (which can hold any quantity we require) to our birth reaction:
```@example dsl_advanced_reaction_metadata
bd_model = @reaction_network begin
    p, 0 --> X, [description="Production reaction", misc=:value]
    d, X --> 0, [description="Degradation reaction"]
end
nothing # hide
```

A reaction's metadata can be accessed using specific functions, e.g. `Catalyst.hasdescription` and `Catalyst.getdescription` can be used to check if a reaction have a description metadata, and to retrieve it, respectively:
```@example dsl_advanced_reaction_metadata
rx = @reaction p, 0 --> X, [description="A production reaction"]
Catalyst.getdescription(rx)
```

A list of all available reaction metadata can be found [in the api](@ref api_rx_metadata).

## [Working with symbolic variables and the DSL](@id dsl_advanced_options_symbolics_and_DSL)
We have previously described how Catalyst represents its models symbolically (enabling e.g. symbolic differentiation of expressions stored in models). While Catalyst utilises this for many internal operation, these symbolic representations can also be accessed and harnessed by the user. Primarily, doing so is much easier during programmatic (as opposed to DSL-based) modelling. Indeed, the section on [programmatic modelling](@ref programmatic_CRN_construction) goes into more details about symbolic representation in models, and how these can be used. It is, however, also ways to utilise these methods during DSL-based modelling. Below we briefly describe two methods for doing so.

### [Using `@unpack` to extract symbolic variables from `ReactionSystem`s](@id dsl_advanced_options_symbolics_and_DSL_unpack)
Let us consider a simple [birth-death process](@ref basic_CRN_library_bd) created using the DSL:
```@example dsl_advanced_programmatic_unpack
using Catalyst # hide
bd_model = @reaction_network begin
    (p,d), 0 <--> X
end
nothing # hide
```
Since we have not explicitly declared `p`, `d`, and `X` using `@parameters` and `@species`, we cannot represent these symbolically (only using `Symbol`s). If we wish to do so, however, we can fetch these into our current scope using the `@unpack` macro:
```@example dsl_advanced_programmatic_unpack
@unpack p, d, X = bd_model
nothing # hide
```
This lists first the quantities we wish to fetch (does not need to be the model's full set of parameters and species), then `=`, followed by the model variable. `p`, `d` and `X` are now symbolic variables in the current scope, just as if they had been declared using `@parameters` or `@species`. We can confirm this:
```@example dsl_advanced_programmatic_unpack
X
```
Next, we can now use these to e.g. designate initial conditions and parameter values for model simulations:
```@example dsl_advanced_programmatic_unpack
using OrdinaryDiffEqDefault, Plots # hide
u0 = [X => 0.1]
tspan = (0.0, 10.0)
ps = [p => 1.0, d => 0.2]
oprob = ODEProblem(bd_model, u0, tspan, ps)
sol = solve(oprob)
plot(sol)
```

!!! warning
    Just like when using `@parameters` and `@species`, `@unpack` will overwrite any variables in the current scope which share name with the imported quantities.

### [Interpolating variables into the DSL](@id dsl_advanced_options_symbolics_and_DSL_interpolation)
Catalyst's DSL allows Julia variables to be interpolated for the network name, within rate constant expressions, or for species/stoichiometries within reactions. Using the lower-level symbolic interface we can then define symbolic variables and parameters outside of `@reaction_network`, which can then be used within expressions in the DSL. 

Interpolation is carried out by pre-appending the interpolating variable with a `$`. For example, here we declare the parameters and species of a birth-death model, and interpolate these into the model:
```@example dsl_advanced_programmatic_interpolation
using Catalyst # hide
t = default_t()
@species X(t)
@parameters p d
bd_model = @reaction_network begin
    ($p, $d), 0 <--> $X
end
```
Additional information (such as default values or metadata) supplied to `p`, `d`, and `X` is carried through to the DSL. However, interpolation for this purpose is of limited value, as such information [can be declared within the DSL](@ref dsl_advanced_options_declaring_species_and_parameters). However, it is possible to interpolate larger algebraic expressions into the DSL, e.g. here
```@example dsl_advanced_programmatic_interpolation
@species X1(t) X2(t) X3(t) E(t)
@parameters d
d_rate = d/(1 + E)
degradation_model = @reaction_network begin
    $d_rate, X1 --> 0
    $d_rate, X2 --> 0
    $d_rate, X3 --> 0
end
```
we declare an expression `d_rate`, which then can be inserted into the DSL via interpolation.

It is also possible to use interpolation in combination with the `@reaction` macro. E.g. the reactions of the above network can be declared individually using
```@example dsl_advanced_programmatic_interpolation
rxs = [
    @reaction $d_rate, $X1 --> 0
    @reaction $d_rate, $X2 --> 0
    @reaction $d_rate, $X3 --> 0
]
nothing # hide 
```

!!! note
    When using interpolation, expressions like `2$spec` won't work; the multiplication symbol must be explicitly included like `2*$spec`.

## [Creating individual reaction using the `@reaction` macro](@id dsl_advanced_options_reaction_macro)
Catalyst exports a macro `@reaction`, which can be used to generate a singular [`Reaction`](@ref) object of the same type which is stored within the [`ReactionSystem`](@ref) structure (which in turn can be generated by `@reaction_network`). Generally, `@reaction` follows [identical rules to those of `@reaction_network`](@ref dsl_description_reactions) for writing and interpreting reactions (however, bi-directional reactions are not permitted). E.g. here we create a simple dimerisation reaction:
```@example dsl_advanced_reaction_macro
using Catalyst # hide
rx_dimerisation = @reaction kD, 2X --> X2
```
Here, `@reaction` is followed by a single line consisting of three parts:
- A rate (at which the reaction occurs).
- Any number of substrates (which are consumed by the reaction).
- Any number of products (which are produced by the reaction).

In the next example, we first create a simple [SIR model](@ref basic_CRN_library_sir). Next, we instead create its individual reaction components using the `@reaction` macro. Finally, we confirm that these are identical to those stored in the initial model (using the [`reactions`](@ref) function).
```@example dsl_advanced_reaction_macro
sir_model = @reaction_network begin
 α, S + I --> 2I
 β, I --> R
end
infection_rx = @reaction α, S + I --> 2I
recovery_rx = @reaction β, I --> R
sir_rxs = [infection_rx, recovery_rx]
issetequal(reactions(sir_model), sir_rxs)
```
One of the primary uses of the `@reaction` macro is to provide some of the convenience of the DSL to [*programmatic modelling](@ref programmatic_CRN_construction). E.g. here we can combine our reactions to create a `ReactionSystem` directly, and also confirm that this is identical to the model created through the DSL:
```@example dsl_advanced_reaction_macro
sir_programmatic = complete(ReactionSystem(sir_rxs, default_t(); name = :sir))
sir_programmatic == sir_model
```

During programmatic modelling, it can be good to keep in mind that already declared symbolic variables can be [*interpolated*](@ref dsl_advanced_options_symbolics_and_DSL_interpolation). E.g. here we create two production reactions both depending on the same Michaelis-Menten function:
```@example dsl_advanced_reaction_macro
t = default_t()
@species X(t)
@parameters v K
mm_term = Catalyst.mm(X, v, K)
rx1 = @reaction $mm_term, 0 --> P1
rx2 = @reaction $mm_term, 0 --> P2
nothing # hide
```

## [Disabling mass action for reactions](@id dsl_advanced_options_disable_ma)

As [described previously](@ref math_models_in_catalyst_rre_odes), Catalyst uses *mass action kinetics* to generate ODEs from reactions. Here, each reaction generates a term for each of its reactants, which consists of the reaction's rate, substrates, and the reactant's stoichiometry. E.g. the following reaction:
```@example dsl_advanced_disable_ma
using Catalyst # hide
rn = @reaction_network begin
  k, X --> ∅
end
```
generates a single term $-k*[X]$:
```@example dsl_advanced_disable_ma
using Latexify
latexify(rn; form = :ode)
```

It is possible to remove the substrate contribution by using any of the following non-filled arrows when declaring the reaction: `<=`, `⇐`, `⟽`, `=>`, `⇒`, `⟾`, `⇔`, `⟺`. This means that the reaction
```@example dsl_advanced_disable_ma
rn = @reaction_network begin
  k, X => ∅
end
latexify(rn; form = :ode)
```
will occur at rate $d[X]/dt = -k$ (which might become a problem since $[X]$ will be degraded at a constant rate even when very small or equal to 0). This functionality allows the user to fully customise the ODEs generated by their models. 

Note, stoichiometric coefficients are still included, i.e. the reaction
```@example dsl_advanced_disable_ma
rn = @reaction_network begin
  k, 2*X ⇒ ∅
end
latexify(rn; form = :ode)
```
has rate $d[X]/dt = -2 k$.

