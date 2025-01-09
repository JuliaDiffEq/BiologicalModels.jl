# [Chemical Reaction Network Theory](@id network_analysis_structural_aspects)
The reaction complex representation can be exploited via [Chemical Reaction
Network Theory](https://en.wikipedia.org/wiki/Chemical_reaction_network_theory)
to provide insight into possible steady state and time-dependent properties of
RRE ODE models and stochastic chemical kinetics models. We'll now illustrate
some of the types of network properties that Catalyst can determine, using the
reaction complex representation in these calculations.

Consider the following reaction network.
```@example s1
rn = @reaction_network begin
    (k1,k2), A + B <--> C
    k3, C --> D+E
    (k4,k5), D+E <--> F
    (k6,k7), 2A <--> B+G
    k8, B+G --> H
    k9, H --> 2A
end
```
with graph
```julia
complexgraph(rn)
```

![network_1](../assets/complex_rn.svg)

# [Linkage classes and sub-networks of the reaction network](@id network_analysis_structural_aspects_linkage)
The preceding reaction complex graph shows that `rn` is composed of two
disconnected sub-graphs, one containing the complexes ``A+B``, ``C``, ``D+E``, and
``F``, the other containing the complexes ``2A``, ``B + G``, and ``H``. These sets,
``\{A+B, C, D+E, F\}`` and ``\{2A, B + G,H\}`` are called the "linkage classes"
of the reaction network. The function [`linkageclasses`](@ref) will calculate
these for a given network, returning a vector of the integer indices of reaction
complexes participating in each set of linkage-classes. Note, indices of
reaction complexes can be determined from the ordering returned by
[`reactioncomplexes`](@ref).
```@example s1
# we must first calculate the reaction complexes -- they are cached in rn
reactioncomplexes(rn)

# now we can calculate the linkage classes
lcs = linkageclasses(rn)
```
It can often be convenient to obtain the disconnected sub-networks as distinct
`ReactionSystem`s, which are returned by the [`subnetworks`](@ref) function:
```@example s1
subnets = subnetworks(rn)

# check the reactions in each subnetwork
reactions.(subnets)
```
The graphs of the reaction complexes in the two sub-networks are then
```julia
  complexgraph(subnets[1])
```

![subnetwork_1](../assets/complex_subnets1.svg)

and,
```julia
 complexgraph(subnets[2])
```

![subnetwork_2](../assets/complex_subnets2.svg)


# [Deficiency of the network](@id network_analysis_structural_aspects_deficiency)
A famous theorem in Chemical Reaction Network Theory, the Deficiency Zero
Theorem [^1], allows us to use knowledge of the net stoichiometry matrix and the
linkage classes of a *mass action* RRE ODE system to draw conclusions about the
system's possible steady states. In this section we'll see how Catalyst can
calculate a network's deficiency.

The rank, ``r``, of a reaction network is defined as the dimension of the
subspace spanned by the net stoichiometry vectors of the reaction-network [^1],
i.e. the span of the columns of the net stoichiometry matrix `N`. It corresponds
to the number of independent species in a chemical reaction network. That is, if
we calculate the linear conservation laws of a network, and use them to
eliminate the dependent species of the network, we will have ``r`` independent
species remaining. For our current example the conservation laws are given by
```@example s1
# first we calculate the conservation laws -- they are cached in rn
conservationlaws(rn)

# then we display them as equations for the dependent variables
conservedequations(rn)
show(stdout, MIME"text/plain"(), ans) # hide
```
Here the parameters `Γ[i]` represent the constants of the three
conservation laws, and we see that there are three dependent species that could
be eliminated. As
```@example s1
numspecies(rn)
```
we find that there are five independent species. Let's check this is correct:
```@example s1
using LinearAlgebra
rank(netstoichmat(rn)) == 5
```
So we know that the rank of our reaction network is five. An extended section discussing how to
utilise conservation law elimination during chemical reaction network modelling can be found
[here](@ref conservation_laws).

The deficiency, ``\delta``, of a reaction network is defined as
```math
\delta = \textrm{(number of complexes)} - \textrm{(number of linkage classes)} - \textrm{(rank)}.
```
For our network this is ``7 - 2 - 5 = 0``, which we can calculate in Catalyst as
```@example s1
# first we calculate the reaction complexes of rn and cache them in rn
reactioncomplexes(rn)

# then we can calculate the deficiency
δ = deficiency(rn)
```
Quoting Feinberg [^1]
> Deficiency zero networks are ones for which the reaction vectors [i.e. net
> stoichiometry vectors] are as independent as the partition of complexes into
> linkage classes will allow.


# [Reversibility of the network](@id network_analysis_structural_aspects_reversibility)
A reaction network is *reversible* if the "arrows" of the reactions are
symmetric so that every reaction is accompanied by its reverse reaction.
Catalyst's API provides the [`isreversible`](@ref) function to determine whether
a reaction network is reversible. As an example, consider
```@example s1
rn = @reaction_network begin
  (k1,k2),A <--> B
  (k3,k4),A + C <--> D
  (k5,k6),D <--> B+E
  (k7,k8),B+E <--> A+C
end

# calculate the set of reaction complexes
reactioncomplexes(rn)

# test if the system is reversible
isreversible(rn)
```
Consider another example,
```@example s1
rn = @reaction_network begin
  (k1,k2),A <--> B
  k3, A + C --> D
  k4, D --> B+E
  k5, B+E --> A+C
end
reactioncomplexes(rn)
isreversible(rn)
```
```julia
complexgraph(rn)
```

![reversibility](../assets/complex_reversibility.svg)

It is evident from the preceding graph that the network is not reversible.
However, it satisfies a weaker property in that there is a path from each
reaction complex back to itself within its associated subgraph. This is known as
*weak reversibility*. One can test a network for weak reversibility by using
the [`isweaklyreversible`](@ref) function:
```@example s1
# need subnetworks from the reaction network first
subnets = subnetworks(rn)
isweaklyreversible(rn, subnets)
```
Every reversible network is also weakly reversible, but not every weakly
reversible network is reversible.

# [Deficiency Zero Theorem](@id network_analysis_structural_aspects_deficiency_zero_theorem)
Knowing the deficiency and weak reversibility of a mass action chemical reaction
network ODE model allows us to make inferences about the corresponding
steady state behavior. Before illustrating how this works for one example, we
need one last definition.

Recall that in the matrix-vector representation for the RRE ODEs, the entries,
``N_{m k}``, of the stoichiometry matrix, ``N``, give the net change in species
``m`` due to reaction ``k``. If we let ``\mathbf{N}_k`` denote the ``k``th
column of this matrix, this vector corresponds to the change in the species
state vector, ``\mathbf{x}(t)``, due to reaction ``k``, i.e. when reaction ``k``
occurs ``\mathbf{x}(t) \to \mathbf{x}(t) + \mathbf{N}_k``. Moreover, by
integrating the ODE
```math
\frac{d\mathbf{x}}{dt} = N \mathbf{v}(\mathbf{x}(t)) = \sum_{k=1}^{K} v_k(\mathbf{x}(t)) \, \mathbf{N}_k
```
we find
```math
\mathbf{x}(t) = \mathbf{x}(0) + \sum_{k=1}^K \left(\int_0^t v_k(\mathbf{x})(s) \, ds\right) \mathbf{N}_k,
```
which demonstrates that ``\mathbf{x}(t) - \mathbf{x}(0)`` is always given by a
linear combination of the stoichiometry vectors, i.e.
```math
\mathbf{x}(t) - \mathbf{x}(0) \in \operatorname{span}\{\mathbf{N}_k \}.
```
In particular, this says that ``\mathbf{x}(t)`` lives in the translation of the
``\operatorname{span}\{\mathbf{N}_k \}`` by ``\mathbf{x}(0)`` which we write as
``(\mathbf{x}(0) + \operatorname{span}\{\mathbf{N}_k\})``. In fact, since the
solution should stay non-negative, if we let $\bar{\mathbb{R}}_+^{M}$ denote the
subset of vectors in $\mathbb{R}^{M}$ with non-negative components, the possible
physical values for the solution, ``\mathbf{x}(t)``, must be in the set
```math
(\mathbf{x}(0) + \operatorname{span}\{\mathbf{N}_k\}) \cap \bar{\mathbb{R}}_+^{M}.
```
This set is called the stoichiometric compatibility class of ``\mathbf{x}(t)``.
The key property of stoichiometric compatibility classes is that they are
invariant under the RRE ODE's dynamics, i.e. a solution will always remain
within the subspace given by the stoichiometric compatibility class. Finally, we
note that the *positive* stoichiometric compatibility class generated by
$\mathbf{x}(0)$ is just ``(\mathbf{x}(0) + \operatorname{span}\{\mathbf{N}_k\})
\cap \mathbb{R}_+^{M}``, where ``\mathbb{R}_+^{M}`` denotes the vectors in
``\mathbb{R}^M`` with strictly positive components.

With these definitions we can now see how knowing the deficiency and weak
reversibility of the network can tell us about its steady state behavior.
Consider the previous example, which we know is weakly reversible. Its
deficiency is
```@example s1
deficiency(rn)
```
We also verify that the system is purely mass action (though it is apparent
from the network's definition):
```@example s1
all(rx -> ismassaction(rx, rn), reactions(rn))
```
We can therefore apply the Deficiency Zero Theorem to draw conclusions about the
system's steady state behavior. The Deficiency Zero Theorem (roughly) says that
a mass action network with deficiency zero satisfies
1. If the network is weakly reversible, then independent of the reaction rate
   constants the RRE ODEs have exactly one equilibrium solution within each
   positive stoichiometric compatibility class. That equilibrium is locally
   asymptotically stable.
2. If the network is not weakly reversible, then the RRE ODEs cannot admit a
   positive equilibrium solution.

See [^1] for a more precise statement, proof, and additional examples.

We can therefore conclude that for any initial condition that is positive, and
hence in some positive stoichiometric compatibility class, `rn` will have
exactly one equilibrium solution which will be positive and locally
asymptotically stable.

As a final example, consider the following network from [^1]
```@example s1
rn = @reaction_network begin
  (k1,k2),A <--> 2B
  (k3,k4), A + C <--> D
  k5, B+E --> C + D
end
reactioncomplexes(rn)
subnets = subnetworks(rn)
isma = all(rx -> ismassaction(rx,rn), reactions(rn))
def = deficiency(rn)
iswr = isweaklyreversible(rn, subnets)
isma,def,iswr
```
which we see is mass action and has deficiency zero, but is not weakly
reversible. As such, we can conclude that for any choice of rate constants the
RRE ODEs cannot have a positive equilibrium solution.

```@docs
satisfiesdeficiencyzero
```

# Deficiency One Theorem
Very analogous to the deficiency zero theorem is the deficiency one theorem. The deficiency one theorem applies to a network with the following properties: 
1. The deficiency of each *linkage class* of the network is at most 1, 
2. The sum of the linkage class deficiencies is the total deficiency of the network, and
3. Each linkage class has at most one terminal linkage class. 

If these conditions are met, then the network will have at most one steady state in each stoichiometric compatibility class for any choice of rate constants and parameters.
Unlike the deficiency zero theorem, networks obeying the deficiency one theorem are not guaranteed to have stable solutions.
```@docs
satisfiesdeficiencyone
```

# [Complex and Detailed Balance](@id complex_and_detailed_balance)
A reaction network's steady state is **complex-balanced** if the total production of each *complex* is zero at the steady state. A reaction network's steady state is **detailed balanced** if every reaction is balanced by its reverse reaction at the steady-state (this corresponds to the usual notion of chemical equilibrium; note that this requires every reaction be reversible). 

Note that detailed balance at a given steady state implies complex balance for that steady state, i.e. detailed balance is a stronger property than complex balance.

Remarkably, having just one positive steady state that is complex (detailed) balance implies that complex (detailed) balance obtains at *every* positive steady state, so we say that a network is complex (detailed) balanced if any one of its steady states are complex (detailed) balanced. Additionally, there will be exactly one steady state in every positive stoichiometric compatibility class, and this steady state is asymptotically stable. (For proofs of these results, please consult Martin Feinberg's *Foundations of Chemical Reaction Network Theory*[^1]). So knowing that a network is complex balanced is really quite powerful.

Let's check that the reaction network defined above is complex balanced by providing a set of rates:
```@example s1
rates = Dict([:k1 => 2.4, :k2 => 4., :k3 => 10., :k4 => 5.5, :k5 => 0.4])
iscomplexbalanced(rn, rates)
```

Complex balance obtains for some sets of rates but not others: 
```@example s1
```

We can do a similar check for detailed balance. Let us make the reaction network 
```@example s1
rn1 = @reaction_network begin
  (k1,k2),A <--> 2B
  (k3,k4), A + C <--> D
  (k5,k6), B + E --> C + D
end
isdetailedbalanced(rn, rates)
```

The reason that the deficiency zero theorem puts such strong restrictions on the steady state properties of the reaction network is because it implies that the reaction network will be complex balanced for any set of rate constants and parameters. The fact that this holds from a purely structural property of the graph, regardless of kinetics, is what makes it so useful. But in some cases it might be desirable to check complex balance on its own, as for higher deficiency networks.

```@docs
iscomplexbalanced
isdetailedbalanced
```
---
## References
[^1]: [Feinberg, M. *Foundations of Chemical Reaction Network Theory*, Applied Mathematical Sciences 202, Springer (2019).](https://link.springer.com/book/10.1007/978-3-030-03858-8?noAccess=true)
