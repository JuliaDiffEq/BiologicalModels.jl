### Reaction Complex Handling ###

# get the species indices and stoichiometry while filtering out constant species.
function filter_constspecs(specs, stoich::AbstractVector{V}, smap) where {V <: Integer}
    isempty(specs) && (return Vector{Int}(), Vector{V}())

    # if any species are constant, go through these manually and add their indices and 
    # stoichiometries to `ids` and `filtered_stoich`
    if any(isconstant, specs)
        ids = Vector{Int}()
        filtered_stoich = Vector{V}()
        for (i, s) in enumerate(specs)
            if !isconstant(s)
                push!(ids, smap[s])
                push!(filtered_stoich, stoich[i])
            end
        end
    else
        ids = map(Base.Fix1(getindex, smap), specs)
        filtered_stoich = copy(stoich)
    end
    ids, filtered_stoich
end

"""
    reactioncomplexmap(rn::ReactionSystem)

Find each [`ReactionComplex`](@ref) within the specified system, constructing a mapping from
the complex to vectors that indicate which reactions it appears in as substrates and
products.

Notes:
- Each [`ReactionComplex`](@ref) is mapped to a vector of pairs, with each pair having the
  form `reactionidx => ± 1`, where `-1` indicates the complex appears as a substrate and
  `+1` as a product in the reaction with integer label `reactionidx`.
- Constant species are ignored as part of a complex. i.e. if species `A` is constant then
  the reaction `A + B --> C + D` is considered to consist of the complexes `B` and `C + D`.
  Likewise `A --> B` would be treated as the same as `0 --> B`.
"""
function reactioncomplexmap(rn::ReactionSystem)
    isempty(get_systems(rn)) ||
        error("reactioncomplexmap does not currently support subsystems.")

    # check if previously calculated and hence cached
    nps = get_networkproperties(rn)
    !isempty(nps.complextorxsmap) && return nps.complextorxsmap
    complextorxsmap = nps.complextorxsmap

    # retrieves system reactions and a map from species to their index in the species vector
    rxs = reactions(rn)
    smap = speciesmap(rn)
    numreactions(rn) > 0 ||
        error("There must be at least one reaction to find reaction complexes.")

    for (i, rx) in enumerate(rxs)
        # Create the `ReactionComplex` corresponding to the reaction's substrates. Adds it 
        # to the reaction complex dictionary (recording it as the substrates of the i'th reaction).
        subids, substoich = filter_constspecs(rx.substrates, rx.substoich, smap)
        subrc = sort!(ReactionComplex(subids, substoich))
        if haskey(complextorxsmap, subrc)
            push!(complextorxsmap[subrc], i => -1)
        else
            complextorxsmap[subrc] = [i => -1]
        end

        # Create the `ReactionComplex` corresponding to the reaction's products. Adds it 
        # to the reaction complex dictionary (recording it as the products of the i'th reaction).
        prodids, prodstoich = filter_constspecs(rx.products, rx.prodstoich, smap)
        prodrc = sort!(ReactionComplex(prodids, prodstoich))
        if haskey(complextorxsmap, prodrc)
            push!(complextorxsmap[prodrc], i => 1)
        else
            complextorxsmap[prodrc] = [i => 1]
        end
    end
    complextorxsmap
end

@doc raw"""
    reactioncomplexes(network::ReactionSystem; sparse=false)

Calculate the reaction complexes and complex incidence matrix for the given
[`ReactionSystem`](@ref).

Notes:
- returns a pair of a vector of [`ReactionComplex`](@ref)s and the complex incidence matrix.
- An empty [`ReactionComplex`](@ref) denotes the null (∅) state (from reactions like ∅ -> A
  or A -> ∅).
- Constant species are ignored in generating a reaction complex. i.e. if A is constant then
  A --> B consists of the complexes ∅ and B.
- The complex incidence matrix, ``B``, is number of complexes by number of reactions with
```math
B_{i j} = \begin{cases}
-1, &\text{if the i'th complex is the substrate of the j'th reaction},\\
1, &\text{if the i'th complex is the product of the j'th reaction},\\
0, &\text{otherwise.}
\end{cases}
```
- Set sparse=true for a sparse matrix representation of the incidence matrix
"""
function reactioncomplexes(rn::ReactionSystem; sparse = false)
    isempty(get_systems(rn)) ||
        error("reactioncomplexes does not currently support subsystems.")
    nps = get_networkproperties(rn)

    # if the complexes have not been cached, or the cached complexes uses a different sparsity
    if isempty(nps.complexes) || (sparse != issparse(nps.complexes))
        # Computes the reaction complex dictionary. Use it to create a sparse/dense matrix.
        complextorxsmap = reactioncomplexmap(rn)
        nps.complexes, nps.incidencemat = if sparse
            reactioncomplexes(SparseMatrixCSC{Int, Int}, rn, complextorxsmap)
        else
            reactioncomplexes(Matrix{Int}, rn, complextorxsmap)
        end
    end
    nps.complexes, nps.incidencemat
end

# creates a *sparse* reaction complex matrix
function reactioncomplexes(::Type{SparseMatrixCSC{Int, Int}}, rn::ReactionSystem,
        complextorxsmap)
    # computes the I, J, and V vectors used for the sparse matrix (read about sparse matrix 
    # representation for more information)
    complexes = collect(keys(complextorxsmap))
    Is = Int[]
    Js = Int[]
    Vs = Int[]
    for (i, c) in enumerate(complexes)
        for (j, σ) in complextorxsmap[c]
            push!(Is, i)
            push!(Js, j)
            push!(Vs, σ)
        end
    end
    B = sparse(Is, Js, Vs, length(complexes), numreactions(rn))
    complexes, B
end

# creates a *dense* reaction complex matrix
function reactioncomplexes(::Type{Matrix{Int}}, rn::ReactionSystem, complextorxsmap)
    complexes = collect(keys(complextorxsmap))
    B = zeros(Int, length(complexes), numreactions(rn))
    for (i, c) in enumerate(complexes)
        for (j, σ) in complextorxsmap[c]
            B[i, j] = σ
        end
    end
    complexes, B
end

"""
    incidencemat(rn::ReactionSystem; sparse=false)

Calculate the incidence matrix of `rn`, see [`reactioncomplexes`](@ref).

Notes:
- Is cached in `rn` so that future calls, assuming the same sparsity, will also be fast.
"""
incidencemat(rn::ReactionSystem; sparse = false) = reactioncomplexes(rn; sparse)[2]

"""
    complexstoichmat(network::ReactionSystem; sparse=false)

Given a [`ReactionSystem`](@ref) and vector of reaction complexes, return a
matrix with positive entries of size number of species by number of complexes,
where the non-zero positive entries in the kth column denote stoichiometric
coefficients of the species participating in the kth reaction complex.

Notes:
- Set sparse=true for a sparse matrix representation
"""
function complexstoichmat(rn::ReactionSystem; sparse = false)
    isempty(get_systems(rn)) ||
        error("complexstoichmat does not currently support subsystems.")
    nps = get_networkproperties(rn)

    # if the complexes stoichiometry matrix has not been cached, or the cached one uses a
    # different sparsity, computes (and caches) it
    if isempty(nps.complexstoichmat) || (sparse != issparse(nps.complexstoichmat))
        nps.complexstoichmat = if sparse
            complexstoichmat(SparseMatrixCSC{Int, Int}, rn, keys(reactioncomplexmap(rn)))
        else
            complexstoichmat(Matrix{Int}, rn, keys(reactioncomplexmap(rn)))
        end
    end
    nps.complexstoichmat
end

# creates a *sparse* reaction complex stoichiometry matrix
function complexstoichmat(::Type{SparseMatrixCSC{Int, Int}}, rn::ReactionSystem, rcs)
    # computes the I, J, and V vectors used for the sparse matrix (read about sparse matrix 
    # representation for more information)
    Is = Int[]
    Js = Int[]
    Vs = Int[]
    for (i, rc) in enumerate(rcs)
        for rcel in rc
            push!(Is, rcel.speciesid)
            push!(Js, i)
            push!(Vs, rcel.speciesstoich)
        end
    end
    Z = sparse(Is, Js, Vs, numspecies(rn), length(rcs))
end

# creates a *dense* reaction complex stoichiometry matrix
function complexstoichmat(::Type{Matrix{Int}}, rn::ReactionSystem, rcs)
    Z = zeros(Int, numspecies(rn), length(rcs))
    for (i, rc) in enumerate(rcs)
        for rcel in rc
            Z[rcel.speciesid, i] = rcel.speciesstoich
        end
    end
    Z
end

@doc raw"""
    complexoutgoingmat(network::ReactionSystem; sparse=false)

Given a [`ReactionSystem`](@ref) and complex incidence matrix, ``B``, return a
matrix of size num of complexes by num of reactions that identifies substrate
complexes.

Notes:
- The complex outgoing matrix, ``\Delta``, is defined by
```math
\Delta_{i j} = \begin{cases}
    = 0,    &\text{if } B_{i j} = 1, \\
    = B_{i j}, &\text{otherwise.}
\end{cases}
```
- Set sparse=true for a sparse matrix representation
"""
function complexoutgoingmat(rn::ReactionSystem; sparse = false)
    isempty(get_systems(rn)) ||
        error("complexoutgoingmat does not currently support subsystems.")
    nps = get_networkproperties(rn)

    # if the outgoing complexes matrix has not been cached, or the cached one uses a
    # different sparsity, computes (and caches) it
    if isempty(nps.complexoutgoingmat) || (sparse != issparse(nps.complexoutgoingmat))
        B = reactioncomplexes(rn, sparse = sparse)[2]
        nps.complexoutgoingmat = if sparse
            complexoutgoingmat(SparseMatrixCSC{Int, Int}, rn, B)
        else
            complexoutgoingmat(Matrix{Int}, rn, B)
        end
    end
    nps.complexoutgoingmat
end

# creates a *sparse* outgoing reaction complex stoichiometry matrix
function complexoutgoingmat(::Type{SparseMatrixCSC{Int, Int}}, rn::ReactionSystem, B)
    # computes the I, J, and V vectors used for the sparse matrix (read about sparse matrix 
    # representation for more information)
    n = size(B, 2)
    rows = rowvals(B)
    vals = nonzeros(B)
    Is = Int[]
    Js = Int[]
    Vs = Int[]

    # allocates space to the vectors (so that it is not done incrementally in the loop)
    sizehint!(Is, div(length(vals), 2))
    sizehint!(Js, div(length(vals), 2))
    sizehint!(Vs, div(length(vals), 2))

    for j in 1:n
        for i in nzrange(B, j)
            if vals[i] != one(eltype(vals))
                push!(Is, rows[i])
                push!(Js, j)
                push!(Vs, vals[i])
            end
        end
    end
    sparse(Is, Js, Vs, size(B, 1), size(B, 2))
end

# creates a *dense* outgoing reaction complex stoichiometry matrix
function complexoutgoingmat(::Type{Matrix{Int}}, rn::ReactionSystem, B)
    Δ = copy(B)
    for (I, b) in pairs(Δ)
        (b == 1) && (Δ[I] = 0)
    end
    Δ
end

"""
    incidencematgraph(rn::ReactionSystem)

Construct a directed simple graph where nodes correspond to reaction complexes and directed
edges to reactions converting between two complexes.

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
incidencematgraph(sir)
```
"""
function incidencematgraph(rn::ReactionSystem)
    nps = get_networkproperties(rn)
    if Graphs.nv(nps.incidencegraph) == 0
        isempty(nps.incidencemat) && reactioncomplexes(rn)
        nps.incidencegraph = incidencematgraph(nps.incidencemat)
    end
    nps.incidencegraph
end

# computes the incidence graph from an *dense* incidence matrix
function incidencematgraph(incidencemat::Matrix{Int})
    @assert all(∈([-1, 0, 1]), incidencemat)
    n = size(incidencemat, 1)  # no. of nodes/complexes
    graph = Graphs.DiGraph(n)

    # Walks through each column (corresponds to reactions). For each, find the input and output
    # complex and add an edge representing these to the incidence graph.
    for col in eachcol(incidencemat)
        src = 0
        dst = 0
        for i in eachindex(col)
            (col[i] == -1) && (src = i)
            (col[i] == 1) && (dst = i)
            (src != 0) && (dst != 0) && break
        end
        Graphs.add_edge!(graph, src, dst)
    end
    return graph
end

# computes the incidence graph from an *sparse* incidence matrix
function incidencematgraph(incidencemat::SparseMatrixCSC{Int, Int})
    @assert all(∈([-1, 0, 1]), incidencemat)
    m, n = size(incidencemat)
    graph = Graphs.DiGraph(m)
    rows = rowvals(incidencemat)
    vals = nonzeros(incidencemat)

    # Loops through the (n) columns. For each column, directly find the index of the input
    # and output complex and add an edge representing these to the incidence graph.
    for j in 1:n
        inds = nzrange(incidencemat, j)
        row = rows[inds]
        val = vals[inds]
        if val[1] == -1
            Graphs.add_edge!(graph, row[1], row[2])
        else
            Graphs.add_edge!(graph, row[2], row[1])
        end
    end
    return graph
end

### Linkage, Deficiency, Reversibility ###

"""
    linkageclasses(rn::ReactionSystem)

Given the incidence graph of a reaction network, return a vector of the
connected components of the graph (i.e. sub-groups of reaction complexes that
are connected in the incidence graph).

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
linkageclasses(sir)
```
gives
```julia
2-element Vector{Vector{Int64}}:
 [1, 2]
 [3, 4]
```
"""
function linkageclasses(rn::ReactionSystem)
    nps = get_networkproperties(rn)
    if isempty(nps.linkageclasses)
        nps.linkageclasses = linkageclasses(incidencematgraph(rn))
    end
    nps.linkageclasses
end

linkageclasses(incidencegraph) = Graphs.connected_components(incidencegraph)

"""
    stronglinkageclasses(rn::ReactionSystem)

    Return the strongly connected components of a reaction network's incidence graph (i.e. sub-groups of reaction complexes such that every complex is reachable from every other one in the sub-group).
"""

function stronglinkageclasses(rn::ReactionSystem)
    nps = get_networkproperties(rn)
    if isempty(nps.stronglinkageclasses)
        nps.stronglinkageclasses = stronglinkageclasses(incidencematgraph(rn))
    end
    nps.stronglinkageclasses
end

stronglinkageclasses(incidencegraph) = Graphs.strongly_connected_components(incidencegraph)

"""
    terminallinkageclasses(rn::ReactionSystem)

    Return the terminal strongly connected components of a reaction network's incidence graph (i.e. sub-groups of reaction complexes that are 1) strongly connected and 2) every outgoing reaction from a complex in the component produces a complex also in the component).
"""

function terminallinkageclasses(rn::ReactionSystem)
    nps = get_networkproperties(rn)
    if isempty(nps.terminallinkageclasses)
        slcs = stronglinkageclasses(rn)
        tslcs = filter(lc -> isterminal(lc, rn), slcs)
        nps.terminallinkageclasses = tslcs
    end
    nps.terminallinkageclasses
end

# Helper function for terminallinkageclasses. Given a linkage class and a reaction network, say 
# whether the linkage class is terminal, i.e. all outgoing reactions from complexes in the linkage 
# class produce a complex also in the linkage class
function isterminal(lc::Vector, rn::ReactionSystem)
    imat = incidencemat(rn)

    for r in 1:size(imat, 2)
        # Find the index of the reactant complex for a given reaction
        s = findfirst(==(-1), @view imat[:, r])

        # If the reactant complex is in the linkage class, check whether the product complex is
        # also in the linkage class. If any of them are not, return false. 
        if s in Set(lc)
            p = findfirst(==(1), @view imat[:, r])
            p in Set(lc) ? continue : return false
        end
    end
    true
end

@doc raw"""
    deficiency(rn::ReactionSystem)

Calculate the deficiency of a reaction network.

Here the deficiency, ``\delta``, of a network with ``n`` reaction complexes,
``\ell`` linkage classes and a rank ``s`` stoichiometric matrix is

```math
\delta = n - \ell - s
```

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
δ = deficiency(sir)
```
"""
function deficiency(rn::ReactionSystem)
    # Precomputes required information. `conservationlaws` caches the conservation laws in `rn`.
    nps = get_networkproperties(rn)
    conservationlaws(rn)
    r = nps.rank
    ig = incidencematgraph(rn)
    lc = linkageclasses(rn)

    # Computes deficiency using its formula. Caches and returns it as output.
    nps.deficiency = Graphs.nv(ig) - length(lc) - r
    nps.deficiency
end

# For a linkage class (set of reaction complexes that form an isolated sub-graph), and some
# additional information of the full network, find the reactions, species, and parameters 
# that constitute the corresponding sub-reaction network.
function subnetworkmapping(linkageclass, allrxs, complextorxsmap, p)
    # Finds the reactions that are part of teh sub-reaction network.
    rxinds = sort!(collect(Set(
        rxidx for rcidx in linkageclass for rxidx in complextorxsmap[rcidx])))
    newrxs = allrxs[rxinds]
    specset = Set(s for rx in newrxs for s in rx.substrates if !isconstant(s))
    for rx in newrxs
        for product in rx.products
            !isconstant(product) && push!(specset, product)
        end
    end
    newspecs = collect(specset)

    # Find the parameters that are part of the sub-reaction network.
    newps = Vector{eltype(p)}()
    for rx in newrxs
        Symbolics.get_variables!(newps, rx.rate, p)
    end

    newrxs, newspecs, newps
end

"""
    subnetworks(rn::ReactionSystem)

Find subnetworks corresponding to each linkage class of the reaction network.

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
subnetworks(sir)
```
"""
function subnetworks(rs::ReactionSystem)
    isempty(get_systems(rs)) || error("subnetworks does not currently support subsystems.")

    # Retrieves required components. `linkageclasses` caches linkage classes in `rs`.
    lcs = linkageclasses(rs)
    rxs = reactions(rs)
    p = parameters(rs)
    t = get_iv(rs)
    spatial_ivs = get_sivs(rs)
    complextorxsmap = [map(first, rcmap) for rcmap in values(reactioncomplexmap(rs))]
    subnetworks = Vector{ReactionSystem}()

    # Loops through each sub-graph of connected reaction complexes. For each, create a 
    # new `ReactionSystem` model and pushes it to the subnetworks vector.
    for i in 1:length(lcs)
        newrxs, newspecs, newps = subnetworkmapping(lcs[i], rxs, complextorxsmap, p)
        newname = Symbol(nameof(rs), "_", i)
        push!(subnetworks,
            ReactionSystem(newrxs, t, newspecs, newps; name = newname, spatial_ivs))
    end
    subnetworks
end

"""
    linkagedeficiencies(network::ReactionSystem)

Calculates the deficiency of each sub-reaction network within `network`.

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
linkage_deficiencies = linkagedeficiencies(sir)
```
"""
function linkagedeficiencies(rs::ReactionSystem)
    lcs = linkageclasses(rs)
    subnets = subnetworks(rs)
    δ = zeros(Int, length(lcs))

    # For each sub-reaction network of the reaction network, compute its deficiency. Returns
    # the full vector of deficiencies for each sub-reaction network.
    for (i, subnet) in enumerate(subnets)
        conservationlaws(subnet)
        nps = get_networkproperties(subnet)
        δ[i] = length(lcs[i]) - 1 - nps.rank
    end
    δ
end

"""
    isreversible(rn::ReactionSystem)

Given a reaction network, returns if the network is reversible or not.

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
isreversible(sir)
```
"""
function isreversible(rn::ReactionSystem)
    ig = incidencematgraph(rn)
    Graphs.reverse(ig) == ig
end

"""
    isweaklyreversible(rn::ReactionSystem, subnetworks)

Determine if the reaction network with the given subnetworks is weakly reversible or not.

For example,
```julia
sir = @reaction_network SIR begin
    β, S + I --> 2I
    ν, I --> R
end
subnets = subnetworks(rn)
isweaklyreversible(rn, subnets)
```
"""
function isweaklyreversible(rn::ReactionSystem, subnets)
    nps = get_networkproperties(rn)
    isempty(nps.incidencemat) && reactioncomplexes(rn)
    sparseig = issparse(nps.incidencemat)

    for subnet in subnets
        subnps = get_networkproperties(subnet)
        isempty(subnps.incidencemat) && reactioncomplexes(subnet; sparse = sparseig)
    end

    # A network is weakly reversible if all of its subnetworks are strongly connected
    all(Graphs.is_strongly_connected ∘ incidencematgraph, subnets)
end

### Conservation Laws ###

# Implements the `conserved` parameter metadata.
struct ConservedParameter end
Symbolics.option_to_metadata_type(::Val{:conserved}) = ConservedParameter

"""
isconserved(p)

Checks if the input parameter (`p`) is a conserved quantity (i.e. have the `conserved`)
metadata.
"""
isconserved(x::Num, args...) = isconserved(Symbolics.unwrap(x), args...)
function isconserved(x, default = false)
    p = Symbolics.getparent(x, nothing)
    p === nothing || (x = p)
    Symbolics.getmetadata(x, ConservedParameter, default)
end

"""
    conservedequations(rn::ReactionSystem)

Calculate symbolic equations from conservation laws, writing dependent variables as
functions of independent variables and the conservation law constants.

Notes:
- Caches the resulting equations in `rn`, so will be fast on subsequent calls.

Examples:
```@repl
rn = @reaction_network begin
    k, A + B --> C
    k2, C --> A + B
    end
conservedequations(rn)
```
gives
```
2-element Vector{Equation}:
 B(t) ~ A(t) + Γ[1]
 C(t) ~ Γ[2] - A(t)
```
"""
function conservedequations(rn::ReactionSystem)
    conservationlaws(rn)
    nps = get_networkproperties(rn)
    nps.conservedeqs
end

"""
    conservationlaw_constants(rn::ReactionSystem)

Calculate symbolic equations from conservation laws, writing the conservation law constants
in terms of the dependent and independent variables.

Notes:
- Caches the resulting equations in `rn`, so will be fast on subsequent calls.

Examples:
```@julia
rn = @reaction_network begin
    k, A + B --> C
    k2, C --> A + B
    end
conservationlaw_constants(rn)
```
gives
```
2-element Vector{Equation}:
 Γ[1] ~ B(t) - A(t)
 Γ[2] ~ A(t) + C(t)
```
"""
function conservationlaw_constants(rn::ReactionSystem)
    conservationlaws(rn)
    nps = get_networkproperties(rn)
    nps.constantdefs
end

"""
    conservationlaws(netstoichmat::AbstractMatrix)::Matrix

Given the net stoichiometry matrix of a reaction system, computes a matrix of
conservation laws, each represented as a row in the output.
"""
function conservationlaws(nsm::T; col_order = nothing) where {T <: AbstractMatrix}
    # compute the left nullspace over the integers
    # the `nullspace` function updates the `col_order`
    N = MT.nullspace(nsm'; col_order)

    # if all coefficients for a conservation law are negative, make positive
    for Nrow in eachcol(N)
        all(r -> r <= 0, Nrow) && (Nrow .*= -1)
    end

    # check we haven't overflowed
    iszero(N' * nsm) || error("Calculation of the conservation law matrix was inaccurate, "
          * "likely due to numerical overflow. Please use a larger integer "
          * "type like Int128 or BigInt for the net stoichiometry matrix.")

    T(N')
end

# used in the subsequent function
function cache_conservationlaw_eqs!(rn::ReactionSystem, N::AbstractMatrix, col_order)
    # Retrieves nullity (the number of conservation laws). `r` is the rank of the netstoichmat.
    nullity = size(N, 1)
    r = numspecies(rn) - nullity

    # Creates vectors of all independent and dependent species (those that are not, and are, 
    # eliminated by the conservation laws). Get vectors both with their indexes and the actual
    # species symbolic variables.
    sps = species(rn)
    indepidxs = col_order[begin:r]
    indepspecs = sps[indepidxs]
    depidxs = col_order[(r + 1):end]
    depspecs = sps[depidxs]

    # declares the conservation law parameters
    constants = MT.unwrap.(MT.scalarize(only(
        @parameters $(CONSERVED_CONSTANT_SYMBOL)[1:nullity] [conserved = true])))

    # Computes the equations for (examples uses simple two-state system, `X1 <--> X2`):
    # - The species eliminated through conservation laws (`conservedeqs`). E.g. `[X2 ~ Γ[1] - X1]`.
    # - The conserved quantity parameters (`constantdefs`). E.g. `[Γ[1] ~ X1 + X2]`.
    conservedeqs = Equation[]
    constantdefs = Equation[]

    # for each conserved quantity
    for (i, depidx) in enumerate(depidxs)
        # finds the coefficient (in the conservation law) of the species that is eliminated
        # by this conservation law
        scaleby = (N[i, depidx] != 1) ? N[i, depidx] : one(eltype(N))
        (scaleby != 0) ||
            error("Error, found a zero in the conservation law matrix where one was not expected.")

        # creates, for this conservation law, the sum of all independent species (weighted by
        # the ratio between the coefficient of the species and the species which is elimianted
        coefs = @view N[i, indepidxs]
        terms = sum(coef / scaleby * sp for (coef, sp) in zip(coefs, indepspecs))

        # computes the two equations corresponding to this conserved quantity
        eq = depspecs[i] ~ constants[i] - terms
        push!(conservedeqs, eq)
        eq = constants[i] ~ depspecs[i] + terms
        push!(constantdefs, eq)
    end

    # cache in the system
    nps = get_networkproperties(rn)
    nps.rank = r
    nps.nullity = nullity
    nps.indepspecs = Set(indepspecs)
    nps.depspecs = Set(depspecs)
    nps.conservedeqs = conservedeqs
    nps.constantdefs = constantdefs

    nothing
end

"""
    conservationlaws(rs::ReactionSystem)

Return the conservation law matrix of the given `ReactionSystem`, calculating it if it is
not already stored within the system, or returning an alias to it.

Notes:
- The first time being called it is calculated and cached in `rn`, subsequent calls should
  be fast.
"""
function conservationlaws(rs::ReactionSystem)
    nps = get_networkproperties(rs)
    !isempty(nps.conservationmat) && (return nps.conservationmat)

    # if the conservation law matrix is not computed, do so and caches the result
    nsm = netstoichmat(rs)
    nps.conservationmat = conservationlaws(nsm; col_order = nps.col_order)
    cache_conservationlaw_eqs!(rs, nps.conservationmat, nps.col_order)
    nps.conservationmat
end

"""
    conservedquantities(state, cons_laws)

Compute conserved quantities for a system with the given conservation laws.
"""
conservedquantities(state, cons_laws) = cons_laws * state

# If u0s are not given while conservation laws are present, throw an error.
# Used in HomotopyContinuation and BifurcationKit extensions.
# Currently only checks if any u0s are given
# (not whether these are enough for computing conserved quantitites, this will yield a less informative error).
function conservationlaw_errorcheck(rs, pre_varmap)
    vars_with_vals = Set(p[1] for p in pre_varmap)
    any(sp -> sp in vars_with_vals, species(rs)) && return
    isempty(conservedequations(Catalyst.flatten(rs))) ||
        error("The system has conservation laws but initial conditions were not provided for some species.")
end

"""
    iscomplexbalanced(rs::ReactionSystem, parametermap)

Constructively compute whether a network will have complex-balanced equilibrium
solutions, following the method in van der Schaft et al., [2015](https://link.springer.com/article/10.1007/s10910-015-0498-2#Sec3). Accepts a dictionary, vector, or tuple of variable-to-value mappings, e.g. [k1 => 1.0, k2 => 2.0,...]. 
"""

function iscomplexbalanced(rs::ReactionSystem, parametermap::Dict)
    if length(parametermap) != numparams(rs)
        error("Incorrect number of parameters specified.")
    end

    pmap = symmap_to_varmap(rs, parametermap)
    pmap = Dict(ModelingToolkit.value(k) => v for (k, v) in pmap)

    sm = speciesmap(rs)
    cm = reactioncomplexmap(rs)
    complexes, D = reactioncomplexes(rs)
    rxns = reactions(rs)
    nc = length(complexes)
    nr = numreactions(rs)
    nm = numspecies(rs)

    if !all(r -> ismassaction(r, rs), rxns)
        error("The supplied ReactionSystem has reactions that are not ismassaction. Testing for being complex balanced is currently only supported for pure mass action networks.")
    end

    rates = [substitute(rate, pmap) for rate in reactionrates(rs)]

    # Construct kinetic matrix, K
    K = zeros(nr, nc)
    for c in 1:nc
        complex = complexes[c]
        for (r, dir) in cm[complex]
            rxn = rxns[r]
            if dir == -1
                K[r, c] = rates[r]
            end
        end
    end

    L = -D * K
    S = netstoichmat(rs)

    # Compute ρ using the matrix-tree theorem
    g = incidencematgraph(rs)
    R = ratematrix(rs, rates)
    ρ = matrixtree(g, R)

    # Determine if 1) ρ is positive and 2) D^T Ln ρ lies in the image of S^T
    if all(>(0), ρ)
        img = D' * log.(ρ)
        if rank(S') == rank(hcat(S', img))
            return true
        else
            return false
        end
    else
        return false
    end
end

function iscomplexbalanced(rs::ReactionSystem, parametermap::Vector{Pair{Symbol, Float64}})
    pdict = Dict(parametermap)
    iscomplexbalanced(rs, pdict)
end

function iscomplexbalanced(rs::ReactionSystem, parametermap::Tuple{Pair{Symbol, Float64}})
    pdict = Dict(parametermap)
    iscomplexbalanced(rs, pdict)
end

function iscomplexbalanced(rs::ReactionSystem, parametermap)
    error("Parameter map must be a dictionary, tuple, or vector of symbol/value pairs.")
end

"""
    ratematrix(rs::ReactionSystem, parametermap)

    Given a reaction system with n complexes, outputs an n-by-n matrix where R_{ij} is the rate constant of the reaction between complex i and complex j. Accepts a dictionary, vector, or tuple of variable-to-value mappings, e.g. [k1 => 1.0, k2 => 2.0,...]. 
"""

function ratematrix(rs::ReactionSystem, rates::Vector{Float64})
    complexes, D = reactioncomplexes(rs)
    n = length(complexes)
    rxns = reactions(rs)
    ratematrix = zeros(n, n)

    for r in 1:length(rxns)
        rxn = rxns[r]
        s = findfirst(==(-1), @view D[:, r])
        p = findfirst(==(1), @view D[:, r])
        ratematrix[s, p] = rates[r]
    end
    ratematrix
end

function ratematrix(rs::ReactionSystem, parametermap::Dict)
    if length(parametermap) != numparams(rs)
        error("Incorrect number of parameters specified.")
    end

    pmap = symmap_to_varmap(rs, parametermap)
    pmap = Dict(ModelingToolkit.value(k) => v for (k, v) in pmap)

    rates = [substitute(rate, pmap) for rate in reactionrates(rs)]
    ratematrix(rs, rates)
end

function ratematrix(rs::ReactionSystem, parametermap::Vector{Pair{Symbol, Float64}})
    pdict = Dict(parametermap)
    ratematrix(rs, pdict)
end

function ratematrix(rs::ReactionSystem, parametermap::Tuple{Pair{Symbol, Float64}})
    pdict = Dict(parametermap)
    ratematrix(rs, pdict)
end

function ratematrix(rs::ReactionSystem, parametermap)
    error("Parameter map must be a dictionary, tuple, or vector of symbol/value pairs.")
end

### BELOW: Helper functions for iscomplexbalanced

function matrixtree(g::SimpleDiGraph, distmx::Matrix)
    n = nv(g)
    if size(distmx) != (n, n)
        error("Size of distance matrix is incorrect")
    end

    π = zeros(n)

    if !Graphs.is_connected(g)
        ccs = Graphs.connected_components(g)
        for cc in ccs
            sg, vmap = Graphs.induced_subgraph(g, cc)
            distmx_s = distmx[cc, cc]
            π_j = matrixtree(sg, distmx_s)
            π[cc] = π_j
        end
        return π
    end

    # generate all spanning trees
    ug = SimpleGraph(SimpleDiGraph(g))
    trees = collect(Combinatorics.combinations(collect(edges(ug)), n - 1))
    trees = SimpleGraph.(trees)
    trees = filter!(t -> isempty(Graphs.cycle_basis(t)), trees)

    # constructed rooted trees for every vertex, compute sum
    for v in 1:n
        rootedTrees = [reverse(Graphs.bfs_tree(t, v, dir = :in)) for t in trees]
        π[v] = sum([treeweight(t, g, distmx) for t in rootedTrees])
    end

    # sum the contributions
    return π
end

function treeweight(t::SimpleDiGraph, g::SimpleDiGraph, distmx::Matrix)
    prod = 1
    for e in edges(t)
        s = Graphs.src(e)
        t = Graphs.dst(e)
        prod *= distmx[s, t]
    end
    prod
end
