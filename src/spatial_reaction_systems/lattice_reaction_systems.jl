### Lattice Reaction Network Structure ###

# Describes a spatial reaction network over a lattice.
# Adding the "<: MT.AbstractTimeDependentSystem" part messes up show, disabling me from creating LRSs. Should be fixed some time.
struct LatticeReactionSystem{Q,R,S,T} # <: MT.AbstractTimeDependentSystem 
    # Input values.
    """The (non-spatial) reaction system within each vertexes."""
    rs::ReactionSystem{Q}
    """The spatial reactions defined between individual vertexes."""
    spatial_reactions::Vector{R}
    """The lattice on which the (discrete) spatial system is defined."""
    lattice::S

    # Derived values.
    """The number of vertexes (compartments)."""
    num_verts::Int64
    """The number of edges."""
    num_edges::Int64
    """The number of species."""
    num_species::Int64

    """List of species that may move spatially."""
    spat_species::Vector{BasicSymbolic{Real}}
    """
    All parameters related to the lattice reaction system
    (both those whose values are tied to vertexes and edges).
    """
    parameters::Vector{BasicSymbolic{Real}}
    """
    Parameters which values are tied to vertexes, 
    e.g. that possibly could have unique values at each vertex of the system.
    """
    vertex_parameters::Vector{BasicSymbolic{Real}}
    """
    Parameters whose values are tied to edges (adjacencies), 
    e.g. that possibly could have unique values at each edge of the system.
    """
    edge_parameters::Vector{BasicSymbolic{Real}}
    """
    An iterator over all the lattice's edges. Currently, the format is always a Vector{Pair{Int64,Int64}}.
    However, in the future, different types could potentially be used for different types of lattice
    (E.g. for a Cartesian grid, we do not technically need to enumerate each edge)
    """
    edge_iterator::T

    function LatticeReactionSystem(rs::ReactionSystem{Q}, spatial_reactions::Vector{R}, lattice::S, 
                                   num_verts::Int64, num_edges::Int64, edge_iterator::T) where {Q,R,S,T}
        # Error checks.
        if !(R <: AbstractSpatialReaction)
            error("The second argument must be a vector of AbstractSpatialReaction subtypes.") 
        end

        # Computes the species which are parts of spatial reactions. Also counts total number of 
        # species types.
        if isempty(spatial_reactions)
            spat_species = Vector{BasicSymbolic{Real}}[]
        else
            spat_species = unique(reduce(vcat, [spatial_species(sr) for sr in spatial_reactions]))
        end
        num_species = length(unique([species(rs); spat_species]))

        # Computes the sets of vertex, edge, and all, parameters.
        rs_edge_parameters = filter(isedgeparameter, parameters(rs))
        if isempty(spatial_reactions)
            srs_edge_parameters = Vector{BasicSymbolic{Real}}[]
        else
            srs_edge_parameters = setdiff(reduce(vcat, [parameters(sr) for sr in spatial_reactions]), parameters(rs))
        end
        edge_parameters = unique([rs_edge_parameters; srs_edge_parameters])
        vertex_parameters = filter(!isedgeparameter, parameters(rs))

        # Ensures the parameter order begins similarly to in the non-spatial ReactionSystem.
        ps = [parameters(rs); setdiff([edge_parameters; vertex_parameters], parameters(rs))]    

        # Checks that all spatial reactions are valid for this reactions system.
        foreach(sr -> check_spatial_reaction_validity(rs, sr; edge_parameters=edge_parameters), spatial_reactions)   

        return new{Q,R,S,T}(rs, spatial_reactions, lattice, num_verts, num_edges, num_species, 
                            spat_species, ps, vertex_parameters, edge_parameters, edge_iterator)
    end
end

# Creates a LatticeReactionSystem from a CartesianGrid lattice (cartesian grid).
function LatticeReactionSystem(rs, srs, lattice_in::CartesianGridRej{S,T}; diagonal_connections=false) where {S,T}
    # Error checks.
    (length(lattice_in.dims) > 3) && error("Grids of higher dimension than 3 is currently not supported.")

    # Ensures that the matrix has a 3d form (used for intermediary computations only,
    # the original is passed to the constructor).
    lattice = CartesianGrid((lattice_in.dims..., fill(1, 3-length(lattice_in.dims))...))

    # Counts vertexes and edges. The `num_edges` count formula counts the number of internal, side, 
    # edge, and corner vertexes (on the grid). The number of edges from each depends on whether diagonal
    # connections are allowed. The formula holds even if l, m, and/or n are 1.
    l,m,n = lattice.dims
    num_verts = l * m * n
    (ni, ns, ne, nc) = diagonal_connections ? (26,17,11,7) : (6,5,4,3)
    num_edges = ni*(l-2)*(m-2)*(n-2) +                            # Edges from internal vertexes.
                ns*(2(l-2)*(m-2) + 2(l-2)*(n-2) + 2(m-2)*(n-2)) + # Edges from side vertexes.
                ne*(4(l-2) + 4(m-2) + 4(n-2)) +                   # Edges from edge vertexes.
                nc*8                                              # Edges from corner vertexes.
    
    # Creates an iterator over all edges.
    # Currently creates a full vector. Future version might be for efficient.
    edge_iterator = Vector{Pair{Int64,Int64}}(undef, num_edges)
    # Loops through, simultaneously, the coordinates of each position in the grid, as well as that 
    # coordinate's (scalar) flat index. For each grid point, loops through all potential neighbours.  
    indices = [(L, M, N) for L in 1:l, M in 1:m, N in 1:n]
    flat_indices = 1:num_verts
    next_vert = 0
    for ((L, M, N), idx) in zip(indices, flat_indices)
        for LL in max(L - 1, 1):min(L + 1, l), 
            MM in max(M - 1, 1):min(M + 1, m), 
            NN in max(N - 1, 1):min(N + 1, n)
    
            # Which (LL,MM,NN) indexes are valid neighbours depends on whether diagonal connects are permitted.
            !diagonal_connections && (count([L==LL, M==MM, N==NN]) == 2) || continue
            diagonal_connections && (L==LL) && (M==MM) && (N==NN) && continue
            
            # Computes the neighbour's flat (scalar) index. Add the edge to edge_iterator.
            neighbour_idx = LL + (MM - 1) * l + (NN - 1) * m * l
            edge_iterator[next_vert += 1] = (idx => neighbour_idx)
        end
    end

    return LatticeReactionSystem(rs, srs, lattice_in, num_verts, num_edges, edge_iterator)
end

# Creates a LatticeReactionSystem from a Boolean Array lattice (masked grid).
function LatticeReactionSystem(rs, srs, lattice_in::Array{Bool, T}; diagonal_connections=false) where {T}  
    # Error checks.
    dims = size(lattice_in)
    (length(dims) > 3) && error("Grids of higher dimension than 3 is currently not supported.")

    # Ensures that the matrix has a 3d form (used for intermediary computations only,
    # the original is passed to the constructor).
    lattice = reshape(lattice_in, [dims...; fill(1, 3-length(dims))]...)
    
    # Counts vertexes (edges have to be counted after the iterator have been created).
    num_verts = count(lattice)


    # Makes a template matrix to store each vertex's index. The matrix is 0 where there is no vertex.
    # The template is used in the next step.
    idx_matrix = fill(0, size(lattice_in))
    cur_vertex_idx = 0
    for flat_idx in 1:length(lattice)
        if lattice[flat_idx]
            idx_matrix[flat_idx] = (cur_vertex_idx += 1)
        end
    end

    # Creates an iterator over all edges. A vector with pairs of each edge's source to its destination.
    edge_iterator = Vector{Pair{Int64,Int64}}()
    # Loops through, the coordinates of each position in the grid. 
    # For each grid point, loops through all potential neighbours and adds edges to edge_iterator.
    l, m, n = size(lattice)
    indices = [(L, M, N) for L in 1:l, M in 1:m, N in 1:n]
    for (L, M, N) in indices
        # Ensures that we are in a valid lattice point.
        lattice[L,M,N] || continue
        for LL in max(L - 1, 1):min(L + 1, l), 
            MM in max(M - 1, 1):min(M + 1, m), 
            NN in max(N - 1, 1):min(N + 1, n)

            # Ensures that the neighbour is a valid lattice point.
            lattice[LL,MM,NN] || continue

            # Which (LL,MM,NN) indexes are valid neighbours depends on whether diagonal connects are permitted.
            !diagonal_connections && (count([L==LL, M==MM, N==NN]) == 2) || continue
            diagonal_connections && (L==LL) && (M==MM) && (N==NN) && continue
            
            # Computes the neighbour's scalar index. Add that connection to `edge_iterator`.
            push!(edge_iterator, idx_matrix[L,M,N] => idx_matrix[LL,MM,NN])
        end
    end
    num_edges = length(edge_iterator)

    return LatticeReactionSystem(rs, srs, lattice_in, num_verts, num_edges, edge_iterator)
end

# Creates a LatticeReactionSystem from a (directed) Graph lattice (graph grid).
function LatticeReactionSystem(rs, srs, lattice::DiGraph)
    num_verts = nv(lattice)
    num_edges = ne(lattice)
    edge_iterator = [e.src => e.dst for e in edges(lattice)]
    return LatticeReactionSystem(rs, srs, lattice, num_verts, num_edges, edge_iterator)
end
# Creates a LatticeReactionSystem from a (undirected) Graph lattice (graph grid).
LatticeReactionSystem(rs, srs, lattice::SimpleGraph) = LatticeReactionSystem(rs, srs, DiGraph(lattice))


### Lattice ReactionSystem Getters ###

# Get all species.
species(lrs::LatticeReactionSystem) = unique([species(lrs.rs); lrs.spat_species])
# Get all species that may be transported.
spatial_species(lrs::LatticeReactionSystem) = lrs.spat_species

# Get all parameters.
ModelingToolkit.parameters(lrs::LatticeReactionSystem) = lrs.parameters
# Get all parameters whose values are tied to vertexes (compartments).
vertex_parameters(lrs::LatticeReactionSystem) = lrs.vertex_parameters
# Get all parameters whose values are tied to edges (adjacencies).
edge_parameters(lrs::LatticeReactionSystem) = lrs.edge_parameters

# Gets the lrs name (same as rs name).
ModelingToolkit.nameof(lrs::LatticeReactionSystem) = nameof(lrs.rs)

# Checks if a lattice reaction system is a pure (linear) transport reaction system.
is_transport_system(lrs::LatticeReactionSystem) = all(sr -> sr isa TransportReaction, lrs.spatial_reactions)

"""
    has_cartesian_lattice(lrs::LatticeReactionSystem)

Returns `true` if `lrs` was created using a cartesian grid lattice (e.g. created via `CartesianGrid(5,5)`). 
Otherwise, returns `false`.
"""
has_cartesian_lattice(lrs::LatticeReactionSystem) = lrs.lattice isa CartesianGridRej{S,T} where {S,T}

"""
    has_masked_lattice(lrs::LatticeReactionSystem)

Returns `true` if `lrs` was created using a masked grid lattice (e.g. created via `[true true; true false]`). 
Otherwise, returns `false`.
"""
has_masked_lattice(lrs::LatticeReactionSystem) = lrs.lattice isa Array{Bool, T} where T

"""
    has_grid_lattice(lrs::LatticeReactionSystem)

Returns `true` if `lrs` was created using a cartesian or masked grid lattice. Otherwise, returns `false`.
"""
has_grid_lattice(lrs::LatticeReactionSystem) = (has_cartesian_lattice(lrs) || has_masked_lattice(lrs))

"""
    has_graph_lattice(lrs::LatticeReactionSystem)

Returns `true` if `lrs` was created using a graph grid lattice (e.g. created via `path_graph(5)`). 
Otherwise, returns `false`.
"""
has_graph_lattice(lrs::LatticeReactionSystem) = lrs.lattice isa SimpleDiGraph

"""
    grid_size(lrs::LatticeReactionSystem)

Returns the size of `lrs`'s lattice (only if it is a cartesian or masked grid lattice). 
E.g. for a lattice `CartesianGrid(4,6)`, `(4,6)` is returned.
"""
function grid_size(lrs::LatticeReactionSystem)
    has_cartesian_lattice(lrs) && (return lrs.lattice.dims)
    has_masked_lattice(lrs) && (return size(lrs.lattice))
    error("Grid size is only defined for LatticeReactionSystems with grid-based lattices (not graph-based).")
end

"""
    grid_dims(lrs::LatticeReactionSystem)

Returns the number of dimensions of `lrs`'s lattice (only if it is a cartesian or masked grid lattice). 
The output is either `1`, `2`, or `3`.
"""
grid_dims(lrs::LatticeReactionSystem) = length(grid_size(lrs))