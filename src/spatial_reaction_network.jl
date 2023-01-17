### Spatial Reaction Structure. ###
# Describing a spatial reaction that involves species from two neighbouring compartments.
# Currently only permit constant rate.
struct SpatialReaction
    """The rate function (excluding mass action terms). Currentl only cosntants supported"""
    rate::Symbol
    """Reaction substrates (source and destination)."""
    substrates::Tuple{Vector{Symbol}, Vector{Symbol}}
    """Reaction products (source and destination)."""
    products::Tuple{Vector{Symbol}, Vector{Symbol}}
    """The stoichiometric coefficients of the reactants (source and destination)."""
    substoich::Tuple{Vector{Int64}, Vector{Int64}}
    """The stoichiometric coefficients of the products (source and destination)."""
    prodstoich::Tuple{Vector{Int64}, Vector{Int64}}
    """The net stoichiometric coefficients of all species changed by the reaction (source and destination)."""
    netstoich::Tuple{Vector{Pair{Symbol,Int64}}, Vector{Pair{Symbol,Int64}}}
    """
    `false` (default) if `rate` should be multiplied by mass action terms to give the rate law.
    `true` if `rate` represents the full reaction rate law.
    Currently only `false`, is supported.
    """
    only_use_rate::Bool
    function SpatialReaction(rate, substrates::Tuple{Vector, Vector}, products::Tuple{Vector, Vector}, substoich::Tuple{Vector{Int64}, Vector{Int64}}, prodstoich::Tuple{Vector{Int64}, Vector{Int64}};
                             only_use_rate = false)
        new(rate, substrates, products, substoich, prodstoich,
            get_netstoich.(substrates, products, substoich, prodstoich), only_use_rate)
    end
end

"""
    DiffusionReaction(rate,species)

Simple function to create a diffusion spatial reaction. 
    Equivalent to SpatialReaction(rate,([species],[]),([],[species]),([1],[]),([],[1]))
"""
DiffusionReaction(rate,species) = SpatialReaction(rate,([species],[]),([],[species]),([1],[]),([],[1]))

"""
    OnewaySpatialReaction(rate, substrates, products, substoich, prodstoich)

Simple function to create a spatial reactions where all substrates are in teh soruce compartment, and all products in the destination.
Equivalent to SpatialReaction(rate,(substrates,[]),([],products),(substoich,[]),([],prodstoich))
"""
OnewaySpatialReaction(rate, substrates::Vector, products::Vector, substoich::Vector{Int64}, prodstoich::Vector{Int64}) = SpatialReaction(rate,(substrates,[]),([],products),(substoich,[]),([],prodstoich))



### Lattice Reaction Network Structure ###
# Couples:
# A reaction network (that is simulated within each compartment).
# A set of spatial reactions (denoting interaction between comaprtments).
# A network of compartments (a meta graph that can contain some additional infro for each compartment).
# The lattice is a DiGraph, normals graphs are converted to DiGraphs (with one edge in each direction).
struct LatticeReactionSystem # <: MT.AbstractTimeDependentSystem # Adding this part messes up show, disabling me from creating LRSs
    """The spatial reactions defined between individual nodes."""
    rs::ReactionSystem
    """The spatial reactions defined between individual nodes."""
    spatial_reactions::Vector{SpatialReaction}
    """The graph on which the lattice is defined."""
    lattice::DiGraph
    """Dependent (state) variables representing amount of each species. Must not contain the
    independent variable."""

    function LatticeReactionSystem(rs, spatial_reactions, lattice::DiGraph)
        return new(rs, spatial_reactions, lattice)
    end
    function LatticeReactionSystem(rs, spatial_reactions, lattice::SimpleGraph)
        return new(rs, spatial_reactions, DiGraph(lattice))
    end
end

### ODEProblem ###
# Creates an ODEProblem from a LatticeReactionSystem.
function DiffEqBase.ODEProblem(lrs::LatticeReactionSystem, u0, tspan,
                               p = DiffEqBase.NullParameters(), args...;
                               kwargs...)
    @unpack rs, spatial_reactions, lattice = lrs

    spatial_params = unique(getfield.(spatial_reactions, :rate))
    pV_in, pE_in = split_parameters(p, spatial_params)
    nV, nE = length.([vertices(lattice), edges(lattice)])
    u_idxs = Dict(reverse.(enumerate(Symbolics.getname.(states(rs)))))
    pV_idxes = Dict(reverse.(enumerate(Symbol.(parameters(rs)))))
    pE_idxes = Dict(reverse.(enumerate(spatial_params)))

    u0 = matrix_form(u0, nV, u_idxs)
    pV = matrix_form(pV_in, nV, pV_idxes)
    pE = matrix_form(pE_in, nE, pE_idxes)

    return ODEProblem(build_f(lrs, u_idxs, pE_idxes), u0, tspan, (pV, pE), args...;
                      kwargs...)
end

# Splits parameters into those for the compartments and those for the connections.
split_parameters(parameters::Tuple, spatial_params) = parameters
function split_parameters(parameters::Vector, spatial_params)
    filter(p -> !in(p[1], spatial_params), parameters),
    filter(p -> in(p[1], spatial_params), parameters)
end;

# Converts species and parameters to matrices form.
matrix_form(input::Matrix, args...) = input
function matrix_form(input::Vector{Pair{Symbol, Vector{Float64}}}, n, index_dict)
    mapreduce(permutedims, vcat, last.(sort(input, by = i -> index_dict[i[1]])))
end
function matrix_form(input::Vector, n, index_dict)
    matrix_form(map(i -> (i[2] isa Vector) ? i[1] => i[2] : i[1] => fill(i[2], n), input),
                n, index_dict)
end;

# Creates a function for simulating the spatial ODE with spatial reactions.
function build_f(lrs::LatticeReactionSystem, u_idxs::Dict{Symbol, Int64},
                 pE_idxes::Dict{Symbol, Int64})
    ofunc = ODEFunction(convert(ODESystem, lrs.rs))

    return function internal___spatial___f(du, u, p, t)
        # Updates for non-spatial reactions.
        for comp_i in 1:size(u, 2)
            ofunc((@view du[:, comp_i]), (@view u[:, comp_i]), p[1], t)
        end

        # Updates for spatial reactions.
        for comp_i in 1:size(u, 2)
            for comp_j::Int64 in (lrs.lattice.fadjlist::Vector{Vector{Int64}})[comp_i],
                sr::SpatialReaction in lrs.spatial_reactions::Vector{SpatialReaction}

                rate = get_rate(sr, p[2], (@view u[:, comp_i]), (@view u[:, comp_j]), u_idxs, pE_idxes)
                for stoich in sr.netstoich[1]
                    du[u_idxs[stoich[1]], comp_i] += rate * stoich[2]
                end
                for stoich in sr.netstoich[2]
                    du[u_idxs[stoich[1]], comp_j] += rate * stoich[2]
                end
            end
        end
    end
end

# Get the rate of a specific reaction.
function get_rate(sr, pE, u_src, u_dst, u_idxs, pE_idxes)
    product = pE[pE_idxes[sr.rate]]
    !isempty(sr.substrates[1]) && for (sub,stoich) in zip(sr.substrates[1], sr.substoich[1])
        product *= prod(u_src[u_idxs[sub]]^stoich / factorial(stoich))
    end
    !isempty(sr.substrates[2]) && for (sub,stoich) in zip(sr.substrates[2], sr.substoich[2])
        product *= prod(u_dst[u_idxs[sub]]^stoich / factorial(stoich))
    end
    return product
end