"""
$(DocStringExtensions.README)
"""
module Catalyst

using DocStringExtensions
using SparseArrays, DiffEqBase, Reexport, Setfield, EnumX
using LaTeXStrings, Latexify, Requires
using LinearAlgebra, Combinatorics
using JumpProcesses: JumpProcesses, JumpProblem,
                     MassActionJump, ConstantRateJump, VariableRateJump,
                     SpatialMassActionJump, CartesianGrid, CartesianGridRej

# ModelingToolkit imports and convenience functions we use
using ModelingToolkit
const MT = ModelingToolkit
using DynamicQuantities #, Unitful # Having Unitful here as well currently gives an error.

@reexport using ModelingToolkit
using Symbolics
using LinearAlgebra
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

import Symbolics: BasicSymbolic
using Symbolics: iscall, sorted_arguments
using ModelingToolkit: Symbolic, value, get_unknowns, get_ps, get_iv, get_systems,
                       get_eqs, get_defaults, toparam, get_var_to_name, get_observed,
                       getvar, has_iv

import ModelingToolkit: get_variables, namespace_expr, namespace_equation, get_variables!,
                        modified_unknowns!, validate, namespace_variables,
                        namespace_parameters, rename, renamespace, getname, flatten,
                        is_alg_equation, is_diff_equation, collect_vars!,
                        eqtype_supports_collect_vars

# internal but needed ModelingToolkit functions
import ModelingToolkit: check_variables,
                        check_parameters, _iszero, _merge, check_units,
                        get_unit, check_equations, iscomplete

import Base: (==), hash, size, getindex, setindex, isless, Sort.defalg, length, show
import MacroTools, Graphs
using MacroTools: striplines
import Graphs: DiGraph, SimpleGraph, SimpleDiGraph, vertices, edges, add_vertices!, nv, ne
import DataStructures: OrderedDict, OrderedSet
import Parameters: @with_kw_noshow
import Symbolics: occursin, wrap
import Symbolics.RewriteHelpers: hasnode, replacenode
import SymbolicUtils: getmetadata, hasmetadata, setmetadata

# globals for the modulate
function default_time_deriv()
    return ModelingToolkit.D_nounits
end
function default_t()
    return ModelingToolkit.t_nounits
end
const DEFAULT_IV = default_t()
const DEFAULT_IV_SYM = Symbol(DEFAULT_IV)
export default_t, default_time_deriv

### Package Constants ###

# Union type of types that can occur in expressions.
const ExprValues = Union{Expr, Symbol, Float64, Int, Bool}

# The symbol used for conserved quantities in conservation law eliminations.
const CONSERVED_CONSTANT_SYMBOL = :Γ

# Declares symbols which may neither be used as parameters nor unknowns.
const forbidden_symbols_skip = Set([:ℯ, :pi, :π, :t, :∅])
const forbidden_symbols_error = union(Set([:im, :nothing, CONSERVED_CONSTANT_SYMBOL]),
    forbidden_symbols_skip)
const forbidden_variables_error = let
    fvars = copy(forbidden_symbols_error)
    delete!(fvars, :t)
    fvars
end

### Package Main ###

# The `Reaction` structure and its functions.
include("reaction.jl")
export isspecies
export Reaction, PhysicalScale

# Union type for `Reaction`s and `Equation`s.
const CatalystEqType = Union{Reaction, Equation}

# The `ReactionSystem` structure and its functions.
include("reactionsystem.jl")
export ReactionSystem, isspatial
export species, nonspecies, reactions, nonreactions, speciesmap, paramsmap
export numspecies, numreactions, setdefaults!
export make_empty_network
export dependants, dependents, substoichmat, prodstoichmat, netstoichmat
export isautonomous
export reactionrates
export isequivalent
export set_default_noise_scaling

# depreciated functions to remove in future releases
export params, numparams

# Conversions of the `ReactionSystem` structure.
include("reactionsystem_conversions.jl")
export ODEProblem, SDEProblem, JumpProblem, NonlinearProblem, DiscreteProblem,
       SteadyStateProblem, JumpInputs
export ismassaction, oderatelaw, jumpratelaw
export symmap_to_varmap

# reaction_network macro
include("expression_utils.jl")
include("dsl.jl")
export @reaction_network, @network_component, @reaction, @species

# Network analysis functionality.
include("network_analysis.jl")
export reactioncomplexmap, reactioncomplexes, incidencemat
export complexstoichmat, laplacianmat, fluxmat, massactionvector, complexoutgoingmat
export incidencematgraph, linkageclasses, stronglinkageclasses,
       terminallinkageclasses, deficiency, subnetworks
export linkagedeficiencies, isreversible, isweaklyreversible
export conservationlaws, conservedquantities, conservedequations, conservationlaw_constants

# registers CRN specific functions using Symbolics.jl
include("registered_functions.jl")
export mm, mmr, hill, hillr, hillar

# functions to query network properties

# for Latex printing of ReactionSystems
include("latexify_recipes.jl")

# for making and saving graphs/plots
include("plotting.jl")

# for creating compounds
include("chemistry_functionality.jl")
export @compound, @compounds
export iscompound, components, coefficients, component_coefficients
export balance_reaction, balance_system

# Functionality for computing the stability of system steady states.
include("steady_state_stability.jl")
export steady_state_stability, steady_state_jac

### Extensions ###

# HomotopyContinuation
function hc_steady_states end
export hc_steady_states

# StructuralIdentifiability
function make_si_ode end
export make_si_ode

# GraphMakie: functionality for plotting species-reaction graphs and complexes
function plot_network end
function plot_complexes end
export plot_network, plot_complexes

### Spatial Reaction Networks ###

# Spatial reactions.
include("spatial_reaction_systems/spatial_reactions.jl")
export TransportReaction, @transport_reaction
export isedgeparameter

# Lattice reaction systems.
include("spatial_reaction_systems/lattice_reaction_systems.jl")
export LatticeReactionSystem
export spatial_species, vertex_parameters, edge_parameters
export CartesianGrid, CartesianGridReJ # (Implemented in JumpProcesses)
export has_cartesian_lattice, has_masked_lattice, has_grid_lattice, has_graph_lattice,
       grid_dims, grid_size
export make_edge_p_values, make_directed_edge_values

# Specific spatial problem types.
include("spatial_reaction_systems/spatial_ODE_systems.jl")
include("spatial_reaction_systems/lattice_jump_systems.jl")

# General spatial modelling utility functions.
include("spatial_reaction_systems/utility.jl")

# Methods for interfacing with from LatticeReactionSystem derived problems, integrators, and solutions.
include("spatial_reaction_systems/lattice_sim_struct_interfacing.jl")
export lat_getp, lat_setp!, lat_getu, lat_setu!, rebuild_lat_internals!

# Functions for plotting of lattice simulations (most of the code is in extensions, not here).
include("spatial_reaction_systems/lattice_simulation_plotting.jl")
export lattice_plot, lattice_animation, lattice_kymograph

### ReactionSystem Serialisation ###
# Has to be at the end (because it uses records of all metadata declared by Catalyst).
include("reactionsystem_serialisation/serialisation_support.jl")
include("reactionsystem_serialisation/serialise_fields.jl")
include("reactionsystem_serialisation/serialise_reactionsystem.jl")
export save_reactionsystem

end # module
