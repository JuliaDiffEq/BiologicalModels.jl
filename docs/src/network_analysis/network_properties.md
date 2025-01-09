# The `NetworkProperties` Object

# [Caching of Network Properties in `ReactionSystems`](@id network_analysis_caching_properties)
When calling many of the network API functions, Catalyst calculates and caches
in `rn` a variety of information. For example the first call to
```julia
rcs,B = reactioncomplexes(rn)
```
calculates, caches, and returns the reaction complexes, `rcs`, and the incidence
matrix, `B`, of `rn`. Subsequent calls simply return `rcs` and `B` from the
cache.

Similarly, the first call to
```julia
N = netstoichmat(rn)
```
calculates, caches and returns the net stoichiometry matrix. Subsequent calls
then simply return the cached value of `N`. Caching such information means users
do not need to manually know which subsets of network properties are needed for
a given calculation (like the deficiency). Generally only
```julia
rcs,B = reactioncomplexes(rn)    # must be called once to cache rcs and B
any_other_network_property(rn)
```
should work to calculate a desired network property, with the API doc strings
indicating when `reactioncomplexes(rn)` must be called at least once before a
given function is used.

Because of the caching of network properties, subsequent calls to most API
functions will be fast, simply returning the previously calculated and cached
values. In some cases it may be desirable to reset the cache and recalculate
these properties. This can be done by calling
```julia
Catalyst.reset_networkproperties!(rn)
```
Network property functions will then recalculate their associated properties and
cache the new values the next time they are called.
