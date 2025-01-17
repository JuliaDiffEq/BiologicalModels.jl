### Preparations ###

# Fetch packages.
using Catalyst, Test
using Catalyst: get_symbolics
using ModelingToolkit: value, get_variables!, collect_vars!, eqtype_supports_collect_vars
using Catalyst: has_physical_scale, get_physical_scale

# Sets the default `t` to use.
t = default_t()

### Reaction Constructor Tests ###

# Checks that `Reaction`s can be successfully created using various complicated inputs.
# Checks that the `Reaction`s have the correct type, and the correct net stoichiometries are generated.
let
    # Declare symbolic variables.
    @parameters k n1 n2::Int32 x [isconstantspecies=true]
    @species X(t) Y(t) Z(t)
    @variables A(t)

    # Tries for different types of rates (should not matter).
    for rate in (k, k*A, 2, 3.0, 4//3)
        # Creates `Reaction`s.
        rx1 = Reaction(rate, [X], [])
        rx2 = Reaction(rate, [x], [Y], [1.5], [1])
        rx3 = Reaction(rate, [x, X], [], [n1 + n2, n2], [])
        rx4 = Reaction(rate, [X, Y], [X, Y, Z], [2//3, 3], [1//3, 1, 2])
        rx5 = Reaction(rate, [X, Y], [X, Y, Z], [2, 3], [1, n1, n2])
        rx6 = Reaction(rate, [X], [x], [n1], [1])

        # Check `Reaction` types.
        @test rx1 isa Reaction{Any,Int64}
        @test rx2 isa Reaction{Any,Float64}
        @test rx3 isa Reaction{Any,Any}
        @test rx4 isa Reaction{Any,Rational{Int64}}
        @test rx5 isa Reaction{Any,Any}
        @test rx6 isa Reaction{Any,Any}

        # Check `Reaction` net stoichiometries.
        issetequal(rx1.netstoich, [X => -1])
        issetequal(rx2.netstoich, [x => -1.5, Y => 1.0])
        issetequal(rx3.netstoich, [x => -n1 - n2, X => -n2])
        issetequal(rx4.netstoich, [X => -1//3, Y => -2//1, Z => 2//1])
        issetequal(rx5.netstoich, [X => -1, Y => n1 - 3, Z => n2])
        issetequal(rx6.netstoich, [X => -n1, x => 1])
    end
end

# Tests that various `Reaction` constructors gives identical inputs.
let
    # Declare symbolic variables.
    @parameters k n1 n2::Int32
    @species X(t) Y(t) Z(t)
    @variables A(t)

    # Tests that the three-argument constructor generates correct result.
    @test Reaction(k*A, [X], [Y, Z]) == Reaction(k*A, [X], [Y, Z], [1], [1, 1])

    # Tests that `[]` and `nothing` can be used interchangeably.
    @test Reaction(k*A, [X, Z], nothing) == Reaction(k*A, [X, Z], [])
    @test Reaction(k*A, nothing, [Y, Z]) == Reaction(k*A, [], [Y, Z])
    @test Reaction(k*A, [X, Z], nothing, [n1 + n2, 2], nothing) == Reaction(k*A, [X, Z], [], [n1 + n2, 2], [])
    @test Reaction(k*A, nothing, [Y, Z], nothing, [n1 + n2, 2]) == Reaction(k*A, [], [Y, Z], [], [n1 + n2, 2])
end

# Tests that various incorrect inputs yields errors.
let
    # Declare symbolic variables.
    @parameters k n1 n2::Int32
    @species X(t) Y(t) Z(t)
    @variables A(t)

    # Neither substrates nor products.
    @test_throws ArgumentError Reaction(k*A, [], [])

    # Substrate vector not of equal length to substrate stoichiometry vector.
    @test_throws ArgumentError Reaction(k*A, [X, X, Z], [], [1, 2], [])

    # Product vector not of equal length to product stoichiometry vector.
    @test_throws ArgumentError Reaction(k*A, [], [X, X, Z], [], [1, 2])

    # Repeated substrates.
    @test_throws ArgumentError Reaction(k*A, [X, X, Z], [])

    # Repeated products.
    @test_throws ArgumentError Reaction(k*A, [], [Y, Z, Z])

    # Non-valid reactants (parameter or variable).
    @test_throws ArgumentError Reaction(k*A, [], [A])
    @test_throws ArgumentError Reaction(k*A, [], [k])
end


### Test Basic Accessors ###

# Tests the `get_variables` function.
let
    # Declare symbolic variables.
    @parameters k1 k2 n1 n2 η1 η2 p
    @species X(t) Y(t) Z(t)
    @variables A(t)

    # Create `Reaction`s.
    rx1 = Reaction(k1, [], [X])
    rx2 = Reaction(k1 + k2, [X], [Y], [1], [n1]; metadata = [:noise_scaling => η1])
    rx3 = Reaction(k1 + k2 + A, [X], [X, Y, Z], [1], [n1 + n2, 2, 1])
    rx4 = Reaction(X + t, [], [Y]; metadata = [:noise_scaling => η1 + η2])

    # Test `get_variables!`.
    @test issetequal(get_variables!([value(p)], rx1), [k1, X, p])
    @test issetequal(get_variables!([value(p)], rx2), [k1, k2, X, Y, n1, η1, p])
    @test issetequal(get_variables!([value(p)], rx3), [k1, k2, A, X, Y, Z, n1, n2, p])
    @test issetequal(get_variables!([value(p)], rx4), [X, t, Y, η1, η2, p])

    # Test `get_symbolics`.
    @test issetequal(get_symbolics(rx1), [k1, X])
    @test issetequal(get_symbolics(rx2), [k1, k2, X, Y, n1, η1])
    @test issetequal(get_symbolics(rx3), [k1, k2, A, X, Y, Z, n1, n2])
    @test issetequal(get_symbolics(rx4), [X, t, Y, η1, η2])
end

### Tests Metadata ###

# Tests creation.
# Tests basic accessor functions.
# Tests that repeated metadata entries are not permitted.
# Tests that attempting to access non-existent metadata throws an error.
let
    @parameters k
    @species X(t) X2(t)

    metadata = [:noise_scaling => 0.0]
    r = Reaction(k, [X], [X2], [2], [1]; metadata=metadata)

    @test Catalyst.getmetadata_dict(r) == [:noise_scaling => 0.0]
    @test hasmetadata(r, :noise_scaling)
    @test !hasmetadata(r, :nonexisting_metadata)
    @test getmetadata(r, :noise_scaling) == 0.0
    @test_throws Exception getmetadata(r, :misc)
    setmetadata(r, :test_metadata, 1234)
    @test getmetadata(r, :test_metadata) == 1234
    setmetadata(r, :test_metadata, 1111)
    @test getmetadata(r, :test_metadata) == 1111

    metadata_repeated = [:noise_scaling => 0.0, :noise_scaling => 1.0, :metadata_entry => "unused"]
    @test_throws Exception Reaction(k, [X], [X2], [2], [1]; metadata=metadata_repeated)
end

# Tests accessors for system without metadata.
let
    @parameters k
    @species X(t) X2(t)

    metadata = Pair{Symbol,Any}[]
    r1 = Reaction(k, [X], [X2], [2], [1])
    r2 = Reaction(k, [X], [X2], [2], [1]; metadata=metadata)

    @test isequal(r1, r2)
    @test Catalyst.getmetadata_dict(r1) == Pair{Symbol,Any}[]
    @test !hasmetadata(r1, :md)
end

# Tests creation.
# Tests basic accessor functions.
# Tests various metadata types.
let
    @parameters k
    @species X(t) X2(t)

    metadata = Pair{Symbol,Any}[]
    push!(metadata, :md_1 => 1.0)
    push!(metadata, :md_2 => false)
    push!(metadata, :md_3 => "Hello world")
    push!(metadata, :md_4 => :sym)
    push!(metadata, :md_5 => X + X2^k -1)
    push!(metadata, :md_6 => (0.1, 2.0))
    r = Reaction(k, [X], [X2], [2], [1]; metadata=metadata)

    @test Catalyst.getmetadata_dict(r) isa Vector{Pair{Symbol,Any}}
    @test hasmetadata(r, :md_1)
    @test hasmetadata(r, :md_2)
    @test hasmetadata(r, :md_3)
    @test hasmetadata(r, :md_4)
    @test hasmetadata(r, :md_5)
    @test hasmetadata(r, :md_6)
    @test !hasmetadata(r, :md_8)

    @test isequal(getmetadata(r, :md_1), 1.0)
    @test isequal(getmetadata(r, :md_2), false)
    @test isequal(getmetadata(r, :md_3), "Hello world")
    @test isequal(getmetadata(r, :md_4), :sym)
    @test isequal(getmetadata(r, :md_5), X + X2^k -1)
    @test isequal(getmetadata(r, :md_6), (0.1, 2.0))
end

# Tests the noise scaling metadata.
let
    @parameters k η
    @species X(t) X2(t)

    r1 = Reaction(k, [X], [X2], [2], [1])
    r2 = Reaction(k, [X], [X2], [2], [1]; metadata=[:noise_scaling => η])

    @test !Catalyst.hasnoisescaling(r1)
    @test Catalyst.hasnoisescaling(r2)
    @test_throws Exception Catalyst.getnoisescaling(r1)
    @test isequal(Catalyst.getnoisescaling(r2), η)
end

# Tests the description metadata.
let
    t = default_t()
    @parameters k η
    @species X(t) X2(t)

    r1 = Reaction(k, [X], [X2], [2], [1])
    r2 = Reaction(k, [X], [X2], [2], [1]; metadata=[:description => "A reaction"])

    @test !Catalyst.hasdescription(r1)
    @test Catalyst.hasdescription(r2)
    @test_throws Exception Catalyst.getdescription(r1)
    @test isequal(Catalyst.getdescription(r2), "A reaction")
end

# Tests the misc metadata.
let
    t = default_t()
    @parameters k η
    @species X(t) X2(t)

    r1 = Reaction(k, [X], [X2], [2], [1])
    r2 = Reaction(k, [X], [X2], [2], [1]; metadata=[:misc => ('M', :M)])

    @test !Catalyst.hasmisc(r1)
    @test Catalyst.hasmisc(r2)
    @test_throws Exception Catalyst.getmisc(r1)
    @test isequal(Catalyst.getmisc(r2), ('M', :M))
end

# tests for collect_vars!
let
    t = default_t()
    @variables E(t) F(t)
    @species A(t) B(t) C(t) D(t) 
    @parameters k1, k2, η

    rx = Reaction(k1*E, [A, B], [C], [k2*D, 3], [F], metadata = [:noise_scaling => η])
    us = Set()
    ps = Set()
    @test eqtype_supports_collect_vars(rx) == true
    collect_vars!(us, ps, rx, t)
    @test us == Set((A, B, C, D, E, F))
    @test ps == Set((k1, k2, η))
end

# tests for PhysicalScales
let
    t = default_t()
    @species A(t) B(t) C(t)
    @parameters k1, k2

    rx = Reaction(k1, [A], [B], [2], [1]; metadata = [:physical_scale => PhysicalScale.Jump])
    @test has_physical_scale(rx)
    @test get_physical_scale(rx) == PhysicalScale.Jump

    rx2 = Reaction(k1, [A], [B], [2], [1])
    @test has_physical_scale(rx2) == false

    rx3 = Reaction(k1, [A], [B], [2], [1]; 
        metadata = [:physical_scale => PhysicalScale.Jump, :noise_scaling => 0.1])
    @test has_physical_scale(rx3)
    @test get_physical_scale(rx3) == PhysicalScale.Jump    
end