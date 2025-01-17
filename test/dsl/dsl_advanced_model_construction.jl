#! format: off

### Prepares Tests ###

# Fetch packages.
using Catalyst, ModelingToolkit

# Set creates the `t` independent variable.
t = default_t()

### Naming Tests ###

# Basic name test.
let
    rn = @reaction_network SIR1 begin
        k1, S + I --> 2I
        k2, I --> R
    end
    @test nameof(rn) == :SIR1
end

# Advanced name tests.
let
    @parameters k
    @species A(t)
    rx = Reaction(k, [A], nothing)
    function rntest(rn, name)
        @test ModelingToolkit.nameof(rn) == name
        @test isequal(species(rn)[1], ModelingToolkit.unwrap(A))
        @test isequal(parameters(rn)[1], ModelingToolkit.unwrap(k))
        @test reactions(rn)[1] == rx
    end

    function emptyrntest(rn, name)
        @test ModelingToolkit.nameof(rn) == name
        @test numreactions(rn) == 0
        @test numspecies(rn) == 0
    end

    rn = @reaction_network name begin
        @parameters k
        k, A --> 0
    end
    rntest(rn, :name)

    name = :blah
    rn = @reaction_network $name begin
        @parameters k
        k, A --> 0
    end
    rntest(rn, :blah)

    rn = @reaction_network begin
        @parameters k
        k, A --> 0
    end
    rntest(rn, ModelingToolkit.nameof(rn))

    function makern(; name)
        @reaction_network $name begin
        @parameters k
            k, A --> 0
        end
    end
    @named testnet = makern()
    rntest(testnet, :testnet)

    rn = @reaction_network name
    emptyrntest(rn, :name)

    rn = @reaction_network $name
    emptyrntest(rn, :blah)
end


### Test Interpolation Within the DSL ###


# Declares parameters and species used across the test.
@parameters α k k1 k2
@species A(t) B(t) C(t) D(t)

# Tests basic interpolation cases.
let
    AA = A
    AAA = A^2 + B
    rn = @reaction_network rn begin
        @parameters k
        @species A(t) B(t) C(t) D(t)
        k*$AAA, C --> D
    end
    rn2 = ReactionSystem([Reaction(k*AAA, [C], [D])], t; name=:rn)
    @test rn == rn2

    rn = @reaction_network rn begin
        @parameters k
        @species A(t) C(t) D(t)
        k, $AA + C --> D
    end
    rn2 = ReactionSystem([Reaction(k, [AA,C], [D])], t; name=:rn)
    @test rn == rn2
end
let
    BB = B; A2 = A
    rn = @reaction_network rn begin
        @parameters k1 k2
        (k1,k2), C + $A2 + $BB + $A2 <--> $BB + $BB
    end
    rn2 = ReactionSystem([Reaction(k1, [C, A, B], [B], [1,2,1],[2]),
                        Reaction(k2, [B], [C, A, B], [2], [1,2,1])],
                        t; name=:rn)
    @test rn == rn2
end
let
    AA = A
    kk1 = k^2*A
    kk2 = k1+k2
    rn = @reaction_network rn begin
        @parameters α k k1 k2
        α+$kk1*$kk2*$AA, 2*$AA + B --> $AA
    end
    rn2 = ReactionSystem([Reaction(α+kk1*kk2*AA, [A, B], [A], [2, 1], [1])], t; name=:rn)
    @test rn == rn2
end

# Creates a reaction network using `eval` and internal function.
let
    ex = quote
        (Ka, Depot --> Central)
        (CL / Vc, Central --> 0)
    end
    # Line number nodes aren't ignored so have to be manually removed
    Base.remove_linenums!(ex)
    exsys = Catalyst.make_reaction_system(ex)
    sys = @eval Catalyst $exsys
    @test sys isa ReactionSystem
end

# Miscellaneous interpolation tests. Unsure what they do here (not related to DSL).
let
    rx = @reaction k*h, A + 2*B --> 3*C + D
    @parameters k h
    @species A(t) B(t) C(t) D(t)
    @test rx == Reaction(k*h,[A,B],[C,D],[1,2],[3,1])

    ex = k*A^2 + B
    V  = A
    rx = @reaction b+$ex, 2*$V + C--> ∅
    @parameters b
    @test rx == Reaction(b+ex, [A,C], nothing, [2,1], nothing)
end

### Tests Reaction Metadata ###

# Tests construction for various types of reaction metadata.
# Tests reaction metadata accessor functions.
let
    # Creates reactions directly.
    t = default_t()
    @parameters k η
    @species X(t) X2(t)

    metadata1 = [:noise_scaling => η]
    r1 = Reaction(k, [X], [X2], [2], [1]; metadata=metadata1)

    metadata2 = Pair{Symbol,Any}[]
    push!(metadata2, :md_1 => 1.0)
    push!(metadata2, :md_2 => false)
    push!(metadata2, :md_3 => "Hello world")
    push!(metadata2, :md_4 => :sym)
    push!(metadata2, :md_5 => X + X2^k -1)
    push!(metadata2, :md_6 => (0.1, 2.0))
    r2 = Reaction(k, [X], [X2], [2], [1]; metadata=metadata2)

    metadata3 = Pair{Symbol,Any}[]
    r3 = Reaction(k, [X], [X2], [2], [1]; metadata=metadata3)

    # Creates reactions using DSL.
    rs = @reaction_network begin
        @parameters η
        k, 2X --> X2, [noise_scaling=η]
        k, 2X --> X2, [md_1=1.0, md_2=false, md_3="Hello world", md_4=:sym, md_5=X+X2^k-1, md_6=(0.1,2.0)]
        k, 2X --> X2
    end

    # Checks DSL reactions are correct.
    rxs = reactions(rs)
    @test isequal([r1, r2, r3], rxs)
    @test isequal(Catalyst.getmetadata_dict(r1), Catalyst.getmetadata_dict(rxs[1]))
    @test isequal(Catalyst.getmetadata_dict(r2), Catalyst.getmetadata_dict(rxs[2]))
    @test isequal(Catalyst.getmetadata_dict(r3), Catalyst.getmetadata_dict(rxs[3]))

    # Checks that accessor functions works on the DSL.
    @test hasmetadata(rxs[1], :noise_scaling)
    @test !hasmetadata(rxs[1], :md_1)
    @test !hasmetadata(rxs[2], :noise_scaling)
    @test hasmetadata(rxs[2], :md_1)
    @test !hasmetadata(rxs[3], :noise_scaling)
    @test !hasmetadata(rxs[3], :md_1)

    @test isequal(getmetadata(rxs[1], :noise_scaling), η)
    @test isequal(getmetadata(rxs[2], :md_1), 1.0)

    # Test that metadata works for @reaction macro.
    rx1 = @reaction k, 2X --> X2, [noise_scaling=$η]
    rx2 = @reaction k, 2X --> X2, [md_1=1.0, md_2=false, md_3="Hello world", md_4=:sym, md_5=X+X2^k-1, md_6=(0.1,2.0)]
    rx3 = @reaction k, 2X --> X2

    @test isequal([rx1, rx2, rx3], rxs)
    @test isequal(Catalyst.getmetadata_dict(rx1), Catalyst.getmetadata_dict(rxs[1]))
    @test isequal(Catalyst.getmetadata_dict(rx2), Catalyst.getmetadata_dict(rxs[2]))
    @test isequal(Catalyst.getmetadata_dict(rx3), Catalyst.getmetadata_dict(rxs[3]))
end

# Checks that repeated metadata throws errors.
let
    @test_throws LoadError @eval @reaction k, 0 --> X, [md1=1.0, md1=2.0]
    @test_throws LoadError @eval @reaction_network begin
        k, 0 --> X, [md1=1.0, md1=1.0]
    end
end

# Tests for nested metadata.
let
    rn1 = @reaction_network reactions begin
        k1, X1 --> Y1, [md1=1.0,md2=2.0]
        k2, X2 --> Y2, [md1=0.0]
        k3, X3 --> Y3, [md3="Hello world"]
        k, X4 --> Y4, [md4=:sym]
        k, X5 --> Y5, [md4=:sym]
        k6, X6 --> Y6, [md6=0.0]
        k6, Y6 --> X6, [md6=2.0]
        k7, X7 --> Y7, [md7="Hi"]
        k7, X7 --> Y8, [md7="Hi"]
        k8, Y7 --> X7, [md7="Hello"]
        k8, Y8 --> X7, [md7="Hi",md8="Yo"]
    end

    rn2 = @reaction_network reactions begin
        (k1,k2,k3), (X1,X2,X3) --> (Y1,Y2,Y3), ([md1=1.0,md2=2.0],[md1=0.0],[md3="Hello world"])
        k, (X4,X5) --> (Y4,Y5), [md4=:sym]
        (k6, k6), X6 <--> Y6, ([md6=0.0],[md6=2.0])
        (k7,k8), X7 <--> (Y7,Y8), ([md7="Hi"],([md7="Hello"],[md7="Hi",md8="Yo"]))
    end

    @test isequal(rn1, rn2)
end

# Tests that `only_use_rate` option works.
let
    rn1 = @reaction_network reactions begin
        k1*X1, X1 + 2Y1 --> Z1
        k2, 4X2 => Z2 + W3
        k3 + X3, Y3 --> Z3
        2*k4 + X4*Y4, 2X2 + 2Y4 => Z4
        k5, 3X5 --> Z5, [unnecessary_metadata=[1,2,3]]
        k6, X6 => Z6, [unnecessary_metadata=true]
    end

    rn2 = @reaction_network reactions begin
        k1*X1, X1 + 2Y1 --> Z1, [only_use_rate=false]
        k2, 4X2 --> Z2 + W3, [only_use_rate=true]
        k3 + X3, Y3 --> Z3
        2*k4 + X4*Y4, 2X2 + 2Y4 --> Z4, [only_use_rate=true]
        k5, 3X5 --> Z5, [only_use_rate=false, unnecessary_metadata=[1,2,3]]
        k6, X6 --> Z6, [only_use_rate=true, unnecessary_metadata=true]
    end

    @test isequal(rn1,rn2)
end

### Other Tests ###

# Test floating point stoichiometry work.
let
    @parameters k
    @species B(t) C(t) D(t)
    rx1 = Reaction(k,[B,C],[B,D], [2.5,1],[3.5, 2.5])
    rx2 = Reaction(2*k, [B], [D], [1], [2.5])
    rx3 = Reaction(2*k, [B], [D], [2.5], [2])
    @named mixedsys = ReactionSystem([rx1,rx2,rx3],t,[B,C,D],[k])
    mixedsys = complete(mixedsys)
    osys = convert(ODESystem, mixedsys; combinatoric_ratelaws=false)
    rn = @reaction_network mixedsys begin
        @parameters k
        k, 2.5*B + C --> 3.5*B + 2.5*D
        2*k, B --> 2.5*D
        2*k, 2.5*B --> 2*D
    end
    @test rn == mixedsys
end

# Test that variables that appear only in rates and aren't ps
# are categorized as species.
let
    rn = @reaction_network begin
        @parameters k k2 n
        @species A(t) B(t) C(t) D(t) H(t)
        π*k*D*hill(B,k2,B*D*H,n), 3*A  --> 2*C
    end
    @parameters k k2 n
    @species A(t) B(t) C(t) D(t) H(t)
    @test issetequal([A,B,C,D,H], species(rn))
    @test issetequal([k,k2,n], parameters(rn))
end

# Test species have the right metadata via the DSL.
let
    rn = @reaction_network begin
        k, 2*A + B --> C
    end
    @test issetequal(unknowns(rn), species(rn))
    @test all(isspecies, species(rn))

    rn2 = @reaction_network begin
        @species A(t) = 1 B(t) = 2 [isbcspecies = true]
        k, A + 2*B --> 2*B
    end
    @unpack A,B = rn2
    D = default_time_deriv()
    eq = D(B) ~ -B
    @named osys = ODESystem([eq], t)
    @named rn2 = extend(osys, rn2)
    rn2 = complete(rn2)
    @test issetequal(unknowns(rn2), species(rn2))
    @test all(isspecies, species(rn))
    @test Catalyst.isbc(ModelingToolkit.value(B))
    @test Catalyst.isbc(ModelingToolkit.value(A)) == false
    osys2 = complete(convert(ODESystem, rn2))
    @test issetequal(unknowns(osys2), unknowns(rn2))
    @test length(equations(osys2)) == 2
end

# Array variables test.
let
    rn = @reaction_network arrtest begin
        @parameters k[1:2] a
        @variables (V(t))[1:2] W(t)
        @species (X(t))[1:2] Y(t)
        k[1]*a+k[2], X[1] + V[1]*X[2] --> V[2]*W*Y + B*C
    end

    @parameters k[1:2] a B
    @variables (V(t))[1:2] W(t)
    @species (X(t))[1:2] Y(t) C(t)
    rx = Reaction(k[1]*a+k[2], [X[1], X[2]], [Y, C], [1, V[1]], [V[2] * W, B])
    @named arrtest = ReactionSystem([rx], t)
    @test arrtest == rn

    rn = @reaction_network twostate begin
        @parameters k[1:2]
        @species (X(t))[1:2]
        (k[1],k[2]), X[1] <--> X[2]
    end
    
    @parameters k[1:2]
    @species (X(t))[1:2]
    rx1 = Reaction(k[1], [X[1]], [X[2]])
    rx2 = Reaction(k[2], [X[2]], [X[1]])
    @named twostate = ReactionSystem([rx1, rx2], t)
    @test twostate == rn
end

############## tests related to hybrid systems ###################
let
    t = default_t()
    D = default_time_deriv()
    @parameters λ k
    @variables V(t)
    @species A(t)
    rx = Reaction(k*V, [], [A])
    eq = D(V) ~ λ*V
    cevents = [[V ~ 2.0] => [V ~ V/2, A ~ A/2]]
    @named hybrid = ReactionSystem([rx, eq], t; continuous_events = cevents)
    rn = @reaction_network hybrid begin
        @parameters λ
        k*V, 0 --> A
        @equations D(V) ~ λ*V
        @continuous_events begin
            [V ~ 2.0] => [V ~ V/2, A ~ A/2]
        end
    end        
    @test hybrid == rn
end