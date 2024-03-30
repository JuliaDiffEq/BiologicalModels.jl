#! format: off

### Prepares Tests ###

# Fetch packages
using Catalyst, JumpProcesses, NonlinearSolve, OrdinaryDiffEq, Plots, StochasticDiffEq
import ModelingToolkit: getp, getu, setp, setu

# Sets rnd number.
using StableRNGs
rng = StableRNG(12345)
seed = rand(rnd, 1:100)

### Basic Tests ###

# Prepares a model and its problems, integrators, and solutions.
begin
    model = @reaction_network begin
        @observables XY ~ X + Y
        (kp,kd), 0 <--> X
        (k1,k2), X <--> Y
    end
    @unpack XY, X, Y, kp, kd, k1, k2 = model

    u0_vals = [X => 4, Y => 5]
    tspan = (0.0, 10.0)
    p_vals = [kp => 1.0, kd => 0.1, k1 => 0.25, k2 => 0.5]

    oprob = ODEProblem(model, u0_vals, tspan, p_vals)
    sprob = SDEProblem(model,u0_vals, tspan, p_vals)
    dprob = DiscreteProblem(model, u0_vals, tspan, p_vals)
    jprob = JumpProblem(model, deepcopy(dprob), Direct(); rng)
    nprob = NonlinearProblem(model, u0_vals, p_vals)
    problems = [oprob, sprob, dprob, jprob, nprob]

    oint = init(oprob, Tsit5(); save_everystep=false)
    sint = init(sprob, ImplicitEM(); save_everystep=false)
    jint = init(jprob, SSAStepper())
    nint = init(nprob, NewtonRaphson(); save_everystep=false)
    integrators = [oint, sint, jint, nint]
    
    osol = solve(oprob, Tsit5())
    ssol = solve(sprob, ImplicitEM(); seed)
    jsol = solve(jprob, SSAStepper(); seed)
    nsol = solve(nprob, NewtonRaphson())
    sols = [osol, ssol, jsol, nsol]
end

# Tests problem indexing and updating.
let 
    for prob in deepcopy(problems)
        # Get u values (including observables).
        @test prob[X] == prob[model.X] == prob[:X] == 4
        @test prob[XY] == prob[model.XY] == prob[:XY] == 9
        @test prob[[XY,Y]] == prob[[model.XY,model.Y]] == prob[[:XY,:Y]] == [9, 5]
        @test_broken prob[(XY,Y)] == prob[(model.XY,model.Y)] == prob[(:XY,:Y)] == (9, 5)
        @test getu(prob, X)(prob) == getu(prob, model.X)(prob) == getu(prob, :X)(prob) == 4
        @test getu(prob, XY)(prob) == getu(prob, model.XY)(prob) == getu(prob, :XY)(prob) == 9 
        @test getu(prob, [XY,Y])(prob) == getu(prob, [model.XY,model.Y])(prob) == getu(prob, [:XY,:Y])(prob) == [9, 5]  
        @test getu(prob, (XY,Y))(prob) == getu(prob, (model.XY,model.Y))(prob) == getu(prob, (:XY,:Y))(prob) == (9, 5)

        # Set u values.
        prob[X] = 20
        @test prob[X] == 20
        prob[model.X] = 30
        @test prob[X] == 30
        prob[:X] = 40
        @test prob[X] == 40
        setu(prob, X)(prob, 50)
        @test prob[X] == 50
        setu(prob, model.X)(prob, 60)
        @test prob[X] == 60
        setu(prob, :X)(prob, 70)
        @test prob[X] == 70

        # Get p values.
        @test prob.ps[kp] == prob.ps[model.kp] == prob.ps[:kp] == 1.0    
        @test prob.ps[[k1,k2]] == prob.ps[[model.k1,model.k2]] == prob.ps[[:k1,:k2]] == [0.25, 0.5]
        @test prob.ps[(k1,k2)] == prob.ps[(model.k1,model.k2)] == prob.ps[(:k1,:k2)] == (0.25, 0.5)
        @test getp(prob, kp)(prob) == getp(prob, model.kp)(prob) == getp(prob, :kp)(prob) == 1.0
        @test getp(prob, [k1,k2])(prob) == getp(prob, [model.k1,model.k2])(prob) == getp(prob, [:k1,:k2])(prob) == [0.25, 0.5]
        @test getp(prob, (k1,k2))(prob) == getp(prob, (model.k1,model.k2))(prob) == getp(prob, (:k1,:k2))(prob) == (0.25, 0.5)
        
        # Set p values.
        prob.ps[kp] = 2.0
        @test prob.ps[kp] == 2.0
        prob.ps[model.kp] = 3.0
        @test prob.ps[kp] == 3.0
        prob.ps[:kp] = 4.0
        @test prob.ps[kp] == 4.0
        setp(prob, kp)(prob, 5.0)
        @test prob.ps[kp] == 5.0
        setp(prob, model.kp)(prob, 6.0)
        @test prob.ps[kp] == 6.0
        setp(prob, :kp)(prob, 7.0)
        @test prob.ps[kp] == 7.0
    end
end

# Test remake function.
let 
    for prob in deepcopy(problems)
        # Remake for all u0s.
        rp = remake(prob; u0 = [X => 1, Y => 2])
        @test rp[[X, Y]] == [1, 2]
        rp = remake(prob; u0 = [model.X => 3, model.Y => 4])
        @test rp[[X, Y]] == [3, 4]
        rp = remake(prob; u0 = [:X => 5, :Y => 6])
        @test rp[[X, Y]] == [5, 6]

        # Remake for a single u0.
        rp = remake(prob; u0 = [Y => 7])
        @test rp[[X, Y]] == [4, 7]
        rp = remake(prob; u0 = [model.Y => 8])
        @test rp[[X, Y]] == [4, 8]
        rp = remake(prob; u0 = [:Y => 9])
        @test rp[[X, Y]] == [4, 9]

        # Remake for all ps.
        rp = remake(prob; p = [kp => 1.0, kd => 2.0, k1 => 3.0, k2 => 4.0])
        @test rp.ps[[kp, kd, k1, k2]] == [1.0, 2.0, 3.0, 4.0]
        rp = remake(prob; p = [model.kp => 5.0, model.kd => 6.0, model.k1 => 7.0, model.k2 => 8.0])
        @test rp.ps[[kp, kd, k1, k2]] == [5.0, 6.0, 7.0, 8.0]
        rp = remake(prob; p = [:kp => 9.0, :kd => 10.0, :k1 => 11.0, :k2 => 12.0])
        @test rp.ps[[kp, kd, k1, k2]] == [9.0, 10.0, 11.0, 12.0]

        # Remake for a single p.
        rp = remake(prob; p = [k2 => 13.0])
        @test rp.ps[[kp, kd, k1, k2]] == [1.0, 0.1, 0.25, 13.0]
        rp = remake(prob; p = [model.k2 => 14.0])
        @test rp.ps[[kp, kd, k1, k2]] == [1.0, 0.1, 0.25, 14.0]
        rp = remake(prob; p = [:k2 => 15.0])
        @test rp.ps[[kp, kd, k1, k2]] == [1.0, 0.1, 0.25, 15.0]
    end
end


# Test integrator indexing.
let 
    @test_broken false # NOTE: Multiple problems for `nint`.
    @test_broken false # NOTE: Multiple problems for `jint`.
    for int in deepcopy([oint, sint])
        # Get u values.
        @test int[X] == int[model.X] == int[:X] == 4
        @test int[XY] == int[model.XY] == int[:XY] == 9
        @test int[[XY,Y]] == int[[model.XY,model.Y]] == int[[:XY,:Y]] == [9, 5]
        @test int[(XY,Y)] == int[(model.XY,model.Y)] == int[(:XY,:Y)] == (9, 5)
        @test getu(int, X)(int) == getu(int, model.X)(int) == getu(int, :X)(int) == 4
        @test getu(int, XY)(int) == getu(int, model.XY)(int) == getu(int, :XY)(int) == 9 
        @test getu(int, [XY,Y])(int) == getu(int, [model.XY,model.Y])(int) == getu(int, [:XY,:Y])(int) == [9, 5]  
        @test getu(int, (XY,Y))(int) == getu(int, (model.XY,model.Y))(int) == getu(int, (:XY,:Y))(int) == (9, 5)

        # Set u values.
        int[X] = 20
        @test int[X] == 20
        int[model.X] = 30
        @test int[X] == 30
        int[:X] = 40
        @test int[X] == 40
        setu(int, X)(int, 50)
        @test int[X] == 50
        setu(int, model.X)(int, 60)
        @test int[X] == 60
        setu(int, :X)(int, 70)
        @test int[X] == 70

        # Get p values.
        @test int.ps[kp] == int.ps[model.kp] == int.ps[:kp] == 1.0    
        @test int.ps[[k1,k2]] == int.ps[[model.k1,model.k2]] == int.ps[[:k1,:k2]] == [0.25, 0.5]
        @test int.ps[(k1,k2)] == int.ps[(model.k1,model.k2)] == int.ps[(:k1,:k2)] == (0.25, 0.5)
        @test getp(int, kp)(int) == getp(int, model.kp)(int) == getp(int, :kp)(int) == 1.0
        @test getp(int, [k1,k2])(int) == getp(int, [model.k1,model.k2])(int) == getp(int, [:k1,:k2])(int) == [0.25, 0.5]
        @test getp(int, (k1,k2))(int) == getp(int, (model.k1,model.k2))(int) == getp(int, (:k1,:k2))(int) == (0.25, 0.5)
        
        # Set p values.
        int.ps[kp] = 2.0
        @test int.ps[kp] == 2.0
        int.ps[model.kp] = 3.0
        @test int.ps[kp] == 3.0
        int.ps[:kp] = 4.0
        @test int.ps[kp] == 4.0
        setp(int, kp)(int, 5.0)
        @test int.ps[kp] == 5.0
        setp(int, model.kp)(int, 6.0)
        @test int.ps[kp] == 6.0
        setp(int, :kp)(int, 7.0)
        @test int.ps[kp] == 7.0
    end
end

# Test solve's save_idxs argument.
let 
    for (prob, solver) in zip(deepcopy([oprob, sprob, jprob]), [Tsit5(), ImplicitEM(), SSAStepper()])
        # Save single variable
        @test_broken solve(prob, solver; seed, save_idxs=X)[X][1] == 4
        @test_broken solve(prob, solver; seed, save_idxs=model.X)[X][1] == 4
        @test_broken solve(prob, solver; seed, save_idxs=:X)[X][1] == 4

        # Save observable.
        @test_broken solve(prob, solver; seed, save_idxs=XY)[XY][1] == 9
        @test_broken solve(prob, solver; seed, save_idxs=model.XY)[XY][1] == 9
        @test_broken solve(prob, solver; seed, save_idxs=:XY)[XY][1] == 9

        # Save vector of stuff.
        @test_broken solve(prob, solver; seed, save_idxs=[XY,Y])[[XY,Y]][1] == [9, 5]
        @test_broken solve(prob, solver; seed, save_idxs=[model.XY,model.Y])[[model.XY,model.Y]][1] == [9, 5]
        @test_broken solve(prob, solver; seed, save_idxs=[:XY,:Y])[[:XY,:Y]][1] == [9, 5]
    end
end

# Tests solution indexing.
let 
    for sol in deepcopy([osol, ssol, jsol])
        # Get u values.
        @test sol[X][1] == sol[model.X][1] == sol[:X][1] == 4
        @test sol[XY][1] == sol[model.XY][1] == sol[:XY][1] == 9
        @test sol[[XY,Y]][1] == sol[[model.XY,model.Y]][1] == sol[[:XY,:Y]][1] == [9, 5]
        @test sol[(XY,Y)][1] == sol[(model.XY,model.Y)][1] == sol[(:XY,:Y)][1] == (9, 5)
        @test getu(sol, X)(sol)[1] == getu(sol, model.X)(sol)[1] == getu(sol, :X)(sol)[1] == 4
        @test getu(sol, XY)(sol)[1] == getu(sol, model.XY)(sol)[1] == getu(sol, :XY)(sol)[1] == 9 
        @test getu(sol, [XY,Y])(sol)[1] == getu(sol, [model.XY,model.Y])(sol)[1] == getu(sol, [:XY,:Y])(sol)[1] == [9, 5]  
        @test getu(sol, (XY,Y))(sol)[1] == getu(sol, (model.XY,model.Y))(sol)[1] == getu(sol, (:XY,:Y))(sol)[1] == (9, 5)       

        # Get u values via idxs and functional call.
        @test osol(0.0; idxs=X) == osol(0.0; idxs=X) == osol(0.0; idxs=X) == 4
        @test osol(0.0; idxs=XY) == osol(0.0; idxs=XY) == osol(0.0; idxs=XY) == 9
        @test_broken osol(0.0; idxs=[model.Y,model.XY]) == osol(0.0; idxs=[model.Y,model.XY]) == osol(0.0; idxs=[model.XY,model.X]) == [9, 5]
        @test_broken osol(0.0; idxs=(:Y,:XY)) == osol(0.0; idxs=(:Y,:XY)) == osol(0.0; idxs=(:XY,:Y)) == (9, 5)

        # Get p values.
        @test sol.ps[kp] == sol.ps[model.kp] == sol.ps[:kp] == 1.0    
        @test sol.ps[[k1,k2]] == sol.ps[[model.k1,model.k2]] == sol.ps[[:k1,:k2]] == [0.25, 0.5]
        @test sol.ps[(k1,k2)] == sol.ps[(model.k1,model.k2)] == sol.ps[(:k1,:k2)] == (0.25, 0.5)
        @test getp(sol, kp)(sol) == getp(sol, model.kp)(sol) == getp(sol, :kp)(sol) == 1.0
        @test getp(sol, [k1,k2])(sol) == getp(sol, [model.k1,model.k2])(sol) == getp(sol, [:k1,:k2])(sol) == [0.25, 0.5]
        @test getp(sol, (k1,k2))(sol) == getp(sol, (model.k1,model.k2))(sol) == getp(sol, (:k1,:k2))(sol) == (0.25, 0.5)
    end

    # Handles nonlinear solution differently.
    let
        sol = deepcopy(nsol)
        # Get u values.
        @test sol[X] == sol[model.X] == sol[:X]
        @test sol[XY] == sol[model.XY][1] == sol[:XY]
        @test sol[[XY,Y]] == sol[[model.XY,model.Y]] == sol[[:XY,:Y]]
        @test_broken sol[(XY,Y)] == sol[(model.XY,model.Y)] == sol[(:XY,:Y)]
        @test getu(sol, X)(sol) == getu(sol, model.X)(sol)[1] == getu(sol, :X)(sol)
        @test getu(sol, XY)(sol) == getu(sol, model.XY)(sol)[1] == getu(sol, :XY)(sol)
        @test getu(sol, [XY,Y])(sol) == getu(sol, [model.XY,model.Y])(sol) == getu(sol, [:XY,:Y])(sol)
        @test_broken getu(sol, (XY,Y))(sol) == getu(sol, (model.XY,model.Y))(sol) == getu(sol, (:XY,:Y))(sol)[1]   

        # Get p values.
        @test sol.ps[kp] == sol.ps[model.kp] == sol.ps[:kp]
        @test sol.ps[[k1,k2]] == sol.ps[[model.k1,model.k2]] == sol.ps[[:k1,:k2]]
        @test sol.ps[(k1,k2)] == sol.ps[(model.k1,model.k2)] == sol.ps[(:k1,:k2)]
        @test getp(sol, kp)(sol) == getp(sol, model.kp)(sol) == getp(sol, :kp)(sol)
        @test getp(sol, [k1,k2])(sol) == getp(sol, [model.k1,model.k2])(sol) == getp(sol, [:k1,:k2])(sol)
        @test getp(sol, (k1,k2))(sol) == getp(sol, (model.k1,model.k2))(sol) == getp(sol, (:k1,:k2))(sol)
    end
end

# Tests plotting.
let 
    @test_broken false # Currently broken for `ssol`.
    for sol in deepcopy([osol, jsol])
        # Single variable.
        @test length(plot(sol; idxs = X).series_list) == 1
        @test length(plot(sol; idxs = XY).series_list) == 1
        @test length(plot(sol; idxs = model.X).series_list) == 1
        @test length(plot(sol; idxs = model.XY).series_list) == 1
        @test length(plot(sol; idxs = :X).series_list) == 1
        @test length(plot(sol; idxs = :XY).series_list) == 1

        # As vector.
        @test length(plot(sol; idxs = [X,Y]).series_list) == 2
        @test length(plot(sol; idxs = [XY,Y]).series_list) == 2
        @test length(plot(sol; idxs = [model.X,model.Y]).series_list) == 2
        @test length(plot(sol; idxs = [model.XY,model.Y]).series_list) == 2
        @test length(plot(sol; idxs = [:X,:Y]).series_list) == 2
        @test length(plot(sol; idxs = [:XY,:Y]).series_list) == 2

        # As tuple.
        @test length(plot(sol; idxs = (X, Y)).series_list) == 1
        @test length(plot(sol; idxs = (XY, Y)).series_list) == 1
        @test length(plot(sol; idxs = (model.X, model.Y)).series_list) == 1
        @test length(plot(sol; idxs = (model.XY, model.Y)).series_list) == 1
        @test length(plot(sol; idxs = (:X, :Y)).series_list) == 1
        @test length(plot(sol; idxs = (:XY, :Y)).series_list) == 1
    end     
end


### Tests For Hierarchical System ###

# TODO