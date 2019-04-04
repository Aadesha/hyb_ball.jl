using Dates; now()
using Revise, MathematicalSystems,Reachability, LinearAlgebra, HybridSystems
using LazySets; Interval(0, 1)
using Reachability: solve
using Plots
using TaylorIntegration


@taylorize function bball_up!(t, x, dx)
    dx[1] = x[2]
    dx[2] = -1.0 - 0.001*(x[1])^2
    return dx
end

@taylorize function bball_down!(t, x, dx)
    dx[1] = x[2]
    dx[2] = -1.0 + 0.001*(x[1])^2
    return dx
end

function bouncing_ball(;use_CBBCS=false)

    automaton = LightAutomaton(2) # two nodes

    inv_up = HPolyhedron([HalfSpace([-1.0, 0.0], 0.0),  # x >= 0 #I guess to combine two HalfSpace this should work
                          HalfSpace([0.0, -1.0], 0.0)]) # v >= 0

    inv_down = HPolyhedron([HalfSpace([-1.0, 0.0], 0.0),  # x >= 0
                         HalfSpace([ 0.0, 1.0], 0.0)])  # v <= 0

    m1 = ConstrainedBlackBoxContinuousSystem(bball_up!, 2, inv_up) #OR ConstrainedBlackBox'Control'System ? or something else?

    m2 = ConstrainedBlackBoxContinuousSystem(bball_down!,2, inv_down)

    modes = [m1,m2]  #modes

    add_transition!(automaton, 1, 2, 1)    #alpha transition
    add_transition!(automaton, 2, 1, 2)    #beta transition



    guard_alpha = HPolyhedron([HalfSpace([1.0, 0.0], 0.0),   #x>=0
                              HalfSpace([-1.0, 0.0], 0.0)])  #x<=0
    guard_beta =  HPolyhedron([HalfSpace([0.0, 1.0], 0.0),   #v<=0
                               HalfSpace([0.0,-1.0], 0.0)])  #v>=0


# what should be the value and dimension of A to map on gaurd ???
    A = [1.0 0.0; 0.0 -0.75]
    t1 = ConstrainedLinearMap(A, guard_alpha)
    t2 = ConstrainedLinearMap(A, guard_beta)

    resetmaps = [t1,t2] #resetmaps

    # switching
    switching = AutonomousSwitching()
    switchings = fill(switching, 2)

    ℋ = HybridSystem(automaton, modes, resetmaps, switchings)

#initial condition in mode one
    X0 = Hyperrectangle(low=[10.0, 0.0], high=[10.2, 0.0])#intial condition in mode two, is it required?(I dont think so)
    initial_condition = [(1,X0)]

    system = InitialValueProblem(ℋ, initial_condition)

    options = Options(:mode=>"reach", :T=>6.0, :plot_vars=>[1, 2], :project_reachset=>false)
    
    return (system, options)
end


problem, options = bouncing_ball(use_CBBCS=true);

options = Options(:mode=>"reach", :T=>12.0, :plot_vars=>[1, 2])

@time sol_TMJets = solve(problem, options, TMJets(:orderT=>5, :orderQ=>2, :abs_tol=>1e-10),LazyDiscretePost(:check_invariant_intersection=>true));


plot(sol_TMJets, use_subindices=false, aspectratio=1, alpha=.5, color=:lightblue)
