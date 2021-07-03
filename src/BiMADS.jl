using Plots
export p_dim, testbi
#TODO BiMADS

"""
points for BiMADS Algo on the BiMADS map
"""
mutable struct B_points
    cost::Vector{Float64}  #f1 and f2 vale on BiMADS map
    x_now::Vector{Float64}
    weight::Int            #weight for each undominated point to determine the δ
    function B_points(cost::Vector{Float64},x::Vector{Float64}, w::Int=0)
        p = new()
        p.cost = cost
        p.x_now=x
        p.weight = w
        return p
    end
end

function DTLZ2n2(x)

    # params
    M = 2; # Number of objectives
    n = 2; # Number of variables
    k = n - M + 1;

    # g(x)
    gx = sum((x[M:n] .- 0.5).^2);

    # functions
    ff = ones(M);
    ff1(x) = (1 + gx ) * prod(cos.(0.5 * pi * x[1:M-1]));

    ff2(x) = (1 + gx) * prod(cos.(0.5 * pi * x[1:M-2])) * sin(0.5 * pi * x[M - 2+ 1]);


    return ff1,ff2

end

function s(x)
    return abs(floor(x + 1 / 2) - x)
end

function f(x::Vector{Float64})
    Taka = 0
    w = 0.9
    for n = 1:100
        Taka -= w^n * s(2^n * x[1])
    end
    return Taka
end

function test1(x)
    f1(x) = (x[1] + 2) .^ 2  - 10.0
    # f2(x) = f(x)
    f2(x) = (x[1] -2) .^ 2 +20
    return f1, f2
end

function test2(x)
f1(x) = (x[1] + 1) .^ 2 +(x[2] - 1) .^2 -10
# f2(x) = (x[1] + 12) .^ 2 +(x[2] - 3) .^2 + 20.
f2(x) = f(x)
return f1,f2
end

function ex005(x)
    return [x[1]^2 - x[2]^2; x[1] / x[2]]
end

# p = DSProblem(2; objective = DTLZ2n2, initial_point = [0.,0.],iteration_limit=50, full_output = false);
p = DSProblem(1; objective = test1, initial_point = [0.],iteration_limit=70, full_output = false);
# p = DSProblem(2; objective = DTLZ2n2, initial_point = [0.,0.],iteration_limit=50, full_output = false);
# SetIterationLimit(p,2)

# SetFunctionEvaluationLimit(p,3)

# cons1(x) = x[1] > -1.
# AddExtremeConstraint(p, cons1)
# cons2(x) = x[1] <1.
# AddExtremeConstraint(p, cons2)

function plot_Bpoint(points::Vector{B_points})
    fig = Plots.scatter()
    for i = 1:length(points)
        fig = Plots.scatter!([points[i].cost[1]],[points[i].cost[2]],legend = false,show = true,)
    end
return fig
end

function get_adj_dis(points::Vector{B_points})
    println("-----------------")
    for i=1:length(points)-1
        display(LinearAlgebra.norm(points[i+1].cost - points[i].cost)^2)
    end
    println("-----------------")
end

"""
    p_split(p::DSProblem)::Tuple{DSProblem,DSProblem}

To split Bi-obj problem into two single obj problems
"""
function p_split(p::DSProblem)::Tuple{DSProblem,DSProblem}
    x = zeros(Float64, p.N)
    f1, f2 = p.objective(x)
    p1::DSProblem = deepcopy(p::DSProblem)
    p2::DSProblem = deepcopy(p::DSProblem)
    SetObjective(p1, f1)
    SetObjective(p2, f2)
    return p1, p2
end

"""
    get_Bpoints(p1::DSProblem,p2::DSProblem,flag::Int)

Get the points to determine the reference point and Pareto front
flag for identifying optimize p1 or p2
"""
function get_Bpoints(p1::DSProblem, p2::DSProblem, flag::Int)
    temp_Bpoints = Vector{B_points}()
    Setup(p1)
    # push B_point for initial point
    flag == 1 &&
        push!(temp_Bpoints, B_points([p1.objective(p1.x), p2.objective(p1.x)],p1.x))

    while _check_stoppingconditions(p1)
        p1.full_output && OutputIterationDetails(p1)
        if OptimizeLoop(p1) == Dominating
            flag == 1 &&
                push!(temp_Bpoints, B_points([p1.x_cost, p2.objective(p1.x)],p1.x))
            flag == 2 &&
                push!(temp_Bpoints, B_points([p2.objective(p1.x), p1.x_cost],p1.x))
        end
    end

    Finish(p1)

    if p1.status.optimization_status != MeshPrecisionLimit && p1.status.optimization_status != PollPrecisionLimit
        error("The $(p1.status.optimization_status_string) is not enough for BiMADS, please increase it")
    end


# display(plot_Bpoint(temp_Bpoints))
    return temp_Bpoints
end

"""
    pareto_front(points::Vector{B_points}, minimize = true)

find the Pareto points, only for this BiMADS.jl, not general methods
"""
function pareto_front(points::Vector{B_points}, minimize = true)
    cmp_strict = minimize ? (<) : (>)
    cmp = minimize ? (<=) : (>=)
    dims = 2 #number of objective functions, needs modification for multi-MADS
    strictly_dominates(u::B_points, v::B_points) =
        all(cmp(u.cost[i], v.cost[i]) for i = 1:dims) &&
        any(cmp_strict(u.cost[i], v.cost[i]) for i = 1:dims)
    undominated(p::B_points) =
        !any(
            strictly_dominates(p2::B_points, p::B_points) for
            p2 ∈ points if p2.cost != p.cost
        )
    filter!(undominated, points)
    sort!(points, by = v -> v.cost, rev = false)
    return points
end

"""
    initial_X_L(p1, p2)::Vector{Float64}

get the initial undominated points from p1 and p2
"""
function initial_X_L(p1, p2)::Vector{B_points}
    x_L = Vector{B_points}()
    # optmization for p1
    append!(x_L, get_Bpoints(p1, p2, 1))
    x_L = reverse(x_L)
    # optmization for p2
    append!(x_L, get_Bpoints(p2, p1, 2))


    # display(x_L)
    # get_adj_dis(x_L)
    # display(plot_Bpoint(x_L))
    return x_L
end

"""
    get_ref(points::Vector{B_points})::Int

get the reference point for L>2
"""
function get_ref(points::Vector{B_points})::Int
    dis(x::Int) = LinearAlgebra.norm(points[x].cost - points[x-1].cost)^2+LinearAlgebra.norm(points[x].cost - points[x+1].cost)^2
    max_index = 2 #point with largest distance to its adjacent points
    max_distance=0.

    for i = 2:length(points)-1
        if dis(i)/(points[i].weight+1) >= max_distance
            max_distance = dis(i)/(points[i].weight+1)
            max_index = i
        end
    end
    return max_index
end


"""
    ReferencePointDetermination(undominated_points::Vector{B_points})::Tuple{Int,Float64,Vector{Float64}}

get the reference point for any case
"""
function ReferencePointDetermination(undominated_points::Vector{B_points})::Tuple{Int,Float64,Vector{Float64}}
    L=length(undominated_points)
    if L==1
        j=1
        ref_point=undominated_points[1].cost
        δ= 1.
    elseif L==2
        j=1 #different definition in paper(j=2) and books(j=1)
        ref_point=[undominated_points[2].cost[1],undominated_points[1].cost[2]]
        δ=(LinearAlgebra.norm(undominated_points[1].cost - undominated_points[2].cost)^2)/(undominated_points[2].weight+1)
    else
        j = get_ref(undominated_points)
        ref_point=[undominated_points[j+1].cost[1],undominated_points[j-1].cost[2]]
        δ=(LinearAlgebra.norm(undominated_points[j].cost - undominated_points[j-1].cost)^2
        +LinearAlgebra.norm(undominated_points[j].cost - undominated_points[j+1].cost)^2)/(undominated_points[j].weight+1)
    end
    undominated_points[j].weight+=1
    # get_adj_dis(undominated_points)
    return j, δ, ref_point
end

"""
    function phi(f1::Function,f2::Function, r::Vector{Float64}, x::Vector{Float64})

Auxiliary Function from BiMADS. Aim to refomulate the Bi-obj problem into
single-obj problem. Return the cost of the refomulated objective function

x is the point to be evaluated, r is the reference point
"""
function phi(f1::Function,f2::Function, r::Vector{Float64}, x::Vector{Float64})
        # @show x
    if (f1(x) <= r[1]) && (f2(x) <= r[2]) #if x dominant r
        # @show -(r[1] - f1(x))^2 * (r[2] - f2(x))^2
        return -(r[1] - f1(x))^2 * (r[2] - f2(x))^2
    else
        # @show (max(f1(x) - r[1], 0))^2 + (max(f2(x) - r[2], 0))^2
        return (max(f1(x) - r[1], 0))^2 + (max(f2(x) - r[2], 0))^2
    end
end

"""
    update_p(p::DSProblem)
reset the reformulated problem for new MADS iteration
Including reset the Mesh/Poll size, IterationLimit, FunctionEvaluationLimit, RuntimeLimit, etc.

"""
function update_p(p_init::DSProblem,p_reform::DSProblem)::DSProblem
    p_reform=deepcopy(p_init)


end

function update_p(p,p1::DSProblem,p2::DSProblem)
    for c in p.stoppingconditions
        i=_get_conditionindexes(p,typeof(c))[1]
        if c isa IterationStoppingCondition
            p.stoppingconditions[i].limit-=p1.status.iteration+p2.status.iteration
            p.stoppingconditions[i].limit<1 && error("The Iteration limit is not enough for BiMADS, please increase it")
        elseif s isa FunctionEvaluationStoppingCondition
            p.stoppingconditions[i].limit-=p1.status.function_evaluations+p2.status.function_evaluations
            p.stoppingconditions[i].limit<1 && error("The Function evaluation limit is not enough for BiMADS, please increase it")
        elseif s isa RuntimeStoppingCondition
            p.stoppingconditions[i].limit -=p1.status.runtime_total+p2.status.runtime_total
            p.stoppingconditions[i].limit<=0 && error("The Runtime limit is not enough for BiMADS, please increase it")
        end
    end



end


"""
    Optimize_Bi!(p::DSProblem)

Run the BiMADS algorithm on problem `p`.

`p` must have had its initial point and objective function set. If extreme
barrier constraints have been set then the initial point must be value for
those constraints.
"""
function Optimize_Bi!(p::DSProblem)
    # term = REPL.Terminals.TTYTerminal("xterm",stdin,stdout,stderr)
    # REPL.Terminals.raw!(term,true)
    # Base.start_reading(stdin)
    println("BiMADS")
    p1::DSProblem, p2::DSProblem = p_split(p)
    f1 = p1.objective
    f2 = p2.objective
    p_reform=DSProblem(1)
    # Initialization
    undominated_points = pareto_front(initial_X_L(p1, p2))
    # display(undominated_points)
# display(plot_Bpoint(undominated_points))
    ref_point=Vector{Float64}()

    update_p(p,p1,p2)

    n_iteration=0
    # n_evaluation=p1.status.function_evaluations+p2.status.function_evaluations,
    # totalTime=p1.status.runtime_total+p2.status.runtime_total]

    j = 0     #initial point for the refomulated sigle-obj problem
    δ = 0.0 # quantify the quality of the coverage of the undominated points
            # here is the coverage of the point with largest distance to its adjacent points

    # Main iteration
    while true
        tt=0
        p_reform=update_p(p,p_reform)
        println("Start of iteration",n_iteration)
        #Reference point determination
        j, δ, ref_point=ReferencePointDetermination(undominated_points)
        println("j=",j,"  ref=",ref_point)
        #Single-objective formulation minimization
        f_reform(x)=phi(f1,f2, ref_point, x)
        SetObjective(p_reform, f_reform)
        SetInitialPoint(p_reform,undominated_points[j].x_now)

# fig1=plot_Bpoint(undominated_points)
# p_reform.full_output=false
        #run MADS for refomulated problem and add new undominated points
        Setup(p_reform)
        while _check_stoppingconditions(p_reform)
            p_reform.full_output && OutputIterationDetails(p_reform)
            if OptimizeLoop(p_reform) == Dominating
                    push!(undominated_points, B_points([p1.objective(p_reform.x), p2.objective(p_reform.x)],p_reform.x))
                    tt+=1
            end
        end
# fig2=plot_Bpoint(undominated_points)
println("Add new points:",tt)
# display(undominated_points)

        #update X_L
        pareto_front(undominated_points)
        println("===============================")
# length(undominated_points)>=474 && break
        n_iteration>=p.stoppingconditions[1].limit ? (break) : (n_iteration+=1)
        # display(undominated_points)
        # fig2=plot_Bpoint(undominated_points)
        # sort!(undominated_points, by = v -> v.cost, rev = false)
        # display(plot(fig1,fig2,fig3))
# sleep(2)
# check()==1 && break
# if !isempty(input)
#     println("dd")
#     display(input)
#     break
# end
    end
println("Total undominated points:", length(undominated_points))
        return undominated_points
end

@time result=Optimize_Bi!(p)

# fig=plot_Bpoint(result)
fig=scatter()
for i in 1:length(result)
    fig=scatter!([result[i].cost[1]],[result[i].cost[2]],legend = false)
end
display(fig)
# savefig(fig, "/Users/zyy/Desktop/XJTLU/MSc_Project/Julia/test_julia/Results/DS_result_$(p.stoppingconditions[1].limit).pdf");





"""
    p_dim(obj::Function,p.N)::Int

The dimension of the given objective function
"""
p_dim(p::DSProblem) = p_dim(p.objective, p.N)
function p_dim(objective::Function, N::Int)::Int
    return length(objective(ones(N)))
end
