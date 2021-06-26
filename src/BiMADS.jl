# Bi-MADS record all dominated points
# include("DirectSearch.jl")

using LinearAlgebra
# using DirectSearch
export p_dim, testbi
#TODO BiMADS

"""
Potential BiMADS points on Pareto Front
"""
mutable struct B_points
    cost::Vector{Float64}
    weight::Int
    function B_points(cost::Vector{Float64}, w::Int = 0)
        p = new()
        p.cost = cost
        p.weight = w
        return p
    end
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
    f1(x) = (x[1] + 2) .^ 2 + (x[2] - 2) .^ 2 - 10.0
    # f2(x) = f(x)
    f2(x) = (x[1] + 45) .^ 2 + (x[2] - 1) .^ 2 -20
    return f1, f2
end

function ex005(x)
    return [x[1]^2 - x[2]^2; x[1] / x[2]]
end

# obj(x) = x[1]^2

p = DSProblem(2; objective = test1, initial_point = [0, 0], full_output = false);


# SetVariableRange(p,1,0.,0.19)
# cons1(x) = x[1] > 0.
# AddExtremeConstraint(p, cons1)
# cons2(x) = x[1] <0.25
# AddExtremeConstraint(p, cons2)

# Optimize!(p)
# OutputIterationDetails(p)
"""
Test
"""
function testbi(p::DSProblem)
    println("recompilesss")
    OutputIterationDetails(p)
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
    SetObjective(p1, f1), SetObjective(p2, f2)
    return p1, p2
end

"""
    function phi()::Float64

Auxiliary Function from BiMADS. Aim to refomulate the Bi-obj problem into
single-obj problem. Return the cost of the refomulated objective function

x is the point to be evaluated, r is the reference point
"""
function phi(p::DSProblem, r::Vector{Float64}, x::Vector{Float64})::Float64
    p1::DSProblem, p2::DSProblem = p_split(p)
    f1::Float64 = p1.objective(x)
    f2::Float64 = p2.objective(x)
    if (f1 <= r[1]) && (f2 <= r[2]) #if x dominant r
        return -(r[1] - f1)^2 * (r[2] - f2)^2
    else
        return (max(f1 - r[1], 0))^2 + (max(f2 - r[2], 0))^2
    end
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
        push!(temp_Bpoints, B_points([p1.objective(p1.x), p2.objective(p1.x)]))

    while _check_stoppingconditions(p1)
        p1.full_output && OutputIterationDetails(p1)
        if OptimizeLoop(p1) == Dominating
            flag == 1 &&
                push!(temp_Bpoints, B_points([p1.x_cost, p2.objective(p1.x)]))
            flag == 2 &&
                push!(temp_Bpoints, B_points([p2.objective(p1.x), p1.x_cost]))
        end
    end

    return temp_Bpoints
end

"""
    pareto_front(points::Vector{B_points}, minimize = true)

find the Pareto points, only for this BiMADS.jl, not general methods
"""
function pareto_front(points::Vector{B_points}, minimize = true)
    display(points)
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
    # ref_point = [0.0, 0.0]
    # initial_point = [0.0, 0.0]
    x_L = Vector{B_points}()
    # optmization for p1
    append!(x_L, get_Bpoints(p1, p2, 1))
    x_L = reverse(x_L)
    # optmization for p2
    append!(x_L, get_Bpoints(p2, p1, 2))
    # pareto_set = pareto_front(pareto_set)

    # max_distance(p1,p2) = LinearAlgebra.norm(pareto_set[i].cost - pareto_set[i-1].cost)
    # +LinearAlgebra.norm(pareto_set[i].cost - pareto_set[i+1].cost)
    # max_index = 2






    # for i = 2:length(pareto_set)-1
    #     dis = LinearAlgebra.norm(pareto_set[i].cost - pareto_set[i-1].cost)
    #     +LinearAlgebra.norm(pareto_set[i].cost - pareto_set[i+1].cost)
    #     if dis >= max_distance
    #         max_distance = dis
    #         max_index = i
    #     end
    # end
    #
    #
    # @show pareto_set
    # @show length(pareto_set)
    # @show max_index
    # @show max_distance


    # find the point with largest gap
    # for i=1:length(all)
    # end

    return x_L
end


"""
    Optimize_Bi!(p::DSProblem)

Run the BiMADS algorithm on problem `p`.

`p` must have had its initial point and objective function set. If extreme
barrier constraints have been set then the initial point must be value for
those constraints.
"""
function Optimize_Bi!(p::DSProblem)
    println("BiMADS")
    p1::DSProblem, p2::DSProblem = p_split(p)
    # Initialization
    undominated_points = pareto_front(initial_X_L(p1, p2))
    display(undominated_points)
    # Main iteration
    δ = 0.0 #measure the coverage
    w = 0 #weight
    # ref_point, initial_point = get_ref_init(p1, p2)

    # p_re=DSProblem()
    # while 1
    #     ref_point, initial_point=get_ref_init(p1,p2)
    #     display(ref_point)
    #
    #     p.full_output && OutputIterationDetails(p_re)
    #     # OptimizeLoop(p_re)
    #     w+=1
    # end


    println("Finish")
end

Optimize_Bi!(p)


"""
    p_dim(obj::Function,p.N)::Int

The dimension of the given objective function
"""
p_dim(p::DSProblem) = p_dim(p.objective, p.N)
function p_dim(objective::Function, N::Int)::Int
    return length(objective(ones(N)))
end
