# Bi-MADS record all dominated points
# include("DirectSearch.jl")


# using DirectSearch
export p_dim, testbi
#TODO BiMADS

"""
Potential BiMADS points on Pareto Front
"""
mutable struct B_points
    f1_cost::Float64
    f2_cost::Float64
end



function test1(x)
    f1(x) = (x[1] + 2) .^ 2 +(x[2] - 2) .^ 2- 10.0
    f2(x) = (x[1] + 2) .^ 2 +(x[2] - 2) .^ 2- 20.0
    return f1, f2
end

function ex005(x)
    return [x[1]^2 - x[2]^2; x[1] / x[2]]
end

# obj(x) = x[1]^2

p = DSProblem(
    2;
    objective = test1,
    initial_point = [0, 0],
    full_output = true,
);


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
    x=zeros(Float64,p.N)
    f1,f2=p.objective(x)
    p1::DSProblem=deepcopy(p::DSProblem)
    p2::DSProblem=deepcopy(p::DSProblem)
    SetObjective(p1,f1),SetObjective(p2,f2)
    return p1,p2
end

"""
    function phi()::Float64

Auxiliary Function from BiMADS. Aim to refomulate the Bi-obj problem into
single-obj problem. Return the cost of the refomulated objective function

x is the point to be evaluated, r is the reference point
"""
function phi(p::DSProblem, r::Vector{Float64},x::Vector{Float64})::Float64
    p1::DSProblem,p2::DSProblem=p_split(p)
    f1::Float64=p1.objective(x)
    f2::Float64=p2.objective(x)
    if (f1<=r[1]) && (f2<=r[2]) #if x dominant r
        return -(r[1]-f1)^2*(r[2]-f2)^2
    else
        return (max(f1-r[1],0))^2+(max(f2-r[2],0))^2
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
    println("BiMADS")
    p1::DSProblem,p2::DSProblem=p_split(p)


    Optimize!(p2)

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
