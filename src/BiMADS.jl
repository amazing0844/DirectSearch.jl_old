# Bi-MADS record all dominated points
include("DirectSearch.jl")


# using DirectSearch
# export
#TODO BiMADS

# title = "Iteration test #$(p.status.iteration):"

function test1(x)
f1 = (x[1]+2).^2 - 10.
# f2 = (x[1]-2).^2 + 20.
f2=1/(1-x[1])
return f1,f2
end

function ex005(x)
    return [x[1]^2 - x[2]^2;
            x[1] / x[2]]
end

obj(x) = x[1]^2 - x[2]^2

p = DSProblem(1; objective=obj,initial_point=[0.17],full_output=true);
# SetVariableRange(p,1,0.,0.19)
# cons1(x) = x[1] > 0.
# AddExtremeConstraint(p, cons1)
# cons2(x) = x[1] <0.25
# AddExtremeConstraint(p, cons2)

Optimize!(p)





"""
    p_dim(obj::Function)::Int

The dimension of the given objective function
"""
p_dim(p::DSProblem) = p_dim(p.objective)
function p_dim(objective::Function)::Int
return length(objective(1.))
end
