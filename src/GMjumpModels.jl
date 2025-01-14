# ==============================================================================
# Modele JuMP pour calculer la relaxation linéaire du 2SPA sur un objectif donne, activant eventuellement une ϵ-contrainte

function computeLinearRelax2SPA(  nbvar::Int,
                                  nbctr::Int,
                                  A::Array{Int,2},
                                  c1::Array{Int,1},
                                  c2::Array{Int,1},
                                  epsilon,
                                  obj::Int
)
  model = Model(default_optimizer)
  set_silent(model)
  @variable(model, 0.0 <= x[1:nbvar] <= 1.0 )
  @constraint(model, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
  if obj == 1
    @objective(model, Min, sum((c1[i])*x[i] for i in 1:nbvar))
    @constraint(model, sum((c2[i])*x[i] for i in 1:nbvar) <= epsilon)
  else
    @objective(model, Min, sum((c2[i])*x[i] for i in 1:nbvar))
    @constraint(model, sum((c1[i])*x[i] for i in 1:nbvar) <= epsilon)
  end
  set_attribute(model, "primal_feasibility_tolerance", 1.)
  optimize!(model)
  return objective_value(model), value.(x)
end

#= Compute the Linear Relaxation on 2SPA, fixing one objective function. Uses an ϵ-constraint enventually. This version forces post-optimization variables
being integer. The user must be aware of possible "communicating vessels" between a posteriori non-integer and integer variables. More precisely, previously non-integer
variables may become non-integer variables after the second call of optimizer.
 =#
function computeLinearRelax2SPAInt(  nbvar::Int, 
  nbctr::Int,
  A::Array{Int,2},
  c1::Array{Int,1},
  c2::Array{Int,1},
  epsilon,
  obj::Int
)
  model = Model(default_optimizer)
  set_silent(model)
  @variable(model, 0.0 <= x[1:nbvar] <= 1.0 )
  @constraint(model, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
  if obj == 1
    @objective(model, Min, sum((c1[i])*x[i] for i in 1:nbvar))
    @constraint(model, sum((c2[i])*x[i] for i in 1:nbvar) <= epsilon)
  else
    @objective(model, Min, sum((c2[i])*x[i] for i in 1:nbvar))
    @constraint(model, sum((c1[i])*x[i] for i in 1:nbvar) <= epsilon)
  end
  optimize!(model)

  # Forcing of non-integer variables to be integer
  xOutput::Vector{Float64} = value.(x)
  [set_binary(model[:x][i]) for i in eachindex(xOutput) if !(isapprox(xOutput[i],0.,atol=1e-3)||isapprox(xOutput[i],1.,atol=1e-3))] # set to binary non integer variables
  optimize!(model) # re-optimization

  return objective_value(model), value.(x)
end