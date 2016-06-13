# Modification of the path (for packages). Should be done in ~/.juliarc.jl with a fixed path instead.
if isempty(findin([abspath("../../../julia")], LOAD_PATH))
    unshift!(LOAD_PATH, abspath("../../../julia"))
end

# Load Dynare package
importall Dynare

# Compile the rbc.mod file -> produce a module with the model definition.

@dynare "rbc1.mod"

importall SteadyState

# First call to the steady state routine (analytical)
@time SteadyState.steady!(model_, oo_)

# First call to the steady state routine (analytical)
@time SteadyState.steady!(model_, oo_)

paramsinit = copy(model_.params);

yinit = deepcopy(oo_.steady_state)

y_init = copy(yinit)
y_init[1] = 1.1*yinit[1]
y_init[2] = 0.9*yinit[2]

# First call to the steady state routine (numerical)
println("First call to the numerical steady state routine")
@time SteadyState.steady!(model_, oo_, yinit)

# Check results
@assert maximum(abs(oo_.steady_state-yinit))<1e-9

y_init = copy(yinit)
yinit[1] = 1.1*yinit[1]
yinit[2] = 0.9*yinit[2]

# Second call to the steady state routine (numerical)
println("Second call to the numerical steady state routine")
@time SteadyState.steady!(model_, oo_, yinit)

# change alpha
println("Change α")
model_.params[4] = max(min(1.0, model_.params[4]*1.2), 0.0)
ys = SteadyState.steady(model_, oo_)
y_init = copy(yinit)
@time SteadyState.steady!(model_, oo_, y_init)
y_init = copy(yinit)
@time SteadyState.steady!(model_, oo_, y_init)
@assert maximum(abs(oo_.steady_state-ys))<1e-9

# change delta
println("Change δ")
model_.params[6] = max(min(1.0, model_.params[6]*1.2), 0.0)
ys = SteadyState.steady(model_, oo_)
y_init = copy(yinit)
@time SteadyState.steady!(model_, oo_, y_init)
y_init = copy(yinit)
@time SteadyState.steady!(model_, oo_, y_init)
@assert maximum(abs(oo_.steady_state-ys))<1e-9

# change beta
println("Change β")
model_.params[1] = max(min(1-1e-6, model_.params[1]*0.99), 0.0)
ys = SteadyState.steady(model_, oo_)
y_init = copy(yinit)
@time SteadyState.steady!(model_, oo_, y_init)
y_init = copy(yinit)
@time SteadyState.steady!(model_, oo_, y_init)
@assert maximum(abs(oo_.steady_state-ys))<1e-9

# change tau
println("Change τ")
model_.params[3] /= 1.5
ys = SteadyState.steady(model_, oo_)
y_init = copy(yinit)
@time SteadyState.steady!(model_, oo_, y_init)
y_init = copy(yinit)
@time SteadyState.steady!(model_, oo_, y_init)
@assert maximum(abs(oo_.steady_state-ys))<1e-9

# change Epsilon
println("Change ϵ")
model_.params[5] *= 1.5
ys = SteadyState.steady(model_, oo_)
y_init = copy(yinit)
@time SteadyState.steady!(model_, oo_, y_init)
y_init = copy(yinit)
@time SteadyState.steady!(model_, oo_, y_init)
@assert maximum(abs(oo_.steady_state-ys))<1e-9
