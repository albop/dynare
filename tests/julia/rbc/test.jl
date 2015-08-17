#clear workspace
workspace()

# Modification of the path (for packages). Should be done in ~/.juliarc.jl with a fixed path instead.
unshift!(LOAD_PATH, abspath("../../../julia"))

# Load Dynare package
importall Dynare

# Compile the rbc.mod file -> produce a module with the model definition.
@dynare "rbc.mod"
