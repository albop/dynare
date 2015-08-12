# Modification of the path (for packages). Should be done in ~/.juliarc.jl with a fixed path instead.
unshift!(LOAD_PATH, abspath("../../../julia"))

# Load Dynare package
using Dynare

# Compile the rbc.mod file -> produce a module with the model definition.
Dynare.dynare("rbc.mod")
