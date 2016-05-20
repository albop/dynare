# Modification of the path (for packages). Should be done in ~/.juliarc.jl with a fixed path instead.
if isempty(findin([abspath("../../../julia")], LOAD_PATH))
    unshift!(LOAD_PATH, abspath("../../../julia"))
end

# Load Dynare package
importall Dynare

# Compile the rbc.mod file -> produce a module with the model definition.


@dynare "rbc2.mod"

print(model_.fname)
