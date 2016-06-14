To compile the Dynare mod files ```rbc1.mod``` and/or ```rbc2.mod``` and produce a julia module, just do

```
include("test1.jl")
```
```
include("test2.jl")
```
or
```
include("test3.jl")
```
in a script or in julia's shell. The first script, evaluates the steady state of the model, using analytical solution or a numerical solver. The two other scripts do nothing except compiling the mod files. The ```test1pfm.jl``` script shows how to solve a perfect foresight model (transition to the steady state or effects of an expected shock).

Note that Julia's packages ```NLsolve.jl``` and ```PyPlot.jl``` are required.
	
