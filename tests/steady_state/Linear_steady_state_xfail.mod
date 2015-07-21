// Test whether a nonlinear equation specfied as a linear model is correctly filtered out
var A;

varexo epsilona;
parameters  rho;

rho =   .42; 

model(linear);
log(A) = rho*log(A(-1)) + epsilona;
end;

initval;
    A = 1.34;
end;

steady;






