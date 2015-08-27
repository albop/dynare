// Example 1 from Collard's guide to Dynare
var h, c, y, k, a, b;
varexo e, u;

parameters beta, rho, alpha, delta, theta, psi, tau;

alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;

phi   = 0.1;

model;
c*theta*exp(h)^(1+psi)=(1-alpha)*y;
exp(k) = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
    *(exp(b(+1))*alpha*y(+1)+(1-delta)*exp(k)));
y = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
exp(k) = exp(b)*(y-c)+(1-delta)*exp(k(-1));
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
end;

steady_state_model;
a = 0;
b = 0;
k_h = ((1-beta*(1-delta))/(beta*alpha))^(1/(alpha-1));
y_h = k_h^alpha;
c_h = y_h - delta*k_h;
h = log((y_h*(1-alpha)/(c_h*theta))^(1/(1+psi)));
k = log(k_h*exp(h));
c = c_h*exp(h);
y = y_h*exp(h);
end;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

stoch_simul;
