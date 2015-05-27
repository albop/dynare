var c k z_forward z_backward;
varexo x z_shock;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;

model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1); // Resource constraint
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam); // Euler equation
z_backward=0.1*1+0.3*z_backward(-1)+0.3*z_backward(-2)+0.3*z_backward(-3)+(x(-4)-1);
z_forward=0.1*1+0.45*z_forward(+1)+0.45*z_forward(+2)+(x(+4)-1);
end;

initval;
c = 1.2;
k = 12;
x = 1;
end;

histval;
x(-1)=1.30;
x(-2)=1.30;
end;

shocks;
var x;
periods 2;
values 0.9;
end;

simul(periods=200,maxit=100);