//Rudebusch Svensson (1999) model with optimized policy rule

var pi y i;
varexo epsilon eta;

parameters alpha_pi1 alpha_pi2 alpha_pi3 alpha_pi4 alpha_y
           beta_y1 beta_y2 beta_r
           gamma_pi1 gamma_pi2 gamma_pi3 gamma_pi4
           gamma_y1 gamma_y2 gamma_i1 gamma_i2 gamma_i3; 


alpha_pi1 = 0.7;
alpha_pi2 = -0.1;
alpha_pi3 = 0.28;
alpha_pi4 = 0.12;
alpha_y   = 0.14;
beta_y1   = 1.16;
beta_y2   = -0.25;
beta_r    = 0.1;
gamma_pi1 = 0.88;
gamma_pi2 = 0.30;
gamma_pi3 = 0.38;
gamma_pi4 = 0.13;
gamma_y1  = 1.30;
gamma_y2  = -0.33;
gamma_i1  = 0.47;
gamma_i2  = -0.06;
gamma_i3  = -0.03;




model(linear);
pi = alpha_pi1*pi(-1) +alpha_pi2*pi(-2) +alpha_pi3*pi(-3) +alpha_pi4*pi(-4) + alpha_y*y(-1) + epsilon;
y  = beta_y1*y(-1) + beta_y2*y(-2) - beta_r*((i(-1)+i(-2)+i(-3)+i(-4))/4 - (pi(-1)+pi(-2)+pi(-3)+pi(-4))/4) + eta;
i =  gamma_pi1*pi + gamma_pi2*pi(-1) + gamma_pi3*pi(-2) + gamma_pi4*pi(-3)
     + gamma_y1*y + gamma_y2*y(-1) + gamma_i1*i(-1) + gamma_i2*i(-2) + gamma_i3*i(-3);
end;          

stoch_simul(irf=0);

A = diag(ones(12,1));
A(3,1) = -0.88;
A(3,2) = -1.30;
B = [ 0.7   0.14  0     -0.1    0    0     0.28  0  0     0.12  0 0;
      0.025 1.16 -0.025 0.025 -0.25 -0.025 0.025 0 -0.025 0.025 0 -0.025;
      0.30  -0.33 0.47  0.38    0   -0.06  0.13  0 -0.03  0     0 0;
      1 0 0 0 0 0 0 0 0 0 0 0;
      0 1 0 0 0 0 0 0 0 0 0 0;
      0 0 1 0 0 0 0 0 0 0 0 0;
      0 0 0 1 0 0 0 0 0 0 0 0;
      0 0 0 0 1 0 0 0 0 0 0 0;
      0 0 0 0 0 1 0 0 0 0 0 0;
      0 0 0 0 0 0 1 0 0 0 0 0;
      0 0 0 0 0 0 0 1 0 0 0 0;
      0 0 0 0 0 0 0 0 1 0 0 0];

eigenvalues = sort(abs(eig(A\B)));
if (max(sort(abs(oo_.dr.eigval)) - eigenvalues(3:end))) > 1e-13;
   error(sprintf('Eigenvalues aren''t correct. Absolute maximum error: %f',max(abs(oo_.dr.eigval - eig(A)))));
end;

