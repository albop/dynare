// Generated by Dolo
// Model basename : two_risky_countries


var D,rI_1,C_2,C_1,I_2,rho_1,d_1,rI_2,alpha_1_1,d_2,rx_2,pS_1,Y_1,W_1,rS_2,rx_1,pS_2,K_2,rho_2,Pf,I_1,chi_2,K_1,w_2,alpha_1_2,Y_2,Rf,w_1,chi_1,rS_1;

varexo eps_1,zeta_1,eps_2,zeta_2,xi;

parameters rho_chi,beta,theta,gam,delta;
rho_chi = 0.5;
beta = 0.96;
gam = 4;
theta = 0.5;
delta = 0.1;

model;


// eq 1 : Production technology
[ name = 'Production technology' ]
Y_1 = K_1(-1)^theta*(1 + eps_1);


// eq 2 : Capital law of motion
[ name = 'Capital law of motion' ]
K_1 = (1 - delta)*K_1(-1) + I_1*chi_1;


// eq 3 : Shock on investment efficiency
[ name = 'Shock on investment efficiency' ]
-1 + chi_1 = -(1 - chi_1(-1))*rho_chi + zeta_1;


// eq 4 : Dividends
[ name = 'Dividends' ]
d_1 = -I_1 + theta*Y_1;


// eq 5 : Wages
[ name = 'Wages' ]
w_1 = (1 - theta)*Y_1;


// eq 6 : Stochastic discount factor
[ name = 'Stochastic discount factor' ]
rho_1 = C_1(-1)^gam*C_1^(-gam)*beta;


// eq 7 : Return on investment
[ name = 'Return on investment' ]
rI_1 = 1 - delta + K_1(-1)^(-1 + theta)*(1 + eps_1)*theta;


// eq 8 : Return on stocks
[ name = 'Return on stocks' ]
rS_1 = d_1/pS_1(-1);


// eq 9 : Production technology
[ name = 'Production technology' ]
Y_2 = K_2(-1)^theta*(1 + eps_2);


// eq 10 : Capital law of motion
[ name = 'Capital law of motion' ]
K_2 = (1 - delta)*K_2(-1) + I_2*chi_2;


// eq 11 : Shock on investment efficiency
[ name = 'Shock on investment efficiency' ]
-1 + chi_2 = -(1 - chi_2(-1))*rho_chi + zeta_2;


// eq 12 : Dividends
[ name = 'Dividends' ]
d_2 = -I_2 + theta*Y_2;


// eq 13 : Wages
[ name = 'Wages' ]
w_2 = (1 - theta)*Y_2;


// eq 14 : Stochastic discount factor
[ name = 'Stochastic discount factor' ]
rho_2 = C_2(-1)^gam*C_2^(-gam)*beta;


// eq 15 : Return on investment
[ name = 'Return on investment' ]
rI_2 = 1 - delta + K_2(-1)^(-1 + theta)*(1 + eps_2)*theta;


// eq 16 : Return on stocks
[ name = 'Return on stocks' ]
rS_2 = d_2/pS_2(-1);


// eq 17 : Risk free return
[ name = 'Risk free return' ]
Rf = 1/Pf(-1);


// eq 18 : Budget constraint - country 1
[ name = 'Budget constraint - country 1' ]
C_1 = -W_1 + (-Rf + rS_1)*alpha_1_1(-1) + (-Rf + rS_2)*alpha_1_2(-1) + Rf*W_1(-1) + d_1 + w_1 + xi;


// eq 19 : Market clearing condition for goods
[ name = 'Market clearing condition for goods' ]
C_1 + C_2 = -I_1 - I_2 + Y_1 + Y_2;


// eq 20 : Excess return - 1
[ name = 'Excess return - 1' ]
rx_1 = -Rf + rS_1;


// eq 21 : Excess return - 2
[ name = 'Excess return - 2' ]
rx_2 = -Rf + rS_2;


// eq 22 : Global stochastic discount factor
[ name = 'Global stochastic discount factor' ]
D = -rho_1 + rho_2;


// eq 23 : Euler condition for investment - country 1
[ name = 'Euler condition for investment - country 1' ]
rI_1(1)*rho_1(1) = 1;


// eq 24 : Euler condition for investment - country 2
[ name = 'Euler condition for investment - country 2' ]
rI_2(1)*rho_2(1) = 1;


// eq 25 : Euler condition for riskfree asset
[ name = 'Euler condition for riskfree asset' ]
rho_1(1) = Pf;


// eq 26 : Euler condition for riskfree asset
[ name = 'Euler condition for riskfree asset' ]
rho_2(1) = Pf;


// eq 27 : Asset pricing equation - stock 1
[ name = 'Asset pricing equation - stock 1' ]
d_1(1)/pS_1 = 1/Pf;


// eq 28 : Asset pricing equation - stock 2
[ name = 'Asset pricing equation - stock 2' ]
d_2(1)/pS_2 = 1/Pf;


// eq 29 : Portolio equation
[ portfolio = 'alpha_1_1' , name = 'Portolio equation' ]
D(1)*rx_1(1) = 0;
//alpha_1_1 = 0;


// eq 30 : Portfolio equation
[ portfolio = 'alpha_1_2' , name = 'Portfolio equation' ]
D(1)*rx_2(1) = 0;
//alpha_1_2 = 0;
end;

initval;
eps_2 = 0;
W_1 = 0;
rS_2 = 1/beta;
Pf = beta;
eps_1 = 0;
chi_1 = 1;
rI_1 = 1/beta;
rho_1 = beta;
rI_2 = 1/beta;
rho_2 = beta;
chi_2 = 1;
Rf = 1/beta;
rS_1 = 1/beta;
K_2 = (-(1 - delta - rS_2)/theta)^(-1/(1 - theta));
K_1 = (-(1 - delta - rS_1)/theta)^(-1/(1 - theta));
Y_2 = K_2^theta;
I_2 = delta*K_2;
d_2 = -I_2 + theta*Y_2;
Y_1 = K_1^theta;
pS_2 = beta*d_2;
I_1 = delta*K_1;
C_2 = -I_2 + Y_2;
C_1 = -I_1 + Y_1;
w_2 = (1 - theta)*Y_2;
w_1 = (1 - theta)*Y_1;
d_1 = -I_1 + theta*Y_1;
pS_1 = beta*d_1;
end;

shocks;
var eps_1 = 0.0100000000000000 ;
var zeta_1 = 0.0100000000000000 ;
var eps_2 = 0.0100000000000000 ;
var zeta_2 = 0.0100000000000000 ;
end;


portfolios_setup('method', 'devereux-sutherland');

stoch_simul(order=2);
