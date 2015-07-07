@#define extended_path_version = 1

var Capital, Output, Labour, Consumption,  Investment, Efficiency, efficiency, residual, marginal_utility;

varexo EfficiencyInnovation;

parameters beta, theta, tau, alpha, psi, delta, rho, effstar, sigma;

/*
** Calibration
*/

beta    =  0.990;
theta   =  0.357;
tau     =  2.000;
alpha   =  0.450;
psi     = -0.200;
delta   =  0.020;
rho     =  0.800;
effstar =  1.000;
sigma   =  0.100;

model;

  efficiency = rho*efficiency(-1) + sigma*EfficiencyInnovation;

  Efficiency = effstar*exp(efficiency);

  [mcp = 'Investment > 0']
  -(((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption + beta*((((Consumption(+1)^theta)*((1-Labour(+1))^(1-theta)))^(1-tau))/Consumption(+1))*(alpha*((Output(+1)/Capital)^(1-psi))+1-delta);

  residual =   (((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption - beta*((((Consumption(+1)^theta)*((1-Labour(+1))^(1-theta)))^(1-tau))/Consumption(+1))*(alpha*((Output(+1)/Capital)^(1-psi))+1-delta);

  ((1-theta)/theta)*(Consumption/(1-Labour)) - (1-alpha)*(Output/Labour)^(1-psi);

  Output = Efficiency*(alpha*(Capital(-1)^psi)+(1-alpha)*(Labour^psi))^(1/psi);

  Output = Consumption + Investment;

  Investment = Capital - (1-delta)*Capital(-1);

  marginal_utility = (((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption;
end;

steady_state_model;
Efficiency = effstar;
y_k = (Efficiency^(-psi)*(1/beta-1+delta)/alpha)^(1/(1-psi));
c_k = y_k - delta;
n_k = (((y_k/Efficiency)^psi-alpha)/(1-alpha))^(1/psi);
y_n = y_k/n_k;
c_n = c_k/n_k;
Labour = y_k*(1-alpha)/(((1-theta)/theta)*c_k*(alpha*n_k^(-psi)+1-alpha)+y_k*(1-alpha));
Capital = Labour/n_k;
Consumption = c_n*Labour;
Output = Efficiency*(alpha*Capital^psi+(1-alpha)*Labour^psi)^(1/psi);
Investment = delta*Capital;
residual = 0;
marginal_utility = (((Consumption^theta)*((1-Labour)^(1-theta)))^(1-tau))/Consumption;
end;

steady;

shocks;
var EfficiencyInnovation;
periods 1;
values -4;
end;

options_.solve_algo = 10;
options_.mcp = 1;
perfect_foresight_setup(periods=100);    

perfect_foresight_solver(stack_solve_algo=7);

rplot Investment;

rplot residual;