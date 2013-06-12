var C_1,C_2,c_11,c_12,c_21,c_22,P_1,P_2,p_1,p_2,y_1,y_2,k_1,k_2,I_1,I_2,i_11,i_12,i_21,i_22,rho_1,rho_2,P_1I,P_2I,w_1,w_2,d_1,d_2,p_1b,p_2b,R_1b,R_2b,theta_1,theta_2,chi_1,chi_2,q,W,Delta,rx;
var p_1S,p_2S,R_1S,R_2S, b_12,rx_1S,rx_2S,s_11,s_12;

varexo epsilon_1_theta,epsilon_2_theta,epsilon_1_chi,epsilon_2_chi,xi;

parameters a,phi,kappa,delta,beta,sigma,omega,a_I,phi_I,phi_L,rho_theta,rho_chi,l_1,l_2;
a = 0.85;
phi = 2;
kappa = 0.4;
delta = 0.1;
beta = 0.96;
sigma = 2;
omega = 0.5;
a_I = 0.85;
phi_I = 2;
phi_L = 2;
rho_theta = 0.75;
rho_chi = 0.79;
l_1 = 1;
l_2 = 1;

model;

// eq 1 : Price index - country 1 - 3
P_1 = (p_1^(1 - phi)*a + p_2^(1 - phi)*(1 - a))^(1/(1 - phi));


// eq 2 : Price index - country 2 - 3
P_2 = (p_1^(1 - phi)*(1 - a) + p_2^(1 - phi)*a)^(1/(1 - phi));


// eq 3
P_1 + P_2 = 1;


// eq 4
q = P_1/P_2;


// eq 5 : Production technology - country 1 - 4
y_1 = k_1^kappa*theta_1;


// eq 6 : Production technology - country 2 - 4
y_2 = k_2^kappa*theta_2;


// eq 7 : Capital law of motion - country 1 - 5
k_1 = (1 - delta)*k_1(-1) + I_1*chi_1;


// eq 8 : Capital law of motion - country 2 - 5
k_2 = (1 - delta)*k_2(-1) + I_2*chi_2;


// eq 9 : Investment price index - country 1 - 7
P_1I = (p_1^(1 - phi_I)*a_I + p_2^(1 - phi_I)*(1 - a_I))^(1/(1 - phi_I));


// eq 10 : Investment price index - country 2 - 7
P_2I = (p_1^(1 - phi_I)*(1 - a_I) + p_2^(1 - phi_I)*a_I)^(1/(1 - phi_I));


// eq 11 : Wage bill - country 1 - 8
w_1 = (1 - kappa)*p_1*y_1;


// eq 12 : Wage bill - country 2 - 8
w_2 = (1 - kappa)*p_2*y_2;


// eq 13 : Dividends - country 1 - 9
d_1 = -I_1*P_1I + (kappa)*p_1*y_1;


// eq 14 : Dividends - country 2 - 9
d_2 = -I_2*P_2I + (kappa)*p_2*y_2;


// eq 15 : FOC for investment - country 1 - 10
1 = ((1 - delta)*P_1I(1)/chi_1(1) + l_1^(1 - kappa)*k_1(1)^(-1 + kappa)*kappa*p_1(1)*theta_1(1))*chi_1*rho_1(1)/P_1I;


// eq 16 : FOC for investment - country 2 - 10
1 = ((1 - delta)*P_2I(1)/chi_2(1) + l_2^(1 - kappa)*k_2(1)^(-1 + kappa)*kappa*p_2(1)*theta_2(1))*chi_2*rho_2(1)/P_2I;


// eq 17 : Pricing kernel - country 1 - 10b
rho_1 = (C_1/C_1(-1))^(-sigma)*beta*P_1(-1)/P_1;


// eq 18 : Pricing kernel - country 2 - 10b
rho_2 = (C_2/C_2(-1))^(-sigma)*beta*P_2(-1)/P_2;


// eq 19 : FOC cost minimization (1) - country 1 - 11
i_11 = (p_1/P_1I)^(-phi_I)*a_I*I_1;


// eq 20 : FOC cost minimization (1) - country 2 - 11
i_22 = (p_2/P_2I)^(-phi_I)*a_I*I_2;


// eq 21 : FOC cost minimization (2) - country 1 - 11
i_12 = (p_2/P_1I)^(-phi_L)*(1 - a_I)*I_1;


// eq 22 : FOC cost minimization (2) - country 2 - 11
i_21 = (p_1/P_2I)^(-phi_L)*(1 - a_I)*I_2;

// eq 23
P_1*C_1 = w_1*l_1 + d_1 + W(-1)* R_1b - W + b_12(-1)*rx + s_11(-1)*rx_1S + s_12(-1)*rx_2S + xi;

// eq 24 : FOC for households (1) - country 1 - 13
c_11 = (p_1/P_1)^(-phi)*a*C_1;


// eq 25 : FOC for households (1) - country 2 - 13
c_22 = (p_2/P_2)^(-phi)*a*C_2;


// eq 26 : FOC for households (2) - country 1 - 13
c_12 = (p_2/P_1)^(-phi)*(1 - a)*C_1;


// eq 27 : FOC for households (2) - country 2 - 13
c_21 = (p_1/P_2)^(-phi)*(1 - a)*C_2;


// eq 37 : Market clearing (goods) - country 1 - 16
c_11 + c_21 + i_11 + i_21 = y_1;


// eq 38 : Market clearing (goods) - country 2 - 16
c_12 + c_22 + i_12 + i_22 = y_2;


// eq 39 : AR(1) process for theta - country 1 - 38
log(theta_1) = rho_theta*log(theta_1(-1)) + epsilon_1_theta;


// eq 40 : AR(1) process for theta - country 2 - 38
log(theta_2) = rho_theta*log(theta_2(-1)) + epsilon_2_theta;


// eq 41 : AR(1) process for chi - country 1 - 38
log(chi_1) = rho_chi*log(chi_1(-1)) + epsilon_1_chi;


// eq 42 : AR(1) process for chi - country 2 - 38
log(chi_2) = rho_chi*log(chi_2(-1)) + epsilon_2_chi;


// eq 28 : FOC for households (portfolios) - country 1 - 14
1 = R_1b(1)*rho_1(1);

// eq 29 : FOC for households (portfolios) - country 1 - 14
rho_1(1) = rho_2(1);


// eq 30 : FOC for households (portfolios) - country 2 - 14
//1 = R_1S(1)*rho_2(1);
p_1S = (p_1S(1)+d_1)*rho_2(1);
//R_1S(+1) = R_1b(+1);


// eq 31 : FOC for households (portfolios) - country 2 - 14
//1 = R_2S(1)*rho_2(1);
p_2S = (p_2S(1)+d_2)*rho_2(1);
//R_2S(+1) = R_1b(+1);


// eq 32 : FOC for households (portfolios) - country 2 - 14
//1 = R_2b(1)*rho_2(1);
p_2b = (p_2+p_2b(1))*rho_2(1);
//R_2b(+1) = R_1b(+1);

// eq 33 : Return on equities - country 1 - 15
R_1S = (d_1 + p_1S)/p_1S(-1);


// eq 34 : Return on equities - country 2 - 15
R_2S = (d_2 + p_2S)/p_2S(-1);


// eq 35 : Return on bonds - country 1 - 15
R_1b = (p_1 + p_1b)/p_1b(-1);


// eq 36 : Return on bonds - country 2 - 15
R_2b = (p_2 + p_2b)/p_2b(-1);

rx = R_2b - R_1b;

rx_1S = R_1S - R_1b;

rx_2S = R_2S - R_1b;

Delta = rho_2 - rho_1;

//[portfolio='b_12']
//Delta(+1)*rx(+1) = 0;
b_12 = 0;

[portfolio='s_11']
Delta(+1)*rx_1S(+1) = 0;
//s_11 = 0;

[portfolio='s_12']
Delta(+1)*rx_2S(+1) = 0;
//s_12 = 0;




end;

initval;
    //x_s11 = 1;
    //x_s21 = 0;
    //x_s22 = 1;
    
    rho_1=beta;
    rho_2=beta;
    R_1S=1/rho_1;
    R_2S=1/rho_2;
    R_1b=1/rho_1;
    R_2b=1/rho_2;

    
    theta_1 = 1;
    theta_2 = 1;
    
    chi_1 = 1;
    chi_2 = 1;

    p_1= 0.5;
    p_2= 0.5;
    P_1= 0.5;
    P_2= 0.5;
    P_1I= 0.5;
    P_2I= 0.5;

    k_1= ( (1/beta - (1-delta))/kappa )^(1/(kappa-1));
    y_1= k_1^kappa;
    w_1= (1-kappa)*y_1*p_1;
    I_1=delta*k_1; 
    d_1= (kappa)*y_1*p_1-I_1*P_1I;
     
    C_1= (w_1 + d_1)/P_1;

    k_2= ( (1/beta - (1-delta))/kappa )^(1/(kappa-1));
    y_2= k_1^kappa;
    w_2= (1-kappa)*y_2*p_2;
    I_2=delta*k_2;
    d_2= (kappa)*y_2*p_2-I_2*P_2I;
    
    C_2= (w_2 + d_2)/P_2;

    c_11= a*C_1;
    c_12= (1-a)*C_1;
    c_21= (1-a)*C_2;
    c_22= a*C_2;
    i_11= a_I*I_1;
    i_12= (1-a_I)*I_1;
    i_21= (1-a_I)*I_2;
    i_22= a_I*I_2;
    
    p_1S= d_1 / (R_1S - 1);
    p_2S= d_2 / (R_1S - 1);
    p_1b= p_1 / (R_1b - 1);
    p_2b= p_2 / (R_2b - 1);
    q= 1;
    W= 0;
end;

shocks;
var epsilon_1_theta = 0.001;
var epsilon_2_theta = 0.001;
var epsilon_1_chi = 0.001;
var epsilon_2_chi = 0.001;

end;

portfolios_setup;

stoch_simul(order=1);
poi;
stoch_simul(order=2);

//close all;
