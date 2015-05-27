// This file deals with the resolution and estimation of a basic DSGE model with
//employment for comparison with the benchmark in Gauss which solves with
//the same particular filter but global methodology.
//
// Juin 2015

var k A c l i y;
varexo e_a;

parameters alp bet tet tau delt rho ;
alp = 0.4;
bet = 0.99;
tet = 0.357 ;
tau =  50 ;
delt = 0.02;
rho = 0.95;

model;
c = ((1 - alp)*tet/(1-tet))*A*(1-l)*((k(-1)/l)^alp) ;
y = A*(k(-1)^alp)*(l^(1-alp)) ;
i = y-c ;
k = (1-delt)*k(-1) + i ;
log(A) = rho*log(A(-1)) + e_a ;
(((c^(tet))*((1-l)^(1-tet)))^(1-tau))/c - bet*((((c(+1)^(tet))*((1-l(+1))^(1-tet)))^(1-tau))/c(+1))*(1 -delt+alp*(A(1)*(k^alp)*(l(1)^(1-alp)))/k)=0 ;
end;

shocks;
var e_a; stderr 0.035;
end;

steady;

estimated_params;
alp, uniform_pdf,,, 0.0001, 0.99;
bet, uniform_pdf,,, 0.0001, 0.99999;
tet, uniform_pdf,,, 0.0001, .999;
tau, uniform_pdf,,, 0.0001, 100;
delt, uniform_pdf,,, 0.0001, 0.05;
rho, uniform_pdf,,, 0.0001, 0.9999;
stderr e_a, uniform_pdf,,, 0.00001, 0.1;
stderr y, uniform_pdf,,, 0.00001, 0.1;
stderr l, uniform_pdf,,, 0.00001, 0.1;
stderr i, uniform_pdf,,, 0.00001, 0.1;
end;

// Risky
estimated_params_init;
alp, 0.4;   
bet, 0.99;
tet, 0.357;
tau, 50;
delt, 0.02;
rho, 0.95;
stderr e_a, .035;
stderr y, .00158;
stderr l, .0011;
stderr i, .000866;
end;

estimated_params_init;
alp, 0.405156;   
bet, 0.9999;
tet, 0.360716;
tau, 30.687089;
delt, 0.017113;
rho, 0.968120;
stderr e_a, .043001;
stderr y, .006791;
stderr l, .001065;
stderr i, .004716;
end;

estimated_params_init;
stderr e_a, 0.043021303067030 ;
stderr y, 0.006670476313281;
stderr l, 0.001074088088005;
stderr i, 0.004799910878108;
alp, 0.405140190096526;
bet, 0.9953935267771;
tet, 0.360731313440020;
tau, 30.687088931783695 ;
delt, 0.017070312172381 ;
rho, 0.968261076182835 ;
end;

// Extreme 
estimated_params_init;
alp, 0.4;   
bet, 0.99;
tet, 0.357;
tau, 50;
delt, 0.02;
rho, 0.95;
stderr e_a, .035;
stderr y, .0175;
stderr l, .00312;
stderr i, .00465;
end;

// Extreme: intermediate solution for SIR
//estimated_params_init;
//stderr e_a, 0.025533519509666;
//stderr y, 0.018507739054788;
//stderr l, 0.003197047316401;
//stderr i, 0.003802850762698;
//alp, 0.363888235523713;
//bet, 0.994999761247253;
//tet, 0.342346704799396;
//tau, 53.199696552640710;
//delt, 0.012206251968057;
//rho, 0.976115559512299 ;
//end ; 

varobs y l i ;

options_.mode_check.neighbourhood_size = .01 ;
options_.mode_check.number_of_points = 250;

//KALMAN
//estimation(datafile=extreme,nograph,order=1,nobs=150,mode_compute=8,mh_replic=0,mode_check);

//SIR 
//estimation(datafile=extreme,order=2,nograph,number_of_particles=1000,nobs=150,mh_replic=0,mode_compute=8,mode_check);

// SISmoothR
//estimation(datafile=extreme,order=2,nograph,number_of_particles=1000,nobs=150,resampling_method=smooth,mode_compute=4,mh_replic=0);
//estimation(datafile=extreme,order=2,nograph,number_of_particles=1000,nobs=150,resampling_method=smooth,mode_compute=8,mode_file=dsge_base2_mode,mh_replic=0,mode_check);
//estimation(datafile=extreme,order=2,nograph,number_of_particles=1000,nobs=150,resampling_method=smooth,mode_compute=4,mode_file=dsge_base2_mode,mh_replic=0,mode_check);

// APF
//estimation(datafile=extreme,order=2,nograph,filter_algorithm=apf,number_of_particles=1000,nobs=150,mh_replic=0,mode_compute=8,mode_check);

// GPF
//estimation(datafile=extreme,order=2,nograph,filter_algorithm=gf,distribution_approximation=montecarlo,number_of_particles=1000,nobs=150,mh_replic=0,mode_compute=8);
//estimation(datafile=extreme,order=2,nograph,filter_algorithm=gf,distribution_approximation=montecarlo,number_of_particles=1000,nobs=150,mode_file=dsge_base2_mode,mh_replic=0,mode_compute=4,mode_check);
//estimation(datafile=risky,order=2,nograph,filter_algorithm=gf,distribution_approximation=montecarlo,number_of_particles=1000,nobs=150,mh_replic=0,mode_compute=8);
//estimation(datafile=risky,order=2,nograph,filter_algorithm=gf,distribution_approximation=montecarlo,number_of_particles=1000,nobs=150,mh_replic=0,mode_compute=4,mode_check);

// GCF
//estimation(datafile=extreme,order=2,nograph,filter_algorithm=gf,nobs=150,mh_replic=0,mode_compute=4);
//estimation(datafile=extreme,order=2,nograph,filter_algorithm=gf,nobs=150,mh_replic=0,mode_compute=8,mode_file=dsge_base2_mode,mode_check);

// GUF
//estimation(datafile=extreme,order=2,nograph,filter_algorithm=gf,proposal_approximation=unscented,distribution_approximation=unscented,nobs=150,mh_replic=0,mode_compute=4);
//estimation(datafile=extreme,order=2,nograph,filter_algorithm=gf,proposal_approximation=unscented,distribution_approximation=unscented,nobs=150,mh_replic=0,mode_compute=8,mode_check);

// GMPF
//estimation(datafile=extreme,nograph,order=2,filter_algorithm=gmf,distribution_approximation=montecarlo,number_of_particles=1000,nobs=150,mh_replic=0,mode_compute=8);
//estimation(datafile=extreme,nograph,order=2,filter_algorithm=gmf,distribution_approximation=montecarlo,number_of_particles=1000,nobs=150,mh_replic=0,mode_file=dsge_base2_mode,mode_compute=8);
//estimation(datafile=extreme,nograph,order=2,filter_algorithm=gmf,distribution_approximation=montecarlo,number_of_particles=1000,nobs=150,mh_replic=0,mode_file=dsge_base2_mode,mode_compute=4,mode_check);

// GMCF
//estimation(datafile=extreme,nograph,order=2,filter_algorithm=gmf,nobs=150,mh_replic=0,mode_file=dsge_base2_mode,mode_compute=8);
//estimation(datafile=extreme,nograph,order=2,filter_algorithm=gmf,nobs=150,mh_replic=0,mode_compute=4,mode_file=dsge_base2_mode,mode_check);
//estimation(datafile=extreme,nograph,order=2,filter_algorithm=gmf,nobs=150,mh_replic=0,mode_compute=0,mode_file=dsge_base2_mode,mode_check);

// P-MH with SIR
%options_.mh_posterior_mode_estimation = 1 ;
%options_.mh_jscale = 1.2e-3 ; 
options_.mh_nblck = 10 ;
//estimation(datafile=extreme,order=2,nobs=150,number_of_particles=50000,mh_replic=0,mode_compute=7);
options_.posterior_sampling_method = 'RWGMH';
options_.rwgmh_scale_shock = (1e-5)*[10 10 1 1 10 10 10 1000 10 10] ;
//estimation(datafile=extreme,order=2,nobs=150,number_of_particles=1000,mh_replic=10000,mode_compute=8);
//estimation(datafile=extreme,order=2,nobs=150,number_of_particles=1000,mh_replic=10000,mode_file=dsge_base2_mode,mode_compute=0);
estimation(datafile=extreme,order=1,mh_replic=5000,mode_compute=0,mode_file=dsge_base2_mode);

//estimation(datafile=extreme,order=2,nobs=150,number_of_particles=50000,mode_file=dsge_base2_mode,mh_replic=60000,mode_compute=7);
//estimation(datafile=extreme,order=2,nobs=150,number_of_particles=10000,load_mh_file,mh_replic=1,mode_compute=0);

// Online
//options_.particle.liu_west_delta = 0.9 ;
//estimation(datafile=extreme,order=2,number_of_particles=10000,nobs=150,mode_compute=11);
