// This file deals with the resolution and estimation of a basic DSGE model with
//employment for comparison with the benchmark in Gauss which solves with
//the same particular filter but global methodology.
//
// December 2014

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
bet, uniform_pdf,,, 0.0001, 0.9999;
tet, uniform_pdf,,, 0.0001, .999;
tau, uniform_pdf,,, 0.0001, 100;
delt, uniform_pdf,,, 0.0001, 0.05;
rho, uniform_pdf,,, 0.0001, 0.9999;
stderr e_a, uniform_pdf,,, 0.00001, 0.1;
stderr y, uniform_pdf,,, 0.00001, 0.1;
stderr l, uniform_pdf,,, 0.00001, 0.1;
stderr i, uniform_pdf,,, 0.00001, 0.1;
end;

//estimated_params_init;
//alp, 0.4;   
//bet, 0.99;
//tet, 0.357;
//tau, 50;
//delt, 0.02;
//rho, 0.95;
//stderr e_a, .035;
//stderr y, .00158;
//stderr l, .0011;
//stderr i, .000866;
//end;

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

varobs y l i ;

options_.mode_check.neighbourhood_size = .05 ;
options_.mode_check.number_of_points = 250;

//estimation(datafile=extreme,order=2,number_of_particles=10000,resampling=generic,mode_compute=7,mh_replic=0,mode_check);

//estimation(datafile=extreme,order=2,number_of_particles=10000,nobs=150,mode_compute=7,mode_check);

//estimation(datafile=extreme,order=2,number_of_particles=50000,nobs=150,resampling_method=smooth,mode_compute=7);
//OK

//estimation(datafile=extreme,order=1,nobs=150,mode_compute=7);
// OK

//estimation(datafile=extreme,order=2,filter_algorithm=gf,nobs=150,mh_replic=0,mode_compute=7,mode_check);
//estimation(datafile=extreme,order=2,filter_algorithm=gf,distribution_approximation=montecarlo,number_of_particles=50000,nobs=150,mh_replic=0,mode_compute=7,mode_check);
// OK

//estimation(datafile=extreme,order=2,filter_algorithm=gmf,nobs=150,mh_replic=0,mode_compute=7,mode_check);
//OK

estimation(datafile=risky,order=2,noconstant,filter_algorithm=apf,number_of_particles=10000,nobs=150,mh_replic=0,mode_compute=7);

//options_.mh_posterior_mode_estimation = 1 ;
//options_.mh_jscale =1.2e-3 ; 
//options_.mh_nblck = 10 ;
//estimation(datafile=extreme,order=2,nobs=150,number_of_particles=50000,mh_replic=0,mode_compute=7);
//estimation(datafile=risky,order=2,nobs=150,number_of_particles=50000,mh_replic=0,mode_compute=7);
//options_.posterior_sampling_method = 'RWGMH';
//estimation(datafile=extreme,order=2,nobs=150,number_of_particles=50000,mode_file=dsge_base2_mode,mh_replic=60000,mode_compute=7);
//estimation(datafile=extreme,order=2,nobs=150,number_of_particles=10000,load_mh_file,mh_replic=1,mode_compute=0);
//OK

//options_.particle.liu_west_delta = 0.99 ;
//estimation(datafile=extreme,order=2,number_of_particles=20000,nobs=150,mode_compute=11);
//OK
