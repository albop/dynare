@#include "fs2000.common.inc"

estimation(mode_compute=6,order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0);

// test the mode file generated with mode_compute=6
estimation(order=1,datafile='../fs2000/fsdat_simul',nobs=192,loglinear,mode_compute=0,mode_file=fs2000_6_mode,mh_replic=10);
