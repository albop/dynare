@#include "fs2000.common.inc"

if exist('fminsearch','file')
  estimation(mode_compute=7,order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0);
end
