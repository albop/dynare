ACLOCAL_AMFLAGS = -I ../../../m4

# libdynare++ must come before gensylv, k_order_perturbation, dynare_simul_
if DO_SOMETHING
SUBDIRS = mjdgges kronecker bytecode libdynare++ gensylv qzcomplex ordschur block_kalman_filter sobol local_state_space_iterations linsolve

if HAVE_MATIO
SUBDIRS += k_order_perturbation dynare_simul_
endif

if HAVE_GSL
if HAVE_MATIO
SUBDIRS += ms_sbvar estimation
endif
endif

if HAVE_SLICOT
SUBDIRS += kalman_steady_state
endif

endif
