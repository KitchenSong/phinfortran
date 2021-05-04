export FFLAGS=-traceback -debug -O2 -static_intel
export MPIFC=mpif90
MKL=/cm/shared/apps/intel-compiler/1.0.080/composer_xe_2013_sp1.0.080/mkl/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group  \
/cm/shared/apps/intel-compiler/1.0.080/composer_xe_2013_sp1.0.080/mkl/lib/intel64/libmkl_intel_lp64.a               \
/cm/shared/apps/intel-compiler/1.0.080/composer_xe_2013_sp1.0.080/mkl/lib/intel64/libmkl_sequential.a               \
/cm/shared/apps/intel-compiler/1.0.080/composer_xe_2013_sp1.0.080/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
export LAPACK=$(MKL)
export LIBS=$(LAPACK)
