.SUFFIXES: .mod .o .f90

FC = ifort
FCFLAGS = -mkl  -fpp -c -qopenmp  -assume nobscc -Dsingle_precision  -I /opt/intel/composer_xe_2015.1.133/mkl/include/fftw -I /opt/intel/composer_xe_2015.1.133/mkl/include/
LDFLAGS = -mkl  -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl   -qopenmp  -I /opt/intel/composer_xe_2015.1.133/mkl/include/ -I /opt/intel/composer_xe_2015.1.133/mkl/include/fftw

SUFFIX = .f90

OBJS =  \
	quadpack.o                 \
	cufft_wrapper.o         \
	m_elsa.o                   \
	m_precision.o              \
	m_string.o                 \
	m_numerical_tools.o        \
	m_crystallography.o        \
	m_electron.o               \
	m_xray_factors.o           \
	global_variables.o     \
	m_user_input.o             \
	m_slicing.o                \
	m_qep.o                    \
	output.o               \
	m_lens.o                   \
	m_tilt.o                   \
	m_absorption.o             \
	m_probe_scan.o             \
	m_potential.o              \
	local_ionization.o     \
	m_multislice.o             \
	s_absorptive_pacbed.o      \
	MS_utilities.o \
	s_qep_pacbed_CPU.o \
	s_absorptive_stem_CPU.o \
	s_qep_tem_CPU.o             \
	s_qep_stem_CPU.o            \
	s_absorptive_tem_CPU.o      \
	muSTEM.o
#F90SRC_MOD = $(OBJS_MOD:.mod=$(SUFFIX))
F90SRC = $(OBJS:.o=$(SUFFIX))

$(OBJS): $(F90SRC) 
	${FC} ${FCFLAGS} -c $*$(SUFFIX)

all: link
link: $(OBJS)
	${FC} ${LDFLAGS}  $(OBJS) -o MuSTEM.x
clean:
	rm *.mod *.o
