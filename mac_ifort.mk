FC=ifort
FLAG_COMPILE_AND_ASSEMBLE= -c
FLAGS=-reentrancy threaded -fstack-protector-strong \
	-fprotect-parens -no-fma -mP2OPT_hpo_vec_divbyzero=F \
	-prec-div -prec-sqrt -fimf-arch-consistency=true \
	-fp-model fast=2 -pc64 -fPIC \
	-fpp -auto -init=zero -init=arrays -msse3 \
	-autodouble -align nocommons -warn noalignments \
#	-heap-arrays \
	# -check all -g -traceback -check bounds -check shape -check stack -warn shape \

	# -ipo -qopt-report-phase vec \
	# -check shape -check stack -warn shape \
	# -fp-stack-check -g -traceback \

include common.mk

##
evpsc: $(objects_evpsc_common) $(objects_evpsc)
	python -u pyscripts/statev_interface.py
	cd umat/src/evpsc_incr/&&make -i test&&./test
	$(FC) $(FLAGS) -o evpsc $(objects_evpsc_common) $(objects_evpsc)

.PHONY: all
all: evpsc

$(objects_evpsc): | $(OBJDIR)
$(OBJDIR):
	mkdir $(OBJDIR)
