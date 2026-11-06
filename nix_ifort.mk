OBJDIR:=objdir
SRCDIR:=src/vpsc_src

FC=ifort
FLAG_COMPILE_AND_ASSEMBLE= -c

FLAGS=-reentrancy threaded -fstack-protector-strong \
	-fprotect-parens -no-fma -mP2OPT_hpo_vec_divbyzero=F \
	-prec-div -prec-sqrt -fimf-arch-consistency=true \
	-fp-model fast=2 -pc64 -fPIC \
	-fpp -auto -init=zero -init=arrays -msse3 \
	-autodouble -align nocommons -warn noalignments  \
	-heap-arrays  # -check all -g -traceback #\
### Below is used for debugging purposes.
#	-ipo -qopt-report-phase vec  -check bounds -g \
#	-check shape -check stack -warn shape \
#	-fp-stack-check -g -traceback \

sed_cmd=sed
include common.mk


## VPSC
vpsc: $(objects_vpsc)
	$(FC) $(FLAGS) -o vpsc $(objects_vpsc)
evpsc: $(objects_evpsc) #vpsc+e
	cd umat/src/evpsc_incr/&&make -i test&&./test
	$(FC) $(FLAGS) -o evpsc $(objects_evpsc)


## archive
DATETIME:=`date +'%Y%m%d_%H%M%S'`

zip: $(objects_evpsc)
	cd $(SRCDIR); zip -r evpsc_$(DATETIME).zip $(src_evpsc_raw); cd -; mv $(SRCDIR)/evpsc_$(DATETIME).zip ./; zip -ur evpsc_$(DATETIME).zip compile_windows.bat
	mkdir /tmp/dum_test
	cp evpsc_$(DATETIME).zip /tmp/dum_test/
	cd /tmp/dum_test/;unzip evpsc_$(DATETIME).zip&&bash compile_windows.bat

tar: $(objects_evpsc)
	cd $(SRCDIR); tar -cvf evpsc_$(DATETIME).tar $(src_evpsc_raw); cd -; mv $(SRCDIR)/evpsc_$(DATETIME).tar ./; tar -rvf evpsc_$(DATETIME).tar compile_windows.bat

.PHONY: all clean zip tar
all: vpsc evpsc #vpsc+e  ## evp zip tar

$(objects_vpsc): | $(OBJDIR)
$(objects_evpsc): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)
clean:
	-rm vpsc evpsc vpsc+e evp $(objects_vpsc) $(objects_evpsc) *.txt *.csv *.tar *.tgz *.out *.OUT *~ ACT_PH?.* RERR.* STATS.* TEX_PH?.* TW_SCH_STATS_PH?_MODE*.* php_ph?.* quartz_ss.* RUN_LOG.* STR_STR.* TW_VAR_STATS_PH*_MODE*.* php_ph?.* sc_convergence.* DIF_STR_PH?.* DIF_WGT_PH?.* CAUCHY.* *.png *.pdf D_slip* P_slip* dis_den* pcys_dis* srho_f* srho_r* tauback*  stats_stress_strain* ELSC.*
	cd umat/src/evpsc_incr/;make clean -i
