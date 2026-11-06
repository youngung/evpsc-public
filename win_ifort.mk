OBJDIR:=objdir
SRCDIR:=src/vpsc_src

FC=ifort
FLAG_COMPILE_AND_ASSEMBLE=-c

FLAGS=/check:all /warn:errors  /warn:all \
	/warn:noexternal /warn:nounused /nologo \
	/traceback /O0 /debug:all /check:stack \
	/check:shape /check:bounds /ftrapuv \
	/ffixed-line-length-none \
	/warn:notruncated_source
#/Wno-error=line-truncation

# FLAGS=  \
#	/fpp  /extend-source  /nologo \
#	/DABQ_FORTRAN  /DABQ_WIN86_64 \
#	/iface:cref /recursive \
#	/Qauto \
#	/Align:array64byte \
#	/Qpc64 \
#	/Qprec-div /Qprec-sqrt \
#	/Qfma- \
#	/fp:precise \
#	/Qimf-arch-consistency:true  \
#	/Qfp-speculation:safe     \
#	/Qprotect-parens          \
#	/Qfp-stack-check          \
#	/reentrancy:threaded      \
#	/F10000000000 -heap-arrays \
#	/traceback \
#	/check:all /check:stack /check:shape /check:bounds \
#	/warn:all /warn:noexternal /warn:nounused
#	/Qinit=zero /Qinit=arrays \
#	/Qinit=snan /Qinit=arrays \
# #	/QxSSE3 /QaxAVX
# #	/warn:all /warn:noexternal /warn:nounused #\
# #	/F3500000 -heap-arrays


include common.mk

## VPSC
vpsc: $(objects_vpsc)
	$(FC) $(FLAGS) -o vpsc $(objects_vpsc)
## EVPSC
evpsc: $(objects_evpsc_common) $(objects_evpsc)
	$(FC) $(FLAGS) -o evpsc $(objects_evpsc_common) $(objects_evpsc)

evpsc_dist: $(objects_evpsc_common)
	$(FC) $(FLAGS) -o evpsc $(objects_evpsc_common)

.PHONY: all clean # zip tar
all: vpsc evpsc #vpsc+e  ## evp zip tar

$(objects_vpsc): | $(OBJDIR)
$(objects_evpsc): | $(OBJDIR)
$(objects_evpsc_common): | $(OBJDIR)

$(OBJDIR):
	mkdir $(OBJDIR)
clean:
	del evpsc $(objects_vpsc) $(objects_evpsc) $(objects_evpsc_common) *.txt *.csv *.tar *.tgz *.out *.OUT *~ ACT_PH?.* RERR.* STATS.* TEX_PH?.* TW_SCH_STATS_PH?_MODE*.* php_ph?.* quartz_ss.* RUN_LOG.* STR_STR.* TW_VAR_STATS_PH*_MODE*.* php_ph?.* sc_convergence.* DIF_STR_PH?.* DIF_WGT_PH?.* CAUCHY.* *.png *.pdf D_slip* P_slip* dis_den* pcys_dis* srho_f* srho_r* tauback*  stats_stress_strain* ELSC.*
	#del evpsc $(objects_vpsc) $(objects_evpsc) *.txt *.csv *.tar *.tgz *.out *.OUT *~ ACT_PH?.* RERR.* STATS.* TEX_PH?.* TW_SCH_STATS_PH?_MODE*.* php_ph?.* quartz_ss.* RUN_LOG.* STR_STR.* TW_VAR_STATS_PH*_MODE*.* php_ph?.* sc_convergence.* DIF_STR_PH?.* DIF_WGT_PH?.* CAUCHY.* *.png *.pdf D_slip* P_slip* dis_den* pcys_dis* srho_f* srho_r* tauback*  stats_stress_strain* ELSC.*
	cd umat/src/evpsc_incr/;make clean -i
