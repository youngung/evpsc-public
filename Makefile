.ONESHELL:
include common.mk

FC=gfortran
FLAG_COMPILE_AND_ASSEMBLE= -c

ENV ?= run

FLAGS= -I./src/vpsc_src/ -finit-integer=-999999 -finit-logical=false -fno-automatic\
	-fno-align-commons -fbackslash \
	-fdefault-double-8 -fdefault-real-8  \
	-Wuninitialized  -O$(opt) -Wno-tabs \
	-g -fbacktrace \
	-fbounds-check -finit-real=snan \
	-Wall -Werror -Wextra -Wno-do-subscript \
	-Wno-compare-reals -Wno-unused  -Wno-zerotrip  -fcheck=array-temps \
#       -Wtabs  -Wno-unused-label  -ffast-math \
	-Wno-unused-parameter -Wno-unused-variable  -Wno-unused-dummy-argument \
	-Wno-conversion \
	-fsanitize=address  -Werror \
	#-static-libasan \
	-ffpe-trap=zero,overflow,underflow,precision,denormal

# FLAGS=-I./src/vpsc_src/ -finit-local-zero -Og -g -fbacktrace -fcheck=bounds -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan -std=legacy -fdefault-real-8 -fdefault-double-8 -finit-derived -finit-integer=-999999

ifeq ($(ENV), run)
	## fortran compiler flags
	FLAGS=-I./src/vpsc_src/ -O3 -g -fdefault-real-8 -std=legacy -fbounds-check #(YJ: Do we want to use '-Wall' flag to remove any potential glitches?)
endif

ifeq ($(ENV), debug)
	# below is Andrea's flags
	FLAGS=-I./src/vpsc_src/ -O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan -std=legacy -fdefault-real-8 -fdefault-double-8 -finit-derived -finit-integer=-999999
endif

ifeq ($(ENV), fulldebug)
	FLAGS=-I./src/vpsc_src/ -O0 -g -fcheck=all -Wall -Wextra -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -fsanitize=address
	FLAGS=-I./src/vpsc_src/ -O0 -g -fcheck=all -Wall -Wextra -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -fsanitize=address -Wno-compare-reals -Wno-conversion -Wno-unused -Wno-zerotrip
endif

ifeq ($(ENV), fulldebug2)
	FLAGS=-I./src/vpsc_src/ -O0 -g -fcheck=all -Wno-compare-reals -Wno-unused -Wall -Wextra -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -fsanitize=address
endif



## archive
DATETIME=`date +'%Y%m%d_%H%M%S'`
hash=`git rev-parse HEAD|cut -c1-8`
tar_filename=evpsc_$(DATETIME)_$(hash).tar
zip_filename=evpsc_$(DATETIME)_$(hash).zip

$(eval tmpdir ?= $(shell mktemp -d))

zip: $(objects_evpsc_common)
	@ rm -fr /tmp/ex02_fcc
	@ echo '** zipping'
	@ zip -qr evpsc_dist_temp.zip $(src_evpsc_raw) src/vpsc_src/parameters.dim pyscripts/build.py
	@ zip -qg evpsc_dist_temp.zip Makefile *.mk
	@ zip -qrg evpsc_dist_temp.zip examples/ex02_FCC/
	@ echo '** copy'
	@ cp evpsc_dist_temp.zip $(tmpdir)
	@ mv evpsc_dist_temp.zip evpsc_dist_$(hash)_$(DATETIME).zip
	@ echo "** unzipping"
	@ (cd $(tmpdir) && unzip -q evpsc_dist_temp.zip)
	@ echo "Unzipped location: ${tmpdir}"

zip_win: $(objects_evpsc_common)
	del -fr c:/temp/dum-zip
	del -fr evpsc_dist_temp.zip
	zip -r evpsc_dist_temp.zip $(src_evpsc_raw) src/vpsc_src/parameters.dim
	zip -g evpsc_dist_temp.zip Makefile *.mk
	zip -rg evpsc_dist_temp.zip examples/ex02_FCC/
	mkdir -p c:/temp/dum-zip
	cp evpsc_dist_temp.zip c:/temp/dum-zip
	mv evpsc_dist_temp.zip evpsc_dist_$(hash)_$(DATETIME).zip
	cd c:/temp/dum-zip; unzip evpsc_dist_temp.zip;

tar: $(objects_evpsc_common)
	del -fr /tmp/dum-tar/
	del -fr evpsc_dist.tar
	tar -cf evpsc_dist.tar $(src_evpsc_raw)
	tar -rf evpsc_dist.tar Makefile *.mk
	tar -rf evpsc_dist.tar examples/ex02_FCC/
	mkdir -p /tmp/dum-tar/
	cp evpsc_dist.tar /tmp/dum-tar/
	mv evpsc_dist.tar evpsc_dist_$(hash)_$(DATETIME).tar
	cd /tmp/dum-tar/; tar -xvf evpsc_dist.tar; make ex02_FCC -f ex.mk

tar_umat:
	tar -cf evpsc_dist.tar
	tar -rf evpsc_dist.tar Makefile *.mk

evpsc: $(objects_evpsc_common) $(files_evpsc_common)
	$(FC) $(FLAGS) -o evpsc $(objects_evpsc_common)

evpscf: $(files_evpsc_common)
	$(FC) $(FLAGS) -o evpsc $(files_evpsc_common)


parsers: objdir/sx_version_tester.o
	$(FC) objdir/sx_version_tester.o fortran_parsers/parser_fnsx.f -o sx_parser

docs:
	make -j all -f make_docs.mk

## distributed version
evpsc_dist: $(objects_evpsc_common)
	$(FC) $(FLAGS) -o evpsc $(objects_evpsc_common)

.PHONY: all clean zip tar evpsc_dist parsers docs evpsc
all: evpsc parsers docs

$(objects_evpsc_common): | $(OBJDIR)
$(OBJDIR):
	mkdir $(OBJDIR)
clean:
	del vpsc evpsc vpsc+e evp    $(objects_evpsc_common) *.txt *.csv *.tar *.tgz *.out *.OUT *~ ACT_PH?.* RERR.* STATS.* TEX_PH?.* TW_SCH_STATS_PH?_MODE*.* php_ph?.* quartz_ss.* RUN_LOG.* STR_STR.* TW_VAR_STATS_PH*_MODE*.* php_ph?.* sc_convergence.* DIF_STR_PH?.* DIF_WGT_PH?.* CAUCHY.* *.png *.pdf D_slip* P_slip* dis_den* pcys_dis* srho_f* srho_r* tauback*  stats_stress_strain* ELSC.*
	rm -fr vpsc evpsc vpsc+e evp $(objects_evpsc_common) *.txt *.csv *.tar *.tgz *.out *.OUT *~ ACT_PH?.* RERR.* STATS.* TEX_PH?.* TW_SCH_STATS_PH?_MODE*.* php_ph?.* quartz_ss.* RUN_LOG.* STR_STR.* TW_VAR_STATS_PH*_MODE*.* php_ph?.* sc_convergence.* DIF_STR_PH?.* DIF_WGT_PH?.* CAUCHY.* *.png *.pdf D_slip* P_slip* dis_den* pcys_dis* srho_f* srho_r* tauback*  stats_stress_strain* ELSC.*
