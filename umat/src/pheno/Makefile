## This make file is to avoid unnecessary compilation of UMAT,
## repetition of which can be quite annoying particularly when
## testing without changing relevant source codes.

## Unless you are a developer, you are advised to make initially, then every once when
## you obtain an updated version of the code.

# Youngung Jeong
# yjeong@changwon.ac.kr


## Use with '-i' option to ignore errors

OBJDIR:=objdir
FC=ifort
FLAG_COMPILE_AND_ASSEMBLE= -c


## List of dependent files (using their relative paths.
dep_common=epl.f Makefile



## One has to find this executable to correctly use the Make procedure.
## This can differ depending on where you install ABAQUS
ABQ=C:\SIMULIA\EstProducts\2020\win_b64\code\bin\ABQLauncher.exe
ifeq ($(OS),Windows_NT)
	FFLAGS= -reentrancy:threaded  -Qprotect-parens -Qfma \
	-Qimf-arch-consistency:true  -fp:fast=2  -extlnk:.o\
	-nologo -Qpc64 -fpp -Qinit:zero -Qinit:arrays -heap-arrays
else
	FFLAGS= -reentrancy threaded -fstack-protector-strong \
#		-heap-arrays -fpp -nologo
#		-fprotect-parens -no-fma -mP2OPT_hpo_vec_divbyzero=F \
#		-prec-div -prec-sqrt -fimf-arch-consistency=true \
#		-fp-model fast=2 -pc64 -fPIC \
#		-fpp -auto -init=zero -init=arrays -msse3 \
#		-autodouble -align nocommons -warn noalignments  \

endif


## UMAT compile on Windows
epl-std.obj: epl.f $(dep_common)
	rm -f epl-std.obj standardU.dll abaqus_v6.env
	$(ABQ) make library=epl.f

iso-std.obj: iso.f $(dep_common)
	rm -f iso-std.obj standardU.dll abaqus_v6.env
	$(ABQ) make library=iso.f

# ## umat compile @Nurion	# cp /apps/commercial/abaqus/6146/6.14-6/SMA/site/abaqus_v6.env .
# umat-std.o:  $(dep_imp) $(dep_common)
#	cp /apps/commercial/abaqus/2020/EstProducts/2020/linux_a64/SMA/site/abaqus_v6.env .
#	rm -f umat-std.o libstandardU.so
#	/apps/commercial/abaqus/Commands/abq2020hf4 make library=umat.f


# ## UMAT debugging @Nurion
# umat-debug.o:  $(dep_imp) $(dep_common)
#	rm -f umat-std.o umat-debug.o libstandardU.so libstandardU-debug.so
#	cp env/abaqus_v6_debug.env abaqus_v6.env
#	/apps/commercial/abaqus/Commands/abq2020hf4 make library=umat.f
#	mv umat-std.o umat-debug.o&&mv libstandardU.so  libstandardU-debug.so

# vumat-xpl.o: $(dep_v) $(dep_common)
#	python prebuild_vumat.py -i 0
#	rm -f vumat-xpl.o libexplicitU.so vumat-xplD.o libexplicitU-D.so
#	cp env/abaqus_v6_debug.env abaqus_v6.env
#	/apps/commercial/abaqus/Commands/abq2020hf4 make library=vumat.f
#	python prebuild_vumat.py -i 1


# ## Diff software to postmortem conduction of diffraction simulation from px.stv generated from FE simulation
# objects_diff=$(addprefix $(OBJDIR)/, \
#	diffraction.o crystal_symmetry.o chg_basis.o dif_planes_fe.o\
#	save_diffraction.o write_texture_fe.o euler.o\
#	)
# $(OBJDIR)/diffraction.o: src/diffraction.f
#	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/diffraction.f -o $(OBJDIR)/diffraction.o
# $(OBJDIR)/crystal_symmetry.o: ../../../src/vpsc_src/crystal_symmetry.f
#	sed -i "s~include 'ABA_PARAM.INC'~implicit none~g" ../../../src/vpsc_src/crystal_symmetry.f
#	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) ../../../src/vpsc_src/crystal_symmetry.f -o $(OBJDIR)/crystal_symmetry.o
#	sed -i "s~implicit none~include 'ABA_PARAM.INC'~g" ../../../src/vpsc_src/crystal_symmetry.f
# $(OBJDIR)/chg_basis.o: ../../../src/vpsc_src/chg_basis.f
#	sed -i "s~include 'ABA_PARAM.INC'~implicit none~g" ../../../src/vpsc_src/chg_basis.f
#	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) ../../../src/vpsc_src/chg_basis.f -o $(OBJDIR)/chg_basis.o
#	sed -i "s~implicit none~include 'ABA_PARAM.INC'~g" ../../../src/vpsc_src/chg_basis.f
# $(OBJDIR)/write_texture_fe.o: ../../../src/vpsc_src/write_texture_fe.f
#	sed -i "s~include 'ABA_PARAM.INC'~implicit none~g" ../../../src/vpsc_src/write_texture_fe.f
#	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) ../../../src/vpsc_src/write_texture_fe.f -o $(OBJDIR)/write_texture_fe.o
#	sed -i "s~implicit none~include 'ABA_PARAM.INC'~g" ../../../src/vpsc_src/write_texture_fe.f
# $(OBJDIR)/euler.o: ../../../src/vpsc_src/euler.f
#	sed -i "s~include 'ABA_PARAM.INC'~implicit none~g" ../../../src/vpsc_src/euler.f
#	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) ../../../src/vpsc_src/euler.f -o $(OBJDIR)/euler.o
#	sed -i "s~implicit none~include 'ABA_PARAM.INC'~g" ../../../src/vpsc_src/euler.f
# $(OBJDIR)/save_diffraction.o: src/save_diffraction.f
#	sed -i "s~include 'ABA_PARAM.INC'~implicit none~g" src/save_diffraction.f
#	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/save_diffraction.f -o $(OBJDIR)/save_diffraction.o
#	sed -i "s~implicit none~include 'ABA_PARAM.INC'~g" src/save_diffraction.f
# $(OBJDIR)/dif_planes_fe.o: src/dif_planes_fe.f
#	sed -i "s~include 'ABA_PARAM.INC'~implicit none~g" src/dif_planes_fe.f
#	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/dif_planes_fe.f -o $(OBJDIR)/dif_planes_fe.o
#	sed -i "s~implicit none~include 'ABA_PARAM.INC'~g" src/dif_planes_fe.f

# ## Diffraction
# # diff: $(objects_diff) Makefile ../../../src/vpsc_src/incr_evpsc_fe.dim ../../../src/vpsc_src/parameters.dim
# #	$(FC) $(FFLAGS) -o diff $(objects_diff)

# diffraction: $(objects_diff) Makefile ../../../src/vpsc_src/incr_evpsc_fe.dim ../../../src/vpsc_src/parameters.dim
#	$(FC) $(FFLAGS) -o diffraction $(objects_diff)

# ## test state variable size
# objects_test=$(addprefix $(OBJDIR)/, test_nstatv_depvar.o statev_io.o)
# $(OBJDIR)/test_nstatv_depvar.o: src/test_nstatv_depvar.f src/call_statev_io.f ../../../src/vpsc_src/incr_evpsc_fe.dim diag.py ../../../src/vpsc_src/parameters.dim
#	sed -i "s~include 'ABA_PARAM.INC'~implicit none~g" src/test_nstatv_depvar.f
#	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/test_nstatv_depvar.f -o $(OBJDIR)/test_nstatv_depvar.o
#	sed -i "s~implicit none~include 'ABA_PARAM.INC'~g" src/test_nstatv_depvar.f
# $(OBJDIR)/statev_io.o: src/statev_io.f src/call_statev_io.f ../../../src/vpsc_src/incr_evpsc_fe.dim ../../../src/vpsc_src/parameters.dim
#	sed -i "s~include 'ABA_PARAM.INC'~implicit none~g" src/statev_io.f
#	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/statev_io.f -o $(OBJDIR)/statev_io.o
#	sed -i "s~implicit none~include 'ABA_PARAM.INC'~g" src/statev_io.f
# test: $(objects_test) ../../../src/vpsc_src/incr_evpsc_fe.dim src/statev_io.f diag.py ../../../src/vpsc_src/parameters.dim
#	$(FC) $(FFLAGS) -o test $(objects_test)

.PHONY: clean all
all: epl-std.obj iso-std.obj

$(objects_diff): | $(OBJDIR)
$(objects_test): | $(OBJDIR)
$(OBJDIR):
	mkdir $(OBJDIR)
clean:
	rm -fr                       $(objects_test) $(objects_diff) test.exe __pycache__ *__genmod.* *.simdir *.dmp diffraction.exe *.SMABulk *.simdir *exception *.lck fort.10 fort.7 eliT_CheckLicense.*.exception program.exe program *.com *.dat *.log *.msg *.odb *.prt *.sim *.sta abaqus.rpy *.odb_* *~ *.cid *.lck *.simlog *.1.SMABulk *.exception  *.stt *.ipm *.mdl abaqus.rpy.? quartz_ss.out *.023 *.OUT *.out *.gif *.png *.mp4 *.dat *.sim *.sta *.rpy* *.SMABulk *.cid *.jnl *.ipm *.log *.res *.obj *.dll *.fil

#	C:/MinGW/msys/1.0/bin/rm -fr $(objects_test) $(objects_diff) test.exe __pycache__ *__genmod.* *.simdir *.dmp diffraction.exe *.SMABulk *.simdir *exception *.lck fort.10 fort.7 eliT_CheckLicense.*.exception program.exe program *.com *.dat *.log *.msg *.odb *.prt *.sim *.sta abaqus.rpy *.odb_* *~ *.cid *.lck *.simlog *.1.SMABulk *.exception  *.stt *.ipm *.mdl abaqus.rpy.? quartz_ss.out *.023 *.OUT *.out *.gif *.png *.mp4 *.dat *.sim *.sta *.rpy* *.SMABulk *.cid *.jnl *.ipm *.log *.res *.obj *.dll *.fil
# del /Q                       $(objects_test) $(objects_diff) test.exe __pycache__ *__genmod.* *.simdir *.dmp diffraction.exe *.SMABulk *.simdir *exception *.lck fort.10 fort.7 eliT_CheckLicense.*.exception program.exe program *.com *.dat *.log *.msg *.odb *.prt *.sim *.sta abaqus.rpy *.odb_* *~ *.cid *.lck *.simlog *.1.SMABulk *.exception  *.stt *.ipm *.mdl abaqus.rpy.? quartz_ss.out *.023 *.OUT *.out *.gif *.png *.mp4 *.dat *.sim *.sta *.rpy* *.SMABulk *.cid *.jnl *.ipm *.log *.res *.obj *.dll *.fil
