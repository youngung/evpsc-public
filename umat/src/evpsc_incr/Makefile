## This make file is to avoid unnecessary compilation of UMAT,
## repetition of which can be quite annoying particularly when
## testing without changing relevant source codes.

## Unless you are a developer, you are advised to make initially, then every once when
## you obtain an updated version of the code.

# Youngung Jeong
# yjeong@changwon.ac.kr


## Use with '-i' option to ignore errors
OBJDIR:=objdir


## One has to find this executable to correctly use the Make procedure.
## This can differ depending on where you install ABAQUS
#ABQ=C:\SIMULIA\EstProducts\2020\win_b64\code\bin\ABQLauncher.exe
ABQ=C:\SIMULIA\EstProducts\2021\win_b64\code\bin\ABQLauncher.exe
ifneq ("$(wildcard $(ABQ))","")
	ABQ=C:\SIMULIA\EstProducts\2021\win_b64\code\bin\ABQLauncher.exe
else
	ABQ=C:\SIMULIA\EstProducts\2020\win_b64\code\bin\ABQLauncher.exe
endif

FC ?= ifort
FLAG_COMPILE_AND_ASSEMBLE= -c


## List of dependent files (using their relative paths.
dep_common=src/load_statev.f src/io_head.f \
	umat-nstv.f \
	src/lib_fe.f \
	src/lib_write_fe.f \
	src/voigt_abq.f \
	src/save_diffraction.f \
	src/calc_devm.f \
	src/rodrigues.f \
	src/twin_orientation_fe.f \
	src/prep_evpsc_call.f \
	src/init_uw_print.f \
	src/update_evpsc_statev.f \
	src/calc_spin.f \
	src/loading_condition_fe.f \
	src/uexternaldb.f \
	../../../src/vpsc_src/statev_io.f \
	../../../src/vpsc_src/call_statev_io.f \
	../../../src/vpsc_src/evpsc_input_fe.f \
	../../../src/vpsc_src/call_evpsc.f \
	../../../src/vpsc_src/call_evpsc_input_fe.f \
	../../../src/vpsc_src/call_elsc_evp_fe.f \
	../../../src/vpsc_src/call_initial_guess.f \
	../../../src/vpsc_src/eshelby_rot.f \
	../../../src/vpsc_src/elsc_evp_fe.f \
	../../../src/vpsc_src/chg_basis.f \
	../../../src/vpsc_src/eshelby.f \
	../../../src/vpsc_src/lu_libs.f \
	../../../src/vpsc_src/funcs.f \
	../../../src/vpsc_src/euler.f \
	../../../src/vpsc_src/crystal_symmetry.f \
	../../../src/vpsc_src/voigt.f \
	../../../src/vpsc_src/neighbours_fe.f \
	../../../src/vpsc_src/library7.sub \
	../../../src/vpsc_src/solve_evp_bc.f \
	../../../src/vpsc_src/initial_state_guess.f \
	../../../src/vpsc_src/grain_rate_and_moduli_evp.f \
	../../../src/vpsc_src/evpsc.f \
	../../../src/vpsc_src/grain_stress_evp.f \
	../../../src/vpsc_src/newton_raphson_evp.f \
	../../../src/vpsc_src/calc_evptilde.f \
	../../../src/vpsc_src/update_macro_property.f \
	../../../src/vpsc_src/update_shape_fe.f \
	../../../src/vpsc_src/data_crystal_fe.f \
	../../../src/vpsc_src/update_crss_voce_fe.f \
	../../../src/vpsc_src/update_crss_rgvb_fe.f \
	../../../src/vpsc_src/data_grain_fe.f \
	../../../src/vpsc_src/scale_gamd0s_evp_fe.f \
	../../../src/vpsc_src/update_fij_fe.f \
	../../../src/vpsc_src/update_orientation_fe.f \
	../../../src/vpsc_src/update_twinning_ptr_fe.f \
	../../../src/vpsc_src/update_twinning_vft_fe.f \
	../../../src/vpsc_src/write_texture_fe.f \
	../../../src/vpsc_src/stat_twinning_fe.f \
	../../../src/vpsc_src/stat_shear_activity.f \
	../../../src/vpsc_src/diff_planes.f \
	../../../src/vpsc_src/incr_evpsc_fe.dim \
	../../../src/vpsc_src/parameters.dim \
	../../../src/vpsc_src/fake_xit.f \
	Makefile

dep_v= vumat.f vumat-io.f src/write_inquiry_vumat.f src/write_inquiry_vumat_call.f src/voigt_exp.f src/loading_condition_fe_explicit.f src/hencky.f src/calc_drot_explicit.f src/update_eacc_explicit.f
dep_imp=umat.f  umat-io.f src/write_inquiry_umat.f  src/uvarm.f src/uvarm_custom.f sdvini.f




## UMAT compile on Windows
umat-std.obj: $(dep_imp) $(dep_common)
	rm -f umat-std.obj standardU.dll abaqus_v6.env
	$(ABQ) make library=umat.f
## VUMAT compile on windows
vumat-xplD.obj: $(dep_v) $(dep_common)
	rm -f vumat-xplD.obj vumat-xpl.obj explicitU.dll explicitU-D.dll
	$(ABQ) make library=vumat.f



## umat compile @Nurion	# cp /apps/commercial/abaqus/6146/6.14-6/SMA/site/abaqus_v6.env .
umat-std.o:  $(dep_imp) $(dep_common)
	cp /apps/commercial/abaqus/2020/EstProducts/2020/linux_a64/SMA/site/abaqus_v6.env .
	rm -f umat-std.o libstandardU.so
	/apps/commercial/abaqus/Commands/abq2020hf4 make library=umat.f


## UMAT debugging @Nurion
umat-debug.o:  $(dep_imp) $(dep_common)
	rm -f umat-std.o umat-debug.o libstandardU.so libstandardU-debug.so
	cp env/abaqus_v6_debug.env abaqus_v6.env
	/apps/commercial/abaqus/Commands/abq2020hf4 make library=umat.f
	mv umat-std.o umat-debug.o&&mv libstandardU.so  libstandardU-debug.so

vumat-xpl.o: $(dep_v) $(dep_common)
	# python prebuild_vumat.py -i 0
	rm -f vumat-xpl.o libexplicitU.so vumat-xplD.o libexplicitU-D.so
	cp env/abaqus_v6_debug.env abaqus_v6.env
	/apps/commercial/abaqus/Commands/abq2020hf4 make library=vumat.f
	# python prebuild_vumat.py -i 1


## Diff software to postmortem conduction of diffraction simulation from px.stv generated from FE simulation
objects_diff=$(addprefix $(OBJDIR)/, \
	diffraction.o crystal_symmetry.o chg_basis.o dif_planes_fe.o \
	save_diffraction.o write_texture_fe.o euler.o fake_xit.o lu_libs.o \
	)
$(OBJDIR)/fake_xit.o:  ../../../src/vpsc_src/fake_xit.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) ../../../src/vpsc_src/fake_xit.f -o $(OBJDIR)/fake_xit.o
$(OBJDIR)/diffraction.o: src/diffraction.f
	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/diffraction.f -o $(OBJDIR)/diffraction.o
$(OBJDIR)/crystal_symmetry.o: ../../../src/vpsc_src/crystal_symmetry.f
	sed -i "s~include 'ABA_PARAM.INC'~implicit none~g" ../../../src/vpsc_src/crystal_symmetry.f
	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) ../../../src/vpsc_src/crystal_symmetry.f -o $(OBJDIR)/crystal_symmetry.o
	sed -i "s~implicit none~include 'ABA_PARAM.INC'~g" ../../../src/vpsc_src/crystal_symmetry.f
$(OBJDIR)/chg_basis.o: ../../../src/vpsc_src/chg_basis.f
	sed -i "s~include 'ABA_PARAM.INC'~implicit none~g" ../../../src/vpsc_src/chg_basis.f
	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) ../../../src/vpsc_src/chg_basis.f -o $(OBJDIR)/chg_basis.o
	sed -i "s~implicit none~include 'ABA_PARAM.INC'~g" ../../../src/vpsc_src/chg_basis.f
$(OBJDIR)/write_texture_fe.o: ../../../src/vpsc_src/write_texture_fe.f
	sed -i "s~include 'ABA_PARAM.INC'~implicit none~g" ../../../src/vpsc_src/write_texture_fe.f
	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) ../../../src/vpsc_src/write_texture_fe.f -o $(OBJDIR)/write_texture_fe.o
	sed -i "s~implicit none~include 'ABA_PARAM.INC'~g" ../../../src/vpsc_src/write_texture_fe.f
$(OBJDIR)/euler.o: ../../../src/vpsc_src/euler.f
	sed -i "s~include 'ABA_PARAM.INC'~implicit none~g" ../../../src/vpsc_src/euler.f
	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) ../../../src/vpsc_src/euler.f -o $(OBJDIR)/euler.o
	sed -i "s~implicit none~include 'ABA_PARAM.INC'~g" ../../../src/vpsc_src/euler.f
$(OBJDIR)/lu_libs.o: ../../../src/vpsc_src/lu_libs.f
	sed -i "s~include 'ABA_PARAM.INC'~implicit none~g" ../../../src/vpsc_src/lu_libs.f
	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) ../../../src/vpsc_src/lu_libs.f -o $(OBJDIR)/lu_libs.o
	sed -i "s~implicit none~include 'ABA_PARAM.INC'~g" ../../../src/vpsc_src/lu_libs.f
$(OBJDIR)/save_diffraction.o: src/save_diffraction.f
	sed -i "s~include 'ABA_PARAM.INC'~implicit none~g" src/save_diffraction.f
	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/save_diffraction.f -o $(OBJDIR)/save_diffraction.o
	sed -i "s~implicit none~include 'ABA_PARAM.INC'~g" src/save_diffraction.f
$(OBJDIR)/dif_planes_fe.o: src/dif_planes_fe.f
	sed -i "s~include 'ABA_PARAM.INC'~implicit none~g" src/dif_planes_fe.f
	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/dif_planes_fe.f -o $(OBJDIR)/dif_planes_fe.o
	sed -i "s~implicit none~include 'ABA_PARAM.INC'~g" src/dif_planes_fe.f
## Diffraction
# diff: $(objects_diff) Makefile ../../../src/vpsc_src/incr_evpsc_fe.dim ../../../src/vpsc_src/parameters.dim
#	$(FC) $(FFLAGS) -o diff $(objects_diff)

diffraction: $(objects_diff) Makefile ../../../src/vpsc_src/incr_evpsc_fe.dim ../../../src/vpsc_src/parameters.dim

	$(FC) $(FFLAGS) -o diffraction $(objects_diff)

## test state variable size
objects_test=$(addprefix $(OBJDIR)/, test_nstatv_depvar.o statev_io.o)
$(OBJDIR)/test_nstatv_depvar.o: src/test_nstatv_depvar.f ../../../src/vpsc_src/call_statev_io.f ../../../src/vpsc_src/incr_evpsc_fe.dim diag.py ../../../src/vpsc_src/parameters.dim
	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/test_nstatv_depvar.f -o $(OBJDIR)/test_nstatv_depvar.o
$(OBJDIR)/statev_io.o: ../../../src/vpsc_src/statev_io.f ../../../src/vpsc_src/call_statev_io.f ../../../src/vpsc_src/incr_evpsc_fe.dim ../../../src/vpsc_src/parameters.dim
	$(FC) $(FFLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) ../../../src/vpsc_src/statev_io.f -o $(OBJDIR)/statev_io.o
test: $(objects_test) ../../../src/vpsc_src/incr_evpsc_fe.dim ../../../src/vpsc_src/statev_io.f diag.py ../../../src/vpsc_src/parameters.dim
	$(FC) $(FFLAGS) -o test $(objects_test)

.PHONY: clean all
all: diffraction test umat-std.obj

$(objects_diff): | $(OBJDIR)
$(objects_test): | $(OBJDIR)
$(OBJDIR):
	mkdir $(OBJDIR)
clean:
	rm -fr                       $(objects_test) $(objects_diff) test.exe __pycache__ *__genmod.* *.simdir *.dmp diffraction.exe *.SMABulk *.simdir *exception *.lck fort.10 fort.7 eliT_CheckLicense.*.exception program.exe program *.com *.dat *.log *.msg *.odb *.prt *.sim *.sta abaqus.rpy *.odb_* *~ *.cid *.lck *.simlog *.1.SMABulk *.exception  *.stt *.ipm *.mdl abaqus.rpy.? quartz_ss.out *.023 *.OUT *.out *.gif *.png *.mp4 *.dat *.sim *.sta *.rpy* *.SMABulk *.cid *.jnl *.ipm *.log *.res *.obj *.dll *.fil test test.exe

#	C:/MinGW/msys/1.0/bin/rm -fr $(objects_test) $(objects_diff) test.exe __pycache__ *__genmod.* *.simdir *.dmp diffraction.exe *.SMABulk *.simdir *exception *.lck fort.10 fort.7 eliT_CheckLicense.*.exception program.exe program *.com *.dat *.log *.msg *.odb *.prt *.sim *.sta abaqus.rpy *.odb_* *~ *.cid *.lck *.simlog *.1.SMABulk *.exception  *.stt *.ipm *.mdl abaqus.rpy.? quartz_ss.out *.023 *.OUT *.out *.gif *.png *.mp4 *.dat *.sim *.sta *.rpy* *.SMABulk *.cid *.jnl *.ipm *.log *.res *.obj *.dll *.fil
# del /Q                       $(objects_test) $(objects_diff) test.exe __pycache__ *__genmod.* *.simdir *.dmp diffraction.exe *.SMABulk *.simdir *exception *.lck fort.10 fort.7 eliT_CheckLicense.*.exception program.exe program *.com *.dat *.log *.msg *.odb *.prt *.sim *.sta abaqus.rpy *.odb_* *~ *.cid *.lck *.simlog *.1.SMABulk *.exception  *.stt *.ipm *.mdl abaqus.rpy.? quartz_ss.out *.023 *.OUT *.out *.gif *.png *.mp4 *.dat *.sim *.sta *.rpy* *.SMABulk *.cid *.jnl *.ipm *.log *.res *.obj *.dll *.fil
