OBJDIR:=objdir
SRCDIR:=src/vpsc_src

included_files=$(addsuffix .f , $(addprefix src/vpsc_src/,\
	call_elsc_evp_fe \
	call_evpsc \
	call_evpsc_input_fe \
	call_initial_guess \
	call_rexdd \
	call_rexdd_srx \
	call_statev_io ))
files_misc=$(addprefix src/vpsc_src/, library7.sub incr_evpsc_fe.dim)
filenames=\
	calc_contribution \
	calc_evptilde \
	calc_rel_err \
	chg_basis \
	cod2grains \
	crossop \
	crystal_symmetry \
	data_crystal_fe \
	data_grain_ebsd \
	data_grain_fe \
	diff_planes \
	elsc_evp_fe \
	eshelby \
	eshelby_rot \
	euler \
	evpsc \
	evpsc_input_fe \
	evpsc_main \
	fake_xit \
	fld_a \
	fld_b \
	funcs \
	gaus \
	grain_info \
	grain_rate_and_moduli_evp \
	grain_stat \
	grain_stress_evp \
	initial_state_guess \
	knuth_shuffle \
	lib \
	lib_write \
	load_condition \
	lu_libs \
	neighbours_fe \
	newton_raphson_evp \
	ph_property_evp_fe \
	phi_delta2vect \
	read_ebsd \
	rexdd \
	rexdd_libs \
	rexdd_nucl \
	rexdd_grow \
	rexdd_init_state \
	rexdd_ori_energy \
	rexdd_vanish \
	rexdd_inherits \
	rexdd_maintain \
	rexdd_disp \
	rexdd_stats \
	rodrigues \
	scale_gamd0s_evp_fe \
	select_points_from_ebsd \
	solve_evp_bc \
	statev_io \
	stat_shear_activity \
	stat_twinning_fe \
	sx_version_tester \
	texture_rotation_fe \
	thermo_emod_cte \
	twin_orientation \
	unload_condition \
	update_crss_dd_fe_new \
	update_crss_rgvb_fe \
	update_crss_voce_fe \
	update_fij_fe \
	update_macro_property \
	update_orientation_fe \
	update_schmid_fe \
	update_shape_fe \
	update_twinning_ptr_fe \
	update_twinning_vft_fe \
	var_dsbar \
	var_vel_grad_evp \
	vector_ang \
	voigt \
	write_activ_fe \
	write_crav \
	write_elsc \
	write_stress_strain_evp \
	write_texture_fe
src_evpsc_raw= \
	$(files_misc)  $(included_files) \
	$(addsuffix .f ,  \
	$(addprefix src/vpsc_src/, \
	$(filenames)\
	) \
	)

# objects_evpsc=$(addprefix $(OBJDIR)/, statev_beta.o)
objects_evpsc_common=$(addprefix $(OBJDIR)/, $(addsuffix .o ,$(filenames)) )
files_evpsc_common=$(addprefix src/vpsc_src/, $(addsuffix .f ,$(filenames)) )



make_files=common.mk mac_ifort.mk win_ifort.mk nix_ifort.mk

## common objects
$(OBJDIR)/calc_contribution.o: $(make_files)  src/vpsc_src/calc_contribution.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/calc_contribution.f -o $(OBJDIR)/calc_contribution.o
$(OBJDIR)/calc_evptilde.o: $(make_files)  src/vpsc_src/calc_evptilde.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/calc_evptilde.f -o $(OBJDIR)/calc_evptilde.o
$(OBJDIR)/calc_rel_err.o: $(make_files)  src/vpsc_src/calc_rel_err.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/calc_rel_err.f -o $(OBJDIR)/calc_rel_err.o
$(OBJDIR)/chg_basis.o: $(make_files)  src/vpsc_src/chg_basis.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/chg_basis.f -o $(OBJDIR)/chg_basis.o
$(OBJDIR)/cod2grains.o: $(make_files)  src/vpsc_src/cod2grains.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/cod2grains.f -o $(OBJDIR)/cod2grains.o
$(OBJDIR)/crossop.o: $(make_files)  src/vpsc_src/crossop.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/crossop.f -o $(OBJDIR)/crossop.o
$(OBJDIR)/crystal_symmetry.o: $(make_files)  src/vpsc_src/crystal_symmetry.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/crystal_symmetry.f -o $(OBJDIR)/crystal_symmetry.o

$(OBJDIR)/data_crystal_fe.o: $(make_files)  src/vpsc_src/data_crystal_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/data_crystal_fe.f -o $(OBJDIR)/data_crystal_fe.o
$(OBJDIR)/data_grain_ebsd.o: $(make_files)  src/vpsc_src/data_grain_ebsd.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/data_grain_ebsd.f -o $(OBJDIR)/data_grain_ebsd.o
$(OBJDIR)/data_grain_fe.o: $(make_files)  src/vpsc_src/data_grain_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/data_grain_fe.f -o $(OBJDIR)/data_grain_fe.o

$(OBJDIR)/diff_planes.o: $(make_files)  src/vpsc_src/diff_planes.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/diff_planes.f -o $(OBJDIR)/diff_planes.o
$(OBJDIR)/elsc_evp_fe.o: $(make_files)  src/vpsc_src/elsc_evp_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/elsc_evp_fe.f -o $(OBJDIR)/elsc_evp_fe.o
$(OBJDIR)/eshelby.o: $(make_files)  src/vpsc_src/eshelby.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/eshelby.f -o $(OBJDIR)/eshelby.o
$(OBJDIR)/eshelby_rot.o: $(make_files)   src/vpsc_src/eshelby_rot.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/eshelby_rot.f -o $(OBJDIR)/eshelby_rot.o
$(OBJDIR)/euler.o: $(make_files)  src/vpsc_src/euler.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/euler.f -o $(OBJDIR)/euler.o
$(OBJDIR)/evpsc.o: $(make_files)  src/vpsc_src/evpsc.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/evpsc.f -o $(OBJDIR)/evpsc.o
$(OBJDIR)/evpsc_input_fe.o: $(make_files)  src/vpsc_src/evpsc_input_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/evpsc_input_fe.f -o $(OBJDIR)/evpsc_input_fe.o
$(OBJDIR)/evpsc_main.o: $(make_files) $(included_files) $(files_misc) src/vpsc_src/evpsc_main.f src/vpsc_src/parameters.dim
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/evpsc_main.f -o $(OBJDIR)/evpsc_main.o

$(OBJDIR)/fake_xit.o: $(make_files)   src/vpsc_src/fake_xit.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/fake_xit.f -o $(OBJDIR)/fake_xit.o
$(OBJDIR)/fld_a.o: $(make_files)  src/vpsc_src/fld_a.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/fld_a.f -o $(OBJDIR)/fld_a.o
$(OBJDIR)/fld_b.o: $(make_files)  src/vpsc_src/fld_b.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/fld_b.f -o $(OBJDIR)/fld_b.o
$(OBJDIR)/funcs.o: $(make_files)  src/vpsc_src/funcs.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/funcs.f -o $(OBJDIR)/funcs.o

$(OBJDIR)/grain_info.o: $(make_files)  src/vpsc_src/grain_info.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/grain_info.f -o $(OBJDIR)/grain_info.o
$(OBJDIR)/grain_rate_and_moduli_evp.o: $(make_files)  src/vpsc_src/grain_rate_and_moduli_evp.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/grain_rate_and_moduli_evp.f -o $(OBJDIR)/grain_rate_and_moduli_evp.o
$(OBJDIR)/grain_stat.o: $(make_files)  src/vpsc_src/grain_stat.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/grain_stat.f -o $(OBJDIR)/grain_stat.o
$(OBJDIR)/grain_stress_evp.o: $(make_files)  src/vpsc_src/grain_stress_evp.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/grain_stress_evp.f  -o $(OBJDIR)/grain_stress_evp.o

$(OBJDIR)/initial_state_guess.o: $(make_files)  src/vpsc_src/initial_state_guess.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/initial_state_guess.f -o $(OBJDIR)/initial_state_guess.o

$(OBJDIR)/knuth_shuffle.o: $(make_files)  src/vpsc_src/knuth_shuffle.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/knuth_shuffle.f -o $(OBJDIR)/knuth_shuffle.o

$(OBJDIR)/lib.o: $(make_files)  src/vpsc_src/lib.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/lib.f -o $(OBJDIR)/lib.o
$(OBJDIR)/lib_write.o: $(make_files)  src/vpsc_src/lib_write.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/lib_write.f -o $(OBJDIR)/lib_write.o
$(OBJDIR)/load_condition.o: $(make_files)  src/vpsc_src/load_condition.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/load_condition.f -o $(OBJDIR)/load_condition.o
$(OBJDIR)/lu_libs.o: $(make_files)  src/vpsc_src/lu_libs.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/lu_libs.f -o $(OBJDIR)/lu_libs.o

$(OBJDIR)/make_temp_fn.o: $(make_files)  src/vpsc_src/make_temp_fn.f src/vpsc_src/vpsc7.dim src/vpsc_src/rgvb.dim
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/make_temp_fn.f -o $(OBJDIR)/make_temp_fn.o

$(OBJDIR)/gaus.o: $(make_files)  src/vpsc_src/gaus.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/gaus.f -o $(OBJDIR)/gaus.o
$(OBJDIR)/gaussian_dispersion.o: $(make_files)  src/vpsc_src/gaussian_dispersion.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/gaussian_dispersion.f -o $(OBJDIR)/gaussian_dispersion.o
$(OBJDIR)/neighbours_fe.o: $(make_files)  src/vpsc_src/neighbours_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/neighbours_fe.f -o $(OBJDIR)/neighbours_fe.o
$(OBJDIR)/newton_raphson_evp.o: $(make_files)  src/vpsc_src/newton_raphson_evp.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/newton_raphson_evp.f -o $(OBJDIR)/newton_raphson_evp.o

$(OBJDIR)/ph_property_evp_fe.o: $(make_files)  src/vpsc_src/ph_property_evp_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/ph_property_evp_fe.f -o $(OBJDIR)/ph_property_evp_fe.o
$(OBJDIR)/phi_delta2vect.o: $(make_files)   src/vpsc_src/phi_delta2vect.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/phi_delta2vect.f -o $(OBJDIR)/phi_delta2vect.o

$(OBJDIR)/read_ebsd.o: $(make_files)  src/vpsc_src/read_ebsd.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/read_ebsd.f -o $(OBJDIR)/read_ebsd.o

## Recrystallization related subroutines
$(OBJDIR)/rexdd.o: $(make_files)  src/vpsc_src/rexdd.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/rexdd.f -o $(OBJDIR)/rexdd.o
$(OBJDIR)/rexdd_nucl.o: $(make_files)  src/vpsc_src/rexdd_nucl.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/rexdd_nucl.f -o $(OBJDIR)/rexdd_nucl.o
$(OBJDIR)/rexdd_grow.o: $(make_files)  src/vpsc_src/rexdd_grow.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/rexdd_grow.f -o $(OBJDIR)/rexdd_grow.o
$(OBJDIR)/rexdd_init_state.o: $(make_files)  src/vpsc_src/rexdd_init_state.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/rexdd_init_state.f -o $(OBJDIR)/rexdd_init_state.o
$(OBJDIR)/rexdd_ori_energy.o: $(make_files)  src/vpsc_src/rexdd_ori_energy.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/rexdd_ori_energy.f -o $(OBJDIR)/rexdd_ori_energy.o
$(OBJDIR)/rexdd_libs.o: $(make_files)  src/vpsc_src/rexdd_libs.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/rexdd_libs.f -o $(OBJDIR)/rexdd_libs.o
$(OBJDIR)/rexdd_vanish.o: $(make_files)  src/vpsc_src/rexdd_vanish.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/rexdd_vanish.f -o $(OBJDIR)/rexdd_vanish.o
$(OBJDIR)/rexdd_inherits.o: $(make_files)  src/vpsc_src/rexdd_inherits.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/rexdd_inherits.f -o $(OBJDIR)/rexdd_inherits.o
$(OBJDIR)/rexdd_maintain.o: $(make_files)  src/vpsc_src/rexdd_maintain.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/rexdd_maintain.f -o $(OBJDIR)/rexdd_maintain.o
$(OBJDIR)/rexdd_disp.o: $(make_files)  src/vpsc_src/rexdd_disp.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/rexdd_disp.f -o $(OBJDIR)/rexdd_disp.o
$(OBJDIR)/rexdd_stats.o: $(make_files)  src/vpsc_src/rexdd_stats.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/rexdd_stats.f -o $(OBJDIR)/rexdd_stats.o


$(OBJDIR)/rodrigues.o: $(make_files)  src/vpsc_src/rodrigues.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/rodrigues.f -o $(OBJDIR)/rodrigues.o

$(OBJDIR)/scale_gamd0s_evp_fe.o: $(make_files)  src/vpsc_src/scale_gamd0s_evp_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/scale_gamd0s_evp_fe.f -o $(OBJDIR)/scale_gamd0s_evp_fe.o
$(OBJDIR)/select_points_from_ebsd.o: $(make_files)  src/vpsc_src/select_points_from_ebsd.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/select_points_from_ebsd.f -o $(OBJDIR)/select_points_from_ebsd.o
$(OBJDIR)/solve_evp_bc.o: $(make_files)  src/vpsc_src/solve_evp_bc.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/solve_evp_bc.f -o $(OBJDIR)/solve_evp_bc.o
$(OBJDIR)/statev_beta.o: $(make_files)  umat/src/evpsc_incr/src/statev_beta.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) umat/src/evpsc_incr/src/statev_beta.f -o $(OBJDIR)/statev_beta.o
$(OBJDIR)/statev_io.o: $(make_files)  src/vpsc_src/statev_io.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/statev_io.f -o $(OBJDIR)/statev_io.o
$(OBJDIR)/stat_shear_activity.o: $(make_files)  src/vpsc_src/stat_shear_activity.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/stat_shear_activity.f -o $(OBJDIR)/stat_shear_activity.o
$(OBJDIR)/stat_twinning_fe.o: $(make_files)  src/vpsc_src/stat_twinning_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/stat_twinning_fe.f -o $(OBJDIR)/stat_twinning_fe.o
$(OBJDIR)/sx_version_tester.o: $(make_files)  src/vpsc_src/sx_version_tester.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/sx_version_tester.f -o $(OBJDIR)/sx_version_tester.o

$(OBJDIR)/texture_rotation_fe.o: $(make_files)  src/vpsc_src/texture_rotation_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/texture_rotation_fe.f -o $(OBJDIR)/texture_rotation_fe.o
$(OBJDIR)/thermo_emod_cte.o: $(make_files)  src/vpsc_src/thermo_emod_cte.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/thermo_emod_cte.f -o $(OBJDIR)/thermo_emod_cte.o
$(OBJDIR)/twin_orientation.o: $(make_files)  src/vpsc_src/twin_orientation.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/twin_orientation.f -o $(OBJDIR)/twin_orientation.o

$(OBJDIR)/unload_condition.o: $(make_files)  src/vpsc_src/unload_condition.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/unload_condition.f -o $(OBJDIR)/unload_condition.o
$(OBJDIR)/update_crss_dd_fe_new.o: $(make_files)  src/vpsc_src/update_crss_dd_fe_new.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/update_crss_dd_fe_new.f -o $(OBJDIR)/update_crss_dd_fe_new.o
$(OBJDIR)/update_crss_rgvb_fe.o: $(make_files)  src/vpsc_src/update_crss_rgvb_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/update_crss_rgvb_fe.f -o $(OBJDIR)/update_crss_rgvb_fe.o
$(OBJDIR)/update_crss_voce_fe.o: $(make_files)  src/vpsc_src/update_crss_voce_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/update_crss_voce_fe.f -o $(OBJDIR)/update_crss_voce_fe.o
$(OBJDIR)/update_fij_fe.o: $(make_files)  src/vpsc_src/update_fij_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/update_fij_fe.f -o $(OBJDIR)/update_fij_fe.o
$(OBJDIR)/update_macro_property.o: $(make_files)  src/vpsc_src/update_macro_property.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/update_macro_property.f -o $(OBJDIR)/update_macro_property.o
$(OBJDIR)/update_orientation_fe.o: $(make_files)  src/vpsc_src/update_orientation_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/update_orientation_fe.f -o $(OBJDIR)/update_orientation_fe.o
$(OBJDIR)/update_schmid_fe.o: $(make_files)  src/vpsc_src/update_schmid_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/update_schmid_fe.f -o $(OBJDIR)/update_schmid_fe.o
$(OBJDIR)/update_shape_fe.o: $(make_files)  src/vpsc_src/update_shape_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/update_shape_fe.f -o $(OBJDIR)/update_shape_fe.o
$(OBJDIR)/update_twinning_ptr_fe.o: $(make_files)  src/vpsc_src/update_twinning_ptr_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/update_twinning_ptr_fe.f -o $(OBJDIR)/update_twinning_ptr_fe.o
$(OBJDIR)/update_twinning_vft_fe.o: $(make_files)  src/vpsc_src/update_twinning_vft_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/update_twinning_vft_fe.f -o $(OBJDIR)/update_twinning_vft_fe.o

$(OBJDIR)/var_dsbar.o: $(make_files)  src/vpsc_src/var_dsbar.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/var_dsbar.f -o $(OBJDIR)/var_dsbar.o
$(OBJDIR)/var_vel_grad_evp.o: $(make_files)  src/vpsc_src/var_vel_grad_evp.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/var_vel_grad_evp.f -o $(OBJDIR)/var_vel_grad_evp.o
$(OBJDIR)/vector_ang.o: $(make_files)  src/vpsc_src/vector_ang.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/vector_ang.f -o $(OBJDIR)/vector_ang.o
$(OBJDIR)/voigt.o: $(make_files)  src/vpsc_src/voigt.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/voigt.f -o $(OBJDIR)/voigt.o

$(OBJDIR)/write_activ_fe.o: $(make_files)  src/vpsc_src/write_activ_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/write_activ_fe.f -o $(OBJDIR)/write_activ_fe.o
$(OBJDIR)/write_elsc.o: $(make_files)  src/vpsc_src/write_elsc.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/write_elsc.f -o $(OBJDIR)/write_elsc.o
$(OBJDIR)/write_crav.o: $(make_files)  src/vpsc_src/write_crav.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/write_crav.f -o $(OBJDIR)/write_crav.o
$(OBJDIR)/write_stress_strain_evp.o: $(make_files)  src/vpsc_src/write_stress_strain_evp.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/write_stress_strain_evp.f -o $(OBJDIR)/write_stress_strain_evp.o
$(OBJDIR)/write_texture_fe.o: $(make_files)  src/vpsc_src/write_texture_fe.f
	$(FC) $(FLAGS) $(FLAG_COMPILE_AND_ASSEMBLE) src/vpsc_src/write_texture_fe.f -o $(OBJDIR)/write_texture_fe.o
