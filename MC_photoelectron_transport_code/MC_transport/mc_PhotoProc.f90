!-------------------------------------------------
!
SUBROUTINE Emit_Photoelectrons_From_This_PFL_Point( i, &
                                                  & jsol, &
                                                  & jchannel, &
                                                  & particle_count, &
                                                  & step_count, &
                                                  & mirror_count )

  USE Photoelectrons
  USE field_line

  implicit none
  
  INTEGER, INTENT(IN) :: i
  INTEGER, INTENT(IN) :: jsol
  INTEGER, INTENT(IN) :: jchannel

  INTEGER, INTENT(INOUT) :: particle_count
  INTEGER(8), INTENT(INOUT) :: step_count
  INTEGER, INTENT(INOUT) :: mirror_count

  REAL(8) neutral_m3
  REAL sigma_cm2
  REAL(8) ioniz_threshold_eV

  SELECT CASE (jchannel)


! O -----------------------------------------------------

     CASE (1)

        neutral_m3 = shared_pfl_point(i)%Nn_O_m3
        sigma_cm2 = sigma_ion_O_4S_cm2(jsol)
        ioniz_threshold_eV = threshold_photoion_O_4S_eV ! 13.05 eV, 950 A

     CASE (2)

        neutral_m3 = shared_pfl_point(i)%Nn_O_m3
        sigma_cm2 = sigma_ion_O_2D_cm2(jsol)
        ioniz_threshold_eV = threshold_photoion_O_2D_eV ! 16.53 eV, 750 A

     CASE (3)

        neutral_m3 = shared_pfl_point(i)%Nn_O_m3
        sigma_cm2 = sigma_ion_O_2P_cm2(jsol)
        ioniz_threshold_eV = threshold_photoion_O_2P_eV ! 17.71 eV, 700 A

     CASE (4)

        neutral_m3 = shared_pfl_point(i)%Nn_O_m3
        sigma_cm2 = sigma_ion_O_4P_cm2(jsol)
        ioniz_threshold_eV = threshold_photoion_O_4P_eV ! 27.55 eV, 450 A

     CASE (5)

        neutral_m3 = shared_pfl_point(i)%Nn_O_m3
        sigma_cm2 = sigma_ion_O_2Pst_cm2(jsol)
        ioniz_threshold_eV = threshold_photoion_O_2Pst_eV ! 40.78 eV, 304 A

! N2 -----------------------------------------------------

     CASE (6)

        neutral_m3 = shared_pfl_point(i)%Nn_N2_m3
        sigma_cm2 = sigma_ion_N2_to_N2_cm2(jsol)
        ioniz_threshold_eV = threshold_photoion_N2_to_N2_eV ! 15.49 eV, 700 A

     CASE (7)

        neutral_m3 = shared_pfl_point(i)%Nn_N2_m3
        sigma_cm2 = sigma_ion_N2_to_N_cm2(jsol)
        ioniz_threshold_eV = threshold_photoion_N2_to_N_eV ! 22.54 eV, 550 A

! O2 -----------------------------------------------------

     CASE (8)

        neutral_m3 = shared_pfl_point(i)%Nn_O2_m3
        sigma_cm2 = sigma_ion_O2_to_O2_cm2(jsol)
        ioniz_threshold_eV = threshold_photoion_O2_to_O2_eV ! 11.81 eV, 1050 A

     CASE (9)

        neutral_m3 = shared_pfl_point(i)%Nn_O2_m3
        sigma_cm2 = sigma_ion_O2_to_O_cm2(jsol)
        ioniz_threshold_eV = threshold_photoion_O2_to_O_eV ! 17.71 eV, 700 A

! He -----------------------------------------------------

     CASE (10)

        neutral_m3 = shared_pfl_point(i)%Nn_He_m3
        sigma_cm2 = sigma_ion_He_cm2(jsol)
        ioniz_threshold_eV = threshold_photoion_He_eV ! 22.54 eV, 550 A

  END SELECT

  CALL START_PFL_PHOTOELECTRONS_WITH_GIVEN_ENERGY( i, &
                                                 & jsol, &
                                                 & neutral_m3, &
                                                 & sigma_cm2, &
                                                 & ioniz_threshold_eV, &
                                                 & particle_count, &
                                                 & step_count, &
                                                 & mirror_count )

END SUBROUTINE Emit_Photoelectrons_From_This_PFL_Point

!----------------------------------------------------------
!
SUBROUTINE START_PFL_PHOTOELECTRONS_WITH_GIVEN_ENERGY( i, &
                                                     & jsol, &
                                                     & my_Nn_m3, &
                                                     & sigma_cm2, &
                                                     & threshold_energy_eV, &
                                                     & particle_count, &
                                                     & step_count, &
                                                     & mirror_count )

  USE field_line
  USE rng_wrapper
  USE PhysicalConstants
  USE Photoelectrons
  USE ParallelOperationValues

  USE MPI

  use, intrinsic :: ieee_arithmetic

  IMPLICIT NONE

  INTEGER ierr

  INTEGER, INTENT(IN) :: i            ! index of a node on the field line where the emitted photoelectrons will be placed

  INTEGER, INTENT(IN) :: jsol         ! index of the solar flux energy bin

  REAL(8), INTENT(IN) :: my_Nn_m3     ! density of neutral species to be ionized by the solar EUV

  REAL, INTENT(IN)    :: sigma_cm2       ! cross section of EUV-neutral ionization interaction [cm^2]

  REAL(8), INTENT(IN) :: threshold_energy_eV    ! initial energy of the photoelectrons

  INTEGER, INTENT(INOUT) :: particle_count
  INTEGER(8), INTENT(INOUT) :: step_count   !###
  INTEGER, INTENT(INOUT) :: mirror_count

  REAL(8) optical_depth

  REAL(8) energy_eV    ! initial energy of the photoelectrons

  REAL(8) EUV_production_rate_m3s
  REAL(8) weight_m3s

  REAL(8) v_ms, vpar_ms, vperp_ms

  INTEGER n

  REAL lambda_mid_A         ! middle of EUVAC bin
  REAL subbin_lambda_min_A  ! minimal wavelength of sub-bin
  REAL photon_lambda_A      ! wavelength of photon macroparticle

! the external procedure checks that the point is not in the shadow, so, calculating the optical depth here is always justified

  optical_depth = 1.0d-22 * DBLE( sigma_tot_He_cm2(jsol) * shared_pfl_point(i)%columnar_content_He_m2 + &
                                & sigma_tot_O_cm2(jsol)  * shared_pfl_point(i)%columnar_content_O_m2 + &
                                & sigma_tot_N2_cm2(jsol) * shared_pfl_point(i)%columnar_content_N2_m2 + &
                                & sigma_tot_O2_cm2(jsol) * shared_pfl_point(i)%columnar_content_O2_m2 )    ! 1e-22 because sigma is in 10^-18 cm^2 while columnar content is in m^-2 

  SELECT CASE (jsol)
     CASE (5,6,8,9,11,14,17,18,20,21,24,26,27,28,33,35,36)
! lines

        energy_eV = solar_flux_energy_bin_eV(jsol) - threshold_energy_eV

        IF (energy_eV.LE.0.0_8) RETURN

! in the EUV production rates below, 1e-9 because solar flux is in 10^9 ph cm^-2 s^-1 and cross sections are in 10^-18 cm^2
! neutral densities are expected in m^-3

        EUV_production_rate_m3s = 1.0d-9 * DBLE(solar_flux_phcm2s1(jsol)) * EXP(-optical_depth) * my_Nn_m3 * DBLE(sigma_cm2)

        weight_m3s = (EUV_production_rate_m3s / DBLE(N_cascades_per_energy_bin)) * shared_pfl_point(i)%fts   ! factor 1/N_of_processes applied in process with rank zero after assembling data from all MPI processes

        if (.not.ieee_is_finite(weight_m3s)) then
           print '("Error in START_PHOTOELECTRONS_WITH_GIVEN_ENERGY :: weight_m3s is ",e12.5," see ",5(e12.5))', weight_m3s, EUV_production_rate_m3s, solar_flux_phcm2s1(jsol), optical_depth, my_Nn_m3, sigma_cm2
           CALL MPI_ABORT(MPI_COMM_WORLD, 11, ierr)
        end if

        v_ms = SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)

        DO n = 1, N_cascades_per_energy_bin

! set initial velocity components
!           vpar_ms = v_ms * COS(pi * RAN(I_random_seed_co_pitch)) ! angle relative to the direction of propagation of guiding center
           vpar_ms = MAX(-v_ms, MIN(v_ms, v_ms * (1.0_8 - 2.0_8 * well_random_number()))) ! isotropic over pitch angle
           vperp_ms = SQRT(MAX(0.0_8, v_ms**2 - vpar_ms**2))
  
           CALL PFL_Run_Photoelectron_Cascade( i, &
                                             & vpar_ms, &
                                             & vperp_ms, &
                                             & weight_m3s, &
                                             & particle_count, &
                                             & step_count, &
                                             & mirror_count )
        END DO

     CASE (3)
! continuum bin #3 which is split

        DO n = 1, N_cascades_per_energy_bin

! in the EUV production rates below, 1e-9 because solar flux is in 10^9 ph cm^-2 s^-1 and cross sections are in 10^-18 cm^2
! neutral densities are expected in m^-3

           lambda_mid_A = 50.0 * jsol + 25.0   ! middle of EUVAC bin
           subbin_lambda_min_A = 50.0 * (REAL(jsol) + REAL(Rank_of_process) / REAL(N_of_processes))   ! minimal wavelength of sub-bin

           photon_lambda_A = subbin_lambda_min_A + 50.0 * (REAL(n) - 0.5) / REAL(N_cascades_per_energy_bin * N_of_processes)  ! wavelength of photon macroparticle

           energy_eV = h_Planck_Js * c_ms * 1.0d10 / (e_Cl * photon_lambda_A) - threshold_energy_eV

           IF (energy_eV.LE.0.0_8) CYCLE

           EUV_production_rate_m3s = 1.0d-9 * DBLE(solar_flux_phcm2s1(jsol)) &
                                   & * &
                                   & (photon_lambda_A / lambda_mid_A) &  ! correction ensuring same values of photon flux an photon energy flux for the whole bin
                                   & * &
                                   & EXP(-optical_depth) * my_Nn_m3 * DBLE(sigma_cm2)

           IF (photon_lambda_A.ge.lambda_mid_A) THEN
              EUV_production_rate_m3s = 1.8*EUV_production_rate_m3s   ! 1.5
           ELSE
              EUV_production_rate_m3s = 0.2*EUV_production_rate_m3s   ! 0.5
           END IF

           weight_m3s = (EUV_production_rate_m3s / DBLE(N_cascades_per_energy_bin)) * shared_pfl_point(i)%fts

           if (.not.ieee_is_finite(weight_m3s)) then
              print '("Error in START_PHOTOELECTRONS_WITH_GIVEN_ENERGY :: weight_m3s is ",e12.5," see ",5(e12.5))', weight_m3s, EUV_production_rate_m3s, solar_flux_phcm2s1(jsol), optical_depth, my_Nn_m3, sigma_cm2
              CALL MPI_ABORT(MPI_COMM_WORLD, 11, ierr)
           end if
           
           v_ms = SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)

! set initial velocity components
!           vpar_ms = v_ms * COS(pi * RAN(I_random_seed_co_pitch)) ! angle relative to the direction of propagation of guiding center
           vpar_ms = MAX(-v_ms, MIN(v_ms, v_ms * (1.0_8 - 2.0_8 * well_random_number()))) ! angle relative to the direction of propagation of guiding center
           vperp_ms = SQRT(MAX(0.0_8, v_ms**2 - vpar_ms**2))
  
           CALL PFL_Run_Photoelectron_Cascade( i, &
                                             & vpar_ms, &
                                             & vperp_ms, &
                                             & weight_m3s, &
                                             & particle_count, &
                                             & step_count, &
                                             & mirror_count )

        END DO

     CASE (1,2,4,7,10,12,13,15,16,19,22,23,25,29:32,34,37)
! other continuum bins

        DO n = 1, N_cascades_per_energy_bin

! in the EUV production rates below, 1e-9 because solar flux is in 10^9 ph cm^-2 s^-1 and cross sections are in 10^-18 cm^2
! neutral densities are expected in m^-3

           lambda_mid_A = 50.0 * jsol + 25.0   ! middle of EUVAC bin
           subbin_lambda_min_A = 50.0 * (REAL(jsol) + REAL(Rank_of_process) / REAL(N_of_processes))   ! minimal wavelength of sub-bin

           photon_lambda_A = subbin_lambda_min_A + 50.0 * (REAL(n) - 0.5) / REAL(N_cascades_per_energy_bin * N_of_processes)  ! wavelength of photon macroparticle

           energy_eV = h_Planck_Js * c_ms * 1.0d10 / (e_Cl * photon_lambda_A) - threshold_energy_eV

           IF (energy_eV.LE.0.0_8) CYCLE

           EUV_production_rate_m3s = 1.0d-9 * DBLE(solar_flux_phcm2s1(jsol)) &
                                   & * &
                                   & (photon_lambda_A / lambda_mid_A) &  ! correction ensuring same values of photon flux and photon energy flux for the whole bin
                                   & * &
                                   & EXP(-optical_depth) * my_Nn_m3 * DBLE(sigma_cm2)

           weight_m3s = (EUV_production_rate_m3s / DBLE(N_cascades_per_energy_bin)) * shared_pfl_point(i)%fts

           if (.not.ieee_is_finite(weight_m3s)) then
              print '("Error in START_PHOTOELECTRONS_WITH_GIVEN_ENERGY :: weight_m3s is ",e12.5," see ",5(e12.5))', weight_m3s, EUV_production_rate_m3s, solar_flux_phcm2s1(jsol), optical_depth, my_Nn_m3, sigma_cm2
              CALL MPI_ABORT(MPI_COMM_WORLD, 11, ierr)
           end if
           
           v_ms = SQRT(2.0_8 * energy_eV * e_Cl / m_e_kg)

! set initial velocity components
!           vpar_ms = v_ms * COS(pi * RAN(I_random_seed_co_pitch)) ! angle relative to the direction of propagation of guiding center
           vpar_ms = MAX(-v_ms, MIN(v_ms, v_ms * (1.0_8 - 2.0_8 * well_random_number()))) ! angle relative to the direction of propagation of guiding center
           vperp_ms = SQRT(MAX(0.0_8, v_ms**2 - vpar_ms**2))
  
           CALL PFL_Run_Photoelectron_Cascade( i, &
                                             & vpar_ms, &
                                             & vperp_ms, &
                                             & weight_m3s, &
                                             & particle_count, &
                                             & step_count, &
                                             & mirror_count )

        END DO

     END SELECT

END SUBROUTINE START_PFL_PHOTOELECTRONS_WITH_GIVEN_ENERGY


!----------------------------------------------------------
!
SUBROUTINE PFL_Run_Photoelectron_Cascade( i_org, &
                                        & org_vpar_ms, &
                                        & org_vperp_ms, &
                                        & org_weight_m3s, &
                                        & particle_count, &
                                        & step_count, &
                                        & mirror_count )

  USE ParallelOperationValues
  USE Photoelectrons
  USE field_line
  USE PhysicalConstants

  USE MPI

  use, intrinsic :: ieee_arithmetic

  IMPLICIT NONE
  
  INTEGER ierr

  INTEGER, INTENT(IN) :: i_org         ! index of point on this field line where the first photoelectron comes from
  REAL(8), INTENT(IN) :: org_vpar_ms
  REAL(8), INTENT(IN) :: org_vperp_ms
  REAL(8), INTENT(IN) :: org_weight_m3s

  INTEGER, INTENT(INOUT) :: particle_count
  INTEGER(8), INTENT(INOUT) :: step_count    !#####
  INTEGER, INTENT(INOUT) :: mirror_count

  INTEGER, PARAMETER :: max_N_of_mirrors = 40

  TYPE electron
     REAL(8) L_m
     REAL(8) vpar_ms
     REAL(8) vperp_ms
     REAL(8) magmom
     INTEGER icell
  end type electron
 
  integer, parameter :: max_N_cascade = 1000       ! maximal expected number of electrons createad by one photoelectron (includes the photoelectron)
  type(electron), allocatable :: cascade_part(:)   ! note: maximal photon energy considered is 6199.25 eV 
                                                   ! which for ionization threshold of 13 eV for O2->O2+ gives 476 ionizations
                                                   ! (13 eV is the minimal energy where cross section for e-O2 collision is given
                                                   !  it is slightly higher than the photoionization threshold value of 12.08 eV)  
  INTEGER ALLOC_ERR

  integer N_cascade             ! actual number of particles in cascade, may vary as particles make new collisions
  integer k

  INTEGER N_of_mirrors

  REAL(8) cp_L_m
  REAL(8) cp_vpar_ms, cp_vperp_ms, cp_magmom
  REAL(8) cp_energy_eV
  INTEGER cell

  REAL(8) L_prev_step_m
  INTEGER cell0

  REAL(8) a_i, b_i    ! interpolation coefficients

  REAL(8) local_geoBT, local_geoBT0

  REAL(8) new_volume
  real(8) recent_dwdt_eVs
  real(8) recent_max_en_coll_frequency_s1

  REAL(8) therm_energy_eV

  real(8) Courant_dt_s
  real(8) Coll_dt_s
  real(8) Coulomb_dt_s
  real(8) photo_delta_t_s

  INTEGER istart, iend, ii

  real(8) dw_eV, alfa

  logical ionization_occurred
  real(8) prod_energy_eV, prod_vpar_ms, prod_vperp_ms

  INTEGER n_vpar_bin, n_vperp_bin

  real(8), allocatable :: local_min_pitch_angle_range(:)
  real(8), allocatable :: local_max_pitch_angle_range(:)
  logical, allocatable :: particle_previously_belonged_to_range(:)
  logical, allocatable :: particle_now_belongs_to_range(:)
  logical particle_collided_with_neutral

  real(8) weight_delta_t_s

  real(8) cp_absv_ms
  real(8) pitch_angle

  real(8) weight_cellm1
  real(8) weight_cell
  real(8) weight_cellp1
  real(8) weight_cellp2
  
  real(8) weight_cellm1_delta_t_s
  real(8) weight_cell_delta_t_s
  real(8) weight_cellp1_delta_t_s
  real(8) weight_cellp2_delta_t_s

  real(8) L_cellm12_m
  real(8) L_cellp12_m
  real(8) L_cellp32_m

  integer nr

! functions

  REAL(8) fl_dwdt_eVs
  REAL(8) fl_Max_frequency_of_en_collisions_s1

  if (N_ranges.gt.0) then
     allocate(local_min_pitch_angle_range(1:N_ranges), stat = alloc_err)
     allocate(local_max_pitch_angle_range(1:N_ranges), stat = alloc_err)
     allocate(particle_previously_belonged_to_range(1:N_ranges), stat = alloc_err)
     allocate(particle_now_belongs_to_range(1:N_ranges), stat = alloc_err)
  end if

  ALLOCATE(cascade_part(1:max_N_cascade), STAT = ALLOC_ERR)

! set initial coordinate, velocities, magnetic moment, and cell containing the particle
  cascade_part(1)%L_m      = shared_pfl_point(i_org)%L_m
  cascade_part(1)%vpar_ms  = org_vpar_ms
  cascade_part(1)%vperp_ms = org_vperp_ms
  cascade_part(1)%magmom =  org_vperp_ms**2 / shared_pfl_point(i_org)%geoB_T
  cascade_part(1)%icell = MIN(i_org, pfl_N_of_points-1)


! the cascade is formed by 
! a [primary] photoelectron and 
! all secondary electrons produced by this [primary] photoelectron and
! all tertiary electrons produced by these secondary electrons and 
! so on

! assume that initially there is only one particle in the cascade
  N_cascade = 1

  k=1  ! particle number in the cascade
  DO WHILE (k.LE.N_cascade)   ! N_cascade can change inside the cycle

     N_of_mirrors = 0

! initiate run-time variables with values stored in the cascade_part array
     cp_L_m   = cascade_part(k)%L_m
     cp_vpar_ms  = cascade_part(k)%vpar_ms
     cp_vperp_ms = cascade_part(k)%vperp_ms
     cp_magmom   = cascade_part(k)%magmom
     cp_energy_eV  = m_e_kg * (cp_vpar_ms**2 + cp_vperp_ms**2) * 0.5_8 / e_Cl

     cell = cascade_part(k)%icell

     if (cell.ge.pfl_N_of_points) then
        print '("proc ",i4," PFL_Run_Photoelectron_Cascade :: error-1 ::: cell too big :: ",i4,"  vs ",i4)', Rank_of_process, cell, pfl_N_of_points
        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
     end if

     cell0 = cell
     L_prev_step_m  = cp_L_m

! prepare interpolation

     a_i = (shared_pfl_point(cell+1)%L_m - cp_L_m) / (shared_pfl_point(cell+1)%L_m - shared_pfl_point(cell)%L_m)
     b_i = (1.0_8 - a_i)

     local_geoBT0 = shared_pfl_point(cell)%geoB_T * a_i + shared_pfl_point(cell+1)%geoB_T * b_i
     new_volume   = shared_pfl_point(cell)%fts   * a_i + shared_pfl_point(cell+1)%fts   * b_i

     recent_dwdt_eVs = fl_dwdt_eVs( cp_energy_eV, &
                                  & MAX(300.0_8, shared_pfl_point(cell)%Te_K  * a_i + shared_pfl_point(cell+1)%Te_K  * b_i), &
                                  & MAX(  0.0_8, shared_pfl_point(cell)%Ne_m3 * a_i + shared_pfl_point(cell+1)%Ne_m3 * b_i) )

     recent_max_en_coll_frequency_s1 = fl_Max_frequency_of_en_collisions_s1( cp_energy_eV, &
                                                                           & MAX(0.0_8, shared_pfl_point(cell)%Nn_He_m3 * a_i  + shared_pfl_point(cell+1)%Nn_He_m3 * b_i), &
                                                                           & MAX(0.0_8, shared_pfl_point(cell)%Nn_O_m3  * a_i  + shared_pfl_point(cell+1)%Nn_O_m3  * b_i), &
                                                                           & MAX(0.0_8, shared_pfl_point(cell)%Nn_N2_m3 * a_i  + shared_pfl_point(cell+1)%Nn_N2_m3 * b_i), &
                                                                           & MAX(0.0_8, shared_pfl_point(cell)%Nn_O2_m3 * a_i  + shared_pfl_point(cell+1)%Nn_O2_m3 * b_i) )

     cp_absv_ms = sqrt(cp_vpar_ms**2 + cp_vperp_ms**2)
     pitch_angle = ACOS(cp_vpar_ms / cp_absv_ms)

     do nr = 1, N_ranges
        call get_local_pitch_angle_range( local_min_pitch_angle_range(nr), &
                                        & local_max_pitch_angle_range(nr), &
                                        & local_geoBT0, &
                                        & min_pitch_angle_range(nr), &
                                        & max_pitch_angle_range(nr), &
                                        & shared_pfl_point(pfl_i_spacecraft)%geoB_T )
     end do

     if (k.eq.1) then
! if the primary particle belongs to an energy range, update source due to EUV
        do nr = 1, N_ranges
           particle_previously_belonged_to_range(nr) = .false.
           if (pitch_angle.lt.local_min_pitch_angle_range(nr)) cycle
           if (pitch_angle.ge.local_max_pitch_angle_range(nr)) cycle
           if (cp_energy_eV.lt.min_energy_eV_range(nr)) cycle
           if (cp_energy_eV.ge.max_energy_eV_range(nr)) cycle
           shared_pfl_point(i_org)%source_EUV_m3s_range(nr) = shared_pfl_point(i_org)%source_EUV_m3s_range(nr) + org_weight_m3s / shared_pfl_point(i_org)%fts
           particle_previously_belonged_to_range(nr) = .true.
        end do
     else
! if the secondary particle belongs to an energy range, update source due to ionization by photoelectrons
        do nr = 1, N_ranges
           particle_previously_belonged_to_range(nr) = .false.
           if (pitch_angle.lt.local_min_pitch_angle_range(nr)) cycle
           if (pitch_angle.ge.local_max_pitch_angle_range(nr)) cycle
           if (cp_energy_eV.lt.min_energy_eV_range(nr)) cycle
           if (cp_energy_eV.ge.max_energy_eV_range(nr)) cycle
! this is NOT the only time when the linear interpolation between two nodes is used
! sources and sinks due t ocollisions with neutrals and interaction with ambient plasma also use the linear interpolation
           shared_pfl_point(cell  )%source_photoe_m3s_range(nr) = shared_pfl_point(cell  )%source_photoe_m3s_range(nr) + (org_weight_m3s / new_volume) * a_i
           shared_pfl_point(cell+1)%source_photoe_m3s_range(nr) = shared_pfl_point(cell+1)%source_photoe_m3s_range(nr) + (org_weight_m3s / new_volume) * b_i
           particle_previously_belonged_to_range(nr) = .true.
        end do
     end if

     DO 

        step_count = step_count + 1_8

! compare particle energy with thermal energy of surrounding plasma ######################################################

        therm_energy_eV = (shared_pfl_point(cell)%Te_K * a_i + shared_pfl_point(cell+1)%Te_K * b_i) * 2.0_8 * k_JK / e_Cl
        IF (cp_energy_eV.LT.therm_energy_eV) EXIT  ! abandon particle

! set timestep ###########################################################################################################

        Courant_dt_s = (shared_pfl_point(cell+1)%L_m - shared_pfl_point(cell)%L_m) / MAX(1.0d5, ABS(cp_vpar_ms))

        Coll_dt_s = Courant_dt_s
        Coulomb_dt_s = Courant_dt_s
        IF (recent_max_en_coll_frequency_s1.GT.0.0_8) Coll_dt_s = 0.2_8 / recent_max_en_coll_frequency_s1   ! was 0.1
        IF (recent_dwdt_eVs.GT.0.0_8) Coulomb_dt_s = cp_energy_eV / recent_dwdt_eVs          
        photo_delta_t_s = MIN(Coulomb_dt_s, Coll_dt_s, Courant_dt_s)

! advance particle ########################################################################################################

        cp_L_m = cp_L_m + photo_delta_t_s * cp_vpar_ms

! check if particle escaped ############################################################################################### update particle and energy losses through the ends

        IF (cp_L_m.GT.shared_pfl_point(pfl_N_of_points)%L_m) EXIT  ! abandon particle
        IF (cp_L_m.LT.shared_pfl_point(1)%L_m) EXIT                !

! update cell ###############################################################################################################

        IF (cp_L_m.GT.L_prev_step_m) THEN
           istart = cell
           iend = pfl_N_of_points
           DO ii = istart+1, iend
              IF (shared_pfl_point(ii)%L_m.GE.cp_L_m) THEN
                 cell = ii-1
                 EXIT
              END IF
           END DO
        ELSE              ! i.e. cp_L_m.le.L_prev_step_m
           istart = cell
           iend = 1
           DO ii = istart, iend, -1
              IF (shared_pfl_point(ii)%L_m.LE.cp_L_m) THEN
                 cell = ii
                 EXIT
              END IF
           END DO
        END IF
! fool proof
        if ((cell.ge.pfl_N_of_points).or.(cell.lt.1)) then
           print '("proc ",i4," Error-2 in PFL_Run_Photoelectron_Cascade :: cell too big/small :: ",i6,"  vs ",i6)', Rank_of_process, cell, pfl_N_of_points
           call MPI_ABORT(MPI_COMM_WORLD, 2, ierr)
        end if
        IF ( (cp_L_m.LT.shared_pfl_point(cell)%L_m).OR.(cp_L_m.GT.shared_pfl_point(cell+1)%L_m) ) THEN
           PRINT '("Process ",i4," Error-3 in PFL_Run_Photoelectron_Cascade :: illegal cell ",3(2x,i5)," Lm: ",3(2x,e14.7)," origin: ",3(2x,e14.7))', &
                & Rank_of_process, istart, cell, iend, shared_pfl_point(cell)%L_m, cp_L_m, shared_pfl_point(cell+1)%L_m, cp_L_m - photo_delta_t_s * cp_vpar_ms, photo_delta_t_s, cp_vpar_ms
           call MPI_ABORT(MPI_COMM_WORLD, 3, ierr)
        END IF

! check if the mirroring occurred ###########################################################################################

! update interpolation coefficients
        a_i = (shared_pfl_point(cell+1)%L_m - cp_L_m) / (shared_pfl_point(cell+1)%L_m - shared_pfl_point(cell)%L_m) 
        b_i = (1.0_8 - a_i)
        if ((a_i.lt.0.0_8).or.(b_i.lt.0.0_8)) then
           print '("proc ",i4," Error-4 in PFL_Run_Photoelectron_Cascade :: ",2(2x,e12.5),2x,i5)', Rank_of_process, a_i, b_i, cell
           call MPI_ABORT(MPI_COMM_WORLD, 4, ierr)
        end if
! modify particle velocity according to the magnetic momentum conservation (apply this before applying collisions)
        local_geoBT = shared_pfl_point(cell)%geoB_T * a_i + shared_pfl_point(cell+1)%geoB_T * b_i
        new_volume   = shared_pfl_point(cell)%fts  * a_i + shared_pfl_point(cell+1)%fts   * b_i
        IF ((cp_magmom * local_geoBT).GE.(2.0_8 * cp_energy_eV * e_Cl / m_e_kg)) THEN
! the particle bounced off the mirror point
           cp_vpar_ms = -cp_vpar_ms
! return it to the initial point of this timestep, restore all related values
           cp_L_m = L_prev_step_m
           local_geoBT = local_geoBT0
           cell = cell0
           a_i = (shared_pfl_point(cell+1)%L_m - cp_L_m) / (shared_pfl_point(cell+1)%L_m - shared_pfl_point(cell)%L_m)
           b_i = (1.0_8 - a_i)   
           new_volume   = shared_pfl_point(cell)%fts * a_i + shared_pfl_point(cell+1)%fts * b_i
! update mirror counters
           mirror_count = mirror_count + 1          ! collective for all particles
           N_of_mirrors = N_of_mirrors + 1          ! individual for this particle
        ELSE
! direction of motion does not change
           cp_vperp_ms = SQRT(cp_magmom * local_geoBT)
           IF (cp_vpar_ms.GE.0.0_8) THEN
              cp_vpar_ms =  SQRT(MAX(0.0_8, 2.0_8 * cp_energy_eV * e_Cl / m_e_kg - cp_magmom * local_geoBT))
           ELSE
              cp_vpar_ms = -SQRT(MAX(0.0_8, 2.0_8 * cp_energy_eV * e_Cl / m_e_kg - cp_magmom * local_geoBT))
           END IF
        END IF

        if (.not.ieee_is_finite(cp_vpar_ms)) then
           print '("Error-5 (NAN) in PFL_Run_Photoelectron_Cascade :: ",6(2x,e14.7))', cp_vpar_ms, cp_vperp_ms, cp_magmom, local_geoBT, a_i, b_i
           call MPI_ABORT(MPI_COMM_WORLD, 5, ierr)
        end if

        if (.not.ieee_is_finite(cp_vperp_ms)) then
           print '("Error-6 (NAN) in PFL_Run_Photoelectron_Cascade :: ",6(2x,e14.7))', cp_vpar_ms, cp_vperp_ms, cp_magmom, local_geoBT, a_i, b_i
           call MPI_ABORT(MPI_COMM_WORLD, 6, ierr)
        end if

! calculate energy decrement due to Coulomb collisions with ambient electrons ###############################################

        recent_dwdt_eVs = fl_dwdt_eVs( cp_energy_eV, &
                                     & MAX(300.0_8, shared_pfl_point(cell)%Te_K  * a_i + shared_pfl_point(cell+1)%Te_K  * b_i), &
                                     & MAX(  0.0_8, shared_pfl_point(cell)%Ne_m3 * a_i + shared_pfl_point(cell+1)%Ne_m3 * b_i) )
        dw_eV = MIN(cp_energy_eV, photo_delta_t_s * recent_dwdt_eVs)
! adjust particle energy and magnetic moment according to the decrement above
        cp_energy_eV = MAX(cp_energy_eV - dw_eV, 0.0_8)
        alfa = SQRT(cp_energy_eV/(cp_energy_eV + dw_eV))
        cp_vpar_ms = cp_vpar_ms  * alfa
        cp_vperp_ms = cp_vperp_ms * alfa

        if (.not.ieee_is_finite(recent_dwdt_eVs+dw_eV+cp_energy_eV+alfa+cp_vpar_ms+cp_vperp_ms)) then
           print '("Error-7 (NAN) in PFL_Run_Photoelectron_Cascade :: ",6(2x,e14.7))', recent_dwdt_eVs, dw_eV, cp_energy_eV, alfa, cp_vpar_ms, cp_vperp_ms
           call MPI_ABORT(MPI_COMM_WORLD, 7, ierr)
        end if

! perform Coulomb scattering #################################################################################################

        CALL FL_Coulomb_Scatter_Electron( MAX(0.0_8, shared_pfl_point(cell)%Ne_m3 * a_i + shared_pfl_point(cell+1)%Ne_m3 * b_i), &
                                        & MAX(300.0_8, shared_pfl_point(cell)%Te_K * a_i + shared_pfl_point(cell+1)%Te_K * b_i), &
                                        & cp_vpar_ms, &
                                        & cp_vperp_ms, &
                                        & photo_delta_t_s ) !, ksi)

! perform collisions with neutrals if necessary ##############################################################################

        CALL FL_Collide_Electron( MAX(0.0_8, shared_pfl_point(cell)%Nn_He_m3 * a_i + shared_pfl_point(cell+1)%Nn_He_m3 * b_i), &
                                & MAX(0.0_8, shared_pfl_point(cell)%Nn_O_m3  * a_i + shared_pfl_point(cell+1)%Nn_O_m3  * b_i), &
                                & MAX(0.0_8, shared_pfl_point(cell)%Nn_N2_m3 * a_i + shared_pfl_point(cell+1)%Nn_N2_m3 * b_i), &
                                & MAX(0.0_8, shared_pfl_point(cell)%Nn_O2_m3 * a_i + shared_pfl_point(cell+1)%Nn_O2_m3 * b_i), &
                                & MAX(300.0_8, shared_pfl_point(cell)%Te_K * a_i + shared_pfl_point(cell+1)%Te_K * b_i), & 
                                & cp_energy_eV, &
                                & cp_vpar_ms, &
                                & cp_vperp_ms, &
                                & recent_max_en_coll_frequency_s1, &
                                & ionization_occurred, &
                                & prod_energy_eV, &
                                & prod_vpar_ms, &
                                & prod_vperp_ms, &
                                & photo_delta_t_s, &
                                & org_weight_m3s / new_volume, &
                                & a_i, &
                                & b_i, &
                                & shared_pfl_point(cell)%rate_iN2_production_m3s, &
                                & shared_pfl_point(cell)%rate_iN_production_m3s, &
                                & shared_pfl_point(cell)%rate_iO2_production_m3s, &
                                & shared_pfl_point(cell)%rate_iO_production_m3s, &
                                & shared_pfl_point(cell)%rate_iHe_production_m3s, &
                                & shared_pfl_point(cell+1)%rate_iN2_production_m3s, &
                                & shared_pfl_point(cell+1)%rate_iN_production_m3s, &
                                & shared_pfl_point(cell+1)%rate_iO2_production_m3s, &
                                & shared_pfl_point(cell+1)%rate_iO_production_m3s, &
                                & shared_pfl_point(cell+1)%rate_iHe_production_m3s, &
                                & particle_collided_with_neutral )

! recalculate the magnetic moment
        cp_magmom = cp_vperp_ms**2 / local_geoBT

        if (.not.ieee_is_finite(cp_magmom)) then
           print '("proc ",i4," Error-8 (NAN) in PFL_Run_Photoelectron_Cascade :: ",4(2x,e14.7))', Rank_of_process, cp_vpar_ms, cp_vperp_ms, cp_magmom, local_geoBT
           call MPI_ABORT(MPI_COMM_WORLD, 8, ierr)
        end if

        cp_absv_ms = sqrt(cp_vpar_ms**2 + cp_vperp_ms**2)
        pitch_angle = ACOS(cp_vpar_ms / cp_absv_ms)

! update pitch angle limits
        do nr = 1, N_ranges
           call get_local_pitch_angle_range( local_min_pitch_angle_range(nr), &
                                           & local_max_pitch_angle_range(nr), &
                                           & local_geoBT, &
                                           & min_pitch_angle_range(nr), &
                                           & max_pitch_angle_range(nr), &
                                           & shared_pfl_point(pfl_i_spacecraft)%geoB_T )
        end do

        do nr = 1, N_ranges
           particle_now_belongs_to_range(nr) = .false.
           if (pitch_angle.lt.local_min_pitch_angle_range(nr)) cycle
           if (pitch_angle.ge.local_max_pitch_angle_range(nr)) cycle
           if (cp_energy_eV.lt.min_energy_eV_range(nr)) cycle
           if (cp_energy_eV.ge.max_energy_eV_range(nr)) cycle
           particle_now_belongs_to_range(nr) = .true.
        end do

! if particle's end position is between nodes i=cell and i=cell+1 its start position can be between nodes i=cell-1 and i=cell+2
! so there are potentially 4 nodes which may receive a contribution from this particle

! find the weights for the 4 nodes

        weight_cellm1 = 0.0_8  ! for source/sink terms
        weight_cell   = 0.0_8
        weight_cellp1 = 0.0_8
        weight_cellp2 = 0.0_8

        weight_cellm1_delta_t_s = 0.0_8 ! for fluxes
        weight_cell_delta_t_s   = 0.0_8
        weight_cellp1_delta_t_s = 0.0_8
        weight_cellp2_delta_t_s = 0.0_8

        if (cell.eq.1) then
           L_cellm12_m = shared_pfl_point(1)%L_m - 0.5_8 * (shared_pfl_point(2)%L_m - shared_pfl_point(1)%L_m)
        else
           L_cellm12_m = 0.5_8 * (shared_pfl_point(cell-1)%L_m + shared_pfl_point(cell)%L_m)
        end if

        L_cellp12_m = 0.5_8 * (shared_pfl_point(cell)%L_m + shared_pfl_point(cell+1)%L_m)

        if (cell.eq.(pfl_N_of_points-1)) then
           L_cellp32_m = shared_pfl_point(pfl_N_of_points)%L_m + 0.5_8 * (shared_pfl_point(pfl_N_of_points)%L_m - shared_pfl_point(pfl_N_of_points-1)%L_m)
        else
           L_cellp32_m = 0.5_8 * (shared_pfl_point(cell+1)%L_m + shared_pfl_point(cell+2)%L_m)
        end if

! it is expected that shared_pfl_point(cell)%L_m <= cp_L_m <= shared_pfl_point(cell+1)%L_m

        if (cp_L_m.gt.L_prev_step_m) then
! particle moves in the positive L-direction (Northward)

           weight_cellm1 = MAX(0.0_8, L_cellm12_m - L_prev_step_m) / (cp_L_m - L_prev_step_m)                              ! zero if L_prev_step_m > L_cellm12_m
           
           weight_cell   = MAX(0.0_8, MIN(cp_L_m,L_cellp12_m) - MAX(L_prev_step_m,L_cellm12_m)) / (cp_L_m - L_prev_step_m) ! zero if L_prev_step_m > L_cellp12_m

           weight_cellp1 = MAX(0.0_8, cp_L_m                  - MAX(L_prev_step_m,L_cellp12_m)) / (cp_L_m - L_prev_step_m) ! zero if cp_L_m < L_cellp12_m
           
           weight_cellm1_delta_t_s = photo_delta_t_s * weight_cellm1
           weight_cell_delta_t_s   = photo_delta_t_s * weight_cell
           weight_cellp1_delta_t_s = photo_delta_t_s * weight_cellp1
           
        else if (cp_L_m.lt.L_prev_step_m) then
! particle moves in the negative L-direction (Southward)

           weight_cell   = MAX(0.0_8, MIN(L_prev_step_m,L_cellp12_m) -                  cp_L_m) / (L_prev_step_m - cp_L_m) ! zero if cp_L_m > L_cellp12_m
              
           weight_cellp1 = MAX(0.0_8, MIN(L_prev_step_m,L_cellp32_m) - MAX(cp_L_m,L_cellp12_m)) / (L_prev_step_m - cp_L_m) ! zero if L_prev_step_m < L_cellp12_m
           
           weight_cellp2 = MAX(0.0_8, L_prev_step_m - L_cellp32_m) / (L_prev_step_m - cp_L_m)                              ! zero if L_prev_step_m < L_cellp32_m

           weight_cell_delta_t_s   = photo_delta_t_s * weight_cell   
           weight_cellp1_delta_t_s = photo_delta_t_s * weight_cellp1
           weight_cellp2_delta_t_s = photo_delta_t_s * weight_cellp2

        else
! exceptional case, particle did not move

           if (cp_L_m.lt.L_cellp12_m) then
              weight_cell = 1.0_8
              weight_cell_delta_t_s   = photo_delta_t_s 
           else if (cp_L_m.gt.L_cellp12_m) then
              weight_cellp1 = 1.0_8
              weight_cellp1_delta_t_s = photo_delta_t_s 
           else
              weight_cell   = 0.5_8
              weight_cellp1 = 0.5_8
              weight_cell_delta_t_s   = 0.5_8 * photo_delta_t_s 
              weight_cellp1_delta_t_s = 0.5_8 * photo_delta_t_s 
           end if

        end if

! update range sources / sinks
        do nr = 1, N_ranges
           if (particle_previously_belonged_to_range(nr)) then
              if (.not.particle_now_belongs_to_range(nr)) then
! particle belonged to the range at the previous step
! but it has left the range
! update loss terms
                 if (particle_collided_with_neutral) then
                    shared_pfl_point(cell  )%sink_neutrals_m3s_range(nr) = shared_pfl_point(cell  )%sink_neutrals_m3s_range(nr) + a_i * (org_weight_m3s/new_volume)
                    shared_pfl_point(cell+1)%sink_neutrals_m3s_range(nr) = shared_pfl_point(cell+1)%sink_neutrals_m3s_range(nr) + b_i * (org_weight_m3s/new_volume)
                 else
                    shared_pfl_point(cell  )%sink_plasma_m3s_range(nr) = shared_pfl_point(cell  )%sink_plasma_m3s_range(nr) + a_i * (org_weight_m3s/new_volume)
                    shared_pfl_point(cell+1)%sink_plasma_m3s_range(nr) = shared_pfl_point(cell+1)%sink_plasma_m3s_range(nr) + b_i * (org_weight_m3s/new_volume)
                 end if
              end if
           else
              if (particle_now_belongs_to_range(nr)) then
! particle did not belong to the range at the previous step
! but now it belongs to the range
! update source terms
                 if (particle_collided_with_neutral) then
                    shared_pfl_point(cell  )%source_neutrals_m3s_range(nr) = shared_pfl_point(cell  )%source_neutrals_m3s_range(nr) + a_i * (org_weight_m3s/new_volume)
                    shared_pfl_point(cell+1)%source_neutrals_m3s_range(nr) = shared_pfl_point(cell+1)%source_neutrals_m3s_range(nr) + b_i * (org_weight_m3s/new_volume)
                 else
                    shared_pfl_point(cell  )%source_plasma_m3s_range(nr) = shared_pfl_point(cell  )%source_plasma_m3s_range(nr) + a_i * (org_weight_m3s/new_volume)
                    shared_pfl_point(cell+1)%source_plasma_m3s_range(nr) = shared_pfl_point(cell+1)%source_plasma_m3s_range(nr) + b_i * (org_weight_m3s/new_volume)
                 end if
              end if
           end if
        end do  !###   do nr = 1, N_ranges

! update range fluxes
        do nr = 1, N_ranges
           if (.not.particle_now_belongs_to_range(nr)) cycle
           if (cell.gt.1)                   shared_pfl_point(cell-1)%Ge_m2s_range(nr) = shared_pfl_point(cell-1)%Ge_m2s_range(nr) + weight_cellm1_delta_t_s * cp_absv_ms * (org_weight_m3s/new_volume)
                                            shared_pfl_point(cell  )%Ge_m2s_range(nr) = shared_pfl_point(cell  )%Ge_m2s_range(nr) + weight_cell_delta_t_s   * cp_absv_ms * (org_weight_m3s/new_volume)
                                            shared_pfl_point(cell+1)%Ge_m2s_range(nr) = shared_pfl_point(cell+1)%Ge_m2s_range(nr) + weight_cellp1_delta_t_s * cp_absv_ms * (org_weight_m3s/new_volume)
           if (cell.lt.(pfl_N_of_points-1)) shared_pfl_point(cell+2)%Ge_m2s_range(nr) = shared_pfl_point(cell+2)%Ge_m2s_range(nr) + weight_cellp2_delta_t_s * cp_absv_ms * (org_weight_m3s/new_volume)

           if (cell.gt.1)                   shared_pfl_point(cell-1)%Ge_par_m2s_range(nr) = shared_pfl_point(cell-1)%Ge_par_m2s_range(nr) + weight_cellm1_delta_t_s * cp_vpar_ms * (org_weight_m3s/new_volume)
                                            shared_pfl_point(cell  )%Ge_par_m2s_range(nr) = shared_pfl_point(cell  )%Ge_par_m2s_range(nr) + weight_cell_delta_t_s   * cp_vpar_ms * (org_weight_m3s/new_volume)
                                            shared_pfl_point(cell+1)%Ge_par_m2s_range(nr) = shared_pfl_point(cell+1)%Ge_par_m2s_range(nr) + weight_cellp1_delta_t_s * cp_vpar_ms * (org_weight_m3s/new_volume)
           if (cell.lt.(pfl_N_of_points-1)) shared_pfl_point(cell+2)%Ge_par_m2s_range(nr) = shared_pfl_point(cell+2)%Ge_par_m2s_range(nr) + weight_cellp2_delta_t_s * cp_vpar_ms * (org_weight_m3s/new_volume)
        end do   !###   do nr = 1, N_ranges

! calculate contributions to the plus/minus fluxes and density #########################################################################

        IF (cp_vpar_ms.GT.0.0_8) THEN

!           shared_pfl_point(cell  )%Ge_plus_m2s = shared_pfl_point(cell  )%Ge_plus_m2s + photo_delta_t_s * cp_vpar_ms * (org_weight_m3s / new_volume) * a_i
!           shared_pfl_point(cell+1)%Ge_plus_m2s = shared_pfl_point(cell+1)%Ge_plus_m2s + photo_delta_t_s * cp_vpar_ms * (org_weight_m3s / new_volume) * b_i

           if (cell.gt.1)                   shared_pfl_point(cell-1)%Ge_plus_m2s = shared_pfl_point(cell-1)%Ge_plus_m2s + weight_cellm1_delta_t_s * cp_vpar_ms * (org_weight_m3s/new_volume)
                                            shared_pfl_point(cell  )%Ge_plus_m2s = shared_pfl_point(cell  )%Ge_plus_m2s + weight_cell_delta_t_s   * cp_vpar_ms * (org_weight_m3s/new_volume)
                                            shared_pfl_point(cell+1)%Ge_plus_m2s = shared_pfl_point(cell+1)%Ge_plus_m2s + weight_cellp1_delta_t_s * cp_vpar_ms * (org_weight_m3s/new_volume)
           if (cell.lt.(pfl_N_of_points-1)) shared_pfl_point(cell+2)%Ge_plus_m2s = shared_pfl_point(cell+2)%Ge_plus_m2s + weight_cellp2_delta_t_s * cp_vpar_ms * (org_weight_m3s/new_volume)

        ELSE

!           shared_pfl_point(cell  )%Ge_minus_m2s = shared_pfl_point(cell  )%Ge_minus_m2s + photo_delta_t_s * cp_vpar_ms * (org_weight_m3s / new_volume) * a_i
!           shared_pfl_point(cell+1)%Ge_minus_m2s = shared_pfl_point(cell+1)%Ge_minus_m2s + photo_delta_t_s * cp_vpar_ms * (org_weight_m3s / new_volume) * b_i

           if (cell.gt.1)                   shared_pfl_point(cell-1)%Ge_minus_m2s = shared_pfl_point(cell-1)%Ge_minus_m2s + weight_cellm1_delta_t_s * cp_vpar_ms * (org_weight_m3s/new_volume)
                                            shared_pfl_point(cell  )%Ge_minus_m2s = shared_pfl_point(cell  )%Ge_minus_m2s + weight_cell_delta_t_s   * cp_vpar_ms * (org_weight_m3s/new_volume)
                                            shared_pfl_point(cell+1)%Ge_minus_m2s = shared_pfl_point(cell+1)%Ge_minus_m2s + weight_cellp1_delta_t_s * cp_vpar_ms * (org_weight_m3s/new_volume)
           if (cell.lt.(pfl_N_of_points-1)) shared_pfl_point(cell+2)%Ge_minus_m2s = shared_pfl_point(cell+2)%Ge_minus_m2s + weight_cellp2_delta_t_s * cp_vpar_ms * (org_weight_m3s/new_volume)

        END IF

!        shared_pfl_point(cell  )%kinetic_ne_m3 = shared_pfl_point(cell  )%kinetic_ne_m3 + photo_delta_t_s * (org_weight_m3s/new_volume) * a_i
!        shared_pfl_point(cell+1)%kinetic_ne_m3 = shared_pfl_point(cell+1)%kinetic_ne_m3 + photo_delta_t_s * (org_weight_m3s/new_volume) * b_i

        if (cell.gt.1)                   shared_pfl_point(cell-1)%kinetic_ne_m3 = shared_pfl_point(cell-1)%kinetic_ne_m3 + weight_cellm1_delta_t_s * (org_weight_m3s/new_volume)
                                         shared_pfl_point(cell  )%kinetic_ne_m3 = shared_pfl_point(cell  )%kinetic_ne_m3 + weight_cell_delta_t_s   * (org_weight_m3s/new_volume)
                                         shared_pfl_point(cell+1)%kinetic_ne_m3 = shared_pfl_point(cell+1)%kinetic_ne_m3 + weight_cellp1_delta_t_s * (org_weight_m3s/new_volume)
        if (cell.lt.(pfl_N_of_points-1)) shared_pfl_point(cell+2)%kinetic_ne_m3 = shared_pfl_point(cell+2)%kinetic_ne_m3 + weight_cellp2_delta_t_s * (org_weight_m3s/new_volume)

! add the particle to the electron velocity distribution function in spacecraft location ##############################################################

        weight_delta_t_s = 0.0_8

        if      (pfl_i_spacecraft.eq.(cell-1)) then
           weight_delta_t_s = weight_cellm1_delta_t_s
        else if (pfl_i_spacecraft.eq.(cell)  ) then
           weight_delta_t_s = weight_cell_delta_t_s
        else if (pfl_i_spacecraft.eq.(cell+1)) then
           weight_delta_t_s = weight_cellp1_delta_t_s
        else if (pfl_i_spacecraft.eq.(cell+2)) then
           weight_delta_t_s = weight_cellp2_delta_t_s
        end if

        if (.not.ieee_is_finite(weight_delta_t_s)) then
           print '("Error-9 in PFL_Run_Photoelectron_Cascade :: NaN weight_delta_t_s :",6(2x,e16.9))', cp_L_m, L_prev_step_m, cp_L_m - L_prev_step_m, L_prev_step_m - cp_L_m
           CALL MPI_ABORT(MPI_COMM_WORLD, 9, ierr)
        end if

! add the particle to the electron velocity distribution function ##############################################################
    
        if (weight_delta_t_s.gt.0.0_8) then
           IF (cp_vpar_ms.GE.(-0.5_8 * dvpar_ms)) THEN
              n_vpar_bin = INT(0.5_8 + cp_vpar_ms / dvpar_ms)
           ELSE
              n_vpar_bin = -INT(0.5_8 - cp_vpar_ms / dvpar_ms)
           END IF
           n_vperp_bin = INT(cp_vperp_ms / dvperp_ms)
        
           IF ((ABS(n_vpar_bin).LE.max_evdf_N_vpar).AND.(n_vperp_bin.LE.max_evdf_N_vperp)) THEN 
              evdf(n_vpar_bin, n_vperp_bin) = evdf(n_vpar_bin, n_vperp_bin) + weight_delta_t_s * (org_weight_m3s/new_volume)
           END IF
        end if

        IF (ionization_occurred) THEN
           IF (N_cascade.LT.max_N_cascade) THEN
! include particle produced in the ionization
              N_cascade = N_cascade+1

              cascade_part(N_cascade)%L_m        = cp_L_m
              cascade_part(N_cascade)%vpar_ms    = prod_vpar_ms
              cascade_part(N_cascade)%vperp_ms   = prod_vperp_ms
              cascade_part(N_cascade)%magmom     = prod_vperp_ms**2 / local_geoBT
              cascade_part(N_cascade)%icell      = cell

           ELSE
!error
              print '("proc ",i4," Error-10 in PFL_Run_Photoelectron_Cascade :: too many particles in cascade")', Rank_of_process
              CALL MPI_ABORT(MPI_COMM_WORLD, 10, ierr)
           END IF
        END IF

! update saved initial position and magnetic field in this position
        L_prev_step_m  = cp_L_m
        local_geoBT0 = local_geoBT
        cell0 = cell

        if (N_ranges.gt.0) particle_previously_belonged_to_range = particle_now_belongs_to_range

        IF (N_of_mirrors.GT.max_N_of_mirrors) EXIT   ! abandon particle which bounces too much

     END DO  ! end of cycle which moves one particle

     k = k+1

  END DO     ! end of cycle over all particles in cascade

  particle_count = particle_count + N_cascade

! cleanup
  DEALLOCATE(cascade_part, STAT = ALLOC_ERR)
  if (N_ranges.gt.0) then 
     deallocate(local_min_pitch_angle_range, stat = alloc_err)
     deallocate(local_max_pitch_angle_range, stat = alloc_err)
     deallocate(particle_previously_belonged_to_range, stat = alloc_err)
     deallocate(particle_now_belongs_to_range, stat = alloc_err)
  end if

END SUBROUTINE PFL_Run_Photoelectron_Cascade

!-------------------------------------
! as in Swartz, Nisbet, and Green, JGR, 76, 8425 (1971)
!
REAL(8) FUNCTION fl_dwdt_eVs(w_eV, Te_K, Ne_m3)

  IMPLICIT NONE

  REAL(8) w_eV
  REAL(8) Te_K
  REAL(8) Ne_m3

  IF (w_eV.LE.(Te_K * 8.618d-5)) THEN
     fl_dwdt_eVs = 0.0_8
     RETURN
  END IF

  fl_dwdt_eVs = MAX(0.0_8, (2.0d-4 * (1.0d-6 * Ne_m3)**0.97 / w_eV**0.44) * (MAX(0.0_8, w_eV - Te_K * 8.618d-5)/(w_eV - 0.53_8 * Te_K * 8.618d-5))**2.36)

END FUNCTION fl_dwdt_eVs

!------------------------------
!
module coefs_for_coulomb_scatter_proc
  use PhysicalConstants, ONLY : pi, e_Cl, m_e_kg, eps_0_Fm, k_JK

  real(8), parameter :: k_b0 = e_Cl**2 / (2.0_8 * pi * eps_0_Fm * 0.5_8 * m_e_kg)
  real(8), parameter :: k_lambdaDe = 2.0_8 * k_JK * eps_0_Fm / (e_Cl**2)

end module coefs_for_coulomb_scatter_proc

!---------------------------------------------------
! based on [K. Nanbu, "Theory of cumulative small-angle collisions in plasmas", Phys.Rev.E, vol.55, pp.4642-4652, 1997]
! and the elastic e-n collision procedure from EDIPIC
! the parallel direction is the x-direction
! the transverse direction is in the y-z plane and only the combined perpendicular velocity is calculated
!
SUBROUTINE FL_Coulomb_Scatter_Electron(Ne_m3, Te_K, vpar_ms, vperp_ms, dt_s ) !, ksi)

  USE PhysicalConstants, ONLY : pi, e_Cl, m_e_kg, eps_0_Fm, k_JK
!  use coefs_for_coulomb_scatter_proc
  USE rng_wrapper

  use, intrinsic :: ieee_arithmetic

  implicit none

  real(8), intent(in) :: Ne_m3, Te_K
  real(8), intent(inout) :: vpar_ms, vperp_ms
  real(8), intent(in)    :: dt_s

!real(8), intent(out) :: ksi

  real(8) g_ms, b0_m, lambdaDe_m, s, A, expA, expmA

  real(8) R

  real(8) CosKsi, SinKsi     ! scattering relative to the initial direction
  real(8) SinPhi             ! azimuthal scattering

  real(8) Vx_ms            ! temporary value to store parallel velocity

! function
  real(8) A_Nanbu

!  k_b0 = e_Cl**2 / (2.0_8 * pi * eps_0_Fm * 0.5_8 * m_e_kg)
!  k_lambdaDe = 2.0_8 * k_JK * eps_0_Fm / (e_Cl**2)

!  g_ms = vpar_ms**2 + vperp_ms**2
!  b0_m = k_b0 / g_ms  ! note, here g_ms is the velocity module squared
!  lambdaDe_m = SQRT(k_lambdaDe * Te_K / Ne_m3)
!  g_ms = SQRT(g_ms)

  g_ms = SQRT(vpar_ms**2 + vperp_ms**2)
  if (g_ms.LE.1.0_8) return
  if (Ne_m3.EQ.0.0_8) return

  b0_m = e_Cl**2 / (2.0_8 * pi * eps_0_Fm * 0.5_8 * m_e_kg * g_ms**2)
  lambdaDe_m = SQRT(2.0_8 * Te_K * k_JK * eps_0_Fm / (Ne_m3 * e_Cl**2))
  s = Ne_m3 * g_ms * pi * b0_m**2 * LOG(lambdaDe_m / b0_m) * dt_s

  if (.not.ieee_is_finite(s)) then
     print '("Error : NaN in FL_Coulomb_Scatter_Electron :: ",9(2x,e12.5))', s,  Ne_m3, g_ms, vpar_ms, vperp_ms, b0_m, lambdaDe_m, Te_K, dt_s
     stop
  end if

  IF (s.LE.0.0_8) RETURN

  A = A_Nanbu(s)

! Calculate the scattering angle relative to the initial direction of electron
  R = well_random_number()

  IF (s.LT.0.1_8) THEN

     CosKsi = 1.0_8 + s * LOG(MAX(R,4.53999297624849d-05))

  ELSE IF (s.GT.6.0_8) THEN

     CosKsi = 2.0_8 * R - 1.0_8

  ELSE

     expa = EXP(A)
     expma = EXP(-A)
     CosKsi = LOG(EXP(-A) + 2.0_8 * R * (expa - expma) / (expa + expma)) / A

  END IF

  CosKsi = MAX(MIN(0.999999999999_8, CosKsi), -0.999999999999_8)   !############ to avoid an unlikely situation when |CosKsi|>1
  SinKsi = SQRT(1.0_8 - CosKsi**2)

! Calculate the azimuthal scattering angle
  SinPhi = SIN(6.28318530718_8 * well_random_number())

! Transform the velocity 
  Vx_ms    = vpar_ms  * CosKsi + SinPhi * vperp_ms * SinKsi

  vperp_ms = sqrt(max(0.0_8, vpar_ms**2 + vperp_ms**2 - Vx_ms**2)) !* alpha
  vpar_ms  = Vx_ms !* alpha

!ksi = acos(CosKsi)

end subroutine FL_Coulomb_Scatter_Electron

!---------------------------------------------------
!
real(8) function A_Nanbu(s)

  implicit none

  real(8) s

  real, parameter, dimension(22) :: stab = (/ 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, &
                                         &     0.1, 0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9, &
                                         &     1.0, 2.0,  3.0,  4.0 /)

  real, parameter, dimension(22) :: atab = (/ 100.5,  50.5,   33.84,  25.5,   20.5,   17.17,  14.79,  13.01,  11.62, &
                                           &   10.51,  5.516,  3.845,  2.987,  2.448,  2.067,  1.779,  1.551,  1.363, &
                                           &    1.207, 0.4105, 0.1496, 0.05496 /)
  integer i
  real ai


  real myA_Nanbu
  logical foundme

  foundme=.false.
  if (real(s).lt.stab(1)) then

!     A_Nanbu = 1.0_8 / s
     myA_Nanbu = 1.0_8 / s
     foundme=.true.

  else if (real(s).gt.stab(22)) then

!     A_Nanbu = 3.0_8 * exp(-s)
     myA_Nanbu = 3.0_8 * exp(-s)
     foundme=.true.

  else

     do i = 1, 21
        if ((real(s).ge.stab(i)).and.(real(s).le.stab(i+1))) then
           ai = (stab(i+1)-real(s))/(stab(i+1)-stab(i))
!           A_Nanbu = dble(atab(i) * ai + atab(i+1) * (1.0 - ai))
           myA_Nanbu = dble(atab(i) * ai + atab(i+1) * (1.0 - ai))
           foundme=.true.
           exit
        end if
     end do

  end if

  if (.not.foundme) then
     print '("error in A_Nanbu for s= ",e12.5)', s
     stop
  end if

  A_Nanbu=myA_Nanbu

end function A_Nanbu

!===================================================================================================
!
SUBROUTINE FL_Collide_Electron( Nn_He_m3, &
                              & Nn_O_m3, &
                              & Nn_N2_m3, &
                              & Nn_O2_m3, &
                              & Te_K, &
                              & energy_eV, &
                              & vpar_ms, &
                              & vperp_ms, &
                              & coll_freq_s1, &
                              & ionization_occurred, &
                              & prod_energy_eV, &
                              & prod_vpar_ms, &
                              & prod_vperp_ms, &
                              & delta_t_s, &
                              & weight_m3s, &
                              & a_i, &
                              & b_i, &
                              & point_i_kinetic_chem_rate_N2_ioniz_m3s, &
                              & point_i_kinetic_chem_rate_N_ioniz_m3s, &
                              & point_i_kinetic_chem_rate_O2_ioniz_m3s, &
                              & point_i_kinetic_chem_rate_O_ioniz_m3s, &
                              & point_i_kinetic_chem_rate_He_ioniz_m3s, &
                              & point_ip1_kinetic_chem_rate_N2_ioniz_m3s, &
                              & point_ip1_kinetic_chem_rate_N_ioniz_m3s, &
                              & point_ip1_kinetic_chem_rate_O2_ioniz_m3s, &
                              & point_ip1_kinetic_chem_rate_O_ioniz_m3s, &
                              & point_ip1_kinetic_chem_rate_He_ioniz_m3s, &
                              & particle_collided_with_neutral )

  use ParallelOperationValues
  use PhysicalConstants, only : k_JK, e_Cl

  use Photoelectrons, only : N_colkind, threshold_eV
  USE rng_wrapper

  IMPLICIT NONE

  REAL(8), INTENT(IN)    :: Nn_He_m3, Nn_O_m3, Nn_N2_m3, Nn_O2_m3, Te_K
  REAL(8), INTENT(INOUT) :: energy_eV, vpar_ms, vperp_ms
  REAL(8), INTENT(OUT)   :: coll_freq_s1
  LOGICAL, INTENT(OUT)   :: ionization_occurred
  REAL(8), INTENT(OUT)   :: prod_energy_eV, prod_vpar_ms, prod_vperp_ms
  REAL(8), INTENT(IN)    :: delta_t_s, weight_m3s

  REAL(8), INTENT(IN) :: a_i, b_i

  REAL(8), INTENT(INOUT) :: point_i_kinetic_chem_rate_N2_ioniz_m3s
  REAL(8), INTENT(INOUT) :: point_i_kinetic_chem_rate_N_ioniz_m3s
  REAL(8), INTENT(INOUT) :: point_i_kinetic_chem_rate_O2_ioniz_m3s
  REAL(8), INTENT(INOUT) :: point_i_kinetic_chem_rate_O_ioniz_m3s
  REAL(8), INTENT(INOUT) :: point_i_kinetic_chem_rate_He_ioniz_m3s

  REAL(8), INTENT(INOUT) :: point_ip1_kinetic_chem_rate_N2_ioniz_m3s
  REAL(8), INTENT(INOUT) :: point_ip1_kinetic_chem_rate_N_ioniz_m3s
  REAL(8), INTENT(INOUT) :: point_ip1_kinetic_chem_rate_O2_ioniz_m3s
  REAL(8), INTENT(INOUT) :: point_ip1_kinetic_chem_rate_O_ioniz_m3s
  REAL(8), INTENT(INOUT) :: point_ip1_kinetic_chem_rate_He_ioniz_m3s

  LOGICAL, INTENT(OUT) :: particle_collided_with_neutral

  REAL(8), PARAMETER :: threshold_Nn_m3 = 1.0d8

  REAL(8) therm_energy_eV

  REAL General_Collision_Probability

  REAL Colprob_kind(1:N_colkind)
  REAL max_f

  INTEGER kk, indx_coll
  
  REAL R

! function
  REAL(8) fl_Max_frequency_of_en_collisions_s1
  REAL(8) fl_Frequency_of_Collision_s1

  particle_collided_with_neutral = .FALSE.

  ionization_occurred = .FALSE.
  coll_freq_s1 = 0.0_8

  IF ( (Nn_He_m3 + Nn_O_m3 + Nn_N2_m3 + Nn_O2_m3).LT.threshold_Nn_m3 ) RETURN

  IF (energy_eV.EQ.0.0_8) RETURN

  indx_coll = 0

  therm_energy_eV = Te_K * 2.0_8 * k_JK / e_Cl
!  therm_energy_eV = (pfl(n)%point(i)%Te_K * a_i + pfl(n)%point(i+1)%Te_K * b_i) * 2.0_8 * k_JK / e_Cl
!  therm_energy_eV = grid_therm_energy_eV(i, j, k, a_i, a_j, a_k)

!print *, therm_energy_eV, energy_eV

! do not scatter particles which are about to be considered as thermalized ones
  IF (energy_eV.LT.therm_energy_eV) RETURN

! check whether a collision occurs
  coll_freq_s1 = fl_Max_frequency_of_en_collisions_s1(energy_eV, Nn_He_m3, Nn_O_m3, Nn_N2_m3, Nn_O2_m3)
  
! fool proof
  IF (coll_freq_s1.EQ.0.0_8) RETURN

  General_Collision_Probability = REAL(1.0_8 - EXP(-coll_freq_s1 * delta_t_s))

  R = well_random_number()

  IF (R.GT.General_Collision_Probability) THEN
! no collision occurred
     RETURN
  END IF

! define which collision occurred
  Colprob_kind(1) = REAL(fl_Frequency_of_Collision_s1(1, energy_eV, Nn_He_m3, Nn_O_m3, Nn_N2_m3, Nn_O2_m3))
  DO kk = 2, N_colkind
     Colprob_kind(kk) = Colprob_kind(kk-1) + REAL(fl_Frequency_of_Collision_s1(kk, energy_eV, Nn_He_m3, Nn_O_m3, Nn_N2_m3, Nn_O2_m3))
  END DO
  max_f = Colprob_kind(N_colkind)
  IF (max_f.GT.0.0) THEN
     Colprob_kind = Colprob_kind / max_f
     Colprob_kind(N_colkind) = 1.0
  else
! weird, the fool proof above should have prevented this
     print '("Collide_Electron :: Error :: trying to process collision while all collision frequencies are zeros:")'
     print '("Coll-Error in proc ",i4," check :: ",6(e14.7))', &
          & Rank_of_process, therm_energy_eV, energy_eV, coll_freq_s1, General_Collision_Probability, R, max_f
     stop
  END IF

  R = well_random_number()
  DO kk = N_colkind, 1, -1 
     IF (R.GT.Colprob_kind(kk)) EXIT 
  END DO
  indx_coll = kk + 1
! fool proof
  indx_coll = MAX(1, MIN(N_colkind, indx_coll))

  particle_collided_with_neutral = .TRUE.

  SELECT CASE (indx_coll)
  CASE (1)       ! e+N2, elastic

     CALL CollideElectron_Elastic(energy_eV, vpar_ms, vperp_ms, 28.0_8)
  
  CASE (24)      ! e+O2, elastic

     CALL CollideElectron_Elastic(energy_eV, vpar_ms, vperp_ms, 32.0_8)

  CASE (37)      ! e+O, elastic

     CALL CollideElectron_Elastic(energy_eV, vpar_ms, vperp_ms, 16.0_8)

  CASE (53)      ! e+He, elastic

     CALL CollideElectron_Elastic(energy_eV, vpar_ms, vperp_ms, 4.0_8)

  CASE (2)     ! e+N2 ionization of N2+

! for N2, shape parameter is 12.7, see Table 1 of [Opal, Peterson, and Beaty, J.Chem.Phys, 55, 4100, 1971]

     IF (energy_eV.LE.threshold_eV(indx_coll)) RETURN
     CALL CollideElectron_Ionization(energy_eV, vpar_ms, vperp_ms, prod_energy_eV, prod_vpar_ms, prod_vperp_ms, threshold_eV(indx_coll), 12.7_8)
     ionization_occurred = .TRUE.

     point_i_kinetic_chem_rate_N2_ioniz_m3s   = point_i_kinetic_chem_rate_N2_ioniz_m3s + a_i * weight_m3s
     point_ip1_kinetic_chem_rate_N2_ioniz_m3s = point_ip1_kinetic_chem_rate_N2_ioniz_m3s + b_i * weight_m3s

  CASE (3)     ! e+N2 ionization of N+

! for N2, shape parameter is 12.7, see Table 1 of [Opal, Peterson, and Beaty, J.Chem.Phys, 55, 4100, 1971]

     IF (energy_eV.LE.threshold_eV(indx_coll)) RETURN
     CALL CollideElectron_Ionization(energy_eV, vpar_ms, vperp_ms, prod_energy_eV, prod_vpar_ms, prod_vperp_ms, threshold_eV(indx_coll), 12.7_8)
     ionization_occurred = .TRUE.

     point_i_kinetic_chem_rate_N_ioniz_m3s   = point_i_kinetic_chem_rate_N_ioniz_m3s + a_i * weight_m3s
     point_ip1_kinetic_chem_rate_N_ioniz_m3s = point_ip1_kinetic_chem_rate_N_ioniz_m3s + b_i * weight_m3s

  CASE (25)     ! e+O2 ionization of O2+

! for O2, shape parameter is 17.5, see Table 1 of [Opal, Peterson, and Beaty, J.Chem.Phys, 55, 4100, 1971]

     IF (energy_eV.LE.threshold_eV(indx_coll)) RETURN
     CALL CollideElectron_Ionization(energy_eV, vpar_ms, vperp_ms, prod_energy_eV, prod_vpar_ms, prod_vperp_ms, threshold_eV(indx_coll), 17.5_8)
     ionization_occurred = .TRUE.

     point_i_kinetic_chem_rate_O2_ioniz_m3s   = point_i_kinetic_chem_rate_O2_ioniz_m3s + a_i * weight_m3s
     point_ip1_kinetic_chem_rate_O2_ioniz_m3s = point_ip1_kinetic_chem_rate_O2_ioniz_m3s + b_i * weight_m3s

  CASE (26)     ! e+O2 ionization of O+

! for O2, shape parameter is 17.5, see Table 1 of [Opal, Peterson, and Beaty, J.Chem.Phys, 55, 4100, 1971]

     IF (energy_eV.LE.threshold_eV(indx_coll)) RETURN
     CALL CollideElectron_Ionization(energy_eV, vpar_ms, vperp_ms, prod_energy_eV, prod_vpar_ms, prod_vperp_ms, threshold_eV(indx_coll), 17.5_8)
     ionization_occurred = .TRUE.

     point_i_kinetic_chem_rate_O_ioniz_m3s   = point_i_kinetic_chem_rate_O_ioniz_m3s + a_i * weight_m3s
     point_ip1_kinetic_chem_rate_O_ioniz_m3s = point_ip1_kinetic_chem_rate_O_ioniz_m3s + b_i * weight_m3s

  CASE (38)     ! e+O ionization

! for O the shape parameter is presently unknown, use the value for Argon ?????

     IF (energy_eV.LE.threshold_eV(indx_coll)) RETURN
     CALL CollideElectron_Ionization(energy_eV, vpar_ms, vperp_ms, prod_energy_eV, prod_vpar_ms, prod_vperp_ms, threshold_eV(indx_coll), 10.0_8)
     ionization_occurred = .TRUE.

     point_i_kinetic_chem_rate_O_ioniz_m3s   = point_i_kinetic_chem_rate_O_ioniz_m3s + a_i * weight_m3s
     point_ip1_kinetic_chem_rate_O_ioniz_m3s = point_ip1_kinetic_chem_rate_O_ioniz_m3s + b_i * weight_m3s

  CASE (54)     ! e+He ionization

! for He, shape parameter is 15.4, see Table 1 of [Opal, Peterson, and Beaty, J.Chem.Phys, 55, 4100, 1971]

     IF (energy_eV.LE.threshold_eV(indx_coll)) RETURN
     CALL CollideElectron_Ionization(energy_eV, vpar_ms, vperp_ms, prod_energy_eV, prod_vpar_ms, prod_vperp_ms, threshold_eV(indx_coll), 15.4_8)
     ionization_occurred = .TRUE.

     point_i_kinetic_chem_rate_He_ioniz_m3s   = point_i_kinetic_chem_rate_He_ioniz_m3s + a_i * weight_m3s
     point_ip1_kinetic_chem_rate_He_ioniz_m3s = point_ip1_kinetic_chem_rate_He_ioniz_m3s + b_i * weight_m3s

  CASE (4:23, 27:36, 39:52, 55:60)    ! e+N2, e+02, e+O, e+He inelastic

     IF (energy_eV.LE.threshold_eV(indx_coll)) RETURN
     CALL CollideElectron_Inelastic(energy_eV, vpar_ms, vperp_ms, threshold_eV(indx_coll))

  END SELECT

END SUBROUTINE FL_Collide_Electron

!---------------------------------------------------
! the procedure is built on the basis of a similar procedure from EDIPIC
! the parallel direction is the x-direction
! the transverse direction is in the y-z plane and only the combined perpendicular velocity is calculated
!
subroutine CollideElectron_Elastic(energy_eV, vpar_ms, vperp_ms, M_neutral_amu)

  use Photoelectrons, only : factor_energy_eV
  use rng_wrapper

  implicit none

  real(8), intent(inout) :: energy_eV, vpar_ms, vperp_ms
  real(8), intent(in)    :: M_neutral_amu

  real(8) R

  real(8) Ksi              ! scattering angle (relative to the initial direction)
  real(8) CosKsi, SinKsi     
  real(8) Phi              ! azimuthal scattering angle 
  real(8) SinPhi  ! CosPhi, 

  real(8) alpha

  real(8) Vx_ms            ! temporary value to store parallel velocity

! Calculate the scattering angle relative to the initial direction of electron
  R = well_random_number()

! #####  CosKsi = (2.0_8 + energy_eV - 2.0_8 * (1.0_8 + energy_eV)**R) / energy_eV #####
! the formula above was in the older code and it was based on Surendra's differential cross section
! below is the corrected expression from Okhrimovsky et al., Phys.Rev.E, 65, 037402 (2002).
  CosKsi = 1.0_8 - 2.0_8 * R / (1.0_8 + 8.0_8 * (energy_eV / 27.21_8) * (1.0_8 - R))

  CosKsi = MAX(MIN(0.999999999999_8, CosKsi), -0.999999999999_8)   !############ to avoid an unlikely situation when |CosKsi|>1
!###  Ksi = ACOS(CosKsi)
!###  SinKsi = SIN(Ksi)
  SinKsi = SQRT(1.0_8 - CosKsi**2)

! Calculate the azimuthal scattering angle
  R = well_random_number()
  Phi = R * 6.28318530718_8
!  CosPhi = COS(Phi)
  SinPhi = SIN(Phi)

! Calculate shrinking factor for the velocity 
! note that (0.001097161_8 / M_neutral_amu) * (1.0_8 - CosKsi) is the relative energy drop (Eq.12 of Vahedi)
! 0.001097161 = 2 * m_e_kg / 1_amu_kg = 2 * 9.109534e-31 / 1.660565e-27 

  alpha = 1.0_8 - (0.001097161_8 / M_neutral_amu) * (1.0_8 - CosKsi)
!  energy_eV = energy_eV * alpha

  alpha = sqrt(alpha)

! Transform the velocity 
  Vx_ms    = vpar_ms  * CosKsi + SinPhi * vperp_ms * SinKsi
  vperp_ms = sqrt(max(0.0_8, vpar_ms**2 + vperp_ms**2 - Vx_ms**2)) * alpha
  vpar_ms  = Vx_ms * alpha

  energy_eV = factor_energy_eV * (vperp_ms**2 + vpar_ms**2) ! factor_energy_eV = 0.5_8 * m_e_kg  / e_Cl

!  Vy_2_ms = (vperp_ms * CosKsi - SinPhi * vpar_ms  * SinKsi)**2
!  Vz_2_ms = (vpar_ms**2 + vperp_ms**2) * (CosPhi * SinKsi)**2

end subroutine CollideElectron_Elastic

!------------------------------------------------
! this procedure is similar to the CollideElectron_Elastic above
! except for the final energy of the particle (drop by the collision threshold energy is introduced)
! and the energy used to define scattering angle (the energy after the collision is used). 
!
subroutine CollideElectron_Inelastic(energy_eV, vpar_ms, vperp_ms, threshold_eV)

  use Photoelectrons, only : factor_energy_eV
  use rng_wrapper

  implicit none

  real(8), intent(inout) :: energy_eV, vpar_ms, vperp_ms
  real(8), intent(in)    :: threshold_eV

  real(8) energy_sc_eV     ! energy of scattered electron after the inelastic collision

  real(8) R

  real(8) Ksi              ! scattering angle (relative to the initial direction)
  real(8) CosKsi, SinKsi     
  real(8) Phi              ! azimuthal scattering angle 
  real(8) SinPhi ! CosPhi,   

  real(8) alpha

  real(8) Vx_ms            ! temporary value to store parallel velocity

! Calculate the energy of the scattered electron 
  energy_sc_eV = max(0.0_8, energy_eV - threshold_eV)

! Calculate the scattering angle relative to the initial direction of electron
  R = well_random_number()

! #####  CosKsi = (2.0_8 + energy_eV - 2.0_8 * (1.0_8 + energy_eV)**R) / energy_eV #####
! the formula above was in the older code and it was based on Surendra's differential cross section
! below is the corrected expression from Okhrimovsky et al., Phys.Rev.E, 65, 037402 (2002).
  CosKsi = 1.0_8 - 2.0_8 * R / (1.0_8 + 8.0_8 * (energy_sc_eV / 27.21_8) * (1.0_8 - R))

  CosKsi = MAX(MIN(0.999999999999_8, CosKsi), -0.999999999999_8)   !############ to avoid an unlikely situation when |CosKsi|>1
!###  Ksi = ACOS(CosKsi)
!###  SinKsi = SIN(Ksi)
  SinKsi = SQRT(1.0_8 - CosKsi**2)

! Calculate the azimuthal scattering angle
  R = well_random_number()
  Phi = R * 6.28318530718_8
!  CosPhi = COS(Phi)
  SinPhi = SIN(Phi)

! Calculate shrinking factor for the velocity 
  alpha = sqrt(energy_sc_eV / energy_eV)
!  energy_eV = energy_sc_eV

! Transform the velocity 
  Vx_ms    = vpar_ms  * CosKsi + SinPhi * vperp_ms * SinKsi
  vperp_ms = sqrt(max(0.0_8, vpar_ms**2 + vperp_ms**2 - Vx_ms**2)) * alpha
  vpar_ms  = Vx_ms * alpha

  energy_eV = factor_energy_eV * (vperp_ms**2 + vpar_ms**2) ! factor_energy_eV = 0.5_8 * m_e_kg  / e_Cl

end subroutine CollideElectron_Inelastic

!-------------------------------------------------------
!
subroutine CollideElectron_Ionization(energy_eV, vpar_ms, vperp_ms, energy_ej_eV, vpar_ej_ms, vperp_ej_ms, threshold_eV, B_of_Einc_eV)

  use Photoelectrons, only : factor_energy_eV
  use rng_wrapper

  implicit none

  real(8), intent(inout) :: energy_eV, vpar_ms, vperp_ms             ! incident/scattered (primary)
  real(8), intent(out)   :: energy_ej_eV, vpar_ej_ms, vperp_ej_ms    ! ejected (secondary)
  real(8), intent(in)    :: threshold_eV

  real(8), intent(in)    :: B_of_Einc_eV

  real(8) R

  real(8) energy_sc_eV

! scattered (primary) electron
  real(8) Ksi_sc, CosKsi_sc, SinKsi_sc     ! pitch scattering angle (relative to the initial direction)
  real(8) Phi_sc, SinPhi_sc                ! azimuthal scattering angle 

! ejected (secondary) electron
  real(8) Ksi_ej, CosKsi_ej, SinKsi_ej     ! pitch scattering angle (relative to the initial direction)
  real(8) Phi_ej, SinPhi_ej                ! azimuthal scattering angle 

  real(8) alpha

  real(8) Vx_ms            ! temporary value to store parallel velocity

!### input parameter now
!###  B_of_Einc_eV = 10.0_8                   ! is constant for Einc < 70 eV, Argon, approximate             
                                          ! the smaller this value the wider the energy spectrum of ejected electrons
                                          ! at high primary electron energies

! Calculate energy of the ejected electron
  R = well_random_number()
  energy_ej_eV = B_of_Einc_eV * TAN(R * ATAN( 0.5_8 * (energy_eV - threshold_eV) / B_of_Einc_eV ))  

! Calculate energy of the scattered (incident or primary) electron
  energy_sc_eV = energy_eV - threshold_eV - energy_ej_eV

! fool-proof check
  if ((energy_ej_eV.lt.0.0_8).or.(energy_sc_eV.lt.0.0_8)) then
     print '("Error in CollideElectron_Ionization :: negative energy of scattered and/or ejected electron :")'
     print '("scattered ",e12.5," eV, ejected ",e12.5," eV")', energy_sc_eV, energy_ej_eV
     print '("incident w/vpar/vperp : ",3(2x,e12.5))', energy_eV, vpar_ms, vperp_ms
     stop
  end if

! Calculate scattering angles for the ejected electron (secondary)
  R = well_random_number()
  CosKsi_ej = 1.0_8 - 2.0_8 * R / (1.0_8 + 8.0_8 * (energy_ej_eV / 27.21_8) * (1.0_8 - R))
  CosKsi_ej = MAX(MIN(0.999999999999_8, CosKsi_ej), -0.999999999999_8)   !############ to avoid an unlikely situation when |CosKsi|>1
!###  Ksi_ej    = ACOS(CosKsi_ej)
!###  SinKsi_ej = SIN(Ksi_ej)
  SinKsi_ej = SQRT(1.0_8 - CosKsi_ej**2)
  R         = well_random_number()
  Phi_ej    = R * 6.28318530718_8
  SinPhi_ej = SIN(Phi_ej)

! Calculate shrinking factor for the velocity 
  alpha = sqrt(energy_ej_eV / energy_eV)

! Transform the velocity 
  Vx_ms       = vpar_ms  * CosKsi_ej + SinPhi_ej * vperp_ms * SinKsi_ej
  vperp_ej_ms = sqrt(max(0.0_8, vpar_ms**2 + vperp_ms**2 - Vx_ms**2)) * alpha
  vpar_ej_ms  = Vx_ms * alpha

  energy_ej_eV = factor_energy_eV * (vperp_ej_ms**2 + vpar_ej_ms**2) ! factor_energy_eV = 0.5_8 * m_e_kg  / e_Cl

! Calculate scattering angles for the scattered electron (primary)
  R = well_random_number()
  CosKsi_sc = 1.0_8 - 2.0_8 * R / (1.0_8 + 8.0_8 * (energy_sc_eV / 27.21_8) * (1.0_8 - R))
  CosKsi_sc = MAX(MIN(0.999999999999_8, CosKsi_sc), -0.999999999999_8)   !############ to avoid an unlikely situation when |CosKsi|>1
!###  Ksi_sc    = ACOS(CosKsi_sc)
!###  SinKsi_sc = SIN(Ksi_sc)
  SinKsi_sc = SQRT(1.0_8 - CosKsi_sc**2)
  R         = well_random_number()
  Phi_sc    = R * 6.28318530718_8
  SinPhi_sc = SIN(Phi_sc)

! Calculate shrinking factor for the velocity 
  alpha = sqrt(energy_sc_eV / energy_eV)

! Transform the velocity 
  Vx_ms    = vpar_ms  * CosKsi_sc + SinPhi_sc * vperp_ms * SinKsi_sc
  vperp_ms = sqrt(max(0.0_8, vpar_ms**2 + vperp_ms**2 - Vx_ms**2)) * alpha
  vpar_ms  = Vx_ms * alpha

  energy_eV = factor_energy_eV * (vperp_ms**2 + vpar_ms**2) ! factor_energy_eV = 0.5_8 * m_e_kg  / e_Cl

end subroutine CollideElectron_Ionization

!----------------------------------------------
!
subroutine get_local_pitch_angle_range( loc_min_angle, &
                                      & loc_max_angle, &
                                      & loc_B, &
                                      & given_min_angle, &
                                      & given_max_angle, &
                                      & given_B )
  
  use PhysicalConstants, only : pi

  implicit none

  real(8), intent(out) :: loc_min_angle
  real(8), intent(out) :: loc_max_angle

  real(8), intent(in) :: loc_B
  real(8), intent(in) :: given_min_angle
  real(8), intent(in) :: given_max_angle
  real(8), intent(in) :: given_B

  real(8) temp

  temp = 1.0_8 - loc_B * (sin(given_max_angle))**2 / given_B

  if (given_max_angle.gt.(0.5_8*pi)) then
     loc_max_angle = acos(-sqrt(max(temp, 0.0_8)))
  else
     loc_max_angle = acos(sqrt(max(temp, 0.0_8)))
  end if

  temp = 1.0_8 - loc_B * (sin(given_min_angle))**2 / given_B

  if (given_max_angle.gt.(0.5_8*pi)) then
     loc_min_angle = acos(-sqrt(max(temp, 0.0_8)))
  else
     loc_min_angle = acos(sqrt(max(temp, 0.0_8)))
  end if

  loc_max_angle = max(loc_max_angle, loc_min_angle)

  return

end subroutine get_local_pitch_angle_range
