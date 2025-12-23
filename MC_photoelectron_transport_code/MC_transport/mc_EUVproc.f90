!---------------------------------
!
SUBROUTINE PREPARE_PHOTOELECTRONS

  USE ParallelOperationValues
  USE PhysicalConstants
  USE photoelectrons
  USE field_line, ONLY : pfl_desired_dL_km

  IMPLICIT NONE

  LOGICAL exists
  CHARACTER(1) buf
  INTEGER ALLOC_ERR

  REAL(8) max_energy_par_eV, max_energy_perp_eV

  INQUIRE (FILE = 'input_photoelectrons.dat', EXIST = exists)
!  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  IF (.NOT.exists) THEN
! file does not exist
     IF (Rank_of_process.eq.0) PRINT '("file input_photoelectrons.dat not found, terminate program")'
     STOP
  END IF

  IF (Rank_of_process.eq.0) PRINT '(2x,"File input_photoelectrons.dat is found. Reading the data file...")'
         
  OPEN (9, FILE = 'input_photoelectrons.dat')

  READ (9, '(A1)') buf ! number of primary particles per solar flux energy bin (>0)
  READ (9, *) N_cascades_per_energy_bin

  READ (9, '(A1)') buf ! desired distance between points on a field line [km]
  READ (9, *) pfl_desired_dL_km

  READ (9, '(A1)') buf ! number of energy bins in electron energy distribution function (EEDF)
  READ (9, *) eedf_N_of_energy_bins

  READ (9, '(A1)') buf ! energy bin size of the EEDF [eV]
  READ (9, *) energy_bin_size_eV

  READ (9, '(A1)') buf ! number of parallel velocity bins in one direction in electron velocity distribution function
  READ (9, *) max_evdf_N_vpar

  READ (9, '(A1)') buf ! number of transverse velocity bins in electron velocity distribution function
  READ (9, *) max_evdf_N_vperp

  READ (9, '(A1)') buf ! maximal energy in the parallel direction [eV]
  READ (9, *) max_energy_par_eV

  READ (9, '(A1)') buf ! maximal energy in the transverse direction [eV]
  READ (9, *) max_energy_perp_eV

!  desired_step_ray_ds = desired_step_ray_ds /  R_Earth_km

  dvpar_ms = SQRT(2.0_8 * max_energy_par_eV * e_Cl / m_e_kg) / max_evdf_N_vpar

  dvperp_ms = SQRT(2.0_8 * max_energy_perp_eV * e_Cl / m_e_kg) / max_evdf_N_vperp

  CLOSE (9, STATUS = 'KEEP')

  factor_energy_eV = 0.5_8 * m_e_kg  / e_Cl

  efactor_eV_to_ms = SQRT(2.0_8 * e_Cl / m_e_kg)

  CALL Initiate_en_collisions       ! ### random number generator seeds ### and thresholds are prepared here

END SUBROUTINE PREPARE_PHOTOELECTRONS

!---------------------------------------------------------------
!
SUBROUTINE UPDATE_TIME(given_year, given_month, given_day_of_month, given_ut_h)

  USE TimeValues

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: given_year, given_month, given_day_of_month
  REAL, INTENT(IN) :: given_ut_h

  INTEGER days_in_the_month(1:12)  ! array with numbers of days in the month
  
  INTEGER m, y

  year = given_year
  month = given_month
  day_of_month = given_day_of_month

  time_ut_h = given_ut_h !time_ut_h_0 + REAL(delta_t_s * T_cntr)/3600.0

! set the number of days in all months except February
  days_in_the_month( 1) = 31  ! January
  days_in_the_month( 3) = 31  ! March
  days_in_the_month( 4) = 30  ! April 
  days_in_the_month( 5) = 31  ! May
  days_in_the_month( 6) = 30  ! June
  days_in_the_month( 7) = 31  ! July
  days_in_the_month( 8) = 31  ! August
  days_in_the_month( 9) = 30  ! September
  days_in_the_month(10) = 31  ! October
  days_in_the_month(11) = 30  ! November
  days_in_the_month(12) = 31  ! December

! set the number of days in the February of the given year
  IF ( (MOD(year,400).EQ.0) .OR. ((MOD(year,4).EQ.0).AND.(MOD(year,100).NE.0)) ) THEN
     days_in_the_month(2) = 29
  ELSE
     days_in_the_month(2) = 28
  END IF

! calculate the Modified Julian Time as a time in days from 0:00 UT on Nov-17, 1858

! start with contributions from days of the given month and hours
  MJD_days = day_of_month - 1 !+ time_ut_h / 24.0

! add contributions from full months of the given year
  DO m = 1, month-1
     MJD_days = MJD_days + days_in_the_month(m)
  END DO

  day_of_year = INT(MJD_days + 1.01)

! add contribution of the period between 11-17-1858 and 01-01-1859
  MJD_days = MJD_days + 14 + 31

! finally, accumulate contribution from full years starting with 01-01-1859
  DO y = 1859, year-1
     IF ( (MOD(y,400).EQ.0) .OR. ((MOD(y,4).EQ.0).AND.(MOD(y,100).NE.0)) )  THEN
! leap year
        MJD_days = MJD_days + 366
     ELSE
! ordinary year
        MJD_days = MJD_days + 365
     END IF
  END DO

!  MJD_days = REAL(MJD_days_0)
!! use the do-loop below to set proper values when restarting
!  DO WHILE (time_ut_h.GE.24.0)
!     time_ut_h = time_ut_h - 24.0
!     MJD_days = MJD_days + 1.0    !### there was a bug here:: MJD_days_0=MJD_days_0+1
!  END DO
     
  MJD_days = MJD_days + time_ut_h / 24.0

! calculate time in Julian centuries from 12:00 UT on Jan-1, 2000
  T0_jcent = (MJD_days - 51544.5) / 36525.0

END SUBROUTINE UPDATE_TIME

!-----------------------------------------------------
!
FUNCTION Get_EUV_production_rates_m3s( columnar_content_He_m2, &
                                     & columnar_content_O_m2, &
                                     & columnar_content_N2_m2, &
                                     & columnar_content_O2_m2, &
                                     & Nn_He_m3, Nn_O_m3, Nn_N2_m3, Nn_O2_m3)

  USE Photoelectrons

  IMPLICIT NONE

  REAL(8) Get_EUV_production_rates_m3s(1:10) ! 1 = He+; 2-6 = O+(4S/2D/2P/4P/2P*) from O; 7,8 = N+,N2+ from N2; 9,10 = O+,O2+ from O2
  
  REAL(8) columnar_content_He_m2
  REAL(8) columnar_content_O_m2
  REAL(8) columnar_content_N2_m2
  REAL(8) columnar_content_O2_m2
  real(8) Nn_He_m3, Nn_O_m3, Nn_N2_m3, Nn_O2_m3

  INTEGER i
  REAL(8) optical_depth
  REAL(8) temp

  Get_EUV_production_rates_m3s = 0.0_8

! special case, the observation point is in the Earth's shadow
  IF ( (columnar_content_He_m2.LT.0.0_8).OR. &
     & (columnar_content_O_m2.LT.0.0_8).OR. &
     & (columnar_content_N2_m2.LT.0.0_8).OR. &
     & (columnar_content_O2_m2.LT.0.0_8) ) RETURN

  DO i = 1, N_solar_bins

     optical_depth = 1.0d-22 * ( sigma_tot_He_cm2(i) * columnar_content_He_m2 + &
                               & sigma_tot_O_cm2(i)  * columnar_content_O_m2 + &
                               & sigma_tot_N2_cm2(i) * columnar_content_N2_m2 + &
                               & sigma_tot_O2_cm2(i) * columnar_content_O2_m2 )    ! 1e-22 because sigma is in 10^-18 cm^2 while columnar content is in m^-2 

     temp = dble(solar_flux_phcm2s1(i)) * EXP(-optical_depth)

     Get_EUV_production_rates_m3s(1) = Get_EUV_production_rates_m3s(1) + temp * sigma_ion_He_cm2(i)

     Get_EUV_production_rates_m3s(2) = Get_EUV_production_rates_m3s(2) + temp * sigma_ion_O_4S_cm2(i)
     Get_EUV_production_rates_m3s(3) = Get_EUV_production_rates_m3s(3) + temp * sigma_ion_O_2D_cm2(i)
     Get_EUV_production_rates_m3s(4) = Get_EUV_production_rates_m3s(4) + temp * sigma_ion_O_2P_cm2(i)
     Get_EUV_production_rates_m3s(5) = Get_EUV_production_rates_m3s(5) + temp * sigma_ion_O_4P_cm2(i)
     Get_EUV_production_rates_m3s(6) = Get_EUV_production_rates_m3s(6) + temp * sigma_ion_O_2Pst_cm2(i)

     Get_EUV_production_rates_m3s(7) = Get_EUV_production_rates_m3s(7) + temp * sigma_ion_N2_to_N_cm2(i)
     Get_EUV_production_rates_m3s(8) = Get_EUV_production_rates_m3s(8) + temp * sigma_ion_N2_to_N2_cm2(i)

     Get_EUV_production_rates_m3s(9)  = Get_EUV_production_rates_m3s(9)  + temp * sigma_ion_O2_to_O_cm2(i)
     Get_EUV_production_rates_m3s(10) = Get_EUV_production_rates_m3s(10) + temp * sigma_ion_O2_to_O2_cm2(i)
  END DO

  Get_EUV_production_rates_m3s(1) = Get_EUV_production_rates_m3s(1) * 1.0d-9 * Nn_He_m3   ! 1e-9 because solar flux is in 10^9 ph cm^-2 s^-1 and cross sections are in 10^-18 cm^2

  Get_EUV_production_rates_m3s(2) = Get_EUV_production_rates_m3s(2) * 1.0d-9 * Nn_O_m3   !
  Get_EUV_production_rates_m3s(3) = Get_EUV_production_rates_m3s(3) * 1.0d-9 * Nn_O_m3   !
  Get_EUV_production_rates_m3s(4) = Get_EUV_production_rates_m3s(4) * 1.0d-9 * Nn_O_m3   !
  Get_EUV_production_rates_m3s(5) = Get_EUV_production_rates_m3s(5) * 1.0d-9 * Nn_O_m3
  Get_EUV_production_rates_m3s(6) = Get_EUV_production_rates_m3s(6) * 1.0d-9 * Nn_O_m3

  Get_EUV_production_rates_m3s(7) = Get_EUV_production_rates_m3s(7) * 1.0d-9 * Nn_N2_m3
  Get_EUV_production_rates_m3s(8) = Get_EUV_production_rates_m3s(8) * 1.0d-9 * Nn_N2_m3

  Get_EUV_production_rates_m3s(9)  = Get_EUV_production_rates_m3s(9)  * 1.0d-9 * Nn_O2_m3
  Get_EUV_production_rates_m3s(10) = Get_EUV_production_rates_m3s(10) * 1.0d-9 * Nn_O2_m3

END FUNCTION Get_EUV_production_rates_m3s
