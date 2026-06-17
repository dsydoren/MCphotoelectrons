
!------------------------------------
!
MODULE ParallelOperationValues

  INTEGER Rank_of_process          ! process rank in the global communicator
  INTEGER N_of_processes           ! number of processes the global communicator

END MODULE  ParallelOperationValues

!------------------------------------
!
MODULE PhysicalConstants  !CurrentProblemValues

! general constants
  REAL(8), PARAMETER :: pi       = 3.14159265358979_8                    ! pi
  REAL(8), PARAMETER :: e_Cl     = 1.602189d-19                          ! Elementary charge [Cl]
  REAL(8), PARAMETER :: m_e_kg   = 9.109534d-31                          ! Electron mass [kg]
  REAL(8), PARAMETER :: eps_0_Fm = 8.854188d-12                          ! The dielectric constant [F/m]
  REAL(8), PARAMETER :: mu_0_Hm  = 1.256637d-6                           ! The magnetic constant [H/m]
  REAL(8), PARAMETER :: k_JK     = 1.3806503d-23                         ! Boltzmann constant [J/K]
!  REAL(8), PARAMETER :: c_ms =299792462.2 != 1.0  / SQRT(eps_0_Fm * mu_0_Hm)          ! The light speed in vacuum [m/s]
!  REAL(8) c_ms !=299792462.2 != 1.0  / SQRT(eps_0_Fm * mu_0_Hm)          ! The light speed in vacuum [m/s]

  REAL(8), PARAMETER :: rad_to_deg = 180.0_8 / pi
  REAL(8), PARAMETER :: deg_to_rad = pi / 180.0_8

  REAL(8), parameter :: h_Planck_Js = 6.62607015d-34
  REAL(8), parameter :: c_ms = 299792458.0_8

! constant Earth parameters
  REAL(8), PARAMETER ::  R_Earth_km  = 6371.0_8 !### changed to 6371.0_8 for consistency with iri ###                  ! Earth surface radius [km]
  REAL(8), PARAMETER ::  g_Earth_ms2 = 9.80665_8                  ! gravity acceleration at the Earth surface [m/s^2]
  REAL(8), PARAMETER ::  M_Earth_Am2 = 8.0d22                     ! Earth magnetic moment [A m^2]

! relative atomic masses of ion species
  REAL(8), PARAMETER :: M_H_amu = 1.0_8                           ! Hydrogen ion mass [a.m.u.]
!  REAL(8), PARAMETER :: M_He_amu = 4.0_8                          ! Helium ion mass [a.m.u.]
  REAL(8), PARAMETER :: M_N_amu = 14.0_8                          ! Nitrogen ion mass [a.m.u.]
  REAL(8), PARAMETER :: M_O_amu = 16.0_8                          ! Oxygen ion mass [a.m.u.]
  REAL(8), PARAMETER :: M_N2_amu = 28.0_8                         ! Nitrogen molecula mass [a.m.u.]
  REAL(8), PARAMETER :: M_NO_amu = 30.0_8                         ! Nitric oxide molecula mass [a.m.u.] 
  REAL(8), PARAMETER :: M_O2_amu = 32.0_8                         ! Oxygen molecula mass [a.m.u.]

  REAL(8), PARAMETER :: amu_kg   = 1.660565d-27                  ! 1 a.m.u. [kg]

END MODULE PhysicalConstants

!---------------------------
!
MODULE TimeValues

  INTEGER year
  INTEGER month
  INTEGER day_of_month
  INTEGER day_of_year
  REAL time_ut_h       ! universal time, hours decimal
  REAL MJD_days        ! Modified Julian Time [days]
  REAL T0_jcent        ! Time from 12:00 UT on Jan-1, 2000 [Julian centuries]

END MODULE TimeValues

!---------------------------
!
MODULE GlobalIndices

  REAL AP(1:7)                  ! array with Ap indices for MSIS
  REAL f10p7_pd, f10p7_81       ! previous day (pd) and 81 day averaged (81) solar radio fluxes for MSIS
  REAL f10p7                    ! same day solar radio flux for EUVAC

END MODULE GlobalIndices

!---------------------------------
!
MODULE CoordinateTransformations

! matrix of transformation from GSE to GEO
  REAL(8) gse_2_geo_11, gse_2_geo_12, gse_2_geo_13
  REAL(8) gse_2_geo_21, gse_2_geo_22, gse_2_geo_23
  REAL(8) gse_2_geo_31, gse_2_geo_32, gse_2_geo_33

END MODULE CoordinateTransformations

!------------------------------
!
MODULE SpacecraftValues

  INTEGER N_of_orbit_points
  INTEGER orbit_number
  INTEGER first_op    ! >=1                 ! only orbit points with ordering number within this range will be processed
  INTEGER last_op     ! <=N_of_orbit_points !

  TYPE spacecraft_data
     INTEGER year
     INTEGER month
     INTEGER day_of_month
     REAL ut_h

     REAL(8) r_km
     REAL(8) theta_geo
     REAL(8) phi_geo

     REAL(8) x_geo_km
     REAL(8) y_geo_km
     REAL(8) z_geo_km

     REAL total_euv_rate_m3s
  END TYPE spacecraft_data

  TYPE(spacecraft_data), ALLOCATABLE :: orbit_point(:)

END MODULE SpacecraftValues

!-----------------------------------------------
!
MODULE Photoelectrons

  INTEGER, PARAMETER :: EUVAC = 0
  INTEGER, PARAMETER :: f76 = 1
  INTEGER, PARAMETER :: f74113 = 2

  INTEGER EUV_spectrum_model
  INTEGER N_split_EUVAC

  INTEGER max_N_cascades_per_EUV_flux
  REAL(8) EUV_energy_flux_onecascade_1e9eVcm2s   ! in units of (10^9 eV/cm^2/s) reasonable values are 0.5-0.2

  REAL(8) collfreq_precision_factor  ! used here :: Coll_dt_s = collfreq_precision_factor / recent_max_en_coll_frequency_s1, do not recommend setting it >0.2

  REAL(8) factor_energy_eV
  REAL(8) efactor_eV_to_ms

  INTEGER N_solar_bins  ! 916  for f76, 853 for f74113

  REAL, ALLOCATABLE :: solar_flux_phcm2s1(:)                 ! solar flux [10^9 photons cm^-2 s^-1]
  REAL, ALLOCATABLE :: solar_flux_energy_bin_eV(:)

! used to calculate photon flux attenuation
  REAL, ALLOCATABLE :: sigma_tot_O_cm2(:)            ! total cross section [10^-18 cm^2] for atomic oxygen          
  REAL, ALLOCATABLE :: sigma_tot_N2_cm2(:)           ! - " -, molecular nitrogen             
  REAL, ALLOCATABLE :: sigma_tot_O2_cm2(:)           ! - " -, molecular oxygen

! used to calculate ionization rates
  REAL, ALLOCATABLE :: sigma_ion_N2_to_N2_cm2(:)     ! ionization cross section [10^-18 cm^2] for molecular nitrogen, N2+
  REAL, ALLOCATABLE :: sigma_ion_N2_to_N_cm2(:)      ! - " -, molecular nitrogen, N+

  REAL, ALLOCATABLE :: sigma_ion_O_4S_cm2(:)         ! - " -, atomic oxygen, O+4S
  REAL, ALLOCATABLE :: sigma_ion_O_2D_cm2(:)         ! - " -, atomic oxygen, O+2D
  REAL, ALLOCATABLE :: sigma_ion_O_2P_cm2(:)         ! - " -, atomic oxygen, O+2P
  REAL, ALLOCATABLE :: sigma_ion_O_4Pst_cm2(:)       ! - " -, atomic oxygen, O+4P*
  REAL, ALLOCATABLE :: sigma_ion_O_2Pst_cm2(:)       ! - " -, atomic oxygen, O+2P*

  REAL, ALLOCATABLE :: sigma_ion_O2_to_O2_cm2(:)     ! - " -, molecular oxygen, O2+
  REAL, ALLOCATABLE :: sigma_ion_O2_to_O_cm2(:)      ! - " -, molecular oxygen, O+

  INTEGER, PARAMETER :: N_photoe_producing_channels = 9          ! number of ionization reactions between an EUV photon and a neutral

  integer, parameter :: N_colkind = 59
  real(8) threshold_eV(1:N_colkind)

  real(8) threshold_photoion_N2_to_N2_eV
  real(8) threshold_photoion_N2_to_N_eV

  real(8) threshold_photoion_O_4S_eV
  real(8) threshold_photoion_O_2D_eV
  real(8) threshold_photoion_O_2P_eV
  real(8) threshold_photoion_O_4Pst_eV
  real(8) threshold_photoion_O_2Pst_eV

  real(8) threshold_photoion_O2_to_O2_eV
  real(8) threshold_photoion_O2_to_O_eV

  INTEGER max_evdf_N_vpar
  INTEGER max_evdf_N_vperp
  REAL(8) dvpar_ms
  REAL(8) dvperp_ms

! in spacecraft location
  REAL(8), ALLOCATABLE :: evdf(:,:)   ! photoelectron distribution function over parallel and transverse velocities

END MODULE Photoelectrons

!---------------------------------------
!
MODULE field_line

  LOGICAL pfl_closed
  INTEGER pfl_N_of_points
  INTEGER pfl_i_spacecraft

  REAL(8) pfl_desired_dL_km

  TYPE field_line_coordinates
     REAL(8) r   ! in units of R_Earth_km
     REAL(8) theta_geo
     REAL(8) phi_geo
     REAL(8) x_geo_km
     REAL(8) y_geo_km
     REAL(8) z_geo_km
  END TYPE field_line_coordinates

  TYPE field_line_data
     REAL(8) L_m
     REAL(8) fts
     REAL(8) columnar_content_O_m2
     REAL(8) columnar_content_N2_m2
     REAL(8) columnar_content_O2_m2
     REAL(8) Nn_O_m3
     REAL(8) Nn_N2_m3
     REAL(8) Nn_O2_m3
     REAL(8) geoB_T
     REAL(8) Te_K
     REAL(8) Ne_m3

     REAL(8) rate_iN_production_m3s
     REAL(8) rate_iO_production_m3s
     REAL(8) rate_iN2_production_m3s
     REAL(8) rate_iO2_production_m3s
     REAL(8) kinetic_ne_m3
     REAL(8) Ge_plus_m2s
     REAL(8) Ge_minus_m2s

     REAL(8), ALLOCATABLE :: Ge_m2s_range(:)
     REAL(8), ALLOCATABLE :: Ge_par_m2s_range(:)

     real(8), allocatable :: source_EUV_m3s_range_channel(:,:)
     real(8), allocatable :: source_photoe_m3s_range_channel(:,:)

     real(8), allocatable :: source_neutrals_m3s_range_coll(:,:)
     real(8), allocatable :: source_plasma_m3s_range(:)
     real(8), allocatable :: source_Coulomb_m3s_range(:)

     real(8), allocatable :: sink_neutrals_m3s_range_coll(:,:)
     real(8), allocatable :: sink_plasma_m3s_range(:)
     real(8), allocatable :: sink_Coulomb_m3s_range(:)

  END TYPE field_line_data

  TYPE(field_line_coordinates), ALLOCATABLE :: pfl_point_coords(:)
  TYPE(field_line_data),        ALLOCATABLE :: shared_pfl_point(:)

  INTEGER N_ranges
  REAL(8), ALLOCATABLE :: min_energy_eV_range(:)
  REAL(8), ALLOCATABLE :: max_energy_eV_range(:)
  REAL(8), ALLOCATABLE :: min_pitch_angle_range(:)
  REAL(8), ALLOCATABLE :: max_pitch_angle_range(:)

  integer N_coll_ranges
  real(8), allocatable :: min_w_eV_colrange(:)
  real(8), allocatable :: max_w_eV_colrange(:)
  real(8), allocatable :: coll_freq_s1_L_colkind_colrange(:,:,:)
  real(8), allocatable :: density_m3_L_colrange(:,:)

END MODULE field_line

!--------------------------------------
!
module combined_collisions

  real(8), parameter :: lowen_dw_eV = 0.05_8
  real(8), parameter :: max_lowen_w_eV = 6.0_8
  real(8), parameter :: max_highen_w_eV = 450.0_8  ! increase if the EUV spectrum includes energies higher than f76 (435.487eV at 28.47A)

  integer N_lowen  ! array size 0:N_lowen

  real(8), allocatable :: lowen_CSV_N2_total_m3s(:)
  real(8), allocatable :: lowen_CSV_O2_total_m3s(:)
  real(8), allocatable :: lowen_CSV_O_total_m3s(:)

  integer min_N_highen  ! 
  integer N_highen      ! array size min_N_highen:N_highen

  real(8), allocatable :: highen_CSV_N2_total_m3s(:)
  real(8), allocatable :: highen_CSV_O2_total_m3s(:)
  real(8), allocatable :: highen_CSV_O_total_m3s(:)

end module combined_collisions
