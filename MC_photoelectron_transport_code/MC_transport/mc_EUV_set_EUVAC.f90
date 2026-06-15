!------------------------------------
!
! all EUVAC parameters are from Appendix J of 
! "Ionospheres. Physics, Plasma Physics, and Chemistry", Second edition, by Robert Schunk and Andrew Nagy
! Cambridge University Press, 2009
!
! the original EUVAC is described in 
! Richards, P. G., Fennelly, J. A., & Torr, D. G. (1994).
! EUVAC: A solar EUV flux model for aeronomic calculations.
! Journal of Geophysical Research: Space Physics, 99 (A5), 8981-8992. doi: https://doi.org/10.1029/94JA00518
!
SUBROUTINE Set_EUVAC

  USE ParallelOperationValues
  USE PhysicalConstants
  USE Photoelectrons
  USE GlobalIndices, ONLY : f10p7, f10p7_81

  IMPLICIT NONE

  REAL lambda_EUVAC_A(37)
  
  REAL factor_f107

  REAL solar_flux_EUVAC_phcm2s1(37)

  REAL sigma_EUVAC_tot_O_cm2(37)
  REAL sigma_EUVAC_tot_O2_cm2(37)
  REAL sigma_EUVAC_tot_N2_cm2(37)

  REAL sigma_EUVAC_ion_N2_to_N2_cm2(37)
  REAL sigma_EUVAC_ion_N2_to_N_cm2(37)

  REAL sigma_EUVAC_ion_O_4S_cm2(37)
  REAL sigma_EUVAC_ion_O_2D_cm2(37)
  REAL sigma_EUVAC_ion_O_2P_cm2(37)
  REAL sigma_EUVAC_ion_O_4Pst_cm2(37)
  REAL sigma_EUVAC_ion_O_2Pst_cm2(37)

  REAL sigma_EUVAC_ion_O2_to_O2_cm2(37)
  REAL sigma_EUVAC_ion_O2_to_O_cm2(37)

  REAL factor_energy_eVA
  REAL d_lambda_A

  REAL, ALLOCATABLE :: lambda_EUV_A(:) !(N_solar_bins)
  INTEGER ALLOC_ERR

  INTEGER i, n, i_split, i_EUVAC

! EUVAC solar flux wavelengths (lines or bin middles)

  lambda_EUVAC_A( 1) = 0.5 * ( 50.0 + 100.0)
  lambda_EUVAC_A( 2) = 0.5 * (100.0 + 150.0) 
  lambda_EUVAC_A( 3) = 0.5 * (150.0 + 200.0) 
  lambda_EUVAC_A( 4) = 0.5 * (200.0 + 250.0) 
  lambda_EUVAC_A( 5) = 256.32 
  lambda_EUVAC_A( 6) = 284.15
  lambda_EUVAC_A( 7) = 0.5 * (250.0 + 300.0)
  lambda_EUVAC_A( 8) = 303.31 
  lambda_EUVAC_A( 9) = 303.78 
  lambda_EUVAC_A(10) = 0.5 * (300.0 + 350.0) 

  lambda_EUVAC_A(11) = 368.07
  lambda_EUVAC_A(12) = 0.5 * (350.0 + 400.0) 
  lambda_EUVAC_A(13) = 0.5 * (400.0 + 450.0) 
  lambda_EUVAC_A(14) = 465.22 
  lambda_EUVAC_A(15) = 0.5 * (450.0 + 500.0) 
  lambda_EUVAC_A(16) = 0.5 * (500.0 + 550.0) 
  lambda_EUVAC_A(17) = 554.37 
  lambda_EUVAC_A(18) = 584.33
  lambda_EUVAC_A(19) = 0.5 * (550.0 + 600.0)
  lambda_EUVAC_A(20) = 609.76

  lambda_EUVAC_A(21) = 629.73
  lambda_EUVAC_A(22) = 0.5 * (600.0 + 650.0)
  lambda_EUVAC_A(23) = 0.5 * (650.0 + 700.0) 
  lambda_EUVAC_A(24) = 703.36 
  lambda_EUVAC_A(25) = 0.5 * (700.0 + 750.0) 
  lambda_EUVAC_A(26) = 765.15 
  lambda_EUVAC_A(27) = 770.41 
  lambda_EUVAC_A(28) = 789.36
  lambda_EUVAC_A(29) = 0.5 * (750.0 + 800.0)
  lambda_EUVAC_A(30) = 0.5 * (800.0 + 850.0) 

  lambda_EUVAC_A(31) = 0.5 * (850.0 + 900.0)
  lambda_EUVAC_A(32) = 0.5 * (900.0 + 950.0) 
  lambda_EUVAC_A(33) = 977.02
  lambda_EUVAC_A(34) = 0.5 * (950.0 + 1000.0)
  lambda_EUVAC_A(35) = 1025.72
  lambda_EUVAC_A(36) = 1031.91
  lambda_EUVAC_A(37) = 0.5 * (1000.0 + 1050.0)

! EUVAC solar flux [in 10^9 ph cm^-2 s^-1]
! correct the solar flux due to the solar activity variation
! use Eqs.(9.19) and (9.20) from the Schunk and Nagy's book

  factor_f107 = 0.5 * (f10p7 + f10p7_81) - 80.0

  solar_flux_EUVAC_phcm2s1( 1) = MAX( 0.0, 1.200 * (1.0 + 1.0017e-2 * factor_f107) )
  solar_flux_EUVAC_phcm2s1( 2) = MAX( 0.0, 0.450 * (1.0 + 7.1250e-3 * factor_f107) )
  solar_flux_EUVAC_phcm2s1( 3) = MAX( 0.0, 4.800 * (1.0 + 1.3375e-2 * factor_f107) )
  solar_flux_EUVAC_phcm2s1( 4) = MAX( 0.0, 3.100 * (1.0 + 1.9450e-2 * factor_f107) )
  solar_flux_EUVAC_phcm2s1( 5) = MAX( 0.0, 0.460 * (1.0 + 2.7750e-3 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1( 6) = MAX( 0.0, 0.210 * (1.0 + 1.3768e-1 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1( 7) = MAX( 0.0, 1.679 * (1.0 + 2.6467e-2 * factor_f107) )
  solar_flux_EUVAC_phcm2s1( 8) = MAX( 0.0, 0.800 * (1.0 + 2.5000e-2 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1( 9) = MAX( 0.0, 6.900 * (1.0 + 3.3333e-3 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1(10) = MAX( 0.0, 0.965 * (1.0 + 2.2450e-2 * factor_f107) )

  solar_flux_EUVAC_phcm2s1(11) = MAX( 0.0, 0.650 * (1.0 + 6.5917e-3 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1(12) = MAX( 0.0, 0.314 * (1.0 + 3.6542e-2 * factor_f107) )
  solar_flux_EUVAC_phcm2s1(13) = MAX( 0.0, 0.383 * (1.0 + 7.4083e-3 * factor_f107) )
  solar_flux_EUVAC_phcm2s1(14) = MAX( 0.0, 0.290 * (1.0 + 7.4917e-3 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1(15) = MAX( 0.0, 0.285 * (1.0 + 2.0225e-2 * factor_f107) )
  solar_flux_EUVAC_phcm2s1(16) = MAX( 0.0, 0.452 * (1.0 + 8.7583e-3 * factor_f107) )
  solar_flux_EUVAC_phcm2s1(17) = MAX( 0.0, 0.720 * (1.0 + 3.2667e-3 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1(18) = MAX( 0.0, 1.270 * (1.0 + 5.1583e-3 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1(19) = MAX( 0.0, 0.357 * (1.0 + 3.6583e-3 * factor_f107) )
  solar_flux_EUVAC_phcm2s1(20) = MAX( 0.0, 0.530 * (1.0 + 1.6175e-2 * factor_f107) )   !### line

  solar_flux_EUVAC_phcm2s1(21) = MAX( 0.0, 1.590 * (1.0 + 3.3250e-3 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1(22) = MAX( 0.0, 0.342 * (1.0 + 1.1800e-2 * factor_f107) )
  solar_flux_EUVAC_phcm2s1(23) = MAX( 0.0, 0.230 * (1.0 + 4.2667e-3 * factor_f107) )
  solar_flux_EUVAC_phcm2s1(24) = MAX( 0.0, 0.360 * (1.0 + 3.0417e-3 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1(25) = MAX( 0.0, 0.141 * (1.0 + 4.7500e-3 * factor_f107) )
  solar_flux_EUVAC_phcm2s1(26) = MAX( 0.0, 0.170 * (1.0 + 3.8500e-3 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1(27) = MAX( 0.0, 0.260 * (1.0 + 1.2808e-2 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1(28) = MAX( 0.0, 0.702 * (1.0 + 3.2750e-3 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1(29) = MAX( 0.0, 0.758 * (1.0 + 4.7667e-3 * factor_f107) )
  solar_flux_EUVAC_phcm2s1(30) = MAX( 0.0, 1.625 * (1.0 + 4.8167e-3 * factor_f107) )

  solar_flux_EUVAC_phcm2s1(31) = MAX( 0.0, 3.537 * (1.0 + 5.6750e-3 * factor_f107) )
  solar_flux_EUVAC_phcm2s1(32) = MAX( 0.0, 3.000 * (1.0 + 4.9833e-3 * factor_f107) )
  solar_flux_EUVAC_phcm2s1(33) = MAX( 0.0, 4.400 * (1.0 + 3.9417e-3 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1(34) = MAX( 0.0, 1.475 * (1.0 + 4.4167e-3 * factor_f107) )
  solar_flux_EUVAC_phcm2s1(35) = MAX( 0.0, 3.500 * (1.0 + 5.1833e-3 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1(36) = MAX( 0.0, 2.100 * (1.0 + 5.2833e-3 * factor_f107) )   !### line
  solar_flux_EUVAC_phcm2s1(37) = MAX( 0.0, 2.467 * (1.0 + 4.3750e-3 * factor_f107) )

! total absorption cross section for atomic oxygen [in 10^-18 cm^2]

  sigma_EUVAC_tot_O_cm2( 1) = 0.730
  sigma_EUVAC_tot_O_cm2( 2) = 1.839
  sigma_EUVAC_tot_O_cm2( 3) = 3.732
  sigma_EUVAC_tot_O_cm2( 4) = 5.202
  sigma_EUVAC_tot_O_cm2( 5) = 6.050
  sigma_EUVAC_tot_O_cm2( 6) = 7.080
  sigma_EUVAC_tot_O_cm2( 7) = 6.461
  sigma_EUVAC_tot_O_cm2( 8) = 7.680
  sigma_EUVAC_tot_O_cm2( 9) = 7.700
  sigma_EUVAC_tot_O_cm2(10) = 8.693

  sigma_EUVAC_tot_O_cm2(11) = 9.840
  sigma_EUVAC_tot_O_cm2(12) = 9.687
  sigma_EUVAC_tot_O_cm2(13) = 11.496
  sigma_EUVAC_tot_O_cm2(14) = 11.930
  sigma_EUVAC_tot_O_cm2(15) = 12.127
  sigma_EUVAC_tot_O_cm2(16) = 12.059
  sigma_EUVAC_tot_O_cm2(17) = 12.590
  sigma_EUVAC_tot_O_cm2(18) = 13.090
  sigma_EUVAC_tot_O_cm2(19) = 13.024
  sigma_EUVAC_tot_O_cm2(20) = 13.400

  sigma_EUVAC_tot_O_cm2(21) = 13.400
  sigma_EUVAC_tot_O_cm2(22) = 13.365
  sigma_EUVAC_tot_O_cm2(23) = 17.245
  sigma_EUVAC_tot_O_cm2(24) = 11.460
  sigma_EUVAC_tot_O_cm2(25) = 10.736
  sigma_EUVAC_tot_O_cm2(26) = 4.000
  sigma_EUVAC_tot_O_cm2(27) = 3.890
  sigma_EUVAC_tot_O_cm2(28) = 3.749
  sigma_EUVAC_tot_O_cm2(29) = 5.091
  sigma_EUVAC_tot_O_cm2(30) = 3.498

  sigma_EUVAC_tot_O_cm2(31) = 4.554
  sigma_EUVAC_tot_O_cm2(32) = 1.315
  sigma_EUVAC_tot_O_cm2(33) = 0.0
  sigma_EUVAC_tot_O_cm2(34) = 0.0  
  sigma_EUVAC_tot_O_cm2(35) = 0.0
  sigma_EUVAC_tot_O_cm2(36) = 0.0
  sigma_EUVAC_tot_O_cm2(37) = 0.0

! total absorption cross section for molecular oxygen [in 10^-18 cm^2]

  sigma_EUVAC_tot_O2_cm2( 1) = 1.316
  sigma_EUVAC_tot_O2_cm2( 2) = 3.806
  sigma_EUVAC_tot_O2_cm2( 3) = 7.509
  sigma_EUVAC_tot_O2_cm2( 4) = 10.900
  sigma_EUVAC_tot_O2_cm2( 5) = 13.370
  sigma_EUVAC_tot_O2_cm2( 6) = 15.790
  sigma_EUVAC_tot_O2_cm2( 7) = 14.387
  sigma_EUVAC_tot_O2_cm2( 8) = 16.800
  sigma_EUVAC_tot_O2_cm2( 9) = 16.810
  sigma_EUVAC_tot_O2_cm2(10) = 17.438

  sigma_EUVAC_tot_O2_cm2(11) = 18.320
  sigma_EUVAC_tot_O2_cm2(12) = 18.118
  sigma_EUVAC_tot_O2_cm2(13) = 20.310
  sigma_EUVAC_tot_O2_cm2(14) = 21.910
  sigma_EUVAC_tot_O2_cm2(15) = 23.101
  sigma_EUVAC_tot_O2_cm2(16) = 24.606
  sigma_EUVAC_tot_O2_cm2(17) = 26.040
  sigma_EUVAC_tot_O2_cm2(18) = 22.720
  sigma_EUVAC_tot_O2_cm2(19) = 26.610
  sigma_EUVAC_tot_O2_cm2(20) = 28.070

  sigma_EUVAC_tot_O2_cm2(21) = 32.060
  sigma_EUVAC_tot_O2_cm2(22) = 26.017
  sigma_EUVAC_tot_O2_cm2(23) = 21.919
  sigma_EUVAC_tot_O2_cm2(24) = 27.440
  sigma_EUVAC_tot_O2_cm2(25) = 28.535
  sigma_EUVAC_tot_O2_cm2(26) = 20.800
  sigma_EUVAC_tot_O2_cm2(27) = 18.910
  sigma_EUVAC_tot_O2_cm2(28) = 26.668
  sigma_EUVAC_tot_O2_cm2(29) = 22.145
  sigma_EUVAC_tot_O2_cm2(30) = 16.631

  sigma_EUVAC_tot_O2_cm2(31) = 8.562
  sigma_EUVAC_tot_O2_cm2(32) = 12.817
  sigma_EUVAC_tot_O2_cm2(33) = 18.730
  sigma_EUVAC_tot_O2_cm2(34) = 21.108
  sigma_EUVAC_tot_O2_cm2(35) = 1.630
  sigma_EUVAC_tot_O2_cm2(36) = 1.050
  sigma_EUVAC_tot_O2_cm2(37) = 1.346

! total absorption cross section for molecular nitrogen [in 10^-18 cm^2]

  sigma_EUVAC_tot_N2_cm2( 1) = 0.720
  sigma_EUVAC_tot_N2_cm2( 2) = 2.261
  sigma_EUVAC_tot_N2_cm2( 3) = 4.958
  sigma_EUVAC_tot_N2_cm2( 4) = 8.392
  sigma_EUVAC_tot_N2_cm2( 5) = 10.210
  sigma_EUVAC_tot_N2_cm2( 6) = 10.900
  sigma_EUVAC_tot_N2_cm2( 7) = 10.493
  sigma_EUVAC_tot_N2_cm2( 8) = 11.670
  sigma_EUVAC_tot_N2_cm2( 9) = 11.700
  sigma_EUVAC_tot_N2_cm2(10) = 13.857

  sigma_EUVAC_tot_N2_cm2(11) = 16.910
  sigma_EUVAC_tot_N2_cm2(12) = 16.395
  sigma_EUVAC_tot_N2_cm2(13) = 21.675
  sigma_EUVAC_tot_N2_cm2(14) = 23.160
  sigma_EUVAC_tot_N2_cm2(15) = 23.471
  sigma_EUVAC_tot_N2_cm2(16) = 24.501
  sigma_EUVAC_tot_N2_cm2(17) = 24.130
  sigma_EUVAC_tot_N2_cm2(18) = 22.400
  sigma_EUVAC_tot_N2_cm2(19) = 22.787
  sigma_EUVAC_tot_N2_cm2(20) = 22.790

  sigma_EUVAC_tot_N2_cm2(21) = 23.370
  sigma_EUVAC_tot_N2_cm2(22) = 23.339
  sigma_EUVAC_tot_N2_cm2(23) = 31.755
  sigma_EUVAC_tot_N2_cm2(24) = 26.540
  sigma_EUVAC_tot_N2_cm2(25) = 24.662
  sigma_EUVAC_tot_N2_cm2(26) = 120.490
  sigma_EUVAC_tot_N2_cm2(27) = 14.180
  sigma_EUVAC_tot_N2_cm2(28) = 16.487
  sigma_EUVAC_tot_N2_cm2(29) = 33.578
  sigma_EUVAC_tot_N2_cm2(30) = 16.992

  sigma_EUVAC_tot_N2_cm2(31) = 20.249
  sigma_EUVAC_tot_N2_cm2(32) = 9.680
  sigma_EUVAC_tot_N2_cm2(33) = 2.240
  sigma_EUVAC_tot_N2_cm2(34) = 50.988
  sigma_EUVAC_tot_N2_cm2(35) = 0.0
  sigma_EUVAC_tot_N2_cm2(36) = 0.0
  sigma_EUVAC_tot_N2_cm2(37) = 0.0

! ionization cross section for molecular nitrogen (N2 -> N2+) [in 10^-18 cm^2]

  sigma_EUVAC_ion_N2_to_N2_cm2( 1) = 0.443
  sigma_EUVAC_ion_N2_to_N2_cm2( 2) = 1.479
  sigma_EUVAC_ion_N2_to_N2_cm2( 3) = 3.153
  sigma_EUVAC_ion_N2_to_N2_cm2( 4) = 5.226
  sigma_EUVAC_ion_N2_to_N2_cm2( 5) = 6.781
  sigma_EUVAC_ion_N2_to_N2_cm2( 6) = 8.100
  sigma_EUVAC_ion_N2_to_N2_cm2( 7) = 7.347
  sigma_EUVAC_ion_N2_to_N2_cm2( 8) = 9.180
  sigma_EUVAC_ion_N2_to_N2_cm2( 9) = 9.210
  sigma_EUVAC_ion_N2_to_N2_cm2(10) = 11.600

  sigma_EUVAC_ion_N2_to_N2_cm2(11) = 15.350
  sigma_EUVAC_ion_N2_to_N2_cm2(12) = 14.669
  sigma_EUVAC_ion_N2_to_N2_cm2(13) = 20.692
  sigma_EUVAC_ion_N2_to_N2_cm2(14) = 22.100
  sigma_EUVAC_ion_N2_to_N2_cm2(15) = 22.772
  sigma_EUVAC_ion_N2_to_N2_cm2(16) = 24.468
  sigma_EUVAC_ion_N2_to_N2_cm2(17) = 24.130
  sigma_EUVAC_ion_N2_to_N2_cm2(18) = 22.400
  sigma_EUVAC_ion_N2_to_N2_cm2(19) = 22.787
  sigma_EUVAC_ion_N2_to_N2_cm2(20) = 22.790

  sigma_EUVAC_ion_N2_to_N2_cm2(21) = 23.370
  sigma_EUVAC_ion_N2_to_N2_cm2(22) = 23.339
  sigma_EUVAC_ion_N2_to_N2_cm2(23) = 29.235
  sigma_EUVAC_ion_N2_to_N2_cm2(24) = 25.480
  sigma_EUVAC_ion_N2_to_N2_cm2(25) = 15.060
  sigma_EUVAC_ion_N2_to_N2_cm2(26) = 65.800
  sigma_EUVAC_ion_N2_to_N2_cm2(27) = 8.500
  sigma_EUVAC_ion_N2_to_N2_cm2(28) = 8.860
  sigma_EUVAC_ion_N2_to_N2_cm2(29) = 14.274
  sigma_EUVAC_ion_N2_to_N2_cm2(30) = 0.0

  sigma_EUVAC_ion_N2_to_N2_cm2(31) = 0.0
  sigma_EUVAC_ion_N2_to_N2_cm2(32) = 0.0
  sigma_EUVAC_ion_N2_to_N2_cm2(33) = 0.0
  sigma_EUVAC_ion_N2_to_N2_cm2(34) = 0.0
  sigma_EUVAC_ion_N2_to_N2_cm2(35) = 0.0
  sigma_EUVAC_ion_N2_to_N2_cm2(36) = 0.0
  sigma_EUVAC_ion_N2_to_N2_cm2(37) = 0.0

! ionization cross section for molecular nitrogen (N2 -> N+) [in 10^-18 cm^2]

  sigma_EUVAC_ion_N2_to_N_cm2( 1) = 0.277
  sigma_EUVAC_ion_N2_to_N_cm2( 2) = 0.782
  sigma_EUVAC_ion_N2_to_N_cm2( 3) = 1.805
  sigma_EUVAC_ion_N2_to_N_cm2( 4) = 3.166
  sigma_EUVAC_ion_N2_to_N_cm2( 5) = 3.420
  sigma_EUVAC_ion_N2_to_N_cm2( 6) = 2.800
  sigma_EUVAC_ion_N2_to_N_cm2( 7) = 3.145
  sigma_EUVAC_ion_N2_to_N_cm2( 8) = 2.490
  sigma_EUVAC_ion_N2_to_N_cm2( 9) = 2.490
  sigma_EUVAC_ion_N2_to_N_cm2(10) = 2.257

  sigma_EUVAC_ion_N2_to_N_cm2(11) = 1.560
  sigma_EUVAC_ion_N2_to_N_cm2(12) = 1.726
  sigma_EUVAC_ion_N2_to_N_cm2(13) = 0.982
  sigma_EUVAC_ion_N2_to_N_cm2(14) = 1.060
  sigma_EUVAC_ion_N2_to_N_cm2(15) = 0.699  
  sigma_EUVAC_ion_N2_to_N_cm2(16) = 0.033
  sigma_EUVAC_ion_N2_to_N_cm2(17) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(18) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(19) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(20) = 0.0

  sigma_EUVAC_ion_N2_to_N_cm2(21) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(22) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(23) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(24) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(25) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(26) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(27) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(28) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(29) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(30) = 0.0

  sigma_EUVAC_ion_N2_to_N_cm2(31) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(32) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(33) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(34) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(35) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(36) = 0.0
  sigma_EUVAC_ion_N2_to_N_cm2(37) = 0.0

! ionization cross section for atomic oxygen (O -> O+4S) [in 10^-18 cm^2]

  sigma_EUVAC_ion_O_4S_cm2( 1) = 0.190
  sigma_EUVAC_ion_O_4S_cm2( 2) = 0.486
  sigma_EUVAC_ion_O_4S_cm2( 3) = 0.952
  sigma_EUVAC_ion_O_4S_cm2( 4) = 1.311
  sigma_EUVAC_ion_O_4S_cm2( 5) = 1.539
  sigma_EUVAC_ion_O_4S_cm2( 6) = 1.770
  sigma_EUVAC_ion_O_4S_cm2( 7) = 1.628
  sigma_EUVAC_ion_O_4S_cm2( 8) = 1.920
  sigma_EUVAC_ion_O_4S_cm2( 9) = 1.925
  sigma_EUVAC_ion_O_4S_cm2(10) = 2.259

  sigma_EUVAC_ion_O_4S_cm2(11) = 2.559
  sigma_EUVAC_ion_O_4S_cm2(12) = 2.523
  sigma_EUVAC_ion_O_4S_cm2(13) = 3.073
  sigma_EUVAC_ion_O_4S_cm2(14) = 3.340
  sigma_EUVAC_ion_O_4S_cm2(15) = 3.394
  sigma_EUVAC_ion_O_4S_cm2(16) = 3.421
  sigma_EUVAC_ion_O_4S_cm2(17) = 3.650
  sigma_EUVAC_ion_O_4S_cm2(18) = 3.920
  sigma_EUVAC_ion_O_4S_cm2(19) = 3.620
  sigma_EUVAC_ion_O_4S_cm2(20) = 3.610

  sigma_EUVAC_ion_O_4S_cm2(21) = 3.880
  sigma_EUVAC_ion_O_4S_cm2(22) = 4.250
  sigma_EUVAC_ion_O_4S_cm2(23) = 5.128
  sigma_EUVAC_ion_O_4S_cm2(24) = 4.890
  sigma_EUVAC_ion_O_4S_cm2(25) = 6.739
  sigma_EUVAC_ion_O_4S_cm2(26) = 4.000
  sigma_EUVAC_ion_O_4S_cm2(27) = 3.890
  sigma_EUVAC_ion_O_4S_cm2(28) = 3.749
  sigma_EUVAC_ion_O_4S_cm2(29) = 5.091
  sigma_EUVAC_ion_O_4S_cm2(30) = 3.498

  sigma_EUVAC_ion_O_4S_cm2(31) = 4.554
  sigma_EUVAC_ion_O_4S_cm2(32) = 1.315
  sigma_EUVAC_ion_O_4S_cm2(33) = 0.0
  sigma_EUVAC_ion_O_4S_cm2(34) = 0.0
  sigma_EUVAC_ion_O_4S_cm2(35) = 0.0
  sigma_EUVAC_ion_O_4S_cm2(36) = 0.0
  sigma_EUVAC_ion_O_4S_cm2(37) = 0.0

! ionization cross section for atomic oxygen (O -> O+2D) [in 10^-18 cm^2]

  sigma_EUVAC_ion_O_2D_cm2( 1) = 0.206
  sigma_EUVAC_ion_O_2D_cm2( 2) = 0.529
  sigma_EUVAC_ion_O_2D_cm2( 3) = 1.171
  sigma_EUVAC_ion_O_2D_cm2( 4) = 1.762
  sigma_EUVAC_ion_O_2D_cm2( 5) = 2.138
  sigma_EUVAC_ion_O_2D_cm2( 6) = 2.620
  sigma_EUVAC_ion_O_2D_cm2( 7) = 2.325
  sigma_EUVAC_ion_O_2D_cm2( 8) = 2.842
  sigma_EUVAC_ion_O_2D_cm2( 9) = 2.849
  sigma_EUVAC_ion_O_2D_cm2(10) = 3.446

  sigma_EUVAC_ion_O_2D_cm2(11) = 3.936
  sigma_EUVAC_ion_O_2D_cm2(12) = 3.883
  sigma_EUVAC_ion_O_2D_cm2(13) = 4.896
  sigma_EUVAC_ion_O_2D_cm2(14) = 5.370
  sigma_EUVAC_ion_O_2D_cm2(15) = 5.459
  sigma_EUVAC_ion_O_2D_cm2(16) = 5.427
  sigma_EUVAC_ion_O_2D_cm2(17) = 5.670
  sigma_EUVAC_ion_O_2D_cm2(18) = 6.020
  sigma_EUVAC_ion_O_2D_cm2(19) = 5.910
  sigma_EUVAC_ion_O_2D_cm2(20) = 6.170

  sigma_EUVAC_ion_O_2D_cm2(21) = 6.290
  sigma_EUVAC_ion_O_2D_cm2(22) = 6.159
  sigma_EUVAC_ion_O_2D_cm2(23) = 11.453
  sigma_EUVAC_ion_O_2D_cm2(24) = 6.570
  sigma_EUVAC_ion_O_2D_cm2(25) = 3.997
  sigma_EUVAC_ion_O_2D_cm2(26) = 0.0
  sigma_EUVAC_ion_O_2D_cm2(27) = 0.0
  sigma_EUVAC_ion_O_2D_cm2(28) = 0.0
  sigma_EUVAC_ion_O_2D_cm2(29) = 0.0
  sigma_EUVAC_ion_O_2D_cm2(30) = 0.0

  sigma_EUVAC_ion_O_2D_cm2(31) = 0.0
  sigma_EUVAC_ion_O_2D_cm2(32) = 0.0
  sigma_EUVAC_ion_O_2D_cm2(33) = 0.0
  sigma_EUVAC_ion_O_2D_cm2(34) = 0.0
  sigma_EUVAC_ion_O_2D_cm2(35) = 0.0
  sigma_EUVAC_ion_O_2D_cm2(36) = 0.0
  sigma_EUVAC_ion_O_2D_cm2(37) = 0.0

! ionization cross section for atomic oxygen (O -> O+2P) [in 10^-18 cm^2]

  sigma_EUVAC_ion_O_2P_cm2( 1) = 0.134
  sigma_EUVAC_ion_O_2P_cm2( 2) = 0.345
  sigma_EUVAC_ion_O_2P_cm2( 3) = 0.768
  sigma_EUVAC_ion_O_2P_cm2( 4) = 1.144
  sigma_EUVAC_ion_O_2P_cm2( 5) = 1.363
  sigma_EUVAC_ion_O_2P_cm2( 6) = 1.630
  sigma_EUVAC_ion_O_2P_cm2( 7) = 1.488
  sigma_EUVAC_ion_O_2P_cm2( 8) = 1.920
  sigma_EUVAC_ion_O_2P_cm2( 9) = 1.925
  sigma_EUVAC_ion_O_2P_cm2(10) = 2.173

  sigma_EUVAC_ion_O_2P_cm2(11) = 2.558
  sigma_EUVAC_ion_O_2P_cm2(12) = 2.422
  sigma_EUVAC_ion_O_2P_cm2(13) = 2.986
  sigma_EUVAC_ion_O_2P_cm2(14) = 3.220
  sigma_EUVAC_ion_O_2P_cm2(15) = 3.274
  sigma_EUVAC_ion_O_2P_cm2(16) = 3.211
  sigma_EUVAC_ion_O_2P_cm2(17) = 3.270
  sigma_EUVAC_ion_O_2P_cm2(18) = 3.150
  sigma_EUVAC_ion_O_2P_cm2(19) = 3.494
  sigma_EUVAC_ion_O_2P_cm2(20) = 3.620

  sigma_EUVAC_ion_O_2P_cm2(21) = 3.230
  sigma_EUVAC_ion_O_2P_cm2(22) = 2.956
  sigma_EUVAC_ion_O_2P_cm2(23) = 0.664
  sigma_EUVAC_ion_O_2P_cm2(24) = 0.0
  sigma_EUVAC_ion_O_2P_cm2(25) = 0.0
  sigma_EUVAC_ion_O_2P_cm2(26) = 0.0
  sigma_EUVAC_ion_O_2P_cm2(27) = 0.0
  sigma_EUVAC_ion_O_2P_cm2(28) = 0.0
  sigma_EUVAC_ion_O_2P_cm2(29) = 0.0
  sigma_EUVAC_ion_O_2P_cm2(30) = 0.0

  sigma_EUVAC_ion_O_2P_cm2(31) = 0.0
  sigma_EUVAC_ion_O_2P_cm2(32) = 0.0
  sigma_EUVAC_ion_O_2P_cm2(33) = 0.0
  sigma_EUVAC_ion_O_2P_cm2(34) = 0.0
  sigma_EUVAC_ion_O_2P_cm2(35) = 0.0
  sigma_EUVAC_ion_O_2P_cm2(36) = 0.0
  sigma_EUVAC_ion_O_2P_cm2(37) = 0.0

! ionization cross section for atomic oxygen (O -> O+4P*) [in 10^-18 cm^2]

  sigma_EUVAC_ion_O_4Pst_cm2( 1) = 0.062
  sigma_EUVAC_ion_O_4Pst_cm2( 2) = 0.163
  sigma_EUVAC_ion_O_4Pst_cm2( 3) = 0.348
  sigma_EUVAC_ion_O_4Pst_cm2( 4) = 0.508
  sigma_EUVAC_ion_O_4Pst_cm2( 5) = 0.598
  sigma_EUVAC_ion_O_4Pst_cm2( 6) = 0.710
  sigma_EUVAC_ion_O_4Pst_cm2( 7) = 0.637
  sigma_EUVAC_ion_O_4Pst_cm2( 8) = 0.691
  sigma_EUVAC_ion_O_4Pst_cm2( 9) = 0.693
  sigma_EUVAC_ion_O_4Pst_cm2(10) = 0.815

  sigma_EUVAC_ion_O_4Pst_cm2(11) = 0.787
  sigma_EUVAC_ion_O_4Pst_cm2(12) = 0.859
  sigma_EUVAC_ion_O_4Pst_cm2(13) = 0.541
  sigma_EUVAC_ion_O_4Pst_cm2(14) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(15) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(16) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(17) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(18) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(19) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(20) = 0.0

  sigma_EUVAC_ion_O_4Pst_cm2(21) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(22) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(23) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(24) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(25) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(26) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(27) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(28) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(29) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(30) = 0.0

  sigma_EUVAC_ion_O_4Pst_cm2(31) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(32) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(33) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(34) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(35) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(36) = 0.0
  sigma_EUVAC_ion_O_4Pst_cm2(37) = 0.0

! ionization cross section for atomic oxygen (O -> O+2P*) [in 10^-18 cm^2]

  sigma_EUVAC_ion_O_2Pst_cm2( 1) = 0.049
  sigma_EUVAC_ion_O_2Pst_cm2( 2) = 0.130
  sigma_EUVAC_ion_O_2Pst_cm2( 3) = 0.278
  sigma_EUVAC_ion_O_2Pst_cm2( 4) = 0.366
  sigma_EUVAC_ion_O_2Pst_cm2( 5) = 0.412
  sigma_EUVAC_ion_O_2Pst_cm2( 6) = 0.350
  sigma_EUVAC_ion_O_2Pst_cm2( 7) = 0.383
  sigma_EUVAC_ion_O_2Pst_cm2( 8) = 0.307
  sigma_EUVAC_ion_O_2Pst_cm2( 9) = 0.308
  sigma_EUVAC_ion_O_2Pst_cm2(10) = 0.0

  sigma_EUVAC_ion_O_2Pst_cm2(11) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(12) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(13) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(14) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(15) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(16) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(17) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(18) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(19) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(20) = 0.0

  sigma_EUVAC_ion_O_2Pst_cm2(21) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(22) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(23) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(24) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(25) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(26) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(27) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(28) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(29) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(30) = 0.0

  sigma_EUVAC_ion_O_2Pst_cm2(31) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(32) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(33) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(34) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(35) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(36) = 0.0
  sigma_EUVAC_ion_O_2Pst_cm2(37) = 0.0

! ionization cross section for molecular oxygen (O2 -> O2+) [in 10^-18 cm^2]

  sigma_EUVAC_ion_O2_to_O2_cm2( 1) = 1.316
  sigma_EUVAC_ion_O2_to_O2_cm2( 2) = 2.346
  sigma_EUVAC_ion_O2_to_O2_cm2( 3) = 4.139
  sigma_EUVAC_ion_O2_to_O2_cm2( 4) = 6.619
  sigma_EUVAC_ion_O2_to_O2_cm2( 5) = 8.460
  sigma_EUVAC_ion_O2_to_O2_cm2( 6) = 9.890
  sigma_EUVAC_ion_O2_to_O2_cm2( 7) = 9.056
  sigma_EUVAC_ion_O2_to_O2_cm2( 8) = 10.860
  sigma_EUVAC_ion_O2_to_O2_cm2( 9) = 10.880
  sigma_EUVAC_ion_O2_to_O2_cm2(10) = 12.229

  sigma_EUVAC_ion_O2_to_O2_cm2(11) = 13.760
  sigma_EUVAC_ion_O2_to_O2_cm2(12) = 13.418
  sigma_EUVAC_ion_O2_to_O2_cm2(13) = 15.490
  sigma_EUVAC_ion_O2_to_O2_cm2(14) = 16.970
  sigma_EUVAC_ion_O2_to_O2_cm2(15) = 17.754
  sigma_EUVAC_ion_O2_to_O2_cm2(16) = 19.469
  sigma_EUVAC_ion_O2_to_O2_cm2(17) = 21.600
  sigma_EUVAC_ion_O2_to_O2_cm2(18) = 18.840
  sigma_EUVAC_ion_O2_to_O2_cm2(19) = 22.789
  sigma_EUVAC_ion_O2_to_O2_cm2(20) = 24.540

  sigma_EUVAC_ion_O2_to_O2_cm2(21) = 30.070
  sigma_EUVAC_ion_O2_to_O2_cm2(22) = 23.974
  sigma_EUVAC_ion_O2_to_O2_cm2(23) = 21.116
  sigma_EUVAC_ion_O2_to_O2_cm2(24) = 23.750
  sigma_EUVAC_ion_O2_to_O2_cm2(25) = 23.805
  sigma_EUVAC_ion_O2_to_O2_cm2(26) = 11.720
  sigma_EUVAC_ion_O2_to_O2_cm2(27) = 8.470
  sigma_EUVAC_ion_O2_to_O2_cm2(28) = 10.191
  sigma_EUVAC_ion_O2_to_O2_cm2(29) = 10.597
  sigma_EUVAC_ion_O2_to_O2_cm2(30) = 6.413

  sigma_EUVAC_ion_O2_to_O2_cm2(31) = 5.494
  sigma_EUVAC_ion_O2_to_O2_cm2(32) = 9.374
  sigma_EUVAC_ion_O2_to_O2_cm2(33) = 15.540
  sigma_EUVAC_ion_O2_to_O2_cm2(34) = 13.940
  sigma_EUVAC_ion_O2_to_O2_cm2(35) = 1.050
  sigma_EUVAC_ion_O2_to_O2_cm2(36) = 0.000
  sigma_EUVAC_ion_O2_to_O2_cm2(37) = 0.259

! ionization cross section for molecular oxygen (O2 -> O+) [in 10^-18 cm^2]

  sigma_EUVAC_ion_O2_to_O_cm2( 1) = 0.0
  sigma_EUVAC_ion_O2_to_O_cm2( 2) = 1.460
  sigma_EUVAC_ion_O2_to_O_cm2( 3) = 3.368
  sigma_EUVAC_ion_O2_to_O_cm2( 4) = 4.281
  sigma_EUVAC_ion_O2_to_O_cm2( 5) = 4.910
  sigma_EUVAC_ion_O2_to_O_cm2( 6) = 5.900
  sigma_EUVAC_ion_O2_to_O_cm2( 7) = 5.332
  sigma_EUVAC_ion_O2_to_O_cm2( 8) = 5.940
  sigma_EUVAC_ion_O2_to_O_cm2( 9) = 5.930
  sigma_EUVAC_ion_O2_to_O_cm2(10) = 5.212

  sigma_EUVAC_ion_O2_to_O_cm2(11) = 4.560
  sigma_EUVAC_ion_O2_to_O_cm2(12) = 4.703
  sigma_EUVAC_ion_O2_to_O_cm2(13) = 4.818
  sigma_EUVAC_ion_O2_to_O_cm2(14) = 4.940
  sigma_EUVAC_ion_O2_to_O_cm2(15) = 5.347
  sigma_EUVAC_ion_O2_to_O_cm2(16) = 5.139
  sigma_EUVAC_ion_O2_to_O_cm2(17) = 4.440
  sigma_EUVAC_ion_O2_to_O_cm2(18) = 3.880
  sigma_EUVAC_ion_O2_to_O_cm2(19) = 3.824
  sigma_EUVAC_ion_O2_to_O_cm2(20) = 1.850

  sigma_EUVAC_ion_O2_to_O_cm2(21) = 1.030
  sigma_EUVAC_ion_O2_to_O_cm2(22) = 0.962
  sigma_EUVAC_ion_O2_to_O_cm2(23) = 0.190
  sigma_EUVAC_ion_O2_to_O_cm2(24) = 0.0
  sigma_EUVAC_ion_O2_to_O_cm2(25) = 0.0
  sigma_EUVAC_ion_O2_to_O_cm2(26) = 0.0
  sigma_EUVAC_ion_O2_to_O_cm2(27) = 0.0
  sigma_EUVAC_ion_O2_to_O_cm2(28) = 0.0
  sigma_EUVAC_ion_O2_to_O_cm2(29) = 0.0
  sigma_EUVAC_ion_O2_to_O_cm2(30) = 0.0

  sigma_EUVAC_ion_O2_to_O_cm2(31) = 0.0
  sigma_EUVAC_ion_O2_to_O_cm2(32) = 0.0
  sigma_EUVAC_ion_O2_to_O_cm2(33) = 0.0
  sigma_EUVAC_ion_O2_to_O_cm2(34) = 0.0
  sigma_EUVAC_ion_O2_to_O_cm2(35) = 0.0
  sigma_EUVAC_ion_O2_to_O_cm2(36) = 0.0
  sigma_EUVAC_ion_O2_to_O_cm2(37) = 0.0

  factor_energy_eVA =  REAL(h_Planck_Js * c_ms * 1.0d10 / e_Cl)   ! 1.0d10 because the wavelength is in Angstrom

  d_lambda_A = 50.0 / N_split_EUVAC  

  ALLOCATE(lambda_EUV_A(N_solar_bins), STAT = ALLOC_ERR)

! the solar spectrum created on the basis of EUVAC has 20*N_split_EUVAC fluxes 
! corresponding to the continuous part of EUVAC, that is the 20 50-Angstrom-wide bins
!
! and 17 fluxes corresponding to the discrete 17 EUVAC lines ::
!
! N_solar_bins = 20 * N_split_EUVAC + 17
!
! in arrays, the lines are placed after all re-binned continuous fluxes
!

! continuous part of the spectrum, re-binned

  i=0
  DO n = 1, 20

     SELECT CASE (n)
        CASE (1)   !  50-100 A
           i_EUVAC = 1
        CASE (2)   ! 100-150 A
           i_EUVAC = 2
        CASE (3)   ! 150-200 A
           i_EUVAC = 3
        CASE (4)   ! 200-250 A
           i_EUVAC = 4
        CASE (5)   ! 250-300 A
           i_EUVAC = 7
        CASE (6)   ! 300-350 A
           i_EUVAC = 10
        CASE (7)   ! 350-400 A
           i_EUVAC = 12
        CASE (8)   ! 400-450 A
           i_EUVAC = 13
        CASE (9)   ! 450-500 A
           i_EUVAC = 15
        CASE (10)  ! 500-550 A
           i_EUVAC = 16
        CASE (11)  ! 550-600 A
           i_EUVAC = 19
        CASE (12)  ! 600-650 A
           i_EUVAC = 22
        CASE (13)  ! 650-700 A
           i_EUVAC = 23
        CASE (14)  ! 700-750 A
           i_EUVAC = 25
        CASE (15)  ! 750-800 A
           i_EUVAC = 29
        CASE (16)  ! 800-850 A
           i_EUVAC = 30
        CASE (17)  ! 850-900 A
           i_EUVAC = 31
        CASE (18)  ! 900-950 A
           i_EUVAC = 32
        CASE (19)  ! 950-1000 A
           i_EUVAC = 34
        CASE (20)  ! 1000-1050 A
           i_EUVAC = 37
     END SELECT

     DO i_split = 1, N_split_EUVAC

        i = i+1

        lambda_EUV_A(i) = lambda_EUVAC_A(i_EUVAC) - 25.0 + (REAL(i_split)-0.5) * d_lambda_A                                       ! modified

        solar_flux_energy_bin_eV(i) = factor_energy_eVA / lambda_EUV_A(i)

        solar_flux_phcm2s1(i) = solar_flux_EUVAC_phcm2s1(i_EUVAC) * (lambda_EUV_A(i) / lambda_EUVAC_A(i_EUVAC)) / N_split_EUVAC   ! modified

        sigma_tot_O_cm2(i)  = sigma_EUVAC_tot_O_cm2(i_EUVAC)
        sigma_tot_N2_cm2(i) = sigma_EUVAC_tot_N2_cm2(i_EUVAC)
        sigma_tot_O2_cm2(i) = sigma_EUVAC_tot_O2_cm2(i_EUVAC)

        sigma_ion_N2_to_N2_cm2(i) = sigma_EUVAC_ion_N2_to_N2_cm2(i_EUVAC)
        sigma_ion_N2_to_N_cm2(i)  = sigma_EUVAC_ion_N2_to_N_cm2(i_EUVAC)

        sigma_ion_O_4S_cm2(i)   = sigma_EUVAC_ion_O_4S_cm2(i_EUVAC)
        sigma_ion_O_2D_cm2(i)   = sigma_EUVAC_ion_O_2D_cm2(i_EUVAC)
        sigma_ion_O_2P_cm2(i)   = sigma_EUVAC_ion_O_2P_cm2(i_EUVAC)
        sigma_ion_O_4Pst_cm2(i) = sigma_EUVAC_ion_O_4Pst_cm2(i_EUVAC)
        sigma_ion_O_2Pst_cm2(i) = sigma_EUVAC_ion_O_2Pst_cm2(i_EUVAC)

        sigma_ion_O2_to_O2_cm2(i) = sigma_EUVAC_ion_O2_to_O2_cm2(i_EUVAC)
        sigma_ion_O2_to_O_cm2(i)  = sigma_EUVAC_ion_O2_to_O_cm2(i_EUVAC)

     END DO   !###   DO i_split = 1, N_split_EUVAC
  END DO   !###   DO n = 1, 20

! lines

  DO i = 20 * N_split_EUVAC + 1, N_solar_bins

     IF (i.EQ.(N_solar_bins-16)) THEN    ! 256.32 A
        i_EUVAC = 5
     ELSE IF (i.EQ.(N_solar_bins-15)) THEN    ! 284.15 A
        i_EUVAC = 6
     ELSE IF (i.EQ.(N_solar_bins-14)) THEN    ! 303.31 A
        i_EUVAC = 8
     ELSE IF (i.EQ.(N_solar_bins-13)) THEN    ! 303.78 A
        i_EUVAC = 9
     ELSE IF (i.EQ.(N_solar_bins-12)) THEN    ! 368.07 A
        i_EUVAC = 11
     ELSE IF (i.EQ.(N_solar_bins-11)) THEN    ! 465.22 A
        i_EUVAC = 14
     ELSE IF (i.EQ.(N_solar_bins-10)) THEN    ! 554.37 A
        i_EUVAC = 17
     ELSE IF (i.EQ.(N_solar_bins- 9)) THEN    ! 584.33 A
        i_EUVAC = 18
     ELSE IF (i.EQ.(N_solar_bins- 8)) THEN    ! 609.76 A
        i_EUVAC = 20
     ELSE IF (i.EQ.(N_solar_bins- 7)) THEN    ! 629.73 A
        i_EUVAC = 21
     ELSE IF (i.EQ.(N_solar_bins- 6)) THEN    ! 703.36 A
        i_EUVAC = 24
     ELSE IF (i.EQ.(N_solar_bins- 5)) THEN    ! 765.15 A
        i_EUVAC = 26
     ELSE IF (i.EQ.(N_solar_bins- 4)) THEN    ! 770.41 A
        i_EUVAC = 27
     ELSE IF (i.EQ.(N_solar_bins- 3)) THEN    ! 789.36 A
        i_EUVAC = 28
     ELSE IF (i.EQ.(N_solar_bins- 2)) THEN    ! 977.02 A
        i_EUVAC = 33
     ELSE IF (i.EQ.(N_solar_bins- 1)) THEN    ! 1025.72 A
        i_EUVAC = 35
     ELSE IF (i.EQ.(N_solar_bins   )) THEN    ! 1031.91 A
        i_EUVAC = 36
     END IF

     lambda_EUV_A(i) = lambda_EUVAC_A(i_EUVAC)                   ! EUVAC value, unchanged

     solar_flux_energy_bin_eV(i) = factor_energy_eVA / lambda_EUV_A(i)

     solar_flux_phcm2s1(i) = solar_flux_EUVAC_phcm2s1(i_EUVAC)   ! EUVAC value, unchanged

     sigma_tot_O_cm2(i)  = sigma_EUVAC_tot_O_cm2(i_EUVAC)
     sigma_tot_N2_cm2(i) = sigma_EUVAC_tot_N2_cm2(i_EUVAC)
     sigma_tot_O2_cm2(i) = sigma_EUVAC_tot_O2_cm2(i_EUVAC)

     sigma_ion_N2_to_N2_cm2(i) = sigma_EUVAC_ion_N2_to_N2_cm2(i_EUVAC)
     sigma_ion_N2_to_N_cm2(i)  = sigma_EUVAC_ion_N2_to_N_cm2(i_EUVAC)

     sigma_ion_O_4S_cm2(i)   = sigma_EUVAC_ion_O_4S_cm2(i_EUVAC)
     sigma_ion_O_2D_cm2(i)   = sigma_EUVAC_ion_O_2D_cm2(i_EUVAC)
     sigma_ion_O_2P_cm2(i)   = sigma_EUVAC_ion_O_2P_cm2(i_EUVAC)
     sigma_ion_O_4Pst_cm2(i) = sigma_EUVAC_ion_O_4Pst_cm2(i_EUVAC)
     sigma_ion_O_2Pst_cm2(i) = sigma_EUVAC_ion_O_2Pst_cm2(i_EUVAC)

     sigma_ion_O2_to_O2_cm2(i) = sigma_EUVAC_ion_O2_to_O2_cm2(i_EUVAC)
     sigma_ion_O2_to_O_cm2(i)  = sigma_EUVAC_ion_O2_to_O_cm2(i_EUVAC)

  END DO

!###########

  IF (Rank_of_process.EQ.0) THEN

     OPEN (10, FILE = 'solar_EUVAC_cont_fluxes_crossections.dat')

     WRITE (10, '("# col  1 is the ordering number of the wavelength bin")')
     WRITE (10, '("# col  2 is the wavelength (A)")')
     WRITE (10, '("# col  3 is the energy of photons with this wavelength (eV)")')

     WRITE (10, '("# col  4 is the photon flux of the bin (10^9 photons/cm^2/s)")')

     WRITE (10, '("# col  5 is the total photon absorption cross-section for  O (10^{-18} cm^2)")')
     WRITE (10, '("# col  6 is the total photon absorption cross-section for N2 (10^{-18} cm^2)")')
     WRITE (10, '("# col  7 is the total photon absorption cross-section for O2 (10^{-18} cm^2)")')

     WRITE (10, '("# col  8 is the ionization cross-section for O+(4S)  from O (10^{-18} cm^2)")')
     WRITE (10, '("# col  9 is the ionization cross-section for O+(2D)  from O (10^{-18} cm^2)")')
     WRITE (10, '("# col 10 is the ionization cross-section for O+(2P)  from O (10^{-18} cm^2)")')
     WRITE (10, '("# col 11 is the ionization cross-section for O+(4P*) from O (10^{-18} cm^2)")')
     WRITE (10, '("# col 12 is the ionization cross-section for O+(2P*) from O (10^{-18} cm^2)")')

     WRITE (10, '("# col 13 is the ionization cross-section for N2+ from N2 (10^{-18} cm^2)")')
     WRITE (10, '("# col 14 is the ionization cross-section for N+  from N2 (10^{-18} cm^2)")')

     WRITE (10, '("# col 15 is the ionization cross-section for O2+ from O2 (10^{-18} cm^2)")')
     WRITE (10, '("# col 16 is the ionization cross-section for O+  from O2 (10^{-18} cm^2)")')

     DO i = 1, N_solar_bins-17

        WRITE (10, '(2x,i4,2(2x,f8.3),2x,e12.5,12(2x,f9.3))') &
             & i, &           ! 1
             & lambda_EUV_A(i), &                  ! 2
             & solar_flux_energy_bin_eV(i), &      ! 3
             & solar_flux_phcm2s1(i), &        ! 4
             & sigma_tot_O_cm2(i), &               ! 5
             & sigma_tot_N2_cm2(i), &              ! 6
             & sigma_tot_O2_cm2(i), &              ! 7
             & sigma_ion_O_4S_cm2(i), &        ! 8
             & sigma_ion_O_2D_cm2(i), &        ! 9
             & sigma_ion_O_2P_cm2(i), &        ! 10
             & sigma_ion_O_4Pst_cm2(i), &      ! 11
             & sigma_ion_O_2Pst_cm2(i), &      ! 12
             & sigma_ion_N2_to_N2_cm2(i), &        ! 13
             & sigma_ion_N2_to_N_cm2(i), &         ! 14
             & sigma_ion_O2_to_O2_cm2(i), &    ! 15
             & sigma_ion_O2_to_O_cm2(i)        ! 16
     END DO
     CLOSE (10, STATUS = 'KEEP')

     PRINT '("### file solar_EUVAC_cont_fluxes_crossections.dat is ready")'

!--------------------------------

     OPEN (10, FILE = 'solar_EUVAC_lines_fluxes_crossections.dat')

     WRITE (10, '("# col  1 is the ordering number of the wavelength bin")')
     WRITE (10, '("# col  2 is the wavelength (A)")')
     WRITE (10, '("# col  3 is the energy of photons with this wavelength (eV)")')

     WRITE (10, '("# col  4 is the photon flux of the bin (10^9 photons/cm^2/s)")')

     WRITE (10, '("# col  5 is the total photon absorption cross-section for  O (10^{-18} cm^2)")')
     WRITE (10, '("# col  6 is the total photon absorption cross-section for N2 (10^{-18} cm^2)")')
     WRITE (10, '("# col  7 is the total photon absorption cross-section for O2 (10^{-18} cm^2)")')

     WRITE (10, '("# col  8 is the ionization cross-section for O+(4S)  from O (10^{-18} cm^2)")')
     WRITE (10, '("# col  9 is the ionization cross-section for O+(2D)  from O (10^{-18} cm^2)")')
     WRITE (10, '("# col 10 is the ionization cross-section for O+(2P)  from O (10^{-18} cm^2)")')
     WRITE (10, '("# col 11 is the ionization cross-section for O+(4P*) from O (10^{-18} cm^2)")')
     WRITE (10, '("# col 12 is the ionization cross-section for O+(2P*) from O (10^{-18} cm^2)")')

     WRITE (10, '("# col 13 is the ionization cross-section for N2+ from N2 (10^{-18} cm^2)")')
     WRITE (10, '("# col 14 is the ionization cross-section for N+  from N2 (10^{-18} cm^2)")')

     WRITE (10, '("# col 15 is the ionization cross-section for O2+ from O2 (10^{-18} cm^2)")')
     WRITE (10, '("# col 16 is the ionization cross-section for O+  from O2 (10^{-18} cm^2)")')

     DO i = N_solar_bins-16, N_solar_bins

        WRITE (10, '(2x,i4,2(2x,f8.3),2x,e12.5,12(2x,f9.3))') &
             & i, &           ! 1
             & lambda_EUV_A(i), &                  ! 2
             & solar_flux_energy_bin_eV(i), &      ! 3
             & solar_flux_phcm2s1(i), &        ! 4
             & sigma_tot_O_cm2(i), &               ! 5
             & sigma_tot_N2_cm2(i), &              ! 6
             & sigma_tot_O2_cm2(i), &              ! 7
             & sigma_ion_O_4S_cm2(i), &        ! 8
             & sigma_ion_O_2D_cm2(i), &        ! 9
             & sigma_ion_O_2P_cm2(i), &        ! 10
             & sigma_ion_O_4Pst_cm2(i), &      ! 11
             & sigma_ion_O_2Pst_cm2(i), &      ! 12
             & sigma_ion_N2_to_N2_cm2(i), &        ! 13
             & sigma_ion_N2_to_N_cm2(i), &         ! 14
             & sigma_ion_O2_to_O2_cm2(i), &    ! 15
             & sigma_ion_O2_to_O_cm2(i)        ! 16
        WRITE (10, '(2x,i4,2(2x,f8.3),2x,e12.5,12(2x,f9.3))') &
             & i, &           ! 1
             & lambda_EUV_A(i), &                  ! 2
             & solar_flux_energy_bin_eV(i), &      ! 3
             & 1.0d-10 * solar_flux_phcm2s1(i), &        ! 4
             & 1.0d-10 * sigma_tot_O_cm2(i), &               ! 5
             & 1.0d-10 * sigma_tot_N2_cm2(i), &              ! 6
             & 1.0d-10 * sigma_tot_O2_cm2(i), &              ! 7
             & 1.0d-10 * sigma_ion_O_4S_cm2(i), &        ! 8
             & 1.0d-10 * sigma_ion_O_2D_cm2(i), &        ! 9
             & 1.0d-10 * sigma_ion_O_2P_cm2(i), &        ! 10
             & 1.0d-10 * sigma_ion_O_4Pst_cm2(i), &      ! 11
             & 1.0d-10 * sigma_ion_O_2Pst_cm2(i), &      ! 12
             & 1.0d-10 * sigma_ion_N2_to_N2_cm2(i), &        ! 13
             & 1.0d-10 * sigma_ion_N2_to_N_cm2(i), &         ! 14
             & 1.0d-10 * sigma_ion_O2_to_O2_cm2(i), &    ! 15
             & 1.0d-10 * sigma_ion_O2_to_O_cm2(i)        ! 16
        WRITE (10, '(" ")')
        WRITE (10, '(" ")')

     END DO
     CLOSE (10, STATUS = 'KEEP')

     PRINT '("### file solar_EUVAC_lines_fluxes_crossections.dat is ready")'

  END IF   !###   IF (Rank_of_process.GT.0) THEN

! cleanup

  DEALLOCATE(lambda_EUV_A, STAT = ALLOC_ERR)

END SUBROUTINE Set_EUVAC
