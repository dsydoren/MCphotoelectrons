
!-------------------------------------------------
! all EUVAC parameters, including the cross-sections, are from Appendix J of 
! "Ionospheres. Physics, Plasma Physics, and Chemistry", Second edition, by Robert Schunk and Andrew Nagy
! Cambridge University Press, 2009
!
SUBROUTINE Initiate_EUVAC_flux_crossections

  USE ParallelOperationValues
  USE PhysicalConstants
  USE Photoelectrons
  USE GlobalIndices

  IMPLICIT NONE

  REAL lambda_euvac_A(N_solar_bins)

  REAL factor_energy_eVA
  INTEGER i

  REAL factor_f107

  REAL subbin_width_A

  INTEGER rank, n

  REAL lambda_mid_A,   energy_mid_eV
  REAL lambda_left_A,  energy_left_eV
  REAL lambda_right_A, energy_right_eV

  REAL subbin_flux_phcm2s1

! EUVAC solar flux wavelengths (lines or bin middles)

  lambda_euvac_A( 1) = 0.5 * ( 50.0 + 100.0)
  lambda_euvac_A( 2) = 0.5 * (100.0 + 150.0) 
  lambda_euvac_A( 3) = 0.5 * (150.0 + 200.0) 
  lambda_euvac_A( 4) = 0.5 * (200.0 + 250.0) 
  lambda_euvac_A( 5) = 256.32 
  lambda_euvac_A( 6) = 284.15
  lambda_euvac_A( 7) = 0.5 * (250.0 + 300.0)
  lambda_euvac_A( 8) = 303.31 
  lambda_euvac_A( 9) = 303.78 
  lambda_euvac_A(10) = 0.5 * (300.0 + 350.0) 

  lambda_euvac_A(11) = 368.07
  lambda_euvac_A(12) = 0.5 * (350.0 + 400.0) 
  lambda_euvac_A(13) = 0.5 * (400.0 + 450.0) 
  lambda_euvac_A(14) = 465.22 
  lambda_euvac_A(15) = 0.5 * (450.0 + 500.0) 
  lambda_euvac_A(16) = 0.5 * (500.0 + 550.0) 
  lambda_euvac_A(17) = 554.37 
  lambda_euvac_A(18) = 584.33
  lambda_euvac_A(19) = 0.5 * (550.0 + 600.0)
  lambda_euvac_A(20) = 609.76

  lambda_euvac_A(21) = 629.73
  lambda_euvac_A(22) = 0.5 * (600.0 + 650.0)
  lambda_euvac_A(23) = 0.5 * (650.0 + 700.0) 
  lambda_euvac_A(24) = 703.36 
  lambda_euvac_A(25) = 0.5 * (700.0 + 750.0) 
  lambda_euvac_A(26) = 765.15 
  lambda_euvac_A(27) = 770.41 
  lambda_euvac_A(28) = 789.36
  lambda_euvac_A(29) = 0.5 * (750.0 + 800.0)
  lambda_euvac_A(30) = 0.5 * (800.0 + 850.0) 

  lambda_euvac_A(31) = 0.5 * (850.0 + 900.0)
  lambda_euvac_A(32) = 0.5 * (900.0 + 950.0) 
  lambda_euvac_A(33) = 977.02
  lambda_euvac_A(34) = 0.5 * (950.0 + 1000.0)
  lambda_euvac_A(35) = 1025.72
  lambda_euvac_A(36) = 1031.91
  lambda_euvac_A(37) = 0.5 * (1000.0 + 1050.0)

! EUVAC solar flux energy bins

  factor_energy_eVA =  REAL(h_Planck_Js * c_ms * 1.0d10 / e_Cl)   ! 1.0d10 because the wavelength is in Angstrom

  DO i = 1, N_solar_bins
     solar_flux_energy_bin_eV(i) = factor_energy_eVA / lambda_euvac_A(i)
  END DO

! EUVAC solar flux [in 10^9 ph cm^-2 s^-1]
! correct the solar flux due to the solar activity variation
! use Eqs.(9.19) and (9.20) from the Schunk and Nagy's book

  factor_f107 = 0.5 * (f10p7 + f10p7_81) - 80.0

  solar_flux_phcm2s1( 1) = MAX( 0.0, 1.200 * (1.0 + 1.0017e-2 * factor_f107) )
  solar_flux_phcm2s1( 2) = MAX( 0.0, 0.450 * (1.0 + 7.1250e-3 * factor_f107) )
  solar_flux_phcm2s1( 3) = MAX( 0.0, 4.800 * (1.0 + 1.3375e-2 * factor_f107) )
  solar_flux_phcm2s1( 4) = MAX( 0.0, 3.100 * (1.0 + 1.9450e-2 * factor_f107) )
  solar_flux_phcm2s1( 5) = MAX( 0.0, 0.460 * (1.0 + 2.7750e-3 * factor_f107) )   !### line
  solar_flux_phcm2s1( 6) = MAX( 0.0, 0.210 * (1.0 + 1.3768e-1 * factor_f107) )   !### line
  solar_flux_phcm2s1( 7) = MAX( 0.0, 1.679 * (1.0 + 2.6467e-2 * factor_f107) )
  solar_flux_phcm2s1( 8) = MAX( 0.0, 0.800 * (1.0 + 2.5000e-2 * factor_f107) )   !### line
  solar_flux_phcm2s1( 9) = MAX( 0.0, 6.900 * (1.0 + 3.3333e-3 * factor_f107) )   !### line
  solar_flux_phcm2s1(10) = MAX( 0.0, 0.965 * (1.0 + 2.2450e-2 * factor_f107) )

  solar_flux_phcm2s1(11) = MAX( 0.0, 0.650 * (1.0 + 6.5917e-3 * factor_f107) )   !### line
  solar_flux_phcm2s1(12) = MAX( 0.0, 0.314 * (1.0 + 3.6542e-2 * factor_f107) )
  solar_flux_phcm2s1(13) = MAX( 0.0, 0.383 * (1.0 + 7.4083e-3 * factor_f107) )
  solar_flux_phcm2s1(14) = MAX( 0.0, 0.290 * (1.0 + 7.4917e-3 * factor_f107) )   !### line
  solar_flux_phcm2s1(15) = MAX( 0.0, 0.285 * (1.0 + 2.0225e-2 * factor_f107) )
  solar_flux_phcm2s1(16) = MAX( 0.0, 0.452 * (1.0 + 8.7583e-3 * factor_f107) )
  solar_flux_phcm2s1(17) = MAX( 0.0, 0.720 * (1.0 + 3.2667e-3 * factor_f107) )   !### line
  solar_flux_phcm2s1(18) = MAX( 0.0, 1.270 * (1.0 + 5.1583e-3 * factor_f107) )   !### line
  solar_flux_phcm2s1(19) = MAX( 0.0, 0.357 * (1.0 + 3.6583e-3 * factor_f107) )
  solar_flux_phcm2s1(20) = MAX( 0.0, 0.530 * (1.0 + 1.6175e-2 * factor_f107) )   !### line

  solar_flux_phcm2s1(21) = MAX( 0.0, 1.590 * (1.0 + 3.3250e-3 * factor_f107) )   !### line
  solar_flux_phcm2s1(22) = MAX( 0.0, 0.342 * (1.0 + 1.1800e-2 * factor_f107) )
  solar_flux_phcm2s1(23) = MAX( 0.0, 0.230 * (1.0 + 4.2667e-3 * factor_f107) )
  solar_flux_phcm2s1(24) = MAX( 0.0, 0.360 * (1.0 + 3.0417e-3 * factor_f107) )   !### line
  solar_flux_phcm2s1(25) = MAX( 0.0, 0.141 * (1.0 + 4.7500e-3 * factor_f107) )
  solar_flux_phcm2s1(26) = MAX( 0.0, 0.170 * (1.0 + 3.8500e-3 * factor_f107) )   !### line
  solar_flux_phcm2s1(27) = MAX( 0.0, 0.260 * (1.0 + 1.2808e-2 * factor_f107) )   !### line
  solar_flux_phcm2s1(28) = MAX( 0.0, 0.702 * (1.0 + 3.2750e-3 * factor_f107) )   !### line
  solar_flux_phcm2s1(29) = MAX( 0.0, 0.758 * (1.0 + 4.7667e-3 * factor_f107) )
  solar_flux_phcm2s1(30) = MAX( 0.0, 1.625 * (1.0 + 4.8167e-3 * factor_f107) )

  solar_flux_phcm2s1(31) = MAX( 0.0, 3.537 * (1.0 + 5.6750e-3 * factor_f107) )
  solar_flux_phcm2s1(32) = MAX( 0.0, 3.000 * (1.0 + 4.9833e-3 * factor_f107) )
  solar_flux_phcm2s1(33) = MAX( 0.0, 4.400 * (1.0 + 3.9417e-3 * factor_f107) )   !### line
  solar_flux_phcm2s1(34) = MAX( 0.0, 1.475 * (1.0 + 4.4167e-3 * factor_f107) )
  solar_flux_phcm2s1(35) = MAX( 0.0, 3.500 * (1.0 + 5.1833e-3 * factor_f107) )   !### line
  solar_flux_phcm2s1(36) = MAX( 0.0, 2.100 * (1.0 + 5.2833e-3 * factor_f107) )   !### line
  solar_flux_phcm2s1(37) = MAX( 0.0, 2.467 * (1.0 + 4.3750e-3 * factor_f107) )

! total absorption cross section for atomic oxygen [in 10^-18 cm^2]

  sigma_tot_O_cm2( 1) = 0.730
  sigma_tot_O_cm2( 2) = 1.839
  sigma_tot_O_cm2( 3) = 3.732
  sigma_tot_O_cm2( 4) = 5.202
  sigma_tot_O_cm2( 5) = 6.050
  sigma_tot_O_cm2( 6) = 7.080
  sigma_tot_O_cm2( 7) = 6.461
  sigma_tot_O_cm2( 8) = 7.680
  sigma_tot_O_cm2( 9) = 7.700
  sigma_tot_O_cm2(10) = 8.693

  sigma_tot_O_cm2(11) = 9.840
  sigma_tot_O_cm2(12) = 9.687
  sigma_tot_O_cm2(13) = 11.496
  sigma_tot_O_cm2(14) = 11.930
  sigma_tot_O_cm2(15) = 12.127
  sigma_tot_O_cm2(16) = 12.059
  sigma_tot_O_cm2(17) = 12.590
  sigma_tot_O_cm2(18) = 13.090
  sigma_tot_O_cm2(19) = 13.024
  sigma_tot_O_cm2(20) = 13.400

  sigma_tot_O_cm2(21) = 13.400
  sigma_tot_O_cm2(22) = 13.365
  sigma_tot_O_cm2(23) = 17.245
  sigma_tot_O_cm2(24) = 11.460
  sigma_tot_O_cm2(25) = 10.736
  sigma_tot_O_cm2(26) = 4.000
  sigma_tot_O_cm2(27) = 3.890
  sigma_tot_O_cm2(28) = 3.749
  sigma_tot_O_cm2(29) = 5.091
  sigma_tot_O_cm2(30) = 3.498

  sigma_tot_O_cm2(31) = 4.554
  sigma_tot_O_cm2(32) = 1.315
  sigma_tot_O_cm2(33) = 0.0
  sigma_tot_O_cm2(34) = 0.0  
  sigma_tot_O_cm2(35) = 0.0
  sigma_tot_O_cm2(36) = 0.0
  sigma_tot_O_cm2(37) = 0.0

! total absorption cross section for molecular oxygen [in 10^-18 cm^2]

  sigma_tot_O2_cm2( 1) = 1.316
  sigma_tot_O2_cm2( 2) = 3.806
  sigma_tot_O2_cm2( 3) = 7.509
  sigma_tot_O2_cm2( 4) = 10.900
  sigma_tot_O2_cm2( 5) = 13.370
  sigma_tot_O2_cm2( 6) = 15.790
  sigma_tot_O2_cm2( 7) = 14.387
  sigma_tot_O2_cm2( 8) = 16.800
  sigma_tot_O2_cm2( 9) = 16.810
  sigma_tot_O2_cm2(10) = 17.438

  sigma_tot_O2_cm2(11) = 18.320
  sigma_tot_O2_cm2(12) = 18.118
  sigma_tot_O2_cm2(13) = 20.310
  sigma_tot_O2_cm2(14) = 21.910
  sigma_tot_O2_cm2(15) = 23.101
  sigma_tot_O2_cm2(16) = 24.606
  sigma_tot_O2_cm2(17) = 26.040
  sigma_tot_O2_cm2(18) = 22.720
  sigma_tot_O2_cm2(19) = 26.610
  sigma_tot_O2_cm2(20) = 28.070

  sigma_tot_O2_cm2(21) = 32.060
  sigma_tot_O2_cm2(22) = 26.017
  sigma_tot_O2_cm2(23) = 21.919
  sigma_tot_O2_cm2(24) = 27.440
  sigma_tot_O2_cm2(25) = 28.535
  sigma_tot_O2_cm2(26) = 20.800
  sigma_tot_O2_cm2(27) = 18.910
  sigma_tot_O2_cm2(28) = 26.668
  sigma_tot_O2_cm2(29) = 22.145
  sigma_tot_O2_cm2(30) = 16.631

  sigma_tot_O2_cm2(31) = 8.562
  sigma_tot_O2_cm2(32) = 12.817
  sigma_tot_O2_cm2(33) = 18.730
  sigma_tot_O2_cm2(34) = 21.108
  sigma_tot_O2_cm2(35) = 1.630
  sigma_tot_O2_cm2(36) = 1.050
  sigma_tot_O2_cm2(37) = 1.346

! total absorption cross section for molecular nitrogen [in 10^-18 cm^2]

  sigma_tot_N2_cm2( 1) = 0.720
  sigma_tot_N2_cm2( 2) = 2.261
  sigma_tot_N2_cm2( 3) = 4.958
  sigma_tot_N2_cm2( 4) = 8.392
  sigma_tot_N2_cm2( 5) = 10.210
  sigma_tot_N2_cm2( 6) = 10.900
  sigma_tot_N2_cm2( 7) = 10.493
  sigma_tot_N2_cm2( 8) = 11.670
  sigma_tot_N2_cm2( 9) = 11.700
  sigma_tot_N2_cm2(10) = 13.857

  sigma_tot_N2_cm2(11) = 16.910
  sigma_tot_N2_cm2(12) = 16.395
  sigma_tot_N2_cm2(13) = 21.675
  sigma_tot_N2_cm2(14) = 23.160
  sigma_tot_N2_cm2(15) = 23.471
  sigma_tot_N2_cm2(16) = 24.501
  sigma_tot_N2_cm2(17) = 24.130
  sigma_tot_N2_cm2(18) = 22.400
  sigma_tot_N2_cm2(19) = 22.787
  sigma_tot_N2_cm2(20) = 22.790

  sigma_tot_N2_cm2(21) = 23.370
  sigma_tot_N2_cm2(22) = 23.339
  sigma_tot_N2_cm2(23) = 31.755
  sigma_tot_N2_cm2(24) = 26.540
  sigma_tot_N2_cm2(25) = 24.662
  sigma_tot_N2_cm2(26) = 120.490
  sigma_tot_N2_cm2(27) = 14.180
  sigma_tot_N2_cm2(28) = 16.487
  sigma_tot_N2_cm2(29) = 33.578
  sigma_tot_N2_cm2(30) = 16.992

  sigma_tot_N2_cm2(31) = 20.249
  sigma_tot_N2_cm2(32) = 9.680
  sigma_tot_N2_cm2(33) = 2.240
  sigma_tot_N2_cm2(34) = 50.988
  sigma_tot_N2_cm2(35) = 0.0
  sigma_tot_N2_cm2(36) = 0.0
  sigma_tot_N2_cm2(37) = 0.0

! total absorption cross section for Helium (same as the ionization) [in 10^-18 cm^2]

  sigma_ion_He_cm2( 1) = 0.1441
  sigma_ion_He_cm2( 2) = 0.4785
  sigma_ion_He_cm2( 3) = 1.1571
  sigma_ion_He_cm2( 4) = 1.6008
  sigma_ion_He_cm2( 5) = 2.1212
  sigma_ion_He_cm2( 6) = 2.5947
  sigma_ion_He_cm2( 7) = 2.3205
  sigma_ion_He_cm2( 8) = 2.9529
  sigma_ion_He_cm2( 9) = 2.9618
  sigma_ion_He_cm2(10) = 3.5437

  sigma_ion_He_cm2(11) = 4.2675
  sigma_ion_He_cm2(12) = 4.1424
  sigma_ion_He_cm2(13) = 5.4466
  sigma_ion_He_cm2(14) = 6.5631
  sigma_ion_He_cm2(15) = 7.2084
  sigma_ion_He_cm2(16) = 0.9581
  sigma_ion_He_cm2(17) = 0.0
  sigma_ion_He_cm2(18) = 0.0
  sigma_ion_He_cm2(19) = 0.0
  sigma_ion_He_cm2(20) = 0.0

  sigma_ion_He_cm2(21) = 0.0
  sigma_ion_He_cm2(22) = 0.0
  sigma_ion_He_cm2(23) = 0.0
  sigma_ion_He_cm2(24) = 0.0
  sigma_ion_He_cm2(25) = 0.0
  sigma_ion_He_cm2(26) = 0.0
  sigma_ion_He_cm2(27) = 0.0
  sigma_ion_He_cm2(28) = 0.0
  sigma_ion_He_cm2(29) = 0.0
  sigma_ion_He_cm2(30) = 0.0

  sigma_ion_He_cm2(31) = 0.0
  sigma_ion_He_cm2(32) = 0.0
  sigma_ion_He_cm2(33) = 0.0
  sigma_ion_He_cm2(34) = 0.0
  sigma_ion_He_cm2(35) = 0.0
  sigma_ion_He_cm2(36) = 0.0
  sigma_ion_He_cm2(37) = 0.0

! ionization cross section for molecular nitrogen (N2 -> N2+) [in 10^-18 cm^2]

  sigma_ion_N2_to_N2_cm2( 1) = 0.443
  sigma_ion_N2_to_N2_cm2( 2) = 1.479
  sigma_ion_N2_to_N2_cm2( 3) = 3.153
  sigma_ion_N2_to_N2_cm2( 4) = 5.226
  sigma_ion_N2_to_N2_cm2( 5) = 6.781
  sigma_ion_N2_to_N2_cm2( 6) = 8.100
  sigma_ion_N2_to_N2_cm2( 7) = 7.347
  sigma_ion_N2_to_N2_cm2( 8) = 9.180
  sigma_ion_N2_to_N2_cm2( 9) = 9.210
  sigma_ion_N2_to_N2_cm2(10) = 11.600

  sigma_ion_N2_to_N2_cm2(11) = 15.350
  sigma_ion_N2_to_N2_cm2(12) = 14.669
  sigma_ion_N2_to_N2_cm2(13) = 20.692
  sigma_ion_N2_to_N2_cm2(14) = 22.100
  sigma_ion_N2_to_N2_cm2(15) = 22.772
  sigma_ion_N2_to_N2_cm2(16) = 24.468
  sigma_ion_N2_to_N2_cm2(17) = 24.130
  sigma_ion_N2_to_N2_cm2(18) = 22.400
  sigma_ion_N2_to_N2_cm2(19) = 22.787
  sigma_ion_N2_to_N2_cm2(20) = 22.790

  sigma_ion_N2_to_N2_cm2(21) = 23.370
  sigma_ion_N2_to_N2_cm2(22) = 23.339
  sigma_ion_N2_to_N2_cm2(23) = 29.235
  sigma_ion_N2_to_N2_cm2(24) = 25.480
  sigma_ion_N2_to_N2_cm2(25) = 15.060
  sigma_ion_N2_to_N2_cm2(26) = 65.800
  sigma_ion_N2_to_N2_cm2(27) = 8.500
  sigma_ion_N2_to_N2_cm2(28) = 8.860
  sigma_ion_N2_to_N2_cm2(29) = 14.274
  sigma_ion_N2_to_N2_cm2(30) = 0.0

  sigma_ion_N2_to_N2_cm2(31) = 0.0
  sigma_ion_N2_to_N2_cm2(32) = 0.0
  sigma_ion_N2_to_N2_cm2(33) = 0.0
  sigma_ion_N2_to_N2_cm2(34) = 0.0
  sigma_ion_N2_to_N2_cm2(35) = 0.0
  sigma_ion_N2_to_N2_cm2(36) = 0.0
  sigma_ion_N2_to_N2_cm2(37) = 0.0

! ionization cross section for molecular nitrogen (N2 -> N+) [in 10^-18 cm^2]

  sigma_ion_N2_to_N_cm2( 1) = 0.277
  sigma_ion_N2_to_N_cm2( 2) = 0.782
  sigma_ion_N2_to_N_cm2( 3) = 1.805
  sigma_ion_N2_to_N_cm2( 4) = 3.166
  sigma_ion_N2_to_N_cm2( 5) = 3.420
  sigma_ion_N2_to_N_cm2( 6) = 2.800
  sigma_ion_N2_to_N_cm2( 7) = 3.145
  sigma_ion_N2_to_N_cm2( 8) = 2.490
  sigma_ion_N2_to_N_cm2( 9) = 2.490
  sigma_ion_N2_to_N_cm2(10) = 2.257

  sigma_ion_N2_to_N_cm2(11) = 1.560
  sigma_ion_N2_to_N_cm2(12) = 1.726
  sigma_ion_N2_to_N_cm2(13) = 0.982
  sigma_ion_N2_to_N_cm2(14) = 1.060
  sigma_ion_N2_to_N_cm2(15) = 0.699  
  sigma_ion_N2_to_N_cm2(16) = 0.033
  sigma_ion_N2_to_N_cm2(17) = 0.0
  sigma_ion_N2_to_N_cm2(18) = 0.0
  sigma_ion_N2_to_N_cm2(19) = 0.0
  sigma_ion_N2_to_N_cm2(20) = 0.0

  sigma_ion_N2_to_N_cm2(21) = 0.0
  sigma_ion_N2_to_N_cm2(22) = 0.0
  sigma_ion_N2_to_N_cm2(23) = 0.0
  sigma_ion_N2_to_N_cm2(24) = 0.0
  sigma_ion_N2_to_N_cm2(25) = 0.0
  sigma_ion_N2_to_N_cm2(26) = 0.0
  sigma_ion_N2_to_N_cm2(27) = 0.0
  sigma_ion_N2_to_N_cm2(28) = 0.0
  sigma_ion_N2_to_N_cm2(29) = 0.0
  sigma_ion_N2_to_N_cm2(30) = 0.0

  sigma_ion_N2_to_N_cm2(31) = 0.0
  sigma_ion_N2_to_N_cm2(32) = 0.0
  sigma_ion_N2_to_N_cm2(33) = 0.0
  sigma_ion_N2_to_N_cm2(34) = 0.0
  sigma_ion_N2_to_N_cm2(35) = 0.0
  sigma_ion_N2_to_N_cm2(36) = 0.0
  sigma_ion_N2_to_N_cm2(37) = 0.0

! ionization cross section for atomic oxygen (O -> O+4S) [in 10^-18 cm^2]

  sigma_ion_O_4S_cm2( 1) = 0.190
  sigma_ion_O_4S_cm2( 2) = 0.486
  sigma_ion_O_4S_cm2( 3) = 0.952
  sigma_ion_O_4S_cm2( 4) = 1.311
  sigma_ion_O_4S_cm2( 5) = 1.539
  sigma_ion_O_4S_cm2( 6) = 1.770
  sigma_ion_O_4S_cm2( 7) = 1.628
  sigma_ion_O_4S_cm2( 8) = 1.920
  sigma_ion_O_4S_cm2( 9) = 1.925
  sigma_ion_O_4S_cm2(10) = 2.259

  sigma_ion_O_4S_cm2(11) = 2.559
  sigma_ion_O_4S_cm2(12) = 2.523
  sigma_ion_O_4S_cm2(13) = 3.073
  sigma_ion_O_4S_cm2(14) = 3.340
  sigma_ion_O_4S_cm2(15) = 3.394
  sigma_ion_O_4S_cm2(16) = 3.421
  sigma_ion_O_4S_cm2(17) = 3.650
  sigma_ion_O_4S_cm2(18) = 3.920
  sigma_ion_O_4S_cm2(19) = 3.620
  sigma_ion_O_4S_cm2(20) = 3.610

  sigma_ion_O_4S_cm2(21) = 3.880
  sigma_ion_O_4S_cm2(22) = 4.250
  sigma_ion_O_4S_cm2(23) = 5.128
  sigma_ion_O_4S_cm2(24) = 4.890
  sigma_ion_O_4S_cm2(25) = 6.739
  sigma_ion_O_4S_cm2(26) = 4.000
  sigma_ion_O_4S_cm2(27) = 3.890
  sigma_ion_O_4S_cm2(28) = 3.749
  sigma_ion_O_4S_cm2(29) = 5.091
  sigma_ion_O_4S_cm2(30) = 3.498

  sigma_ion_O_4S_cm2(31) = 4.554
  sigma_ion_O_4S_cm2(32) = 1.315
  sigma_ion_O_4S_cm2(33) = 0.0
  sigma_ion_O_4S_cm2(34) = 0.0
  sigma_ion_O_4S_cm2(35) = 0.0
  sigma_ion_O_4S_cm2(36) = 0.0
  sigma_ion_O_4S_cm2(37) = 0.0

! ionization cross section for atomic oxygen (O -> O+2D) [in 10^-18 cm^2]

  sigma_ion_O_2D_cm2( 1) = 0.206
  sigma_ion_O_2D_cm2( 2) = 0.529
  sigma_ion_O_2D_cm2( 3) = 1.171
  sigma_ion_O_2D_cm2( 4) = 1.762
  sigma_ion_O_2D_cm2( 5) = 2.138
  sigma_ion_O_2D_cm2( 6) = 2.620
  sigma_ion_O_2D_cm2( 7) = 2.325
  sigma_ion_O_2D_cm2( 8) = 2.842
  sigma_ion_O_2D_cm2( 9) = 2.849
  sigma_ion_O_2D_cm2(10) = 3.446

  sigma_ion_O_2D_cm2(11) = 3.936
  sigma_ion_O_2D_cm2(12) = 3.883
  sigma_ion_O_2D_cm2(13) = 4.896
  sigma_ion_O_2D_cm2(14) = 5.370
  sigma_ion_O_2D_cm2(15) = 5.459
  sigma_ion_O_2D_cm2(16) = 5.427
  sigma_ion_O_2D_cm2(17) = 5.670
  sigma_ion_O_2D_cm2(18) = 6.020
  sigma_ion_O_2D_cm2(19) = 5.910
  sigma_ion_O_2D_cm2(20) = 6.170

  sigma_ion_O_2D_cm2(21) = 6.290
  sigma_ion_O_2D_cm2(22) = 6.159
  sigma_ion_O_2D_cm2(23) = 11.453
  sigma_ion_O_2D_cm2(24) = 6.570
  sigma_ion_O_2D_cm2(25) = 3.997
  sigma_ion_O_2D_cm2(26) = 0.0
  sigma_ion_O_2D_cm2(27) = 0.0
  sigma_ion_O_2D_cm2(28) = 0.0
  sigma_ion_O_2D_cm2(29) = 0.0
  sigma_ion_O_2D_cm2(30) = 0.0

  sigma_ion_O_2D_cm2(31) = 0.0
  sigma_ion_O_2D_cm2(32) = 0.0
  sigma_ion_O_2D_cm2(33) = 0.0
  sigma_ion_O_2D_cm2(34) = 0.0
  sigma_ion_O_2D_cm2(35) = 0.0
  sigma_ion_O_2D_cm2(36) = 0.0
  sigma_ion_O_2D_cm2(37) = 0.0

! ionization cross section for atomic oxygen (O -> O+2P) [in 10^-18 cm^2]

  sigma_ion_O_2P_cm2( 1) = 0.134
  sigma_ion_O_2P_cm2( 2) = 0.345
  sigma_ion_O_2P_cm2( 3) = 0.768
  sigma_ion_O_2P_cm2( 4) = 1.144
  sigma_ion_O_2P_cm2( 5) = 1.363
  sigma_ion_O_2P_cm2( 6) = 1.630
  sigma_ion_O_2P_cm2( 7) = 1.488
  sigma_ion_O_2P_cm2( 8) = 1.920
  sigma_ion_O_2P_cm2( 9) = 1.925
  sigma_ion_O_2P_cm2(10) = 2.173

  sigma_ion_O_2P_cm2(11) = 2.558
  sigma_ion_O_2P_cm2(12) = 2.422
  sigma_ion_O_2P_cm2(13) = 2.986
  sigma_ion_O_2P_cm2(14) = 3.220
  sigma_ion_O_2P_cm2(15) = 3.274
  sigma_ion_O_2P_cm2(16) = 3.211
  sigma_ion_O_2P_cm2(17) = 3.270
  sigma_ion_O_2P_cm2(18) = 3.150
  sigma_ion_O_2P_cm2(19) = 3.494
  sigma_ion_O_2P_cm2(20) = 3.620

  sigma_ion_O_2P_cm2(21) = 3.230
  sigma_ion_O_2P_cm2(22) = 2.956
  sigma_ion_O_2P_cm2(23) = 0.664
  sigma_ion_O_2P_cm2(24) = 0.0
  sigma_ion_O_2P_cm2(25) = 0.0
  sigma_ion_O_2P_cm2(26) = 0.0
  sigma_ion_O_2P_cm2(27) = 0.0
  sigma_ion_O_2P_cm2(28) = 0.0
  sigma_ion_O_2P_cm2(29) = 0.0
  sigma_ion_O_2P_cm2(30) = 0.0

  sigma_ion_O_2P_cm2(31) = 0.0
  sigma_ion_O_2P_cm2(32) = 0.0
  sigma_ion_O_2P_cm2(33) = 0.0
  sigma_ion_O_2P_cm2(34) = 0.0
  sigma_ion_O_2P_cm2(35) = 0.0
  sigma_ion_O_2P_cm2(36) = 0.0
  sigma_ion_O_2P_cm2(37) = 0.0

! ionization cross section for atomic oxygen (O -> O+4P) [in 10^-18 cm^2]

  sigma_ion_O_4P_cm2( 1) = 0.062
  sigma_ion_O_4P_cm2( 2) = 0.163
  sigma_ion_O_4P_cm2( 3) = 0.348
  sigma_ion_O_4P_cm2( 4) = 0.508
  sigma_ion_O_4P_cm2( 5) = 0.598
  sigma_ion_O_4P_cm2( 6) = 0.710
  sigma_ion_O_4P_cm2( 7) = 0.637
  sigma_ion_O_4P_cm2( 8) = 0.691
  sigma_ion_O_4P_cm2( 9) = 0.693
  sigma_ion_O_4P_cm2(10) = 0.815

  sigma_ion_O_4P_cm2(11) = 0.787
  sigma_ion_O_4P_cm2(12) = 0.859
  sigma_ion_O_4P_cm2(13) = 0.541
  sigma_ion_O_4P_cm2(14) = 0.0
  sigma_ion_O_4P_cm2(15) = 0.0
  sigma_ion_O_4P_cm2(16) = 0.0
  sigma_ion_O_4P_cm2(17) = 0.0
  sigma_ion_O_4P_cm2(18) = 0.0
  sigma_ion_O_4P_cm2(19) = 0.0
  sigma_ion_O_4P_cm2(20) = 0.0

  sigma_ion_O_4P_cm2(21) = 0.0
  sigma_ion_O_4P_cm2(22) = 0.0
  sigma_ion_O_4P_cm2(23) = 0.0
  sigma_ion_O_4P_cm2(24) = 0.0
  sigma_ion_O_4P_cm2(25) = 0.0
  sigma_ion_O_4P_cm2(26) = 0.0
  sigma_ion_O_4P_cm2(27) = 0.0
  sigma_ion_O_4P_cm2(28) = 0.0
  sigma_ion_O_4P_cm2(29) = 0.0
  sigma_ion_O_4P_cm2(30) = 0.0

  sigma_ion_O_4P_cm2(31) = 0.0
  sigma_ion_O_4P_cm2(32) = 0.0
  sigma_ion_O_4P_cm2(33) = 0.0
  sigma_ion_O_4P_cm2(34) = 0.0
  sigma_ion_O_4P_cm2(35) = 0.0
  sigma_ion_O_4P_cm2(36) = 0.0
  sigma_ion_O_4P_cm2(37) = 0.0

! ionization cross section for atomic oxygen (O -> O+2P*) [in 10^-18 cm^2]

  sigma_ion_O_2Pst_cm2( 1) = 0.049
  sigma_ion_O_2Pst_cm2( 2) = 0.130
  sigma_ion_O_2Pst_cm2( 3) = 0.278
  sigma_ion_O_2Pst_cm2( 4) = 0.366
  sigma_ion_O_2Pst_cm2( 5) = 0.412
  sigma_ion_O_2Pst_cm2( 6) = 0.350
  sigma_ion_O_2Pst_cm2( 7) = 0.383
  sigma_ion_O_2Pst_cm2( 8) = 0.307
  sigma_ion_O_2Pst_cm2( 9) = 0.308
  sigma_ion_O_2Pst_cm2(10) = 0.0

  sigma_ion_O_2Pst_cm2(11) = 0.0
  sigma_ion_O_2Pst_cm2(12) = 0.0
  sigma_ion_O_2Pst_cm2(13) = 0.0
  sigma_ion_O_2Pst_cm2(14) = 0.0
  sigma_ion_O_2Pst_cm2(15) = 0.0
  sigma_ion_O_2Pst_cm2(16) = 0.0
  sigma_ion_O_2Pst_cm2(17) = 0.0
  sigma_ion_O_2Pst_cm2(18) = 0.0
  sigma_ion_O_2Pst_cm2(19) = 0.0
  sigma_ion_O_2Pst_cm2(20) = 0.0

  sigma_ion_O_2Pst_cm2(21) = 0.0
  sigma_ion_O_2Pst_cm2(22) = 0.0
  sigma_ion_O_2Pst_cm2(23) = 0.0
  sigma_ion_O_2Pst_cm2(24) = 0.0
  sigma_ion_O_2Pst_cm2(25) = 0.0
  sigma_ion_O_2Pst_cm2(26) = 0.0
  sigma_ion_O_2Pst_cm2(27) = 0.0
  sigma_ion_O_2Pst_cm2(28) = 0.0
  sigma_ion_O_2Pst_cm2(29) = 0.0
  sigma_ion_O_2Pst_cm2(30) = 0.0

  sigma_ion_O_2Pst_cm2(31) = 0.0
  sigma_ion_O_2Pst_cm2(32) = 0.0
  sigma_ion_O_2Pst_cm2(33) = 0.0
  sigma_ion_O_2Pst_cm2(34) = 0.0
  sigma_ion_O_2Pst_cm2(35) = 0.0
  sigma_ion_O_2Pst_cm2(36) = 0.0
  sigma_ion_O_2Pst_cm2(37) = 0.0

! ionization cross section for molecular oxygen (O2 -> O2+) [in 10^-18 cm^2]

  sigma_ion_O2_to_O2_cm2( 1) = 1.316
  sigma_ion_O2_to_O2_cm2( 2) = 2.346
  sigma_ion_O2_to_O2_cm2( 3) = 4.139
  sigma_ion_O2_to_O2_cm2( 4) = 6.619
  sigma_ion_O2_to_O2_cm2( 5) = 8.460
  sigma_ion_O2_to_O2_cm2( 6) = 9.890
  sigma_ion_O2_to_O2_cm2( 7) = 9.056
  sigma_ion_O2_to_O2_cm2( 8) = 10.860
  sigma_ion_O2_to_O2_cm2( 9) = 10.880
  sigma_ion_O2_to_O2_cm2(10) = 12.229

  sigma_ion_O2_to_O2_cm2(11) = 13.760
  sigma_ion_O2_to_O2_cm2(12) = 13.418
  sigma_ion_O2_to_O2_cm2(13) = 15.490
  sigma_ion_O2_to_O2_cm2(14) = 16.970
  sigma_ion_O2_to_O2_cm2(15) = 17.754
  sigma_ion_O2_to_O2_cm2(16) = 19.469
  sigma_ion_O2_to_O2_cm2(17) = 21.600
  sigma_ion_O2_to_O2_cm2(18) = 18.840
  sigma_ion_O2_to_O2_cm2(19) = 22.789
  sigma_ion_O2_to_O2_cm2(20) = 24.540

  sigma_ion_O2_to_O2_cm2(21) = 30.070
  sigma_ion_O2_to_O2_cm2(22) = 23.974
  sigma_ion_O2_to_O2_cm2(23) = 21.116
  sigma_ion_O2_to_O2_cm2(24) = 23.750
  sigma_ion_O2_to_O2_cm2(25) = 23.805
  sigma_ion_O2_to_O2_cm2(26) = 11.720
  sigma_ion_O2_to_O2_cm2(27) = 8.470
  sigma_ion_O2_to_O2_cm2(28) = 10.191
  sigma_ion_O2_to_O2_cm2(29) = 10.597
  sigma_ion_O2_to_O2_cm2(30) = 6.413

  sigma_ion_O2_to_O2_cm2(31) = 5.494
  sigma_ion_O2_to_O2_cm2(32) = 9.374
  sigma_ion_O2_to_O2_cm2(33) = 15.540
  sigma_ion_O2_to_O2_cm2(34) = 13.940
  sigma_ion_O2_to_O2_cm2(35) = 1.050
  sigma_ion_O2_to_O2_cm2(36) = 0.000
  sigma_ion_O2_to_O2_cm2(37) = 0.259

! ionization cross section for molecular oxygen (O2 -> O+) [in 10^-18 cm^2]

  sigma_ion_O2_to_O_cm2( 1) = 0.0
  sigma_ion_O2_to_O_cm2( 2) = 1.460
  sigma_ion_O2_to_O_cm2( 3) = 3.368
  sigma_ion_O2_to_O_cm2( 4) = 4.281
  sigma_ion_O2_to_O_cm2( 5) = 4.910
  sigma_ion_O2_to_O_cm2( 6) = 5.900
  sigma_ion_O2_to_O_cm2( 7) = 5.332
  sigma_ion_O2_to_O_cm2( 8) = 5.940
  sigma_ion_O2_to_O_cm2( 9) = 5.930
  sigma_ion_O2_to_O_cm2(10) = 5.212

  sigma_ion_O2_to_O_cm2(11) = 4.560
  sigma_ion_O2_to_O_cm2(12) = 4.703
  sigma_ion_O2_to_O_cm2(13) = 4.818
  sigma_ion_O2_to_O_cm2(14) = 4.940
  sigma_ion_O2_to_O_cm2(15) = 5.347
  sigma_ion_O2_to_O_cm2(16) = 5.139
  sigma_ion_O2_to_O_cm2(17) = 4.440
  sigma_ion_O2_to_O_cm2(18) = 3.880
  sigma_ion_O2_to_O_cm2(19) = 3.824
  sigma_ion_O2_to_O_cm2(20) = 1.850

  sigma_ion_O2_to_O_cm2(21) = 1.030
  sigma_ion_O2_to_O_cm2(22) = 0.962
  sigma_ion_O2_to_O_cm2(23) = 0.190
  sigma_ion_O2_to_O_cm2(24) = 0.0
  sigma_ion_O2_to_O_cm2(25) = 0.0
  sigma_ion_O2_to_O_cm2(26) = 0.0
  sigma_ion_O2_to_O_cm2(27) = 0.0
  sigma_ion_O2_to_O_cm2(28) = 0.0
  sigma_ion_O2_to_O_cm2(29) = 0.0
  sigma_ion_O2_to_O_cm2(30) = 0.0

  sigma_ion_O2_to_O_cm2(31) = 0.0
  sigma_ion_O2_to_O_cm2(32) = 0.0
  sigma_ion_O2_to_O_cm2(33) = 0.0
  sigma_ion_O2_to_O_cm2(34) = 0.0
  sigma_ion_O2_to_O_cm2(35) = 0.0
  sigma_ion_O2_to_O_cm2(36) = 0.0
  sigma_ion_O2_to_O_cm2(37) = 0.0

! ionization cross section for Helium (He -> He+) [in 10^-18 cm^2]

  sigma_ion_He_cm2( 1) = 0.1441
  sigma_ion_He_cm2( 2) = 0.4785
  sigma_ion_He_cm2( 3) = 1.1571
  sigma_ion_He_cm2( 4) = 1.6008
  sigma_ion_He_cm2( 5) = 2.1212
  sigma_ion_He_cm2( 6) = 2.5947
  sigma_ion_He_cm2( 7) = 2.3205
  sigma_ion_He_cm2( 8) = 2.9529
  sigma_ion_He_cm2( 9) = 2.9618
  sigma_ion_He_cm2(10) = 3.5437

  sigma_ion_He_cm2(11) = 4.2675
  sigma_ion_He_cm2(12) = 4.1424
  sigma_ion_He_cm2(13) = 5.4466
  sigma_ion_He_cm2(14) = 6.5631
  sigma_ion_He_cm2(15) = 7.2084
  sigma_ion_He_cm2(16) = 0.9581
  sigma_ion_He_cm2(17) = 0.0
  sigma_ion_He_cm2(18) = 0.0
  sigma_ion_He_cm2(19) = 0.0
  sigma_ion_He_cm2(20) = 0.0

  sigma_ion_He_cm2(21) = 0.0
  sigma_ion_He_cm2(22) = 0.0
  sigma_ion_He_cm2(23) = 0.0
  sigma_ion_He_cm2(24) = 0.0
  sigma_ion_He_cm2(25) = 0.0
  sigma_ion_He_cm2(26) = 0.0
  sigma_ion_He_cm2(27) = 0.0
  sigma_ion_He_cm2(28) = 0.0
  sigma_ion_He_cm2(29) = 0.0
  sigma_ion_He_cm2(30) = 0.0

  sigma_ion_He_cm2(31) = 0.0
  sigma_ion_He_cm2(32) = 0.0
  sigma_ion_He_cm2(33) = 0.0
  sigma_ion_He_cm2(34) = 0.0
  sigma_ion_He_cm2(35) = 0.0
  sigma_ion_He_cm2(36) = 0.0
  sigma_ion_He_cm2(37) = 0.0

  IF (Rank_of_process.NE.0) RETURN

! the spectrum is saved by the zero rank process only

!###########

  OPEN (10, FILE = 'solar_EUV_spectrum_EUVAC_based_cont.dat')

  WRITE (10, '("# col  1 is the ordering number of the continuous EUVAC bin")')
  WRITE (10, '("# col  2 is the rank of MPI process using the wavelength")')
  WRITE (10, '("# col  3 is the ordering number of primary photoelectron particle emitted by the MPI process")')
  WRITE (10, '("# col  4 is the wavelength of the middle of spectrum sub-bin (A)")')
  WRITE (10, '("# col  5 is the energy of photons with the wavelength of the middle of spectrum sub-bin (eV)")')
  WRITE (10, '("# col  6 is the photon flux of spectrum sub-bin (10^9 photons/cm^2/s)")')
  WRITE (10, '("# col  7 is the average differential (over wavelength) photon flux of spectrum sub-bin (10^9 photons/cm^2/s/A)")')
  WRITE (10, '("# col  8 is the average differential (over energy)     photon flux of spectrum sub-bin (10^9 photons/cm^2/s/eV)")')
  WRITE (10, '("# col  9 is the photon energy flux of spectrum sub-bin (10^9 eV/cm^2/s)")')
  WRITE (10, '("# col 10 is the average differential (over wavelength) photon energy flux of spectrum sub-bin (10^9 eV/cm^2/s/A)")')
  WRITE (10, '("# col 11 is the average differential (over energy) photon energy flux of spectrum sub-bin (10^9 eV/cm^2/s/eV)")')

  subbin_width_A = 50.0 / REAL(N_cascades_per_energy_bin * N_of_processes)

  DO i = 1, N_solar_bins

! skip the lines
     SELECT CASE (i)
        CASE(5,6,8,9,11,14,17,18,20,21,24,26,27,28,33,35,36)
           CYCLE
     END SELECT

     DO rank = 0, N_of_processes-1
        
        DO n = 1, N_cascades_per_energy_bin

           lambda_mid_A = lambda_euvac_A(i) - 25.0 + &                                                !
                        & 50.0 * (REAL(rank) / REAL(N_of_processes)) + &                              ! middle wavelength of spectrum sub-bin (A)
                        & 50.0 * (REAL(n) - 0.5) / REAL(N_cascades_per_energy_bin * N_of_processes)   !

           energy_mid_eV = factor_energy_eVA / lambda_mid_A

           lambda_left_A  = lambda_mid_A - 25.0 / REAL(N_cascades_per_energy_bin * N_of_processes)
           lambda_right_A = lambda_mid_A + 25.0 / REAL(N_cascades_per_energy_bin * N_of_processes)

           energy_left_eV  = factor_energy_eVA / lambda_left_A
           energy_right_eV = factor_energy_eVA / lambda_right_A

           subbin_flux_phcm2s1 = (lambda_mid_A / lambda_euvac_A(i)) * solar_flux_phcm2s1(i) / REAL(N_cascades_per_energy_bin * N_of_processes) ! in units of 10^9 photons/cm^2/s

           IF (i.EQ.3) THEN
              IF (lambda_mid_A.GE.lambda_euvac_A(i)) THEN
                 subbin_flux_phcm2s1 = 1.8 * subbin_flux_phcm2s1   ! 1.5
              ELSE
                 subbin_flux_phcm2s1 = 0.2 * subbin_flux_phcm2s1   ! 0.5
              END IF
           END IF

           WRITE (10, '(3(2x,i4),2x,f7.2,2x,f8.3,6(2x,e12.5))') &
                & i, &           ! 1
                & rank, &            ! 2
                & n, &           ! 3
                & lambda_mid_A, &               ! 4
                & energy_mid_eV, &              ! 5
                & subbin_flux_phcm2s1, &                                            ! 6
                & subbin_flux_phcm2s1 / subbin_width_A, &                           ! 7
                & subbin_flux_phcm2s1 / (energy_right_eV - energy_left_eV), &       ! 8
                & subbin_flux_phcm2s1 * energy_mid_eV, &                                         ! 9
                & subbin_flux_phcm2s1 * energy_mid_eV / subbin_width_A, &                        ! 10
                & subbin_flux_phcm2s1 * energy_mid_eV / (energy_right_eV - energy_left_eV)       ! 11

        END DO
     END DO
  END DO
  CLOSE (10, STATUS = 'KEEP')
  PRINT '("### file  solar_EUV_spectrum_EUVAC_based_cont.dat is ready")'

!###########

  OPEN (10, FILE = 'solar_EUV_spectrum_EUVAC_lines.dat') 

  WRITE (10, '("# col  1 is the ordering number of the discrete EUVAC line")')
  WRITE (10, '("# col  2 is the wavelength of the line (A)")')
  WRITE (10, '("# col  3 is the energy of photons with the wavelength of the line (eV)")')
  WRITE (10, '("# col  4 is the photon flux of the line (10^9 photons/cm^2/s)")')

  DO i = 1, N_solar_bins

! skip the continuous bins
     SELECT CASE (i)
        CASE(1:4,7,10,12,13,15,16,19,22,23,25,29:32,34,37)
           CYCLE
     END SELECT

     WRITE (10, '(2x,i2,2x,f7.2,2x,f8.3,2x,e12.5)') i, lambda_euvac_A(i), solar_flux_energy_bin_eV(i), 1.0d-10 * solar_flux_phcm2s1(i)
     WRITE (10, '(2x,i2,2x,f7.2,2x,f8.3,2x,e12.5)') i, lambda_euvac_A(i), solar_flux_energy_bin_eV(i),           solar_flux_phcm2s1(i)
     WRITE (10, '(" ")')
     WRITE (10, '(" ")')

  END DO

  CLOSE (10, STATUS = 'KEEP')
  PRINT '("### file  solar_EUV_spectrum_EUVAC_lines.dat is ready")'

!###########

  OPEN (10, FILE = 'solar_EUV_spectrum_EUVAC_cont.dat')
  
  WRITE (10, '("# col  1 is the ordering number of the continuous EUVAC bin")')
  WRITE (10, '("# col  2 is the wavelength of the middle of the bin (A)")')
  WRITE (10, '("# col  3 is the energy of photons with the middle wavelength of the bin (eV)")')
  WRITE (10, '("# col  4 is the photon flux of the bin (10^9 photons/cm^2/s)")')

  DO i = 1, N_solar_bins
! skip the lines
     SELECT CASE (i)
        CASE(5,6,8,9,11,14,17,18,20,21,24,26,27,28,33,35,36)
           CYCLE
     END SELECT

     WRITE (10, '(2x,i2,2x,f7.2,2x,f8.3,2x,e12.5)') &
          & i, &
          & lambda_euvac_A(i), &
          & solar_flux_energy_bin_eV(i), &
          & solar_flux_phcm2s1(i)
  END DO
  CLOSE (10, STATUS = 'KEEP')
  PRINT '("### file  solar_EUV_spectrum_EUVAC_cont.dat is ready")'

!###########
! this file represents the continuous part of the EUVAC spectrum as a continuous spectrum
! with the differential flux equal to the average differential EUVAC flux within the corresponding continuous bin

  OPEN (10, FILE = 'solar_EUV_spectrum_EUVAC_cont_box.dat')

  WRITE (10, '("# col 1 is the wavelength (A)")')
  WRITE (10, '("# col 2 is the energy of photons corresponding to the wavelength (eV)")')
  WRITE (10, '("# col 3 is the average differential (over wavelength) photon flux of spectrum sub-bin (10^9 photons/cm^2/s/A)")')
  WRITE (10, '("# col 4 is the average differential (over energy)     photon flux of spectrum sub-bin (10^9 photons/cm^2/s/eV)")')
  WRITE (10, '("# col 5 is the average differential (over wavelength) photon energy flux of spectrum sub-bin (10^9 eV/cm^2/s/A)")')
  WRITE (10, '("# col 6 is the average differential (over energy) photon energy flux of spectrum sub-bin (10^9 eV/cm^2/s/eV)")')

  DO i = 1, N_solar_bins
! skip the lines
     SELECT CASE (i)
        CASE(5,6,8,9,11,14,17,18,20,21,24,26,27,28,33,35,36)
           CYCLE
     END SELECT

     lambda_left_A  = lambda_euvac_A(i) - 25.0
     lambda_right_A = lambda_euvac_A(i) + 25.0

     energy_left_eV  = factor_energy_eVA / lambda_left_A
     energy_mid_eV   = factor_energy_eVA / lambda_euvac_A(i)
     energy_right_eV = factor_energy_eVA / lambda_right_A

     WRITE (10, '(2x,f7.2,2x,f8.3,4(2x,e12.5))') &
          & lambda_left_A, &
          & energy_left_eV, &
          & 1.0d-10 * solar_flux_phcm2s1(i) / 50.0, &
          & 1.0d-10 * solar_flux_phcm2s1(i) / ABS(energy_right_eV - energy_left_eV), &
          & 1.0d-10 * solar_flux_phcm2s1(i) * solar_flux_energy_bin_eV(i) / 50.0, &
          & 1.0d-10 * solar_flux_phcm2s1(i) * solar_flux_energy_bin_eV(i) / ABS(energy_right_eV - energy_left_eV)

     WRITE (10, '(2x,f7.2,2x,f8.3,4(2x,e12.5))') &
          & lambda_left_A, &
          & energy_left_eV, &
          & solar_flux_phcm2s1(i) / 50.0, &
          & solar_flux_phcm2s1(i) / ABS(energy_right_eV - energy_left_eV), &
          & solar_flux_phcm2s1(i) * solar_flux_energy_bin_eV(i) / 50.0, &
          & solar_flux_phcm2s1(i) * solar_flux_energy_bin_eV(i) / ABS(energy_right_eV - energy_left_eV)

     WRITE (10, '(2x,f7.2,2x,f8.3,4(2x,e12.5))') &
          & lambda_euvac_A(i), &
          & energy_mid_eV, &
          & solar_flux_phcm2s1(i) / 50.0, &
          & solar_flux_phcm2s1(i) / ABS(energy_right_eV - energy_left_eV), &
          & solar_flux_phcm2s1(i) * solar_flux_energy_bin_eV(i) / 50.0, &
          & solar_flux_phcm2s1(i) * solar_flux_energy_bin_eV(i) / ABS(energy_right_eV - energy_left_eV)

     WRITE (10, '(2x,f7.2,2x,f8.3,4(2x,e12.5))') &
          & lambda_right_A, &
          & energy_right_eV, &
          & solar_flux_phcm2s1(i) / 50.0, &
          & solar_flux_phcm2s1(i) / ABS(energy_right_eV - energy_left_eV), &
          & solar_flux_phcm2s1(i) * solar_flux_energy_bin_eV(i) / 50.0, &
          & solar_flux_phcm2s1(i) * solar_flux_energy_bin_eV(i) / ABS(energy_right_eV - energy_left_eV)

     WRITE (10, '(2x,f7.2,2x,f8.3,4(2x,e12.5))') &
          & lambda_right_A, &
          & energy_right_eV, &
          & 1.0d-10 * solar_flux_phcm2s1(i) / 50.0, &
          & 1.0d-10 * solar_flux_phcm2s1(i) / ABS(energy_right_eV - energy_left_eV), &
          & 1.0d-10 * solar_flux_phcm2s1(i) * solar_flux_energy_bin_eV(i) / 50.0, &
          & 1.0d-10 * solar_flux_phcm2s1(i) * solar_flux_energy_bin_eV(i) / ABS(energy_right_eV - energy_left_eV)

     WRITE (10, '(" ")')
     WRITE (10, '(" ")')
  END DO
  CLOSE (10, STATUS = 'KEEP')
  PRINT '("### file  solar_EUV_spectrum_EUVAC_cont_box.dat is ready")'

!###########
! this file draws the cross-sections for the continuous part of the spectrum as steps constant within one continuous bin

  OPEN (10, FILE = 'solar_EUV_crossections_cont.dat')

  WRITE (10, '("# col  1 is the ordering number of the continuous EUVAC bin")')
  WRITE (10, '("# col  2 is the wavelength (A)")')
  WRITE (10, '("# col  3 is the energy of photons with this wavelength (eV)")')

  WRITE (10, '("# col  4 is the total photon absorption cross-section for He (10^{-18} cm^2)")')
  WRITE (10, '("# col  5 is the total photon absorption cross-section for  O (10^{-18} cm^2)")')
  WRITE (10, '("# col  6 is the total photon absorption cross-section for N2 (10^{-18} cm^2)")')
  WRITE (10, '("# col  7 is the total photon absorption cross-section for O2 (10^{-18} cm^2)")')

  WRITE (10, '("# col  8 is the ionization cross-section for He+ from He (10^{-18} cm^2)")')

  WRITE (10, '("# col  9 is the ionization cross-section for O+(4S)  from O (10^{-18} cm^2)")')
  WRITE (10, '("# col 10 is the ionization cross-section for O+(2D)  from O (10^{-18} cm^2)")')
  WRITE (10, '("# col 11 is the ionization cross-section for O+(2P)  from O (10^{-18} cm^2)")')
  WRITE (10, '("# col 12 is the ionization cross-section for O+(4P)  from O (10^{-18} cm^2)")')
  WRITE (10, '("# col 13 is the ionization cross-section for O+(2P*) from O (10^{-18} cm^2)")')

  WRITE (10, '("# col 14 is the ionization cross-section for N2+ from N2 (10^{-18} cm^2)")')
  WRITE (10, '("# col 15 is the ionization cross-section for N+  from N2 (10^{-18} cm^2)")')

  WRITE (10, '("# col 16 is the ionization cross-section for O2+ from O2 (10^{-18} cm^2)")')
  WRITE (10, '("# col 17 is the ionization cross-section for O+  from O2 (10^{-18} cm^2)")')

  DO i = 1, N_solar_bins
! skip the lines
     SELECT CASE (i)
        CASE(5,6,8,9,11,14,17,18,20,21,24,26,27,28,33,35,36)
           CYCLE
     END SELECT

     lambda_left_A  = lambda_euvac_A(i) - 25.0
     lambda_right_A = lambda_euvac_A(i) + 25.0

     energy_left_eV  = factor_energy_eVA / lambda_left_A
     energy_mid_eV   = factor_energy_eVA / lambda_euvac_A(i)
     energy_right_eV = factor_energy_eVA / lambda_right_A

     WRITE (10, '(2x,i4,2(2x,f8.3),14(2x,e12.5))') &
          & i, &              ! 1
          & lambda_left_A, &       ! 2
          & energy_left_eV, &          ! 3
          & sigma_tot_He_cm2(i), &          ! 4
          & sigma_tot_O_cm2(i), &           ! 5
          & sigma_tot_N2_cm2(i), &          ! 6
          & sigma_tot_O2_cm2(i), &          ! 7
          & sigma_ion_He_cm2(i), &     ! 8
          & sigma_ion_O_4S_cm2(i), &        ! 9
          & sigma_ion_O_2D_cm2(i), &        ! 10
          & sigma_ion_O_2P_cm2(i), &        ! 11
          & sigma_ion_O_4P_cm2(i), &        ! 12
          & sigma_ion_O_2Pst_cm2(i), &      ! 13
          & sigma_ion_N2_to_N2_cm2(i), &         ! 14
          & sigma_ion_N2_to_N_cm2(i), &          ! 15
          & sigma_ion_O2_to_O2_cm2(i), &    ! 16
          & sigma_ion_O2_to_O_cm2(i)        ! 17

     WRITE (10, '(2x,i4,2(2x,f8.3),14(2x,e12.5))') &
          & i, &           ! 1
          & lambda_euvac_A(i), &   ! 2
          & energy_mid_eV, &           ! 3
          & sigma_tot_He_cm2(i), &          ! 4
          & sigma_tot_O_cm2(i), &           ! 5
          & sigma_tot_N2_cm2(i), &          ! 6
          & sigma_tot_O2_cm2(i), &          ! 7
          & sigma_ion_He_cm2(i), &     ! 8
          & sigma_ion_O_4S_cm2(i), &        ! 9
          & sigma_ion_O_2D_cm2(i), &        ! 10
          & sigma_ion_O_2P_cm2(i), &        ! 11
          & sigma_ion_O_4P_cm2(i), &        ! 12
          & sigma_ion_O_2Pst_cm2(i), &      ! 13
          & sigma_ion_N2_to_N2_cm2(i), &         ! 14
          & sigma_ion_N2_to_N_cm2(i), &          ! 15
          & sigma_ion_O2_to_O2_cm2(i), &    ! 16
          & sigma_ion_O2_to_O_cm2(i)        ! 17

     WRITE (10, '(2x,i4,2(2x,f8.3),14(2x,e12.5))') &
          & i, &           ! 1
          & lambda_right_A, &   ! 2
          & energy_right_eV, &         ! 3
          & sigma_tot_He_cm2(i), &          ! 4
          & sigma_tot_O_cm2(i), &           ! 5
          & sigma_tot_N2_cm2(i), &          ! 6
          & sigma_tot_O2_cm2(i), &          ! 7
          & sigma_ion_He_cm2(i), &     ! 8
          & sigma_ion_O_4S_cm2(i), &        ! 9
          & sigma_ion_O_2D_cm2(i), &        ! 10
          & sigma_ion_O_2P_cm2(i), &        ! 11
          & sigma_ion_O_4P_cm2(i), &        ! 12
          & sigma_ion_O_2Pst_cm2(i), &      ! 13
          & sigma_ion_N2_to_N2_cm2(i), &         ! 14
          & sigma_ion_N2_to_N_cm2(i), &          ! 15
          & sigma_ion_O2_to_O2_cm2(i), &    ! 16
          & sigma_ion_O2_to_O_cm2(i)        ! 17

  END DO
  CLOSE (10, STATUS = 'KEEP')
  PRINT '("### file solar_EUV_crossections_cont.dat is ready")'

!###########
  OPEN (10, FILE = 'solar_EUV_crossections_lines.dat')

  WRITE (10, '("# col  1 is the ordering number of the discrete EUVAC bin")')
  WRITE (10, '("# col  2 is the wavelength (A)")')
  WRITE (10, '("# col  3 is the energy of photons with this wavelength (eV)")')

  WRITE (10, '("# col  4 is the total photon absorption cross-section for He (10^{-18} cm^2)")')
  WRITE (10, '("# col  5 is the total photon absorption cross-section for  O (10^{-18} cm^2)")')
  WRITE (10, '("# col  6 is the total photon absorption cross-section for N2 (10^{-18} cm^2)")')
  WRITE (10, '("# col  7 is the total photon absorption cross-section for O2 (10^{-18} cm^2)")')

  WRITE (10, '("# col  8 is the ionization cross-section for He+ from He (10^{-18} cm^2)")')

  WRITE (10, '("# col  9 is the ionization cross-section for O+(4S)  from O (10^{-18} cm^2)")')
  WRITE (10, '("# col 10 is the ionization cross-section for O+(2D)  from O (10^{-18} cm^2)")')
  WRITE (10, '("# col 11 is the ionization cross-section for O+(2P)  from O (10^{-18} cm^2)")')
  WRITE (10, '("# col 12 is the ionization cross-section for O+(4P)  from O (10^{-18} cm^2)")')
  WRITE (10, '("# col 13 is the ionization cross-section for O+(2P*) from O (10^{-18} cm^2)")')

  WRITE (10, '("# col 14 is the ionization cross-section for N2+ from N2 (10^{-18} cm^2)")')
  WRITE (10, '("# col 15 is the ionization cross-section for N+  from N2 (10^{-18} cm^2)")')

  WRITE (10, '("# col 16 is the ionization cross-section for O2+ from O2 (10^{-18} cm^2)")')
  WRITE (10, '("# col 17 is the ionization cross-section for O+  from O2 (10^{-18} cm^2)")')

  DO i = 1, N_solar_bins

! skip the continuous bins
     SELECT CASE (i)
        CASE(1:4,7,10,12,13,15,16,19,22,23,25,29:32,34,37)
           CYCLE
     END SELECT

     WRITE (10, '(2x,i4,2(2x,f8.3),14(2x,e12.5))') &
          & i, &           ! 1
          & lambda_euvac_A(i), &
          & solar_flux_energy_bin_eV(i), &      ! 2
          & sigma_tot_He_cm2(i), &          ! 4
          & sigma_tot_O_cm2(i), &           ! 5
          & sigma_tot_N2_cm2(i), &          ! 6
          & sigma_tot_O2_cm2(i), &          ! 7
          & sigma_ion_He_cm2(i), &     ! 8
          & sigma_ion_O_4S_cm2(i), &        ! 9
          & sigma_ion_O_2D_cm2(i), &        ! 10
          & sigma_ion_O_2P_cm2(i), &        ! 11
          & sigma_ion_O_4P_cm2(i), &        ! 12
          & sigma_ion_O_2Pst_cm2(i), &      ! 13
          & sigma_ion_N2_to_N2_cm2(i), &         ! 14
          & sigma_ion_N2_to_N_cm2(i), &          ! 15
          & sigma_ion_O2_to_O2_cm2(i), &    ! 16
          & sigma_ion_O2_to_O_cm2(i)        ! 17
  END DO
  CLOSE (10, STATUS = 'KEEP')
  PRINT '("### file solar_EUV_crossections_lines.dat is ready")'

END SUBROUTINE Initiate_EUVAC_flux_crossections
