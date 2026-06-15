!---------------------------------------------------------------
! the wavelengths and fluxes are from file f76ref.dat downloaded from
!
! ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_UV/AE_E/REF_SPEC/f74113.dat
!
! this link can be found here: https://www.ngdc.noaa.gov/stp/solar/solaruv.html
!
! The data were edited: fluxes equal to 0 and 0.1 (in original units which are 1e10/m2/s) were removed,
! fluxes with wavelength shorter than 28.47 A and longer than 1025 A were removed
!
SUBROUTINE Use_f74113(lambda_EUV_A)

  USE Photoelectrons, ONLY : N_solar_bins, solar_flux_phcm2s1, solar_flux_energy_bin_eV
  USE GlobalIndices, ONLY : f10p7, f10p7_81
  USE PhysicalConstants

  IMPLICIT NONE

  REAL, INTENT(INOUT) :: lambda_EUV_A(N_solar_bins)

  REAL factor_f107

  INTEGER n1, n2
  REAL scale_factor
  INTEGER i

  REAL factor_energy_eVA
  
! solar EUV wavelengths (A) and flux (10^10 photons/m^2/s)

! lines shorter than 23.7 A are omitted because the photoionization cross-sectionsi in Fennelly's paper are given for 
! wavelenths like this or longer

! below 50A
 lambda_EUV_A(  1) =   28.47 ; solar_flux_phcm2s1(  1) =    0.8  !   # C VI
 lambda_EUV_A(  2) =   28.79 ; solar_flux_phcm2s1(  2) =    3.2  !   # N VI
 lambda_EUV_A(  3) =   29.52 ; solar_flux_phcm2s1(  3) =    2.8  !   # N VI
 lambda_EUV_A(  4) =   30.02 ; solar_flux_phcm2s1(  4) =    0.7  !   #UNSPEC'D
 lambda_EUV_A(  5) =   30.43 ; solar_flux_phcm2s1(  5) =    1.0  !   # S XIV
 lambda_EUV_A(  6) =   33.74 ; solar_flux_phcm2s1(  6) =    1.8  !   # C VI
 lambda_EUV_A(  7) =   40.95 ; solar_flux_phcm2s1(  7) =    0.8  !   #SI XII
 lambda_EUV_A(  8) =   43.76 ; solar_flux_phcm2s1(  8) =    2.4  !   #SI XI
 lambda_EUV_A(  9) =   44.02 ; solar_flux_phcm2s1(  9) =    1.0  !   #SI XII
 lambda_EUV_A( 10) =   44.16 ; solar_flux_phcm2s1( 10) =    1.2  !   #SI XII
 lambda_EUV_A( 11) =   45.66 ; solar_flux_phcm2s1( 11) =    0.6  !   #SI XII
 lambda_EUV_A( 12) =   46.40 ; solar_flux_phcm2s1( 12) =    3.0  !   #SI XI
 lambda_EUV_A( 13) =   46.67 ; solar_flux_phcm2s1( 13) =    3.2  !   #UNSPEC'D
 lambda_EUV_A( 14) =   47.87 ; solar_flux_phcm2s1( 14) =    3.6  !   #UNSPEC'D
 lambda_EUV_A( 15) =   49.22 ; solar_flux_phcm2s1( 15) =    4.8  !   #SI XI

! 50-100A, EUVAC bin 1
 lambda_EUV_A( 16) =   50.36 ; solar_flux_phcm2s1( 16) =    0.2  !   #FE XVI
 lambda_EUV_A( 17) =   50.52 ; solar_flux_phcm2s1( 17) =    6.0  !   #SI X
 lambda_EUV_A( 18) =   50.69 ; solar_flux_phcm2s1( 18) =    6.0  !   #SI X
 lambda_EUV_A( 19) =   52.30 ; solar_flux_phcm2s1( 19) =    3.2  !   #SI XI
 lambda_EUV_A( 20) =   52.91 ; solar_flux_phcm2s1( 20) =    0.4  !   #FE XV
 lambda_EUV_A( 21) =   54.15 ; solar_flux_phcm2s1( 21) =    5.3  !   #BLEND
 lambda_EUV_A( 22) =   54.42 ; solar_flux_phcm2s1( 22) =    2.3  !   #UNSPEC'D
 lambda_EUV_A( 23) =   54.70 ; solar_flux_phcm2s1( 23) =    0.2  !   #FE XVI
 lambda_EUV_A( 24) =   55.06 ; solar_flux_phcm2s1( 24) =    2.8  !   #MG IX
 lambda_EUV_A( 25) =   55.34 ; solar_flux_phcm2s1( 25) =    8.9  !   #SI IX
 lambda_EUV_A( 26) =   56.08 ; solar_flux_phcm2s1( 26) =    2.0  !   # S IX
 lambda_EUV_A( 27) =   56.92 ; solar_flux_phcm2s1( 27) =    4.8  !   #UNSPEC'D
 lambda_EUV_A( 28) =   57.36 ; solar_flux_phcm2s1( 28) =    4.0  !   #UNSPEC'D
 lambda_EUV_A( 29) =   57.56 ; solar_flux_phcm2s1( 29) =    3.2  !   #UNSPEC'D
 lambda_EUV_A( 30) =   57.88 ; solar_flux_phcm2s1( 30) =    2.9  !   #MG X
 lambda_EUV_A( 31) =   58.96 ; solar_flux_phcm2s1( 31) =    1.5  !   #FE XIV
 lambda_EUV_A( 32) =   59.62 ; solar_flux_phcm2s1( 32) =    1.4  !   #FE XIV
 lambda_EUV_A( 33) =   60.30 ; solar_flux_phcm2s1( 33) =    2.0  !   #UNSPEC'D
 lambda_EUV_A( 34) =   60.85 ; solar_flux_phcm2s1( 34) =    2.9  !   #UNSPEC'D
 lambda_EUV_A( 35) =   61.07 ; solar_flux_phcm2s1( 35) =    5.2  !   #SI VIII
 lambda_EUV_A( 36) =   61.63 ; solar_flux_phcm2s1( 36) =    2.4  !   # S VIII
 lambda_EUV_A( 37) =   61.90 ; solar_flux_phcm2s1( 37) =    4.0  !   #BLEND
 lambda_EUV_A( 38) =   62.35 ; solar_flux_phcm2s1( 38) =    1.7  !   #FE XII
 lambda_EUV_A( 39) =   62.77 ; solar_flux_phcm2s1( 39) =    3.0  !   #MG IX
 lambda_EUV_A( 40) =   63.16 ; solar_flux_phcm2s1( 40) =    2.8  !   #MG X
 lambda_EUV_A( 41) =   63.30 ; solar_flux_phcm2s1( 41) =    4.5  !   #MG X
 lambda_EUV_A( 42) =   63.65 ; solar_flux_phcm2s1( 42) =    3.3  !   #UNSPEC'D
 lambda_EUV_A( 43) =   63.72 ; solar_flux_phcm2s1( 43) =    0.2  !   #FE XVI
 lambda_EUV_A( 44) =   64.11 ; solar_flux_phcm2s1( 44) =    2.2  !   #FE XIII
 lambda_EUV_A( 45) =   64.60 ; solar_flux_phcm2s1( 45) =    2.0  !   #UNSPEC'D
 lambda_EUV_A( 46) =   65.21 ; solar_flux_phcm2s1( 46) =    2.4  !   #UNSPEC'D
 lambda_EUV_A( 47) =   65.71 ; solar_flux_phcm2s1( 47) =    3.3  !   #UNSPEC'D
 lambda_EUV_A( 48) =   65.85 ; solar_flux_phcm2s1( 48) =    2.5  !   #MG X
 lambda_EUV_A( 49) =   66.30 ; solar_flux_phcm2s1( 49) =    6.0  !   #FE XII
 lambda_EUV_A( 50) =   67.14 ; solar_flux_phcm2s1( 50) =    2.6  !   #MG IX
 lambda_EUV_A( 51) =   67.35 ; solar_flux_phcm2s1( 51) =    3.1  !   #FE XII
 lambda_EUV_A( 52) =   68.35 ; solar_flux_phcm2s1( 52) =    1.8  !   #UNSPEC'D
 lambda_EUV_A( 53) =   69.65 ; solar_flux_phcm2s1( 53) =   11.0  !   #UNSPEC'D
 lambda_EUV_A( 54) =   70.00 ; solar_flux_phcm2s1( 54) =    0.4  !   #FE XV
 lambda_EUV_A( 55) =   70.54 ; solar_flux_phcm2s1( 55) =    2.6  !   #UNSPEC'D
 lambda_EUV_A( 56) =   70.75 ; solar_flux_phcm2s1( 56) =    2.4  !   #UNSPEC'D
 lambda_EUV_A( 57) =   71.00 ; solar_flux_phcm2s1( 57) =    3.5  !   #UNSPEC'D
 lambda_EUV_A( 58) =   71.94 ; solar_flux_phcm2s1( 58) =    3.1  !   #FE XIV
 lambda_EUV_A( 59) =   72.31 ; solar_flux_phcm2s1( 59) =    5.1  !   #BLEND
 lambda_EUV_A( 60) =   72.63 ; solar_flux_phcm2s1( 60) =    1.9  !   #FE XI
 lambda_EUV_A( 61) =   72.80 ; solar_flux_phcm2s1( 61) =    1.7  !   #UNSPEC'D
 lambda_EUV_A( 62) =   72.95 ; solar_flux_phcm2s1( 62) =    2.7  !   #UNSPEC'D
 lambda_EUV_A( 63) =   73.55 ; solar_flux_phcm2s1( 63) =    2.0  !   #NE VIII
 lambda_EUV_A( 64) =   74.21 ; solar_flux_phcm2s1( 64) =    2.0  !   #UNSPEC'D
 lambda_EUV_A( 65) =   74.44 ; solar_flux_phcm2s1( 65) =    1.0  !   #UNSPEC'D
 lambda_EUV_A( 66) =   74.83 ; solar_flux_phcm2s1( 66) =    3.2  !   #BLEND
 lambda_EUV_A( 67) =   75.03 ; solar_flux_phcm2s1( 67) =    4.1  !   #MG VIII
 lambda_EUV_A( 68) =   75.29 ; solar_flux_phcm2s1( 68) =    2.0  !   #UNSPEC'D
 lambda_EUV_A( 69) =   75.46 ; solar_flux_phcm2s1( 69) =    3.0  !   #UNSPEC'D
 lambda_EUV_A( 70) =   75.73 ; solar_flux_phcm2s1( 70) =    2.0  !   #UNSPEC'D
 lambda_EUV_A( 71) =   76.01 ; solar_flux_phcm2s1( 71) =    5.6  !   #FE XIII
 lambda_EUV_A( 72) =   76.48 ; solar_flux_phcm2s1( 72) =    2.0  !   #FE XIII
 lambda_EUV_A( 73) =   76.83 ; solar_flux_phcm2s1( 73) =    3.3  !   #BLEND
 lambda_EUV_A( 74) =   76.94 ; solar_flux_phcm2s1( 74) =    2.6  !   #UNSPEC'D
 lambda_EUV_A( 75) =   77.30 ; solar_flux_phcm2s1( 75) =    2.3  !   #UNSPEC'D
 lambda_EUV_A( 76) =   77.74 ; solar_flux_phcm2s1( 76) =    3.4  !   #MG IX
 lambda_EUV_A( 77) =   78.56 ; solar_flux_phcm2s1( 77) =    2.4  !   #UNSPEC'D
 lambda_EUV_A( 78) =   78.70 ; solar_flux_phcm2s1( 78) =    2.9  !   #NI XI
 lambda_EUV_A( 79) =   79.08 ; solar_flux_phcm2s1( 79) =    1.5  !   #UNSPEC'D
 lambda_EUV_A( 80) =   79.48 ; solar_flux_phcm2s1( 80) =    2.8  !   #FE XII
 lambda_EUV_A( 81) =   79.76 ; solar_flux_phcm2s1( 81) =    1.8  !   #UNSPEC'D
 lambda_EUV_A( 82) =   80.00 ; solar_flux_phcm2s1( 82) =    2.3  !   #FE XII
 lambda_EUV_A( 83) =   80.21 ; solar_flux_phcm2s1( 83) =    2.3  !   #UNSPEC'D
 lambda_EUV_A( 84) =   80.55 ; solar_flux_phcm2s1( 84) =    3.6  !   #FE XII
 lambda_EUV_A( 85) =   80.94 ; solar_flux_phcm2s1( 85) =    1.6  !   #UNSPEC'D
 lambda_EUV_A( 86) =   81.16 ; solar_flux_phcm2s1( 86) =    2.1  !   #UNSPEC'D
 lambda_EUV_A( 87) =   81.58 ; solar_flux_phcm2s1( 87) =    2.4  !   #UNSPEC'D
 lambda_EUV_A( 88) =   81.94 ; solar_flux_phcm2s1( 88) =    3.4  !   #UNSPEC'D
 lambda_EUV_A( 89) =   82.43 ; solar_flux_phcm2s1( 89) =    4.8  !   #FE IX
 lambda_EUV_A( 90) =   82.67 ; solar_flux_phcm2s1( 90) =    6.7  !   #UNSPEC'D
 lambda_EUV_A( 91) =   83.25 ; solar_flux_phcm2s1( 91) =    3.5  !   #UNSPEC'D
 lambda_EUV_A( 92) =   83.42 ; solar_flux_phcm2s1( 92) =    4.1  !   #MG VII
 lambda_EUV_A( 93) =   83.67 ; solar_flux_phcm2s1( 93) =    3.5  !   #MG VII
 lambda_EUV_A( 94) =   84.00 ; solar_flux_phcm2s1( 94) =    4.7  !   #MG VII
 lambda_EUV_A( 95) =   84.26 ; solar_flux_phcm2s1( 95) =    3.0  !   #UNSPEC'D
 lambda_EUV_A( 96) =   84.50 ; solar_flux_phcm2s1( 96) =    4.7  !   #UNSPEC'D
 lambda_EUV_A( 97) =   84.72 ; solar_flux_phcm2s1( 97) =    3.3  !   #UNSPEC'D
 lambda_EUV_A( 98) =   84.86 ; solar_flux_phcm2s1( 98) =    3.1  !   #UNSPEC'D
 lambda_EUV_A( 99) =   85.16 ; solar_flux_phcm2s1( 99) =    2.4  !   #UNSPEC'D
 lambda_EUV_A(100) =   85.50 ; solar_flux_phcm2s1(100) =    4.7  !   #UNSPEC'D
 lambda_EUV_A(101) =   85.69 ; solar_flux_phcm2s1(101) =    2.2  !   #UNSPEC'D
 lambda_EUV_A(102) =   85.87 ; solar_flux_phcm2s1(102) =    2.3  !   #UNSPEC'D
 lambda_EUV_A(103) =   86.23 ; solar_flux_phcm2s1(103) =    1.7  !   #UNSPEC'D
 lambda_EUV_A(104) =   86.40 ; solar_flux_phcm2s1(104) =    1.4  !   #UNSPEC'D
 lambda_EUV_A(105) =   86.77 ; solar_flux_phcm2s1(105) =    5.9  !   #FE XI
 lambda_EUV_A(106) =   86.98 ; solar_flux_phcm2s1(106) =    3.6  !   #FE XI
 lambda_EUV_A(107) =   87.30 ; solar_flux_phcm2s1(107) =    1.9  !   #UNSPEC'D
 lambda_EUV_A(108) =   87.61 ; solar_flux_phcm2s1(108) =    1.6  !   #UNSPEC'D
 lambda_EUV_A(109) =   88.11 ; solar_flux_phcm2s1(109) =    4.9  !   #BLEND
 lambda_EUV_A(110) =   88.42 ; solar_flux_phcm2s1(110) =    1.5  !   #UNSPEC'D
 lambda_EUV_A(111) =   88.64 ; solar_flux_phcm2s1(111) =    1.9  !   #UNSPEC'D
 lambda_EUV_A(112) =   88.90 ; solar_flux_phcm2s1(112) =    4.5  !   #FE XI
 lambda_EUV_A(113) =   89.14 ; solar_flux_phcm2s1(113) =    3.2  !   #FE XI
 lambda_EUV_A(114) =   89.70 ; solar_flux_phcm2s1(114) =    3.5  !   #FE XI
 lambda_EUV_A(115) =   90.14 ; solar_flux_phcm2s1(115) =    3.7  !   #FE XI
 lambda_EUV_A(116) =   90.45 ; solar_flux_phcm2s1(116) =    2.3  !   #FE XI
 lambda_EUV_A(117) =   90.71 ; solar_flux_phcm2s1(117) =    2.4  !   #UNSPEC'D
 lambda_EUV_A(118) =   91.00 ; solar_flux_phcm2s1(118) =    3.1  !   #UNSPEC'D
 lambda_EUV_A(119) =   91.48 ; solar_flux_phcm2s1(119) =    1.5  !   #UNSPEC'D
 lambda_EUV_A(120) =   91.69 ; solar_flux_phcm2s1(120) =    4.6  !   #NI X
 lambda_EUV_A(121) =   91.81 ; solar_flux_phcm2s1(121) =    3.8  !   #UNSPEC'D
 lambda_EUV_A(122) =   92.09 ; solar_flux_phcm2s1(122) =    3.0  !   #UNSPEC'D
 lambda_EUV_A(123) =   92.55 ; solar_flux_phcm2s1(123) =    2.4  !   #UNSPEC'D
 lambda_EUV_A(124) =   92.81 ; solar_flux_phcm2s1(124) =    3.0  !   #UNSPEC'D
 lambda_EUV_A(125) =   93.61 ; solar_flux_phcm2s1(125) =    4.5  !   #UNSPEC'D
 lambda_EUV_A(126) =   94.07 ; solar_flux_phcm2s1(126) =    6.8  !   #FE X
 lambda_EUV_A(127) =   94.39 ; solar_flux_phcm2s1(127) =    1.4  !   #UNSPEC'D
 lambda_EUV_A(128) =   94.81 ; solar_flux_phcm2s1(128) =    1.1  !   #UNSPEC'D
 lambda_EUV_A(129) =   95.37 ; solar_flux_phcm2s1(129) =    5.0  !   #FE X
 lambda_EUV_A(130) =   95.51 ; solar_flux_phcm2s1(130) =    2.3  !   #UNSPEC'D
 lambda_EUV_A(131) =   95.81 ; solar_flux_phcm2s1(131) =    2.3  !   #UNSPEC'D
 lambda_EUV_A(132) =   96.05 ; solar_flux_phcm2s1(132) =    9.2  !   #FE X
 lambda_EUV_A(133) =   96.49 ; solar_flux_phcm2s1(133) =    1.5  !   #UNSPEC'D
 lambda_EUV_A(134) =   96.83 ; solar_flux_phcm2s1(134) =    2.6  !   #FE X
 lambda_EUV_A(135) =   97.12 ; solar_flux_phcm2s1(135) =    5.0  !   #FE X
 lambda_EUV_A(136) =   97.51 ; solar_flux_phcm2s1(136) =    2.4  !   #BLEND
 lambda_EUV_A(137) =   97.87 ; solar_flux_phcm2s1(137) =    2.2  !   #FE X
 lambda_EUV_A(138) =   98.12 ; solar_flux_phcm2s1(138) =    4.5  !   #BLEND
 lambda_EUV_A(139) =   98.23 ; solar_flux_phcm2s1(139) =    5.2  !   #BLEND
 lambda_EUV_A(140) =   98.50 ; solar_flux_phcm2s1(140) =    2.1  !   #UNSPEC'D
 lambda_EUV_A(141) =   98.88 ; solar_flux_phcm2s1(141) =    1.2  !   #UNSPEC'D
 lambda_EUV_A(142) =   99.44 ; solar_flux_phcm2s1(142) =    1.3  !   #UNSPEC'D
 lambda_EUV_A(143) =   99.71 ; solar_flux_phcm2s1(143) =    1.6  !   #UNSPEC'D
 lambda_EUV_A(144) =   99.99 ; solar_flux_phcm2s1(144) =    2.2  !   #UNSPEC'D

! 100-150A, EUVAC bin 2
 lambda_EUV_A(145) =  100.54 ; solar_flux_phcm2s1(145) =    6.7  !   #UNSPEC'D
 lambda_EUV_A(146) =  100.96 ; solar_flux_phcm2s1(146) =    2.1  !   #UNSPEC'D
 lambda_EUV_A(147) =  101.57 ; solar_flux_phcm2s1(147) =    2.7  !   #UNSPEC'D
 lambda_EUV_A(148) =  102.15 ; solar_flux_phcm2s1(148) =    3.6  !   #UNSPEC'D
 lambda_EUV_A(149) =  103.01 ; solar_flux_phcm2s1(149) =    1.6  !   #NE VIII
 lambda_EUV_A(150) =  103.17 ; solar_flux_phcm2s1(150) =    3.4  !   #UNSPEC'D
 lambda_EUV_A(151) =  103.58 ; solar_flux_phcm2s1(151) =    6.0  !   #FE IX
 lambda_EUV_A(152) =  103.94 ; solar_flux_phcm2s1(152) =    4.9  !   #UNSPEC'D
 lambda_EUV_A(153) =  104.23 ; solar_flux_phcm2s1(153) =    1.3  !   #UNSPEC'D
 lambda_EUV_A(154) =  104.76 ; solar_flux_phcm2s1(154) =    1.6  !   #UNSPEC'D
 lambda_EUV_A(155) =  105.23 ; solar_flux_phcm2s1(155) =    5.0  !   #FE IX
 lambda_EUV_A(156) =  106.25 ; solar_flux_phcm2s1(156) =    1.8  !   #NE VII
 lambda_EUV_A(157) =  106.57 ; solar_flux_phcm2s1(157) =    0.9  !   #UNSPEC'D
 lambda_EUV_A(158) =  106.93 ; solar_flux_phcm2s1(158) =    1.5  !   #UNSPEC'D
 lambda_EUV_A(159) =  108.05 ; solar_flux_phcm2s1(159) =    1.5  !   #FE VIII
 lambda_EUV_A(160) =  108.46 ; solar_flux_phcm2s1(160) =    1.6  !   #UNSPEC'D
 lambda_EUV_A(161) =  109.50 ; solar_flux_phcm2s1(161) =    1.6  !   #UNSPEC'D
 lambda_EUV_A(162) =  110.76 ; solar_flux_phcm2s1(162) =    1.2  !   #UNSPEC'D
 lambda_EUV_A(163) =  111.25 ; solar_flux_phcm2s1(163) =    3.8  !   #UNSPEC'D
 lambda_EUV_A(164) =  113.80 ; solar_flux_phcm2s1(164) =    2.4  !   #UNSPEC'D
 lambda_EUV_A(165) =  114.09 ; solar_flux_phcm2s1(165) =    2.1  !   #UNSPEC'D
 lambda_EUV_A(166) =  115.82 ; solar_flux_phcm2s1(166) =    1.9  !   #UNSPEC'D
 lambda_EUV_A(167) =  116.75 ; solar_flux_phcm2s1(167) =    3.1  !   #UNSPEC'D
 lambda_EUV_A(168) =  117.20 ; solar_flux_phcm2s1(168) =    2.1  !   #UNSPEC'D
 lambda_EUV_A(169) =  121.79 ; solar_flux_phcm2s1(169) =    0.9  !   #UNSPEC'D
 lambda_EUV_A(170) =  122.70 ; solar_flux_phcm2s1(170) =    3.7  !   #NE VI
 lambda_EUV_A(171) =  123.50 ; solar_flux_phcm2s1(171) =    2.1  !   #UNSPEC'D
 lambda_EUV_A(172) =  127.65 ; solar_flux_phcm2s1(172) =    3.6  !   #UNSPEC'D
 lambda_EUV_A(173) =  129.87 ; solar_flux_phcm2s1(173) =    3.7  !   # O VI
 lambda_EUV_A(174) =  131.02 ; solar_flux_phcm2s1(174) =    4.2  !   #BLEND
 lambda_EUV_A(175) =  131.21 ; solar_flux_phcm2s1(175) =    3.8  !   #BLEND
 lambda_EUV_A(176) =  141.20 ; solar_flux_phcm2s1(176) =    9.0  !   #UNSPEC'D
 lambda_EUV_A(177) =  144.27 ; solar_flux_phcm2s1(177) =    0.9  !   #NI X
 lambda_EUV_A(178) =  145.04 ; solar_flux_phcm2s1(178) =   11.9  !   #NI X
 lambda_EUV_A(179) =  148.40 ; solar_flux_phcm2s1(179) =   40.0  !   #NI XI

! 150-200A, EUVAC bin 3
 lambda_EUV_A(180) =  150.10 ; solar_flux_phcm2s1(180) =   19.2  !   # O VI
 lambda_EUV_A(181) =  152.15 ; solar_flux_phcm2s1(181) =   19.2  !   #NI XII
 lambda_EUV_A(182) =  154.18 ; solar_flux_phcm2s1(182) =   12.8  !   #NI XII
 lambda_EUV_A(183) =  157.73 ; solar_flux_phcm2s1(183) =    7.5  !   #NI XIII
 lambda_EUV_A(184) =  158.37 ; solar_flux_phcm2s1(184) =   12.0  !   #NI X
 lambda_EUV_A(185) =  159.98 ; solar_flux_phcm2s1(185) =   11.0  !   #NI X
 lambda_EUV_A(186) =  160.37 ; solar_flux_phcm2s1(186) =    6.0  !   #UNSPEC'D
 lambda_EUV_A(187) =  164.15 ; solar_flux_phcm2s1(187) =    3.2  !   #NI XIV
 lambda_EUV_A(188) =  167.50 ; solar_flux_phcm2s1(188) =   19.3  !   #FE VIII
 lambda_EUV_A(189) =  168.17 ; solar_flux_phcm2s1(189) =   35.6  !   #FE VIII
 lambda_EUV_A(190) =  168.55 ; solar_flux_phcm2s1(190) =   20.2  !   #FE VIII
 lambda_EUV_A(191) =  168.92 ; solar_flux_phcm2s1(191) =   12.4  !   #FE VIII
 lambda_EUV_A(192) =  169.70 ; solar_flux_phcm2s1(192) =   18.7  !   #UNSPEC'D
 lambda_EUV_A(193) =  171.08 ; solar_flux_phcm2s1(193) =  300.0  !   #FE IX
 lambda_EUV_A(194) =  172.17 ; solar_flux_phcm2s1(194) =    9.6  !   # O V
 lambda_EUV_A(195) =  173.08 ; solar_flux_phcm2s1(195) =   17.1  !   #BLEND
 lambda_EUV_A(196) =  174.58 ; solar_flux_phcm2s1(196) =  256.0  !   #FE X
 lambda_EUV_A(197) =  175.26 ; solar_flux_phcm2s1(197) =   26.5  !   #BLEND
 lambda_EUV_A(198) =  177.24 ; solar_flux_phcm2s1(198) =  136.0  !   #FE X
 lambda_EUV_A(199) =  178.05 ; solar_flux_phcm2s1(199) =   18.2  !   #FE XI
 lambda_EUV_A(200) =  179.27 ; solar_flux_phcm2s1(200) =    0.4  !   #NI XV
 lambda_EUV_A(201) =  179.75 ; solar_flux_phcm2s1(201) =   18.4  !   #FE XI
 lambda_EUV_A(202) =  180.41 ; solar_flux_phcm2s1(202) =  220.0  !   #FE XI
 lambda_EUV_A(203) =  181.14 ; solar_flux_phcm2s1(203) =   25.0  !   #FE XI
 lambda_EUV_A(204) =  182.17 ; solar_flux_phcm2s1(204) =   66.0  !   #FE XI
 lambda_EUV_A(205) =  183.45 ; solar_flux_phcm2s1(205) =   16.9  !   #CA XIV
 lambda_EUV_A(206) =  184.53 ; solar_flux_phcm2s1(206) =  103.7  !   #FE X
 lambda_EUV_A(207) =  184.80 ; solar_flux_phcm2s1(207) =    6.9  !   #FE XI
 lambda_EUV_A(208) =  185.21 ; solar_flux_phcm2s1(208) =   44.5  !   #FE VIII B
 lambda_EUV_A(209) =  186.60 ; solar_flux_phcm2s1(209) =   25.8  !   #CA XIV
 lambda_EUV_A(210) =  186.87 ; solar_flux_phcm2s1(210) =   70.0  !   # S XI   B
 lambda_EUV_A(211) =  187.95 ; solar_flux_phcm2s1(211) =    1.9  !   #AR XIV
 lambda_EUV_A(212) =  188.31 ; solar_flux_phcm2s1(212) =  191.0  !   #FE XI
 lambda_EUV_A(213) =  188.23 ; solar_flux_phcm2s1(213) =   15.4  !   #FE XII
 lambda_EUV_A(214) =  190.02 ; solar_flux_phcm2s1(214) =   48.0  !   #FE X
 lambda_EUV_A(215) =  191.04 ; solar_flux_phcm2s1(215) =   19.9  !   #FE XII
 lambda_EUV_A(216) =  191.34 ; solar_flux_phcm2s1(216) =   22.5  !   #FE XIII B
 lambda_EUV_A(217) =  192.40 ; solar_flux_phcm2s1(217) =   64.0  !   #FE XII
 lambda_EUV_A(218) =  192.82 ; solar_flux_phcm2s1(218) =   76.8  !   #FE XI
 lambda_EUV_A(219) =  193.52 ; solar_flux_phcm2s1(219) =  100.0  !   #FE XII
 lambda_EUV_A(220) =  195.13 ; solar_flux_phcm2s1(220) =  174.0  !   #FE XII
 lambda_EUV_A(221) =  196.52 ; solar_flux_phcm2s1(221) =   45.3  !   #FE XIII
 lambda_EUV_A(222) =  196.65 ; solar_flux_phcm2s1(222) =    9.9  !   #FE XII
 lambda_EUV_A(223) =  197.44 ; solar_flux_phcm2s1(223) =   18.2  !   #FE XIII
 lambda_EUV_A(224) =  198.58 ; solar_flux_phcm2s1(224) =   23.3  !   #FE XII  B

! 200-250A, EUVAC bin 4
 lambda_EUV_A(225) =  200.02 ; solar_flux_phcm2s1(225) =   37.9  !   #FE XIII
 lambda_EUV_A(226) =  201.13 ; solar_flux_phcm2s1(226) =   66.0  !   #FE XIII B
 lambda_EUV_A(227) =  202.05 ; solar_flux_phcm2s1(227) =   96.0  !   #FE XIII
 lambda_EUV_A(228) =  202.64 ; solar_flux_phcm2s1(228) =   20.7  !   #BLEND
 lambda_EUV_A(229) =  203.81 ; solar_flux_phcm2s1(229) =   39.0  !   #FE XIII
 lambda_EUV_A(230) =  204.25 ; solar_flux_phcm2s1(230) =   15.0  !   #FE XIII
 lambda_EUV_A(231) =  204.94 ; solar_flux_phcm2s1(231) =    9.9  !   #FE XIII
 lambda_EUV_A(232) =  206.26 ; solar_flux_phcm2s1(232) =    0.9  !   #UNSPEC'D
 lambda_EUV_A(233) =  206.38 ; solar_flux_phcm2s1(233) =    0.9  !   #UNSPEC'D
 lambda_EUV_A(234) =  207.46 ; solar_flux_phcm2s1(234) =    0.9  !   #UNSPEC'D
 lambda_EUV_A(235) =  208.33 ; solar_flux_phcm2s1(235) =    1.9  !   # S X    B
 lambda_EUV_A(236) =  209.63 ; solar_flux_phcm2s1(236) =    1.9  !   #FE XIII
 lambda_EUV_A(237) =  209.78 ; solar_flux_phcm2s1(237) =    0.9  !   #UNSPEC'D
 lambda_EUV_A(238) =  211.32 ; solar_flux_phcm2s1(238) =   54.0  !   #FE XIV
 lambda_EUV_A(239) =  212.14 ; solar_flux_phcm2s1(239) =   11.8  !   # S XII
 lambda_EUV_A(240) =  213.78 ; solar_flux_phcm2s1(240) =    9.9  !   #FEXIII
 lambda_EUV_A(241) =  214.75 ; solar_flux_phcm2s1(241) =    7.6  !   #SI VIII
 lambda_EUV_A(242) =  215.16 ; solar_flux_phcm2s1(242) =   21.7  !   # S XII
 lambda_EUV_A(243) =  216.88 ; solar_flux_phcm2s1(243) =   32.8  !   #SI VIII B
 lambda_EUV_A(244) =  217.00 ; solar_flux_phcm2s1(244) =   51.0  !   #UNSPEC'D
 lambda_EUV_A(245) =  218.19 ; solar_flux_phcm2s1(245) =   65.0  !   # S XII
 lambda_EUV_A(246) =  219.13 ; solar_flux_phcm2s1(246) =   12.0  !   #FE XIV
 lambda_EUV_A(247) =  220.08 ; solar_flux_phcm2s1(247) =   17.0  !   #FE XIV
 lambda_EUV_A(248) =  221.44 ; solar_flux_phcm2s1(248) =   20.1  !   # S XII
 lambda_EUV_A(249) =  221.82 ; solar_flux_phcm2s1(249) =    1.9  !   #FE XIII
 lambda_EUV_A(250) =  224.74 ; solar_flux_phcm2s1(250) =   37.4  !   # S IX   B
 lambda_EUV_A(251) =  225.12 ; solar_flux_phcm2s1(251) =   51.0  !   #BLEND
 lambda_EUV_A(252) =  227.01 ; solar_flux_phcm2s1(252) =   90.3  !   #SI IX
 lambda_EUV_A(253) =  227.19 ; solar_flux_phcm2s1(253) =    1.9  !   #FE XV
 lambda_EUV_A(254) =  227.47 ; solar_flux_phcm2s1(254) =   29.8  !   # S XII
 lambda_EUV_A(255) =  228.70 ; solar_flux_phcm2s1(255) =   18.4  !   # S X
 lambda_EUV_A(256) =  230.65 ; solar_flux_phcm2s1(256) =   10.5  !   #HE II
 lambda_EUV_A(257) =  231.55 ; solar_flux_phcm2s1(257) =   12.8  !   #HE II
 lambda_EUV_A(258) =  232.60 ; solar_flux_phcm2s1(258) =   18.7  !   #HE II
 lambda_EUV_A(259) =  233.84 ; solar_flux_phcm2s1(259) =   10.0  !   #FE XV
 lambda_EUV_A(260) =  234.38 ; solar_flux_phcm2s1(260) =   82.6  !   #BLEND
 lambda_EUV_A(261) =  237.33 ; solar_flux_phcm2s1(261) =   55.0  !   #HE II
 lambda_EUV_A(262) =  239.87 ; solar_flux_phcm2s1(262) =   16.9  !   # S XI
 lambda_EUV_A(263) =  240.71 ; solar_flux_phcm2s1(263) =   54.0  !   #FE XIII
 lambda_EUV_A(264) =  241.74 ; solar_flux_phcm2s1(264) =  138.0  !   #FE IX
 lambda_EUV_A(265) =  243.03 ; solar_flux_phcm2s1(265) =   98.0  !   #HE II
 lambda_EUV_A(266) =  243.86 ; solar_flux_phcm2s1(266) =   50.0  !   #UNSPEC'D
 lambda_EUV_A(267) =  244.92 ; solar_flux_phcm2s1(267) =   67.5  !   #FE IX
 lambda_EUV_A(268) =  245.94 ; solar_flux_phcm2s1(268) =   14.9  !   #BLEND
 lambda_EUV_A(269) =  246.24 ; solar_flux_phcm2s1(269) =   50.0  !   #BLEND
 lambda_EUV_A(270) =  246.91 ; solar_flux_phcm2s1(270) =   18.7  !   #BLEND
 lambda_EUV_A(271) =  247.18 ; solar_flux_phcm2s1(271) =   19.9  !   #BLEND
 lambda_EUV_A(272) =  249.18 ; solar_flux_phcm2s1(272) =   20.2  !   #NI XVII

! 250-300A, EUVAC bin 7
 lambda_EUV_A(273) =  251.10 ; solar_flux_phcm2s1(273) =    9.6  !   #FE XVI
 lambda_EUV_A(274) =  251.95 ; solar_flux_phcm2s1(274) =   73.0  !   #FE XIII
 lambda_EUV_A(275) =  252.19 ; solar_flux_phcm2s1(275) =   35.0  !   #FE XIV
 lambda_EUV_A(276) =  253.78 ; solar_flux_phcm2s1(276) =   29.0  !   #SI X
 lambda_EUV_A(277) =  256.32 ; solar_flux_phcm2s1(277) =  400.0  !   #HE II       ! EUVAC bin/line 5
 lambda_EUV_A(278) =  256.38 ; solar_flux_phcm2s1(278) =   60.0  !   #SI X
 lambda_EUV_A(279) =  256.64 ; solar_flux_phcm2s1(279) =   30.0  !   # S XIII
 lambda_EUV_A(280) =  257.16 ; solar_flux_phcm2s1(280) =  200.0  !   # S X    B
 lambda_EUV_A(281) =  257.39 ; solar_flux_phcm2s1(281) =  200.0  !   #FE XIV  B
 lambda_EUV_A(282) =  258.36 ; solar_flux_phcm2s1(282) =  140.0  !   #SI X
 lambda_EUV_A(283) =  259.52 ; solar_flux_phcm2s1(283) =   56.0  !   # S X
 lambda_EUV_A(284) =  261.05 ; solar_flux_phcm2s1(284) =   62.0  !   #SI X
 lambda_EUV_A(285) =  262.99 ; solar_flux_phcm2s1(285) =   10.0  !   #FE XVI
 lambda_EUV_A(286) =  264.24 ; solar_flux_phcm2s1(286) =   68.0  !   # S X
 lambda_EUV_A(287) =  264.80 ; solar_flux_phcm2s1(287) =   64.0  !   #FE XIV
 lambda_EUV_A(288) =  270.51 ; solar_flux_phcm2s1(288) =   30.0  !   #FE XIV
 lambda_EUV_A(289) =  271.99 ; solar_flux_phcm2s1(289) =   60.0  !   #SI X
 lambda_EUV_A(290) =  272.64 ; solar_flux_phcm2s1(290) =   25.0  !   #SI VII
 lambda_EUV_A(291) =  274.19 ; solar_flux_phcm2s1(291) =   90.0  !   #FE XIV
 lambda_EUV_A(292) =  275.35 ; solar_flux_phcm2s1(292) =   30.0  !   #SI VII
 lambda_EUV_A(293) =  275.67 ; solar_flux_phcm2s1(293) =   25.0  !   #SI VII
 lambda_EUV_A(294) =  276.15 ; solar_flux_phcm2s1(294) =    5.6  !   #MG VII
 lambda_EUV_A(295) =  276.84 ; solar_flux_phcm2s1(295) =   10.0  !   #BLEND
 lambda_EUV_A(296) =  277.00 ; solar_flux_phcm2s1(296) =   40.0  !   #BLEND
 lambda_EUV_A(297) =  277.27 ; solar_flux_phcm2s1(297) =   60.0  !   #BLEND
 lambda_EUV_A(298) =  278.40 ; solar_flux_phcm2s1(298) =   50.0  !   #MG VII  B
 lambda_EUV_A(299) =  281.41 ; solar_flux_phcm2s1(299) =   15.0  !   # S XI
 lambda_EUV_A(300) =  284.15 ; solar_flux_phcm2s1(300) =  210.0  !   #FE XV        ! EUVAC bin/line 6
 lambda_EUV_A(301) =  285.70 ; solar_flux_phcm2s1(301) =   24.0  !   # S XI   G
 lambda_EUV_A(302) =  289.17 ; solar_flux_phcm2s1(302) =   13.0  !   #FE XIV
 lambda_EUV_A(303) =  290.69 ; solar_flux_phcm2s1(303) =   40.0  !   #SI IX
 lambda_EUV_A(304) =  291.70 ; solar_flux_phcm2s1(304) =   18.0  !   # S XI   G
 lambda_EUV_A(305) =  292.78 ; solar_flux_phcm2s1(305) =   63.4  !   #SI IX
 lambda_EUV_A(306) =  296.19 ; solar_flux_phcm2s1(306) =   94.0  !   #SI IX
 lambda_EUV_A(307) =  299.50 ; solar_flux_phcm2s1(307) =    9.8  !   # S XII

! 300-350A, EUVAC bin 10

 lambda_EUV_A(308) =  303.31 ; solar_flux_phcm2s1(308) =  800.0  !   #SI XI        ! EUVAC bin/line 8

 lambda_EUV_A(309) =  303.78 ; solar_flux_phcm2s1(309) = 6900.0  !   #HE II        ! EUVAC bin/line 9

 lambda_EUV_A(310) =  315.02 ; solar_flux_phcm2s1(310) =  120.0  !   #MG VIII
 lambda_EUV_A(311) =  319.01 ; solar_flux_phcm2s1(311) =   25.0  !   #NI XV
 lambda_EUV_A(312) =  320.56 ; solar_flux_phcm2s1(312) =   20.0  !   #NI XVIIIB
 lambda_EUV_A(313) =  316.20 ; solar_flux_phcm2s1(313) =  100.0  !   #SI VIII
 lambda_EUV_A(314) =  319.83 ; solar_flux_phcm2s1(314) =  130.0  !   #SI VIII
 lambda_EUV_A(315) =  335.41 ; solar_flux_phcm2s1(315) =  140.0  !   #FE XVI
 lambda_EUV_A(316) =  345.13 ; solar_flux_phcm2s1(316) =  100.0  !   #SI IX
 lambda_EUV_A(317) =  345.74 ; solar_flux_phcm2s1(317) =   80.0  !   #FE X
 lambda_EUV_A(318) =  347.39 ; solar_flux_phcm2s1(318) =  150.0  !   #SI X
 lambda_EUV_A(319) =  349.85 ; solar_flux_phcm2s1(319) =  100.0  !   #SI IX

! 350-400A, EUVAC bin 12
 lambda_EUV_A(320) =  356.01 ; solar_flux_phcm2s1(320) =  110.0  !   #SI X
 lambda_EUV_A(321) =  360.80 ; solar_flux_phcm2s1(321) =   70.0  !   #FE XVI
 lambda_EUV_A(322) =  364.48 ; solar_flux_phcm2s1(322) =  120.0  !   #FE XII

 lambda_EUV_A(323) =  368.07 ; solar_flux_phcm2s1(323) =  650.0  !   #MG IX        ! EUVAC bin/line 11

 lambda_EUV_A(324) =  399.82 ; solar_flux_phcm2s1(324) =   14.0  !   #NE VI

! 400-450A, EUVAC bin 13
 lambda_EUV_A(325) =  401.14 ; solar_flux_phcm2s1(325) =   31.0  !   #NE VI
 lambda_EUV_A(326) =  401.94 ; solar_flux_phcm2s1(326) =   82.0  !   #NE VI
 lambda_EUV_A(327) =  403.26 ; solar_flux_phcm2s1(327) =   48.0  !   #NE VI
 lambda_EUV_A(328) =  417.24 ; solar_flux_phcm2s1(328) =   27.0  !   #FE XV
 lambda_EUV_A(329) =  423.00 ; solar_flux_phcm2s1(329) =    0.3  !   #HE I (C)
 lambda_EUV_A(330) =  424.00 ; solar_flux_phcm2s1(330) =    0.3  !   #HE I (C)
 lambda_EUV_A(331) =  425.00 ; solar_flux_phcm2s1(331) =    0.3  !   #HE I (C)
 lambda_EUV_A(332) =  426.00 ; solar_flux_phcm2s1(332) =    0.3  !   #HE I (C)
 lambda_EUV_A(333) =  427.00 ; solar_flux_phcm2s1(333) =    0.3  !   #HE I (C)
 lambda_EUV_A(334) =  428.00 ; solar_flux_phcm2s1(334) =    0.3  !   #HE I (C)
 lambda_EUV_A(335) =  429.00 ; solar_flux_phcm2s1(335) =    0.3  !   #HE I (C)
 lambda_EUV_A(336) =  430.00 ; solar_flux_phcm2s1(336) =    0.3  !   #HE I (C)
 lambda_EUV_A(337) =  430.47 ; solar_flux_phcm2s1(337) =   74.0  !   #MG VIII
 lambda_EUV_A(338) =  431.00 ; solar_flux_phcm2s1(338) =    0.3  !   #HE I (C)
 lambda_EUV_A(339) =  432.00 ; solar_flux_phcm2s1(339) =    0.3  !   #HE I (C)
 lambda_EUV_A(340) =  433.00 ; solar_flux_phcm2s1(340) =    0.3  !   #HE I (C)
 lambda_EUV_A(341) =  434.00 ; solar_flux_phcm2s1(341) =    0.3  !   #HE I (C)
 lambda_EUV_A(342) =  435.00 ; solar_flux_phcm2s1(342) =    0.3  !   #HE I (C)
 lambda_EUV_A(343) =  436.00 ; solar_flux_phcm2s1(343) =    0.3  !   #HE I (C)
 lambda_EUV_A(344) =  436.70 ; solar_flux_phcm2s1(344) =  110.0  !   #MG  VIII
 lambda_EUV_A(345) =  437.00 ; solar_flux_phcm2s1(345) =    0.3  !   #HE I (C)
 lambda_EUV_A(346) =  438.00 ; solar_flux_phcm2s1(346) =    0.3  !   #HE I (C)
 lambda_EUV_A(347) =  439.00 ; solar_flux_phcm2s1(347) =    0.3  !   #HE I (C)
 lambda_EUV_A(348) =  440.00 ; solar_flux_phcm2s1(348) =    0.3  !   #HE I (C)
 lambda_EUV_A(349) =  441.00 ; solar_flux_phcm2s1(349) =    0.3  !   #HE I (C)
 lambda_EUV_A(350) =  442.00 ; solar_flux_phcm2s1(350) =    0.3  !   #HE I (C)
 lambda_EUV_A(351) =  443.00 ; solar_flux_phcm2s1(351) =    0.3  !   #HE I (C)
 lambda_EUV_A(352) =  444.00 ; solar_flux_phcm2s1(352) =    0.7  !   #HE I (C)
 lambda_EUV_A(353) =  445.00 ; solar_flux_phcm2s1(353) =    0.7  !   #HE I (C)
 lambda_EUV_A(354) =  446.00 ; solar_flux_phcm2s1(354) =    0.7  !   #HE I (C)
 lambda_EUV_A(355) =  447.00 ; solar_flux_phcm2s1(355) =    0.7  !   #HE I (C)
 lambda_EUV_A(356) =  448.00 ; solar_flux_phcm2s1(356) =    0.7  !   #HE I (C)
 lambda_EUV_A(357) =  449.00 ; solar_flux_phcm2s1(357) =    0.7  !   #HE I (C)

! 450-500A, EUVAC bin 15
 lambda_EUV_A(358) =  450.00 ; solar_flux_phcm2s1(358) =    0.7  !   #HE I (C)
 lambda_EUV_A(359) =  451.00 ; solar_flux_phcm2s1(359) =    0.7  !   #HE I (C)
 lambda_EUV_A(360) =  452.00 ; solar_flux_phcm2s1(360) =    0.7  !   #HE I (C)
 lambda_EUV_A(361) =  453.00 ; solar_flux_phcm2s1(361) =    1.0  !   #HE I (C)
 lambda_EUV_A(362) =  454.00 ; solar_flux_phcm2s1(362) =    1.0  !   #HE I (C)
 lambda_EUV_A(363) =  455.00 ; solar_flux_phcm2s1(363) =    1.0  !   #HE I (C)
 lambda_EUV_A(364) =  456.00 ; solar_flux_phcm2s1(364) =    1.0  !   #HE I (C)
 lambda_EUV_A(365) =  457.00 ; solar_flux_phcm2s1(365) =    1.0  !   #HE I (C)
 lambda_EUV_A(366) =  458.00 ; solar_flux_phcm2s1(366) =    1.0  !   #HE I (C)
 lambda_EUV_A(367) =  459.00 ; solar_flux_phcm2s1(367) =    1.0  !   #HE I (C)
 lambda_EUV_A(368) =  460.00 ; solar_flux_phcm2s1(368) =    1.3  !   #HE I (C)
 lambda_EUV_A(369) =  461.00 ; solar_flux_phcm2s1(369) =    1.3  !   #HE I (C)
 lambda_EUV_A(370) =  462.00 ; solar_flux_phcm2s1(370) =    1.3  !   #HE I (C)
 lambda_EUV_A(371) =  463.00 ; solar_flux_phcm2s1(371) =    1.3  !   #HE I (C)
 lambda_EUV_A(372) =  464.00 ; solar_flux_phcm2s1(372) =    1.3  !   #HE I (C)
 lambda_EUV_A(373) =  465.00 ; solar_flux_phcm2s1(373) =    1.6  !   #HE I (C)

 lambda_EUV_A(374) =  465.22 ; solar_flux_phcm2s1(374) =  290.0  !   #NE VII      ! EUVAC bin/line 14

 lambda_EUV_A(375) =  466.00 ; solar_flux_phcm2s1(375) =    1.6  !   #HE I (C)
 lambda_EUV_A(376) =  467.00 ; solar_flux_phcm2s1(376) =    1.6  !   #HE I (C)
 lambda_EUV_A(377) =  468.00 ; solar_flux_phcm2s1(377) =    2.0  !   #HE I (C)
 lambda_EUV_A(378) =  469.00 ; solar_flux_phcm2s1(378) =    2.0  !   #HE I (C)
 lambda_EUV_A(379) =  470.00 ; solar_flux_phcm2s1(379) =    2.0  !   #HE I (C)
 lambda_EUV_A(380) =  471.00 ; solar_flux_phcm2s1(380) =    2.0  !   #HE I (C)
 lambda_EUV_A(381) =  472.00 ; solar_flux_phcm2s1(381) =    2.3  !   #HE I (C)
 lambda_EUV_A(382) =  473.00 ; solar_flux_phcm2s1(382) =    2.3  !   #HE I (C)
 lambda_EUV_A(383) =  474.00 ; solar_flux_phcm2s1(383) =    2.6  !   #HE I (C)
 lambda_EUV_A(384) =  475.00 ; solar_flux_phcm2s1(384) =    2.6  !   #HE I (C)
 lambda_EUV_A(385) =  476.00 ; solar_flux_phcm2s1(385) =    2.6  !   #HE I (C)
 lambda_EUV_A(386) =  477.00 ; solar_flux_phcm2s1(386) =    3.0  !   #HE I (C)
 lambda_EUV_A(387) =  478.00 ; solar_flux_phcm2s1(387) =    3.0  !   #HE I (C)
 lambda_EUV_A(388) =  479.00 ; solar_flux_phcm2s1(388) =    3.3  !   #HE I (C)
 lambda_EUV_A(389) =  480.00 ; solar_flux_phcm2s1(389) =    3.3  !   #HE I (C)
 lambda_EUV_A(390) =  481.00 ; solar_flux_phcm2s1(390) =    3.6  !   #HE I (C)
 lambda_EUV_A(391) =  482.00 ; solar_flux_phcm2s1(391) =    3.6  !   #HE I (C)
 lambda_EUV_A(392) =  483.00 ; solar_flux_phcm2s1(392) =    4.0  !   #HE I (C)
 lambda_EUV_A(393) =  484.00 ; solar_flux_phcm2s1(393) =    4.3  !   #HE I (C)
 lambda_EUV_A(394) =  485.00 ; solar_flux_phcm2s1(394) =    4.3  !   #HE I (C)
 lambda_EUV_A(395) =  486.00 ; solar_flux_phcm2s1(395) =    4.6  !   #HE I (C)
 lambda_EUV_A(396) =  487.00 ; solar_flux_phcm2s1(396) =    4.9  !   #HE I (C)
 lambda_EUV_A(397) =  488.00 ; solar_flux_phcm2s1(397) =    5.3  !   #HE I (C)
 lambda_EUV_A(398) =  489.00 ; solar_flux_phcm2s1(398) =    5.3  !   #HE I (C)
 lambda_EUV_A(399) =  489.50 ; solar_flux_phcm2s1(399) =   10.0  !   #NE III
 lambda_EUV_A(400) =  490.00 ; solar_flux_phcm2s1(400) =    5.6  !   #HE I (C)
 lambda_EUV_A(401) =  491.00 ; solar_flux_phcm2s1(401) =    5.9  !   #HE I (C)
 lambda_EUV_A(402) =  492.00 ; solar_flux_phcm2s1(402) =    6.3  !   #HE I (C)
 lambda_EUV_A(403) =  493.00 ; solar_flux_phcm2s1(403) =    6.6  !   #HE I (C)
 lambda_EUV_A(404) =  494.00 ; solar_flux_phcm2s1(404) =    6.9  !   #HE I (C)
 lambda_EUV_A(405) =  495.00 ; solar_flux_phcm2s1(405) =    7.6  !   #HE I (C)
 lambda_EUV_A(406) =  496.00 ; solar_flux_phcm2s1(406) =    7.9  !   #HE I (C)
 lambda_EUV_A(407) =  497.00 ; solar_flux_phcm2s1(407) =    8.3  !   #HE I (C)
 lambda_EUV_A(408) =  498.00 ; solar_flux_phcm2s1(408) =    8.6  !   #HE I (C)
 lambda_EUV_A(409) =  499.00 ; solar_flux_phcm2s1(409) =    9.2  !   #HE I (C)
 lambda_EUV_A(410) =  499.37 ; solar_flux_phcm2s1(410) =  100.0  !   #SI XII

! 500-550A, EUVAC bin 16
 lambda_EUV_A(411) =  500.00 ; solar_flux_phcm2s1(411) =    9.6  !   #HE I (C)
 lambda_EUV_A(412) =  501.00 ; solar_flux_phcm2s1(412) =   10.2  !   #HE I (C)
 lambda_EUV_A(413) =  502.00 ; solar_flux_phcm2s1(413) =   10.6  !   #HE I (C)
 lambda_EUV_A(414) =  503.00 ; solar_flux_phcm2s1(414) =   11.2  !   #HE I (C)
 lambda_EUV_A(415) =  504.00 ; solar_flux_phcm2s1(415) =   11.9  !   #HE I (C)
 lambda_EUV_A(416) =  507.93 ; solar_flux_phcm2s1(416) =  110.0  !   # O III  G
 lambda_EUV_A(417) =  515.60 ; solar_flux_phcm2s1(417) =   25.0  !   #HE I
 lambda_EUV_A(418) =  520.66 ; solar_flux_phcm2s1(418) =   46.0  !   #SI XII
 lambda_EUV_A(419) =  525.80 ; solar_flux_phcm2s1(419) =   65.0  !   # O III
 lambda_EUV_A(420) =  537.02 ; solar_flux_phcm2s1(420) =  120.0  !   #HE I

! 550-600A, EUVAC bin 19
 lambda_EUV_A(421) =  554.37 ; solar_flux_phcm2s1(421) =  720.0  !   # O IV   G   ! EUVAC bin/line 17

 lambda_EUV_A(422) =  558.60 ; solar_flux_phcm2s1(422) =   45.0  !   #NE VI
 lambda_EUV_A(423) =  562.80 ; solar_flux_phcm2s1(423) =   70.0  !   #NE VI

 lambda_EUV_A(424) =  584.33 ; solar_flux_phcm2s1(424) = 1270.0  !   #HE I        ! EUVAC bin/line 18

 lambda_EUV_A(425) =  599.60 ; solar_flux_phcm2s1(425) =  140.0  !   # O III

! 600-650A, EUVAC bin 22

 lambda_EUV_A(426) =  609.76 ; solar_flux_phcm2s1(426) =  530.0  !   #MG X        ! EUVAC bin/line 20

 lambda_EUV_A(427) =  616.60 ; solar_flux_phcm2s1(427) =   15.7  !   # O II
 lambda_EUV_A(428) =  624.00 ; solar_flux_phcm2s1(428) =    0.2  !   # H LY(C)
 lambda_EUV_A(429) =  624.93 ; solar_flux_phcm2s1(429) =  240.0  !   #MG X
 lambda_EUV_A(430) =  625.00 ; solar_flux_phcm2s1(430) =    0.2  !   # H LY(C)
 lambda_EUV_A(431) =  626.00 ; solar_flux_phcm2s1(431) =    0.2  !   # H LY(C)
 lambda_EUV_A(432) =  627.00 ; solar_flux_phcm2s1(432) =    0.2  !   # H LY(C)
 lambda_EUV_A(433) =  628.00 ; solar_flux_phcm2s1(433) =    0.2  !   # H LY(C)
 lambda_EUV_A(434) =  629.00 ; solar_flux_phcm2s1(434) =    0.2  !   # H LY(C)

 lambda_EUV_A(435) =  629.73 ; solar_flux_phcm2s1(435) = 1590.0  !   # O V        ! EUVAC bin/line 21

 lambda_EUV_A(436) =  630.00 ; solar_flux_phcm2s1(436) =    0.2  !   # H LY(C)
 lambda_EUV_A(437) =  640.41 ; solar_flux_phcm2s1(437) =   10.9  !   # S II
 lambda_EUV_A(438) =  640.93 ; solar_flux_phcm2s1(438) =   12.9  !   # S II
 lambda_EUV_A(439) =  641.81 ; solar_flux_phcm2s1(439) =   16.9  !   # S II
 lambda_EUV_A(440) =  644.00 ; solar_flux_phcm2s1(440) =    0.2  !   # H LY(C)
 lambda_EUV_A(441) =  645.00 ; solar_flux_phcm2s1(441) =    0.2  !   # H LY(C)
 lambda_EUV_A(442) =  646.00 ; solar_flux_phcm2s1(442) =    0.2  !   # H LY(C)
 lambda_EUV_A(443) =  647.00 ; solar_flux_phcm2s1(443) =    0.2  !   # H LY(C)
 lambda_EUV_A(444) =  648.00 ; solar_flux_phcm2s1(444) =    0.2  !   # H LY(C)
 lambda_EUV_A(445) =  649.00 ; solar_flux_phcm2s1(445) =    0.2  !   # H LY(C)

! 650-700A, EUVAC bin 23
 lambda_EUV_A(446) =  650.00 ; solar_flux_phcm2s1(446) =    0.2  !   # H LY(C)
 lambda_EUV_A(447) =  651.00 ; solar_flux_phcm2s1(447) =    0.2  !   # H LY(C)
 lambda_EUV_A(448) =  652.00 ; solar_flux_phcm2s1(448) =    0.2  !   # H LY(C)
 lambda_EUV_A(449) =  653.00 ; solar_flux_phcm2s1(449) =    0.2  !   # H LY(C)
 lambda_EUV_A(450) =  654.00 ; solar_flux_phcm2s1(450) =    0.2  !   # H LY(C)
 lambda_EUV_A(451) =  655.00 ; solar_flux_phcm2s1(451) =    0.2  !   # H LY(C)
 lambda_EUV_A(452) =  656.00 ; solar_flux_phcm2s1(452) =    0.2  !   # H LY(C)
 lambda_EUV_A(453) =  657.00 ; solar_flux_phcm2s1(453) =    0.3  !   # H LY(C)
 lambda_EUV_A(454) =  657.30 ; solar_flux_phcm2s1(454) =   11.0  !   # S IV
 lambda_EUV_A(455) =  658.00 ; solar_flux_phcm2s1(455) =    0.3  !   # H LY(C)
 lambda_EUV_A(456) =  659.00 ; solar_flux_phcm2s1(456) =    0.3  !   # H LY(C)
 lambda_EUV_A(457) =  660.00 ; solar_flux_phcm2s1(457) =    0.3  !   # H LY(C)
 lambda_EUV_A(458) =  661.00 ; solar_flux_phcm2s1(458) =    0.3  !   # H LY(C)
 lambda_EUV_A(459) =  661.40 ; solar_flux_phcm2s1(459) =   11.0  !   # S IV
 lambda_EUV_A(460) =  662.00 ; solar_flux_phcm2s1(460) =    0.3  !   # H LY(C)
 lambda_EUV_A(461) =  663.00 ; solar_flux_phcm2s1(461) =    0.3  !   # H LY(C)
 lambda_EUV_A(462) =  664.00 ; solar_flux_phcm2s1(462) =    0.3  !   # H LY(C)
 lambda_EUV_A(463) =  665.00 ; solar_flux_phcm2s1(463) =    0.3  !   # H LY(C)
 lambda_EUV_A(464) =  666.00 ; solar_flux_phcm2s1(464) =    0.3  !   # H LY(C)
 lambda_EUV_A(465) =  667.00 ; solar_flux_phcm2s1(465) =    0.4  !   # H LY(C)
 lambda_EUV_A(466) =  668.00 ; solar_flux_phcm2s1(466) =    0.4  !   # H LY(C)
 lambda_EUV_A(467) =  669.00 ; solar_flux_phcm2s1(467) =    0.4  !   # H LY(C)
 lambda_EUV_A(468) =  670.00 ; solar_flux_phcm2s1(468) =    0.4  !   # H LY(C)
 lambda_EUV_A(469) =  671.00 ; solar_flux_phcm2s1(469) =    0.4  !   # H LY(C)
 lambda_EUV_A(470) =  672.00 ; solar_flux_phcm2s1(470) =    0.4  !   # H LY(C)
 lambda_EUV_A(471) =  673.00 ; solar_flux_phcm2s1(471) =    0.4  !   # H LY(C)
 lambda_EUV_A(472) =  674.00 ; solar_flux_phcm2s1(472) =    0.4  !   # H LY(C)
 lambda_EUV_A(473) =  675.00 ; solar_flux_phcm2s1(473) =    0.5  !   # H LY(C)
 lambda_EUV_A(474) =  676.00 ; solar_flux_phcm2s1(474) =    0.5  !   # H LY(C)
 lambda_EUV_A(475) =  677.00 ; solar_flux_phcm2s1(475) =    0.5  !   # H LY(C)
 lambda_EUV_A(476) =  678.00 ; solar_flux_phcm2s1(476) =    0.5  !   # H LY(C)
 lambda_EUV_A(477) =  679.00 ; solar_flux_phcm2s1(477) =    0.5  !   # H LY(C)
 lambda_EUV_A(478) =  680.00 ; solar_flux_phcm2s1(478) =    0.5  !   # H LY(C)
 lambda_EUV_A(479) =  681.00 ; solar_flux_phcm2s1(479) =    0.5  !   # H LY(C)
 lambda_EUV_A(480) =  682.00 ; solar_flux_phcm2s1(480) =    0.5  !   # H LY(C)
 lambda_EUV_A(481) =  683.00 ; solar_flux_phcm2s1(481) =    0.5  !   # H LY(C)
 lambda_EUV_A(482) =  684.00 ; solar_flux_phcm2s1(482) =    0.5  !   # H LY(C)
 lambda_EUV_A(483) =  685.00 ; solar_flux_phcm2s1(483) =    0.5  !   # H LY(C)
 lambda_EUV_A(484) =  685.71 ; solar_flux_phcm2s1(484) =   90.0  !   # N III  G
 lambda_EUV_A(485) =  686.00 ; solar_flux_phcm2s1(485) =    0.5  !   # H LY(C)
 lambda_EUV_A(486) =  687.00 ; solar_flux_phcm2s1(486) =    0.6  !   # H LY(C)
 lambda_EUV_A(487) =  688.00 ; solar_flux_phcm2s1(487) =    0.6  !   # H LY(C)
 lambda_EUV_A(488) =  689.00 ; solar_flux_phcm2s1(488) =    0.6  !   # H LY(C)
 lambda_EUV_A(489) =  690.00 ; solar_flux_phcm2s1(489) =    0.6  !   # H LY(C)
 lambda_EUV_A(490) =  691.00 ; solar_flux_phcm2s1(490) =    0.6  !   # H LY(C)
 lambda_EUV_A(491) =  692.00 ; solar_flux_phcm2s1(491) =    0.7  !   # H LY(C)
 lambda_EUV_A(492) =  693.00 ; solar_flux_phcm2s1(492) =    0.7  !   # H LY(C)
 lambda_EUV_A(493) =  694.00 ; solar_flux_phcm2s1(493) =    0.7  !   # H LY(C)
 lambda_EUV_A(494) =  695.00 ; solar_flux_phcm2s1(494) =    0.7  !   # H LY(C)
 lambda_EUV_A(495) =  696.00 ; solar_flux_phcm2s1(495) =    0.7  !   # H LY(C)
 lambda_EUV_A(496) =  697.00 ; solar_flux_phcm2s1(496) =    0.7  !   # H LY(C)
 lambda_EUV_A(497) =  698.00 ; solar_flux_phcm2s1(497) =    0.7  !   # H LY(C)
 lambda_EUV_A(498) =  699.00 ; solar_flux_phcm2s1(498) =    0.7  !   # H LY(C)

! 700-750A, EUVAC bin 25
 lambda_EUV_A(499) =  700.00 ; solar_flux_phcm2s1(499) =    0.8  !   # H LY(C)
 lambda_EUV_A(500) =  701.00 ; solar_flux_phcm2s1(500) =    0.8  !   # H LY(C)
 lambda_EUV_A(501) =  702.00 ; solar_flux_phcm2s1(501) =    0.8  !   # H LY(C)
 lambda_EUV_A(502) =  703.00 ; solar_flux_phcm2s1(502) =    0.9  !   # H LY(C)

 lambda_EUV_A(503) =  703.36 ; solar_flux_phcm2s1(503) =  360.0  !   # O III  G   ! EUVAC bin/line 24

 lambda_EUV_A(504) =  704.00 ; solar_flux_phcm2s1(504) =    0.9  !   # H LY(C)
 lambda_EUV_A(505) =  705.00 ; solar_flux_phcm2s1(505) =    0.9  !   # H LY(C)
 lambda_EUV_A(506) =  706.00 ; solar_flux_phcm2s1(506) =    0.9  !   # H LY(C)
 lambda_EUV_A(507) =  707.00 ; solar_flux_phcm2s1(507) =    0.9  !   # H LY(C)
 lambda_EUV_A(508) =  708.00 ; solar_flux_phcm2s1(508) =    0.9  !   # H LY(C)
 lambda_EUV_A(509) =  709.00 ; solar_flux_phcm2s1(509) =    0.9  !   # H LY(C)
 lambda_EUV_A(510) =  710.00 ; solar_flux_phcm2s1(510) =    1.0  !   # H LY(C)
 lambda_EUV_A(511) =  711.00 ; solar_flux_phcm2s1(511) =    1.0  !   # H LY(C)
 lambda_EUV_A(512) =  712.00 ; solar_flux_phcm2s1(512) =    1.0  !   # H LY(C)
 lambda_EUV_A(513) =  713.00 ; solar_flux_phcm2s1(513) =    1.1  !   # H LY(C)
 lambda_EUV_A(514) =  714.00 ; solar_flux_phcm2s1(514) =    1.1  !   # H LY(C)
 lambda_EUV_A(515) =  715.00 ; solar_flux_phcm2s1(515) =    1.2  !   # H LY(C)
 lambda_EUV_A(516) =  716.00 ; solar_flux_phcm2s1(516) =    1.2  !   # H LY(C)
 lambda_EUV_A(517) =  717.00 ; solar_flux_phcm2s1(517) =    1.2  !   # H LY(C)
 lambda_EUV_A(518) =  718.00 ; solar_flux_phcm2s1(518) =    1.3  !   # H LY(C)
 lambda_EUV_A(519) =  718.50 ; solar_flux_phcm2s1(519) =   50.0  !   # O II
 lambda_EUV_A(520) =  719.00 ; solar_flux_phcm2s1(520) =    1.3  !   # H LY(C)
 lambda_EUV_A(521) =  720.00 ; solar_flux_phcm2s1(521) =    1.3  !   # H LY(C)
 lambda_EUV_A(522) =  721.00 ; solar_flux_phcm2s1(522) =    1.3  !   # H LY(C)
 lambda_EUV_A(523) =  722.00 ; solar_flux_phcm2s1(523) =    1.4  !   # H LY(C)
 lambda_EUV_A(524) =  723.00 ; solar_flux_phcm2s1(524) =    1.4  !   # H LY(C)
 lambda_EUV_A(525) =  724.00 ; solar_flux_phcm2s1(525) =    1.5  !   # H LY(C)
 lambda_EUV_A(526) =  725.00 ; solar_flux_phcm2s1(526) =    1.5  !   # H LY(C)
 lambda_EUV_A(527) =  726.00 ; solar_flux_phcm2s1(527) =    1.5  !   # H LY(C)
 lambda_EUV_A(528) =  727.00 ; solar_flux_phcm2s1(528) =    1.5  !   # H LY(C)
 lambda_EUV_A(529) =  728.00 ; solar_flux_phcm2s1(529) =    1.6  !   # H LY(C)
 lambda_EUV_A(530) =  729.00 ; solar_flux_phcm2s1(530) =    1.6  !   # H LY(C)
 lambda_EUV_A(531) =  730.00 ; solar_flux_phcm2s1(531) =    1.7  !   # H LY(C)
 lambda_EUV_A(532) =  731.00 ; solar_flux_phcm2s1(532) =    1.7  !   # H LY(C)
 lambda_EUV_A(533) =  732.00 ; solar_flux_phcm2s1(533) =    1.7  !   # H LY(C)
 lambda_EUV_A(534) =  733.00 ; solar_flux_phcm2s1(534) =    1.9  !   # H LY(C)
 lambda_EUV_A(535) =  734.00 ; solar_flux_phcm2s1(535) =    1.9  !   # H LY(C)
 lambda_EUV_A(536) =  735.00 ; solar_flux_phcm2s1(536) =    1.9  !   # H LY(C)
 lambda_EUV_A(537) =  736.00 ; solar_flux_phcm2s1(537) =    2.0  !   # H LY(C)
 lambda_EUV_A(538) =  737.00 ; solar_flux_phcm2s1(538) =    2.0  !   # H LY(C)
 lambda_EUV_A(539) =  738.00 ; solar_flux_phcm2s1(539) =    2.1  !   # H LY(C)
 lambda_EUV_A(540) =  739.00 ; solar_flux_phcm2s1(540) =    2.1  !   # H LY(C)
 lambda_EUV_A(541) =  740.00 ; solar_flux_phcm2s1(541) =    2.1  !   # H LY(C)
 lambda_EUV_A(542) =  741.00 ; solar_flux_phcm2s1(542) =    2.2  !   # H LY(C)
 lambda_EUV_A(543) =  742.00 ; solar_flux_phcm2s1(543) =    2.3  !   # H LY(C)
 lambda_EUV_A(544) =  743.00 ; solar_flux_phcm2s1(544) =    2.3  !   # H LY(C)
 lambda_EUV_A(545) =  744.00 ; solar_flux_phcm2s1(545) =    2.4  !   # H LY(C)
 lambda_EUV_A(546) =  745.00 ; solar_flux_phcm2s1(546) =    2.4  !   # H LY(C)
 lambda_EUV_A(547) =  746.00 ; solar_flux_phcm2s1(547) =    2.5  !   # H LY(C)
 lambda_EUV_A(548) =  747.00 ; solar_flux_phcm2s1(548) =    2.6  !   # H LY(C)
 lambda_EUV_A(549) =  748.00 ; solar_flux_phcm2s1(549) =    2.7  !   # H LY(C)
 lambda_EUV_A(550) =  749.00 ; solar_flux_phcm2s1(550) =    2.7  !   # H LY(C)

! 750-800A, EUVAC bin 29
 lambda_EUV_A(551) =  750.00 ; solar_flux_phcm2s1(551) =    2.8  !   # H LY(C)
 lambda_EUV_A(552) =  751.00 ; solar_flux_phcm2s1(552) =    2.9  !   # H LY(C)
 lambda_EUV_A(553) =  752.00 ; solar_flux_phcm2s1(553) =    2.9  !   # H LY(C)
 lambda_EUV_A(554) =  753.00 ; solar_flux_phcm2s1(554) =    3.0  !   # H LY(C)
 lambda_EUV_A(555) =  754.00 ; solar_flux_phcm2s1(555) =    3.1  !   # H LY(C)
 lambda_EUV_A(556) =  755.00 ; solar_flux_phcm2s1(556) =    3.1  !   # H LY(C)
 lambda_EUV_A(557) =  756.00 ; solar_flux_phcm2s1(557) =    3.2  !   # H LY(C)
 lambda_EUV_A(558) =  757.00 ; solar_flux_phcm2s1(558) =    3.3  !   # H LY(C)
 lambda_EUV_A(559) =  758.00 ; solar_flux_phcm2s1(559) =    3.4  !   # H LY(C)
 lambda_EUV_A(560) =  758.68 ; solar_flux_phcm2s1(560) =   30.0  !   # O V
 lambda_EUV_A(561) =  759.00 ; solar_flux_phcm2s1(561) =    3.6  !   # H LY(C)
 lambda_EUV_A(562) =  759.44 ; solar_flux_phcm2s1(562) =   23.0  !   # O V
 lambda_EUV_A(563) =  760.00 ; solar_flux_phcm2s1(563) =    3.6  !   # H LY(C)
 lambda_EUV_A(564) =  760.30 ; solar_flux_phcm2s1(564) =   80.0  !   # O V
 lambda_EUV_A(565) =  761.00 ; solar_flux_phcm2s1(565) =    3.7  !   # H LY(C)
 lambda_EUV_A(566) =  761.13 ; solar_flux_phcm2s1(566) =   20.0  !   # O V
 lambda_EUV_A(567) =  762.00 ; solar_flux_phcm2s1(567) =    3.8  !   # H LY(C)
 lambda_EUV_A(568) =  762.00 ; solar_flux_phcm2s1(568) =   30.0  !   # O V
 lambda_EUV_A(569) =  763.00 ; solar_flux_phcm2s1(569) =    3.9  !   # H LY(C)
 lambda_EUV_A(570) =  764.00 ; solar_flux_phcm2s1(570) =    4.0  !   # H LY(C)
 lambda_EUV_A(571) =  765.00 ; solar_flux_phcm2s1(571) =    4.1  !   # H LY(C)

 lambda_EUV_A(572) =  765.15 ; solar_flux_phcm2s1(572) =  170.0  !   # N IV         ! EUVAC bin/line 26

 lambda_EUV_A(573) =  766.00 ; solar_flux_phcm2s1(573) =    4.2  !   # H LY(C)
 lambda_EUV_A(574) =  767.00 ; solar_flux_phcm2s1(574) =    4.3  !   # H LY(C)
 lambda_EUV_A(575) =  768.00 ; solar_flux_phcm2s1(575) =    4.4  !   # H LY(C)
 lambda_EUV_A(576) =  769.00 ; solar_flux_phcm2s1(576) =    4.5  !   # H LY(C)
 lambda_EUV_A(577) =  770.00 ; solar_flux_phcm2s1(577) =    4.6  !   # H LY(C)

 lambda_EUV_A(578) =  770.41 ; solar_flux_phcm2s1(578) =  260.0  !   #NE VIII       ! EUVAC bin/line 27(1) ???

 lambda_EUV_A(579) =  771.00 ; solar_flux_phcm2s1(579) =    4.8  !   # H LY(C)
 lambda_EUV_A(580) =  772.00 ; solar_flux_phcm2s1(580) =    4.9  !   # H LY(C)
 lambda_EUV_A(581) =  773.00 ; solar_flux_phcm2s1(581) =    5.0  !   # H LY(C)
 lambda_EUV_A(582) =  774.00 ; solar_flux_phcm2s1(582) =    5.2  !   # H LY(C)
 lambda_EUV_A(583) =  775.00 ; solar_flux_phcm2s1(583) =    5.2  !   # H LY(C)
 lambda_EUV_A(584) =  776.00 ; solar_flux_phcm2s1(584) =   11.0  !   # N II
 lambda_EUV_A(585) =  776.00 ; solar_flux_phcm2s1(585) =    5.4  !   # H LY(C)
 lambda_EUV_A(586) =  777.00 ; solar_flux_phcm2s1(586) =    5.6  !   # H LY(C)
 lambda_EUV_A(587) =  778.00 ; solar_flux_phcm2s1(587) =    5.7  !   # H LY(C)
 lambda_EUV_A(588) =  779.00 ; solar_flux_phcm2s1(588) =    5.8  !   # H LY(C)
 lambda_EUV_A(589) =  780.00 ; solar_flux_phcm2s1(589) =    6.0  !   # H LY(C)

 lambda_EUV_A(590) =  780.32 ; solar_flux_phcm2s1(590) =  140.0  !   #NE VIII       ! EUVAC bin/line 27(2) ???

 lambda_EUV_A(591) =  781.00 ; solar_flux_phcm2s1(591) =    6.1  !   # H LY(C)
 lambda_EUV_A(592) =  782.00 ; solar_flux_phcm2s1(592) =    6.3  !   # H LY(C)
 lambda_EUV_A(593) =  783.00 ; solar_flux_phcm2s1(593) =    6.4  !   # H LY(C)
 lambda_EUV_A(594) =  784.00 ; solar_flux_phcm2s1(594) =    6.6  !   # H LY(C)
 lambda_EUV_A(595) =  785.00 ; solar_flux_phcm2s1(595) =    6.8  !   # H LY(C)
 lambda_EUV_A(596) =  786.00 ; solar_flux_phcm2s1(596) =    6.9  !   # H LY(C)
 lambda_EUV_A(597) =  786.47 ; solar_flux_phcm2s1(597) =  130.0  !   # S V
 lambda_EUV_A(598) =  787.00 ; solar_flux_phcm2s1(598) =    7.2  !   # H LY(C)

 lambda_EUV_A(599) =  787.71 ; solar_flux_phcm2s1(599) =  250.0  !   # O IV         ! EUVAC bin/line 28(1) ???

 lambda_EUV_A(600) =  788.00 ; solar_flux_phcm2s1(600) =    7.3  !   # H LY(C)
 lambda_EUV_A(601) =  789.00 ; solar_flux_phcm2s1(601) =    7.5  !   # H LY(C)
 lambda_EUV_A(602) =  790.00 ; solar_flux_phcm2s1(602) =    7.6  !   # H LY(C)

 lambda_EUV_A(603) =  790.15 ; solar_flux_phcm2s1(603) =  430.0  !   # O IV   G     ! EUVAC bin/line 28(2) ???

 lambda_EUV_A(604) =  791.00 ; solar_flux_phcm2s1(604) =    7.9  !   # H LY(C)
 lambda_EUV_A(605) =  792.00 ; solar_flux_phcm2s1(605) =    8.1  !   # H LY(C)
 lambda_EUV_A(606) =  793.00 ; solar_flux_phcm2s1(606) =    8.2  !   # H LY(C)
 lambda_EUV_A(607) =  794.00 ; solar_flux_phcm2s1(607) =    8.5  !   # H LY(C)
 lambda_EUV_A(608) =  795.00 ; solar_flux_phcm2s1(608) =    8.7  !   # H LY(C)
 lambda_EUV_A(609) =  796.00 ; solar_flux_phcm2s1(609) =    8.9  !   # H LY(C)
 lambda_EUV_A(610) =  797.00 ; solar_flux_phcm2s1(610) =    9.1  !   # H LY(C)
 lambda_EUV_A(611) =  798.00 ; solar_flux_phcm2s1(611) =    9.4  !   # H LY(C)
 lambda_EUV_A(612) =  799.00 ; solar_flux_phcm2s1(612) =    9.6  !   # H LY(C)

! 800-850A, EUVAC bin 30
 lambda_EUV_A(613) =  800.00 ; solar_flux_phcm2s1(613) =    9.8  !   # H LY(C)
 lambda_EUV_A(614) =  801.00 ; solar_flux_phcm2s1(614) =   10.2  !   # H LY(C)
 lambda_EUV_A(615) =  802.00 ; solar_flux_phcm2s1(615) =   10.4  !   # H LY(C)
 lambda_EUV_A(616) =  803.00 ; solar_flux_phcm2s1(616) =   10.7  !   # H LY(C)
 lambda_EUV_A(617) =  804.00 ; solar_flux_phcm2s1(617) =   11.0  !   # H LY(C)
 lambda_EUV_A(618) =  805.00 ; solar_flux_phcm2s1(618) =   11.2  !   # H LY(C)
 lambda_EUV_A(619) =  806.00 ; solar_flux_phcm2s1(619) =   11.5  !   # H LY(C)
 lambda_EUV_A(620) =  807.00 ; solar_flux_phcm2s1(620) =   11.8  !   # H LY(C)
 lambda_EUV_A(621) =  808.00 ; solar_flux_phcm2s1(621) =   12.1  !   # H LY(C)
 lambda_EUV_A(622) =  809.00 ; solar_flux_phcm2s1(622) =   12.4  !   # H LY(C)
 lambda_EUV_A(623) =  810.00 ; solar_flux_phcm2s1(623) =   12.7  !   # H LY(C)
 lambda_EUV_A(624) =  811.00 ; solar_flux_phcm2s1(624) =   13.1  !   # H LY(C)
 lambda_EUV_A(625) =  812.00 ; solar_flux_phcm2s1(625) =   13.4  !   # H LY(C)
 lambda_EUV_A(626) =  813.00 ; solar_flux_phcm2s1(626) =   13.7  !   # H LY(C)
 lambda_EUV_A(627) =  814.00 ; solar_flux_phcm2s1(627) =   14.1  !   # H LY(C)
 lambda_EUV_A(628) =  815.00 ; solar_flux_phcm2s1(628) =   14.5  !   # H LY(C)
 lambda_EUV_A(629) =  816.00 ; solar_flux_phcm2s1(629) =   14.8  !   # H LY(C)
 lambda_EUV_A(630) =  817.00 ; solar_flux_phcm2s1(630) =   15.2  !   # H LY(C)
 lambda_EUV_A(631) =  818.00 ; solar_flux_phcm2s1(631) =   15.5  !   # H LY(C)
 lambda_EUV_A(632) =  819.00 ; solar_flux_phcm2s1(632) =   16.0  !   # H LY(C)
 lambda_EUV_A(633) =  820.00 ; solar_flux_phcm2s1(633) =   16.3  !   # H LY(C)
 lambda_EUV_A(634) =  821.00 ; solar_flux_phcm2s1(634) =   16.8  !   # H LY(C)
 lambda_EUV_A(635) =  822.00 ; solar_flux_phcm2s1(635) =   17.2  !   # H LY(C)
 lambda_EUV_A(636) =  823.00 ; solar_flux_phcm2s1(636) =   17.7  !   # H LY(C)
 lambda_EUV_A(637) =  824.00 ; solar_flux_phcm2s1(637) =   18.2  !   # H LY(C)
 lambda_EUV_A(638) =  825.00 ; solar_flux_phcm2s1(638) =   18.6  !   # H LY(C)
 lambda_EUV_A(639) =  826.00 ; solar_flux_phcm2s1(639) =   19.1  !   # H LY(C)
 lambda_EUV_A(640) =  827.00 ; solar_flux_phcm2s1(640) =   19.6  !   # H LY(C)
 lambda_EUV_A(641) =  828.00 ; solar_flux_phcm2s1(641) =   20.0  !   # H LY(C)
 lambda_EUV_A(642) =  829.00 ; solar_flux_phcm2s1(642) =   20.6  !   # H LY(C)
 lambda_EUV_A(643) =  830.00 ; solar_flux_phcm2s1(643) =   21.1  !   # H LY(C)
 lambda_EUV_A(644) =  831.00 ; solar_flux_phcm2s1(644) =   21.6  !   # H LY(C)
 lambda_EUV_A(645) =  832.00 ; solar_flux_phcm2s1(645) =   22.2  !   # H LY(C)
 lambda_EUV_A(646) =  833.00 ; solar_flux_phcm2s1(646) =   22.8  !   # H LY(C)
 lambda_EUV_A(647) =  834.00 ; solar_flux_phcm2s1(647) =   23.3  !   # H LY(C)
 lambda_EUV_A(648) =  834.20 ; solar_flux_phcm2s1(648) =  620.0  !   # OII,IIIG
 lambda_EUV_A(649) =  835.00 ; solar_flux_phcm2s1(649) =   24.0  !   # H LY(C)
 lambda_EUV_A(650) =  836.00 ; solar_flux_phcm2s1(650) =   24.5  !   # H LY(C)
 lambda_EUV_A(651) =  837.00 ; solar_flux_phcm2s1(651) =   25.1  !   # H LY(C)
 lambda_EUV_A(652) =  838.00 ; solar_flux_phcm2s1(652) =   25.8  !   # H LY(C)
 lambda_EUV_A(653) =  839.00 ; solar_flux_phcm2s1(653) =   26.5  !   # H LY(C)
 lambda_EUV_A(654) =  840.00 ; solar_flux_phcm2s1(654) =   27.1  !   # H LY(C)
 lambda_EUV_A(655) =  841.00 ; solar_flux_phcm2s1(655) =   27.9  !   # H LY(C)
 lambda_EUV_A(656) =  842.00 ; solar_flux_phcm2s1(656) =   28.5  !   # H LY(C)
 lambda_EUV_A(657) =  843.00 ; solar_flux_phcm2s1(657) =   29.3  !   # H LY(C)
 lambda_EUV_A(658) =  844.00 ; solar_flux_phcm2s1(658) =   30.0  !   # H LY(C)
 lambda_EUV_A(659) =  845.00 ; solar_flux_phcm2s1(659) =   30.8  !   # H LY(C)
 lambda_EUV_A(660) =  846.00 ; solar_flux_phcm2s1(660) =   31.6  !   # H LY(C)
 lambda_EUV_A(661) =  847.00 ; solar_flux_phcm2s1(661) =   32.4  !   # H LY(C)
 lambda_EUV_A(662) =  848.00 ; solar_flux_phcm2s1(662) =   33.2  !   # H LY(C)
 lambda_EUV_A(663) =  849.00 ; solar_flux_phcm2s1(663) =   34.0  !   # H LY(C)

! 850-900A, EUVAC bin 31
 lambda_EUV_A(664) =  850.00 ; solar_flux_phcm2s1(664) =   34.9  !   # H LY(C)
 lambda_EUV_A(665) =  851.00 ; solar_flux_phcm2s1(665) =   35.8  !   # H LY(C)
 lambda_EUV_A(666) =  852.00 ; solar_flux_phcm2s1(666) =   36.7  !   # H LY(C)
 lambda_EUV_A(667) =  853.00 ; solar_flux_phcm2s1(667) =   37.7  !   # H LY(C)
 lambda_EUV_A(668) =  854.00 ; solar_flux_phcm2s1(668) =   38.6  !   # H LY(C)
 lambda_EUV_A(669) =  855.00 ; solar_flux_phcm2s1(669) =   39.6  !   # H LY(C)
 lambda_EUV_A(670) =  856.00 ; solar_flux_phcm2s1(670) =   40.6  !   # H LY(C)
 lambda_EUV_A(671) =  857.00 ; solar_flux_phcm2s1(671) =   41.7  !   # H LY(C)
 lambda_EUV_A(672) =  858.00 ; solar_flux_phcm2s1(672) =   42.7  !   # H LY(C)
 lambda_EUV_A(673) =  859.00 ; solar_flux_phcm2s1(673) =   43.8  !   # H LY(C)
 lambda_EUV_A(674) =  860.00 ; solar_flux_phcm2s1(674) =   44.9  !   # H LY(C)
 lambda_EUV_A(675) =  861.00 ; solar_flux_phcm2s1(675) =   46.1  !   # H LY(C)
 lambda_EUV_A(676) =  862.00 ; solar_flux_phcm2s1(676) =   47.2  !   # H LY(C)
 lambda_EUV_A(677) =  863.00 ; solar_flux_phcm2s1(677) =   48.4  !   # H LY(C)
 lambda_EUV_A(678) =  864.00 ; solar_flux_phcm2s1(678) =   49.6  !   # H LY(C)
 lambda_EUV_A(679) =  865.00 ; solar_flux_phcm2s1(679) =   50.9  !   # H LY(C)
 lambda_EUV_A(680) =  866.00 ; solar_flux_phcm2s1(680) =   52.2  !   # H LY(C)
 lambda_EUV_A(681) =  867.00 ; solar_flux_phcm2s1(681) =   53.5  !   # H LY(C)
 lambda_EUV_A(682) =  868.00 ; solar_flux_phcm2s1(682) =   54.9  !   # H LY(C)
 lambda_EUV_A(683) =  869.00 ; solar_flux_phcm2s1(683) =   56.3  !   # H LY(C)
 lambda_EUV_A(684) =  870.00 ; solar_flux_phcm2s1(684) =   57.7  !   # H LY(C)
 lambda_EUV_A(685) =  871.00 ; solar_flux_phcm2s1(685) =   59.2  !   # H LY(C)
 lambda_EUV_A(686) =  872.00 ; solar_flux_phcm2s1(686) =   60.7  !   # H LY(C)
 lambda_EUV_A(687) =  873.00 ; solar_flux_phcm2s1(687) =   62.3  !   # H LY(C)
 lambda_EUV_A(688) =  874.00 ; solar_flux_phcm2s1(688) =   63.9  !   # H LY(C)
 lambda_EUV_A(689) =  875.00 ; solar_flux_phcm2s1(689) =   65.5  !   # H LY(C)
 lambda_EUV_A(690) =  876.00 ; solar_flux_phcm2s1(690) =   67.1  !   # H LY(C)
 lambda_EUV_A(691) =  877.00 ; solar_flux_phcm2s1(691) =   68.8  !   # H LY(C)
 lambda_EUV_A(692) =  878.00 ; solar_flux_phcm2s1(692) =   70.6  !   # H LY(C)
 lambda_EUV_A(693) =  879.00 ; solar_flux_phcm2s1(693) =   72.4  !   # H LY(C)
 lambda_EUV_A(694) =  880.00 ; solar_flux_phcm2s1(694) =   74.2  !   # H LY(C)
 lambda_EUV_A(695) =  881.00 ; solar_flux_phcm2s1(695) =   76.1  !   # H LY(C)
 lambda_EUV_A(696) =  882.00 ; solar_flux_phcm2s1(696) =   78.1  !   # H LY(C)
 lambda_EUV_A(697) =  883.00 ; solar_flux_phcm2s1(697) =   80.1  !   # H LY(C)
 lambda_EUV_A(698) =  884.00 ; solar_flux_phcm2s1(698) =   82.1  !   # H LY(C)
 lambda_EUV_A(699) =  885.00 ; solar_flux_phcm2s1(699) =   84.2  !   # H LY(C)
 lambda_EUV_A(700) =  886.00 ; solar_flux_phcm2s1(700) =   86.3  !   # H LY(C)
 lambda_EUV_A(701) =  887.00 ; solar_flux_phcm2s1(701) =   88.5  !   # H LY(C)
 lambda_EUV_A(702) =  888.00 ; solar_flux_phcm2s1(702) =   90.8  !   # H LY(C)
 lambda_EUV_A(703) =  889.00 ; solar_flux_phcm2s1(703) =   93.1  !   # H LY(C)
 lambda_EUV_A(704) =  890.00 ; solar_flux_phcm2s1(704) =   95.5  !   # H LY(C)
 lambda_EUV_A(705) =  891.00 ; solar_flux_phcm2s1(705) =   97.9  !   # H LY(C)
 lambda_EUV_A(706) =  892.00 ; solar_flux_phcm2s1(706) =  100.4  !   # H LY(C)
 lambda_EUV_A(707) =  893.00 ; solar_flux_phcm2s1(707) =  102.9  !   # H LY(C)
 lambda_EUV_A(708) =  894.00 ; solar_flux_phcm2s1(708) =  105.5  !   # H LY(C)
 lambda_EUV_A(709) =  895.00 ; solar_flux_phcm2s1(709) =  108.2  !   # H LY(C)
 lambda_EUV_A(710) =  896.00 ; solar_flux_phcm2s1(710) =  111.0  !   # H LY(C)
 lambda_EUV_A(711) =  897.00 ; solar_flux_phcm2s1(711) =  113.8  !   # H LY(C)
 lambda_EUV_A(712) =  898.00 ; solar_flux_phcm2s1(712) =  116.6  !   # H LY(C)
 lambda_EUV_A(713) =  899.00 ; solar_flux_phcm2s1(713) =  119.6  !   # H LY(C)

! 900-950A, EUVAC bin 32
 lambda_EUV_A(714) =  900.00 ; solar_flux_phcm2s1(714) =  122.7  !   # H LY(C)
 lambda_EUV_A(715) =  901.00 ; solar_flux_phcm2s1(715) =  125.8  !   # H LY(C)
 lambda_EUV_A(716) =  902.00 ; solar_flux_phcm2s1(716) =  129.0  !   # H LY(C)
 lambda_EUV_A(717) =  903.00 ; solar_flux_phcm2s1(717) =  132.3  !   # H LY(C)
 lambda_EUV_A(718) =  904.00 ; solar_flux_phcm2s1(718) =  135.7  !   # H LY(C)
 lambda_EUV_A(719) =  904.10 ; solar_flux_phcm2s1(719) =  110.0  !   # C II   G
 lambda_EUV_A(720) =  905.00 ; solar_flux_phcm2s1(720) =  139.1  !   # H LY(C)
 lambda_EUV_A(721) =  906.00 ; solar_flux_phcm2s1(721) =  142.6  !   # H LY(C)
 lambda_EUV_A(722) =  907.00 ; solar_flux_phcm2s1(722) =  146.3  !   # H LY(C)
 lambda_EUV_A(723) =  908.00 ; solar_flux_phcm2s1(723) =  150.0  !   # H LY(C)
 lambda_EUV_A(724) =  909.00 ; solar_flux_phcm2s1(724) =  153.8  !   # H LY(C)
 lambda_EUV_A(725) =  910.00 ; solar_flux_phcm2s1(725) =  157.7  !   # H LY(C)
 lambda_EUV_A(726) =  911.00 ; solar_flux_phcm2s1(726) =  161.7  !   # H LY(C)
 lambda_EUV_A(727) =  912.00 ; solar_flux_phcm2s1(727) =  165.8  !   # H LY(C)
 lambda_EUV_A(728) =  913.00 ; solar_flux_phcm2s1(728) =    2.3  !   # C I (C)
 lambda_EUV_A(729) =  914.00 ; solar_flux_phcm2s1(729) =    2.3  !   # C I (C)
 lambda_EUV_A(730) =  915.00 ; solar_flux_phcm2s1(730) =    2.3  !   # C I (C)
 lambda_EUV_A(731) =  916.00 ; solar_flux_phcm2s1(731) =    2.4  !   # C I (C)
 lambda_EUV_A(732) =  917.00 ; solar_flux_phcm2s1(732) =    2.4  !   # C I (C)
 lambda_EUV_A(733) =  918.00 ; solar_flux_phcm2s1(733) =    2.5  !   # C I (C)
 lambda_EUV_A(734) =  919.00 ; solar_flux_phcm2s1(734) =    2.5  !   # C I (C)
 lambda_EUV_A(735) =  920.00 ; solar_flux_phcm2s1(735) =    2.6  !   # C I (C)
 lambda_EUV_A(736) =  920.96 ; solar_flux_phcm2s1(736) =   56.0  !   # H LY-9
 lambda_EUV_A(737) =  921.00 ; solar_flux_phcm2s1(737) =    2.6  !   # C I (C)
 lambda_EUV_A(738) =  922.00 ; solar_flux_phcm2s1(738) =    2.6  !   # C I (C)
 lambda_EUV_A(739) =  923.00 ; solar_flux_phcm2s1(739) =    2.8  !   # C I (C)
 lambda_EUV_A(740) =  923.15 ; solar_flux_phcm2s1(740) =   80.0  !   # H LY-8 B
 lambda_EUV_A(741) =  924.00 ; solar_flux_phcm2s1(741) =    2.8  !   # C I (C)
 lambda_EUV_A(742) =  925.00 ; solar_flux_phcm2s1(742) =    2.9  !   # C I (C)
 lambda_EUV_A(743) =  926.00 ; solar_flux_phcm2s1(743) =    2.9  !   # C I (C)
 lambda_EUV_A(744) =  926.20 ; solar_flux_phcm2s1(744) =  110.0  !   # H LY-7
 lambda_EUV_A(745) =  927.00 ; solar_flux_phcm2s1(745) =    3.0  !   # C I (C)
 lambda_EUV_A(746) =  928.00 ; solar_flux_phcm2s1(746) =    3.0  !   # C I (C)
 lambda_EUV_A(747) =  929.00 ; solar_flux_phcm2s1(747) =    3.1  !   # C I (C)
 lambda_EUV_A(748) =  930.00 ; solar_flux_phcm2s1(748) =    3.1  !   # C I (C)
 lambda_EUV_A(749) =  930.75 ; solar_flux_phcm2s1(749) =  130.0  !   # H LY-6 B
 lambda_EUV_A(750) =  931.00 ; solar_flux_phcm2s1(750) =    3.2  !   # C I (C)
 lambda_EUV_A(751) =  932.00 ; solar_flux_phcm2s1(751) =    3.2  !   # C I (C)
 lambda_EUV_A(752) =  933.00 ; solar_flux_phcm2s1(752) =    3.3  !   # C I (C)
 lambda_EUV_A(753) =  933.38 ; solar_flux_phcm2s1(753) =  103.0  !   # S VI
 lambda_EUV_A(754) =  934.00 ; solar_flux_phcm2s1(754) =    3.3  !   # C I (C)
 lambda_EUV_A(755) =  935.00 ; solar_flux_phcm2s1(755) =    3.4  !   # C I (C)
 lambda_EUV_A(756) =  936.00 ; solar_flux_phcm2s1(756) =    3.4  !   # C I (C)
 lambda_EUV_A(757) =  937.00 ; solar_flux_phcm2s1(757) =    3.5  !   # C I (C)
 lambda_EUV_A(758) =  937.80 ; solar_flux_phcm2s1(758) =  180.0  !   # H LY-5
 lambda_EUV_A(759) =  938.00 ; solar_flux_phcm2s1(759) =    3.6  !   # C I (C)
 lambda_EUV_A(760) =  939.00 ; solar_flux_phcm2s1(760) =    3.6  !   # C I (C)
 lambda_EUV_A(761) =  940.00 ; solar_flux_phcm2s1(761) =    3.7  !   # C I (C)
 lambda_EUV_A(762) =  941.00 ; solar_flux_phcm2s1(762) =    3.7  !   # C I (C)
 lambda_EUV_A(763) =  942.00 ; solar_flux_phcm2s1(763) =    3.9  !   # C I (C)
 lambda_EUV_A(764) =  943.00 ; solar_flux_phcm2s1(764) =    4.0  !   # C I (C)
 lambda_EUV_A(765) =  944.00 ; solar_flux_phcm2s1(765) =    4.0  !   # C I (C)
 lambda_EUV_A(766) =  944.52 ; solar_flux_phcm2s1(766) =   68.0  !   # S VI
 lambda_EUV_A(767) =  945.00 ; solar_flux_phcm2s1(767) =    4.1  !   # C I (C)
 lambda_EUV_A(768) =  946.00 ; solar_flux_phcm2s1(768) =    4.2  !   # C I (C)
 lambda_EUV_A(769) =  947.00 ; solar_flux_phcm2s1(769) =    4.2  !   # C I (C)
 lambda_EUV_A(770) =  948.00 ; solar_flux_phcm2s1(770) =    4.3  !   # C I (C)
 lambda_EUV_A(771) =  949.00 ; solar_flux_phcm2s1(771) =    4.4  !   # C I (C)
 lambda_EUV_A(772) =  949.74 ; solar_flux_phcm2s1(772) =  300.0  !   # H LY-4

! 950-1000A, EUVAC bin 34
 lambda_EUV_A(773) =  950.00 ; solar_flux_phcm2s1(773) =    4.4  !   # C I (C)
 lambda_EUV_A(774) =  951.00 ; solar_flux_phcm2s1(774) =    4.5  !   # C I (C)
 lambda_EUV_A(775) =  952.00 ; solar_flux_phcm2s1(775) =    4.6  !   # C I (C)
 lambda_EUV_A(776) =  953.00 ; solar_flux_phcm2s1(776) =    4.7  !   # C I (C)
 lambda_EUV_A(777) =  954.00 ; solar_flux_phcm2s1(777) =    4.7  !   # C I (C)
 lambda_EUV_A(778) =  955.00 ; solar_flux_phcm2s1(778) =    4.8  !   # C I (C)
 lambda_EUV_A(779) =  956.00 ; solar_flux_phcm2s1(779) =    5.0  !   # C I (C)
 lambda_EUV_A(780) =  957.00 ; solar_flux_phcm2s1(780) =    5.1  !   # C I (C)
 lambda_EUV_A(781) =  958.00 ; solar_flux_phcm2s1(781) =    5.2  !   # C I (C)
 lambda_EUV_A(782) =  959.00 ; solar_flux_phcm2s1(782) =    5.3  !   # C I (C)
 lambda_EUV_A(783) =  960.00 ; solar_flux_phcm2s1(783) =    5.3  !   # C I (C)
 lambda_EUV_A(784) =  961.00 ; solar_flux_phcm2s1(784) =    5.4  !   # C I (C)
 lambda_EUV_A(785) =  962.00 ; solar_flux_phcm2s1(785) =    5.5  !   # C I (C)
 lambda_EUV_A(786) =  963.00 ; solar_flux_phcm2s1(786) =    5.6  !   # C I (C)
 lambda_EUV_A(787) =  964.00 ; solar_flux_phcm2s1(787) =    5.7  !   # C I (C)
 lambda_EUV_A(788) =  965.00 ; solar_flux_phcm2s1(788) =    5.8  !   # C I (C)
 lambda_EUV_A(789) =  966.00 ; solar_flux_phcm2s1(789) =    5.9  !   # C I (C)
 lambda_EUV_A(790) =  967.00 ; solar_flux_phcm2s1(790) =    6.1  !   # C I (C)
 lambda_EUV_A(791) =  968.00 ; solar_flux_phcm2s1(791) =    6.2  !   # C I (C)
 lambda_EUV_A(792) =  969.00 ; solar_flux_phcm2s1(792) =    6.3  !   # C I (C)
 lambda_EUV_A(793) =  970.00 ; solar_flux_phcm2s1(793) =    6.4  !   # C I (C)
 lambda_EUV_A(794) =  971.00 ; solar_flux_phcm2s1(794) =    6.5  !   # C I (C)
 lambda_EUV_A(795) =  972.00 ; solar_flux_phcm2s1(795) =    6.6  !   # C I (C)
 lambda_EUV_A(796) =  972.54 ; solar_flux_phcm2s1(796) =  600.0  !   # H LY-3
 lambda_EUV_A(797) =  973.00 ; solar_flux_phcm2s1(797) =    6.7  !   # C I (C)
 lambda_EUV_A(798) =  974.00 ; solar_flux_phcm2s1(798) =    6.8  !   # C I (C)
 lambda_EUV_A(799) =  975.00 ; solar_flux_phcm2s1(799) =    7.0  !   # C I (C)
 lambda_EUV_A(800) =  976.00 ; solar_flux_phcm2s1(800) =    7.2  !   # C I (C)
 lambda_EUV_A(801) =  977.00 ; solar_flux_phcm2s1(801) =    7.3  !   # C I (C)

 lambda_EUV_A(802) =  977.02 ; solar_flux_phcm2s1(802) = 4400.0  !   # C III        ! EUVAC bin/line 33

 lambda_EUV_A(803) =  978.00 ; solar_flux_phcm2s1(803) =    7.4  !   # C I (C)
 lambda_EUV_A(804) =  979.00 ; solar_flux_phcm2s1(804) =    7.5  !   # C I (C)
 lambda_EUV_A(805) =  980.00 ; solar_flux_phcm2s1(805) =    7.7  !   # C I (C)
 lambda_EUV_A(806) =  981.00 ; solar_flux_phcm2s1(806) =    7.8  !   # C I (C)
 lambda_EUV_A(807) =  982.00 ; solar_flux_phcm2s1(807) =    7.9  !   # C I (C)
 lambda_EUV_A(808) =  983.00 ; solar_flux_phcm2s1(808) =    8.0  !   # C I (C)
 lambda_EUV_A(809) =  984.00 ; solar_flux_phcm2s1(809) =    8.3  !   # C I (C)
 lambda_EUV_A(810) =  985.00 ; solar_flux_phcm2s1(810) =    8.4  !   # C I (C)
 lambda_EUV_A(811) =  986.00 ; solar_flux_phcm2s1(811) =    8.5  !   # C I (C)
 lambda_EUV_A(812) =  987.00 ; solar_flux_phcm2s1(812) =    8.7  !   # C I (C)
 lambda_EUV_A(813) =  988.00 ; solar_flux_phcm2s1(813) =    8.8  !   # C I (C)
 lambda_EUV_A(814) =  989.00 ; solar_flux_phcm2s1(814) =    9.0  !   # C I (C)
 lambda_EUV_A(815) =  989.79 ; solar_flux_phcm2s1(815) =  170.0  !   # N III
 lambda_EUV_A(816) =  990.00 ; solar_flux_phcm2s1(816) =    9.1  !   # C I (C)
 lambda_EUV_A(817) =  991.00 ; solar_flux_phcm2s1(817) =    9.4  !   # C I (C)
 lambda_EUV_A(818) =  991.55 ; solar_flux_phcm2s1(818) =  340.0  !   # N III
 lambda_EUV_A(819) =  992.00 ; solar_flux_phcm2s1(819) =    9.5  !   # C I (C)
 lambda_EUV_A(820) =  993.00 ; solar_flux_phcm2s1(820) =    9.7  !   # C I (C)
 lambda_EUV_A(821) =  994.00 ; solar_flux_phcm2s1(821) =    9.8  !   # C I (C)
 lambda_EUV_A(822) =  995.00 ; solar_flux_phcm2s1(822) =   10.0  !   # C I (C)
 lambda_EUV_A(823) =  996.00 ; solar_flux_phcm2s1(823) =   10.2  !   # C I (C)
 lambda_EUV_A(824) =  997.00 ; solar_flux_phcm2s1(824) =   10.3  !   # C I (C)
 lambda_EUV_A(825) =  998.00 ; solar_flux_phcm2s1(825) =   10.6  !   # C I (C)
 lambda_EUV_A(826) =  999.00 ; solar_flux_phcm2s1(826) =   10.8  !   # C I (C)

! 1000-1050A, EUVAC bin 37 (actual range 1000-1025A)
 lambda_EUV_A(827) = 1000.00 ; solar_flux_phcm2s1(827) =   11.0  !   # C I (C)
 lambda_EUV_A(828) = 1001.00 ; solar_flux_phcm2s1(828) =   11.1  !   # C I (C)
 lambda_EUV_A(829) = 1002.00 ; solar_flux_phcm2s1(829) =   11.3  !   # C I (C)
 lambda_EUV_A(830) = 1003.00 ; solar_flux_phcm2s1(830) =   11.6  !   # C I (C)
 lambda_EUV_A(831) = 1004.00 ; solar_flux_phcm2s1(831) =   11.8  !   # C I (C)
 lambda_EUV_A(832) = 1005.00 ; solar_flux_phcm2s1(832) =   12.0  !   # C I (C)
 lambda_EUV_A(833) = 1006.00 ; solar_flux_phcm2s1(833) =   12.2  !   # C I (C)
 lambda_EUV_A(834) = 1007.00 ; solar_flux_phcm2s1(834) =   12.4  !   # C I (C)
 lambda_EUV_A(835) = 1008.00 ; solar_flux_phcm2s1(835) =   12.7  !   # C I (C)
 lambda_EUV_A(836) = 1009.00 ; solar_flux_phcm2s1(836) =   12.9  !   # C I (C)
 lambda_EUV_A(837) = 1010.00 ; solar_flux_phcm2s1(837) =   13.1  !   # C I (C)
 lambda_EUV_A(838) = 1010.20 ; solar_flux_phcm2s1(838) =   80.0  !   # C II
 lambda_EUV_A(839) = 1011.00 ; solar_flux_phcm2s1(839) =   13.3  !   # C I (C)
 lambda_EUV_A(840) = 1012.00 ; solar_flux_phcm2s1(840) =   13.6  !   # C I (C)
 lambda_EUV_A(841) = 1013.00 ; solar_flux_phcm2s1(841) =   13.9  !   # C I (C)
 lambda_EUV_A(842) = 1014.00 ; solar_flux_phcm2s1(842) =   14.1  !   # C I (C)
 lambda_EUV_A(843) = 1015.00 ; solar_flux_phcm2s1(843) =   14.4  !   # C I (C)
 lambda_EUV_A(844) = 1016.00 ; solar_flux_phcm2s1(844) =   14.6  !   # C I (C)
 lambda_EUV_A(845) = 1017.00 ; solar_flux_phcm2s1(845) =   14.9  !   # C I (C)
 lambda_EUV_A(846) = 1018.00 ; solar_flux_phcm2s1(846) =   15.2  !   # C I (C)
 lambda_EUV_A(847) = 1019.00 ; solar_flux_phcm2s1(847) =   15.4  !   # C I (C)
 lambda_EUV_A(848) = 1020.00 ; solar_flux_phcm2s1(848) =   15.7  !   # C I (C)
 lambda_EUV_A(849) = 1021.00 ; solar_flux_phcm2s1(849) =   16.0  !   # C I (C)
 lambda_EUV_A(850) = 1022.00 ; solar_flux_phcm2s1(850) =   16.3  !   # C I (C)
 lambda_EUV_A(851) = 1023.00 ; solar_flux_phcm2s1(851) =   16.6  !   # C I (C)
 lambda_EUV_A(852) = 1024.00 ; solar_flux_phcm2s1(852) =   16.8  !   # C I (C)
 lambda_EUV_A(853) = 1025.00 ; solar_flux_phcm2s1(853) =   17.2  !   # C I (C)

! EUVAC lines 35 (1025.72 A or 12.0874 eV) and 36 (1031.91 A or 12.015 eV) have not enough energy for ionization:
! the lowest photoionization threshold is 12.0876 eV for O2+ from O2 
! this is why these lines (and all wavelengths longer than 1025 A or 12.096 eV) are omitted

! the flux above is in units of [10^10 1/m^2/s]
! to convert it to units of [10^9 1/cm^2/s] use 
! 10^10 1/m^2/s = 10*10^9 1/(10^4 cm^2)/s = 0.001*(10^9 1/cm^2/s) ::

 solar_flux_phcm2s1 = 0.001 * solar_flux_phcm2s1  ! now the flux is in units of [10^9 1/cm^2/s]

 solar_flux_phcm2s1(  1:179) = 3.0 * solar_flux_phcm2s1(  1:179) ! below 150A
 solar_flux_phcm2s1(180:272) = 2.0 * solar_flux_phcm2s1(180:272) ! between 150A and 250 A

! correct the solar flux due to the solar activity variation same way it is done in the EUVAC model
! use Eqs.(9.19) and (9.20) from the Schunk and Nagy's book

 factor_f107 = 0.5 * (f10p7 + f10p7_81) - 80.0

! below 50A  (same as 50-100A)
 n1 = 1
 n2 = 15
 scale_factor = 1.0 + 1.0017e-2 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 50-100A, EUVAC bin 1
 n1 = 16
 n2 = 144
 scale_factor = 1.0 + 1.0017e-2 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 100-150A, EUVAC bin 2
 n1 = 145
 n2 = 179
 scale_factor = 1.0 + 7.1250e-3 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 150-200A, EUVAC bin 3
 n1 = 180
 n2 = 224
 scale_factor = 1.0 + 1.3375e-2 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 200-250A, EUVAC bin 4
 n1 = 225
 n2 = 272
 scale_factor = 1.0 + 1.9450e-2 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 250-300A, EUVAC bin 7
 n1 = 273
 n2 = 307
 DO i = n1, n2
    SELECT CASE (i)
       CASE (277)
! 256.32A, EUVAC bin/line 5
          scale_factor = 1.0 + 2.7750e-3 * factor_f107
       CASE (300)
! 284.15A, EUVAC bin/line 6
          scale_factor = 1.0 + 1.3768e-1 * factor_f107
       CASE DEFAULT
! bin 7
          scale_factor = 1.0 + 2.6467e-2 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 300-350A, EUVAC bin 10
 n1 = 308
 n2 = 319
 DO i = n1, n2
    SELECT CASE (i)
       CASE (308)
! 303.31A, EUVAC bin/line 8
          scale_factor = 1.0 + 2.5000e-2 * factor_f107
       CASE (309)
! 303.78A, EUVAC bin/line 9
          scale_factor = 1.0 + 3.3333e-3 * factor_f107
       CASE DEFAULT
! bin 10
          scale_factor = 1.0 + 2.2450e-2 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 350-400A, EUVAC bin 12
 n1 = 320
 n2 = 324
 DO i = n1, n2
    SELECT CASE (i)
       CASE (323)
! 368.07A, EUVAC bin/line 11
          scale_factor = 1.0 + 6.5917e-3 * factor_f107
       CASE DEFAULT
! bin 12
          scale_factor = 1.0 + 3.6542e-2 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 400-450A, EUVAC bin 13
 n1 = 325
 n2 = 357
 scale_factor = 1.0 + 7.4083e-3 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 450-500A, EUVAC bin 15
 n1 = 358
 n2 = 410
 DO i = n1, n2
    SELECT CASE (i)
       CASE (374)
! 465.22A, EUVAC bin/line 14
          scale_factor = 1.0 + 7.4917e-3 * factor_f107
       CASE DEFAULT
! bin 15
          scale_factor = 1.0 + 2.0225e-2 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 500-550A, EUVAC bin 16
 n1 = 411
 n2 = 420
 scale_factor = 1.0 + 8.7583e-3 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 550-600A, EUVAC bin 19
 n1 = 421
 n2 = 425
 DO i = n1, n2
    SELECT CASE (i)
       CASE (421)
! 554.37A, EUVAC bin/line 17
          scale_factor = 1.0 + 3.2667e-3 * factor_f107
       CASE (424)
! 584.33A, EUVAC bin/line 18
          scale_factor = 1.0 + 5.1583e-3 * factor_f107
       CASE DEFAULT
! bin 19
          scale_factor = 1.0 + 3.6583e-3 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 600-650A, EUVAC bin 22
 n1 = 426
 n2 = 445
 DO i = n1, n2
    SELECT CASE (i)
       CASE (426)
! 609.76A, EUVAC bin/line 20
          scale_factor = 1.0 + 1.6175e-2 * factor_f107
       CASE (435)
! 629.73A, EUVAC bin/line 21
          scale_factor = 1.0 + 3.3250e-3 * factor_f107
       CASE DEFAULT
! bin 22
          scale_factor = 1.0 + 1.1800e-2 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 650-700A, EUVAC bin 23
 n1 = 446
 n2 = 498
 scale_factor = 1.0 + 4.2667e-3 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 700-750A, EUVAC bin 25
 n1 = 499
 n2 = 550
 DO i = n1, n2
    SELECT CASE (i)
       CASE (503)
! 703.36A, EUVAC bin/line 24
          scale_factor = 1.0 + 3.0417e-3 * factor_f107
       CASE DEFAULT
! bin 25
          scale_factor = 1.0 + 4.7500e-3 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 750-800A, EUVAC bin 29
 n1 = 551
 n2 = 612
 DO i = n1, n2
    SELECT CASE (i)
       CASE (572)
! 765.15A, EUVAC bin/line 26
          scale_factor = 1.0 + 3.8500e-3 * factor_f107
       CASE (578,590)
! 578 is 770.41A, EUVAC bin/line 27, identified as NE VIII
! 590 is intense line with wavelength 780.32, also identified as NE VIII but not mentioned in EUVAC
          scale_factor = 1.0 + 1.2808e-2 * factor_f107
       CASE (599,603)
! 789.36A, EUVAC bin/line 28, identified as O IV
! apparently, EUVAC wavelength and flux are averaged over multiplet [Heroux and Hinteregger, JGR, 83, 5305, 1978]
! 599 is line with wavelength 787.71 identified as O IV
! 603 is line with wavelength 790.15 identified as O IV
          scale_factor = 1.0 + 3.2750e-3 * factor_f107
       CASE DEFAULT
! bin 29
          scale_factor = 1.0 + 4.7667e-3 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 800-850A, EUVAC bin 30
 n1 = 613
 n2 = 663
 scale_factor = 1.0 + 4.8167e-3 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 850-900A, EUVAC bin 31
 n1 = 664
 n2 = 713
 scale_factor = 1.0 + 5.6750e-3 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 900-950A, EUVAC bin 32
 n1 = 714
 n2 = 772
 scale_factor = 1.0 + 4.9833e-3 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 950-1000A, EUVAC bin 34
 n1 = 773
 n2 = 826
 DO i = n1, n2
    SELECT CASE (i)
       CASE (802)
! 977.02A, EUVAC bin/line 33
          scale_factor = 1.0 + 3.9417e-3 * factor_f107
       CASE DEFAULT
! bin 34
          scale_factor = 1.0 + 4.4167e-3 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 1000-1050A, EUVAC bin 37 (actual range 1000-1026A)
 n1 = 827
 n2 = 853
 scale_factor = 1.0 + 4.3750e-3 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! solar flux energy bins

  factor_energy_eVA =  REAL(h_Planck_Js * c_ms * 1.0d10 / e_Cl)   ! 1.0d10 because the wavelength is in Angstrom

  DO i = 1, N_solar_bins
     solar_flux_energy_bin_eV(i) = factor_energy_eVA / lambda_EUV_A(i)
  END DO


END SUBROUTINE Use_f74113
