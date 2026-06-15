!---------------------------------------------------------------
! the wavelengths and fluxes are from file f76ref.dat downloaded from
!
! ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_UV/AE_E/REF_SPEC/f76ref.dat
!
! this link can be found here: https://www.ngdc.noaa.gov/stp/solar/solaruv.html
!
! The data were edited: zero fluxes and fluxes with wavelength shorter than 28.47 A and longer than 1026 A were removed
!
SUBROUTINE Use_f76(lambda_EUV_A)

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

! below 50A
 lambda_EUV_A(  1) =    28.47 ; solar_flux_phcm2s1(  1) =     .5
 lambda_EUV_A(  2) =    28.79 ; solar_flux_phcm2s1(  2) =    2.5
 lambda_EUV_A(  3) =    29.52 ; solar_flux_phcm2s1(  3) =    2.2
 lambda_EUV_A(  4) =    30.02 ; solar_flux_phcm2s1(  4) =     .9
 lambda_EUV_A(  5) =    30.43 ; solar_flux_phcm2s1(  5) =     .6
 lambda_EUV_A(  6) =    33.74 ; solar_flux_phcm2s1(  6) =    1.1
 lambda_EUV_A(  7) =    40.95 ; solar_flux_phcm2s1(  7) =     .6
 lambda_EUV_A(  8) =    43.76 ; solar_flux_phcm2s1(  8) =    2.1
 lambda_EUV_A(  9) =    44.02 ; solar_flux_phcm2s1(  9) =     .8
 lambda_EUV_A( 10) =    44.16 ; solar_flux_phcm2s1( 10) =     .9
 lambda_EUV_A( 11) =    45.66 ; solar_flux_phcm2s1( 11) =     .5
 lambda_EUV_A( 12) =    46.40 ; solar_flux_phcm2s1( 12) =    2.7
 lambda_EUV_A( 13) =    46.67 ; solar_flux_phcm2s1( 13) =    4.0
 lambda_EUV_A( 14) =    47.87 ; solar_flux_phcm2s1( 14) =    4.5
 lambda_EUV_A( 15) =    49.22 ; solar_flux_phcm2s1( 15) =    4.3

! 50-100A, EUVAC bin 1
 lambda_EUV_A( 16) =    50.52 ; solar_flux_phcm2s1( 16) =    5.6
 lambda_EUV_A( 17) =    50.69 ; solar_flux_phcm2s1( 17) =    5.6
 lambda_EUV_A( 18) =    52.30 ; solar_flux_phcm2s1( 18) =    2.8
 lambda_EUV_A( 19) =    52.91 ; solar_flux_phcm2s1( 19) =     .1
 lambda_EUV_A( 20) =    54.15 ; solar_flux_phcm2s1( 20) =    6.6
 lambda_EUV_A( 21) =    54.42 ; solar_flux_phcm2s1( 21) =    2.9
 lambda_EUV_A( 22) =    55.06 ; solar_flux_phcm2s1( 22) =    3.2
 lambda_EUV_A( 23) =    55.34 ; solar_flux_phcm2s1( 23) =    8.6
 lambda_EUV_A( 24) =    56.08 ; solar_flux_phcm2s1( 24) =    1.9
 lambda_EUV_A( 25) =    56.92 ; solar_flux_phcm2s1( 25) =    6.0
 lambda_EUV_A( 26) =    57.36 ; solar_flux_phcm2s1( 26) =    5.0
 lambda_EUV_A( 27) =    57.56 ; solar_flux_phcm2s1( 27) =    4.0
 lambda_EUV_A( 28) =    57.88 ; solar_flux_phcm2s1( 28) =    3.5
 lambda_EUV_A( 29) =    58.96 ; solar_flux_phcm2s1( 29) =     .6
 lambda_EUV_A( 30) =    59.62 ; solar_flux_phcm2s1( 30) =     .6
 lambda_EUV_A( 31) =    60.30 ; solar_flux_phcm2s1( 31) =    2.5
 lambda_EUV_A( 32) =    60.85 ; solar_flux_phcm2s1( 32) =    3.6
 lambda_EUV_A( 33) =    61.07 ; solar_flux_phcm2s1( 33) =    5.8
 lambda_EUV_A( 34) =    61.63 ; solar_flux_phcm2s1( 34) =    2.9
 lambda_EUV_A( 35) =    61.90 ; solar_flux_phcm2s1( 35) =    5.0
 lambda_EUV_A( 36) =    62.30 ; solar_flux_phcm2s1( 36) =     .1
 lambda_EUV_A( 37) =    62.35 ; solar_flux_phcm2s1( 37) =    1.1
 lambda_EUV_A( 38) =    62.77 ; solar_flux_phcm2s1( 38) =    3.4
 lambda_EUV_A( 39) =    63.16 ; solar_flux_phcm2s1( 39) =    3.3
 lambda_EUV_A( 40) =    63.30 ; solar_flux_phcm2s1( 40) =    5.4
 lambda_EUV_A( 41) =    63.65 ; solar_flux_phcm2s1( 41) =    4.1
 lambda_EUV_A( 42) =    64.11 ; solar_flux_phcm2s1( 42) =    1.1
 lambda_EUV_A( 43) =    64.60 ; solar_flux_phcm2s1( 43) =    2.5
 lambda_EUV_A( 44) =    65.21 ; solar_flux_phcm2s1( 44) =    3.0
 lambda_EUV_A( 45) =    65.71 ; solar_flux_phcm2s1( 45) =    4.1
 lambda_EUV_A( 46) =    65.85 ; solar_flux_phcm2s1( 46) =    3.0
 lambda_EUV_A( 47) =    66.30 ; solar_flux_phcm2s1( 47) =    3.8
 lambda_EUV_A( 48) =    67.14 ; solar_flux_phcm2s1( 48) =    3.0
 lambda_EUV_A( 49) =    67.35 ; solar_flux_phcm2s1( 49) =    2.0
 lambda_EUV_A( 50) =    68.35 ; solar_flux_phcm2s1( 50) =    2.3
 lambda_EUV_A( 51) =    69.65 ; solar_flux_phcm2s1( 51) =   13.8
 lambda_EUV_A( 52) =    70.00 ; solar_flux_phcm2s1( 52) =     .1
 lambda_EUV_A( 53) =    70.54 ; solar_flux_phcm2s1( 53) =    3.3
 lambda_EUV_A( 54) =    70.75 ; solar_flux_phcm2s1( 54) =    3.0
 lambda_EUV_A( 55) =    71.00 ; solar_flux_phcm2s1( 55) =    4.4
 lambda_EUV_A( 56) =    71.94 ; solar_flux_phcm2s1( 56) =    1.2
 lambda_EUV_A( 57) =    72.31 ; solar_flux_phcm2s1( 57) =    6.4
 lambda_EUV_A( 58) =    72.63 ; solar_flux_phcm2s1( 58) =    1.6
 lambda_EUV_A( 59) =    72.80 ; solar_flux_phcm2s1( 59) =    2.1
 lambda_EUV_A( 60) =    72.95 ; solar_flux_phcm2s1( 60) =    3.4
 lambda_EUV_A( 61) =    73.55 ; solar_flux_phcm2s1( 61) =    1.9
 lambda_EUV_A( 62) =    74.21 ; solar_flux_phcm2s1( 62) =    2.5
 lambda_EUV_A( 63) =    74.44 ; solar_flux_phcm2s1( 63) =    1.3
 lambda_EUV_A( 64) =    74.83 ; solar_flux_phcm2s1( 64) =    4.0
 lambda_EUV_A( 65) =    75.03 ; solar_flux_phcm2s1( 65) =    4.6
 lambda_EUV_A( 66) =    75.29 ; solar_flux_phcm2s1( 66) =    2.5
 lambda_EUV_A( 67) =    75.46 ; solar_flux_phcm2s1( 67) =    3.8
 lambda_EUV_A( 68) =    75.73 ; solar_flux_phcm2s1( 68) =    2.5
 lambda_EUV_A( 69) =    76.01 ; solar_flux_phcm2s1( 69) =    2.7
 lambda_EUV_A( 70) =    76.48 ; solar_flux_phcm2s1( 70) =    1.0
 lambda_EUV_A( 71) =    76.83 ; solar_flux_phcm2s1( 71) =    4.1
 lambda_EUV_A( 72) =    76.94 ; solar_flux_phcm2s1( 72) =    3.3
 lambda_EUV_A( 73) =    77.30 ; solar_flux_phcm2s1( 73) =    2.9
 lambda_EUV_A( 74) =    77.74 ; solar_flux_phcm2s1( 74) =    3.9
 lambda_EUV_A( 75) =    78.56 ; solar_flux_phcm2s1( 75) =    3.0
 lambda_EUV_A( 76) =    78.70 ; solar_flux_phcm2s1( 76) =    2.8
 lambda_EUV_A( 77) =    79.08 ; solar_flux_phcm2s1( 77) =    1.9
 lambda_EUV_A( 78) =    79.48 ; solar_flux_phcm2s1( 78) =    1.8
 lambda_EUV_A( 79) =    79.76 ; solar_flux_phcm2s1( 79) =    2.3
 lambda_EUV_A( 80) =    80.00 ; solar_flux_phcm2s1( 80) =    1.4
 lambda_EUV_A( 81) =    80.21 ; solar_flux_phcm2s1( 81) =    2.9
 lambda_EUV_A( 82) =    80.55 ; solar_flux_phcm2s1( 82) =    2.3
 lambda_EUV_A( 83) =    80.94 ; solar_flux_phcm2s1( 83) =    2.0
 lambda_EUV_A( 84) =    81.16 ; solar_flux_phcm2s1( 84) =    2.6
 lambda_EUV_A( 85) =    81.58 ; solar_flux_phcm2s1( 85) =    3.0
 lambda_EUV_A( 86) =    81.94 ; solar_flux_phcm2s1( 86) =    4.2
 lambda_EUV_A( 87) =    82.43 ; solar_flux_phcm2s1( 87) =    4.9
 lambda_EUV_A( 88) =    82.67 ; solar_flux_phcm2s1( 88) =    8.4
 lambda_EUV_A( 89) =    83.25 ; solar_flux_phcm2s1( 89) =    4.4
 lambda_EUV_A( 90) =    83.42 ; solar_flux_phcm2s1( 90) =    4.5
 lambda_EUV_A( 91) =    83.67 ; solar_flux_phcm2s1( 91) =    3.8
 lambda_EUV_A( 92) =    84.00 ; solar_flux_phcm2s1( 92) =    5.1
 lambda_EUV_A( 93) =    84.26 ; solar_flux_phcm2s1( 93) =    3.8
 lambda_EUV_A( 94) =    84.50 ; solar_flux_phcm2s1( 94) =    5.9
 lambda_EUV_A( 95) =    84.72 ; solar_flux_phcm2s1( 95) =    4.1
 lambda_EUV_A( 96) =    84.86 ; solar_flux_phcm2s1( 96) =    3.9
 lambda_EUV_A( 97) =    85.16 ; solar_flux_phcm2s1( 97) =    3.0
 lambda_EUV_A( 98) =    85.50 ; solar_flux_phcm2s1( 98) =    5.9
 lambda_EUV_A( 99) =    85.69 ; solar_flux_phcm2s1( 99) =    2.7
 lambda_EUV_A(100) =    85.87 ; solar_flux_phcm2s1(100) =    2.9
 lambda_EUV_A(101) =    86.23 ; solar_flux_phcm2s1(101) =    2.1
 lambda_EUV_A(102) =    86.40 ; solar_flux_phcm2s1(102) =    1.7
 lambda_EUV_A(103) =    86.77 ; solar_flux_phcm2s1(103) =    5.0
 lambda_EUV_A(104) =    86.98 ; solar_flux_phcm2s1(104) =    3.1
 lambda_EUV_A(105) =    87.30 ; solar_flux_phcm2s1(105) =    2.4
 lambda_EUV_A(106) =    87.61 ; solar_flux_phcm2s1(106) =    2.0
 lambda_EUV_A(107) =    88.10 ; solar_flux_phcm2s1(107) =     .1
 lambda_EUV_A(108) =    88.11 ; solar_flux_phcm2s1(108) =    6.1
 lambda_EUV_A(109) =    88.14 ; solar_flux_phcm2s1(109) =     .1
 lambda_EUV_A(110) =    88.42 ; solar_flux_phcm2s1(110) =    1.9
 lambda_EUV_A(111) =    88.64 ; solar_flux_phcm2s1(111) =    2.4
 lambda_EUV_A(112) =    88.90 ; solar_flux_phcm2s1(112) =    3.8
 lambda_EUV_A(113) =    89.14 ; solar_flux_phcm2s1(113) =    2.7
 lambda_EUV_A(114) =    89.70 ; solar_flux_phcm2s1(114) =    3.0
 lambda_EUV_A(115) =    90.14 ; solar_flux_phcm2s1(115) =    3.1
 lambda_EUV_A(116) =    90.45 ; solar_flux_phcm2s1(116) =    2.0
 lambda_EUV_A(117) =    90.71 ; solar_flux_phcm2s1(117) =    3.0
 lambda_EUV_A(118) =    91.00 ; solar_flux_phcm2s1(118) =    3.9
 lambda_EUV_A(119) =    91.48 ; solar_flux_phcm2s1(119) =    1.9
 lambda_EUV_A(120) =    91.69 ; solar_flux_phcm2s1(120) =    5.4
 lambda_EUV_A(121) =    91.81 ; solar_flux_phcm2s1(121) =    4.8
 lambda_EUV_A(122) =    92.09 ; solar_flux_phcm2s1(122) =    3.8
 lambda_EUV_A(123) =    92.55 ; solar_flux_phcm2s1(123) =    3.0
 lambda_EUV_A(124) =    92.81 ; solar_flux_phcm2s1(124) =    3.8
 lambda_EUV_A(125) =    93.61 ; solar_flux_phcm2s1(125) =    5.6
 lambda_EUV_A(126) =    94.07 ; solar_flux_phcm2s1(126) =    7.1
 lambda_EUV_A(127) =    94.25 ; solar_flux_phcm2s1(127) =     .1
 lambda_EUV_A(128) =    94.39 ; solar_flux_phcm2s1(128) =    1.7
 lambda_EUV_A(129) =    94.81 ; solar_flux_phcm2s1(129) =    1.4
 lambda_EUV_A(130) =    94.90 ; solar_flux_phcm2s1(130) =     .1
 lambda_EUV_A(131) =    95.37 ; solar_flux_phcm2s1(131) =    5.2
 lambda_EUV_A(132) =    95.51 ; solar_flux_phcm2s1(132) =    2.9
 lambda_EUV_A(133) =    95.81 ; solar_flux_phcm2s1(133) =    2.9
 lambda_EUV_A(134) =    96.05 ; solar_flux_phcm2s1(134) =    9.7
 lambda_EUV_A(135) =    96.49 ; solar_flux_phcm2s1(135) =    1.9
 lambda_EUV_A(136) =    96.83 ; solar_flux_phcm2s1(136) =    2.7
 lambda_EUV_A(137) =    97.12 ; solar_flux_phcm2s1(137) =    5.2
 lambda_EUV_A(138) =    97.51 ; solar_flux_phcm2s1(138) =    3.0
 lambda_EUV_A(139) =    97.87 ; solar_flux_phcm2s1(139) =    2.3
 lambda_EUV_A(140) =    98.12 ; solar_flux_phcm2s1(140) =    5.6
 lambda_EUV_A(141) =    98.23 ; solar_flux_phcm2s1(141) =    6.5
 lambda_EUV_A(142) =    98.50 ; solar_flux_phcm2s1(142) =    2.6
 lambda_EUV_A(143) =    98.88 ; solar_flux_phcm2s1(143) =    1.5
 lambda_EUV_A(144) =    99.44 ; solar_flux_phcm2s1(144) =    1.6
 lambda_EUV_A(145) =    99.71 ; solar_flux_phcm2s1(145) =    2.0
 lambda_EUV_A(146) =    99.99 ; solar_flux_phcm2s1(146) =    2.7

! 100-150A, EUVAC bin 2
 lambda_EUV_A(147) =   100.54 ; solar_flux_phcm2s1(147) =    8.4
 lambda_EUV_A(148) =   100.96 ; solar_flux_phcm2s1(148) =    2.6
 lambda_EUV_A(149) =   101.57 ; solar_flux_phcm2s1(149) =    3.4
 lambda_EUV_A(150) =   102.15 ; solar_flux_phcm2s1(150) =    4.5
 lambda_EUV_A(151) =   103.01 ; solar_flux_phcm2s1(151) =    1.5
 lambda_EUV_A(152) =   103.15 ; solar_flux_phcm2s1(152) =     .1
 lambda_EUV_A(153) =   103.17 ; solar_flux_phcm2s1(153) =    4.2
 lambda_EUV_A(154) =   103.58 ; solar_flux_phcm2s1(154) =    6.1
 lambda_EUV_A(155) =   103.94 ; solar_flux_phcm2s1(155) =    6.1
 lambda_EUV_A(156) =   104.23 ; solar_flux_phcm2s1(156) =    1.6
 lambda_EUV_A(157) =   104.76 ; solar_flux_phcm2s1(157) =    2.0
 lambda_EUV_A(158) =   105.23 ; solar_flux_phcm2s1(158) =    5.1
 lambda_EUV_A(159) =   106.25 ; solar_flux_phcm2s1(159) =    2.0
 lambda_EUV_A(160) =   106.57 ; solar_flux_phcm2s1(160) =    1.1
 lambda_EUV_A(161) =   106.93 ; solar_flux_phcm2s1(161) =    1.9
 lambda_EUV_A(162) =   108.05 ; solar_flux_phcm2s1(162) =    1.5
 lambda_EUV_A(163) =   108.46 ; solar_flux_phcm2s1(163) =    2.0
 lambda_EUV_A(164) =   109.50 ; solar_flux_phcm2s1(164) =    2.0
 lambda_EUV_A(165) =   109.98 ; solar_flux_phcm2s1(165) =     .1
 lambda_EUV_A(166) =   110.56 ; solar_flux_phcm2s1(166) =     .1
 lambda_EUV_A(167) =   110.62 ; solar_flux_phcm2s1(167) =     .1
 lambda_EUV_A(168) =   110.76 ; solar_flux_phcm2s1(168) =    1.5
 lambda_EUV_A(169) =   111.16 ; solar_flux_phcm2s1(169) =     .1
 lambda_EUV_A(170) =   111.25 ; solar_flux_phcm2s1(170) =    4.8
 lambda_EUV_A(171) =   113.80 ; solar_flux_phcm2s1(171) =    3.0
 lambda_EUV_A(172) =   114.09 ; solar_flux_phcm2s1(172) =    2.6
 lambda_EUV_A(173) =   114.24 ; solar_flux_phcm2s1(173) =     .1
 lambda_EUV_A(174) =   115.39 ; solar_flux_phcm2s1(174) =     .1
 lambda_EUV_A(175) =   115.82 ; solar_flux_phcm2s1(175) =    2.4
 lambda_EUV_A(176) =   116.75 ; solar_flux_phcm2s1(176) =    3.9
 lambda_EUV_A(177) =   117.20 ; solar_flux_phcm2s1(177) =    2.6
 lambda_EUV_A(178) =   120.40 ; solar_flux_phcm2s1(178) =     .1
 lambda_EUV_A(179) =   121.15 ; solar_flux_phcm2s1(179) =     .1
 lambda_EUV_A(180) =   121.79 ; solar_flux_phcm2s1(180) =    1.1
 lambda_EUV_A(181) =   122.70 ; solar_flux_phcm2s1(181) =    4.1
 lambda_EUV_A(182) =   123.50 ; solar_flux_phcm2s1(182) =    2.6
 lambda_EUV_A(183) =   127.65 ; solar_flux_phcm2s1(183) =    4.5
 lambda_EUV_A(184) =   129.87 ; solar_flux_phcm2s1(184) =    3.4
 lambda_EUV_A(185) =   130.30 ; solar_flux_phcm2s1(185) =     .1
 lambda_EUV_A(186) =   131.02 ; solar_flux_phcm2s1(186) =    5.2
 lambda_EUV_A(187) =   131.21 ; solar_flux_phcm2s1(187) =    4.8
 lambda_EUV_A(188) =   136.21 ; solar_flux_phcm2s1(188) =     .1
 lambda_EUV_A(189) =   136.28 ; solar_flux_phcm2s1(189) =     .1
 lambda_EUV_A(190) =   136.34 ; solar_flux_phcm2s1(190) =     .1
 lambda_EUV_A(191) =   136.45 ; solar_flux_phcm2s1(191) =     .1
 lambda_EUV_A(192) =   136.48 ; solar_flux_phcm2s1(192) =     .1
 lambda_EUV_A(193) =   141.20 ; solar_flux_phcm2s1(193) =   11.3
 lambda_EUV_A(194) =   144.27 ; solar_flux_phcm2s1(194) =    1.1
 lambda_EUV_A(195) =   145.04 ; solar_flux_phcm2s1(195) =   14.1
 lambda_EUV_A(196) =   148.40 ; solar_flux_phcm2s1(196) =   38.2

! 150-200A, EUVAC bin 3
 lambda_EUV_A(197) =   150.10 ; solar_flux_phcm2s1(197) =   17.8
 lambda_EUV_A(198) =   152.15 ; solar_flux_phcm2s1(198) =   17.5
 lambda_EUV_A(199) =   154.18 ; solar_flux_phcm2s1(199) =   11.6
 lambda_EUV_A(200) =   157.73 ; solar_flux_phcm2s1(200) =    6.5
 lambda_EUV_A(201) =   158.37 ; solar_flux_phcm2s1(201) =   14.2
 lambda_EUV_A(202) =   159.98 ; solar_flux_phcm2s1(202) =   13.0
 lambda_EUV_A(203) =   160.37 ; solar_flux_phcm2s1(203) =    7.5
 lambda_EUV_A(204) =   164.15 ; solar_flux_phcm2s1(204) =    2.6
 lambda_EUV_A(205) =   167.50 ; solar_flux_phcm2s1(205) =   19.5
 lambda_EUV_A(206) =   168.17 ; solar_flux_phcm2s1(206) =   36.0
 lambda_EUV_A(207) =   168.55 ; solar_flux_phcm2s1(207) =   20.5
 lambda_EUV_A(208) =   168.92 ; solar_flux_phcm2s1(208) =   12.6
 lambda_EUV_A(209) =   169.70 ; solar_flux_phcm2s1(209) =   23.4
 lambda_EUV_A(210) =   171.08 ; solar_flux_phcm2s1(210) =  307.5
 lambda_EUV_A(211) =   172.17 ; solar_flux_phcm2s1(211) =   11.2
 lambda_EUV_A(212) =   173.08 ; solar_flux_phcm2s1(212) =   21.4
 lambda_EUV_A(213) =   174.58 ; solar_flux_phcm2s1(213) =  268.8
 lambda_EUV_A(214) =   175.26 ; solar_flux_phcm2s1(214) =   33.1
 lambda_EUV_A(215) =   177.24 ; solar_flux_phcm2s1(215) =  142.8
 lambda_EUV_A(216) =   178.05 ; solar_flux_phcm2s1(216) =   15.5
 lambda_EUV_A(217) =   179.27 ; solar_flux_phcm2s1(217) =     .3
 lambda_EUV_A(218) =   179.75 ; solar_flux_phcm2s1(218) =   15.6
 lambda_EUV_A(219) =   180.41 ; solar_flux_phcm2s1(219) =  138.0
 lambda_EUV_A(220) =   181.14 ; solar_flux_phcm2s1(220) =   21.3
 lambda_EUV_A(221) =   182.17 ; solar_flux_phcm2s1(221) =   39.2
 lambda_EUV_A(222) =   183.45 ; solar_flux_phcm2s1(222) =   10.0
 lambda_EUV_A(223) =   184.53 ; solar_flux_phcm2s1(223) =   78.2
 lambda_EUV_A(224) =   184.80 ; solar_flux_phcm2s1(224) =    4.0
 lambda_EUV_A(225) =   185.21 ; solar_flux_phcm2s1(225) =   32.5
 lambda_EUV_A(226) =   186.60 ; solar_flux_phcm2s1(226) =   16.5
 lambda_EUV_A(227) =   186.87 ; solar_flux_phcm2s1(227) =   30.0
 lambda_EUV_A(228) =   187.95 ; solar_flux_phcm2s1(228) =    1.5
 lambda_EUV_A(229) =   188.23 ; solar_flux_phcm2s1(229) =    9.7
 lambda_EUV_A(230) =   188.31 ; solar_flux_phcm2s1(230) =  120.0
 lambda_EUV_A(231) =   190.02 ; solar_flux_phcm2s1(231) =   50.4
 lambda_EUV_A(232) =   191.04 ; solar_flux_phcm2s1(232) =   10.5
 lambda_EUV_A(233) =   191.34 ; solar_flux_phcm2s1(233) =    9.0
 lambda_EUV_A(234) =   192.40 ; solar_flux_phcm2s1(234) =   32.0
 lambda_EUV_A(235) =   192.82 ; solar_flux_phcm2s1(235) =   49.0
 lambda_EUV_A(236) =   193.52 ; solar_flux_phcm2s1(236) =   57.0
 lambda_EUV_A(237) =   195.13 ; solar_flux_phcm2s1(237) =  100.0
 lambda_EUV_A(238) =   196.52 ; solar_flux_phcm2s1(238) =   17.0
 lambda_EUV_A(239) =   196.65 ; solar_flux_phcm2s1(239) =    5.0
 lambda_EUV_A(240) =   197.44 ; solar_flux_phcm2s1(240) =    7.0
 lambda_EUV_A(241) =   198.58 ; solar_flux_phcm2s1(241) =   12.5

! 200-250A, EUVAC bin 4
 lambda_EUV_A(242) =   200.02 ; solar_flux_phcm2s1(242) =   15.0
 lambda_EUV_A(243) =   201.13 ; solar_flux_phcm2s1(243) =   25.0
 lambda_EUV_A(244) =   202.05 ; solar_flux_phcm2s1(244) =   38.0
 lambda_EUV_A(245) =   202.64 ; solar_flux_phcm2s1(245) =   21.0
 lambda_EUV_A(246) =   203.81 ; solar_flux_phcm2s1(246) =   16.0
 lambda_EUV_A(247) =   204.25 ; solar_flux_phcm2s1(247) =    6.0
 lambda_EUV_A(248) =   204.94 ; solar_flux_phcm2s1(248) =    4.0
 lambda_EUV_A(249) =   206.26 ; solar_flux_phcm2s1(249) =    1.1
 lambda_EUV_A(250) =   206.38 ; solar_flux_phcm2s1(250) =    1.1
 lambda_EUV_A(251) =   207.46 ; solar_flux_phcm2s1(251) =    1.1
 lambda_EUV_A(252) =   208.33 ; solar_flux_phcm2s1(252) =    1.8
 lambda_EUV_A(253) =   209.63 ; solar_flux_phcm2s1(253) =     .9
 lambda_EUV_A(254) =   209.78 ; solar_flux_phcm2s1(254) =    1.1
 lambda_EUV_A(255) =   211.32 ; solar_flux_phcm2s1(255) =   21.6
 lambda_EUV_A(256) =   212.14 ; solar_flux_phcm2s1(256) =   10.5
 lambda_EUV_A(257) =   213.78 ; solar_flux_phcm2s1(257) =    4.9
 lambda_EUV_A(258) =   214.75 ; solar_flux_phcm2s1(258) =    8.6
 lambda_EUV_A(259) =   215.16 ; solar_flux_phcm2s1(259) =   19.3
 lambda_EUV_A(260) =   216.88 ; solar_flux_phcm2s1(260) =   36.9
 lambda_EUV_A(261) =   217.00 ; solar_flux_phcm2s1(261) =   63.8
 lambda_EUV_A(262) =   218.19 ; solar_flux_phcm2s1(262) =   57.7
 lambda_EUV_A(263) =   219.13 ; solar_flux_phcm2s1(263) =    4.8
 lambda_EUV_A(264) =   220.08 ; solar_flux_phcm2s1(264) =    6.8
 lambda_EUV_A(265) =   221.44 ; solar_flux_phcm2s1(265) =   17.8
 lambda_EUV_A(266) =   221.82 ; solar_flux_phcm2s1(266) =     .9
 lambda_EUV_A(267) =   224.74 ; solar_flux_phcm2s1(267) =   36.3
 lambda_EUV_A(268) =   225.12 ; solar_flux_phcm2s1(268) =   63.8
 lambda_EUV_A(269) =   227.01 ; solar_flux_phcm2s1(269) =   86.9
 lambda_EUV_A(270) =   227.19 ; solar_flux_phcm2s1(270) =     .6
 lambda_EUV_A(271) =   227.47 ; solar_flux_phcm2s1(271) =   26.4
 lambda_EUV_A(272) =   228.70 ; solar_flux_phcm2s1(272) =   17.7
 lambda_EUV_A(273) =   230.65 ; solar_flux_phcm2s1(273) =   11.8
 lambda_EUV_A(274) =   231.55 ; solar_flux_phcm2s1(274) =   14.4
 lambda_EUV_A(275) =   232.60 ; solar_flux_phcm2s1(275) =   21.0
 lambda_EUV_A(276) =   233.84 ; solar_flux_phcm2s1(276) =    3.3
 lambda_EUV_A(277) =   234.38 ; solar_flux_phcm2s1(277) =  103.3
 lambda_EUV_A(278) =   237.33 ; solar_flux_phcm2s1(278) =   61.9
 lambda_EUV_A(279) =   239.87 ; solar_flux_phcm2s1(279) =   16.1
 lambda_EUV_A(280) =   240.71 ; solar_flux_phcm2s1(280) =   26.5
 lambda_EUV_A(281) =   241.74 ; solar_flux_phcm2s1(281) =  141.4
 lambda_EUV_A(282) =   243.03 ; solar_flux_phcm2s1(282) =  110.3
 lambda_EUV_A(283) =   243.86 ; solar_flux_phcm2s1(283) =   62.5
 lambda_EUV_A(284) =   244.92 ; solar_flux_phcm2s1(284) =   69.2
 lambda_EUV_A(285) =   245.94 ; solar_flux_phcm2s1(285) =   18.6
 lambda_EUV_A(286) =   246.24 ; solar_flux_phcm2s1(286) =   62.5
 lambda_EUV_A(287) =   246.91 ; solar_flux_phcm2s1(287) =   23.4
 lambda_EUV_A(288) =   247.18 ; solar_flux_phcm2s1(288) =   24.9
 lambda_EUV_A(289) =   249.18 ; solar_flux_phcm2s1(289) =    6.6

! 250-300A, EUVAC bin 7
 lambda_EUV_A(290) =   251.10 ; solar_flux_phcm2s1(290) =    2.4
 lambda_EUV_A(291) =   251.95 ; solar_flux_phcm2s1(291) =   29.0
 lambda_EUV_A(292) =   252.19 ; solar_flux_phcm2s1(292) =   14.0
 lambda_EUV_A(293) =   253.78 ; solar_flux_phcm2s1(293) =   27.3
 lambda_EUV_A(294) =   256.32 ; solar_flux_phcm2s1(294) =  450.0   ! EUVAC bin/line 5
 lambda_EUV_A(295) =   256.38 ; solar_flux_phcm2s1(295) =   56.4
 lambda_EUV_A(296) =   256.64 ; solar_flux_phcm2s1(296) =   24.4
 lambda_EUV_A(297) =   257.16 ; solar_flux_phcm2s1(297) =  192.5
 lambda_EUV_A(298) =   257.39 ; solar_flux_phcm2s1(298) =   80.0
 lambda_EUV_A(299) =   258.36 ; solar_flux_phcm2s1(299) =  131.6
 lambda_EUV_A(300) =   259.52 ; solar_flux_phcm2s1(300) =   53.9
 lambda_EUV_A(301) =   261.05 ; solar_flux_phcm2s1(301) =   58.3
 lambda_EUV_A(302) =   262.99 ; solar_flux_phcm2s1(302) =    2.5
 lambda_EUV_A(303) =   264.24 ; solar_flux_phcm2s1(303) =   65.4
 lambda_EUV_A(304) =   264.80 ; solar_flux_phcm2s1(304) =   25.6
 lambda_EUV_A(305) =   270.51 ; solar_flux_phcm2s1(305) =   12.0
 lambda_EUV_A(306) =   271.99 ; solar_flux_phcm2s1(306) =   56.4
 lambda_EUV_A(307) =   272.64 ; solar_flux_phcm2s1(307) =   27.5
 lambda_EUV_A(308) =   274.19 ; solar_flux_phcm2s1(308) =   36.0
 lambda_EUV_A(309) =   275.35 ; solar_flux_phcm2s1(309) =   33.0
 lambda_EUV_A(310) =   275.67 ; solar_flux_phcm2s1(310) =   27.5
 lambda_EUV_A(311) =   276.15 ; solar_flux_phcm2s1(311) =    6.1
 lambda_EUV_A(312) =   276.84 ; solar_flux_phcm2s1(312) =   12.5
 lambda_EUV_A(313) =   277.00 ; solar_flux_phcm2s1(313) =   50.0
 lambda_EUV_A(314) =   277.27 ; solar_flux_phcm2s1(314) =   75.0
 lambda_EUV_A(315) =   278.40 ; solar_flux_phcm2s1(315) =   54.4
 lambda_EUV_A(316) =   281.41 ; solar_flux_phcm2s1(316) =   14.3
 lambda_EUV_A(317) =   284.15 ; solar_flux_phcm2s1(317) =   77.3    ! EUVAC bin/line 6
 lambda_EUV_A(318) =   285.70 ; solar_flux_phcm2s1(318) =   22.9
 lambda_EUV_A(319) =   289.17 ; solar_flux_phcm2s1(319) =    5.2
 lambda_EUV_A(320) =   290.69 ; solar_flux_phcm2s1(320) =   38.5
 lambda_EUV_A(321) =   291.70 ; solar_flux_phcm2s1(321) =   17.2
 lambda_EUV_A(322) =   292.78 ; solar_flux_phcm2s1(322) =   61.0
 lambda_EUV_A(323) =   296.19 ; solar_flux_phcm2s1(323) =   90.5
 lambda_EUV_A(324) =   299.50 ; solar_flux_phcm2s1(324) =    8.7

! 300-350A, EUVAC bin 10
 lambda_EUV_A(325) =   303.31 ; solar_flux_phcm2s1(325) =  600.0    ! EUVAC bin/line 8
 lambda_EUV_A(326) =   303.78 ; solar_flux_phcm2s1(326) = 7762.5    ! EUVAC bin/line 9
 lambda_EUV_A(327) =   315.02 ; solar_flux_phcm2s1(327) =  133.5
 lambda_EUV_A(328) =   316.20 ; solar_flux_phcm2s1(328) =  112.5
 lambda_EUV_A(329) =   319.01 ; solar_flux_phcm2s1(329) =   17.5
 lambda_EUV_A(330) =   319.83 ; solar_flux_phcm2s1(330) =  146.3
 lambda_EUV_A(331) =   320.56 ; solar_flux_phcm2s1(331) =    4.7
 lambda_EUV_A(332) =   335.41 ; solar_flux_phcm2s1(332) =   35.0
 lambda_EUV_A(333) =   345.13 ; solar_flux_phcm2s1(333) =   96.3
 lambda_EUV_A(334) =   345.74 ; solar_flux_phcm2s1(334) =   84.0
 lambda_EUV_A(335) =   347.39 ; solar_flux_phcm2s1(335) =  141.0
 lambda_EUV_A(336) =   349.85 ; solar_flux_phcm2s1(336) =   96.3

! 350-400A, EUVAC bin 12
 lambda_EUV_A(337) =   356.01 ; solar_flux_phcm2s1(337) =  103.4
 lambda_EUV_A(338) =   360.80 ; solar_flux_phcm2s1(338) =   17.5
 lambda_EUV_A(339) =   364.48 ; solar_flux_phcm2s1(339) =   75.6
 lambda_EUV_A(340) =   368.07 ; solar_flux_phcm2s1(340) =  739.4    ! EUVAC bin/line 11
 lambda_EUV_A(341) =   399.82 ; solar_flux_phcm2s1(341) =   15.6

! 400-450A, EUVAC bin 13
 lambda_EUV_A(342) =   401.14 ; solar_flux_phcm2s1(342) =   34.5
 lambda_EUV_A(343) =   401.94 ; solar_flux_phcm2s1(343) =   91.2
 lambda_EUV_A(344) =   403.26 ; solar_flux_phcm2s1(344) =   53.4
 lambda_EUV_A(345) =   417.24 ; solar_flux_phcm2s1(345) =    8.8
 lambda_EUV_A(346) =   423.00 ; solar_flux_phcm2s1(346) =     .4
 lambda_EUV_A(347) =   424.00 ; solar_flux_phcm2s1(347) =     .4
 lambda_EUV_A(348) =   425.00 ; solar_flux_phcm2s1(348) =     .4
 lambda_EUV_A(349) =   426.00 ; solar_flux_phcm2s1(349) =     .4
 lambda_EUV_A(350) =   427.00 ; solar_flux_phcm2s1(350) =     .4
 lambda_EUV_A(351) =   428.00 ; solar_flux_phcm2s1(351) =     .4
 lambda_EUV_A(352) =   429.00 ; solar_flux_phcm2s1(352) =     .4
 lambda_EUV_A(353) =   430.00 ; solar_flux_phcm2s1(353) =     .4
 lambda_EUV_A(354) =   430.47 ; solar_flux_phcm2s1(354) =   82.3
 lambda_EUV_A(355) =   431.00 ; solar_flux_phcm2s1(355) =     .4
 lambda_EUV_A(356) =   432.00 ; solar_flux_phcm2s1(356) =     .4
 lambda_EUV_A(357) =   433.00 ; solar_flux_phcm2s1(357) =     .4
 lambda_EUV_A(358) =   434.00 ; solar_flux_phcm2s1(358) =     .4
 lambda_EUV_A(359) =   435.00 ; solar_flux_phcm2s1(359) =     .4
 lambda_EUV_A(360) =   436.00 ; solar_flux_phcm2s1(360) =     .4
 lambda_EUV_A(361) =   436.70 ; solar_flux_phcm2s1(361) =  122.4
 lambda_EUV_A(362) =   437.00 ; solar_flux_phcm2s1(362) =     .4
 lambda_EUV_A(363) =   438.00 ; solar_flux_phcm2s1(363) =     .4
 lambda_EUV_A(364) =   439.00 ; solar_flux_phcm2s1(364) =     .4
 lambda_EUV_A(365) =   440.00 ; solar_flux_phcm2s1(365) =     .4
 lambda_EUV_A(366) =   441.00 ; solar_flux_phcm2s1(366) =     .4
 lambda_EUV_A(367) =   442.00 ; solar_flux_phcm2s1(367) =     .4
 lambda_EUV_A(368) =   443.00 ; solar_flux_phcm2s1(368) =     .4
 lambda_EUV_A(369) =   444.00 ; solar_flux_phcm2s1(369) =     .9
 lambda_EUV_A(370) =   445.00 ; solar_flux_phcm2s1(370) =     .9
 lambda_EUV_A(371) =   446.00 ; solar_flux_phcm2s1(371) =     .9
 lambda_EUV_A(372) =   447.00 ; solar_flux_phcm2s1(372) =     .9
 lambda_EUV_A(373) =   448.00 ; solar_flux_phcm2s1(373) =     .9
 lambda_EUV_A(374) =   449.00 ; solar_flux_phcm2s1(374) =     .9

! 450-500A, EUVAC bin 15
 lambda_EUV_A(375) =   450.00 ; solar_flux_phcm2s1(375) =     .9
 lambda_EUV_A(376) =   451.00 ; solar_flux_phcm2s1(376) =     .9
 lambda_EUV_A(377) =   452.00 ; solar_flux_phcm2s1(377) =     .9
 lambda_EUV_A(378) =   453.00 ; solar_flux_phcm2s1(378) =    1.3
 lambda_EUV_A(379) =   454.00 ; solar_flux_phcm2s1(379) =    1.3
 lambda_EUV_A(380) =   455.00 ; solar_flux_phcm2s1(380) =    1.3
 lambda_EUV_A(381) =   456.00 ; solar_flux_phcm2s1(381) =    1.3
 lambda_EUV_A(382) =   457.00 ; solar_flux_phcm2s1(382) =    1.3
 lambda_EUV_A(383) =   458.00 ; solar_flux_phcm2s1(383) =    1.3
 lambda_EUV_A(384) =   459.00 ; solar_flux_phcm2s1(384) =    1.3
 lambda_EUV_A(385) =   460.00 ; solar_flux_phcm2s1(385) =    1.6
 lambda_EUV_A(386) =   461.00 ; solar_flux_phcm2s1(386) =    1.6
 lambda_EUV_A(387) =   462.00 ; solar_flux_phcm2s1(387) =    1.6
 lambda_EUV_A(388) =   463.00 ; solar_flux_phcm2s1(388) =    1.6
 lambda_EUV_A(389) =   464.00 ; solar_flux_phcm2s1(389) =    1.6
 lambda_EUV_A(390) =   465.00 ; solar_flux_phcm2s1(390) =    2.0
 lambda_EUV_A(391) =   465.22 ; solar_flux_phcm2s1(391) =  329.9    ! EUVAC bin/line 14
 lambda_EUV_A(392) =   466.00 ; solar_flux_phcm2s1(392) =    2.0
 lambda_EUV_A(393) =   467.00 ; solar_flux_phcm2s1(393) =    2.0
 lambda_EUV_A(394) =   468.00 ; solar_flux_phcm2s1(394) =    2.5
 lambda_EUV_A(395) =   469.00 ; solar_flux_phcm2s1(395) =    2.5
 lambda_EUV_A(396) =   470.00 ; solar_flux_phcm2s1(396) =    2.5
 lambda_EUV_A(397) =   471.00 ; solar_flux_phcm2s1(397) =    2.5
 lambda_EUV_A(398) =   472.00 ; solar_flux_phcm2s1(398) =    2.9
 lambda_EUV_A(399) =   473.00 ; solar_flux_phcm2s1(399) =    2.9
 lambda_EUV_A(400) =   474.00 ; solar_flux_phcm2s1(400) =    3.3
 lambda_EUV_A(401) =   475.00 ; solar_flux_phcm2s1(401) =    3.3
 lambda_EUV_A(402) =   476.00 ; solar_flux_phcm2s1(402) =    3.3
 lambda_EUV_A(403) =   477.00 ; solar_flux_phcm2s1(403) =    3.8
 lambda_EUV_A(404) =   478.00 ; solar_flux_phcm2s1(404) =    3.8
 lambda_EUV_A(405) =   479.00 ; solar_flux_phcm2s1(405) =    4.1
 lambda_EUV_A(406) =   480.00 ; solar_flux_phcm2s1(406) =    4.1
 lambda_EUV_A(407) =   481.00 ; solar_flux_phcm2s1(407) =    4.5
 lambda_EUV_A(408) =   482.00 ; solar_flux_phcm2s1(408) =    4.5
 lambda_EUV_A(409) =   483.00 ; solar_flux_phcm2s1(409) =    5.0
 lambda_EUV_A(410) =   484.00 ; solar_flux_phcm2s1(410) =    5.4
 lambda_EUV_A(411) =   485.00 ; solar_flux_phcm2s1(411) =    5.4
 lambda_EUV_A(412) =   486.00 ; solar_flux_phcm2s1(412) =    5.8
 lambda_EUV_A(413) =   487.00 ; solar_flux_phcm2s1(413) =    6.1
 lambda_EUV_A(414) =   488.00 ; solar_flux_phcm2s1(414) =    6.6
 lambda_EUV_A(415) =   489.00 ; solar_flux_phcm2s1(415) =    6.6
 lambda_EUV_A(416) =   489.50 ; solar_flux_phcm2s1(416) =   11.0
 lambda_EUV_A(417) =   490.00 ; solar_flux_phcm2s1(417) =    7.0
 lambda_EUV_A(418) =   491.00 ; solar_flux_phcm2s1(418) =    7.4
 lambda_EUV_A(419) =   492.00 ; solar_flux_phcm2s1(419) =    7.9
 lambda_EUV_A(420) =   493.00 ; solar_flux_phcm2s1(420) =    8.3
 lambda_EUV_A(421) =   494.00 ; solar_flux_phcm2s1(421) =    8.6
 lambda_EUV_A(422) =   495.00 ; solar_flux_phcm2s1(422) =    9.5
 lambda_EUV_A(423) =   496.00 ; solar_flux_phcm2s1(423) =    9.9
 lambda_EUV_A(424) =   497.00 ; solar_flux_phcm2s1(424) =   10.4
 lambda_EUV_A(425) =   498.00 ; solar_flux_phcm2s1(425) =   10.8
 lambda_EUV_A(426) =   499.00 ; solar_flux_phcm2s1(426) =   11.5
 lambda_EUV_A(427) =   499.37 ; solar_flux_phcm2s1(427) =   77.5

! 500-550A, EUVAC bin 16
 lambda_EUV_A(428) =   500.00 ; solar_flux_phcm2s1(428) =   12.0
 lambda_EUV_A(429) =   501.00 ; solar_flux_phcm2s1(429) =   12.7
 lambda_EUV_A(430) =   502.00 ; solar_flux_phcm2s1(430) =   13.3
 lambda_EUV_A(431) =   503.00 ; solar_flux_phcm2s1(431) =   14.0
 lambda_EUV_A(432) =   504.00 ; solar_flux_phcm2s1(432) =   14.9
 lambda_EUV_A(433) =   507.93 ; solar_flux_phcm2s1(433) =  119.6
 lambda_EUV_A(434) =   515.60 ; solar_flux_phcm2s1(434) =   31.3
 lambda_EUV_A(435) =   520.66 ; solar_flux_phcm2s1(435) =   35.6
 lambda_EUV_A(436) =   525.80 ; solar_flux_phcm2s1(436) =   70.7
 lambda_EUV_A(437) =   537.02 ; solar_flux_phcm2s1(437) =  150.0
 lambda_EUV_A(438) =   542.80 ; solar_flux_phcm2s1(438) =   23.2

! 550-600A, EUVAC bin 19
 lambda_EUV_A(439) =   550.00 ; solar_flux_phcm2s1(439) =   23.2
 lambda_EUV_A(440) =   554.37 ; solar_flux_phcm2s1(440) =  799.2    ! EUVAC bin/line 17
 lambda_EUV_A(441) =   558.60 ; solar_flux_phcm2s1(441) =   61.5
 lambda_EUV_A(442) =   562.80 ; solar_flux_phcm2s1(442) =   81.9
 lambda_EUV_A(443) =   568.50 ; solar_flux_phcm2s1(443) =   51.7
 lambda_EUV_A(444) =   572.30 ; solar_flux_phcm2s1(444) =   68.9
 lambda_EUV_A(445) =   580.40 ; solar_flux_phcm2s1(445) =   11.1
 lambda_EUV_A(446) =   584.33 ; solar_flux_phcm2s1(446) = 1587.5    ! EUVAC bin/line 18
 lambda_EUV_A(447) =   592.40 ; solar_flux_phcm2s1(447) =   19.2
 lambda_EUV_A(448) =   599.60 ; solar_flux_phcm2s1(448) =  190.0

! 600-650A, EUVAC bin 22
 lambda_EUV_A(449) =   608.00 ; solar_flux_phcm2s1(449) =     .1
 lambda_EUV_A(450) =   609.00 ; solar_flux_phcm2s1(450) =     .1
 lambda_EUV_A(451) =   609.76 ; solar_flux_phcm2s1(451) =  633.3    ! EUVAC bin/line 20
 lambda_EUV_A(452) =   610.00 ; solar_flux_phcm2s1(452) =     .1
 lambda_EUV_A(453) =   611.00 ; solar_flux_phcm2s1(453) =     .1
 lambda_EUV_A(454) =   612.00 ; solar_flux_phcm2s1(454) =     .1
 lambda_EUV_A(455) =   613.00 ; solar_flux_phcm2s1(455) =     .1
 lambda_EUV_A(456) =   614.00 ; solar_flux_phcm2s1(456) =     .1
 lambda_EUV_A(457) =   615.00 ; solar_flux_phcm2s1(457) =     .1
 lambda_EUV_A(458) =   616.00 ; solar_flux_phcm2s1(458) =     .1
 lambda_EUV_A(459) =   616.60 ; solar_flux_phcm2s1(459) =   16.7
 lambda_EUV_A(460) =   617.00 ; solar_flux_phcm2s1(460) =     .1
 lambda_EUV_A(461) =   618.00 ; solar_flux_phcm2s1(461) =     .1
 lambda_EUV_A(462) =   619.00 ; solar_flux_phcm2s1(462) =     .1
 lambda_EUV_A(463) =   620.00 ; solar_flux_phcm2s1(463) =     .1
 lambda_EUV_A(464) =   621.00 ; solar_flux_phcm2s1(464) =     .1
 lambda_EUV_A(465) =   622.00 ; solar_flux_phcm2s1(465) =     .1
 lambda_EUV_A(466) =   623.00 ; solar_flux_phcm2s1(466) =     .1
 lambda_EUV_A(467) =   624.00 ; solar_flux_phcm2s1(467) =     .3
 lambda_EUV_A(468) =   624.93 ; solar_flux_phcm2s1(468) =  286.8
 lambda_EUV_A(469) =   625.00 ; solar_flux_phcm2s1(469) =     .3
 lambda_EUV_A(470) =   626.00 ; solar_flux_phcm2s1(470) =     .3
 lambda_EUV_A(471) =   627.00 ; solar_flux_phcm2s1(471) =     .3
 lambda_EUV_A(472) =   628.00 ; solar_flux_phcm2s1(472) =     .3
 lambda_EUV_A(473) =   629.00 ; solar_flux_phcm2s1(473) =     .3
 lambda_EUV_A(474) =   629.73 ; solar_flux_phcm2s1(474) = 1848.4    ! EUVAC bin/line 21
 lambda_EUV_A(475) =   630.00 ; solar_flux_phcm2s1(475) =     .3
 lambda_EUV_A(476) =   631.00 ; solar_flux_phcm2s1(476) =     .1
 lambda_EUV_A(477) =   632.00 ; solar_flux_phcm2s1(477) =     .1
 lambda_EUV_A(478) =   633.00 ; solar_flux_phcm2s1(478) =     .1
 lambda_EUV_A(479) =   634.00 ; solar_flux_phcm2s1(479) =     .1
 lambda_EUV_A(480) =   635.00 ; solar_flux_phcm2s1(480) =     .1
 lambda_EUV_A(481) =   636.00 ; solar_flux_phcm2s1(481) =     .1
 lambda_EUV_A(482) =   637.00 ; solar_flux_phcm2s1(482) =     .1
 lambda_EUV_A(483) =   638.00 ; solar_flux_phcm2s1(483) =     .1
 lambda_EUV_A(484) =   638.50 ; solar_flux_phcm2s1(484) =   24.9
 lambda_EUV_A(485) =   639.00 ; solar_flux_phcm2s1(485) =     .1
 lambda_EUV_A(486) =   640.00 ; solar_flux_phcm2s1(486) =     .1
 lambda_EUV_A(487) =   640.41 ; solar_flux_phcm2s1(487) =   11.6
 lambda_EUV_A(488) =   640.93 ; solar_flux_phcm2s1(488) =   13.7
 lambda_EUV_A(489) =   641.00 ; solar_flux_phcm2s1(489) =     .1
 lambda_EUV_A(490) =   641.81 ; solar_flux_phcm2s1(490) =   18.0
 lambda_EUV_A(491) =   642.00 ; solar_flux_phcm2s1(491) =     .1
 lambda_EUV_A(492) =   643.00 ; solar_flux_phcm2s1(492) =     .1
 lambda_EUV_A(493) =   644.00 ; solar_flux_phcm2s1(493) =     .3
 lambda_EUV_A(494) =   644.10 ; solar_flux_phcm2s1(494) =   21.4
 lambda_EUV_A(495) =   645.00 ; solar_flux_phcm2s1(495) =     .3
 lambda_EUV_A(496) =   646.00 ; solar_flux_phcm2s1(496) =     .3
 lambda_EUV_A(497) =   647.00 ; solar_flux_phcm2s1(497) =     .3
 lambda_EUV_A(498) =   648.00 ; solar_flux_phcm2s1(498) =     .3
 lambda_EUV_A(499) =   649.00 ; solar_flux_phcm2s1(499) =     .3

! 650-700A, EUVAC bin 23
 lambda_EUV_A(500) =   650.00 ; solar_flux_phcm2s1(500) =     .3
 lambda_EUV_A(501) =   650.30 ; solar_flux_phcm2s1(501) =   15.6
 lambda_EUV_A(502) =   651.00 ; solar_flux_phcm2s1(502) =     .3
 lambda_EUV_A(503) =   652.00 ; solar_flux_phcm2s1(503) =     .3
 lambda_EUV_A(504) =   653.00 ; solar_flux_phcm2s1(504) =     .3
 lambda_EUV_A(505) =   654.00 ; solar_flux_phcm2s1(505) =     .3
 lambda_EUV_A(506) =   655.00 ; solar_flux_phcm2s1(506) =     .3
 lambda_EUV_A(507) =   656.00 ; solar_flux_phcm2s1(507) =     .3
 lambda_EUV_A(508) =   657.00 ; solar_flux_phcm2s1(508) =     .4
 lambda_EUV_A(509) =   657.30 ; solar_flux_phcm2s1(509) =   12.1
 lambda_EUV_A(510) =   658.00 ; solar_flux_phcm2s1(510) =     .4
 lambda_EUV_A(511) =   659.00 ; solar_flux_phcm2s1(511) =     .4
 lambda_EUV_A(512) =   660.00 ; solar_flux_phcm2s1(512) =     .4
 lambda_EUV_A(513) =   661.00 ; solar_flux_phcm2s1(513) =     .4
 lambda_EUV_A(514) =   661.40 ; solar_flux_phcm2s1(514) =   12.1
 lambda_EUV_A(515) =   662.00 ; solar_flux_phcm2s1(515) =     .4
 lambda_EUV_A(516) =   663.00 ; solar_flux_phcm2s1(516) =     .4
 lambda_EUV_A(517) =   664.00 ; solar_flux_phcm2s1(517) =     .4
 lambda_EUV_A(518) =   665.00 ; solar_flux_phcm2s1(518) =     .4
 lambda_EUV_A(519) =   666.00 ; solar_flux_phcm2s1(519) =     .4
 lambda_EUV_A(520) =   667.00 ; solar_flux_phcm2s1(520) =     .5
 lambda_EUV_A(521) =   668.00 ; solar_flux_phcm2s1(521) =     .5
 lambda_EUV_A(522) =   669.00 ; solar_flux_phcm2s1(522) =     .5
 lambda_EUV_A(523) =   670.00 ; solar_flux_phcm2s1(523) =     .5
 lambda_EUV_A(524) =   671.00 ; solar_flux_phcm2s1(524) =     .5
 lambda_EUV_A(525) =   671.50 ; solar_flux_phcm2s1(525) =   10.6
 lambda_EUV_A(526) =   672.00 ; solar_flux_phcm2s1(526) =     .5
 lambda_EUV_A(527) =   673.00 ; solar_flux_phcm2s1(527) =     .5
 lambda_EUV_A(528) =   674.00 ; solar_flux_phcm2s1(528) =     .5
 lambda_EUV_A(529) =   675.00 ; solar_flux_phcm2s1(529) =     .6
 lambda_EUV_A(530) =   676.00 ; solar_flux_phcm2s1(530) =     .6
 lambda_EUV_A(531) =   677.00 ; solar_flux_phcm2s1(531) =     .6
 lambda_EUV_A(532) =   678.00 ; solar_flux_phcm2s1(532) =     .6
 lambda_EUV_A(533) =   679.00 ; solar_flux_phcm2s1(533) =     .6
 lambda_EUV_A(534) =   680.00 ; solar_flux_phcm2s1(534) =     .6
 lambda_EUV_A(535) =   681.00 ; solar_flux_phcm2s1(535) =     .6
 lambda_EUV_A(536) =   681.70 ; solar_flux_phcm2s1(536) =   37.9
 lambda_EUV_A(537) =   682.00 ; solar_flux_phcm2s1(537) =     .6
 lambda_EUV_A(538) =   683.00 ; solar_flux_phcm2s1(538) =     .6
 lambda_EUV_A(539) =   684.00 ; solar_flux_phcm2s1(539) =     .6
 lambda_EUV_A(540) =   685.00 ; solar_flux_phcm2s1(540) =     .6
 lambda_EUV_A(541) =   685.71 ; solar_flux_phcm2s1(541) =  101.3
 lambda_EUV_A(542) =   686.00 ; solar_flux_phcm2s1(542) =     .6
 lambda_EUV_A(543) =   687.00 ; solar_flux_phcm2s1(543) =     .8
 lambda_EUV_A(544) =   688.00 ; solar_flux_phcm2s1(544) =     .8
 lambda_EUV_A(545) =   689.00 ; solar_flux_phcm2s1(545) =     .8
 lambda_EUV_A(546) =   690.00 ; solar_flux_phcm2s1(546) =     .8
 lambda_EUV_A(547) =   690.80 ; solar_flux_phcm2s1(547) =   20.9
 lambda_EUV_A(548) =   691.00 ; solar_flux_phcm2s1(548) =     .8
 lambda_EUV_A(549) =   692.00 ; solar_flux_phcm2s1(549) =     .9
 lambda_EUV_A(550) =   693.00 ; solar_flux_phcm2s1(550) =     .9
 lambda_EUV_A(551) =   694.00 ; solar_flux_phcm2s1(551) =     .9
 lambda_EUV_A(552) =   694.30 ; solar_flux_phcm2s1(552) =   22.6
 lambda_EUV_A(553) =   695.00 ; solar_flux_phcm2s1(553) =     .9
 lambda_EUV_A(554) =   696.00 ; solar_flux_phcm2s1(554) =     .9
 lambda_EUV_A(555) =   697.00 ; solar_flux_phcm2s1(555) =     .9
 lambda_EUV_A(556) =   698.00 ; solar_flux_phcm2s1(556) =     .9
 lambda_EUV_A(557) =   699.00 ; solar_flux_phcm2s1(557) =     .9

! 700-750A, EUVAC bin 25
 lambda_EUV_A(558) =   700.00 ; solar_flux_phcm2s1(558) =    1.0
 lambda_EUV_A(559) =   701.00 ; solar_flux_phcm2s1(559) =    1.0
 lambda_EUV_A(560) =   702.00 ; solar_flux_phcm2s1(560) =    1.0
 lambda_EUV_A(561) =   703.00 ; solar_flux_phcm2s1(561) =    1.1
 lambda_EUV_A(562) =   703.36 ; solar_flux_phcm2s1(562) =  391.5    ! EUVAC bin/line 24
 lambda_EUV_A(563) =   704.00 ; solar_flux_phcm2s1(563) =    1.1
 lambda_EUV_A(564) =   705.00 ; solar_flux_phcm2s1(564) =    1.1
 lambda_EUV_A(565) =   706.00 ; solar_flux_phcm2s1(565) =    1.1
 lambda_EUV_A(566) =   707.00 ; solar_flux_phcm2s1(566) =    1.1
 lambda_EUV_A(567) =   708.00 ; solar_flux_phcm2s1(567) =    1.1
 lambda_EUV_A(568) =   709.00 ; solar_flux_phcm2s1(568) =    1.1
 lambda_EUV_A(569) =   710.00 ; solar_flux_phcm2s1(569) =    1.3
 lambda_EUV_A(570) =   711.00 ; solar_flux_phcm2s1(570) =    1.3
 lambda_EUV_A(571) =   712.00 ; solar_flux_phcm2s1(571) =    1.3
 lambda_EUV_A(572) =   712.70 ; solar_flux_phcm2s1(572) =   12.4
 lambda_EUV_A(573) =   713.00 ; solar_flux_phcm2s1(573) =    1.4
 lambda_EUV_A(574) =   714.00 ; solar_flux_phcm2s1(574) =    1.4
 lambda_EUV_A(575) =   715.00 ; solar_flux_phcm2s1(575) =    1.5
 lambda_EUV_A(576) =   716.00 ; solar_flux_phcm2s1(576) =    1.5
 lambda_EUV_A(577) =   717.00 ; solar_flux_phcm2s1(577) =    1.5
 lambda_EUV_A(578) =   718.00 ; solar_flux_phcm2s1(578) =    1.7
 lambda_EUV_A(579) =   718.50 ; solar_flux_phcm2s1(579) =   53.1
 lambda_EUV_A(580) =   719.00 ; solar_flux_phcm2s1(580) =    1.7
 lambda_EUV_A(581) =   720.00 ; solar_flux_phcm2s1(581) =    1.7
 lambda_EUV_A(582) =   721.00 ; solar_flux_phcm2s1(582) =    1.7
 lambda_EUV_A(583) =   722.00 ; solar_flux_phcm2s1(583) =    1.8
 lambda_EUV_A(584) =   723.00 ; solar_flux_phcm2s1(584) =    1.8
 lambda_EUV_A(585) =   724.00 ; solar_flux_phcm2s1(585) =    1.9
 lambda_EUV_A(586) =   725.00 ; solar_flux_phcm2s1(586) =    1.9
 lambda_EUV_A(587) =   726.00 ; solar_flux_phcm2s1(587) =    1.9
 lambda_EUV_A(588) =   727.00 ; solar_flux_phcm2s1(588) =    1.9
 lambda_EUV_A(589) =   728.00 ; solar_flux_phcm2s1(589) =    2.0
 lambda_EUV_A(590) =   729.00 ; solar_flux_phcm2s1(590) =    2.0
 lambda_EUV_A(591) =   730.00 ; solar_flux_phcm2s1(591) =    2.2
 lambda_EUV_A(592) =   731.00 ; solar_flux_phcm2s1(592) =    2.2
 lambda_EUV_A(593) =   732.00 ; solar_flux_phcm2s1(593) =    2.2
 lambda_EUV_A(594) =   733.00 ; solar_flux_phcm2s1(594) =    2.4
 lambda_EUV_A(595) =   734.00 ; solar_flux_phcm2s1(595) =    2.4
 lambda_EUV_A(596) =   735.00 ; solar_flux_phcm2s1(596) =    2.4
 lambda_EUV_A(597) =   736.00 ; solar_flux_phcm2s1(597) =    2.5
 lambda_EUV_A(598) =   737.00 ; solar_flux_phcm2s1(598) =    2.5
 lambda_EUV_A(599) =   738.00 ; solar_flux_phcm2s1(599) =    2.7
 lambda_EUV_A(600) =   739.00 ; solar_flux_phcm2s1(600) =    2.7
 lambda_EUV_A(601) =   740.00 ; solar_flux_phcm2s1(601) =    2.7
 lambda_EUV_A(602) =   741.00 ; solar_flux_phcm2s1(602) =    2.8
 lambda_EUV_A(603) =   742.00 ; solar_flux_phcm2s1(603) =    2.9
 lambda_EUV_A(604) =   743.00 ; solar_flux_phcm2s1(604) =    2.9
 lambda_EUV_A(605) =   744.00 ; solar_flux_phcm2s1(605) =    3.0
 lambda_EUV_A(606) =   745.00 ; solar_flux_phcm2s1(606) =    3.0
 lambda_EUV_A(607) =   746.00 ; solar_flux_phcm2s1(607) =    3.2
 lambda_EUV_A(608) =   747.00 ; solar_flux_phcm2s1(608) =    3.3
 lambda_EUV_A(609) =   748.00 ; solar_flux_phcm2s1(609) =    3.4
 lambda_EUV_A(610) =   749.00 ; solar_flux_phcm2s1(610) =    3.4

! 750-800A, EUVAC bin 29
 lambda_EUV_A(611) =   750.00 ; solar_flux_phcm2s1(611) =    3.5
 lambda_EUV_A(612) =   750.01 ; solar_flux_phcm2s1(612) =   40.5
 lambda_EUV_A(613) =   751.00 ; solar_flux_phcm2s1(613) =    3.7
 lambda_EUV_A(614) =   752.00 ; solar_flux_phcm2s1(614) =    3.7
 lambda_EUV_A(615) =   753.00 ; solar_flux_phcm2s1(615) =    3.8
 lambda_EUV_A(616) =   754.00 ; solar_flux_phcm2s1(616) =    3.9
 lambda_EUV_A(617) =   755.00 ; solar_flux_phcm2s1(617) =    3.9
 lambda_EUV_A(618) =   756.00 ; solar_flux_phcm2s1(618) =    4.1
 lambda_EUV_A(619) =   757.00 ; solar_flux_phcm2s1(619) =    4.2
 lambda_EUV_A(620) =   758.00 ; solar_flux_phcm2s1(620) =    4.3
 lambda_EUV_A(621) =   758.68 ; solar_flux_phcm2s1(621) =   34.9
 lambda_EUV_A(622) =   759.00 ; solar_flux_phcm2s1(622) =    4.6
 lambda_EUV_A(623) =   759.44 ; solar_flux_phcm2s1(623) =   26.7
 lambda_EUV_A(624) =   760.00 ; solar_flux_phcm2s1(624) =    4.6
 lambda_EUV_A(625) =   760.30 ; solar_flux_phcm2s1(625) =   93.0
 lambda_EUV_A(626) =   761.00 ; solar_flux_phcm2s1(626) =    4.7
 lambda_EUV_A(627) =   761.13 ; solar_flux_phcm2s1(627) =   23.2
 lambda_EUV_A(628) =   762.00 ; solar_flux_phcm2s1(628) =    4.8
 lambda_EUV_A(629) =   762.00 ; solar_flux_phcm2s1(629) =   34.9
 lambda_EUV_A(630) =   763.00 ; solar_flux_phcm2s1(630) =    4.9
 lambda_EUV_A(631) =   764.00 ; solar_flux_phcm2s1(631) =    5.1
 lambda_EUV_A(632) =   765.00 ; solar_flux_phcm2s1(632) =    5.2
 lambda_EUV_A(633) =   765.15 ; solar_flux_phcm2s1(633) =  199.7    ! EUVAC bin/line 26
 lambda_EUV_A(634) =   766.00 ; solar_flux_phcm2s1(634) =    5.3
 lambda_EUV_A(635) =   767.00 ; solar_flux_phcm2s1(635) =    5.4
 lambda_EUV_A(636) =   768.00 ; solar_flux_phcm2s1(636) =    5.6
 lambda_EUV_A(637) =   769.00 ; solar_flux_phcm2s1(637) =    5.7
 lambda_EUV_A(638) =   770.00 ; solar_flux_phcm2s1(638) =    5.8
 lambda_EUV_A(639) =   770.41 ; solar_flux_phcm2s1(639) =  242.5    ! EUVAC bin/line 27(1)
 lambda_EUV_A(640) =   771.00 ; solar_flux_phcm2s1(640) =    6.1
 lambda_EUV_A(641) =   772.00 ; solar_flux_phcm2s1(641) =    6.2
 lambda_EUV_A(642) =   773.00 ; solar_flux_phcm2s1(642) =    6.3
 lambda_EUV_A(643) =   774.00 ; solar_flux_phcm2s1(643) =    6.6
 lambda_EUV_A(644) =   775.00 ; solar_flux_phcm2s1(644) =    6.6
 lambda_EUV_A(645) =   776.00 ; solar_flux_phcm2s1(645) =   11.8
 lambda_EUV_A(646) =   776.00 ; solar_flux_phcm2s1(646) =    6.8
 lambda_EUV_A(647) =   777.00 ; solar_flux_phcm2s1(647) =    7.1
 lambda_EUV_A(648) =   778.00 ; solar_flux_phcm2s1(648) =    7.2
 lambda_EUV_A(649) =   779.00 ; solar_flux_phcm2s1(649) =    7.3
 lambda_EUV_A(650) =   780.00 ; solar_flux_phcm2s1(650) =    7.6
 lambda_EUV_A(651) =   780.32 ; solar_flux_phcm2s1(651) =  130.6    ! EUVAC bin/line 27(2) ???
 lambda_EUV_A(652) =   781.00 ; solar_flux_phcm2s1(652) =    7.7
 lambda_EUV_A(653) =   782.00 ; solar_flux_phcm2s1(653) =    8.0
 lambda_EUV_A(654) =   783.00 ; solar_flux_phcm2s1(654) =    8.1
 lambda_EUV_A(655) =   784.00 ; solar_flux_phcm2s1(655) =    8.3
 lambda_EUV_A(656) =   785.00 ; solar_flux_phcm2s1(656) =    8.6
 lambda_EUV_A(657) =   786.00 ; solar_flux_phcm2s1(657) =    8.7
 lambda_EUV_A(658) =   786.47 ; solar_flux_phcm2s1(658) =  146.3
 lambda_EUV_A(659) =   787.00 ; solar_flux_phcm2s1(659) =    9.1
 lambda_EUV_A(660) =   787.71 ; solar_flux_phcm2s1(660) =  277.5    ! EUVAC bin/line 28(1) ???
 lambda_EUV_A(661) =   788.00 ; solar_flux_phcm2s1(661) =    9.2
 lambda_EUV_A(662) =   789.00 ; solar_flux_phcm2s1(662) =    9.5
 lambda_EUV_A(663) =   790.00 ; solar_flux_phcm2s1(663) =    9.6
 lambda_EUV_A(664) =   790.15 ; solar_flux_phcm2s1(664) =  477.3    ! EUVAC bin/line 28(2) ???
 lambda_EUV_A(665) =   791.00 ; solar_flux_phcm2s1(665) =   10.0
 lambda_EUV_A(666) =   792.00 ; solar_flux_phcm2s1(666) =   10.2
 lambda_EUV_A(667) =   793.00 ; solar_flux_phcm2s1(667) =   10.3
 lambda_EUV_A(668) =   794.00 ; solar_flux_phcm2s1(668) =   10.7
 lambda_EUV_A(669) =   795.00 ; solar_flux_phcm2s1(669) =   11.0
 lambda_EUV_A(670) =   796.00 ; solar_flux_phcm2s1(670) =   11.2
 lambda_EUV_A(671) =   797.00 ; solar_flux_phcm2s1(671) =   11.5
 lambda_EUV_A(672) =   798.00 ; solar_flux_phcm2s1(672) =   11.9
 lambda_EUV_A(673) =   799.00 ; solar_flux_phcm2s1(673) =   12.1

! 800-850A, EUVAC bin 30
 lambda_EUV_A(674) =   800.00 ; solar_flux_phcm2s1(674) =   12.4
 lambda_EUV_A(675) =   801.00 ; solar_flux_phcm2s1(675) =   12.9
 lambda_EUV_A(676) =   802.00 ; solar_flux_phcm2s1(676) =   13.1
 lambda_EUV_A(677) =   803.00 ; solar_flux_phcm2s1(677) =   13.5
 lambda_EUV_A(678) =   804.00 ; solar_flux_phcm2s1(678) =   13.9
 lambda_EUV_A(679) =   805.00 ; solar_flux_phcm2s1(679) =   14.1
 lambda_EUV_A(680) =   806.00 ; solar_flux_phcm2s1(680) =   14.5
 lambda_EUV_A(681) =   807.00 ; solar_flux_phcm2s1(681) =   14.9
 lambda_EUV_A(682) =   808.00 ; solar_flux_phcm2s1(682) =   15.3
 lambda_EUV_A(683) =   809.00 ; solar_flux_phcm2s1(683) =   15.6
 lambda_EUV_A(684) =   810.00 ; solar_flux_phcm2s1(684) =   16.0
 lambda_EUV_A(685) =   811.00 ; solar_flux_phcm2s1(685) =   16.5
 lambda_EUV_A(686) =   812.00 ; solar_flux_phcm2s1(686) =   16.9
 lambda_EUV_A(687) =   813.00 ; solar_flux_phcm2s1(687) =   17.3
 lambda_EUV_A(688) =   814.00 ; solar_flux_phcm2s1(688) =   17.8
 lambda_EUV_A(689) =   815.00 ; solar_flux_phcm2s1(689) =   18.3
 lambda_EUV_A(690) =   816.00 ; solar_flux_phcm2s1(690) =   18.6
 lambda_EUV_A(691) =   817.00 ; solar_flux_phcm2s1(691) =   19.1
 lambda_EUV_A(692) =   818.00 ; solar_flux_phcm2s1(692) =   19.5
 lambda_EUV_A(693) =   819.00 ; solar_flux_phcm2s1(693) =   20.2
 lambda_EUV_A(694) =   820.00 ; solar_flux_phcm2s1(694) =   20.5
 lambda_EUV_A(695) =   821.00 ; solar_flux_phcm2s1(695) =   21.2
 lambda_EUV_A(696) =   822.00 ; solar_flux_phcm2s1(696) =   21.7
 lambda_EUV_A(697) =   823.00 ; solar_flux_phcm2s1(697) =   22.3
 lambda_EUV_A(698) =   824.00 ; solar_flux_phcm2s1(698) =   22.9
 lambda_EUV_A(699) =   825.00 ; solar_flux_phcm2s1(699) =   23.4
 lambda_EUV_A(700) =   826.00 ; solar_flux_phcm2s1(700) =   24.0
 lambda_EUV_A(701) =   827.00 ; solar_flux_phcm2s1(701) =   24.7
 lambda_EUV_A(702) =   828.00 ; solar_flux_phcm2s1(702) =   25.2
 lambda_EUV_A(703) =   829.00 ; solar_flux_phcm2s1(703) =   25.9
 lambda_EUV_A(704) =   830.00 ; solar_flux_phcm2s1(704) =   26.6
 lambda_EUV_A(705) =   831.00 ; solar_flux_phcm2s1(705) =   27.2
 lambda_EUV_A(706) =   832.00 ; solar_flux_phcm2s1(706) =   27.9
 lambda_EUV_A(707) =   833.00 ; solar_flux_phcm2s1(707) =   28.7
 lambda_EUV_A(708) =   834.00 ; solar_flux_phcm2s1(708) =   29.3
 lambda_EUV_A(709) =   834.20 ; solar_flux_phcm2s1(709) =  666.5
 lambda_EUV_A(710) =   835.00 ; solar_flux_phcm2s1(710) =   30.2
 lambda_EUV_A(711) =   836.00 ; solar_flux_phcm2s1(711) =   30.8
 lambda_EUV_A(712) =   837.00 ; solar_flux_phcm2s1(712) =   31.6
 lambda_EUV_A(713) =   838.00 ; solar_flux_phcm2s1(713) =   32.4
 lambda_EUV_A(714) =   839.00 ; solar_flux_phcm2s1(714) =   33.3
 lambda_EUV_A(715) =   840.00 ; solar_flux_phcm2s1(715) =   34.1
 lambda_EUV_A(716) =   841.00 ; solar_flux_phcm2s1(716) =   35.1
 lambda_EUV_A(717) =   842.00 ; solar_flux_phcm2s1(717) =   35.8
 lambda_EUV_A(718) =   843.00 ; solar_flux_phcm2s1(718) =   36.8
 lambda_EUV_A(719) =   844.00 ; solar_flux_phcm2s1(719) =   37.7
 lambda_EUV_A(720) =   845.00 ; solar_flux_phcm2s1(720) =   38.7
 lambda_EUV_A(721) =   846.00 ; solar_flux_phcm2s1(721) =   39.7
 lambda_EUV_A(722) =   847.00 ; solar_flux_phcm2s1(722) =   40.7
 lambda_EUV_A(723) =   848.00 ; solar_flux_phcm2s1(723) =   41.7
 lambda_EUV_A(724) =   849.00 ; solar_flux_phcm2s1(724) =   42.7

! 850-900A, EUVAC bin 31
 lambda_EUV_A(725) =   850.00 ; solar_flux_phcm2s1(725) =   43.8
 lambda_EUV_A(726) =   851.00 ; solar_flux_phcm2s1(726) =   45.0
 lambda_EUV_A(727) =   852.00 ; solar_flux_phcm2s1(727) =   46.1
 lambda_EUV_A(728) =   853.00 ; solar_flux_phcm2s1(728) =   47.4
 lambda_EUV_A(729) =   854.00 ; solar_flux_phcm2s1(729) =   48.5
 lambda_EUV_A(730) =   855.00 ; solar_flux_phcm2s1(730) =   49.7
 lambda_EUV_A(731) =   856.00 ; solar_flux_phcm2s1(731) =   51.0
 lambda_EUV_A(732) =   857.00 ; solar_flux_phcm2s1(732) =   52.4
 lambda_EUV_A(733) =   858.00 ; solar_flux_phcm2s1(733) =   53.6
 lambda_EUV_A(734) =   859.00 ; solar_flux_phcm2s1(734) =   55.0
 lambda_EUV_A(735) =   860.00 ; solar_flux_phcm2s1(735) =   56.4
 lambda_EUV_A(736) =   861.00 ; solar_flux_phcm2s1(736) =   57.9
 lambda_EUV_A(737) =   862.00 ; solar_flux_phcm2s1(737) =   59.2
 lambda_EUV_A(738) =   863.00 ; solar_flux_phcm2s1(738) =   60.7
 lambda_EUV_A(739) =   864.00 ; solar_flux_phcm2s1(739) =   62.2
 lambda_EUV_A(740) =   865.00 ; solar_flux_phcm2s1(740) =   63.9
 lambda_EUV_A(741) =   866.00 ; solar_flux_phcm2s1(741) =   65.5
 lambda_EUV_A(742) =   867.00 ; solar_flux_phcm2s1(742) =   67.1
 lambda_EUV_A(743) =   868.00 ; solar_flux_phcm2s1(743) =   68.9
 lambda_EUV_A(744) =   869.00 ; solar_flux_phcm2s1(744) =   70.6
 lambda_EUV_A(745) =   870.00 ; solar_flux_phcm2s1(745) =   72.4
 lambda_EUV_A(746) =   871.00 ; solar_flux_phcm2s1(746) =   74.2
 lambda_EUV_A(747) =   872.00 ; solar_flux_phcm2s1(747) =   76.1
 lambda_EUV_A(748) =   873.00 ; solar_flux_phcm2s1(748) =   78.1
 lambda_EUV_A(749) =   874.00 ; solar_flux_phcm2s1(749) =   80.1
 lambda_EUV_A(750) =   875.00 ; solar_flux_phcm2s1(750) =   82.1
 lambda_EUV_A(751) =   876.00 ; solar_flux_phcm2s1(751) =   84.1
 lambda_EUV_A(752) =   877.00 ; solar_flux_phcm2s1(752) =   86.2
 lambda_EUV_A(753) =   878.00 ; solar_flux_phcm2s1(753) =   88.5
 lambda_EUV_A(754) =   879.00 ; solar_flux_phcm2s1(754) =   90.7
 lambda_EUV_A(755) =   880.00 ; solar_flux_phcm2s1(755) =   93.0
 lambda_EUV_A(756) =   881.00 ; solar_flux_phcm2s1(756) =   95.4
 lambda_EUV_A(757) =   882.00 ; solar_flux_phcm2s1(757) =   97.9
 lambda_EUV_A(758) =   883.00 ; solar_flux_phcm2s1(758) =  100.4
 lambda_EUV_A(759) =   884.00 ; solar_flux_phcm2s1(759) =  102.9
 lambda_EUV_A(760) =   885.00 ; solar_flux_phcm2s1(760) =  105.5
 lambda_EUV_A(761) =   886.00 ; solar_flux_phcm2s1(761) =  108.1
 lambda_EUV_A(762) =   887.00 ; solar_flux_phcm2s1(762) =  110.8
 lambda_EUV_A(763) =   888.00 ; solar_flux_phcm2s1(763) =  113.7
 lambda_EUV_A(764) =   889.00 ; solar_flux_phcm2s1(764) =  116.6
 lambda_EUV_A(765) =   890.00 ; solar_flux_phcm2s1(765) =  119.6
 lambda_EUV_A(766) =   891.00 ; solar_flux_phcm2s1(766) =  122.6
 lambda_EUV_A(767) =   892.00 ; solar_flux_phcm2s1(767) =  125.7
 lambda_EUV_A(768) =   893.00 ; solar_flux_phcm2s1(768) =  128.8
 lambda_EUV_A(769) =   894.00 ; solar_flux_phcm2s1(769) =  132.1
 lambda_EUV_A(770) =   895.00 ; solar_flux_phcm2s1(770) =  135.4
 lambda_EUV_A(771) =   896.00 ; solar_flux_phcm2s1(771) =  138.9
 lambda_EUV_A(772) =   897.00 ; solar_flux_phcm2s1(772) =  142.4
 lambda_EUV_A(773) =   898.00 ; solar_flux_phcm2s1(773) =  145.9
 lambda_EUV_A(774) =   899.00 ; solar_flux_phcm2s1(774) =  149.7

! 900-950A, EUVAC bin 32
 lambda_EUV_A(775) =   900.00 ; solar_flux_phcm2s1(775) =  153.5
 lambda_EUV_A(776) =   901.00 ; solar_flux_phcm2s1(776) =  157.4
 lambda_EUV_A(777) =   902.00 ; solar_flux_phcm2s1(777) =  161.4
 lambda_EUV_A(778) =   903.00 ; solar_flux_phcm2s1(778) =  165.5
 lambda_EUV_A(779) =   904.00 ; solar_flux_phcm2s1(779) =  169.7
 lambda_EUV_A(780) =   904.10 ; solar_flux_phcm2s1(780) =  116.9
 lambda_EUV_A(781) =   905.00 ; solar_flux_phcm2s1(781) =  174.0
 lambda_EUV_A(782) =   906.00 ; solar_flux_phcm2s1(782) =  178.3
 lambda_EUV_A(783) =   907.00 ; solar_flux_phcm2s1(783) =  182.9
 lambda_EUV_A(784) =   908.00 ; solar_flux_phcm2s1(784) =  187.6
 lambda_EUV_A(785) =   909.00 ; solar_flux_phcm2s1(785) =  192.3
 lambda_EUV_A(786) =   910.00 ; solar_flux_phcm2s1(786) =  197.2
 lambda_EUV_A(787) =   911.00 ; solar_flux_phcm2s1(787) =  202.1
 lambda_EUV_A(788) =   912.00 ; solar_flux_phcm2s1(788) =  207.3
 lambda_EUV_A(789) =   913.00 ; solar_flux_phcm2s1(789) =    2.4
 lambda_EUV_A(790) =   914.00 ; solar_flux_phcm2s1(790) =    2.4
 lambda_EUV_A(791) =   915.00 ; solar_flux_phcm2s1(791) =    2.4
 lambda_EUV_A(792) =   916.00 ; solar_flux_phcm2s1(792) =    2.5
 lambda_EUV_A(793) =   917.00 ; solar_flux_phcm2s1(793) =    2.5
 lambda_EUV_A(794) =   918.00 ; solar_flux_phcm2s1(794) =    2.7
 lambda_EUV_A(795) =   919.00 ; solar_flux_phcm2s1(795) =    2.7
 lambda_EUV_A(796) =   920.00 ; solar_flux_phcm2s1(796) =    2.8
 lambda_EUV_A(797) =   920.96 ; solar_flux_phcm2s1(797) =   71.1
 lambda_EUV_A(798) =   921.00 ; solar_flux_phcm2s1(798) =    2.8
 lambda_EUV_A(799) =   922.00 ; solar_flux_phcm2s1(799) =    2.8
 lambda_EUV_A(800) =   923.00 ; solar_flux_phcm2s1(800) =    3.0
 lambda_EUV_A(801) =   923.15 ; solar_flux_phcm2s1(801) =  101.5
 lambda_EUV_A(802) =   924.00 ; solar_flux_phcm2s1(802) =    3.0
 lambda_EUV_A(803) =   925.00 ; solar_flux_phcm2s1(803) =    3.1
 lambda_EUV_A(804) =   926.00 ; solar_flux_phcm2s1(804) =    3.1
 lambda_EUV_A(805) =   926.20 ; solar_flux_phcm2s1(805) =  139.4
 lambda_EUV_A(806) =   927.00 ; solar_flux_phcm2s1(806) =    3.2
 lambda_EUV_A(807) =   928.00 ; solar_flux_phcm2s1(807) =    3.2
 lambda_EUV_A(808) =   929.00 ; solar_flux_phcm2s1(808) =    3.3
 lambda_EUV_A(809) =   930.00 ; solar_flux_phcm2s1(809) =    3.3
 lambda_EUV_A(810) =   930.75 ; solar_flux_phcm2s1(810) =  164.6
 lambda_EUV_A(811) =   931.00 ; solar_flux_phcm2s1(811) =    3.4
 lambda_EUV_A(812) =   932.00 ; solar_flux_phcm2s1(812) =    3.4
 lambda_EUV_A(813) =   933.00 ; solar_flux_phcm2s1(813) =    3.5
 lambda_EUV_A(814) =   933.38 ; solar_flux_phcm2s1(814) =  115.9
 lambda_EUV_A(815) =   934.00 ; solar_flux_phcm2s1(815) =    3.5
 lambda_EUV_A(816) =   935.00 ; solar_flux_phcm2s1(816) =    3.6
 lambda_EUV_A(817) =   936.00 ; solar_flux_phcm2s1(817) =    3.6
 lambda_EUV_A(818) =   937.00 ; solar_flux_phcm2s1(818) =    3.7
 lambda_EUV_A(819) =   937.80 ; solar_flux_phcm2s1(819) =  227.7
 lambda_EUV_A(820) =   938.00 ; solar_flux_phcm2s1(820) =    3.8
 lambda_EUV_A(821) =   939.00 ; solar_flux_phcm2s1(821) =    3.8
 lambda_EUV_A(822) =   940.00 ; solar_flux_phcm2s1(822) =    3.9
 lambda_EUV_A(823) =   941.00 ; solar_flux_phcm2s1(823) =    3.9
 lambda_EUV_A(824) =   942.00 ; solar_flux_phcm2s1(824) =    4.1
 lambda_EUV_A(825) =   943.00 ; solar_flux_phcm2s1(825) =    4.3
 lambda_EUV_A(826) =   944.00 ; solar_flux_phcm2s1(826) =    4.3
 lambda_EUV_A(827) =   944.52 ; solar_flux_phcm2s1(827) =   76.5
 lambda_EUV_A(828) =   945.00 ; solar_flux_phcm2s1(828) =    4.4
 lambda_EUV_A(829) =   946.00 ; solar_flux_phcm2s1(829) =    4.5
 lambda_EUV_A(830) =   947.00 ; solar_flux_phcm2s1(830) =    4.5
 lambda_EUV_A(831) =   948.00 ; solar_flux_phcm2s1(831) =    4.6
 lambda_EUV_A(832) =   949.00 ; solar_flux_phcm2s1(832) =    4.7
 lambda_EUV_A(833) =   949.74 ; solar_flux_phcm2s1(833) =  378.7

! 950-1000A, EUVAC bin 34
 lambda_EUV_A(834) =   950.00 ; solar_flux_phcm2s1(834) =    4.7
 lambda_EUV_A(835) =   951.00 ; solar_flux_phcm2s1(835) =    4.8
 lambda_EUV_A(836) =   952.00 ; solar_flux_phcm2s1(836) =    4.9
 lambda_EUV_A(837) =   953.00 ; solar_flux_phcm2s1(837) =    5.0
 lambda_EUV_A(838) =   954.00 ; solar_flux_phcm2s1(838) =    5.0
 lambda_EUV_A(839) =   955.00 ; solar_flux_phcm2s1(839) =    5.1
 lambda_EUV_A(840) =   956.00 ; solar_flux_phcm2s1(840) =    5.3
 lambda_EUV_A(841) =   957.00 ; solar_flux_phcm2s1(841) =    5.4
 lambda_EUV_A(842) =   958.00 ; solar_flux_phcm2s1(842) =    5.5
 lambda_EUV_A(843) =   959.00 ; solar_flux_phcm2s1(843) =    5.6
 lambda_EUV_A(844) =   960.00 ; solar_flux_phcm2s1(844) =    5.6
 lambda_EUV_A(845) =   961.00 ; solar_flux_phcm2s1(845) =    5.7
 lambda_EUV_A(846) =   962.00 ; solar_flux_phcm2s1(846) =    5.8
 lambda_EUV_A(847) =   963.00 ; solar_flux_phcm2s1(847) =    5.9
 lambda_EUV_A(848) =   964.00 ; solar_flux_phcm2s1(848) =    6.1
 lambda_EUV_A(849) =   965.00 ; solar_flux_phcm2s1(849) =    6.2
 lambda_EUV_A(850) =   966.00 ; solar_flux_phcm2s1(850) =    6.3
 lambda_EUV_A(851) =   967.00 ; solar_flux_phcm2s1(851) =    6.5
 lambda_EUV_A(852) =   968.00 ; solar_flux_phcm2s1(852) =    6.6
 lambda_EUV_A(853) =   969.00 ; solar_flux_phcm2s1(853) =    6.7
 lambda_EUV_A(854) =   970.00 ; solar_flux_phcm2s1(854) =    6.8
 lambda_EUV_A(855) =   971.00 ; solar_flux_phcm2s1(855) =    6.9
 lambda_EUV_A(856) =   972.00 ; solar_flux_phcm2s1(856) =    7.0
 lambda_EUV_A(857) =   972.54 ; solar_flux_phcm2s1(857) =  754.5
 lambda_EUV_A(858) =   973.00 ; solar_flux_phcm2s1(858) =    7.1
 lambda_EUV_A(859) =   974.00 ; solar_flux_phcm2s1(859) =    7.2
 lambda_EUV_A(860) =   975.00 ; solar_flux_phcm2s1(860) =    7.4
 lambda_EUV_A(861) =   976.00 ; solar_flux_phcm2s1(861) =    7.6
 lambda_EUV_A(862) =   977.00 ; solar_flux_phcm2s1(862) =    7.8
 lambda_EUV_A(863) =   977.02 ; solar_flux_phcm2s1(863) = 4840.0    ! EUVAC bin/line 33
 lambda_EUV_A(864) =   978.00 ; solar_flux_phcm2s1(864) =    7.9
 lambda_EUV_A(865) =   979.00 ; solar_flux_phcm2s1(865) =    8.0
 lambda_EUV_A(866) =   980.00 ; solar_flux_phcm2s1(866) =    8.2
 lambda_EUV_A(867) =   981.00 ; solar_flux_phcm2s1(867) =    8.3
 lambda_EUV_A(868) =   982.00 ; solar_flux_phcm2s1(868) =    8.4
 lambda_EUV_A(869) =   983.00 ; solar_flux_phcm2s1(869) =    8.5
 lambda_EUV_A(870) =   984.00 ; solar_flux_phcm2s1(870) =    8.8
 lambda_EUV_A(871) =   985.00 ; solar_flux_phcm2s1(871) =    8.9
 lambda_EUV_A(872) =   986.00 ; solar_flux_phcm2s1(872) =    9.0
 lambda_EUV_A(873) =   987.00 ; solar_flux_phcm2s1(873) =    9.2
 lambda_EUV_A(874) =   988.00 ; solar_flux_phcm2s1(874) =    9.3
 lambda_EUV_A(875) =   989.00 ; solar_flux_phcm2s1(875) =    9.6
 lambda_EUV_A(876) =   989.79 ; solar_flux_phcm2s1(876) =  191.3
 lambda_EUV_A(877) =   990.00 ; solar_flux_phcm2s1(877) =    9.7
 lambda_EUV_A(878) =   991.00 ; solar_flux_phcm2s1(878) =   10.0
 lambda_EUV_A(879) =   991.55 ; solar_flux_phcm2s1(879) =  382.5
 lambda_EUV_A(880) =   992.00 ; solar_flux_phcm2s1(880) =   10.1
 lambda_EUV_A(881) =   993.00 ; solar_flux_phcm2s1(881) =   10.3
 lambda_EUV_A(882) =   994.00 ; solar_flux_phcm2s1(882) =   10.4
 lambda_EUV_A(883) =   995.00 ; solar_flux_phcm2s1(883) =   10.6
 lambda_EUV_A(884) =   996.00 ; solar_flux_phcm2s1(884) =   10.8
 lambda_EUV_A(885) =   997.00 ; solar_flux_phcm2s1(885) =   10.9
 lambda_EUV_A(886) =   998.00 ; solar_flux_phcm2s1(886) =   11.3
 lambda_EUV_A(887) =   999.00 ; solar_flux_phcm2s1(887) =   11.5

! 1000-1050A, EUVAC bin 37 (actual range 1000-1026A)
 lambda_EUV_A(888) =  1000.00 ; solar_flux_phcm2s1(888) =   11.7
 lambda_EUV_A(889) =  1001.00 ; solar_flux_phcm2s1(889) =   11.8
 lambda_EUV_A(890) =  1002.00 ; solar_flux_phcm2s1(890) =   12.0
 lambda_EUV_A(891) =  1003.00 ; solar_flux_phcm2s1(891) =   12.3
 lambda_EUV_A(892) =  1004.00 ; solar_flux_phcm2s1(892) =   12.5
 lambda_EUV_A(893) =  1005.00 ; solar_flux_phcm2s1(893) =   12.8
 lambda_EUV_A(894) =  1006.00 ; solar_flux_phcm2s1(894) =   13.0
 lambda_EUV_A(895) =  1007.00 ; solar_flux_phcm2s1(895) =   13.2
 lambda_EUV_A(896) =  1008.00 ; solar_flux_phcm2s1(896) =   13.5
 lambda_EUV_A(897) =  1009.00 ; solar_flux_phcm2s1(897) =   13.7
 lambda_EUV_A(898) =  1010.00 ; solar_flux_phcm2s1(898) =   13.9
 lambda_EUV_A(899) =  1010.20 ; solar_flux_phcm2s1(899) =   85.0
 lambda_EUV_A(900) =  1011.00 ; solar_flux_phcm2s1(900) =   14.1
 lambda_EUV_A(901) =  1012.00 ; solar_flux_phcm2s1(901) =   14.4
 lambda_EUV_A(902) =  1013.00 ; solar_flux_phcm2s1(902) =   14.8
 lambda_EUV_A(903) =  1014.00 ; solar_flux_phcm2s1(903) =   15.0
 lambda_EUV_A(904) =  1015.00 ; solar_flux_phcm2s1(904) =   15.3
 lambda_EUV_A(905) =  1016.00 ; solar_flux_phcm2s1(905) =   15.5
 lambda_EUV_A(906) =  1017.00 ; solar_flux_phcm2s1(906) =   15.8
 lambda_EUV_A(907) =  1018.00 ; solar_flux_phcm2s1(907) =   16.1
 lambda_EUV_A(908) =  1019.00 ; solar_flux_phcm2s1(908) =   16.4
 lambda_EUV_A(909) =  1020.00 ; solar_flux_phcm2s1(909) =   16.7
 lambda_EUV_A(910) =  1021.00 ; solar_flux_phcm2s1(910) =   17.0
 lambda_EUV_A(911) =  1022.00 ; solar_flux_phcm2s1(911) =   17.3
 lambda_EUV_A(912) =  1023.00 ; solar_flux_phcm2s1(912) =   17.6
 lambda_EUV_A(913) =  1024.00 ; solar_flux_phcm2s1(913) =   17.8
 lambda_EUV_A(914) =  1025.00 ; solar_flux_phcm2s1(914) =   18.3
 lambda_EUV_A(915) =  1025.72 ; solar_flux_phcm2s1(915) = 4375.0    ! EUVAC bin/line 35
 lambda_EUV_A(916) =  1026.00 ; solar_flux_phcm2s1(916) =   18.6

! the flux above is in units of [10^10 1/m^2/s]
! to convert it to units of [10^9 1/cm^2/s] use 
! 10^10 1/m^2/s = 10*10^9 1/(10^4 cm^2)/s = 0.001*(10^9 1/cm^2/s) ::

 solar_flux_phcm2s1 = 0.001 * solar_flux_phcm2s1  ! now the flux is in units of [10^9 1/cm^2/s]

 solar_flux_phcm2s1(  1:196) = 3.0 * solar_flux_phcm2s1(  1:196) ! below 150A
 solar_flux_phcm2s1(197:289) = 2.0 * solar_flux_phcm2s1(197:289) ! between 150A and 250 A

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
 n2 = 146
 scale_factor = 1.0 + 1.0017e-2 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 100-150A, EUVAC bin 2
 n1 = 147
 n2 = 196
 scale_factor = 1.0 + 7.1250e-3 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 150-200A, EUVAC bin 3
 n1 = 197
 n2 = 241
 scale_factor = 1.0 + 1.3375e-2 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 200-250A, EUVAC bin 4
 n1 = 242
 n2 = 289
 scale_factor = 1.0 + 1.9450e-2 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 250-300A, EUVAC bin 7
 n1 = 290
 n2 = 324
 DO i = n1, n2
    SELECT CASE (i)
       CASE (294)
! 256.32A, EUVAC bin/line 5
          scale_factor = 1.0 + 2.7750e-3 * factor_f107
       CASE (317)
! 284.15A, EUVAC bin/line 6
          scale_factor = 1.0 + 1.3768e-1 * factor_f107
       CASE DEFAULT
! bin 7
          scale_factor = 1.0 + 2.6467e-2 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 300-350A, EUVAC bin 10
 n1 = 325
 n2 = 336
 DO i = n1, n2
    SELECT CASE (i)
       CASE (325)
! 303.31A, EUVAC bin/line 8
          scale_factor = 1.0 + 2.5000e-2 * factor_f107
       CASE (326)
! 303.78A, EUVAC bin/line 9
          scale_factor = 1.0 + 3.3333e-3 * factor_f107
       CASE DEFAULT
! bin 10
          scale_factor = 1.0 + 2.2450e-2 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 350-400A, EUVAC bin 12
 n1 = 337
 n2 = 341
 DO i = n1, n2
    SELECT CASE (i)
       CASE (340)
! 368.07A, EUVAC bin/line 11
          scale_factor = 1.0 + 6.5917e-3 * factor_f107
       CASE DEFAULT
! bin 12
          scale_factor = 1.0 + 3.6542e-2 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 400-450A, EUVAC bin 13
 n1 = 342
 n2 = 374
 scale_factor = 1.0 + 7.4083e-3 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 450-500A, EUVAC bin 15
 n1 = 375
 n2 = 427
 DO i = n1, n2
    SELECT CASE (i)
       CASE (391)
! 465.22A, EUVAC bin/line 14
          scale_factor = 1.0 + 7.4917e-3 * factor_f107
       CASE DEFAULT
! bin 15
          scale_factor = 1.0 + 2.0225e-2 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 500-550A, EUVAC bin 16
 n1 = 428
 n2 = 438
 scale_factor = 1.0 + 8.7583e-3 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 550-600A, EUVAC bin 19
 n1 = 439
 n2 = 448
 DO i = n1, n2
    SELECT CASE (i)
       CASE (440)
! 554.37A, EUVAC bin/line 17
          scale_factor = 1.0 + 3.2667e-3 * factor_f107
       CASE (446)
! 584.33A, EUVAC bin/line 18
          scale_factor = 1.0 + 5.1583e-3 * factor_f107
       CASE DEFAULT
! bin 19
          scale_factor = 1.0 + 3.6583e-3 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 600-650A, EUVAC bin 22
 n1 = 449
 n2 = 499
 DO i = n1, n2
    SELECT CASE (i)
       CASE (451)
! 609.76A, EUVAC bin/line 20
          scale_factor = 1.0 + 1.6175e-2 * factor_f107
       CASE (474)
! 629.73A, EUVAC bin/line 21
          scale_factor = 1.0 + 3.3250e-3 * factor_f107
       CASE DEFAULT
! bin 22
          scale_factor = 1.0 + 1.1800e-2 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 650-700A, EUVAC bin 23
 n1 = 500
 n2 = 557
 scale_factor = 1.0 + 4.2667e-3 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 700-750A, EUVAC bin 25
 n1 = 558
 n2 = 610
 DO i = n1, n2
    SELECT CASE (i)
       CASE (562)
! 703.36A, EUVAC bin/line 24
          scale_factor = 1.0 + 3.0417e-3 * factor_f107
       CASE DEFAULT
! bin 25
          scale_factor = 1.0 + 4.7500e-3 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 750-800A, EUVAC bin 29
 n1 = 611
 n2 = 673
 DO i = n1, n2
    SELECT CASE (i)
       CASE (633)
! 765.15A, EUVAC bin/line 26
          scale_factor = 1.0 + 3.8500e-3 * factor_f107
       CASE (639,651)
! 639 is 770.41A, EUVAC bin/line 27, identified as NE VIII
! 651 is intense line with wavelength 780.32, also identified as NE VIII but not mentioned in EUVAC
          scale_factor = 1.0 + 1.2808e-2 * factor_f107
       CASE (660,664)
! 789.36A, EUVAC bin/line 28, identified as O IV
! apparently, EUVAC wavelength and flux are averaged over multiplet [Heroux and Hinteregger, JGR, 83, 5305, 1978]
! 660 is line with wavelength 787.71 identified as O IV
! 664 is line with wavelength 790.15 identified as O IV
          scale_factor = 1.0 + 3.2750e-3 * factor_f107
       CASE DEFAULT
! bin 29
          scale_factor = 1.0 + 4.7667e-3 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 800-850A, EUVAC bin 30
 n1 = 674
 n2 = 724
 scale_factor = 1.0 + 4.8167e-3 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 850-900A, EUVAC bin 31
 n1 = 725
 n2 = 774
 scale_factor = 1.0 + 5.6750e-3 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 900-950A, EUVAC bin 32
 n1 = 775
 n2 = 833
 scale_factor = 1.0 + 4.9833e-3 * factor_f107
 DO i = n1, n2
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 950-1000A, EUVAC bin 34
 n1 = 834
 n2 = 887
 DO i = n1, n2
    SELECT CASE (i)
       CASE (863)
! 977.02A, EUVAC bin/line 33
          scale_factor = 1.0 + 3.9417e-3 * factor_f107
       CASE DEFAULT
! bin 34
          scale_factor = 1.0 + 4.4167e-3 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! 1000-1050A, EUVAC bin 37 (actual range 1000-1026A)
 n1 = 888
 n2 = 916
 DO i = n1, n2
    SELECT CASE (i)
       CASE (915)
! 1025.72A, EUVAC bin/line 35
          scale_factor = 1.0 + 5.1833e-3 * factor_f107
       CASE DEFAULT
!  bin 37
          scale_factor = 1.0 + 4.3750e-3 * factor_f107
    END SELECT
    solar_flux_phcm2s1(i) = solar_flux_phcm2s1(i) * scale_factor
 END DO

! EUVAC bin/line 36 with wavelength 1031.91A has energy below any ionization threshold which is why it is not considered

! solar flux energy bins

  factor_energy_eVA =  REAL(h_Planck_Js * c_ms * 1.0d10 / e_Cl)   ! 1.0d10 because the wavelength is in Angstrom

  DO i = 1, N_solar_bins
     solar_flux_energy_bin_eV(i) = factor_energy_eVA / lambda_EUV_A(i)
  END DO


END SUBROUTINE Use_f76
