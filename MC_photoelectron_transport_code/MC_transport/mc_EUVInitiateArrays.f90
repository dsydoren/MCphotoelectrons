!-------------------------------------------------
! all cross-sections below are from
! Fennelly, J., and Torr, D. (1992). 
! "Photoionization and photoabsorption cross sections of O, N2, O2, and N for aeronomic calculations". 
! Atomic Data and Nuclear Data Tables, 51 (2), 321-363. doi: https://doi.org/10.1016/0092-640X(92)90004-2
!
SUBROUTINE Initiate_EUV_fluxes_crossections

  USE ParallelOperationValues
  USE PhysicalConstants
  USE Photoelectrons

  IMPLICIT NONE

  INTEGER ALLOC_ERR

! only need variables declared below if use f76 or f74113 spectrum

  REAL, ALLOCATABLE :: lambda_EUV_A(:)                      ! wavelengths of solar EUV fluxes

  integer, parameter :: N_given_xs_points = 1571

  REAL, ALLOCATABLE :: given_lambda_A(:)                    ! wavelengths corresponding to the cross-sections listed below

  REAL, ALLOCATABLE :: given_sigma_tot_O_cm2(:)             ! total cross section [10^-18 cm^2] for atomic oxygen          
  REAL, ALLOCATABLE :: given_sigma_tot_N2_cm2(:)            ! - " -, molecular nitrogen             
  REAL, ALLOCATABLE :: given_sigma_tot_O2_cm2(:)            ! - " -, molecular oxygen

! used to calculate ionization rates
  REAL, ALLOCATABLE :: given_sigma_ion_N2_to_N2_cm2(:)      ! ionization cross section [10^-18 cm^2] for molecular nitrogen, N2+
  REAL, ALLOCATABLE :: given_sigma_ion_N2_to_N_cm2(:)       ! - " -, molecular nitrogen, N+

  REAL, ALLOCATABLE :: given_sigma_ion_O_4S_cm2(:)          ! - " -, atomic oxygen, O+4S
  REAL, ALLOCATABLE :: given_sigma_ion_O_2D_cm2(:)          ! - " -, atomic oxygen, O+2D
  REAL, ALLOCATABLE :: given_sigma_ion_O_2P_cm2(:)          ! - " -, atomic oxygen, O+2P
  REAL, ALLOCATABLE :: given_sigma_ion_O_4Pst_cm2(:)        ! - " -, atomic oxygen, O+4P*
  REAL, ALLOCATABLE :: given_sigma_ion_O_2Pst_cm2(:)        ! - " -, atomic oxygen, O+2P*

  REAL, ALLOCATABLE :: given_sigma_ion_O2_to_O2_cm2(:)      ! - " -, molecular oxygen, O2+
  REAL, ALLOCATABLE :: given_sigma_ion_O2_to_O_cm2(:)       ! - " -, molecular oxygen, O+

  INTEGER i

  INTEGER j
  REAL aj, ajp1

  SELECT CASE (EUV_spectrum_model)
     CASE (EUVAC)
        N_solar_bins = 20 * N_split_EUVAC + 17
     CASE (f76)
        N_solar_bins = 916
     CASE (f74113)
        N_solar_bins = 853
  END SELECT

  ALLOCATE(solar_flux_phcm2s1(      1:N_solar_bins), STAT = ALLOC_ERR)
  ALLOCATE(solar_flux_energy_bin_eV(1:N_solar_bins), STAT = ALLOC_ERR)

  ALLOCATE(sigma_tot_O_cm2( 1:N_solar_bins), STAT = ALLOC_ERR)
  ALLOCATE(sigma_tot_N2_cm2(1:N_solar_bins), STAT = ALLOC_ERR)
  ALLOCATE(sigma_tot_O2_cm2(1:N_solar_bins), STAT = ALLOC_ERR)

  ALLOCATE(sigma_ion_N2_to_N2_cm2(1:N_solar_bins), STAT = ALLOC_ERR)
  ALLOCATE(sigma_ion_N2_to_N_cm2( 1:N_solar_bins), STAT = ALLOC_ERR)

  ALLOCATE(sigma_ion_O_4S_cm2(  1:N_solar_bins), STAT = ALLOC_ERR)
  ALLOCATE(sigma_ion_O_2D_cm2(  1:N_solar_bins), STAT = ALLOC_ERR)
  ALLOCATE(sigma_ion_O_2P_cm2(  1:N_solar_bins), STAT = ALLOC_ERR)
  ALLOCATE(sigma_ion_O_4Pst_cm2(1:N_solar_bins), STAT = ALLOC_ERR)
  ALLOCATE(sigma_ion_O_2Pst_cm2(1:N_solar_bins), STAT = ALLOC_ERR)

  ALLOCATE(sigma_ion_O2_to_O2_cm2(1:N_solar_bins), STAT = ALLOC_ERR)
  ALLOCATE(sigma_ion_O2_to_O_cm2( 1:N_solar_bins), STAT = ALLOC_ERR)

  IF (EUV_spectrum_model.EQ.EUVAC) THEN
     CALL Set_EUVAC   ! sets fluxes, energies, and cross-sections
     RETURN
  END IF

  ALLOCATE(lambda_EUV_A(1:N_solar_bins), STAT = ALLOC_ERR)

  ALLOCATE(given_lambda_A(1:N_given_xs_points), STAT = ALLOC_ERR)                   !       

  ALLOCATE(given_sigma_tot_O_cm2(1:N_given_xs_points), STAT = ALLOC_ERR)            ! total cross section [10^-18 cm^2] for atomic oxygen          
  ALLOCATE(given_sigma_tot_N2_cm2(1:N_given_xs_points), STAT = ALLOC_ERR)           ! - " -, molecular nitrogen             
  ALLOCATE(given_sigma_tot_O2_cm2(1:N_given_xs_points), STAT = ALLOC_ERR)           ! - " -, molecular oxygen

! used to calculate ionization rates
  ALLOCATE(given_sigma_ion_N2_to_N2_cm2(1:N_given_xs_points), STAT = ALLOC_ERR)     ! ionization cross section [10^-18 cm^2] for molecular nitrogen, N2+
  ALLOCATE(given_sigma_ion_N2_to_N_cm2(1:N_given_xs_points), STAT = ALLOC_ERR)      ! - " -, molecular nitrogen, N+

  ALLOCATE(given_sigma_ion_O_4S_cm2(1:N_given_xs_points), STAT = ALLOC_ERR)         ! - " -, atomic oxygen, O+4S
  ALLOCATE(given_sigma_ion_O_2D_cm2(1:N_given_xs_points), STAT = ALLOC_ERR)         ! - " -, atomic oxygen, O+2D
  ALLOCATE(given_sigma_ion_O_2P_cm2(1:N_given_xs_points), STAT = ALLOC_ERR)         ! - " -, atomic oxygen, O+2P
  ALLOCATE(given_sigma_ion_O_4Pst_cm2(1:N_given_xs_points), STAT = ALLOC_ERR)       ! - " -, atomic oxygen, O+4P*
  ALLOCATE(given_sigma_ion_O_2Pst_cm2(1:N_given_xs_points), STAT = ALLOC_ERR)       ! - " -, atomic oxygen, O+2P*

  ALLOCATE(given_sigma_ion_O2_to_O2_cm2(1:N_given_xs_points), STAT = ALLOC_ERR)     ! - " -, molecular oxygen, O2+
  ALLOCATE(given_sigma_ion_O2_to_O_cm2(1:N_given_xs_points), STAT = ALLOC_ERR)      ! - " -, molecular oxygen, O+

  SELECT CASE (EUV_spectrum_model)
     CASE (f76)
        CALL Use_f76(lambda_EUV_A)
     CASE (f74113)
        CALL Use_f74113(lambda_EUV_A)
  END SELECT

! total absorption cross section for atomic oxygen [in 10^-18 cm^2]

 given_lambda_A(   1) =   23.70
 given_lambda_A(   2) =   28.47
 given_lambda_A(   3) =   28.79
 given_lambda_A(   4) =   29.52
 given_lambda_A(   5) =   30.00
 given_lambda_A(   6) =   30.27
 given_lambda_A(   7) =   30.43
 given_lambda_A(   8) =   31.00
 given_lambda_A(   9) =   31.62
 given_lambda_A(  10) =   33.74
 given_lambda_A(  11) =   35.40
 given_lambda_A(  12) =   40.95
 given_lambda_A(  13) =   41.30
 given_lambda_A(  14) =   43.76
 given_lambda_A(  15) =   44.02
 given_lambda_A(  16) =   44.16
 given_lambda_A(  17) =   44.70
 given_lambda_A(  18) =   45.30
 given_lambda_A(  19) =   45.66
 given_lambda_A(  20) =   45.90
 given_lambda_A(  21) =   46.40
 given_lambda_A(  22) =   46.67
 given_lambda_A(  23) =   47.70
 given_lambda_A(  24) =   47.87
 given_lambda_A(  25) =   49.22
 given_lambda_A(  26) =   49.60
 given_lambda_A(  27) =   50.36
 given_lambda_A(  28) =   50.52
 given_lambda_A(  29) =   50.69
 given_lambda_A(  30) =   51.70
 given_lambda_A(  31) =   52.30
 given_lambda_A(  32) =   52.91
 given_lambda_A(  33) =   53.90
 given_lambda_A(  34) =   54.15
 given_lambda_A(  35) =   54.42
 given_lambda_A(  36) =   54.70
 given_lambda_A(  37) =   55.06
 given_lambda_A(  38) =   55.34
 given_lambda_A(  39) =   56.08
 given_lambda_A(  40) =   56.40
 given_lambda_A(  41) =   56.92
 given_lambda_A(  42) =   57.36
 given_lambda_A(  43) =   57.56
 given_lambda_A(  44) =   57.88
 given_lambda_A(  45) =   58.96
 given_lambda_A(  46) =   59.62
 given_lambda_A(  47) =   60.30
 given_lambda_A(  48) =   60.85
 given_lambda_A(  49) =   61.07
 given_lambda_A(  50) =   61.63
 given_lambda_A(  51) =   61.90
 given_lambda_A(  52) =   62.30
 given_lambda_A(  53) =   62.77
 given_lambda_A(  54) =   62.92
 given_lambda_A(  55) =   63.16
 given_lambda_A(  56) =   63.30
 given_lambda_A(  57) =   63.65
 given_lambda_A(  58) =   64.11
 given_lambda_A(  59) =   64.38
 given_lambda_A(  60) =   64.60
 given_lambda_A(  61) =   65.21
 given_lambda_A(  62) =   65.71
 given_lambda_A(  63) =   65.85
 given_lambda_A(  64) =   66.30
 given_lambda_A(  65) =   67.14
 given_lambda_A(  66) =   67.35
 given_lambda_A(  67) =   67.60
 given_lambda_A(  68) =   68.35
 given_lambda_A(  69) =   68.90
 given_lambda_A(  70) =   69.65
 given_lambda_A(  71) =   70.00
 given_lambda_A(  72) =   70.54
 given_lambda_A(  73) =   70.75
 given_lambda_A(  74) =   71.00
 given_lambda_A(  75) =   71.94
 given_lambda_A(  76) =   72.19
 given_lambda_A(  77) =   72.31
 given_lambda_A(  78) =   72.63
 given_lambda_A(  79) =   72.80
 given_lambda_A(  80) =   73.55
 given_lambda_A(  81) =   74.21
 given_lambda_A(  82) =   74.44
 given_lambda_A(  83) =   74.83
 given_lambda_A(  84) =   75.03
 given_lambda_A(  85) =   75.29
 given_lambda_A(  86) =   75.46
 given_lambda_A(  87) =   75.73
 given_lambda_A(  88) =   76.01
 given_lambda_A(  89) =   76.48
 given_lambda_A(  90) =   76.83
 given_lambda_A(  91) =   76.94
 given_lambda_A(  92) =   77.30
 given_lambda_A(  93) =   77.50
 given_lambda_A(  94) =   77.74
 given_lambda_A(  95) =   78.56
 given_lambda_A(  96) =   78.70
 given_lambda_A(  97) =   79.08
 given_lambda_A(  98) =   79.48
 given_lambda_A(  99) =   79.76
 given_lambda_A( 100) =   80.00
 given_lambda_A( 101) =   80.21
 given_lambda_A( 102) =   80.55
 given_lambda_A( 103) =   80.94
 given_lambda_A( 104) =   81.16
 given_lambda_A( 105) =   81.58
 given_lambda_A( 106) =   81.94
 given_lambda_A( 107) =   82.10
 given_lambda_A( 108) =   82.43
 given_lambda_A( 109) =   82.67
 given_lambda_A( 110) =   83.25
 given_lambda_A( 111) =   83.42
 given_lambda_A( 112) =   83.67
 given_lambda_A( 113) =   84.00
 given_lambda_A( 114) =   84.26
 given_lambda_A( 115) =   84.50
 given_lambda_A( 116) =   84.72
 given_lambda_A( 117) =   84.86
 given_lambda_A( 118) =   85.16
 given_lambda_A( 119) =   85.50
 given_lambda_A( 120) =   85.69
 given_lambda_A( 121) =   85.87
 given_lambda_A( 122) =   86.23
 given_lambda_A( 123) =   86.40
 given_lambda_A( 124) =   86.77
 given_lambda_A( 125) =   86.98
 given_lambda_A( 126) =   87.30
 given_lambda_A( 127) =   87.61
 given_lambda_A( 128) =   88.10
 given_lambda_A( 129) =   88.42
 given_lambda_A( 130) =   88.60
 given_lambda_A( 131) =   88.90
 given_lambda_A( 132) =   89.14
 given_lambda_A( 133) =   89.70
 given_lambda_A( 134) =   90.14
 given_lambda_A( 135) =   90.45
 given_lambda_A( 136) =   90.71
 given_lambda_A( 137) =   91.00
 given_lambda_A( 138) =   91.48
 given_lambda_A( 139) =   91.69
 given_lambda_A( 140) =   91.81
 given_lambda_A( 141) =   92.09
 given_lambda_A( 142) =   92.55
 given_lambda_A( 143) =   92.81
 given_lambda_A( 144) =   93.61
 given_lambda_A( 145) =   94.07
 given_lambda_A( 146) =   94.25
 given_lambda_A( 147) =   94.39
 given_lambda_A( 148) =   94.81
 given_lambda_A( 149) =   95.37
 given_lambda_A( 150) =   95.51
 given_lambda_A( 151) =   95.81
 given_lambda_A( 152) =   96.05
 given_lambda_A( 153) =   96.49
 given_lambda_A( 154) =   96.83
 given_lambda_A( 155) =   97.12
 given_lambda_A( 156) =   97.51
 given_lambda_A( 157) =   97.87
 given_lambda_A( 158) =   98.12
 given_lambda_A( 159) =   98.23
 given_lambda_A( 160) =   98.50
 given_lambda_A( 161) =   98.88
 given_lambda_A( 162) =   99.44
 given_lambda_A( 163) =   99.71
 given_lambda_A( 164) =   99.99
 given_lambda_A( 165) =  100.54
 given_lambda_A( 166) =  100.96
 given_lambda_A( 167) =  101.57
 given_lambda_A( 168) =  102.15
 given_lambda_A( 169) =  103.01
 given_lambda_A( 170) =  103.15
 given_lambda_A( 171) =  103.30
 given_lambda_A( 172) =  103.58
 given_lambda_A( 173) =  103.94
 given_lambda_A( 174) =  104.23
 given_lambda_A( 175) =  104.76
 given_lambda_A( 176) =  105.23
 given_lambda_A( 177) =  106.25
 given_lambda_A( 178) =  106.57
 given_lambda_A( 179) =  106.93
 given_lambda_A( 180) =  107.80
 given_lambda_A( 181) =  108.05
 given_lambda_A( 182) =  108.46
 given_lambda_A( 183) =  109.50
 given_lambda_A( 184) =  109.98
 given_lambda_A( 185) =  110.56
 given_lambda_A( 186) =  110.76
 given_lambda_A( 187) =  111.16
 given_lambda_A( 188) =  112.70
 given_lambda_A( 189) =  113.80
 given_lambda_A( 190) =  114.09
 given_lambda_A( 191) =  114.24
 given_lambda_A( 192) =  115.39
 given_lambda_A( 193) =  115.80
 given_lambda_A( 194) =  116.75
 given_lambda_A( 195) =  117.20
 given_lambda_A( 196) =  118.10
 given_lambda_A( 197) =  120.00
 given_lambda_A( 198) =  120.40
 given_lambda_A( 199) =  121.15
 given_lambda_A( 200) =  121.79
 given_lambda_A( 201) =  122.70
 given_lambda_A( 202) =  123.50
 given_lambda_A( 203) =  124.00
 given_lambda_A( 204) =  125.00
 given_lambda_A( 205) =  127.65
 given_lambda_A( 206) =  129.87
 given_lambda_A( 207) =  130.00
 given_lambda_A( 208) =  130.30
 given_lambda_A( 209) =  130.50
 given_lambda_A( 210) =  131.02
 given_lambda_A( 211) =  131.21
 given_lambda_A( 212) =  135.00
 given_lambda_A( 213) =  136.21
 given_lambda_A( 214) =  136.45
 given_lambda_A( 215) =  137.80
 given_lambda_A( 216) =  140.00
 given_lambda_A( 217) =  141.20
 given_lambda_A( 218) =  144.21
 given_lambda_A( 219) =  145.00
 given_lambda_A( 220) =  145.90
 given_lambda_A( 221) =  148.38
 given_lambda_A( 222) =  150.00
 given_lambda_A( 223) =  150.10
 given_lambda_A( 224) =  152.15
 given_lambda_A( 225) =  152.84
 given_lambda_A( 226) =  154.18
 given_lambda_A( 227) =  155.00
 given_lambda_A( 228) =  157.73
 given_lambda_A( 229) =  158.37
 given_lambda_A( 230) =  159.94
 given_lambda_A( 231) =  160.37
 given_lambda_A( 232) =  164.13
 given_lambda_A( 233) =  165.00
 given_lambda_A( 234) =  165.30
 given_lambda_A( 235) =  167.50
 given_lambda_A( 236) =  168.17
 given_lambda_A( 237) =  168.55
 given_lambda_A( 238) =  168.92
 given_lambda_A( 239) =  169.70
 given_lambda_A( 240) =  170.00
 given_lambda_A( 241) =  171.06
 given_lambda_A( 242) =  172.12
 given_lambda_A( 243) =  172.92
 given_lambda_A( 244) =  173.08
 given_lambda_A( 245) =  174.53
 given_lambda_A( 246) =  175.00
 given_lambda_A( 247) =  175.24
 given_lambda_A( 248) =  175.47
 given_lambda_A( 249) =  177.10
 given_lambda_A( 250) =  177.22
 given_lambda_A( 251) =  178.02
 given_lambda_A( 252) =  179.27
 given_lambda_A( 253) =  179.74
 given_lambda_A( 254) =  180.00
 given_lambda_A( 255) =  180.40
 given_lambda_A( 256) =  180.71
 given_lambda_A( 257) =  181.14
 given_lambda_A( 258) =  182.16
 given_lambda_A( 259) =  182.39
 given_lambda_A( 260) =  183.45
 given_lambda_A( 261) =  183.91
 given_lambda_A( 262) =  184.10
 given_lambda_A( 263) =  184.52
 given_lambda_A( 264) =  184.76
 given_lambda_A( 265) =  185.00
 given_lambda_A( 266) =  185.21
 given_lambda_A( 267) =  186.60
 given_lambda_A( 268) =  186.87
 given_lambda_A( 269) =  187.95
 given_lambda_A( 270) =  188.23
 given_lambda_A( 271) =  188.70
 given_lambda_A( 272) =  190.00
 given_lambda_A( 273) =  190.70
 given_lambda_A( 274) =  191.04
 given_lambda_A( 275) =  191.29
 given_lambda_A( 276) =  192.38
 given_lambda_A( 277) =  192.80
 given_lambda_A( 278) =  193.50
 given_lambda_A( 279) =  195.00
 given_lambda_A( 280) =  195.13
 given_lambda_A( 281) =  196.52
 given_lambda_A( 282) =  196.63
 given_lambda_A( 283) =  197.41
 given_lambda_A( 284) =  198.53
 given_lambda_A( 285) =  200.00
 given_lambda_A( 286) =  201.10
 given_lambda_A( 287) =  202.05
 given_lambda_A( 288) =  202.64
 given_lambda_A( 289) =  203.78
 given_lambda_A( 290) =  204.25
 given_lambda_A( 291) =  204.89
 given_lambda_A( 292) =  206.26
 given_lambda_A( 293) =  206.38
 given_lambda_A( 294) =  206.60
 given_lambda_A( 295) =  207.46
 given_lambda_A( 296) =  208.33
 given_lambda_A( 297) =  209.63
 given_lambda_A( 298) =  209.78
 given_lambda_A( 299) =  209.93
 given_lambda_A( 300) =  211.32
 given_lambda_A( 301) =  212.14
 given_lambda_A( 302) =  213.78
 given_lambda_A( 303) =  214.75
 given_lambda_A( 304) =  215.00
 given_lambda_A( 305) =  215.16
 given_lambda_A( 306) =  216.88
 given_lambda_A( 307) =  217.00
 given_lambda_A( 308) =  218.19
 given_lambda_A( 309) =  219.09
 given_lambda_A( 310) =  220.00
 given_lambda_A( 311) =  221.26
 given_lambda_A( 312) =  221.40
 given_lambda_A( 313) =  221.82
 given_lambda_A( 314) =  223.26
 given_lambda_A( 315) =  223.72
 given_lambda_A( 316) =  224.74
 given_lambda_A( 317) =  225.00
 given_lambda_A( 318) =  225.12
 given_lambda_A( 319) =  227.01
 given_lambda_A( 320) =  227.19
 given_lambda_A( 321) =  227.47
 given_lambda_A( 322) =  228.70
 given_lambda_A( 323) =  229.60
 given_lambda_A( 324) =  230.00
 given_lambda_A( 325) =  230.65
 given_lambda_A( 326) =  231.55
 given_lambda_A( 327) =  232.60
 given_lambda_A( 328) =  233.64
 given_lambda_A( 329) =  234.38
 given_lambda_A( 330) =  235.00
 given_lambda_A( 331) =  235.55
 given_lambda_A( 332) =  237.33
 given_lambda_A( 333) =  238.40
 given_lambda_A( 334) =  239.87
 given_lambda_A( 335) =  240.00
 given_lambda_A( 336) =  240.71
 given_lambda_A( 337) =  241.74
 given_lambda_A( 338) =  243.03
 given_lambda_A( 339) =  243.86
 given_lambda_A( 340) =  244.92
 given_lambda_A( 341) =  245.94
 given_lambda_A( 342) =  246.24
 given_lambda_A( 343) =  246.91
 given_lambda_A( 344) =  247.18
 given_lambda_A( 345) =  248.00
 given_lambda_A( 346) =  249.18
 given_lambda_A( 347) =  250.00
 given_lambda_A( 348) =  251.10
 given_lambda_A( 349) =  251.95
 given_lambda_A( 350) =  252.17
 given_lambda_A( 351) =  253.78
 given_lambda_A( 352) =  254.40
 given_lambda_A( 353) =  255.00
 given_lambda_A( 354) =  256.32
 given_lambda_A( 355) =  256.64
 given_lambda_A( 356) =  257.16
 given_lambda_A( 357) =  257.39
 given_lambda_A( 358) =  258.30
 given_lambda_A( 359) =  259.50
 given_lambda_A( 360) =  260.00
 given_lambda_A( 361) =  261.05
 given_lambda_A( 362) =  262.99
 given_lambda_A( 363) =  264.24
 given_lambda_A( 364) =  264.80
 given_lambda_A( 365) =  265.00
 given_lambda_A( 366) =  269.50
 given_lambda_A( 367) =  270.00
 given_lambda_A( 368) =  270.50
 given_lambda_A( 369) =  271.99
 given_lambda_A( 370) =  272.64
 given_lambda_A( 371) =  274.19
 given_lambda_A( 372) =  275.00
 given_lambda_A( 373) =  275.35
 given_lambda_A( 374) =  275.67
 given_lambda_A( 375) =  276.15
 given_lambda_A( 376) =  276.77
 given_lambda_A( 377) =  277.00
 given_lambda_A( 378) =  277.27
 given_lambda_A( 379) =  278.40
 given_lambda_A( 380) =  280.00
 given_lambda_A( 381) =  280.90
 given_lambda_A( 382) =  281.41
 given_lambda_A( 383) =  284.15
 given_lambda_A( 384) =  285.00
 given_lambda_A( 385) =  285.70
 given_lambda_A( 386) =  285.85
 given_lambda_A( 387) =  288.36
 given_lambda_A( 388) =  289.17
 given_lambda_A( 389) =  290.00
 given_lambda_A( 390) =  290.69
 given_lambda_A( 391) =  291.63
 given_lambda_A( 392) =  292.00
 given_lambda_A( 393) =  292.78
 given_lambda_A( 394) =  295.00
 given_lambda_A( 395) =  295.57
 given_lambda_A( 396) =  296.17
 given_lambda_A( 397) =  299.50
 given_lambda_A( 398) =  300.00
 given_lambda_A( 399) =  300.80
 given_lambda_A( 400) =  303.31
 given_lambda_A( 401) =  303.78
 given_lambda_A( 402) =  310.00
 given_lambda_A( 403) =  315.00
 given_lambda_A( 404) =  316.10
 given_lambda_A( 405) =  316.20
 given_lambda_A( 406) =  319.01
 given_lambda_A( 407) =  319.83
 given_lambda_A( 408) =  320.00
 given_lambda_A( 409) =  320.56
 given_lambda_A( 410) =  325.00
 given_lambda_A( 411) =  330.00
 given_lambda_A( 412) =  335.00
 given_lambda_A( 413) =  335.39
 given_lambda_A( 414) =  340.00
 given_lambda_A( 415) =  345.00
 given_lambda_A( 416) =  345.13
 given_lambda_A( 417) =  345.74
 given_lambda_A( 418) =  347.39
 given_lambda_A( 419) =  349.85
 given_lambda_A( 420) =  350.00
 given_lambda_A( 421) =  353.86
 given_lambda_A( 422) =  356.01
 given_lambda_A( 423) =  360.00
 given_lambda_A( 424) =  360.76
 given_lambda_A( 425) =  364.48
 given_lambda_A( 426) =  364.80
 given_lambda_A( 427) =  368.07
 given_lambda_A( 428) =  370.00
 given_lambda_A( 429) =  374.74
 given_lambda_A( 430) =  380.00
 given_lambda_A( 431) =  390.00
 given_lambda_A( 432) =  399.82
 given_lambda_A( 433) =  400.00
 given_lambda_A( 434) =  401.14
 given_lambda_A( 435) =  401.70
 given_lambda_A( 436) =  401.94
 given_lambda_A( 437) =  403.26
 given_lambda_A( 438) =  405.00
 given_lambda_A( 439) =  406.00
 given_lambda_A( 440) =  407.00
 given_lambda_A( 441) =  408.00
 given_lambda_A( 442) =  409.00
 given_lambda_A( 443) =  410.00
 given_lambda_A( 444) =  411.00
 given_lambda_A( 445) =  412.00
 given_lambda_A( 446) =  413.00
 given_lambda_A( 447) =  414.00
 given_lambda_A( 448) =  415.00
 given_lambda_A( 449) =  416.00
 given_lambda_A( 450) =  417.00
 given_lambda_A( 451) =  417.24
 given_lambda_A( 452) =  417.71
 given_lambda_A( 453) =  418.00
 given_lambda_A( 454) =  419.00
 given_lambda_A( 455) =  420.00
 given_lambda_A( 456) =  421.00
 given_lambda_A( 457) =  422.00
 given_lambda_A( 458) =  423.00
 given_lambda_A( 459) =  424.00
 given_lambda_A( 460) =  425.00
 given_lambda_A( 461) =  426.00
 given_lambda_A( 462) =  427.00
 given_lambda_A( 463) =  428.00
 given_lambda_A( 464) =  429.00
 given_lambda_A( 465) =  430.00
 given_lambda_A( 466) =  430.47
 given_lambda_A( 467) =  431.00
 given_lambda_A( 468) =  432.00
 given_lambda_A( 469) =  433.00
 given_lambda_A( 470) =  434.00
 given_lambda_A( 471) =  435.00
 given_lambda_A( 472) =  436.00
 given_lambda_A( 473) =  436.10
 given_lambda_A( 474) =  436.67
 given_lambda_A( 475) =  437.00
 given_lambda_A( 476) =  438.00
 given_lambda_A( 477) =  439.00
 given_lambda_A( 478) =  440.00
 given_lambda_A( 479) =  440.50
 given_lambda_A( 480) =  441.00
 given_lambda_A( 481) =  442.00
 given_lambda_A( 482) =  442.30
 given_lambda_A( 483) =  442.70
 given_lambda_A( 484) =  443.00
 given_lambda_A( 485) =  443.60
 given_lambda_A( 486) =  444.00
 given_lambda_A( 487) =  445.00
 given_lambda_A( 488) =  445.10
 given_lambda_A( 489) =  446.00
 given_lambda_A( 490) =  446.20
 given_lambda_A( 491) =  446.70
 given_lambda_A( 492) =  447.00
 given_lambda_A( 493) =  447.60
 given_lambda_A( 494) =  448.00
 given_lambda_A( 495) =  448.60
 given_lambda_A( 496) =  449.00
 given_lambda_A( 497) =  450.00
 given_lambda_A( 498) =  450.20
 given_lambda_A( 499) =  451.00
 given_lambda_A( 500) =  452.00
 given_lambda_A( 501) =  452.60
 given_lambda_A( 502) =  453.00
 given_lambda_A( 503) =  454.00
 given_lambda_A( 504) =  454.30
 given_lambda_A( 505) =  454.70
 given_lambda_A( 506) =  455.00
 given_lambda_A( 507) =  455.40
 given_lambda_A( 508) =  456.00
 given_lambda_A( 509) =  456.50
 given_lambda_A( 510) =  457.00
 given_lambda_A( 511) =  458.00
 given_lambda_A( 512) =  459.00
 given_lambda_A( 513) =  459.20
 given_lambda_A( 514) =  460.00
 given_lambda_A( 515) =  461.00
 given_lambda_A( 516) =  462.00
 given_lambda_A( 517) =  463.00
 given_lambda_A( 518) =  464.00
 given_lambda_A( 519) =  465.00
 given_lambda_A( 520) =  465.22
 given_lambda_A( 521) =  466.00
 given_lambda_A( 522) =  467.00
 given_lambda_A( 523) =  468.00
 given_lambda_A( 524) =  469.00
 given_lambda_A( 525) =  469.80
 given_lambda_A( 526) =  470.00
 given_lambda_A( 527) =  471.00
 given_lambda_A( 528) =  472.00
 given_lambda_A( 529) =  473.00
 given_lambda_A( 530) =  473.20
 given_lambda_A( 531) =  473.80
 given_lambda_A( 532) =  474.00
 given_lambda_A( 533) =  475.00
 given_lambda_A( 534) =  475.90
 given_lambda_A( 535) =  476.00
 given_lambda_A( 536) =  476.40
 given_lambda_A( 537) =  477.00
 given_lambda_A( 538) =  478.00
 given_lambda_A( 539) =  478.50
 given_lambda_A( 540) =  479.00
 given_lambda_A( 541) =  480.00
 given_lambda_A( 542) =  480.80
 given_lambda_A( 543) =  481.00
 given_lambda_A( 544) =  481.70
 given_lambda_A( 545) =  482.00
 given_lambda_A( 546) =  482.10
 given_lambda_A( 547) =  482.60
 given_lambda_A( 548) =  483.00
 given_lambda_A( 549) =  483.30
 given_lambda_A( 550) =  484.00
 given_lambda_A( 551) =  484.70
 given_lambda_A( 552) =  485.00
 given_lambda_A( 553) =  486.00
 given_lambda_A( 554) =  486.60
 given_lambda_A( 555) =  487.00
 given_lambda_A( 556) =  488.00
 given_lambda_A( 557) =  489.00
 given_lambda_A( 558) =  489.50
 given_lambda_A( 559) =  490.00
 given_lambda_A( 560) =  491.00
 given_lambda_A( 561) =  491.70
 given_lambda_A( 562) =  492.00
 given_lambda_A( 563) =  493.00
 given_lambda_A( 564) =  494.00
 given_lambda_A( 565) =  495.00
 given_lambda_A( 566) =  496.00
 given_lambda_A( 567) =  497.00
 given_lambda_A( 568) =  498.00
 given_lambda_A( 569) =  499.00
 given_lambda_A( 570) =  499.27
 given_lambda_A( 571) =  499.37
 given_lambda_A( 572) =  500.00
 given_lambda_A( 573) =  501.00
 given_lambda_A( 574) =  501.10
 given_lambda_A( 575) =  502.00
 given_lambda_A( 576) =  503.00
 given_lambda_A( 577) =  504.00
 given_lambda_A( 578) =  505.00
 given_lambda_A( 579) =  507.90
 given_lambda_A( 580) =  510.00
 given_lambda_A( 581) =  515.60
 given_lambda_A( 582) =  520.00
 given_lambda_A( 583) =  520.66
 given_lambda_A( 584) =  521.10
 given_lambda_A( 585) =  525.80
 given_lambda_A( 586) =  530.00
 given_lambda_A( 587) =  537.02
 given_lambda_A( 588) =  540.00
 given_lambda_A( 589) =  542.80
 given_lambda_A( 590) =  544.70
 given_lambda_A( 591) =  548.90
 given_lambda_A( 592) =  550.00
 given_lambda_A( 593) =  551.40
 given_lambda_A( 594) =  554.37
 given_lambda_A( 595) =  554.51
 given_lambda_A( 596) =  555.60
 given_lambda_A( 597) =  558.60
 given_lambda_A( 598) =  560.00
 given_lambda_A( 599) =  562.80
 given_lambda_A( 600) =  565.00
 given_lambda_A( 601) =  568.50
 given_lambda_A( 602) =  570.00
 given_lambda_A( 603) =  572.30
 given_lambda_A( 604) =  575.00
 given_lambda_A( 605) =  580.00
 given_lambda_A( 606) =  580.40
 given_lambda_A( 607) =  584.00
 given_lambda_A( 608) =  584.33
 given_lambda_A( 609) =  585.00
 given_lambda_A( 610) =  585.80
 given_lambda_A( 611) =  588.90
 given_lambda_A( 612) =  590.00
 given_lambda_A( 613) =  592.40
 given_lambda_A( 614) =  594.10
 given_lambda_A( 615) =  595.00
 given_lambda_A( 616) =  596.70
 given_lambda_A( 617) =  599.60
 given_lambda_A( 618) =  600.00
 given_lambda_A( 619) =  604.30
 given_lambda_A( 620) =  608.00
 given_lambda_A( 621) =  608.40
 given_lambda_A( 622) =  608.90
 given_lambda_A( 623) =  609.70
 given_lambda_A( 624) =  610.00
 given_lambda_A( 625) =  610.70
 given_lambda_A( 626) =  610.90
 given_lambda_A( 627) =  612.00
 given_lambda_A( 628) =  612.40
 given_lambda_A( 629) =  612.70
 given_lambda_A( 630) =  613.00
 given_lambda_A( 631) =  614.00
 given_lambda_A( 632) =  615.00
 given_lambda_A( 633) =  615.20
 given_lambda_A( 634) =  616.00
 given_lambda_A( 635) =  616.30
 given_lambda_A( 636) =  616.60
 given_lambda_A( 637) =  616.90
 given_lambda_A( 638) =  617.70
 given_lambda_A( 639) =  618.00
 given_lambda_A( 640) =  618.20
 given_lambda_A( 641) =  618.70
 given_lambda_A( 642) =  618.90
 given_lambda_A( 643) =  619.80
 given_lambda_A( 644) =  620.00
 given_lambda_A( 645) =  620.50
 given_lambda_A( 646) =  621.00
 given_lambda_A( 647) =  621.90
 given_lambda_A( 648) =  623.00
 given_lambda_A( 649) =  624.00
 given_lambda_A( 650) =  624.50
 given_lambda_A( 651) =  624.93
 given_lambda_A( 652) =  625.90
 given_lambda_A( 653) =  626.60
 given_lambda_A( 654) =  627.00
 given_lambda_A( 655) =  628.00
 given_lambda_A( 656) =  629.00
 given_lambda_A( 657) =  629.60
 given_lambda_A( 658) =  629.70
 given_lambda_A( 659) =  630.00
 given_lambda_A( 660) =  630.30
 given_lambda_A( 661) =  631.00
 given_lambda_A( 662) =  631.60
 given_lambda_A( 663) =  632.00
 given_lambda_A( 664) =  633.00
 given_lambda_A( 665) =  634.00
 given_lambda_A( 666) =  634.40
 given_lambda_A( 667) =  634.70
 given_lambda_A( 668) =  635.00
 given_lambda_A( 669) =  635.80
 given_lambda_A( 670) =  636.00
 given_lambda_A( 671) =  636.30
 given_lambda_A( 672) =  637.00
 given_lambda_A( 673) =  637.30
 given_lambda_A( 674) =  638.00
 given_lambda_A( 675) =  638.50
 given_lambda_A( 676) =  638.70
 given_lambda_A( 677) =  639.00
 given_lambda_A( 678) =  640.00
 given_lambda_A( 679) =  640.41
 given_lambda_A( 680) =  640.93
 given_lambda_A( 681) =  641.81
 given_lambda_A( 682) =  642.00
 given_lambda_A( 683) =  642.50
 given_lambda_A( 684) =  643.00
 given_lambda_A( 685) =  644.00
 given_lambda_A( 686) =  644.20
 given_lambda_A( 687) =  644.40
 given_lambda_A( 688) =  645.00
 given_lambda_A( 689) =  645.90
 given_lambda_A( 690) =  646.60
 given_lambda_A( 691) =  646.70
 given_lambda_A( 692) =  647.00
 given_lambda_A( 693) =  647.50
 given_lambda_A( 694) =  648.00
 given_lambda_A( 695) =  649.00
 given_lambda_A( 696) =  649.40
 given_lambda_A( 697) =  650.00
 given_lambda_A( 698) =  650.30
 given_lambda_A( 699) =  651.00
 given_lambda_A( 700) =  651.80
 given_lambda_A( 701) =  652.00
 given_lambda_A( 702) =  653.00
 given_lambda_A( 703) =  654.00
 given_lambda_A( 704) =  655.00
 given_lambda_A( 705) =  656.00
 given_lambda_A( 706) =  657.00
 given_lambda_A( 707) =  657.30
 given_lambda_A( 708) =  658.00
 given_lambda_A( 709) =  659.00
 given_lambda_A( 710) =  660.00
 given_lambda_A( 711) =  660.30
 given_lambda_A( 712) =  661.00
 given_lambda_A( 713) =  661.40
 given_lambda_A( 714) =  661.90
 given_lambda_A( 715) =  663.00
 given_lambda_A( 716) =  664.00
 given_lambda_A( 717) =  664.60
 given_lambda_A( 718) =  664.90
 given_lambda_A( 719) =  665.30
 given_lambda_A( 720) =  665.80
 given_lambda_A( 721) =  666.00
 given_lambda_A( 722) =  667.00
 given_lambda_A( 723) =  667.30
 given_lambda_A( 724) =  668.00
 given_lambda_A( 725) =  669.00
 given_lambda_A( 726) =  669.60
 given_lambda_A( 727) =  670.00
 given_lambda_A( 728) =  670.40
 given_lambda_A( 729) =  671.00
 given_lambda_A( 730) =  671.50
 given_lambda_A( 731) =  671.90
 given_lambda_A( 732) =  672.90
 given_lambda_A( 733) =  673.60
 given_lambda_A( 734) =  673.80
 given_lambda_A( 735) =  674.00
 given_lambda_A( 736) =  674.40
 given_lambda_A( 737) =  675.00
 given_lambda_A( 738) =  675.20
 given_lambda_A( 739) =  675.70
 given_lambda_A( 740) =  676.00
 given_lambda_A( 741) =  676.20
 given_lambda_A( 742) =  676.60
 given_lambda_A( 743) =  677.00
 given_lambda_A( 744) =  677.50
 given_lambda_A( 745) =  677.90
 given_lambda_A( 746) =  678.30
 given_lambda_A( 747) =  678.80
 given_lambda_A( 748) =  679.00
 given_lambda_A( 749) =  679.20
 given_lambda_A( 750) =  679.90
 given_lambda_A( 751) =  680.20
 given_lambda_A( 752) =  680.40
 given_lambda_A( 753) =  680.70
 given_lambda_A( 754) =  681.00
 given_lambda_A( 755) =  681.30
 given_lambda_A( 756) =  681.40
 given_lambda_A( 757) =  681.60
 given_lambda_A( 758) =  681.70
 given_lambda_A( 759) =  682.00
 given_lambda_A( 760) =  682.30
 given_lambda_A( 761) =  682.80
 given_lambda_A( 762) =  683.00
 given_lambda_A( 763) =  683.30
 given_lambda_A( 764) =  683.80
 given_lambda_A( 765) =  683.90
 given_lambda_A( 766) =  684.30
 given_lambda_A( 767) =  684.50
 given_lambda_A( 768) =  684.80
 given_lambda_A( 769) =  684.90
 given_lambda_A( 770) =  685.50
 given_lambda_A( 771) =  685.71
 given_lambda_A( 772) =  686.00
 given_lambda_A( 773) =  686.30
 given_lambda_A( 774) =  686.60
 given_lambda_A( 775) =  686.70
 given_lambda_A( 776) =  687.00
 given_lambda_A( 777) =  687.30
 given_lambda_A( 778) =  687.90
 given_lambda_A( 779) =  688.40
 given_lambda_A( 780) =  688.70
 given_lambda_A( 781) =  689.00
 given_lambda_A( 782) =  690.00
 given_lambda_A( 783) =  690.60
 given_lambda_A( 784) =  690.80
 given_lambda_A( 785) =  691.00
 given_lambda_A( 786) =  691.20
 given_lambda_A( 787) =  691.40
 given_lambda_A( 788) =  692.00
 given_lambda_A( 789) =  692.20
 given_lambda_A( 790) =  692.40
 given_lambda_A( 791) =  692.70
 given_lambda_A( 792) =  693.00
 given_lambda_A( 793) =  693.80
 given_lambda_A( 794) =  694.00
 given_lambda_A( 795) =  694.30
 given_lambda_A( 796) =  694.90
 given_lambda_A( 797) =  695.20
 given_lambda_A( 798) =  696.00
 given_lambda_A( 799) =  696.50
 given_lambda_A( 800) =  697.00
 given_lambda_A( 801) =  697.30
 given_lambda_A( 802) =  697.50
 given_lambda_A( 803) =  697.70
 given_lambda_A( 804) =  698.00
 given_lambda_A( 805) =  698.30
 given_lambda_A( 806) =  698.90
 given_lambda_A( 807) =  699.40
 given_lambda_A( 808) =  699.60
 given_lambda_A( 809) =  699.70
 given_lambda_A( 810) =  700.00
 given_lambda_A( 811) =  700.30
 given_lambda_A( 812) =  700.40
 given_lambda_A( 813) =  700.80
 given_lambda_A( 814) =  701.00
 given_lambda_A( 815) =  701.20
 given_lambda_A( 816) =  701.60
 given_lambda_A( 817) =  701.80
 given_lambda_A( 818) =  702.00
 given_lambda_A( 819) =  703.00
 given_lambda_A( 820) =  703.36
 given_lambda_A( 821) =  703.60
 given_lambda_A( 822) =  704.00
 given_lambda_A( 823) =  704.50
 given_lambda_A( 824) =  705.00
 given_lambda_A( 825) =  705.30
 given_lambda_A( 826) =  705.90
 given_lambda_A( 827) =  706.60
 given_lambda_A( 828) =  707.00
 given_lambda_A( 829) =  707.60
 given_lambda_A( 830) =  707.90
 given_lambda_A( 831) =  708.90
 given_lambda_A( 832) =  709.20
 given_lambda_A( 833) =  709.70
 given_lambda_A( 834) =  710.00
 given_lambda_A( 835) =  710.80
 given_lambda_A( 836) =  711.00
 given_lambda_A( 837) =  711.90
 given_lambda_A( 838) =  712.50
 given_lambda_A( 839) =  712.70
 given_lambda_A( 840) =  712.90
 given_lambda_A( 841) =  713.50
 given_lambda_A( 842) =  713.90
 given_lambda_A( 843) =  714.70
 given_lambda_A( 844) =  715.00
 given_lambda_A( 845) =  715.60
 given_lambda_A( 846) =  715.70
 given_lambda_A( 847) =  716.00
 given_lambda_A( 848) =  716.50
 given_lambda_A( 849) =  717.00
 given_lambda_A( 850) =  717.20
 given_lambda_A( 851) =  717.60
 given_lambda_A( 852) =  718.00
 given_lambda_A( 853) =  718.50
 given_lambda_A( 854) =  719.00
 given_lambda_A( 855) =  719.20
 given_lambda_A( 856) =  719.40
 given_lambda_A( 857) =  720.00
 given_lambda_A( 858) =  720.40
 given_lambda_A( 859) =  720.90
 given_lambda_A( 860) =  721.30
 given_lambda_A( 861) =  721.40
 given_lambda_A( 862) =  721.80
 given_lambda_A( 863) =  722.00
 given_lambda_A( 864) =  722.50
 given_lambda_A( 865) =  722.90
 given_lambda_A( 866) =  723.40
 given_lambda_A( 867) =  723.90
 given_lambda_A( 868) =  724.20
 given_lambda_A( 869) =  724.80
 given_lambda_A( 870) =  724.90
 given_lambda_A( 871) =  725.50
 given_lambda_A( 872) =  725.70
 given_lambda_A( 873) =  726.00
 given_lambda_A( 874) =  726.40
 given_lambda_A( 875) =  727.00
 given_lambda_A( 876) =  727.30
 given_lambda_A( 877) =  727.50
 given_lambda_A( 878) =  728.00
 given_lambda_A( 879) =  728.30
 given_lambda_A( 880) =  728.70
 given_lambda_A( 881) =  729.00
 given_lambda_A( 882) =  729.40
 given_lambda_A( 883) =  729.80
 given_lambda_A( 884) =  730.00
 given_lambda_A( 885) =  730.60
 given_lambda_A( 886) =  730.90
 given_lambda_A( 887) =  731.50
 given_lambda_A( 888) =  731.80
 given_lambda_A( 889) =  732.00
 given_lambda_A( 890) =  732.20
 given_lambda_A( 891) =  732.50
 given_lambda_A( 892) =  733.00
 given_lambda_A( 893) =  733.30
 given_lambda_A( 894) =  734.00
 given_lambda_A( 895) =  734.50
 given_lambda_A( 896) =  735.00
 given_lambda_A( 897) =  735.30
 given_lambda_A( 898) =  735.90
 given_lambda_A( 899) =  736.50
 given_lambda_A( 900) =  737.00
 given_lambda_A( 901) =  737.20
 given_lambda_A( 902) =  737.50
 given_lambda_A( 903) =  738.00
 given_lambda_A( 904) =  738.40
 given_lambda_A( 905) =  739.00
 given_lambda_A( 906) =  739.20
 given_lambda_A( 907) =  740.00
 given_lambda_A( 908) =  740.20
 given_lambda_A( 909) =  741.00
 given_lambda_A( 910) =  741.20
 given_lambda_A( 911) =  742.00
 given_lambda_A( 912) =  742.20
 given_lambda_A( 913) =  742.40
 given_lambda_A( 914) =  743.00
 given_lambda_A( 915) =  743.20
 given_lambda_A( 916) =  743.50
 given_lambda_A( 917) =  743.70
 given_lambda_A( 918) =  744.00
 given_lambda_A( 919) =  744.50
 given_lambda_A( 920) =  744.90
 given_lambda_A( 921) =  746.00
 given_lambda_A( 922) =  746.40
 given_lambda_A( 923) =  746.70
 given_lambda_A( 924) =  747.00
 given_lambda_A( 925) =  747.50
 given_lambda_A( 926) =  748.00
 given_lambda_A( 927) =  748.50
 given_lambda_A( 928) =  749.00
 given_lambda_A( 929) =  749.20
 given_lambda_A( 930) =  750.00
 given_lambda_A( 931) =  750.70
 given_lambda_A( 932) =  751.00
 given_lambda_A( 933) =  751.60
 given_lambda_A( 934) =  752.00
 given_lambda_A( 935) =  752.30
 given_lambda_A( 936) =  752.90
 given_lambda_A( 937) =  754.00
 given_lambda_A( 938) =  754.40
 given_lambda_A( 939) =  755.00
 given_lambda_A( 940) =  755.20
 given_lambda_A( 941) =  755.50
 given_lambda_A( 942) =  755.80
 given_lambda_A( 943) =  756.00
 given_lambda_A( 944) =  756.20
 given_lambda_A( 945) =  756.50
 given_lambda_A( 946) =  757.00
 given_lambda_A( 947) =  757.20
 given_lambda_A( 948) =  757.50
 given_lambda_A( 949) =  757.70
 given_lambda_A( 950) =  758.00
 given_lambda_A( 951) =  758.30
 given_lambda_A( 952) =  758.40
 given_lambda_A( 953) =  758.68
 given_lambda_A( 954) =  759.00
 given_lambda_A( 955) =  759.40
 given_lambda_A( 956) =  760.00
 given_lambda_A( 957) =  760.20
 given_lambda_A( 958) =  760.40
 given_lambda_A( 959) =  760.70
 given_lambda_A( 960) =  761.00
 given_lambda_A( 961) =  761.50
 given_lambda_A( 962) =  762.00
 given_lambda_A( 963) =  762.20
 given_lambda_A( 964) =  762.50
 given_lambda_A( 965) =  763.00
 given_lambda_A( 966) =  763.20
 given_lambda_A( 967) =  763.50
 given_lambda_A( 968) =  763.70
 given_lambda_A( 969) =  764.00
 given_lambda_A( 970) =  764.40
 given_lambda_A( 971) =  764.70
 given_lambda_A( 972) =  765.00
 given_lambda_A( 973) =  765.30
 given_lambda_A( 974) =  765.40
 given_lambda_A( 975) =  765.70
 given_lambda_A( 976) =  766.00
 given_lambda_A( 977) =  766.50
 given_lambda_A( 978) =  766.70
 given_lambda_A( 979) =  767.00
 given_lambda_A( 980) =  767.30
 given_lambda_A( 981) =  767.50
 given_lambda_A( 982) =  767.70
 given_lambda_A( 983) =  767.90
 given_lambda_A( 984) =  768.30
 given_lambda_A( 985) =  768.40
 given_lambda_A( 986) =  768.70
 given_lambda_A( 987) =  769.00
 given_lambda_A( 988) =  769.20
 given_lambda_A( 989) =  769.50
 given_lambda_A( 990) =  770.00
 given_lambda_A( 991) =  770.20
 given_lambda_A( 992) =  770.40
 given_lambda_A( 993) =  770.80
 given_lambda_A( 994) =  771.00
 given_lambda_A( 995) =  771.50
 given_lambda_A( 996) =  772.00
 given_lambda_A( 997) =  772.40
 given_lambda_A( 998) =  773.00
 given_lambda_A( 999) =  773.50
 given_lambda_A(1000) =  774.00
 given_lambda_A(1001) =  774.50
 given_lambda_A(1002) =  775.00
 given_lambda_A(1003) =  775.50
 given_lambda_A(1004) =  775.70
 given_lambda_A(1005) =  776.00
 given_lambda_A(1006) =  776.50
 given_lambda_A(1007) =  777.00
 given_lambda_A(1008) =  777.50
 given_lambda_A(1009) =  778.00
 given_lambda_A(1010) =  778.50
 given_lambda_A(1011) =  778.70
 given_lambda_A(1012) =  779.00
 given_lambda_A(1013) =  779.50
 given_lambda_A(1014) =  779.80
 given_lambda_A(1015) =  779.90
 given_lambda_A(1016) =  780.30
 given_lambda_A(1017) =  780.50
 given_lambda_A(1018) =  781.00
 given_lambda_A(1019) =  781.20
 given_lambda_A(1020) =  781.50
 given_lambda_A(1021) =  782.00
 given_lambda_A(1022) =  782.50
 given_lambda_A(1023) =  782.90
 given_lambda_A(1024) =  783.20
 given_lambda_A(1025) =  783.50
 given_lambda_A(1026) =  783.80
 given_lambda_A(1027) =  784.00
 given_lambda_A(1028) =  784.40
 given_lambda_A(1029) =  784.80
 given_lambda_A(1030) =  785.00
 given_lambda_A(1031) =  785.50
 given_lambda_A(1032) =  786.00
 given_lambda_A(1033) =  786.20
 given_lambda_A(1034) =  786.40
 given_lambda_A(1035) =  787.00
 given_lambda_A(1036) =  787.50
 given_lambda_A(1037) =  787.70
 given_lambda_A(1038) =  788.00
 given_lambda_A(1039) =  788.50
 given_lambda_A(1040) =  789.00
 given_lambda_A(1041) =  789.50
 given_lambda_A(1042) =  790.00
 given_lambda_A(1043) =  790.50
 given_lambda_A(1044) =  790.80
 given_lambda_A(1045) =  791.00
 given_lambda_A(1046) =  791.30
 given_lambda_A(1047) =  791.40
 given_lambda_A(1048) =  791.80
 given_lambda_A(1049) =  792.40
 given_lambda_A(1050) =  792.80
 given_lambda_A(1051) =  792.92
 given_lambda_A(1052) =  793.14
 given_lambda_A(1053) =  793.50
 given_lambda_A(1054) =  794.00
 given_lambda_A(1055) =  794.50
 given_lambda_A(1056) =  795.00
 given_lambda_A(1057) =  795.20
 given_lambda_A(1058) =  796.00
 given_lambda_A(1059) =  797.00
 given_lambda_A(1060) =  797.70
 given_lambda_A(1061) =  798.00
 given_lambda_A(1062) =  799.00
 given_lambda_A(1063) =  799.50
 given_lambda_A(1064) =  800.00
 given_lambda_A(1065) =  800.50
 given_lambda_A(1066) =  801.00
 given_lambda_A(1067) =  801.50
 given_lambda_A(1068) =  802.00
 given_lambda_A(1069) =  802.60
 given_lambda_A(1070) =  803.00
 given_lambda_A(1071) =  803.50
 given_lambda_A(1072) =  804.00
 given_lambda_A(1073) =  804.27
 given_lambda_A(1074) =  804.38
 given_lambda_A(1075) =  804.50
 given_lambda_A(1076) =  804.78
 given_lambda_A(1077) =  805.00
 given_lambda_A(1078) =  805.29
 given_lambda_A(1079) =  805.44
 given_lambda_A(1080) =  805.74
 given_lambda_A(1081) =  806.00
 given_lambda_A(1082) =  806.23
 given_lambda_A(1083) =  806.42
 given_lambda_A(1084) =  807.00
 given_lambda_A(1085) =  808.00
 given_lambda_A(1086) =  808.20
 given_lambda_A(1087) =  808.50
 given_lambda_A(1088) =  808.80
 given_lambda_A(1089) =  809.00
 given_lambda_A(1090) =  809.30
 given_lambda_A(1091) =  809.50
 given_lambda_A(1092) =  810.00
 given_lambda_A(1093) =  810.50
 given_lambda_A(1094) =  810.66
 given_lambda_A(1095) =  810.85
 given_lambda_A(1096) =  811.00
 given_lambda_A(1097) =  811.26
 given_lambda_A(1098) =  811.49
 given_lambda_A(1099) =  811.61
 given_lambda_A(1100) =  811.80
 given_lambda_A(1101) =  812.00
 given_lambda_A(1102) =  812.27
 given_lambda_A(1103) =  812.50
 given_lambda_A(1104) =  812.80
 given_lambda_A(1105) =  813.00
 given_lambda_A(1106) =  813.50
 given_lambda_A(1107) =  813.70
 given_lambda_A(1108) =  814.00
 given_lambda_A(1109) =  814.50
 given_lambda_A(1110) =  814.90
 given_lambda_A(1111) =  815.50
 given_lambda_A(1112) =  816.00
 given_lambda_A(1113) =  816.42
 given_lambda_A(1114) =  816.77
 given_lambda_A(1115) =  817.00
 given_lambda_A(1116) =  817.19
 given_lambda_A(1117) =  817.50
 given_lambda_A(1118) =  817.78
 given_lambda_A(1119) =  818.00
 given_lambda_A(1120) =  818.20
 given_lambda_A(1121) =  818.34
 given_lambda_A(1122) =  818.50
 given_lambda_A(1123) =  819.00
 given_lambda_A(1124) =  819.50
 given_lambda_A(1125) =  819.80
 given_lambda_A(1126) =  820.00
 given_lambda_A(1127) =  820.50
 given_lambda_A(1128) =  821.00
 given_lambda_A(1129) =  821.30
 given_lambda_A(1130) =  821.50
 given_lambda_A(1131) =  822.00
 given_lambda_A(1132) =  823.00
 given_lambda_A(1133) =  823.20
 given_lambda_A(1134) =  824.00
 given_lambda_A(1135) =  824.50
 given_lambda_A(1136) =  824.90
 given_lambda_A(1137) =  825.30
 given_lambda_A(1138) =  825.50
 given_lambda_A(1139) =  826.00
 given_lambda_A(1140) =  826.50
 given_lambda_A(1141) =  826.60
 given_lambda_A(1142) =  827.00
 given_lambda_A(1143) =  827.50
 given_lambda_A(1144) =  827.80
 given_lambda_A(1145) =  828.00
 given_lambda_A(1146) =  828.30
 given_lambda_A(1147) =  829.00
 given_lambda_A(1148) =  829.40
 given_lambda_A(1149) =  829.60
 given_lambda_A(1150) =  829.80
 given_lambda_A(1151) =  830.00
 given_lambda_A(1152) =  831.00
 given_lambda_A(1153) =  832.00
 given_lambda_A(1154) =  832.50
 given_lambda_A(1155) =  832.80
 given_lambda_A(1156) =  832.90
 given_lambda_A(1157) =  833.50
 given_lambda_A(1158) =  833.50
 given_lambda_A(1159) =  834.00
 given_lambda_A(1160) =  834.20
 given_lambda_A(1161) =  834.50
 given_lambda_A(1162) =  835.00
 given_lambda_A(1163) =  835.20
 given_lambda_A(1164) =  835.40
 given_lambda_A(1165) =  836.00
 given_lambda_A(1166) =  836.30
 given_lambda_A(1167) =  837.00
 given_lambda_A(1168) =  837.50
 given_lambda_A(1169) =  837.80
 given_lambda_A(1170) =  838.00
 given_lambda_A(1171) =  838.60
 given_lambda_A(1172) =  838.90
 given_lambda_A(1173) =  840.00
 given_lambda_A(1174) =  840.50
 given_lambda_A(1175) =  840.70
 given_lambda_A(1176) =  841.00
 given_lambda_A(1177) =  841.50
 given_lambda_A(1178) =  842.00
 given_lambda_A(1179) =  842.50
 given_lambda_A(1180) =  843.00
 given_lambda_A(1181) =  843.50
 given_lambda_A(1182) =  843.80
 given_lambda_A(1183) =  844.00
 given_lambda_A(1184) =  844.50
 given_lambda_A(1185) =  845.00
 given_lambda_A(1186) =  845.50
 given_lambda_A(1187) =  845.90
 given_lambda_A(1188) =  847.00
 given_lambda_A(1189) =  847.60
 given_lambda_A(1190) =  848.00
 given_lambda_A(1191) =  848.50
 given_lambda_A(1192) =  849.00
 given_lambda_A(1193) =  849.20
 given_lambda_A(1194) =  849.50
 given_lambda_A(1195) =  850.00
 given_lambda_A(1196) =  850.60
 given_lambda_A(1197) =  851.00
 given_lambda_A(1198) =  851.50
 given_lambda_A(1199) =  851.80
 given_lambda_A(1200) =  852.00
 given_lambda_A(1201) =  852.19
 given_lambda_A(1202) =  852.50
 given_lambda_A(1203) =  853.00
 given_lambda_A(1204) =  853.20
 given_lambda_A(1205) =  853.50
 given_lambda_A(1206) =  854.00
 given_lambda_A(1207) =  854.50
 given_lambda_A(1208) =  855.00
 given_lambda_A(1209) =  855.50
 given_lambda_A(1210) =  856.00
 given_lambda_A(1211) =  856.20
 given_lambda_A(1212) =  856.50
 given_lambda_A(1213) =  857.00
 given_lambda_A(1214) =  857.30
 given_lambda_A(1215) =  857.50
 given_lambda_A(1216) =  858.00
 given_lambda_A(1217) =  858.50
 given_lambda_A(1218) =  859.00
 given_lambda_A(1219) =  859.20
 given_lambda_A(1220) =  860.00
 given_lambda_A(1221) =  860.40
 given_lambda_A(1222) =  861.00
 given_lambda_A(1223) =  861.50
 given_lambda_A(1224) =  862.00
 given_lambda_A(1225) =  863.00
 given_lambda_A(1226) =  863.20
 given_lambda_A(1227) =  864.00
 given_lambda_A(1228) =  864.60
 given_lambda_A(1229) =  865.00
 given_lambda_A(1230) =  865.20
 given_lambda_A(1231) =  865.40
 given_lambda_A(1232) =  866.00
 given_lambda_A(1233) =  866.50
 given_lambda_A(1234) =  867.00
 given_lambda_A(1235) =  867.50
 given_lambda_A(1236) =  868.00
 given_lambda_A(1237) =  868.50
 given_lambda_A(1238) =  869.00
 given_lambda_A(1239) =  869.50
 given_lambda_A(1240) =  870.00
 given_lambda_A(1241) =  870.80
 given_lambda_A(1242) =  871.00
 given_lambda_A(1243) =  871.40
 given_lambda_A(1244) =  872.00
 given_lambda_A(1245) =  872.50
 given_lambda_A(1246) =  873.00
 given_lambda_A(1247) =  873.50
 given_lambda_A(1248) =  874.00
 given_lambda_A(1249) =  874.50
 given_lambda_A(1250) =  875.00
 given_lambda_A(1251) =  875.20
 given_lambda_A(1252) =  875.50
 given_lambda_A(1253) =  876.00
 given_lambda_A(1254) =  876.20
 given_lambda_A(1255) =  876.50
 given_lambda_A(1256) =  877.00
 given_lambda_A(1257) =  877.50
 given_lambda_A(1258) =  877.72
 given_lambda_A(1259) =  878.50
 given_lambda_A(1260) =  878.92
 given_lambda_A(1261) =  879.20
 given_lambda_A(1262) =  879.50
 given_lambda_A(1263) =  879.80
 given_lambda_A(1264) =  880.00
 given_lambda_A(1265) =  880.50
 given_lambda_A(1266) =  881.00
 given_lambda_A(1267) =  881.50
 given_lambda_A(1268) =  882.00
 given_lambda_A(1269) =  882.50
 given_lambda_A(1270) =  882.70
 given_lambda_A(1271) =  883.00
 given_lambda_A(1272) =  883.30
 given_lambda_A(1273) =  883.50
 given_lambda_A(1274) =  884.00
 given_lambda_A(1275) =  884.50
 given_lambda_A(1276) =  885.00
 given_lambda_A(1277) =  885.50
 given_lambda_A(1278) =  885.80
 given_lambda_A(1279) =  886.00
 given_lambda_A(1280) =  886.50
 given_lambda_A(1281) =  886.90
 given_lambda_A(1282) =  887.50
 given_lambda_A(1283) =  888.00
 given_lambda_A(1284) =  888.50
 given_lambda_A(1285) =  889.00
 given_lambda_A(1286) =  889.50
 given_lambda_A(1287) =  890.00
 given_lambda_A(1288) =  890.50
 given_lambda_A(1289) =  891.00
 given_lambda_A(1290) =  891.50
 given_lambda_A(1291) =  892.00
 given_lambda_A(1292) =  892.50
 given_lambda_A(1293) =  893.00
 given_lambda_A(1294) =  893.50
 given_lambda_A(1295) =  894.00
 given_lambda_A(1296) =  894.50
 given_lambda_A(1297) =  894.70
 given_lambda_A(1298) =  895.00
 given_lambda_A(1299) =  895.50
 given_lambda_A(1300) =  895.80
 given_lambda_A(1301) =  896.00
 given_lambda_A(1302) =  896.50
 given_lambda_A(1303) =  897.00
 given_lambda_A(1304) =  897.20
 given_lambda_A(1305) =  897.50
 given_lambda_A(1306) =  898.00
 given_lambda_A(1307) =  898.50
 given_lambda_A(1308) =  898.50
 given_lambda_A(1309) =  899.00
 given_lambda_A(1310) =  899.50
 given_lambda_A(1311) =  900.00
 given_lambda_A(1312) =  900.20
 given_lambda_A(1313) =  900.50
 given_lambda_A(1314) =  901.00
 given_lambda_A(1315) =  901.30
 given_lambda_A(1316) =  901.50
 given_lambda_A(1317) =  902.00
 given_lambda_A(1318) =  903.00
 given_lambda_A(1319) =  903.50
 given_lambda_A(1320) =  903.80
 given_lambda_A(1321) =  904.00
 given_lambda_A(1322) =  904.50
 given_lambda_A(1323) =  905.00
 given_lambda_A(1324) =  905.50
 given_lambda_A(1325) =  906.00
 given_lambda_A(1326) =  906.40
 given_lambda_A(1327) =  907.00
 given_lambda_A(1328) =  907.40
 given_lambda_A(1329) =  908.00
 given_lambda_A(1330) =  908.50
 given_lambda_A(1331) =  909.00
 given_lambda_A(1332) =  909.50
 given_lambda_A(1333) =  909.80
 given_lambda_A(1334) =  910.00
 given_lambda_A(1335) =  910.40
 given_lambda_A(1336) =  911.00
 given_lambda_A(1337) =  911.50
 given_lambda_A(1338) =  911.76
 given_lambda_A(1339) =  912.00
 given_lambda_A(1340) =  912.32
 given_lambda_A(1341) =  912.50
 given_lambda_A(1342) =  913.00
 given_lambda_A(1343) =  913.30
 given_lambda_A(1344) =  913.50
 given_lambda_A(1345) =  914.00
 given_lambda_A(1346) =  914.50
 given_lambda_A(1347) =  914.70
 given_lambda_A(1348) =  915.00
 given_lambda_A(1349) =  915.50
 given_lambda_A(1350) =  916.00
 given_lambda_A(1351) =  916.50
 given_lambda_A(1352) =  917.00
 given_lambda_A(1353) =  917.20
 given_lambda_A(1354) =  917.50
 given_lambda_A(1355) =  918.00
 given_lambda_A(1356) =  918.50
 given_lambda_A(1357) =  918.90
 given_lambda_A(1358) =  919.50
 given_lambda_A(1359) =  919.90
 given_lambda_A(1360) =  920.40
 given_lambda_A(1361) =  920.96
 given_lambda_A(1362) =  921.50
 given_lambda_A(1363) =  922.00
 given_lambda_A(1364) =  922.50
 given_lambda_A(1365) =  922.80
 given_lambda_A(1366) =  923.00
 given_lambda_A(1367) =  923.50
 given_lambda_A(1368) =  923.70
 given_lambda_A(1369) =  924.00
 given_lambda_A(1370) =  924.30
 given_lambda_A(1371) =  924.50
 given_lambda_A(1372) =  925.00
 given_lambda_A(1373) =  925.50
 given_lambda_A(1374) =  926.00
 given_lambda_A(1375) =  926.20
 given_lambda_A(1376) =  926.40
 given_lambda_A(1377) =  927.00
 given_lambda_A(1378) =  927.60
 given_lambda_A(1379) =  928.00
 given_lambda_A(1380) =  928.50
 given_lambda_A(1381) =  929.00
 given_lambda_A(1382) =  930.00
 given_lambda_A(1383) =  930.50
 given_lambda_A(1384) =  930.75
 given_lambda_A(1385) =  931.00
 given_lambda_A(1386) =  931.50
 given_lambda_A(1387) =  931.90
 given_lambda_A(1388) =  932.40
 given_lambda_A(1389) =  933.00
 given_lambda_A(1390) =  933.38
 given_lambda_A(1391) =  933.50
 given_lambda_A(1392) =  934.00
 given_lambda_A(1393) =  934.50
 given_lambda_A(1394) =  935.00
 given_lambda_A(1395) =  935.30
 given_lambda_A(1396) =  935.50
 given_lambda_A(1397) =  936.00
 given_lambda_A(1398) =  936.50
 given_lambda_A(1399) =  937.00
 given_lambda_A(1400) =  937.50
 given_lambda_A(1401) =  937.80
 given_lambda_A(1402) =  937.90
 given_lambda_A(1403) =  938.50
 given_lambda_A(1404) =  939.00
 given_lambda_A(1405) =  939.30
 given_lambda_A(1406) =  939.50
 given_lambda_A(1407) =  940.00
 given_lambda_A(1408) =  940.50
 given_lambda_A(1409) =  941.00
 given_lambda_A(1410) =  941.50
 given_lambda_A(1411) =  942.00
 given_lambda_A(1412) =  942.40
 given_lambda_A(1413) =  943.00
 given_lambda_A(1414) =  943.30
 given_lambda_A(1415) =  943.50
 given_lambda_A(1416) =  944.00
 given_lambda_A(1417) =  944.50
 given_lambda_A(1418) =  945.00
 given_lambda_A(1419) =  946.00
 given_lambda_A(1420) =  947.00
 given_lambda_A(1421) =  947.70
 given_lambda_A(1422) =  948.00
 given_lambda_A(1423) =  948.50
 given_lambda_A(1424) =  949.00
 given_lambda_A(1425) =  949.50
 given_lambda_A(1426) =  949.74
 given_lambda_A(1427) =  950.00
 given_lambda_A(1428) =  950.30
 given_lambda_A(1429) =  950.50
 given_lambda_A(1430) =  951.00
 given_lambda_A(1431) =  952.00
 given_lambda_A(1432) =  953.00
 given_lambda_A(1433) =  954.00
 given_lambda_A(1434) =  955.00
 given_lambda_A(1435) =  955.90
 given_lambda_A(1436) =  956.50
 given_lambda_A(1437) =  956.70
 given_lambda_A(1438) =  957.00
 given_lambda_A(1439) =  957.50
 given_lambda_A(1440) =  958.00
 given_lambda_A(1441) =  958.20
 given_lambda_A(1442) =  958.50
 given_lambda_A(1443) =  958.80
 given_lambda_A(1444) =  959.00
 given_lambda_A(1445) =  959.50
 given_lambda_A(1446) =  960.00
 given_lambda_A(1447) =  960.50
 given_lambda_A(1448) =  961.00
 given_lambda_A(1449) =  961.50
 given_lambda_A(1450) =  961.90
 given_lambda_A(1451) =  962.50
 given_lambda_A(1452) =  962.80
 given_lambda_A(1453) =  963.00
 given_lambda_A(1454) =  964.00
 given_lambda_A(1455) =  965.00
 given_lambda_A(1456) =  965.50
 given_lambda_A(1457) =  966.00
 given_lambda_A(1458) =  966.50
 given_lambda_A(1459) =  967.00
 given_lambda_A(1460) =  967.50
 given_lambda_A(1461) =  968.00
 given_lambda_A(1462) =  969.00
 given_lambda_A(1463) =  969.50
 given_lambda_A(1464) =  970.00
 given_lambda_A(1465) =  970.40
 given_lambda_A(1466) =  971.00
 given_lambda_A(1467) =  971.50
 given_lambda_A(1468) =  972.00
 given_lambda_A(1469) =  972.50
 given_lambda_A(1470) =  972.90
 given_lambda_A(1471) =  973.50
 given_lambda_A(1472) =  974.00
 given_lambda_A(1473) =  974.50
 given_lambda_A(1474) =  975.00
 given_lambda_A(1475) =  975.30
 given_lambda_A(1476) =  975.50
 given_lambda_A(1477) =  976.00
 given_lambda_A(1478) =  976.50
 given_lambda_A(1479) =  977.00
 given_lambda_A(1480) =  977.50
 given_lambda_A(1481) =  978.00
 given_lambda_A(1482) =  978.50
 given_lambda_A(1483) =  979.00
 given_lambda_A(1484) =  979.50
 given_lambda_A(1485) =  980.00
 given_lambda_A(1486) =  980.50
 given_lambda_A(1487) =  981.00
 given_lambda_A(1488) =  981.50
 given_lambda_A(1489) =  982.00
 given_lambda_A(1490) =  982.50
 given_lambda_A(1491) =  983.00
 given_lambda_A(1492) =  983.30
 given_lambda_A(1493) =  983.50
 given_lambda_A(1494) =  984.00
 given_lambda_A(1495) =  984.50
 given_lambda_A(1496) =  985.00
 given_lambda_A(1497) =  985.20
 given_lambda_A(1498) =  985.50
 given_lambda_A(1499) =  985.90
 given_lambda_A(1500) =  986.30
 given_lambda_A(1501) =  987.00
 given_lambda_A(1502) =  988.00
 given_lambda_A(1503) =  988.50
 given_lambda_A(1504) =  989.00
 given_lambda_A(1505) =  989.60
 given_lambda_A(1506) =  989.79
 given_lambda_A(1507) =  990.00
 given_lambda_A(1508) =  991.00
 given_lambda_A(1509) =  991.50
 given_lambda_A(1510) =  992.00
 given_lambda_A(1511) =  992.90
 given_lambda_A(1512) =  993.20
 given_lambda_A(1513) =  993.50
 given_lambda_A(1514) =  994.00
 given_lambda_A(1515) =  995.00
 given_lambda_A(1516) =  996.00
 given_lambda_A(1517) =  997.00
 given_lambda_A(1518) =  997.20
 given_lambda_A(1519) =  998.00
 given_lambda_A(1520) =  999.00
 given_lambda_A(1521) = 1000.00
 given_lambda_A(1522) = 1001.00
 given_lambda_A(1523) = 1002.00
 given_lambda_A(1524) = 1003.00
 given_lambda_A(1525) = 1004.00
 given_lambda_A(1526) = 1004.30
 given_lambda_A(1527) = 1004.60
 given_lambda_A(1528) = 1005.00
 given_lambda_A(1529) = 1006.00
 given_lambda_A(1530) = 1006.80
 given_lambda_A(1531) = 1007.00
 given_lambda_A(1532) = 1007.90
 given_lambda_A(1533) = 1009.00
 given_lambda_A(1534) = 1009.40
 given_lambda_A(1535) = 1010.00
 given_lambda_A(1536) = 1010.20
 given_lambda_A(1537) = 1011.00
 given_lambda_A(1538) = 1011.40
 given_lambda_A(1539) = 1012.00
 given_lambda_A(1540) = 1012.30
 given_lambda_A(1541) = 1013.00
 given_lambda_A(1542) = 1013.50
 given_lambda_A(1543) = 1013.90
 given_lambda_A(1544) = 1015.00
 given_lambda_A(1545) = 1015.80
 given_lambda_A(1546) = 1016.00
 given_lambda_A(1547) = 1016.40
 given_lambda_A(1548) = 1016.90
 given_lambda_A(1549) = 1017.20
 given_lambda_A(1550) = 1017.80
 given_lambda_A(1551) = 1018.00
 given_lambda_A(1552) = 1018.30
 given_lambda_A(1553) = 1018.80
 given_lambda_A(1554) = 1019.00
 given_lambda_A(1555) = 1019.40
 given_lambda_A(1556) = 1020.00
 given_lambda_A(1557) = 1020.40
 given_lambda_A(1558) = 1020.80
 given_lambda_A(1559) = 1021.00
 given_lambda_A(1560) = 1021.60
 given_lambda_A(1561) = 1022.00
 given_lambda_A(1562) = 1022.40
 given_lambda_A(1563) = 1023.00
 given_lambda_A(1564) = 1023.40
 given_lambda_A(1565) = 1024.00
 given_lambda_A(1566) = 1024.60
 given_lambda_A(1567) = 1025.00
 given_lambda_A(1568) = 1025.30
 given_lambda_A(1569) = 1025.70
 given_lambda_A(1570) = 1026.00
 given_lambda_A(1571) = 1027.00






  given_sigma_ion_O_4S_cm2(   1) =  0.01 ; given_sigma_ion_O_2D_cm2(   1) =    0.01 ; given_sigma_ion_O_2P_cm2(   1) =    0.01 ; given_sigma_ion_O_4Pst_cm2(   1) =    0.00 ; given_sigma_ion_O_2Pst_cm2(   1) =    0.00
  given_sigma_ion_O_4S_cm2(   2) =  0.02 ; given_sigma_ion_O_2D_cm2(   2) =    0.02 ; given_sigma_ion_O_2P_cm2(   2) =    0.01 ; given_sigma_ion_O_4Pst_cm2(   2) =    0.01 ; given_sigma_ion_O_2Pst_cm2(   2) =    0.00
  given_sigma_ion_O_4S_cm2(   3) =  0.02 ; given_sigma_ion_O_2D_cm2(   3) =    0.02 ; given_sigma_ion_O_2P_cm2(   3) =    0.01 ; given_sigma_ion_O_4Pst_cm2(   3) =    0.01 ; given_sigma_ion_O_2Pst_cm2(   3) =    0.00
  given_sigma_ion_O_4S_cm2(   4) =  0.02 ; given_sigma_ion_O_2D_cm2(   4) =    0.02 ; given_sigma_ion_O_2P_cm2(   4) =    0.02 ; given_sigma_ion_O_4Pst_cm2(   4) =    0.01 ; given_sigma_ion_O_2Pst_cm2(   4) =    0.00
  given_sigma_ion_O_4S_cm2(   5) =  0.02 ; given_sigma_ion_O_2D_cm2(   5) =    0.02 ; given_sigma_ion_O_2P_cm2(   5) =    0.02 ; given_sigma_ion_O_4Pst_cm2(   5) =    0.01 ; given_sigma_ion_O_2Pst_cm2(   5) =    0.01
  given_sigma_ion_O_4S_cm2(   6) =  0.02 ; given_sigma_ion_O_2D_cm2(   6) =    0.02 ; given_sigma_ion_O_2P_cm2(   6) =    0.02 ; given_sigma_ion_O_4Pst_cm2(   6) =    0.01 ; given_sigma_ion_O_2Pst_cm2(   6) =    0.01
  given_sigma_ion_O_4S_cm2(   7) =  0.02 ; given_sigma_ion_O_2D_cm2(   7) =    0.03 ; given_sigma_ion_O_2P_cm2(   7) =    0.02 ; given_sigma_ion_O_4Pst_cm2(   7) =    0.01 ; given_sigma_ion_O_2Pst_cm2(   7) =    0.01
  given_sigma_ion_O_4S_cm2(   8) =  0.02 ; given_sigma_ion_O_2D_cm2(   8) =    0.03 ; given_sigma_ion_O_2P_cm2(   8) =    0.02 ; given_sigma_ion_O_4Pst_cm2(   8) =    0.01 ; given_sigma_ion_O_2Pst_cm2(   8) =    0.01
  given_sigma_ion_O_4S_cm2(   9) =  0.02 ; given_sigma_ion_O_2D_cm2(   9) =    0.03 ; given_sigma_ion_O_2P_cm2(   9) =    0.02 ; given_sigma_ion_O_4Pst_cm2(   9) =    0.01 ; given_sigma_ion_O_2Pst_cm2(   9) =    0.01
  given_sigma_ion_O_4S_cm2(  10) =  0.03 ; given_sigma_ion_O_2D_cm2(  10) =    0.03 ; given_sigma_ion_O_2P_cm2(  10) =    0.02 ; given_sigma_ion_O_4Pst_cm2(  10) =    0.01 ; given_sigma_ion_O_2Pst_cm2(  10) =    0.01
  given_sigma_ion_O_4S_cm2(  11) =  0.03 ; given_sigma_ion_O_2D_cm2(  11) =    0.03 ; given_sigma_ion_O_2P_cm2(  11) =    0.02 ; given_sigma_ion_O_4Pst_cm2(  11) =    0.01 ; given_sigma_ion_O_2Pst_cm2(  11) =    0.01
  given_sigma_ion_O_4S_cm2(  12) =  0.04 ; given_sigma_ion_O_2D_cm2(  12) =    0.05 ; given_sigma_ion_O_2P_cm2(  12) =    0.03 ; given_sigma_ion_O_4Pst_cm2(  12) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  12) =    0.01
  given_sigma_ion_O_4S_cm2(  13) =  0.04 ; given_sigma_ion_O_2D_cm2(  13) =    0.05 ; given_sigma_ion_O_2P_cm2(  13) =    0.03 ; given_sigma_ion_O_4Pst_cm2(  13) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  13) =    0.01
  given_sigma_ion_O_4S_cm2(  14) =  0.05 ; given_sigma_ion_O_2D_cm2(  14) =    0.05 ; given_sigma_ion_O_2P_cm2(  14) =    0.03 ; given_sigma_ion_O_4Pst_cm2(  14) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  14) =    0.01
  given_sigma_ion_O_4S_cm2(  15) =  0.05 ; given_sigma_ion_O_2D_cm2(  15) =    0.05 ; given_sigma_ion_O_2P_cm2(  15) =    0.03 ; given_sigma_ion_O_4Pst_cm2(  15) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  15) =    0.01
  given_sigma_ion_O_4S_cm2(  16) =  0.05 ; given_sigma_ion_O_2D_cm2(  16) =    0.05 ; given_sigma_ion_O_2P_cm2(  16) =    0.03 ; given_sigma_ion_O_4Pst_cm2(  16) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  16) =    0.01
  given_sigma_ion_O_4S_cm2(  17) =  0.05 ; given_sigma_ion_O_2D_cm2(  17) =    0.05 ; given_sigma_ion_O_2P_cm2(  17) =    0.03 ; given_sigma_ion_O_4Pst_cm2(  17) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  17) =    0.01
  given_sigma_ion_O_4S_cm2(  18) =  0.05 ; given_sigma_ion_O_2D_cm2(  18) =    0.05 ; given_sigma_ion_O_2P_cm2(  18) =    0.03 ; given_sigma_ion_O_4Pst_cm2(  18) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  18) =    0.01
  given_sigma_ion_O_4S_cm2(  19) =  0.05 ; given_sigma_ion_O_2D_cm2(  19) =    0.05 ; given_sigma_ion_O_2P_cm2(  19) =    0.03 ; given_sigma_ion_O_4Pst_cm2(  19) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  19) =    0.01
  given_sigma_ion_O_4S_cm2(  20) =  0.05 ; given_sigma_ion_O_2D_cm2(  20) =    0.05 ; given_sigma_ion_O_2P_cm2(  20) =    0.03 ; given_sigma_ion_O_4Pst_cm2(  20) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  20) =    0.01
  given_sigma_ion_O_4S_cm2(  21) =  0.05 ; given_sigma_ion_O_2D_cm2(  21) =    0.05 ; given_sigma_ion_O_2P_cm2(  21) =    0.03 ; given_sigma_ion_O_4Pst_cm2(  21) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  21) =    0.01
  given_sigma_ion_O_4S_cm2(  22) =  0.05 ; given_sigma_ion_O_2D_cm2(  22) =    0.05 ; given_sigma_ion_O_2P_cm2(  22) =    0.04 ; given_sigma_ion_O_4Pst_cm2(  22) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  22) =    0.01
  given_sigma_ion_O_4S_cm2(  23) =  0.05 ; given_sigma_ion_O_2D_cm2(  23) =    0.05 ; given_sigma_ion_O_2P_cm2(  23) =    0.04 ; given_sigma_ion_O_4Pst_cm2(  23) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  23) =    0.01
  given_sigma_ion_O_4S_cm2(  24) =  0.05 ; given_sigma_ion_O_2D_cm2(  24) =    0.05 ; given_sigma_ion_O_2P_cm2(  24) =    0.04 ; given_sigma_ion_O_4Pst_cm2(  24) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  24) =    0.01
  given_sigma_ion_O_4S_cm2(  25) =  0.06 ; given_sigma_ion_O_2D_cm2(  25) =    0.06 ; given_sigma_ion_O_2P_cm2(  25) =    0.04 ; given_sigma_ion_O_4Pst_cm2(  25) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  25) =    0.01
  given_sigma_ion_O_4S_cm2(  26) =  0.05 ; given_sigma_ion_O_2D_cm2(  26) =    0.06 ; given_sigma_ion_O_2P_cm2(  26) =    0.04 ; given_sigma_ion_O_4Pst_cm2(  26) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  26) =    0.01
  given_sigma_ion_O_4S_cm2(  27) =  0.05 ; given_sigma_ion_O_2D_cm2(  27) =    0.06 ; given_sigma_ion_O_2P_cm2(  27) =    0.04 ; given_sigma_ion_O_4Pst_cm2(  27) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  27) =    0.01
  given_sigma_ion_O_4S_cm2(  28) =  0.05 ; given_sigma_ion_O_2D_cm2(  28) =    0.06 ; given_sigma_ion_O_2P_cm2(  28) =    0.04 ; given_sigma_ion_O_4Pst_cm2(  28) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  28) =    0.01
  given_sigma_ion_O_4S_cm2(  29) =  0.06 ; given_sigma_ion_O_2D_cm2(  29) =    0.06 ; given_sigma_ion_O_2P_cm2(  29) =    0.04 ; given_sigma_ion_O_4Pst_cm2(  29) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  29) =    0.01
  given_sigma_ion_O_4S_cm2(  30) =  0.06 ; given_sigma_ion_O_2D_cm2(  30) =    0.07 ; given_sigma_ion_O_2P_cm2(  30) =    0.04 ; given_sigma_ion_O_4Pst_cm2(  30) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  30) =    0.02
  given_sigma_ion_O_4S_cm2(  31) =  0.06 ; given_sigma_ion_O_2D_cm2(  31) =    0.07 ; given_sigma_ion_O_2P_cm2(  31) =    0.04 ; given_sigma_ion_O_4Pst_cm2(  31) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  31) =    0.02
  given_sigma_ion_O_4S_cm2(  32) =  0.06 ; given_sigma_ion_O_2D_cm2(  32) =    0.07 ; given_sigma_ion_O_2P_cm2(  32) =    0.04 ; given_sigma_ion_O_4Pst_cm2(  32) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  32) =    0.02
  given_sigma_ion_O_4S_cm2(  33) =  0.06 ; given_sigma_ion_O_2D_cm2(  33) =    0.07 ; given_sigma_ion_O_2P_cm2(  33) =    0.04 ; given_sigma_ion_O_4Pst_cm2(  33) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  33) =    0.02
  given_sigma_ion_O_4S_cm2(  34) =  0.06 ; given_sigma_ion_O_2D_cm2(  34) =    0.08 ; given_sigma_ion_O_2P_cm2(  34) =    0.04 ; given_sigma_ion_O_4Pst_cm2(  34) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  34) =    0.02
  given_sigma_ion_O_4S_cm2(  35) =  0.07 ; given_sigma_ion_O_2D_cm2(  35) =    0.08 ; given_sigma_ion_O_2P_cm2(  35) =    0.05 ; given_sigma_ion_O_4Pst_cm2(  35) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  35) =    0.02
  given_sigma_ion_O_4S_cm2(  36) =  0.06 ; given_sigma_ion_O_2D_cm2(  36) =    0.08 ; given_sigma_ion_O_2P_cm2(  36) =    0.05 ; given_sigma_ion_O_4Pst_cm2(  36) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  36) =    0.02
  given_sigma_ion_O_4S_cm2(  37) =  0.07 ; given_sigma_ion_O_2D_cm2(  37) =    0.08 ; given_sigma_ion_O_2P_cm2(  37) =    0.05 ; given_sigma_ion_O_4Pst_cm2(  37) =    0.02 ; given_sigma_ion_O_2Pst_cm2(  37) =    0.02
  given_sigma_ion_O_4S_cm2(  38) =  0.07 ; given_sigma_ion_O_2D_cm2(  38) =    0.08 ; given_sigma_ion_O_2P_cm2(  38) =    0.05 ; given_sigma_ion_O_4Pst_cm2(  38) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  38) =    0.02
  given_sigma_ion_O_4S_cm2(  39) =  0.07 ; given_sigma_ion_O_2D_cm2(  39) =    0.08 ; given_sigma_ion_O_2P_cm2(  39) =    0.05 ; given_sigma_ion_O_4Pst_cm2(  39) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  39) =    0.02
  given_sigma_ion_O_4S_cm2(  40) =  0.07 ; given_sigma_ion_O_2D_cm2(  40) =    0.08 ; given_sigma_ion_O_2P_cm2(  40) =    0.05 ; given_sigma_ion_O_4Pst_cm2(  40) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  40) =    0.02
  given_sigma_ion_O_4S_cm2(  41) =  0.07 ; given_sigma_ion_O_2D_cm2(  41) =    0.09 ; given_sigma_ion_O_2P_cm2(  41) =    0.05 ; given_sigma_ion_O_4Pst_cm2(  41) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  41) =    0.02
  given_sigma_ion_O_4S_cm2(  42) =  0.07 ; given_sigma_ion_O_2D_cm2(  42) =    0.09 ; given_sigma_ion_O_2P_cm2(  42) =    0.05 ; given_sigma_ion_O_4Pst_cm2(  42) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  42) =    0.02
  given_sigma_ion_O_4S_cm2(  43) =  0.07 ; given_sigma_ion_O_2D_cm2(  43) =    0.09 ; given_sigma_ion_O_2P_cm2(  43) =    0.06 ; given_sigma_ion_O_4Pst_cm2(  43) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  43) =    0.02
  given_sigma_ion_O_4S_cm2(  44) =  0.08 ; given_sigma_ion_O_2D_cm2(  44) =    0.09 ; given_sigma_ion_O_2P_cm2(  44) =    0.05 ; given_sigma_ion_O_4Pst_cm2(  44) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  44) =    0.02
  given_sigma_ion_O_4S_cm2(  45) =  0.08 ; given_sigma_ion_O_2D_cm2(  45) =    0.10 ; given_sigma_ion_O_2P_cm2(  45) =    0.06 ; given_sigma_ion_O_4Pst_cm2(  45) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  45) =    0.02
  given_sigma_ion_O_4S_cm2(  46) =  0.08 ; given_sigma_ion_O_2D_cm2(  46) =    0.10 ; given_sigma_ion_O_2P_cm2(  46) =    0.06 ; given_sigma_ion_O_4Pst_cm2(  46) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  46) =    0.02
  given_sigma_ion_O_4S_cm2(  47) =  0.09 ; given_sigma_ion_O_2D_cm2(  47) =    0.10 ; given_sigma_ion_O_2P_cm2(  47) =    0.06 ; given_sigma_ion_O_4Pst_cm2(  47) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  47) =    0.02
  given_sigma_ion_O_4S_cm2(  48) =  0.09 ; given_sigma_ion_O_2D_cm2(  48) =    0.10 ; given_sigma_ion_O_2P_cm2(  48) =    0.06 ; given_sigma_ion_O_4Pst_cm2(  48) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  48) =    0.02
  given_sigma_ion_O_4S_cm2(  49) =  0.09 ; given_sigma_ion_O_2D_cm2(  49) =    0.11 ; given_sigma_ion_O_2P_cm2(  49) =    0.06 ; given_sigma_ion_O_4Pst_cm2(  49) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  49) =    0.02
  given_sigma_ion_O_4S_cm2(  50) =  0.09 ; given_sigma_ion_O_2D_cm2(  50) =    0.11 ; given_sigma_ion_O_2P_cm2(  50) =    0.07 ; given_sigma_ion_O_4Pst_cm2(  50) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  50) =    0.03
  given_sigma_ion_O_4S_cm2(  51) =  0.10 ; given_sigma_ion_O_2D_cm2(  51) =    0.11 ; given_sigma_ion_O_2P_cm2(  51) =    0.07 ; given_sigma_ion_O_4Pst_cm2(  51) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  51) =    0.03
  given_sigma_ion_O_4S_cm2(  52) =  0.10 ; given_sigma_ion_O_2D_cm2(  52) =    0.11 ; given_sigma_ion_O_2P_cm2(  52) =    0.07 ; given_sigma_ion_O_4Pst_cm2(  52) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  52) =    0.03
  given_sigma_ion_O_4S_cm2(  53) =  0.10 ; given_sigma_ion_O_2D_cm2(  53) =    0.11 ; given_sigma_ion_O_2P_cm2(  53) =    0.07 ; given_sigma_ion_O_4Pst_cm2(  53) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  53) =    0.03
  given_sigma_ion_O_4S_cm2(  54) =  0.10 ; given_sigma_ion_O_2D_cm2(  54) =    0.11 ; given_sigma_ion_O_2P_cm2(  54) =    0.07 ; given_sigma_ion_O_4Pst_cm2(  54) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  54) =    0.03
  given_sigma_ion_O_4S_cm2(  55) =  0.10 ; given_sigma_ion_O_2D_cm2(  55) =    0.11 ; given_sigma_ion_O_2P_cm2(  55) =    0.07 ; given_sigma_ion_O_4Pst_cm2(  55) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  55) =    0.03
  given_sigma_ion_O_4S_cm2(  56) =  0.10 ; given_sigma_ion_O_2D_cm2(  56) =    0.12 ; given_sigma_ion_O_2P_cm2(  56) =    0.07 ; given_sigma_ion_O_4Pst_cm2(  56) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  56) =    0.03
  given_sigma_ion_O_4S_cm2(  57) =  0.10 ; given_sigma_ion_O_2D_cm2(  57) =    0.12 ; given_sigma_ion_O_2P_cm2(  57) =    0.08 ; given_sigma_ion_O_4Pst_cm2(  57) =    0.04 ; given_sigma_ion_O_2Pst_cm2(  57) =    0.03
  given_sigma_ion_O_4S_cm2(  58) =  0.11 ; given_sigma_ion_O_2D_cm2(  58) =    0.12 ; given_sigma_ion_O_2P_cm2(  58) =    0.08 ; given_sigma_ion_O_4Pst_cm2(  58) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  58) =    0.03
  given_sigma_ion_O_4S_cm2(  59) =  0.11 ; given_sigma_ion_O_2D_cm2(  59) =    0.12 ; given_sigma_ion_O_2P_cm2(  59) =    0.07 ; given_sigma_ion_O_4Pst_cm2(  59) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  59) =    0.03
  given_sigma_ion_O_4S_cm2(  60) =  0.11 ; given_sigma_ion_O_2D_cm2(  60) =    0.12 ; given_sigma_ion_O_2P_cm2(  60) =    0.07 ; given_sigma_ion_O_4Pst_cm2(  60) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  60) =    0.03
  given_sigma_ion_O_4S_cm2(  61) =  0.11 ; given_sigma_ion_O_2D_cm2(  61) =    0.12 ; given_sigma_ion_O_2P_cm2(  61) =    0.08 ; given_sigma_ion_O_4Pst_cm2(  61) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  61) =    0.03
  given_sigma_ion_O_4S_cm2(  62) =  0.12 ; given_sigma_ion_O_2D_cm2(  62) =    0.13 ; given_sigma_ion_O_2P_cm2(  62) =    0.08 ; given_sigma_ion_O_4Pst_cm2(  62) =    0.03 ; given_sigma_ion_O_2Pst_cm2(  62) =    0.03
  given_sigma_ion_O_4S_cm2(  63) =  0.12 ; given_sigma_ion_O_2D_cm2(  63) =    0.13 ; given_sigma_ion_O_2P_cm2(  63) =    0.08 ; given_sigma_ion_O_4Pst_cm2(  63) =    0.04 ; given_sigma_ion_O_2Pst_cm2(  63) =    0.03
  given_sigma_ion_O_4S_cm2(  64) =  0.12 ; given_sigma_ion_O_2D_cm2(  64) =    0.13 ; given_sigma_ion_O_2P_cm2(  64) =    0.08 ; given_sigma_ion_O_4Pst_cm2(  64) =    0.04 ; given_sigma_ion_O_2Pst_cm2(  64) =    0.03
  given_sigma_ion_O_4S_cm2(  65) =  0.12 ; given_sigma_ion_O_2D_cm2(  65) =    0.14 ; given_sigma_ion_O_2P_cm2(  65) =    0.09 ; given_sigma_ion_O_4Pst_cm2(  65) =    0.04 ; given_sigma_ion_O_2Pst_cm2(  65) =    0.03
  given_sigma_ion_O_4S_cm2(  66) =  0.12 ; given_sigma_ion_O_2D_cm2(  66) =    0.13 ; given_sigma_ion_O_2P_cm2(  66) =    0.09 ; given_sigma_ion_O_4Pst_cm2(  66) =    0.04 ; given_sigma_ion_O_2Pst_cm2(  66) =    0.03
  given_sigma_ion_O_4S_cm2(  67) =  0.12 ; given_sigma_ion_O_2D_cm2(  67) =    0.14 ; given_sigma_ion_O_2P_cm2(  67) =    0.09 ; given_sigma_ion_O_4Pst_cm2(  67) =    0.04 ; given_sigma_ion_O_2Pst_cm2(  67) =    0.03
  given_sigma_ion_O_4S_cm2(  68) =  0.13 ; given_sigma_ion_O_2D_cm2(  68) =    0.14 ; given_sigma_ion_O_2P_cm2(  68) =    0.09 ; given_sigma_ion_O_4Pst_cm2(  68) =    0.04 ; given_sigma_ion_O_2Pst_cm2(  68) =    0.04
  given_sigma_ion_O_4S_cm2(  69) =  0.13 ; given_sigma_ion_O_2D_cm2(  69) =    0.15 ; given_sigma_ion_O_2P_cm2(  69) =    0.09 ; given_sigma_ion_O_4Pst_cm2(  69) =    0.04 ; given_sigma_ion_O_2Pst_cm2(  69) =    0.03
  given_sigma_ion_O_4S_cm2(  70) =  0.13 ; given_sigma_ion_O_2D_cm2(  70) =    0.15 ; given_sigma_ion_O_2P_cm2(  70) =    0.10 ; given_sigma_ion_O_4Pst_cm2(  70) =    0.04 ; given_sigma_ion_O_2Pst_cm2(  70) =    0.03
  given_sigma_ion_O_4S_cm2(  71) =  0.14 ; given_sigma_ion_O_2D_cm2(  71) =    0.15 ; given_sigma_ion_O_2P_cm2(  71) =    0.09 ; given_sigma_ion_O_4Pst_cm2(  71) =    0.05 ; given_sigma_ion_O_2Pst_cm2(  71) =    0.03
  given_sigma_ion_O_4S_cm2(  72) =  0.14 ; given_sigma_ion_O_2D_cm2(  72) =    0.16 ; given_sigma_ion_O_2P_cm2(  72) =    0.10 ; given_sigma_ion_O_4Pst_cm2(  72) =    0.05 ; given_sigma_ion_O_2Pst_cm2(  72) =    0.03
  given_sigma_ion_O_4S_cm2(  73) =  0.14 ; given_sigma_ion_O_2D_cm2(  73) =    0.15 ; given_sigma_ion_O_2P_cm2(  73) =    0.10 ; given_sigma_ion_O_4Pst_cm2(  73) =    0.05 ; given_sigma_ion_O_2Pst_cm2(  73) =    0.03
  given_sigma_ion_O_4S_cm2(  74) =  0.15 ; given_sigma_ion_O_2D_cm2(  74) =    0.16 ; given_sigma_ion_O_2P_cm2(  74) =    0.10 ; given_sigma_ion_O_4Pst_cm2(  74) =    0.05 ; given_sigma_ion_O_2Pst_cm2(  74) =    0.03
  given_sigma_ion_O_4S_cm2(  75) =  0.15 ; given_sigma_ion_O_2D_cm2(  75) =    0.16 ; given_sigma_ion_O_2P_cm2(  75) =    0.11 ; given_sigma_ion_O_4Pst_cm2(  75) =    0.05 ; given_sigma_ion_O_2Pst_cm2(  75) =    0.04
  given_sigma_ion_O_4S_cm2(  76) =  0.15 ; given_sigma_ion_O_2D_cm2(  76) =    0.16 ; given_sigma_ion_O_2P_cm2(  76) =    0.11 ; given_sigma_ion_O_4Pst_cm2(  76) =    0.05 ; given_sigma_ion_O_2Pst_cm2(  76) =    0.04
  given_sigma_ion_O_4S_cm2(  77) =  0.15 ; given_sigma_ion_O_2D_cm2(  77) =    0.16 ; given_sigma_ion_O_2P_cm2(  77) =    0.11 ; given_sigma_ion_O_4Pst_cm2(  77) =    0.05 ; given_sigma_ion_O_2Pst_cm2(  77) =    0.04
  given_sigma_ion_O_4S_cm2(  78) =  0.15 ; given_sigma_ion_O_2D_cm2(  78) =    0.16 ; given_sigma_ion_O_2P_cm2(  78) =    0.11 ; given_sigma_ion_O_4Pst_cm2(  78) =    0.05 ; given_sigma_ion_O_2Pst_cm2(  78) =    0.04
  given_sigma_ion_O_4S_cm2(  79) =  0.16 ; given_sigma_ion_O_2D_cm2(  79) =    0.16 ; given_sigma_ion_O_2P_cm2(  79) =    0.11 ; given_sigma_ion_O_4Pst_cm2(  79) =    0.05 ; given_sigma_ion_O_2Pst_cm2(  79) =    0.04
  given_sigma_ion_O_4S_cm2(  80) =  0.15 ; given_sigma_ion_O_2D_cm2(  80) =    0.17 ; given_sigma_ion_O_2P_cm2(  80) =    0.11 ; given_sigma_ion_O_4Pst_cm2(  80) =    0.05 ; given_sigma_ion_O_2Pst_cm2(  80) =    0.04
  given_sigma_ion_O_4S_cm2(  81) =  0.16 ; given_sigma_ion_O_2D_cm2(  81) =    0.17 ; given_sigma_ion_O_2P_cm2(  81) =    0.11 ; given_sigma_ion_O_4Pst_cm2(  81) =    0.05 ; given_sigma_ion_O_2Pst_cm2(  81) =    0.04
  given_sigma_ion_O_4S_cm2(  82) =  0.17 ; given_sigma_ion_O_2D_cm2(  82) =    0.18 ; given_sigma_ion_O_2P_cm2(  82) =    0.12 ; given_sigma_ion_O_4Pst_cm2(  82) =    0.06 ; given_sigma_ion_O_2Pst_cm2(  82) =    0.04
  given_sigma_ion_O_4S_cm2(  83) =  0.16 ; given_sigma_ion_O_2D_cm2(  83) =    0.18 ; given_sigma_ion_O_2P_cm2(  83) =    0.12 ; given_sigma_ion_O_4Pst_cm2(  83) =    0.06 ; given_sigma_ion_O_2Pst_cm2(  83) =    0.04
  given_sigma_ion_O_4S_cm2(  84) =  0.17 ; given_sigma_ion_O_2D_cm2(  84) =    0.18 ; given_sigma_ion_O_2P_cm2(  84) =    0.12 ; given_sigma_ion_O_4Pst_cm2(  84) =    0.06 ; given_sigma_ion_O_2Pst_cm2(  84) =    0.04
  given_sigma_ion_O_4S_cm2(  85) =  0.17 ; given_sigma_ion_O_2D_cm2(  85) =    0.18 ; given_sigma_ion_O_2P_cm2(  85) =    0.12 ; given_sigma_ion_O_4Pst_cm2(  85) =    0.06 ; given_sigma_ion_O_2Pst_cm2(  85) =    0.04
  given_sigma_ion_O_4S_cm2(  86) =  0.17 ; given_sigma_ion_O_2D_cm2(  86) =    0.18 ; given_sigma_ion_O_2P_cm2(  86) =    0.12 ; given_sigma_ion_O_4Pst_cm2(  86) =    0.06 ; given_sigma_ion_O_2Pst_cm2(  86) =    0.04
  given_sigma_ion_O_4S_cm2(  87) =  0.17 ; given_sigma_ion_O_2D_cm2(  87) =    0.18 ; given_sigma_ion_O_2P_cm2(  87) =    0.12 ; given_sigma_ion_O_4Pst_cm2(  87) =    0.05 ; given_sigma_ion_O_2Pst_cm2(  87) =    0.05
  given_sigma_ion_O_4S_cm2(  88) =  0.17 ; given_sigma_ion_O_2D_cm2(  88) =    0.19 ; given_sigma_ion_O_2P_cm2(  88) =    0.12 ; given_sigma_ion_O_4Pst_cm2(  88) =    0.05 ; given_sigma_ion_O_2Pst_cm2(  88) =    0.05
  given_sigma_ion_O_4S_cm2(  89) =  0.17 ; given_sigma_ion_O_2D_cm2(  89) =    0.18 ; given_sigma_ion_O_2P_cm2(  89) =    0.12 ; given_sigma_ion_O_4Pst_cm2(  89) =    0.06 ; given_sigma_ion_O_2Pst_cm2(  89) =    0.05
  given_sigma_ion_O_4S_cm2(  90) =  0.18 ; given_sigma_ion_O_2D_cm2(  90) =    0.19 ; given_sigma_ion_O_2P_cm2(  90) =    0.13 ; given_sigma_ion_O_4Pst_cm2(  90) =    0.06 ; given_sigma_ion_O_2Pst_cm2(  90) =    0.05
  given_sigma_ion_O_4S_cm2(  91) =  0.18 ; given_sigma_ion_O_2D_cm2(  91) =    0.19 ; given_sigma_ion_O_2P_cm2(  91) =    0.13 ; given_sigma_ion_O_4Pst_cm2(  91) =    0.06 ; given_sigma_ion_O_2Pst_cm2(  91) =    0.05
  given_sigma_ion_O_4S_cm2(  92) =  0.18 ; given_sigma_ion_O_2D_cm2(  92) =    0.19 ; given_sigma_ion_O_2P_cm2(  92) =    0.12 ; given_sigma_ion_O_4Pst_cm2(  92) =    0.06 ; given_sigma_ion_O_2Pst_cm2(  92) =    0.05
  given_sigma_ion_O_4S_cm2(  93) =  0.18 ; given_sigma_ion_O_2D_cm2(  93) =    0.20 ; given_sigma_ion_O_2P_cm2(  93) =    0.12 ; given_sigma_ion_O_4Pst_cm2(  93) =    0.06 ; given_sigma_ion_O_2Pst_cm2(  93) =    0.05
  given_sigma_ion_O_4S_cm2(  94) =  0.18 ; given_sigma_ion_O_2D_cm2(  94) =    0.20 ; given_sigma_ion_O_2P_cm2(  94) =    0.13 ; given_sigma_ion_O_4Pst_cm2(  94) =    0.06 ; given_sigma_ion_O_2Pst_cm2(  94) =    0.04
  given_sigma_ion_O_4S_cm2(  95) =  0.18 ; given_sigma_ion_O_2D_cm2(  95) =    0.20 ; given_sigma_ion_O_2P_cm2(  95) =    0.13 ; given_sigma_ion_O_4Pst_cm2(  95) =    0.06 ; given_sigma_ion_O_2Pst_cm2(  95) =    0.05
  given_sigma_ion_O_4S_cm2(  96) =  0.19 ; given_sigma_ion_O_2D_cm2(  96) =    0.20 ; given_sigma_ion_O_2P_cm2(  96) =    0.13 ; given_sigma_ion_O_4Pst_cm2(  96) =    0.06 ; given_sigma_ion_O_2Pst_cm2(  96) =    0.05
  given_sigma_ion_O_4S_cm2(  97) =  0.19 ; given_sigma_ion_O_2D_cm2(  97) =    0.21 ; given_sigma_ion_O_2P_cm2(  97) =    0.14 ; given_sigma_ion_O_4Pst_cm2(  97) =    0.06 ; given_sigma_ion_O_2Pst_cm2(  97) =    0.05
  given_sigma_ion_O_4S_cm2(  98) =  0.20 ; given_sigma_ion_O_2D_cm2(  98) =    0.21 ; given_sigma_ion_O_2P_cm2(  98) =    0.14 ; given_sigma_ion_O_4Pst_cm2(  98) =    0.07 ; given_sigma_ion_O_2Pst_cm2(  98) =    0.05
  given_sigma_ion_O_4S_cm2(  99) =  0.20 ; given_sigma_ion_O_2D_cm2(  99) =    0.21 ; given_sigma_ion_O_2P_cm2(  99) =    0.14 ; given_sigma_ion_O_4Pst_cm2(  99) =    0.07 ; given_sigma_ion_O_2Pst_cm2(  99) =    0.05
  given_sigma_ion_O_4S_cm2( 100) =  0.20 ; given_sigma_ion_O_2D_cm2( 100) =    0.21 ; given_sigma_ion_O_2P_cm2( 100) =    0.14 ; given_sigma_ion_O_4Pst_cm2( 100) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 100) =    0.05
  given_sigma_ion_O_4S_cm2( 101) =  0.20 ; given_sigma_ion_O_2D_cm2( 101) =    0.22 ; given_sigma_ion_O_2P_cm2( 101) =    0.14 ; given_sigma_ion_O_4Pst_cm2( 101) =    0.06 ; given_sigma_ion_O_2Pst_cm2( 101) =    0.05
  given_sigma_ion_O_4S_cm2( 102) =  0.20 ; given_sigma_ion_O_2D_cm2( 102) =    0.21 ; given_sigma_ion_O_2P_cm2( 102) =    0.14 ; given_sigma_ion_O_4Pst_cm2( 102) =    0.06 ; given_sigma_ion_O_2Pst_cm2( 102) =    0.05
  given_sigma_ion_O_4S_cm2( 103) =  0.21 ; given_sigma_ion_O_2D_cm2( 103) =    0.22 ; given_sigma_ion_O_2P_cm2( 103) =    0.14 ; given_sigma_ion_O_4Pst_cm2( 103) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 103) =    0.05
  given_sigma_ion_O_4S_cm2( 104) =  0.21 ; given_sigma_ion_O_2D_cm2( 104) =    0.22 ; given_sigma_ion_O_2P_cm2( 104) =    0.15 ; given_sigma_ion_O_4Pst_cm2( 104) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 104) =    0.06
  given_sigma_ion_O_4S_cm2( 105) =  0.21 ; given_sigma_ion_O_2D_cm2( 105) =    0.23 ; given_sigma_ion_O_2P_cm2( 105) =    0.15 ; given_sigma_ion_O_4Pst_cm2( 105) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 105) =    0.06
  given_sigma_ion_O_4S_cm2( 106) =  0.21 ; given_sigma_ion_O_2D_cm2( 106) =    0.23 ; given_sigma_ion_O_2P_cm2( 106) =    0.15 ; given_sigma_ion_O_4Pst_cm2( 106) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 106) =    0.06
  given_sigma_ion_O_4S_cm2( 107) =  0.21 ; given_sigma_ion_O_2D_cm2( 107) =    0.23 ; given_sigma_ion_O_2P_cm2( 107) =    0.15 ; given_sigma_ion_O_4Pst_cm2( 107) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 107) =    0.06
  given_sigma_ion_O_4S_cm2( 108) =  0.22 ; given_sigma_ion_O_2D_cm2( 108) =    0.23 ; given_sigma_ion_O_2P_cm2( 108) =    0.15 ; given_sigma_ion_O_4Pst_cm2( 108) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 108) =    0.06
  given_sigma_ion_O_4S_cm2( 109) =  0.22 ; given_sigma_ion_O_2D_cm2( 109) =    0.23 ; given_sigma_ion_O_2P_cm2( 109) =    0.15 ; given_sigma_ion_O_4Pst_cm2( 109) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 109) =    0.06
  given_sigma_ion_O_4S_cm2( 110) =  0.22 ; given_sigma_ion_O_2D_cm2( 110) =    0.23 ; given_sigma_ion_O_2P_cm2( 110) =    0.16 ; given_sigma_ion_O_4Pst_cm2( 110) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 110) =    0.05
  given_sigma_ion_O_4S_cm2( 111) =  0.22 ; given_sigma_ion_O_2D_cm2( 111) =    0.24 ; given_sigma_ion_O_2P_cm2( 111) =    0.16 ; given_sigma_ion_O_4Pst_cm2( 111) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 111) =    0.05
  given_sigma_ion_O_4S_cm2( 112) =  0.22 ; given_sigma_ion_O_2D_cm2( 112) =    0.24 ; given_sigma_ion_O_2P_cm2( 112) =    0.16 ; given_sigma_ion_O_4Pst_cm2( 112) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 112) =    0.05
  given_sigma_ion_O_4S_cm2( 113) =  0.23 ; given_sigma_ion_O_2D_cm2( 113) =    0.24 ; given_sigma_ion_O_2P_cm2( 113) =    0.16 ; given_sigma_ion_O_4Pst_cm2( 113) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 113) =    0.05
  given_sigma_ion_O_4S_cm2( 114) =  0.23 ; given_sigma_ion_O_2D_cm2( 114) =    0.24 ; given_sigma_ion_O_2P_cm2( 114) =    0.16 ; given_sigma_ion_O_4Pst_cm2( 114) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 114) =    0.06
  given_sigma_ion_O_4S_cm2( 115) =  0.23 ; given_sigma_ion_O_2D_cm2( 115) =    0.24 ; given_sigma_ion_O_2P_cm2( 115) =    0.16 ; given_sigma_ion_O_4Pst_cm2( 115) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 115) =    0.06
  given_sigma_ion_O_4S_cm2( 116) =  0.23 ; given_sigma_ion_O_2D_cm2( 116) =    0.25 ; given_sigma_ion_O_2P_cm2( 116) =    0.16 ; given_sigma_ion_O_4Pst_cm2( 116) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 116) =    0.06
  given_sigma_ion_O_4S_cm2( 117) =  0.23 ; given_sigma_ion_O_2D_cm2( 117) =    0.25 ; given_sigma_ion_O_2P_cm2( 117) =    0.16 ; given_sigma_ion_O_4Pst_cm2( 117) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 117) =    0.06
  given_sigma_ion_O_4S_cm2( 118) =  0.23 ; given_sigma_ion_O_2D_cm2( 118) =    0.25 ; given_sigma_ion_O_2P_cm2( 118) =    0.16 ; given_sigma_ion_O_4Pst_cm2( 118) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 118) =    0.06
  given_sigma_ion_O_4S_cm2( 119) =  0.23 ; given_sigma_ion_O_2D_cm2( 119) =    0.25 ; given_sigma_ion_O_2P_cm2( 119) =    0.16 ; given_sigma_ion_O_4Pst_cm2( 119) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 119) =    0.06
  given_sigma_ion_O_4S_cm2( 120) =  0.24 ; given_sigma_ion_O_2D_cm2( 120) =    0.25 ; given_sigma_ion_O_2P_cm2( 120) =    0.17 ; given_sigma_ion_O_4Pst_cm2( 120) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 120) =    0.06
  given_sigma_ion_O_4S_cm2( 121) =  0.24 ; given_sigma_ion_O_2D_cm2( 121) =    0.25 ; given_sigma_ion_O_2P_cm2( 121) =    0.17 ; given_sigma_ion_O_4Pst_cm2( 121) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 121) =    0.06
  given_sigma_ion_O_4S_cm2( 122) =  0.24 ; given_sigma_ion_O_2D_cm2( 122) =    0.25 ; given_sigma_ion_O_2P_cm2( 122) =    0.17 ; given_sigma_ion_O_4Pst_cm2( 122) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 122) =    0.06
  given_sigma_ion_O_4S_cm2( 123) =  0.24 ; given_sigma_ion_O_2D_cm2( 123) =    0.26 ; given_sigma_ion_O_2P_cm2( 123) =    0.17 ; given_sigma_ion_O_4Pst_cm2( 123) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 123) =    0.06
  given_sigma_ion_O_4S_cm2( 124) =  0.24 ; given_sigma_ion_O_2D_cm2( 124) =    0.26 ; given_sigma_ion_O_2P_cm2( 124) =    0.17 ; given_sigma_ion_O_4Pst_cm2( 124) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 124) =    0.06
  given_sigma_ion_O_4S_cm2( 125) =  0.24 ; given_sigma_ion_O_2D_cm2( 125) =    0.26 ; given_sigma_ion_O_2P_cm2( 125) =    0.17 ; given_sigma_ion_O_4Pst_cm2( 125) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 125) =    0.06
  given_sigma_ion_O_4S_cm2( 126) =  0.24 ; given_sigma_ion_O_2D_cm2( 126) =    0.26 ; given_sigma_ion_O_2P_cm2( 126) =    0.17 ; given_sigma_ion_O_4Pst_cm2( 126) =    0.07 ; given_sigma_ion_O_2Pst_cm2( 126) =    0.07
  given_sigma_ion_O_4S_cm2( 127) =  0.25 ; given_sigma_ion_O_2D_cm2( 127) =    0.26 ; given_sigma_ion_O_2P_cm2( 127) =    0.17 ; given_sigma_ion_O_4Pst_cm2( 127) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 127) =    0.07
  given_sigma_ion_O_4S_cm2( 128) =  0.25 ; given_sigma_ion_O_2D_cm2( 128) =    0.27 ; given_sigma_ion_O_2P_cm2( 128) =    0.17 ; given_sigma_ion_O_4Pst_cm2( 128) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 128) =    0.07
  given_sigma_ion_O_4S_cm2( 129) =  0.25 ; given_sigma_ion_O_2D_cm2( 129) =    0.27 ; given_sigma_ion_O_2P_cm2( 129) =    0.18 ; given_sigma_ion_O_4Pst_cm2( 129) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 129) =    0.07
  given_sigma_ion_O_4S_cm2( 130) =  0.25 ; given_sigma_ion_O_2D_cm2( 130) =    0.27 ; given_sigma_ion_O_2P_cm2( 130) =    0.18 ; given_sigma_ion_O_4Pst_cm2( 130) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 130) =    0.07
  given_sigma_ion_O_4S_cm2( 131) =  0.25 ; given_sigma_ion_O_2D_cm2( 131) =    0.27 ; given_sigma_ion_O_2P_cm2( 131) =    0.18 ; given_sigma_ion_O_4Pst_cm2( 131) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 131) =    0.06
  given_sigma_ion_O_4S_cm2( 132) =  0.26 ; given_sigma_ion_O_2D_cm2( 132) =    0.27 ; given_sigma_ion_O_2P_cm2( 132) =    0.18 ; given_sigma_ion_O_4Pst_cm2( 132) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 132) =    0.06
  given_sigma_ion_O_4S_cm2( 133) =  0.26 ; given_sigma_ion_O_2D_cm2( 133) =    0.28 ; given_sigma_ion_O_2P_cm2( 133) =    0.18 ; given_sigma_ion_O_4Pst_cm2( 133) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 133) =    0.07
  given_sigma_ion_O_4S_cm2( 134) =  0.26 ; given_sigma_ion_O_2D_cm2( 134) =    0.28 ; given_sigma_ion_O_2P_cm2( 134) =    0.18 ; given_sigma_ion_O_4Pst_cm2( 134) =    0.09 ; given_sigma_ion_O_2Pst_cm2( 134) =    0.07
  given_sigma_ion_O_4S_cm2( 135) =  0.27 ; given_sigma_ion_O_2D_cm2( 135) =    0.28 ; given_sigma_ion_O_2P_cm2( 135) =    0.19 ; given_sigma_ion_O_4Pst_cm2( 135) =    0.09 ; given_sigma_ion_O_2Pst_cm2( 135) =    0.07
  given_sigma_ion_O_4S_cm2( 136) =  0.26 ; given_sigma_ion_O_2D_cm2( 136) =    0.28 ; given_sigma_ion_O_2P_cm2( 136) =    0.19 ; given_sigma_ion_O_4Pst_cm2( 136) =    0.09 ; given_sigma_ion_O_2Pst_cm2( 136) =    0.07
  given_sigma_ion_O_4S_cm2( 137) =  0.27 ; given_sigma_ion_O_2D_cm2( 137) =    0.29 ; given_sigma_ion_O_2P_cm2( 137) =    0.19 ; given_sigma_ion_O_4Pst_cm2( 137) =    0.09 ; given_sigma_ion_O_2Pst_cm2( 137) =    0.07
  given_sigma_ion_O_4S_cm2( 138) =  0.27 ; given_sigma_ion_O_2D_cm2( 138) =    0.29 ; given_sigma_ion_O_2P_cm2( 138) =    0.19 ; given_sigma_ion_O_4Pst_cm2( 138) =    0.08 ; given_sigma_ion_O_2Pst_cm2( 138) =    0.06
  given_sigma_ion_O_4S_cm2( 139) =  0.27 ; given_sigma_ion_O_2D_cm2( 139) =    0.29 ; given_sigma_ion_O_2P_cm2( 139) =    0.19 ; given_sigma_ion_O_4Pst_cm2( 139) =    0.09 ; given_sigma_ion_O_2Pst_cm2( 139) =    0.06
  given_sigma_ion_O_4S_cm2( 140) =  0.27 ; given_sigma_ion_O_2D_cm2( 140) =    0.29 ; given_sigma_ion_O_2P_cm2( 140) =    0.19 ; given_sigma_ion_O_4Pst_cm2( 140) =    0.09 ; given_sigma_ion_O_2Pst_cm2( 140) =    0.06
  given_sigma_ion_O_4S_cm2( 141) =  0.28 ; given_sigma_ion_O_2D_cm2( 141) =    0.30 ; given_sigma_ion_O_2P_cm2( 141) =    0.19 ; given_sigma_ion_O_4Pst_cm2( 141) =    0.09 ; given_sigma_ion_O_2Pst_cm2( 141) =    0.07
  given_sigma_ion_O_4S_cm2( 142) =  0.28 ; given_sigma_ion_O_2D_cm2( 142) =    0.30 ; given_sigma_ion_O_2P_cm2( 142) =    0.20 ; given_sigma_ion_O_4Pst_cm2( 142) =    0.09 ; given_sigma_ion_O_2Pst_cm2( 142) =    0.07
  given_sigma_ion_O_4S_cm2( 143) =  0.28 ; given_sigma_ion_O_2D_cm2( 143) =    0.30 ; given_sigma_ion_O_2P_cm2( 143) =    0.20 ; given_sigma_ion_O_4Pst_cm2( 143) =    0.09 ; given_sigma_ion_O_2Pst_cm2( 143) =    0.08
  given_sigma_ion_O_4S_cm2( 144) =  0.29 ; given_sigma_ion_O_2D_cm2( 144) =    0.31 ; given_sigma_ion_O_2P_cm2( 144) =    0.20 ; given_sigma_ion_O_4Pst_cm2( 144) =    0.09 ; given_sigma_ion_O_2Pst_cm2( 144) =    0.08
  given_sigma_ion_O_4S_cm2( 145) =  0.29 ; given_sigma_ion_O_2D_cm2( 145) =    0.31 ; given_sigma_ion_O_2P_cm2( 145) =    0.20 ; given_sigma_ion_O_4Pst_cm2( 145) =    0.10 ; given_sigma_ion_O_2Pst_cm2( 145) =    0.08
  given_sigma_ion_O_4S_cm2( 146) =  0.29 ; given_sigma_ion_O_2D_cm2( 146) =    0.31 ; given_sigma_ion_O_2P_cm2( 146) =    0.20 ; given_sigma_ion_O_4Pst_cm2( 146) =    0.10 ; given_sigma_ion_O_2Pst_cm2( 146) =    0.07
  given_sigma_ion_O_4S_cm2( 147) =  0.29 ; given_sigma_ion_O_2D_cm2( 147) =    0.31 ; given_sigma_ion_O_2P_cm2( 147) =    0.21 ; given_sigma_ion_O_4Pst_cm2( 147) =    0.10 ; given_sigma_ion_O_2Pst_cm2( 147) =    0.07
  given_sigma_ion_O_4S_cm2( 148) =  0.30 ; given_sigma_ion_O_2D_cm2( 148) =    0.32 ; given_sigma_ion_O_2P_cm2( 148) =    0.21 ; given_sigma_ion_O_4Pst_cm2( 148) =    0.10 ; given_sigma_ion_O_2Pst_cm2( 148) =    0.08
  given_sigma_ion_O_4S_cm2( 149) =  0.30 ; given_sigma_ion_O_2D_cm2( 149) =    0.32 ; given_sigma_ion_O_2P_cm2( 149) =    0.21 ; given_sigma_ion_O_4Pst_cm2( 149) =    0.09 ; given_sigma_ion_O_2Pst_cm2( 149) =    0.08
  given_sigma_ion_O_4S_cm2( 150) =  0.30 ; given_sigma_ion_O_2D_cm2( 150) =    0.31 ; given_sigma_ion_O_2P_cm2( 150) =    0.21 ; given_sigma_ion_O_4Pst_cm2( 150) =    0.09 ; given_sigma_ion_O_2Pst_cm2( 150) =    0.08
  given_sigma_ion_O_4S_cm2( 151) =  0.30 ; given_sigma_ion_O_2D_cm2( 151) =    0.32 ; given_sigma_ion_O_2P_cm2( 151) =    0.21 ; given_sigma_ion_O_4Pst_cm2( 151) =    0.09 ; given_sigma_ion_O_2Pst_cm2( 151) =    0.08
  given_sigma_ion_O_4S_cm2( 152) =  0.30 ; given_sigma_ion_O_2D_cm2( 152) =    0.32 ; given_sigma_ion_O_2P_cm2( 152) =    0.21 ; given_sigma_ion_O_4Pst_cm2( 152) =    0.10 ; given_sigma_ion_O_2Pst_cm2( 152) =    0.08
  given_sigma_ion_O_4S_cm2( 153) =  0.31 ; given_sigma_ion_O_2D_cm2( 153) =    0.33 ; given_sigma_ion_O_2P_cm2( 153) =    0.21 ; given_sigma_ion_O_4Pst_cm2( 153) =    0.10 ; given_sigma_ion_O_2Pst_cm2( 153) =    0.08
  given_sigma_ion_O_4S_cm2( 154) =  0.31 ; given_sigma_ion_O_2D_cm2( 154) =    0.33 ; given_sigma_ion_O_2P_cm2( 154) =    0.22 ; given_sigma_ion_O_4Pst_cm2( 154) =    0.10 ; given_sigma_ion_O_2Pst_cm2( 154) =    0.07
  given_sigma_ion_O_4S_cm2( 155) =  0.31 ; given_sigma_ion_O_2D_cm2( 155) =    0.33 ; given_sigma_ion_O_2P_cm2( 155) =    0.22 ; given_sigma_ion_O_4Pst_cm2( 155) =    0.10 ; given_sigma_ion_O_2Pst_cm2( 155) =    0.07
  given_sigma_ion_O_4S_cm2( 156) =  0.31 ; given_sigma_ion_O_2D_cm2( 156) =    0.33 ; given_sigma_ion_O_2P_cm2( 156) =    0.22 ; given_sigma_ion_O_4Pst_cm2( 156) =    0.09 ; given_sigma_ion_O_2Pst_cm2( 156) =    0.08
  given_sigma_ion_O_4S_cm2( 157) =  0.32 ; given_sigma_ion_O_2D_cm2( 157) =    0.34 ; given_sigma_ion_O_2P_cm2( 157) =    0.22 ; given_sigma_ion_O_4Pst_cm2( 157) =    0.11 ; given_sigma_ion_O_2Pst_cm2( 157) =    0.08
  given_sigma_ion_O_4S_cm2( 158) =  0.32 ; given_sigma_ion_O_2D_cm2( 158) =    0.34 ; given_sigma_ion_O_2P_cm2( 158) =    0.22 ; given_sigma_ion_O_4Pst_cm2( 158) =    0.11 ; given_sigma_ion_O_2Pst_cm2( 158) =    0.08
  given_sigma_ion_O_4S_cm2( 159) =  0.32 ; given_sigma_ion_O_2D_cm2( 159) =    0.34 ; given_sigma_ion_O_2P_cm2( 159) =    0.22 ; given_sigma_ion_O_4Pst_cm2( 159) =    0.11 ; given_sigma_ion_O_2Pst_cm2( 159) =    0.08
  given_sigma_ion_O_4S_cm2( 160) =  0.32 ; given_sigma_ion_O_2D_cm2( 160) =    0.34 ; given_sigma_ion_O_2P_cm2( 160) =    0.22 ; given_sigma_ion_O_4Pst_cm2( 160) =    0.11 ; given_sigma_ion_O_2Pst_cm2( 160) =    0.09
  given_sigma_ion_O_4S_cm2( 161) =  0.32 ; given_sigma_ion_O_2D_cm2( 161) =    0.34 ; given_sigma_ion_O_2P_cm2( 161) =    0.23 ; given_sigma_ion_O_4Pst_cm2( 161) =    0.11 ; given_sigma_ion_O_2Pst_cm2( 161) =    0.09
  given_sigma_ion_O_4S_cm2( 162) =  0.33 ; given_sigma_ion_O_2D_cm2( 162) =    0.35 ; given_sigma_ion_O_2P_cm2( 162) =    0.23 ; given_sigma_ion_O_4Pst_cm2( 162) =    0.10 ; given_sigma_ion_O_2Pst_cm2( 162) =    0.08
  given_sigma_ion_O_4S_cm2( 163) =  0.33 ; given_sigma_ion_O_2D_cm2( 163) =    0.35 ; given_sigma_ion_O_2P_cm2( 163) =    0.23 ; given_sigma_ion_O_4Pst_cm2( 163) =    0.10 ; given_sigma_ion_O_2Pst_cm2( 163) =    0.08
  given_sigma_ion_O_4S_cm2( 164) =  0.33 ; given_sigma_ion_O_2D_cm2( 164) =    0.35 ; given_sigma_ion_O_2P_cm2( 164) =    0.23 ; given_sigma_ion_O_4Pst_cm2( 164) =    0.11 ; given_sigma_ion_O_2Pst_cm2( 164) =    0.09
  given_sigma_ion_O_4S_cm2( 165) =  0.33 ; given_sigma_ion_O_2D_cm2( 165) =    0.36 ; given_sigma_ion_O_2P_cm2( 165) =    0.23 ; given_sigma_ion_O_4Pst_cm2( 165) =    0.11 ; given_sigma_ion_O_2Pst_cm2( 165) =    0.09
  given_sigma_ion_O_4S_cm2( 166) =  0.34 ; given_sigma_ion_O_2D_cm2( 166) =    0.36 ; given_sigma_ion_O_2P_cm2( 166) =    0.24 ; given_sigma_ion_O_4Pst_cm2( 166) =    0.11 ; given_sigma_ion_O_2Pst_cm2( 166) =    0.09
  given_sigma_ion_O_4S_cm2( 167) =  0.34 ; given_sigma_ion_O_2D_cm2( 167) =    0.36 ; given_sigma_ion_O_2P_cm2( 167) =    0.24 ; given_sigma_ion_O_4Pst_cm2( 167) =    0.10 ; given_sigma_ion_O_2Pst_cm2( 167) =    0.09
  given_sigma_ion_O_4S_cm2( 168) =  0.34 ; given_sigma_ion_O_2D_cm2( 168) =    0.37 ; given_sigma_ion_O_2P_cm2( 168) =    0.24 ; given_sigma_ion_O_4Pst_cm2( 168) =    0.11 ; given_sigma_ion_O_2Pst_cm2( 168) =    0.08
  given_sigma_ion_O_4S_cm2( 169) =  0.35 ; given_sigma_ion_O_2D_cm2( 169) =    0.37 ; given_sigma_ion_O_2P_cm2( 169) =    0.24 ; given_sigma_ion_O_4Pst_cm2( 169) =    0.12 ; given_sigma_ion_O_2Pst_cm2( 169) =    0.09
  given_sigma_ion_O_4S_cm2( 170) =  0.35 ; given_sigma_ion_O_2D_cm2( 170) =    0.37 ; given_sigma_ion_O_2P_cm2( 170) =    0.25 ; given_sigma_ion_O_4Pst_cm2( 170) =    0.12 ; given_sigma_ion_O_2Pst_cm2( 170) =    0.09
  given_sigma_ion_O_4S_cm2( 171) =  0.35 ; given_sigma_ion_O_2D_cm2( 171) =    0.37 ; given_sigma_ion_O_2P_cm2( 171) =    0.25 ; given_sigma_ion_O_4Pst_cm2( 171) =    0.11 ; given_sigma_ion_O_2Pst_cm2( 171) =    0.09
  given_sigma_ion_O_4S_cm2( 172) =  0.35 ; given_sigma_ion_O_2D_cm2( 172) =    0.38 ; given_sigma_ion_O_2P_cm2( 172) =    0.25 ; given_sigma_ion_O_4Pst_cm2( 172) =    0.11 ; given_sigma_ion_O_2Pst_cm2( 172) =    0.09
  given_sigma_ion_O_4S_cm2( 173) =  0.36 ; given_sigma_ion_O_2D_cm2( 173) =    0.38 ; given_sigma_ion_O_2P_cm2( 173) =    0.25 ; given_sigma_ion_O_4Pst_cm2( 173) =    0.12 ; given_sigma_ion_O_2Pst_cm2( 173) =    0.09
  given_sigma_ion_O_4S_cm2( 174) =  0.36 ; given_sigma_ion_O_2D_cm2( 174) =    0.38 ; given_sigma_ion_O_2P_cm2( 174) =    0.25 ; given_sigma_ion_O_4Pst_cm2( 174) =    0.12 ; given_sigma_ion_O_2Pst_cm2( 174) =    0.10
  given_sigma_ion_O_4S_cm2( 175) =  0.36 ; given_sigma_ion_O_2D_cm2( 175) =    0.39 ; given_sigma_ion_O_2P_cm2( 175) =    0.25 ; given_sigma_ion_O_4Pst_cm2( 175) =    0.12 ; given_sigma_ion_O_2Pst_cm2( 175) =    0.08
  given_sigma_ion_O_4S_cm2( 176) =  0.36 ; given_sigma_ion_O_2D_cm2( 176) =    0.39 ; given_sigma_ion_O_2P_cm2( 176) =    0.26 ; given_sigma_ion_O_4Pst_cm2( 176) =    0.12 ; given_sigma_ion_O_2Pst_cm2( 176) =    0.10
  given_sigma_ion_O_4S_cm2( 177) =  0.37 ; given_sigma_ion_O_2D_cm2( 177) =    0.40 ; given_sigma_ion_O_2P_cm2( 177) =    0.26 ; given_sigma_ion_O_4Pst_cm2( 177) =    0.12 ; given_sigma_ion_O_2Pst_cm2( 177) =    0.10
  given_sigma_ion_O_4S_cm2( 178) =  0.37 ; given_sigma_ion_O_2D_cm2( 178) =    0.40 ; given_sigma_ion_O_2P_cm2( 178) =    0.26 ; given_sigma_ion_O_4Pst_cm2( 178) =    0.12 ; given_sigma_ion_O_2Pst_cm2( 178) =    0.10
  given_sigma_ion_O_4S_cm2( 179) =  0.38 ; given_sigma_ion_O_2D_cm2( 179) =    0.40 ; given_sigma_ion_O_2P_cm2( 179) =    0.26 ; given_sigma_ion_O_4Pst_cm2( 179) =    0.13 ; given_sigma_ion_O_2Pst_cm2( 179) =    0.09
  given_sigma_ion_O_4S_cm2( 180) =  0.38 ; given_sigma_ion_O_2D_cm2( 180) =    0.41 ; given_sigma_ion_O_2P_cm2( 180) =    0.27 ; given_sigma_ion_O_4Pst_cm2( 180) =    0.13 ; given_sigma_ion_O_2Pst_cm2( 180) =    0.10
  given_sigma_ion_O_4S_cm2( 181) =  0.38 ; given_sigma_ion_O_2D_cm2( 181) =    0.41 ; given_sigma_ion_O_2P_cm2( 181) =    0.27 ; given_sigma_ion_O_4Pst_cm2( 181) =    0.13 ; given_sigma_ion_O_2Pst_cm2( 181) =    0.10
  given_sigma_ion_O_4S_cm2( 182) =  0.39 ; given_sigma_ion_O_2D_cm2( 182) =    0.41 ; given_sigma_ion_O_2P_cm2( 182) =    0.27 ; given_sigma_ion_O_4Pst_cm2( 182) =    0.13 ; given_sigma_ion_O_2Pst_cm2( 182) =    0.10
  given_sigma_ion_O_4S_cm2( 183) =  0.39 ; given_sigma_ion_O_2D_cm2( 183) =    0.42 ; given_sigma_ion_O_2P_cm2( 183) =    0.28 ; given_sigma_ion_O_4Pst_cm2( 183) =    0.13 ; given_sigma_ion_O_2Pst_cm2( 183) =    0.11
  given_sigma_ion_O_4S_cm2( 184) =  0.40 ; given_sigma_ion_O_2D_cm2( 184) =    0.42 ; given_sigma_ion_O_2P_cm2( 184) =    0.28 ; given_sigma_ion_O_4Pst_cm2( 184) =    0.13 ; given_sigma_ion_O_2Pst_cm2( 184) =    0.11
  given_sigma_ion_O_4S_cm2( 185) =  0.40 ; given_sigma_ion_O_2D_cm2( 185) =    0.43 ; given_sigma_ion_O_2P_cm2( 185) =    0.28 ; given_sigma_ion_O_4Pst_cm2( 185) =    0.12 ; given_sigma_ion_O_2Pst_cm2( 185) =    0.11
  given_sigma_ion_O_4S_cm2( 186) =  0.40 ; given_sigma_ion_O_2D_cm2( 186) =    0.43 ; given_sigma_ion_O_2P_cm2( 186) =    0.28 ; given_sigma_ion_O_4Pst_cm2( 186) =    0.12 ; given_sigma_ion_O_2Pst_cm2( 186) =    0.11
  given_sigma_ion_O_4S_cm2( 187) =  0.41 ; given_sigma_ion_O_2D_cm2( 187) =    0.43 ; given_sigma_ion_O_2P_cm2( 187) =    0.28 ; given_sigma_ion_O_4Pst_cm2( 187) =    0.14 ; given_sigma_ion_O_2Pst_cm2( 187) =    0.09
  given_sigma_ion_O_4S_cm2( 188) =  0.42 ; given_sigma_ion_O_2D_cm2( 188) =    0.44 ; given_sigma_ion_O_2P_cm2( 188) =    0.29 ; given_sigma_ion_O_4Pst_cm2( 188) =    0.13 ; given_sigma_ion_O_2Pst_cm2( 188) =    0.11
  given_sigma_ion_O_4S_cm2( 189) =  0.42 ; given_sigma_ion_O_2D_cm2( 189) =    0.45 ; given_sigma_ion_O_2P_cm2( 189) =    0.30 ; given_sigma_ion_O_4Pst_cm2( 189) =    0.13 ; given_sigma_ion_O_2Pst_cm2( 189) =    0.11
  given_sigma_ion_O_4S_cm2( 190) =  0.43 ; given_sigma_ion_O_2D_cm2( 190) =    0.46 ; given_sigma_ion_O_2P_cm2( 190) =    0.30 ; given_sigma_ion_O_4Pst_cm2( 190) =    0.14 ; given_sigma_ion_O_2Pst_cm2( 190) =    0.11
  given_sigma_ion_O_4S_cm2( 191) =  0.43 ; given_sigma_ion_O_2D_cm2( 191) =    0.46 ; given_sigma_ion_O_2P_cm2( 191) =    0.30 ; given_sigma_ion_O_4Pst_cm2( 191) =    0.14 ; given_sigma_ion_O_2Pst_cm2( 191) =    0.11
  given_sigma_ion_O_4S_cm2( 192) =  0.44 ; given_sigma_ion_O_2D_cm2( 192) =    0.46 ; given_sigma_ion_O_2P_cm2( 192) =    0.30 ; given_sigma_ion_O_4Pst_cm2( 192) =    0.13 ; given_sigma_ion_O_2Pst_cm2( 192) =    0.10
  given_sigma_ion_O_4S_cm2( 193) =  0.44 ; given_sigma_ion_O_2D_cm2( 193) =    0.47 ; given_sigma_ion_O_2P_cm2( 193) =    0.31 ; given_sigma_ion_O_4Pst_cm2( 193) =    0.15 ; given_sigma_ion_O_2Pst_cm2( 193) =    0.12
  given_sigma_ion_O_4S_cm2( 194) =  0.44 ; given_sigma_ion_O_2D_cm2( 194) =    0.47 ; given_sigma_ion_O_2P_cm2( 194) =    0.31 ; given_sigma_ion_O_4Pst_cm2( 194) =    0.15 ; given_sigma_ion_O_2Pst_cm2( 194) =    0.12
  given_sigma_ion_O_4S_cm2( 195) =  0.45 ; given_sigma_ion_O_2D_cm2( 195) =    0.48 ; given_sigma_ion_O_2P_cm2( 195) =    0.31 ; given_sigma_ion_O_4Pst_cm2( 195) =    0.13 ; given_sigma_ion_O_2Pst_cm2( 195) =    0.12
  given_sigma_ion_O_4S_cm2( 196) =  0.45 ; given_sigma_ion_O_2D_cm2( 196) =    0.48 ; given_sigma_ion_O_2P_cm2( 196) =    0.32 ; given_sigma_ion_O_4Pst_cm2( 196) =    0.14 ; given_sigma_ion_O_2Pst_cm2( 196) =    0.12
  given_sigma_ion_O_4S_cm2( 197) =  0.47 ; given_sigma_ion_O_2D_cm2( 197) =    0.50 ; given_sigma_ion_O_2P_cm2( 197) =    0.33 ; given_sigma_ion_O_4Pst_cm2( 197) =    0.15 ; given_sigma_ion_O_2Pst_cm2( 197) =    0.12
  given_sigma_ion_O_4S_cm2( 198) =  0.47 ; given_sigma_ion_O_2D_cm2( 198) =    0.50 ; given_sigma_ion_O_2P_cm2( 198) =    0.33 ; given_sigma_ion_O_4Pst_cm2( 198) =    0.16 ; given_sigma_ion_O_2Pst_cm2( 198) =    0.13
  given_sigma_ion_O_4S_cm2( 199) =  0.48 ; given_sigma_ion_O_2D_cm2( 199) =    0.51 ; given_sigma_ion_O_2P_cm2( 199) =    0.33 ; given_sigma_ion_O_4Pst_cm2( 199) =    0.16 ; given_sigma_ion_O_2Pst_cm2( 199) =    0.13
  given_sigma_ion_O_4S_cm2( 200) =  0.48 ; given_sigma_ion_O_2D_cm2( 200) =    0.51 ; given_sigma_ion_O_2P_cm2( 200) =    0.34 ; given_sigma_ion_O_4Pst_cm2( 200) =    0.16 ; given_sigma_ion_O_2Pst_cm2( 200) =    0.13
  given_sigma_ion_O_4S_cm2( 201) =  0.49 ; given_sigma_ion_O_2D_cm2( 201) =    0.52 ; given_sigma_ion_O_2P_cm2( 201) =    0.34 ; given_sigma_ion_O_4Pst_cm2( 201) =    0.16 ; given_sigma_ion_O_2Pst_cm2( 201) =    0.13
  given_sigma_ion_O_4S_cm2( 202) =  0.49 ; given_sigma_ion_O_2D_cm2( 202) =    0.52 ; given_sigma_ion_O_2P_cm2( 202) =    0.34 ; given_sigma_ion_O_4Pst_cm2( 202) =    0.16 ; given_sigma_ion_O_2Pst_cm2( 202) =    0.13
  given_sigma_ion_O_4S_cm2( 203) =  0.50 ; given_sigma_ion_O_2D_cm2( 203) =    0.53 ; given_sigma_ion_O_2P_cm2( 203) =    0.35 ; given_sigma_ion_O_4Pst_cm2( 203) =    0.17 ; given_sigma_ion_O_2Pst_cm2( 203) =    0.13
  given_sigma_ion_O_4S_cm2( 204) =  0.50 ; given_sigma_ion_O_2D_cm2( 204) =    0.54 ; given_sigma_ion_O_2P_cm2( 204) =    0.35 ; given_sigma_ion_O_4Pst_cm2( 204) =    0.17 ; given_sigma_ion_O_2Pst_cm2( 204) =    0.13
  given_sigma_ion_O_4S_cm2( 205) =  0.53 ; given_sigma_ion_O_2D_cm2( 205) =    0.56 ; given_sigma_ion_O_2P_cm2( 205) =    0.37 ; given_sigma_ion_O_4Pst_cm2( 205) =    0.18 ; given_sigma_ion_O_2Pst_cm2( 205) =    0.14
  given_sigma_ion_O_4S_cm2( 206) =  0.55 ; given_sigma_ion_O_2D_cm2( 206) =    0.58 ; given_sigma_ion_O_2P_cm2( 206) =    0.38 ; given_sigma_ion_O_4Pst_cm2( 206) =    0.18 ; given_sigma_ion_O_2Pst_cm2( 206) =    0.15
  given_sigma_ion_O_4S_cm2( 207) =  0.55 ; given_sigma_ion_O_2D_cm2( 207) =    0.59 ; given_sigma_ion_O_2P_cm2( 207) =    0.38 ; given_sigma_ion_O_4Pst_cm2( 207) =    0.18 ; given_sigma_ion_O_2Pst_cm2( 207) =    0.15
  given_sigma_ion_O_4S_cm2( 208) =  0.55 ; given_sigma_ion_O_2D_cm2( 208) =    0.59 ; given_sigma_ion_O_2P_cm2( 208) =    0.39 ; given_sigma_ion_O_4Pst_cm2( 208) =    0.18 ; given_sigma_ion_O_2Pst_cm2( 208) =    0.15
  given_sigma_ion_O_4S_cm2( 209) =  0.55 ; given_sigma_ion_O_2D_cm2( 209) =    0.59 ; given_sigma_ion_O_2P_cm2( 209) =    0.39 ; given_sigma_ion_O_4Pst_cm2( 209) =    0.18 ; given_sigma_ion_O_2Pst_cm2( 209) =    0.15
  given_sigma_ion_O_4S_cm2( 210) =  0.56 ; given_sigma_ion_O_2D_cm2( 210) =    0.59 ; given_sigma_ion_O_2P_cm2( 210) =    0.39 ; given_sigma_ion_O_4Pst_cm2( 210) =    0.19 ; given_sigma_ion_O_2Pst_cm2( 210) =    0.15
  given_sigma_ion_O_4S_cm2( 211) =  0.56 ; given_sigma_ion_O_2D_cm2( 211) =    0.60 ; given_sigma_ion_O_2P_cm2( 211) =    0.39 ; given_sigma_ion_O_4Pst_cm2( 211) =    0.19 ; given_sigma_ion_O_2Pst_cm2( 211) =    0.15
  given_sigma_ion_O_4S_cm2( 212) =  0.58 ; given_sigma_ion_O_2D_cm2( 212) =    0.63 ; given_sigma_ion_O_2P_cm2( 212) =    0.42 ; given_sigma_ion_O_4Pst_cm2( 212) =    0.20 ; given_sigma_ion_O_2Pst_cm2( 212) =    0.16
  given_sigma_ion_O_4S_cm2( 213) =  0.58 ; given_sigma_ion_O_2D_cm2( 213) =    0.64 ; given_sigma_ion_O_2P_cm2( 213) =    0.42 ; given_sigma_ion_O_4Pst_cm2( 213) =    0.20 ; given_sigma_ion_O_2Pst_cm2( 213) =    0.16
  given_sigma_ion_O_4S_cm2( 214) =  0.59 ; given_sigma_ion_O_2D_cm2( 214) =    0.65 ; given_sigma_ion_O_2P_cm2( 214) =    0.42 ; given_sigma_ion_O_4Pst_cm2( 214) =    0.20 ; given_sigma_ion_O_2Pst_cm2( 214) =    0.16
  given_sigma_ion_O_4S_cm2( 215) =  0.60 ; given_sigma_ion_O_2D_cm2( 215) =    0.66 ; given_sigma_ion_O_2P_cm2( 215) =    0.43 ; given_sigma_ion_O_4Pst_cm2( 215) =    0.21 ; given_sigma_ion_O_2Pst_cm2( 215) =    0.16
  given_sigma_ion_O_4S_cm2( 216) =  0.60 ; given_sigma_ion_O_2D_cm2( 216) =    0.66 ; given_sigma_ion_O_2P_cm2( 216) =    0.43 ; given_sigma_ion_O_4Pst_cm2( 216) =    0.21 ; given_sigma_ion_O_2Pst_cm2( 216) =    0.16
  given_sigma_ion_O_4S_cm2( 217) =  0.59 ; given_sigma_ion_O_2D_cm2( 217) =    0.66 ; given_sigma_ion_O_2P_cm2( 217) =    0.43 ; given_sigma_ion_O_4Pst_cm2( 217) =    0.20 ; given_sigma_ion_O_2Pst_cm2( 217) =    0.16
  given_sigma_ion_O_4S_cm2( 218) =  0.59 ; given_sigma_ion_O_2D_cm2( 218) =    0.65 ; given_sigma_ion_O_2P_cm2( 218) =    0.43 ; given_sigma_ion_O_4Pst_cm2( 218) =    0.20 ; given_sigma_ion_O_2Pst_cm2( 218) =    0.16
  given_sigma_ion_O_4S_cm2( 219) =  0.59 ; given_sigma_ion_O_2D_cm2( 219) =    0.67 ; given_sigma_ion_O_2P_cm2( 219) =    0.43 ; given_sigma_ion_O_4Pst_cm2( 219) =    0.20 ; given_sigma_ion_O_2Pst_cm2( 219) =    0.16
  given_sigma_ion_O_4S_cm2( 220) =  0.59 ; given_sigma_ion_O_2D_cm2( 220) =    0.67 ; given_sigma_ion_O_2P_cm2( 220) =    0.43 ; given_sigma_ion_O_4Pst_cm2( 220) =    0.20 ; given_sigma_ion_O_2Pst_cm2( 220) =    0.16
  given_sigma_ion_O_4S_cm2( 221) =  0.59 ; given_sigma_ion_O_2D_cm2( 221) =    0.65 ; given_sigma_ion_O_2P_cm2( 221) =    0.42 ; given_sigma_ion_O_4Pst_cm2( 221) =    0.20 ; given_sigma_ion_O_2Pst_cm2( 221) =    0.16
  given_sigma_ion_O_4S_cm2( 222) =  0.59 ; given_sigma_ion_O_2D_cm2( 222) =    0.66 ; given_sigma_ion_O_2P_cm2( 222) =    0.42 ; given_sigma_ion_O_4Pst_cm2( 222) =    0.20 ; given_sigma_ion_O_2Pst_cm2( 222) =    0.16
  given_sigma_ion_O_4S_cm2( 223) =  0.58 ; given_sigma_ion_O_2D_cm2( 223) =    0.67 ; given_sigma_ion_O_2P_cm2( 223) =    0.42 ; given_sigma_ion_O_4Pst_cm2( 223) =    0.20 ; given_sigma_ion_O_2Pst_cm2( 223) =    0.16
  given_sigma_ion_O_4S_cm2( 224) =  0.58 ; given_sigma_ion_O_2D_cm2( 224) =    0.66 ; given_sigma_ion_O_2P_cm2( 224) =    0.42 ; given_sigma_ion_O_4Pst_cm2( 224) =    0.20 ; given_sigma_ion_O_2Pst_cm2( 224) =    0.16
  given_sigma_ion_O_4S_cm2( 225) =  0.56 ; given_sigma_ion_O_2D_cm2( 225) =    0.66 ; given_sigma_ion_O_2P_cm2( 225) =    0.42 ; given_sigma_ion_O_4Pst_cm2( 225) =    0.20 ; given_sigma_ion_O_2Pst_cm2( 225) =    0.16
  given_sigma_ion_O_4S_cm2( 226) =  0.56 ; given_sigma_ion_O_2D_cm2( 226) =    0.66 ; given_sigma_ion_O_2P_cm2( 226) =    0.42 ; given_sigma_ion_O_4Pst_cm2( 226) =    0.20 ; given_sigma_ion_O_2Pst_cm2( 226) =    0.16
  given_sigma_ion_O_4S_cm2( 227) =  0.56 ; given_sigma_ion_O_2D_cm2( 227) =    0.66 ; given_sigma_ion_O_2P_cm2( 227) =    0.42 ; given_sigma_ion_O_4Pst_cm2( 227) =    0.20 ; given_sigma_ion_O_2Pst_cm2( 227) =    0.16
  given_sigma_ion_O_4S_cm2( 228) =  0.61 ; given_sigma_ion_O_2D_cm2( 228) =    0.72 ; given_sigma_ion_O_2P_cm2( 228) =    0.46 ; given_sigma_ion_O_4Pst_cm2( 228) =    0.22 ; given_sigma_ion_O_2Pst_cm2( 228) =    0.17
  given_sigma_ion_O_4S_cm2( 229) =  0.62 ; given_sigma_ion_O_2D_cm2( 229) =    0.73 ; given_sigma_ion_O_2P_cm2( 229) =    0.49 ; given_sigma_ion_O_4Pst_cm2( 229) =    0.22 ; given_sigma_ion_O_2Pst_cm2( 229) =    0.18
  given_sigma_ion_O_4S_cm2( 230) =  0.65 ; given_sigma_ion_O_2D_cm2( 230) =    0.77 ; given_sigma_ion_O_2P_cm2( 230) =    0.51 ; given_sigma_ion_O_4Pst_cm2( 230) =    0.23 ; given_sigma_ion_O_2Pst_cm2( 230) =    0.19
  given_sigma_ion_O_4S_cm2( 231) =  0.66 ; given_sigma_ion_O_2D_cm2( 231) =    0.77 ; given_sigma_ion_O_2P_cm2( 231) =    0.52 ; given_sigma_ion_O_4Pst_cm2( 231) =    0.23 ; given_sigma_ion_O_2Pst_cm2( 231) =    0.19
  given_sigma_ion_O_4S_cm2( 232) =  0.73 ; given_sigma_ion_O_2D_cm2( 232) =    0.85 ; given_sigma_ion_O_2P_cm2( 232) =    0.57 ; given_sigma_ion_O_4Pst_cm2( 232) =    0.26 ; given_sigma_ion_O_2Pst_cm2( 232) =    0.21
  given_sigma_ion_O_4S_cm2( 233) =  0.74 ; given_sigma_ion_O_2D_cm2( 233) =    0.87 ; given_sigma_ion_O_2P_cm2( 233) =    0.58 ; given_sigma_ion_O_4Pst_cm2( 233) =    0.26 ; given_sigma_ion_O_2Pst_cm2( 233) =    0.21
  given_sigma_ion_O_4S_cm2( 234) =  0.75 ; given_sigma_ion_O_2D_cm2( 234) =    0.88 ; given_sigma_ion_O_2P_cm2( 234) =    0.59 ; given_sigma_ion_O_4Pst_cm2( 234) =    0.27 ; given_sigma_ion_O_2Pst_cm2( 234) =    0.21
  given_sigma_ion_O_4S_cm2( 235) =  0.79 ; given_sigma_ion_O_2D_cm2( 235) =    0.93 ; given_sigma_ion_O_2P_cm2( 235) =    0.62 ; given_sigma_ion_O_4Pst_cm2( 235) =    0.28 ; given_sigma_ion_O_2Pst_cm2( 235) =    0.22
  given_sigma_ion_O_4S_cm2( 236) =  0.80 ; given_sigma_ion_O_2D_cm2( 236) =    0.94 ; given_sigma_ion_O_2P_cm2( 236) =    0.63 ; given_sigma_ion_O_4Pst_cm2( 236) =    0.29 ; given_sigma_ion_O_2Pst_cm2( 236) =    0.23
  given_sigma_ion_O_4S_cm2( 237) =  0.81 ; given_sigma_ion_O_2D_cm2( 237) =    0.95 ; given_sigma_ion_O_2P_cm2( 237) =    0.63 ; given_sigma_ion_O_4Pst_cm2( 237) =    0.29 ; given_sigma_ion_O_2Pst_cm2( 237) =    0.23
  given_sigma_ion_O_4S_cm2( 238) =  0.81 ; given_sigma_ion_O_2D_cm2( 238) =    0.96 ; given_sigma_ion_O_2P_cm2( 238) =    0.64 ; given_sigma_ion_O_4Pst_cm2( 238) =    0.29 ; given_sigma_ion_O_2Pst_cm2( 238) =    0.23
  given_sigma_ion_O_4S_cm2( 239) =  0.83 ; given_sigma_ion_O_2D_cm2( 239) =    0.97 ; given_sigma_ion_O_2P_cm2( 239) =    0.65 ; given_sigma_ion_O_4Pst_cm2( 239) =    0.30 ; given_sigma_ion_O_2Pst_cm2( 239) =    0.24
  given_sigma_ion_O_4S_cm2( 240) =  0.83 ; given_sigma_ion_O_2D_cm2( 240) =    0.98 ; given_sigma_ion_O_2P_cm2( 240) =    0.65 ; given_sigma_ion_O_4Pst_cm2( 240) =    0.30 ; given_sigma_ion_O_2Pst_cm2( 240) =    0.24
  given_sigma_ion_O_4S_cm2( 241) =  0.85 ; given_sigma_ion_O_2D_cm2( 241) =    1.00 ; given_sigma_ion_O_2P_cm2( 241) =    0.67 ; given_sigma_ion_O_4Pst_cm2( 241) =    0.30 ; given_sigma_ion_O_2Pst_cm2( 241) =    0.24
  given_sigma_ion_O_4S_cm2( 242) =  0.87 ; given_sigma_ion_O_2D_cm2( 242) =    1.03 ; given_sigma_ion_O_2P_cm2( 242) =    0.68 ; given_sigma_ion_O_4Pst_cm2( 242) =    0.31 ; given_sigma_ion_O_2Pst_cm2( 242) =    0.25
  given_sigma_ion_O_4S_cm2( 243) =  0.85 ; given_sigma_ion_O_2D_cm2( 243) =    1.04 ; given_sigma_ion_O_2P_cm2( 243) =    0.70 ; given_sigma_ion_O_4Pst_cm2( 243) =    0.32 ; given_sigma_ion_O_2Pst_cm2( 243) =    0.25
  given_sigma_ion_O_4S_cm2( 244) =  0.86 ; given_sigma_ion_O_2D_cm2( 244) =    1.05 ; given_sigma_ion_O_2P_cm2( 244) =    0.70 ; given_sigma_ion_O_4Pst_cm2( 244) =    0.32 ; given_sigma_ion_O_2Pst_cm2( 244) =    0.25
  given_sigma_ion_O_4S_cm2( 245) =  0.91 ; given_sigma_ion_O_2D_cm2( 245) =    1.08 ; given_sigma_ion_O_2P_cm2( 245) =    0.72 ; given_sigma_ion_O_4Pst_cm2( 245) =    0.33 ; given_sigma_ion_O_2Pst_cm2( 245) =    0.26
  given_sigma_ion_O_4S_cm2( 246) =  0.90 ; given_sigma_ion_O_2D_cm2( 246) =    1.09 ; given_sigma_ion_O_2P_cm2( 246) =    0.72 ; given_sigma_ion_O_4Pst_cm2( 246) =    0.33 ; given_sigma_ion_O_2Pst_cm2( 246) =    0.26
  given_sigma_ion_O_4S_cm2( 247) =  0.89 ; given_sigma_ion_O_2D_cm2( 247) =    1.09 ; given_sigma_ion_O_2P_cm2( 247) =    0.73 ; given_sigma_ion_O_4Pst_cm2( 247) =    0.33 ; given_sigma_ion_O_2Pst_cm2( 247) =    0.26
  given_sigma_ion_O_4S_cm2( 248) =  0.90 ; given_sigma_ion_O_2D_cm2( 248) =    1.10 ; given_sigma_ion_O_2P_cm2( 248) =    0.73 ; given_sigma_ion_O_4Pst_cm2( 248) =    0.33 ; given_sigma_ion_O_2Pst_cm2( 248) =    0.27
  given_sigma_ion_O_4S_cm2( 249) =  0.93 ; given_sigma_ion_O_2D_cm2( 249) =    1.13 ; given_sigma_ion_O_2P_cm2( 249) =    0.75 ; given_sigma_ion_O_4Pst_cm2( 249) =    0.34 ; given_sigma_ion_O_2Pst_cm2( 249) =    0.27
  given_sigma_ion_O_4S_cm2( 250) =  0.93 ; given_sigma_ion_O_2D_cm2( 250) =    1.13 ; given_sigma_ion_O_2P_cm2( 250) =    0.76 ; given_sigma_ion_O_4Pst_cm2( 250) =    0.34 ; given_sigma_ion_O_2Pst_cm2( 250) =    0.27
  given_sigma_ion_O_4S_cm2( 251) =  0.93 ; given_sigma_ion_O_2D_cm2( 251) =    1.14 ; given_sigma_ion_O_2P_cm2( 251) =    0.76 ; given_sigma_ion_O_4Pst_cm2( 251) =    0.35 ; given_sigma_ion_O_2Pst_cm2( 251) =    0.28
  given_sigma_ion_O_4S_cm2( 252) =  0.95 ; given_sigma_ion_O_2D_cm2( 252) =    1.16 ; given_sigma_ion_O_2P_cm2( 252) =    0.77 ; given_sigma_ion_O_4Pst_cm2( 252) =    0.35 ; given_sigma_ion_O_2Pst_cm2( 252) =    0.28
  given_sigma_ion_O_4S_cm2( 253) =  0.95 ; given_sigma_ion_O_2D_cm2( 253) =    1.16 ; given_sigma_ion_O_2P_cm2( 253) =    0.77 ; given_sigma_ion_O_4Pst_cm2( 253) =    0.35 ; given_sigma_ion_O_2Pst_cm2( 253) =    0.28
  given_sigma_ion_O_4S_cm2( 254) =  0.95 ; given_sigma_ion_O_2D_cm2( 254) =    1.16 ; given_sigma_ion_O_2P_cm2( 254) =    0.78 ; given_sigma_ion_O_4Pst_cm2( 254) =    0.35 ; given_sigma_ion_O_2Pst_cm2( 254) =    0.28
  given_sigma_ion_O_4S_cm2( 255) =  0.96 ; given_sigma_ion_O_2D_cm2( 255) =    1.17 ; given_sigma_ion_O_2P_cm2( 255) =    0.78 ; given_sigma_ion_O_4Pst_cm2( 255) =    0.35 ; given_sigma_ion_O_2Pst_cm2( 255) =    0.28
  given_sigma_ion_O_4S_cm2( 256) =  0.96 ; given_sigma_ion_O_2D_cm2( 256) =    1.21 ; given_sigma_ion_O_2P_cm2( 256) =    0.78 ; given_sigma_ion_O_4Pst_cm2( 256) =    0.36 ; given_sigma_ion_O_2Pst_cm2( 256) =    0.28
  given_sigma_ion_O_4S_cm2( 257) =  0.96 ; given_sigma_ion_O_2D_cm2( 257) =    1.18 ; given_sigma_ion_O_2P_cm2( 257) =    0.79 ; given_sigma_ion_O_4Pst_cm2( 257) =    0.36 ; given_sigma_ion_O_2Pst_cm2( 257) =    0.29
  given_sigma_ion_O_4S_cm2( 258) =  0.97 ; given_sigma_ion_O_2D_cm2( 258) =    1.23 ; given_sigma_ion_O_2P_cm2( 258) =    0.79 ; given_sigma_ion_O_4Pst_cm2( 258) =    0.36 ; given_sigma_ion_O_2Pst_cm2( 258) =    0.29
  given_sigma_ion_O_4S_cm2( 259) =  0.98 ; given_sigma_ion_O_2D_cm2( 259) =    1.23 ; given_sigma_ion_O_2P_cm2( 259) =    0.79 ; given_sigma_ion_O_4Pst_cm2( 259) =    0.36 ; given_sigma_ion_O_2Pst_cm2( 259) =    0.29
  given_sigma_ion_O_4S_cm2( 260) =  0.98 ; given_sigma_ion_O_2D_cm2( 260) =    1.24 ; given_sigma_ion_O_2P_cm2( 260) =    0.80 ; given_sigma_ion_O_4Pst_cm2( 260) =    0.36 ; given_sigma_ion_O_2Pst_cm2( 260) =    0.29
  given_sigma_ion_O_4S_cm2( 261) =  0.99 ; given_sigma_ion_O_2D_cm2( 261) =    1.25 ; given_sigma_ion_O_2P_cm2( 261) =    0.81 ; given_sigma_ion_O_4Pst_cm2( 261) =    0.37 ; given_sigma_ion_O_2Pst_cm2( 261) =    0.29
  given_sigma_ion_O_4S_cm2( 262) =  0.99 ; given_sigma_ion_O_2D_cm2( 262) =    1.21 ; given_sigma_ion_O_2P_cm2( 262) =    0.81 ; given_sigma_ion_O_4Pst_cm2( 262) =    0.37 ; given_sigma_ion_O_2Pst_cm2( 262) =    0.29
  given_sigma_ion_O_4S_cm2( 263) =  0.99 ; given_sigma_ion_O_2D_cm2( 263) =    1.25 ; given_sigma_ion_O_2P_cm2( 263) =    0.81 ; given_sigma_ion_O_4Pst_cm2( 263) =    0.37 ; given_sigma_ion_O_2Pst_cm2( 263) =    0.29
  given_sigma_ion_O_4S_cm2( 264) =  1.00 ; given_sigma_ion_O_2D_cm2( 264) =    1.26 ; given_sigma_ion_O_2P_cm2( 264) =    0.81 ; given_sigma_ion_O_4Pst_cm2( 264) =    0.37 ; given_sigma_ion_O_2Pst_cm2( 264) =    0.30
  given_sigma_ion_O_4S_cm2( 265) =  1.00 ; given_sigma_ion_O_2D_cm2( 265) =    1.26 ; given_sigma_ion_O_2P_cm2( 265) =    0.81 ; given_sigma_ion_O_4Pst_cm2( 265) =    0.37 ; given_sigma_ion_O_2Pst_cm2( 265) =    0.30
  given_sigma_ion_O_4S_cm2( 266) =  1.00 ; given_sigma_ion_O_2D_cm2( 266) =    1.26 ; given_sigma_ion_O_2P_cm2( 266) =    0.82 ; given_sigma_ion_O_4Pst_cm2( 266) =    0.37 ; given_sigma_ion_O_2Pst_cm2( 266) =    0.30
  given_sigma_ion_O_4S_cm2( 267) =  1.01 ; given_sigma_ion_O_2D_cm2( 267) =    1.28 ; given_sigma_ion_O_2P_cm2( 267) =    0.83 ; given_sigma_ion_O_4Pst_cm2( 267) =    0.38 ; given_sigma_ion_O_2Pst_cm2( 267) =    0.30
  given_sigma_ion_O_4S_cm2( 268) =  1.02 ; given_sigma_ion_O_2D_cm2( 268) =    1.28 ; given_sigma_ion_O_2P_cm2( 268) =    0.83 ; given_sigma_ion_O_4Pst_cm2( 268) =    0.38 ; given_sigma_ion_O_2Pst_cm2( 268) =    0.30
  given_sigma_ion_O_4S_cm2( 269) =  1.03 ; given_sigma_ion_O_2D_cm2( 269) =    1.29 ; given_sigma_ion_O_2P_cm2( 269) =    0.84 ; given_sigma_ion_O_4Pst_cm2( 269) =    0.38 ; given_sigma_ion_O_2Pst_cm2( 269) =    0.30
  given_sigma_ion_O_4S_cm2( 270) =  1.03 ; given_sigma_ion_O_2D_cm2( 270) =    1.30 ; given_sigma_ion_O_2P_cm2( 270) =    0.84 ; given_sigma_ion_O_4Pst_cm2( 270) =    0.38 ; given_sigma_ion_O_2Pst_cm2( 270) =    0.30
  given_sigma_ion_O_4S_cm2( 271) =  1.03 ; given_sigma_ion_O_2D_cm2( 271) =    1.30 ; given_sigma_ion_O_2P_cm2( 271) =    0.84 ; given_sigma_ion_O_4Pst_cm2( 271) =    0.38 ; given_sigma_ion_O_2Pst_cm2( 271) =    0.31
  given_sigma_ion_O_4S_cm2( 272) =  1.05 ; given_sigma_ion_O_2D_cm2( 272) =    1.32 ; given_sigma_ion_O_2P_cm2( 272) =    0.85 ; given_sigma_ion_O_4Pst_cm2( 272) =    0.39 ; given_sigma_ion_O_2Pst_cm2( 272) =    0.31
  given_sigma_ion_O_4S_cm2( 273) =  1.05 ; given_sigma_ion_O_2D_cm2( 273) =    1.32 ; given_sigma_ion_O_2P_cm2( 273) =    0.86 ; given_sigma_ion_O_4Pst_cm2( 273) =    0.39 ; given_sigma_ion_O_2Pst_cm2( 273) =    0.31
  given_sigma_ion_O_4S_cm2( 274) =  1.05 ; given_sigma_ion_O_2D_cm2( 274) =    1.33 ; given_sigma_ion_O_2P_cm2( 274) =    0.86 ; given_sigma_ion_O_4Pst_cm2( 274) =    0.39 ; given_sigma_ion_O_2Pst_cm2( 274) =    0.31
  given_sigma_ion_O_4S_cm2( 275) =  1.06 ; given_sigma_ion_O_2D_cm2( 275) =    1.33 ; given_sigma_ion_O_2P_cm2( 275) =    0.86 ; given_sigma_ion_O_4Pst_cm2( 275) =    0.39 ; given_sigma_ion_O_2Pst_cm2( 275) =    0.31
  given_sigma_ion_O_4S_cm2( 276) =  1.07 ; given_sigma_ion_O_2D_cm2( 276) =    1.34 ; given_sigma_ion_O_2P_cm2( 276) =    0.87 ; given_sigma_ion_O_4Pst_cm2( 276) =    0.40 ; given_sigma_ion_O_2Pst_cm2( 276) =    0.32
  given_sigma_ion_O_4S_cm2( 277) =  1.07 ; given_sigma_ion_O_2D_cm2( 277) =    1.35 ; given_sigma_ion_O_2P_cm2( 277) =    0.87 ; given_sigma_ion_O_4Pst_cm2( 277) =    0.40 ; given_sigma_ion_O_2Pst_cm2( 277) =    0.32
  given_sigma_ion_O_4S_cm2( 278) =  1.08 ; given_sigma_ion_O_2D_cm2( 278) =    1.36 ; given_sigma_ion_O_2P_cm2( 278) =    0.88 ; given_sigma_ion_O_4Pst_cm2( 278) =    0.40 ; given_sigma_ion_O_2Pst_cm2( 278) =    0.32
  given_sigma_ion_O_4S_cm2( 279) =  1.09 ; given_sigma_ion_O_2D_cm2( 279) =    1.37 ; given_sigma_ion_O_2P_cm2( 279) =    0.89 ; given_sigma_ion_O_4Pst_cm2( 279) =    0.40 ; given_sigma_ion_O_2Pst_cm2( 279) =    0.32
  given_sigma_ion_O_4S_cm2( 280) =  1.09 ; given_sigma_ion_O_2D_cm2( 280) =    1.38 ; given_sigma_ion_O_2P_cm2( 280) =    0.89 ; given_sigma_ion_O_4Pst_cm2( 280) =    0.40 ; given_sigma_ion_O_2Pst_cm2( 280) =    0.32
  given_sigma_ion_O_4S_cm2( 281) =  1.11 ; given_sigma_ion_O_2D_cm2( 281) =    1.39 ; given_sigma_ion_O_2P_cm2( 281) =    0.90 ; given_sigma_ion_O_4Pst_cm2( 281) =    0.41 ; given_sigma_ion_O_2Pst_cm2( 281) =    0.33
  given_sigma_ion_O_4S_cm2( 282) =  1.11 ; given_sigma_ion_O_2D_cm2( 282) =    1.39 ; given_sigma_ion_O_2P_cm2( 282) =    0.90 ; given_sigma_ion_O_4Pst_cm2( 282) =    0.41 ; given_sigma_ion_O_2Pst_cm2( 282) =    0.33
  given_sigma_ion_O_4S_cm2( 283) =  1.11 ; given_sigma_ion_O_2D_cm2( 283) =    1.40 ; given_sigma_ion_O_2P_cm2( 283) =    0.91 ; given_sigma_ion_O_4Pst_cm2( 283) =    0.41 ; given_sigma_ion_O_2Pst_cm2( 283) =    0.33
  given_sigma_ion_O_4S_cm2( 284) =  1.12 ; given_sigma_ion_O_2D_cm2( 284) =    1.42 ; given_sigma_ion_O_2P_cm2( 284) =    0.92 ; given_sigma_ion_O_4Pst_cm2( 284) =    0.42 ; given_sigma_ion_O_2Pst_cm2( 284) =    0.33
  given_sigma_ion_O_4S_cm2( 285) =  1.14 ; given_sigma_ion_O_2D_cm2( 285) =    1.43 ; given_sigma_ion_O_2P_cm2( 285) =    0.93 ; given_sigma_ion_O_4Pst_cm2( 285) =    0.42 ; given_sigma_ion_O_2Pst_cm2( 285) =    0.34
  given_sigma_ion_O_4S_cm2( 286) =  1.11 ; given_sigma_ion_O_2D_cm2( 286) =    1.45 ; given_sigma_ion_O_2P_cm2( 286) =    0.94 ; given_sigma_ion_O_4Pst_cm2( 286) =    0.43 ; given_sigma_ion_O_2Pst_cm2( 286) =    0.34
  given_sigma_ion_O_4S_cm2( 287) =  1.11 ; given_sigma_ion_O_2D_cm2( 287) =    1.46 ; given_sigma_ion_O_2P_cm2( 287) =    0.94 ; given_sigma_ion_O_4Pst_cm2( 287) =    0.43 ; given_sigma_ion_O_2Pst_cm2( 287) =    0.34
  given_sigma_ion_O_4S_cm2( 288) =  1.12 ; given_sigma_ion_O_2D_cm2( 288) =    1.46 ; given_sigma_ion_O_2P_cm2( 288) =    0.95 ; given_sigma_ion_O_4Pst_cm2( 288) =    0.43 ; given_sigma_ion_O_2Pst_cm2( 288) =    0.34
  given_sigma_ion_O_4S_cm2( 289) =  1.13 ; given_sigma_ion_O_2D_cm2( 289) =    1.48 ; given_sigma_ion_O_2P_cm2( 289) =    0.96 ; given_sigma_ion_O_4Pst_cm2( 289) =    0.43 ; given_sigma_ion_O_2Pst_cm2( 289) =    0.35
  given_sigma_ion_O_4S_cm2( 290) =  1.13 ; given_sigma_ion_O_2D_cm2( 290) =    1.48 ; given_sigma_ion_O_2P_cm2( 290) =    0.96 ; given_sigma_ion_O_4Pst_cm2( 290) =    0.44 ; given_sigma_ion_O_2Pst_cm2( 290) =    0.35
  given_sigma_ion_O_4S_cm2( 291) =  1.14 ; given_sigma_ion_O_2D_cm2( 291) =    1.49 ; given_sigma_ion_O_2P_cm2( 291) =    0.96 ; given_sigma_ion_O_4Pst_cm2( 291) =    0.44 ; given_sigma_ion_O_2Pst_cm2( 291) =    0.35
  given_sigma_ion_O_4S_cm2( 292) =  1.15 ; given_sigma_ion_O_2D_cm2( 292) =    1.51 ; given_sigma_ion_O_2P_cm2( 292) =    0.97 ; given_sigma_ion_O_4Pst_cm2( 292) =    0.44 ; given_sigma_ion_O_2Pst_cm2( 292) =    0.35
  given_sigma_ion_O_4S_cm2( 293) =  1.15 ; given_sigma_ion_O_2D_cm2( 293) =    1.51 ; given_sigma_ion_O_2P_cm2( 293) =    0.98 ; given_sigma_ion_O_4Pst_cm2( 293) =    0.44 ; given_sigma_ion_O_2Pst_cm2( 293) =    0.35
  given_sigma_ion_O_4S_cm2( 294) =  1.15 ; given_sigma_ion_O_2D_cm2( 294) =    1.51 ; given_sigma_ion_O_2P_cm2( 294) =    0.98 ; given_sigma_ion_O_4Pst_cm2( 294) =    0.44 ; given_sigma_ion_O_2Pst_cm2( 294) =    0.36
  given_sigma_ion_O_4S_cm2( 295) =  1.16 ; given_sigma_ion_O_2D_cm2( 295) =    1.52 ; given_sigma_ion_O_2P_cm2( 295) =    0.98 ; given_sigma_ion_O_4Pst_cm2( 295) =    0.45 ; given_sigma_ion_O_2Pst_cm2( 295) =    0.36
  given_sigma_ion_O_4S_cm2( 296) =  1.17 ; given_sigma_ion_O_2D_cm2( 296) =    1.53 ; given_sigma_ion_O_2P_cm2( 296) =    0.99 ; given_sigma_ion_O_4Pst_cm2( 296) =    0.45 ; given_sigma_ion_O_2Pst_cm2( 296) =    0.36
  given_sigma_ion_O_4S_cm2( 297) =  1.18 ; given_sigma_ion_O_2D_cm2( 297) =    1.54 ; given_sigma_ion_O_2P_cm2( 297) =    1.00 ; given_sigma_ion_O_4Pst_cm2( 297) =    0.45 ; given_sigma_ion_O_2Pst_cm2( 297) =    0.36
  given_sigma_ion_O_4S_cm2( 298) =  1.18 ; given_sigma_ion_O_2D_cm2( 298) =    1.54 ; given_sigma_ion_O_2P_cm2( 298) =    1.00 ; given_sigma_ion_O_4Pst_cm2( 298) =    0.45 ; given_sigma_ion_O_2Pst_cm2( 298) =    0.36
  given_sigma_ion_O_4S_cm2( 299) =  1.18 ; given_sigma_ion_O_2D_cm2( 299) =    1.55 ; given_sigma_ion_O_2P_cm2( 299) =    1.00 ; given_sigma_ion_O_4Pst_cm2( 299) =    0.45 ; given_sigma_ion_O_2Pst_cm2( 299) =    0.36
  given_sigma_ion_O_4S_cm2( 300) =  1.19 ; given_sigma_ion_O_2D_cm2( 300) =    1.56 ; given_sigma_ion_O_2P_cm2( 300) =    1.01 ; given_sigma_ion_O_4Pst_cm2( 300) =    0.46 ; given_sigma_ion_O_2Pst_cm2( 300) =    0.37
  given_sigma_ion_O_4S_cm2( 301) =  1.20 ; given_sigma_ion_O_2D_cm2( 301) =    1.57 ; given_sigma_ion_O_2P_cm2( 301) =    1.02 ; given_sigma_ion_O_4Pst_cm2( 301) =    0.46 ; given_sigma_ion_O_2Pst_cm2( 301) =    0.32
  given_sigma_ion_O_4S_cm2( 302) =  1.21 ; given_sigma_ion_O_2D_cm2( 302) =    1.59 ; given_sigma_ion_O_2P_cm2( 302) =    1.03 ; given_sigma_ion_O_4Pst_cm2( 302) =    0.47 ; given_sigma_ion_O_2Pst_cm2( 302) =    0.36
  given_sigma_ion_O_4S_cm2( 303) =  1.22 ; given_sigma_ion_O_2D_cm2( 303) =    1.60 ; given_sigma_ion_O_2P_cm2( 303) =    1.03 ; given_sigma_ion_O_4Pst_cm2( 303) =    0.47 ; given_sigma_ion_O_2Pst_cm2( 303) =    0.38
  given_sigma_ion_O_4S_cm2( 304) =  1.22 ; given_sigma_ion_O_2D_cm2( 304) =    1.60 ; given_sigma_ion_O_2P_cm2( 304) =    1.04 ; given_sigma_ion_O_4Pst_cm2( 304) =    0.47 ; given_sigma_ion_O_2Pst_cm2( 304) =    0.38
  given_sigma_ion_O_4S_cm2( 305) =  1.23 ; given_sigma_ion_O_2D_cm2( 305) =    1.60 ; given_sigma_ion_O_2P_cm2( 305) =    1.04 ; given_sigma_ion_O_4Pst_cm2( 305) =    0.47 ; given_sigma_ion_O_2Pst_cm2( 305) =    0.38
  given_sigma_ion_O_4S_cm2( 306) =  1.24 ; given_sigma_ion_O_2D_cm2( 306) =    1.62 ; given_sigma_ion_O_2P_cm2( 306) =    1.05 ; given_sigma_ion_O_4Pst_cm2( 306) =    0.48 ; given_sigma_ion_O_2Pst_cm2( 306) =    0.33
  given_sigma_ion_O_4S_cm2( 307) =  1.24 ; given_sigma_ion_O_2D_cm2( 307) =    1.62 ; given_sigma_ion_O_2P_cm2( 307) =    1.05 ; given_sigma_ion_O_4Pst_cm2( 307) =    0.48 ; given_sigma_ion_O_2Pst_cm2( 307) =    0.33
  given_sigma_ion_O_4S_cm2( 308) =  1.25 ; given_sigma_ion_O_2D_cm2( 308) =    1.64 ; given_sigma_ion_O_2P_cm2( 308) =    1.06 ; given_sigma_ion_O_4Pst_cm2( 308) =    0.48 ; given_sigma_ion_O_2Pst_cm2( 308) =    0.34
  given_sigma_ion_O_4S_cm2( 309) =  1.26 ; given_sigma_ion_O_2D_cm2( 309) =    1.65 ; given_sigma_ion_O_2P_cm2( 309) =    1.07 ; given_sigma_ion_O_4Pst_cm2( 309) =    0.48 ; given_sigma_ion_O_2Pst_cm2( 309) =    0.34
  given_sigma_ion_O_4S_cm2( 310) =  1.27 ; given_sigma_ion_O_2D_cm2( 310) =    1.66 ; given_sigma_ion_O_2P_cm2( 310) =    1.07 ; given_sigma_ion_O_4Pst_cm2( 310) =    0.49 ; given_sigma_ion_O_2Pst_cm2( 310) =    0.34
  given_sigma_ion_O_4S_cm2( 311) =  1.28 ; given_sigma_ion_O_2D_cm2( 311) =    1.72 ; given_sigma_ion_O_2P_cm2( 311) =    1.08 ; given_sigma_ion_O_4Pst_cm2( 311) =    0.49 ; given_sigma_ion_O_2Pst_cm2( 311) =    0.34
  given_sigma_ion_O_4S_cm2( 312) =  1.28 ; given_sigma_ion_O_2D_cm2( 312) =    1.72 ; given_sigma_ion_O_2P_cm2( 312) =    1.08 ; given_sigma_ion_O_4Pst_cm2( 312) =    0.49 ; given_sigma_ion_O_2Pst_cm2( 312) =    0.34
  given_sigma_ion_O_4S_cm2( 313) =  1.28 ; given_sigma_ion_O_2D_cm2( 313) =    1.72 ; given_sigma_ion_O_2P_cm2( 313) =    1.08 ; given_sigma_ion_O_4Pst_cm2( 313) =    0.49 ; given_sigma_ion_O_2Pst_cm2( 313) =    0.35
  given_sigma_ion_O_4S_cm2( 314) =  1.29 ; given_sigma_ion_O_2D_cm2( 314) =    1.69 ; given_sigma_ion_O_2P_cm2( 314) =    1.09 ; given_sigma_ion_O_4Pst_cm2( 314) =    0.50 ; given_sigma_ion_O_2Pst_cm2( 314) =    0.35
  given_sigma_ion_O_4S_cm2( 315) =  1.30 ; given_sigma_ion_O_2D_cm2( 315) =    1.75 ; given_sigma_ion_O_2P_cm2( 315) =    1.10 ; given_sigma_ion_O_4Pst_cm2( 315) =    0.50 ; given_sigma_ion_O_2Pst_cm2( 315) =    0.35
  given_sigma_ion_O_4S_cm2( 316) =  1.31 ; given_sigma_ion_O_2D_cm2( 316) =    1.76 ; given_sigma_ion_O_2P_cm2( 316) =    1.11 ; given_sigma_ion_O_4Pst_cm2( 316) =    0.50 ; given_sigma_ion_O_2Pst_cm2( 316) =    0.35
  given_sigma_ion_O_4S_cm2( 317) =  1.31 ; given_sigma_ion_O_2D_cm2( 317) =    1.76 ; given_sigma_ion_O_2P_cm2( 317) =    1.11 ; given_sigma_ion_O_4Pst_cm2( 317) =    0.50 ; given_sigma_ion_O_2Pst_cm2( 317) =    0.35
  given_sigma_ion_O_4S_cm2( 318) =  1.31 ; given_sigma_ion_O_2D_cm2( 318) =    1.76 ; given_sigma_ion_O_2P_cm2( 318) =    1.11 ; given_sigma_ion_O_4Pst_cm2( 318) =    0.50 ; given_sigma_ion_O_2Pst_cm2( 318) =    0.35
  given_sigma_ion_O_4S_cm2( 319) =  1.32 ; given_sigma_ion_O_2D_cm2( 319) =    1.78 ; given_sigma_ion_O_2P_cm2( 319) =    1.12 ; given_sigma_ion_O_4Pst_cm2( 319) =    0.51 ; given_sigma_ion_O_2Pst_cm2( 319) =    0.36
  given_sigma_ion_O_4S_cm2( 320) =  1.33 ; given_sigma_ion_O_2D_cm2( 320) =    1.79 ; given_sigma_ion_O_2P_cm2( 320) =    1.12 ; given_sigma_ion_O_4Pst_cm2( 320) =    0.51 ; given_sigma_ion_O_2Pst_cm2( 320) =    0.36
  given_sigma_ion_O_4S_cm2( 321) =  1.33 ; given_sigma_ion_O_2D_cm2( 321) =    1.79 ; given_sigma_ion_O_2P_cm2( 321) =    1.12 ; given_sigma_ion_O_4Pst_cm2( 321) =    0.51 ; given_sigma_ion_O_2Pst_cm2( 321) =    0.36
  given_sigma_ion_O_4S_cm2( 322) =  1.34 ; given_sigma_ion_O_2D_cm2( 322) =    1.80 ; given_sigma_ion_O_2P_cm2( 322) =    1.13 ; given_sigma_ion_O_4Pst_cm2( 322) =    0.51 ; given_sigma_ion_O_2Pst_cm2( 322) =    0.36
  given_sigma_ion_O_4S_cm2( 323) =  1.35 ; given_sigma_ion_O_2D_cm2( 323) =    1.81 ; given_sigma_ion_O_2P_cm2( 323) =    1.14 ; given_sigma_ion_O_4Pst_cm2( 323) =    0.52 ; given_sigma_ion_O_2Pst_cm2( 323) =    0.36
  given_sigma_ion_O_4S_cm2( 324) =  1.35 ; given_sigma_ion_O_2D_cm2( 324) =    1.82 ; given_sigma_ion_O_2P_cm2( 324) =    1.14 ; given_sigma_ion_O_4Pst_cm2( 324) =    0.52 ; given_sigma_ion_O_2Pst_cm2( 324) =    0.36
  given_sigma_ion_O_4S_cm2( 325) =  1.36 ; given_sigma_ion_O_2D_cm2( 325) =    1.82 ; given_sigma_ion_O_2P_cm2( 325) =    1.15 ; given_sigma_ion_O_4Pst_cm2( 325) =    0.52 ; given_sigma_ion_O_2Pst_cm2( 325) =    0.36
  given_sigma_ion_O_4S_cm2( 326) =  1.36 ; given_sigma_ion_O_2D_cm2( 326) =    1.83 ; given_sigma_ion_O_2P_cm2( 326) =    1.21 ; given_sigma_ion_O_4Pst_cm2( 326) =    0.52 ; given_sigma_ion_O_2Pst_cm2( 326) =    0.37
  given_sigma_ion_O_4S_cm2( 327) =  1.37 ; given_sigma_ion_O_2D_cm2( 327) =    1.85 ; given_sigma_ion_O_2P_cm2( 327) =    1.21 ; given_sigma_ion_O_4Pst_cm2( 327) =    0.53 ; given_sigma_ion_O_2Pst_cm2( 327) =    0.37
  given_sigma_ion_O_4S_cm2( 328) =  1.38 ; given_sigma_ion_O_2D_cm2( 328) =    1.86 ; given_sigma_ion_O_2P_cm2( 328) =    1.22 ; given_sigma_ion_O_4Pst_cm2( 328) =    0.53 ; given_sigma_ion_O_2Pst_cm2( 328) =    0.37
  given_sigma_ion_O_4S_cm2( 329) =  1.39 ; given_sigma_ion_O_2D_cm2( 329) =    1.87 ; given_sigma_ion_O_2P_cm2( 329) =    1.23 ; given_sigma_ion_O_4Pst_cm2( 329) =    0.53 ; given_sigma_ion_O_2Pst_cm2( 329) =    0.37
  given_sigma_ion_O_4S_cm2( 330) =  1.39 ; given_sigma_ion_O_2D_cm2( 330) =    1.87 ; given_sigma_ion_O_2P_cm2( 330) =    1.23 ; given_sigma_ion_O_4Pst_cm2( 330) =    0.54 ; given_sigma_ion_O_2Pst_cm2( 330) =    0.37
  given_sigma_ion_O_4S_cm2( 331) =  1.40 ; given_sigma_ion_O_2D_cm2( 331) =    1.88 ; given_sigma_ion_O_2P_cm2( 331) =    1.24 ; given_sigma_ion_O_4Pst_cm2( 331) =    0.54 ; given_sigma_ion_O_2Pst_cm2( 331) =    0.38
  given_sigma_ion_O_4S_cm2( 332) =  1.41 ; given_sigma_ion_O_2D_cm2( 332) =    1.90 ; given_sigma_ion_O_2P_cm2( 332) =    1.25 ; given_sigma_ion_O_4Pst_cm2( 332) =    0.54 ; given_sigma_ion_O_2Pst_cm2( 332) =    0.38
  given_sigma_ion_O_4S_cm2( 333) =  1.42 ; given_sigma_ion_O_2D_cm2( 333) =    1.91 ; given_sigma_ion_O_2P_cm2( 333) =    1.26 ; given_sigma_ion_O_4Pst_cm2( 333) =    0.55 ; given_sigma_ion_O_2Pst_cm2( 333) =    0.38
  given_sigma_ion_O_4S_cm2( 334) =  1.43 ; given_sigma_ion_O_2D_cm2( 334) =    1.93 ; given_sigma_ion_O_2P_cm2( 334) =    1.27 ; given_sigma_ion_O_4Pst_cm2( 334) =    0.55 ; given_sigma_ion_O_2Pst_cm2( 334) =    0.39
  given_sigma_ion_O_4S_cm2( 335) =  1.43 ; given_sigma_ion_O_2D_cm2( 335) =    1.93 ; given_sigma_ion_O_2P_cm2( 335) =    1.27 ; given_sigma_ion_O_4Pst_cm2( 335) =    0.55 ; given_sigma_ion_O_2Pst_cm2( 335) =    0.39
  given_sigma_ion_O_4S_cm2( 336) =  1.44 ; given_sigma_ion_O_2D_cm2( 336) =    1.94 ; given_sigma_ion_O_2P_cm2( 336) =    1.27 ; given_sigma_ion_O_4Pst_cm2( 336) =    0.55 ; given_sigma_ion_O_2Pst_cm2( 336) =    0.39
  given_sigma_ion_O_4S_cm2( 337) =  1.42 ; given_sigma_ion_O_2D_cm2( 337) =    1.95 ; given_sigma_ion_O_2P_cm2( 337) =    1.28 ; given_sigma_ion_O_4Pst_cm2( 337) =    0.56 ; given_sigma_ion_O_2Pst_cm2( 337) =    0.39
  given_sigma_ion_O_4S_cm2( 338) =  1.40 ; given_sigma_ion_O_2D_cm2( 338) =    1.96 ; given_sigma_ion_O_2P_cm2( 338) =    1.29 ; given_sigma_ion_O_4Pst_cm2( 338) =    0.56 ; given_sigma_ion_O_2Pst_cm2( 338) =    0.39
  given_sigma_ion_O_4S_cm2( 339) =  1.41 ; given_sigma_ion_O_2D_cm2( 339) =    1.97 ; given_sigma_ion_O_2P_cm2( 339) =    1.30 ; given_sigma_ion_O_4Pst_cm2( 339) =    0.56 ; given_sigma_ion_O_2Pst_cm2( 339) =    0.39
  given_sigma_ion_O_4S_cm2( 340) =  1.42 ; given_sigma_ion_O_2D_cm2( 340) =    1.98 ; given_sigma_ion_O_2P_cm2( 340) =    1.30 ; given_sigma_ion_O_4Pst_cm2( 340) =    0.57 ; given_sigma_ion_O_2Pst_cm2( 340) =    0.40
  given_sigma_ion_O_4S_cm2( 341) =  1.43 ; given_sigma_ion_O_2D_cm2( 341) =    2.00 ; given_sigma_ion_O_2P_cm2( 341) =    1.31 ; given_sigma_ion_O_4Pst_cm2( 341) =    0.57 ; given_sigma_ion_O_2Pst_cm2( 341) =    0.40
  given_sigma_ion_O_4S_cm2( 342) =  1.49 ; given_sigma_ion_O_2D_cm2( 342) =    2.00 ; given_sigma_ion_O_2P_cm2( 342) =    1.31 ; given_sigma_ion_O_4Pst_cm2( 342) =    0.57 ; given_sigma_ion_O_2Pst_cm2( 342) =    0.40
  given_sigma_ion_O_4S_cm2( 343) =  1.49 ; given_sigma_ion_O_2D_cm2( 343) =    2.01 ; given_sigma_ion_O_2P_cm2( 343) =    1.32 ; given_sigma_ion_O_4Pst_cm2( 343) =    0.57 ; given_sigma_ion_O_2Pst_cm2( 343) =    0.40
  given_sigma_ion_O_4S_cm2( 344) =  1.49 ; given_sigma_ion_O_2D_cm2( 344) =    2.01 ; given_sigma_ion_O_2P_cm2( 344) =    1.32 ; given_sigma_ion_O_4Pst_cm2( 344) =    0.57 ; given_sigma_ion_O_2Pst_cm2( 344) =    0.40
  given_sigma_ion_O_4S_cm2( 345) =  1.50 ; given_sigma_ion_O_2D_cm2( 345) =    2.04 ; given_sigma_ion_O_2P_cm2( 345) =    1.33 ; given_sigma_ion_O_4Pst_cm2( 345) =    0.58 ; given_sigma_ion_O_2Pst_cm2( 345) =    0.40
  given_sigma_ion_O_4S_cm2( 346) =  1.51 ; given_sigma_ion_O_2D_cm2( 346) =    2.09 ; given_sigma_ion_O_2P_cm2( 346) =    1.34 ; given_sigma_ion_O_4Pst_cm2( 346) =    0.58 ; given_sigma_ion_O_2Pst_cm2( 346) =    0.41
  given_sigma_ion_O_4S_cm2( 347) =  1.52 ; given_sigma_ion_O_2D_cm2( 347) =    2.10 ; given_sigma_ion_O_2P_cm2( 347) =    1.34 ; given_sigma_ion_O_4Pst_cm2( 347) =    0.58 ; given_sigma_ion_O_2Pst_cm2( 347) =    0.41
  given_sigma_ion_O_4S_cm2( 348) =  1.53 ; given_sigma_ion_O_2D_cm2( 348) =    2.12 ; given_sigma_ion_O_2P_cm2( 348) =    1.35 ; given_sigma_ion_O_4Pst_cm2( 348) =    0.59 ; given_sigma_ion_O_2Pst_cm2( 348) =    0.41
  given_sigma_ion_O_4S_cm2( 349) =  1.48 ; given_sigma_ion_O_2D_cm2( 349) =    2.13 ; given_sigma_ion_O_2P_cm2( 349) =    1.36 ; given_sigma_ion_O_4Pst_cm2( 349) =    0.59 ; given_sigma_ion_O_2Pst_cm2( 349) =    0.41
  given_sigma_ion_O_4S_cm2( 350) =  1.54 ; given_sigma_ion_O_2D_cm2( 350) =    2.13 ; given_sigma_ion_O_2P_cm2( 350) =    1.36 ; given_sigma_ion_O_4Pst_cm2( 350) =    0.59 ; given_sigma_ion_O_2Pst_cm2( 350) =    0.41
  given_sigma_ion_O_4S_cm2( 351) =  1.55 ; given_sigma_ion_O_2D_cm2( 351) =    2.15 ; given_sigma_ion_O_2P_cm2( 351) =    1.37 ; given_sigma_ion_O_4Pst_cm2( 351) =    0.60 ; given_sigma_ion_O_2Pst_cm2( 351) =    0.42
  given_sigma_ion_O_4S_cm2( 352) =  1.56 ; given_sigma_ion_O_2D_cm2( 352) =    2.16 ; given_sigma_ion_O_2P_cm2( 352) =    1.38 ; given_sigma_ion_O_4Pst_cm2( 352) =    0.60 ; given_sigma_ion_O_2Pst_cm2( 352) =    0.42
  given_sigma_ion_O_4S_cm2( 353) =  1.56 ; given_sigma_ion_O_2D_cm2( 353) =    2.17 ; given_sigma_ion_O_2P_cm2( 353) =    1.38 ; given_sigma_ion_O_4Pst_cm2( 353) =    0.60 ; given_sigma_ion_O_2Pst_cm2( 353) =    0.42
  given_sigma_ion_O_4S_cm2( 354) =  1.57 ; given_sigma_ion_O_2D_cm2( 354) =    2.18 ; given_sigma_ion_O_2P_cm2( 354) =    1.39 ; given_sigma_ion_O_4Pst_cm2( 354) =    0.61 ; given_sigma_ion_O_2Pst_cm2( 354) =    0.42
  given_sigma_ion_O_4S_cm2( 355) =  1.58 ; given_sigma_ion_O_2D_cm2( 355) =    2.18 ; given_sigma_ion_O_2P_cm2( 355) =    1.39 ; given_sigma_ion_O_4Pst_cm2( 355) =    0.61 ; given_sigma_ion_O_2Pst_cm2( 355) =    0.42
  given_sigma_ion_O_4S_cm2( 356) =  1.58 ; given_sigma_ion_O_2D_cm2( 356) =    2.19 ; given_sigma_ion_O_2P_cm2( 356) =    1.40 ; given_sigma_ion_O_4Pst_cm2( 356) =    0.61 ; given_sigma_ion_O_2Pst_cm2( 356) =    0.39
  given_sigma_ion_O_4S_cm2( 357) =  1.58 ; given_sigma_ion_O_2D_cm2( 357) =    2.19 ; given_sigma_ion_O_2P_cm2( 357) =    1.40 ; given_sigma_ion_O_4Pst_cm2( 357) =    0.61 ; given_sigma_ion_O_2Pst_cm2( 357) =    0.37
  given_sigma_ion_O_4S_cm2( 358) =  1.59 ; given_sigma_ion_O_2D_cm2( 358) =    2.20 ; given_sigma_ion_O_2P_cm2( 358) =    1.40 ; given_sigma_ion_O_4Pst_cm2( 358) =    0.61 ; given_sigma_ion_O_2Pst_cm2( 358) =    0.37
  given_sigma_ion_O_4S_cm2( 359) =  1.53 ; given_sigma_ion_O_2D_cm2( 359) =    2.21 ; given_sigma_ion_O_2P_cm2( 359) =    1.41 ; given_sigma_ion_O_4Pst_cm2( 359) =    0.61 ; given_sigma_ion_O_2Pst_cm2( 359) =    0.37
  given_sigma_ion_O_4S_cm2( 360) =  1.54 ; given_sigma_ion_O_2D_cm2( 360) =    2.21 ; given_sigma_ion_O_2P_cm2( 360) =    1.41 ; given_sigma_ion_O_4Pst_cm2( 360) =    0.62 ; given_sigma_ion_O_2Pst_cm2( 360) =    0.37
  given_sigma_ion_O_4S_cm2( 361) =  1.55 ; given_sigma_ion_O_2D_cm2( 361) =    2.23 ; given_sigma_ion_O_2P_cm2( 361) =    1.42 ; given_sigma_ion_O_4Pst_cm2( 361) =    0.62 ; given_sigma_ion_O_2Pst_cm2( 361) =    0.37
  given_sigma_ion_O_4S_cm2( 362) =  1.56 ; given_sigma_ion_O_2D_cm2( 362) =    2.25 ; given_sigma_ion_O_2P_cm2( 362) =    1.44 ; given_sigma_ion_O_4Pst_cm2( 362) =    0.63 ; given_sigma_ion_O_2Pst_cm2( 362) =    0.38
  given_sigma_ion_O_4S_cm2( 363) =  1.57 ; given_sigma_ion_O_2D_cm2( 363) =    2.27 ; given_sigma_ion_O_2P_cm2( 363) =    1.45 ; given_sigma_ion_O_4Pst_cm2( 363) =    0.63 ; given_sigma_ion_O_2Pst_cm2( 363) =    0.38
  given_sigma_ion_O_4S_cm2( 364) =  1.58 ; given_sigma_ion_O_2D_cm2( 364) =    2.27 ; given_sigma_ion_O_2P_cm2( 364) =    1.45 ; given_sigma_ion_O_4Pst_cm2( 364) =    0.63 ; given_sigma_ion_O_2Pst_cm2( 364) =    0.38
  given_sigma_ion_O_4S_cm2( 365) =  1.58 ; given_sigma_ion_O_2D_cm2( 365) =    2.28 ; given_sigma_ion_O_2P_cm2( 365) =    1.45 ; given_sigma_ion_O_4Pst_cm2( 365) =    0.63 ; given_sigma_ion_O_2Pst_cm2( 365) =    0.38
  given_sigma_ion_O_4S_cm2( 366) =  1.62 ; given_sigma_ion_O_2D_cm2( 366) =    2.33 ; given_sigma_ion_O_2P_cm2( 366) =    1.49 ; given_sigma_ion_O_4Pst_cm2( 366) =    0.65 ; given_sigma_ion_O_2Pst_cm2( 366) =    0.39
  given_sigma_ion_O_4S_cm2( 367) =  1.63 ; given_sigma_ion_O_2D_cm2( 367) =    2.34 ; given_sigma_ion_O_2P_cm2( 367) =    1.50 ; given_sigma_ion_O_4Pst_cm2( 367) =    0.65 ; given_sigma_ion_O_2Pst_cm2( 367) =    0.39
  given_sigma_ion_O_4S_cm2( 368) =  1.63 ; given_sigma_ion_O_2D_cm2( 368) =    2.35 ; given_sigma_ion_O_2P_cm2( 368) =    1.50 ; given_sigma_ion_O_4Pst_cm2( 368) =    0.65 ; given_sigma_ion_O_2Pst_cm2( 368) =    0.39
  given_sigma_ion_O_4S_cm2( 369) =  1.64 ; given_sigma_ion_O_2D_cm2( 369) =    2.37 ; given_sigma_ion_O_2P_cm2( 369) =    1.51 ; given_sigma_ion_O_4Pst_cm2( 369) =    0.66 ; given_sigma_ion_O_2Pst_cm2( 369) =    0.39
  given_sigma_ion_O_4S_cm2( 370) =  1.65 ; given_sigma_ion_O_2D_cm2( 370) =    2.37 ; given_sigma_ion_O_2P_cm2( 370) =    1.52 ; given_sigma_ion_O_4Pst_cm2( 370) =    0.66 ; given_sigma_ion_O_2Pst_cm2( 370) =    0.40
  given_sigma_ion_O_4S_cm2( 371) =  1.66 ; given_sigma_ion_O_2D_cm2( 371) =    2.39 ; given_sigma_ion_O_2P_cm2( 371) =    1.53 ; given_sigma_ion_O_4Pst_cm2( 371) =    0.66 ; given_sigma_ion_O_2Pst_cm2( 371) =    0.40
  given_sigma_ion_O_4S_cm2( 372) =  1.67 ; given_sigma_ion_O_2D_cm2( 372) =    2.40 ; given_sigma_ion_O_2P_cm2( 372) =    1.54 ; given_sigma_ion_O_4Pst_cm2( 372) =    0.67 ; given_sigma_ion_O_2Pst_cm2( 372) =    0.40
  given_sigma_ion_O_4S_cm2( 373) =  1.67 ; given_sigma_ion_O_2D_cm2( 373) =    2.41 ; given_sigma_ion_O_2P_cm2( 373) =    1.54 ; given_sigma_ion_O_4Pst_cm2( 373) =    0.67 ; given_sigma_ion_O_2Pst_cm2( 373) =    0.40
  given_sigma_ion_O_4S_cm2( 374) =  1.67 ; given_sigma_ion_O_2D_cm2( 374) =    2.41 ; given_sigma_ion_O_2P_cm2( 374) =    1.54 ; given_sigma_ion_O_4Pst_cm2( 374) =    0.67 ; given_sigma_ion_O_2Pst_cm2( 374) =    0.40
  given_sigma_ion_O_4S_cm2( 375) =  1.68 ; given_sigma_ion_O_2D_cm2( 375) =    2.42 ; given_sigma_ion_O_2P_cm2( 375) =    1.54 ; given_sigma_ion_O_4Pst_cm2( 375) =    0.67 ; given_sigma_ion_O_2Pst_cm2( 375) =    0.40
  given_sigma_ion_O_4S_cm2( 376) =  1.68 ; given_sigma_ion_O_2D_cm2( 376) =    2.43 ; given_sigma_ion_O_2P_cm2( 376) =    1.55 ; given_sigma_ion_O_4Pst_cm2( 376) =    0.67 ; given_sigma_ion_O_2Pst_cm2( 376) =    0.40
  given_sigma_ion_O_4S_cm2( 377) =  1.69 ; given_sigma_ion_O_2D_cm2( 377) =    2.43 ; given_sigma_ion_O_2P_cm2( 377) =    1.55 ; given_sigma_ion_O_4Pst_cm2( 377) =    0.67 ; given_sigma_ion_O_2Pst_cm2( 377) =    0.40
  given_sigma_ion_O_4S_cm2( 378) =  1.69 ; given_sigma_ion_O_2D_cm2( 378) =    2.43 ; given_sigma_ion_O_2P_cm2( 378) =    1.55 ; given_sigma_ion_O_4Pst_cm2( 378) =    0.68 ; given_sigma_ion_O_2Pst_cm2( 378) =    0.41
  given_sigma_ion_O_4S_cm2( 379) =  1.70 ; given_sigma_ion_O_2D_cm2( 379) =    2.45 ; given_sigma_ion_O_2P_cm2( 379) =    1.56 ; given_sigma_ion_O_4Pst_cm2( 379) =    0.68 ; given_sigma_ion_O_2Pst_cm2( 379) =    0.41
  given_sigma_ion_O_4S_cm2( 380) =  1.71 ; given_sigma_ion_O_2D_cm2( 380) =    2.47 ; given_sigma_ion_O_2P_cm2( 380) =    1.58 ; given_sigma_ion_O_4Pst_cm2( 380) =    0.69 ; given_sigma_ion_O_2Pst_cm2( 380) =    0.37
  given_sigma_ion_O_4S_cm2( 381) =  1.72 ; given_sigma_ion_O_2D_cm2( 381) =    2.48 ; given_sigma_ion_O_2P_cm2( 381) =    1.59 ; given_sigma_ion_O_4Pst_cm2( 381) =    0.69 ; given_sigma_ion_O_2Pst_cm2( 381) =    0.36
  given_sigma_ion_O_4S_cm2( 382) =  1.73 ; given_sigma_ion_O_2D_cm2( 382) =    2.49 ; given_sigma_ion_O_2P_cm2( 382) =    1.59 ; given_sigma_ion_O_4Pst_cm2( 382) =    0.69 ; given_sigma_ion_O_2Pst_cm2( 382) =    0.35
  given_sigma_ion_O_4S_cm2( 383) =  1.77 ; given_sigma_ion_O_2D_cm2( 383) =    2.62 ; given_sigma_ion_O_2P_cm2( 383) =    1.63 ; given_sigma_ion_O_4Pst_cm2( 383) =    0.71 ; given_sigma_ion_O_2Pst_cm2( 383) =    0.35
  given_sigma_ion_O_4S_cm2( 384) =  1.78 ; given_sigma_ion_O_2D_cm2( 384) =    2.64 ; given_sigma_ion_O_2P_cm2( 384) =    1.64 ; given_sigma_ion_O_4Pst_cm2( 384) =    0.71 ; given_sigma_ion_O_2Pst_cm2( 384) =    0.36
  given_sigma_ion_O_4S_cm2( 385) =  1.79 ; given_sigma_ion_O_2D_cm2( 385) =    2.65 ; given_sigma_ion_O_2P_cm2( 385) =    1.65 ; given_sigma_ion_O_4Pst_cm2( 385) =    0.72 ; given_sigma_ion_O_2Pst_cm2( 385) =    0.36
  given_sigma_ion_O_4S_cm2( 386) =  1.79 ; given_sigma_ion_O_2D_cm2( 386) =    2.65 ; given_sigma_ion_O_2P_cm2( 386) =    1.65 ; given_sigma_ion_O_4Pst_cm2( 386) =    0.72 ; given_sigma_ion_O_2Pst_cm2( 386) =    0.36
  given_sigma_ion_O_4S_cm2( 387) =  1.83 ; given_sigma_ion_O_2D_cm2( 387) =    2.70 ; given_sigma_ion_O_2P_cm2( 387) =    1.75 ; given_sigma_ion_O_4Pst_cm2( 387) =    0.73 ; given_sigma_ion_O_2Pst_cm2( 387) =    0.37
  given_sigma_ion_O_4S_cm2( 388) =  1.84 ; given_sigma_ion_O_2D_cm2( 388) =    2.72 ; given_sigma_ion_O_2P_cm2( 388) =    1.77 ; given_sigma_ion_O_4Pst_cm2( 388) =    0.66 ; given_sigma_ion_O_2Pst_cm2( 388) =    0.37
  given_sigma_ion_O_4S_cm2( 389) =  1.85 ; given_sigma_ion_O_2D_cm2( 389) =    2.74 ; given_sigma_ion_O_2P_cm2( 389) =    1.78 ; given_sigma_ion_O_4Pst_cm2( 389) =    0.71 ; given_sigma_ion_O_2Pst_cm2( 389) =    0.37
  given_sigma_ion_O_4S_cm2( 390) =  1.85 ; given_sigma_ion_O_2D_cm2( 390) =    2.74 ; given_sigma_ion_O_2P_cm2( 390) =    1.78 ; given_sigma_ion_O_4Pst_cm2( 390) =    0.74 ; given_sigma_ion_O_2Pst_cm2( 390) =    0.37
  given_sigma_ion_O_4S_cm2( 391) =  1.86 ; given_sigma_ion_O_2D_cm2( 391) =    2.75 ; given_sigma_ion_O_2P_cm2( 391) =    1.78 ; given_sigma_ion_O_4Pst_cm2( 391) =    0.74 ; given_sigma_ion_O_2Pst_cm2( 391) =    0.37
  given_sigma_ion_O_4S_cm2( 392) =  1.86 ; given_sigma_ion_O_2D_cm2( 392) =    2.75 ; given_sigma_ion_O_2P_cm2( 392) =    1.78 ; given_sigma_ion_O_4Pst_cm2( 392) =    0.74 ; given_sigma_ion_O_2Pst_cm2( 392) =    0.37
  given_sigma_ion_O_4S_cm2( 393) =  1.86 ; given_sigma_ion_O_2D_cm2( 393) =    2.75 ; given_sigma_ion_O_2P_cm2( 393) =    1.79 ; given_sigma_ion_O_4Pst_cm2( 393) =    0.74 ; given_sigma_ion_O_2Pst_cm2( 393) =    0.37
  given_sigma_ion_O_4S_cm2( 394) =  1.87 ; given_sigma_ion_O_2D_cm2( 394) =    2.77 ; given_sigma_ion_O_2P_cm2( 394) =    1.79 ; given_sigma_ion_O_4Pst_cm2( 394) =    0.69 ; given_sigma_ion_O_2Pst_cm2( 394) =    0.37
  given_sigma_ion_O_4S_cm2( 395) =  1.87 ; given_sigma_ion_O_2D_cm2( 395) =    2.77 ; given_sigma_ion_O_2P_cm2( 395) =    1.80 ; given_sigma_ion_O_4Pst_cm2( 395) =    0.67 ; given_sigma_ion_O_2Pst_cm2( 395) =    0.37
  given_sigma_ion_O_4S_cm2( 396) =  1.87 ; given_sigma_ion_O_2D_cm2( 396) =    2.77 ; given_sigma_ion_O_2P_cm2( 396) =    1.80 ; given_sigma_ion_O_4Pst_cm2( 396) =    0.67 ; given_sigma_ion_O_2Pst_cm2( 396) =    0.37
  given_sigma_ion_O_4S_cm2( 397) =  1.89 ; given_sigma_ion_O_2D_cm2( 397) =    2.79 ; given_sigma_ion_O_2P_cm2( 397) =    1.81 ; given_sigma_ion_O_4Pst_cm2( 397) =    0.68 ; given_sigma_ion_O_2Pst_cm2( 397) =    0.30
  given_sigma_ion_O_4S_cm2( 398) =  1.89 ; given_sigma_ion_O_2D_cm2( 398) =    2.79 ; given_sigma_ion_O_2P_cm2( 398) =    1.81 ; given_sigma_ion_O_4Pst_cm2( 398) =    0.68 ; given_sigma_ion_O_2Pst_cm2( 398) =    0.30
  given_sigma_ion_O_4S_cm2( 399) =  1.90 ; given_sigma_ion_O_2D_cm2( 399) =    2.81 ; given_sigma_ion_O_2P_cm2( 399) =    1.82 ; given_sigma_ion_O_4Pst_cm2( 399) =    0.68 ; given_sigma_ion_O_2Pst_cm2( 399) =    0.30
  given_sigma_ion_O_4S_cm2( 400) =  1.92 ; given_sigma_ion_O_2D_cm2( 400) =    2.84 ; given_sigma_ion_O_2P_cm2( 400) =    1.84 ; given_sigma_ion_O_4Pst_cm2( 400) =    0.69 ; given_sigma_ion_O_2Pst_cm2( 400) =    0.31
  given_sigma_ion_O_4S_cm2( 401) =  1.92 ; given_sigma_ion_O_2D_cm2( 401) =    2.85 ; given_sigma_ion_O_2P_cm2( 401) =    1.85 ; given_sigma_ion_O_4Pst_cm2( 401) =    0.69 ; given_sigma_ion_O_2Pst_cm2( 401) =    0.31
  given_sigma_ion_O_4S_cm2( 402) =  2.02 ; given_sigma_ion_O_2D_cm2( 402) =    3.01 ; given_sigma_ion_O_2P_cm2( 402) =    1.94 ; given_sigma_ion_O_4Pst_cm2( 402) =    0.75 ; given_sigma_ion_O_2Pst_cm2( 402) =    0.14
  given_sigma_ion_O_4S_cm2( 403) =  2.10 ; given_sigma_ion_O_2D_cm2( 403) =    3.15 ; given_sigma_ion_O_2P_cm2( 403) =    2.02 ; given_sigma_ion_O_4Pst_cm2( 403) =    0.81 ; given_sigma_ion_O_2Pst_cm2( 403) =    0.00
  given_sigma_ion_O_4S_cm2( 404) =  2.11 ; given_sigma_ion_O_2D_cm2( 404) =    3.16 ; given_sigma_ion_O_2P_cm2( 404) =    2.03 ; given_sigma_ion_O_4Pst_cm2( 404) =    0.81 ; given_sigma_ion_O_2Pst_cm2( 404) =    0.00
  given_sigma_ion_O_4S_cm2( 405) =  2.11 ; given_sigma_ion_O_2D_cm2( 405) =    3.17 ; given_sigma_ion_O_2P_cm2( 405) =    2.03 ; given_sigma_ion_O_4Pst_cm2( 405) =    0.81 ; given_sigma_ion_O_2Pst_cm2( 405) =    0.00
  given_sigma_ion_O_4S_cm2( 406) =  2.14 ; given_sigma_ion_O_2D_cm2( 406) =    3.20 ; given_sigma_ion_O_2P_cm2( 406) =    2.05 ; given_sigma_ion_O_4Pst_cm2( 406) =    0.82 ; given_sigma_ion_O_2Pst_cm2( 406) =    0.00
  given_sigma_ion_O_4S_cm2( 407) =  2.14 ; given_sigma_ion_O_2D_cm2( 407) =    3.22 ; given_sigma_ion_O_2P_cm2( 407) =    2.06 ; given_sigma_ion_O_4Pst_cm2( 407) =    0.82 ; given_sigma_ion_O_2Pst_cm2( 407) =    0.00
  given_sigma_ion_O_4S_cm2( 408) =  2.14 ; given_sigma_ion_O_2D_cm2( 408) =    3.22 ; given_sigma_ion_O_2P_cm2( 408) =    2.06 ; given_sigma_ion_O_4Pst_cm2( 408) =    0.82 ; given_sigma_ion_O_2Pst_cm2( 408) =    0.00
  given_sigma_ion_O_4S_cm2( 409) =  2.15 ; given_sigma_ion_O_2D_cm2( 409) =    3.23 ; given_sigma_ion_O_2P_cm2( 409) =    2.07 ; given_sigma_ion_O_4Pst_cm2( 409) =    0.82 ; given_sigma_ion_O_2Pst_cm2( 409) =    0.00
  given_sigma_ion_O_4S_cm2( 410) =  2.19 ; given_sigma_ion_O_2D_cm2( 410) =    3.31 ; given_sigma_ion_O_2P_cm2( 410) =    2.11 ; given_sigma_ion_O_4Pst_cm2( 410) =    0.81 ; given_sigma_ion_O_2Pst_cm2( 410) =    0.00
  given_sigma_ion_O_4S_cm2( 411) =  2.24 ; given_sigma_ion_O_2D_cm2( 411) =    3.41 ; given_sigma_ion_O_2P_cm2( 411) =    2.15 ; given_sigma_ion_O_4Pst_cm2( 411) =    0.80 ; given_sigma_ion_O_2Pst_cm2( 411) =    0.00
  given_sigma_ion_O_4S_cm2( 412) =  2.28 ; given_sigma_ion_O_2D_cm2( 412) =    3.50 ; given_sigma_ion_O_2P_cm2( 412) =    2.19 ; given_sigma_ion_O_4Pst_cm2( 412) =    0.79 ; given_sigma_ion_O_2Pst_cm2( 412) =    0.00
  given_sigma_ion_O_4S_cm2( 413) =  2.28 ; given_sigma_ion_O_2D_cm2( 413) =    3.51 ; given_sigma_ion_O_2P_cm2( 413) =    2.19 ; given_sigma_ion_O_4Pst_cm2( 413) =    0.79 ; given_sigma_ion_O_2Pst_cm2( 413) =    0.00
  given_sigma_ion_O_4S_cm2( 414) =  2.32 ; given_sigma_ion_O_2D_cm2( 414) =    3.57 ; given_sigma_ion_O_2P_cm2( 414) =    2.23 ; given_sigma_ion_O_4Pst_cm2( 414) =    0.80 ; given_sigma_ion_O_2Pst_cm2( 414) =    0.00
  given_sigma_ion_O_4S_cm2( 415) =  2.36 ; given_sigma_ion_O_2D_cm2( 415) =    3.63 ; given_sigma_ion_O_2P_cm2( 415) =    2.27 ; given_sigma_ion_O_4Pst_cm2( 415) =    0.82 ; given_sigma_ion_O_2Pst_cm2( 415) =    0.00
  given_sigma_ion_O_4S_cm2( 416) =  2.36 ; given_sigma_ion_O_2D_cm2( 416) =    3.64 ; given_sigma_ion_O_2P_cm2( 416) =    2.27 ; given_sigma_ion_O_4Pst_cm2( 416) =    0.82 ; given_sigma_ion_O_2Pst_cm2( 416) =    0.00
  given_sigma_ion_O_4S_cm2( 417) =  2.37 ; given_sigma_ion_O_2D_cm2( 417) =    3.64 ; given_sigma_ion_O_2P_cm2( 417) =    2.28 ; given_sigma_ion_O_4Pst_cm2( 417) =    0.82 ; given_sigma_ion_O_2Pst_cm2( 417) =    0.00
  given_sigma_ion_O_4S_cm2( 418) =  2.38 ; given_sigma_ion_O_2D_cm2( 418) =    3.67 ; given_sigma_ion_O_2P_cm2( 418) =    2.29 ; given_sigma_ion_O_4Pst_cm2( 418) =    0.82 ; given_sigma_ion_O_2Pst_cm2( 418) =    0.00
  given_sigma_ion_O_4S_cm2( 419) =  2.40 ; given_sigma_ion_O_2D_cm2( 419) =    3.70 ; given_sigma_ion_O_2P_cm2( 419) =    2.31 ; given_sigma_ion_O_4Pst_cm2( 419) =    0.83 ; given_sigma_ion_O_2Pst_cm2( 419) =    0.00
  given_sigma_ion_O_4S_cm2( 420) =  2.40 ; given_sigma_ion_O_2D_cm2( 420) =    3.70 ; given_sigma_ion_O_2P_cm2( 420) =    2.31 ; given_sigma_ion_O_4Pst_cm2( 420) =    0.83 ; given_sigma_ion_O_2Pst_cm2( 420) =    0.00
  given_sigma_ion_O_4S_cm2( 421) =  2.44 ; given_sigma_ion_O_2D_cm2( 421) =    3.75 ; given_sigma_ion_O_2P_cm2( 421) =    2.34 ; given_sigma_ion_O_4Pst_cm2( 421) =    0.84 ; given_sigma_ion_O_2Pst_cm2( 421) =    0.00
  given_sigma_ion_O_4S_cm2( 422) =  2.46 ; given_sigma_ion_O_2D_cm2( 422) =    3.78 ; given_sigma_ion_O_2P_cm2( 422) =    2.36 ; given_sigma_ion_O_4Pst_cm2( 422) =    0.85 ; given_sigma_ion_O_2Pst_cm2( 422) =    0.00
  given_sigma_ion_O_4S_cm2( 423) =  2.49 ; given_sigma_ion_O_2D_cm2( 423) =    3.83 ; given_sigma_ion_O_2P_cm2( 423) =    2.39 ; given_sigma_ion_O_4Pst_cm2( 423) =    0.86 ; given_sigma_ion_O_2Pst_cm2( 423) =    0.00
  given_sigma_ion_O_4S_cm2( 424) =  2.49 ; given_sigma_ion_O_2D_cm2( 424) =    3.84 ; given_sigma_ion_O_2P_cm2( 424) =    2.40 ; given_sigma_ion_O_4Pst_cm2( 424) =    0.86 ; given_sigma_ion_O_2Pst_cm2( 424) =    0.00
  given_sigma_ion_O_4S_cm2( 425) =  2.53 ; given_sigma_ion_O_2D_cm2( 425) =    3.89 ; given_sigma_ion_O_2P_cm2( 425) =    2.43 ; given_sigma_ion_O_4Pst_cm2( 425) =    0.87 ; given_sigma_ion_O_2Pst_cm2( 425) =    0.00
  given_sigma_ion_O_4S_cm2( 426) =  2.53 ; given_sigma_ion_O_2D_cm2( 426) =    3.89 ; given_sigma_ion_O_2P_cm2( 426) =    2.43 ; given_sigma_ion_O_4Pst_cm2( 426) =    0.88 ; given_sigma_ion_O_2Pst_cm2( 426) =    0.00
  given_sigma_ion_O_4S_cm2( 427) =  2.56 ; given_sigma_ion_O_2D_cm2( 427) =    3.93 ; given_sigma_ion_O_2P_cm2( 427) =    2.46 ; given_sigma_ion_O_4Pst_cm2( 427) =    0.79 ; given_sigma_ion_O_2Pst_cm2( 427) =    0.00
  given_sigma_ion_O_4S_cm2( 428) =  2.57 ; given_sigma_ion_O_2D_cm2( 428) =    3.97 ; given_sigma_ion_O_2P_cm2( 428) =    2.47 ; given_sigma_ion_O_4Pst_cm2( 428) =    0.79 ; given_sigma_ion_O_2Pst_cm2( 428) =    0.00
  given_sigma_ion_O_4S_cm2( 429) =  2.59 ; given_sigma_ion_O_2D_cm2( 429) =    4.00 ; given_sigma_ion_O_2P_cm2( 429) =    2.49 ; given_sigma_ion_O_4Pst_cm2( 429) =    0.78 ; given_sigma_ion_O_2Pst_cm2( 429) =    0.00
  given_sigma_ion_O_4S_cm2( 430) =  2.65 ; given_sigma_ion_O_2D_cm2( 430) =    4.12 ; given_sigma_ion_O_2P_cm2( 430) =    2.55 ; given_sigma_ion_O_4Pst_cm2( 430) =    0.78 ; given_sigma_ion_O_2Pst_cm2( 430) =    0.00
  given_sigma_ion_O_4S_cm2( 431) =  2.73 ; given_sigma_ion_O_2D_cm2( 431) =    4.27 ; given_sigma_ion_O_2P_cm2( 431) =    2.63 ; given_sigma_ion_O_4Pst_cm2( 431) =    0.77 ; given_sigma_ion_O_2Pst_cm2( 431) =    0.00
  given_sigma_ion_O_4S_cm2( 432) =  3.06 ; given_sigma_ion_O_2D_cm2( 432) =    4.82 ; given_sigma_ion_O_2P_cm2( 432) =    2.94 ; given_sigma_ion_O_4Pst_cm2( 432) =    0.83 ; given_sigma_ion_O_2Pst_cm2( 432) =    0.00
  given_sigma_ion_O_4S_cm2( 433) =  3.07 ; given_sigma_ion_O_2D_cm2( 433) =    4.83 ; given_sigma_ion_O_2P_cm2( 433) =    2.95 ; given_sigma_ion_O_4Pst_cm2( 433) =    0.83 ; given_sigma_ion_O_2Pst_cm2( 433) =    0.00
  given_sigma_ion_O_4S_cm2( 434) =  3.04 ; given_sigma_ion_O_2D_cm2( 434) =    4.80 ; given_sigma_ion_O_2P_cm2( 434) =    2.93 ; given_sigma_ion_O_4Pst_cm2( 434) =    0.82 ; given_sigma_ion_O_2Pst_cm2( 434) =    0.00
  given_sigma_ion_O_4S_cm2( 435) =  3.03 ; given_sigma_ion_O_2D_cm2( 435) =    4.78 ; given_sigma_ion_O_2P_cm2( 435) =    2.92 ; given_sigma_ion_O_4Pst_cm2( 435) =    0.82 ; given_sigma_ion_O_2Pst_cm2( 435) =    0.00
  given_sigma_ion_O_4S_cm2( 436) =  3.03 ; given_sigma_ion_O_2D_cm2( 436) =    4.77 ; given_sigma_ion_O_2P_cm2( 436) =    2.91 ; given_sigma_ion_O_4Pst_cm2( 436) =    0.82 ; given_sigma_ion_O_2Pst_cm2( 436) =    0.00
  given_sigma_ion_O_4S_cm2( 437) =  3.00 ; given_sigma_ion_O_2D_cm2( 437) =    4.73 ; given_sigma_ion_O_2P_cm2( 437) =    2.88 ; given_sigma_ion_O_4Pst_cm2( 437) =    0.81 ; given_sigma_ion_O_2Pst_cm2( 437) =    0.00
  given_sigma_ion_O_4S_cm2( 438) =  2.96 ; given_sigma_ion_O_2D_cm2( 438) =    4.67 ; given_sigma_ion_O_2P_cm2( 438) =    2.85 ; given_sigma_ion_O_4Pst_cm2( 438) =    0.80 ; given_sigma_ion_O_2Pst_cm2( 438) =    0.00
  given_sigma_ion_O_4S_cm2( 439) =  2.94 ; given_sigma_ion_O_2D_cm2( 439) =    4.64 ; given_sigma_ion_O_2P_cm2( 439) =    2.94 ; given_sigma_ion_O_4Pst_cm2( 439) =    0.79 ; given_sigma_ion_O_2Pst_cm2( 439) =    0.00
  given_sigma_ion_O_4S_cm2( 440) =  2.92 ; given_sigma_ion_O_2D_cm2( 440) =    4.61 ; given_sigma_ion_O_2P_cm2( 440) =    2.92 ; given_sigma_ion_O_4Pst_cm2( 440) =    0.79 ; given_sigma_ion_O_2Pst_cm2( 440) =    0.00
  given_sigma_ion_O_4S_cm2( 441) =  2.90 ; given_sigma_ion_O_2D_cm2( 441) =    4.58 ; given_sigma_ion_O_2P_cm2( 441) =    2.90 ; given_sigma_ion_O_4Pst_cm2( 441) =    0.78 ; given_sigma_ion_O_2Pst_cm2( 441) =    0.00
  given_sigma_ion_O_4S_cm2( 442) =  2.88 ; given_sigma_ion_O_2D_cm2( 442) =    4.54 ; given_sigma_ion_O_2P_cm2( 442) =    2.88 ; given_sigma_ion_O_4Pst_cm2( 442) =    0.78 ; given_sigma_ion_O_2Pst_cm2( 442) =    0.00
  given_sigma_ion_O_4S_cm2( 443) =  2.86 ; given_sigma_ion_O_2D_cm2( 443) =    4.51 ; given_sigma_ion_O_2P_cm2( 443) =    2.86 ; given_sigma_ion_O_4Pst_cm2( 443) =    0.77 ; given_sigma_ion_O_2Pst_cm2( 443) =    0.00
  given_sigma_ion_O_4S_cm2( 444) =  2.87 ; given_sigma_ion_O_2D_cm2( 444) =    4.41 ; given_sigma_ion_O_2P_cm2( 444) =    2.87 ; given_sigma_ion_O_4Pst_cm2( 444) =    0.77 ; given_sigma_ion_O_2Pst_cm2( 444) =    0.00
  given_sigma_ion_O_4S_cm2( 445) =  2.87 ; given_sigma_ion_O_2D_cm2( 445) =    4.53 ; given_sigma_ion_O_2P_cm2( 445) =    2.87 ; given_sigma_ion_O_4Pst_cm2( 445) =    0.77 ; given_sigma_ion_O_2Pst_cm2( 445) =    0.00
  given_sigma_ion_O_4S_cm2( 446) =  2.88 ; given_sigma_ion_O_2D_cm2( 446) =    4.54 ; given_sigma_ion_O_2P_cm2( 446) =    2.88 ; given_sigma_ion_O_4Pst_cm2( 446) =    0.78 ; given_sigma_ion_O_2Pst_cm2( 446) =    0.00
  given_sigma_ion_O_4S_cm2( 447) =  2.89 ; given_sigma_ion_O_2D_cm2( 447) =    4.56 ; given_sigma_ion_O_2P_cm2( 447) =    2.89 ; given_sigma_ion_O_4Pst_cm2( 447) =    0.78 ; given_sigma_ion_O_2Pst_cm2( 447) =    0.00
  given_sigma_ion_O_4S_cm2( 448) =  2.90 ; given_sigma_ion_O_2D_cm2( 448) =    4.57 ; given_sigma_ion_O_2P_cm2( 448) =    2.90 ; given_sigma_ion_O_4Pst_cm2( 448) =    0.78 ; given_sigma_ion_O_2Pst_cm2( 448) =    0.00
  given_sigma_ion_O_4S_cm2( 449) =  2.90 ; given_sigma_ion_O_2D_cm2( 449) =    4.58 ; given_sigma_ion_O_2P_cm2( 449) =    2.90 ; given_sigma_ion_O_4Pst_cm2( 449) =    0.78 ; given_sigma_ion_O_2Pst_cm2( 449) =    0.00
  given_sigma_ion_O_4S_cm2( 450) =  2.91 ; given_sigma_ion_O_2D_cm2( 450) =    4.59 ; given_sigma_ion_O_2P_cm2( 450) =    2.91 ; given_sigma_ion_O_4Pst_cm2( 450) =    0.78 ; given_sigma_ion_O_2Pst_cm2( 450) =    0.00
  given_sigma_ion_O_4S_cm2( 451) =  2.91 ; given_sigma_ion_O_2D_cm2( 451) =    4.60 ; given_sigma_ion_O_2P_cm2( 451) =    2.91 ; given_sigma_ion_O_4Pst_cm2( 451) =    0.78 ; given_sigma_ion_O_2Pst_cm2( 451) =    0.00
  given_sigma_ion_O_4S_cm2( 452) =  2.92 ; given_sigma_ion_O_2D_cm2( 452) =    4.60 ; given_sigma_ion_O_2P_cm2( 452) =    2.92 ; given_sigma_ion_O_4Pst_cm2( 452) =    0.79 ; given_sigma_ion_O_2Pst_cm2( 452) =    0.00
  given_sigma_ion_O_4S_cm2( 453) =  2.92 ; given_sigma_ion_O_2D_cm2( 453) =    4.61 ; given_sigma_ion_O_2P_cm2( 453) =    2.92 ; given_sigma_ion_O_4Pst_cm2( 453) =    0.79 ; given_sigma_ion_O_2Pst_cm2( 453) =    0.00
  given_sigma_ion_O_4S_cm2( 454) =  2.93 ; given_sigma_ion_O_2D_cm2( 454) =    4.73 ; given_sigma_ion_O_2P_cm2( 454) =    2.93 ; given_sigma_ion_O_4Pst_cm2( 454) =    0.79 ; given_sigma_ion_O_2Pst_cm2( 454) =    0.00
  given_sigma_ion_O_4S_cm2( 455) =  2.94 ; given_sigma_ion_O_2D_cm2( 455) =    4.63 ; given_sigma_ion_O_2P_cm2( 455) =    2.94 ; given_sigma_ion_O_4Pst_cm2( 455) =    0.79 ; given_sigma_ion_O_2Pst_cm2( 455) =    0.00
  given_sigma_ion_O_4S_cm2( 456) =  2.94 ; given_sigma_ion_O_2D_cm2( 456) =    4.64 ; given_sigma_ion_O_2P_cm2( 456) =    2.94 ; given_sigma_ion_O_4Pst_cm2( 456) =    0.79 ; given_sigma_ion_O_2Pst_cm2( 456) =    0.00
  given_sigma_ion_O_4S_cm2( 457) =  2.95 ; given_sigma_ion_O_2D_cm2( 457) =    4.76 ; given_sigma_ion_O_2P_cm2( 457) =    2.95 ; given_sigma_ion_O_4Pst_cm2( 457) =    0.79 ; given_sigma_ion_O_2Pst_cm2( 457) =    0.00
  given_sigma_ion_O_4S_cm2( 458) =  2.95 ; given_sigma_ion_O_2D_cm2( 458) =    4.77 ; given_sigma_ion_O_2P_cm2( 458) =    2.95 ; given_sigma_ion_O_4Pst_cm2( 458) =    0.80 ; given_sigma_ion_O_2Pst_cm2( 458) =    0.00
  given_sigma_ion_O_4S_cm2( 459) =  2.96 ; given_sigma_ion_O_2D_cm2( 459) =    4.78 ; given_sigma_ion_O_2P_cm2( 459) =    2.96 ; given_sigma_ion_O_4Pst_cm2( 459) =    0.80 ; given_sigma_ion_O_2Pst_cm2( 459) =    0.00
  given_sigma_ion_O_4S_cm2( 460) =  2.96 ; given_sigma_ion_O_2D_cm2( 460) =    4.79 ; given_sigma_ion_O_2P_cm2( 460) =    2.96 ; given_sigma_ion_O_4Pst_cm2( 460) =    0.80 ; given_sigma_ion_O_2Pst_cm2( 460) =    0.00
  given_sigma_ion_O_4S_cm2( 461) =  2.97 ; given_sigma_ion_O_2D_cm2( 461) =    4.80 ; given_sigma_ion_O_2P_cm2( 461) =    2.97 ; given_sigma_ion_O_4Pst_cm2( 461) =    0.80 ; given_sigma_ion_O_2Pst_cm2( 461) =    0.00
  given_sigma_ion_O_4S_cm2( 462) =  2.97 ; given_sigma_ion_O_2D_cm2( 462) =    4.80 ; given_sigma_ion_O_2P_cm2( 462) =    2.97 ; given_sigma_ion_O_4Pst_cm2( 462) =    0.80 ; given_sigma_ion_O_2Pst_cm2( 462) =    0.00
  given_sigma_ion_O_4S_cm2( 463) =  2.98 ; given_sigma_ion_O_2D_cm2( 463) =    4.81 ; given_sigma_ion_O_2P_cm2( 463) =    2.98 ; given_sigma_ion_O_4Pst_cm2( 463) =    0.69 ; given_sigma_ion_O_2Pst_cm2( 463) =    0.00
  given_sigma_ion_O_4S_cm2( 464) =  2.98 ; given_sigma_ion_O_2D_cm2( 464) =    4.82 ; given_sigma_ion_O_2P_cm2( 464) =    2.98 ; given_sigma_ion_O_4Pst_cm2( 464) =    0.69 ; given_sigma_ion_O_2Pst_cm2( 464) =    0.00
  given_sigma_ion_O_4S_cm2( 465) =  2.99 ; given_sigma_ion_O_2D_cm2( 465) =    4.83 ; given_sigma_ion_O_2P_cm2( 465) =    2.99 ; given_sigma_ion_O_4Pst_cm2( 465) =    0.69 ; given_sigma_ion_O_2Pst_cm2( 465) =    0.00
  given_sigma_ion_O_4S_cm2( 466) =  2.99 ; given_sigma_ion_O_2D_cm2( 466) =    4.83 ; given_sigma_ion_O_2P_cm2( 466) =    2.99 ; given_sigma_ion_O_4Pst_cm2( 466) =    0.69 ; given_sigma_ion_O_2Pst_cm2( 466) =    0.00
  given_sigma_ion_O_4S_cm2( 467) =  2.99 ; given_sigma_ion_O_2D_cm2( 467) =    4.82 ; given_sigma_ion_O_2P_cm2( 467) =    2.99 ; given_sigma_ion_O_4Pst_cm2( 467) =    0.69 ; given_sigma_ion_O_2Pst_cm2( 467) =    0.00
  given_sigma_ion_O_4S_cm2( 468) =  2.98 ; given_sigma_ion_O_2D_cm2( 468) =    4.81 ; given_sigma_ion_O_2P_cm2( 468) =    2.98 ; given_sigma_ion_O_4Pst_cm2( 468) =    0.69 ; given_sigma_ion_O_2Pst_cm2( 468) =    0.00
  given_sigma_ion_O_4S_cm2( 469) =  2.98 ; given_sigma_ion_O_2D_cm2( 469) =    4.81 ; given_sigma_ion_O_2P_cm2( 469) =    2.98 ; given_sigma_ion_O_4Pst_cm2( 469) =    0.69 ; given_sigma_ion_O_2Pst_cm2( 469) =    0.00
  given_sigma_ion_O_4S_cm2( 470) =  2.97 ; given_sigma_ion_O_2D_cm2( 470) =    4.80 ; given_sigma_ion_O_2P_cm2( 470) =    2.97 ; given_sigma_ion_O_4Pst_cm2( 470) =    0.69 ; given_sigma_ion_O_2Pst_cm2( 470) =    0.00
  given_sigma_ion_O_4S_cm2( 471) =  2.97 ; given_sigma_ion_O_2D_cm2( 471) =    4.79 ; given_sigma_ion_O_2P_cm2( 471) =    2.97 ; given_sigma_ion_O_4Pst_cm2( 471) =    0.68 ; given_sigma_ion_O_2Pst_cm2( 471) =    0.00
  given_sigma_ion_O_4S_cm2( 472) =  3.19 ; given_sigma_ion_O_2D_cm2( 472) =    5.13 ; given_sigma_ion_O_2P_cm2( 472) =    3.08 ; given_sigma_ion_O_4Pst_cm2( 472) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 472) =    0.00
  given_sigma_ion_O_4S_cm2( 473) =  3.19 ; given_sigma_ion_O_2D_cm2( 473) =    5.13 ; given_sigma_ion_O_2P_cm2( 473) =    3.08 ; given_sigma_ion_O_4Pst_cm2( 473) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 473) =    0.00
  given_sigma_ion_O_4S_cm2( 474) =  3.19 ; given_sigma_ion_O_2D_cm2( 474) =    5.12 ; given_sigma_ion_O_2P_cm2( 474) =    3.07 ; given_sigma_ion_O_4Pst_cm2( 474) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 474) =    0.00
  given_sigma_ion_O_4S_cm2( 475) =  3.20 ; given_sigma_ion_O_2D_cm2( 475) =    5.15 ; given_sigma_ion_O_2P_cm2( 475) =    3.09 ; given_sigma_ion_O_4Pst_cm2( 475) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 475) =    0.00
  given_sigma_ion_O_4S_cm2( 476) =  3.25 ; given_sigma_ion_O_2D_cm2( 476) =    5.22 ; given_sigma_ion_O_2P_cm2( 476) =    3.13 ; given_sigma_ion_O_4Pst_cm2( 476) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 476) =    0.00
  given_sigma_ion_O_4S_cm2( 477) =  3.29 ; given_sigma_ion_O_2D_cm2( 477) =    5.29 ; given_sigma_ion_O_2P_cm2( 477) =    3.17 ; given_sigma_ion_O_4Pst_cm2( 477) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 477) =    0.00
  given_sigma_ion_O_4S_cm2( 478) =  3.34 ; given_sigma_ion_O_2D_cm2( 478) =    5.36 ; given_sigma_ion_O_2P_cm2( 478) =    3.22 ; given_sigma_ion_O_4Pst_cm2( 478) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 478) =    0.00
  given_sigma_ion_O_4S_cm2( 479) =  3.36 ; given_sigma_ion_O_2D_cm2( 479) =    5.40 ; given_sigma_ion_O_2P_cm2( 479) =    3.24 ; given_sigma_ion_O_4Pst_cm2( 479) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 479) =    0.00
  given_sigma_ion_O_4S_cm2( 480) =  3.27 ; given_sigma_ion_O_2D_cm2( 480) =    5.26 ; given_sigma_ion_O_2P_cm2( 480) =    3.16 ; given_sigma_ion_O_4Pst_cm2( 480) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 480) =    0.00
  given_sigma_ion_O_4S_cm2( 481) =  3.10 ; given_sigma_ion_O_2D_cm2( 481) =    4.99 ; given_sigma_ion_O_2P_cm2( 481) =    2.99 ; given_sigma_ion_O_4Pst_cm2( 481) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 481) =    0.00
  given_sigma_ion_O_4S_cm2( 482) =  3.05 ; given_sigma_ion_O_2D_cm2( 482) =    4.90 ; given_sigma_ion_O_2P_cm2( 482) =    2.94 ; given_sigma_ion_O_4Pst_cm2( 482) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 482) =    0.00
  given_sigma_ion_O_4S_cm2( 483) =  3.61 ; given_sigma_ion_O_2D_cm2( 483) =    5.80 ; given_sigma_ion_O_2P_cm2( 483) =    3.48 ; given_sigma_ion_O_4Pst_cm2( 483) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 483) =    0.00
  given_sigma_ion_O_4S_cm2( 484) =  3.49 ; given_sigma_ion_O_2D_cm2( 484) =    5.61 ; given_sigma_ion_O_2P_cm2( 484) =    3.37 ; given_sigma_ion_O_4Pst_cm2( 484) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 484) =    0.00
  given_sigma_ion_O_4S_cm2( 485) =  3.25 ; given_sigma_ion_O_2D_cm2( 485) =    5.22 ; given_sigma_ion_O_2P_cm2( 485) =    3.13 ; given_sigma_ion_O_4Pst_cm2( 485) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 485) =    0.00
  given_sigma_ion_O_4S_cm2( 486) =  3.21 ; given_sigma_ion_O_2D_cm2( 486) =    5.16 ; given_sigma_ion_O_2P_cm2( 486) =    3.10 ; given_sigma_ion_O_4Pst_cm2( 486) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 486) =    0.00
  given_sigma_ion_O_4S_cm2( 487) =  3.12 ; given_sigma_ion_O_2D_cm2( 487) =    5.01 ; given_sigma_ion_O_2P_cm2( 487) =    3.01 ; given_sigma_ion_O_4Pst_cm2( 487) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 487) =    0.00
  given_sigma_ion_O_4S_cm2( 488) =  3.11 ; given_sigma_ion_O_2D_cm2( 488) =    4.99 ; given_sigma_ion_O_2P_cm2( 488) =    3.00 ; given_sigma_ion_O_4Pst_cm2( 488) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 488) =    0.00
  given_sigma_ion_O_4S_cm2( 489) =  2.97 ; given_sigma_ion_O_2D_cm2( 489) =    4.77 ; given_sigma_ion_O_2P_cm2( 489) =    2.86 ; given_sigma_ion_O_4Pst_cm2( 489) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 489) =    0.00
  given_sigma_ion_O_4S_cm2( 490) =  2.94 ; given_sigma_ion_O_2D_cm2( 490) =    4.72 ; given_sigma_ion_O_2P_cm2( 490) =    2.84 ; given_sigma_ion_O_4Pst_cm2( 490) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 490) =    0.00
  given_sigma_ion_O_4S_cm2( 491) =  3.92 ; given_sigma_ion_O_2D_cm2( 491) =    6.30 ; given_sigma_ion_O_2P_cm2( 491) =    3.78 ; given_sigma_ion_O_4Pst_cm2( 491) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 491) =    0.00
  given_sigma_ion_O_4S_cm2( 492) =  3.74 ; given_sigma_ion_O_2D_cm2( 492) =    6.02 ; given_sigma_ion_O_2P_cm2( 492) =    3.61 ; given_sigma_ion_O_4Pst_cm2( 492) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 492) =    0.00
  given_sigma_ion_O_4S_cm2( 493) =  3.39 ; given_sigma_ion_O_2D_cm2( 493) =    5.45 ; given_sigma_ion_O_2P_cm2( 493) =    3.25 ; given_sigma_ion_O_4Pst_cm2( 493) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 493) =    0.00
  given_sigma_ion_O_4S_cm2( 494) =  3.28 ; given_sigma_ion_O_2D_cm2( 494) =    5.27 ; given_sigma_ion_O_2P_cm2( 494) =    3.16 ; given_sigma_ion_O_4Pst_cm2( 494) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 494) =    0.00
  given_sigma_ion_O_4S_cm2( 495) =  3.11 ; given_sigma_ion_O_2D_cm2( 495) =    4.99 ; given_sigma_ion_O_2P_cm2( 495) =    3.00 ; given_sigma_ion_O_4Pst_cm2( 495) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 495) =    0.00
  given_sigma_ion_O_4S_cm2( 496) =  3.12 ; given_sigma_ion_O_2D_cm2( 496) =    5.02 ; given_sigma_ion_O_2P_cm2( 496) =    3.01 ; given_sigma_ion_O_4Pst_cm2( 496) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 496) =    0.00
  given_sigma_ion_O_4S_cm2( 497) =  3.16 ; given_sigma_ion_O_2D_cm2( 497) =    5.07 ; given_sigma_ion_O_2P_cm2( 497) =    3.04 ; given_sigma_ion_O_4Pst_cm2( 497) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 497) =    0.00
  given_sigma_ion_O_4S_cm2( 498) =  3.16 ; given_sigma_ion_O_2D_cm2( 498) =    5.09 ; given_sigma_ion_O_2P_cm2( 498) =    3.05 ; given_sigma_ion_O_4Pst_cm2( 498) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 498) =    0.00
  given_sigma_ion_O_4S_cm2( 499) =  3.14 ; given_sigma_ion_O_2D_cm2( 499) =    5.04 ; given_sigma_ion_O_2P_cm2( 499) =    3.02 ; given_sigma_ion_O_4Pst_cm2( 499) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 499) =    0.00
  given_sigma_ion_O_4S_cm2( 500) =  3.10 ; given_sigma_ion_O_2D_cm2( 500) =    4.98 ; given_sigma_ion_O_2P_cm2( 500) =    2.99 ; given_sigma_ion_O_4Pst_cm2( 500) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 500) =    0.00
  given_sigma_ion_O_4S_cm2( 501) =  3.08 ; given_sigma_ion_O_2D_cm2( 501) =    4.95 ; given_sigma_ion_O_2P_cm2( 501) =    2.97 ; given_sigma_ion_O_4Pst_cm2( 501) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 501) =    0.00
  given_sigma_ion_O_4S_cm2( 502) =  3.01 ; given_sigma_ion_O_2D_cm2( 502) =    4.84 ; given_sigma_ion_O_2P_cm2( 502) =    2.91 ; given_sigma_ion_O_4Pst_cm2( 502) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 502) =    0.00
  given_sigma_ion_O_4S_cm2( 503) =  2.85 ; given_sigma_ion_O_2D_cm2( 503) =    4.58 ; given_sigma_ion_O_2P_cm2( 503) =    2.75 ; given_sigma_ion_O_4Pst_cm2( 503) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 503) =    0.00
  given_sigma_ion_O_4S_cm2( 504) =  2.80 ; given_sigma_ion_O_2D_cm2( 504) =    4.47 ; given_sigma_ion_O_2P_cm2( 504) =    2.70 ; given_sigma_ion_O_4Pst_cm2( 504) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 504) =    0.00
  given_sigma_ion_O_4S_cm2( 505) =  3.08 ; given_sigma_ion_O_2D_cm2( 505) =    4.87 ; given_sigma_ion_O_2P_cm2( 505) =    2.97 ; given_sigma_ion_O_4Pst_cm2( 505) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 505) =    0.00
  given_sigma_ion_O_4S_cm2( 506) =  3.57 ; given_sigma_ion_O_2D_cm2( 506) =    5.61 ; given_sigma_ion_O_2P_cm2( 506) =    3.44 ; given_sigma_ion_O_4Pst_cm2( 506) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 506) =    0.00
  given_sigma_ion_O_4S_cm2( 507) =  4.23 ; given_sigma_ion_O_2D_cm2( 507) =    6.70 ; given_sigma_ion_O_2P_cm2( 507) =    4.08 ; given_sigma_ion_O_4Pst_cm2( 507) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 507) =    0.00
  given_sigma_ion_O_4S_cm2( 508) =  3.88 ; given_sigma_ion_O_2D_cm2( 508) =    6.23 ; given_sigma_ion_O_2P_cm2( 508) =    3.74 ; given_sigma_ion_O_4Pst_cm2( 508) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 508) =    0.00
  given_sigma_ion_O_4S_cm2( 509) =  3.52 ; given_sigma_ion_O_2D_cm2( 509) =    5.76 ; given_sigma_ion_O_2P_cm2( 509) =    3.46 ; given_sigma_ion_O_4Pst_cm2( 509) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 509) =    0.00
  given_sigma_ion_O_4S_cm2( 510) =  3.41 ; given_sigma_ion_O_2D_cm2( 510) =    5.68 ; given_sigma_ion_O_2P_cm2( 510) =    3.41 ; given_sigma_ion_O_4Pst_cm2( 510) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 510) =    0.00
  given_sigma_ion_O_4S_cm2( 511) =  3.31 ; given_sigma_ion_O_2D_cm2( 511) =    5.51 ; given_sigma_ion_O_2P_cm2( 511) =    3.31 ; given_sigma_ion_O_4Pst_cm2( 511) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 511) =    0.00
  given_sigma_ion_O_4S_cm2( 512) =  3.21 ; given_sigma_ion_O_2D_cm2( 512) =    5.34 ; given_sigma_ion_O_2P_cm2( 512) =    3.21 ; given_sigma_ion_O_4Pst_cm2( 512) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 512) =    0.00
  given_sigma_ion_O_4S_cm2( 513) =  3.21 ; given_sigma_ion_O_2D_cm2( 513) =    5.31 ; given_sigma_ion_O_2P_cm2( 513) =    3.19 ; given_sigma_ion_O_4Pst_cm2( 513) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 513) =    0.00
  given_sigma_ion_O_4S_cm2( 514) =  3.31 ; given_sigma_ion_O_2D_cm2( 514) =    5.32 ; given_sigma_ion_O_2P_cm2( 514) =    3.19 ; given_sigma_ion_O_4Pst_cm2( 514) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 514) =    0.00
  given_sigma_ion_O_4S_cm2( 515) =  3.31 ; given_sigma_ion_O_2D_cm2( 515) =    5.33 ; given_sigma_ion_O_2P_cm2( 515) =    3.20 ; given_sigma_ion_O_4Pst_cm2( 515) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 515) =    0.00
  given_sigma_ion_O_4S_cm2( 516) =  3.32 ; given_sigma_ion_O_2D_cm2( 516) =    5.34 ; given_sigma_ion_O_2P_cm2( 516) =    3.20 ; given_sigma_ion_O_4Pst_cm2( 516) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 516) =    0.00
  given_sigma_ion_O_4S_cm2( 517) =  3.33 ; given_sigma_ion_O_2D_cm2( 517) =    5.35 ; given_sigma_ion_O_2P_cm2( 517) =    3.21 ; given_sigma_ion_O_4Pst_cm2( 517) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 517) =    0.00
  given_sigma_ion_O_4S_cm2( 518) =  3.33 ; given_sigma_ion_O_2D_cm2( 518) =    5.36 ; given_sigma_ion_O_2P_cm2( 518) =    3.21 ; given_sigma_ion_O_4Pst_cm2( 518) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 518) =    0.00
  given_sigma_ion_O_4S_cm2( 519) =  3.34 ; given_sigma_ion_O_2D_cm2( 519) =    5.37 ; given_sigma_ion_O_2P_cm2( 519) =    3.22 ; given_sigma_ion_O_4Pst_cm2( 519) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 519) =    0.00
  given_sigma_ion_O_4S_cm2( 520) =  3.34 ; given_sigma_ion_O_2D_cm2( 520) =    5.37 ; given_sigma_ion_O_2P_cm2( 520) =    3.22 ; given_sigma_ion_O_4Pst_cm2( 520) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 520) =    0.00
  given_sigma_ion_O_4S_cm2( 521) =  3.34 ; given_sigma_ion_O_2D_cm2( 521) =    5.38 ; given_sigma_ion_O_2P_cm2( 521) =    3.23 ; given_sigma_ion_O_4Pst_cm2( 521) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 521) =    0.00
  given_sigma_ion_O_4S_cm2( 522) =  3.35 ; given_sigma_ion_O_2D_cm2( 522) =    5.39 ; given_sigma_ion_O_2P_cm2( 522) =    3.23 ; given_sigma_ion_O_4Pst_cm2( 522) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 522) =    0.00
  given_sigma_ion_O_4S_cm2( 523) =  3.36 ; given_sigma_ion_O_2D_cm2( 523) =    5.39 ; given_sigma_ion_O_2P_cm2( 523) =    3.24 ; given_sigma_ion_O_4Pst_cm2( 523) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 523) =    0.00
  given_sigma_ion_O_4S_cm2( 524) =  3.36 ; given_sigma_ion_O_2D_cm2( 524) =    5.40 ; given_sigma_ion_O_2P_cm2( 524) =    3.24 ; given_sigma_ion_O_4Pst_cm2( 524) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 524) =    0.00
  given_sigma_ion_O_4S_cm2( 525) =  3.37 ; given_sigma_ion_O_2D_cm2( 525) =    5.41 ; given_sigma_ion_O_2P_cm2( 525) =    3.25 ; given_sigma_ion_O_4Pst_cm2( 525) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 525) =    0.00
  given_sigma_ion_O_4S_cm2( 526) =  3.37 ; given_sigma_ion_O_2D_cm2( 526) =    5.41 ; given_sigma_ion_O_2P_cm2( 526) =    3.25 ; given_sigma_ion_O_4Pst_cm2( 526) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 526) =    0.00
  given_sigma_ion_O_4S_cm2( 527) =  3.37 ; given_sigma_ion_O_2D_cm2( 527) =    5.42 ; given_sigma_ion_O_2P_cm2( 527) =    3.25 ; given_sigma_ion_O_4Pst_cm2( 527) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 527) =    0.00
  given_sigma_ion_O_4S_cm2( 528) =  3.38 ; given_sigma_ion_O_2D_cm2( 528) =    5.43 ; given_sigma_ion_O_2P_cm2( 528) =    3.26 ; given_sigma_ion_O_4Pst_cm2( 528) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 528) =    0.00
  given_sigma_ion_O_4S_cm2( 529) =  3.39 ; given_sigma_ion_O_2D_cm2( 529) =    5.44 ; given_sigma_ion_O_2P_cm2( 529) =    3.27 ; given_sigma_ion_O_4Pst_cm2( 529) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 529) =    0.00
  given_sigma_ion_O_4S_cm2( 530) =  3.39 ; given_sigma_ion_O_2D_cm2( 530) =    5.45 ; given_sigma_ion_O_2P_cm2( 530) =    3.27 ; given_sigma_ion_O_4Pst_cm2( 530) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 530) =    0.00
  given_sigma_ion_O_4S_cm2( 531) =  3.36 ; given_sigma_ion_O_2D_cm2( 531) =    5.40 ; given_sigma_ion_O_2P_cm2( 531) =    3.24 ; given_sigma_ion_O_4Pst_cm2( 531) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 531) =    0.00
  given_sigma_ion_O_4S_cm2( 532) =  3.35 ; given_sigma_ion_O_2D_cm2( 532) =    5.38 ; given_sigma_ion_O_2P_cm2( 532) =    3.23 ; given_sigma_ion_O_4Pst_cm2( 532) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 532) =    0.00
  given_sigma_ion_O_4S_cm2( 533) =  3.30 ; given_sigma_ion_O_2D_cm2( 533) =    5.31 ; given_sigma_ion_O_2P_cm2( 533) =    3.19 ; given_sigma_ion_O_4Pst_cm2( 533) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 533) =    0.00
  given_sigma_ion_O_4S_cm2( 534) =  3.22 ; given_sigma_ion_O_2D_cm2( 534) =    5.17 ; given_sigma_ion_O_2P_cm2( 534) =    3.11 ; given_sigma_ion_O_4Pst_cm2( 534) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 534) =    0.00
  given_sigma_ion_O_4S_cm2( 535) =  3.21 ; given_sigma_ion_O_2D_cm2( 535) =    5.16 ; given_sigma_ion_O_2P_cm2( 535) =    3.09 ; given_sigma_ion_O_4Pst_cm2( 535) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 535) =    0.00
  given_sigma_ion_O_4S_cm2( 536) =  3.16 ; given_sigma_ion_O_2D_cm2( 536) =    5.09 ; given_sigma_ion_O_2P_cm2( 536) =    3.05 ; given_sigma_ion_O_4Pst_cm2( 536) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 536) =    0.00
  given_sigma_ion_O_4S_cm2( 537) =  3.25 ; given_sigma_ion_O_2D_cm2( 537) =    5.23 ; given_sigma_ion_O_2P_cm2( 537) =    3.14 ; given_sigma_ion_O_4Pst_cm2( 537) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 537) =    0.00
  given_sigma_ion_O_4S_cm2( 538) =  3.40 ; given_sigma_ion_O_2D_cm2( 538) =    5.46 ; given_sigma_ion_O_2P_cm2( 538) =    3.28 ; given_sigma_ion_O_4Pst_cm2( 538) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 538) =    0.00
  given_sigma_ion_O_4S_cm2( 539) =  3.47 ; given_sigma_ion_O_2D_cm2( 539) =    5.58 ; given_sigma_ion_O_2P_cm2( 539) =    3.35 ; given_sigma_ion_O_4Pst_cm2( 539) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 539) =    0.00
  given_sigma_ion_O_4S_cm2( 540) =  3.45 ; given_sigma_ion_O_2D_cm2( 540) =    5.55 ; given_sigma_ion_O_2P_cm2( 540) =    3.33 ; given_sigma_ion_O_4Pst_cm2( 540) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 540) =    0.00
  given_sigma_ion_O_4S_cm2( 541) =  3.42 ; given_sigma_ion_O_2D_cm2( 541) =    5.49 ; given_sigma_ion_O_2P_cm2( 541) =    3.29 ; given_sigma_ion_O_4Pst_cm2( 541) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 541) =    0.00
  given_sigma_ion_O_4S_cm2( 542) =  4.51 ; given_sigma_ion_O_2D_cm2( 542) =    7.24 ; given_sigma_ion_O_2P_cm2( 542) =    4.35 ; given_sigma_ion_O_4Pst_cm2( 542) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 542) =    0.00
  given_sigma_ion_O_4S_cm2( 543) =  4.36 ; given_sigma_ion_O_2D_cm2( 543) =    7.00 ; given_sigma_ion_O_2P_cm2( 543) =    4.20 ; given_sigma_ion_O_4Pst_cm2( 543) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 543) =    0.00
  given_sigma_ion_O_4S_cm2( 544) =  3.84 ; given_sigma_ion_O_2D_cm2( 544) =    6.16 ; given_sigma_ion_O_2P_cm2( 544) =    3.70 ; given_sigma_ion_O_4Pst_cm2( 544) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 544) =    0.00
  given_sigma_ion_O_4S_cm2( 545) =  3.78 ; given_sigma_ion_O_2D_cm2( 545) =    6.08 ; given_sigma_ion_O_2P_cm2( 545) =    3.65 ; given_sigma_ion_O_4Pst_cm2( 545) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 545) =    0.00
  given_sigma_ion_O_4S_cm2( 546) =  3.76 ; given_sigma_ion_O_2D_cm2( 546) =    6.05 ; given_sigma_ion_O_2P_cm2( 546) =    3.63 ; given_sigma_ion_O_4Pst_cm2( 546) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 546) =    0.00
  given_sigma_ion_O_4S_cm2( 547) =  3.67 ; given_sigma_ion_O_2D_cm2( 547) =    5.89 ; given_sigma_ion_O_2P_cm2( 547) =    3.54 ; given_sigma_ion_O_4Pst_cm2( 547) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 547) =    0.00
  given_sigma_ion_O_4S_cm2( 548) =  3.65 ; given_sigma_ion_O_2D_cm2( 548) =    5.87 ; given_sigma_ion_O_2P_cm2( 548) =    3.52 ; given_sigma_ion_O_4Pst_cm2( 548) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 548) =    0.00
  given_sigma_ion_O_4S_cm2( 549) =  3.64 ; given_sigma_ion_O_2D_cm2( 549) =    5.85 ; given_sigma_ion_O_2P_cm2( 549) =    3.51 ; given_sigma_ion_O_4Pst_cm2( 549) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 549) =    0.00
  given_sigma_ion_O_4S_cm2( 550) =  3.63 ; given_sigma_ion_O_2D_cm2( 550) =    5.83 ; given_sigma_ion_O_2P_cm2( 550) =    3.50 ; given_sigma_ion_O_4Pst_cm2( 550) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 550) =    0.00
  given_sigma_ion_O_4S_cm2( 551) =  3.61 ; given_sigma_ion_O_2D_cm2( 551) =    5.80 ; given_sigma_ion_O_2P_cm2( 551) =    3.48 ; given_sigma_ion_O_4Pst_cm2( 551) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 551) =    0.00
  given_sigma_ion_O_4S_cm2( 552) =  3.60 ; given_sigma_ion_O_2D_cm2( 552) =    5.79 ; given_sigma_ion_O_2P_cm2( 552) =    3.47 ; given_sigma_ion_O_4Pst_cm2( 552) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 552) =    0.00
  given_sigma_ion_O_4S_cm2( 553) =  3.56 ; given_sigma_ion_O_2D_cm2( 553) =    5.73 ; given_sigma_ion_O_2P_cm2( 553) =    3.44 ; given_sigma_ion_O_4Pst_cm2( 553) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 553) =    0.00
  given_sigma_ion_O_4S_cm2( 554) =  3.54 ; given_sigma_ion_O_2D_cm2( 554) =    5.69 ; given_sigma_ion_O_2P_cm2( 554) =    3.42 ; given_sigma_ion_O_4Pst_cm2( 554) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 554) =    0.00
  given_sigma_ion_O_4S_cm2( 555) =  3.52 ; given_sigma_ion_O_2D_cm2( 555) =    5.66 ; given_sigma_ion_O_2P_cm2( 555) =    3.39 ; given_sigma_ion_O_4Pst_cm2( 555) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 555) =    0.00
  given_sigma_ion_O_4S_cm2( 556) =  3.45 ; given_sigma_ion_O_2D_cm2( 556) =    5.57 ; given_sigma_ion_O_2P_cm2( 556) =    3.34 ; given_sigma_ion_O_4Pst_cm2( 556) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 556) =    0.00
  given_sigma_ion_O_4S_cm2( 557) =  3.41 ; given_sigma_ion_O_2D_cm2( 557) =    5.49 ; given_sigma_ion_O_2P_cm2( 557) =    3.29 ; given_sigma_ion_O_4Pst_cm2( 557) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 557) =    0.00
  given_sigma_ion_O_4S_cm2( 558) =  3.39 ; given_sigma_ion_O_2D_cm2( 558) =    5.44 ; given_sigma_ion_O_2P_cm2( 558) =    3.27 ; given_sigma_ion_O_4Pst_cm2( 558) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 558) =    0.00
  given_sigma_ion_O_4S_cm2( 559) =  3.36 ; given_sigma_ion_O_2D_cm2( 559) =    5.40 ; given_sigma_ion_O_2P_cm2( 559) =    3.24 ; given_sigma_ion_O_4Pst_cm2( 559) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 559) =    0.00
  given_sigma_ion_O_4S_cm2( 560) =  3.40 ; given_sigma_ion_O_2D_cm2( 560) =    5.47 ; given_sigma_ion_O_2P_cm2( 560) =    3.28 ; given_sigma_ion_O_4Pst_cm2( 560) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 560) =    0.00
  given_sigma_ion_O_4S_cm2( 561) =  3.44 ; given_sigma_ion_O_2D_cm2( 561) =    5.52 ; given_sigma_ion_O_2P_cm2( 561) =    3.31 ; given_sigma_ion_O_4Pst_cm2( 561) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 561) =    0.00
  given_sigma_ion_O_4S_cm2( 562) =  3.43 ; given_sigma_ion_O_2D_cm2( 562) =    5.52 ; given_sigma_ion_O_2P_cm2( 562) =    3.31 ; given_sigma_ion_O_4Pst_cm2( 562) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 562) =    0.00
  given_sigma_ion_O_4S_cm2( 563) =  3.42 ; given_sigma_ion_O_2D_cm2( 563) =    5.50 ; given_sigma_ion_O_2P_cm2( 563) =    3.30 ; given_sigma_ion_O_4Pst_cm2( 563) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 563) =    0.00
  given_sigma_ion_O_4S_cm2( 564) =  3.41 ; given_sigma_ion_O_2D_cm2( 564) =    5.48 ; given_sigma_ion_O_2P_cm2( 564) =    3.29 ; given_sigma_ion_O_4Pst_cm2( 564) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 564) =    0.00
  given_sigma_ion_O_4S_cm2( 565) =  3.39 ; given_sigma_ion_O_2D_cm2( 565) =    5.46 ; given_sigma_ion_O_2P_cm2( 565) =    3.27 ; given_sigma_ion_O_4Pst_cm2( 565) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 565) =    0.00
  given_sigma_ion_O_4S_cm2( 566) =  3.38 ; given_sigma_ion_O_2D_cm2( 566) =    5.44 ; given_sigma_ion_O_2P_cm2( 566) =    3.26 ; given_sigma_ion_O_4Pst_cm2( 566) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 566) =    0.00
  given_sigma_ion_O_4S_cm2( 567) =  3.37 ; given_sigma_ion_O_2D_cm2( 567) =    5.42 ; given_sigma_ion_O_2P_cm2( 567) =    3.25 ; given_sigma_ion_O_4Pst_cm2( 567) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 567) =    0.00
  given_sigma_ion_O_4S_cm2( 568) =  3.36 ; given_sigma_ion_O_2D_cm2( 568) =    5.40 ; given_sigma_ion_O_2P_cm2( 568) =    3.24 ; given_sigma_ion_O_4Pst_cm2( 568) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 568) =    0.00
  given_sigma_ion_O_4S_cm2( 569) =  3.34 ; given_sigma_ion_O_2D_cm2( 569) =    5.38 ; given_sigma_ion_O_2P_cm2( 569) =    3.23 ; given_sigma_ion_O_4Pst_cm2( 569) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 569) =    0.00
  given_sigma_ion_O_4S_cm2( 570) =  3.34 ; given_sigma_ion_O_2D_cm2( 570) =    5.37 ; given_sigma_ion_O_2P_cm2( 570) =    3.22 ; given_sigma_ion_O_4Pst_cm2( 570) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 570) =    0.00
  given_sigma_ion_O_4S_cm2( 571) =  3.34 ; given_sigma_ion_O_2D_cm2( 571) =    5.37 ; given_sigma_ion_O_2P_cm2( 571) =    3.22 ; given_sigma_ion_O_4Pst_cm2( 571) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 571) =    0.00
  given_sigma_ion_O_4S_cm2( 572) =  3.33 ; given_sigma_ion_O_2D_cm2( 572) =    5.35 ; given_sigma_ion_O_2P_cm2( 572) =    3.21 ; given_sigma_ion_O_4Pst_cm2( 572) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 572) =    0.00
  given_sigma_ion_O_4S_cm2( 573) =  3.30 ; given_sigma_ion_O_2D_cm2( 573) =    5.30 ; given_sigma_ion_O_2P_cm2( 573) =    3.18 ; given_sigma_ion_O_4Pst_cm2( 573) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 573) =    0.00
  given_sigma_ion_O_4S_cm2( 574) =  3.29 ; given_sigma_ion_O_2D_cm2( 574) =    5.29 ; given_sigma_ion_O_2P_cm2( 574) =    3.18 ; given_sigma_ion_O_4Pst_cm2( 574) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 574) =    0.00
  given_sigma_ion_O_4S_cm2( 575) =  3.30 ; given_sigma_ion_O_2D_cm2( 575) =    5.30 ; given_sigma_ion_O_2P_cm2( 575) =    3.18 ; given_sigma_ion_O_4Pst_cm2( 575) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 575) =    0.00
  given_sigma_ion_O_4S_cm2( 576) =  3.30 ; given_sigma_ion_O_2D_cm2( 576) =    5.31 ; given_sigma_ion_O_2P_cm2( 576) =    3.18 ; given_sigma_ion_O_4Pst_cm2( 576) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 576) =    0.00
  given_sigma_ion_O_4S_cm2( 577) =  3.31 ; given_sigma_ion_O_2D_cm2( 577) =    5.31 ; given_sigma_ion_O_2P_cm2( 577) =    3.19 ; given_sigma_ion_O_4Pst_cm2( 577) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 577) =    0.00
  given_sigma_ion_O_4S_cm2( 578) =  3.31 ; given_sigma_ion_O_2D_cm2( 578) =    5.32 ; given_sigma_ion_O_2P_cm2( 578) =    3.19 ; given_sigma_ion_O_4Pst_cm2( 578) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 578) =    0.00
  given_sigma_ion_O_4S_cm2( 579) =  3.32 ; given_sigma_ion_O_2D_cm2( 579) =    5.34 ; given_sigma_ion_O_2P_cm2( 579) =    3.20 ; given_sigma_ion_O_4Pst_cm2( 579) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 579) =    0.00
  given_sigma_ion_O_4S_cm2( 580) =  3.33 ; given_sigma_ion_O_2D_cm2( 580) =    5.35 ; given_sigma_ion_O_2P_cm2( 580) =    3.21 ; given_sigma_ion_O_4Pst_cm2( 580) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 580) =    0.00
  given_sigma_ion_O_4S_cm2( 581) =  3.35 ; given_sigma_ion_O_2D_cm2( 581) =    5.38 ; given_sigma_ion_O_2P_cm2( 581) =    3.23 ; given_sigma_ion_O_4Pst_cm2( 581) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 581) =    0.00
  given_sigma_ion_O_4S_cm2( 582) =  3.36 ; given_sigma_ion_O_2D_cm2( 582) =    5.40 ; given_sigma_ion_O_2P_cm2( 582) =    3.24 ; given_sigma_ion_O_4Pst_cm2( 582) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 582) =    0.00
  given_sigma_ion_O_4S_cm2( 583) =  3.36 ; given_sigma_ion_O_2D_cm2( 583) =    5.40 ; given_sigma_ion_O_2P_cm2( 583) =    3.24 ; given_sigma_ion_O_4Pst_cm2( 583) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 583) =    0.00
  given_sigma_ion_O_4S_cm2( 584) =  3.36 ; given_sigma_ion_O_2D_cm2( 584) =    5.40 ; given_sigma_ion_O_2P_cm2( 584) =    3.24 ; given_sigma_ion_O_4Pst_cm2( 584) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 584) =    0.00
  given_sigma_ion_O_4S_cm2( 585) =  3.38 ; given_sigma_ion_O_2D_cm2( 585) =    5.43 ; given_sigma_ion_O_2P_cm2( 585) =    3.26 ; given_sigma_ion_O_4Pst_cm2( 585) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 585) =    0.00
  given_sigma_ion_O_4S_cm2( 586) =  3.43 ; given_sigma_ion_O_2D_cm2( 586) =    5.45 ; given_sigma_ion_O_2P_cm2( 586) =    3.22 ; given_sigma_ion_O_4Pst_cm2( 586) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 586) =    0.00
  given_sigma_ion_O_4S_cm2( 587) =  3.55 ; given_sigma_ion_O_2D_cm2( 587) =    5.51 ; given_sigma_ion_O_2P_cm2( 587) =    3.18 ; given_sigma_ion_O_4Pst_cm2( 587) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 587) =    0.00
  given_sigma_ion_O_4S_cm2( 588) =  3.57 ; given_sigma_ion_O_2D_cm2( 588) =    5.53 ; given_sigma_ion_O_2P_cm2( 588) =    3.20 ; given_sigma_ion_O_4Pst_cm2( 588) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 588) =    0.00
  given_sigma_ion_O_4S_cm2( 589) =  3.58 ; given_sigma_ion_O_2D_cm2( 589) =    5.56 ; given_sigma_ion_O_2P_cm2( 589) =    3.21 ; given_sigma_ion_O_4Pst_cm2( 589) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 589) =    0.00
  given_sigma_ion_O_4S_cm2( 590) =  3.59 ; given_sigma_ion_O_2D_cm2( 590) =    5.58 ; given_sigma_ion_O_2P_cm2( 590) =    3.22 ; given_sigma_ion_O_4Pst_cm2( 590) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 590) =    0.00
  given_sigma_ion_O_4S_cm2( 591) =  3.62 ; given_sigma_ion_O_2D_cm2( 591) =    5.62 ; given_sigma_ion_O_2P_cm2( 591) =    3.24 ; given_sigma_ion_O_4Pst_cm2( 591) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 591) =    0.00
  given_sigma_ion_O_4S_cm2( 592) =  3.63 ; given_sigma_ion_O_2D_cm2( 592) =    5.63 ; given_sigma_ion_O_2P_cm2( 592) =    3.25 ; given_sigma_ion_O_4Pst_cm2( 592) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 592) =    0.00
  given_sigma_ion_O_4S_cm2( 593) =  3.63 ; given_sigma_ion_O_2D_cm2( 593) =    5.64 ; given_sigma_ion_O_2P_cm2( 593) =    3.26 ; given_sigma_ion_O_4Pst_cm2( 593) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 593) =    0.00
  given_sigma_ion_O_4S_cm2( 594) =  3.65 ; given_sigma_ion_O_2D_cm2( 594) =    5.66 ; given_sigma_ion_O_2P_cm2( 594) =    3.27 ; given_sigma_ion_O_4Pst_cm2( 594) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 594) =    0.00
  given_sigma_ion_O_4S_cm2( 595) =  3.65 ; given_sigma_ion_O_2D_cm2( 595) =    5.67 ; given_sigma_ion_O_2P_cm2( 595) =    3.27 ; given_sigma_ion_O_4Pst_cm2( 595) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 595) =    0.00
  given_sigma_ion_O_4S_cm2( 596) =  3.66 ; given_sigma_ion_O_2D_cm2( 596) =    5.68 ; given_sigma_ion_O_2P_cm2( 596) =    3.28 ; given_sigma_ion_O_4Pst_cm2( 596) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 596) =    0.00
  given_sigma_ion_O_4S_cm2( 597) =  3.67 ; given_sigma_ion_O_2D_cm2( 597) =    5.70 ; given_sigma_ion_O_2P_cm2( 597) =    3.29 ; given_sigma_ion_O_4Pst_cm2( 597) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 597) =    0.00
  given_sigma_ion_O_4S_cm2( 598) =  3.68 ; given_sigma_ion_O_2D_cm2( 598) =    5.71 ; given_sigma_ion_O_2P_cm2( 598) =    3.30 ; given_sigma_ion_O_4Pst_cm2( 598) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 598) =    0.00
  given_sigma_ion_O_4S_cm2( 599) =  3.70 ; given_sigma_ion_O_2D_cm2( 599) =    5.74 ; given_sigma_ion_O_2P_cm2( 599) =    3.32 ; given_sigma_ion_O_4Pst_cm2( 599) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 599) =    0.00
  given_sigma_ion_O_4S_cm2( 600) =  3.71 ; given_sigma_ion_O_2D_cm2( 600) =    5.76 ; given_sigma_ion_O_2P_cm2( 600) =    3.33 ; given_sigma_ion_O_4Pst_cm2( 600) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 600) =    0.00
  given_sigma_ion_O_4S_cm2( 601) =  3.73 ; given_sigma_ion_O_2D_cm2( 601) =    5.79 ; given_sigma_ion_O_2P_cm2( 601) =    3.35 ; given_sigma_ion_O_4Pst_cm2( 601) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 601) =    0.00
  given_sigma_ion_O_4S_cm2( 602) =  3.74 ; given_sigma_ion_O_2D_cm2( 602) =    5.80 ; given_sigma_ion_O_2P_cm2( 602) =    3.35 ; given_sigma_ion_O_4Pst_cm2( 602) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 602) =    0.00
  given_sigma_ion_O_4S_cm2( 603) =  3.75 ; given_sigma_ion_O_2D_cm2( 603) =    5.82 ; given_sigma_ion_O_2P_cm2( 603) =    3.36 ; given_sigma_ion_O_4Pst_cm2( 603) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 603) =    0.00
  given_sigma_ion_O_4S_cm2( 604) =  3.76 ; given_sigma_ion_O_2D_cm2( 604) =    5.78 ; given_sigma_ion_O_2P_cm2( 604) =    3.41 ; given_sigma_ion_O_4Pst_cm2( 604) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 604) =    0.00
  given_sigma_ion_O_4S_cm2( 605) =  3.77 ; given_sigma_ion_O_2D_cm2( 605) =    5.72 ; given_sigma_ion_O_2P_cm2( 605) =    3.51 ; given_sigma_ion_O_4Pst_cm2( 605) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 605) =    0.00
  given_sigma_ion_O_4S_cm2( 606) =  3.79 ; given_sigma_ion_O_2D_cm2( 606) =    5.75 ; given_sigma_ion_O_2P_cm2( 606) =    3.47 ; given_sigma_ion_O_4Pst_cm2( 606) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 606) =    0.00
  given_sigma_ion_O_4S_cm2( 607) =  3.92 ; given_sigma_ion_O_2D_cm2( 607) =    6.02 ; given_sigma_ion_O_2P_cm2( 607) =    3.14 ; given_sigma_ion_O_4Pst_cm2( 607) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 607) =    0.00
  given_sigma_ion_O_4S_cm2( 608) =  3.92 ; given_sigma_ion_O_2D_cm2( 608) =    6.02 ; given_sigma_ion_O_2P_cm2( 608) =    3.15 ; given_sigma_ion_O_4Pst_cm2( 608) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 608) =    0.00
  given_sigma_ion_O_4S_cm2( 609) =  3.90 ; given_sigma_ion_O_2D_cm2( 609) =    6.03 ; given_sigma_ion_O_2P_cm2( 609) =    3.18 ; given_sigma_ion_O_4Pst_cm2( 609) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 609) =    0.00
  given_sigma_ion_O_4S_cm2( 610) =  3.88 ; given_sigma_ion_O_2D_cm2( 610) =    6.03 ; given_sigma_ion_O_2P_cm2( 610) =    3.21 ; given_sigma_ion_O_4Pst_cm2( 610) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 610) =    0.00
  given_sigma_ion_O_4S_cm2( 611) =  3.79 ; given_sigma_ion_O_2D_cm2( 611) =    6.06 ; given_sigma_ion_O_2P_cm2( 611) =    3.32 ; given_sigma_ion_O_4Pst_cm2( 611) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 611) =    0.00
  given_sigma_ion_O_4S_cm2( 612) =  3.76 ; given_sigma_ion_O_2D_cm2( 612) =    6.07 ; given_sigma_ion_O_2P_cm2( 612) =    3.37 ; given_sigma_ion_O_4Pst_cm2( 612) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 612) =    0.00
  given_sigma_ion_O_4S_cm2( 613) =  3.69 ; given_sigma_ion_O_2D_cm2( 613) =    6.08 ; given_sigma_ion_O_2P_cm2( 613) =    3.45 ; given_sigma_ion_O_4Pst_cm2( 613) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 613) =    0.00
  given_sigma_ion_O_4S_cm2( 614) =  3.64 ; given_sigma_ion_O_2D_cm2( 614) =    6.09 ; given_sigma_ion_O_2P_cm2( 614) =    3.51 ; given_sigma_ion_O_4Pst_cm2( 614) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 614) =    0.00
  given_sigma_ion_O_4S_cm2( 615) =  3.61 ; given_sigma_ion_O_2D_cm2( 615) =    6.10 ; given_sigma_ion_O_2P_cm2( 615) =    3.54 ; given_sigma_ion_O_4Pst_cm2( 615) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 615) =    0.00
  given_sigma_ion_O_4S_cm2( 616) =  3.56 ; given_sigma_ion_O_2D_cm2( 616) =    6.10 ; given_sigma_ion_O_2P_cm2( 616) =    3.61 ; given_sigma_ion_O_4Pst_cm2( 616) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 616) =    0.00
  given_sigma_ion_O_4S_cm2( 617) =  3.47 ; given_sigma_ion_O_2D_cm2( 617) =    6.12 ; given_sigma_ion_O_2P_cm2( 617) =    3.71 ; given_sigma_ion_O_4Pst_cm2( 617) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 617) =    0.00
  given_sigma_ion_O_4S_cm2( 618) =  3.46 ; given_sigma_ion_O_2D_cm2( 618) =    6.12 ; given_sigma_ion_O_2P_cm2( 618) =    3.72 ; given_sigma_ion_O_4Pst_cm2( 618) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 618) =    0.00
  given_sigma_ion_O_4S_cm2( 619) =  3.53 ; given_sigma_ion_O_2D_cm2( 619) =    6.14 ; given_sigma_ion_O_2P_cm2( 619) =    3.68 ; given_sigma_ion_O_4Pst_cm2( 619) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 619) =    0.00
  given_sigma_ion_O_4S_cm2( 620) =  3.59 ; given_sigma_ion_O_2D_cm2( 620) =    6.15 ; given_sigma_ion_O_2P_cm2( 620) =    3.64 ; given_sigma_ion_O_4Pst_cm2( 620) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 620) =    0.00
  given_sigma_ion_O_4S_cm2( 621) =  3.59 ; given_sigma_ion_O_2D_cm2( 621) =    6.16 ; given_sigma_ion_O_2P_cm2( 621) =    3.64 ; given_sigma_ion_O_4Pst_cm2( 621) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 621) =    0.00
  given_sigma_ion_O_4S_cm2( 622) =  3.60 ; given_sigma_ion_O_2D_cm2( 622) =    6.16 ; given_sigma_ion_O_2P_cm2( 622) =    3.63 ; given_sigma_ion_O_4Pst_cm2( 622) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 622) =    0.00
  given_sigma_ion_O_4S_cm2( 623) =  3.61 ; given_sigma_ion_O_2D_cm2( 623) =    6.16 ; given_sigma_ion_O_2P_cm2( 623) =    3.62 ; given_sigma_ion_O_4Pst_cm2( 623) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 623) =    0.00
  given_sigma_ion_O_4S_cm2( 624) =  3.62 ; given_sigma_ion_O_2D_cm2( 624) =    6.16 ; given_sigma_ion_O_2P_cm2( 624) =    3.62 ; given_sigma_ion_O_4Pst_cm2( 624) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 624) =    0.00
  given_sigma_ion_O_4S_cm2( 625) =  3.63 ; given_sigma_ion_O_2D_cm2( 625) =    6.16 ; given_sigma_ion_O_2P_cm2( 625) =    3.61 ; given_sigma_ion_O_4Pst_cm2( 625) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 625) =    0.00
  given_sigma_ion_O_4S_cm2( 626) =  3.63 ; given_sigma_ion_O_2D_cm2( 626) =    6.16 ; given_sigma_ion_O_2P_cm2( 626) =    3.61 ; given_sigma_ion_O_4Pst_cm2( 626) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 626) =    0.00
  given_sigma_ion_O_4S_cm2( 627) =  3.64 ; given_sigma_ion_O_2D_cm2( 627) =    6.16 ; given_sigma_ion_O_2P_cm2( 627) =    3.59 ; given_sigma_ion_O_4Pst_cm2( 627) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 627) =    0.00
  given_sigma_ion_O_4S_cm2( 628) =  3.65 ; given_sigma_ion_O_2D_cm2( 628) =    6.16 ; given_sigma_ion_O_2P_cm2( 628) =    3.59 ; given_sigma_ion_O_4Pst_cm2( 628) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 628) =    0.00
  given_sigma_ion_O_4S_cm2( 629) =  3.65 ; given_sigma_ion_O_2D_cm2( 629) =    6.16 ; given_sigma_ion_O_2P_cm2( 629) =    3.58 ; given_sigma_ion_O_4Pst_cm2( 629) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 629) =    0.00
  given_sigma_ion_O_4S_cm2( 630) =  3.66 ; given_sigma_ion_O_2D_cm2( 630) =    6.16 ; given_sigma_ion_O_2P_cm2( 630) =    3.58 ; given_sigma_ion_O_4Pst_cm2( 630) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 630) =    0.00
  given_sigma_ion_O_4S_cm2( 631) =  3.67 ; given_sigma_ion_O_2D_cm2( 631) =    6.16 ; given_sigma_ion_O_2P_cm2( 631) =    3.56 ; given_sigma_ion_O_4Pst_cm2( 631) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 631) =    0.00
  given_sigma_ion_O_4S_cm2( 632) =  3.68 ; given_sigma_ion_O_2D_cm2( 632) =    6.16 ; given_sigma_ion_O_2P_cm2( 632) =    3.55 ; given_sigma_ion_O_4Pst_cm2( 632) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 632) =    0.00
  given_sigma_ion_O_4S_cm2( 633) =  3.69 ; given_sigma_ion_O_2D_cm2( 633) =    6.16 ; given_sigma_ion_O_2P_cm2( 633) =    3.55 ; given_sigma_ion_O_4Pst_cm2( 633) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 633) =    0.00
  given_sigma_ion_O_4S_cm2( 634) =  3.70 ; given_sigma_ion_O_2D_cm2( 634) =    6.16 ; given_sigma_ion_O_2P_cm2( 634) =    3.54 ; given_sigma_ion_O_4Pst_cm2( 634) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 634) =    0.00
  given_sigma_ion_O_4S_cm2( 635) =  3.70 ; given_sigma_ion_O_2D_cm2( 635) =    6.16 ; given_sigma_ion_O_2P_cm2( 635) =    3.53 ; given_sigma_ion_O_4Pst_cm2( 635) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 635) =    0.00
  given_sigma_ion_O_4S_cm2( 636) =  3.71 ; given_sigma_ion_O_2D_cm2( 636) =    6.16 ; given_sigma_ion_O_2P_cm2( 636) =    3.53 ; given_sigma_ion_O_4Pst_cm2( 636) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 636) =    0.00
  given_sigma_ion_O_4S_cm2( 637) =  3.71 ; given_sigma_ion_O_2D_cm2( 637) =    6.16 ; given_sigma_ion_O_2P_cm2( 637) =    3.53 ; given_sigma_ion_O_4Pst_cm2( 637) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 637) =    0.00
  given_sigma_ion_O_4S_cm2( 638) =  3.72 ; given_sigma_ion_O_2D_cm2( 638) =    6.16 ; given_sigma_ion_O_2P_cm2( 638) =    3.51 ; given_sigma_ion_O_4Pst_cm2( 638) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 638) =    0.00
  given_sigma_ion_O_4S_cm2( 639) =  3.73 ; given_sigma_ion_O_2D_cm2( 639) =    6.16 ; given_sigma_ion_O_2P_cm2( 639) =    3.51 ; given_sigma_ion_O_4Pst_cm2( 639) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 639) =    0.00
  given_sigma_ion_O_4S_cm2( 640) =  3.73 ; given_sigma_ion_O_2D_cm2( 640) =    6.16 ; given_sigma_ion_O_2P_cm2( 640) =    3.51 ; given_sigma_ion_O_4Pst_cm2( 640) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 640) =    0.00
  given_sigma_ion_O_4S_cm2( 641) =  3.73 ; given_sigma_ion_O_2D_cm2( 641) =    6.16 ; given_sigma_ion_O_2P_cm2( 641) =    3.50 ; given_sigma_ion_O_4Pst_cm2( 641) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 641) =    0.00
  given_sigma_ion_O_4S_cm2( 642) =  3.74 ; given_sigma_ion_O_2D_cm2( 642) =    6.16 ; given_sigma_ion_O_2P_cm2( 642) =    3.50 ; given_sigma_ion_O_4Pst_cm2( 642) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 642) =    0.00
  given_sigma_ion_O_4S_cm2( 643) =  3.75 ; given_sigma_ion_O_2D_cm2( 643) =    6.16 ; given_sigma_ion_O_2P_cm2( 643) =    3.49 ; given_sigma_ion_O_4Pst_cm2( 643) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 643) =    0.00
  given_sigma_ion_O_4S_cm2( 644) =  3.75 ; given_sigma_ion_O_2D_cm2( 644) =    6.16 ; given_sigma_ion_O_2P_cm2( 644) =    3.48 ; given_sigma_ion_O_4Pst_cm2( 644) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 644) =    0.00
  given_sigma_ion_O_4S_cm2( 645) =  3.76 ; given_sigma_ion_O_2D_cm2( 645) =    6.17 ; given_sigma_ion_O_2P_cm2( 645) =    3.47 ; given_sigma_ion_O_4Pst_cm2( 645) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 645) =    0.00
  given_sigma_ion_O_4S_cm2( 646) =  3.77 ; given_sigma_ion_O_2D_cm2( 646) =    6.18 ; given_sigma_ion_O_2P_cm2( 646) =    3.46 ; given_sigma_ion_O_4Pst_cm2( 646) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 646) =    0.00
  given_sigma_ion_O_4S_cm2( 647) =  3.78 ; given_sigma_ion_O_2D_cm2( 647) =    6.19 ; given_sigma_ion_O_2P_cm2( 647) =    3.43 ; given_sigma_ion_O_4Pst_cm2( 647) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 647) =    0.00
  given_sigma_ion_O_4S_cm2( 648) =  3.79 ; given_sigma_ion_O_2D_cm2( 648) =    6.20 ; given_sigma_ion_O_2P_cm2( 648) =    3.40 ; given_sigma_ion_O_4Pst_cm2( 648) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 648) =    0.00
  given_sigma_ion_O_4S_cm2( 649) =  3.81 ; given_sigma_ion_O_2D_cm2( 649) =    6.22 ; given_sigma_ion_O_2P_cm2( 649) =    3.38 ; given_sigma_ion_O_4Pst_cm2( 649) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 649) =    0.00
  given_sigma_ion_O_4S_cm2( 650) =  3.81 ; given_sigma_ion_O_2D_cm2( 650) =    6.22 ; given_sigma_ion_O_2P_cm2( 650) =    3.36 ; given_sigma_ion_O_4Pst_cm2( 650) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 650) =    0.00
  given_sigma_ion_O_4S_cm2( 651) =  3.82 ; given_sigma_ion_O_2D_cm2( 651) =    6.23 ; given_sigma_ion_O_2P_cm2( 651) =    3.35 ; given_sigma_ion_O_4Pst_cm2( 651) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 651) =    0.00
  given_sigma_ion_O_4S_cm2( 652) =  3.83 ; given_sigma_ion_O_2D_cm2( 652) =    6.24 ; given_sigma_ion_O_2P_cm2( 652) =    3.33 ; given_sigma_ion_O_4Pst_cm2( 652) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 652) =    0.00
  given_sigma_ion_O_4S_cm2( 653) =  3.84 ; given_sigma_ion_O_2D_cm2( 653) =    6.25 ; given_sigma_ion_O_2P_cm2( 653) =    3.31 ; given_sigma_ion_O_4Pst_cm2( 653) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 653) =    0.00
  given_sigma_ion_O_4S_cm2( 654) =  3.85 ; given_sigma_ion_O_2D_cm2( 654) =    6.26 ; given_sigma_ion_O_2P_cm2( 654) =    3.30 ; given_sigma_ion_O_4Pst_cm2( 654) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 654) =    0.00
  given_sigma_ion_O_4S_cm2( 655) =  3.86 ; given_sigma_ion_O_2D_cm2( 655) =    6.27 ; given_sigma_ion_O_2P_cm2( 655) =    3.27 ; given_sigma_ion_O_4Pst_cm2( 655) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 655) =    0.00
  given_sigma_ion_O_4S_cm2( 656) =  3.87 ; given_sigma_ion_O_2D_cm2( 656) =    6.28 ; given_sigma_ion_O_2P_cm2( 656) =    3.24 ; given_sigma_ion_O_4Pst_cm2( 656) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 656) =    0.00
  given_sigma_ion_O_4S_cm2( 657) =  3.88 ; given_sigma_ion_O_2D_cm2( 657) =    6.29 ; given_sigma_ion_O_2P_cm2( 657) =    3.23 ; given_sigma_ion_O_4Pst_cm2( 657) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 657) =    0.00
  given_sigma_ion_O_4S_cm2( 658) =  3.88 ; given_sigma_ion_O_2D_cm2( 658) =    6.29 ; given_sigma_ion_O_2P_cm2( 658) =    3.22 ; given_sigma_ion_O_4Pst_cm2( 658) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 658) =    0.00
  given_sigma_ion_O_4S_cm2( 659) =  3.89 ; given_sigma_ion_O_2D_cm2( 659) =    6.30 ; given_sigma_ion_O_2P_cm2( 659) =    3.22 ; given_sigma_ion_O_4Pst_cm2( 659) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 659) =    0.00
  given_sigma_ion_O_4S_cm2( 660) =  3.95 ; given_sigma_ion_O_2D_cm2( 660) =    6.29 ; given_sigma_ion_O_2P_cm2( 660) =    3.15 ; given_sigma_ion_O_4Pst_cm2( 660) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 660) =    0.00
  given_sigma_ion_O_4S_cm2( 661) =  4.11 ; given_sigma_ion_O_2D_cm2( 661) =    6.27 ; given_sigma_ion_O_2P_cm2( 661) =    3.01 ; given_sigma_ion_O_4Pst_cm2( 661) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 661) =    0.00
  given_sigma_ion_O_4S_cm2( 662) =  4.25 ; given_sigma_ion_O_2D_cm2( 662) =    6.25 ; given_sigma_ion_O_2P_cm2( 662) =    2.89 ; given_sigma_ion_O_4Pst_cm2( 662) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 662) =    0.00
  given_sigma_ion_O_4S_cm2( 663) =  4.34 ; given_sigma_ion_O_2D_cm2( 663) =    6.24 ; given_sigma_ion_O_2P_cm2( 663) =    2.81 ; given_sigma_ion_O_4Pst_cm2( 663) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 663) =    0.00
  given_sigma_ion_O_4S_cm2( 664) =  4.56 ; given_sigma_ion_O_2D_cm2( 664) =    6.20 ; given_sigma_ion_O_2P_cm2( 664) =    2.61 ; given_sigma_ion_O_4Pst_cm2( 664) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 664) =    0.00
  given_sigma_ion_O_4S_cm2( 665) =  4.78 ; given_sigma_ion_O_2D_cm2( 665) =    6.17 ; given_sigma_ion_O_2P_cm2( 665) =    2.40 ; given_sigma_ion_O_4Pst_cm2( 665) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 665) =    0.00
  given_sigma_ion_O_4S_cm2( 666) =  4.87 ; given_sigma_ion_O_2D_cm2( 666) =    6.16 ; given_sigma_ion_O_2P_cm2( 666) =    2.32 ; given_sigma_ion_O_4Pst_cm2( 666) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 666) =    0.00
  given_sigma_ion_O_4S_cm2( 667) =  4.94 ; given_sigma_ion_O_2D_cm2( 667) =    6.15 ; given_sigma_ion_O_2P_cm2( 667) =    2.26 ; given_sigma_ion_O_4Pst_cm2( 667) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 667) =    0.00
  given_sigma_ion_O_4S_cm2( 668) =  5.01 ; given_sigma_ion_O_2D_cm2( 668) =    6.14 ; given_sigma_ion_O_2P_cm2( 668) =    2.20 ; given_sigma_ion_O_4Pst_cm2( 668) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 668) =    0.00
  given_sigma_ion_O_4S_cm2( 669) =  5.18 ; given_sigma_ion_O_2D_cm2( 669) =    6.12 ; given_sigma_ion_O_2P_cm2( 669) =    2.04 ; given_sigma_ion_O_4Pst_cm2( 669) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 669) =    0.00
  given_sigma_ion_O_4S_cm2( 670) =  5.23 ; given_sigma_ion_O_2D_cm2( 670) =    6.11 ; given_sigma_ion_O_2P_cm2( 670) =    2.00 ; given_sigma_ion_O_4Pst_cm2( 670) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 670) =    0.00
  given_sigma_ion_O_4S_cm2( 671) =  5.30 ; given_sigma_ion_O_2D_cm2( 671) =    6.10 ; given_sigma_ion_O_2P_cm2( 671) =    1.94 ; given_sigma_ion_O_4Pst_cm2( 671) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 671) =    0.00
  given_sigma_ion_O_4S_cm2( 672) =  5.45 ; given_sigma_ion_O_2D_cm2( 672) =    6.08 ; given_sigma_ion_O_2P_cm2( 672) =    1.80 ; given_sigma_ion_O_4Pst_cm2( 672) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 672) =    0.00
  given_sigma_ion_O_4S_cm2( 673) =  5.52 ; given_sigma_ion_O_2D_cm2( 673) =    6.07 ; given_sigma_ion_O_2P_cm2( 673) =    1.74 ; given_sigma_ion_O_4Pst_cm2( 673) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 673) =    0.00
  given_sigma_ion_O_4S_cm2( 674) =  5.67 ; given_sigma_ion_O_2D_cm2( 674) =    6.05 ; given_sigma_ion_O_2P_cm2( 674) =    1.60 ; given_sigma_ion_O_4Pst_cm2( 674) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 674) =    0.00
  given_sigma_ion_O_4S_cm2( 675) =  5.79 ; given_sigma_ion_O_2D_cm2( 675) =    6.03 ; given_sigma_ion_O_2P_cm2( 675) =    1.50 ; given_sigma_ion_O_4Pst_cm2( 675) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 675) =    0.00
  given_sigma_ion_O_4S_cm2( 676) =  5.83 ; given_sigma_ion_O_2D_cm2( 676) =    6.03 ; given_sigma_ion_O_2P_cm2( 676) =    1.46 ; given_sigma_ion_O_4Pst_cm2( 676) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 676) =    0.00
  given_sigma_ion_O_4S_cm2( 677) =  5.90 ; given_sigma_ion_O_2D_cm2( 677) =    6.02 ; given_sigma_ion_O_2P_cm2( 677) =    1.40 ; given_sigma_ion_O_4Pst_cm2( 677) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 677) =    0.00
  given_sigma_ion_O_4S_cm2( 678) =  6.12 ; given_sigma_ion_O_2D_cm2( 678) =    5.99 ; given_sigma_ion_O_2P_cm2( 678) =    1.20 ; given_sigma_ion_O_4Pst_cm2( 678) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 678) =    0.00
  given_sigma_ion_O_4S_cm2( 679) =  6.01 ; given_sigma_ion_O_2D_cm2( 679) =    5.97 ; given_sigma_ion_O_2P_cm2( 679) =    1.30 ; given_sigma_ion_O_4Pst_cm2( 679) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 679) =    0.00
  given_sigma_ion_O_4S_cm2( 680) =  5.88 ; given_sigma_ion_O_2D_cm2( 680) =    5.96 ; given_sigma_ion_O_2P_cm2( 680) =    1.43 ; given_sigma_ion_O_4Pst_cm2( 680) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 680) =    0.00
  given_sigma_ion_O_4S_cm2( 681) =  5.66 ; given_sigma_ion_O_2D_cm2( 681) =    5.94 ; given_sigma_ion_O_2P_cm2( 681) =    1.65 ; given_sigma_ion_O_4Pst_cm2( 681) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 681) =    0.00
  given_sigma_ion_O_4S_cm2( 682) =  5.61 ; given_sigma_ion_O_2D_cm2( 682) =    5.93 ; given_sigma_ion_O_2P_cm2( 682) =    1.69 ; given_sigma_ion_O_4Pst_cm2( 682) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 682) =    0.00
  given_sigma_ion_O_4S_cm2( 683) =  5.49 ; given_sigma_ion_O_2D_cm2( 683) =    5.92 ; given_sigma_ion_O_2P_cm2( 683) =    1.82 ; given_sigma_ion_O_4Pst_cm2( 683) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 683) =    0.00
  given_sigma_ion_O_4S_cm2( 684) =  5.36 ; given_sigma_ion_O_2D_cm2( 684) =    5.90 ; given_sigma_ion_O_2P_cm2( 684) =    1.94 ; given_sigma_ion_O_4Pst_cm2( 684) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 684) =    0.00
  given_sigma_ion_O_4S_cm2( 685) =  5.11 ; given_sigma_ion_O_2D_cm2( 685) =    5.88 ; given_sigma_ion_O_2P_cm2( 685) =    2.19 ; given_sigma_ion_O_4Pst_cm2( 685) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 685) =    0.00
  given_sigma_ion_O_4S_cm2( 686) =  5.06 ; given_sigma_ion_O_2D_cm2( 686) =    5.87 ; given_sigma_ion_O_2P_cm2( 686) =    2.24 ; given_sigma_ion_O_4Pst_cm2( 686) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 686) =    0.00
  given_sigma_ion_O_4S_cm2( 687) =  5.01 ; given_sigma_ion_O_2D_cm2( 687) =    5.87 ; given_sigma_ion_O_2P_cm2( 687) =    2.29 ; given_sigma_ion_O_4Pst_cm2( 687) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 687) =    0.00
  given_sigma_ion_O_4S_cm2( 688) =  4.87 ; given_sigma_ion_O_2D_cm2( 688) =    5.85 ; given_sigma_ion_O_2P_cm2( 688) =    2.43 ; given_sigma_ion_O_4Pst_cm2( 688) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 688) =    0.00
  given_sigma_ion_O_4S_cm2( 689) =  4.64 ; given_sigma_ion_O_2D_cm2( 689) =    5.83 ; given_sigma_ion_O_2P_cm2( 689) =    2.65 ; given_sigma_ion_O_4Pst_cm2( 689) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 689) =    0.00
  given_sigma_ion_O_4S_cm2( 690) =  4.47 ; given_sigma_ion_O_2D_cm2( 690) =    5.81 ; given_sigma_ion_O_2P_cm2( 690) =    2.82 ; given_sigma_ion_O_4Pst_cm2( 690) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 690) =    0.00
  given_sigma_ion_O_4S_cm2( 691) =  4.45 ; given_sigma_ion_O_2D_cm2( 691) =    5.81 ; given_sigma_ion_O_2P_cm2( 691) =    2.85 ; given_sigma_ion_O_4Pst_cm2( 691) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 691) =    0.00
  given_sigma_ion_O_4S_cm2( 692) =  4.37 ; given_sigma_ion_O_2D_cm2( 692) =    5.80 ; given_sigma_ion_O_2P_cm2( 692) =    2.92 ; given_sigma_ion_O_4Pst_cm2( 692) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 692) =    0.00
  given_sigma_ion_O_4S_cm2( 693) =  4.25 ; given_sigma_ion_O_2D_cm2( 693) =    5.79 ; given_sigma_ion_O_2P_cm2( 693) =    3.04 ; given_sigma_ion_O_4Pst_cm2( 693) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 693) =    0.00
  given_sigma_ion_O_4S_cm2( 694) =  4.13 ; given_sigma_ion_O_2D_cm2( 694) =    5.77 ; given_sigma_ion_O_2P_cm2( 694) =    3.16 ; given_sigma_ion_O_4Pst_cm2( 694) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 694) =    0.00
  given_sigma_ion_O_4S_cm2( 695) =  3.88 ; given_sigma_ion_O_2D_cm2( 695) =    5.75 ; given_sigma_ion_O_2P_cm2( 695) =    3.40 ; given_sigma_ion_O_4Pst_cm2( 695) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 695) =    0.00
  given_sigma_ion_O_4S_cm2( 696) =  3.79 ; given_sigma_ion_O_2D_cm2( 696) =    5.74 ; given_sigma_ion_O_2P_cm2( 696) =    3.50 ; given_sigma_ion_O_4Pst_cm2( 696) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 696) =    0.00
  given_sigma_ion_O_4S_cm2( 697) =  3.64 ; given_sigma_ion_O_2D_cm2( 697) =    5.72 ; given_sigma_ion_O_2P_cm2( 697) =    3.64 ; given_sigma_ion_O_4Pst_cm2( 697) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 697) =    0.00
  given_sigma_ion_O_4S_cm2( 698) =  3.63 ; given_sigma_ion_O_2D_cm2( 698) =    5.72 ; given_sigma_ion_O_2P_cm2( 698) =    3.64 ; given_sigma_ion_O_4Pst_cm2( 698) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 698) =    0.00
  given_sigma_ion_O_4S_cm2( 699) =  3.60 ; given_sigma_ion_O_2D_cm2( 699) =    5.72 ; given_sigma_ion_O_2P_cm2( 699) =    3.64 ; given_sigma_ion_O_4Pst_cm2( 699) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 699) =    0.00
  given_sigma_ion_O_4S_cm2( 700) =  3.57 ; given_sigma_ion_O_2D_cm2( 700) =    5.71 ; given_sigma_ion_O_2P_cm2( 700) =    3.64 ; given_sigma_ion_O_4Pst_cm2( 700) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 700) =    0.00
  given_sigma_ion_O_4S_cm2( 701) =  3.57 ; given_sigma_ion_O_2D_cm2( 701) =    5.71 ; given_sigma_ion_O_2P_cm2( 701) =    3.64 ; given_sigma_ion_O_4Pst_cm2( 701) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 701) =    0.00
  given_sigma_ion_O_4S_cm2( 702) =  3.53 ; given_sigma_ion_O_2D_cm2( 702) =    5.71 ; given_sigma_ion_O_2P_cm2( 702) =    3.65 ; given_sigma_ion_O_4Pst_cm2( 702) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 702) =    0.00
  given_sigma_ion_O_4S_cm2( 703) =  3.49 ; given_sigma_ion_O_2D_cm2( 703) =    5.70 ; given_sigma_ion_O_2P_cm2( 703) =    3.65 ; given_sigma_ion_O_4Pst_cm2( 703) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 703) =    0.00
  given_sigma_ion_O_4S_cm2( 704) =  3.46 ; given_sigma_ion_O_2D_cm2( 704) =    5.70 ; given_sigma_ion_O_2P_cm2( 704) =    3.65 ; given_sigma_ion_O_4Pst_cm2( 704) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 704) =    0.00
  given_sigma_ion_O_4S_cm2( 705) =  3.42 ; given_sigma_ion_O_2D_cm2( 705) =    5.69 ; given_sigma_ion_O_2P_cm2( 705) =    3.65 ; given_sigma_ion_O_4Pst_cm2( 705) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 705) =    0.00
  given_sigma_ion_O_4S_cm2( 706) =  3.38 ; given_sigma_ion_O_2D_cm2( 706) =    5.69 ; given_sigma_ion_O_2P_cm2( 706) =    3.65 ; given_sigma_ion_O_4Pst_cm2( 706) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 706) =    0.00
  given_sigma_ion_O_4S_cm2( 707) =  3.37 ; given_sigma_ion_O_2D_cm2( 707) =    5.68 ; given_sigma_ion_O_2P_cm2( 707) =    3.65 ; given_sigma_ion_O_4Pst_cm2( 707) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 707) =    0.00
  given_sigma_ion_O_4S_cm2( 708) =  3.35 ; given_sigma_ion_O_2D_cm2( 708) =    5.68 ; given_sigma_ion_O_2P_cm2( 708) =    3.65 ; given_sigma_ion_O_4Pst_cm2( 708) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 708) =    0.00
  given_sigma_ion_O_4S_cm2( 709) =  3.31 ; given_sigma_ion_O_2D_cm2( 709) =    5.68 ; given_sigma_ion_O_2P_cm2( 709) =    3.65 ; given_sigma_ion_O_4Pst_cm2( 709) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 709) =    0.00
  given_sigma_ion_O_4S_cm2( 710) =  3.28 ; given_sigma_ion_O_2D_cm2( 710) =    5.67 ; given_sigma_ion_O_2P_cm2( 710) =    3.65 ; given_sigma_ion_O_4Pst_cm2( 710) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 710) =    0.00
  given_sigma_ion_O_4S_cm2( 711) =  3.31 ; given_sigma_ion_O_2D_cm2( 711) =    5.65 ; given_sigma_ion_O_2P_cm2( 711) =    3.60 ; given_sigma_ion_O_4Pst_cm2( 711) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 711) =    0.00
  given_sigma_ion_O_4S_cm2( 712) =  3.40 ; given_sigma_ion_O_2D_cm2( 712) =    5.59 ; given_sigma_ion_O_2P_cm2( 712) =    3.47 ; given_sigma_ion_O_4Pst_cm2( 712) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 712) =    0.00
  given_sigma_ion_O_4S_cm2( 713) =  3.44 ; given_sigma_ion_O_2D_cm2( 713) =    5.56 ; given_sigma_ion_O_2P_cm2( 713) =    3.40 ; given_sigma_ion_O_4Pst_cm2( 713) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 713) =    0.00
  given_sigma_ion_O_4S_cm2( 714) =  3.50 ; given_sigma_ion_O_2D_cm2( 714) =    5.53 ; given_sigma_ion_O_2P_cm2( 714) =    3.31 ; given_sigma_ion_O_4Pst_cm2( 714) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 714) =    0.00
  given_sigma_ion_O_4S_cm2( 715) =  3.63 ; given_sigma_ion_O_2D_cm2( 715) =    5.44 ; given_sigma_ion_O_2P_cm2( 715) =    3.11 ; given_sigma_ion_O_4Pst_cm2( 715) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 715) =    0.00
  given_sigma_ion_O_4S_cm2( 716) =  3.74 ; given_sigma_ion_O_2D_cm2( 716) =    5.37 ; given_sigma_ion_O_2P_cm2( 716) =    2.94 ; given_sigma_ion_O_4Pst_cm2( 716) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 716) =    0.00
  given_sigma_ion_O_4S_cm2( 717) =  3.81 ; given_sigma_ion_O_2D_cm2( 717) =    5.32 ; given_sigma_ion_O_2P_cm2( 717) =    2.84 ; given_sigma_ion_O_4Pst_cm2( 717) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 717) =    0.00
  given_sigma_ion_O_4S_cm2( 718) =  3.84 ; given_sigma_ion_O_2D_cm2( 718) =    5.30 ; given_sigma_ion_O_2P_cm2( 718) =    2.78 ; given_sigma_ion_O_4Pst_cm2( 718) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 718) =    0.00
  given_sigma_ion_O_4S_cm2( 719) =  4.16 ; given_sigma_ion_O_2D_cm2( 719) =    5.78 ; given_sigma_ion_O_2P_cm2( 719) =    1.93 ; given_sigma_ion_O_4Pst_cm2( 719) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 719) =    0.00
  given_sigma_ion_O_4S_cm2( 720) =  4.67 ; given_sigma_ion_O_2D_cm2( 720) =    6.58 ; given_sigma_ion_O_2P_cm2( 720) =    0.55 ; given_sigma_ion_O_4Pst_cm2( 720) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 720) =    0.00
  given_sigma_ion_O_4S_cm2( 721) =  4.87 ; given_sigma_ion_O_2D_cm2( 721) =    6.89 ; given_sigma_ion_O_2P_cm2( 721) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 721) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 721) =    0.00
  given_sigma_ion_O_4S_cm2( 722) =  4.77 ; given_sigma_ion_O_2D_cm2( 722) =    6.85 ; given_sigma_ion_O_2P_cm2( 722) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 722) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 722) =    0.00
  given_sigma_ion_O_4S_cm2( 723) =  4.74 ; given_sigma_ion_O_2D_cm2( 723) =    6.83 ; given_sigma_ion_O_2P_cm2( 723) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 723) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 723) =    0.00
  given_sigma_ion_O_4S_cm2( 724) =  4.66 ; given_sigma_ion_O_2D_cm2( 724) =    6.80 ; given_sigma_ion_O_2P_cm2( 724) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 724) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 724) =    0.00
  given_sigma_ion_O_4S_cm2( 725) =  4.56 ; given_sigma_ion_O_2D_cm2( 725) =    6.75 ; given_sigma_ion_O_2P_cm2( 725) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 725) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 725) =    0.00
  given_sigma_ion_O_4S_cm2( 726) =  4.50 ; given_sigma_ion_O_2D_cm2( 726) =    6.72 ; given_sigma_ion_O_2P_cm2( 726) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 726) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 726) =    0.00
  given_sigma_ion_O_4S_cm2( 727) =  4.46 ; given_sigma_ion_O_2D_cm2( 727) =    6.70 ; given_sigma_ion_O_2P_cm2( 727) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 727) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 727) =    0.00
  given_sigma_ion_O_4S_cm2( 728) =  4.42 ; given_sigma_ion_O_2D_cm2( 728) =    6.67 ; given_sigma_ion_O_2P_cm2( 728) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 728) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 728) =    0.00
  given_sigma_ion_O_4S_cm2( 729) =  4.36 ; given_sigma_ion_O_2D_cm2( 729) =    6.64 ; given_sigma_ion_O_2P_cm2( 729) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 729) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 729) =    0.00
  given_sigma_ion_O_4S_cm2( 730) =  4.31 ; given_sigma_ion_O_2D_cm2( 730) =    6.62 ; given_sigma_ion_O_2P_cm2( 730) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 730) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 730) =    0.00
  given_sigma_ion_O_4S_cm2( 731) =  4.27 ; given_sigma_ion_O_2D_cm2( 731) =    6.59 ; given_sigma_ion_O_2P_cm2( 731) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 731) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 731) =    0.00
  given_sigma_ion_O_4S_cm2( 732) =  4.17 ; given_sigma_ion_O_2D_cm2( 732) =    6.54 ; given_sigma_ion_O_2P_cm2( 732) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 732) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 732) =    0.00
  given_sigma_ion_O_4S_cm2( 733) =  4.10 ; given_sigma_ion_O_2D_cm2( 733) =    6.50 ; given_sigma_ion_O_2P_cm2( 733) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 733) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 733) =    0.00
  given_sigma_ion_O_4S_cm2( 734) =  4.08 ; given_sigma_ion_O_2D_cm2( 734) =    6.49 ; given_sigma_ion_O_2P_cm2( 734) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 734) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 734) =    0.00
  given_sigma_ion_O_4S_cm2( 735) =  4.06 ; given_sigma_ion_O_2D_cm2( 735) =    6.48 ; given_sigma_ion_O_2P_cm2( 735) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 735) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 735) =    0.00
  given_sigma_ion_O_4S_cm2( 736) =  4.02 ; given_sigma_ion_O_2D_cm2( 736) =    6.45 ; given_sigma_ion_O_2P_cm2( 736) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 736) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 736) =    0.00
  given_sigma_ion_O_4S_cm2( 737) =  3.96 ; given_sigma_ion_O_2D_cm2( 737) =    6.42 ; given_sigma_ion_O_2P_cm2( 737) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 737) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 737) =    0.00
  given_sigma_ion_O_4S_cm2( 738) =  3.94 ; given_sigma_ion_O_2D_cm2( 738) =    6.41 ; given_sigma_ion_O_2P_cm2( 738) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 738) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 738) =    0.00
  given_sigma_ion_O_4S_cm2( 739) =  3.90 ; given_sigma_ion_O_2D_cm2( 739) =    6.38 ; given_sigma_ion_O_2P_cm2( 739) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 739) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 739) =    0.00
  given_sigma_ion_O_4S_cm2( 740) =  3.87 ; given_sigma_ion_O_2D_cm2( 740) =    6.36 ; given_sigma_ion_O_2P_cm2( 740) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 740) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 740) =    0.00
  given_sigma_ion_O_4S_cm2( 741) =  3.85 ; given_sigma_ion_O_2D_cm2( 741) =    6.35 ; given_sigma_ion_O_2P_cm2( 741) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 741) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 741) =    0.00
  given_sigma_ion_O_4S_cm2( 742) =  3.81 ; given_sigma_ion_O_2D_cm2( 742) =    6.32 ; given_sigma_ion_O_2P_cm2( 742) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 742) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 742) =    0.00
  given_sigma_ion_O_4S_cm2( 743) =  3.77 ; given_sigma_ion_O_2D_cm2( 743) =    6.30 ; given_sigma_ion_O_2P_cm2( 743) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 743) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 743) =    0.00
  given_sigma_ion_O_4S_cm2( 744) =  3.73 ; given_sigma_ion_O_2D_cm2( 744) =    6.27 ; given_sigma_ion_O_2P_cm2( 744) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 744) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 744) =    0.00
  given_sigma_ion_O_4S_cm2( 745) =  3.82 ; given_sigma_ion_O_2D_cm2( 745) =    6.47 ; given_sigma_ion_O_2P_cm2( 745) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 745) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 745) =    0.00
  given_sigma_ion_O_4S_cm2( 746) =  3.92 ; given_sigma_ion_O_2D_cm2( 746) =    6.68 ; given_sigma_ion_O_2P_cm2( 746) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 746) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 746) =    0.00
  given_sigma_ion_O_4S_cm2( 747) =  3.96 ; given_sigma_ion_O_2D_cm2( 747) =    8.04 ; given_sigma_ion_O_2P_cm2( 747) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 747) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 747) =    0.00
  given_sigma_ion_O_4S_cm2( 748) =  5.48 ; given_sigma_ion_O_2D_cm2( 748) =   12.19 ; given_sigma_ion_O_2P_cm2( 748) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 748) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 748) =    0.00
  given_sigma_ion_O_4S_cm2( 749) =  4.47 ; given_sigma_ion_O_2D_cm2( 749) =   10.93 ; given_sigma_ion_O_2P_cm2( 749) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 749) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 749) =    0.00
  given_sigma_ion_O_4S_cm2( 750) =  4.43 ; given_sigma_ion_O_2D_cm2( 750) =   10.09 ; given_sigma_ion_O_2P_cm2( 750) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 750) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 750) =    0.00
  given_sigma_ion_O_4S_cm2( 751) =  4.41 ; given_sigma_ion_O_2D_cm2( 751) =    9.74 ; given_sigma_ion_O_2P_cm2( 751) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 751) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 751) =    0.00
  given_sigma_ion_O_4S_cm2( 752) =  4.39 ; given_sigma_ion_O_2D_cm2( 752) =    9.51 ; given_sigma_ion_O_2P_cm2( 752) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 752) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 752) =    0.00
  given_sigma_ion_O_4S_cm2( 753) =  4.61 ; given_sigma_ion_O_2D_cm2( 753) =    9.68 ; given_sigma_ion_O_2P_cm2( 753) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 753) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 753) =    0.00
  given_sigma_ion_O_4S_cm2( 754) =  4.83 ; given_sigma_ion_O_2D_cm2( 754) =    9.85 ; given_sigma_ion_O_2P_cm2( 754) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 754) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 754) =    0.00
  given_sigma_ion_O_4S_cm2( 755) =  5.05 ; given_sigma_ion_O_2D_cm2( 755) =   10.02 ; given_sigma_ion_O_2P_cm2( 755) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 755) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 755) =    0.00
  given_sigma_ion_O_4S_cm2( 756) =  5.13 ; given_sigma_ion_O_2D_cm2( 756) =   10.07 ; given_sigma_ion_O_2P_cm2( 756) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 756) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 756) =    0.00
  given_sigma_ion_O_4S_cm2( 757) =  4.93 ; given_sigma_ion_O_2D_cm2( 757) =    9.49 ; given_sigma_ion_O_2P_cm2( 757) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 757) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 757) =    0.00
  given_sigma_ion_O_4S_cm2( 758) =  4.82 ; given_sigma_ion_O_2D_cm2( 758) =    9.20 ; given_sigma_ion_O_2P_cm2( 758) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 758) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 758) =    0.00
  given_sigma_ion_O_4S_cm2( 759) =  4.50 ; given_sigma_ion_O_2D_cm2( 759) =    8.34 ; given_sigma_ion_O_2P_cm2( 759) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 759) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 759) =    0.00
  given_sigma_ion_O_4S_cm2( 760) =  4.16 ; given_sigma_ion_O_2D_cm2( 760) =    7.50 ; given_sigma_ion_O_2P_cm2( 760) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 760) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 760) =    0.00
  given_sigma_ion_O_4S_cm2( 761) =  3.56 ; given_sigma_ion_O_2D_cm2( 761) =    6.13 ; given_sigma_ion_O_2P_cm2( 761) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 761) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 761) =    0.00
  given_sigma_ion_O_4S_cm2( 762) =  3.59 ; given_sigma_ion_O_2D_cm2( 762) =    6.06 ; given_sigma_ion_O_2P_cm2( 762) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 762) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 762) =    0.00
  given_sigma_ion_O_4S_cm2( 763) =  3.60 ; given_sigma_ion_O_2D_cm2( 763) =    5.92 ; given_sigma_ion_O_2P_cm2( 763) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 763) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 763) =    0.00
  given_sigma_ion_O_4S_cm2( 764) =  3.62 ; given_sigma_ion_O_2D_cm2( 764) =    5.68 ; given_sigma_ion_O_2P_cm2( 764) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 764) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 764) =    0.00
  given_sigma_ion_O_4S_cm2( 765) =  3.62 ; given_sigma_ion_O_2D_cm2( 765) =    5.63 ; given_sigma_ion_O_2P_cm2( 765) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 765) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 765) =    0.00
  given_sigma_ion_O_4S_cm2( 766) =  3.88 ; given_sigma_ion_O_2D_cm2( 766) =    5.82 ; given_sigma_ion_O_2P_cm2( 766) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 766) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 766) =    0.00
  given_sigma_ion_O_4S_cm2( 767) =  3.92 ; given_sigma_ion_O_2D_cm2( 767) =    6.28 ; given_sigma_ion_O_2P_cm2( 767) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 767) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 767) =    0.00
  given_sigma_ion_O_4S_cm2( 768) =  4.03 ; given_sigma_ion_O_2D_cm2( 768) =    7.17 ; given_sigma_ion_O_2P_cm2( 768) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 768) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 768) =    0.00
  given_sigma_ion_O_4S_cm2( 769) =  4.58 ; given_sigma_ion_O_2D_cm2( 769) =    8.55 ; given_sigma_ion_O_2P_cm2( 769) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 769) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 769) =    0.00
  given_sigma_ion_O_4S_cm2( 770) =  6.94 ; given_sigma_ion_O_2D_cm2( 770) =   17.74 ; given_sigma_ion_O_2P_cm2( 770) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 770) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 770) =    0.00
  given_sigma_ion_O_4S_cm2( 771) =  6.94 ; given_sigma_ion_O_2D_cm2( 771) =   18.55 ; given_sigma_ion_O_2P_cm2( 771) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 771) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 771) =    0.00
  given_sigma_ion_O_4S_cm2( 772) =  6.28 ; given_sigma_ion_O_2D_cm2( 772) =   16.29 ; given_sigma_ion_O_2P_cm2( 772) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 772) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 772) =    0.00
  given_sigma_ion_O_4S_cm2( 773) =  5.55 ; given_sigma_ion_O_2D_cm2( 773) =   13.98 ; given_sigma_ion_O_2P_cm2( 773) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 773) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 773) =    0.00
  given_sigma_ion_O_4S_cm2( 774) =  4.79 ; given_sigma_ion_O_2D_cm2( 774) =   11.72 ; given_sigma_ion_O_2P_cm2( 774) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 774) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 774) =    0.00
  given_sigma_ion_O_4S_cm2( 775) =  4.53 ; given_sigma_ion_O_2D_cm2( 775) =   10.97 ; given_sigma_ion_O_2P_cm2( 775) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 775) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 775) =    0.00
  given_sigma_ion_O_4S_cm2( 776) =  3.98 ; given_sigma_ion_O_2D_cm2( 776) =    9.37 ; given_sigma_ion_O_2P_cm2( 776) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 776) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 776) =    0.00
  given_sigma_ion_O_4S_cm2( 777) =  3.41 ; given_sigma_ion_O_2D_cm2( 777) =    7.79 ; given_sigma_ion_O_2P_cm2( 777) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 777) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 777) =    0.00
  given_sigma_ion_O_4S_cm2( 778) =  3.95 ; given_sigma_ion_O_2D_cm2( 778) =    8.52 ; given_sigma_ion_O_2P_cm2( 778) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 778) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 778) =    0.00
  given_sigma_ion_O_4S_cm2( 779) =  4.06 ; given_sigma_ion_O_2D_cm2( 779) =    8.35 ; given_sigma_ion_O_2P_cm2( 779) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 779) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 779) =    0.00
  given_sigma_ion_O_4S_cm2( 780) =  3.97 ; given_sigma_ion_O_2D_cm2( 780) =    7.95 ; given_sigma_ion_O_2P_cm2( 780) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 780) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 780) =    0.00
  given_sigma_ion_O_4S_cm2( 781) =  3.88 ; given_sigma_ion_O_2D_cm2( 781) =    7.56 ; given_sigma_ion_O_2P_cm2( 781) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 781) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 781) =    0.00
  given_sigma_ion_O_4S_cm2( 782) =  3.52 ; given_sigma_ion_O_2D_cm2( 782) =    6.28 ; given_sigma_ion_O_2P_cm2( 782) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 782) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 782) =    0.00
  given_sigma_ion_O_4S_cm2( 783) =  3.52 ; given_sigma_ion_O_2D_cm2( 783) =    5.95 ; given_sigma_ion_O_2P_cm2( 783) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 783) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 783) =    0.00
  given_sigma_ion_O_4S_cm2( 784) =  3.56 ; given_sigma_ion_O_2D_cm2( 784) =    5.92 ; given_sigma_ion_O_2P_cm2( 784) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 784) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 784) =    0.00
  given_sigma_ion_O_4S_cm2( 785) =  3.60 ; given_sigma_ion_O_2D_cm2( 785) =    5.88 ; given_sigma_ion_O_2P_cm2( 785) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 785) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 785) =    0.00
  given_sigma_ion_O_4S_cm2( 786) =  3.64 ; given_sigma_ion_O_2D_cm2( 786) =    5.85 ; given_sigma_ion_O_2P_cm2( 786) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 786) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 786) =    0.00
  given_sigma_ion_O_4S_cm2( 787) =  3.68 ; given_sigma_ion_O_2D_cm2( 787) =    5.81 ; given_sigma_ion_O_2P_cm2( 787) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 787) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 787) =    0.00
  given_sigma_ion_O_4S_cm2( 788) =  3.80 ; given_sigma_ion_O_2D_cm2( 788) =    5.70 ; given_sigma_ion_O_2P_cm2( 788) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 788) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 788) =    0.00
  given_sigma_ion_O_4S_cm2( 789) =  3.81 ; given_sigma_ion_O_2D_cm2( 789) =    5.70 ; given_sigma_ion_O_2P_cm2( 789) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 789) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 789) =    0.00
  given_sigma_ion_O_4S_cm2( 790) =  3.85 ; given_sigma_ion_O_2D_cm2( 790) =    5.72 ; given_sigma_ion_O_2P_cm2( 790) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 790) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 790) =    0.00
  given_sigma_ion_O_4S_cm2( 791) =  3.84 ; given_sigma_ion_O_2D_cm2( 791) =    5.67 ; given_sigma_ion_O_2P_cm2( 791) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 791) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 791) =    0.00
  given_sigma_ion_O_4S_cm2( 792) =  3.79 ; given_sigma_ion_O_2D_cm2( 792) =    5.57 ; given_sigma_ion_O_2P_cm2( 792) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 792) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 792) =    0.00
  given_sigma_ion_O_4S_cm2( 793) =  3.68 ; given_sigma_ion_O_2D_cm2( 793) =    5.32 ; given_sigma_ion_O_2P_cm2( 793) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 793) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 793) =    0.00
  given_sigma_ion_O_4S_cm2( 794) =  3.69 ; given_sigma_ion_O_2D_cm2( 794) =    5.31 ; given_sigma_ion_O_2P_cm2( 794) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 794) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 794) =    0.00
  given_sigma_ion_O_4S_cm2( 795) =  3.71 ; given_sigma_ion_O_2D_cm2( 795) =    5.31 ; given_sigma_ion_O_2P_cm2( 795) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 795) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 795) =    0.00
  given_sigma_ion_O_4S_cm2( 796) =  3.75 ; given_sigma_ion_O_2D_cm2( 796) =    5.30 ; given_sigma_ion_O_2P_cm2( 796) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 796) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 796) =    0.00
  given_sigma_ion_O_4S_cm2( 797) =  3.77 ; given_sigma_ion_O_2D_cm2( 797) =    5.29 ; given_sigma_ion_O_2P_cm2( 797) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 797) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 797) =    0.00
  given_sigma_ion_O_4S_cm2( 798) =  3.74 ; given_sigma_ion_O_2D_cm2( 798) =    5.17 ; given_sigma_ion_O_2P_cm2( 798) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 798) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 798) =    0.00
  given_sigma_ion_O_4S_cm2( 799) =  3.96 ; given_sigma_ion_O_2D_cm2( 799) =    4.84 ; given_sigma_ion_O_2P_cm2( 799) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 799) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 799) =    0.00
  given_sigma_ion_O_4S_cm2( 800) =  3.24 ; given_sigma_ion_O_2D_cm2( 800) =    6.12 ; given_sigma_ion_O_2P_cm2( 800) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 800) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 800) =    0.00
  given_sigma_ion_O_4S_cm2( 801) =  2.75 ; given_sigma_ion_O_2D_cm2( 801) =    6.95 ; given_sigma_ion_O_2P_cm2( 801) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 801) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 801) =    0.00
  given_sigma_ion_O_4S_cm2( 802) =  9.67 ; given_sigma_ion_O_2D_cm2( 802) =   30.33 ; given_sigma_ion_O_2P_cm2( 802) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 802) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 802) =    0.00
  given_sigma_ion_O_4S_cm2( 803) =  7.42 ; given_sigma_ion_O_2D_cm2( 803) =   29.67 ; given_sigma_ion_O_2P_cm2( 803) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 803) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 803) =    0.00
  given_sigma_ion_O_4S_cm2( 804) =  7.65 ; given_sigma_ion_O_2D_cm2( 804) =   25.06 ; given_sigma_ion_O_2P_cm2( 804) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 804) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 804) =    0.00
  given_sigma_ion_O_4S_cm2( 805) =  7.59 ; given_sigma_ion_O_2D_cm2( 805) =   20.75 ; given_sigma_ion_O_2P_cm2( 805) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 805) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 805) =    0.00
  given_sigma_ion_O_4S_cm2( 806) =  6.58 ; given_sigma_ion_O_2D_cm2( 806) =   13.02 ; given_sigma_ion_O_2P_cm2( 806) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 806) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 806) =    0.00
  given_sigma_ion_O_4S_cm2( 807) =  6.34 ; given_sigma_ion_O_2D_cm2( 807) =    9.83 ; given_sigma_ion_O_2P_cm2( 807) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 807) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 807) =    0.00
  given_sigma_ion_O_4S_cm2( 808) =  6.14 ; given_sigma_ion_O_2D_cm2( 808) =    8.66 ; given_sigma_ion_O_2P_cm2( 808) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 808) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 808) =    0.00
  given_sigma_ion_O_4S_cm2( 809) =  6.01 ; given_sigma_ion_O_2D_cm2( 809) =    8.10 ; given_sigma_ion_O_2P_cm2( 809) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 809) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 809) =    0.00
  given_sigma_ion_O_4S_cm2( 810) =  5.54 ; given_sigma_ion_O_2D_cm2( 810) =    6.51 ; given_sigma_ion_O_2P_cm2( 810) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 810) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 810) =    0.00
  given_sigma_ion_O_4S_cm2( 811) =  4.56 ; given_sigma_ion_O_2D_cm2( 811) =    5.42 ; given_sigma_ion_O_2P_cm2( 811) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 811) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 811) =    0.00
  given_sigma_ion_O_4S_cm2( 812) =  4.24 ; given_sigma_ion_O_2D_cm2( 812) =    5.06 ; given_sigma_ion_O_2P_cm2( 812) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 812) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 812) =    0.00
  given_sigma_ion_O_4S_cm2( 813) =  3.82 ; given_sigma_ion_O_2D_cm2( 813) =    4.63 ; given_sigma_ion_O_2P_cm2( 813) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 813) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 813) =    0.00
  given_sigma_ion_O_4S_cm2( 814) =  3.61 ; given_sigma_ion_O_2D_cm2( 814) =    4.41 ; given_sigma_ion_O_2P_cm2( 814) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 814) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 814) =    0.00
  given_sigma_ion_O_4S_cm2( 815) =  3.40 ; given_sigma_ion_O_2D_cm2( 815) =    4.20 ; given_sigma_ion_O_2P_cm2( 815) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 815) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 815) =    0.00
  given_sigma_ion_O_4S_cm2( 816) =  5.24 ; given_sigma_ion_O_2D_cm2( 816) =    6.56 ; given_sigma_ion_O_2P_cm2( 816) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 816) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 816) =    0.00
  given_sigma_ion_O_4S_cm2( 817) =  6.14 ; given_sigma_ion_O_2D_cm2( 817) =    7.76 ; given_sigma_ion_O_2P_cm2( 817) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 817) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 817) =    0.00
  given_sigma_ion_O_4S_cm2( 818) =  5.98 ; given_sigma_ion_O_2D_cm2( 818) =    7.61 ; given_sigma_ion_O_2P_cm2( 818) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 818) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 818) =    0.00
  given_sigma_ion_O_4S_cm2( 819) =  5.17 ; given_sigma_ion_O_2D_cm2( 819) =    6.85 ; given_sigma_ion_O_2P_cm2( 819) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 819) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 819) =    0.00
  given_sigma_ion_O_4S_cm2( 820) =  4.89 ; given_sigma_ion_O_2D_cm2( 820) =    6.57 ; given_sigma_ion_O_2P_cm2( 820) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 820) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 820) =    0.00
  given_sigma_ion_O_4S_cm2( 821) =  4.70 ; given_sigma_ion_O_2D_cm2( 821) =    6.38 ; given_sigma_ion_O_2P_cm2( 821) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 821) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 821) =    0.00
  given_sigma_ion_O_4S_cm2( 822) =  4.39 ; given_sigma_ion_O_2D_cm2( 822) =    6.06 ; given_sigma_ion_O_2P_cm2( 822) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 822) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 822) =    0.00
  given_sigma_ion_O_4S_cm2( 823) =  4.15 ; given_sigma_ion_O_2D_cm2( 823) =    5.84 ; given_sigma_ion_O_2P_cm2( 823) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 823) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 823) =    0.00
  given_sigma_ion_O_4S_cm2( 824) =  3.94 ; given_sigma_ion_O_2D_cm2( 824) =    5.66 ; given_sigma_ion_O_2P_cm2( 824) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 824) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 824) =    0.00
  given_sigma_ion_O_4S_cm2( 825) =  3.96 ; given_sigma_ion_O_2D_cm2( 825) =    5.67 ; given_sigma_ion_O_2P_cm2( 825) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 825) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 825) =    0.00
  given_sigma_ion_O_4S_cm2( 826) =  4.02 ; given_sigma_ion_O_2D_cm2( 826) =    5.69 ; given_sigma_ion_O_2P_cm2( 826) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 826) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 826) =    0.00
  given_sigma_ion_O_4S_cm2( 827) =  4.07 ; given_sigma_ion_O_2D_cm2( 827) =    5.71 ; given_sigma_ion_O_2P_cm2( 827) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 827) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 827) =    0.00
  given_sigma_ion_O_4S_cm2( 828) =  4.10 ; given_sigma_ion_O_2D_cm2( 828) =    5.72 ; given_sigma_ion_O_2P_cm2( 828) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 828) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 828) =    0.00
  given_sigma_ion_O_4S_cm2( 829) =  4.15 ; given_sigma_ion_O_2D_cm2( 829) =    5.73 ; given_sigma_ion_O_2P_cm2( 829) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 829) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 829) =    0.00
  given_sigma_ion_O_4S_cm2( 830) =  4.15 ; given_sigma_ion_O_2D_cm2( 830) =    5.70 ; given_sigma_ion_O_2P_cm2( 830) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 830) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 830) =    0.00
  given_sigma_ion_O_4S_cm2( 831) =  4.15 ; given_sigma_ion_O_2D_cm2( 831) =    5.60 ; given_sigma_ion_O_2P_cm2( 831) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 831) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 831) =    0.00
  given_sigma_ion_O_4S_cm2( 832) =  4.14 ; given_sigma_ion_O_2D_cm2( 832) =    5.57 ; given_sigma_ion_O_2P_cm2( 832) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 832) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 832) =    0.00
  given_sigma_ion_O_4S_cm2( 833) =  4.01 ; given_sigma_ion_O_2D_cm2( 833) =    5.35 ; given_sigma_ion_O_2P_cm2( 833) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 833) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 833) =    0.00
  given_sigma_ion_O_4S_cm2( 834) =  3.93 ; given_sigma_ion_O_2D_cm2( 834) =    5.22 ; given_sigma_ion_O_2P_cm2( 834) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 834) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 834) =    0.00
  given_sigma_ion_O_4S_cm2( 835) =  4.03 ; given_sigma_ion_O_2D_cm2( 835) =    5.27 ; given_sigma_ion_O_2P_cm2( 835) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 835) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 835) =    0.00
  given_sigma_ion_O_4S_cm2( 836) =  4.03 ; given_sigma_ion_O_2D_cm2( 836) =    5.26 ; given_sigma_ion_O_2P_cm2( 836) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 836) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 836) =    0.00
  given_sigma_ion_O_4S_cm2( 837) =  4.06 ; given_sigma_ion_O_2D_cm2( 837) =    5.22 ; given_sigma_ion_O_2P_cm2( 837) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 837) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 837) =    0.00
  given_sigma_ion_O_4S_cm2( 838) =  4.07 ; given_sigma_ion_O_2D_cm2( 838) =    5.19 ; given_sigma_ion_O_2P_cm2( 838) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 838) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 838) =    0.00
  given_sigma_ion_O_4S_cm2( 839) =  4.08 ; given_sigma_ion_O_2D_cm2( 839) =    5.17 ; given_sigma_ion_O_2P_cm2( 839) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 839) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 839) =    0.00
  given_sigma_ion_O_4S_cm2( 840) =  4.08 ; given_sigma_ion_O_2D_cm2( 840) =    5.16 ; given_sigma_ion_O_2P_cm2( 840) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 840) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 840) =    0.00
  given_sigma_ion_O_4S_cm2( 841) =  4.09 ; given_sigma_ion_O_2D_cm2( 841) =    5.12 ; given_sigma_ion_O_2P_cm2( 841) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 841) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 841) =    0.00
  given_sigma_ion_O_4S_cm2( 842) =  4.10 ; given_sigma_ion_O_2D_cm2( 842) =    5.09 ; given_sigma_ion_O_2P_cm2( 842) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 842) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 842) =    0.00
  given_sigma_ion_O_4S_cm2( 843) =  4.13 ; given_sigma_ion_O_2D_cm2( 843) =    5.07 ; given_sigma_ion_O_2P_cm2( 843) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 843) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 843) =    0.00
  given_sigma_ion_O_4S_cm2( 844) =  4.15 ; given_sigma_ion_O_2D_cm2( 844) =    5.07 ; given_sigma_ion_O_2P_cm2( 844) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 844) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 844) =    0.00
  given_sigma_ion_O_4S_cm2( 845) =  4.18 ; given_sigma_ion_O_2D_cm2( 845) =    5.06 ; given_sigma_ion_O_2P_cm2( 845) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 845) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 845) =    0.00
  given_sigma_ion_O_4S_cm2( 846) =  4.18 ; given_sigma_ion_O_2D_cm2( 846) =    5.06 ; given_sigma_ion_O_2P_cm2( 846) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 846) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 846) =    0.00
  given_sigma_ion_O_4S_cm2( 847) =  4.16 ; given_sigma_ion_O_2D_cm2( 847) =    5.01 ; given_sigma_ion_O_2P_cm2( 847) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 847) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 847) =    0.00
  given_sigma_ion_O_4S_cm2( 848) =  4.13 ; given_sigma_ion_O_2D_cm2( 848) =    4.92 ; given_sigma_ion_O_2P_cm2( 848) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 848) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 848) =    0.00
  given_sigma_ion_O_4S_cm2( 849) =  4.09 ; given_sigma_ion_O_2D_cm2( 849) =    4.84 ; given_sigma_ion_O_2P_cm2( 849) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 849) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 849) =    0.00
  given_sigma_ion_O_4S_cm2( 850) =  4.08 ; given_sigma_ion_O_2D_cm2( 850) =    4.81 ; given_sigma_ion_O_2P_cm2( 850) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 850) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 850) =    0.00
  given_sigma_ion_O_4S_cm2( 851) =  4.05 ; given_sigma_ion_O_2D_cm2( 851) =    4.74 ; given_sigma_ion_O_2P_cm2( 851) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 851) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 851) =    0.00
  given_sigma_ion_O_4S_cm2( 852) =  4.09 ; given_sigma_ion_O_2D_cm2( 852) =    4.76 ; given_sigma_ion_O_2P_cm2( 852) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 852) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 852) =    0.00
  given_sigma_ion_O_4S_cm2( 853) =  4.13 ; given_sigma_ion_O_2D_cm2( 853) =    4.77 ; given_sigma_ion_O_2P_cm2( 853) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 853) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 853) =    0.00
  given_sigma_ion_O_4S_cm2( 854) =  4.18 ; given_sigma_ion_O_2D_cm2( 854) =    4.79 ; given_sigma_ion_O_2P_cm2( 854) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 854) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 854) =    0.00
  given_sigma_ion_O_4S_cm2( 855) =  4.20 ; given_sigma_ion_O_2D_cm2( 855) =    4.80 ; given_sigma_ion_O_2P_cm2( 855) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 855) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 855) =    0.00
  given_sigma_ion_O_4S_cm2( 856) =  4.21 ; given_sigma_ion_O_2D_cm2( 856) =    4.79 ; given_sigma_ion_O_2P_cm2( 856) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 856) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 856) =    0.00
  given_sigma_ion_O_4S_cm2( 857) =  4.23 ; given_sigma_ion_O_2D_cm2( 857) =    4.77 ; given_sigma_ion_O_2P_cm2( 857) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 857) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 857) =    0.00
  given_sigma_ion_O_4S_cm2( 858) =  4.10 ; given_sigma_ion_O_2D_cm2( 858) =    4.86 ; given_sigma_ion_O_2P_cm2( 858) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 858) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 858) =    0.00
  given_sigma_ion_O_4S_cm2( 859) =  3.95 ; given_sigma_ion_O_2D_cm2( 859) =    4.96 ; given_sigma_ion_O_2P_cm2( 859) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 859) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 859) =    0.00
  given_sigma_ion_O_4S_cm2( 860) =  3.85 ; given_sigma_ion_O_2D_cm2( 860) =    5.09 ; given_sigma_ion_O_2P_cm2( 860) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 860) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 860) =    0.00
  given_sigma_ion_O_4S_cm2( 861) =  3.83 ; given_sigma_ion_O_2D_cm2( 861) =    5.12 ; given_sigma_ion_O_2P_cm2( 861) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 861) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 861) =    0.00
  given_sigma_ion_O_4S_cm2( 862) =  3.73 ; given_sigma_ion_O_2D_cm2( 862) =    5.24 ; given_sigma_ion_O_2P_cm2( 862) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 862) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 862) =    0.00
  given_sigma_ion_O_4S_cm2( 863) =  3.68 ; given_sigma_ion_O_2D_cm2( 863) =    5.30 ; given_sigma_ion_O_2P_cm2( 863) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 863) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 863) =    0.00
  given_sigma_ion_O_4S_cm2( 864) =  3.56 ; given_sigma_ion_O_2D_cm2( 864) =    5.46 ; given_sigma_ion_O_2P_cm2( 864) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 864) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 864) =    0.00
  given_sigma_ion_O_4S_cm2( 865) =  3.41 ; given_sigma_ion_O_2D_cm2( 865) =    5.49 ; given_sigma_ion_O_2P_cm2( 865) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 865) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 865) =    0.00
  given_sigma_ion_O_4S_cm2( 866) =  3.61 ; given_sigma_ion_O_2D_cm2( 866) =    5.14 ; given_sigma_ion_O_2P_cm2( 866) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 866) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 866) =    0.00
  given_sigma_ion_O_4S_cm2( 867) =  3.68 ; given_sigma_ion_O_2D_cm2( 867) =    4.92 ; given_sigma_ion_O_2P_cm2( 867) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 867) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 867) =    0.00
  given_sigma_ion_O_4S_cm2( 868) =  9.72 ; given_sigma_ion_O_2D_cm2( 868) =   16.39 ; given_sigma_ion_O_2P_cm2( 868) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 868) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 868) =    0.00
  given_sigma_ion_O_4S_cm2( 869) = 12.23 ; given_sigma_ion_O_2D_cm2( 869) =   49.27 ; given_sigma_ion_O_2P_cm2( 869) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 869) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 869) =    0.00
  given_sigma_ion_O_4S_cm2( 870) =  9.23 ; given_sigma_ion_O_2D_cm2( 870) =   45.07 ; given_sigma_ion_O_2P_cm2( 870) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 870) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 870) =    0.00
  given_sigma_ion_O_4S_cm2( 871) =  6.59 ; given_sigma_ion_O_2D_cm2( 871) =   27.31 ; given_sigma_ion_O_2P_cm2( 871) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 871) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 871) =    0.00
  given_sigma_ion_O_4S_cm2( 872) =  5.48 ; given_sigma_ion_O_2D_cm2( 872) =   21.62 ; given_sigma_ion_O_2P_cm2( 872) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 872) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 872) =    0.00
  given_sigma_ion_O_4S_cm2( 873) =  5.16 ; given_sigma_ion_O_2D_cm2( 873) =   18.90 ; given_sigma_ion_O_2P_cm2( 873) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 873) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 873) =    0.00
  given_sigma_ion_O_4S_cm2( 874) =  4.62 ; given_sigma_ion_O_2D_cm2( 874) =   15.40 ; given_sigma_ion_O_2P_cm2( 874) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 874) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 874) =    0.00
  given_sigma_ion_O_4S_cm2( 875) =  3.56 ; given_sigma_ion_O_2D_cm2( 875) =   10.40 ; given_sigma_ion_O_2P_cm2( 875) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 875) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 875) =    0.00
  given_sigma_ion_O_4S_cm2( 876) =  2.92 ; given_sigma_ion_O_2D_cm2( 876) =    8.00 ; given_sigma_ion_O_2P_cm2( 876) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 876) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 876) =    0.00
  given_sigma_ion_O_4S_cm2( 877) =  2.45 ; given_sigma_ion_O_2D_cm2( 877) =    6.45 ; given_sigma_ion_O_2P_cm2( 877) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 877) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 877) =    0.00
  given_sigma_ion_O_4S_cm2( 878) =  2.66 ; given_sigma_ion_O_2D_cm2( 878) =    6.33 ; given_sigma_ion_O_2P_cm2( 878) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 878) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 878) =    0.00
  given_sigma_ion_O_4S_cm2( 879) =  2.78 ; given_sigma_ion_O_2D_cm2( 879) =    6.26 ; given_sigma_ion_O_2P_cm2( 879) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 879) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 879) =    0.00
  given_sigma_ion_O_4S_cm2( 880) =  2.95 ; given_sigma_ion_O_2D_cm2( 880) =    6.16 ; given_sigma_ion_O_2P_cm2( 880) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 880) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 880) =    0.00
  given_sigma_ion_O_4S_cm2( 881) =  3.08 ; given_sigma_ion_O_2D_cm2( 881) =    6.09 ; given_sigma_ion_O_2P_cm2( 881) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 881) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 881) =    0.00
  given_sigma_ion_O_4S_cm2( 882) =  3.16 ; given_sigma_ion_O_2D_cm2( 882) =    5.82 ; given_sigma_ion_O_2P_cm2( 882) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 882) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 882) =    0.00
  given_sigma_ion_O_4S_cm2( 883) =  3.24 ; given_sigma_ion_O_2D_cm2( 883) =    5.55 ; given_sigma_ion_O_2P_cm2( 883) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 883) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 883) =    0.00
  given_sigma_ion_O_4S_cm2( 884) =  3.27 ; given_sigma_ion_O_2D_cm2( 884) =    5.42 ; given_sigma_ion_O_2P_cm2( 884) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 884) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 884) =    0.00
  given_sigma_ion_O_4S_cm2( 885) =  3.37 ; given_sigma_ion_O_2D_cm2( 885) =    5.03 ; given_sigma_ion_O_2P_cm2( 885) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 885) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 885) =    0.00
  given_sigma_ion_O_4S_cm2( 886) =  3.44 ; given_sigma_ion_O_2D_cm2( 886) =    4.89 ; given_sigma_ion_O_2P_cm2( 886) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 886) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 886) =    0.00
  given_sigma_ion_O_4S_cm2( 887) =  3.58 ; given_sigma_ion_O_2D_cm2( 887) =    4.61 ; given_sigma_ion_O_2P_cm2( 887) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 887) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 887) =    0.00
  given_sigma_ion_O_4S_cm2( 888) =  3.65 ; given_sigma_ion_O_2D_cm2( 888) =    4.47 ; given_sigma_ion_O_2P_cm2( 888) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 888) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 888) =    0.00
  given_sigma_ion_O_4S_cm2( 889) =  8.08 ; given_sigma_ion_O_2D_cm2( 889) =    0.00 ; given_sigma_ion_O_2P_cm2( 889) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 889) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 889) =    0.00
  given_sigma_ion_O_4S_cm2( 890) =  8.03 ; given_sigma_ion_O_2D_cm2( 890) =    0.00 ; given_sigma_ion_O_2P_cm2( 890) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 890) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 890) =    0.00
  given_sigma_ion_O_4S_cm2( 891) =  6.97 ; given_sigma_ion_O_2D_cm2( 891) =    0.00 ; given_sigma_ion_O_2P_cm2( 891) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 891) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 891) =    0.00
  given_sigma_ion_O_4S_cm2( 892) =  5.21 ; given_sigma_ion_O_2D_cm2( 892) =    0.00 ; given_sigma_ion_O_2P_cm2( 892) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 892) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 892) =    0.00
  given_sigma_ion_O_4S_cm2( 893) =  4.15 ; given_sigma_ion_O_2D_cm2( 893) =    0.00 ; given_sigma_ion_O_2P_cm2( 893) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 893) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 893) =    0.00
  given_sigma_ion_O_4S_cm2( 894) =  8.62 ; given_sigma_ion_O_2D_cm2( 894) =    0.00 ; given_sigma_ion_O_2P_cm2( 894) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 894) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 894) =    0.00
  given_sigma_ion_O_4S_cm2( 895) = 11.81 ; given_sigma_ion_O_2D_cm2( 895) =    0.00 ; given_sigma_ion_O_2P_cm2( 895) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 895) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 895) =    0.00
  given_sigma_ion_O_4S_cm2( 896) = 15.00 ; given_sigma_ion_O_2D_cm2( 896) =    0.00 ; given_sigma_ion_O_2P_cm2( 896) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 896) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 896) =    0.00
  given_sigma_ion_O_4S_cm2( 897) = 15.34 ; given_sigma_ion_O_2D_cm2( 897) =    0.00 ; given_sigma_ion_O_2P_cm2( 897) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 897) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 897) =    0.00
  given_sigma_ion_O_4S_cm2( 898) = 16.01 ; given_sigma_ion_O_2D_cm2( 898) =    0.00 ; given_sigma_ion_O_2P_cm2( 898) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 898) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 898) =    0.00
  given_sigma_ion_O_4S_cm2( 899) = 16.69 ; given_sigma_ion_O_2D_cm2( 899) =    0.00 ; given_sigma_ion_O_2P_cm2( 899) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 899) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 899) =    0.00
  given_sigma_ion_O_4S_cm2( 900) = 10.40 ; given_sigma_ion_O_2D_cm2( 900) =    0.00 ; given_sigma_ion_O_2P_cm2( 900) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 900) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 900) =    0.00
  given_sigma_ion_O_4S_cm2( 901) = 12.30 ; given_sigma_ion_O_2D_cm2( 901) =    0.00 ; given_sigma_ion_O_2P_cm2( 901) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 901) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 901) =    0.00
  given_sigma_ion_O_4S_cm2( 902) = 15.15 ; given_sigma_ion_O_2D_cm2( 902) =    0.00 ; given_sigma_ion_O_2P_cm2( 902) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 902) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 902) =    0.00
  given_sigma_ion_O_4S_cm2( 903) = 12.50 ; given_sigma_ion_O_2D_cm2( 903) =    0.00 ; given_sigma_ion_O_2P_cm2( 903) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 903) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 903) =    0.00
  given_sigma_ion_O_4S_cm2( 904) =  8.90 ; given_sigma_ion_O_2D_cm2( 904) =    0.00 ; given_sigma_ion_O_2P_cm2( 904) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 904) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 904) =    0.00
  given_sigma_ion_O_4S_cm2( 905) = 14.37 ; given_sigma_ion_O_2D_cm2( 905) =    0.00 ; given_sigma_ion_O_2P_cm2( 905) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 905) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 905) =    0.00
  given_sigma_ion_O_4S_cm2( 906) = 16.20 ; given_sigma_ion_O_2D_cm2( 906) =    0.00 ; given_sigma_ion_O_2P_cm2( 906) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 906) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 906) =    0.00
  given_sigma_ion_O_4S_cm2( 907) =  9.08 ; given_sigma_ion_O_2D_cm2( 907) =    0.00 ; given_sigma_ion_O_2P_cm2( 907) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 907) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 907) =    0.00
  given_sigma_ion_O_4S_cm2( 908) =  7.30 ; given_sigma_ion_O_2D_cm2( 908) =    0.00 ; given_sigma_ion_O_2P_cm2( 908) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 908) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 908) =    0.00
  given_sigma_ion_O_4S_cm2( 909) = 18.41 ; given_sigma_ion_O_2D_cm2( 909) =    0.00 ; given_sigma_ion_O_2P_cm2( 909) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 909) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 909) =    0.00
  given_sigma_ion_O_4S_cm2( 910) = 18.91 ; given_sigma_ion_O_2D_cm2( 910) =    0.00 ; given_sigma_ion_O_2P_cm2( 910) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 910) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 910) =    0.00
  given_sigma_ion_O_4S_cm2( 911) = 11.77 ; given_sigma_ion_O_2D_cm2( 911) =    0.00 ; given_sigma_ion_O_2P_cm2( 911) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 911) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 911) =    0.00
  given_sigma_ion_O_4S_cm2( 912) =  9.98 ; given_sigma_ion_O_2D_cm2( 912) =    0.00 ; given_sigma_ion_O_2P_cm2( 912) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 912) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 912) =    0.00
  given_sigma_ion_O_4S_cm2( 913) =  8.20 ; given_sigma_ion_O_2D_cm2( 913) =    0.00 ; given_sigma_ion_O_2P_cm2( 913) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 913) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 913) =    0.00
  given_sigma_ion_O_4S_cm2( 914) =  6.33 ; given_sigma_ion_O_2D_cm2( 914) =    0.00 ; given_sigma_ion_O_2P_cm2( 914) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 914) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 914) =    0.00
  given_sigma_ion_O_4S_cm2( 915) =  5.71 ; given_sigma_ion_O_2D_cm2( 915) =    0.00 ; given_sigma_ion_O_2P_cm2( 915) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 915) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 915) =    0.00
  given_sigma_ion_O_4S_cm2( 916) = 10.08 ; given_sigma_ion_O_2D_cm2( 916) =    0.00 ; given_sigma_ion_O_2P_cm2( 916) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 916) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 916) =    0.00
  given_sigma_ion_O_4S_cm2( 917) = 14.75 ; given_sigma_ion_O_2D_cm2( 917) =    0.00 ; given_sigma_ion_O_2P_cm2( 917) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 917) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 917) =    0.00
  given_sigma_ion_O_4S_cm2( 918) = 21.76 ; given_sigma_ion_O_2D_cm2( 918) =    0.00 ; given_sigma_ion_O_2P_cm2( 918) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 918) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 918) =    0.00
  given_sigma_ion_O_4S_cm2( 919) = 17.10 ; given_sigma_ion_O_2D_cm2( 919) =    0.00 ; given_sigma_ion_O_2P_cm2( 919) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 919) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 919) =    0.00
  given_sigma_ion_O_4S_cm2( 920) = 10.10 ; given_sigma_ion_O_2D_cm2( 920) =    0.00 ; given_sigma_ion_O_2P_cm2( 920) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 920) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 920) =    0.00
  given_sigma_ion_O_4S_cm2( 921) =  5.42 ; given_sigma_ion_O_2D_cm2( 921) =    0.00 ; given_sigma_ion_O_2P_cm2( 921) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 921) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 921) =    0.00
  given_sigma_ion_O_4S_cm2( 922) =  4.85 ; given_sigma_ion_O_2D_cm2( 922) =    0.00 ; given_sigma_ion_O_2P_cm2( 922) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 922) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 922) =    0.00
  given_sigma_ion_O_4S_cm2( 923) =  4.70 ; given_sigma_ion_O_2D_cm2( 923) =    0.00 ; given_sigma_ion_O_2P_cm2( 923) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 923) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 923) =    0.00
  given_sigma_ion_O_4S_cm2( 924) =  5.49 ; given_sigma_ion_O_2D_cm2( 924) =    0.00 ; given_sigma_ion_O_2P_cm2( 924) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 924) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 924) =    0.00
  given_sigma_ion_O_4S_cm2( 925) =  6.80 ; given_sigma_ion_O_2D_cm2( 925) =    0.00 ; given_sigma_ion_O_2P_cm2( 925) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 925) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 925) =    0.00
  given_sigma_ion_O_4S_cm2( 926) = 20.10 ; given_sigma_ion_O_2D_cm2( 926) =    0.00 ; given_sigma_ion_O_2P_cm2( 926) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 926) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 926) =    0.00
  given_sigma_ion_O_4S_cm2( 927) = 33.40 ; given_sigma_ion_O_2D_cm2( 927) =    0.00 ; given_sigma_ion_O_2P_cm2( 927) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 927) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 927) =    0.00
  given_sigma_ion_O_4S_cm2( 928) = 20.33 ; given_sigma_ion_O_2D_cm2( 928) =    0.00 ; given_sigma_ion_O_2P_cm2( 928) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 928) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 928) =    0.00
  given_sigma_ion_O_4S_cm2( 929) = 15.10 ; given_sigma_ion_O_2D_cm2( 929) =    0.00 ; given_sigma_ion_O_2P_cm2( 929) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 929) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 929) =    0.00
  given_sigma_ion_O_4S_cm2( 930) =  9.82 ; given_sigma_ion_O_2D_cm2( 930) =    0.00 ; given_sigma_ion_O_2P_cm2( 930) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 930) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 930) =    0.00
  given_sigma_ion_O_4S_cm2( 931) =  5.20 ; given_sigma_ion_O_2D_cm2( 931) =    0.00 ; given_sigma_ion_O_2P_cm2( 931) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 931) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 931) =    0.00
  given_sigma_ion_O_4S_cm2( 932) =  5.00 ; given_sigma_ion_O_2D_cm2( 932) =    0.00 ; given_sigma_ion_O_2P_cm2( 932) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 932) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 932) =    0.00
  given_sigma_ion_O_4S_cm2( 933) =  4.61 ; given_sigma_ion_O_2D_cm2( 933) =    0.00 ; given_sigma_ion_O_2P_cm2( 933) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 933) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 933) =    0.00
  given_sigma_ion_O_4S_cm2( 934) =  4.35 ; given_sigma_ion_O_2D_cm2( 934) =    0.00 ; given_sigma_ion_O_2P_cm2( 934) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 934) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 934) =    0.00
  given_sigma_ion_O_4S_cm2( 935) =  4.15 ; given_sigma_ion_O_2D_cm2( 935) =    0.00 ; given_sigma_ion_O_2P_cm2( 935) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 935) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 935) =    0.00
  given_sigma_ion_O_4S_cm2( 936) =  4.36 ; given_sigma_ion_O_2D_cm2( 936) =    0.00 ; given_sigma_ion_O_2P_cm2( 936) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 936) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 936) =    0.00
  given_sigma_ion_O_4S_cm2( 937) =  4.76 ; given_sigma_ion_O_2D_cm2( 937) =    0.00 ; given_sigma_ion_O_2P_cm2( 937) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 937) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 937) =    0.00
  given_sigma_ion_O_4S_cm2( 938) =  4.90 ; given_sigma_ion_O_2D_cm2( 938) =    0.00 ; given_sigma_ion_O_2P_cm2( 938) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 938) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 938) =    0.00
  given_sigma_ion_O_4S_cm2( 939) = 22.34 ; given_sigma_ion_O_2D_cm2( 939) =    0.00 ; given_sigma_ion_O_2P_cm2( 939) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 939) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 939) =    0.00
  given_sigma_ion_O_4S_cm2( 940) = 28.16 ; given_sigma_ion_O_2D_cm2( 940) =    0.00 ; given_sigma_ion_O_2P_cm2( 940) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 940) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 940) =    0.00
  given_sigma_ion_O_4S_cm2( 941) = 36.88 ; given_sigma_ion_O_2D_cm2( 941) =    0.00 ; given_sigma_ion_O_2P_cm2( 941) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 941) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 941) =    0.00
  given_sigma_ion_O_4S_cm2( 942) = 45.60 ; given_sigma_ion_O_2D_cm2( 942) =    0.00 ; given_sigma_ion_O_2P_cm2( 942) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 942) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 942) =    0.00
  given_sigma_ion_O_4S_cm2( 943) = 40.19 ; given_sigma_ion_O_2D_cm2( 943) =    0.00 ; given_sigma_ion_O_2P_cm2( 943) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 943) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 943) =    0.00
  given_sigma_ion_O_4S_cm2( 944) = 34.77 ; given_sigma_ion_O_2D_cm2( 944) =    0.00 ; given_sigma_ion_O_2P_cm2( 944) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 944) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 944) =    0.00
  given_sigma_ion_O_4S_cm2( 945) = 26.65 ; given_sigma_ion_O_2D_cm2( 945) =    0.00 ; given_sigma_ion_O_2P_cm2( 945) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 945) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 945) =    0.00
  given_sigma_ion_O_4S_cm2( 946) = 13.11 ; given_sigma_ion_O_2D_cm2( 946) =    0.00 ; given_sigma_ion_O_2P_cm2( 946) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 946) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 946) =    0.00
  given_sigma_ion_O_4S_cm2( 947) =  7.70 ; given_sigma_ion_O_2D_cm2( 947) =    0.00 ; given_sigma_ion_O_2P_cm2( 947) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 947) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 947) =    0.00
  given_sigma_ion_O_4S_cm2( 948) =  9.57 ; given_sigma_ion_O_2D_cm2( 948) =    0.00 ; given_sigma_ion_O_2P_cm2( 948) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 948) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 948) =    0.00
  given_sigma_ion_O_4S_cm2( 949) = 10.83 ; given_sigma_ion_O_2D_cm2( 949) =    0.00 ; given_sigma_ion_O_2P_cm2( 949) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 949) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 949) =    0.00
  given_sigma_ion_O_4S_cm2( 950) = 12.70 ; given_sigma_ion_O_2D_cm2( 950) =    0.00 ; given_sigma_ion_O_2P_cm2( 950) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 950) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 950) =    0.00
  given_sigma_ion_O_4S_cm2( 951) = 11.38 ; given_sigma_ion_O_2D_cm2( 951) =    0.00 ; given_sigma_ion_O_2P_cm2( 951) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 951) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 951) =    0.00
  given_sigma_ion_O_4S_cm2( 952) = 10.94 ; given_sigma_ion_O_2D_cm2( 952) =    0.00 ; given_sigma_ion_O_2P_cm2( 952) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 952) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 952) =    0.00
  given_sigma_ion_O_4S_cm2( 953) =  9.71 ; given_sigma_ion_O_2D_cm2( 953) =    0.00 ; given_sigma_ion_O_2P_cm2( 953) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 953) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 953) =    0.00
  given_sigma_ion_O_4S_cm2( 954) =  8.30 ; given_sigma_ion_O_2D_cm2( 954) =    0.00 ; given_sigma_ion_O_2P_cm2( 954) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 954) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 954) =    0.00
  given_sigma_ion_O_4S_cm2( 955) =  6.54 ; given_sigma_ion_O_2D_cm2( 955) =    0.00 ; given_sigma_ion_O_2P_cm2( 955) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 955) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 955) =    0.00
  given_sigma_ion_O_4S_cm2( 956) =  3.90 ; given_sigma_ion_O_2D_cm2( 956) =    0.00 ; given_sigma_ion_O_2P_cm2( 956) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 956) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 956) =    0.00
  given_sigma_ion_O_4S_cm2( 957) =  4.05 ; given_sigma_ion_O_2D_cm2( 957) =    0.00 ; given_sigma_ion_O_2P_cm2( 957) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 957) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 957) =    0.00
  given_sigma_ion_O_4S_cm2( 958) =  4.06 ; given_sigma_ion_O_2D_cm2( 958) =    0.00 ; given_sigma_ion_O_2P_cm2( 958) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 958) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 958) =    0.00
  given_sigma_ion_O_4S_cm2( 959) =  4.05 ; given_sigma_ion_O_2D_cm2( 959) =    0.00 ; given_sigma_ion_O_2P_cm2( 959) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 959) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 959) =    0.00
  given_sigma_ion_O_4S_cm2( 960) =  4.05 ; given_sigma_ion_O_2D_cm2( 960) =    0.00 ; given_sigma_ion_O_2P_cm2( 960) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 960) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 960) =    0.00
  given_sigma_ion_O_4S_cm2( 961) =  4.04 ; given_sigma_ion_O_2D_cm2( 961) =    0.00 ; given_sigma_ion_O_2P_cm2( 961) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 961) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 961) =    0.00
  given_sigma_ion_O_4S_cm2( 962) =  4.04 ; given_sigma_ion_O_2D_cm2( 962) =    0.00 ; given_sigma_ion_O_2P_cm2( 962) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 962) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 962) =    0.00
  given_sigma_ion_O_4S_cm2( 963) =  4.04 ; given_sigma_ion_O_2D_cm2( 963) =    0.00 ; given_sigma_ion_O_2P_cm2( 963) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 963) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 963) =    0.00
  given_sigma_ion_O_4S_cm2( 964) =  4.03 ; given_sigma_ion_O_2D_cm2( 964) =    0.00 ; given_sigma_ion_O_2P_cm2( 964) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 964) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 964) =    0.00
  given_sigma_ion_O_4S_cm2( 965) =  4.03 ; given_sigma_ion_O_2D_cm2( 965) =    0.00 ; given_sigma_ion_O_2P_cm2( 965) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 965) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 965) =    0.00
  given_sigma_ion_O_4S_cm2( 966) =  4.02 ; given_sigma_ion_O_2D_cm2( 966) =    0.00 ; given_sigma_ion_O_2P_cm2( 966) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 966) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 966) =    0.00
  given_sigma_ion_O_4S_cm2( 967) =  4.02 ; given_sigma_ion_O_2D_cm2( 967) =    0.00 ; given_sigma_ion_O_2P_cm2( 967) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 967) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 967) =    0.00
  given_sigma_ion_O_4S_cm2( 968) =  4.02 ; given_sigma_ion_O_2D_cm2( 968) =    0.00 ; given_sigma_ion_O_2P_cm2( 968) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 968) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 968) =    0.00
  given_sigma_ion_O_4S_cm2( 969) =  4.01 ; given_sigma_ion_O_2D_cm2( 969) =    0.00 ; given_sigma_ion_O_2P_cm2( 969) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 969) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 969) =    0.00
  given_sigma_ion_O_4S_cm2( 970) =  4.01 ; given_sigma_ion_O_2D_cm2( 970) =    0.00 ; given_sigma_ion_O_2P_cm2( 970) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 970) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 970) =    0.00
  given_sigma_ion_O_4S_cm2( 971) =  4.00 ; given_sigma_ion_O_2D_cm2( 971) =    0.00 ; given_sigma_ion_O_2P_cm2( 971) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 971) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 971) =    0.00
  given_sigma_ion_O_4S_cm2( 972) =  4.00 ; given_sigma_ion_O_2D_cm2( 972) =    0.00 ; given_sigma_ion_O_2P_cm2( 972) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 972) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 972) =    0.00
  given_sigma_ion_O_4S_cm2( 973) =  3.99 ; given_sigma_ion_O_2D_cm2( 973) =    0.00 ; given_sigma_ion_O_2P_cm2( 973) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 973) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 973) =    0.00
  given_sigma_ion_O_4S_cm2( 974) =  3.99 ; given_sigma_ion_O_2D_cm2( 974) =    0.00 ; given_sigma_ion_O_2P_cm2( 974) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 974) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 974) =    0.00
  given_sigma_ion_O_4S_cm2( 975) =  3.99 ; given_sigma_ion_O_2D_cm2( 975) =    0.00 ; given_sigma_ion_O_2P_cm2( 975) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 975) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 975) =    0.00
  given_sigma_ion_O_4S_cm2( 976) =  3.98 ; given_sigma_ion_O_2D_cm2( 976) =    0.00 ; given_sigma_ion_O_2P_cm2( 976) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 976) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 976) =    0.00
  given_sigma_ion_O_4S_cm2( 977) =  3.97 ; given_sigma_ion_O_2D_cm2( 977) =    0.00 ; given_sigma_ion_O_2P_cm2( 977) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 977) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 977) =    0.00
  given_sigma_ion_O_4S_cm2( 978) =  3.97 ; given_sigma_ion_O_2D_cm2( 978) =    0.00 ; given_sigma_ion_O_2P_cm2( 978) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 978) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 978) =    0.00
  given_sigma_ion_O_4S_cm2( 979) =  3.96 ; given_sigma_ion_O_2D_cm2( 979) =    0.00 ; given_sigma_ion_O_2P_cm2( 979) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 979) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 979) =    0.00
  given_sigma_ion_O_4S_cm2( 980) =  3.95 ; given_sigma_ion_O_2D_cm2( 980) =    0.00 ; given_sigma_ion_O_2P_cm2( 980) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 980) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 980) =    0.00
  given_sigma_ion_O_4S_cm2( 981) =  3.95 ; given_sigma_ion_O_2D_cm2( 981) =    0.00 ; given_sigma_ion_O_2P_cm2( 981) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 981) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 981) =    0.00
  given_sigma_ion_O_4S_cm2( 982) =  3.95 ; given_sigma_ion_O_2D_cm2( 982) =    0.00 ; given_sigma_ion_O_2P_cm2( 982) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 982) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 982) =    0.00
  given_sigma_ion_O_4S_cm2( 983) =  3.94 ; given_sigma_ion_O_2D_cm2( 983) =    0.00 ; given_sigma_ion_O_2P_cm2( 983) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 983) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 983) =    0.00
  given_sigma_ion_O_4S_cm2( 984) =  3.93 ; given_sigma_ion_O_2D_cm2( 984) =    0.00 ; given_sigma_ion_O_2P_cm2( 984) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 984) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 984) =    0.00
  given_sigma_ion_O_4S_cm2( 985) =  3.93 ; given_sigma_ion_O_2D_cm2( 985) =    0.00 ; given_sigma_ion_O_2P_cm2( 985) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 985) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 985) =    0.00
  given_sigma_ion_O_4S_cm2( 986) =  3.93 ; given_sigma_ion_O_2D_cm2( 986) =    0.00 ; given_sigma_ion_O_2P_cm2( 986) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 986) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 986) =    0.00
  given_sigma_ion_O_4S_cm2( 987) =  3.92 ; given_sigma_ion_O_2D_cm2( 987) =    0.00 ; given_sigma_ion_O_2P_cm2( 987) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 987) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 987) =    0.00
  given_sigma_ion_O_4S_cm2( 988) =  3.92 ; given_sigma_ion_O_2D_cm2( 988) =    0.00 ; given_sigma_ion_O_2P_cm2( 988) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 988) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 988) =    0.00
  given_sigma_ion_O_4S_cm2( 989) =  3.91 ; given_sigma_ion_O_2D_cm2( 989) =    0.00 ; given_sigma_ion_O_2P_cm2( 989) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 989) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 989) =    0.00
  given_sigma_ion_O_4S_cm2( 990) =  3.90 ; given_sigma_ion_O_2D_cm2( 990) =    0.00 ; given_sigma_ion_O_2P_cm2( 990) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 990) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 990) =    0.00
  given_sigma_ion_O_4S_cm2( 991) =  3.90 ; given_sigma_ion_O_2D_cm2( 991) =    0.00 ; given_sigma_ion_O_2P_cm2( 991) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 991) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 991) =    0.00
  given_sigma_ion_O_4S_cm2( 992) =  3.89 ; given_sigma_ion_O_2D_cm2( 992) =    0.00 ; given_sigma_ion_O_2P_cm2( 992) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 992) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 992) =    0.00
  given_sigma_ion_O_4S_cm2( 993) =  3.88 ; given_sigma_ion_O_2D_cm2( 993) =    0.00 ; given_sigma_ion_O_2P_cm2( 993) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 993) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 993) =    0.00
  given_sigma_ion_O_4S_cm2( 994) =  3.88 ; given_sigma_ion_O_2D_cm2( 994) =    0.00 ; given_sigma_ion_O_2P_cm2( 994) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 994) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 994) =    0.00
  given_sigma_ion_O_4S_cm2( 995) =  3.87 ; given_sigma_ion_O_2D_cm2( 995) =    0.00 ; given_sigma_ion_O_2P_cm2( 995) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 995) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 995) =    0.00
  given_sigma_ion_O_4S_cm2( 996) =  3.86 ; given_sigma_ion_O_2D_cm2( 996) =    0.00 ; given_sigma_ion_O_2P_cm2( 996) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 996) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 996) =    0.00
  given_sigma_ion_O_4S_cm2( 997) =  3.85 ; given_sigma_ion_O_2D_cm2( 997) =    0.00 ; given_sigma_ion_O_2P_cm2( 997) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 997) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 997) =    0.00
  given_sigma_ion_O_4S_cm2( 998) =  3.84 ; given_sigma_ion_O_2D_cm2( 998) =    0.00 ; given_sigma_ion_O_2P_cm2( 998) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 998) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 998) =    0.00
  given_sigma_ion_O_4S_cm2( 999) =  3.83 ; given_sigma_ion_O_2D_cm2( 999) =    0.00 ; given_sigma_ion_O_2P_cm2( 999) =    0.00 ; given_sigma_ion_O_4Pst_cm2( 999) =    0.00 ; given_sigma_ion_O_2Pst_cm2( 999) =    0.00
  given_sigma_ion_O_4S_cm2(1000) =  3.82 ; given_sigma_ion_O_2D_cm2(1000) =    0.00 ; given_sigma_ion_O_2P_cm2(1000) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1000) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1000) =    0.00
  given_sigma_ion_O_4S_cm2(1001) =  3.81 ; given_sigma_ion_O_2D_cm2(1001) =    0.00 ; given_sigma_ion_O_2P_cm2(1001) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1001) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1001) =    0.00
  given_sigma_ion_O_4S_cm2(1002) =  3.80 ; given_sigma_ion_O_2D_cm2(1002) =    0.00 ; given_sigma_ion_O_2P_cm2(1002) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1002) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1002) =    0.00
  given_sigma_ion_O_4S_cm2(1003) =  3.79 ; given_sigma_ion_O_2D_cm2(1003) =    0.00 ; given_sigma_ion_O_2P_cm2(1003) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1003) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1003) =    0.00
  given_sigma_ion_O_4S_cm2(1004) =  3.79 ; given_sigma_ion_O_2D_cm2(1004) =    0.00 ; given_sigma_ion_O_2P_cm2(1004) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1004) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1004) =    0.00
  given_sigma_ion_O_4S_cm2(1005) =  3.78 ; given_sigma_ion_O_2D_cm2(1005) =    0.00 ; given_sigma_ion_O_2P_cm2(1005) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1005) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1005) =    0.00
  given_sigma_ion_O_4S_cm2(1006) =  3.77 ; given_sigma_ion_O_2D_cm2(1006) =    0.00 ; given_sigma_ion_O_2P_cm2(1006) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1006) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1006) =    0.00
  given_sigma_ion_O_4S_cm2(1007) =  3.76 ; given_sigma_ion_O_2D_cm2(1007) =    0.00 ; given_sigma_ion_O_2P_cm2(1007) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1007) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1007) =    0.00
  given_sigma_ion_O_4S_cm2(1008) =  3.75 ; given_sigma_ion_O_2D_cm2(1008) =    0.00 ; given_sigma_ion_O_2P_cm2(1008) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1008) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1008) =    0.00
  given_sigma_ion_O_4S_cm2(1009) =  3.74 ; given_sigma_ion_O_2D_cm2(1009) =    0.00 ; given_sigma_ion_O_2P_cm2(1009) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1009) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1009) =    0.00
  given_sigma_ion_O_4S_cm2(1010) =  3.73 ; given_sigma_ion_O_2D_cm2(1010) =    0.00 ; given_sigma_ion_O_2P_cm2(1010) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1010) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1010) =    0.00
  given_sigma_ion_O_4S_cm2(1011) =  3.73 ; given_sigma_ion_O_2D_cm2(1011) =    0.00 ; given_sigma_ion_O_2P_cm2(1011) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1011) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1011) =    0.00
  given_sigma_ion_O_4S_cm2(1012) =  3.72 ; given_sigma_ion_O_2D_cm2(1012) =    0.00 ; given_sigma_ion_O_2P_cm2(1012) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1012) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1012) =    0.00
  given_sigma_ion_O_4S_cm2(1013) =  3.71 ; given_sigma_ion_O_2D_cm2(1013) =    0.00 ; given_sigma_ion_O_2P_cm2(1013) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1013) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1013) =    0.00
  given_sigma_ion_O_4S_cm2(1014) =  3.70 ; given_sigma_ion_O_2D_cm2(1014) =    0.00 ; given_sigma_ion_O_2P_cm2(1014) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1014) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1014) =    0.00
  given_sigma_ion_O_4S_cm2(1015) =  3.70 ; given_sigma_ion_O_2D_cm2(1015) =    0.00 ; given_sigma_ion_O_2P_cm2(1015) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1015) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1015) =    0.00
  given_sigma_ion_O_4S_cm2(1016) =  3.70 ; given_sigma_ion_O_2D_cm2(1016) =    0.00 ; given_sigma_ion_O_2P_cm2(1016) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1016) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1016) =    0.00
  given_sigma_ion_O_4S_cm2(1017) =  3.70 ; given_sigma_ion_O_2D_cm2(1017) =    0.00 ; given_sigma_ion_O_2P_cm2(1017) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1017) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1017) =    0.00
  given_sigma_ion_O_4S_cm2(1018) =  3.70 ; given_sigma_ion_O_2D_cm2(1018) =    0.00 ; given_sigma_ion_O_2P_cm2(1018) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1018) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1018) =    0.00
  given_sigma_ion_O_4S_cm2(1019) =  3.70 ; given_sigma_ion_O_2D_cm2(1019) =    0.00 ; given_sigma_ion_O_2P_cm2(1019) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1019) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1019) =    0.00
  given_sigma_ion_O_4S_cm2(1020) =  3.70 ; given_sigma_ion_O_2D_cm2(1020) =    0.00 ; given_sigma_ion_O_2P_cm2(1020) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1020) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1020) =    0.00
  given_sigma_ion_O_4S_cm2(1021) =  3.70 ; given_sigma_ion_O_2D_cm2(1021) =    0.00 ; given_sigma_ion_O_2P_cm2(1021) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1021) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1021) =    0.00
  given_sigma_ion_O_4S_cm2(1022) =  3.70 ; given_sigma_ion_O_2D_cm2(1022) =    0.00 ; given_sigma_ion_O_2P_cm2(1022) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1022) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1022) =    0.00
  given_sigma_ion_O_4S_cm2(1023) =  3.70 ; given_sigma_ion_O_2D_cm2(1023) =    0.00 ; given_sigma_ion_O_2P_cm2(1023) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1023) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1023) =    0.00
  given_sigma_ion_O_4S_cm2(1024) =  3.70 ; given_sigma_ion_O_2D_cm2(1024) =    0.00 ; given_sigma_ion_O_2P_cm2(1024) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1024) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1024) =    0.00
  given_sigma_ion_O_4S_cm2(1025) =  3.70 ; given_sigma_ion_O_2D_cm2(1025) =    0.00 ; given_sigma_ion_O_2P_cm2(1025) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1025) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1025) =    0.00
  given_sigma_ion_O_4S_cm2(1026) =  3.70 ; given_sigma_ion_O_2D_cm2(1026) =    0.00 ; given_sigma_ion_O_2P_cm2(1026) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1026) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1026) =    0.00
  given_sigma_ion_O_4S_cm2(1027) =  3.70 ; given_sigma_ion_O_2D_cm2(1027) =    0.00 ; given_sigma_ion_O_2P_cm2(1027) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1027) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1027) =    0.00
  given_sigma_ion_O_4S_cm2(1028) =  3.70 ; given_sigma_ion_O_2D_cm2(1028) =    0.00 ; given_sigma_ion_O_2P_cm2(1028) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1028) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1028) =    0.00
  given_sigma_ion_O_4S_cm2(1029) =  3.70 ; given_sigma_ion_O_2D_cm2(1029) =    0.00 ; given_sigma_ion_O_2P_cm2(1029) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1029) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1029) =    0.00
  given_sigma_ion_O_4S_cm2(1030) =  3.70 ; given_sigma_ion_O_2D_cm2(1030) =    0.00 ; given_sigma_ion_O_2P_cm2(1030) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1030) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1030) =    0.00
  given_sigma_ion_O_4S_cm2(1031) =  3.71 ; given_sigma_ion_O_2D_cm2(1031) =    0.00 ; given_sigma_ion_O_2P_cm2(1031) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1031) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1031) =    0.00
  given_sigma_ion_O_4S_cm2(1032) =  3.71 ; given_sigma_ion_O_2D_cm2(1032) =    0.00 ; given_sigma_ion_O_2P_cm2(1032) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1032) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1032) =    0.00
  given_sigma_ion_O_4S_cm2(1033) =  3.72 ; given_sigma_ion_O_2D_cm2(1033) =    0.00 ; given_sigma_ion_O_2P_cm2(1033) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1033) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1033) =    0.00
  given_sigma_ion_O_4S_cm2(1034) =  3.72 ; given_sigma_ion_O_2D_cm2(1034) =    0.00 ; given_sigma_ion_O_2P_cm2(1034) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1034) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1034) =    0.00
  given_sigma_ion_O_4S_cm2(1035) =  3.72 ; given_sigma_ion_O_2D_cm2(1035) =    0.00 ; given_sigma_ion_O_2P_cm2(1035) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1035) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1035) =    0.00
  given_sigma_ion_O_4S_cm2(1036) =  3.73 ; given_sigma_ion_O_2D_cm2(1036) =    0.00 ; given_sigma_ion_O_2P_cm2(1036) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1036) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1036) =    0.00
  given_sigma_ion_O_4S_cm2(1037) =  3.73 ; given_sigma_ion_O_2D_cm2(1037) =    0.00 ; given_sigma_ion_O_2P_cm2(1037) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1037) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1037) =    0.00
  given_sigma_ion_O_4S_cm2(1038) =  3.74 ; given_sigma_ion_O_2D_cm2(1038) =    0.00 ; given_sigma_ion_O_2P_cm2(1038) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1038) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1038) =    0.00
  given_sigma_ion_O_4S_cm2(1039) =  3.74 ; given_sigma_ion_O_2D_cm2(1039) =    0.00 ; given_sigma_ion_O_2P_cm2(1039) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1039) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1039) =    0.00
  given_sigma_ion_O_4S_cm2(1040) =  3.75 ; given_sigma_ion_O_2D_cm2(1040) =    0.00 ; given_sigma_ion_O_2P_cm2(1040) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1040) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1040) =    0.00
  given_sigma_ion_O_4S_cm2(1041) =  3.76 ; given_sigma_ion_O_2D_cm2(1041) =    0.00 ; given_sigma_ion_O_2P_cm2(1041) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1041) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1041) =    0.00
  given_sigma_ion_O_4S_cm2(1042) =  3.76 ; given_sigma_ion_O_2D_cm2(1042) =    0.00 ; given_sigma_ion_O_2P_cm2(1042) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1042) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1042) =    0.00
  given_sigma_ion_O_4S_cm2(1043) =  3.77 ; given_sigma_ion_O_2D_cm2(1043) =    0.00 ; given_sigma_ion_O_2P_cm2(1043) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1043) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1043) =    0.00
  given_sigma_ion_O_4S_cm2(1044) =  3.77 ; given_sigma_ion_O_2D_cm2(1044) =    0.00 ; given_sigma_ion_O_2P_cm2(1044) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1044) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1044) =    0.00
  given_sigma_ion_O_4S_cm2(1045) =  3.77 ; given_sigma_ion_O_2D_cm2(1045) =    0.00 ; given_sigma_ion_O_2P_cm2(1045) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1045) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1045) =    0.00
  given_sigma_ion_O_4S_cm2(1046) =  3.78 ; given_sigma_ion_O_2D_cm2(1046) =    0.00 ; given_sigma_ion_O_2P_cm2(1046) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1046) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1046) =    0.00
  given_sigma_ion_O_4S_cm2(1047) =  3.78 ; given_sigma_ion_O_2D_cm2(1047) =    0.00 ; given_sigma_ion_O_2P_cm2(1047) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1047) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1047) =    0.00
  given_sigma_ion_O_4S_cm2(1048) =  4.23 ; given_sigma_ion_O_2D_cm2(1048) =    0.00 ; given_sigma_ion_O_2P_cm2(1048) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1048) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1048) =    0.00
  given_sigma_ion_O_4S_cm2(1049) =  3.78 ; given_sigma_ion_O_2D_cm2(1049) =    0.00 ; given_sigma_ion_O_2P_cm2(1049) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1049) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1049) =    0.00
  given_sigma_ion_O_4S_cm2(1050) =  4.86 ; given_sigma_ion_O_2D_cm2(1050) =    0.00 ; given_sigma_ion_O_2P_cm2(1050) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1050) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1050) =    0.00
  given_sigma_ion_O_4S_cm2(1051) = 19.97 ; given_sigma_ion_O_2D_cm2(1051) =    0.00 ; given_sigma_ion_O_2P_cm2(1051) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1051) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1051) =    0.00
  given_sigma_ion_O_4S_cm2(1052) =  2.70 ; given_sigma_ion_O_2D_cm2(1052) =    0.00 ; given_sigma_ion_O_2P_cm2(1052) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1052) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1052) =    0.00
  given_sigma_ion_O_4S_cm2(1053) =  2.74 ; given_sigma_ion_O_2D_cm2(1053) =    0.00 ; given_sigma_ion_O_2P_cm2(1053) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1053) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1053) =    0.00
  given_sigma_ion_O_4S_cm2(1054) =  2.80 ; given_sigma_ion_O_2D_cm2(1054) =    0.00 ; given_sigma_ion_O_2P_cm2(1054) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1054) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1054) =    0.00
  given_sigma_ion_O_4S_cm2(1055) =  2.86 ; given_sigma_ion_O_2D_cm2(1055) =    0.00 ; given_sigma_ion_O_2P_cm2(1055) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1055) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1055) =    0.00
  given_sigma_ion_O_4S_cm2(1056) =  2.92 ; given_sigma_ion_O_2D_cm2(1056) =    0.00 ; given_sigma_ion_O_2P_cm2(1056) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1056) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1056) =    0.00
  given_sigma_ion_O_4S_cm2(1057) =  2.94 ; given_sigma_ion_O_2D_cm2(1057) =    0.00 ; given_sigma_ion_O_2P_cm2(1057) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1057) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1057) =    0.00
  given_sigma_ion_O_4S_cm2(1058) =  3.03 ; given_sigma_ion_O_2D_cm2(1058) =    0.00 ; given_sigma_ion_O_2P_cm2(1058) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1058) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1058) =    0.00
  given_sigma_ion_O_4S_cm2(1059) =  3.15 ; given_sigma_ion_O_2D_cm2(1059) =    0.00 ; given_sigma_ion_O_2P_cm2(1059) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1059) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1059) =    0.00
  given_sigma_ion_O_4S_cm2(1060) =  3.23 ; given_sigma_ion_O_2D_cm2(1060) =    0.00 ; given_sigma_ion_O_2P_cm2(1060) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1060) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1060) =    0.00
  given_sigma_ion_O_4S_cm2(1061) =  3.26 ; given_sigma_ion_O_2D_cm2(1061) =    0.00 ; given_sigma_ion_O_2P_cm2(1061) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1061) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1061) =    0.00
  given_sigma_ion_O_4S_cm2(1062) =  3.38 ; given_sigma_ion_O_2D_cm2(1062) =    0.00 ; given_sigma_ion_O_2P_cm2(1062) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1062) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1062) =    0.00
  given_sigma_ion_O_4S_cm2(1063) =  3.44 ; given_sigma_ion_O_2D_cm2(1063) =    0.00 ; given_sigma_ion_O_2P_cm2(1063) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1063) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1063) =    0.00
  given_sigma_ion_O_4S_cm2(1064) =  3.50 ; given_sigma_ion_O_2D_cm2(1064) =    0.00 ; given_sigma_ion_O_2P_cm2(1064) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1064) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1064) =    0.00
  given_sigma_ion_O_4S_cm2(1065) =  3.55 ; given_sigma_ion_O_2D_cm2(1065) =    0.00 ; given_sigma_ion_O_2P_cm2(1065) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1065) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1065) =    0.00
  given_sigma_ion_O_4S_cm2(1066) =  3.61 ; given_sigma_ion_O_2D_cm2(1066) =    0.00 ; given_sigma_ion_O_2P_cm2(1066) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1066) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1066) =    0.00
  given_sigma_ion_O_4S_cm2(1067) =  3.67 ; given_sigma_ion_O_2D_cm2(1067) =    0.00 ; given_sigma_ion_O_2P_cm2(1067) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1067) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1067) =    0.00
  given_sigma_ion_O_4S_cm2(1068) =  3.73 ; given_sigma_ion_O_2D_cm2(1068) =    0.00 ; given_sigma_ion_O_2P_cm2(1068) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1068) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1068) =    0.00
  given_sigma_ion_O_4S_cm2(1069) =  3.80 ; given_sigma_ion_O_2D_cm2(1069) =    0.00 ; given_sigma_ion_O_2P_cm2(1069) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1069) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1069) =    0.00
  given_sigma_ion_O_4S_cm2(1070) =  3.84 ; given_sigma_ion_O_2D_cm2(1070) =    0.00 ; given_sigma_ion_O_2P_cm2(1070) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1070) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1070) =    0.00
  given_sigma_ion_O_4S_cm2(1071) =  3.90 ; given_sigma_ion_O_2D_cm2(1071) =    0.00 ; given_sigma_ion_O_2P_cm2(1071) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1071) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1071) =    0.00
  given_sigma_ion_O_4S_cm2(1072) =  3.96 ; given_sigma_ion_O_2D_cm2(1072) =    0.00 ; given_sigma_ion_O_2P_cm2(1072) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1072) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1072) =    0.00
  given_sigma_ion_O_4S_cm2(1073) = 36.41 ; given_sigma_ion_O_2D_cm2(1073) =    0.00 ; given_sigma_ion_O_2P_cm2(1073) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1073) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1073) =    0.00
  given_sigma_ion_O_4S_cm2(1074) =  9.14 ; given_sigma_ion_O_2D_cm2(1074) =    0.00 ; given_sigma_ion_O_2P_cm2(1074) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1074) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1074) =    0.00
  given_sigma_ion_O_4S_cm2(1075) =  5.34 ; given_sigma_ion_O_2D_cm2(1075) =    0.00 ; given_sigma_ion_O_2P_cm2(1075) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1075) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1075) =    0.00
  given_sigma_ion_O_4S_cm2(1076) = 48.94 ; given_sigma_ion_O_2D_cm2(1076) =    0.00 ; given_sigma_ion_O_2P_cm2(1076) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1076) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1076) =    0.00
  given_sigma_ion_O_4S_cm2(1077) =  7.20 ; given_sigma_ion_O_2D_cm2(1077) =    0.00 ; given_sigma_ion_O_2P_cm2(1077) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1077) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1077) =    0.00
  given_sigma_ion_O_4S_cm2(1078) = 11.73 ; given_sigma_ion_O_2D_cm2(1078) =    0.00 ; given_sigma_ion_O_2P_cm2(1078) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1078) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1078) =    0.00
  given_sigma_ion_O_4S_cm2(1079) =  3.64 ; given_sigma_ion_O_2D_cm2(1079) =    0.00 ; given_sigma_ion_O_2P_cm2(1079) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1079) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1079) =    0.00
  given_sigma_ion_O_4S_cm2(1080) = 15.37 ; given_sigma_ion_O_2D_cm2(1080) =    0.00 ; given_sigma_ion_O_2P_cm2(1080) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1080) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1080) =    0.00
  given_sigma_ion_O_4S_cm2(1081) =  3.64 ; given_sigma_ion_O_2D_cm2(1081) =    0.00 ; given_sigma_ion_O_2P_cm2(1081) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1081) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1081) =    0.00
  given_sigma_ion_O_4S_cm2(1082) =  5.91 ; given_sigma_ion_O_2D_cm2(1082) =    0.00 ; given_sigma_ion_O_2P_cm2(1082) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1082) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1082) =    0.00
  given_sigma_ion_O_4S_cm2(1083) =  1.62 ; given_sigma_ion_O_2D_cm2(1083) =    0.00 ; given_sigma_ion_O_2P_cm2(1083) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1083) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1083) =    0.00
  given_sigma_ion_O_4S_cm2(1084) =  2.06 ; given_sigma_ion_O_2D_cm2(1084) =    0.00 ; given_sigma_ion_O_2P_cm2(1084) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1084) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1084) =    0.00
  given_sigma_ion_O_4S_cm2(1085) =  2.13 ; given_sigma_ion_O_2D_cm2(1085) =    0.00 ; given_sigma_ion_O_2P_cm2(1085) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1085) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1085) =    0.00
  given_sigma_ion_O_4S_cm2(1086) =  2.14 ; given_sigma_ion_O_2D_cm2(1086) =    0.00 ; given_sigma_ion_O_2P_cm2(1086) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1086) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1086) =    0.00
  given_sigma_ion_O_4S_cm2(1087) =  2.17 ; given_sigma_ion_O_2D_cm2(1087) =    0.00 ; given_sigma_ion_O_2P_cm2(1087) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1087) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1087) =    0.00
  given_sigma_ion_O_4S_cm2(1088) =  2.19 ; given_sigma_ion_O_2D_cm2(1088) =    0.00 ; given_sigma_ion_O_2P_cm2(1088) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1088) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1088) =    0.00
  given_sigma_ion_O_4S_cm2(1089) =  2.20 ; given_sigma_ion_O_2D_cm2(1089) =    0.00 ; given_sigma_ion_O_2P_cm2(1089) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1089) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1089) =    0.00
  given_sigma_ion_O_4S_cm2(1090) =  2.22 ; given_sigma_ion_O_2D_cm2(1090) =    0.00 ; given_sigma_ion_O_2P_cm2(1090) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1090) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1090) =    0.00
  given_sigma_ion_O_4S_cm2(1091) =  2.24 ; given_sigma_ion_O_2D_cm2(1091) =    0.00 ; given_sigma_ion_O_2P_cm2(1091) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1091) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1091) =    0.00
  given_sigma_ion_O_4S_cm2(1092) =  2.27 ; given_sigma_ion_O_2D_cm2(1092) =    0.00 ; given_sigma_ion_O_2P_cm2(1092) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1092) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1092) =    0.00
  given_sigma_ion_O_4S_cm2(1093) =  2.31 ; given_sigma_ion_O_2D_cm2(1093) =    0.00 ; given_sigma_ion_O_2P_cm2(1093) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1093) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1093) =    0.00
  given_sigma_ion_O_4S_cm2(1094) = 15.23 ; given_sigma_ion_O_2D_cm2(1094) =    0.00 ; given_sigma_ion_O_2P_cm2(1094) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1094) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1094) =    0.00
  given_sigma_ion_O_4S_cm2(1095) =  2.82 ; given_sigma_ion_O_2D_cm2(1095) =    0.00 ; given_sigma_ion_O_2P_cm2(1095) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1095) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1095) =    0.00
  given_sigma_ion_O_4S_cm2(1096) = 32.62 ; given_sigma_ion_O_2D_cm2(1096) =    0.00 ; given_sigma_ion_O_2P_cm2(1096) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1096) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1096) =    0.00
  given_sigma_ion_O_4S_cm2(1097) =  3.39 ; given_sigma_ion_O_2D_cm2(1097) =    0.00 ; given_sigma_ion_O_2P_cm2(1097) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1097) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1097) =    0.00
  given_sigma_ion_O_4S_cm2(1098) =  7.90 ; given_sigma_ion_O_2D_cm2(1098) =    0.00 ; given_sigma_ion_O_2P_cm2(1098) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1098) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1098) =    0.00
  given_sigma_ion_O_4S_cm2(1099) =  4.74 ; given_sigma_ion_O_2D_cm2(1099) =    0.00 ; given_sigma_ion_O_2P_cm2(1099) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1099) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1099) =    0.00
  given_sigma_ion_O_4S_cm2(1100) =  2.80 ; given_sigma_ion_O_2D_cm2(1100) =    0.00 ; given_sigma_ion_O_2P_cm2(1100) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1100) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1100) =    0.00
  given_sigma_ion_O_4S_cm2(1101) =  3.67 ; given_sigma_ion_O_2D_cm2(1101) =    0.00 ; given_sigma_ion_O_2P_cm2(1101) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1101) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1101) =    0.00
  given_sigma_ion_O_4S_cm2(1102) =  3.95 ; given_sigma_ion_O_2D_cm2(1102) =    0.00 ; given_sigma_ion_O_2P_cm2(1102) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1102) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1102) =    0.00
  given_sigma_ion_O_4S_cm2(1103) =  1.97 ; given_sigma_ion_O_2D_cm2(1103) =    0.00 ; given_sigma_ion_O_2P_cm2(1103) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1103) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1103) =    0.00
  given_sigma_ion_O_4S_cm2(1104) =  2.21 ; given_sigma_ion_O_2D_cm2(1104) =    0.00 ; given_sigma_ion_O_2P_cm2(1104) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1104) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1104) =    0.00
  given_sigma_ion_O_4S_cm2(1105) =  2.38 ; given_sigma_ion_O_2D_cm2(1105) =    0.00 ; given_sigma_ion_O_2P_cm2(1105) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1105) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1105) =    0.00
  given_sigma_ion_O_4S_cm2(1106) =  2.79 ; given_sigma_ion_O_2D_cm2(1106) =    0.00 ; given_sigma_ion_O_2P_cm2(1106) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1106) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1106) =    0.00
  given_sigma_ion_O_4S_cm2(1107) =  2.95 ; given_sigma_ion_O_2D_cm2(1107) =    0.00 ; given_sigma_ion_O_2P_cm2(1107) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1107) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1107) =    0.00
  given_sigma_ion_O_4S_cm2(1108) =  3.19 ; given_sigma_ion_O_2D_cm2(1108) =    0.00 ; given_sigma_ion_O_2P_cm2(1108) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1108) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1108) =    0.00
  given_sigma_ion_O_4S_cm2(1109) =  3.60 ; given_sigma_ion_O_2D_cm2(1109) =    0.00 ; given_sigma_ion_O_2P_cm2(1109) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1109) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1109) =    0.00
  given_sigma_ion_O_4S_cm2(1110) =  3.93 ; given_sigma_ion_O_2D_cm2(1110) =    0.00 ; given_sigma_ion_O_2P_cm2(1110) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1110) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1110) =    0.00
  given_sigma_ion_O_4S_cm2(1111) =  4.42 ; given_sigma_ion_O_2D_cm2(1111) =    0.00 ; given_sigma_ion_O_2P_cm2(1111) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1111) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1111) =    0.00
  given_sigma_ion_O_4S_cm2(1112) =  4.83 ; given_sigma_ion_O_2D_cm2(1112) =    0.00 ; given_sigma_ion_O_2P_cm2(1112) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1112) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1112) =    0.00
  given_sigma_ion_O_4S_cm2(1113) =  5.17 ; given_sigma_ion_O_2D_cm2(1113) =    0.00 ; given_sigma_ion_O_2P_cm2(1113) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1113) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1113) =    0.00
  given_sigma_ion_O_4S_cm2(1114) = 65.14 ; given_sigma_ion_O_2D_cm2(1114) =    0.00 ; given_sigma_ion_O_2P_cm2(1114) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1114) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1114) =    0.00
  given_sigma_ion_O_4S_cm2(1115) =  4.14 ; given_sigma_ion_O_2D_cm2(1115) =    0.00 ; given_sigma_ion_O_2P_cm2(1115) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1115) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1115) =    0.00
  given_sigma_ion_O_4S_cm2(1116) =  1.86 ; given_sigma_ion_O_2D_cm2(1116) =    0.00 ; given_sigma_ion_O_2P_cm2(1116) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1116) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1116) =    0.00
  given_sigma_ion_O_4S_cm2(1117) =  3.10 ; given_sigma_ion_O_2D_cm2(1117) =    0.00 ; given_sigma_ion_O_2P_cm2(1117) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1117) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1117) =    0.00
  given_sigma_ion_O_4S_cm2(1118) = 22.02 ; given_sigma_ion_O_2D_cm2(1118) =    0.00 ; given_sigma_ion_O_2P_cm2(1118) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1118) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1118) =    0.00
  given_sigma_ion_O_4S_cm2(1119) =  2.59 ; given_sigma_ion_O_2D_cm2(1119) =    0.00 ; given_sigma_ion_O_2P_cm2(1119) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1119) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1119) =    0.00
  given_sigma_ion_O_4S_cm2(1120) =  7.40 ; given_sigma_ion_O_2D_cm2(1120) =    0.00 ; given_sigma_ion_O_2P_cm2(1120) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1120) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1120) =    0.00
  given_sigma_ion_O_4S_cm2(1121) =  1.96 ; given_sigma_ion_O_2D_cm2(1121) =    0.00 ; given_sigma_ion_O_2P_cm2(1121) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1121) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1121) =    0.00
  given_sigma_ion_O_4S_cm2(1122) =  1.55 ; given_sigma_ion_O_2D_cm2(1122) =    0.00 ; given_sigma_ion_O_2P_cm2(1122) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1122) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1122) =    0.00
  given_sigma_ion_O_4S_cm2(1123) =  1.68 ; given_sigma_ion_O_2D_cm2(1123) =    0.00 ; given_sigma_ion_O_2P_cm2(1123) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1123) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1123) =    0.00
  given_sigma_ion_O_4S_cm2(1124) =  1.82 ; given_sigma_ion_O_2D_cm2(1124) =    0.00 ; given_sigma_ion_O_2P_cm2(1124) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1124) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1124) =    0.00
  given_sigma_ion_O_4S_cm2(1125) =  1.90 ; given_sigma_ion_O_2D_cm2(1125) =    0.00 ; given_sigma_ion_O_2P_cm2(1125) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1125) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1125) =    0.00
  given_sigma_ion_O_4S_cm2(1126) =  1.95 ; given_sigma_ion_O_2D_cm2(1126) =    0.00 ; given_sigma_ion_O_2P_cm2(1126) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1126) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1126) =    0.00
  given_sigma_ion_O_4S_cm2(1127) =  2.09 ; given_sigma_ion_O_2D_cm2(1127) =    0.00 ; given_sigma_ion_O_2P_cm2(1127) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1127) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1127) =    0.00
  given_sigma_ion_O_4S_cm2(1128) =  2.22 ; given_sigma_ion_O_2D_cm2(1128) =    0.00 ; given_sigma_ion_O_2P_cm2(1128) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1128) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1128) =    0.00
  given_sigma_ion_O_4S_cm2(1129) =  2.30 ; given_sigma_ion_O_2D_cm2(1129) =    0.00 ; given_sigma_ion_O_2P_cm2(1129) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1129) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1129) =    0.00
  given_sigma_ion_O_4S_cm2(1130) =  2.36 ; given_sigma_ion_O_2D_cm2(1130) =    0.00 ; given_sigma_ion_O_2P_cm2(1130) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1130) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1130) =    0.00
  given_sigma_ion_O_4S_cm2(1131) =  2.49 ; given_sigma_ion_O_2D_cm2(1131) =    0.00 ; given_sigma_ion_O_2P_cm2(1131) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1131) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1131) =    0.00
  given_sigma_ion_O_4S_cm2(1132) =  2.76 ; given_sigma_ion_O_2D_cm2(1132) =    0.00 ; given_sigma_ion_O_2P_cm2(1132) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1132) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1132) =    0.00
  given_sigma_ion_O_4S_cm2(1133) =  2.82 ; given_sigma_ion_O_2D_cm2(1133) =    0.00 ; given_sigma_ion_O_2P_cm2(1133) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1133) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1133) =    0.00
  given_sigma_ion_O_4S_cm2(1134) =  3.03 ; given_sigma_ion_O_2D_cm2(1134) =    0.00 ; given_sigma_ion_O_2P_cm2(1134) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1134) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1134) =    0.00
  given_sigma_ion_O_4S_cm2(1135) =  3.17 ; given_sigma_ion_O_2D_cm2(1135) =    0.00 ; given_sigma_ion_O_2P_cm2(1135) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1135) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1135) =    0.00
  given_sigma_ion_O_4S_cm2(1136) =  3.27 ; given_sigma_ion_O_2D_cm2(1136) =    0.00 ; given_sigma_ion_O_2P_cm2(1136) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1136) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1136) =    0.00
  given_sigma_ion_O_4S_cm2(1137) =  3.30 ; given_sigma_ion_O_2D_cm2(1137) =    0.00 ; given_sigma_ion_O_2P_cm2(1137) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1137) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1137) =    0.00
  given_sigma_ion_O_4S_cm2(1138) =  3.30 ; given_sigma_ion_O_2D_cm2(1138) =    0.00 ; given_sigma_ion_O_2P_cm2(1138) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1138) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1138) =    0.00
  given_sigma_ion_O_4S_cm2(1139) =  3.30 ; given_sigma_ion_O_2D_cm2(1139) =    0.00 ; given_sigma_ion_O_2P_cm2(1139) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1139) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1139) =    0.00
  given_sigma_ion_O_4S_cm2(1140) =  3.30 ; given_sigma_ion_O_2D_cm2(1140) =    0.00 ; given_sigma_ion_O_2P_cm2(1140) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1140) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1140) =    0.00
  given_sigma_ion_O_4S_cm2(1141) =  3.30 ; given_sigma_ion_O_2D_cm2(1141) =    0.00 ; given_sigma_ion_O_2P_cm2(1141) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1141) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1141) =    0.00
  given_sigma_ion_O_4S_cm2(1142) =  3.30 ; given_sigma_ion_O_2D_cm2(1142) =    0.00 ; given_sigma_ion_O_2P_cm2(1142) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1142) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1142) =    0.00
  given_sigma_ion_O_4S_cm2(1143) =  3.30 ; given_sigma_ion_O_2D_cm2(1143) =    0.00 ; given_sigma_ion_O_2P_cm2(1143) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1143) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1143) =    0.00
  given_sigma_ion_O_4S_cm2(1144) =  3.30 ; given_sigma_ion_O_2D_cm2(1144) =    0.00 ; given_sigma_ion_O_2P_cm2(1144) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1144) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1144) =    0.00
  given_sigma_ion_O_4S_cm2(1145) =  3.30 ; given_sigma_ion_O_2D_cm2(1145) =    0.00 ; given_sigma_ion_O_2P_cm2(1145) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1145) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1145) =    0.00
  given_sigma_ion_O_4S_cm2(1146) =  3.30 ; given_sigma_ion_O_2D_cm2(1146) =    0.00 ; given_sigma_ion_O_2P_cm2(1146) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1146) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1146) =    0.00
  given_sigma_ion_O_4S_cm2(1147) =  3.30 ; given_sigma_ion_O_2D_cm2(1147) =    0.00 ; given_sigma_ion_O_2P_cm2(1147) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1147) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1147) =    0.00
  given_sigma_ion_O_4S_cm2(1148) =  3.30 ; given_sigma_ion_O_2D_cm2(1148) =    0.00 ; given_sigma_ion_O_2P_cm2(1148) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1148) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1148) =    0.00
  given_sigma_ion_O_4S_cm2(1149) =  3.30 ; given_sigma_ion_O_2D_cm2(1149) =    0.00 ; given_sigma_ion_O_2P_cm2(1149) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1149) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1149) =    0.00
  given_sigma_ion_O_4S_cm2(1150) =  3.30 ; given_sigma_ion_O_2D_cm2(1150) =    0.00 ; given_sigma_ion_O_2P_cm2(1150) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1150) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1150) =    0.00
  given_sigma_ion_O_4S_cm2(1151) =  3.30 ; given_sigma_ion_O_2D_cm2(1151) =    0.00 ; given_sigma_ion_O_2P_cm2(1151) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1151) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1151) =    0.00
  given_sigma_ion_O_4S_cm2(1152) =  3.30 ; given_sigma_ion_O_2D_cm2(1152) =    0.00 ; given_sigma_ion_O_2P_cm2(1152) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1152) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1152) =    0.00
  given_sigma_ion_O_4S_cm2(1153) =  3.30 ; given_sigma_ion_O_2D_cm2(1153) =    0.00 ; given_sigma_ion_O_2P_cm2(1153) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1153) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1153) =    0.00
  given_sigma_ion_O_4S_cm2(1154) =  3.30 ; given_sigma_ion_O_2D_cm2(1154) =    0.00 ; given_sigma_ion_O_2P_cm2(1154) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1154) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1154) =    0.00
  given_sigma_ion_O_4S_cm2(1155) =  3.30 ; given_sigma_ion_O_2D_cm2(1155) =    0.00 ; given_sigma_ion_O_2P_cm2(1155) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1155) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1155) =    0.00
  given_sigma_ion_O_4S_cm2(1156) =  3.30 ; given_sigma_ion_O_2D_cm2(1156) =    0.00 ; given_sigma_ion_O_2P_cm2(1156) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1156) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1156) =    0.00
  given_sigma_ion_O_4S_cm2(1157) =  3.30 ; given_sigma_ion_O_2D_cm2(1157) =    0.00 ; given_sigma_ion_O_2P_cm2(1157) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1157) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1157) =    0.00
  given_sigma_ion_O_4S_cm2(1158) =  3.30 ; given_sigma_ion_O_2D_cm2(1158) =    0.00 ; given_sigma_ion_O_2P_cm2(1158) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1158) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1158) =    0.00
  given_sigma_ion_O_4S_cm2(1159) =  3.30 ; given_sigma_ion_O_2D_cm2(1159) =    0.00 ; given_sigma_ion_O_2P_cm2(1159) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1159) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1159) =    0.00
  given_sigma_ion_O_4S_cm2(1160) =  3.30 ; given_sigma_ion_O_2D_cm2(1160) =    0.00 ; given_sigma_ion_O_2P_cm2(1160) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1160) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1160) =    0.00
  given_sigma_ion_O_4S_cm2(1161) =  3.30 ; given_sigma_ion_O_2D_cm2(1161) =    0.00 ; given_sigma_ion_O_2P_cm2(1161) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1161) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1161) =    0.00
  given_sigma_ion_O_4S_cm2(1162) =  3.30 ; given_sigma_ion_O_2D_cm2(1162) =    0.00 ; given_sigma_ion_O_2P_cm2(1162) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1162) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1162) =    0.00
  given_sigma_ion_O_4S_cm2(1163) =  3.30 ; given_sigma_ion_O_2D_cm2(1163) =    0.00 ; given_sigma_ion_O_2P_cm2(1163) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1163) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1163) =    0.00
  given_sigma_ion_O_4S_cm2(1164) =  3.30 ; given_sigma_ion_O_2D_cm2(1164) =    0.00 ; given_sigma_ion_O_2P_cm2(1164) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1164) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1164) =    0.00
  given_sigma_ion_O_4S_cm2(1165) =  3.30 ; given_sigma_ion_O_2D_cm2(1165) =    0.00 ; given_sigma_ion_O_2P_cm2(1165) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1165) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1165) =    0.00
  given_sigma_ion_O_4S_cm2(1166) =  3.30 ; given_sigma_ion_O_2D_cm2(1166) =    0.00 ; given_sigma_ion_O_2P_cm2(1166) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1166) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1166) =    0.00
  given_sigma_ion_O_4S_cm2(1167) =  3.30 ; given_sigma_ion_O_2D_cm2(1167) =    0.00 ; given_sigma_ion_O_2P_cm2(1167) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1167) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1167) =    0.00
  given_sigma_ion_O_4S_cm2(1168) =  3.30 ; given_sigma_ion_O_2D_cm2(1168) =    0.00 ; given_sigma_ion_O_2P_cm2(1168) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1168) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1168) =    0.00
  given_sigma_ion_O_4S_cm2(1169) =  3.30 ; given_sigma_ion_O_2D_cm2(1169) =    0.00 ; given_sigma_ion_O_2P_cm2(1169) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1169) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1169) =    0.00
  given_sigma_ion_O_4S_cm2(1170) =  3.30 ; given_sigma_ion_O_2D_cm2(1170) =    0.00 ; given_sigma_ion_O_2P_cm2(1170) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1170) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1170) =    0.00
  given_sigma_ion_O_4S_cm2(1171) =  3.30 ; given_sigma_ion_O_2D_cm2(1171) =    0.00 ; given_sigma_ion_O_2P_cm2(1171) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1171) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1171) =    0.00
  given_sigma_ion_O_4S_cm2(1172) =  3.30 ; given_sigma_ion_O_2D_cm2(1172) =    0.00 ; given_sigma_ion_O_2P_cm2(1172) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1172) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1172) =    0.00
  given_sigma_ion_O_4S_cm2(1173) =  3.30 ; given_sigma_ion_O_2D_cm2(1173) =    0.00 ; given_sigma_ion_O_2P_cm2(1173) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1173) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1173) =    0.00
  given_sigma_ion_O_4S_cm2(1174) =  3.30 ; given_sigma_ion_O_2D_cm2(1174) =    0.00 ; given_sigma_ion_O_2P_cm2(1174) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1174) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1174) =    0.00
  given_sigma_ion_O_4S_cm2(1175) =  3.29 ; given_sigma_ion_O_2D_cm2(1175) =    0.00 ; given_sigma_ion_O_2P_cm2(1175) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1175) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1175) =    0.00
  given_sigma_ion_O_4S_cm2(1176) =  3.29 ; given_sigma_ion_O_2D_cm2(1176) =    0.00 ; given_sigma_ion_O_2P_cm2(1176) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1176) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1176) =    0.00
  given_sigma_ion_O_4S_cm2(1177) =  3.29 ; given_sigma_ion_O_2D_cm2(1177) =    0.00 ; given_sigma_ion_O_2P_cm2(1177) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1177) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1177) =    0.00
  given_sigma_ion_O_4S_cm2(1178) =  3.28 ; given_sigma_ion_O_2D_cm2(1178) =    0.00 ; given_sigma_ion_O_2P_cm2(1178) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1178) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1178) =    0.00
  given_sigma_ion_O_4S_cm2(1179) =  3.28 ; given_sigma_ion_O_2D_cm2(1179) =    0.00 ; given_sigma_ion_O_2P_cm2(1179) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1179) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1179) =    0.00
  given_sigma_ion_O_4S_cm2(1180) =  3.27 ; given_sigma_ion_O_2D_cm2(1180) =    0.00 ; given_sigma_ion_O_2P_cm2(1180) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1180) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1180) =    0.00
  given_sigma_ion_O_4S_cm2(1181) =  3.27 ; given_sigma_ion_O_2D_cm2(1181) =    0.00 ; given_sigma_ion_O_2P_cm2(1181) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1181) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1181) =    0.00
  given_sigma_ion_O_4S_cm2(1182) =  3.26 ; given_sigma_ion_O_2D_cm2(1182) =    0.00 ; given_sigma_ion_O_2P_cm2(1182) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1182) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1182) =    0.00
  given_sigma_ion_O_4S_cm2(1183) =  3.26 ; given_sigma_ion_O_2D_cm2(1183) =    0.00 ; given_sigma_ion_O_2P_cm2(1183) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1183) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1183) =    0.00
  given_sigma_ion_O_4S_cm2(1184) =  3.26 ; given_sigma_ion_O_2D_cm2(1184) =    0.00 ; given_sigma_ion_O_2P_cm2(1184) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1184) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1184) =    0.00
  given_sigma_ion_O_4S_cm2(1185) =  3.25 ; given_sigma_ion_O_2D_cm2(1185) =    0.00 ; given_sigma_ion_O_2P_cm2(1185) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1185) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1185) =    0.00
  given_sigma_ion_O_4S_cm2(1186) =  3.25 ; given_sigma_ion_O_2D_cm2(1186) =    0.00 ; given_sigma_ion_O_2P_cm2(1186) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1186) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1186) =    0.00
  given_sigma_ion_O_4S_cm2(1187) =  3.24 ; given_sigma_ion_O_2D_cm2(1187) =    0.00 ; given_sigma_ion_O_2P_cm2(1187) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1187) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1187) =    0.00
  given_sigma_ion_O_4S_cm2(1188) =  3.23 ; given_sigma_ion_O_2D_cm2(1188) =    0.00 ; given_sigma_ion_O_2P_cm2(1188) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1188) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1188) =    0.00
  given_sigma_ion_O_4S_cm2(1189) =  3.22 ; given_sigma_ion_O_2D_cm2(1189) =    0.00 ; given_sigma_ion_O_2P_cm2(1189) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1189) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1189) =    0.00
  given_sigma_ion_O_4S_cm2(1190) =  3.22 ; given_sigma_ion_O_2D_cm2(1190) =    0.00 ; given_sigma_ion_O_2P_cm2(1190) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1190) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1190) =    0.00
  given_sigma_ion_O_4S_cm2(1191) =  3.22 ; given_sigma_ion_O_2D_cm2(1191) =    0.00 ; given_sigma_ion_O_2P_cm2(1191) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1191) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1191) =    0.00
  given_sigma_ion_O_4S_cm2(1192) =  3.21 ; given_sigma_ion_O_2D_cm2(1192) =    0.00 ; given_sigma_ion_O_2P_cm2(1192) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1192) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1192) =    0.00
  given_sigma_ion_O_4S_cm2(1193) =  3.21 ; given_sigma_ion_O_2D_cm2(1193) =    0.00 ; given_sigma_ion_O_2P_cm2(1193) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1193) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1193) =    0.00
  given_sigma_ion_O_4S_cm2(1194) =  3.21 ; given_sigma_ion_O_2D_cm2(1194) =    0.00 ; given_sigma_ion_O_2P_cm2(1194) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1194) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1194) =    0.00
  given_sigma_ion_O_4S_cm2(1195) =  3.20 ; given_sigma_ion_O_2D_cm2(1195) =    0.00 ; given_sigma_ion_O_2P_cm2(1195) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1195) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1195) =    0.00
  given_sigma_ion_O_4S_cm2(1196) =  3.19 ; given_sigma_ion_O_2D_cm2(1196) =    0.00 ; given_sigma_ion_O_2P_cm2(1196) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1196) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1196) =    0.00
  given_sigma_ion_O_4S_cm2(1197) =  3.18 ; given_sigma_ion_O_2D_cm2(1197) =    0.00 ; given_sigma_ion_O_2P_cm2(1197) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1197) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1197) =    0.00
  given_sigma_ion_O_4S_cm2(1198) =  3.17 ; given_sigma_ion_O_2D_cm2(1198) =    0.00 ; given_sigma_ion_O_2P_cm2(1198) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1198) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1198) =    0.00
  given_sigma_ion_O_4S_cm2(1199) =  3.16 ; given_sigma_ion_O_2D_cm2(1199) =    0.00 ; given_sigma_ion_O_2P_cm2(1199) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1199) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1199) =    0.00
  given_sigma_ion_O_4S_cm2(1200) =  3.16 ; given_sigma_ion_O_2D_cm2(1200) =    0.00 ; given_sigma_ion_O_2P_cm2(1200) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1200) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1200) =    0.00
  given_sigma_ion_O_4S_cm2(1201) =  3.16 ; given_sigma_ion_O_2D_cm2(1201) =    0.00 ; given_sigma_ion_O_2P_cm2(1201) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1201) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1201) =    0.00
  given_sigma_ion_O_4S_cm2(1202) =  3.15 ; given_sigma_ion_O_2D_cm2(1202) =    0.00 ; given_sigma_ion_O_2P_cm2(1202) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1202) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1202) =    0.00
  given_sigma_ion_O_4S_cm2(1203) =  3.14 ; given_sigma_ion_O_2D_cm2(1203) =    0.00 ; given_sigma_ion_O_2P_cm2(1203) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1203) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1203) =    0.00
  given_sigma_ion_O_4S_cm2(1204) =  3.14 ; given_sigma_ion_O_2D_cm2(1204) =    0.00 ; given_sigma_ion_O_2P_cm2(1204) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1204) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1204) =    0.00
  given_sigma_ion_O_4S_cm2(1205) =  3.13 ; given_sigma_ion_O_2D_cm2(1205) =    0.00 ; given_sigma_ion_O_2P_cm2(1205) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1205) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1205) =    0.00
  given_sigma_ion_O_4S_cm2(1206) =  3.12 ; given_sigma_ion_O_2D_cm2(1206) =    0.00 ; given_sigma_ion_O_2P_cm2(1206) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1206) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1206) =    0.00
  given_sigma_ion_O_4S_cm2(1207) =  3.11 ; given_sigma_ion_O_2D_cm2(1207) =    0.00 ; given_sigma_ion_O_2P_cm2(1207) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1207) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1207) =    0.00
  given_sigma_ion_O_4S_cm2(1208) =  3.10 ; given_sigma_ion_O_2D_cm2(1208) =    0.00 ; given_sigma_ion_O_2P_cm2(1208) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1208) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1208) =    0.00
  given_sigma_ion_O_4S_cm2(1209) =  3.10 ; given_sigma_ion_O_2D_cm2(1209) =    0.00 ; given_sigma_ion_O_2P_cm2(1209) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1209) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1209) =    0.00
  given_sigma_ion_O_4S_cm2(1210) =  3.10 ; given_sigma_ion_O_2D_cm2(1210) =    0.00 ; given_sigma_ion_O_2P_cm2(1210) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1210) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1210) =    0.00
  given_sigma_ion_O_4S_cm2(1211) =  3.10 ; given_sigma_ion_O_2D_cm2(1211) =    0.00 ; given_sigma_ion_O_2P_cm2(1211) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1211) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1211) =    0.00
  given_sigma_ion_O_4S_cm2(1212) =  3.10 ; given_sigma_ion_O_2D_cm2(1212) =    0.00 ; given_sigma_ion_O_2P_cm2(1212) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1212) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1212) =    0.00
  given_sigma_ion_O_4S_cm2(1213) =  3.10 ; given_sigma_ion_O_2D_cm2(1213) =    0.00 ; given_sigma_ion_O_2P_cm2(1213) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1213) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1213) =    0.00
  given_sigma_ion_O_4S_cm2(1214) =  3.10 ; given_sigma_ion_O_2D_cm2(1214) =    0.00 ; given_sigma_ion_O_2P_cm2(1214) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1214) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1214) =    0.00
  given_sigma_ion_O_4S_cm2(1215) =  3.10 ; given_sigma_ion_O_2D_cm2(1215) =    0.00 ; given_sigma_ion_O_2P_cm2(1215) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1215) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1215) =    0.00
  given_sigma_ion_O_4S_cm2(1216) =  3.10 ; given_sigma_ion_O_2D_cm2(1216) =    0.00 ; given_sigma_ion_O_2P_cm2(1216) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1216) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1216) =    0.00
  given_sigma_ion_O_4S_cm2(1217) =  3.10 ; given_sigma_ion_O_2D_cm2(1217) =    0.00 ; given_sigma_ion_O_2P_cm2(1217) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1217) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1217) =    0.00
  given_sigma_ion_O_4S_cm2(1218) =  3.10 ; given_sigma_ion_O_2D_cm2(1218) =    0.00 ; given_sigma_ion_O_2P_cm2(1218) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1218) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1218) =    0.00
  given_sigma_ion_O_4S_cm2(1219) =  3.10 ; given_sigma_ion_O_2D_cm2(1219) =    0.00 ; given_sigma_ion_O_2P_cm2(1219) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1219) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1219) =    0.00
  given_sigma_ion_O_4S_cm2(1220) =  3.10 ; given_sigma_ion_O_2D_cm2(1220) =    0.00 ; given_sigma_ion_O_2P_cm2(1220) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1220) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1220) =    0.00
  given_sigma_ion_O_4S_cm2(1221) =  3.10 ; given_sigma_ion_O_2D_cm2(1221) =    0.00 ; given_sigma_ion_O_2P_cm2(1221) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1221) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1221) =    0.00
  given_sigma_ion_O_4S_cm2(1222) =  3.10 ; given_sigma_ion_O_2D_cm2(1222) =    0.00 ; given_sigma_ion_O_2P_cm2(1222) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1222) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1222) =    0.00
  given_sigma_ion_O_4S_cm2(1223) =  3.10 ; given_sigma_ion_O_2D_cm2(1223) =    0.00 ; given_sigma_ion_O_2P_cm2(1223) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1223) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1223) =    0.00
  given_sigma_ion_O_4S_cm2(1224) =  3.10 ; given_sigma_ion_O_2D_cm2(1224) =    0.00 ; given_sigma_ion_O_2P_cm2(1224) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1224) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1224) =    0.00
  given_sigma_ion_O_4S_cm2(1225) =  3.10 ; given_sigma_ion_O_2D_cm2(1225) =    0.00 ; given_sigma_ion_O_2P_cm2(1225) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1225) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1225) =    0.00
  given_sigma_ion_O_4S_cm2(1226) =  3.10 ; given_sigma_ion_O_2D_cm2(1226) =    0.00 ; given_sigma_ion_O_2P_cm2(1226) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1226) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1226) =    0.00
  given_sigma_ion_O_4S_cm2(1227) =  3.10 ; given_sigma_ion_O_2D_cm2(1227) =    0.00 ; given_sigma_ion_O_2P_cm2(1227) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1227) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1227) =    0.00
  given_sigma_ion_O_4S_cm2(1228) =  3.10 ; given_sigma_ion_O_2D_cm2(1228) =    0.00 ; given_sigma_ion_O_2P_cm2(1228) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1228) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1228) =    0.00
  given_sigma_ion_O_4S_cm2(1229) =  3.10 ; given_sigma_ion_O_2D_cm2(1229) =    0.00 ; given_sigma_ion_O_2P_cm2(1229) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1229) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1229) =    0.00
  given_sigma_ion_O_4S_cm2(1230) =  3.10 ; given_sigma_ion_O_2D_cm2(1230) =    0.00 ; given_sigma_ion_O_2P_cm2(1230) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1230) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1230) =    0.00
  given_sigma_ion_O_4S_cm2(1231) =  3.10 ; given_sigma_ion_O_2D_cm2(1231) =    0.00 ; given_sigma_ion_O_2P_cm2(1231) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1231) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1231) =    0.00
  given_sigma_ion_O_4S_cm2(1232) =  3.10 ; given_sigma_ion_O_2D_cm2(1232) =    0.00 ; given_sigma_ion_O_2P_cm2(1232) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1232) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1232) =    0.00
  given_sigma_ion_O_4S_cm2(1233) =  3.10 ; given_sigma_ion_O_2D_cm2(1233) =    0.00 ; given_sigma_ion_O_2P_cm2(1233) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1233) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1233) =    0.00
  given_sigma_ion_O_4S_cm2(1234) =  3.10 ; given_sigma_ion_O_2D_cm2(1234) =    0.00 ; given_sigma_ion_O_2P_cm2(1234) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1234) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1234) =    0.00
  given_sigma_ion_O_4S_cm2(1235) =  3.10 ; given_sigma_ion_O_2D_cm2(1235) =    0.00 ; given_sigma_ion_O_2P_cm2(1235) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1235) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1235) =    0.00
  given_sigma_ion_O_4S_cm2(1236) =  3.10 ; given_sigma_ion_O_2D_cm2(1236) =    0.00 ; given_sigma_ion_O_2P_cm2(1236) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1236) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1236) =    0.00
  given_sigma_ion_O_4S_cm2(1237) =  3.10 ; given_sigma_ion_O_2D_cm2(1237) =    0.00 ; given_sigma_ion_O_2P_cm2(1237) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1237) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1237) =    0.00
  given_sigma_ion_O_4S_cm2(1238) =  3.10 ; given_sigma_ion_O_2D_cm2(1238) =    0.00 ; given_sigma_ion_O_2P_cm2(1238) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1238) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1238) =    0.00
  given_sigma_ion_O_4S_cm2(1239) =  3.10 ; given_sigma_ion_O_2D_cm2(1239) =    0.00 ; given_sigma_ion_O_2P_cm2(1239) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1239) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1239) =    0.00
  given_sigma_ion_O_4S_cm2(1240) =  3.10 ; given_sigma_ion_O_2D_cm2(1240) =    0.00 ; given_sigma_ion_O_2P_cm2(1240) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1240) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1240) =    0.00
  given_sigma_ion_O_4S_cm2(1241) =  3.77 ; given_sigma_ion_O_2D_cm2(1241) =    0.00 ; given_sigma_ion_O_2P_cm2(1241) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1241) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1241) =    0.00
  given_sigma_ion_O_4S_cm2(1242) =  3.94 ; given_sigma_ion_O_2D_cm2(1242) =    0.00 ; given_sigma_ion_O_2P_cm2(1242) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1242) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1242) =    0.00
  given_sigma_ion_O_4S_cm2(1243) =  4.27 ; given_sigma_ion_O_2D_cm2(1243) =    0.00 ; given_sigma_ion_O_2P_cm2(1243) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1243) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1243) =    0.00
  given_sigma_ion_O_4S_cm2(1244) =  4.77 ; given_sigma_ion_O_2D_cm2(1244) =    0.00 ; given_sigma_ion_O_2P_cm2(1244) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1244) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1244) =    0.00
  given_sigma_ion_O_4S_cm2(1245) =  5.19 ; given_sigma_ion_O_2D_cm2(1245) =    0.00 ; given_sigma_ion_O_2P_cm2(1245) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1245) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1245) =    0.00
  given_sigma_ion_O_4S_cm2(1246) =  5.61 ; given_sigma_ion_O_2D_cm2(1246) =    0.00 ; given_sigma_ion_O_2P_cm2(1246) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1246) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1246) =    0.00
  given_sigma_ion_O_4S_cm2(1247) =  6.02 ; given_sigma_ion_O_2D_cm2(1247) =    0.00 ; given_sigma_ion_O_2P_cm2(1247) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1247) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1247) =    0.00
  given_sigma_ion_O_4S_cm2(1248) =  6.44 ; given_sigma_ion_O_2D_cm2(1248) =    0.00 ; given_sigma_ion_O_2P_cm2(1248) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1248) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1248) =    0.00
  given_sigma_ion_O_4S_cm2(1249) =  6.86 ; given_sigma_ion_O_2D_cm2(1249) =    0.00 ; given_sigma_ion_O_2P_cm2(1249) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1249) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1249) =    0.00
  given_sigma_ion_O_4S_cm2(1250) =  7.28 ; given_sigma_ion_O_2D_cm2(1250) =    0.00 ; given_sigma_ion_O_2P_cm2(1250) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1250) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1250) =    0.00
  given_sigma_ion_O_4S_cm2(1251) =  7.44 ; given_sigma_ion_O_2D_cm2(1251) =    0.00 ; given_sigma_ion_O_2P_cm2(1251) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1251) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1251) =    0.00
  given_sigma_ion_O_4S_cm2(1252) =  7.70 ; given_sigma_ion_O_2D_cm2(1252) =    0.00 ; given_sigma_ion_O_2P_cm2(1252) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1252) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1252) =    0.00
  given_sigma_ion_O_4S_cm2(1253) =  8.11 ; given_sigma_ion_O_2D_cm2(1253) =    0.00 ; given_sigma_ion_O_2P_cm2(1253) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1253) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1253) =    0.00
  given_sigma_ion_O_4S_cm2(1254) =  8.28 ; given_sigma_ion_O_2D_cm2(1254) =    0.00 ; given_sigma_ion_O_2P_cm2(1254) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1254) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1254) =    0.00
  given_sigma_ion_O_4S_cm2(1255) =  8.53 ; given_sigma_ion_O_2D_cm2(1255) =    0.00 ; given_sigma_ion_O_2P_cm2(1255) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1255) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1255) =    0.00
  given_sigma_ion_O_4S_cm2(1256) =  8.95 ; given_sigma_ion_O_2D_cm2(1256) =    0.00 ; given_sigma_ion_O_2P_cm2(1256) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1256) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1256) =    0.00
  given_sigma_ion_O_4S_cm2(1257) =  9.37 ; given_sigma_ion_O_2D_cm2(1257) =    0.00 ; given_sigma_ion_O_2P_cm2(1257) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1257) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1257) =    0.00
  given_sigma_ion_O_4S_cm2(1258) = 11.82 ; given_sigma_ion_O_2D_cm2(1258) =    0.00 ; given_sigma_ion_O_2P_cm2(1258) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1258) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1258) =    0.00
  given_sigma_ion_O_4S_cm2(1259) = 34.88 ; given_sigma_ion_O_2D_cm2(1259) =    0.00 ; given_sigma_ion_O_2P_cm2(1259) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1259) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1259) =    0.00
  given_sigma_ion_O_4S_cm2(1260) =  9.93 ; given_sigma_ion_O_2D_cm2(1260) =    0.00 ; given_sigma_ion_O_2P_cm2(1260) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1260) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1260) =    0.00
  given_sigma_ion_O_4S_cm2(1261) =  4.73 ; given_sigma_ion_O_2D_cm2(1261) =    0.00 ; given_sigma_ion_O_2P_cm2(1261) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1261) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1261) =    0.00
  given_sigma_ion_O_4S_cm2(1262) =  6.81 ; given_sigma_ion_O_2D_cm2(1262) =    0.00 ; given_sigma_ion_O_2P_cm2(1262) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1262) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1262) =    0.00
  given_sigma_ion_O_4S_cm2(1263) =  5.20 ; given_sigma_ion_O_2D_cm2(1263) =    0.00 ; given_sigma_ion_O_2P_cm2(1263) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1263) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1263) =    0.00
  given_sigma_ion_O_4S_cm2(1264) =  5.11 ; given_sigma_ion_O_2D_cm2(1264) =    0.00 ; given_sigma_ion_O_2P_cm2(1264) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1264) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1264) =    0.00
  given_sigma_ion_O_4S_cm2(1265) =  4.88 ; given_sigma_ion_O_2D_cm2(1265) =    0.00 ; given_sigma_ion_O_2P_cm2(1265) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1265) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1265) =    0.00
  given_sigma_ion_O_4S_cm2(1266) =  4.66 ; given_sigma_ion_O_2D_cm2(1266) =    0.00 ; given_sigma_ion_O_2P_cm2(1266) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1266) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1266) =    0.00
  given_sigma_ion_O_4S_cm2(1267) =  4.43 ; given_sigma_ion_O_2D_cm2(1267) =    0.00 ; given_sigma_ion_O_2P_cm2(1267) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1267) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1267) =    0.00
  given_sigma_ion_O_4S_cm2(1268) =  4.21 ; given_sigma_ion_O_2D_cm2(1268) =    0.00 ; given_sigma_ion_O_2P_cm2(1268) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1268) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1268) =    0.00
  given_sigma_ion_O_4S_cm2(1269) =  3.98 ; given_sigma_ion_O_2D_cm2(1269) =    0.00 ; given_sigma_ion_O_2P_cm2(1269) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1269) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1269) =    0.00
  given_sigma_ion_O_4S_cm2(1270) =  3.89 ; given_sigma_ion_O_2D_cm2(1270) =    0.00 ; given_sigma_ion_O_2P_cm2(1270) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1270) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1270) =    0.00
  given_sigma_ion_O_4S_cm2(1271) =  3.75 ; given_sigma_ion_O_2D_cm2(1271) =    0.00 ; given_sigma_ion_O_2P_cm2(1271) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1271) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1271) =    0.00
  given_sigma_ion_O_4S_cm2(1272) =  3.62 ; given_sigma_ion_O_2D_cm2(1272) =    0.00 ; given_sigma_ion_O_2P_cm2(1272) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1272) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1272) =    0.00
  given_sigma_ion_O_4S_cm2(1273) =  3.53 ; given_sigma_ion_O_2D_cm2(1273) =    0.00 ; given_sigma_ion_O_2P_cm2(1273) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1273) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1273) =    0.00
  given_sigma_ion_O_4S_cm2(1274) =  3.30 ; given_sigma_ion_O_2D_cm2(1274) =    0.00 ; given_sigma_ion_O_2P_cm2(1274) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1274) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1274) =    0.00
  given_sigma_ion_O_4S_cm2(1275) =  3.08 ; given_sigma_ion_O_2D_cm2(1275) =    0.00 ; given_sigma_ion_O_2P_cm2(1275) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1275) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1275) =    0.00
  given_sigma_ion_O_4S_cm2(1276) =  2.85 ; given_sigma_ion_O_2D_cm2(1276) =    0.00 ; given_sigma_ion_O_2P_cm2(1276) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1276) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1276) =    0.00
  given_sigma_ion_O_4S_cm2(1277) =  2.84 ; given_sigma_ion_O_2D_cm2(1277) =    0.00 ; given_sigma_ion_O_2P_cm2(1277) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1277) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1277) =    0.00
  given_sigma_ion_O_4S_cm2(1278) =  2.83 ; given_sigma_ion_O_2D_cm2(1278) =    0.00 ; given_sigma_ion_O_2P_cm2(1278) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1278) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1278) =    0.00
  given_sigma_ion_O_4S_cm2(1279) =  2.82 ; given_sigma_ion_O_2D_cm2(1279) =    0.00 ; given_sigma_ion_O_2P_cm2(1279) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1279) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1279) =    0.00
  given_sigma_ion_O_4S_cm2(1280) =  2.81 ; given_sigma_ion_O_2D_cm2(1280) =    0.00 ; given_sigma_ion_O_2P_cm2(1280) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1280) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1280) =    0.00
  given_sigma_ion_O_4S_cm2(1281) =  2.79 ; given_sigma_ion_O_2D_cm2(1281) =    0.00 ; given_sigma_ion_O_2P_cm2(1281) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1281) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1281) =    0.00
  given_sigma_ion_O_4S_cm2(1282) =  2.78 ; given_sigma_ion_O_2D_cm2(1282) =    0.00 ; given_sigma_ion_O_2P_cm2(1282) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1282) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1282) =    0.00
  given_sigma_ion_O_4S_cm2(1283) =  2.76 ; given_sigma_ion_O_2D_cm2(1283) =    0.00 ; given_sigma_ion_O_2P_cm2(1283) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1283) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1283) =    0.00
  given_sigma_ion_O_4S_cm2(1284) =  2.75 ; given_sigma_ion_O_2D_cm2(1284) =    0.00 ; given_sigma_ion_O_2P_cm2(1284) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1284) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1284) =    0.00
  given_sigma_ion_O_4S_cm2(1285) =  2.73 ; given_sigma_ion_O_2D_cm2(1285) =    0.00 ; given_sigma_ion_O_2P_cm2(1285) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1285) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1285) =    0.00
  given_sigma_ion_O_4S_cm2(1286) =  2.72 ; given_sigma_ion_O_2D_cm2(1286) =    0.00 ; given_sigma_ion_O_2P_cm2(1286) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1286) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1286) =    0.00
  given_sigma_ion_O_4S_cm2(1287) =  2.70 ; given_sigma_ion_O_2D_cm2(1287) =    0.00 ; given_sigma_ion_O_2P_cm2(1287) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1287) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1287) =    0.00
  given_sigma_ion_O_4S_cm2(1288) =  2.71 ; given_sigma_ion_O_2D_cm2(1288) =    0.00 ; given_sigma_ion_O_2P_cm2(1288) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1288) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1288) =    0.00
  given_sigma_ion_O_4S_cm2(1289) =  2.71 ; given_sigma_ion_O_2D_cm2(1289) =    0.00 ; given_sigma_ion_O_2P_cm2(1289) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1289) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1289) =    0.00
  given_sigma_ion_O_4S_cm2(1290) =  2.72 ; given_sigma_ion_O_2D_cm2(1290) =    0.00 ; given_sigma_ion_O_2P_cm2(1290) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1290) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1290) =    0.00
  given_sigma_ion_O_4S_cm2(1291) =  2.72 ; given_sigma_ion_O_2D_cm2(1291) =    0.00 ; given_sigma_ion_O_2P_cm2(1291) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1291) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1291) =    0.00
  given_sigma_ion_O_4S_cm2(1292) =  2.73 ; given_sigma_ion_O_2D_cm2(1292) =    0.00 ; given_sigma_ion_O_2P_cm2(1292) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1292) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1292) =    0.00
  given_sigma_ion_O_4S_cm2(1293) =  2.73 ; given_sigma_ion_O_2D_cm2(1293) =    0.00 ; given_sigma_ion_O_2P_cm2(1293) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1293) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1293) =    0.00
  given_sigma_ion_O_4S_cm2(1294) =  2.74 ; given_sigma_ion_O_2D_cm2(1294) =    0.00 ; given_sigma_ion_O_2P_cm2(1294) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1294) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1294) =    0.00
  given_sigma_ion_O_4S_cm2(1295) =  2.74 ; given_sigma_ion_O_2D_cm2(1295) =    0.00 ; given_sigma_ion_O_2P_cm2(1295) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1295) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1295) =    0.00
  given_sigma_ion_O_4S_cm2(1296) =  2.75 ; given_sigma_ion_O_2D_cm2(1296) =    0.00 ; given_sigma_ion_O_2P_cm2(1296) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1296) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1296) =    0.00
  given_sigma_ion_O_4S_cm2(1297) =  2.75 ; given_sigma_ion_O_2D_cm2(1297) =    0.00 ; given_sigma_ion_O_2P_cm2(1297) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1297) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1297) =    0.00
  given_sigma_ion_O_4S_cm2(1298) =  2.75 ; given_sigma_ion_O_2D_cm2(1298) =    0.00 ; given_sigma_ion_O_2P_cm2(1298) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1298) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1298) =    0.00
  given_sigma_ion_O_4S_cm2(1299) =  2.76 ; given_sigma_ion_O_2D_cm2(1299) =    0.00 ; given_sigma_ion_O_2P_cm2(1299) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1299) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1299) =    0.00
  given_sigma_ion_O_4S_cm2(1300) =  2.77 ; given_sigma_ion_O_2D_cm2(1300) =    0.00 ; given_sigma_ion_O_2P_cm2(1300) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1300) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1300) =    0.00
  given_sigma_ion_O_4S_cm2(1301) =  2.77 ; given_sigma_ion_O_2D_cm2(1301) =    0.00 ; given_sigma_ion_O_2P_cm2(1301) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1301) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1301) =    0.00
  given_sigma_ion_O_4S_cm2(1302) =  2.78 ; given_sigma_ion_O_2D_cm2(1302) =    0.00 ; given_sigma_ion_O_2P_cm2(1302) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1302) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1302) =    0.00
  given_sigma_ion_O_4S_cm2(1303) =  2.79 ; given_sigma_ion_O_2D_cm2(1303) =    0.00 ; given_sigma_ion_O_2P_cm2(1303) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1303) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1303) =    0.00
  given_sigma_ion_O_4S_cm2(1304) =  2.79 ; given_sigma_ion_O_2D_cm2(1304) =    0.00 ; given_sigma_ion_O_2P_cm2(1304) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1304) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1304) =    0.00
  given_sigma_ion_O_4S_cm2(1305) =  2.80 ; given_sigma_ion_O_2D_cm2(1305) =    0.00 ; given_sigma_ion_O_2P_cm2(1305) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1305) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1305) =    0.00
  given_sigma_ion_O_4S_cm2(1306) =  2.81 ; given_sigma_ion_O_2D_cm2(1306) =    0.00 ; given_sigma_ion_O_2P_cm2(1306) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1306) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1306) =    0.00
  given_sigma_ion_O_4S_cm2(1307) =  2.82 ; given_sigma_ion_O_2D_cm2(1307) =    0.00 ; given_sigma_ion_O_2P_cm2(1307) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1307) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1307) =    0.00
  given_sigma_ion_O_4S_cm2(1308) =  2.82 ; given_sigma_ion_O_2D_cm2(1308) =    0.00 ; given_sigma_ion_O_2P_cm2(1308) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1308) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1308) =    0.00
  given_sigma_ion_O_4S_cm2(1309) =  2.83 ; given_sigma_ion_O_2D_cm2(1309) =    0.00 ; given_sigma_ion_O_2P_cm2(1309) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1309) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1309) =    0.00
  given_sigma_ion_O_4S_cm2(1310) =  2.84 ; given_sigma_ion_O_2D_cm2(1310) =    0.00 ; given_sigma_ion_O_2P_cm2(1310) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1310) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1310) =    0.00
  given_sigma_ion_O_4S_cm2(1311) =  2.85 ; given_sigma_ion_O_2D_cm2(1311) =    0.00 ; given_sigma_ion_O_2P_cm2(1311) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1311) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1311) =    0.00
  given_sigma_ion_O_4S_cm2(1312) =  2.85 ; given_sigma_ion_O_2D_cm2(1312) =    0.00 ; given_sigma_ion_O_2P_cm2(1312) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1312) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1312) =    0.00
  given_sigma_ion_O_4S_cm2(1313) =  2.84 ; given_sigma_ion_O_2D_cm2(1313) =    0.00 ; given_sigma_ion_O_2P_cm2(1313) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1313) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1313) =    0.00
  given_sigma_ion_O_4S_cm2(1314) =  2.83 ; given_sigma_ion_O_2D_cm2(1314) =    0.00 ; given_sigma_ion_O_2P_cm2(1314) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1314) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1314) =    0.00
  given_sigma_ion_O_4S_cm2(1315) =  2.82 ; given_sigma_ion_O_2D_cm2(1315) =    0.00 ; given_sigma_ion_O_2P_cm2(1315) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1315) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1315) =    0.00
  given_sigma_ion_O_4S_cm2(1316) =  2.82 ; given_sigma_ion_O_2D_cm2(1316) =    0.00 ; given_sigma_ion_O_2P_cm2(1316) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1316) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1316) =    0.00
  given_sigma_ion_O_4S_cm2(1317) =  2.81 ; given_sigma_ion_O_2D_cm2(1317) =    0.00 ; given_sigma_ion_O_2P_cm2(1317) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1317) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1317) =    0.00
  given_sigma_ion_O_4S_cm2(1318) =  2.79 ; given_sigma_ion_O_2D_cm2(1318) =    0.00 ; given_sigma_ion_O_2P_cm2(1318) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1318) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1318) =    0.00
  given_sigma_ion_O_4S_cm2(1319) =  2.78 ; given_sigma_ion_O_2D_cm2(1319) =    0.00 ; given_sigma_ion_O_2P_cm2(1319) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1319) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1319) =    0.00
  given_sigma_ion_O_4S_cm2(1320) =  2.77 ; given_sigma_ion_O_2D_cm2(1320) =    0.00 ; given_sigma_ion_O_2P_cm2(1320) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1320) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1320) =    0.00
  given_sigma_ion_O_4S_cm2(1321) =  2.77 ; given_sigma_ion_O_2D_cm2(1321) =    0.00 ; given_sigma_ion_O_2P_cm2(1321) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1321) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1321) =    0.00
  given_sigma_ion_O_4S_cm2(1322) =  2.76 ; given_sigma_ion_O_2D_cm2(1322) =    0.00 ; given_sigma_ion_O_2P_cm2(1322) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1322) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1322) =    0.00
  given_sigma_ion_O_4S_cm2(1323) =  2.75 ; given_sigma_ion_O_2D_cm2(1323) =    0.00 ; given_sigma_ion_O_2P_cm2(1323) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1323) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1323) =    0.00
  given_sigma_ion_O_4S_cm2(1324) =  2.75 ; given_sigma_ion_O_2D_cm2(1324) =    0.00 ; given_sigma_ion_O_2P_cm2(1324) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1324) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1324) =    0.00
  given_sigma_ion_O_4S_cm2(1325) =  2.75 ; given_sigma_ion_O_2D_cm2(1325) =    0.00 ; given_sigma_ion_O_2P_cm2(1325) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1325) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1325) =    0.00
  given_sigma_ion_O_4S_cm2(1326) =  2.75 ; given_sigma_ion_O_2D_cm2(1326) =    0.00 ; given_sigma_ion_O_2P_cm2(1326) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1326) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1326) =    0.00
  given_sigma_ion_O_4S_cm2(1327) =  2.75 ; given_sigma_ion_O_2D_cm2(1327) =    0.00 ; given_sigma_ion_O_2P_cm2(1327) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1327) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1327) =    0.00
  given_sigma_ion_O_4S_cm2(1328) =  2.75 ; given_sigma_ion_O_2D_cm2(1328) =    0.00 ; given_sigma_ion_O_2P_cm2(1328) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1328) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1328) =    0.00
  given_sigma_ion_O_4S_cm2(1329) =  2.68 ; given_sigma_ion_O_2D_cm2(1329) =    0.00 ; given_sigma_ion_O_2P_cm2(1329) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1329) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1329) =    0.00
  given_sigma_ion_O_4S_cm2(1330) =  2.62 ; given_sigma_ion_O_2D_cm2(1330) =    0.00 ; given_sigma_ion_O_2P_cm2(1330) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1330) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1330) =    0.00
  given_sigma_ion_O_4S_cm2(1331) =  2.55 ; given_sigma_ion_O_2D_cm2(1331) =    0.00 ; given_sigma_ion_O_2P_cm2(1331) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1331) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1331) =    0.00
  given_sigma_ion_O_4S_cm2(1332) =  2.49 ; given_sigma_ion_O_2D_cm2(1332) =    0.00 ; given_sigma_ion_O_2P_cm2(1332) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1332) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1332) =    0.00
  given_sigma_ion_O_4S_cm2(1333) =  2.45 ; given_sigma_ion_O_2D_cm2(1333) =    0.00 ; given_sigma_ion_O_2P_cm2(1333) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1333) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1333) =    0.00
  given_sigma_ion_O_4S_cm2(1334) =  1.15 ; given_sigma_ion_O_2D_cm2(1334) =    0.00 ; given_sigma_ion_O_2P_cm2(1334) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1334) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1334) =    0.00
  given_sigma_ion_O_4S_cm2(1335) =  0.75 ; given_sigma_ion_O_2D_cm2(1335) =    0.00 ; given_sigma_ion_O_2P_cm2(1335) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1335) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1335) =    0.00
  given_sigma_ion_O_4S_cm2(1336) =  0.49 ; given_sigma_ion_O_2D_cm2(1336) =    0.00 ; given_sigma_ion_O_2P_cm2(1336) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1336) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1336) =    0.00
  given_sigma_ion_O_4S_cm2(1337) =  0.30 ; given_sigma_ion_O_2D_cm2(1337) =    0.00 ; given_sigma_ion_O_2P_cm2(1337) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1337) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1337) =    0.00
  given_sigma_ion_O_4S_cm2(1338) =  0.20 ; given_sigma_ion_O_2D_cm2(1338) =    0.00 ; given_sigma_ion_O_2P_cm2(1338) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1338) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1338) =    0.00
  given_sigma_ion_O_4S_cm2(1339) =  0.14 ; given_sigma_ion_O_2D_cm2(1339) =    0.00 ; given_sigma_ion_O_2P_cm2(1339) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1339) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1339) =    0.00
  given_sigma_ion_O_4S_cm2(1340) =  0.06 ; given_sigma_ion_O_2D_cm2(1340) =    0.00 ; given_sigma_ion_O_2P_cm2(1340) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1340) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1340) =    0.00
  given_sigma_ion_O_4S_cm2(1341) =  0.04 ; given_sigma_ion_O_2D_cm2(1341) =    0.00 ; given_sigma_ion_O_2P_cm2(1341) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1341) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1341) =    0.00
  given_sigma_ion_O_4S_cm2(1342) =  0.00 ; given_sigma_ion_O_2D_cm2(1342) =    0.00 ; given_sigma_ion_O_2P_cm2(1342) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1342) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1342) =    0.00
  given_sigma_ion_O_4S_cm2(1343) =  0.00 ; given_sigma_ion_O_2D_cm2(1343) =    0.00 ; given_sigma_ion_O_2P_cm2(1343) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1343) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1343) =    0.00
  given_sigma_ion_O_4S_cm2(1344) =  0.00 ; given_sigma_ion_O_2D_cm2(1344) =    0.00 ; given_sigma_ion_O_2P_cm2(1344) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1344) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1344) =    0.00
  given_sigma_ion_O_4S_cm2(1345) =  0.00 ; given_sigma_ion_O_2D_cm2(1345) =    0.00 ; given_sigma_ion_O_2P_cm2(1345) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1345) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1345) =    0.00
  given_sigma_ion_O_4S_cm2(1346) =  0.00 ; given_sigma_ion_O_2D_cm2(1346) =    0.00 ; given_sigma_ion_O_2P_cm2(1346) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1346) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1346) =    0.00
  given_sigma_ion_O_4S_cm2(1347) =  0.00 ; given_sigma_ion_O_2D_cm2(1347) =    0.00 ; given_sigma_ion_O_2P_cm2(1347) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1347) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1347) =    0.00
  given_sigma_ion_O_4S_cm2(1348) =  0.00 ; given_sigma_ion_O_2D_cm2(1348) =    0.00 ; given_sigma_ion_O_2P_cm2(1348) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1348) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1348) =    0.00
  given_sigma_ion_O_4S_cm2(1349) =  0.00 ; given_sigma_ion_O_2D_cm2(1349) =    0.00 ; given_sigma_ion_O_2P_cm2(1349) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1349) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1349) =    0.00
  given_sigma_ion_O_4S_cm2(1350) =  0.00 ; given_sigma_ion_O_2D_cm2(1350) =    0.00 ; given_sigma_ion_O_2P_cm2(1350) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1350) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1350) =    0.00
  given_sigma_ion_O_4S_cm2(1351) =  0.00 ; given_sigma_ion_O_2D_cm2(1351) =    0.00 ; given_sigma_ion_O_2P_cm2(1351) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1351) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1351) =    0.00
  given_sigma_ion_O_4S_cm2(1352) =  0.00 ; given_sigma_ion_O_2D_cm2(1352) =    0.00 ; given_sigma_ion_O_2P_cm2(1352) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1352) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1352) =    0.00
  given_sigma_ion_O_4S_cm2(1353) =  0.00 ; given_sigma_ion_O_2D_cm2(1353) =    0.00 ; given_sigma_ion_O_2P_cm2(1353) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1353) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1353) =    0.00
  given_sigma_ion_O_4S_cm2(1354) =  0.00 ; given_sigma_ion_O_2D_cm2(1354) =    0.00 ; given_sigma_ion_O_2P_cm2(1354) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1354) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1354) =    0.00
  given_sigma_ion_O_4S_cm2(1355) =  0.00 ; given_sigma_ion_O_2D_cm2(1355) =    0.00 ; given_sigma_ion_O_2P_cm2(1355) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1355) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1355) =    0.00
  given_sigma_ion_O_4S_cm2(1356) =  0.00 ; given_sigma_ion_O_2D_cm2(1356) =    0.00 ; given_sigma_ion_O_2P_cm2(1356) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1356) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1356) =    0.00
  given_sigma_ion_O_4S_cm2(1357) =  0.00 ; given_sigma_ion_O_2D_cm2(1357) =    0.00 ; given_sigma_ion_O_2P_cm2(1357) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1357) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1357) =    0.00
  given_sigma_ion_O_4S_cm2(1358) =  0.00 ; given_sigma_ion_O_2D_cm2(1358) =    0.00 ; given_sigma_ion_O_2P_cm2(1358) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1358) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1358) =    0.00
  given_sigma_ion_O_4S_cm2(1359) =  0.00 ; given_sigma_ion_O_2D_cm2(1359) =    0.00 ; given_sigma_ion_O_2P_cm2(1359) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1359) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1359) =    0.00
  given_sigma_ion_O_4S_cm2(1360) =  0.00 ; given_sigma_ion_O_2D_cm2(1360) =    0.00 ; given_sigma_ion_O_2P_cm2(1360) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1360) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1360) =    0.00
  given_sigma_ion_O_4S_cm2(1361) =  0.00 ; given_sigma_ion_O_2D_cm2(1361) =    0.00 ; given_sigma_ion_O_2P_cm2(1361) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1361) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1361) =    0.00
  given_sigma_ion_O_4S_cm2(1362) =  0.00 ; given_sigma_ion_O_2D_cm2(1362) =    0.00 ; given_sigma_ion_O_2P_cm2(1362) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1362) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1362) =    0.00
  given_sigma_ion_O_4S_cm2(1363) =  0.00 ; given_sigma_ion_O_2D_cm2(1363) =    0.00 ; given_sigma_ion_O_2P_cm2(1363) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1363) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1363) =    0.00
  given_sigma_ion_O_4S_cm2(1364) =  0.00 ; given_sigma_ion_O_2D_cm2(1364) =    0.00 ; given_sigma_ion_O_2P_cm2(1364) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1364) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1364) =    0.00
  given_sigma_ion_O_4S_cm2(1365) =  0.00 ; given_sigma_ion_O_2D_cm2(1365) =    0.00 ; given_sigma_ion_O_2P_cm2(1365) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1365) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1365) =    0.00
  given_sigma_ion_O_4S_cm2(1366) =  0.00 ; given_sigma_ion_O_2D_cm2(1366) =    0.00 ; given_sigma_ion_O_2P_cm2(1366) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1366) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1366) =    0.00
  given_sigma_ion_O_4S_cm2(1367) =  0.00 ; given_sigma_ion_O_2D_cm2(1367) =    0.00 ; given_sigma_ion_O_2P_cm2(1367) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1367) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1367) =    0.00
  given_sigma_ion_O_4S_cm2(1368) =  0.00 ; given_sigma_ion_O_2D_cm2(1368) =    0.00 ; given_sigma_ion_O_2P_cm2(1368) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1368) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1368) =    0.00
  given_sigma_ion_O_4S_cm2(1369) =  0.00 ; given_sigma_ion_O_2D_cm2(1369) =    0.00 ; given_sigma_ion_O_2P_cm2(1369) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1369) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1369) =    0.00
  given_sigma_ion_O_4S_cm2(1370) =  0.00 ; given_sigma_ion_O_2D_cm2(1370) =    0.00 ; given_sigma_ion_O_2P_cm2(1370) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1370) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1370) =    0.00
  given_sigma_ion_O_4S_cm2(1371) =  0.00 ; given_sigma_ion_O_2D_cm2(1371) =    0.00 ; given_sigma_ion_O_2P_cm2(1371) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1371) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1371) =    0.00
  given_sigma_ion_O_4S_cm2(1372) =  0.00 ; given_sigma_ion_O_2D_cm2(1372) =    0.00 ; given_sigma_ion_O_2P_cm2(1372) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1372) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1372) =    0.00
  given_sigma_ion_O_4S_cm2(1373) =  0.00 ; given_sigma_ion_O_2D_cm2(1373) =    0.00 ; given_sigma_ion_O_2P_cm2(1373) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1373) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1373) =    0.00
  given_sigma_ion_O_4S_cm2(1374) =  0.00 ; given_sigma_ion_O_2D_cm2(1374) =    0.00 ; given_sigma_ion_O_2P_cm2(1374) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1374) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1374) =    0.00
  given_sigma_ion_O_4S_cm2(1375) =  0.00 ; given_sigma_ion_O_2D_cm2(1375) =    0.00 ; given_sigma_ion_O_2P_cm2(1375) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1375) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1375) =    0.00
  given_sigma_ion_O_4S_cm2(1376) =  0.00 ; given_sigma_ion_O_2D_cm2(1376) =    0.00 ; given_sigma_ion_O_2P_cm2(1376) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1376) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1376) =    0.00
  given_sigma_ion_O_4S_cm2(1377) =  0.00 ; given_sigma_ion_O_2D_cm2(1377) =    0.00 ; given_sigma_ion_O_2P_cm2(1377) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1377) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1377) =    0.00
  given_sigma_ion_O_4S_cm2(1378) =  0.00 ; given_sigma_ion_O_2D_cm2(1378) =    0.00 ; given_sigma_ion_O_2P_cm2(1378) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1378) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1378) =    0.00
  given_sigma_ion_O_4S_cm2(1379) =  0.00 ; given_sigma_ion_O_2D_cm2(1379) =    0.00 ; given_sigma_ion_O_2P_cm2(1379) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1379) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1379) =    0.00
  given_sigma_ion_O_4S_cm2(1380) =  0.00 ; given_sigma_ion_O_2D_cm2(1380) =    0.00 ; given_sigma_ion_O_2P_cm2(1380) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1380) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1380) =    0.00
  given_sigma_ion_O_4S_cm2(1381) =  0.00 ; given_sigma_ion_O_2D_cm2(1381) =    0.00 ; given_sigma_ion_O_2P_cm2(1381) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1381) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1381) =    0.00
  given_sigma_ion_O_4S_cm2(1382) =  0.00 ; given_sigma_ion_O_2D_cm2(1382) =    0.00 ; given_sigma_ion_O_2P_cm2(1382) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1382) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1382) =    0.00
  given_sigma_ion_O_4S_cm2(1383) =  0.00 ; given_sigma_ion_O_2D_cm2(1383) =    0.00 ; given_sigma_ion_O_2P_cm2(1383) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1383) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1383) =    0.00
  given_sigma_ion_O_4S_cm2(1384) =  0.00 ; given_sigma_ion_O_2D_cm2(1384) =    0.00 ; given_sigma_ion_O_2P_cm2(1384) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1384) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1384) =    0.00
  given_sigma_ion_O_4S_cm2(1385) =  0.00 ; given_sigma_ion_O_2D_cm2(1385) =    0.00 ; given_sigma_ion_O_2P_cm2(1385) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1385) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1385) =    0.00
  given_sigma_ion_O_4S_cm2(1386) =  0.00 ; given_sigma_ion_O_2D_cm2(1386) =    0.00 ; given_sigma_ion_O_2P_cm2(1386) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1386) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1386) =    0.00
  given_sigma_ion_O_4S_cm2(1387) =  0.00 ; given_sigma_ion_O_2D_cm2(1387) =    0.00 ; given_sigma_ion_O_2P_cm2(1387) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1387) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1387) =    0.00
  given_sigma_ion_O_4S_cm2(1388) =  0.00 ; given_sigma_ion_O_2D_cm2(1388) =    0.00 ; given_sigma_ion_O_2P_cm2(1388) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1388) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1388) =    0.00
  given_sigma_ion_O_4S_cm2(1389) =  0.00 ; given_sigma_ion_O_2D_cm2(1389) =    0.00 ; given_sigma_ion_O_2P_cm2(1389) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1389) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1389) =    0.00
  given_sigma_ion_O_4S_cm2(1390) =  0.00 ; given_sigma_ion_O_2D_cm2(1390) =    0.00 ; given_sigma_ion_O_2P_cm2(1390) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1390) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1390) =    0.00
  given_sigma_ion_O_4S_cm2(1391) =  0.00 ; given_sigma_ion_O_2D_cm2(1391) =    0.00 ; given_sigma_ion_O_2P_cm2(1391) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1391) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1391) =    0.00
  given_sigma_ion_O_4S_cm2(1392) =  0.00 ; given_sigma_ion_O_2D_cm2(1392) =    0.00 ; given_sigma_ion_O_2P_cm2(1392) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1392) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1392) =    0.00
  given_sigma_ion_O_4S_cm2(1393) =  0.00 ; given_sigma_ion_O_2D_cm2(1393) =    0.00 ; given_sigma_ion_O_2P_cm2(1393) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1393) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1393) =    0.00
  given_sigma_ion_O_4S_cm2(1394) =  0.00 ; given_sigma_ion_O_2D_cm2(1394) =    0.00 ; given_sigma_ion_O_2P_cm2(1394) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1394) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1394) =    0.00
  given_sigma_ion_O_4S_cm2(1395) =  0.00 ; given_sigma_ion_O_2D_cm2(1395) =    0.00 ; given_sigma_ion_O_2P_cm2(1395) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1395) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1395) =    0.00
  given_sigma_ion_O_4S_cm2(1396) =  0.00 ; given_sigma_ion_O_2D_cm2(1396) =    0.00 ; given_sigma_ion_O_2P_cm2(1396) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1396) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1396) =    0.00
  given_sigma_ion_O_4S_cm2(1397) =  0.00 ; given_sigma_ion_O_2D_cm2(1397) =    0.00 ; given_sigma_ion_O_2P_cm2(1397) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1397) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1397) =    0.00
  given_sigma_ion_O_4S_cm2(1398) =  0.00 ; given_sigma_ion_O_2D_cm2(1398) =    0.00 ; given_sigma_ion_O_2P_cm2(1398) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1398) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1398) =    0.00
  given_sigma_ion_O_4S_cm2(1399) =  0.00 ; given_sigma_ion_O_2D_cm2(1399) =    0.00 ; given_sigma_ion_O_2P_cm2(1399) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1399) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1399) =    0.00
  given_sigma_ion_O_4S_cm2(1400) =  0.00 ; given_sigma_ion_O_2D_cm2(1400) =    0.00 ; given_sigma_ion_O_2P_cm2(1400) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1400) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1400) =    0.00
  given_sigma_ion_O_4S_cm2(1401) =  0.00 ; given_sigma_ion_O_2D_cm2(1401) =    0.00 ; given_sigma_ion_O_2P_cm2(1401) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1401) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1401) =    0.00
  given_sigma_ion_O_4S_cm2(1402) =  0.00 ; given_sigma_ion_O_2D_cm2(1402) =    0.00 ; given_sigma_ion_O_2P_cm2(1402) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1402) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1402) =    0.00
  given_sigma_ion_O_4S_cm2(1403) =  0.00 ; given_sigma_ion_O_2D_cm2(1403) =    0.00 ; given_sigma_ion_O_2P_cm2(1403) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1403) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1403) =    0.00
  given_sigma_ion_O_4S_cm2(1404) =  0.00 ; given_sigma_ion_O_2D_cm2(1404) =    0.00 ; given_sigma_ion_O_2P_cm2(1404) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1404) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1404) =    0.00
  given_sigma_ion_O_4S_cm2(1405) =  0.00 ; given_sigma_ion_O_2D_cm2(1405) =    0.00 ; given_sigma_ion_O_2P_cm2(1405) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1405) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1405) =    0.00
  given_sigma_ion_O_4S_cm2(1406) =  0.00 ; given_sigma_ion_O_2D_cm2(1406) =    0.00 ; given_sigma_ion_O_2P_cm2(1406) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1406) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1406) =    0.00
  given_sigma_ion_O_4S_cm2(1407) =  0.00 ; given_sigma_ion_O_2D_cm2(1407) =    0.00 ; given_sigma_ion_O_2P_cm2(1407) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1407) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1407) =    0.00
  given_sigma_ion_O_4S_cm2(1408) =  0.00 ; given_sigma_ion_O_2D_cm2(1408) =    0.00 ; given_sigma_ion_O_2P_cm2(1408) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1408) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1408) =    0.00
  given_sigma_ion_O_4S_cm2(1409) =  0.00 ; given_sigma_ion_O_2D_cm2(1409) =    0.00 ; given_sigma_ion_O_2P_cm2(1409) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1409) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1409) =    0.00
  given_sigma_ion_O_4S_cm2(1410) =  0.00 ; given_sigma_ion_O_2D_cm2(1410) =    0.00 ; given_sigma_ion_O_2P_cm2(1410) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1410) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1410) =    0.00
  given_sigma_ion_O_4S_cm2(1411) =  0.00 ; given_sigma_ion_O_2D_cm2(1411) =    0.00 ; given_sigma_ion_O_2P_cm2(1411) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1411) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1411) =    0.00
  given_sigma_ion_O_4S_cm2(1412) =  0.00 ; given_sigma_ion_O_2D_cm2(1412) =    0.00 ; given_sigma_ion_O_2P_cm2(1412) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1412) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1412) =    0.00
  given_sigma_ion_O_4S_cm2(1413) =  0.00 ; given_sigma_ion_O_2D_cm2(1413) =    0.00 ; given_sigma_ion_O_2P_cm2(1413) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1413) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1413) =    0.00
  given_sigma_ion_O_4S_cm2(1414) =  0.00 ; given_sigma_ion_O_2D_cm2(1414) =    0.00 ; given_sigma_ion_O_2P_cm2(1414) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1414) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1414) =    0.00
  given_sigma_ion_O_4S_cm2(1415) =  0.00 ; given_sigma_ion_O_2D_cm2(1415) =    0.00 ; given_sigma_ion_O_2P_cm2(1415) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1415) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1415) =    0.00
  given_sigma_ion_O_4S_cm2(1416) =  0.00 ; given_sigma_ion_O_2D_cm2(1416) =    0.00 ; given_sigma_ion_O_2P_cm2(1416) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1416) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1416) =    0.00
  given_sigma_ion_O_4S_cm2(1417) =  0.00 ; given_sigma_ion_O_2D_cm2(1417) =    0.00 ; given_sigma_ion_O_2P_cm2(1417) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1417) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1417) =    0.00
  given_sigma_ion_O_4S_cm2(1418) =  0.00 ; given_sigma_ion_O_2D_cm2(1418) =    0.00 ; given_sigma_ion_O_2P_cm2(1418) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1418) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1418) =    0.00
  given_sigma_ion_O_4S_cm2(1419) =  0.00 ; given_sigma_ion_O_2D_cm2(1419) =    0.00 ; given_sigma_ion_O_2P_cm2(1419) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1419) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1419) =    0.00
  given_sigma_ion_O_4S_cm2(1420) =  0.00 ; given_sigma_ion_O_2D_cm2(1420) =    0.00 ; given_sigma_ion_O_2P_cm2(1420) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1420) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1420) =    0.00
  given_sigma_ion_O_4S_cm2(1421) =  0.00 ; given_sigma_ion_O_2D_cm2(1421) =    0.00 ; given_sigma_ion_O_2P_cm2(1421) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1421) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1421) =    0.00
  given_sigma_ion_O_4S_cm2(1422) =  0.00 ; given_sigma_ion_O_2D_cm2(1422) =    0.00 ; given_sigma_ion_O_2P_cm2(1422) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1422) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1422) =    0.00
  given_sigma_ion_O_4S_cm2(1423) =  0.00 ; given_sigma_ion_O_2D_cm2(1423) =    0.00 ; given_sigma_ion_O_2P_cm2(1423) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1423) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1423) =    0.00
  given_sigma_ion_O_4S_cm2(1424) =  0.00 ; given_sigma_ion_O_2D_cm2(1424) =    0.00 ; given_sigma_ion_O_2P_cm2(1424) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1424) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1424) =    0.00
  given_sigma_ion_O_4S_cm2(1425) =  0.00 ; given_sigma_ion_O_2D_cm2(1425) =    0.00 ; given_sigma_ion_O_2P_cm2(1425) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1425) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1425) =    0.00
  given_sigma_ion_O_4S_cm2(1426) =  0.00 ; given_sigma_ion_O_2D_cm2(1426) =    0.00 ; given_sigma_ion_O_2P_cm2(1426) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1426) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1426) =    0.00
  given_sigma_ion_O_4S_cm2(1427) =  0.00 ; given_sigma_ion_O_2D_cm2(1427) =    0.00 ; given_sigma_ion_O_2P_cm2(1427) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1427) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1427) =    0.00
  given_sigma_ion_O_4S_cm2(1428) =  0.00 ; given_sigma_ion_O_2D_cm2(1428) =    0.00 ; given_sigma_ion_O_2P_cm2(1428) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1428) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1428) =    0.00
  given_sigma_ion_O_4S_cm2(1429) =  0.00 ; given_sigma_ion_O_2D_cm2(1429) =    0.00 ; given_sigma_ion_O_2P_cm2(1429) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1429) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1429) =    0.00
  given_sigma_ion_O_4S_cm2(1430) =  0.00 ; given_sigma_ion_O_2D_cm2(1430) =    0.00 ; given_sigma_ion_O_2P_cm2(1430) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1430) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1430) =    0.00
  given_sigma_ion_O_4S_cm2(1431) =  0.00 ; given_sigma_ion_O_2D_cm2(1431) =    0.00 ; given_sigma_ion_O_2P_cm2(1431) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1431) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1431) =    0.00
  given_sigma_ion_O_4S_cm2(1432) =  0.00 ; given_sigma_ion_O_2D_cm2(1432) =    0.00 ; given_sigma_ion_O_2P_cm2(1432) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1432) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1432) =    0.00
  given_sigma_ion_O_4S_cm2(1433) =  0.00 ; given_sigma_ion_O_2D_cm2(1433) =    0.00 ; given_sigma_ion_O_2P_cm2(1433) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1433) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1433) =    0.00
  given_sigma_ion_O_4S_cm2(1434) =  0.00 ; given_sigma_ion_O_2D_cm2(1434) =    0.00 ; given_sigma_ion_O_2P_cm2(1434) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1434) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1434) =    0.00
  given_sigma_ion_O_4S_cm2(1435) =  0.00 ; given_sigma_ion_O_2D_cm2(1435) =    0.00 ; given_sigma_ion_O_2P_cm2(1435) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1435) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1435) =    0.00
  given_sigma_ion_O_4S_cm2(1436) =  0.00 ; given_sigma_ion_O_2D_cm2(1436) =    0.00 ; given_sigma_ion_O_2P_cm2(1436) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1436) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1436) =    0.00
  given_sigma_ion_O_4S_cm2(1437) =  0.00 ; given_sigma_ion_O_2D_cm2(1437) =    0.00 ; given_sigma_ion_O_2P_cm2(1437) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1437) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1437) =    0.00
  given_sigma_ion_O_4S_cm2(1438) =  0.00 ; given_sigma_ion_O_2D_cm2(1438) =    0.00 ; given_sigma_ion_O_2P_cm2(1438) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1438) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1438) =    0.00
  given_sigma_ion_O_4S_cm2(1439) =  0.00 ; given_sigma_ion_O_2D_cm2(1439) =    0.00 ; given_sigma_ion_O_2P_cm2(1439) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1439) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1439) =    0.00
  given_sigma_ion_O_4S_cm2(1440) =  0.00 ; given_sigma_ion_O_2D_cm2(1440) =    0.00 ; given_sigma_ion_O_2P_cm2(1440) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1440) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1440) =    0.00
  given_sigma_ion_O_4S_cm2(1441) =  0.00 ; given_sigma_ion_O_2D_cm2(1441) =    0.00 ; given_sigma_ion_O_2P_cm2(1441) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1441) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1441) =    0.00
  given_sigma_ion_O_4S_cm2(1442) =  0.00 ; given_sigma_ion_O_2D_cm2(1442) =    0.00 ; given_sigma_ion_O_2P_cm2(1442) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1442) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1442) =    0.00
  given_sigma_ion_O_4S_cm2(1443) =  0.00 ; given_sigma_ion_O_2D_cm2(1443) =    0.00 ; given_sigma_ion_O_2P_cm2(1443) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1443) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1443) =    0.00
  given_sigma_ion_O_4S_cm2(1444) =  0.00 ; given_sigma_ion_O_2D_cm2(1444) =    0.00 ; given_sigma_ion_O_2P_cm2(1444) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1444) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1444) =    0.00
  given_sigma_ion_O_4S_cm2(1445) =  0.00 ; given_sigma_ion_O_2D_cm2(1445) =    0.00 ; given_sigma_ion_O_2P_cm2(1445) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1445) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1445) =    0.00
  given_sigma_ion_O_4S_cm2(1446) =  0.00 ; given_sigma_ion_O_2D_cm2(1446) =    0.00 ; given_sigma_ion_O_2P_cm2(1446) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1446) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1446) =    0.00
  given_sigma_ion_O_4S_cm2(1447) =  0.00 ; given_sigma_ion_O_2D_cm2(1447) =    0.00 ; given_sigma_ion_O_2P_cm2(1447) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1447) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1447) =    0.00
  given_sigma_ion_O_4S_cm2(1448) =  0.00 ; given_sigma_ion_O_2D_cm2(1448) =    0.00 ; given_sigma_ion_O_2P_cm2(1448) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1448) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1448) =    0.00
  given_sigma_ion_O_4S_cm2(1449) =  0.00 ; given_sigma_ion_O_2D_cm2(1449) =    0.00 ; given_sigma_ion_O_2P_cm2(1449) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1449) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1449) =    0.00
  given_sigma_ion_O_4S_cm2(1450) =  0.00 ; given_sigma_ion_O_2D_cm2(1450) =    0.00 ; given_sigma_ion_O_2P_cm2(1450) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1450) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1450) =    0.00
  given_sigma_ion_O_4S_cm2(1451) =  0.00 ; given_sigma_ion_O_2D_cm2(1451) =    0.00 ; given_sigma_ion_O_2P_cm2(1451) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1451) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1451) =    0.00
  given_sigma_ion_O_4S_cm2(1452) =  0.00 ; given_sigma_ion_O_2D_cm2(1452) =    0.00 ; given_sigma_ion_O_2P_cm2(1452) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1452) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1452) =    0.00
  given_sigma_ion_O_4S_cm2(1453) =  0.00 ; given_sigma_ion_O_2D_cm2(1453) =    0.00 ; given_sigma_ion_O_2P_cm2(1453) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1453) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1453) =    0.00
  given_sigma_ion_O_4S_cm2(1454) =  0.00 ; given_sigma_ion_O_2D_cm2(1454) =    0.00 ; given_sigma_ion_O_2P_cm2(1454) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1454) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1454) =    0.00
  given_sigma_ion_O_4S_cm2(1455) =  0.00 ; given_sigma_ion_O_2D_cm2(1455) =    0.00 ; given_sigma_ion_O_2P_cm2(1455) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1455) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1455) =    0.00
  given_sigma_ion_O_4S_cm2(1456) =  0.00 ; given_sigma_ion_O_2D_cm2(1456) =    0.00 ; given_sigma_ion_O_2P_cm2(1456) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1456) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1456) =    0.00
  given_sigma_ion_O_4S_cm2(1457) =  0.00 ; given_sigma_ion_O_2D_cm2(1457) =    0.00 ; given_sigma_ion_O_2P_cm2(1457) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1457) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1457) =    0.00
  given_sigma_ion_O_4S_cm2(1458) =  0.00 ; given_sigma_ion_O_2D_cm2(1458) =    0.00 ; given_sigma_ion_O_2P_cm2(1458) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1458) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1458) =    0.00
  given_sigma_ion_O_4S_cm2(1459) =  0.00 ; given_sigma_ion_O_2D_cm2(1459) =    0.00 ; given_sigma_ion_O_2P_cm2(1459) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1459) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1459) =    0.00
  given_sigma_ion_O_4S_cm2(1460) =  0.00 ; given_sigma_ion_O_2D_cm2(1460) =    0.00 ; given_sigma_ion_O_2P_cm2(1460) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1460) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1460) =    0.00
  given_sigma_ion_O_4S_cm2(1461) =  0.00 ; given_sigma_ion_O_2D_cm2(1461) =    0.00 ; given_sigma_ion_O_2P_cm2(1461) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1461) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1461) =    0.00
  given_sigma_ion_O_4S_cm2(1462) =  0.00 ; given_sigma_ion_O_2D_cm2(1462) =    0.00 ; given_sigma_ion_O_2P_cm2(1462) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1462) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1462) =    0.00
  given_sigma_ion_O_4S_cm2(1463) =  0.00 ; given_sigma_ion_O_2D_cm2(1463) =    0.00 ; given_sigma_ion_O_2P_cm2(1463) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1463) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1463) =    0.00
  given_sigma_ion_O_4S_cm2(1464) =  0.00 ; given_sigma_ion_O_2D_cm2(1464) =    0.00 ; given_sigma_ion_O_2P_cm2(1464) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1464) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1464) =    0.00
  given_sigma_ion_O_4S_cm2(1465) =  0.00 ; given_sigma_ion_O_2D_cm2(1465) =    0.00 ; given_sigma_ion_O_2P_cm2(1465) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1465) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1465) =    0.00
  given_sigma_ion_O_4S_cm2(1466) =  0.00 ; given_sigma_ion_O_2D_cm2(1466) =    0.00 ; given_sigma_ion_O_2P_cm2(1466) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1466) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1466) =    0.00
  given_sigma_ion_O_4S_cm2(1467) =  0.00 ; given_sigma_ion_O_2D_cm2(1467) =    0.00 ; given_sigma_ion_O_2P_cm2(1467) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1467) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1467) =    0.00
  given_sigma_ion_O_4S_cm2(1468) =  0.00 ; given_sigma_ion_O_2D_cm2(1468) =    0.00 ; given_sigma_ion_O_2P_cm2(1468) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1468) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1468) =    0.00
  given_sigma_ion_O_4S_cm2(1469) =  0.00 ; given_sigma_ion_O_2D_cm2(1469) =    0.00 ; given_sigma_ion_O_2P_cm2(1469) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1469) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1469) =    0.00
  given_sigma_ion_O_4S_cm2(1470) =  0.00 ; given_sigma_ion_O_2D_cm2(1470) =    0.00 ; given_sigma_ion_O_2P_cm2(1470) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1470) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1470) =    0.00
  given_sigma_ion_O_4S_cm2(1471) =  0.00 ; given_sigma_ion_O_2D_cm2(1471) =    0.00 ; given_sigma_ion_O_2P_cm2(1471) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1471) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1471) =    0.00
  given_sigma_ion_O_4S_cm2(1472) =  0.00 ; given_sigma_ion_O_2D_cm2(1472) =    0.00 ; given_sigma_ion_O_2P_cm2(1472) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1472) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1472) =    0.00
  given_sigma_ion_O_4S_cm2(1473) =  0.00 ; given_sigma_ion_O_2D_cm2(1473) =    0.00 ; given_sigma_ion_O_2P_cm2(1473) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1473) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1473) =    0.00
  given_sigma_ion_O_4S_cm2(1474) =  0.00 ; given_sigma_ion_O_2D_cm2(1474) =    0.00 ; given_sigma_ion_O_2P_cm2(1474) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1474) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1474) =    0.00
  given_sigma_ion_O_4S_cm2(1475) =  0.00 ; given_sigma_ion_O_2D_cm2(1475) =    0.00 ; given_sigma_ion_O_2P_cm2(1475) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1475) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1475) =    0.00
  given_sigma_ion_O_4S_cm2(1476) =  0.00 ; given_sigma_ion_O_2D_cm2(1476) =    0.00 ; given_sigma_ion_O_2P_cm2(1476) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1476) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1476) =    0.00
  given_sigma_ion_O_4S_cm2(1477) =  0.00 ; given_sigma_ion_O_2D_cm2(1477) =    0.00 ; given_sigma_ion_O_2P_cm2(1477) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1477) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1477) =    0.00
  given_sigma_ion_O_4S_cm2(1478) =  0.00 ; given_sigma_ion_O_2D_cm2(1478) =    0.00 ; given_sigma_ion_O_2P_cm2(1478) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1478) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1478) =    0.00
  given_sigma_ion_O_4S_cm2(1479) =  0.00 ; given_sigma_ion_O_2D_cm2(1479) =    0.00 ; given_sigma_ion_O_2P_cm2(1479) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1479) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1479) =    0.00
  given_sigma_ion_O_4S_cm2(1480) =  0.00 ; given_sigma_ion_O_2D_cm2(1480) =    0.00 ; given_sigma_ion_O_2P_cm2(1480) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1480) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1480) =    0.00
  given_sigma_ion_O_4S_cm2(1481) =  0.00 ; given_sigma_ion_O_2D_cm2(1481) =    0.00 ; given_sigma_ion_O_2P_cm2(1481) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1481) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1481) =    0.00
  given_sigma_ion_O_4S_cm2(1482) =  0.00 ; given_sigma_ion_O_2D_cm2(1482) =    0.00 ; given_sigma_ion_O_2P_cm2(1482) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1482) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1482) =    0.00
  given_sigma_ion_O_4S_cm2(1483) =  0.00 ; given_sigma_ion_O_2D_cm2(1483) =    0.00 ; given_sigma_ion_O_2P_cm2(1483) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1483) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1483) =    0.00
  given_sigma_ion_O_4S_cm2(1484) =  0.00 ; given_sigma_ion_O_2D_cm2(1484) =    0.00 ; given_sigma_ion_O_2P_cm2(1484) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1484) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1484) =    0.00
  given_sigma_ion_O_4S_cm2(1485) =  0.00 ; given_sigma_ion_O_2D_cm2(1485) =    0.00 ; given_sigma_ion_O_2P_cm2(1485) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1485) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1485) =    0.00
  given_sigma_ion_O_4S_cm2(1486) =  0.00 ; given_sigma_ion_O_2D_cm2(1486) =    0.00 ; given_sigma_ion_O_2P_cm2(1486) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1486) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1486) =    0.00
  given_sigma_ion_O_4S_cm2(1487) =  0.00 ; given_sigma_ion_O_2D_cm2(1487) =    0.00 ; given_sigma_ion_O_2P_cm2(1487) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1487) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1487) =    0.00
  given_sigma_ion_O_4S_cm2(1488) =  0.00 ; given_sigma_ion_O_2D_cm2(1488) =    0.00 ; given_sigma_ion_O_2P_cm2(1488) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1488) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1488) =    0.00
  given_sigma_ion_O_4S_cm2(1489) =  0.00 ; given_sigma_ion_O_2D_cm2(1489) =    0.00 ; given_sigma_ion_O_2P_cm2(1489) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1489) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1489) =    0.00
  given_sigma_ion_O_4S_cm2(1490) =  0.00 ; given_sigma_ion_O_2D_cm2(1490) =    0.00 ; given_sigma_ion_O_2P_cm2(1490) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1490) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1490) =    0.00
  given_sigma_ion_O_4S_cm2(1491) =  0.00 ; given_sigma_ion_O_2D_cm2(1491) =    0.00 ; given_sigma_ion_O_2P_cm2(1491) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1491) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1491) =    0.00
  given_sigma_ion_O_4S_cm2(1492) =  0.00 ; given_sigma_ion_O_2D_cm2(1492) =    0.00 ; given_sigma_ion_O_2P_cm2(1492) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1492) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1492) =    0.00
  given_sigma_ion_O_4S_cm2(1493) =  0.00 ; given_sigma_ion_O_2D_cm2(1493) =    0.00 ; given_sigma_ion_O_2P_cm2(1493) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1493) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1493) =    0.00
  given_sigma_ion_O_4S_cm2(1494) =  0.00 ; given_sigma_ion_O_2D_cm2(1494) =    0.00 ; given_sigma_ion_O_2P_cm2(1494) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1494) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1494) =    0.00
  given_sigma_ion_O_4S_cm2(1495) =  0.00 ; given_sigma_ion_O_2D_cm2(1495) =    0.00 ; given_sigma_ion_O_2P_cm2(1495) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1495) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1495) =    0.00
  given_sigma_ion_O_4S_cm2(1496) =  0.00 ; given_sigma_ion_O_2D_cm2(1496) =    0.00 ; given_sigma_ion_O_2P_cm2(1496) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1496) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1496) =    0.00
  given_sigma_ion_O_4S_cm2(1497) =  0.00 ; given_sigma_ion_O_2D_cm2(1497) =    0.00 ; given_sigma_ion_O_2P_cm2(1497) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1497) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1497) =    0.00
  given_sigma_ion_O_4S_cm2(1498) =  0.00 ; given_sigma_ion_O_2D_cm2(1498) =    0.00 ; given_sigma_ion_O_2P_cm2(1498) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1498) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1498) =    0.00
  given_sigma_ion_O_4S_cm2(1499) =  0.00 ; given_sigma_ion_O_2D_cm2(1499) =    0.00 ; given_sigma_ion_O_2P_cm2(1499) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1499) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1499) =    0.00
  given_sigma_ion_O_4S_cm2(1500) =  0.00 ; given_sigma_ion_O_2D_cm2(1500) =    0.00 ; given_sigma_ion_O_2P_cm2(1500) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1500) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1500) =    0.00
  given_sigma_ion_O_4S_cm2(1501) =  0.00 ; given_sigma_ion_O_2D_cm2(1501) =    0.00 ; given_sigma_ion_O_2P_cm2(1501) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1501) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1501) =    0.00
  given_sigma_ion_O_4S_cm2(1502) =  0.00 ; given_sigma_ion_O_2D_cm2(1502) =    0.00 ; given_sigma_ion_O_2P_cm2(1502) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1502) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1502) =    0.00
  given_sigma_ion_O_4S_cm2(1503) =  0.00 ; given_sigma_ion_O_2D_cm2(1503) =    0.00 ; given_sigma_ion_O_2P_cm2(1503) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1503) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1503) =    0.00
  given_sigma_ion_O_4S_cm2(1504) =  0.00 ; given_sigma_ion_O_2D_cm2(1504) =    0.00 ; given_sigma_ion_O_2P_cm2(1504) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1504) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1504) =    0.00
  given_sigma_ion_O_4S_cm2(1505) =  0.00 ; given_sigma_ion_O_2D_cm2(1505) =    0.00 ; given_sigma_ion_O_2P_cm2(1505) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1505) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1505) =    0.00
  given_sigma_ion_O_4S_cm2(1506) =  0.00 ; given_sigma_ion_O_2D_cm2(1506) =    0.00 ; given_sigma_ion_O_2P_cm2(1506) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1506) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1506) =    0.00
  given_sigma_ion_O_4S_cm2(1507) =  0.00 ; given_sigma_ion_O_2D_cm2(1507) =    0.00 ; given_sigma_ion_O_2P_cm2(1507) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1507) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1507) =    0.00
  given_sigma_ion_O_4S_cm2(1508) =  0.00 ; given_sigma_ion_O_2D_cm2(1508) =    0.00 ; given_sigma_ion_O_2P_cm2(1508) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1508) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1508) =    0.00
  given_sigma_ion_O_4S_cm2(1509) =  0.00 ; given_sigma_ion_O_2D_cm2(1509) =    0.00 ; given_sigma_ion_O_2P_cm2(1509) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1509) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1509) =    0.00
  given_sigma_ion_O_4S_cm2(1510) =  0.00 ; given_sigma_ion_O_2D_cm2(1510) =    0.00 ; given_sigma_ion_O_2P_cm2(1510) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1510) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1510) =    0.00
  given_sigma_ion_O_4S_cm2(1511) =  0.00 ; given_sigma_ion_O_2D_cm2(1511) =    0.00 ; given_sigma_ion_O_2P_cm2(1511) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1511) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1511) =    0.00
  given_sigma_ion_O_4S_cm2(1512) =  0.00 ; given_sigma_ion_O_2D_cm2(1512) =    0.00 ; given_sigma_ion_O_2P_cm2(1512) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1512) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1512) =    0.00
  given_sigma_ion_O_4S_cm2(1513) =  0.00 ; given_sigma_ion_O_2D_cm2(1513) =    0.00 ; given_sigma_ion_O_2P_cm2(1513) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1513) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1513) =    0.00
  given_sigma_ion_O_4S_cm2(1514) =  0.00 ; given_sigma_ion_O_2D_cm2(1514) =    0.00 ; given_sigma_ion_O_2P_cm2(1514) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1514) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1514) =    0.00
  given_sigma_ion_O_4S_cm2(1515) =  0.00 ; given_sigma_ion_O_2D_cm2(1515) =    0.00 ; given_sigma_ion_O_2P_cm2(1515) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1515) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1515) =    0.00
  given_sigma_ion_O_4S_cm2(1516) =  0.00 ; given_sigma_ion_O_2D_cm2(1516) =    0.00 ; given_sigma_ion_O_2P_cm2(1516) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1516) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1516) =    0.00
  given_sigma_ion_O_4S_cm2(1517) =  0.00 ; given_sigma_ion_O_2D_cm2(1517) =    0.00 ; given_sigma_ion_O_2P_cm2(1517) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1517) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1517) =    0.00
  given_sigma_ion_O_4S_cm2(1518) =  0.00 ; given_sigma_ion_O_2D_cm2(1518) =    0.00 ; given_sigma_ion_O_2P_cm2(1518) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1518) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1518) =    0.00
  given_sigma_ion_O_4S_cm2(1519) =  0.00 ; given_sigma_ion_O_2D_cm2(1519) =    0.00 ; given_sigma_ion_O_2P_cm2(1519) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1519) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1519) =    0.00
  given_sigma_ion_O_4S_cm2(1520) =  0.00 ; given_sigma_ion_O_2D_cm2(1520) =    0.00 ; given_sigma_ion_O_2P_cm2(1520) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1520) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1520) =    0.00
  given_sigma_ion_O_4S_cm2(1521) =  0.00 ; given_sigma_ion_O_2D_cm2(1521) =    0.00 ; given_sigma_ion_O_2P_cm2(1521) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1521) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1521) =    0.00
  given_sigma_ion_O_4S_cm2(1522) =  0.00 ; given_sigma_ion_O_2D_cm2(1522) =    0.00 ; given_sigma_ion_O_2P_cm2(1522) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1522) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1522) =    0.00
  given_sigma_ion_O_4S_cm2(1523) =  0.00 ; given_sigma_ion_O_2D_cm2(1523) =    0.00 ; given_sigma_ion_O_2P_cm2(1523) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1523) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1523) =    0.00
  given_sigma_ion_O_4S_cm2(1524) =  0.00 ; given_sigma_ion_O_2D_cm2(1524) =    0.00 ; given_sigma_ion_O_2P_cm2(1524) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1524) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1524) =    0.00
  given_sigma_ion_O_4S_cm2(1525) =  0.00 ; given_sigma_ion_O_2D_cm2(1525) =    0.00 ; given_sigma_ion_O_2P_cm2(1525) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1525) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1525) =    0.00
  given_sigma_ion_O_4S_cm2(1526) =  0.00 ; given_sigma_ion_O_2D_cm2(1526) =    0.00 ; given_sigma_ion_O_2P_cm2(1526) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1526) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1526) =    0.00
  given_sigma_ion_O_4S_cm2(1527) =  0.00 ; given_sigma_ion_O_2D_cm2(1527) =    0.00 ; given_sigma_ion_O_2P_cm2(1527) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1527) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1527) =    0.00
  given_sigma_ion_O_4S_cm2(1528) =  0.00 ; given_sigma_ion_O_2D_cm2(1528) =    0.00 ; given_sigma_ion_O_2P_cm2(1528) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1528) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1528) =    0.00
  given_sigma_ion_O_4S_cm2(1529) =  0.00 ; given_sigma_ion_O_2D_cm2(1529) =    0.00 ; given_sigma_ion_O_2P_cm2(1529) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1529) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1529) =    0.00
  given_sigma_ion_O_4S_cm2(1530) =  0.00 ; given_sigma_ion_O_2D_cm2(1530) =    0.00 ; given_sigma_ion_O_2P_cm2(1530) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1530) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1530) =    0.00
  given_sigma_ion_O_4S_cm2(1531) =  0.00 ; given_sigma_ion_O_2D_cm2(1531) =    0.00 ; given_sigma_ion_O_2P_cm2(1531) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1531) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1531) =    0.00
  given_sigma_ion_O_4S_cm2(1532) =  0.00 ; given_sigma_ion_O_2D_cm2(1532) =    0.00 ; given_sigma_ion_O_2P_cm2(1532) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1532) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1532) =    0.00
  given_sigma_ion_O_4S_cm2(1533) =  0.00 ; given_sigma_ion_O_2D_cm2(1533) =    0.00 ; given_sigma_ion_O_2P_cm2(1533) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1533) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1533) =    0.00
  given_sigma_ion_O_4S_cm2(1534) =  0.00 ; given_sigma_ion_O_2D_cm2(1534) =    0.00 ; given_sigma_ion_O_2P_cm2(1534) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1534) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1534) =    0.00
  given_sigma_ion_O_4S_cm2(1535) =  0.00 ; given_sigma_ion_O_2D_cm2(1535) =    0.00 ; given_sigma_ion_O_2P_cm2(1535) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1535) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1535) =    0.00
  given_sigma_ion_O_4S_cm2(1536) =  0.00 ; given_sigma_ion_O_2D_cm2(1536) =    0.00 ; given_sigma_ion_O_2P_cm2(1536) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1536) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1536) =    0.00
  given_sigma_ion_O_4S_cm2(1537) =  0.00 ; given_sigma_ion_O_2D_cm2(1537) =    0.00 ; given_sigma_ion_O_2P_cm2(1537) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1537) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1537) =    0.00
  given_sigma_ion_O_4S_cm2(1538) =  0.00 ; given_sigma_ion_O_2D_cm2(1538) =    0.00 ; given_sigma_ion_O_2P_cm2(1538) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1538) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1538) =    0.00
  given_sigma_ion_O_4S_cm2(1539) =  0.00 ; given_sigma_ion_O_2D_cm2(1539) =    0.00 ; given_sigma_ion_O_2P_cm2(1539) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1539) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1539) =    0.00
  given_sigma_ion_O_4S_cm2(1540) =  0.00 ; given_sigma_ion_O_2D_cm2(1540) =    0.00 ; given_sigma_ion_O_2P_cm2(1540) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1540) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1540) =    0.00
  given_sigma_ion_O_4S_cm2(1541) =  0.00 ; given_sigma_ion_O_2D_cm2(1541) =    0.00 ; given_sigma_ion_O_2P_cm2(1541) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1541) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1541) =    0.00
  given_sigma_ion_O_4S_cm2(1542) =  0.00 ; given_sigma_ion_O_2D_cm2(1542) =    0.00 ; given_sigma_ion_O_2P_cm2(1542) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1542) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1542) =    0.00
  given_sigma_ion_O_4S_cm2(1543) =  0.00 ; given_sigma_ion_O_2D_cm2(1543) =    0.00 ; given_sigma_ion_O_2P_cm2(1543) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1543) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1543) =    0.00
  given_sigma_ion_O_4S_cm2(1544) =  0.00 ; given_sigma_ion_O_2D_cm2(1544) =    0.00 ; given_sigma_ion_O_2P_cm2(1544) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1544) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1544) =    0.00
  given_sigma_ion_O_4S_cm2(1545) =  0.00 ; given_sigma_ion_O_2D_cm2(1545) =    0.00 ; given_sigma_ion_O_2P_cm2(1545) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1545) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1545) =    0.00
  given_sigma_ion_O_4S_cm2(1546) =  0.00 ; given_sigma_ion_O_2D_cm2(1546) =    0.00 ; given_sigma_ion_O_2P_cm2(1546) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1546) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1546) =    0.00
  given_sigma_ion_O_4S_cm2(1547) =  0.00 ; given_sigma_ion_O_2D_cm2(1547) =    0.00 ; given_sigma_ion_O_2P_cm2(1547) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1547) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1547) =    0.00
  given_sigma_ion_O_4S_cm2(1548) =  0.00 ; given_sigma_ion_O_2D_cm2(1548) =    0.00 ; given_sigma_ion_O_2P_cm2(1548) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1548) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1548) =    0.00
  given_sigma_ion_O_4S_cm2(1549) =  0.00 ; given_sigma_ion_O_2D_cm2(1549) =    0.00 ; given_sigma_ion_O_2P_cm2(1549) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1549) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1549) =    0.00
  given_sigma_ion_O_4S_cm2(1550) =  0.00 ; given_sigma_ion_O_2D_cm2(1550) =    0.00 ; given_sigma_ion_O_2P_cm2(1550) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1550) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1550) =    0.00
  given_sigma_ion_O_4S_cm2(1551) =  0.00 ; given_sigma_ion_O_2D_cm2(1551) =    0.00 ; given_sigma_ion_O_2P_cm2(1551) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1551) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1551) =    0.00
  given_sigma_ion_O_4S_cm2(1552) =  0.00 ; given_sigma_ion_O_2D_cm2(1552) =    0.00 ; given_sigma_ion_O_2P_cm2(1552) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1552) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1552) =    0.00
  given_sigma_ion_O_4S_cm2(1553) =  0.00 ; given_sigma_ion_O_2D_cm2(1553) =    0.00 ; given_sigma_ion_O_2P_cm2(1553) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1553) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1553) =    0.00
  given_sigma_ion_O_4S_cm2(1554) =  0.00 ; given_sigma_ion_O_2D_cm2(1554) =    0.00 ; given_sigma_ion_O_2P_cm2(1554) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1554) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1554) =    0.00
  given_sigma_ion_O_4S_cm2(1555) =  0.00 ; given_sigma_ion_O_2D_cm2(1555) =    0.00 ; given_sigma_ion_O_2P_cm2(1555) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1555) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1555) =    0.00
  given_sigma_ion_O_4S_cm2(1556) =  0.00 ; given_sigma_ion_O_2D_cm2(1556) =    0.00 ; given_sigma_ion_O_2P_cm2(1556) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1556) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1556) =    0.00
  given_sigma_ion_O_4S_cm2(1557) =  0.00 ; given_sigma_ion_O_2D_cm2(1557) =    0.00 ; given_sigma_ion_O_2P_cm2(1557) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1557) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1557) =    0.00
  given_sigma_ion_O_4S_cm2(1558) =  0.00 ; given_sigma_ion_O_2D_cm2(1558) =    0.00 ; given_sigma_ion_O_2P_cm2(1558) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1558) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1558) =    0.00
  given_sigma_ion_O_4S_cm2(1559) =  0.00 ; given_sigma_ion_O_2D_cm2(1559) =    0.00 ; given_sigma_ion_O_2P_cm2(1559) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1559) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1559) =    0.00
  given_sigma_ion_O_4S_cm2(1560) =  0.00 ; given_sigma_ion_O_2D_cm2(1560) =    0.00 ; given_sigma_ion_O_2P_cm2(1560) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1560) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1560) =    0.00
  given_sigma_ion_O_4S_cm2(1561) =  0.00 ; given_sigma_ion_O_2D_cm2(1561) =    0.00 ; given_sigma_ion_O_2P_cm2(1561) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1561) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1561) =    0.00
  given_sigma_ion_O_4S_cm2(1562) =  0.00 ; given_sigma_ion_O_2D_cm2(1562) =    0.00 ; given_sigma_ion_O_2P_cm2(1562) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1562) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1562) =    0.00
  given_sigma_ion_O_4S_cm2(1563) =  0.00 ; given_sigma_ion_O_2D_cm2(1563) =    0.00 ; given_sigma_ion_O_2P_cm2(1563) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1563) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1563) =    0.00
  given_sigma_ion_O_4S_cm2(1564) =  0.00 ; given_sigma_ion_O_2D_cm2(1564) =    0.00 ; given_sigma_ion_O_2P_cm2(1564) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1564) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1564) =    0.00
  given_sigma_ion_O_4S_cm2(1565) =  0.00 ; given_sigma_ion_O_2D_cm2(1565) =    0.00 ; given_sigma_ion_O_2P_cm2(1565) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1565) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1565) =    0.00
  given_sigma_ion_O_4S_cm2(1566) =  0.00 ; given_sigma_ion_O_2D_cm2(1566) =    0.00 ; given_sigma_ion_O_2P_cm2(1566) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1566) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1566) =    0.00
  given_sigma_ion_O_4S_cm2(1567) =  0.00 ; given_sigma_ion_O_2D_cm2(1567) =    0.00 ; given_sigma_ion_O_2P_cm2(1567) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1567) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1567) =    0.00
  given_sigma_ion_O_4S_cm2(1568) =  0.00 ; given_sigma_ion_O_2D_cm2(1568) =    0.00 ; given_sigma_ion_O_2P_cm2(1568) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1568) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1568) =    0.00
  given_sigma_ion_O_4S_cm2(1569) =  0.00 ; given_sigma_ion_O_2D_cm2(1569) =    0.00 ; given_sigma_ion_O_2P_cm2(1569) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1569) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1569) =    0.00
  given_sigma_ion_O_4S_cm2(1570) =  0.00 ; given_sigma_ion_O_2D_cm2(1570) =    0.00 ; given_sigma_ion_O_2P_cm2(1570) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1570) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1570) =    0.00
  given_sigma_ion_O_4S_cm2(1571) =  0.00 ; given_sigma_ion_O_2D_cm2(1571) =    0.00 ; given_sigma_ion_O_2P_cm2(1571) =    0.00 ; given_sigma_ion_O_4Pst_cm2(1571) =    0.00 ; given_sigma_ion_O_2Pst_cm2(1571) =    0.00


 given_sigma_tot_O_cm2(   1) =    0.04
 given_sigma_tot_O_cm2(   2) =    0.06
 given_sigma_tot_O_cm2(   3) =    0.07
 given_sigma_tot_O_cm2(   4) =    0.07
 given_sigma_tot_O_cm2(   5) =    0.07
 given_sigma_tot_O_cm2(   6) =    0.07
 given_sigma_tot_O_cm2(   7) =    0.08
 given_sigma_tot_O_cm2(   8) =    0.08
 given_sigma_tot_O_cm2(   9) =    0.08
 given_sigma_tot_O_cm2(  10) =    0.10
 given_sigma_tot_O_cm2(  11) =    0.11
 given_sigma_tot_O_cm2(  12) =    0.14
 given_sigma_tot_O_cm2(  13) =    0.14
 given_sigma_tot_O_cm2(  14) =    0.16
 given_sigma_tot_O_cm2(  15) =    0.16
 given_sigma_tot_O_cm2(  16) =    0.16
 given_sigma_tot_O_cm2(  17) =    0.18
 given_sigma_tot_O_cm2(  18) =    0.18
 given_sigma_tot_O_cm2(  19) =    0.18
 given_sigma_tot_O_cm2(  20) =    0.18
 given_sigma_tot_O_cm2(  21) =    0.18
 given_sigma_tot_O_cm2(  22) =    0.18
 given_sigma_tot_O_cm2(  23) =    0.18
 given_sigma_tot_O_cm2(  24) =    0.18
 given_sigma_tot_O_cm2(  25) =    0.19
 given_sigma_tot_O_cm2(  26) =    0.20
 given_sigma_tot_O_cm2(  27) =    0.20
 given_sigma_tot_O_cm2(  28) =    0.21
 given_sigma_tot_O_cm2(  29) =    0.21
 given_sigma_tot_O_cm2(  30) =    0.22
 given_sigma_tot_O_cm2(  31) =    0.22
 given_sigma_tot_O_cm2(  32) =    0.23
 given_sigma_tot_O_cm2(  33) =    0.24
 given_sigma_tot_O_cm2(  34) =    0.24
 given_sigma_tot_O_cm2(  35) =    0.25
 given_sigma_tot_O_cm2(  36) =    0.25
 given_sigma_tot_O_cm2(  37) =    0.26
 given_sigma_tot_O_cm2(  38) =    0.26
 given_sigma_tot_O_cm2(  39) =    0.27
 given_sigma_tot_O_cm2(  40) =    0.28
 given_sigma_tot_O_cm2(  41) =    0.28
 given_sigma_tot_O_cm2(  42) =    0.29
 given_sigma_tot_O_cm2(  43) =    0.29
 given_sigma_tot_O_cm2(  44) =    0.30
 given_sigma_tot_O_cm2(  45) =    0.31
 given_sigma_tot_O_cm2(  46) =    0.32
 given_sigma_tot_O_cm2(  47) =    0.34
 given_sigma_tot_O_cm2(  48) =    0.35
 given_sigma_tot_O_cm2(  49) =    0.35
 given_sigma_tot_O_cm2(  50) =    0.36
 given_sigma_tot_O_cm2(  51) =    0.37
 given_sigma_tot_O_cm2(  52) =    0.38
 given_sigma_tot_O_cm2(  53) =    0.38
 given_sigma_tot_O_cm2(  54) =    0.39
 given_sigma_tot_O_cm2(  55) =    0.39
 given_sigma_tot_O_cm2(  56) =    0.39
 given_sigma_tot_O_cm2(  57) =    0.40
 given_sigma_tot_O_cm2(  58) =    0.41
 given_sigma_tot_O_cm2(  59) =    0.41
 given_sigma_tot_O_cm2(  60) =    0.42
 given_sigma_tot_O_cm2(  61) =    0.43
 given_sigma_tot_O_cm2(  62) =    0.44
 given_sigma_tot_O_cm2(  63) =    0.44
 given_sigma_tot_O_cm2(  64) =    0.45
 given_sigma_tot_O_cm2(  65) =    0.47
 given_sigma_tot_O_cm2(  66) =    0.48
 given_sigma_tot_O_cm2(  67) =    0.48
 given_sigma_tot_O_cm2(  68) =    0.50
 given_sigma_tot_O_cm2(  69) =    0.51
 given_sigma_tot_O_cm2(  70) =    0.52
 given_sigma_tot_O_cm2(  71) =    0.53
 given_sigma_tot_O_cm2(  72) =    0.54
 given_sigma_tot_O_cm2(  73) =    0.55
 given_sigma_tot_O_cm2(  74) =    0.55
 given_sigma_tot_O_cm2(  75) =    0.57
 given_sigma_tot_O_cm2(  76) =    0.58
 given_sigma_tot_O_cm2(  77) =    0.58
 given_sigma_tot_O_cm2(  78) =    0.59
 given_sigma_tot_O_cm2(  79) =    0.59
 given_sigma_tot_O_cm2(  80) =    0.61
 given_sigma_tot_O_cm2(  81) =    0.62
 given_sigma_tot_O_cm2(  82) =    0.63
 given_sigma_tot_O_cm2(  83) =    0.64
 given_sigma_tot_O_cm2(  84) =    0.64
 given_sigma_tot_O_cm2(  85) =    0.65
 given_sigma_tot_O_cm2(  86) =    0.65
 given_sigma_tot_O_cm2(  87) =    0.66
 given_sigma_tot_O_cm2(  88) =    0.66
 given_sigma_tot_O_cm2(  89) =    0.67
 given_sigma_tot_O_cm2(  90) =    0.68
 given_sigma_tot_O_cm2(  91) =    0.68
 given_sigma_tot_O_cm2(  92) =    0.69
 given_sigma_tot_O_cm2(  93) =    0.70
 given_sigma_tot_O_cm2(  94) =    0.70
 given_sigma_tot_O_cm2(  95) =    0.72
 given_sigma_tot_O_cm2(  96) =    0.73
 given_sigma_tot_O_cm2(  97) =    0.74
 given_sigma_tot_O_cm2(  98) =    0.75
 given_sigma_tot_O_cm2(  99) =    0.76
 given_sigma_tot_O_cm2( 100) =    0.76
 given_sigma_tot_O_cm2( 101) =    0.77
 given_sigma_tot_O_cm2( 102) =    0.78
 given_sigma_tot_O_cm2( 103) =    0.79
 given_sigma_tot_O_cm2( 104) =    0.79
 given_sigma_tot_O_cm2( 105) =    0.81
 given_sigma_tot_O_cm2( 106) =    0.82
 given_sigma_tot_O_cm2( 107) =    0.82
 given_sigma_tot_O_cm2( 108) =    0.83
 given_sigma_tot_O_cm2( 109) =    0.84
 given_sigma_tot_O_cm2( 110) =    0.85
 given_sigma_tot_O_cm2( 111) =    0.85
 given_sigma_tot_O_cm2( 112) =    0.86
 given_sigma_tot_O_cm2( 113) =    0.86
 given_sigma_tot_O_cm2( 114) =    0.87
 given_sigma_tot_O_cm2( 115) =    0.87
 given_sigma_tot_O_cm2( 116) =    0.88
 given_sigma_tot_O_cm2( 117) =    0.88
 given_sigma_tot_O_cm2( 118) =    0.89
 given_sigma_tot_O_cm2( 119) =    0.90
 given_sigma_tot_O_cm2( 120) =    0.90
 given_sigma_tot_O_cm2( 121) =    0.90
 given_sigma_tot_O_cm2( 122) =    0.91
 given_sigma_tot_O_cm2( 123) =    0.92
 given_sigma_tot_O_cm2( 124) =    0.92
 given_sigma_tot_O_cm2( 125) =    0.93
 given_sigma_tot_O_cm2( 126) =    0.93
 given_sigma_tot_O_cm2( 127) =    0.94
 given_sigma_tot_O_cm2( 128) =    0.95
 given_sigma_tot_O_cm2( 129) =    0.96
 given_sigma_tot_O_cm2( 130) =    0.96
 given_sigma_tot_O_cm2( 131) =    0.97
 given_sigma_tot_O_cm2( 132) =    0.98
 given_sigma_tot_O_cm2( 133) =    0.99
 given_sigma_tot_O_cm2( 134) =    1.00
 given_sigma_tot_O_cm2( 135) =    1.01
 given_sigma_tot_O_cm2( 136) =    1.02
 given_sigma_tot_O_cm2( 137) =    1.03
 given_sigma_tot_O_cm2( 138) =    1.04
 given_sigma_tot_O_cm2( 139) =    1.04
 given_sigma_tot_O_cm2( 140) =    1.05
 given_sigma_tot_O_cm2( 141) =    1.06
 given_sigma_tot_O_cm2( 142) =    1.07
 given_sigma_tot_O_cm2( 143) =    1.07
 given_sigma_tot_O_cm2( 144) =    1.10
 given_sigma_tot_O_cm2( 145) =    1.11
 given_sigma_tot_O_cm2( 146) =    1.11
 given_sigma_tot_O_cm2( 147) =    1.12
 given_sigma_tot_O_cm2( 148) =    1.13
 given_sigma_tot_O_cm2( 149) =    1.14
 given_sigma_tot_O_cm2( 150) =    1.15
 given_sigma_tot_O_cm2( 151) =    1.15
 given_sigma_tot_O_cm2( 152) =    1.16
 given_sigma_tot_O_cm2( 153) =    1.17
 given_sigma_tot_O_cm2( 154) =    1.18
 given_sigma_tot_O_cm2( 155) =    1.18
 given_sigma_tot_O_cm2( 156) =    1.19
 given_sigma_tot_O_cm2( 157) =    1.20
 given_sigma_tot_O_cm2( 158) =    1.21
 given_sigma_tot_O_cm2( 159) =    1.21
 given_sigma_tot_O_cm2( 160) =    1.22
 given_sigma_tot_O_cm2( 161) =    1.23
 given_sigma_tot_O_cm2( 162) =    1.24
 given_sigma_tot_O_cm2( 163) =    1.25
 given_sigma_tot_O_cm2( 164) =    1.25
 given_sigma_tot_O_cm2( 165) =    1.27
 given_sigma_tot_O_cm2( 166) =    1.28
 given_sigma_tot_O_cm2( 167) =    1.29
 given_sigma_tot_O_cm2( 168) =    1.31
 given_sigma_tot_O_cm2( 169) =    1.33
 given_sigma_tot_O_cm2( 170) =    1.33
 given_sigma_tot_O_cm2( 171) =    1.33
 given_sigma_tot_O_cm2( 172) =    1.34
 given_sigma_tot_O_cm2( 173) =    1.35
 given_sigma_tot_O_cm2( 174) =    1.36
 given_sigma_tot_O_cm2( 175) =    1.37
 given_sigma_tot_O_cm2( 176) =    1.38
 given_sigma_tot_O_cm2( 177) =    1.41
 given_sigma_tot_O_cm2( 178) =    1.42
 given_sigma_tot_O_cm2( 179) =    1.43
 given_sigma_tot_O_cm2( 180) =    1.45
 given_sigma_tot_O_cm2( 181) =    1.46
 given_sigma_tot_O_cm2( 182) =    1.47
 given_sigma_tot_O_cm2( 183) =    1.49
 given_sigma_tot_O_cm2( 184) =    1.51
 given_sigma_tot_O_cm2( 185) =    1.52
 given_sigma_tot_O_cm2( 186) =    1.53
 given_sigma_tot_O_cm2( 187) =    1.54
 given_sigma_tot_O_cm2( 188) =    1.58
 given_sigma_tot_O_cm2( 189) =    1.60
 given_sigma_tot_O_cm2( 190) =    1.61
 given_sigma_tot_O_cm2( 191) =    1.61
 given_sigma_tot_O_cm2( 192) =    1.64
 given_sigma_tot_O_cm2( 193) =    1.65
 given_sigma_tot_O_cm2( 194) =    1.68
 given_sigma_tot_O_cm2( 195) =    1.69
 given_sigma_tot_O_cm2( 196) =    1.71
 given_sigma_tot_O_cm2( 197) =    1.76
 given_sigma_tot_O_cm2( 198) =    1.77
 given_sigma_tot_O_cm2( 199) =    1.79
 given_sigma_tot_O_cm2( 200) =    1.80
 given_sigma_tot_O_cm2( 201) =    1.82
 given_sigma_tot_O_cm2( 202) =    1.84
 given_sigma_tot_O_cm2( 203) =    1.86
 given_sigma_tot_O_cm2( 204) =    1.89
 given_sigma_tot_O_cm2( 205) =    1.97
 given_sigma_tot_O_cm2( 206) =    2.04
 given_sigma_tot_O_cm2( 207) =    2.04
 given_sigma_tot_O_cm2( 208) =    2.05
 given_sigma_tot_O_cm2( 209) =    2.06
 given_sigma_tot_O_cm2( 210) =    2.07
 given_sigma_tot_O_cm2( 211) =    2.08
 given_sigma_tot_O_cm2( 212) =    2.20
 given_sigma_tot_O_cm2( 213) =    2.23
 given_sigma_tot_O_cm2( 214) =    2.24
 given_sigma_tot_O_cm2( 215) =    2.28
 given_sigma_tot_O_cm2( 216) =    2.27
 given_sigma_tot_O_cm2( 217) =    2.26
 given_sigma_tot_O_cm2( 218) =    2.24
 given_sigma_tot_O_cm2( 219) =    2.24
 given_sigma_tot_O_cm2( 220) =    2.23
 given_sigma_tot_O_cm2( 221) =    2.22
 given_sigma_tot_O_cm2( 222) =    2.21
 given_sigma_tot_O_cm2( 223) =    2.21
 given_sigma_tot_O_cm2( 224) =    2.19
 given_sigma_tot_O_cm2( 225) =    2.19
 given_sigma_tot_O_cm2( 226) =    2.18
 given_sigma_tot_O_cm2( 227) =    2.18
 given_sigma_tot_O_cm2( 228) =    2.36
 given_sigma_tot_O_cm2( 229) =    2.41
 given_sigma_tot_O_cm2( 230) =    2.51
 given_sigma_tot_O_cm2( 231) =    2.54
 given_sigma_tot_O_cm2( 232) =    2.80
 given_sigma_tot_O_cm2( 233) =    2.85
 given_sigma_tot_O_cm2( 234) =    2.87
 given_sigma_tot_O_cm2( 235) =    3.02
 given_sigma_tot_O_cm2( 236) =    3.07
 given_sigma_tot_O_cm2( 237) =    3.09
 given_sigma_tot_O_cm2( 238) =    3.12
 given_sigma_tot_O_cm2( 239) =    3.17
 given_sigma_tot_O_cm2( 240) =    3.19
 given_sigma_tot_O_cm2( 241) =    3.26
 given_sigma_tot_O_cm2( 242) =    3.33
 given_sigma_tot_O_cm2( 243) =    3.38
 given_sigma_tot_O_cm2( 244) =    3.39
 given_sigma_tot_O_cm2( 245) =    3.49
 given_sigma_tot_O_cm2( 246) =    3.52
 given_sigma_tot_O_cm2( 247) =    3.53
 given_sigma_tot_O_cm2( 248) =    3.55
 given_sigma_tot_O_cm2( 249) =    3.66
 given_sigma_tot_O_cm2( 250) =    3.66
 given_sigma_tot_O_cm2( 251) =    3.69
 given_sigma_tot_O_cm2( 252) =    3.73
 given_sigma_tot_O_cm2( 253) =    3.74
 given_sigma_tot_O_cm2( 254) =    3.75
 given_sigma_tot_O_cm2( 255) =    3.77
 given_sigma_tot_O_cm2( 256) =    3.78
 given_sigma_tot_O_cm2( 257) =    3.79
 given_sigma_tot_O_cm2( 258) =    3.83
 given_sigma_tot_O_cm2( 259) =    3.83
 given_sigma_tot_O_cm2( 260) =    3.87
 given_sigma_tot_O_cm2( 261) =    3.88
 given_sigma_tot_O_cm2( 262) =    3.89
 given_sigma_tot_O_cm2( 263) =    3.90
 given_sigma_tot_O_cm2( 264) =    3.91
 given_sigma_tot_O_cm2( 265) =    3.92
 given_sigma_tot_O_cm2( 266) =    3.93
 given_sigma_tot_O_cm2( 267) =    3.97
 given_sigma_tot_O_cm2( 268) =    3.98
 given_sigma_tot_O_cm2( 269) =    4.02
 given_sigma_tot_O_cm2( 270) =    4.03
 given_sigma_tot_O_cm2( 271) =    4.04
 given_sigma_tot_O_cm2( 272) =    4.09
 given_sigma_tot_O_cm2( 273) =    4.11
 given_sigma_tot_O_cm2( 274) =    4.12
 given_sigma_tot_O_cm2( 275) =    4.13
 given_sigma_tot_O_cm2( 276) =    4.16
 given_sigma_tot_O_cm2( 277) =    4.18
 given_sigma_tot_O_cm2( 278) =    4.20
 given_sigma_tot_O_cm2( 279) =    4.25
 given_sigma_tot_O_cm2( 280) =    4.25
 given_sigma_tot_O_cm2( 281) =    4.30
 given_sigma_tot_O_cm2( 282) =    4.30
 given_sigma_tot_O_cm2( 283) =    4.33
 given_sigma_tot_O_cm2( 284) =    4.37
 given_sigma_tot_O_cm2( 285) =    4.41
 given_sigma_tot_O_cm2( 286) =    4.45
 given_sigma_tot_O_cm2( 287) =    4.48
 given_sigma_tot_O_cm2( 288) =    4.50
 given_sigma_tot_O_cm2( 289) =    4.54
 given_sigma_tot_O_cm2( 290) =    4.55
 given_sigma_tot_O_cm2( 291) =    4.57
 given_sigma_tot_O_cm2( 292) =    4.62
 given_sigma_tot_O_cm2( 293) =    4.62
 given_sigma_tot_O_cm2( 294) =    4.63
 given_sigma_tot_O_cm2( 295) =    4.65
 given_sigma_tot_O_cm2( 296) =    4.68
 given_sigma_tot_O_cm2( 297) =    4.72
 given_sigma_tot_O_cm2( 298) =    4.72
 given_sigma_tot_O_cm2( 299) =    4.72
 given_sigma_tot_O_cm2( 300) =    4.76
 given_sigma_tot_O_cm2( 301) =    4.79
 given_sigma_tot_O_cm2( 302) =    4.84
 given_sigma_tot_O_cm2( 303) =    4.86
 given_sigma_tot_O_cm2( 304) =    4.87
 given_sigma_tot_O_cm2( 305) =    4.88
 given_sigma_tot_O_cm2( 306) =    4.93
 given_sigma_tot_O_cm2( 307) =    4.93
 given_sigma_tot_O_cm2( 308) =    4.96
 given_sigma_tot_O_cm2( 309) =    4.99
 given_sigma_tot_O_cm2( 310) =    5.01
 given_sigma_tot_O_cm2( 311) =    5.05
 given_sigma_tot_O_cm2( 312) =    5.05
 given_sigma_tot_O_cm2( 313) =    5.07
 given_sigma_tot_O_cm2( 314) =    5.11
 given_sigma_tot_O_cm2( 315) =    5.12
 given_sigma_tot_O_cm2( 316) =    5.15
 given_sigma_tot_O_cm2( 317) =    5.16
 given_sigma_tot_O_cm2( 318) =    5.16
 given_sigma_tot_O_cm2( 319) =    5.21
 given_sigma_tot_O_cm2( 320) =    5.22
 given_sigma_tot_O_cm2( 321) =    5.23
 given_sigma_tot_O_cm2( 322) =    5.26
 given_sigma_tot_O_cm2( 323) =    5.29
 given_sigma_tot_O_cm2( 324) =    5.30
 given_sigma_tot_O_cm2( 325) =    5.32
 given_sigma_tot_O_cm2( 326) =    5.34
 given_sigma_tot_O_cm2( 327) =    5.37
 given_sigma_tot_O_cm2( 328) =    5.40
 given_sigma_tot_O_cm2( 329) =    5.42
 given_sigma_tot_O_cm2( 330) =    5.44
 given_sigma_tot_O_cm2( 331) =    5.45
 given_sigma_tot_O_cm2( 332) =    5.50
 given_sigma_tot_O_cm2( 333) =    5.53
 given_sigma_tot_O_cm2( 334) =    5.57
 given_sigma_tot_O_cm2( 335) =    5.57
 given_sigma_tot_O_cm2( 336) =    5.59
 given_sigma_tot_O_cm2( 337) =    5.62
 given_sigma_tot_O_cm2( 338) =    5.66
 given_sigma_tot_O_cm2( 339) =    5.68
 given_sigma_tot_O_cm2( 340) =    5.71
 given_sigma_tot_O_cm2( 341) =    5.73
 given_sigma_tot_O_cm2( 342) =    5.74
 given_sigma_tot_O_cm2( 343) =    5.76
 given_sigma_tot_O_cm2( 344) =    5.77
 given_sigma_tot_O_cm2( 345) =    5.79
 given_sigma_tot_O_cm2( 346) =    5.81
 given_sigma_tot_O_cm2( 347) =    5.84
 given_sigma_tot_O_cm2( 348) =    5.88
 given_sigma_tot_O_cm2( 349) =    5.91
 given_sigma_tot_O_cm2( 350) =    5.92
 given_sigma_tot_O_cm2( 351) =    5.98
 given_sigma_tot_O_cm2( 352) =    6.00
 given_sigma_tot_O_cm2( 353) =    6.02
 given_sigma_tot_O_cm2( 354) =    6.05
 given_sigma_tot_O_cm2( 355) =    6.06
 given_sigma_tot_O_cm2( 356) =    6.07
 given_sigma_tot_O_cm2( 357) =    6.08
 given_sigma_tot_O_cm2( 358) =    6.10
 given_sigma_tot_O_cm2( 359) =    6.14
 given_sigma_tot_O_cm2( 360) =    6.15
 given_sigma_tot_O_cm2( 361) =    6.19
 given_sigma_tot_O_cm2( 362) =    6.25
 given_sigma_tot_O_cm2( 363) =    6.30
 given_sigma_tot_O_cm2( 364) =    6.32
 given_sigma_tot_O_cm2( 365) =    6.33
 given_sigma_tot_O_cm2( 366) =    6.48
 given_sigma_tot_O_cm2( 367) =    6.50
 given_sigma_tot_O_cm2( 368) =    6.52
 given_sigma_tot_O_cm2( 369) =    6.57
 given_sigma_tot_O_cm2( 370) =    6.59
 given_sigma_tot_O_cm2( 371) =    6.65
 given_sigma_tot_O_cm2( 372) =    6.68
 given_sigma_tot_O_cm2( 373) =    6.69
 given_sigma_tot_O_cm2( 374) =    6.70
 given_sigma_tot_O_cm2( 375) =    6.72
 given_sigma_tot_O_cm2( 376) =    6.74
 given_sigma_tot_O_cm2( 377) =    6.74
 given_sigma_tot_O_cm2( 378) =    6.75
 given_sigma_tot_O_cm2( 379) =    6.79
 given_sigma_tot_O_cm2( 380) =    6.85
 given_sigma_tot_O_cm2( 381) =    6.90
 given_sigma_tot_O_cm2( 382) =    6.93
 given_sigma_tot_O_cm2( 383) =    7.08
 given_sigma_tot_O_cm2( 384) =    7.13
 given_sigma_tot_O_cm2( 385) =    7.16
 given_sigma_tot_O_cm2( 386) =    7.17
 given_sigma_tot_O_cm2( 387) =    7.31
 given_sigma_tot_O_cm2( 388) =    7.35
 given_sigma_tot_O_cm2( 389) =    7.40
 given_sigma_tot_O_cm2( 390) =    7.41
 given_sigma_tot_O_cm2( 391) =    7.42
 given_sigma_tot_O_cm2( 392) =    7.43
 given_sigma_tot_O_cm2( 393) =    7.44
 given_sigma_tot_O_cm2( 394) =    7.48
 given_sigma_tot_O_cm2( 395) =    7.48
 given_sigma_tot_O_cm2( 396) =    7.49
 given_sigma_tot_O_cm2( 397) =    7.54
 given_sigma_tot_O_cm2( 398) =    7.55
 given_sigma_tot_O_cm2( 399) =    7.58
 given_sigma_tot_O_cm2( 400) =    7.68
 given_sigma_tot_O_cm2( 401) =    7.70
 given_sigma_tot_O_cm2( 402) =    7.90
 given_sigma_tot_O_cm2( 403) =    8.07
 given_sigma_tot_O_cm2( 404) =    8.11
 given_sigma_tot_O_cm2( 405) =    8.12
 given_sigma_tot_O_cm2( 406) =    8.22
 given_sigma_tot_O_cm2( 407) =    8.24
 given_sigma_tot_O_cm2( 408) =    8.25
 given_sigma_tot_O_cm2( 409) =    8.27
 given_sigma_tot_O_cm2( 410) =    8.43
 given_sigma_tot_O_cm2( 411) =    8.60
 given_sigma_tot_O_cm2( 412) =    8.76
 given_sigma_tot_O_cm2( 413) =    8.77
 given_sigma_tot_O_cm2( 414) =    8.92
 given_sigma_tot_O_cm2( 415) =    9.09
 given_sigma_tot_O_cm2( 416) =    9.09
 given_sigma_tot_O_cm2( 417) =    9.11
 given_sigma_tot_O_cm2( 418) =    9.16
 given_sigma_tot_O_cm2( 419) =    9.25
 given_sigma_tot_O_cm2( 420) =    9.25
 given_sigma_tot_O_cm2( 421) =    9.37
 given_sigma_tot_O_cm2( 422) =    9.44
 given_sigma_tot_O_cm2( 423) =    9.57
 given_sigma_tot_O_cm2( 424) =    9.60
 given_sigma_tot_O_cm2( 425) =    9.72
 given_sigma_tot_O_cm2( 426) =    9.73
 given_sigma_tot_O_cm2( 427) =    9.84
 given_sigma_tot_O_cm2( 428) =    9.90
 given_sigma_tot_O_cm2( 429) =    9.95
 given_sigma_tot_O_cm2( 430) =   10.20
 given_sigma_tot_O_cm2( 431) =   10.50
 given_sigma_tot_O_cm2( 432) =   11.78
 given_sigma_tot_O_cm2( 433) =   11.80
 given_sigma_tot_O_cm2( 434) =   11.71
 given_sigma_tot_O_cm2( 435) =   11.66
 given_sigma_tot_O_cm2( 436) =   11.64
 given_sigma_tot_O_cm2( 437) =   11.54
 given_sigma_tot_O_cm2( 438) =   11.40
 given_sigma_tot_O_cm2( 439) =   11.32
 given_sigma_tot_O_cm2( 440) =   11.24
 given_sigma_tot_O_cm2( 441) =   11.16
 given_sigma_tot_O_cm2( 442) =   11.08
 given_sigma_tot_O_cm2( 443) =   11.00
 given_sigma_tot_O_cm2( 444) =   11.03
 given_sigma_tot_O_cm2( 445) =   11.06
 given_sigma_tot_O_cm2( 446) =   11.09
 given_sigma_tot_O_cm2( 447) =   11.11
 given_sigma_tot_O_cm2( 448) =   11.14
 given_sigma_tot_O_cm2( 449) =   11.17
 given_sigma_tot_O_cm2( 450) =   11.20
 given_sigma_tot_O_cm2( 451) =   11.21
 given_sigma_tot_O_cm2( 452) =   11.23
 given_sigma_tot_O_cm2( 453) =   11.24
 given_sigma_tot_O_cm2( 454) =   11.27
 given_sigma_tot_O_cm2( 455) =   11.30
 given_sigma_tot_O_cm2( 456) =   11.32
 given_sigma_tot_O_cm2( 457) =   11.34
 given_sigma_tot_O_cm2( 458) =   11.36
 given_sigma_tot_O_cm2( 459) =   11.38
 given_sigma_tot_O_cm2( 460) =   11.40
 given_sigma_tot_O_cm2( 461) =   11.42
 given_sigma_tot_O_cm2( 462) =   11.44
 given_sigma_tot_O_cm2( 463) =   11.46
 given_sigma_tot_O_cm2( 464) =   11.48
 given_sigma_tot_O_cm2( 465) =   11.50
 given_sigma_tot_O_cm2( 466) =   11.49
 given_sigma_tot_O_cm2( 467) =   11.48
 given_sigma_tot_O_cm2( 468) =   11.46
 given_sigma_tot_O_cm2( 469) =   11.45
 given_sigma_tot_O_cm2( 470) =   11.43
 given_sigma_tot_O_cm2( 471) =   11.41
 given_sigma_tot_O_cm2( 472) =   11.39
 given_sigma_tot_O_cm2( 473) =   11.39
 given_sigma_tot_O_cm2( 474) =   11.38
 given_sigma_tot_O_cm2( 475) =   11.43
 given_sigma_tot_O_cm2( 476) =   11.60
 given_sigma_tot_O_cm2( 477) =   11.76
 given_sigma_tot_O_cm2( 478) =   11.92
 given_sigma_tot_O_cm2( 479) =   12.00
 given_sigma_tot_O_cm2( 480) =   11.69
 given_sigma_tot_O_cm2( 481) =   11.08
 given_sigma_tot_O_cm2( 482) =   10.90
 given_sigma_tot_O_cm2( 483) =   12.90
 given_sigma_tot_O_cm2( 484) =   12.47
 given_sigma_tot_O_cm2( 485) =   11.60
 given_sigma_tot_O_cm2( 486) =   11.47
 given_sigma_tot_O_cm2( 487) =   11.13
 given_sigma_tot_O_cm2( 488) =   11.10
 given_sigma_tot_O_cm2( 489) =   10.61
 given_sigma_tot_O_cm2( 490) =   10.50
 given_sigma_tot_O_cm2( 491) =   14.00
 given_sigma_tot_O_cm2( 492) =   13.37
 given_sigma_tot_O_cm2( 493) =   12.10
 given_sigma_tot_O_cm2( 494) =   11.70
 given_sigma_tot_O_cm2( 495) =   11.10
 given_sigma_tot_O_cm2( 496) =   11.15
 given_sigma_tot_O_cm2( 497) =   11.27
 given_sigma_tot_O_cm2( 498) =   11.30
 given_sigma_tot_O_cm2( 499) =   11.20
 given_sigma_tot_O_cm2( 500) =   11.08
 given_sigma_tot_O_cm2( 501) =   11.00
 given_sigma_tot_O_cm2( 502) =   10.76
 given_sigma_tot_O_cm2( 503) =   10.18
 given_sigma_tot_O_cm2( 504) =   10.00
 given_sigma_tot_O_cm2( 505) =   11.00
 given_sigma_tot_O_cm2( 506) =   12.76
 given_sigma_tot_O_cm2( 507) =   15.10
 given_sigma_tot_O_cm2( 508) =   13.85
 given_sigma_tot_O_cm2( 509) =   12.80
 given_sigma_tot_O_cm2( 510) =   12.61
 given_sigma_tot_O_cm2( 511) =   12.24
 given_sigma_tot_O_cm2( 512) =   11.87
 given_sigma_tot_O_cm2( 513) =   11.80
 given_sigma_tot_O_cm2( 514) =   11.82
 given_sigma_tot_O_cm2( 515) =   11.84
 given_sigma_tot_O_cm2( 516) =   11.86
 given_sigma_tot_O_cm2( 517) =   11.88
 given_sigma_tot_O_cm2( 518) =   11.90
 given_sigma_tot_O_cm2( 519) =   11.92
 given_sigma_tot_O_cm2( 520) =   11.93
 given_sigma_tot_O_cm2( 521) =   11.95
 given_sigma_tot_O_cm2( 522) =   11.97
 given_sigma_tot_O_cm2( 523) =   11.99
 given_sigma_tot_O_cm2( 524) =   12.01
 given_sigma_tot_O_cm2( 525) =   12.03
 given_sigma_tot_O_cm2( 526) =   12.03
 given_sigma_tot_O_cm2( 527) =   12.05
 given_sigma_tot_O_cm2( 528) =   12.07
 given_sigma_tot_O_cm2( 529) =   12.10
 given_sigma_tot_O_cm2( 530) =   12.10
 given_sigma_tot_O_cm2( 531) =   12.00
 given_sigma_tot_O_cm2( 532) =   11.97
 given_sigma_tot_O_cm2( 533) =   11.80
 given_sigma_tot_O_cm2( 534) =   11.50
 given_sigma_tot_O_cm2( 535) =   11.46
 given_sigma_tot_O_cm2( 536) =   11.30
 given_sigma_tot_O_cm2( 537) =   11.61
 given_sigma_tot_O_cm2( 538) =   12.14
 given_sigma_tot_O_cm2( 539) =   12.40
 given_sigma_tot_O_cm2( 540) =   12.33
 given_sigma_tot_O_cm2( 541) =   12.20
 given_sigma_tot_O_cm2( 542) =   16.10
 given_sigma_tot_O_cm2( 543) =   15.57
 given_sigma_tot_O_cm2( 544) =   13.70
 given_sigma_tot_O_cm2( 545) =   13.50
 given_sigma_tot_O_cm2( 546) =   13.43
 given_sigma_tot_O_cm2( 547) =   13.10
 given_sigma_tot_O_cm2( 548) =   13.04
 given_sigma_tot_O_cm2( 549) =   13.00
 given_sigma_tot_O_cm2( 550) =   12.95
 given_sigma_tot_O_cm2( 551) =   12.90
 given_sigma_tot_O_cm2( 552) =   12.86
 given_sigma_tot_O_cm2( 553) =   12.73
 given_sigma_tot_O_cm2( 554) =   12.65
 given_sigma_tot_O_cm2( 555) =   12.57
 given_sigma_tot_O_cm2( 556) =   12.38
 given_sigma_tot_O_cm2( 557) =   12.19
 given_sigma_tot_O_cm2( 558) =   12.10
 given_sigma_tot_O_cm2( 559) =   12.00
 given_sigma_tot_O_cm2( 560) =   12.16
 given_sigma_tot_O_cm2( 561) =   12.27
 given_sigma_tot_O_cm2( 562) =   12.26
 given_sigma_tot_O_cm2( 563) =   12.21
 given_sigma_tot_O_cm2( 564) =   12.17
 given_sigma_tot_O_cm2( 565) =   12.12
 given_sigma_tot_O_cm2( 566) =   12.08
 given_sigma_tot_O_cm2( 567) =   12.03
 given_sigma_tot_O_cm2( 568) =   11.99
 given_sigma_tot_O_cm2( 569) =   11.94
 given_sigma_tot_O_cm2( 570) =   11.93
 given_sigma_tot_O_cm2( 571) =   11.93
 given_sigma_tot_O_cm2( 572) =   11.90
 given_sigma_tot_O_cm2( 573) =   11.77
 given_sigma_tot_O_cm2( 574) =   11.76
 given_sigma_tot_O_cm2( 575) =   11.77
 given_sigma_tot_O_cm2( 576) =   11.79
 given_sigma_tot_O_cm2( 577) =   11.81
 given_sigma_tot_O_cm2( 578) =   11.82
 given_sigma_tot_O_cm2( 579) =   11.87
 given_sigma_tot_O_cm2( 580) =   11.90
 given_sigma_tot_O_cm2( 581) =   11.96
 given_sigma_tot_O_cm2( 582) =   12.00
 given_sigma_tot_O_cm2( 583) =   12.01
 given_sigma_tot_O_cm2( 584) =   12.01
 given_sigma_tot_O_cm2( 585) =   12.06
 given_sigma_tot_O_cm2( 586) =   12.10
 given_sigma_tot_O_cm2( 587) =   12.24
 given_sigma_tot_O_cm2( 588) =   12.30
 given_sigma_tot_O_cm2( 589) =   12.36
 given_sigma_tot_O_cm2( 590) =   12.39
 given_sigma_tot_O_cm2( 591) =   12.48
 given_sigma_tot_O_cm2( 592) =   12.50
 given_sigma_tot_O_cm2( 593) =   12.53
 given_sigma_tot_O_cm2( 594) =   12.59
 given_sigma_tot_O_cm2( 595) =   12.59
 given_sigma_tot_O_cm2( 596) =   12.61
 given_sigma_tot_O_cm2( 597) =   12.67
 given_sigma_tot_O_cm2( 598) =   12.70
 given_sigma_tot_O_cm2( 599) =   12.76
 given_sigma_tot_O_cm2( 600) =   12.80
 given_sigma_tot_O_cm2( 601) =   12.87
 given_sigma_tot_O_cm2( 602) =   12.90
 given_sigma_tot_O_cm2( 603) =   12.92
 given_sigma_tot_O_cm2( 604) =   12.95
 given_sigma_tot_O_cm2( 605) =   13.00
 given_sigma_tot_O_cm2( 606) =   13.01
 given_sigma_tot_O_cm2( 607) =   13.08
 given_sigma_tot_O_cm2( 608) =   13.09
 given_sigma_tot_O_cm2( 609) =   13.10
 given_sigma_tot_O_cm2( 610) =   13.12
 given_sigma_tot_O_cm2( 611) =   13.18
 given_sigma_tot_O_cm2( 612) =   13.20
 given_sigma_tot_O_cm2( 613) =   13.22
 given_sigma_tot_O_cm2( 614) =   13.24
 given_sigma_tot_O_cm2( 615) =   13.25
 given_sigma_tot_O_cm2( 616) =   13.27
 given_sigma_tot_O_cm2( 617) =   13.30
 given_sigma_tot_O_cm2( 618) =   13.30
 given_sigma_tot_O_cm2( 619) =   13.34
 given_sigma_tot_O_cm2( 620) =   13.38
 given_sigma_tot_O_cm2( 621) =   13.38
 given_sigma_tot_O_cm2( 622) =   13.39
 given_sigma_tot_O_cm2( 623) =   13.40
 given_sigma_tot_O_cm2( 624) =   13.40
 given_sigma_tot_O_cm2( 625) =   13.40
 given_sigma_tot_O_cm2( 626) =   13.40
 given_sigma_tot_O_cm2( 627) =   13.40
 given_sigma_tot_O_cm2( 628) =   13.40
 given_sigma_tot_O_cm2( 629) =   13.40
 given_sigma_tot_O_cm2( 630) =   13.40
 given_sigma_tot_O_cm2( 631) =   13.40
 given_sigma_tot_O_cm2( 632) =   13.40
 given_sigma_tot_O_cm2( 633) =   13.40
 given_sigma_tot_O_cm2( 634) =   13.40
 given_sigma_tot_O_cm2( 635) =   13.40
 given_sigma_tot_O_cm2( 636) =   13.40
 given_sigma_tot_O_cm2( 637) =   13.40
 given_sigma_tot_O_cm2( 638) =   13.40
 given_sigma_tot_O_cm2( 639) =   13.40
 given_sigma_tot_O_cm2( 640) =   13.40
 given_sigma_tot_O_cm2( 641) =   13.40
 given_sigma_tot_O_cm2( 642) =   13.40
 given_sigma_tot_O_cm2( 643) =   13.40
 given_sigma_tot_O_cm2( 644) =   13.40
 given_sigma_tot_O_cm2( 645) =   13.40
 given_sigma_tot_O_cm2( 646) =   13.40
 given_sigma_tot_O_cm2( 647) =   13.40
 given_sigma_tot_O_cm2( 648) =   13.40
 given_sigma_tot_O_cm2( 649) =   13.40
 given_sigma_tot_O_cm2( 650) =   13.40
 given_sigma_tot_O_cm2( 651) =   13.40
 given_sigma_tot_O_cm2( 652) =   13.40
 given_sigma_tot_O_cm2( 653) =   13.40
 given_sigma_tot_O_cm2( 654) =   13.40
 given_sigma_tot_O_cm2( 655) =   13.40
 given_sigma_tot_O_cm2( 656) =   13.40
 given_sigma_tot_O_cm2( 657) =   13.40
 given_sigma_tot_O_cm2( 658) =   13.40
 given_sigma_tot_O_cm2( 659) =   13.40
 given_sigma_tot_O_cm2( 660) =   13.40
 given_sigma_tot_O_cm2( 661) =   13.39
 given_sigma_tot_O_cm2( 662) =   13.38
 given_sigma_tot_O_cm2( 663) =   13.38
 given_sigma_tot_O_cm2( 664) =   13.37
 given_sigma_tot_O_cm2( 665) =   13.36
 given_sigma_tot_O_cm2( 666) =   13.36
 given_sigma_tot_O_cm2( 667) =   13.35
 given_sigma_tot_O_cm2( 668) =   13.35
 given_sigma_tot_O_cm2( 669) =   13.34
 given_sigma_tot_O_cm2( 670) =   13.34
 given_sigma_tot_O_cm2( 671) =   13.34
 given_sigma_tot_O_cm2( 672) =   13.33
 given_sigma_tot_O_cm2( 673) =   13.33
 given_sigma_tot_O_cm2( 674) =   13.32
 given_sigma_tot_O_cm2( 675) =   13.32
 given_sigma_tot_O_cm2( 676) =   13.31
 given_sigma_tot_O_cm2( 677) =   13.31
 given_sigma_tot_O_cm2( 678) =   13.30
 given_sigma_tot_O_cm2( 679) =   13.29
 given_sigma_tot_O_cm2( 680) =   13.27
 given_sigma_tot_O_cm2( 681) =   13.25
 given_sigma_tot_O_cm2( 682) =   13.24
 given_sigma_tot_O_cm2( 683) =   13.23
 given_sigma_tot_O_cm2( 684) =   13.21
 given_sigma_tot_O_cm2( 685) =   13.18
 given_sigma_tot_O_cm2( 686) =   13.17
 given_sigma_tot_O_cm2( 687) =   13.17
 given_sigma_tot_O_cm2( 688) =   13.15
 given_sigma_tot_O_cm2( 689) =   13.12
 given_sigma_tot_O_cm2( 690) =   13.10
 given_sigma_tot_O_cm2( 691) =   13.10
 given_sigma_tot_O_cm2( 692) =   13.09
 given_sigma_tot_O_cm2( 693) =   13.07
 given_sigma_tot_O_cm2( 694) =   13.06
 given_sigma_tot_O_cm2( 695) =   13.03
 given_sigma_tot_O_cm2( 696) =   13.02
 given_sigma_tot_O_cm2( 697) =   13.00
 given_sigma_tot_O_cm2( 698) =   12.99
 given_sigma_tot_O_cm2( 699) =   12.96
 given_sigma_tot_O_cm2( 700) =   12.93
 given_sigma_tot_O_cm2( 701) =   12.92
 given_sigma_tot_O_cm2( 702) =   12.88
 given_sigma_tot_O_cm2( 703) =   12.84
 given_sigma_tot_O_cm2( 704) =   12.80
 given_sigma_tot_O_cm2( 705) =   12.76
 given_sigma_tot_O_cm2( 706) =   12.72
 given_sigma_tot_O_cm2( 707) =   12.71
 given_sigma_tot_O_cm2( 708) =   12.68
 given_sigma_tot_O_cm2( 709) =   12.64
 given_sigma_tot_O_cm2( 710) =   12.60
 given_sigma_tot_O_cm2( 711) =   12.57
 given_sigma_tot_O_cm2( 712) =   12.49
 given_sigma_tot_O_cm2( 713) =   12.44
 given_sigma_tot_O_cm2( 714) =   12.38
 given_sigma_tot_O_cm2( 715) =   12.26
 given_sigma_tot_O_cm2( 716) =   12.15
 given_sigma_tot_O_cm2( 717) =   12.08
 given_sigma_tot_O_cm2( 718) =   12.05
 given_sigma_tot_O_cm2( 719) =   12.00
 given_sigma_tot_O_cm2( 720) =   11.92
 given_sigma_tot_O_cm2( 721) =   11.89
 given_sigma_tot_O_cm2( 722) =   11.72
 given_sigma_tot_O_cm2( 723) =   11.67
 given_sigma_tot_O_cm2( 724) =   11.56
 given_sigma_tot_O_cm2( 725) =   11.39
 given_sigma_tot_O_cm2( 726) =   11.30
 given_sigma_tot_O_cm2( 727) =   11.23
 given_sigma_tot_O_cm2( 728) =   11.16
 given_sigma_tot_O_cm2( 729) =   11.07
 given_sigma_tot_O_cm2( 730) =   10.98
 given_sigma_tot_O_cm2( 731) =   10.92
 given_sigma_tot_O_cm2( 732) =   10.75
 given_sigma_tot_O_cm2( 733) =   10.64
 given_sigma_tot_O_cm2( 734) =   10.61
 given_sigma_tot_O_cm2( 735) =   10.57
 given_sigma_tot_O_cm2( 736) =   10.51
 given_sigma_tot_O_cm2( 737) =   10.41
 given_sigma_tot_O_cm2( 738) =   10.38
 given_sigma_tot_O_cm2( 739) =   10.30
 given_sigma_tot_O_cm2( 740) =   10.25
 given_sigma_tot_O_cm2( 741) =   10.21
 given_sigma_tot_O_cm2( 742) =   10.15
 given_sigma_tot_O_cm2( 743) =   10.08
 given_sigma_tot_O_cm2( 744) =   10.00
 given_sigma_tot_O_cm2( 745) =   10.30
 given_sigma_tot_O_cm2( 746) =   10.60
 given_sigma_tot_O_cm2( 747) =   12.00
 given_sigma_tot_O_cm2( 748) =   17.67
 given_sigma_tot_O_cm2( 749) =   15.40
 given_sigma_tot_O_cm2( 750) =   14.52
 given_sigma_tot_O_cm2( 751) =   14.15
 given_sigma_tot_O_cm2( 752) =   13.90
 given_sigma_tot_O_cm2( 753) =   14.29
 given_sigma_tot_O_cm2( 754) =   14.68
 given_sigma_tot_O_cm2( 755) =   15.07
 given_sigma_tot_O_cm2( 756) =   15.20
 given_sigma_tot_O_cm2( 757) =   14.41
 given_sigma_tot_O_cm2( 758) =   14.02
 given_sigma_tot_O_cm2( 759) =   12.84
 given_sigma_tot_O_cm2( 760) =   11.66
 given_sigma_tot_O_cm2( 761) =    9.69
 given_sigma_tot_O_cm2( 762) =    9.65
 given_sigma_tot_O_cm2( 763) =    9.52
 given_sigma_tot_O_cm2( 764) =    9.29
 given_sigma_tot_O_cm2( 765) =    9.25
 given_sigma_tot_O_cm2( 766) =    9.70
 given_sigma_tot_O_cm2( 767) =   10.20
 given_sigma_tot_O_cm2( 768) =   11.20
 given_sigma_tot_O_cm2( 769) =   13.13
 given_sigma_tot_O_cm2( 770) =   24.68
 given_sigma_tot_O_cm2( 771) =   25.49
 given_sigma_tot_O_cm2( 772) =   22.56
 given_sigma_tot_O_cm2( 773) =   19.54
 given_sigma_tot_O_cm2( 774) =   16.51
 given_sigma_tot_O_cm2( 775) =   15.50
 given_sigma_tot_O_cm2( 776) =   13.35
 given_sigma_tot_O_cm2( 777) =   11.20
 given_sigma_tot_O_cm2( 778) =   12.48
 given_sigma_tot_O_cm2( 779) =   12.41
 given_sigma_tot_O_cm2( 780) =   11.92
 given_sigma_tot_O_cm2( 781) =   11.43
 given_sigma_tot_O_cm2( 782) =    9.80
 given_sigma_tot_O_cm2( 783) =    9.47
 given_sigma_tot_O_cm2( 784) =    9.48
 given_sigma_tot_O_cm2( 785) =    9.48
 given_sigma_tot_O_cm2( 786) =    9.49
 given_sigma_tot_O_cm2( 787) =    9.49
 given_sigma_tot_O_cm2( 788) =    9.51
 given_sigma_tot_O_cm2( 789) =    9.51
 given_sigma_tot_O_cm2( 790) =    9.57
 given_sigma_tot_O_cm2( 791) =    9.51
 given_sigma_tot_O_cm2( 792) =    9.37
 given_sigma_tot_O_cm2( 793) =    9.00
 given_sigma_tot_O_cm2( 794) =    9.01
 given_sigma_tot_O_cm2( 795) =    9.02
 given_sigma_tot_O_cm2( 796) =    9.04
 given_sigma_tot_O_cm2( 797) =    9.06
 given_sigma_tot_O_cm2( 798) =    8.91
 given_sigma_tot_O_cm2( 799) =    8.80
 given_sigma_tot_O_cm2( 800) =    9.36
 given_sigma_tot_O_cm2( 801) =    9.70
 given_sigma_tot_O_cm2( 802) =   40.00
 given_sigma_tot_O_cm2( 803) =   37.09
 given_sigma_tot_O_cm2( 804) =   32.71
 given_sigma_tot_O_cm2( 805) =   28.34
 given_sigma_tot_O_cm2( 806) =   19.60
 given_sigma_tot_O_cm2( 807) =   16.17
 given_sigma_tot_O_cm2( 808) =   14.79
 given_sigma_tot_O_cm2( 809) =   14.11
 given_sigma_tot_O_cm2( 810) =   12.05
 given_sigma_tot_O_cm2( 811) =    9.99
 given_sigma_tot_O_cm2( 812) =    9.30
 given_sigma_tot_O_cm2( 813) =    8.45
 given_sigma_tot_O_cm2( 814) =    8.03
 given_sigma_tot_O_cm2( 815) =    7.60
 given_sigma_tot_O_cm2( 816) =   11.80
 given_sigma_tot_O_cm2( 817) =   13.90
 given_sigma_tot_O_cm2( 818) =   13.59
 given_sigma_tot_O_cm2( 819) =   12.02
 given_sigma_tot_O_cm2( 820) =   11.46
 given_sigma_tot_O_cm2( 821) =   11.08
 given_sigma_tot_O_cm2( 822) =   10.46
 given_sigma_tot_O_cm2( 823) =    9.99
 given_sigma_tot_O_cm2( 824) =    9.60
 given_sigma_tot_O_cm2( 825) =    9.64
 given_sigma_tot_O_cm2( 826) =    9.71
 given_sigma_tot_O_cm2( 827) =    9.78
 given_sigma_tot_O_cm2( 828) =    9.82
 given_sigma_tot_O_cm2( 829) =    9.88
 given_sigma_tot_O_cm2( 830) =    9.85
 given_sigma_tot_O_cm2( 831) =    9.74
 given_sigma_tot_O_cm2( 832) =    9.71
 given_sigma_tot_O_cm2( 833) =    9.36
 given_sigma_tot_O_cm2( 834) =    9.15
 given_sigma_tot_O_cm2( 835) =    9.30
 given_sigma_tot_O_cm2( 836) =    9.30
 given_sigma_tot_O_cm2( 837) =    9.25
 given_sigma_tot_O_cm2( 838) =    9.26
 given_sigma_tot_O_cm2( 839) =    9.25
 given_sigma_tot_O_cm2( 840) =    9.24
 given_sigma_tot_O_cm2( 841) =    9.21
 given_sigma_tot_O_cm2( 842) =    9.19
 given_sigma_tot_O_cm2( 843) =    9.20
 given_sigma_tot_O_cm2( 844) =    9.21
 given_sigma_tot_O_cm2( 845) =    9.24
 given_sigma_tot_O_cm2( 846) =    9.24
 given_sigma_tot_O_cm2( 847) =    9.17
 given_sigma_tot_O_cm2( 848) =    9.05
 given_sigma_tot_O_cm2( 849) =    8.93
 given_sigma_tot_O_cm2( 850) =    8.88
 given_sigma_tot_O_cm2( 851) =    8.79
 given_sigma_tot_O_cm2( 852) =    8.84
 given_sigma_tot_O_cm2( 853) =    8.91
 given_sigma_tot_O_cm2( 854) =    8.97
 given_sigma_tot_O_cm2( 855) =    9.00
 given_sigma_tot_O_cm2( 856) =    9.00
 given_sigma_tot_O_cm2( 857) =    9.00
 given_sigma_tot_O_cm2( 858) =    8.96
 given_sigma_tot_O_cm2( 859) =    8.91
 given_sigma_tot_O_cm2( 860) =    8.94
 given_sigma_tot_O_cm2( 861) =    8.94
 given_sigma_tot_O_cm2( 862) =    8.97
 given_sigma_tot_O_cm2( 863) =    8.99
 given_sigma_tot_O_cm2( 864) =    9.02
 given_sigma_tot_O_cm2( 865) =    8.90
 given_sigma_tot_O_cm2( 866) =    8.75
 given_sigma_tot_O_cm2( 867) =    8.60
 given_sigma_tot_O_cm2( 868) =   26.10
 given_sigma_tot_O_cm2( 869) =   61.50
 given_sigma_tot_O_cm2( 870) =   54.30
 given_sigma_tot_O_cm2( 871) =   33.90
 given_sigma_tot_O_cm2( 872) =   27.10
 given_sigma_tot_O_cm2( 873) =   24.07
 given_sigma_tot_O_cm2( 874) =   20.02
 given_sigma_tot_O_cm2( 875) =   13.96
 given_sigma_tot_O_cm2( 876) =   10.92
 given_sigma_tot_O_cm2( 877) =    8.90
 given_sigma_tot_O_cm2( 878) =    8.99
 given_sigma_tot_O_cm2( 879) =    9.04
 given_sigma_tot_O_cm2( 880) =    9.12
 given_sigma_tot_O_cm2( 881) =    9.17
 given_sigma_tot_O_cm2( 882) =    8.98
 given_sigma_tot_O_cm2( 883) =    8.78
 given_sigma_tot_O_cm2( 884) =    8.69
 given_sigma_tot_O_cm2( 885) =    8.40
 given_sigma_tot_O_cm2( 886) =    8.33
 given_sigma_tot_O_cm2( 887) =    8.19
 given_sigma_tot_O_cm2( 888) =    8.12
 given_sigma_tot_O_cm2( 889) =    8.08
 given_sigma_tot_O_cm2( 890) =    8.03
 given_sigma_tot_O_cm2( 891) =    6.97
 given_sigma_tot_O_cm2( 892) =    5.21
 given_sigma_tot_O_cm2( 893) =    4.15
 given_sigma_tot_O_cm2( 894) =    8.62
 given_sigma_tot_O_cm2( 895) =   11.81
 given_sigma_tot_O_cm2( 896) =   15.00
 given_sigma_tot_O_cm2( 897) =   15.34
 given_sigma_tot_O_cm2( 898) =   16.01
 given_sigma_tot_O_cm2( 899) =   16.69
 given_sigma_tot_O_cm2( 900) =   10.40
 given_sigma_tot_O_cm2( 901) =   12.30
 given_sigma_tot_O_cm2( 902) =   15.15
 given_sigma_tot_O_cm2( 903) =   12.50
 given_sigma_tot_O_cm2( 904) =    8.90
 given_sigma_tot_O_cm2( 905) =   14.37
 given_sigma_tot_O_cm2( 906) =   16.20
 given_sigma_tot_O_cm2( 907) =    9.08
 given_sigma_tot_O_cm2( 908) =    7.30
 given_sigma_tot_O_cm2( 909) =   18.41
 given_sigma_tot_O_cm2( 910) =   18.91
 given_sigma_tot_O_cm2( 911) =   11.77
 given_sigma_tot_O_cm2( 912) =    9.98
 given_sigma_tot_O_cm2( 913) =    8.20
 given_sigma_tot_O_cm2( 914) =    6.33
 given_sigma_tot_O_cm2( 915) =    5.71
 given_sigma_tot_O_cm2( 916) =   10.08
 given_sigma_tot_O_cm2( 917) =   14.75
 given_sigma_tot_O_cm2( 918) =   21.76
 given_sigma_tot_O_cm2( 919) =   17.10
 given_sigma_tot_O_cm2( 920) =   10.10
 given_sigma_tot_O_cm2( 921) =    5.42
 given_sigma_tot_O_cm2( 922) =    4.85
 given_sigma_tot_O_cm2( 923) =    4.70
 given_sigma_tot_O_cm2( 924) =    5.49
 given_sigma_tot_O_cm2( 925) =    6.80
 given_sigma_tot_O_cm2( 926) =   20.10
 given_sigma_tot_O_cm2( 927) =   33.40
 given_sigma_tot_O_cm2( 928) =   20.33
 given_sigma_tot_O_cm2( 929) =   15.10
 given_sigma_tot_O_cm2( 930) =    9.82
 given_sigma_tot_O_cm2( 931) =    5.20
 given_sigma_tot_O_cm2( 932) =    5.00
 given_sigma_tot_O_cm2( 933) =    4.61
 given_sigma_tot_O_cm2( 934) =    4.35
 given_sigma_tot_O_cm2( 935) =    4.15
 given_sigma_tot_O_cm2( 936) =    4.36
 given_sigma_tot_O_cm2( 937) =    4.76
 given_sigma_tot_O_cm2( 938) =    4.90
 given_sigma_tot_O_cm2( 939) =   22.34
 given_sigma_tot_O_cm2( 940) =   28.16
 given_sigma_tot_O_cm2( 941) =   36.88
 given_sigma_tot_O_cm2( 942) =   45.60
 given_sigma_tot_O_cm2( 943) =   40.19
 given_sigma_tot_O_cm2( 944) =   34.77
 given_sigma_tot_O_cm2( 945) =   26.65
 given_sigma_tot_O_cm2( 946) =   13.11
 given_sigma_tot_O_cm2( 947) =    7.70
 given_sigma_tot_O_cm2( 948) =    9.57
 given_sigma_tot_O_cm2( 949) =   10.83
 given_sigma_tot_O_cm2( 950) =   12.70
 given_sigma_tot_O_cm2( 951) =   11.38
 given_sigma_tot_O_cm2( 952) =   10.94
 given_sigma_tot_O_cm2( 953) =    9.71
 given_sigma_tot_O_cm2( 954) =    8.30
 given_sigma_tot_O_cm2( 955) =    6.54
 given_sigma_tot_O_cm2( 956) =    3.90
 given_sigma_tot_O_cm2( 957) =    4.05
 given_sigma_tot_O_cm2( 958) =    4.06
 given_sigma_tot_O_cm2( 959) =    4.05
 given_sigma_tot_O_cm2( 960) =    4.05
 given_sigma_tot_O_cm2( 961) =    4.04
 given_sigma_tot_O_cm2( 962) =    4.04
 given_sigma_tot_O_cm2( 963) =    4.04
 given_sigma_tot_O_cm2( 964) =    4.03
 given_sigma_tot_O_cm2( 965) =    4.03
 given_sigma_tot_O_cm2( 966) =    4.02
 given_sigma_tot_O_cm2( 967) =    4.02
 given_sigma_tot_O_cm2( 968) =    4.02
 given_sigma_tot_O_cm2( 969) =    4.01
 given_sigma_tot_O_cm2( 970) =    4.01
 given_sigma_tot_O_cm2( 971) =    4.00
 given_sigma_tot_O_cm2( 972) =    4.00
 given_sigma_tot_O_cm2( 973) =    3.99
 given_sigma_tot_O_cm2( 974) =    3.99
 given_sigma_tot_O_cm2( 975) =    3.99
 given_sigma_tot_O_cm2( 976) =    3.98
 given_sigma_tot_O_cm2( 977) =    3.97
 given_sigma_tot_O_cm2( 978) =    3.97
 given_sigma_tot_O_cm2( 979) =    3.96
 given_sigma_tot_O_cm2( 980) =    3.95
 given_sigma_tot_O_cm2( 981) =    3.95
 given_sigma_tot_O_cm2( 982) =    3.95
 given_sigma_tot_O_cm2( 983) =    3.94
 given_sigma_tot_O_cm2( 984) =    3.93
 given_sigma_tot_O_cm2( 985) =    3.93
 given_sigma_tot_O_cm2( 986) =    3.93
 given_sigma_tot_O_cm2( 987) =    3.92
 given_sigma_tot_O_cm2( 988) =    3.92
 given_sigma_tot_O_cm2( 989) =    3.91
 given_sigma_tot_O_cm2( 990) =    3.90
 given_sigma_tot_O_cm2( 991) =    3.90
 given_sigma_tot_O_cm2( 992) =    3.89
 given_sigma_tot_O_cm2( 993) =    3.88
 given_sigma_tot_O_cm2( 994) =    3.88
 given_sigma_tot_O_cm2( 995) =    3.87
 given_sigma_tot_O_cm2( 996) =    3.86
 given_sigma_tot_O_cm2( 997) =    3.85
 given_sigma_tot_O_cm2( 998) =    3.84
 given_sigma_tot_O_cm2( 999) =    3.83
 given_sigma_tot_O_cm2(1000) =    3.82
 given_sigma_tot_O_cm2(1001) =    3.81
 given_sigma_tot_O_cm2(1002) =    3.80
 given_sigma_tot_O_cm2(1003) =    3.79
 given_sigma_tot_O_cm2(1004) =    3.79
 given_sigma_tot_O_cm2(1005) =    3.78
 given_sigma_tot_O_cm2(1006) =    3.77
 given_sigma_tot_O_cm2(1007) =    3.76
 given_sigma_tot_O_cm2(1008) =    3.75
 given_sigma_tot_O_cm2(1009) =    3.74
 given_sigma_tot_O_cm2(1010) =    3.73
 given_sigma_tot_O_cm2(1011) =    3.73
 given_sigma_tot_O_cm2(1012) =    3.72
 given_sigma_tot_O_cm2(1013) =    3.71
 given_sigma_tot_O_cm2(1014) =    3.70
 given_sigma_tot_O_cm2(1015) =    3.70
 given_sigma_tot_O_cm2(1016) =    3.70
 given_sigma_tot_O_cm2(1017) =    3.70
 given_sigma_tot_O_cm2(1018) =    3.70
 given_sigma_tot_O_cm2(1019) =    3.70
 given_sigma_tot_O_cm2(1020) =    3.70
 given_sigma_tot_O_cm2(1021) =    3.70
 given_sigma_tot_O_cm2(1022) =    3.70
 given_sigma_tot_O_cm2(1023) =    3.70
 given_sigma_tot_O_cm2(1024) =    3.70
 given_sigma_tot_O_cm2(1025) =    3.70
 given_sigma_tot_O_cm2(1026) =    3.70
 given_sigma_tot_O_cm2(1027) =    3.70
 given_sigma_tot_O_cm2(1028) =    3.70
 given_sigma_tot_O_cm2(1029) =    3.70
 given_sigma_tot_O_cm2(1030) =    3.70
 given_sigma_tot_O_cm2(1031) =    3.71
 given_sigma_tot_O_cm2(1032) =    3.71
 given_sigma_tot_O_cm2(1033) =    3.72
 given_sigma_tot_O_cm2(1034) =    3.72
 given_sigma_tot_O_cm2(1035) =    3.72
 given_sigma_tot_O_cm2(1036) =    3.73
 given_sigma_tot_O_cm2(1037) =    3.73
 given_sigma_tot_O_cm2(1038) =    3.74
 given_sigma_tot_O_cm2(1039) =    3.74
 given_sigma_tot_O_cm2(1040) =    3.75
 given_sigma_tot_O_cm2(1041) =    3.76
 given_sigma_tot_O_cm2(1042) =    3.76
 given_sigma_tot_O_cm2(1043) =    3.77
 given_sigma_tot_O_cm2(1044) =    3.77
 given_sigma_tot_O_cm2(1045) =    3.77
 given_sigma_tot_O_cm2(1046) =    3.78
 given_sigma_tot_O_cm2(1047) =    3.78
 given_sigma_tot_O_cm2(1048) =    4.23
 given_sigma_tot_O_cm2(1049) =    3.78
 given_sigma_tot_O_cm2(1050) =    4.86
 given_sigma_tot_O_cm2(1051) =   19.97
 given_sigma_tot_O_cm2(1052) =    2.70
 given_sigma_tot_O_cm2(1053) =    2.74
 given_sigma_tot_O_cm2(1054) =    2.80
 given_sigma_tot_O_cm2(1055) =    2.86
 given_sigma_tot_O_cm2(1056) =    2.92
 given_sigma_tot_O_cm2(1057) =    2.94
 given_sigma_tot_O_cm2(1058) =    3.03
 given_sigma_tot_O_cm2(1059) =    3.15
 given_sigma_tot_O_cm2(1060) =    3.23
 given_sigma_tot_O_cm2(1061) =    3.26
 given_sigma_tot_O_cm2(1062) =    3.38
 given_sigma_tot_O_cm2(1063) =    3.44
 given_sigma_tot_O_cm2(1064) =    3.50
 given_sigma_tot_O_cm2(1065) =    3.55
 given_sigma_tot_O_cm2(1066) =    3.61
 given_sigma_tot_O_cm2(1067) =    3.67
 given_sigma_tot_O_cm2(1068) =    3.73
 given_sigma_tot_O_cm2(1069) =    3.80
 given_sigma_tot_O_cm2(1070) =    3.84
 given_sigma_tot_O_cm2(1071) =    3.90
 given_sigma_tot_O_cm2(1072) =    3.96
 given_sigma_tot_O_cm2(1073) =   36.41
 given_sigma_tot_O_cm2(1074) =    9.14
 given_sigma_tot_O_cm2(1075) =    5.34
 given_sigma_tot_O_cm2(1076) =   48.94
 given_sigma_tot_O_cm2(1077) =    7.20
 given_sigma_tot_O_cm2(1078) =   11.73
 given_sigma_tot_O_cm2(1079) =    3.64
 given_sigma_tot_O_cm2(1080) =   15.37
 given_sigma_tot_O_cm2(1081) =    3.64
 given_sigma_tot_O_cm2(1082) =    5.91
 given_sigma_tot_O_cm2(1083) =    1.62
 given_sigma_tot_O_cm2(1084) =    2.06
 given_sigma_tot_O_cm2(1085) =    2.13
 given_sigma_tot_O_cm2(1086) =    2.14
 given_sigma_tot_O_cm2(1087) =    2.17
 given_sigma_tot_O_cm2(1088) =    2.19
 given_sigma_tot_O_cm2(1089) =    2.20
 given_sigma_tot_O_cm2(1090) =    2.22
 given_sigma_tot_O_cm2(1091) =    2.24
 given_sigma_tot_O_cm2(1092) =    2.27
 given_sigma_tot_O_cm2(1093) =    2.31
 given_sigma_tot_O_cm2(1094) =   15.23
 given_sigma_tot_O_cm2(1095) =    2.82
 given_sigma_tot_O_cm2(1096) =   32.62
 given_sigma_tot_O_cm2(1097) =    3.39
 given_sigma_tot_O_cm2(1098) =    7.90
 given_sigma_tot_O_cm2(1099) =    4.74
 given_sigma_tot_O_cm2(1100) =    2.80
 given_sigma_tot_O_cm2(1101) =    3.67
 given_sigma_tot_O_cm2(1102) =    3.95
 given_sigma_tot_O_cm2(1103) =    1.97
 given_sigma_tot_O_cm2(1104) =    2.21
 given_sigma_tot_O_cm2(1105) =    2.38
 given_sigma_tot_O_cm2(1106) =    2.79
 given_sigma_tot_O_cm2(1107) =    2.95
 given_sigma_tot_O_cm2(1108) =    3.19
 given_sigma_tot_O_cm2(1109) =    3.60
 given_sigma_tot_O_cm2(1110) =    3.93
 given_sigma_tot_O_cm2(1111) =    4.42
 given_sigma_tot_O_cm2(1112) =    4.83
 given_sigma_tot_O_cm2(1113) =    5.17
 given_sigma_tot_O_cm2(1114) =   65.14
 given_sigma_tot_O_cm2(1115) =    4.14
 given_sigma_tot_O_cm2(1116) =    1.86
 given_sigma_tot_O_cm2(1117) =    3.10
 given_sigma_tot_O_cm2(1118) =   22.02
 given_sigma_tot_O_cm2(1119) =    2.59
 given_sigma_tot_O_cm2(1120) =    7.40
 given_sigma_tot_O_cm2(1121) =    1.96
 given_sigma_tot_O_cm2(1122) =    1.55
 given_sigma_tot_O_cm2(1123) =    1.68
 given_sigma_tot_O_cm2(1124) =    1.82
 given_sigma_tot_O_cm2(1125) =    1.90
 given_sigma_tot_O_cm2(1126) =    1.95
 given_sigma_tot_O_cm2(1127) =    2.09
 given_sigma_tot_O_cm2(1128) =    2.22
 given_sigma_tot_O_cm2(1129) =    2.30
 given_sigma_tot_O_cm2(1130) =    2.36
 given_sigma_tot_O_cm2(1131) =    2.49
 given_sigma_tot_O_cm2(1132) =    2.76
 given_sigma_tot_O_cm2(1133) =    2.82
 given_sigma_tot_O_cm2(1134) =    3.03
 given_sigma_tot_O_cm2(1135) =    3.17
 given_sigma_tot_O_cm2(1136) =    3.27
 given_sigma_tot_O_cm2(1137) =    3.30
 given_sigma_tot_O_cm2(1138) =    3.30
 given_sigma_tot_O_cm2(1139) =    3.30
 given_sigma_tot_O_cm2(1140) =    3.30
 given_sigma_tot_O_cm2(1141) =    3.30
 given_sigma_tot_O_cm2(1142) =    3.30
 given_sigma_tot_O_cm2(1143) =    3.30
 given_sigma_tot_O_cm2(1144) =    3.30
 given_sigma_tot_O_cm2(1145) =    3.30
 given_sigma_tot_O_cm2(1146) =    3.30
 given_sigma_tot_O_cm2(1147) =    3.30
 given_sigma_tot_O_cm2(1148) =    3.30
 given_sigma_tot_O_cm2(1149) =    3.30
 given_sigma_tot_O_cm2(1150) =    3.30
 given_sigma_tot_O_cm2(1151) =    3.30
 given_sigma_tot_O_cm2(1152) =    3.30
 given_sigma_tot_O_cm2(1153) =    3.30
 given_sigma_tot_O_cm2(1154) =    3.30
 given_sigma_tot_O_cm2(1155) =    3.30
 given_sigma_tot_O_cm2(1156) =    3.30
 given_sigma_tot_O_cm2(1157) =    3.30
 given_sigma_tot_O_cm2(1158) =    3.30
 given_sigma_tot_O_cm2(1159) =    3.30
 given_sigma_tot_O_cm2(1160) =    3.30
 given_sigma_tot_O_cm2(1161) =    3.30
 given_sigma_tot_O_cm2(1162) =    3.30
 given_sigma_tot_O_cm2(1163) =    3.30
 given_sigma_tot_O_cm2(1164) =    3.30
 given_sigma_tot_O_cm2(1165) =    3.30
 given_sigma_tot_O_cm2(1166) =    3.30
 given_sigma_tot_O_cm2(1167) =    3.30
 given_sigma_tot_O_cm2(1168) =    3.30
 given_sigma_tot_O_cm2(1169) =    3.30
 given_sigma_tot_O_cm2(1170) =    3.30
 given_sigma_tot_O_cm2(1171) =    3.30
 given_sigma_tot_O_cm2(1172) =    3.30
 given_sigma_tot_O_cm2(1173) =    3.30
 given_sigma_tot_O_cm2(1174) =    3.30
 given_sigma_tot_O_cm2(1175) =    3.29
 given_sigma_tot_O_cm2(1176) =    3.29
 given_sigma_tot_O_cm2(1177) =    3.29
 given_sigma_tot_O_cm2(1178) =    3.28
 given_sigma_tot_O_cm2(1179) =    3.28
 given_sigma_tot_O_cm2(1180) =    3.27
 given_sigma_tot_O_cm2(1181) =    3.27
 given_sigma_tot_O_cm2(1182) =    3.26
 given_sigma_tot_O_cm2(1183) =    3.26
 given_sigma_tot_O_cm2(1184) =    3.26
 given_sigma_tot_O_cm2(1185) =    3.25
 given_sigma_tot_O_cm2(1186) =    3.25
 given_sigma_tot_O_cm2(1187) =    3.24
 given_sigma_tot_O_cm2(1188) =    3.23
 given_sigma_tot_O_cm2(1189) =    3.22
 given_sigma_tot_O_cm2(1190) =    3.22
 given_sigma_tot_O_cm2(1191) =    3.22
 given_sigma_tot_O_cm2(1192) =    3.21
 given_sigma_tot_O_cm2(1193) =    3.21
 given_sigma_tot_O_cm2(1194) =    3.21
 given_sigma_tot_O_cm2(1195) =    3.20
 given_sigma_tot_O_cm2(1196) =    3.19
 given_sigma_tot_O_cm2(1197) =    3.18
 given_sigma_tot_O_cm2(1198) =    3.17
 given_sigma_tot_O_cm2(1199) =    3.16
 given_sigma_tot_O_cm2(1200) =    3.16
 given_sigma_tot_O_cm2(1201) =    3.16
 given_sigma_tot_O_cm2(1202) =    3.15
 given_sigma_tot_O_cm2(1203) =    3.14
 given_sigma_tot_O_cm2(1204) =    3.14
 given_sigma_tot_O_cm2(1205) =    3.13
 given_sigma_tot_O_cm2(1206) =    3.12
 given_sigma_tot_O_cm2(1207) =    3.11
 given_sigma_tot_O_cm2(1208) =    3.10
 given_sigma_tot_O_cm2(1209) =    3.10
 given_sigma_tot_O_cm2(1210) =    3.10
 given_sigma_tot_O_cm2(1211) =    3.10
 given_sigma_tot_O_cm2(1212) =    3.10
 given_sigma_tot_O_cm2(1213) =    3.10
 given_sigma_tot_O_cm2(1214) =    3.10
 given_sigma_tot_O_cm2(1215) =    3.10
 given_sigma_tot_O_cm2(1216) =    3.10
 given_sigma_tot_O_cm2(1217) =    3.10
 given_sigma_tot_O_cm2(1218) =    3.10
 given_sigma_tot_O_cm2(1219) =    3.10
 given_sigma_tot_O_cm2(1220) =    3.10
 given_sigma_tot_O_cm2(1221) =    3.10
 given_sigma_tot_O_cm2(1222) =    3.10
 given_sigma_tot_O_cm2(1223) =    3.10
 given_sigma_tot_O_cm2(1224) =    3.10
 given_sigma_tot_O_cm2(1225) =    3.10
 given_sigma_tot_O_cm2(1226) =    3.10
 given_sigma_tot_O_cm2(1227) =    3.10
 given_sigma_tot_O_cm2(1228) =    3.10
 given_sigma_tot_O_cm2(1229) =    3.10
 given_sigma_tot_O_cm2(1230) =    3.10
 given_sigma_tot_O_cm2(1231) =    3.10
 given_sigma_tot_O_cm2(1232) =    3.10
 given_sigma_tot_O_cm2(1233) =    3.10
 given_sigma_tot_O_cm2(1234) =    3.10
 given_sigma_tot_O_cm2(1235) =    3.10
 given_sigma_tot_O_cm2(1236) =    3.10
 given_sigma_tot_O_cm2(1237) =    3.10
 given_sigma_tot_O_cm2(1238) =    3.10
 given_sigma_tot_O_cm2(1239) =    3.10
 given_sigma_tot_O_cm2(1240) =    3.10
 given_sigma_tot_O_cm2(1241) =    3.77
 given_sigma_tot_O_cm2(1242) =    3.94
 given_sigma_tot_O_cm2(1243) =    4.27
 given_sigma_tot_O_cm2(1244) =    4.77
 given_sigma_tot_O_cm2(1245) =    5.19
 given_sigma_tot_O_cm2(1246) =    5.61
 given_sigma_tot_O_cm2(1247) =    6.02
 given_sigma_tot_O_cm2(1248) =    6.44
 given_sigma_tot_O_cm2(1249) =    6.86
 given_sigma_tot_O_cm2(1250) =    7.28
 given_sigma_tot_O_cm2(1251) =    7.44
 given_sigma_tot_O_cm2(1252) =    7.70
 given_sigma_tot_O_cm2(1253) =    8.11
 given_sigma_tot_O_cm2(1254) =    8.28
 given_sigma_tot_O_cm2(1255) =    8.53
 given_sigma_tot_O_cm2(1256) =    8.95
 given_sigma_tot_O_cm2(1257) =    9.37
 given_sigma_tot_O_cm2(1258) =   11.82
 given_sigma_tot_O_cm2(1259) =   34.88
 given_sigma_tot_O_cm2(1260) =    9.93
 given_sigma_tot_O_cm2(1261) =    4.73
 given_sigma_tot_O_cm2(1262) =    6.81
 given_sigma_tot_O_cm2(1263) =    5.20
 given_sigma_tot_O_cm2(1264) =    5.11
 given_sigma_tot_O_cm2(1265) =    4.88
 given_sigma_tot_O_cm2(1266) =    4.66
 given_sigma_tot_O_cm2(1267) =    4.43
 given_sigma_tot_O_cm2(1268) =    4.21
 given_sigma_tot_O_cm2(1269) =    3.98
 given_sigma_tot_O_cm2(1270) =    3.89
 given_sigma_tot_O_cm2(1271) =    3.75
 given_sigma_tot_O_cm2(1272) =    3.62
 given_sigma_tot_O_cm2(1273) =    3.53
 given_sigma_tot_O_cm2(1274) =    3.30
 given_sigma_tot_O_cm2(1275) =    3.08
 given_sigma_tot_O_cm2(1276) =    2.85
 given_sigma_tot_O_cm2(1277) =    2.84
 given_sigma_tot_O_cm2(1278) =    2.83
 given_sigma_tot_O_cm2(1279) =    2.82
 given_sigma_tot_O_cm2(1280) =    2.81
 given_sigma_tot_O_cm2(1281) =    2.79
 given_sigma_tot_O_cm2(1282) =    2.78
 given_sigma_tot_O_cm2(1283) =    2.76
 given_sigma_tot_O_cm2(1284) =    2.75
 given_sigma_tot_O_cm2(1285) =    2.73
 given_sigma_tot_O_cm2(1286) =    2.72
 given_sigma_tot_O_cm2(1287) =    2.70
 given_sigma_tot_O_cm2(1288) =    2.71
 given_sigma_tot_O_cm2(1289) =    2.71
 given_sigma_tot_O_cm2(1290) =    2.72
 given_sigma_tot_O_cm2(1291) =    2.72
 given_sigma_tot_O_cm2(1292) =    2.73
 given_sigma_tot_O_cm2(1293) =    2.73
 given_sigma_tot_O_cm2(1294) =    2.74
 given_sigma_tot_O_cm2(1295) =    2.74
 given_sigma_tot_O_cm2(1296) =    2.75
 given_sigma_tot_O_cm2(1297) =    2.75
 given_sigma_tot_O_cm2(1298) =    2.75
 given_sigma_tot_O_cm2(1299) =    2.76
 given_sigma_tot_O_cm2(1300) =    2.77
 given_sigma_tot_O_cm2(1301) =    2.77
 given_sigma_tot_O_cm2(1302) =    2.78
 given_sigma_tot_O_cm2(1303) =    2.79
 given_sigma_tot_O_cm2(1304) =    2.79
 given_sigma_tot_O_cm2(1305) =    2.80
 given_sigma_tot_O_cm2(1306) =    2.81
 given_sigma_tot_O_cm2(1307) =    2.82
 given_sigma_tot_O_cm2(1308) =    2.82
 given_sigma_tot_O_cm2(1309) =    2.83
 given_sigma_tot_O_cm2(1310) =    2.84
 given_sigma_tot_O_cm2(1311) =    2.85
 given_sigma_tot_O_cm2(1312) =    2.85
 given_sigma_tot_O_cm2(1313) =    2.84
 given_sigma_tot_O_cm2(1314) =    2.83
 given_sigma_tot_O_cm2(1315) =    2.82
 given_sigma_tot_O_cm2(1316) =    2.82
 given_sigma_tot_O_cm2(1317) =    2.81
 given_sigma_tot_O_cm2(1318) =    2.79
 given_sigma_tot_O_cm2(1319) =    2.78
 given_sigma_tot_O_cm2(1320) =    2.77
 given_sigma_tot_O_cm2(1321) =    2.77
 given_sigma_tot_O_cm2(1322) =    2.76
 given_sigma_tot_O_cm2(1323) =    2.75
 given_sigma_tot_O_cm2(1324) =    2.75
 given_sigma_tot_O_cm2(1325) =    2.75
 given_sigma_tot_O_cm2(1326) =    2.75
 given_sigma_tot_O_cm2(1327) =    2.75
 given_sigma_tot_O_cm2(1328) =    2.75
 given_sigma_tot_O_cm2(1329) =    2.68
 given_sigma_tot_O_cm2(1330) =    2.62
 given_sigma_tot_O_cm2(1331) =    2.55
 given_sigma_tot_O_cm2(1332) =    2.49
 given_sigma_tot_O_cm2(1333) =    2.45
 given_sigma_tot_O_cm2(1334) =    1.15
 given_sigma_tot_O_cm2(1335) =    0.75
 given_sigma_tot_O_cm2(1336) =    0.49
 given_sigma_tot_O_cm2(1337) =    0.30
 given_sigma_tot_O_cm2(1338) =    0.20
 given_sigma_tot_O_cm2(1339) =    0.14
 given_sigma_tot_O_cm2(1340) =    0.06
 given_sigma_tot_O_cm2(1341) =    0.04
 given_sigma_tot_O_cm2(1342) =    0.00
 given_sigma_tot_O_cm2(1343) =    0.00
 given_sigma_tot_O_cm2(1344) =    0.00
 given_sigma_tot_O_cm2(1345) =    0.00
 given_sigma_tot_O_cm2(1346) =    0.00
 given_sigma_tot_O_cm2(1347) =    0.00
 given_sigma_tot_O_cm2(1348) =    0.00
 given_sigma_tot_O_cm2(1349) =    0.00
 given_sigma_tot_O_cm2(1350) =    0.00
 given_sigma_tot_O_cm2(1351) =    0.00
 given_sigma_tot_O_cm2(1352) =    0.00
 given_sigma_tot_O_cm2(1353) =    0.00
 given_sigma_tot_O_cm2(1354) =    0.00
 given_sigma_tot_O_cm2(1355) =    0.00
 given_sigma_tot_O_cm2(1356) =    0.00
 given_sigma_tot_O_cm2(1357) =    0.00
 given_sigma_tot_O_cm2(1358) =    0.00
 given_sigma_tot_O_cm2(1359) =    0.00
 given_sigma_tot_O_cm2(1360) =    0.00
 given_sigma_tot_O_cm2(1361) =    0.00
 given_sigma_tot_O_cm2(1362) =    0.00
 given_sigma_tot_O_cm2(1363) =    0.00
 given_sigma_tot_O_cm2(1364) =    0.00
 given_sigma_tot_O_cm2(1365) =    0.00
 given_sigma_tot_O_cm2(1366) =    0.00
 given_sigma_tot_O_cm2(1367) =    0.00
 given_sigma_tot_O_cm2(1368) =    0.00
 given_sigma_tot_O_cm2(1369) =    0.00
 given_sigma_tot_O_cm2(1370) =    0.00
 given_sigma_tot_O_cm2(1371) =    0.00
 given_sigma_tot_O_cm2(1372) =    0.00
 given_sigma_tot_O_cm2(1373) =    0.00
 given_sigma_tot_O_cm2(1374) =    0.00
 given_sigma_tot_O_cm2(1375) =    0.00
 given_sigma_tot_O_cm2(1376) =    0.00
 given_sigma_tot_O_cm2(1377) =    0.00
 given_sigma_tot_O_cm2(1378) =    0.00
 given_sigma_tot_O_cm2(1379) =    0.00
 given_sigma_tot_O_cm2(1380) =    0.00
 given_sigma_tot_O_cm2(1381) =    0.00
 given_sigma_tot_O_cm2(1382) =    0.00
 given_sigma_tot_O_cm2(1383) =    0.00
 given_sigma_tot_O_cm2(1384) =    0.00
 given_sigma_tot_O_cm2(1385) =    0.00
 given_sigma_tot_O_cm2(1386) =    0.00
 given_sigma_tot_O_cm2(1387) =    0.00
 given_sigma_tot_O_cm2(1388) =    0.00
 given_sigma_tot_O_cm2(1389) =    0.00
 given_sigma_tot_O_cm2(1390) =    0.00
 given_sigma_tot_O_cm2(1391) =    0.00
 given_sigma_tot_O_cm2(1392) =    0.00
 given_sigma_tot_O_cm2(1393) =    0.00
 given_sigma_tot_O_cm2(1394) =    0.00
 given_sigma_tot_O_cm2(1395) =    0.00
 given_sigma_tot_O_cm2(1396) =    0.00
 given_sigma_tot_O_cm2(1397) =    0.00
 given_sigma_tot_O_cm2(1398) =    0.00
 given_sigma_tot_O_cm2(1399) =    0.00
 given_sigma_tot_O_cm2(1400) =    0.00
 given_sigma_tot_O_cm2(1401) =    0.00
 given_sigma_tot_O_cm2(1402) =    0.00
 given_sigma_tot_O_cm2(1403) =    0.00
 given_sigma_tot_O_cm2(1404) =    0.00
 given_sigma_tot_O_cm2(1405) =    0.00
 given_sigma_tot_O_cm2(1406) =    0.00
 given_sigma_tot_O_cm2(1407) =    0.00
 given_sigma_tot_O_cm2(1408) =    0.00
 given_sigma_tot_O_cm2(1409) =    0.00
 given_sigma_tot_O_cm2(1410) =    0.00
 given_sigma_tot_O_cm2(1411) =    0.00
 given_sigma_tot_O_cm2(1412) =    0.00
 given_sigma_tot_O_cm2(1413) =    0.00
 given_sigma_tot_O_cm2(1414) =    0.00
 given_sigma_tot_O_cm2(1415) =    0.00
 given_sigma_tot_O_cm2(1416) =    0.00
 given_sigma_tot_O_cm2(1417) =    0.00
 given_sigma_tot_O_cm2(1418) =    0.00
 given_sigma_tot_O_cm2(1419) =    0.00
 given_sigma_tot_O_cm2(1420) =    0.00
 given_sigma_tot_O_cm2(1421) =    0.00
 given_sigma_tot_O_cm2(1422) =    0.00
 given_sigma_tot_O_cm2(1423) =    0.00
 given_sigma_tot_O_cm2(1424) =    0.00
 given_sigma_tot_O_cm2(1425) =    0.00
 given_sigma_tot_O_cm2(1426) =    0.00
 given_sigma_tot_O_cm2(1427) =    0.00
 given_sigma_tot_O_cm2(1428) =    0.00
 given_sigma_tot_O_cm2(1429) =    0.00
 given_sigma_tot_O_cm2(1430) =    0.00
 given_sigma_tot_O_cm2(1431) =    0.00
 given_sigma_tot_O_cm2(1432) =    0.00
 given_sigma_tot_O_cm2(1433) =    0.00
 given_sigma_tot_O_cm2(1434) =    0.00
 given_sigma_tot_O_cm2(1435) =    0.00
 given_sigma_tot_O_cm2(1436) =    0.00
 given_sigma_tot_O_cm2(1437) =    0.00
 given_sigma_tot_O_cm2(1438) =    0.00
 given_sigma_tot_O_cm2(1439) =    0.00
 given_sigma_tot_O_cm2(1440) =    0.00
 given_sigma_tot_O_cm2(1441) =    0.00
 given_sigma_tot_O_cm2(1442) =    0.00
 given_sigma_tot_O_cm2(1443) =    0.00
 given_sigma_tot_O_cm2(1444) =    0.00
 given_sigma_tot_O_cm2(1445) =    0.00
 given_sigma_tot_O_cm2(1446) =    0.00
 given_sigma_tot_O_cm2(1447) =    0.00
 given_sigma_tot_O_cm2(1448) =    0.00
 given_sigma_tot_O_cm2(1449) =    0.00
 given_sigma_tot_O_cm2(1450) =    0.00
 given_sigma_tot_O_cm2(1451) =    0.00
 given_sigma_tot_O_cm2(1452) =    0.00
 given_sigma_tot_O_cm2(1453) =    0.00
 given_sigma_tot_O_cm2(1454) =    0.00
 given_sigma_tot_O_cm2(1455) =    0.00
 given_sigma_tot_O_cm2(1456) =    0.00
 given_sigma_tot_O_cm2(1457) =    0.00
 given_sigma_tot_O_cm2(1458) =    0.00
 given_sigma_tot_O_cm2(1459) =    0.00
 given_sigma_tot_O_cm2(1460) =    0.00
 given_sigma_tot_O_cm2(1461) =    0.00
 given_sigma_tot_O_cm2(1462) =    0.00
 given_sigma_tot_O_cm2(1463) =    0.00
 given_sigma_tot_O_cm2(1464) =    0.00
 given_sigma_tot_O_cm2(1465) =    0.00
 given_sigma_tot_O_cm2(1466) =    0.00
 given_sigma_tot_O_cm2(1467) =    0.00
 given_sigma_tot_O_cm2(1468) =    0.00
 given_sigma_tot_O_cm2(1469) =    0.00
 given_sigma_tot_O_cm2(1470) =    0.00
 given_sigma_tot_O_cm2(1471) =    0.00
 given_sigma_tot_O_cm2(1472) =    0.00
 given_sigma_tot_O_cm2(1473) =    0.00
 given_sigma_tot_O_cm2(1474) =    0.00
 given_sigma_tot_O_cm2(1475) =    0.00
 given_sigma_tot_O_cm2(1476) =    0.00
 given_sigma_tot_O_cm2(1477) =    0.00
 given_sigma_tot_O_cm2(1478) =    0.00
 given_sigma_tot_O_cm2(1479) =    0.00
 given_sigma_tot_O_cm2(1480) =    0.00
 given_sigma_tot_O_cm2(1481) =    0.00
 given_sigma_tot_O_cm2(1482) =    0.00
 given_sigma_tot_O_cm2(1483) =    0.00
 given_sigma_tot_O_cm2(1484) =    0.00
 given_sigma_tot_O_cm2(1485) =    0.00
 given_sigma_tot_O_cm2(1486) =    0.00
 given_sigma_tot_O_cm2(1487) =    0.00
 given_sigma_tot_O_cm2(1488) =    0.00
 given_sigma_tot_O_cm2(1489) =    0.00
 given_sigma_tot_O_cm2(1490) =    0.00
 given_sigma_tot_O_cm2(1491) =    0.00
 given_sigma_tot_O_cm2(1492) =    0.00
 given_sigma_tot_O_cm2(1493) =    0.00
 given_sigma_tot_O_cm2(1494) =    0.00
 given_sigma_tot_O_cm2(1495) =    0.00
 given_sigma_tot_O_cm2(1496) =    0.00
 given_sigma_tot_O_cm2(1497) =    0.00
 given_sigma_tot_O_cm2(1498) =    0.00
 given_sigma_tot_O_cm2(1499) =    0.00
 given_sigma_tot_O_cm2(1500) =    0.00
 given_sigma_tot_O_cm2(1501) =    0.00
 given_sigma_tot_O_cm2(1502) =    0.00
 given_sigma_tot_O_cm2(1503) =    0.00
 given_sigma_tot_O_cm2(1504) =    0.00
 given_sigma_tot_O_cm2(1505) =    0.00
 given_sigma_tot_O_cm2(1506) =    0.00
 given_sigma_tot_O_cm2(1507) =    0.00
 given_sigma_tot_O_cm2(1508) =    0.00
 given_sigma_tot_O_cm2(1509) =    0.00
 given_sigma_tot_O_cm2(1510) =    0.00
 given_sigma_tot_O_cm2(1511) =    0.00
 given_sigma_tot_O_cm2(1512) =    0.00
 given_sigma_tot_O_cm2(1513) =    0.00
 given_sigma_tot_O_cm2(1514) =    0.00
 given_sigma_tot_O_cm2(1515) =    0.00
 given_sigma_tot_O_cm2(1516) =    0.00
 given_sigma_tot_O_cm2(1517) =    0.00
 given_sigma_tot_O_cm2(1518) =    0.00
 given_sigma_tot_O_cm2(1519) =    0.00
 given_sigma_tot_O_cm2(1520) =    0.00
 given_sigma_tot_O_cm2(1521) =    0.00
 given_sigma_tot_O_cm2(1522) =    0.00
 given_sigma_tot_O_cm2(1523) =    0.00
 given_sigma_tot_O_cm2(1524) =    0.00
 given_sigma_tot_O_cm2(1525) =    0.00
 given_sigma_tot_O_cm2(1526) =    0.00
 given_sigma_tot_O_cm2(1527) =    0.00
 given_sigma_tot_O_cm2(1528) =    0.00
 given_sigma_tot_O_cm2(1529) =    0.00
 given_sigma_tot_O_cm2(1530) =    0.00
 given_sigma_tot_O_cm2(1531) =    0.00
 given_sigma_tot_O_cm2(1532) =    0.00
 given_sigma_tot_O_cm2(1533) =    0.00
 given_sigma_tot_O_cm2(1534) =    0.00
 given_sigma_tot_O_cm2(1535) =    0.00
 given_sigma_tot_O_cm2(1536) =    0.00
 given_sigma_tot_O_cm2(1537) =    0.00
 given_sigma_tot_O_cm2(1538) =    0.00
 given_sigma_tot_O_cm2(1539) =    0.00
 given_sigma_tot_O_cm2(1540) =    0.00
 given_sigma_tot_O_cm2(1541) =    0.00
 given_sigma_tot_O_cm2(1542) =    0.00
 given_sigma_tot_O_cm2(1543) =    0.00
 given_sigma_tot_O_cm2(1544) =    0.00
 given_sigma_tot_O_cm2(1545) =    0.00
 given_sigma_tot_O_cm2(1546) =    0.00
 given_sigma_tot_O_cm2(1547) =    0.00
 given_sigma_tot_O_cm2(1548) =    0.00
 given_sigma_tot_O_cm2(1549) =    0.00
 given_sigma_tot_O_cm2(1550) =    0.00
 given_sigma_tot_O_cm2(1551) =    0.00
 given_sigma_tot_O_cm2(1552) =    0.00
 given_sigma_tot_O_cm2(1553) =    0.00
 given_sigma_tot_O_cm2(1554) =    0.00
 given_sigma_tot_O_cm2(1555) =    0.00
 given_sigma_tot_O_cm2(1556) =    0.00
 given_sigma_tot_O_cm2(1557) =    0.00
 given_sigma_tot_O_cm2(1558) =    0.00
 given_sigma_tot_O_cm2(1559) =    0.00
 given_sigma_tot_O_cm2(1560) =    0.00
 given_sigma_tot_O_cm2(1561) =    0.00
 given_sigma_tot_O_cm2(1562) =    0.00
 given_sigma_tot_O_cm2(1563) =    0.00
 given_sigma_tot_O_cm2(1564) =    0.00
 given_sigma_tot_O_cm2(1565) =    0.00
 given_sigma_tot_O_cm2(1566) =    0.00
 given_sigma_tot_O_cm2(1567) =    0.00
 given_sigma_tot_O_cm2(1568) =    0.00
 given_sigma_tot_O_cm2(1569) =    0.00
 given_sigma_tot_O_cm2(1570) =    0.00
 given_sigma_tot_O_cm2(1571) =    0.00







 given_sigma_ion_N2_to_N2_cm2(   1) =    0.80 ;  given_sigma_ion_N2_to_N_cm2(   1) =    0.00
 given_sigma_ion_N2_to_N2_cm2(   2) =    0.96 ;  given_sigma_ion_N2_to_N_cm2(   2) =    0.00
 given_sigma_ion_N2_to_N2_cm2(   3) =    0.97 ;  given_sigma_ion_N2_to_N_cm2(   3) =    0.00
 given_sigma_ion_N2_to_N2_cm2(   4) =    0.99 ;  given_sigma_ion_N2_to_N_cm2(   4) =    0.00
 given_sigma_ion_N2_to_N2_cm2(   5) =    0.61 ;  given_sigma_ion_N2_to_N_cm2(   5) =    0.40
 given_sigma_ion_N2_to_N2_cm2(   6) =    0.61 ;  given_sigma_ion_N2_to_N_cm2(   6) =    0.41
 given_sigma_ion_N2_to_N2_cm2(   7) =    0.55 ;  given_sigma_ion_N2_to_N_cm2(   7) =    0.36
 given_sigma_ion_N2_to_N2_cm2(   8) =    0.31 ;  given_sigma_ion_N2_to_N_cm2(   8) =    0.20
 given_sigma_ion_N2_to_N2_cm2(   9) =    0.05 ;  given_sigma_ion_N2_to_N_cm2(   9) =    0.03
 given_sigma_ion_N2_to_N2_cm2(  10) =    0.06 ;  given_sigma_ion_N2_to_N_cm2(  10) =    0.04
 given_sigma_ion_N2_to_N2_cm2(  11) =    0.06 ;  given_sigma_ion_N2_to_N_cm2(  11) =    0.04
 given_sigma_ion_N2_to_N2_cm2(  12) =    0.09 ;  given_sigma_ion_N2_to_N_cm2(  12) =    0.06
 given_sigma_ion_N2_to_N2_cm2(  13) =    0.09 ;  given_sigma_ion_N2_to_N_cm2(  13) =    0.06
 given_sigma_ion_N2_to_N2_cm2(  14) =    0.10 ;  given_sigma_ion_N2_to_N_cm2(  14) =    0.07
 given_sigma_ion_N2_to_N2_cm2(  15) =    0.10 ;  given_sigma_ion_N2_to_N_cm2(  15) =    0.07
 given_sigma_ion_N2_to_N2_cm2(  16) =    0.10 ;  given_sigma_ion_N2_to_N_cm2(  16) =    0.07
 given_sigma_ion_N2_to_N2_cm2(  17) =    0.11 ;  given_sigma_ion_N2_to_N_cm2(  17) =    0.07
 given_sigma_ion_N2_to_N2_cm2(  18) =    0.11 ;  given_sigma_ion_N2_to_N_cm2(  18) =    0.07
 given_sigma_ion_N2_to_N2_cm2(  19) =    0.11 ;  given_sigma_ion_N2_to_N_cm2(  19) =    0.08
 given_sigma_ion_N2_to_N2_cm2(  20) =    0.11 ;  given_sigma_ion_N2_to_N_cm2(  20) =    0.08
 given_sigma_ion_N2_to_N2_cm2(  21) =    0.12 ;  given_sigma_ion_N2_to_N_cm2(  21) =    0.08
 given_sigma_ion_N2_to_N2_cm2(  22) =    0.12 ;  given_sigma_ion_N2_to_N_cm2(  22) =    0.08
 given_sigma_ion_N2_to_N2_cm2(  23) =    0.13 ;  given_sigma_ion_N2_to_N_cm2(  23) =    0.08
 given_sigma_ion_N2_to_N2_cm2(  24) =    0.13 ;  given_sigma_ion_N2_to_N_cm2(  24) =    0.08
 given_sigma_ion_N2_to_N2_cm2(  25) =    0.14 ;  given_sigma_ion_N2_to_N_cm2(  25) =    0.09
 given_sigma_ion_N2_to_N2_cm2(  26) =    0.14 ;  given_sigma_ion_N2_to_N_cm2(  26) =    0.09
 given_sigma_ion_N2_to_N2_cm2(  27) =    0.14 ;  given_sigma_ion_N2_to_N_cm2(  27) =    0.09
 given_sigma_ion_N2_to_N2_cm2(  28) =    0.14 ;  given_sigma_ion_N2_to_N_cm2(  28) =    0.09
 given_sigma_ion_N2_to_N2_cm2(  29) =    0.15 ;  given_sigma_ion_N2_to_N_cm2(  29) =    0.10
 given_sigma_ion_N2_to_N2_cm2(  30) =    0.15 ;  given_sigma_ion_N2_to_N_cm2(  30) =    0.10
 given_sigma_ion_N2_to_N2_cm2(  31) =    0.16 ;  given_sigma_ion_N2_to_N_cm2(  31) =    0.10
 given_sigma_ion_N2_to_N2_cm2(  32) =    0.16 ;  given_sigma_ion_N2_to_N_cm2(  32) =    0.11
 given_sigma_ion_N2_to_N2_cm2(  33) =    0.17 ;  given_sigma_ion_N2_to_N_cm2(  33) =    0.11
 given_sigma_ion_N2_to_N2_cm2(  34) =    0.17 ;  given_sigma_ion_N2_to_N_cm2(  34) =    0.11
 given_sigma_ion_N2_to_N2_cm2(  35) =    0.17 ;  given_sigma_ion_N2_to_N_cm2(  35) =    0.11
 given_sigma_ion_N2_to_N2_cm2(  36) =    0.18 ;  given_sigma_ion_N2_to_N_cm2(  36) =    0.12
 given_sigma_ion_N2_to_N2_cm2(  37) =    0.18 ;  given_sigma_ion_N2_to_N_cm2(  37) =    0.12
 given_sigma_ion_N2_to_N2_cm2(  38) =    0.18 ;  given_sigma_ion_N2_to_N_cm2(  38) =    0.12
 given_sigma_ion_N2_to_N2_cm2(  39) =    0.19 ;  given_sigma_ion_N2_to_N_cm2(  39) =    0.12
 given_sigma_ion_N2_to_N2_cm2(  40) =    0.19 ;  given_sigma_ion_N2_to_N_cm2(  40) =    0.12
 given_sigma_ion_N2_to_N2_cm2(  41) =    0.19 ;  given_sigma_ion_N2_to_N_cm2(  41) =    0.13
 given_sigma_ion_N2_to_N2_cm2(  42) =    0.20 ;  given_sigma_ion_N2_to_N_cm2(  42) =    0.13
 given_sigma_ion_N2_to_N2_cm2(  43) =    0.20 ;  given_sigma_ion_N2_to_N_cm2(  43) =    0.13
 given_sigma_ion_N2_to_N2_cm2(  44) =    0.20 ;  given_sigma_ion_N2_to_N_cm2(  44) =    0.13
 given_sigma_ion_N2_to_N2_cm2(  45) =    0.21 ;  given_sigma_ion_N2_to_N_cm2(  45) =    0.14
 given_sigma_ion_N2_to_N2_cm2(  46) =    0.22 ;  given_sigma_ion_N2_to_N_cm2(  46) =    0.14
 given_sigma_ion_N2_to_N2_cm2(  47) =    0.22 ;  given_sigma_ion_N2_to_N_cm2(  47) =    0.15
 given_sigma_ion_N2_to_N2_cm2(  48) =    0.23 ;  given_sigma_ion_N2_to_N_cm2(  48) =    0.15
 given_sigma_ion_N2_to_N2_cm2(  49) =    0.23 ;  given_sigma_ion_N2_to_N_cm2(  49) =    0.15
 given_sigma_ion_N2_to_N2_cm2(  50) =    0.24 ;  given_sigma_ion_N2_to_N_cm2(  50) =    0.15
 given_sigma_ion_N2_to_N2_cm2(  51) =    0.24 ;  given_sigma_ion_N2_to_N_cm2(  51) =    0.15
 given_sigma_ion_N2_to_N2_cm2(  52) =    0.24 ;  given_sigma_ion_N2_to_N_cm2(  52) =    0.16
 given_sigma_ion_N2_to_N2_cm2(  53) =    0.25 ;  given_sigma_ion_N2_to_N_cm2(  53) =    0.16
 given_sigma_ion_N2_to_N2_cm2(  54) =    0.25 ;  given_sigma_ion_N2_to_N_cm2(  54) =    0.16
 given_sigma_ion_N2_to_N2_cm2(  55) =    0.25 ;  given_sigma_ion_N2_to_N_cm2(  55) =    0.16
 given_sigma_ion_N2_to_N2_cm2(  56) =    0.25 ;  given_sigma_ion_N2_to_N_cm2(  56) =    0.16
 given_sigma_ion_N2_to_N2_cm2(  57) =    0.25 ;  given_sigma_ion_N2_to_N_cm2(  57) =    0.16
 given_sigma_ion_N2_to_N2_cm2(  58) =    0.26 ;  given_sigma_ion_N2_to_N_cm2(  58) =    0.17
 given_sigma_ion_N2_to_N2_cm2(  59) =    0.26 ;  given_sigma_ion_N2_to_N_cm2(  59) =    0.17
 given_sigma_ion_N2_to_N2_cm2(  60) =    0.26 ;  given_sigma_ion_N2_to_N_cm2(  60) =    0.17
 given_sigma_ion_N2_to_N2_cm2(  61) =    0.27 ;  given_sigma_ion_N2_to_N_cm2(  61) =    0.17
 given_sigma_ion_N2_to_N2_cm2(  62) =    0.27 ;  given_sigma_ion_N2_to_N_cm2(  62) =    0.17
 given_sigma_ion_N2_to_N2_cm2(  63) =    0.27 ;  given_sigma_ion_N2_to_N_cm2(  63) =    0.17
 given_sigma_ion_N2_to_N2_cm2(  64) =    0.27 ;  given_sigma_ion_N2_to_N_cm2(  64) =    0.18
 given_sigma_ion_N2_to_N2_cm2(  65) =    0.28 ;  given_sigma_ion_N2_to_N_cm2(  65) =    0.18
 given_sigma_ion_N2_to_N2_cm2(  66) =    0.28 ;  given_sigma_ion_N2_to_N_cm2(  66) =    0.18
 given_sigma_ion_N2_to_N2_cm2(  67) =    0.29 ;  given_sigma_ion_N2_to_N_cm2(  67) =    0.18
 given_sigma_ion_N2_to_N2_cm2(  68) =    0.29 ;  given_sigma_ion_N2_to_N_cm2(  68) =    0.19
 given_sigma_ion_N2_to_N2_cm2(  69) =    0.30 ;  given_sigma_ion_N2_to_N_cm2(  69) =    0.19
 given_sigma_ion_N2_to_N2_cm2(  70) =    0.31 ;  given_sigma_ion_N2_to_N_cm2(  70) =    0.20
 given_sigma_ion_N2_to_N2_cm2(  71) =    0.32 ;  given_sigma_ion_N2_to_N_cm2(  71) =    0.20
 given_sigma_ion_N2_to_N2_cm2(  72) =    0.32 ;  given_sigma_ion_N2_to_N_cm2(  72) =    0.21
 given_sigma_ion_N2_to_N2_cm2(  73) =    0.33 ;  given_sigma_ion_N2_to_N_cm2(  73) =    0.21
 given_sigma_ion_N2_to_N2_cm2(  74) =    0.33 ;  given_sigma_ion_N2_to_N_cm2(  74) =    0.21
 given_sigma_ion_N2_to_N2_cm2(  75) =    0.34 ;  given_sigma_ion_N2_to_N_cm2(  75) =    0.22
 given_sigma_ion_N2_to_N2_cm2(  76) =    0.34 ;  given_sigma_ion_N2_to_N_cm2(  76) =    0.22
 given_sigma_ion_N2_to_N2_cm2(  77) =    0.35 ;  given_sigma_ion_N2_to_N_cm2(  77) =    0.22
 given_sigma_ion_N2_to_N2_cm2(  78) =    0.35 ;  given_sigma_ion_N2_to_N_cm2(  78) =    0.22
 given_sigma_ion_N2_to_N2_cm2(  79) =    0.35 ;  given_sigma_ion_N2_to_N_cm2(  79) =    0.22
 given_sigma_ion_N2_to_N2_cm2(  80) =    0.36 ;  given_sigma_ion_N2_to_N_cm2(  80) =    0.23
 given_sigma_ion_N2_to_N2_cm2(  81) =    0.37 ;  given_sigma_ion_N2_to_N_cm2(  81) =    0.23
 given_sigma_ion_N2_to_N2_cm2(  82) =    0.37 ;  given_sigma_ion_N2_to_N_cm2(  82) =    0.23
 given_sigma_ion_N2_to_N2_cm2(  83) =    0.37 ;  given_sigma_ion_N2_to_N_cm2(  83) =    0.24
 given_sigma_ion_N2_to_N2_cm2(  84) =    0.38 ;  given_sigma_ion_N2_to_N_cm2(  84) =    0.24
 given_sigma_ion_N2_to_N2_cm2(  85) =    0.38 ;  given_sigma_ion_N2_to_N_cm2(  85) =    0.24
 given_sigma_ion_N2_to_N2_cm2(  86) =    0.38 ;  given_sigma_ion_N2_to_N_cm2(  86) =    0.24
 given_sigma_ion_N2_to_N2_cm2(  87) =    0.38 ;  given_sigma_ion_N2_to_N_cm2(  87) =    0.24
 given_sigma_ion_N2_to_N2_cm2(  88) =    0.39 ;  given_sigma_ion_N2_to_N_cm2(  88) =    0.25
 given_sigma_ion_N2_to_N2_cm2(  89) =    0.39 ;  given_sigma_ion_N2_to_N_cm2(  89) =    0.25
 given_sigma_ion_N2_to_N2_cm2(  90) =    0.40 ;  given_sigma_ion_N2_to_N_cm2(  90) =    0.25
 given_sigma_ion_N2_to_N2_cm2(  91) =    0.40 ;  given_sigma_ion_N2_to_N_cm2(  91) =    0.25
 given_sigma_ion_N2_to_N2_cm2(  92) =    0.40 ;  given_sigma_ion_N2_to_N_cm2(  92) =    0.25
 given_sigma_ion_N2_to_N2_cm2(  93) =    0.40 ;  given_sigma_ion_N2_to_N_cm2(  93) =    0.26
 given_sigma_ion_N2_to_N2_cm2(  94) =    0.41 ;  given_sigma_ion_N2_to_N_cm2(  94) =    0.26
 given_sigma_ion_N2_to_N2_cm2(  95) =    0.42 ;  given_sigma_ion_N2_to_N_cm2(  95) =    0.26
 given_sigma_ion_N2_to_N2_cm2(  96) =    0.42 ;  given_sigma_ion_N2_to_N_cm2(  96) =    0.26
 given_sigma_ion_N2_to_N2_cm2(  97) =    0.42 ;  given_sigma_ion_N2_to_N_cm2(  97) =    0.27
 given_sigma_ion_N2_to_N2_cm2(  98) =    0.43 ;  given_sigma_ion_N2_to_N_cm2(  98) =    0.27
 given_sigma_ion_N2_to_N2_cm2(  99) =    0.43 ;  given_sigma_ion_N2_to_N_cm2(  99) =    0.27
 given_sigma_ion_N2_to_N2_cm2( 100) =    0.43 ;  given_sigma_ion_N2_to_N_cm2( 100) =    0.27
 given_sigma_ion_N2_to_N2_cm2( 101) =    0.43 ;  given_sigma_ion_N2_to_N_cm2( 101) =    0.27
 given_sigma_ion_N2_to_N2_cm2( 102) =    0.44 ;  given_sigma_ion_N2_to_N_cm2( 102) =    0.28
 given_sigma_ion_N2_to_N2_cm2( 103) =    0.44 ;  given_sigma_ion_N2_to_N_cm2( 103) =    0.28
 given_sigma_ion_N2_to_N2_cm2( 104) =    0.45 ;  given_sigma_ion_N2_to_N_cm2( 104) =    0.28
 given_sigma_ion_N2_to_N2_cm2( 105) =    0.45 ;  given_sigma_ion_N2_to_N_cm2( 105) =    0.28
 given_sigma_ion_N2_to_N2_cm2( 106) =    0.45 ;  given_sigma_ion_N2_to_N_cm2( 106) =    0.28
 given_sigma_ion_N2_to_N2_cm2( 107) =    0.46 ;  given_sigma_ion_N2_to_N_cm2( 107) =    0.29
 given_sigma_ion_N2_to_N2_cm2( 108) =    0.46 ;  given_sigma_ion_N2_to_N_cm2( 108) =    0.29
 given_sigma_ion_N2_to_N2_cm2( 109) =    0.47 ;  given_sigma_ion_N2_to_N_cm2( 109) =    0.29
 given_sigma_ion_N2_to_N2_cm2( 110) =    0.48 ;  given_sigma_ion_N2_to_N_cm2( 110) =    0.30
 given_sigma_ion_N2_to_N2_cm2( 111) =    0.48 ;  given_sigma_ion_N2_to_N_cm2( 111) =    0.30
 given_sigma_ion_N2_to_N2_cm2( 112) =    0.49 ;  given_sigma_ion_N2_to_N_cm2( 112) =    0.30
 given_sigma_ion_N2_to_N2_cm2( 113) =    0.49 ;  given_sigma_ion_N2_to_N_cm2( 113) =    0.31
 given_sigma_ion_N2_to_N2_cm2( 114) =    0.50 ;  given_sigma_ion_N2_to_N_cm2( 114) =    0.31
 given_sigma_ion_N2_to_N2_cm2( 115) =    0.50 ;  given_sigma_ion_N2_to_N_cm2( 115) =    0.31
 given_sigma_ion_N2_to_N2_cm2( 116) =    0.51 ;  given_sigma_ion_N2_to_N_cm2( 116) =    0.32
 given_sigma_ion_N2_to_N2_cm2( 117) =    0.51 ;  given_sigma_ion_N2_to_N_cm2( 117) =    0.32
 given_sigma_ion_N2_to_N2_cm2( 118) =    0.52 ;  given_sigma_ion_N2_to_N_cm2( 118) =    0.32
 given_sigma_ion_N2_to_N2_cm2( 119) =    0.52 ;  given_sigma_ion_N2_to_N_cm2( 119) =    0.32
 given_sigma_ion_N2_to_N2_cm2( 120) =    0.53 ;  given_sigma_ion_N2_to_N_cm2( 120) =    0.33
 given_sigma_ion_N2_to_N2_cm2( 121) =    0.53 ;  given_sigma_ion_N2_to_N_cm2( 121) =    0.33
 given_sigma_ion_N2_to_N2_cm2( 122) =    0.54 ;  given_sigma_ion_N2_to_N_cm2( 122) =    0.33
 given_sigma_ion_N2_to_N2_cm2( 123) =    0.54 ;  given_sigma_ion_N2_to_N_cm2( 123) =    0.34
 given_sigma_ion_N2_to_N2_cm2( 124) =    0.55 ;  given_sigma_ion_N2_to_N_cm2( 124) =    0.34
 given_sigma_ion_N2_to_N2_cm2( 125) =    0.55 ;  given_sigma_ion_N2_to_N_cm2( 125) =    0.34
 given_sigma_ion_N2_to_N2_cm2( 126) =    0.56 ;  given_sigma_ion_N2_to_N_cm2( 126) =    0.35
 given_sigma_ion_N2_to_N2_cm2( 127) =    0.56 ;  given_sigma_ion_N2_to_N_cm2( 127) =    0.35
 given_sigma_ion_N2_to_N2_cm2( 128) =    0.57 ;  given_sigma_ion_N2_to_N_cm2( 128) =    0.35
 given_sigma_ion_N2_to_N2_cm2( 129) =    0.58 ;  given_sigma_ion_N2_to_N_cm2( 129) =    0.36
 given_sigma_ion_N2_to_N2_cm2( 130) =    0.58 ;  given_sigma_ion_N2_to_N_cm2( 130) =    0.36
 given_sigma_ion_N2_to_N2_cm2( 131) =    0.59 ;  given_sigma_ion_N2_to_N_cm2( 131) =    0.36
 given_sigma_ion_N2_to_N2_cm2( 132) =    0.59 ;  given_sigma_ion_N2_to_N_cm2( 132) =    0.37
 given_sigma_ion_N2_to_N2_cm2( 133) =    0.60 ;  given_sigma_ion_N2_to_N_cm2( 133) =    0.37
 given_sigma_ion_N2_to_N2_cm2( 134) =    0.61 ;  given_sigma_ion_N2_to_N_cm2( 134) =    0.38
 given_sigma_ion_N2_to_N2_cm2( 135) =    0.62 ;  given_sigma_ion_N2_to_N_cm2( 135) =    0.38
 given_sigma_ion_N2_to_N2_cm2( 136) =    0.62 ;  given_sigma_ion_N2_to_N_cm2( 136) =    0.38
 given_sigma_ion_N2_to_N2_cm2( 137) =    0.63 ;  given_sigma_ion_N2_to_N_cm2( 137) =    0.39
 given_sigma_ion_N2_to_N2_cm2( 138) =    0.64 ;  given_sigma_ion_N2_to_N_cm2( 138) =    0.39
 given_sigma_ion_N2_to_N2_cm2( 139) =    0.64 ;  given_sigma_ion_N2_to_N_cm2( 139) =    0.40
 given_sigma_ion_N2_to_N2_cm2( 140) =    0.64 ;  given_sigma_ion_N2_to_N_cm2( 140) =    0.40
 given_sigma_ion_N2_to_N2_cm2( 141) =    0.65 ;  given_sigma_ion_N2_to_N_cm2( 141) =    0.40
 given_sigma_ion_N2_to_N2_cm2( 142) =    0.66 ;  given_sigma_ion_N2_to_N_cm2( 142) =    0.41
 given_sigma_ion_N2_to_N2_cm2( 143) =    0.66 ;  given_sigma_ion_N2_to_N_cm2( 143) =    0.41
 given_sigma_ion_N2_to_N2_cm2( 144) =    0.68 ;  given_sigma_ion_N2_to_N_cm2( 144) =    0.42
 given_sigma_ion_N2_to_N2_cm2( 145) =    0.69 ;  given_sigma_ion_N2_to_N_cm2( 145) =    0.42
 given_sigma_ion_N2_to_N2_cm2( 146) =    0.69 ;  given_sigma_ion_N2_to_N_cm2( 146) =    0.42
 given_sigma_ion_N2_to_N2_cm2( 147) =    0.69 ;  given_sigma_ion_N2_to_N_cm2( 147) =    0.43
 given_sigma_ion_N2_to_N2_cm2( 148) =    0.70 ;  given_sigma_ion_N2_to_N_cm2( 148) =    0.43
 given_sigma_ion_N2_to_N2_cm2( 149) =    0.71 ;  given_sigma_ion_N2_to_N_cm2( 149) =    0.44
 given_sigma_ion_N2_to_N2_cm2( 150) =    0.72 ;  given_sigma_ion_N2_to_N_cm2( 150) =    0.44
 given_sigma_ion_N2_to_N2_cm2( 151) =    0.72 ;  given_sigma_ion_N2_to_N_cm2( 151) =    0.44
 given_sigma_ion_N2_to_N2_cm2( 152) =    0.73 ;  given_sigma_ion_N2_to_N_cm2( 152) =    0.44
 given_sigma_ion_N2_to_N2_cm2( 153) =    0.74 ;  given_sigma_ion_N2_to_N_cm2( 153) =    0.45
 given_sigma_ion_N2_to_N2_cm2( 154) =    0.74 ;  given_sigma_ion_N2_to_N_cm2( 154) =    0.45
 given_sigma_ion_N2_to_N2_cm2( 155) =    0.75 ;  given_sigma_ion_N2_to_N_cm2( 155) =    0.46
 given_sigma_ion_N2_to_N2_cm2( 156) =    0.76 ;  given_sigma_ion_N2_to_N_cm2( 156) =    0.46
 given_sigma_ion_N2_to_N2_cm2( 157) =    0.76 ;  given_sigma_ion_N2_to_N_cm2( 157) =    0.46
 given_sigma_ion_N2_to_N2_cm2( 158) =    0.77 ;  given_sigma_ion_N2_to_N_cm2( 158) =    0.47
 given_sigma_ion_N2_to_N2_cm2( 159) =    0.77 ;  given_sigma_ion_N2_to_N_cm2( 159) =    0.47
 given_sigma_ion_N2_to_N2_cm2( 160) =    0.78 ;  given_sigma_ion_N2_to_N_cm2( 160) =    0.47
 given_sigma_ion_N2_to_N2_cm2( 161) =    0.78 ;  given_sigma_ion_N2_to_N_cm2( 161) =    0.48
 given_sigma_ion_N2_to_N2_cm2( 162) =    0.79 ;  given_sigma_ion_N2_to_N_cm2( 162) =    0.48
 given_sigma_ion_N2_to_N2_cm2( 163) =    0.80 ;  given_sigma_ion_N2_to_N_cm2( 163) =    0.49
 given_sigma_ion_N2_to_N2_cm2( 164) =    0.80 ;  given_sigma_ion_N2_to_N_cm2( 164) =    0.49
 given_sigma_ion_N2_to_N2_cm2( 165) =    0.82 ;  given_sigma_ion_N2_to_N_cm2( 165) =    0.49
 given_sigma_ion_N2_to_N2_cm2( 166) =    0.82 ;  given_sigma_ion_N2_to_N_cm2( 166) =    0.50
 given_sigma_ion_N2_to_N2_cm2( 167) =    0.84 ;  given_sigma_ion_N2_to_N_cm2( 167) =    0.51
 given_sigma_ion_N2_to_N2_cm2( 168) =    0.85 ;  given_sigma_ion_N2_to_N_cm2( 168) =    0.51
 given_sigma_ion_N2_to_N2_cm2( 169) =    0.86 ;  given_sigma_ion_N2_to_N_cm2( 169) =    0.52
 given_sigma_ion_N2_to_N2_cm2( 170) =    0.87 ;  given_sigma_ion_N2_to_N_cm2( 170) =    0.52
 given_sigma_ion_N2_to_N2_cm2( 171) =    0.87 ;  given_sigma_ion_N2_to_N_cm2( 171) =    0.53
 given_sigma_ion_N2_to_N2_cm2( 172) =    0.88 ;  given_sigma_ion_N2_to_N_cm2( 172) =    0.53
 given_sigma_ion_N2_to_N2_cm2( 173) =    0.88 ;  given_sigma_ion_N2_to_N_cm2( 173) =    0.53
 given_sigma_ion_N2_to_N2_cm2( 174) =    0.89 ;  given_sigma_ion_N2_to_N_cm2( 174) =    0.54
 given_sigma_ion_N2_to_N2_cm2( 175) =    0.90 ;  given_sigma_ion_N2_to_N_cm2( 175) =    0.54
 given_sigma_ion_N2_to_N2_cm2( 176) =    0.91 ;  given_sigma_ion_N2_to_N_cm2( 176) =    0.55
 given_sigma_ion_N2_to_N2_cm2( 177) =    0.93 ;  given_sigma_ion_N2_to_N_cm2( 177) =    0.56
 given_sigma_ion_N2_to_N2_cm2( 178) =    0.93 ;  given_sigma_ion_N2_to_N_cm2( 178) =    0.56
 given_sigma_ion_N2_to_N2_cm2( 179) =    0.94 ;  given_sigma_ion_N2_to_N_cm2( 179) =    0.57
 given_sigma_ion_N2_to_N2_cm2( 180) =    0.96 ;  given_sigma_ion_N2_to_N_cm2( 180) =    0.57
 given_sigma_ion_N2_to_N2_cm2( 181) =    0.96 ;  given_sigma_ion_N2_to_N_cm2( 181) =    0.58
 given_sigma_ion_N2_to_N2_cm2( 182) =    0.97 ;  given_sigma_ion_N2_to_N_cm2( 182) =    0.58
 given_sigma_ion_N2_to_N2_cm2( 183) =    0.99 ;  given_sigma_ion_N2_to_N_cm2( 183) =    0.59
 given_sigma_ion_N2_to_N2_cm2( 184) =    1.00 ;  given_sigma_ion_N2_to_N_cm2( 184) =    0.60
 given_sigma_ion_N2_to_N2_cm2( 185) =    1.01 ;  given_sigma_ion_N2_to_N_cm2( 185) =    0.60
 given_sigma_ion_N2_to_N2_cm2( 186) =    1.02 ;  given_sigma_ion_N2_to_N_cm2( 186) =    0.61
 given_sigma_ion_N2_to_N2_cm2( 187) =    1.03 ;  given_sigma_ion_N2_to_N_cm2( 187) =    0.61
 given_sigma_ion_N2_to_N2_cm2( 188) =    1.06 ;  given_sigma_ion_N2_to_N_cm2( 188) =    0.63
 given_sigma_ion_N2_to_N2_cm2( 189) =    1.08 ;  given_sigma_ion_N2_to_N_cm2( 189) =    0.64
 given_sigma_ion_N2_to_N2_cm2( 190) =    1.08 ;  given_sigma_ion_N2_to_N_cm2( 190) =    0.64
 given_sigma_ion_N2_to_N2_cm2( 191) =    1.09 ;  given_sigma_ion_N2_to_N_cm2( 191) =    0.64
 given_sigma_ion_N2_to_N2_cm2( 192) =    1.11 ;  given_sigma_ion_N2_to_N_cm2( 192) =    0.66
 given_sigma_ion_N2_to_N2_cm2( 193) =    1.12 ;  given_sigma_ion_N2_to_N_cm2( 193) =    0.66
 given_sigma_ion_N2_to_N2_cm2( 194) =    1.15 ;  given_sigma_ion_N2_to_N_cm2( 194) =    0.66
 given_sigma_ion_N2_to_N2_cm2( 195) =    1.17 ;  given_sigma_ion_N2_to_N_cm2( 195) =    0.67
 given_sigma_ion_N2_to_N2_cm2( 196) =    1.21 ;  given_sigma_ion_N2_to_N_cm2( 196) =    0.67
 given_sigma_ion_N2_to_N2_cm2( 197) =    1.28 ;  given_sigma_ion_N2_to_N_cm2( 197) =    0.67
 given_sigma_ion_N2_to_N2_cm2( 198) =    1.29 ;  given_sigma_ion_N2_to_N_cm2( 198) =    0.68
 given_sigma_ion_N2_to_N2_cm2( 199) =    1.31 ;  given_sigma_ion_N2_to_N_cm2( 199) =    0.69
 given_sigma_ion_N2_to_N2_cm2( 200) =    1.33 ;  given_sigma_ion_N2_to_N_cm2( 200) =    0.70
 given_sigma_ion_N2_to_N2_cm2( 201) =    1.35 ;  given_sigma_ion_N2_to_N_cm2( 201) =    0.71
 given_sigma_ion_N2_to_N2_cm2( 202) =    1.37 ;  given_sigma_ion_N2_to_N_cm2( 202) =    0.72
 given_sigma_ion_N2_to_N2_cm2( 203) =    1.38 ;  given_sigma_ion_N2_to_N_cm2( 203) =    0.73
 given_sigma_ion_N2_to_N2_cm2( 204) =    1.41 ;  given_sigma_ion_N2_to_N_cm2( 204) =    0.74
 given_sigma_ion_N2_to_N2_cm2( 205) =    1.48 ;  given_sigma_ion_N2_to_N_cm2( 205) =    0.79
 given_sigma_ion_N2_to_N2_cm2( 206) =    1.54 ;  given_sigma_ion_N2_to_N_cm2( 206) =    0.83
 given_sigma_ion_N2_to_N2_cm2( 207) =    1.54 ;  given_sigma_ion_N2_to_N_cm2( 207) =    0.83
 given_sigma_ion_N2_to_N2_cm2( 208) =    1.55 ;  given_sigma_ion_N2_to_N_cm2( 208) =    0.83
 given_sigma_ion_N2_to_N2_cm2( 209) =    1.55 ;  given_sigma_ion_N2_to_N_cm2( 209) =    0.83
 given_sigma_ion_N2_to_N2_cm2( 210) =    1.57 ;  given_sigma_ion_N2_to_N_cm2( 210) =    0.84
 given_sigma_ion_N2_to_N2_cm2( 211) =    1.58 ;  given_sigma_ion_N2_to_N_cm2( 211) =    0.84
 given_sigma_ion_N2_to_N2_cm2( 212) =    1.69 ;  given_sigma_ion_N2_to_N_cm2( 212) =    0.87
 given_sigma_ion_N2_to_N2_cm2( 213) =    1.73 ;  given_sigma_ion_N2_to_N_cm2( 213) =    0.88
 given_sigma_ion_N2_to_N2_cm2( 214) =    1.73 ;  given_sigma_ion_N2_to_N_cm2( 214) =    0.89
 given_sigma_ion_N2_to_N2_cm2( 215) =    1.77 ;  given_sigma_ion_N2_to_N_cm2( 215) =    0.90
 given_sigma_ion_N2_to_N2_cm2( 216) =    1.84 ;  given_sigma_ion_N2_to_N_cm2( 216) =    0.93
 given_sigma_ion_N2_to_N2_cm2( 217) =    1.88 ;  given_sigma_ion_N2_to_N_cm2( 217) =    0.94
 given_sigma_ion_N2_to_N2_cm2( 218) =    1.97 ;  given_sigma_ion_N2_to_N_cm2( 218) =    0.97
 given_sigma_ion_N2_to_N2_cm2( 219) =    2.00 ;  given_sigma_ion_N2_to_N_cm2( 219) =    0.98
 given_sigma_ion_N2_to_N2_cm2( 220) =    2.03 ;  given_sigma_ion_N2_to_N_cm2( 220) =    1.00
 given_sigma_ion_N2_to_N2_cm2( 221) =    2.09 ;  given_sigma_ion_N2_to_N_cm2( 221) =    1.04
 given_sigma_ion_N2_to_N2_cm2( 222) =    2.14 ;  given_sigma_ion_N2_to_N_cm2( 222) =    1.07
 given_sigma_ion_N2_to_N2_cm2( 223) =    2.14 ;  given_sigma_ion_N2_to_N_cm2( 223) =    1.07
 given_sigma_ion_N2_to_N2_cm2( 224) =    2.20 ;  given_sigma_ion_N2_to_N_cm2( 224) =    1.12
 given_sigma_ion_N2_to_N2_cm2( 225) =    2.22 ;  given_sigma_ion_N2_to_N_cm2( 225) =    1.13
 given_sigma_ion_N2_to_N2_cm2( 226) =    2.26 ;  given_sigma_ion_N2_to_N_cm2( 226) =    1.16
 given_sigma_ion_N2_to_N2_cm2( 227) =    2.28 ;  given_sigma_ion_N2_to_N_cm2( 227) =    1.18
 given_sigma_ion_N2_to_N2_cm2( 228) =    2.35 ;  given_sigma_ion_N2_to_N_cm2( 228) =    1.24
 given_sigma_ion_N2_to_N2_cm2( 229) =    2.37 ;  given_sigma_ion_N2_to_N_cm2( 229) =    1.26
 given_sigma_ion_N2_to_N2_cm2( 230) =    2.41 ;  given_sigma_ion_N2_to_N_cm2( 230) =    1.30
 given_sigma_ion_N2_to_N2_cm2( 231) =    2.42 ;  given_sigma_ion_N2_to_N_cm2( 231) =    1.31
 given_sigma_ion_N2_to_N2_cm2( 232) =    2.54 ;  given_sigma_ion_N2_to_N_cm2( 232) =    1.38
 given_sigma_ion_N2_to_N2_cm2( 233) =    2.57 ;  given_sigma_ion_N2_to_N_cm2( 233) =    1.40
 given_sigma_ion_N2_to_N2_cm2( 234) =    2.58 ;  given_sigma_ion_N2_to_N_cm2( 234) =    1.41
 given_sigma_ion_N2_to_N2_cm2( 235) =    2.66 ;  given_sigma_ion_N2_to_N_cm2( 235) =    1.45
 given_sigma_ion_N2_to_N2_cm2( 236) =    2.68 ;  given_sigma_ion_N2_to_N_cm2( 236) =    1.46
 given_sigma_ion_N2_to_N2_cm2( 237) =    2.69 ;  given_sigma_ion_N2_to_N_cm2( 237) =    1.47
 given_sigma_ion_N2_to_N2_cm2( 238) =    2.70 ;  given_sigma_ion_N2_to_N_cm2( 238) =    1.48
 given_sigma_ion_N2_to_N2_cm2( 239) =    2.73 ;  given_sigma_ion_N2_to_N_cm2( 239) =    1.49
 given_sigma_ion_N2_to_N2_cm2( 240) =    2.74 ;  given_sigma_ion_N2_to_N_cm2( 240) =    1.50
 given_sigma_ion_N2_to_N2_cm2( 241) =    2.78 ;  given_sigma_ion_N2_to_N_cm2( 241) =    1.52
 given_sigma_ion_N2_to_N2_cm2( 242) =    2.81 ;  given_sigma_ion_N2_to_N_cm2( 242) =    1.54
 given_sigma_ion_N2_to_N2_cm2( 243) =    2.84 ;  given_sigma_ion_N2_to_N_cm2( 243) =    1.56
 given_sigma_ion_N2_to_N2_cm2( 244) =    2.84 ;  given_sigma_ion_N2_to_N_cm2( 244) =    1.56
 given_sigma_ion_N2_to_N2_cm2( 245) =    2.89 ;  given_sigma_ion_N2_to_N_cm2( 245) =    1.59
 given_sigma_ion_N2_to_N2_cm2( 246) =    2.91 ;  given_sigma_ion_N2_to_N_cm2( 246) =    1.60
 given_sigma_ion_N2_to_N2_cm2( 247) =    2.92 ;  given_sigma_ion_N2_to_N_cm2( 247) =    1.61
 given_sigma_ion_N2_to_N2_cm2( 248) =    2.93 ;  given_sigma_ion_N2_to_N_cm2( 248) =    1.61
 given_sigma_ion_N2_to_N2_cm2( 249) =    2.98 ;  given_sigma_ion_N2_to_N_cm2( 249) =    1.65
 given_sigma_ion_N2_to_N2_cm2( 250) =    2.99 ;  given_sigma_ion_N2_to_N_cm2( 250) =    1.66
 given_sigma_ion_N2_to_N2_cm2( 251) =    3.01 ;  given_sigma_ion_N2_to_N_cm2( 251) =    1.68
 given_sigma_ion_N2_to_N2_cm2( 252) =    3.06 ;  given_sigma_ion_N2_to_N_cm2( 252) =    1.71
 given_sigma_ion_N2_to_N2_cm2( 253) =    3.07 ;  given_sigma_ion_N2_to_N_cm2( 253) =    1.72
 given_sigma_ion_N2_to_N2_cm2( 254) =    3.08 ;  given_sigma_ion_N2_to_N_cm2( 254) =    1.73
 given_sigma_ion_N2_to_N2_cm2( 255) =    3.09 ;  given_sigma_ion_N2_to_N_cm2( 255) =    1.74
 given_sigma_ion_N2_to_N2_cm2( 256) =    3.11 ;  given_sigma_ion_N2_to_N_cm2( 256) =    1.75
 given_sigma_ion_N2_to_N2_cm2( 257) =    3.12 ;  given_sigma_ion_N2_to_N_cm2( 257) =    1.77
 given_sigma_ion_N2_to_N2_cm2( 258) =    3.16 ;  given_sigma_ion_N2_to_N_cm2( 258) =    1.80
 given_sigma_ion_N2_to_N2_cm2( 259) =    3.17 ;  given_sigma_ion_N2_to_N_cm2( 259) =    1.81
 given_sigma_ion_N2_to_N2_cm2( 260) =    3.20 ;  given_sigma_ion_N2_to_N_cm2( 260) =    1.85
 given_sigma_ion_N2_to_N2_cm2( 261) =    3.22 ;  given_sigma_ion_N2_to_N_cm2( 261) =    1.86
 given_sigma_ion_N2_to_N2_cm2( 262) =    3.23 ;  given_sigma_ion_N2_to_N_cm2( 262) =    1.87
 given_sigma_ion_N2_to_N2_cm2( 263) =    3.24 ;  given_sigma_ion_N2_to_N_cm2( 263) =    1.88
 given_sigma_ion_N2_to_N2_cm2( 264) =    3.25 ;  given_sigma_ion_N2_to_N_cm2( 264) =    1.89
 given_sigma_ion_N2_to_N2_cm2( 265) =    3.26 ;  given_sigma_ion_N2_to_N_cm2( 265) =    1.90
 given_sigma_ion_N2_to_N2_cm2( 266) =    3.27 ;  given_sigma_ion_N2_to_N_cm2( 266) =    1.91
 given_sigma_ion_N2_to_N2_cm2( 267) =    3.32 ;  given_sigma_ion_N2_to_N_cm2( 267) =    1.95
 given_sigma_ion_N2_to_N2_cm2( 268) =    3.33 ;  given_sigma_ion_N2_to_N_cm2( 268) =    1.96
 given_sigma_ion_N2_to_N2_cm2( 269) =    3.37 ;  given_sigma_ion_N2_to_N_cm2( 269) =    1.99
 given_sigma_ion_N2_to_N2_cm2( 270) =    3.38 ;  given_sigma_ion_N2_to_N_cm2( 270) =    2.00
 given_sigma_ion_N2_to_N2_cm2( 271) =    3.40 ;  given_sigma_ion_N2_to_N_cm2( 271) =    2.01
 given_sigma_ion_N2_to_N2_cm2( 272) =    3.45 ;  given_sigma_ion_N2_to_N_cm2( 272) =    2.05
 given_sigma_ion_N2_to_N2_cm2( 273) =    3.48 ;  given_sigma_ion_N2_to_N_cm2( 273) =    2.07
 given_sigma_ion_N2_to_N2_cm2( 274) =    3.49 ;  given_sigma_ion_N2_to_N_cm2( 274) =    2.08
 given_sigma_ion_N2_to_N2_cm2( 275) =    3.50 ;  given_sigma_ion_N2_to_N_cm2( 275) =    2.08
 given_sigma_ion_N2_to_N2_cm2( 276) =    3.54 ;  given_sigma_ion_N2_to_N_cm2( 276) =    2.11
 given_sigma_ion_N2_to_N2_cm2( 277) =    3.56 ;  given_sigma_ion_N2_to_N_cm2( 277) =    2.12
 given_sigma_ion_N2_to_N2_cm2( 278) =    3.59 ;  given_sigma_ion_N2_to_N_cm2( 278) =    2.14
 given_sigma_ion_N2_to_N2_cm2( 279) =    3.65 ;  given_sigma_ion_N2_to_N_cm2( 279) =    2.18
 given_sigma_ion_N2_to_N2_cm2( 280) =    3.65 ;  given_sigma_ion_N2_to_N_cm2( 280) =    2.19
 given_sigma_ion_N2_to_N2_cm2( 281) =    3.72 ;  given_sigma_ion_N2_to_N_cm2( 281) =    2.23
 given_sigma_ion_N2_to_N2_cm2( 282) =    3.72 ;  given_sigma_ion_N2_to_N_cm2( 282) =    2.23
 given_sigma_ion_N2_to_N2_cm2( 283) =    3.76 ;  given_sigma_ion_N2_to_N_cm2( 283) =    2.25
 given_sigma_ion_N2_to_N2_cm2( 284) =    3.80 ;  given_sigma_ion_N2_to_N_cm2( 284) =    2.29
 given_sigma_ion_N2_to_N2_cm2( 285) =    3.87 ;  given_sigma_ion_N2_to_N_cm2( 285) =    2.33
 given_sigma_ion_N2_to_N2_cm2( 286) =    3.92 ;  given_sigma_ion_N2_to_N_cm2( 286) =    2.36
 given_sigma_ion_N2_to_N2_cm2( 287) =    3.96 ;  given_sigma_ion_N2_to_N_cm2( 287) =    2.39
 given_sigma_ion_N2_to_N2_cm2( 288) =    3.99 ;  given_sigma_ion_N2_to_N_cm2( 288) =    2.41
 given_sigma_ion_N2_to_N2_cm2( 289) =    4.04 ;  given_sigma_ion_N2_to_N_cm2( 289) =    2.44
 given_sigma_ion_N2_to_N2_cm2( 290) =    4.07 ;  given_sigma_ion_N2_to_N_cm2( 290) =    2.46
 given_sigma_ion_N2_to_N2_cm2( 291) =    4.09 ;  given_sigma_ion_N2_to_N_cm2( 291) =    2.48
 given_sigma_ion_N2_to_N2_cm2( 292) =    4.16 ;  given_sigma_ion_N2_to_N_cm2( 292) =    2.52
 given_sigma_ion_N2_to_N2_cm2( 293) =    4.17 ;  given_sigma_ion_N2_to_N_cm2( 293) =    2.52
 given_sigma_ion_N2_to_N2_cm2( 294) =    4.18 ;  given_sigma_ion_N2_to_N_cm2( 294) =    2.53
 given_sigma_ion_N2_to_N2_cm2( 295) =    4.22 ;  given_sigma_ion_N2_to_N_cm2( 295) =    2.56
 given_sigma_ion_N2_to_N2_cm2( 296) =    4.27 ;  given_sigma_ion_N2_to_N_cm2( 296) =    2.59
 given_sigma_ion_N2_to_N2_cm2( 297) =    4.33 ;  given_sigma_ion_N2_to_N_cm2( 297) =    2.63
 given_sigma_ion_N2_to_N2_cm2( 298) =    4.34 ;  given_sigma_ion_N2_to_N_cm2( 298) =    2.63
 given_sigma_ion_N2_to_N2_cm2( 299) =    4.35 ;  given_sigma_ion_N2_to_N_cm2( 299) =    2.64
 given_sigma_ion_N2_to_N2_cm2( 300) =    4.41 ;  given_sigma_ion_N2_to_N_cm2( 300) =    2.69
 given_sigma_ion_N2_to_N2_cm2( 301) =    4.45 ;  given_sigma_ion_N2_to_N_cm2( 301) =    2.72
 given_sigma_ion_N2_to_N2_cm2( 302) =    4.52 ;  given_sigma_ion_N2_to_N_cm2( 302) =    2.78
 given_sigma_ion_N2_to_N2_cm2( 303) =    4.57 ;  given_sigma_ion_N2_to_N_cm2( 303) =    2.81
 given_sigma_ion_N2_to_N2_cm2( 304) =    4.58 ;  given_sigma_ion_N2_to_N_cm2( 304) =    2.82
 given_sigma_ion_N2_to_N2_cm2( 305) =    4.59 ;  given_sigma_ion_N2_to_N_cm2( 305) =    2.83
 given_sigma_ion_N2_to_N2_cm2( 306) =    4.68 ;  given_sigma_ion_N2_to_N_cm2( 306) =    2.89
 given_sigma_ion_N2_to_N2_cm2( 307) =    4.68 ;  given_sigma_ion_N2_to_N_cm2( 307) =    2.90
 given_sigma_ion_N2_to_N2_cm2( 308) =    4.75 ;  given_sigma_ion_N2_to_N_cm2( 308) =    2.94
 given_sigma_ion_N2_to_N2_cm2( 309) =    4.79 ;  given_sigma_ion_N2_to_N_cm2( 309) =    2.98
 given_sigma_ion_N2_to_N2_cm2( 310) =    4.84 ;  given_sigma_ion_N2_to_N_cm2( 310) =    3.01
 given_sigma_ion_N2_to_N2_cm2( 311) =    4.91 ;  given_sigma_ion_N2_to_N_cm2( 311) =    3.05
 given_sigma_ion_N2_to_N2_cm2( 312) =    4.91 ;  given_sigma_ion_N2_to_N_cm2( 312) =    3.06
 given_sigma_ion_N2_to_N2_cm2( 313) =    4.93 ;  given_sigma_ion_N2_to_N_cm2( 313) =    3.07
 given_sigma_ion_N2_to_N2_cm2( 314) =    5.01 ;  given_sigma_ion_N2_to_N_cm2( 314) =    3.12
 given_sigma_ion_N2_to_N2_cm2( 315) =    5.03 ;  given_sigma_ion_N2_to_N_cm2( 315) =    3.14
 given_sigma_ion_N2_to_N2_cm2( 316) =    5.09 ;  given_sigma_ion_N2_to_N_cm2( 316) =    3.17
 given_sigma_ion_N2_to_N2_cm2( 317) =    5.10 ;  given_sigma_ion_N2_to_N_cm2( 317) =    3.18
 given_sigma_ion_N2_to_N2_cm2( 318) =    5.11 ;  given_sigma_ion_N2_to_N_cm2( 318) =    3.18
 given_sigma_ion_N2_to_N2_cm2( 319) =    5.20 ;  given_sigma_ion_N2_to_N_cm2( 319) =    3.25
 given_sigma_ion_N2_to_N2_cm2( 320) =    5.21 ;  given_sigma_ion_N2_to_N_cm2( 320) =    3.25
 given_sigma_ion_N2_to_N2_cm2( 321) =    5.23 ;  given_sigma_ion_N2_to_N_cm2( 321) =    3.26
 given_sigma_ion_N2_to_N2_cm2( 322) =    5.29 ;  given_sigma_ion_N2_to_N_cm2( 322) =    3.31
 given_sigma_ion_N2_to_N2_cm2( 323) =    5.34 ;  given_sigma_ion_N2_to_N_cm2( 323) =    3.34
 given_sigma_ion_N2_to_N2_cm2( 324) =    5.36 ;  given_sigma_ion_N2_to_N_cm2( 324) =    3.35
 given_sigma_ion_N2_to_N2_cm2( 325) =    5.40 ;  given_sigma_ion_N2_to_N_cm2( 325) =    3.37
 given_sigma_ion_N2_to_N2_cm2( 326) =    5.45 ;  given_sigma_ion_N2_to_N_cm2( 326) =    3.39
 given_sigma_ion_N2_to_N2_cm2( 327) =    5.51 ;  given_sigma_ion_N2_to_N_cm2( 327) =    3.42
 given_sigma_ion_N2_to_N2_cm2( 328) =    5.57 ;  given_sigma_ion_N2_to_N_cm2( 328) =    3.45
 given_sigma_ion_N2_to_N2_cm2( 329) =    5.61 ;  given_sigma_ion_N2_to_N_cm2( 329) =    3.46
 given_sigma_ion_N2_to_N2_cm2( 330) =    5.64 ;  given_sigma_ion_N2_to_N_cm2( 330) =    3.48
 given_sigma_ion_N2_to_N2_cm2( 331) =    5.67 ;  given_sigma_ion_N2_to_N_cm2( 331) =    3.49
 given_sigma_ion_N2_to_N2_cm2( 332) =    5.77 ;  given_sigma_ion_N2_to_N_cm2( 332) =    3.53
 given_sigma_ion_N2_to_N2_cm2( 333) =    5.82 ;  given_sigma_ion_N2_to_N_cm2( 333) =    3.56
 given_sigma_ion_N2_to_N2_cm2( 334) =    5.90 ;  given_sigma_ion_N2_to_N_cm2( 334) =    3.59
 given_sigma_ion_N2_to_N2_cm2( 335) =    5.91 ;  given_sigma_ion_N2_to_N_cm2( 335) =    3.59
 given_sigma_ion_N2_to_N2_cm2( 336) =    5.95 ;  given_sigma_ion_N2_to_N_cm2( 336) =    3.60
 given_sigma_ion_N2_to_N2_cm2( 337) =    6.00 ;  given_sigma_ion_N2_to_N_cm2( 337) =    3.60
 given_sigma_ion_N2_to_N2_cm2( 338) =    6.07 ;  given_sigma_ion_N2_to_N_cm2( 338) =    3.61
 given_sigma_ion_N2_to_N2_cm2( 339) =    6.11 ;  given_sigma_ion_N2_to_N_cm2( 339) =    3.62
 given_sigma_ion_N2_to_N2_cm2( 340) =    6.17 ;  given_sigma_ion_N2_to_N_cm2( 340) =    3.63
 given_sigma_ion_N2_to_N2_cm2( 341) =    6.22 ;  given_sigma_ion_N2_to_N_cm2( 341) =    3.62
 given_sigma_ion_N2_to_N2_cm2( 342) =    6.24 ;  given_sigma_ion_N2_to_N_cm2( 342) =    3.62
 given_sigma_ion_N2_to_N2_cm2( 343) =    6.28 ;  given_sigma_ion_N2_to_N_cm2( 343) =    3.61
 given_sigma_ion_N2_to_N2_cm2( 344) =    6.29 ;  given_sigma_ion_N2_to_N_cm2( 344) =    3.60
 given_sigma_ion_N2_to_N2_cm2( 345) =    6.34 ;  given_sigma_ion_N2_to_N_cm2( 345) =    3.59
 given_sigma_ion_N2_to_N2_cm2( 346) =    6.40 ;  given_sigma_ion_N2_to_N_cm2( 346) =    3.58
 given_sigma_ion_N2_to_N2_cm2( 347) =    6.45 ;  given_sigma_ion_N2_to_N_cm2( 347) =    3.57
 given_sigma_ion_N2_to_N2_cm2( 348) =    6.51 ;  given_sigma_ion_N2_to_N_cm2( 348) =    3.55
 given_sigma_ion_N2_to_N2_cm2( 349) =    6.55 ;  given_sigma_ion_N2_to_N_cm2( 349) =    3.53
 given_sigma_ion_N2_to_N2_cm2( 350) =    6.56 ;  given_sigma_ion_N2_to_N_cm2( 350) =    3.53
 given_sigma_ion_N2_to_N2_cm2( 351) =    6.65 ;  given_sigma_ion_N2_to_N_cm2( 351) =    3.49
 given_sigma_ion_N2_to_N2_cm2( 352) =    6.68 ;  given_sigma_ion_N2_to_N_cm2( 352) =    3.48
 given_sigma_ion_N2_to_N2_cm2( 353) =    6.71 ;  given_sigma_ion_N2_to_N_cm2( 353) =    3.47
 given_sigma_ion_N2_to_N2_cm2( 354) =    6.78 ;  given_sigma_ion_N2_to_N_cm2( 354) =    3.42
 given_sigma_ion_N2_to_N2_cm2( 355) =    6.80 ;  given_sigma_ion_N2_to_N_cm2( 355) =    3.41
 given_sigma_ion_N2_to_N2_cm2( 356) =    6.83 ;  given_sigma_ion_N2_to_N_cm2( 356) =    3.39
 given_sigma_ion_N2_to_N2_cm2( 357) =    6.84 ;  given_sigma_ion_N2_to_N_cm2( 357) =    3.38
 given_sigma_ion_N2_to_N2_cm2( 358) =    6.89 ;  given_sigma_ion_N2_to_N_cm2( 358) =    3.35
 given_sigma_ion_N2_to_N2_cm2( 359) =    6.96 ;  given_sigma_ion_N2_to_N_cm2( 359) =    3.31
 given_sigma_ion_N2_to_N2_cm2( 360) =    6.99 ;  given_sigma_ion_N2_to_N_cm2( 360) =    3.29
 given_sigma_ion_N2_to_N2_cm2( 361) =    7.03 ;  given_sigma_ion_N2_to_N_cm2( 361) =    3.27
 given_sigma_ion_N2_to_N2_cm2( 362) =    7.12 ;  given_sigma_ion_N2_to_N_cm2( 362) =    3.24
 given_sigma_ion_N2_to_N2_cm2( 363) =    7.17 ;  given_sigma_ion_N2_to_N_cm2( 363) =    3.21
 given_sigma_ion_N2_to_N2_cm2( 364) =    7.19 ;  given_sigma_ion_N2_to_N_cm2( 364) =    3.20
 given_sigma_ion_N2_to_N2_cm2( 365) =    7.20 ;  given_sigma_ion_N2_to_N_cm2( 365) =    3.20
 given_sigma_ion_N2_to_N2_cm2( 366) =    7.42 ;  given_sigma_ion_N2_to_N_cm2( 366) =    3.07
 given_sigma_ion_N2_to_N2_cm2( 367) =    7.44 ;  given_sigma_ion_N2_to_N_cm2( 367) =    3.06
 given_sigma_ion_N2_to_N2_cm2( 368) =    7.46 ;  given_sigma_ion_N2_to_N_cm2( 368) =    3.05
 given_sigma_ion_N2_to_N2_cm2( 369) =    7.54 ;  given_sigma_ion_N2_to_N_cm2( 369) =    3.01
 given_sigma_ion_N2_to_N2_cm2( 370) =    7.57 ;  given_sigma_ion_N2_to_N_cm2( 370) =    2.99
 given_sigma_ion_N2_to_N2_cm2( 371) =    7.64 ;  given_sigma_ion_N2_to_N_cm2( 371) =    2.95
 given_sigma_ion_N2_to_N2_cm2( 372) =    7.68 ;  given_sigma_ion_N2_to_N_cm2( 372) =    2.93
 given_sigma_ion_N2_to_N2_cm2( 373) =    7.69 ;  given_sigma_ion_N2_to_N_cm2( 373) =    2.93
 given_sigma_ion_N2_to_N2_cm2( 374) =    7.71 ;  given_sigma_ion_N2_to_N_cm2( 374) =    2.93
 given_sigma_ion_N2_to_N2_cm2( 375) =    7.73 ;  given_sigma_ion_N2_to_N_cm2( 375) =    2.92
 given_sigma_ion_N2_to_N2_cm2( 376) =    7.75 ;  given_sigma_ion_N2_to_N_cm2( 376) =    2.92
 given_sigma_ion_N2_to_N2_cm2( 377) =    7.76 ;  given_sigma_ion_N2_to_N_cm2( 377) =    2.92
 given_sigma_ion_N2_to_N2_cm2( 378) =    7.77 ;  given_sigma_ion_N2_to_N_cm2( 378) =    2.92
 given_sigma_ion_N2_to_N2_cm2( 379) =    7.82 ;  given_sigma_ion_N2_to_N_cm2( 379) =    2.91
 given_sigma_ion_N2_to_N2_cm2( 380) =    7.88 ;  given_sigma_ion_N2_to_N_cm2( 380) =    2.90
 given_sigma_ion_N2_to_N2_cm2( 381) =    7.93 ;  given_sigma_ion_N2_to_N_cm2( 381) =    2.88
 given_sigma_ion_N2_to_N2_cm2( 382) =    7.95 ;  given_sigma_ion_N2_to_N_cm2( 382) =    2.87
 given_sigma_ion_N2_to_N2_cm2( 383) =    8.10 ;  given_sigma_ion_N2_to_N_cm2( 383) =    2.80
 given_sigma_ion_N2_to_N2_cm2( 384) =    8.14 ;  given_sigma_ion_N2_to_N_cm2( 384) =    2.78
 given_sigma_ion_N2_to_N2_cm2( 385) =    8.18 ;  given_sigma_ion_N2_to_N_cm2( 385) =    2.77
 given_sigma_ion_N2_to_N2_cm2( 386) =    8.18 ;  given_sigma_ion_N2_to_N_cm2( 386) =    2.77
 given_sigma_ion_N2_to_N2_cm2( 387) =    8.31 ;  given_sigma_ion_N2_to_N_cm2( 387) =    2.73
 given_sigma_ion_N2_to_N2_cm2( 388) =    8.36 ;  given_sigma_ion_N2_to_N_cm2( 388) =    2.71
 given_sigma_ion_N2_to_N2_cm2( 389) =    8.40 ;  given_sigma_ion_N2_to_N_cm2( 389) =    2.70
 given_sigma_ion_N2_to_N2_cm2( 390) =    8.43 ;  given_sigma_ion_N2_to_N_cm2( 390) =    2.69
 given_sigma_ion_N2_to_N2_cm2( 391) =    8.48 ;  given_sigma_ion_N2_to_N_cm2( 391) =    2.67
 given_sigma_ion_N2_to_N2_cm2( 392) =    8.50 ;  given_sigma_ion_N2_to_N_cm2( 392) =    2.67
 given_sigma_ion_N2_to_N2_cm2( 393) =    8.54 ;  given_sigma_ion_N2_to_N_cm2( 393) =    2.66
 given_sigma_ion_N2_to_N2_cm2( 394) =    8.65 ;  given_sigma_ion_N2_to_N_cm2( 394) =    2.62
 given_sigma_ion_N2_to_N2_cm2( 395) =    8.68 ;  given_sigma_ion_N2_to_N_cm2( 395) =    2.62
 given_sigma_ion_N2_to_N2_cm2( 396) =    8.71 ;  given_sigma_ion_N2_to_N_cm2( 396) =    2.61
 given_sigma_ion_N2_to_N2_cm2( 397) =    8.88 ;  given_sigma_ion_N2_to_N_cm2( 397) =    2.59
 given_sigma_ion_N2_to_N2_cm2( 398) =    8.91 ;  given_sigma_ion_N2_to_N_cm2( 398) =    2.59
 given_sigma_ion_N2_to_N2_cm2( 399) =    9.04 ;  given_sigma_ion_N2_to_N_cm2( 399) =    2.51
 given_sigma_ion_N2_to_N2_cm2( 400) =    9.18 ;  given_sigma_ion_N2_to_N_cm2( 400) =    2.49
 given_sigma_ion_N2_to_N2_cm2( 401) =    9.21 ;  given_sigma_ion_N2_to_N_cm2( 401) =    2.49
 given_sigma_ion_N2_to_N2_cm2( 402) =    9.57 ;  given_sigma_ion_N2_to_N_cm2( 402) =    2.45
 given_sigma_ion_N2_to_N2_cm2( 403) =    9.92 ;  given_sigma_ion_N2_to_N_cm2( 403) =    2.45
 given_sigma_ion_N2_to_N2_cm2( 404) =   10.00 ;  given_sigma_ion_N2_to_N_cm2( 404) =    2.45
 given_sigma_ion_N2_to_N2_cm2( 405) =   10.01 ;  given_sigma_ion_N2_to_N_cm2( 405) =    2.45
 given_sigma_ion_N2_to_N2_cm2( 406) =   10.22 ;  given_sigma_ion_N2_to_N_cm2( 406) =    2.44
 given_sigma_ion_N2_to_N2_cm2( 407) =   10.29 ;  given_sigma_ion_N2_to_N_cm2( 407) =    2.43
 given_sigma_ion_N2_to_N2_cm2( 408) =   10.30 ;  given_sigma_ion_N2_to_N_cm2( 408) =    2.43
 given_sigma_ion_N2_to_N2_cm2( 409) =   10.34 ;  given_sigma_ion_N2_to_N_cm2( 409) =    2.43
 given_sigma_ion_N2_to_N2_cm2( 410) =   10.68 ;  given_sigma_ion_N2_to_N_cm2( 410) =    2.45
 given_sigma_ion_N2_to_N2_cm2( 411) =   11.11 ;  given_sigma_ion_N2_to_N_cm2( 411) =    2.44
 given_sigma_ion_N2_to_N2_cm2( 412) =   11.60 ;  given_sigma_ion_N2_to_N_cm2( 412) =    2.35
 given_sigma_ion_N2_to_N2_cm2( 413) =   11.65 ;  given_sigma_ion_N2_to_N_cm2( 413) =    2.34
 given_sigma_ion_N2_to_N2_cm2( 414) =   12.20 ;  given_sigma_ion_N2_to_N_cm2( 414) =    2.18
 given_sigma_ion_N2_to_N2_cm2( 415) =   12.70 ;  given_sigma_ion_N2_to_N_cm2( 415) =    2.12
 given_sigma_ion_N2_to_N2_cm2( 416) =   12.72 ;  given_sigma_ion_N2_to_N_cm2( 416) =    2.12
 given_sigma_ion_N2_to_N2_cm2( 417) =   12.79 ;  given_sigma_ion_N2_to_N_cm2( 417) =    2.10
 given_sigma_ion_N2_to_N2_cm2( 418) =   12.99 ;  given_sigma_ion_N2_to_N_cm2( 418) =    2.05
 given_sigma_ion_N2_to_N2_cm2( 419) =   13.28 ;  given_sigma_ion_N2_to_N_cm2( 419) =    1.98
 given_sigma_ion_N2_to_N2_cm2( 420) =   13.30 ;  given_sigma_ion_N2_to_N_cm2( 420) =    1.98
 given_sigma_ion_N2_to_N2_cm2( 421) =   13.72 ;  given_sigma_ion_N2_to_N_cm2( 421) =    1.91
 given_sigma_ion_N2_to_N2_cm2( 422) =   13.96 ;  given_sigma_ion_N2_to_N_cm2( 422) =    1.86
 given_sigma_ion_N2_to_N2_cm2( 423) =   14.40 ;  given_sigma_ion_N2_to_N_cm2( 423) =    1.78
 given_sigma_ion_N2_to_N2_cm2( 424) =   14.49 ;  given_sigma_ion_N2_to_N_cm2( 424) =    1.76
 given_sigma_ion_N2_to_N2_cm2( 425) =   14.92 ;  given_sigma_ion_N2_to_N_cm2( 425) =    1.66
 given_sigma_ion_N2_to_N2_cm2( 426) =   14.96 ;  given_sigma_ion_N2_to_N_cm2( 426) =    1.65
 given_sigma_ion_N2_to_N2_cm2( 427) =   15.35 ;  given_sigma_ion_N2_to_N_cm2( 427) =    1.56
 given_sigma_ion_N2_to_N2_cm2( 428) =   15.58 ;  given_sigma_ion_N2_to_N_cm2( 428) =    1.50
 given_sigma_ion_N2_to_N2_cm2( 429) =   16.15 ;  given_sigma_ion_N2_to_N_cm2( 429) =    1.38
 given_sigma_ion_N2_to_N2_cm2( 430) =   16.80 ;  given_sigma_ion_N2_to_N_cm2( 430) =    1.23
 given_sigma_ion_N2_to_N2_cm2( 431) =   17.90 ;  given_sigma_ion_N2_to_N_cm2( 431) =    1.15
 given_sigma_ion_N2_to_N2_cm2( 432) =   18.98 ;  given_sigma_ion_N2_to_N_cm2( 432) =    1.07
 given_sigma_ion_N2_to_N2_cm2( 433) =   19.00 ;  given_sigma_ion_N2_to_N_cm2( 433) =    1.07
 given_sigma_ion_N2_to_N2_cm2( 434) =   19.11 ;  given_sigma_ion_N2_to_N_cm2( 434) =    1.06
 given_sigma_ion_N2_to_N2_cm2( 435) =   19.17 ;  given_sigma_ion_N2_to_N_cm2( 435) =    1.06
 given_sigma_ion_N2_to_N2_cm2( 436) =   19.19 ;  given_sigma_ion_N2_to_N_cm2( 436) =    1.06
 given_sigma_ion_N2_to_N2_cm2( 437) =   19.33 ;  given_sigma_ion_N2_to_N_cm2( 437) =    1.05
 given_sigma_ion_N2_to_N2_cm2( 438) =   19.50 ;  given_sigma_ion_N2_to_N_cm2( 438) =    1.05
 given_sigma_ion_N2_to_N2_cm2( 439) =   19.60 ;  given_sigma_ion_N2_to_N_cm2( 439) =    1.04
 given_sigma_ion_N2_to_N2_cm2( 440) =   19.70 ;  given_sigma_ion_N2_to_N_cm2( 440) =    1.04
 given_sigma_ion_N2_to_N2_cm2( 441) =   19.80 ;  given_sigma_ion_N2_to_N_cm2( 441) =    1.03
 given_sigma_ion_N2_to_N2_cm2( 442) =   19.90 ;  given_sigma_ion_N2_to_N_cm2( 442) =    1.03
 given_sigma_ion_N2_to_N2_cm2( 443) =   20.00 ;  given_sigma_ion_N2_to_N_cm2( 443) =    1.02
 given_sigma_ion_N2_to_N2_cm2( 444) =   20.10 ;  given_sigma_ion_N2_to_N_cm2( 444) =    1.00
 given_sigma_ion_N2_to_N2_cm2( 445) =   20.20 ;  given_sigma_ion_N2_to_N_cm2( 445) =    0.99
 given_sigma_ion_N2_to_N2_cm2( 446) =   20.30 ;  given_sigma_ion_N2_to_N_cm2( 446) =    0.97
 given_sigma_ion_N2_to_N2_cm2( 447) =   20.40 ;  given_sigma_ion_N2_to_N_cm2( 447) =    0.95
 given_sigma_ion_N2_to_N2_cm2( 448) =   20.50 ;  given_sigma_ion_N2_to_N_cm2( 448) =    0.94
 given_sigma_ion_N2_to_N2_cm2( 449) =   20.60 ;  given_sigma_ion_N2_to_N_cm2( 449) =    0.92
 given_sigma_ion_N2_to_N2_cm2( 450) =   20.70 ;  given_sigma_ion_N2_to_N_cm2( 450) =    0.90
 given_sigma_ion_N2_to_N2_cm2( 451) =   20.72 ;  given_sigma_ion_N2_to_N_cm2( 451) =    0.90
 given_sigma_ion_N2_to_N2_cm2( 452) =   20.77 ;  given_sigma_ion_N2_to_N_cm2( 452) =    0.89
 given_sigma_ion_N2_to_N2_cm2( 453) =   20.80 ;  given_sigma_ion_N2_to_N_cm2( 453) =    0.89
 given_sigma_ion_N2_to_N2_cm2( 454) =   20.90 ;  given_sigma_ion_N2_to_N_cm2( 454) =    0.87
 given_sigma_ion_N2_to_N2_cm2( 455) =   21.00 ;  given_sigma_ion_N2_to_N_cm2( 455) =    0.85
 given_sigma_ion_N2_to_N2_cm2( 456) =   21.07 ;  given_sigma_ion_N2_to_N_cm2( 456) =    0.86
 given_sigma_ion_N2_to_N2_cm2( 457) =   21.14 ;  given_sigma_ion_N2_to_N_cm2( 457) =    0.86
 given_sigma_ion_N2_to_N2_cm2( 458) =   21.22 ;  given_sigma_ion_N2_to_N_cm2( 458) =    0.87
 given_sigma_ion_N2_to_N2_cm2( 459) =   21.29 ;  given_sigma_ion_N2_to_N_cm2( 459) =    0.88
 given_sigma_ion_N2_to_N2_cm2( 460) =   21.36 ;  given_sigma_ion_N2_to_N_cm2( 460) =    0.88
 given_sigma_ion_N2_to_N2_cm2( 461) =   21.43 ;  given_sigma_ion_N2_to_N_cm2( 461) =    0.89
 given_sigma_ion_N2_to_N2_cm2( 462) =   21.50 ;  given_sigma_ion_N2_to_N_cm2( 462) =    0.90
 given_sigma_ion_N2_to_N2_cm2( 463) =   21.58 ;  given_sigma_ion_N2_to_N_cm2( 463) =    0.91
 given_sigma_ion_N2_to_N2_cm2( 464) =   21.65 ;  given_sigma_ion_N2_to_N_cm2( 464) =    0.91
 given_sigma_ion_N2_to_N2_cm2( 465) =   21.72 ;  given_sigma_ion_N2_to_N_cm2( 465) =    0.92
 given_sigma_ion_N2_to_N2_cm2( 466) =   21.74 ;  given_sigma_ion_N2_to_N_cm2( 466) =    0.92
 given_sigma_ion_N2_to_N2_cm2( 467) =   21.76 ;  given_sigma_ion_N2_to_N_cm2( 467) =    0.92
 given_sigma_ion_N2_to_N2_cm2( 468) =   21.81 ;  given_sigma_ion_N2_to_N_cm2( 468) =    0.93
 given_sigma_ion_N2_to_N2_cm2( 469) =   21.85 ;  given_sigma_ion_N2_to_N_cm2( 469) =    0.93
 given_sigma_ion_N2_to_N2_cm2( 470) =   21.89 ;  given_sigma_ion_N2_to_N_cm2( 470) =    0.93
 given_sigma_ion_N2_to_N2_cm2( 471) =   21.94 ;  given_sigma_ion_N2_to_N_cm2( 471) =    0.93
 given_sigma_ion_N2_to_N2_cm2( 472) =   21.98 ;  given_sigma_ion_N2_to_N_cm2( 472) =    0.94
 given_sigma_ion_N2_to_N2_cm2( 473) =   21.98 ;  given_sigma_ion_N2_to_N_cm2( 473) =    0.94
 given_sigma_ion_N2_to_N2_cm2( 474) =   22.01 ;  given_sigma_ion_N2_to_N_cm2( 474) =    0.94
 given_sigma_ion_N2_to_N2_cm2( 475) =   22.02 ;  given_sigma_ion_N2_to_N_cm2( 475) =    0.94
 given_sigma_ion_N2_to_N2_cm2( 476) =   22.06 ;  given_sigma_ion_N2_to_N_cm2( 476) =    0.94
 given_sigma_ion_N2_to_N2_cm2( 477) =   22.11 ;  given_sigma_ion_N2_to_N_cm2( 477) =    0.95
 given_sigma_ion_N2_to_N2_cm2( 478) =   22.15 ;  given_sigma_ion_N2_to_N_cm2( 478) =    0.95
 given_sigma_ion_N2_to_N2_cm2( 479) =   22.15 ;  given_sigma_ion_N2_to_N_cm2( 479) =    0.95
 given_sigma_ion_N2_to_N2_cm2( 480) =   22.16 ;  given_sigma_ion_N2_to_N_cm2( 480) =    0.94
 given_sigma_ion_N2_to_N2_cm2( 481) =   22.16 ;  given_sigma_ion_N2_to_N_cm2( 481) =    0.94
 given_sigma_ion_N2_to_N2_cm2( 482) =   22.16 ;  given_sigma_ion_N2_to_N_cm2( 482) =    0.94
 given_sigma_ion_N2_to_N2_cm2( 483) =   22.16 ;  given_sigma_ion_N2_to_N_cm2( 483) =    0.94
 given_sigma_ion_N2_to_N2_cm2( 484) =   22.17 ;  given_sigma_ion_N2_to_N_cm2( 484) =    0.93
 given_sigma_ion_N2_to_N2_cm2( 485) =   22.17 ;  given_sigma_ion_N2_to_N_cm2( 485) =    0.93
 given_sigma_ion_N2_to_N2_cm2( 486) =   22.17 ;  given_sigma_ion_N2_to_N_cm2( 486) =    0.93
 given_sigma_ion_N2_to_N2_cm2( 487) =   22.18 ;  given_sigma_ion_N2_to_N_cm2( 487) =    0.93
 given_sigma_ion_N2_to_N2_cm2( 488) =   22.18 ;  given_sigma_ion_N2_to_N_cm2( 488) =    0.92
 given_sigma_ion_N2_to_N2_cm2( 489) =   22.18 ;  given_sigma_ion_N2_to_N_cm2( 489) =    0.92
 given_sigma_ion_N2_to_N2_cm2( 490) =   22.18 ;  given_sigma_ion_N2_to_N_cm2( 490) =    0.92
 given_sigma_ion_N2_to_N2_cm2( 491) =   22.18 ;  given_sigma_ion_N2_to_N_cm2( 491) =    0.92
 given_sigma_ion_N2_to_N2_cm2( 492) =   22.19 ;  given_sigma_ion_N2_to_N_cm2( 492) =    0.92
 given_sigma_ion_N2_to_N2_cm2( 493) =   22.19 ;  given_sigma_ion_N2_to_N_cm2( 493) =    0.91
 given_sigma_ion_N2_to_N2_cm2( 494) =   22.19 ;  given_sigma_ion_N2_to_N_cm2( 494) =    0.91
 given_sigma_ion_N2_to_N2_cm2( 495) =   22.19 ;  given_sigma_ion_N2_to_N_cm2( 495) =    0.91
 given_sigma_ion_N2_to_N2_cm2( 496) =   22.19 ;  given_sigma_ion_N2_to_N_cm2( 496) =    0.90
 given_sigma_ion_N2_to_N2_cm2( 497) =   22.20 ;  given_sigma_ion_N2_to_N_cm2( 497) =    0.90
 given_sigma_ion_N2_to_N2_cm2( 498) =   22.20 ;  given_sigma_ion_N2_to_N_cm2( 498) =    0.90
 given_sigma_ion_N2_to_N2_cm2( 499) =   22.18 ;  given_sigma_ion_N2_to_N_cm2( 499) =    0.91
 given_sigma_ion_N2_to_N2_cm2( 500) =   22.17 ;  given_sigma_ion_N2_to_N_cm2( 500) =    0.92
 given_sigma_ion_N2_to_N2_cm2( 501) =   22.16 ;  given_sigma_ion_N2_to_N_cm2( 501) =    0.93
 given_sigma_ion_N2_to_N2_cm2( 502) =   22.15 ;  given_sigma_ion_N2_to_N_cm2( 502) =    0.93
 given_sigma_ion_N2_to_N2_cm2( 503) =   22.14 ;  given_sigma_ion_N2_to_N_cm2( 503) =    0.94
 given_sigma_ion_N2_to_N2_cm2( 504) =   22.14 ;  given_sigma_ion_N2_to_N_cm2( 504) =    0.94
 given_sigma_ion_N2_to_N2_cm2( 505) =   22.13 ;  given_sigma_ion_N2_to_N_cm2( 505) =    0.95
 given_sigma_ion_N2_to_N2_cm2( 506) =   22.12 ;  given_sigma_ion_N2_to_N_cm2( 506) =    0.95
 given_sigma_ion_N2_to_N2_cm2( 507) =   22.12 ;  given_sigma_ion_N2_to_N_cm2( 507) =    0.95
 given_sigma_ion_N2_to_N2_cm2( 508) =   22.11 ;  given_sigma_ion_N2_to_N_cm2( 508) =    0.96
 given_sigma_ion_N2_to_N2_cm2( 509) =   22.10 ;  given_sigma_ion_N2_to_N_cm2( 509) =    0.97
 given_sigma_ion_N2_to_N2_cm2( 510) =   22.09 ;  given_sigma_ion_N2_to_N_cm2( 510) =    0.97
 given_sigma_ion_N2_to_N2_cm2( 511) =   22.08 ;  given_sigma_ion_N2_to_N_cm2( 511) =    0.98
 given_sigma_ion_N2_to_N2_cm2( 512) =   22.06 ;  given_sigma_ion_N2_to_N_cm2( 512) =    0.99
 given_sigma_ion_N2_to_N2_cm2( 513) =   22.06 ;  given_sigma_ion_N2_to_N_cm2( 513) =    0.99
 given_sigma_ion_N2_to_N2_cm2( 514) =   22.05 ;  given_sigma_ion_N2_to_N_cm2( 514) =    1.00
 given_sigma_ion_N2_to_N2_cm2( 515) =   22.06 ;  given_sigma_ion_N2_to_N_cm2( 515) =    1.01
 given_sigma_ion_N2_to_N2_cm2( 516) =   22.07 ;  given_sigma_ion_N2_to_N_cm2( 516) =    1.02
 given_sigma_ion_N2_to_N2_cm2( 517) =   22.08 ;  given_sigma_ion_N2_to_N_cm2( 517) =    1.04
 given_sigma_ion_N2_to_N2_cm2( 518) =   22.09 ;  given_sigma_ion_N2_to_N_cm2( 518) =    1.05
 given_sigma_ion_N2_to_N2_cm2( 519) =   22.10 ;  given_sigma_ion_N2_to_N_cm2( 519) =    1.06
 given_sigma_ion_N2_to_N2_cm2( 520) =   22.10 ;  given_sigma_ion_N2_to_N_cm2( 520) =    1.06
 given_sigma_ion_N2_to_N2_cm2( 521) =   22.11 ;  given_sigma_ion_N2_to_N_cm2( 521) =    1.07
 given_sigma_ion_N2_to_N2_cm2( 522) =   22.12 ;  given_sigma_ion_N2_to_N_cm2( 522) =    1.08
 given_sigma_ion_N2_to_N2_cm2( 523) =   22.13 ;  given_sigma_ion_N2_to_N_cm2( 523) =    1.10
 given_sigma_ion_N2_to_N2_cm2( 524) =   22.14 ;  given_sigma_ion_N2_to_N_cm2( 524) =    1.11
 given_sigma_ion_N2_to_N2_cm2( 525) =   22.15 ;  given_sigma_ion_N2_to_N_cm2( 525) =    1.12
 given_sigma_ion_N2_to_N2_cm2( 526) =   22.15 ;  given_sigma_ion_N2_to_N_cm2( 526) =    1.12
 given_sigma_ion_N2_to_N2_cm2( 527) =   22.17 ;  given_sigma_ion_N2_to_N_cm2( 527) =    1.13
 given_sigma_ion_N2_to_N2_cm2( 528) =   22.19 ;  given_sigma_ion_N2_to_N_cm2( 528) =    1.15
 given_sigma_ion_N2_to_N2_cm2( 529) =   22.21 ;  given_sigma_ion_N2_to_N_cm2( 529) =    1.16
 given_sigma_ion_N2_to_N2_cm2( 530) =   22.21 ;  given_sigma_ion_N2_to_N_cm2( 530) =    1.16
 given_sigma_ion_N2_to_N2_cm2( 531) =   22.23 ;  given_sigma_ion_N2_to_N_cm2( 531) =    1.17
 given_sigma_ion_N2_to_N2_cm2( 532) =   22.23 ;  given_sigma_ion_N2_to_N_cm2( 532) =    1.17
 given_sigma_ion_N2_to_N2_cm2( 533) =   22.25 ;  given_sigma_ion_N2_to_N_cm2( 533) =    1.18
 given_sigma_ion_N2_to_N2_cm2( 534) =   22.27 ;  given_sigma_ion_N2_to_N_cm2( 534) =    1.20
 given_sigma_ion_N2_to_N2_cm2( 535) =   22.27 ;  given_sigma_ion_N2_to_N_cm2( 535) =    1.20
 given_sigma_ion_N2_to_N2_cm2( 536) =   22.28 ;  given_sigma_ion_N2_to_N_cm2( 536) =    1.20
 given_sigma_ion_N2_to_N2_cm2( 537) =   22.29 ;  given_sigma_ion_N2_to_N_cm2( 537) =    1.21
 given_sigma_ion_N2_to_N2_cm2( 538) =   22.31 ;  given_sigma_ion_N2_to_N_cm2( 538) =    1.22
 given_sigma_ion_N2_to_N2_cm2( 539) =   22.32 ;  given_sigma_ion_N2_to_N_cm2( 539) =    1.23
 given_sigma_ion_N2_to_N2_cm2( 540) =   22.33 ;  given_sigma_ion_N2_to_N_cm2( 540) =    1.24
 given_sigma_ion_N2_to_N2_cm2( 541) =   22.35 ;  given_sigma_ion_N2_to_N_cm2( 541) =    1.25
 given_sigma_ion_N2_to_N2_cm2( 542) =   22.37 ;  given_sigma_ion_N2_to_N_cm2( 542) =    1.22
 given_sigma_ion_N2_to_N2_cm2( 543) =   22.38 ;  given_sigma_ion_N2_to_N_cm2( 543) =    1.21
 given_sigma_ion_N2_to_N2_cm2( 544) =   22.39 ;  given_sigma_ion_N2_to_N_cm2( 544) =    1.19
 given_sigma_ion_N2_to_N2_cm2( 545) =   22.40 ;  given_sigma_ion_N2_to_N_cm2( 545) =    1.18
 given_sigma_ion_N2_to_N2_cm2( 546) =   22.40 ;  given_sigma_ion_N2_to_N_cm2( 546) =    1.18
 given_sigma_ion_N2_to_N2_cm2( 547) =   22.42 ;  given_sigma_ion_N2_to_N_cm2( 547) =    1.16
 given_sigma_ion_N2_to_N2_cm2( 548) =   22.43 ;  given_sigma_ion_N2_to_N_cm2( 548) =    1.14
 given_sigma_ion_N2_to_N2_cm2( 549) =   22.43 ;  given_sigma_ion_N2_to_N_cm2( 549) =    1.13
 given_sigma_ion_N2_to_N2_cm2( 550) =   22.45 ;  given_sigma_ion_N2_to_N_cm2( 550) =    1.11
 given_sigma_ion_N2_to_N2_cm2( 551) =   22.47 ;  given_sigma_ion_N2_to_N_cm2( 551) =    1.09
 given_sigma_ion_N2_to_N2_cm2( 552) =   22.48 ;  given_sigma_ion_N2_to_N_cm2( 552) =    1.07
 given_sigma_ion_N2_to_N2_cm2( 553) =   22.50 ;  given_sigma_ion_N2_to_N_cm2( 553) =    1.04
 given_sigma_ion_N2_to_N2_cm2( 554) =   22.52 ;  given_sigma_ion_N2_to_N_cm2( 554) =    1.02
 given_sigma_ion_N2_to_N2_cm2( 555) =   22.53 ;  given_sigma_ion_N2_to_N_cm2( 555) =    1.00
 given_sigma_ion_N2_to_N2_cm2( 556) =   22.55 ;  given_sigma_ion_N2_to_N_cm2( 556) =    0.97
 given_sigma_ion_N2_to_N2_cm2( 557) =   22.58 ;  given_sigma_ion_N2_to_N_cm2( 557) =    0.93
 given_sigma_ion_N2_to_N2_cm2( 558) =   22.59 ;  given_sigma_ion_N2_to_N_cm2( 558) =    0.92
 given_sigma_ion_N2_to_N2_cm2( 559) =   22.60 ;  given_sigma_ion_N2_to_N_cm2( 559) =    0.90
 given_sigma_ion_N2_to_N2_cm2( 560) =   22.65 ;  given_sigma_ion_N2_to_N_cm2( 560) =    0.85
 given_sigma_ion_N2_to_N2_cm2( 561) =   22.68 ;  given_sigma_ion_N2_to_N_cm2( 561) =    0.81
 given_sigma_ion_N2_to_N2_cm2( 562) =   22.70 ;  given_sigma_ion_N2_to_N_cm2( 562) =    0.80
 given_sigma_ion_N2_to_N2_cm2( 563) =   22.75 ;  given_sigma_ion_N2_to_N_cm2( 563) =    0.75
 given_sigma_ion_N2_to_N2_cm2( 564) =   22.80 ;  given_sigma_ion_N2_to_N_cm2( 564) =    0.70
 given_sigma_ion_N2_to_N2_cm2( 565) =   22.85 ;  given_sigma_ion_N2_to_N_cm2( 565) =    0.65
 given_sigma_ion_N2_to_N2_cm2( 566) =   22.91 ;  given_sigma_ion_N2_to_N_cm2( 566) =    0.59
 given_sigma_ion_N2_to_N2_cm2( 567) =   22.97 ;  given_sigma_ion_N2_to_N_cm2( 567) =    0.53
 given_sigma_ion_N2_to_N2_cm2( 568) =   23.03 ;  given_sigma_ion_N2_to_N_cm2( 568) =    0.47
 given_sigma_ion_N2_to_N2_cm2( 569) =   23.09 ;  given_sigma_ion_N2_to_N_cm2( 569) =    0.41
 given_sigma_ion_N2_to_N2_cm2( 570) =   23.11 ;  given_sigma_ion_N2_to_N_cm2( 570) =    0.39
 given_sigma_ion_N2_to_N2_cm2( 571) =   23.11 ;  given_sigma_ion_N2_to_N_cm2( 571) =    0.39
 given_sigma_ion_N2_to_N2_cm2( 572) =   23.15 ;  given_sigma_ion_N2_to_N_cm2( 572) =    0.35
 given_sigma_ion_N2_to_N2_cm2( 573) =   23.22 ;  given_sigma_ion_N2_to_N_cm2( 573) =    0.30
 given_sigma_ion_N2_to_N2_cm2( 574) =   23.23 ;  given_sigma_ion_N2_to_N_cm2( 574) =    0.29
 given_sigma_ion_N2_to_N2_cm2( 575) =   23.29 ;  given_sigma_ion_N2_to_N_cm2( 575) =    0.24
 given_sigma_ion_N2_to_N2_cm2( 576) =   23.36 ;  given_sigma_ion_N2_to_N_cm2( 576) =    0.19
 given_sigma_ion_N2_to_N2_cm2( 577) =   23.43 ;  given_sigma_ion_N2_to_N_cm2( 577) =    0.13
 given_sigma_ion_N2_to_N2_cm2( 578) =   23.50 ;  given_sigma_ion_N2_to_N_cm2( 578) =    0.08
 given_sigma_ion_N2_to_N2_cm2( 579) =   23.67 ;  given_sigma_ion_N2_to_N_cm2( 579) =    0.05
 given_sigma_ion_N2_to_N2_cm2( 580) =   23.80 ;  given_sigma_ion_N2_to_N_cm2( 580) =    0.03
 given_sigma_ion_N2_to_N2_cm2( 581) =   24.25 ;  given_sigma_ion_N2_to_N_cm2( 581) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 582) =   24.58 ;  given_sigma_ion_N2_to_N_cm2( 582) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 583) =   24.61 ;  given_sigma_ion_N2_to_N_cm2( 583) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 584) =   24.63 ;  given_sigma_ion_N2_to_N_cm2( 584) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 585) =   24.86 ;  given_sigma_ion_N2_to_N_cm2( 585) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 586) =   25.07 ;  given_sigma_ion_N2_to_N_cm2( 586) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 587) =   25.23 ;  given_sigma_ion_N2_to_N_cm2( 587) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 588) =   25.30 ;  given_sigma_ion_N2_to_N_cm2( 588) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 589) =   25.13 ;  given_sigma_ion_N2_to_N_cm2( 589) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 590) =   25.02 ;  given_sigma_ion_N2_to_N_cm2( 590) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 591) =   24.77 ;  given_sigma_ion_N2_to_N_cm2( 591) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 592) =   24.70 ;  given_sigma_ion_N2_to_N_cm2( 592) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 593) =   24.52 ;  given_sigma_ion_N2_to_N_cm2( 593) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 594) =   24.13 ;  given_sigma_ion_N2_to_N_cm2( 594) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 595) =   24.11 ;  given_sigma_ion_N2_to_N_cm2( 595) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 596) =   23.97 ;  given_sigma_ion_N2_to_N_cm2( 596) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 597) =   23.58 ;  given_sigma_ion_N2_to_N_cm2( 597) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 598) =   23.40 ;  given_sigma_ion_N2_to_N_cm2( 598) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 599) =   23.15 ;  given_sigma_ion_N2_to_N_cm2( 599) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 600) =   22.95 ;  given_sigma_ion_N2_to_N_cm2( 600) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 601) =   22.64 ;  given_sigma_ion_N2_to_N_cm2( 601) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 602) =   22.50 ;  given_sigma_ion_N2_to_N_cm2( 602) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 603) =   22.48 ;  given_sigma_ion_N2_to_N_cm2( 603) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 604) =   22.45 ;  given_sigma_ion_N2_to_N_cm2( 604) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 605) =   22.40 ;  given_sigma_ion_N2_to_N_cm2( 605) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 606) =   22.40 ;  given_sigma_ion_N2_to_N_cm2( 606) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 607) =   22.40 ;  given_sigma_ion_N2_to_N_cm2( 607) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 608) =   22.40 ;  given_sigma_ion_N2_to_N_cm2( 608) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 609) =   22.40 ;  given_sigma_ion_N2_to_N_cm2( 609) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 610) =   22.40 ;  given_sigma_ion_N2_to_N_cm2( 610) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 611) =   22.40 ;  given_sigma_ion_N2_to_N_cm2( 611) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 612) =   22.40 ;  given_sigma_ion_N2_to_N_cm2( 612) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 613) =   22.44 ;  given_sigma_ion_N2_to_N_cm2( 613) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 614) =   22.47 ;  given_sigma_ion_N2_to_N_cm2( 614) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 615) =   22.49 ;  given_sigma_ion_N2_to_N_cm2( 615) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 616) =   22.52 ;  given_sigma_ion_N2_to_N_cm2( 616) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 617) =   22.57 ;  given_sigma_ion_N2_to_N_cm2( 617) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 618) =   22.58 ;  given_sigma_ion_N2_to_N_cm2( 618) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 619) =   22.67 ;  given_sigma_ion_N2_to_N_cm2( 619) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 620) =   22.76 ;  given_sigma_ion_N2_to_N_cm2( 620) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 621) =   22.76 ;  given_sigma_ion_N2_to_N_cm2( 621) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 622) =   22.78 ;  given_sigma_ion_N2_to_N_cm2( 622) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 623) =   22.79 ;  given_sigma_ion_N2_to_N_cm2( 623) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 624) =   22.80 ;  given_sigma_ion_N2_to_N_cm2( 624) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 625) =   22.82 ;  given_sigma_ion_N2_to_N_cm2( 625) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 626) =   22.83 ;  given_sigma_ion_N2_to_N_cm2( 626) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 627) =   22.86 ;  given_sigma_ion_N2_to_N_cm2( 627) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 628) =   22.87 ;  given_sigma_ion_N2_to_N_cm2( 628) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 629) =   22.88 ;  given_sigma_ion_N2_to_N_cm2( 629) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 630) =   22.89 ;  given_sigma_ion_N2_to_N_cm2( 630) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 631) =   22.92 ;  given_sigma_ion_N2_to_N_cm2( 631) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 632) =   22.95 ;  given_sigma_ion_N2_to_N_cm2( 632) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 633) =   22.96 ;  given_sigma_ion_N2_to_N_cm2( 633) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 634) =   22.98 ;  given_sigma_ion_N2_to_N_cm2( 634) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 635) =   22.99 ;  given_sigma_ion_N2_to_N_cm2( 635) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 636) =   23.00 ;  given_sigma_ion_N2_to_N_cm2( 636) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 637) =   23.01 ;  given_sigma_ion_N2_to_N_cm2( 637) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 638) =   23.03 ;  given_sigma_ion_N2_to_N_cm2( 638) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 639) =   23.04 ;  given_sigma_ion_N2_to_N_cm2( 639) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 640) =   23.05 ;  given_sigma_ion_N2_to_N_cm2( 640) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 641) =   23.06 ;  given_sigma_ion_N2_to_N_cm2( 641) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 642) =   23.07 ;  given_sigma_ion_N2_to_N_cm2( 642) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 643) =   23.09 ;  given_sigma_ion_N2_to_N_cm2( 643) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 644) =   23.10 ;  given_sigma_ion_N2_to_N_cm2( 644) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 645) =   23.11 ;  given_sigma_ion_N2_to_N_cm2( 645) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 646) =   23.13 ;  given_sigma_ion_N2_to_N_cm2( 646) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 647) =   23.15 ;  given_sigma_ion_N2_to_N_cm2( 647) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 648) =   23.18 ;  given_sigma_ion_N2_to_N_cm2( 648) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 649) =   23.21 ;  given_sigma_ion_N2_to_N_cm2( 649) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 650) =   23.23 ;  given_sigma_ion_N2_to_N_cm2( 650) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 651) =   23.24 ;  given_sigma_ion_N2_to_N_cm2( 651) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 652) =   23.27 ;  given_sigma_ion_N2_to_N_cm2( 652) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 653) =   23.28 ;  given_sigma_ion_N2_to_N_cm2( 653) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 654) =   23.30 ;  given_sigma_ion_N2_to_N_cm2( 654) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 655) =   23.32 ;  given_sigma_ion_N2_to_N_cm2( 655) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 656) =   23.35 ;  given_sigma_ion_N2_to_N_cm2( 656) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 657) =   23.37 ;  given_sigma_ion_N2_to_N_cm2( 657) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 658) =   23.37 ;  given_sigma_ion_N2_to_N_cm2( 658) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 659) =   23.38 ;  given_sigma_ion_N2_to_N_cm2( 659) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 660) =   23.39 ;  given_sigma_ion_N2_to_N_cm2( 660) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 661) =   23.41 ;  given_sigma_ion_N2_to_N_cm2( 661) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 662) =   23.42 ;  given_sigma_ion_N2_to_N_cm2( 662) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 663) =   23.44 ;  given_sigma_ion_N2_to_N_cm2( 663) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 664) =   23.46 ;  given_sigma_ion_N2_to_N_cm2( 664) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 665) =   23.49 ;  given_sigma_ion_N2_to_N_cm2( 665) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 666) =   23.50 ;  given_sigma_ion_N2_to_N_cm2( 666) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 667) =   23.51 ;  given_sigma_ion_N2_to_N_cm2( 667) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 668) =   23.52 ;  given_sigma_ion_N2_to_N_cm2( 668) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 669) =   23.54 ;  given_sigma_ion_N2_to_N_cm2( 669) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 670) =   23.55 ;  given_sigma_ion_N2_to_N_cm2( 670) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 671) =   23.56 ;  given_sigma_ion_N2_to_N_cm2( 671) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 672) =   23.58 ;  given_sigma_ion_N2_to_N_cm2( 672) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 673) =   23.58 ;  given_sigma_ion_N2_to_N_cm2( 673) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 674) =   23.60 ;  given_sigma_ion_N2_to_N_cm2( 674) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 675) =   23.62 ;  given_sigma_ion_N2_to_N_cm2( 675) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 676) =   23.62 ;  given_sigma_ion_N2_to_N_cm2( 676) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 677) =   23.63 ;  given_sigma_ion_N2_to_N_cm2( 677) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 678) =   23.66 ;  given_sigma_ion_N2_to_N_cm2( 678) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 679) =   23.67 ;  given_sigma_ion_N2_to_N_cm2( 679) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 680) =   23.69 ;  given_sigma_ion_N2_to_N_cm2( 680) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 681) =   23.71 ;  given_sigma_ion_N2_to_N_cm2( 681) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 682) =   23.72 ;  given_sigma_ion_N2_to_N_cm2( 682) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 683) =   23.73 ;  given_sigma_ion_N2_to_N_cm2( 683) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 684) =   23.75 ;  given_sigma_ion_N2_to_N_cm2( 684) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 685) =   23.78 ;  given_sigma_ion_N2_to_N_cm2( 685) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 686) =   23.78 ;  given_sigma_ion_N2_to_N_cm2( 686) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 687) =   23.79 ;  given_sigma_ion_N2_to_N_cm2( 687) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 688) =   23.81 ;  given_sigma_ion_N2_to_N_cm2( 688) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 689) =   23.83 ;  given_sigma_ion_N2_to_N_cm2( 689) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 690) =   23.85 ;  given_sigma_ion_N2_to_N_cm2( 690) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 691) =   23.85 ;  given_sigma_ion_N2_to_N_cm2( 691) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 692) =   23.86 ;  given_sigma_ion_N2_to_N_cm2( 692) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 693) =   23.88 ;  given_sigma_ion_N2_to_N_cm2( 693) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 694) =   23.89 ;  given_sigma_ion_N2_to_N_cm2( 694) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 695) =   23.92 ;  given_sigma_ion_N2_to_N_cm2( 695) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 696) =   23.93 ;  given_sigma_ion_N2_to_N_cm2( 696) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 697) =   23.95 ;  given_sigma_ion_N2_to_N_cm2( 697) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 698) =   23.96 ;  given_sigma_ion_N2_to_N_cm2( 698) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 699) =   23.98 ;  given_sigma_ion_N2_to_N_cm2( 699) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 700) =   24.00 ;  given_sigma_ion_N2_to_N_cm2( 700) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 701) =   24.00 ;  given_sigma_ion_N2_to_N_cm2( 701) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 702) =   24.03 ;  given_sigma_ion_N2_to_N_cm2( 702) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 703) =   24.05 ;  given_sigma_ion_N2_to_N_cm2( 703) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 704) =   24.08 ;  given_sigma_ion_N2_to_N_cm2( 704) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 705) =   24.10 ;  given_sigma_ion_N2_to_N_cm2( 705) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 706) =   24.13 ;  given_sigma_ion_N2_to_N_cm2( 706) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 707) =   24.13 ;  given_sigma_ion_N2_to_N_cm2( 707) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 708) =   24.15 ;  given_sigma_ion_N2_to_N_cm2( 708) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 709) =   24.18 ;  given_sigma_ion_N2_to_N_cm2( 709) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 710) =   24.20 ;  given_sigma_ion_N2_to_N_cm2( 710) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 711) =   24.30 ;  given_sigma_ion_N2_to_N_cm2( 711) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 712) =   23.00 ;  given_sigma_ion_N2_to_N_cm2( 712) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 713) =   23.36 ;  given_sigma_ion_N2_to_N_cm2( 713) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 714) =   23.80 ;  given_sigma_ion_N2_to_N_cm2( 714) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 715) =   21.00 ;  given_sigma_ion_N2_to_N_cm2( 715) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 716) =   24.30 ;  given_sigma_ion_N2_to_N_cm2( 716) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 717) =   25.20 ;  given_sigma_ion_N2_to_N_cm2( 717) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 718) =   23.80 ;  given_sigma_ion_N2_to_N_cm2( 718) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 719) =   25.60 ;  given_sigma_ion_N2_to_N_cm2( 719) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 720) =   25.60 ;  given_sigma_ion_N2_to_N_cm2( 720) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 721) =   25.60 ;  given_sigma_ion_N2_to_N_cm2( 721) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 722) =   28.20 ;  given_sigma_ion_N2_to_N_cm2( 722) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 723) =   38.10 ;  given_sigma_ion_N2_to_N_cm2( 723) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 724) =   21.20 ;  given_sigma_ion_N2_to_N_cm2( 724) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 725) =   46.10 ;  given_sigma_ion_N2_to_N_cm2( 725) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 726) =   30.86 ;  given_sigma_ion_N2_to_N_cm2( 726) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 727) =   20.70 ;  given_sigma_ion_N2_to_N_cm2( 727) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 728) =   21.40 ;  given_sigma_ion_N2_to_N_cm2( 728) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 729) =   26.40 ;  given_sigma_ion_N2_to_N_cm2( 729) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 730) =   33.30 ;  given_sigma_ion_N2_to_N_cm2( 730) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 731) =   27.50 ;  given_sigma_ion_N2_to_N_cm2( 731) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 732) =   21.40 ;  given_sigma_ion_N2_to_N_cm2( 732) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 733) =   20.83 ;  given_sigma_ion_N2_to_N_cm2( 733) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 734) =   20.50 ;  given_sigma_ion_N2_to_N_cm2( 734) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 735) =   14.80 ;  given_sigma_ion_N2_to_N_cm2( 735) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 736) =   27.60 ;  given_sigma_ion_N2_to_N_cm2( 736) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 737) =   28.70 ;  given_sigma_ion_N2_to_N_cm2( 737) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 738) =   61.40 ;  given_sigma_ion_N2_to_N_cm2( 738) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 739) =   38.71 ;  given_sigma_ion_N2_to_N_cm2( 739) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 740) =   25.10 ;  given_sigma_ion_N2_to_N_cm2( 740) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 741) =   24.08 ;  given_sigma_ion_N2_to_N_cm2( 741) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 742) =   22.04 ;  given_sigma_ion_N2_to_N_cm2( 742) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 743) =   20.00 ;  given_sigma_ion_N2_to_N_cm2( 743) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 744) =   20.85 ;  given_sigma_ion_N2_to_N_cm2( 744) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 745) =   21.53 ;  given_sigma_ion_N2_to_N_cm2( 745) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 746) =   22.03 ;  given_sigma_ion_N2_to_N_cm2( 746) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 747) =   22.58 ;  given_sigma_ion_N2_to_N_cm2( 747) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 748) =   22.80 ;  given_sigma_ion_N2_to_N_cm2( 748) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 749) =   19.68 ;  given_sigma_ion_N2_to_N_cm2( 749) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 750) =    8.76 ;  given_sigma_ion_N2_to_N_cm2( 750) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 751) =   32.80 ;  given_sigma_ion_N2_to_N_cm2( 751) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 752) =   31.95 ;  given_sigma_ion_N2_to_N_cm2( 752) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 753) =   30.67 ;  given_sigma_ion_N2_to_N_cm2( 753) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 754) =   29.40 ;  given_sigma_ion_N2_to_N_cm2( 754) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 755) =   23.50 ;  given_sigma_ion_N2_to_N_cm2( 755) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 756) =   30.49 ;  given_sigma_ion_N2_to_N_cm2( 756) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 757) =   44.46 ;  given_sigma_ion_N2_to_N_cm2( 757) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 758) =   51.44 ;  given_sigma_ion_N2_to_N_cm2( 758) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 759) =   72.40 ;  given_sigma_ion_N2_to_N_cm2( 759) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 760) =   59.38 ;  given_sigma_ion_N2_to_N_cm2( 760) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 761) =   37.68 ;  given_sigma_ion_N2_to_N_cm2( 761) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 762) =   29.00 ;  given_sigma_ion_N2_to_N_cm2( 762) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 763) =   23.30 ;  given_sigma_ion_N2_to_N_cm2( 763) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 764) =   22.97 ;  given_sigma_ion_N2_to_N_cm2( 764) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 765) =   22.90 ;  given_sigma_ion_N2_to_N_cm2( 765) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 766) =   22.63 ;  given_sigma_ion_N2_to_N_cm2( 766) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 767) =   22.50 ;  given_sigma_ion_N2_to_N_cm2( 767) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 768) =   22.30 ;  given_sigma_ion_N2_to_N_cm2( 768) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 769) =   24.35 ;  given_sigma_ion_N2_to_N_cm2( 769) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 770) =   23.96 ;  given_sigma_ion_N2_to_N_cm2( 770) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 771) =   22.94 ;  given_sigma_ion_N2_to_N_cm2( 771) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 772) =   22.50 ;  given_sigma_ion_N2_to_N_cm2( 772) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 773) =   22.50 ;  given_sigma_ion_N2_to_N_cm2( 773) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 774) =   21.71 ;  given_sigma_ion_N2_to_N_cm2( 774) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 775) =   21.45 ;  given_sigma_ion_N2_to_N_cm2( 775) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 776) =   20.67 ;  given_sigma_ion_N2_to_N_cm2( 776) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 777) =   19.88 ;  given_sigma_ion_N2_to_N_cm2( 777) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 778) =   18.31 ;  given_sigma_ion_N2_to_N_cm2( 778) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 779) =   17.00 ;  given_sigma_ion_N2_to_N_cm2( 779) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 780) =   16.40 ;  given_sigma_ion_N2_to_N_cm2( 780) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 781) =   15.80 ;  given_sigma_ion_N2_to_N_cm2( 781) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 782) =   13.70 ;  given_sigma_ion_N2_to_N_cm2( 782) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 783) =    9.08 ;  given_sigma_ion_N2_to_N_cm2( 783) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 784) =    7.54 ;  given_sigma_ion_N2_to_N_cm2( 784) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 785) =    6.00 ;  given_sigma_ion_N2_to_N_cm2( 785) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 786) =   12.20 ;  given_sigma_ion_N2_to_N_cm2( 786) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 787) =   17.90 ;  given_sigma_ion_N2_to_N_cm2( 787) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 788) =   18.10 ;  given_sigma_ion_N2_to_N_cm2( 788) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 789) =   18.43 ;  given_sigma_ion_N2_to_N_cm2( 789) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 790) =   20.03 ;  given_sigma_ion_N2_to_N_cm2( 790) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 791) =   24.31 ;  given_sigma_ion_N2_to_N_cm2( 791) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 792) =   28.60 ;  given_sigma_ion_N2_to_N_cm2( 792) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 793) =   35.52 ;  given_sigma_ion_N2_to_N_cm2( 793) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 794) =   37.50 ;  given_sigma_ion_N2_to_N_cm2( 794) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 795) =   46.59 ;  given_sigma_ion_N2_to_N_cm2( 795) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 796) =   64.77 ;  given_sigma_ion_N2_to_N_cm2( 796) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 797) =   62.98 ;  given_sigma_ion_N2_to_N_cm2( 797) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 798) =   43.70 ;  given_sigma_ion_N2_to_N_cm2( 798) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 799) =   32.90 ;  given_sigma_ion_N2_to_N_cm2( 799) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 800) =   22.10 ;  given_sigma_ion_N2_to_N_cm2( 800) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 801) =   22.25 ;  given_sigma_ion_N2_to_N_cm2( 801) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 802) =   22.35 ;  given_sigma_ion_N2_to_N_cm2( 802) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 803) =   22.45 ;  given_sigma_ion_N2_to_N_cm2( 803) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 804) =   22.60 ;  given_sigma_ion_N2_to_N_cm2( 804) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 805) =   22.51 ;  given_sigma_ion_N2_to_N_cm2( 805) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 806) =   22.33 ;  given_sigma_ion_N2_to_N_cm2( 806) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 807) =   20.80 ;  given_sigma_ion_N2_to_N_cm2( 807) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 808) =   20.80 ;  given_sigma_ion_N2_to_N_cm2( 808) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 809) =   20.80 ;  given_sigma_ion_N2_to_N_cm2( 809) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 810) =   22.40 ;  given_sigma_ion_N2_to_N_cm2( 810) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 811) =   20.80 ;  given_sigma_ion_N2_to_N_cm2( 811) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 812) =   21.27 ;  given_sigma_ion_N2_to_N_cm2( 812) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 813) =   23.16 ;  given_sigma_ion_N2_to_N_cm2( 813) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 814) =   24.10 ;  given_sigma_ion_N2_to_N_cm2( 814) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 815) =   25.00 ;  given_sigma_ion_N2_to_N_cm2( 815) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 816) =   26.80 ;  given_sigma_ion_N2_to_N_cm2( 816) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 817) =   27.70 ;  given_sigma_ion_N2_to_N_cm2( 817) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 818) =   28.60 ;  given_sigma_ion_N2_to_N_cm2( 818) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 819) =   27.10 ;  given_sigma_ion_N2_to_N_cm2( 819) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 820) =   25.48 ;  given_sigma_ion_N2_to_N_cm2( 820) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 821) =   24.40 ;  given_sigma_ion_N2_to_N_cm2( 821) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 822) =   22.60 ;  given_sigma_ion_N2_to_N_cm2( 822) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 823) =   25.00 ;  given_sigma_ion_N2_to_N_cm2( 823) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 824) =   22.30 ;  given_sigma_ion_N2_to_N_cm2( 824) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 825) =   21.91 ;  given_sigma_ion_N2_to_N_cm2( 825) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 826) =   21.13 ;  given_sigma_ion_N2_to_N_cm2( 826) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 827) =   20.16 ;  given_sigma_ion_N2_to_N_cm2( 827) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 828) =   19.60 ;  given_sigma_ion_N2_to_N_cm2( 828) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 829) =   19.24 ;  given_sigma_ion_N2_to_N_cm2( 829) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 830) =   19.06 ;  given_sigma_ion_N2_to_N_cm2( 830) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 831) =   18.19 ;  given_sigma_ion_N2_to_N_cm2( 831) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 832) =   18.04 ;  given_sigma_ion_N2_to_N_cm2( 832) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 833) =   17.89 ;  given_sigma_ion_N2_to_N_cm2( 833) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 834) =   17.80 ;  given_sigma_ion_N2_to_N_cm2( 834) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 835) =   17.96 ;  given_sigma_ion_N2_to_N_cm2( 835) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 836) =   18.00 ;  given_sigma_ion_N2_to_N_cm2( 836) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 837) =   18.00 ;  given_sigma_ion_N2_to_N_cm2( 837) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 838) =   17.30 ;  given_sigma_ion_N2_to_N_cm2( 838) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 839) =   17.02 ;  given_sigma_ion_N2_to_N_cm2( 839) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 840) =   16.74 ;  given_sigma_ion_N2_to_N_cm2( 840) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 841) =   15.40 ;  given_sigma_ion_N2_to_N_cm2( 841) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 842) =   14.44 ;  given_sigma_ion_N2_to_N_cm2( 842) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 843) =    7.41 ;  given_sigma_ion_N2_to_N_cm2( 843) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 844) =    4.50 ;  given_sigma_ion_N2_to_N_cm2( 844) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 845) =    3.60 ;  given_sigma_ion_N2_to_N_cm2( 845) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 846) =    3.45 ;  given_sigma_ion_N2_to_N_cm2( 846) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 847) =    3.00 ;  given_sigma_ion_N2_to_N_cm2( 847) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 848) =    6.75 ;  given_sigma_ion_N2_to_N_cm2( 848) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 849) =   10.50 ;  given_sigma_ion_N2_to_N_cm2( 849) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 850) =   10.60 ;  given_sigma_ion_N2_to_N_cm2( 850) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 851) =   10.80 ;  given_sigma_ion_N2_to_N_cm2( 851) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 852) =   11.00 ;  given_sigma_ion_N2_to_N_cm2( 852) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 853) =   17.55 ;  given_sigma_ion_N2_to_N_cm2( 853) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 854) =   24.10 ;  given_sigma_ion_N2_to_N_cm2( 854) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 855) =   22.08 ;  given_sigma_ion_N2_to_N_cm2( 855) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 856) =   20.06 ;  given_sigma_ion_N2_to_N_cm2( 856) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 857) =   14.00 ;  given_sigma_ion_N2_to_N_cm2( 857) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 858) =   13.48 ;  given_sigma_ion_N2_to_N_cm2( 858) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 859) =   12.83 ;  given_sigma_ion_N2_to_N_cm2( 859) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 860) =   21.17 ;  given_sigma_ion_N2_to_N_cm2( 860) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 861) =   24.00 ;  given_sigma_ion_N2_to_N_cm2( 861) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 862) =   18.53 ;  given_sigma_ion_N2_to_N_cm2( 862) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 863) =   15.80 ;  given_sigma_ion_N2_to_N_cm2( 863) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 864) =   25.70 ;  given_sigma_ion_N2_to_N_cm2( 864) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 865) =   33.62 ;  given_sigma_ion_N2_to_N_cm2( 865) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 866) =   34.60 ;  given_sigma_ion_N2_to_N_cm2( 866) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 867) =   33.35 ;  given_sigma_ion_N2_to_N_cm2( 867) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 868) =   31.00 ;  given_sigma_ion_N2_to_N_cm2( 868) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 869) =   24.70 ;  given_sigma_ion_N2_to_N_cm2( 869) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 870) =   23.65 ;  given_sigma_ion_N2_to_N_cm2( 870) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 871) =   19.15 ;  given_sigma_ion_N2_to_N_cm2( 871) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 872) =   17.77 ;  given_sigma_ion_N2_to_N_cm2( 872) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 873) =   15.70 ;  given_sigma_ion_N2_to_N_cm2( 873) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 874) =   16.10 ;  given_sigma_ion_N2_to_N_cm2( 874) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 875) =   16.70 ;  given_sigma_ion_N2_to_N_cm2( 875) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 876) =   13.82 ;  given_sigma_ion_N2_to_N_cm2( 876) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 877) =   11.90 ;  given_sigma_ion_N2_to_N_cm2( 877) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 878) =    7.10 ;  given_sigma_ion_N2_to_N_cm2( 878) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 879) =    6.86 ;  given_sigma_ion_N2_to_N_cm2( 879) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 880) =    6.54 ;  given_sigma_ion_N2_to_N_cm2( 880) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 881) =    6.30 ;  given_sigma_ion_N2_to_N_cm2( 881) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 882) =   10.54 ;  given_sigma_ion_N2_to_N_cm2( 882) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 883) =   14.78 ;  given_sigma_ion_N2_to_N_cm2( 883) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 884) =   16.90 ;  given_sigma_ion_N2_to_N_cm2( 884) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 885) =   14.02 ;  given_sigma_ion_N2_to_N_cm2( 885) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 886) =   12.58 ;  given_sigma_ion_N2_to_N_cm2( 886) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 887) =    9.04 ;  given_sigma_ion_N2_to_N_cm2( 887) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 888) =   11.48 ;  given_sigma_ion_N2_to_N_cm2( 888) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 889) =   13.10 ;  given_sigma_ion_N2_to_N_cm2( 889) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 890) =   12.24 ;  given_sigma_ion_N2_to_N_cm2( 890) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 891) =   10.95 ;  given_sigma_ion_N2_to_N_cm2( 891) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 892) =    8.80 ;  given_sigma_ion_N2_to_N_cm2( 892) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 893) =    9.22 ;  given_sigma_ion_N2_to_N_cm2( 893) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 894) =   10.20 ;  given_sigma_ion_N2_to_N_cm2( 894) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 895) =   15.10 ;  given_sigma_ion_N2_to_N_cm2( 895) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 896) =   11.30 ;  given_sigma_ion_N2_to_N_cm2( 896) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 897) =   10.55 ;  given_sigma_ion_N2_to_N_cm2( 897) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 898) =    9.05 ;  given_sigma_ion_N2_to_N_cm2( 898) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 899) =   14.80 ;  given_sigma_ion_N2_to_N_cm2( 899) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 900) =    9.20 ;  given_sigma_ion_N2_to_N_cm2( 900) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 901) =   11.96 ;  given_sigma_ion_N2_to_N_cm2( 901) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 902) =   16.10 ;  given_sigma_ion_N2_to_N_cm2( 902) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 903) =   23.70 ;  given_sigma_ion_N2_to_N_cm2( 903) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 904) =   18.18 ;  given_sigma_ion_N2_to_N_cm2( 904) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 905) =    9.90 ;  given_sigma_ion_N2_to_N_cm2( 905) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 906) =    9.92 ;  given_sigma_ion_N2_to_N_cm2( 906) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 907) =   10.00 ;  given_sigma_ion_N2_to_N_cm2( 907) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 908) =    9.50 ;  given_sigma_ion_N2_to_N_cm2( 908) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 909) =    7.50 ;  given_sigma_ion_N2_to_N_cm2( 909) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 910) =    8.02 ;  given_sigma_ion_N2_to_N_cm2( 910) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 911) =   10.10 ;  given_sigma_ion_N2_to_N_cm2( 911) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 912) =    9.24 ;  given_sigma_ion_N2_to_N_cm2( 912) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 913) =    8.38 ;  given_sigma_ion_N2_to_N_cm2( 913) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 914) =    5.80 ;  given_sigma_ion_N2_to_N_cm2( 914) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 915) =    8.72 ;  given_sigma_ion_N2_to_N_cm2( 915) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 916) =   13.10 ;  given_sigma_ion_N2_to_N_cm2( 916) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 917) =   11.42 ;  given_sigma_ion_N2_to_N_cm2( 917) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 918) =    8.90 ;  given_sigma_ion_N2_to_N_cm2( 918) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 919) =   17.60 ;  given_sigma_ion_N2_to_N_cm2( 919) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 920) =   12.56 ;  given_sigma_ion_N2_to_N_cm2( 920) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 921) =    6.30 ;  given_sigma_ion_N2_to_N_cm2( 921) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 922) =    6.14 ;  given_sigma_ion_N2_to_N_cm2( 922) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 923) =    6.02 ;  given_sigma_ion_N2_to_N_cm2( 923) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 924) =    5.90 ;  given_sigma_ion_N2_to_N_cm2( 924) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 925) =   21.90 ;  given_sigma_ion_N2_to_N_cm2( 925) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 926) =   11.30 ;  given_sigma_ion_N2_to_N_cm2( 926) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 927) =    8.85 ;  given_sigma_ion_N2_to_N_cm2( 927) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 928) =    6.40 ;  given_sigma_ion_N2_to_N_cm2( 928) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 929) =    5.12 ;  given_sigma_ion_N2_to_N_cm2( 929) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 930) =   10.00 ;  given_sigma_ion_N2_to_N_cm2( 930) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 931) =    5.85 ;  given_sigma_ion_N2_to_N_cm2( 931) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 932) =    4.10 ;  given_sigma_ion_N2_to_N_cm2( 932) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 933) =    4.76 ;  given_sigma_ion_N2_to_N_cm2( 933) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 934) =    5.20 ;  given_sigma_ion_N2_to_N_cm2( 934) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 935) =    4.99 ;  given_sigma_ion_N2_to_N_cm2( 935) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 936) =    4.57 ;  given_sigma_ion_N2_to_N_cm2( 936) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 937) =   13.80 ;  given_sigma_ion_N2_to_N_cm2( 937) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 938) =   46.12 ;  given_sigma_ion_N2_to_N_cm2( 938) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 939) =   28.50 ;  given_sigma_ion_N2_to_N_cm2( 939) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 940) =   24.74 ;  given_sigma_ion_N2_to_N_cm2( 940) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 941) =   19.10 ;  given_sigma_ion_N2_to_N_cm2( 941) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 942) =   16.22 ;  given_sigma_ion_N2_to_N_cm2( 942) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 943) =   14.30 ;  given_sigma_ion_N2_to_N_cm2( 943) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 944) =   12.22 ;  given_sigma_ion_N2_to_N_cm2( 944) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 945) =    9.10 ;  given_sigma_ion_N2_to_N_cm2( 945) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 946) =   24.20 ;  given_sigma_ion_N2_to_N_cm2( 946) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 947) =   20.24 ;  given_sigma_ion_N2_to_N_cm2( 947) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 948) =   14.30 ;  given_sigma_ion_N2_to_N_cm2( 948) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 949) =   12.70 ;  given_sigma_ion_N2_to_N_cm2( 949) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 950) =   10.30 ;  given_sigma_ion_N2_to_N_cm2( 950) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 951) =   37.60 ;  given_sigma_ion_N2_to_N_cm2( 951) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 952) =   32.20 ;  given_sigma_ion_N2_to_N_cm2( 952) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 953) =   19.15 ;  given_sigma_ion_N2_to_N_cm2( 953) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 954) =   14.50 ;  given_sigma_ion_N2_to_N_cm2( 954) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 955) =   10.20 ;  given_sigma_ion_N2_to_N_cm2( 955) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 956) =   12.30 ;  given_sigma_ion_N2_to_N_cm2( 956) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 957) =   15.80 ;  given_sigma_ion_N2_to_N_cm2( 957) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 958) =   13.27 ;  given_sigma_ion_N2_to_N_cm2( 958) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 959) =   17.00 ;  given_sigma_ion_N2_to_N_cm2( 959) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 960) =   24.50 ;  given_sigma_ion_N2_to_N_cm2( 960) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 961) =   18.80 ;  given_sigma_ion_N2_to_N_cm2( 961) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 962) =   11.70 ;  given_sigma_ion_N2_to_N_cm2( 962) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 963) =   10.46 ;  given_sigma_ion_N2_to_N_cm2( 963) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 964) =    8.60 ;  given_sigma_ion_N2_to_N_cm2( 964) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 965) =   10.60 ;  given_sigma_ion_N2_to_N_cm2( 965) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 966) =   18.73 ;  given_sigma_ion_N2_to_N_cm2( 966) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 967) =   21.40 ;  given_sigma_ion_N2_to_N_cm2( 967) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 968) =   38.30 ;  given_sigma_ion_N2_to_N_cm2( 968) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 969) =   22.20 ;  given_sigma_ion_N2_to_N_cm2( 969) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 970) =    9.20 ;  given_sigma_ion_N2_to_N_cm2( 970) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 971) =   11.40 ;  given_sigma_ion_N2_to_N_cm2( 971) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 972) =  106.90 ;  given_sigma_ion_N2_to_N_cm2( 972) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 973) =   52.90 ;  given_sigma_ion_N2_to_N_cm2( 973) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 974) =   44.30 ;  given_sigma_ion_N2_to_N_cm2( 974) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 975) =   27.02 ;  given_sigma_ion_N2_to_N_cm2( 975) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 976) =   14.00 ;  given_sigma_ion_N2_to_N_cm2( 976) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 977) =   11.40 ;  given_sigma_ion_N2_to_N_cm2( 977) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 978) =   10.72 ;  given_sigma_ion_N2_to_N_cm2( 978) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 979) =    9.70 ;  given_sigma_ion_N2_to_N_cm2( 979) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 980) =   10.06 ;  given_sigma_ion_N2_to_N_cm2( 980) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 981) =   10.30 ;  given_sigma_ion_N2_to_N_cm2( 981) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 982) =   10.86 ;  given_sigma_ion_N2_to_N_cm2( 982) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 983) =   11.42 ;  given_sigma_ion_N2_to_N_cm2( 983) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 984) =   14.58 ;  given_sigma_ion_N2_to_N_cm2( 984) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 985) =   15.54 ;  given_sigma_ion_N2_to_N_cm2( 985) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 986) =   28.50 ;  given_sigma_ion_N2_to_N_cm2( 986) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 987) =   16.00 ;  given_sigma_ion_N2_to_N_cm2( 987) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 988) =   14.72 ;  given_sigma_ion_N2_to_N_cm2( 988) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 989) =   12.80 ;  given_sigma_ion_N2_to_N_cm2( 989) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 990) =    9.40 ;  given_sigma_ion_N2_to_N_cm2( 990) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 991) =    8.96 ;  given_sigma_ion_N2_to_N_cm2( 991) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 992) =    8.52 ;  given_sigma_ion_N2_to_N_cm2( 992) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 993) =   10.16 ;  given_sigma_ion_N2_to_N_cm2( 993) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 994) =   11.40 ;  given_sigma_ion_N2_to_N_cm2( 994) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 995) =   51.30 ;  given_sigma_ion_N2_to_N_cm2( 995) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 996) =   26.00 ;  given_sigma_ion_N2_to_N_cm2( 996) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 997) =   14.80 ;  given_sigma_ion_N2_to_N_cm2( 997) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 998) =   10.30 ;  given_sigma_ion_N2_to_N_cm2( 998) =    0.00
 given_sigma_ion_N2_to_N2_cm2( 999) =    8.30 ;  given_sigma_ion_N2_to_N_cm2( 999) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1000) =    8.00 ;  given_sigma_ion_N2_to_N_cm2(1000) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1001) =   13.20 ;  given_sigma_ion_N2_to_N_cm2(1001) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1002) =   12.50 ;  given_sigma_ion_N2_to_N_cm2(1002) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1003) =   48.20 ;  given_sigma_ion_N2_to_N_cm2(1003) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1004) =  126.90 ;  given_sigma_ion_N2_to_N_cm2(1004) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1005) =   63.30 ;  given_sigma_ion_N2_to_N_cm2(1005) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1006) =   20.30 ;  given_sigma_ion_N2_to_N_cm2(1006) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1007) =   14.05 ;  given_sigma_ion_N2_to_N_cm2(1007) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1008) =    7.80 ;  given_sigma_ion_N2_to_N_cm2(1008) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1009) =    7.20 ;  given_sigma_ion_N2_to_N_cm2(1009) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1010) =    9.50 ;  given_sigma_ion_N2_to_N_cm2(1010) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1011) =    9.14 ;  given_sigma_ion_N2_to_N_cm2(1011) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1012) =    8.60 ;  given_sigma_ion_N2_to_N_cm2(1012) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1013) =    8.10 ;  given_sigma_ion_N2_to_N_cm2(1013) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1014) =    7.70 ;  given_sigma_ion_N2_to_N_cm2(1014) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1015) =    7.35 ;  given_sigma_ion_N2_to_N_cm2(1015) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1016) =    8.50 ;  given_sigma_ion_N2_to_N_cm2(1016) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1017) =    9.50 ;  given_sigma_ion_N2_to_N_cm2(1017) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1018) =   41.70 ;  given_sigma_ion_N2_to_N_cm2(1018) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1019) =   70.70 ;  given_sigma_ion_N2_to_N_cm2(1019) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1020) =   32.30 ;  given_sigma_ion_N2_to_N_cm2(1020) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1021) =   13.40 ;  given_sigma_ion_N2_to_N_cm2(1021) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1022) =    6.50 ;  given_sigma_ion_N2_to_N_cm2(1022) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1023) =    7.58 ;  given_sigma_ion_N2_to_N_cm2(1023) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1024) =   29.96 ;  given_sigma_ion_N2_to_N_cm2(1024) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1025) =   63.20 ;  given_sigma_ion_N2_to_N_cm2(1025) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1026) =  125.30 ;  given_sigma_ion_N2_to_N_cm2(1026) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1027) =   48.40 ;  given_sigma_ion_N2_to_N_cm2(1027) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1028) =   18.80 ;  given_sigma_ion_N2_to_N_cm2(1028) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1029) =    8.40 ;  given_sigma_ion_N2_to_N_cm2(1029) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1030) =    6.40 ;  given_sigma_ion_N2_to_N_cm2(1030) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1031) =    9.70 ;  given_sigma_ion_N2_to_N_cm2(1031) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1032) =   22.80 ;  given_sigma_ion_N2_to_N_cm2(1032) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1033) =   18.04 ;  given_sigma_ion_N2_to_N_cm2(1033) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1034) =   13.28 ;  given_sigma_ion_N2_to_N_cm2(1034) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1035) =    7.80 ;  given_sigma_ion_N2_to_N_cm2(1035) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1036) =    7.30 ;  given_sigma_ion_N2_to_N_cm2(1036) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1037) =    6.94 ;  given_sigma_ion_N2_to_N_cm2(1037) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1038) =    6.40 ;  given_sigma_ion_N2_to_N_cm2(1038) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1039) =    5.80 ;  given_sigma_ion_N2_to_N_cm2(1039) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1040) =    7.20 ;  given_sigma_ion_N2_to_N_cm2(1040) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1041) =   11.30 ;  given_sigma_ion_N2_to_N_cm2(1041) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1042) =   11.00 ;  given_sigma_ion_N2_to_N_cm2(1042) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1043) =    9.60 ;  given_sigma_ion_N2_to_N_cm2(1043) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1044) =  101.90 ;  given_sigma_ion_N2_to_N_cm2(1044) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1045) =   41.00 ;  given_sigma_ion_N2_to_N_cm2(1045) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1046) =   37.94 ;  given_sigma_ion_N2_to_N_cm2(1046) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1047) =   36.92 ;  given_sigma_ion_N2_to_N_cm2(1047) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1048) =   44.90 ;  given_sigma_ion_N2_to_N_cm2(1048) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1049) =   38.90 ;  given_sigma_ion_N2_to_N_cm2(1049) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1050) =   21.62 ;  given_sigma_ion_N2_to_N_cm2(1050) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1051) =   15.91 ;  given_sigma_ion_N2_to_N_cm2(1051) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1052) =   10.56 ;  given_sigma_ion_N2_to_N_cm2(1052) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1053) =    6.60 ;  given_sigma_ion_N2_to_N_cm2(1053) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1054) =    6.20 ;  given_sigma_ion_N2_to_N_cm2(1054) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1055) =    6.30 ;  given_sigma_ion_N2_to_N_cm2(1055) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1056) =   12.30 ;  given_sigma_ion_N2_to_N_cm2(1056) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1057) =    5.92 ;  given_sigma_ion_N2_to_N_cm2(1057) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1058) =    4.50 ;  given_sigma_ion_N2_to_N_cm2(1058) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1059) =    1.00 ;  given_sigma_ion_N2_to_N_cm2(1059) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1060) =    0.30 ;  given_sigma_ion_N2_to_N_cm2(1060) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1061) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1061) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1062) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1062) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1063) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1063) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1064) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1064) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1065) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1065) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1066) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1066) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1067) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1067) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1068) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1068) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1069) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1069) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1070) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1070) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1071) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1071) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1072) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1072) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1073) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1073) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1074) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1074) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1075) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1075) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1076) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1076) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1077) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1077) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1078) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1078) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1079) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1079) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1080) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1080) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1081) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1081) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1082) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1082) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1083) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1083) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1084) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1084) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1085) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1085) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1086) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1086) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1087) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1087) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1088) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1088) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1089) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1089) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1090) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1090) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1091) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1091) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1092) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1092) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1093) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1093) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1094) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1094) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1095) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1095) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1096) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1096) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1097) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1097) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1098) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1098) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1099) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1099) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1100) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1100) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1101) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1101) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1102) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1102) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1103) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1103) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1104) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1104) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1105) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1105) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1106) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1106) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1107) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1107) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1108) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1108) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1109) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1109) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1110) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1110) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1111) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1111) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1112) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1112) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1113) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1113) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1114) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1114) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1115) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1115) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1116) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1116) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1117) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1117) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1118) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1118) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1119) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1119) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1120) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1120) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1121) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1121) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1122) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1122) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1123) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1123) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1124) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1124) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1125) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1125) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1126) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1126) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1127) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1127) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1128) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1128) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1129) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1129) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1130) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1130) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1131) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1131) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1132) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1132) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1133) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1133) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1134) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1134) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1135) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1135) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1136) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1136) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1137) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1137) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1138) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1138) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1139) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1139) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1140) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1140) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1141) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1141) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1142) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1142) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1143) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1143) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1144) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1144) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1145) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1145) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1146) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1146) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1147) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1147) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1148) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1148) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1149) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1149) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1150) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1150) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1151) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1151) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1152) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1152) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1153) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1153) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1154) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1154) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1155) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1155) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1156) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1156) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1157) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1157) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1158) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1158) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1159) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1159) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1160) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1160) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1161) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1161) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1162) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1162) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1163) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1163) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1164) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1164) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1165) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1165) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1166) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1166) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1167) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1167) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1168) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1168) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1169) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1169) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1170) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1170) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1171) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1171) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1172) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1172) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1173) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1173) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1174) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1174) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1175) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1175) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1176) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1176) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1177) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1177) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1178) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1178) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1179) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1179) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1180) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1180) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1181) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1181) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1182) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1182) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1183) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1183) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1184) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1184) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1185) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1185) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1186) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1186) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1187) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1187) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1188) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1188) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1189) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1189) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1190) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1190) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1191) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1191) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1192) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1192) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1193) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1193) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1194) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1194) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1195) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1195) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1196) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1196) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1197) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1197) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1198) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1198) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1199) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1199) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1200) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1200) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1201) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1201) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1202) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1202) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1203) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1203) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1204) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1204) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1205) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1205) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1206) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1206) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1207) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1207) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1208) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1208) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1209) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1209) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1210) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1210) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1211) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1211) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1212) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1212) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1213) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1213) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1214) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1214) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1215) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1215) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1216) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1216) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1217) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1217) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1218) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1218) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1219) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1219) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1220) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1220) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1221) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1221) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1222) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1222) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1223) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1223) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1224) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1224) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1225) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1225) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1226) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1226) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1227) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1227) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1228) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1228) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1229) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1229) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1230) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1230) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1231) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1231) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1232) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1232) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1233) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1233) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1234) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1234) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1235) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1235) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1236) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1236) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1237) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1237) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1238) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1238) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1239) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1239) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1240) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1240) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1241) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1241) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1242) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1242) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1243) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1243) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1244) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1244) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1245) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1245) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1246) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1246) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1247) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1247) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1248) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1248) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1249) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1249) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1250) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1250) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1251) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1251) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1252) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1252) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1253) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1253) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1254) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1254) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1255) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1255) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1256) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1256) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1257) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1257) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1258) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1258) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1259) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1259) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1260) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1260) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1261) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1261) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1262) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1262) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1263) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1263) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1264) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1264) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1265) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1265) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1266) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1266) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1267) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1267) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1268) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1268) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1269) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1269) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1270) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1270) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1271) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1271) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1272) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1272) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1273) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1273) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1274) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1274) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1275) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1275) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1276) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1276) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1277) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1277) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1278) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1278) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1279) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1279) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1280) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1280) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1281) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1281) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1282) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1282) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1283) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1283) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1284) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1284) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1285) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1285) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1286) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1286) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1287) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1287) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1288) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1288) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1289) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1289) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1290) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1290) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1291) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1291) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1292) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1292) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1293) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1293) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1294) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1294) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1295) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1295) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1296) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1296) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1297) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1297) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1298) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1298) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1299) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1299) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1300) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1300) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1301) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1301) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1302) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1302) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1303) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1303) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1304) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1304) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1305) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1305) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1306) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1306) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1307) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1307) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1308) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1308) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1309) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1309) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1310) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1310) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1311) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1311) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1312) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1312) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1313) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1313) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1314) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1314) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1315) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1315) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1316) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1316) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1317) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1317) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1318) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1318) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1319) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1319) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1320) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1320) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1321) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1321) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1322) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1322) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1323) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1323) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1324) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1324) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1325) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1325) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1326) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1326) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1327) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1327) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1328) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1328) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1329) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1329) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1330) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1330) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1331) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1331) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1332) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1332) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1333) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1333) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1334) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1334) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1335) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1335) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1336) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1336) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1337) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1337) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1338) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1338) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1339) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1339) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1340) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1340) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1341) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1341) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1342) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1342) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1343) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1343) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1344) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1344) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1345) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1345) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1346) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1346) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1347) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1347) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1348) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1348) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1349) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1349) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1350) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1350) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1351) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1351) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1352) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1352) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1353) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1353) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1354) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1354) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1355) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1355) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1356) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1356) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1357) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1357) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1358) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1358) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1359) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1359) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1360) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1360) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1361) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1361) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1362) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1362) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1363) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1363) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1364) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1364) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1365) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1365) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1366) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1366) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1367) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1367) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1368) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1368) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1369) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1369) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1370) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1370) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1371) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1371) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1372) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1372) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1373) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1373) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1374) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1374) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1375) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1375) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1376) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1376) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1377) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1377) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1378) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1378) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1379) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1379) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1380) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1380) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1381) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1381) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1382) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1382) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1383) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1383) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1384) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1384) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1385) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1385) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1386) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1386) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1387) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1387) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1388) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1388) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1389) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1389) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1390) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1390) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1391) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1391) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1392) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1392) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1393) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1393) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1394) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1394) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1395) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1395) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1396) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1396) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1397) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1397) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1398) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1398) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1399) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1399) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1400) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1400) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1401) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1401) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1402) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1402) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1403) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1403) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1404) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1404) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1405) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1405) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1406) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1406) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1407) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1407) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1408) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1408) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1409) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1409) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1410) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1410) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1411) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1411) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1412) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1412) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1413) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1413) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1414) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1414) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1415) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1415) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1416) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1416) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1417) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1417) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1418) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1418) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1419) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1419) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1420) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1420) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1421) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1421) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1422) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1422) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1423) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1423) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1424) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1424) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1425) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1425) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1426) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1426) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1427) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1427) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1428) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1428) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1429) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1429) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1430) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1430) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1431) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1431) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1432) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1432) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1433) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1433) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1434) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1434) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1435) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1435) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1436) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1436) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1437) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1437) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1438) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1438) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1439) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1439) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1440) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1440) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1441) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1441) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1442) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1442) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1443) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1443) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1444) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1444) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1445) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1445) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1446) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1446) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1447) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1447) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1448) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1448) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1449) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1449) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1450) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1450) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1451) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1451) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1452) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1452) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1453) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1453) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1454) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1454) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1455) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1455) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1456) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1456) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1457) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1457) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1458) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1458) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1459) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1459) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1460) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1460) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1461) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1461) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1462) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1462) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1463) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1463) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1464) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1464) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1465) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1465) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1466) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1466) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1467) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1467) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1468) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1468) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1469) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1469) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1470) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1470) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1471) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1471) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1472) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1472) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1473) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1473) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1474) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1474) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1475) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1475) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1476) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1476) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1477) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1477) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1478) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1478) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1479) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1479) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1480) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1480) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1481) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1481) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1482) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1482) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1483) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1483) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1484) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1484) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1485) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1485) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1486) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1486) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1487) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1487) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1488) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1488) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1489) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1489) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1490) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1490) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1491) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1491) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1492) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1492) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1493) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1493) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1494) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1494) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1495) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1495) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1496) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1496) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1497) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1497) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1498) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1498) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1499) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1499) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1500) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1500) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1501) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1501) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1502) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1502) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1503) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1503) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1504) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1504) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1505) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1505) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1506) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1506) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1507) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1507) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1508) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1508) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1509) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1509) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1510) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1510) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1511) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1511) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1512) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1512) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1513) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1513) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1514) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1514) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1515) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1515) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1516) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1516) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1517) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1517) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1518) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1518) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1519) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1519) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1520) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1520) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1521) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1521) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1522) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1522) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1523) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1523) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1524) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1524) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1525) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1525) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1526) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1526) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1527) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1527) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1528) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1528) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1529) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1529) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1530) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1530) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1531) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1531) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1532) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1532) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1533) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1533) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1534) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1534) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1535) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1535) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1536) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1536) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1537) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1537) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1538) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1538) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1539) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1539) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1540) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1540) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1541) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1541) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1542) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1542) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1543) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1543) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1544) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1544) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1545) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1545) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1546) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1546) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1547) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1547) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1548) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1548) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1549) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1549) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1550) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1550) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1551) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1551) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1552) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1552) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1553) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1553) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1554) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1554) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1555) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1555) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1556) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1556) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1557) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1557) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1558) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1558) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1559) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1559) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1560) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1560) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1561) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1561) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1562) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1562) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1563) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1563) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1564) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1564) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1565) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1565) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1566) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1566) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1567) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1567) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1568) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1568) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1569) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1569) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1570) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1570) =    0.00
 given_sigma_ion_N2_to_N2_cm2(1571) =    0.00 ;  given_sigma_ion_N2_to_N_cm2(1571) =    0.00










 given_sigma_tot_N2_cm2(   1) =     0.80
 given_sigma_tot_N2_cm2(   2) =     0.96
 given_sigma_tot_N2_cm2(   3) =     0.97
 given_sigma_tot_N2_cm2(   4) =     0.99
 given_sigma_tot_N2_cm2(   5) =     1.01
 given_sigma_tot_N2_cm2(   6) =     1.02
 given_sigma_tot_N2_cm2(   7) =     0.91
 given_sigma_tot_N2_cm2(   8) =     0.51
 given_sigma_tot_N2_cm2(   9) =     0.08
 given_sigma_tot_N2_cm2(  10) =     0.10
 given_sigma_tot_N2_cm2(  11) =     0.11
 given_sigma_tot_N2_cm2(  12) =     0.15
 given_sigma_tot_N2_cm2(  13) =     0.15
 given_sigma_tot_N2_cm2(  14) =     0.17
 given_sigma_tot_N2_cm2(  15) =     0.17
 given_sigma_tot_N2_cm2(  16) =     0.17
 given_sigma_tot_N2_cm2(  17) =     0.18
 given_sigma_tot_N2_cm2(  18) =     0.18
 given_sigma_tot_N2_cm2(  19) =     0.19
 given_sigma_tot_N2_cm2(  20) =     0.19
 given_sigma_tot_N2_cm2(  21) =     0.20
 given_sigma_tot_N2_cm2(  22) =     0.20
 given_sigma_tot_N2_cm2(  23) =     0.21
 given_sigma_tot_N2_cm2(  24) =     0.21
 given_sigma_tot_N2_cm2(  25) =     0.23
 given_sigma_tot_N2_cm2(  26) =     0.23
 given_sigma_tot_N2_cm2(  27) =     0.24
 given_sigma_tot_N2_cm2(  28) =     0.24
 given_sigma_tot_N2_cm2(  29) =     0.24
 given_sigma_tot_N2_cm2(  30) =     0.25
 given_sigma_tot_N2_cm2(  31) =     0.26
 given_sigma_tot_N2_cm2(  32) =     0.27
 given_sigma_tot_N2_cm2(  33) =     0.28
 given_sigma_tot_N2_cm2(  34) =     0.28
 given_sigma_tot_N2_cm2(  35) =     0.29
 given_sigma_tot_N2_cm2(  36) =     0.29
 given_sigma_tot_N2_cm2(  37) =     0.30
 given_sigma_tot_N2_cm2(  38) =     0.30
 given_sigma_tot_N2_cm2(  39) =     0.31
 given_sigma_tot_N2_cm2(  40) =     0.31
 given_sigma_tot_N2_cm2(  41) =     0.32
 given_sigma_tot_N2_cm2(  42) =     0.33
 given_sigma_tot_N2_cm2(  43) =     0.33
 given_sigma_tot_N2_cm2(  44) =     0.34
 given_sigma_tot_N2_cm2(  45) =     0.35
 given_sigma_tot_N2_cm2(  46) =     0.36
 given_sigma_tot_N2_cm2(  47) =     0.37
 given_sigma_tot_N2_cm2(  48) =     0.38
 given_sigma_tot_N2_cm2(  49) =     0.38
 given_sigma_tot_N2_cm2(  50) =     0.39
 given_sigma_tot_N2_cm2(  51) =     0.39
 given_sigma_tot_N2_cm2(  52) =     0.40
 given_sigma_tot_N2_cm2(  53) =     0.40
 given_sigma_tot_N2_cm2(  54) =     0.41
 given_sigma_tot_N2_cm2(  55) =     0.41
 given_sigma_tot_N2_cm2(  56) =     0.41
 given_sigma_tot_N2_cm2(  57) =     0.42
 given_sigma_tot_N2_cm2(  58) =     0.42
 given_sigma_tot_N2_cm2(  59) =     0.43
 given_sigma_tot_N2_cm2(  60) =     0.43
 given_sigma_tot_N2_cm2(  61) =     0.44
 given_sigma_tot_N2_cm2(  62) =     0.44
 given_sigma_tot_N2_cm2(  63) =     0.45
 given_sigma_tot_N2_cm2(  64) =     0.45
 given_sigma_tot_N2_cm2(  65) =     0.46
 given_sigma_tot_N2_cm2(  66) =     0.46
 given_sigma_tot_N2_cm2(  67) =     0.47
 given_sigma_tot_N2_cm2(  68) =     0.48
 given_sigma_tot_N2_cm2(  69) =     0.50
 given_sigma_tot_N2_cm2(  70) =     0.51
 given_sigma_tot_N2_cm2(  71) =     0.52
 given_sigma_tot_N2_cm2(  72) =     0.53
 given_sigma_tot_N2_cm2(  73) =     0.53
 given_sigma_tot_N2_cm2(  74) =     0.54
 given_sigma_tot_N2_cm2(  75) =     0.56
 given_sigma_tot_N2_cm2(  76) =     0.56
 given_sigma_tot_N2_cm2(  77) =     0.57
 given_sigma_tot_N2_cm2(  78) =     0.57
 given_sigma_tot_N2_cm2(  79) =     0.57
 given_sigma_tot_N2_cm2(  80) =     0.59
 given_sigma_tot_N2_cm2(  81) =     0.60
 given_sigma_tot_N2_cm2(  82) =     0.60
 given_sigma_tot_N2_cm2(  83) =     0.61
 given_sigma_tot_N2_cm2(  84) =     0.62
 given_sigma_tot_N2_cm2(  85) =     0.62
 given_sigma_tot_N2_cm2(  86) =     0.62
 given_sigma_tot_N2_cm2(  87) =     0.63
 given_sigma_tot_N2_cm2(  88) =     0.63
 given_sigma_tot_N2_cm2(  89) =     0.64
 given_sigma_tot_N2_cm2(  90) =     0.65
 given_sigma_tot_N2_cm2(  91) =     0.65
 given_sigma_tot_N2_cm2(  92) =     0.66
 given_sigma_tot_N2_cm2(  93) =     0.66
 given_sigma_tot_N2_cm2(  94) =     0.66
 given_sigma_tot_N2_cm2(  95) =     0.68
 given_sigma_tot_N2_cm2(  96) =     0.68
 given_sigma_tot_N2_cm2(  97) =     0.69
 given_sigma_tot_N2_cm2(  98) =     0.69
 given_sigma_tot_N2_cm2(  99) =     0.70
 given_sigma_tot_N2_cm2( 100) =     0.70
 given_sigma_tot_N2_cm2( 101) =     0.71
 given_sigma_tot_N2_cm2( 102) =     0.71
 given_sigma_tot_N2_cm2( 103) =     0.72
 given_sigma_tot_N2_cm2( 104) =     0.73
 given_sigma_tot_N2_cm2( 105) =     0.73
 given_sigma_tot_N2_cm2( 106) =     0.74
 given_sigma_tot_N2_cm2( 107) =     0.74
 given_sigma_tot_N2_cm2( 108) =     0.75
 given_sigma_tot_N2_cm2( 109) =     0.76
 given_sigma_tot_N2_cm2( 110) =     0.78
 given_sigma_tot_N2_cm2( 111) =     0.78
 given_sigma_tot_N2_cm2( 112) =     0.79
 given_sigma_tot_N2_cm2( 113) =     0.80
 given_sigma_tot_N2_cm2( 114) =     0.81
 given_sigma_tot_N2_cm2( 115) =     0.82
 given_sigma_tot_N2_cm2( 116) =     0.82
 given_sigma_tot_N2_cm2( 117) =     0.83
 given_sigma_tot_N2_cm2( 118) =     0.84
 given_sigma_tot_N2_cm2( 119) =     0.85
 given_sigma_tot_N2_cm2( 120) =     0.85
 given_sigma_tot_N2_cm2( 121) =     0.86
 given_sigma_tot_N2_cm2( 122) =     0.87
 given_sigma_tot_N2_cm2( 123) =     0.87
 given_sigma_tot_N2_cm2( 124) =     0.89
 given_sigma_tot_N2_cm2( 125) =     0.89
 given_sigma_tot_N2_cm2( 126) =     0.90
 given_sigma_tot_N2_cm2( 127) =     0.91
 given_sigma_tot_N2_cm2( 128) =     0.93
 given_sigma_tot_N2_cm2( 129) =     0.94
 given_sigma_tot_N2_cm2( 130) =     0.94
 given_sigma_tot_N2_cm2( 131) =     0.95
 given_sigma_tot_N2_cm2( 132) =     0.96
 given_sigma_tot_N2_cm2( 133) =     0.98
 given_sigma_tot_N2_cm2( 134) =     0.99
 given_sigma_tot_N2_cm2( 135) =     1.00
 given_sigma_tot_N2_cm2( 136) =     1.01
 given_sigma_tot_N2_cm2( 137) =     1.02
 given_sigma_tot_N2_cm2( 138) =     1.03
 given_sigma_tot_N2_cm2( 139) =     1.04
 given_sigma_tot_N2_cm2( 140) =     1.04
 given_sigma_tot_N2_cm2( 141) =     1.05
 given_sigma_tot_N2_cm2( 142) =     1.06
 given_sigma_tot_N2_cm2( 143) =     1.07
 given_sigma_tot_N2_cm2( 144) =     1.10
 given_sigma_tot_N2_cm2( 145) =     1.11
 given_sigma_tot_N2_cm2( 146) =     1.12
 given_sigma_tot_N2_cm2( 147) =     1.12
 given_sigma_tot_N2_cm2( 148) =     1.13
 given_sigma_tot_N2_cm2( 149) =     1.15
 given_sigma_tot_N2_cm2( 150) =     1.16
 given_sigma_tot_N2_cm2( 151) =     1.16
 given_sigma_tot_N2_cm2( 152) =     1.17
 given_sigma_tot_N2_cm2( 153) =     1.19
 given_sigma_tot_N2_cm2( 154) =     1.20
 given_sigma_tot_N2_cm2( 155) =     1.20
 given_sigma_tot_N2_cm2( 156) =     1.22
 given_sigma_tot_N2_cm2( 157) =     1.23
 given_sigma_tot_N2_cm2( 158) =     1.24
 given_sigma_tot_N2_cm2( 159) =     1.24
 given_sigma_tot_N2_cm2( 160) =     1.25
 given_sigma_tot_N2_cm2( 161) =     1.26
 given_sigma_tot_N2_cm2( 162) =     1.28
 given_sigma_tot_N2_cm2( 163) =     1.28
 given_sigma_tot_N2_cm2( 164) =     1.29
 given_sigma_tot_N2_cm2( 165) =     1.31
 given_sigma_tot_N2_cm2( 166) =     1.32
 given_sigma_tot_N2_cm2( 167) =     1.34
 given_sigma_tot_N2_cm2( 168) =     1.36
 given_sigma_tot_N2_cm2( 169) =     1.39
 given_sigma_tot_N2_cm2( 170) =     1.39
 given_sigma_tot_N2_cm2( 171) =     1.39
 given_sigma_tot_N2_cm2( 172) =     1.40
 given_sigma_tot_N2_cm2( 173) =     1.41
 given_sigma_tot_N2_cm2( 174) =     1.42
 given_sigma_tot_N2_cm2( 175) =     1.44
 given_sigma_tot_N2_cm2( 176) =     1.45
 given_sigma_tot_N2_cm2( 177) =     1.49
 given_sigma_tot_N2_cm2( 178) =     1.50
 given_sigma_tot_N2_cm2( 179) =     1.51
 given_sigma_tot_N2_cm2( 180) =     1.53
 given_sigma_tot_N2_cm2( 181) =     1.54
 given_sigma_tot_N2_cm2( 182) =     1.55
 given_sigma_tot_N2_cm2( 183) =     1.59
 given_sigma_tot_N2_cm2( 184) =     1.60
 given_sigma_tot_N2_cm2( 185) =     1.62
 given_sigma_tot_N2_cm2( 186) =     1.62
 given_sigma_tot_N2_cm2( 187) =     1.64
 given_sigma_tot_N2_cm2( 188) =     1.68
 given_sigma_tot_N2_cm2( 189) =     1.72
 given_sigma_tot_N2_cm2( 190) =     1.73
 given_sigma_tot_N2_cm2( 191) =     1.73
 given_sigma_tot_N2_cm2( 192) =     1.77
 given_sigma_tot_N2_cm2( 193) =     1.78
 given_sigma_tot_N2_cm2( 194) =     1.82
 given_sigma_tot_N2_cm2( 195) =     1.84
 given_sigma_tot_N2_cm2( 196) =     1.87
 given_sigma_tot_N2_cm2( 197) =     1.95
 given_sigma_tot_N2_cm2( 198) =     1.97
 given_sigma_tot_N2_cm2( 199) =     2.00
 given_sigma_tot_N2_cm2( 200) =     2.02
 given_sigma_tot_N2_cm2( 201) =     2.06
 given_sigma_tot_N2_cm2( 202) =     2.09
 given_sigma_tot_N2_cm2( 203) =     2.11
 given_sigma_tot_N2_cm2( 204) =     2.15
 given_sigma_tot_N2_cm2( 205) =     2.27
 given_sigma_tot_N2_cm2( 206) =     2.36
 given_sigma_tot_N2_cm2( 207) =     2.37
 given_sigma_tot_N2_cm2( 208) =     2.38
 given_sigma_tot_N2_cm2( 209) =     2.39
 given_sigma_tot_N2_cm2( 210) =     2.41
 given_sigma_tot_N2_cm2( 211) =     2.42
 given_sigma_tot_N2_cm2( 212) =     2.56
 given_sigma_tot_N2_cm2( 213) =     2.61
 given_sigma_tot_N2_cm2( 214) =     2.62
 given_sigma_tot_N2_cm2( 215) =     2.68
 given_sigma_tot_N2_cm2( 216) =     2.77
 given_sigma_tot_N2_cm2( 217) =     2.82
 given_sigma_tot_N2_cm2( 218) =     2.95
 given_sigma_tot_N2_cm2( 219) =     2.98
 given_sigma_tot_N2_cm2( 220) =     3.02
 given_sigma_tot_N2_cm2( 221) =     3.14
 given_sigma_tot_N2_cm2( 222) =     3.21
 given_sigma_tot_N2_cm2( 223) =     3.22
 given_sigma_tot_N2_cm2( 224) =     3.32
 given_sigma_tot_N2_cm2( 225) =     3.35
 given_sigma_tot_N2_cm2( 226) =     3.42
 given_sigma_tot_N2_cm2( 227) =     3.46
 given_sigma_tot_N2_cm2( 228) =     3.60
 given_sigma_tot_N2_cm2( 229) =     3.63
 given_sigma_tot_N2_cm2( 230) =     3.71
 given_sigma_tot_N2_cm2( 231) =     3.73
 given_sigma_tot_N2_cm2( 232) =     3.92
 given_sigma_tot_N2_cm2( 233) =     3.97
 given_sigma_tot_N2_cm2( 234) =     3.99
 given_sigma_tot_N2_cm2( 235) =     4.11
 given_sigma_tot_N2_cm2( 236) =     4.14
 given_sigma_tot_N2_cm2( 237) =     4.16
 given_sigma_tot_N2_cm2( 238) =     4.18
 given_sigma_tot_N2_cm2( 239) =     4.22
 given_sigma_tot_N2_cm2( 240) =     4.24
 given_sigma_tot_N2_cm2( 241) =     4.30
 given_sigma_tot_N2_cm2( 242) =     4.35
 given_sigma_tot_N2_cm2( 243) =     4.40
 given_sigma_tot_N2_cm2( 244) =     4.41
 given_sigma_tot_N2_cm2( 245) =     4.48
 given_sigma_tot_N2_cm2( 246) =     4.51
 given_sigma_tot_N2_cm2( 247) =     4.52
 given_sigma_tot_N2_cm2( 248) =     4.54
 given_sigma_tot_N2_cm2( 249) =     4.64
 given_sigma_tot_N2_cm2( 250) =     4.64
 given_sigma_tot_N2_cm2( 251) =     4.69
 given_sigma_tot_N2_cm2( 252) =     4.77
 given_sigma_tot_N2_cm2( 253) =     4.79
 given_sigma_tot_N2_cm2( 254) =     4.81
 given_sigma_tot_N2_cm2( 255) =     4.84
 given_sigma_tot_N2_cm2( 256) =     4.86
 given_sigma_tot_N2_cm2( 257) =     4.89
 given_sigma_tot_N2_cm2( 258) =     4.96
 given_sigma_tot_N2_cm2( 259) =     4.98
 given_sigma_tot_N2_cm2( 260) =     5.05
 given_sigma_tot_N2_cm2( 261) =     5.08
 given_sigma_tot_N2_cm2( 262) =     5.10
 given_sigma_tot_N2_cm2( 263) =     5.13
 given_sigma_tot_N2_cm2( 264) =     5.14
 given_sigma_tot_N2_cm2( 265) =     5.16
 given_sigma_tot_N2_cm2( 266) =     5.15
 given_sigma_tot_N2_cm2( 267) =     5.27
 given_sigma_tot_N2_cm2( 268) =     5.29
 given_sigma_tot_N2_cm2( 269) =     5.36
 given_sigma_tot_N2_cm2( 270) =     5.38
 given_sigma_tot_N2_cm2( 271) =     5.41
 given_sigma_tot_N2_cm2( 272) =     5.50
 given_sigma_tot_N2_cm2( 273) =     5.55
 given_sigma_tot_N2_cm2( 274) =     5.57
 given_sigma_tot_N2_cm2( 275) =     5.59
 given_sigma_tot_N2_cm2( 276) =     5.66
 given_sigma_tot_N2_cm2( 277) =     5.68
 given_sigma_tot_N2_cm2( 278) =     5.73
 given_sigma_tot_N2_cm2( 279) =     5.83
 given_sigma_tot_N2_cm2( 280) =     5.84
 given_sigma_tot_N2_cm2( 281) =     5.94
 given_sigma_tot_N2_cm2( 282) =     5.95
 given_sigma_tot_N2_cm2( 283) =     6.01
 given_sigma_tot_N2_cm2( 284) =     6.09
 given_sigma_tot_N2_cm2( 285) =     6.20
 given_sigma_tot_N2_cm2( 286) =     6.28
 given_sigma_tot_N2_cm2( 287) =     6.36
 given_sigma_tot_N2_cm2( 288) =     6.40
 given_sigma_tot_N2_cm2( 289) =     6.49
 given_sigma_tot_N2_cm2( 290) =     6.52
 given_sigma_tot_N2_cm2( 291) =     6.57
 given_sigma_tot_N2_cm2( 292) =     6.66
 given_sigma_tot_N2_cm2( 293) =     6.69
 given_sigma_tot_N2_cm2( 294) =     6.71
 given_sigma_tot_N2_cm2( 295) =     6.78
 given_sigma_tot_N2_cm2( 296) =     6.85
 given_sigma_tot_N2_cm2( 297) =     6.96
 given_sigma_tot_N2_cm2( 298) =     6.97
 given_sigma_tot_N2_cm2( 299) =     6.98
 given_sigma_tot_N2_cm2( 300) =     7.10
 given_sigma_tot_N2_cm2( 301) =     7.17
 given_sigma_tot_N2_cm2( 302) =     7.30
 given_sigma_tot_N2_cm2( 303) =     7.38
 given_sigma_tot_N2_cm2( 304) =     7.40
 given_sigma_tot_N2_cm2( 305) =     7.41
 given_sigma_tot_N2_cm2( 306) =     7.57
 given_sigma_tot_N2_cm2( 307) =     7.58
 given_sigma_tot_N2_cm2( 308) =     7.69
 given_sigma_tot_N2_cm2( 309) =     7.77
 given_sigma_tot_N2_cm2( 310) =     7.85
 given_sigma_tot_N2_cm2( 311) =     7.96
 given_sigma_tot_N2_cm2( 312) =     7.97
 given_sigma_tot_N2_cm2( 313) =     8.01
 given_sigma_tot_N2_cm2( 314) =     8.13
 given_sigma_tot_N2_cm2( 315) =     8.17
 given_sigma_tot_N2_cm2( 316) =     8.26
 given_sigma_tot_N2_cm2( 317) =     8.28
 given_sigma_tot_N2_cm2( 318) =     8.29
 given_sigma_tot_N2_cm2( 319) =     8.45
 given_sigma_tot_N2_cm2( 320) =     8.47
 given_sigma_tot_N2_cm2( 321) =     8.49
 given_sigma_tot_N2_cm2( 322) =     8.60
 given_sigma_tot_N2_cm2( 323) =     8.68
 given_sigma_tot_N2_cm2( 324) =     8.71
 given_sigma_tot_N2_cm2( 325) =     8.76
 given_sigma_tot_N2_cm2( 326) =     8.84
 given_sigma_tot_N2_cm2( 327) =     8.92
 given_sigma_tot_N2_cm2( 328) =     9.02
 given_sigma_tot_N2_cm2( 329) =     9.07
 given_sigma_tot_N2_cm2( 330) =     9.12
 given_sigma_tot_N2_cm2( 331) =     9.16
 given_sigma_tot_N2_cm2( 332) =     9.30
 given_sigma_tot_N2_cm2( 333) =     9.38
 given_sigma_tot_N2_cm2( 334) =     9.49
 given_sigma_tot_N2_cm2( 335) =     9.50
 given_sigma_tot_N2_cm2( 336) =     9.54
 given_sigma_tot_N2_cm2( 337) =     9.60
 given_sigma_tot_N2_cm2( 338) =     9.68
 given_sigma_tot_N2_cm2( 339) =     9.73
 given_sigma_tot_N2_cm2( 340) =     9.80
 given_sigma_tot_N2_cm2( 341) =     9.84
 given_sigma_tot_N2_cm2( 342) =     9.85
 given_sigma_tot_N2_cm2( 343) =     9.88
 given_sigma_tot_N2_cm2( 344) =     9.90
 given_sigma_tot_N2_cm2( 345) =     9.93
 given_sigma_tot_N2_cm2( 346) =     9.98
 given_sigma_tot_N2_cm2( 347) =    10.02
 given_sigma_tot_N2_cm2( 348) =    10.06
 given_sigma_tot_N2_cm2( 349) =    10.08
 given_sigma_tot_N2_cm2( 350) =    10.09
 given_sigma_tot_N2_cm2( 351) =    10.14
 given_sigma_tot_N2_cm2( 352) =    10.16
 given_sigma_tot_N2_cm2( 353) =    10.18
 given_sigma_tot_N2_cm2( 354) =    10.21
 given_sigma_tot_N2_cm2( 355) =    10.21
 given_sigma_tot_N2_cm2( 356) =    10.22
 given_sigma_tot_N2_cm2( 357) =    10.23
 given_sigma_tot_N2_cm2( 358) =    10.25
 given_sigma_tot_N2_cm2( 359) =    10.27
 given_sigma_tot_N2_cm2( 360) =    10.28
 given_sigma_tot_N2_cm2( 361) =    10.31
 given_sigma_tot_N2_cm2( 362) =    10.35
 given_sigma_tot_N2_cm2( 363) =    10.38
 given_sigma_tot_N2_cm2( 364) =    10.40
 given_sigma_tot_N2_cm2( 365) =    10.40
 given_sigma_tot_N2_cm2( 366) =    10.49
 given_sigma_tot_N2_cm2( 367) =    10.50
 given_sigma_tot_N2_cm2( 368) =    10.51
 given_sigma_tot_N2_cm2( 369) =    10.54
 given_sigma_tot_N2_cm2( 370) =    10.56
 given_sigma_tot_N2_cm2( 371) =    10.59
 given_sigma_tot_N2_cm2( 372) =    10.61
 given_sigma_tot_N2_cm2( 373) =    10.62
 given_sigma_tot_N2_cm2( 374) =    10.63
 given_sigma_tot_N2_cm2( 375) =    10.65
 given_sigma_tot_N2_cm2( 376) =    10.67
 given_sigma_tot_N2_cm2( 377) =    10.68
 given_sigma_tot_N2_cm2( 378) =    10.69
 given_sigma_tot_N2_cm2( 379) =    10.73
 given_sigma_tot_N2_cm2( 380) =    10.78
 given_sigma_tot_N2_cm2( 381) =    10.81
 given_sigma_tot_N2_cm2( 382) =    10.82
 given_sigma_tot_N2_cm2( 383) =    10.90
 given_sigma_tot_N2_cm2( 384) =    10.92
 given_sigma_tot_N2_cm2( 385) =    10.95
 given_sigma_tot_N2_cm2( 386) =    10.95
 given_sigma_tot_N2_cm2( 387) =    11.04
 given_sigma_tot_N2_cm2( 388) =    11.07
 given_sigma_tot_N2_cm2( 389) =    11.10
 given_sigma_tot_N2_cm2( 390) =    11.12
 given_sigma_tot_N2_cm2( 391) =    11.16
 given_sigma_tot_N2_cm2( 392) =    11.17
 given_sigma_tot_N2_cm2( 393) =    11.19
 given_sigma_tot_N2_cm2( 394) =    11.27
 given_sigma_tot_N2_cm2( 395) =    11.30
 given_sigma_tot_N2_cm2( 396) =    11.32
 given_sigma_tot_N2_cm2( 397) =    11.48
 given_sigma_tot_N2_cm2( 398) =    11.50
 given_sigma_tot_N2_cm2( 399) =    11.54
 given_sigma_tot_N2_cm2( 400) =    11.67
 given_sigma_tot_N2_cm2( 401) =    11.70
 given_sigma_tot_N2_cm2( 402) =    12.02
 given_sigma_tot_N2_cm2( 403) =    12.37
 given_sigma_tot_N2_cm2( 404) =    12.45
 given_sigma_tot_N2_cm2( 405) =    12.46
 given_sigma_tot_N2_cm2( 406) =    12.66
 given_sigma_tot_N2_cm2( 407) =    12.72
 given_sigma_tot_N2_cm2( 408) =    12.73
 given_sigma_tot_N2_cm2( 409) =    12.77
 given_sigma_tot_N2_cm2( 410) =    13.13
 given_sigma_tot_N2_cm2( 411) =    13.55
 given_sigma_tot_N2_cm2( 412) =    13.95
 given_sigma_tot_N2_cm2( 413) =    13.98
 given_sigma_tot_N2_cm2( 414) =    14.38
 given_sigma_tot_N2_cm2( 415) =    14.82
 given_sigma_tot_N2_cm2( 416) =    14.83
 given_sigma_tot_N2_cm2( 417) =    14.89
 given_sigma_tot_N2_cm2( 418) =    15.04
 given_sigma_tot_N2_cm2( 419) =    15.27
 given_sigma_tot_N2_cm2( 420) =    15.28
 given_sigma_tot_N2_cm2( 421) =    15.63
 given_sigma_tot_N2_cm2( 422) =    15.82
 given_sigma_tot_N2_cm2( 423) =    16.18
 given_sigma_tot_N2_cm2( 424) =    16.25
 given_sigma_tot_N2_cm2( 425) =    16.58
 given_sigma_tot_N2_cm2( 426) =    16.61
 given_sigma_tot_N2_cm2( 427) =    16.91
 given_sigma_tot_N2_cm2( 428) =    17.08
 given_sigma_tot_N2_cm2( 429) =    17.53
 given_sigma_tot_N2_cm2( 430) =    18.03
 given_sigma_tot_N2_cm2( 431) =    19.05
 given_sigma_tot_N2_cm2( 432) =    20.05
 given_sigma_tot_N2_cm2( 433) =    20.07
 given_sigma_tot_N2_cm2( 434) =    20.18
 given_sigma_tot_N2_cm2( 435) =    20.23
 given_sigma_tot_N2_cm2( 436) =    20.25
 given_sigma_tot_N2_cm2( 437) =    20.38
 given_sigma_tot_N2_cm2( 438) =    20.55
 given_sigma_tot_N2_cm2( 439) =    20.64
 given_sigma_tot_N2_cm2( 440) =    20.74
 given_sigma_tot_N2_cm2( 441) =    20.83
 given_sigma_tot_N2_cm2( 442) =    20.93
 given_sigma_tot_N2_cm2( 443) =    21.02
 given_sigma_tot_N2_cm2( 444) =    21.10
 given_sigma_tot_N2_cm2( 445) =    21.19
 given_sigma_tot_N2_cm2( 446) =    21.27
 given_sigma_tot_N2_cm2( 447) =    21.35
 given_sigma_tot_N2_cm2( 448) =    21.44
 given_sigma_tot_N2_cm2( 449) =    21.52
 given_sigma_tot_N2_cm2( 450) =    21.60
 given_sigma_tot_N2_cm2( 451) =    21.62
 given_sigma_tot_N2_cm2( 452) =    21.66
 given_sigma_tot_N2_cm2( 453) =    21.68
 given_sigma_tot_N2_cm2( 454) =    21.77
 given_sigma_tot_N2_cm2( 455) =    21.85
 given_sigma_tot_N2_cm2( 456) =    21.93
 given_sigma_tot_N2_cm2( 457) =    22.01
 given_sigma_tot_N2_cm2( 458) =    22.09
 given_sigma_tot_N2_cm2( 459) =    22.17
 given_sigma_tot_N2_cm2( 460) =    22.25
 given_sigma_tot_N2_cm2( 461) =    22.32
 given_sigma_tot_N2_cm2( 462) =    22.40
 given_sigma_tot_N2_cm2( 463) =    22.48
 given_sigma_tot_N2_cm2( 464) =    22.56
 given_sigma_tot_N2_cm2( 465) =    22.64
 given_sigma_tot_N2_cm2( 466) =    22.66
 given_sigma_tot_N2_cm2( 467) =    22.69
 given_sigma_tot_N2_cm2( 468) =    22.73
 given_sigma_tot_N2_cm2( 469) =    22.78
 given_sigma_tot_N2_cm2( 470) =    22.82
 given_sigma_tot_N2_cm2( 471) =    22.87
 given_sigma_tot_N2_cm2( 472) =    22.92
 given_sigma_tot_N2_cm2( 473) =    22.92
 given_sigma_tot_N2_cm2( 474) =    22.95
 given_sigma_tot_N2_cm2( 475) =    22.96
 given_sigma_tot_N2_cm2( 476) =    23.01
 given_sigma_tot_N2_cm2( 477) =    23.05
 given_sigma_tot_N2_cm2( 478) =    23.10
 given_sigma_tot_N2_cm2( 479) =    23.10
 given_sigma_tot_N2_cm2( 480) =    23.10
 given_sigma_tot_N2_cm2( 481) =    23.10
 given_sigma_tot_N2_cm2( 482) =    23.10
 given_sigma_tot_N2_cm2( 483) =    23.10
 given_sigma_tot_N2_cm2( 484) =    23.10
 given_sigma_tot_N2_cm2( 485) =    23.10
 given_sigma_tot_N2_cm2( 486) =    23.10
 given_sigma_tot_N2_cm2( 487) =    23.10
 given_sigma_tot_N2_cm2( 488) =    23.10
 given_sigma_tot_N2_cm2( 489) =    23.10
 given_sigma_tot_N2_cm2( 490) =    23.10
 given_sigma_tot_N2_cm2( 491) =    23.10
 given_sigma_tot_N2_cm2( 492) =    23.10
 given_sigma_tot_N2_cm2( 493) =    23.10
 given_sigma_tot_N2_cm2( 494) =    23.10
 given_sigma_tot_N2_cm2( 495) =    23.10
 given_sigma_tot_N2_cm2( 496) =    23.10
 given_sigma_tot_N2_cm2( 497) =    23.10
 given_sigma_tot_N2_cm2( 498) =    23.10
 given_sigma_tot_N2_cm2( 499) =    23.10
 given_sigma_tot_N2_cm2( 500) =    23.09
 given_sigma_tot_N2_cm2( 501) =    23.09
 given_sigma_tot_N2_cm2( 502) =    23.09
 given_sigma_tot_N2_cm2( 503) =    23.08
 given_sigma_tot_N2_cm2( 504) =    23.08
 given_sigma_tot_N2_cm2( 505) =    23.08
 given_sigma_tot_N2_cm2( 506) =    23.08
 given_sigma_tot_N2_cm2( 507) =    23.07
 given_sigma_tot_N2_cm2( 508) =    23.07
 given_sigma_tot_N2_cm2( 509) =    23.07
 given_sigma_tot_N2_cm2( 510) =    23.07
 given_sigma_tot_N2_cm2( 511) =    23.06
 given_sigma_tot_N2_cm2( 512) =    23.06
 given_sigma_tot_N2_cm2( 513) =    23.05
 given_sigma_tot_N2_cm2( 514) =    23.05
 given_sigma_tot_N2_cm2( 515) =    23.07
 given_sigma_tot_N2_cm2( 516) =    23.09
 given_sigma_tot_N2_cm2( 517) =    23.12
 given_sigma_tot_N2_cm2( 518) =    23.14
 given_sigma_tot_N2_cm2( 519) =    23.16
 given_sigma_tot_N2_cm2( 520) =    23.16
 given_sigma_tot_N2_cm2( 521) =    23.18
 given_sigma_tot_N2_cm2( 522) =    23.20
 given_sigma_tot_N2_cm2( 523) =    23.23
 given_sigma_tot_N2_cm2( 524) =    23.25
 given_sigma_tot_N2_cm2( 525) =    23.27
 given_sigma_tot_N2_cm2( 526) =    23.27
 given_sigma_tot_N2_cm2( 527) =    23.30
 given_sigma_tot_N2_cm2( 528) =    23.34
 given_sigma_tot_N2_cm2( 529) =    23.37
 given_sigma_tot_N2_cm2( 530) =    23.38
 given_sigma_tot_N2_cm2( 531) =    23.40
 given_sigma_tot_N2_cm2( 532) =    23.40
 given_sigma_tot_N2_cm2( 533) =    23.44
 given_sigma_tot_N2_cm2( 534) =    23.46
 given_sigma_tot_N2_cm2( 535) =    23.47
 given_sigma_tot_N2_cm2( 536) =    23.48
 given_sigma_tot_N2_cm2( 537) =    23.50
 given_sigma_tot_N2_cm2( 538) =    23.53
 given_sigma_tot_N2_cm2( 539) =    23.55
 given_sigma_tot_N2_cm2( 540) =    23.57
 given_sigma_tot_N2_cm2( 541) =    23.60
 given_sigma_tot_N2_cm2( 542) =    23.59
 given_sigma_tot_N2_cm2( 543) =    23.59
 given_sigma_tot_N2_cm2( 544) =    23.58
 given_sigma_tot_N2_cm2( 545) =    23.58
 given_sigma_tot_N2_cm2( 546) =    23.58
 given_sigma_tot_N2_cm2( 547) =    23.57
 given_sigma_tot_N2_cm2( 548) =    23.57
 given_sigma_tot_N2_cm2( 549) =    23.57
 given_sigma_tot_N2_cm2( 550) =    23.56
 given_sigma_tot_N2_cm2( 551) =    23.55
 given_sigma_tot_N2_cm2( 552) =    23.55
 given_sigma_tot_N2_cm2( 553) =    23.54
 given_sigma_tot_N2_cm2( 554) =    23.53
 given_sigma_tot_N2_cm2( 555) =    23.53
 given_sigma_tot_N2_cm2( 556) =    23.52
 given_sigma_tot_N2_cm2( 557) =    23.51
 given_sigma_tot_N2_cm2( 558) =    23.50
 given_sigma_tot_N2_cm2( 559) =    23.50
 given_sigma_tot_N2_cm2( 560) =    23.50
 given_sigma_tot_N2_cm2( 561) =    23.50
 given_sigma_tot_N2_cm2( 562) =    23.50
 given_sigma_tot_N2_cm2( 563) =    23.50
 given_sigma_tot_N2_cm2( 564) =    23.50
 given_sigma_tot_N2_cm2( 565) =    23.50
 given_sigma_tot_N2_cm2( 566) =    23.50
 given_sigma_tot_N2_cm2( 567) =    23.50
 given_sigma_tot_N2_cm2( 568) =    23.50
 given_sigma_tot_N2_cm2( 569) =    23.50
 given_sigma_tot_N2_cm2( 570) =    23.50
 given_sigma_tot_N2_cm2( 571) =    23.50
 given_sigma_tot_N2_cm2( 572) =    23.50
 given_sigma_tot_N2_cm2( 573) =    23.52
 given_sigma_tot_N2_cm2( 574) =    23.52
 given_sigma_tot_N2_cm2( 575) =    23.53
 given_sigma_tot_N2_cm2( 576) =    23.55
 given_sigma_tot_N2_cm2( 577) =    23.56
 given_sigma_tot_N2_cm2( 578) =    23.58
 given_sigma_tot_N2_cm2( 579) =    23.73
 given_sigma_tot_N2_cm2( 580) =    23.83
 given_sigma_tot_N2_cm2( 581) =    24.25
 given_sigma_tot_N2_cm2( 582) =    24.58
 given_sigma_tot_N2_cm2( 583) =    24.61
 given_sigma_tot_N2_cm2( 584) =    24.63
 given_sigma_tot_N2_cm2( 585) =    24.86
 given_sigma_tot_N2_cm2( 586) =    25.07
 given_sigma_tot_N2_cm2( 587) =    25.23
 given_sigma_tot_N2_cm2( 588) =    25.30
 given_sigma_tot_N2_cm2( 589) =    25.13
 given_sigma_tot_N2_cm2( 590) =    25.02
 given_sigma_tot_N2_cm2( 591) =    24.77
 given_sigma_tot_N2_cm2( 592) =    24.70
 given_sigma_tot_N2_cm2( 593) =    24.52
 given_sigma_tot_N2_cm2( 594) =    24.13
 given_sigma_tot_N2_cm2( 595) =    24.11
 given_sigma_tot_N2_cm2( 596) =    23.97
 given_sigma_tot_N2_cm2( 597) =    23.58
 given_sigma_tot_N2_cm2( 598) =    23.40
 given_sigma_tot_N2_cm2( 599) =    23.15
 given_sigma_tot_N2_cm2( 600) =    22.95
 given_sigma_tot_N2_cm2( 601) =    22.64
 given_sigma_tot_N2_cm2( 602) =    22.50
 given_sigma_tot_N2_cm2( 603) =    22.48
 given_sigma_tot_N2_cm2( 604) =    22.45
 given_sigma_tot_N2_cm2( 605) =    22.40
 given_sigma_tot_N2_cm2( 606) =    22.40
 given_sigma_tot_N2_cm2( 607) =    22.40
 given_sigma_tot_N2_cm2( 608) =    22.40
 given_sigma_tot_N2_cm2( 609) =    22.40
 given_sigma_tot_N2_cm2( 610) =    22.40
 given_sigma_tot_N2_cm2( 611) =    22.40
 given_sigma_tot_N2_cm2( 612) =    22.40
 given_sigma_tot_N2_cm2( 613) =    22.44
 given_sigma_tot_N2_cm2( 614) =    22.47
 given_sigma_tot_N2_cm2( 615) =    22.49
 given_sigma_tot_N2_cm2( 616) =    22.52
 given_sigma_tot_N2_cm2( 617) =    22.57
 given_sigma_tot_N2_cm2( 618) =    22.58
 given_sigma_tot_N2_cm2( 619) =    22.67
 given_sigma_tot_N2_cm2( 620) =    22.76
 given_sigma_tot_N2_cm2( 621) =    22.76
 given_sigma_tot_N2_cm2( 622) =    22.78
 given_sigma_tot_N2_cm2( 623) =    22.79
 given_sigma_tot_N2_cm2( 624) =    22.80
 given_sigma_tot_N2_cm2( 625) =    22.82
 given_sigma_tot_N2_cm2( 626) =    22.83
 given_sigma_tot_N2_cm2( 627) =    22.86
 given_sigma_tot_N2_cm2( 628) =    22.87
 given_sigma_tot_N2_cm2( 629) =    22.88
 given_sigma_tot_N2_cm2( 630) =    22.89
 given_sigma_tot_N2_cm2( 631) =    22.92
 given_sigma_tot_N2_cm2( 632) =    22.95
 given_sigma_tot_N2_cm2( 633) =    22.96
 given_sigma_tot_N2_cm2( 634) =    22.98
 given_sigma_tot_N2_cm2( 635) =    22.99
 given_sigma_tot_N2_cm2( 636) =    23.00
 given_sigma_tot_N2_cm2( 637) =    23.01
 given_sigma_tot_N2_cm2( 638) =    23.03
 given_sigma_tot_N2_cm2( 639) =    23.04
 given_sigma_tot_N2_cm2( 640) =    23.05
 given_sigma_tot_N2_cm2( 641) =    23.06
 given_sigma_tot_N2_cm2( 642) =    23.07
 given_sigma_tot_N2_cm2( 643) =    23.09
 given_sigma_tot_N2_cm2( 644) =    23.10
 given_sigma_tot_N2_cm2( 645) =    23.11
 given_sigma_tot_N2_cm2( 646) =    23.13
 given_sigma_tot_N2_cm2( 647) =    23.15
 given_sigma_tot_N2_cm2( 648) =    23.18
 given_sigma_tot_N2_cm2( 649) =    23.21
 given_sigma_tot_N2_cm2( 650) =    23.23
 given_sigma_tot_N2_cm2( 651) =    23.24
 given_sigma_tot_N2_cm2( 652) =    23.27
 given_sigma_tot_N2_cm2( 653) =    23.28
 given_sigma_tot_N2_cm2( 654) =    23.30
 given_sigma_tot_N2_cm2( 655) =    23.32
 given_sigma_tot_N2_cm2( 656) =    23.35
 given_sigma_tot_N2_cm2( 657) =    23.37
 given_sigma_tot_N2_cm2( 658) =    23.37
 given_sigma_tot_N2_cm2( 659) =    23.38
 given_sigma_tot_N2_cm2( 660) =    23.39
 given_sigma_tot_N2_cm2( 661) =    23.41
 given_sigma_tot_N2_cm2( 662) =    23.42
 given_sigma_tot_N2_cm2( 663) =    23.44
 given_sigma_tot_N2_cm2( 664) =    23.46
 given_sigma_tot_N2_cm2( 665) =    23.49
 given_sigma_tot_N2_cm2( 666) =    23.50
 given_sigma_tot_N2_cm2( 667) =    23.51
 given_sigma_tot_N2_cm2( 668) =    23.52
 given_sigma_tot_N2_cm2( 669) =    23.54
 given_sigma_tot_N2_cm2( 670) =    23.55
 given_sigma_tot_N2_cm2( 671) =    23.56
 given_sigma_tot_N2_cm2( 672) =    23.58
 given_sigma_tot_N2_cm2( 673) =    23.58
 given_sigma_tot_N2_cm2( 674) =    23.60
 given_sigma_tot_N2_cm2( 675) =    23.62
 given_sigma_tot_N2_cm2( 676) =    23.62
 given_sigma_tot_N2_cm2( 677) =    23.63
 given_sigma_tot_N2_cm2( 678) =    23.66
 given_sigma_tot_N2_cm2( 679) =    23.67
 given_sigma_tot_N2_cm2( 680) =    23.69
 given_sigma_tot_N2_cm2( 681) =    23.71
 given_sigma_tot_N2_cm2( 682) =    23.72
 given_sigma_tot_N2_cm2( 683) =    23.73
 given_sigma_tot_N2_cm2( 684) =    23.75
 given_sigma_tot_N2_cm2( 685) =    23.78
 given_sigma_tot_N2_cm2( 686) =    23.78
 given_sigma_tot_N2_cm2( 687) =    23.79
 given_sigma_tot_N2_cm2( 688) =    23.81
 given_sigma_tot_N2_cm2( 689) =    23.83
 given_sigma_tot_N2_cm2( 690) =    23.85
 given_sigma_tot_N2_cm2( 691) =    23.85
 given_sigma_tot_N2_cm2( 692) =    23.86
 given_sigma_tot_N2_cm2( 693) =    23.88
 given_sigma_tot_N2_cm2( 694) =    23.89
 given_sigma_tot_N2_cm2( 695) =    23.92
 given_sigma_tot_N2_cm2( 696) =    23.93
 given_sigma_tot_N2_cm2( 697) =    23.95
 given_sigma_tot_N2_cm2( 698) =    23.96
 given_sigma_tot_N2_cm2( 699) =    23.98
 given_sigma_tot_N2_cm2( 700) =    24.00
 given_sigma_tot_N2_cm2( 701) =    24.00
 given_sigma_tot_N2_cm2( 702) =    24.03
 given_sigma_tot_N2_cm2( 703) =    24.05
 given_sigma_tot_N2_cm2( 704) =    24.08
 given_sigma_tot_N2_cm2( 705) =    24.10
 given_sigma_tot_N2_cm2( 706) =    24.13
 given_sigma_tot_N2_cm2( 707) =    24.13
 given_sigma_tot_N2_cm2( 708) =    24.15
 given_sigma_tot_N2_cm2( 709) =    24.18
 given_sigma_tot_N2_cm2( 710) =    24.20
 given_sigma_tot_N2_cm2( 711) =    24.30
 given_sigma_tot_N2_cm2( 712) =    25.22
 given_sigma_tot_N2_cm2( 713) =    25.74
 given_sigma_tot_N2_cm2( 714) =    26.40
 given_sigma_tot_N2_cm2( 715) =    26.69
 given_sigma_tot_N2_cm2( 716) =    26.94
 given_sigma_tot_N2_cm2( 717) =    27.10
 given_sigma_tot_N2_cm2( 718) =    25.60
 given_sigma_tot_N2_cm2( 719) =    27.47
 given_sigma_tot_N2_cm2( 720) =    27.42
 given_sigma_tot_N2_cm2( 721) =    27.40
 given_sigma_tot_N2_cm2( 722) =    30.00
 given_sigma_tot_N2_cm2( 723) =    40.40
 given_sigma_tot_N2_cm2( 724) =    22.40
 given_sigma_tot_N2_cm2( 725) =    48.40
 given_sigma_tot_N2_cm2( 726) =    32.32
 given_sigma_tot_N2_cm2( 727) =    21.60
 given_sigma_tot_N2_cm2( 728) =    22.00
 given_sigma_tot_N2_cm2( 729) =    28.71
 given_sigma_tot_N2_cm2( 730) =    34.30
 given_sigma_tot_N2_cm2( 731) =    28.40
 given_sigma_tot_N2_cm2( 732) =    23.00
 given_sigma_tot_N2_cm2( 733) =    23.00
 given_sigma_tot_N2_cm2( 734) =    23.00
 given_sigma_tot_N2_cm2( 735) =    30.41
 given_sigma_tot_N2_cm2( 736) =    45.24
 given_sigma_tot_N2_cm2( 737) =    67.49
 given_sigma_tot_N2_cm2( 738) =    74.90
 given_sigma_tot_N2_cm2( 739) =    65.77
 given_sigma_tot_N2_cm2( 740) =    60.30
 given_sigma_tot_N2_cm2( 741) =    56.65
 given_sigma_tot_N2_cm2( 742) =    49.35
 given_sigma_tot_N2_cm2( 743) =    42.05
 given_sigma_tot_N2_cm2( 744) =    32.93
 given_sigma_tot_N2_cm2( 745) =    25.62
 given_sigma_tot_N2_cm2( 746) =    25.42
 given_sigma_tot_N2_cm2( 747) =    28.13
 given_sigma_tot_N2_cm2( 748) =    29.21
 given_sigma_tot_N2_cm2( 749) =    30.29
 given_sigma_tot_N2_cm2( 750) =    34.08
 given_sigma_tot_N2_cm2( 751) =    35.70
 given_sigma_tot_N2_cm2( 752) =    34.67
 given_sigma_tot_N2_cm2( 753) =    33.14
 given_sigma_tot_N2_cm2( 754) =    31.60
 given_sigma_tot_N2_cm2( 755) =    25.30
 given_sigma_tot_N2_cm2( 756) =    32.80
 given_sigma_tot_N2_cm2( 757) =    47.80
 given_sigma_tot_N2_cm2( 758) =    55.30
 given_sigma_tot_N2_cm2( 759) =    77.80
 given_sigma_tot_N2_cm2( 760) =    65.57
 given_sigma_tot_N2_cm2( 761) =    45.18
 given_sigma_tot_N2_cm2( 762) =    37.03
 given_sigma_tot_N2_cm2( 763) =    24.80
 given_sigma_tot_N2_cm2( 764) =    24.43
 given_sigma_tot_N2_cm2( 765) =    24.36
 given_sigma_tot_N2_cm2( 766) =    24.07
 given_sigma_tot_N2_cm2( 767) =    23.92
 given_sigma_tot_N2_cm2( 768) =    23.70
 given_sigma_tot_N2_cm2( 769) =    23.70
 given_sigma_tot_N2_cm2( 770) =    23.70
 given_sigma_tot_N2_cm2( 771) =    23.70
 given_sigma_tot_N2_cm2( 772) =    23.70
 given_sigma_tot_N2_cm2( 773) =    23.70
 given_sigma_tot_N2_cm2( 774) =    23.69
 given_sigma_tot_N2_cm2( 775) =    23.68
 given_sigma_tot_N2_cm2( 776) =    23.67
 given_sigma_tot_N2_cm2( 777) =    23.65
 given_sigma_tot_N2_cm2( 778) =    23.62
 given_sigma_tot_N2_cm2( 779) =    23.60
 given_sigma_tot_N2_cm2( 780) =    23.60
 given_sigma_tot_N2_cm2( 781) =    23.60
 given_sigma_tot_N2_cm2( 782) =    23.33
 given_sigma_tot_N2_cm2( 783) =    23.16
 given_sigma_tot_N2_cm2( 784) =    23.11
 given_sigma_tot_N2_cm2( 785) =    23.05
 given_sigma_tot_N2_cm2( 786) =    23.00
 given_sigma_tot_N2_cm2( 787) =    23.00
 given_sigma_tot_N2_cm2( 788) =    23.80
 given_sigma_tot_N2_cm2( 789) =    24.07
 given_sigma_tot_N2_cm2( 790) =    25.13
 given_sigma_tot_N2_cm2( 791) =    27.91
 given_sigma_tot_N2_cm2( 792) =    30.69
 given_sigma_tot_N2_cm2( 793) =    38.10
 given_sigma_tot_N2_cm2( 794) =    39.96
 given_sigma_tot_N2_cm2( 795) =    42.74
 given_sigma_tot_N2_cm2( 796) =    48.30
 given_sigma_tot_N2_cm2( 797) =    76.10
 given_sigma_tot_N2_cm2( 798) =    61.30
 given_sigma_tot_N2_cm2( 799) =    52.05
 given_sigma_tot_N2_cm2( 800) =    42.80
 given_sigma_tot_N2_cm2( 801) =    37.25
 given_sigma_tot_N2_cm2( 802) =    33.55
 given_sigma_tot_N2_cm2( 803) =    29.85
 given_sigma_tot_N2_cm2( 804) =    24.30
 given_sigma_tot_N2_cm2( 805) =    23.83
 given_sigma_tot_N2_cm2( 806) =    22.89
 given_sigma_tot_N2_cm2( 807) =    22.10
 given_sigma_tot_N2_cm2( 808) =    22.10
 given_sigma_tot_N2_cm2( 809) =    22.10
 given_sigma_tot_N2_cm2( 810) =    23.80
 given_sigma_tot_N2_cm2( 811) =    22.10
 given_sigma_tot_N2_cm2( 812) =    22.59
 given_sigma_tot_N2_cm2( 813) =    24.53
 given_sigma_tot_N2_cm2( 814) =    25.50
 given_sigma_tot_N2_cm2( 815) =    26.42
 given_sigma_tot_N2_cm2( 816) =    28.26
 given_sigma_tot_N2_cm2( 817) =    29.18
 given_sigma_tot_N2_cm2( 818) =    30.10
 given_sigma_tot_N2_cm2( 819) =    28.30
 given_sigma_tot_N2_cm2( 820) =    26.54
 given_sigma_tot_N2_cm2( 821) =    25.36
 given_sigma_tot_N2_cm2( 822) =    23.40
 given_sigma_tot_N2_cm2( 823) =    25.80
 given_sigma_tot_N2_cm2( 824) =    25.79
 given_sigma_tot_N2_cm2( 825) =    25.79
 given_sigma_tot_N2_cm2( 826) =    25.78
 given_sigma_tot_N2_cm2( 827) =    25.78
 given_sigma_tot_N2_cm2( 828) =    25.77
 given_sigma_tot_N2_cm2( 829) =    25.77
 given_sigma_tot_N2_cm2( 830) =    25.76
 given_sigma_tot_N2_cm2( 831) =    25.75
 given_sigma_tot_N2_cm2( 832) =    25.75
 given_sigma_tot_N2_cm2( 833) =    25.74
 given_sigma_tot_N2_cm2( 834) =    25.74
 given_sigma_tot_N2_cm2( 835) =    25.73
 given_sigma_tot_N2_cm2( 836) =    25.73
 given_sigma_tot_N2_cm2( 837) =    25.72
 given_sigma_tot_N2_cm2( 838) =    25.71
 given_sigma_tot_N2_cm2( 839) =    25.71
 given_sigma_tot_N2_cm2( 840) =    25.71
 given_sigma_tot_N2_cm2( 841) =    25.70
 given_sigma_tot_N2_cm2( 842) =    25.70
 given_sigma_tot_N2_cm2( 843) =    20.57
 given_sigma_tot_N2_cm2( 844) =    18.55
 given_sigma_tot_N2_cm2( 845) =    14.50
 given_sigma_tot_N2_cm2( 846) =    14.94
 given_sigma_tot_N2_cm2( 847) =    16.26
 given_sigma_tot_N2_cm2( 848) =    18.46
 given_sigma_tot_N2_cm2( 849) =    20.66
 given_sigma_tot_N2_cm2( 850) =    21.54
 given_sigma_tot_N2_cm2( 851) =    23.30
 given_sigma_tot_N2_cm2( 852) =    25.06
 given_sigma_tot_N2_cm2( 853) =    23.20
 given_sigma_tot_N2_cm2( 854) =    27.40
 given_sigma_tot_N2_cm2( 855) =    27.05
 given_sigma_tot_N2_cm2( 856) =    26.71
 given_sigma_tot_N2_cm2( 857) =    25.66
 given_sigma_tot_N2_cm2( 858) =    24.97
 given_sigma_tot_N2_cm2( 859) =    24.10
 given_sigma_tot_N2_cm2( 860) =    26.80
 given_sigma_tot_N2_cm2( 861) =    27.48
 given_sigma_tot_N2_cm2( 862) =    30.19
 given_sigma_tot_N2_cm2( 863) =    31.54
 given_sigma_tot_N2_cm2( 864) =    34.92
 given_sigma_tot_N2_cm2( 865) =    37.62
 given_sigma_tot_N2_cm2( 866) =    69.30
 given_sigma_tot_N2_cm2( 867) =    59.51
 given_sigma_tot_N2_cm2( 868) =    53.63
 given_sigma_tot_N2_cm2( 869) =    41.88
 given_sigma_tot_N2_cm2( 870) =    39.92
 given_sigma_tot_N2_cm2( 871) =    30.60
 given_sigma_tot_N2_cm2( 872) =    30.06
 given_sigma_tot_N2_cm2( 873) =    29.24
 given_sigma_tot_N2_cm2( 874) =    28.16
 given_sigma_tot_N2_cm2( 875) =    26.53
 given_sigma_tot_N2_cm2( 876) =    25.71
 given_sigma_tot_N2_cm2( 877) =    25.17
 given_sigma_tot_N2_cm2( 878) =    24.17
 given_sigma_tot_N2_cm2( 879) =    23.63
 given_sigma_tot_N2_cm2( 880) =    22.90
 given_sigma_tot_N2_cm2( 881) =    22.94
 given_sigma_tot_N2_cm2( 882) =    23.00
 given_sigma_tot_N2_cm2( 883) =    23.05
 given_sigma_tot_N2_cm2( 884) =    23.08
 given_sigma_tot_N2_cm2( 885) =    23.16
 given_sigma_tot_N2_cm2( 886) =    23.20
 given_sigma_tot_N2_cm2( 887) =    22.68
 given_sigma_tot_N2_cm2( 888) =    22.42
 given_sigma_tot_N2_cm2( 889) =    22.24
 given_sigma_tot_N2_cm2( 890) =    22.07
 given_sigma_tot_N2_cm2( 891) =    21.81
 given_sigma_tot_N2_cm2( 892) =    21.37
 given_sigma_tot_N2_cm2( 893) =    21.11
 given_sigma_tot_N2_cm2( 894) =    20.50
 given_sigma_tot_N2_cm2( 895) =    30.50
 given_sigma_tot_N2_cm2( 896) =    25.00
 given_sigma_tot_N2_cm2( 897) =    23.65
 given_sigma_tot_N2_cm2( 898) =    20.95
 given_sigma_tot_N2_cm2( 899) =    31.80
 given_sigma_tot_N2_cm2( 900) =    26.40
 given_sigma_tot_N2_cm2( 901) =    24.04
 given_sigma_tot_N2_cm2( 902) =    20.50
 given_sigma_tot_N2_cm2( 903) =    18.20
 given_sigma_tot_N2_cm2( 904) =    21.84
 given_sigma_tot_N2_cm2( 905) =    27.30
 given_sigma_tot_N2_cm2( 906) =    26.56
 given_sigma_tot_N2_cm2( 907) =    23.60
 given_sigma_tot_N2_cm2( 908) =    22.34
 given_sigma_tot_N2_cm2( 909) =    17.30
 given_sigma_tot_N2_cm2( 910) =    17.48
 given_sigma_tot_N2_cm2( 911) =    18.20
 given_sigma_tot_N2_cm2( 912) =    17.20
 given_sigma_tot_N2_cm2( 913) =    16.20
 given_sigma_tot_N2_cm2( 914) =    13.20
 given_sigma_tot_N2_cm2( 915) =    17.92
 given_sigma_tot_N2_cm2( 916) =    25.00
 given_sigma_tot_N2_cm2( 917) =    27.72
 given_sigma_tot_N2_cm2( 918) =    31.80
 given_sigma_tot_N2_cm2( 919) =    46.40
 given_sigma_tot_N2_cm2( 920) =    34.12
 given_sigma_tot_N2_cm2( 921) =    15.90
 given_sigma_tot_N2_cm2( 922) =    18.62
 given_sigma_tot_N2_cm2( 923) =    20.66
 given_sigma_tot_N2_cm2( 924) =    22.79
 given_sigma_tot_N2_cm2( 925) =    51.12
 given_sigma_tot_N2_cm2( 926) =    38.60
 given_sigma_tot_N2_cm2( 927) =    28.60
 given_sigma_tot_N2_cm2( 928) =    18.60
 given_sigma_tot_N2_cm2( 929) =    21.52
 given_sigma_tot_N2_cm2( 930) =    33.20
 given_sigma_tot_N2_cm2( 931) =    25.85
 given_sigma_tot_N2_cm2( 932) =    22.70
 given_sigma_tot_N2_cm2( 933) =    20.00
 given_sigma_tot_N2_cm2( 934) =    18.20
 given_sigma_tot_N2_cm2( 935) =    16.82
 given_sigma_tot_N2_cm2( 936) =    14.06
 given_sigma_tot_N2_cm2( 937) =    27.30
 given_sigma_tot_N2_cm2( 938) =    61.94
 given_sigma_tot_N2_cm2( 939) =    30.00
 given_sigma_tot_N2_cm2( 940) =    44.00
 given_sigma_tot_N2_cm2( 941) =    25.00
 given_sigma_tot_N2_cm2( 942) =    22.00
 given_sigma_tot_N2_cm2( 943) =    20.00
 given_sigma_tot_N2_cm2( 944) =    17.00
 given_sigma_tot_N2_cm2( 945) =    12.50
 given_sigma_tot_N2_cm2( 946) =    50.00
 given_sigma_tot_N2_cm2( 947) =    69.20
 given_sigma_tot_N2_cm2( 948) =    42.68
 given_sigma_tot_N2_cm2( 949) =    25.00
 given_sigma_tot_N2_cm2( 950) =    16.80
 given_sigma_tot_N2_cm2( 951) =    66.70
 given_sigma_tot_N2_cm2( 952) =    62.85
 given_sigma_tot_N2_cm2( 953) =    27.86
 given_sigma_tot_N2_cm2( 954) =    19.04
 given_sigma_tot_N2_cm2( 955) =    11.90
 given_sigma_tot_N2_cm2( 956) =    23.82
 given_sigma_tot_N2_cm2( 957) =    27.80
 given_sigma_tot_N2_cm2( 958) =    27.80
 given_sigma_tot_N2_cm2( 959) =    36.90
 given_sigma_tot_N2_cm2( 960) =    46.00
 given_sigma_tot_N2_cm2( 961) =    48.00
 given_sigma_tot_N2_cm2( 962) =    36.40
 given_sigma_tot_N2_cm2( 963) =    25.00
 given_sigma_tot_N2_cm2( 964) =    17.30
 given_sigma_tot_N2_cm2( 965) =    16.00
 given_sigma_tot_N2_cm2( 966) =    25.00
 given_sigma_tot_N2_cm2( 967) =    25.00
 given_sigma_tot_N2_cm2( 968) =    68.10
 given_sigma_tot_N2_cm2( 969) =    18.20
 given_sigma_tot_N2_cm2( 970) =    13.30
 given_sigma_tot_N2_cm2( 971) =    18.00
 given_sigma_tot_N2_cm2( 972) =   150.00
 given_sigma_tot_N2_cm2( 973) =    91.00
 given_sigma_tot_N2_cm2( 974) =    70.59
 given_sigma_tot_N2_cm2( 975) =    31.80
 given_sigma_tot_N2_cm2( 976) =    14.50
 given_sigma_tot_N2_cm2( 977) =    13.95
 given_sigma_tot_N2_cm2( 978) =    13.74
 given_sigma_tot_N2_cm2( 979) =    13.41
 given_sigma_tot_N2_cm2( 980) =    13.30
 given_sigma_tot_N2_cm2( 981) =    13.30
 given_sigma_tot_N2_cm2( 982) =    13.30
 given_sigma_tot_N2_cm2( 983) =    40.90
 given_sigma_tot_N2_cm2( 984) =    22.30
 given_sigma_tot_N2_cm2( 985) =    25.30
 given_sigma_tot_N2_cm2( 986) =    42.30
 given_sigma_tot_N2_cm2( 987) =    25.00
 given_sigma_tot_N2_cm2( 988) =    32.70
 given_sigma_tot_N2_cm2( 989) =    42.30
 given_sigma_tot_N2_cm2( 990) =    26.85
 given_sigma_tot_N2_cm2( 991) =    20.67
 given_sigma_tot_N2_cm2( 992) =    14.49
 given_sigma_tot_N2_cm2( 993) =    19.56
 given_sigma_tot_N2_cm2( 994) =    25.00
 given_sigma_tot_N2_cm2( 995) =   141.80
 given_sigma_tot_N2_cm2( 996) =    50.00
 given_sigma_tot_N2_cm2( 997) =    30.00
 given_sigma_tot_N2_cm2( 998) =    21.50
 given_sigma_tot_N2_cm2( 999) =    18.00
 given_sigma_tot_N2_cm2(1000) =    14.50
 given_sigma_tot_N2_cm2(1001) =    33.10
 given_sigma_tot_N2_cm2(1002) =    40.90
 given_sigma_tot_N2_cm2(1003) =    31.40
 given_sigma_tot_N2_cm2(1004) =   143.81
 given_sigma_tot_N2_cm2(1005) =   125.00
 given_sigma_tot_N2_cm2(1006) =    59.10
 given_sigma_tot_N2_cm2(1007) =    24.50
 given_sigma_tot_N2_cm2(1008) =    20.20
 given_sigma_tot_N2_cm2(1009) =    15.90
 given_sigma_tot_N2_cm2(1010) =    38.61
 given_sigma_tot_N2_cm2(1011) =    47.70
 given_sigma_tot_N2_cm2(1012) =    40.54
 given_sigma_tot_N2_cm2(1013) =    28.60
 given_sigma_tot_N2_cm2(1014) =    14.90
 given_sigma_tot_N2_cm2(1015) =    14.70
 given_sigma_tot_N2_cm2(1016) =    25.15
 given_sigma_tot_N2_cm2(1017) =    32.25
 given_sigma_tot_N2_cm2(1018) =    50.00
 given_sigma_tot_N2_cm2(1019) =   118.20
 given_sigma_tot_N2_cm2(1020) =    75.00
 given_sigma_tot_N2_cm2(1021) =    34.00
 given_sigma_tot_N2_cm2(1022) =    15.00
 given_sigma_tot_N2_cm2(1023) =    17.56
 given_sigma_tot_N2_cm2(1024) =    70.92
 given_sigma_tot_N2_cm2(1025) =   150.00
 given_sigma_tot_N2_cm2(1026) =   234.10
 given_sigma_tot_N2_cm2(1027) =   224.36
 given_sigma_tot_N2_cm2(1028) =   204.87
 given_sigma_tot_N2_cm2(1029) =   155.00
 given_sigma_tot_N2_cm2(1030) =   125.00
 given_sigma_tot_N2_cm2(1031) =    50.00
 given_sigma_tot_N2_cm2(1032) =    45.50
 given_sigma_tot_N2_cm2(1033) =    56.80
 given_sigma_tot_N2_cm2(1034) =    35.60
 given_sigma_tot_N2_cm2(1035) =    10.50
 given_sigma_tot_N2_cm2(1036) =     8.48
 given_sigma_tot_N2_cm2(1037) =     7.67
 given_sigma_tot_N2_cm2(1038) =     9.04
 given_sigma_tot_N2_cm2(1039) =    11.32
 given_sigma_tot_N2_cm2(1040) =    13.60
 given_sigma_tot_N2_cm2(1041) =    45.50
 given_sigma_tot_N2_cm2(1042) =    25.58
 given_sigma_tot_N2_cm2(1043) =    18.60
 given_sigma_tot_N2_cm2(1044) =   125.00
 given_sigma_tot_N2_cm2(1045) =    58.10
 given_sigma_tot_N2_cm2(1046) =    53.24
 given_sigma_tot_N2_cm2(1047) =    51.62
 given_sigma_tot_N2_cm2(1048) =   125.00
 given_sigma_tot_N2_cm2(1049) =    71.80
 given_sigma_tot_N2_cm2(1050) =    43.60
 given_sigma_tot_N2_cm2(1051) =    35.14
 given_sigma_tot_N2_cm2(1052) =    25.69
 given_sigma_tot_N2_cm2(1053) =    15.90
 given_sigma_tot_N2_cm2(1054) =    25.00
 given_sigma_tot_N2_cm2(1055) =    39.75
 given_sigma_tot_N2_cm2(1056) =    54.50
 given_sigma_tot_N2_cm2(1057) =    46.22
 given_sigma_tot_N2_cm2(1058) =    13.10
 given_sigma_tot_N2_cm2(1059) =     8.00
 given_sigma_tot_N2_cm2(1060) =     8.00
 given_sigma_tot_N2_cm2(1061) =     8.00
 given_sigma_tot_N2_cm2(1062) =     8.00
 given_sigma_tot_N2_cm2(1063) =     8.00
 given_sigma_tot_N2_cm2(1064) =     8.00
 given_sigma_tot_N2_cm2(1065) =    16.70
 given_sigma_tot_N2_cm2(1066) =    40.00
 given_sigma_tot_N2_cm2(1067) =   104.00
 given_sigma_tot_N2_cm2(1068) =    40.00
 given_sigma_tot_N2_cm2(1069) =    16.60
 given_sigma_tot_N2_cm2(1070) =     1.00
 given_sigma_tot_N2_cm2(1071) =     3.85
 given_sigma_tot_N2_cm2(1072) =     6.70
 given_sigma_tot_N2_cm2(1073) =     7.82
 given_sigma_tot_N2_cm2(1074) =     8.30
 given_sigma_tot_N2_cm2(1075) =     8.80
 given_sigma_tot_N2_cm2(1076) =     9.96
 given_sigma_tot_N2_cm2(1077) =    10.90
 given_sigma_tot_N2_cm2(1078) =    12.14
 given_sigma_tot_N2_cm2(1079) =    12.75
 given_sigma_tot_N2_cm2(1080) =    56.86
 given_sigma_tot_N2_cm2(1081) =    18.00
 given_sigma_tot_N2_cm2(1082) =    95.95
 given_sigma_tot_N2_cm2(1083) =   159.70
 given_sigma_tot_N2_cm2(1084) =    46.70
 given_sigma_tot_N2_cm2(1085) =     6.70
 given_sigma_tot_N2_cm2(1086) =     4.42
 given_sigma_tot_N2_cm2(1087) =     1.00
 given_sigma_tot_N2_cm2(1088) =    33.30
 given_sigma_tot_N2_cm2(1089) =    33.30
 given_sigma_tot_N2_cm2(1090) =    13.92
 given_sigma_tot_N2_cm2(1091) =     1.00
 given_sigma_tot_N2_cm2(1092) =     8.00
 given_sigma_tot_N2_cm2(1093) =     4.50
 given_sigma_tot_N2_cm2(1094) =     3.35
 given_sigma_tot_N2_cm2(1095) =     2.05
 given_sigma_tot_N2_cm2(1096) =     1.00
 given_sigma_tot_N2_cm2(1097) =     1.96
 given_sigma_tot_N2_cm2(1098) =     2.81
 given_sigma_tot_N2_cm2(1099) =     3.26
 given_sigma_tot_N2_cm2(1100) =     3.96
 given_sigma_tot_N2_cm2(1101) =     4.70
 given_sigma_tot_N2_cm2(1102) =    15.50
 given_sigma_tot_N2_cm2(1103) =    24.70
 given_sigma_tot_N2_cm2(1104) =   100.00
 given_sigma_tot_N2_cm2(1105) =    33.40
 given_sigma_tot_N2_cm2(1106) =    50.00
 given_sigma_tot_N2_cm2(1107) =    62.00
 given_sigma_tot_N2_cm2(1108) =    80.00
 given_sigma_tot_N2_cm2(1109) =    14.70
 given_sigma_tot_N2_cm2(1110) =    11.50
 given_sigma_tot_N2_cm2(1111) =     6.70
 given_sigma_tot_N2_cm2(1112) =    14.00
 given_sigma_tot_N2_cm2(1113) =    24.92
 given_sigma_tot_N2_cm2(1114) =    33.92
 given_sigma_tot_N2_cm2(1115) =    40.00
 given_sigma_tot_N2_cm2(1116) =    44.94
 given_sigma_tot_N2_cm2(1117) =    53.00
 given_sigma_tot_N2_cm2(1118) =    60.20
 given_sigma_tot_N2_cm2(1119) =    66.00
 given_sigma_tot_N2_cm2(1120) =    71.20
 given_sigma_tot_N2_cm2(1121) =    74.84
 given_sigma_tot_N2_cm2(1122) =    79.00
 given_sigma_tot_N2_cm2(1123) =    92.00
 given_sigma_tot_N2_cm2(1124) =    12.00
 given_sigma_tot_N2_cm2(1125) =     6.60
 given_sigma_tot_N2_cm2(1126) =     3.00
 given_sigma_tot_N2_cm2(1127) =     7.30
 given_sigma_tot_N2_cm2(1128) =    26.70
 given_sigma_tot_N2_cm2(1129) =    18.66
 given_sigma_tot_N2_cm2(1130) =    13.30
 given_sigma_tot_N2_cm2(1131) =     3.00
 given_sigma_tot_N2_cm2(1132) =     3.00
 given_sigma_tot_N2_cm2(1133) =     3.00
 given_sigma_tot_N2_cm2(1134) =     3.00
 given_sigma_tot_N2_cm2(1135) =     8.00
 given_sigma_tot_N2_cm2(1136) =    12.24
 given_sigma_tot_N2_cm2(1137) =     8.92
 given_sigma_tot_N2_cm2(1138) =     6.00
 given_sigma_tot_N2_cm2(1139) =    33.30
 given_sigma_tot_N2_cm2(1140) =    25.30
 given_sigma_tot_N2_cm2(1141) =    18.10
 given_sigma_tot_N2_cm2(1142) =    13.30
 given_sigma_tot_N2_cm2(1143) =     3.00
 given_sigma_tot_N2_cm2(1144) =     3.00
 given_sigma_tot_N2_cm2(1145) =     3.00
 given_sigma_tot_N2_cm2(1146) =     3.00
 given_sigma_tot_N2_cm2(1147) =     3.00
 given_sigma_tot_N2_cm2(1148) =     3.00
 given_sigma_tot_N2_cm2(1149) =     3.00
 given_sigma_tot_N2_cm2(1150) =     3.00
 given_sigma_tot_N2_cm2(1151) =     3.00
 given_sigma_tot_N2_cm2(1152) =     3.00
 given_sigma_tot_N2_cm2(1153) =     3.00
 given_sigma_tot_N2_cm2(1154) =     3.00
 given_sigma_tot_N2_cm2(1155) =     3.00
 given_sigma_tot_N2_cm2(1156) =     3.00
 given_sigma_tot_N2_cm2(1157) =     5.30
 given_sigma_tot_N2_cm2(1158) =     8.50
 given_sigma_tot_N2_cm2(1159) =    13.30
 given_sigma_tot_N2_cm2(1160) =    13.30
 given_sigma_tot_N2_cm2(1161) =    13.30
 given_sigma_tot_N2_cm2(1162) =    83.30
 given_sigma_tot_N2_cm2(1163) =   180.00
 given_sigma_tot_N2_cm2(1164) =   135.75
 given_sigma_tot_N2_cm2(1165) =     3.00
 given_sigma_tot_N2_cm2(1166) =     3.00
 given_sigma_tot_N2_cm2(1167) =     3.00
 given_sigma_tot_N2_cm2(1168) =    20.00
 given_sigma_tot_N2_cm2(1169) =    15.98
 given_sigma_tot_N2_cm2(1170) =    13.30
 given_sigma_tot_N2_cm2(1171) =    41.32
 given_sigma_tot_N2_cm2(1172) =    55.33
 given_sigma_tot_N2_cm2(1173) =     3.00
 given_sigma_tot_N2_cm2(1174) =     1.00
 given_sigma_tot_N2_cm2(1175) =    50.00
 given_sigma_tot_N2_cm2(1176) =    26.70
 given_sigma_tot_N2_cm2(1177) =    20.00
 given_sigma_tot_N2_cm2(1178) =    96.70
 given_sigma_tot_N2_cm2(1179) =    13.30
 given_sigma_tot_N2_cm2(1180) =     6.70
 given_sigma_tot_N2_cm2(1181) =     3.00
 given_sigma_tot_N2_cm2(1182) =     1.80
 given_sigma_tot_N2_cm2(1183) =     1.00
 given_sigma_tot_N2_cm2(1184) =    26.70
 given_sigma_tot_N2_cm2(1185) =    13.30
 given_sigma_tot_N2_cm2(1186) =     3.00
 given_sigma_tot_N2_cm2(1187) =     1.40
 given_sigma_tot_N2_cm2(1188) =     1.00
 given_sigma_tot_N2_cm2(1189) =     1.00
 given_sigma_tot_N2_cm2(1190) =     1.00
 given_sigma_tot_N2_cm2(1191) =     1.00
 given_sigma_tot_N2_cm2(1192) =     1.00
 given_sigma_tot_N2_cm2(1193) =     7.28
 given_sigma_tot_N2_cm2(1194) =    16.70
 given_sigma_tot_N2_cm2(1195) =    46.70
 given_sigma_tot_N2_cm2(1196) =    26.66
 given_sigma_tot_N2_cm2(1197) =    13.30
 given_sigma_tot_N2_cm2(1198) =     6.70
 given_sigma_tot_N2_cm2(1199) =     4.48
 given_sigma_tot_N2_cm2(1200) =     3.00
 given_sigma_tot_N2_cm2(1201) =     2.24
 given_sigma_tot_N2_cm2(1202) =     1.00
 given_sigma_tot_N2_cm2(1203) =     3.00
 given_sigma_tot_N2_cm2(1204) =    36.48
 given_sigma_tot_N2_cm2(1205) =    86.70
 given_sigma_tot_N2_cm2(1206) =    66.50
 given_sigma_tot_N2_cm2(1207) =    46.70
 given_sigma_tot_N2_cm2(1208) =    20.00
 given_sigma_tot_N2_cm2(1209) =    13.30
 given_sigma_tot_N2_cm2(1210) =     3.00
 given_sigma_tot_N2_cm2(1211) =   173.30
 given_sigma_tot_N2_cm2(1212) =    20.00
 given_sigma_tot_N2_cm2(1213) =     8.70
 given_sigma_tot_N2_cm2(1214) =    12.13
 given_sigma_tot_N2_cm2(1215) =    16.70
 given_sigma_tot_N2_cm2(1216) =    13.30
 given_sigma_tot_N2_cm2(1217) =     6.70
 given_sigma_tot_N2_cm2(1218) =     3.00
 given_sigma_tot_N2_cm2(1219) =     3.00
 given_sigma_tot_N2_cm2(1220) =     3.00
 given_sigma_tot_N2_cm2(1221) =     1.00
 given_sigma_tot_N2_cm2(1222) =     6.00
 given_sigma_tot_N2_cm2(1223) =     1.00
 given_sigma_tot_N2_cm2(1224) =     1.00
 given_sigma_tot_N2_cm2(1225) =     1.00
 given_sigma_tot_N2_cm2(1226) =    33.40
 given_sigma_tot_N2_cm2(1227) =     6.00
 given_sigma_tot_N2_cm2(1228) =     3.00
 given_sigma_tot_N2_cm2(1229) =     1.00
 given_sigma_tot_N2_cm2(1230) =    86.70
 given_sigma_tot_N2_cm2(1231) =   106.70
 given_sigma_tot_N2_cm2(1232) =     3.00
 given_sigma_tot_N2_cm2(1233) =     1.00
 given_sigma_tot_N2_cm2(1234) =    65.00
 given_sigma_tot_N2_cm2(1235) =    32.00
 given_sigma_tot_N2_cm2(1236) =    14.00
 given_sigma_tot_N2_cm2(1237) =    11.00
 given_sigma_tot_N2_cm2(1238) =    10.00
 given_sigma_tot_N2_cm2(1239) =     6.00
 given_sigma_tot_N2_cm2(1240) =     3.00
 given_sigma_tot_N2_cm2(1241) =   120.00
 given_sigma_tot_N2_cm2(1242) =    25.00
 given_sigma_tot_N2_cm2(1243) =    26.00
 given_sigma_tot_N2_cm2(1244) =    80.00
 given_sigma_tot_N2_cm2(1245) =    33.00
 given_sigma_tot_N2_cm2(1246) =    12.00
 given_sigma_tot_N2_cm2(1247) =    10.00
 given_sigma_tot_N2_cm2(1248) =     6.00
 given_sigma_tot_N2_cm2(1249) =     3.00
 given_sigma_tot_N2_cm2(1250) =     3.00
 given_sigma_tot_N2_cm2(1251) =     3.00
 given_sigma_tot_N2_cm2(1252) =     3.00
 given_sigma_tot_N2_cm2(1253) =   107.00
 given_sigma_tot_N2_cm2(1254) =    60.00
 given_sigma_tot_N2_cm2(1255) =    54.00
 given_sigma_tot_N2_cm2(1256) =    30.00
 given_sigma_tot_N2_cm2(1257) =    24.00
 given_sigma_tot_N2_cm2(1258) =    23.12
 given_sigma_tot_N2_cm2(1259) =    20.00
 given_sigma_tot_N2_cm2(1260) =    13.28
 given_sigma_tot_N2_cm2(1261) =     9.60
 given_sigma_tot_N2_cm2(1262) =     6.00
 given_sigma_tot_N2_cm2(1263) =     4.20
 given_sigma_tot_N2_cm2(1264) =     3.00
 given_sigma_tot_N2_cm2(1265) =     5.00
 given_sigma_tot_N2_cm2(1266) =   110.00
 given_sigma_tot_N2_cm2(1267) =    67.00
 given_sigma_tot_N2_cm2(1268) =    51.00
 given_sigma_tot_N2_cm2(1269) =     1.00
 given_sigma_tot_N2_cm2(1270) =    54.00
 given_sigma_tot_N2_cm2(1271) =    27.00
 given_sigma_tot_N2_cm2(1272) =    16.80
 given_sigma_tot_N2_cm2(1273) =    10.00
 given_sigma_tot_N2_cm2(1274) =     1.00
 given_sigma_tot_N2_cm2(1275) =     3.00
 given_sigma_tot_N2_cm2(1276) =     1.00
 given_sigma_tot_N2_cm2(1277) =     3.00
 given_sigma_tot_N2_cm2(1278) =     6.60
 given_sigma_tot_N2_cm2(1279) =     9.00
 given_sigma_tot_N2_cm2(1280) =    15.00
 given_sigma_tot_N2_cm2(1281) =    80.00
 given_sigma_tot_N2_cm2(1282) =    18.00
 given_sigma_tot_N2_cm2(1283) =    30.00
 given_sigma_tot_N2_cm2(1284) =    12.00
 given_sigma_tot_N2_cm2(1285) =    20.00
 given_sigma_tot_N2_cm2(1286) =    14.00
 given_sigma_tot_N2_cm2(1287) =    12.00
 given_sigma_tot_N2_cm2(1288) =     3.00
 given_sigma_tot_N2_cm2(1289) =    65.00
 given_sigma_tot_N2_cm2(1290) =    25.00
 given_sigma_tot_N2_cm2(1291) =    15.00
 given_sigma_tot_N2_cm2(1292) =    15.00
 given_sigma_tot_N2_cm2(1293) =     1.00
 given_sigma_tot_N2_cm2(1294) =    14.00
 given_sigma_tot_N2_cm2(1295) =    10.00
 given_sigma_tot_N2_cm2(1296) =     3.00
 given_sigma_tot_N2_cm2(1297) =     3.00
 given_sigma_tot_N2_cm2(1298) =     3.00
 given_sigma_tot_N2_cm2(1299) =     3.00
 given_sigma_tot_N2_cm2(1300) =     3.00
 given_sigma_tot_N2_cm2(1301) =     3.00
 given_sigma_tot_N2_cm2(1302) =    25.00
 given_sigma_tot_N2_cm2(1303) =    11.00
 given_sigma_tot_N2_cm2(1304) =   104.00
 given_sigma_tot_N2_cm2(1305) =    64.00
 given_sigma_tot_N2_cm2(1306) =     1.00
 given_sigma_tot_N2_cm2(1307) =     3.00
 given_sigma_tot_N2_cm2(1308) =     4.20
 given_sigma_tot_N2_cm2(1309) =     6.00
 given_sigma_tot_N2_cm2(1310) =    10.00
 given_sigma_tot_N2_cm2(1311) =    12.00
 given_sigma_tot_N2_cm2(1312) =    11.20
 given_sigma_tot_N2_cm2(1313) =    10.00
 given_sigma_tot_N2_cm2(1314) =     3.00
 given_sigma_tot_N2_cm2(1315) =    31.00
 given_sigma_tot_N2_cm2(1316) =    27.00
 given_sigma_tot_N2_cm2(1317) =     6.00
 given_sigma_tot_N2_cm2(1318) =    10.00
 given_sigma_tot_N2_cm2(1319) =     1.00
 given_sigma_tot_N2_cm2(1320) =     5.20
 given_sigma_tot_N2_cm2(1321) =     8.00
 given_sigma_tot_N2_cm2(1322) =    14.00
 given_sigma_tot_N2_cm2(1323) =    12.00
 given_sigma_tot_N2_cm2(1324) =    10.00
 given_sigma_tot_N2_cm2(1325) =     3.00
 given_sigma_tot_N2_cm2(1326) =     1.40
 given_sigma_tot_N2_cm2(1327) =     1.00
 given_sigma_tot_N2_cm2(1328) =     1.00
 given_sigma_tot_N2_cm2(1329) =    15.00
 given_sigma_tot_N2_cm2(1330) =    15.00
 given_sigma_tot_N2_cm2(1331) =     9.00
 given_sigma_tot_N2_cm2(1332) =     6.00
 given_sigma_tot_N2_cm2(1333) =     4.20
 given_sigma_tot_N2_cm2(1334) =     3.00
 given_sigma_tot_N2_cm2(1335) =     1.00
 given_sigma_tot_N2_cm2(1336) =    14.00
 given_sigma_tot_N2_cm2(1337) =     6.00
 given_sigma_tot_N2_cm2(1338) =     3.40
 given_sigma_tot_N2_cm2(1339) =     1.00
 given_sigma_tot_N2_cm2(1340) =     1.00
 given_sigma_tot_N2_cm2(1341) =     1.00
 given_sigma_tot_N2_cm2(1342) =    82.00
 given_sigma_tot_N2_cm2(1343) =    53.80
 given_sigma_tot_N2_cm2(1344) =    35.00
 given_sigma_tot_N2_cm2(1345) =    15.00
 given_sigma_tot_N2_cm2(1346) =     3.00
 given_sigma_tot_N2_cm2(1347) =     3.00
 given_sigma_tot_N2_cm2(1348) =     3.00
 given_sigma_tot_N2_cm2(1349) =     1.00
 given_sigma_tot_N2_cm2(1350) =     1.00
 given_sigma_tot_N2_cm2(1351) =    42.00
 given_sigma_tot_N2_cm2(1352) =     1.00
 given_sigma_tot_N2_cm2(1353) =     1.80
 given_sigma_tot_N2_cm2(1354) =     3.00
 given_sigma_tot_N2_cm2(1355) =     6.00
 given_sigma_tot_N2_cm2(1356) =     3.00
 given_sigma_tot_N2_cm2(1357) =     3.00
 given_sigma_tot_N2_cm2(1358) =     1.00
 given_sigma_tot_N2_cm2(1359) =    97.00
 given_sigma_tot_N2_cm2(1360) =    50.75
 given_sigma_tot_N2_cm2(1361) =    27.64
 given_sigma_tot_N2_cm2(1362) =    17.00
 given_sigma_tot_N2_cm2(1363) =    11.00
 given_sigma_tot_N2_cm2(1364) =     3.00
 given_sigma_tot_N2_cm2(1365) =    58.00
 given_sigma_tot_N2_cm2(1366) =    17.00
 given_sigma_tot_N2_cm2(1367) =    18.50
 given_sigma_tot_N2_cm2(1368) =    19.10
 given_sigma_tot_N2_cm2(1369) =    20.00
 given_sigma_tot_N2_cm2(1370) =    15.20
 given_sigma_tot_N2_cm2(1371) =    12.00
 given_sigma_tot_N2_cm2(1372) =    22.00
 given_sigma_tot_N2_cm2(1373) =     3.00
 given_sigma_tot_N2_cm2(1374) =    14.00
 given_sigma_tot_N2_cm2(1375) =     9.60
 given_sigma_tot_N2_cm2(1376) =     5.20
 given_sigma_tot_N2_cm2(1377) =    18.00
 given_sigma_tot_N2_cm2(1378) =    13.20
 given_sigma_tot_N2_cm2(1379) =    10.00
 given_sigma_tot_N2_cm2(1380) =     2.00
 given_sigma_tot_N2_cm2(1381) =     2.00
 given_sigma_tot_N2_cm2(1382) =    29.00
 given_sigma_tot_N2_cm2(1383) =     3.00
 given_sigma_tot_N2_cm2(1384) =     4.50
 given_sigma_tot_N2_cm2(1385) =     6.00
 given_sigma_tot_N2_cm2(1386) =     6.00
 given_sigma_tot_N2_cm2(1387) =    15.00
 given_sigma_tot_N2_cm2(1388) =     8.20
 given_sigma_tot_N2_cm2(1389) =     7.00
 given_sigma_tot_N2_cm2(1390) =     4.72
 given_sigma_tot_N2_cm2(1391) =     4.00
 given_sigma_tot_N2_cm2(1392) =    10.00
 given_sigma_tot_N2_cm2(1393) =    11.30
 given_sigma_tot_N2_cm2(1394) =     5.00
 given_sigma_tot_N2_cm2(1395) =    12.00
 given_sigma_tot_N2_cm2(1396) =     7.00
 given_sigma_tot_N2_cm2(1397) =     6.00
 given_sigma_tot_N2_cm2(1398) =     6.00
 given_sigma_tot_N2_cm2(1399) =     1.00
 given_sigma_tot_N2_cm2(1400) =     7.00
 given_sigma_tot_N2_cm2(1401) =    20.50
 given_sigma_tot_N2_cm2(1402) =    25.00
 given_sigma_tot_N2_cm2(1403) =     1.00
 given_sigma_tot_N2_cm2(1404) =   148.00
 given_sigma_tot_N2_cm2(1405) =   106.00
 given_sigma_tot_N2_cm2(1406) =    78.00
 given_sigma_tot_N2_cm2(1407) =    27.00
 given_sigma_tot_N2_cm2(1408) =    29.00
 given_sigma_tot_N2_cm2(1409) =    18.00
 given_sigma_tot_N2_cm2(1410) =    15.00
 given_sigma_tot_N2_cm2(1411) =    22.00
 given_sigma_tot_N2_cm2(1412) =     1.00
 given_sigma_tot_N2_cm2(1413) =    27.00
 given_sigma_tot_N2_cm2(1414) =    22.80
 given_sigma_tot_N2_cm2(1415) =    20.00
 given_sigma_tot_N2_cm2(1416) =    11.00
 given_sigma_tot_N2_cm2(1417) =    11.00
 given_sigma_tot_N2_cm2(1418) =     1.00
 given_sigma_tot_N2_cm2(1419) =     1.00
 given_sigma_tot_N2_cm2(1420) =     3.00
 given_sigma_tot_N2_cm2(1421) =     3.00
 given_sigma_tot_N2_cm2(1422) =     3.00
 given_sigma_tot_N2_cm2(1423) =     1.00
 given_sigma_tot_N2_cm2(1424) =     1.00
 given_sigma_tot_N2_cm2(1425) =    25.00
 given_sigma_tot_N2_cm2(1426) =    13.48
 given_sigma_tot_N2_cm2(1427) =     1.00
 given_sigma_tot_N2_cm2(1428) =    17.20
 given_sigma_tot_N2_cm2(1429) =    28.00
 given_sigma_tot_N2_cm2(1430) =     2.00
 given_sigma_tot_N2_cm2(1431) =     2.00
 given_sigma_tot_N2_cm2(1432) =    12.00
 given_sigma_tot_N2_cm2(1433) =     3.00
 given_sigma_tot_N2_cm2(1434) =     6.00
 given_sigma_tot_N2_cm2(1435) =    18.60
 given_sigma_tot_N2_cm2(1436) =    12.00
 given_sigma_tot_N2_cm2(1437) =     8.40
 given_sigma_tot_N2_cm2(1438) =     3.00
 given_sigma_tot_N2_cm2(1439) =     6.00
 given_sigma_tot_N2_cm2(1440) =    18.00
 given_sigma_tot_N2_cm2(1441) =   134.00
 given_sigma_tot_N2_cm2(1442) =    53.00
 given_sigma_tot_N2_cm2(1443) =    62.00
 given_sigma_tot_N2_cm2(1444) =    47.00
 given_sigma_tot_N2_cm2(1445) =    21.00
 given_sigma_tot_N2_cm2(1446) =     1.00
 given_sigma_tot_N2_cm2(1447) =    98.00
 given_sigma_tot_N2_cm2(1448) =    51.00
 given_sigma_tot_N2_cm2(1449) =    30.00
 given_sigma_tot_N2_cm2(1450) =    20.40
 given_sigma_tot_N2_cm2(1451) =    10.00
 given_sigma_tot_N2_cm2(1452) =     5.80
 given_sigma_tot_N2_cm2(1453) =     3.00
 given_sigma_tot_N2_cm2(1454) =     3.00
 given_sigma_tot_N2_cm2(1455) =     3.00
 given_sigma_tot_N2_cm2(1456) =     1.00
 given_sigma_tot_N2_cm2(1457) =    67.00
 given_sigma_tot_N2_cm2(1458) =    60.00
 given_sigma_tot_N2_cm2(1459) =    31.00
 given_sigma_tot_N2_cm2(1460) =    20.00
 given_sigma_tot_N2_cm2(1461) =    27.00
 given_sigma_tot_N2_cm2(1462) =    10.00
 given_sigma_tot_N2_cm2(1463) =     3.00
 given_sigma_tot_N2_cm2(1464) =     6.00
 given_sigma_tot_N2_cm2(1465) =     2.00
 given_sigma_tot_N2_cm2(1466) =     1.00
 given_sigma_tot_N2_cm2(1467) =     6.00
 given_sigma_tot_N2_cm2(1468) =    27.00
 given_sigma_tot_N2_cm2(1469) =   120.00
 given_sigma_tot_N2_cm2(1470) =    94.40
 given_sigma_tot_N2_cm2(1471) =    31.00
 given_sigma_tot_N2_cm2(1472) =    20.00
 given_sigma_tot_N2_cm2(1473) =     3.00
 given_sigma_tot_N2_cm2(1474) =     5.00
 given_sigma_tot_N2_cm2(1475) =     3.80
 given_sigma_tot_N2_cm2(1476) =     3.00
 given_sigma_tot_N2_cm2(1477) =     3.00
 given_sigma_tot_N2_cm2(1478) =     1.00
 given_sigma_tot_N2_cm2(1479) =     2.00
 given_sigma_tot_N2_cm2(1480) =     8.00
 given_sigma_tot_N2_cm2(1481) =     5.00
 given_sigma_tot_N2_cm2(1482) =     2.00
 given_sigma_tot_N2_cm2(1483) =    88.00
 given_sigma_tot_N2_cm2(1484) =    33.00
 given_sigma_tot_N2_cm2(1485) =    27.00
 given_sigma_tot_N2_cm2(1486) =    22.00
 given_sigma_tot_N2_cm2(1487) =    26.00
 given_sigma_tot_N2_cm2(1488) =    11.00
 given_sigma_tot_N2_cm2(1489) =    14.00
 given_sigma_tot_N2_cm2(1490) =    25.00
 given_sigma_tot_N2_cm2(1491) =    23.00
 given_sigma_tot_N2_cm2(1492) =    18.20
 given_sigma_tot_N2_cm2(1493) =    15.00
 given_sigma_tot_N2_cm2(1494) =    17.00
 given_sigma_tot_N2_cm2(1495) =    15.00
 given_sigma_tot_N2_cm2(1496) =    15.00
 given_sigma_tot_N2_cm2(1497) =    13.00
 given_sigma_tot_N2_cm2(1498) =    10.00
 given_sigma_tot_N2_cm2(1499) =    22.00
 given_sigma_tot_N2_cm2(1500) =     0.00
 given_sigma_tot_N2_cm2(1501) =     0.00
 given_sigma_tot_N2_cm2(1502) =     0.00
 given_sigma_tot_N2_cm2(1503) =     0.00
 given_sigma_tot_N2_cm2(1504) =     0.00
 given_sigma_tot_N2_cm2(1505) =     0.00
 given_sigma_tot_N2_cm2(1506) =     0.00
 given_sigma_tot_N2_cm2(1507) =     0.00
 given_sigma_tot_N2_cm2(1508) =     0.00
 given_sigma_tot_N2_cm2(1509) =     0.00
 given_sigma_tot_N2_cm2(1510) =     0.00
 given_sigma_tot_N2_cm2(1511) =     0.00
 given_sigma_tot_N2_cm2(1512) =     0.00
 given_sigma_tot_N2_cm2(1513) =     0.00
 given_sigma_tot_N2_cm2(1514) =     0.00
 given_sigma_tot_N2_cm2(1515) =     0.00
 given_sigma_tot_N2_cm2(1516) =     0.00
 given_sigma_tot_N2_cm2(1517) =     0.00
 given_sigma_tot_N2_cm2(1518) =     0.00
 given_sigma_tot_N2_cm2(1519) =     0.00
 given_sigma_tot_N2_cm2(1520) =     0.00
 given_sigma_tot_N2_cm2(1521) =     0.00
 given_sigma_tot_N2_cm2(1522) =     0.00
 given_sigma_tot_N2_cm2(1523) =     0.00
 given_sigma_tot_N2_cm2(1524) =     0.00
 given_sigma_tot_N2_cm2(1525) =     0.00
 given_sigma_tot_N2_cm2(1526) =     0.00
 given_sigma_tot_N2_cm2(1527) =     0.00
 given_sigma_tot_N2_cm2(1528) =     0.00
 given_sigma_tot_N2_cm2(1529) =     0.00
 given_sigma_tot_N2_cm2(1530) =     0.00
 given_sigma_tot_N2_cm2(1531) =     0.00
 given_sigma_tot_N2_cm2(1532) =     0.00
 given_sigma_tot_N2_cm2(1533) =     0.00
 given_sigma_tot_N2_cm2(1534) =     0.00
 given_sigma_tot_N2_cm2(1535) =     0.00
 given_sigma_tot_N2_cm2(1536) =     0.00
 given_sigma_tot_N2_cm2(1537) =     0.00
 given_sigma_tot_N2_cm2(1538) =     0.00
 given_sigma_tot_N2_cm2(1539) =     0.00
 given_sigma_tot_N2_cm2(1540) =     0.00
 given_sigma_tot_N2_cm2(1541) =     0.00
 given_sigma_tot_N2_cm2(1542) =     0.00
 given_sigma_tot_N2_cm2(1543) =     0.00
 given_sigma_tot_N2_cm2(1544) =     0.00
 given_sigma_tot_N2_cm2(1545) =     0.00
 given_sigma_tot_N2_cm2(1546) =     0.00
 given_sigma_tot_N2_cm2(1547) =     0.00
 given_sigma_tot_N2_cm2(1548) =     0.00
 given_sigma_tot_N2_cm2(1549) =     0.00
 given_sigma_tot_N2_cm2(1550) =     0.00
 given_sigma_tot_N2_cm2(1551) =     0.00
 given_sigma_tot_N2_cm2(1552) =     0.00
 given_sigma_tot_N2_cm2(1553) =     0.00
 given_sigma_tot_N2_cm2(1554) =     0.00
 given_sigma_tot_N2_cm2(1555) =     0.00
 given_sigma_tot_N2_cm2(1556) =     0.00
 given_sigma_tot_N2_cm2(1557) =     0.00
 given_sigma_tot_N2_cm2(1558) =     0.00
 given_sigma_tot_N2_cm2(1559) =     0.00
 given_sigma_tot_N2_cm2(1560) =     0.00
 given_sigma_tot_N2_cm2(1561) =     0.00
 given_sigma_tot_N2_cm2(1562) =     0.00
 given_sigma_tot_N2_cm2(1563) =     0.00
 given_sigma_tot_N2_cm2(1564) =     0.00
 given_sigma_tot_N2_cm2(1565) =     0.00
 given_sigma_tot_N2_cm2(1566) =     0.00
 given_sigma_tot_N2_cm2(1567) =     0.00
 given_sigma_tot_N2_cm2(1568) =     0.00
 given_sigma_tot_N2_cm2(1569) =     0.00
 given_sigma_tot_N2_cm2(1570) =     0.00
 given_sigma_tot_N2_cm2(1571) =     0.00












  given_sigma_ion_O2_to_O2_cm2(   1) =     0.07 ; given_sigma_ion_O2_to_O_cm2(   1) =    0.00
  given_sigma_ion_O2_to_O2_cm2(   2) =     0.13 ; given_sigma_ion_O2_to_O_cm2(   2) =    0.00
  given_sigma_ion_O2_to_O2_cm2(   3) =     0.13 ; given_sigma_ion_O2_to_O_cm2(   3) =    0.00
  given_sigma_ion_O2_to_O2_cm2(   4) =     0.14 ; given_sigma_ion_O2_to_O_cm2(   4) =    0.00
  given_sigma_ion_O2_to_O2_cm2(   5) =     0.14 ; given_sigma_ion_O2_to_O_cm2(   5) =    0.00
  given_sigma_ion_O2_to_O2_cm2(   6) =     0.15 ; given_sigma_ion_O2_to_O_cm2(   6) =    0.00
  given_sigma_ion_O2_to_O2_cm2(   7) =     0.15 ; given_sigma_ion_O2_to_O_cm2(   7) =    0.00
  given_sigma_ion_O2_to_O2_cm2(   8) =     0.16 ; given_sigma_ion_O2_to_O_cm2(   8) =    0.00
  given_sigma_ion_O2_to_O2_cm2(   9) =     0.16 ; given_sigma_ion_O2_to_O_cm2(   9) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  10) =     0.19 ; given_sigma_ion_O2_to_O_cm2(  10) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  11) =     0.21 ; given_sigma_ion_O2_to_O_cm2(  11) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  12) =     0.28 ; given_sigma_ion_O2_to_O_cm2(  12) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  13) =     0.28 ; given_sigma_ion_O2_to_O_cm2(  13) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  14) =     0.31 ; given_sigma_ion_O2_to_O_cm2(  14) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  15) =     0.31 ; given_sigma_ion_O2_to_O_cm2(  15) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  16) =     0.31 ; given_sigma_ion_O2_to_O_cm2(  16) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  17) =     0.32 ; given_sigma_ion_O2_to_O_cm2(  17) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  18) =     0.33 ; given_sigma_ion_O2_to_O_cm2(  18) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  19) =     0.34 ; given_sigma_ion_O2_to_O_cm2(  19) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  20) =     0.35 ; given_sigma_ion_O2_to_O_cm2(  20) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  21) =     0.36 ; given_sigma_ion_O2_to_O_cm2(  21) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  22) =     0.36 ; given_sigma_ion_O2_to_O_cm2(  22) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  23) =     0.38 ; given_sigma_ion_O2_to_O_cm2(  23) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  24) =     0.39 ; given_sigma_ion_O2_to_O_cm2(  24) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  25) =     0.42 ; given_sigma_ion_O2_to_O_cm2(  25) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  26) =     0.43 ; given_sigma_ion_O2_to_O_cm2(  26) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  27) =     0.44 ; given_sigma_ion_O2_to_O_cm2(  27) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  28) =     0.44 ; given_sigma_ion_O2_to_O_cm2(  28) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  29) =     0.45 ; given_sigma_ion_O2_to_O_cm2(  29) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  30) =     0.47 ; given_sigma_ion_O2_to_O_cm2(  30) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  31) =     0.48 ; given_sigma_ion_O2_to_O_cm2(  31) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  32) =     0.50 ; given_sigma_ion_O2_to_O_cm2(  32) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  33) =     0.52 ; given_sigma_ion_O2_to_O_cm2(  33) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  34) =     0.53 ; given_sigma_ion_O2_to_O_cm2(  34) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  35) =     0.54 ; given_sigma_ion_O2_to_O_cm2(  35) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  36) =     0.54 ; given_sigma_ion_O2_to_O_cm2(  36) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  37) =     0.55 ; given_sigma_ion_O2_to_O_cm2(  37) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  38) =     0.56 ; given_sigma_ion_O2_to_O_cm2(  38) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  39) =     0.58 ; given_sigma_ion_O2_to_O_cm2(  39) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  40) =     0.59 ; given_sigma_ion_O2_to_O_cm2(  40) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  41) =     0.60 ; given_sigma_ion_O2_to_O_cm2(  41) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  42) =     0.61 ; given_sigma_ion_O2_to_O_cm2(  42) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  43) =     0.62 ; given_sigma_ion_O2_to_O_cm2(  43) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  44) =     0.62 ; given_sigma_ion_O2_to_O_cm2(  44) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  45) =     0.65 ; given_sigma_ion_O2_to_O_cm2(  45) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  46) =     0.67 ; given_sigma_ion_O2_to_O_cm2(  46) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  47) =     0.68 ; given_sigma_ion_O2_to_O_cm2(  47) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  48) =     0.70 ; given_sigma_ion_O2_to_O_cm2(  48) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  49) =     0.70 ; given_sigma_ion_O2_to_O_cm2(  49) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  50) =     0.72 ; given_sigma_ion_O2_to_O_cm2(  50) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  51) =     0.73 ; given_sigma_ion_O2_to_O_cm2(  51) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  52) =     0.74 ; given_sigma_ion_O2_to_O_cm2(  52) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  53) =     0.75 ; given_sigma_ion_O2_to_O_cm2(  53) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  54) =     0.75 ; given_sigma_ion_O2_to_O_cm2(  54) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  55) =     0.76 ; given_sigma_ion_O2_to_O_cm2(  55) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  56) =     0.76 ; given_sigma_ion_O2_to_O_cm2(  56) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  57) =     0.77 ; given_sigma_ion_O2_to_O_cm2(  57) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  58) =     0.78 ; given_sigma_ion_O2_to_O_cm2(  58) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  59) =     0.79 ; given_sigma_ion_O2_to_O_cm2(  59) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  60) =     0.80 ; given_sigma_ion_O2_to_O_cm2(  60) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  61) =     0.82 ; given_sigma_ion_O2_to_O_cm2(  61) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  62) =     0.84 ; given_sigma_ion_O2_to_O_cm2(  62) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  63) =     0.84 ; given_sigma_ion_O2_to_O_cm2(  63) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  64) =     0.86 ; given_sigma_ion_O2_to_O_cm2(  64) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  65) =     0.89 ; given_sigma_ion_O2_to_O_cm2(  65) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  66) =     0.89 ; given_sigma_ion_O2_to_O_cm2(  66) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  67) =     0.90 ; given_sigma_ion_O2_to_O_cm2(  67) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  68) =     0.93 ; given_sigma_ion_O2_to_O_cm2(  68) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  69) =     0.94 ; given_sigma_ion_O2_to_O_cm2(  69) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  70) =     0.96 ; given_sigma_ion_O2_to_O_cm2(  70) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  71) =     0.98 ; given_sigma_ion_O2_to_O_cm2(  71) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  72) =     0.99 ; given_sigma_ion_O2_to_O_cm2(  72) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  73) =     1.00 ; given_sigma_ion_O2_to_O_cm2(  73) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  74) =     1.00 ; given_sigma_ion_O2_to_O_cm2(  74) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  75) =     1.03 ; given_sigma_ion_O2_to_O_cm2(  75) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  76) =     1.04 ; given_sigma_ion_O2_to_O_cm2(  76) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  77) =     1.04 ; given_sigma_ion_O2_to_O_cm2(  77) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  78) =     1.06 ; given_sigma_ion_O2_to_O_cm2(  78) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  79) =     1.06 ; given_sigma_ion_O2_to_O_cm2(  79) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  80) =     1.09 ; given_sigma_ion_O2_to_O_cm2(  80) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  81) =     1.11 ; given_sigma_ion_O2_to_O_cm2(  81) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  82) =     1.12 ; given_sigma_ion_O2_to_O_cm2(  82) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  83) =     1.14 ; given_sigma_ion_O2_to_O_cm2(  83) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  84) =     1.14 ; given_sigma_ion_O2_to_O_cm2(  84) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  85) =     1.15 ; given_sigma_ion_O2_to_O_cm2(  85) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  86) =     1.16 ; given_sigma_ion_O2_to_O_cm2(  86) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  87) =     1.17 ; given_sigma_ion_O2_to_O_cm2(  87) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  88) =     1.18 ; given_sigma_ion_O2_to_O_cm2(  88) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  89) =     1.19 ; given_sigma_ion_O2_to_O_cm2(  89) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  90) =     1.21 ; given_sigma_ion_O2_to_O_cm2(  90) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  91) =     1.21 ; given_sigma_ion_O2_to_O_cm2(  91) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  92) =     1.22 ; given_sigma_ion_O2_to_O_cm2(  92) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  93) =     1.23 ; given_sigma_ion_O2_to_O_cm2(  93) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  94) =     1.24 ; given_sigma_ion_O2_to_O_cm2(  94) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  95) =     1.27 ; given_sigma_ion_O2_to_O_cm2(  95) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  96) =     1.28 ; given_sigma_ion_O2_to_O_cm2(  96) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  97) =     1.29 ; given_sigma_ion_O2_to_O_cm2(  97) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  98) =     1.30 ; given_sigma_ion_O2_to_O_cm2(  98) =    0.00
  given_sigma_ion_O2_to_O2_cm2(  99) =     1.31 ; given_sigma_ion_O2_to_O_cm2(  99) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 100) =     1.32 ; given_sigma_ion_O2_to_O_cm2( 100) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 101) =     1.33 ; given_sigma_ion_O2_to_O_cm2( 101) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 102) =     1.34 ; given_sigma_ion_O2_to_O_cm2( 102) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 103) =     1.36 ; given_sigma_ion_O2_to_O_cm2( 103) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 104) =     1.36 ; given_sigma_ion_O2_to_O_cm2( 104) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 105) =     1.38 ; given_sigma_ion_O2_to_O_cm2( 105) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 106) =     1.39 ; given_sigma_ion_O2_to_O_cm2( 106) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 107) =     1.40 ; given_sigma_ion_O2_to_O_cm2( 107) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 108) =     1.41 ; given_sigma_ion_O2_to_O_cm2( 108) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 109) =     1.43 ; given_sigma_ion_O2_to_O_cm2( 109) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 110) =     1.45 ; given_sigma_ion_O2_to_O_cm2( 110) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 111) =     1.46 ; given_sigma_ion_O2_to_O_cm2( 111) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 112) =     1.47 ; given_sigma_ion_O2_to_O_cm2( 112) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 113) =     1.49 ; given_sigma_ion_O2_to_O_cm2( 113) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 114) =     1.50 ; given_sigma_ion_O2_to_O_cm2( 114) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 115) =     1.51 ; given_sigma_ion_O2_to_O_cm2( 115) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 116) =     1.53 ; given_sigma_ion_O2_to_O_cm2( 116) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 117) =     1.53 ; given_sigma_ion_O2_to_O_cm2( 117) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 118) =     1.55 ; given_sigma_ion_O2_to_O_cm2( 118) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 119) =     1.56 ; given_sigma_ion_O2_to_O_cm2( 119) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 120) =     1.57 ; given_sigma_ion_O2_to_O_cm2( 120) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 121) =     1.58 ; given_sigma_ion_O2_to_O_cm2( 121) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 122) =     1.60 ; given_sigma_ion_O2_to_O_cm2( 122) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 123) =     1.61 ; given_sigma_ion_O2_to_O_cm2( 123) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 124) =     1.62 ; given_sigma_ion_O2_to_O_cm2( 124) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 125) =     1.64 ; given_sigma_ion_O2_to_O_cm2( 125) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 126) =     1.65 ; given_sigma_ion_O2_to_O_cm2( 126) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 127) =     1.67 ; given_sigma_ion_O2_to_O_cm2( 127) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 128) =     1.69 ; given_sigma_ion_O2_to_O_cm2( 128) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 129) =     1.71 ; given_sigma_ion_O2_to_O_cm2( 129) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 130) =     1.71 ; given_sigma_ion_O2_to_O_cm2( 130) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 131) =     1.73 ; given_sigma_ion_O2_to_O_cm2( 131) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 132) =     1.74 ; given_sigma_ion_O2_to_O_cm2( 132) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 133) =     1.77 ; given_sigma_ion_O2_to_O_cm2( 133) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 134) =     1.79 ; given_sigma_ion_O2_to_O_cm2( 134) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 135) =     1.80 ; given_sigma_ion_O2_to_O_cm2( 135) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 136) =     1.82 ; given_sigma_ion_O2_to_O_cm2( 136) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 137) =     1.83 ; given_sigma_ion_O2_to_O_cm2( 137) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 138) =     1.85 ; given_sigma_ion_O2_to_O_cm2( 138) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 139) =     1.86 ; given_sigma_ion_O2_to_O_cm2( 139) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 140) =     1.87 ; given_sigma_ion_O2_to_O_cm2( 140) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 141) =     1.88 ; given_sigma_ion_O2_to_O_cm2( 141) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 142) =     1.91 ; given_sigma_ion_O2_to_O_cm2( 142) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 143) =     1.92 ; given_sigma_ion_O2_to_O_cm2( 143) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 144) =     1.96 ; given_sigma_ion_O2_to_O_cm2( 144) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 145) =     1.98 ; given_sigma_ion_O2_to_O_cm2( 145) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 146) =     1.99 ; given_sigma_ion_O2_to_O_cm2( 146) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 147) =     2.00 ; given_sigma_ion_O2_to_O_cm2( 147) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 148) =     2.02 ; given_sigma_ion_O2_to_O_cm2( 148) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 149) =     2.04 ; given_sigma_ion_O2_to_O_cm2( 149) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 150) =     2.05 ; given_sigma_ion_O2_to_O_cm2( 150) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 151) =     2.06 ; given_sigma_ion_O2_to_O_cm2( 151) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 152) =     2.08 ; given_sigma_ion_O2_to_O_cm2( 152) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 153) =     2.10 ; given_sigma_ion_O2_to_O_cm2( 153) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 154) =     2.11 ; given_sigma_ion_O2_to_O_cm2( 154) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 155) =     2.13 ; given_sigma_ion_O2_to_O_cm2( 155) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 156) =     2.15 ; given_sigma_ion_O2_to_O_cm2( 156) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 157) =     2.16 ; given_sigma_ion_O2_to_O_cm2( 157) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 158) =     2.18 ; given_sigma_ion_O2_to_O_cm2( 158) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 159) =     2.18 ; given_sigma_ion_O2_to_O_cm2( 159) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 160) =     2.20 ; given_sigma_ion_O2_to_O_cm2( 160) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 161) =     2.21 ; given_sigma_ion_O2_to_O_cm2( 161) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 162) =     2.24 ; given_sigma_ion_O2_to_O_cm2( 162) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 163) =     2.25 ; given_sigma_ion_O2_to_O_cm2( 163) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 164) =     2.27 ; given_sigma_ion_O2_to_O_cm2( 164) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 165) =     2.29 ; given_sigma_ion_O2_to_O_cm2( 165) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 166) =     2.31 ; given_sigma_ion_O2_to_O_cm2( 166) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 167) =     2.34 ; given_sigma_ion_O2_to_O_cm2( 167) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 168) =     2.37 ; given_sigma_ion_O2_to_O_cm2( 168) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 169) =     2.41 ; given_sigma_ion_O2_to_O_cm2( 169) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 170) =     2.42 ; given_sigma_ion_O2_to_O_cm2( 170) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 171) =     2.43 ; given_sigma_ion_O2_to_O_cm2( 171) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 172) =     2.44 ; given_sigma_ion_O2_to_O_cm2( 172) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 173) =     2.46 ; given_sigma_ion_O2_to_O_cm2( 173) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 174) =     2.47 ; given_sigma_ion_O2_to_O_cm2( 174) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 175) =     2.50 ; given_sigma_ion_O2_to_O_cm2( 175) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 176) =     2.52 ; given_sigma_ion_O2_to_O_cm2( 176) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 177) =     2.57 ; given_sigma_ion_O2_to_O_cm2( 177) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 178) =     2.59 ; given_sigma_ion_O2_to_O_cm2( 178) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 179) =     2.60 ; given_sigma_ion_O2_to_O_cm2( 179) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 180) =     2.65 ; given_sigma_ion_O2_to_O_cm2( 180) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 181) =     2.66 ; given_sigma_ion_O2_to_O_cm2( 181) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 182) =     2.68 ; given_sigma_ion_O2_to_O_cm2( 182) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 183) =     2.73 ; given_sigma_ion_O2_to_O_cm2( 183) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 184) =     2.75 ; given_sigma_ion_O2_to_O_cm2( 184) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 185) =     2.78 ; given_sigma_ion_O2_to_O_cm2( 185) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 186) =     2.79 ; given_sigma_ion_O2_to_O_cm2( 186) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 187) =     2.81 ; given_sigma_ion_O2_to_O_cm2( 187) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 188) =     2.89 ; given_sigma_ion_O2_to_O_cm2( 188) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 189) =     2.94 ; given_sigma_ion_O2_to_O_cm2( 189) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 190) =     2.95 ; given_sigma_ion_O2_to_O_cm2( 190) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 191) =     2.96 ; given_sigma_ion_O2_to_O_cm2( 191) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 192) =     3.02 ; given_sigma_ion_O2_to_O_cm2( 192) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 193) =     3.04 ; given_sigma_ion_O2_to_O_cm2( 193) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 194) =     3.08 ; given_sigma_ion_O2_to_O_cm2( 194) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 195) =     3.10 ; given_sigma_ion_O2_to_O_cm2( 195) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 196) =     3.15 ; given_sigma_ion_O2_to_O_cm2( 196) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 197) =     1.13 ; given_sigma_ion_O2_to_O_cm2( 197) =    2.11
  given_sigma_ion_O2_to_O2_cm2( 198) =     1.15 ; given_sigma_ion_O2_to_O_cm2( 198) =    2.12
  given_sigma_ion_O2_to_O2_cm2( 199) =     1.18 ; given_sigma_ion_O2_to_O_cm2( 199) =    2.14
  given_sigma_ion_O2_to_O2_cm2( 200) =     1.21 ; given_sigma_ion_O2_to_O_cm2( 200) =    2.16
  given_sigma_ion_O2_to_O2_cm2( 201) =     1.24 ; given_sigma_ion_O2_to_O_cm2( 201) =    2.19
  given_sigma_ion_O2_to_O2_cm2( 202) =     1.28 ; given_sigma_ion_O2_to_O_cm2( 202) =    2.21
  given_sigma_ion_O2_to_O2_cm2( 203) =     1.30 ; given_sigma_ion_O2_to_O_cm2( 203) =    2.22
  given_sigma_ion_O2_to_O2_cm2( 204) =     1.35 ; given_sigma_ion_O2_to_O_cm2( 204) =    2.25
  given_sigma_ion_O2_to_O2_cm2( 205) =     1.47 ; given_sigma_ion_O2_to_O_cm2( 205) =    2.32
  given_sigma_ion_O2_to_O2_cm2( 206) =     1.57 ; given_sigma_ion_O2_to_O_cm2( 206) =    2.37
  given_sigma_ion_O2_to_O2_cm2( 207) =     1.58 ; given_sigma_ion_O2_to_O_cm2( 207) =    2.37
  given_sigma_ion_O2_to_O2_cm2( 208) =     1.60 ; given_sigma_ion_O2_to_O_cm2( 208) =    2.38
  given_sigma_ion_O2_to_O2_cm2( 209) =     1.61 ; given_sigma_ion_O2_to_O_cm2( 209) =    2.38
  given_sigma_ion_O2_to_O2_cm2( 210) =     1.63 ; given_sigma_ion_O2_to_O_cm2( 210) =    2.39
  given_sigma_ion_O2_to_O2_cm2( 211) =     1.64 ; given_sigma_ion_O2_to_O_cm2( 211) =    2.39
  given_sigma_ion_O2_to_O2_cm2( 212) =     1.84 ; given_sigma_ion_O2_to_O_cm2( 212) =    2.46
  given_sigma_ion_O2_to_O2_cm2( 213) =     1.91 ; given_sigma_ion_O2_to_O_cm2( 213) =    2.48
  given_sigma_ion_O2_to_O2_cm2( 214) =     1.92 ; given_sigma_ion_O2_to_O_cm2( 214) =    2.48
  given_sigma_ion_O2_to_O2_cm2( 215) =     2.00 ; given_sigma_ion_O2_to_O_cm2( 215) =    2.50
  given_sigma_ion_O2_to_O2_cm2( 216) =     2.12 ; given_sigma_ion_O2_to_O_cm2( 216) =    2.53
  given_sigma_ion_O2_to_O2_cm2( 217) =     2.17 ; given_sigma_ion_O2_to_O_cm2( 217) =    2.56
  given_sigma_ion_O2_to_O2_cm2( 218) =     2.30 ; given_sigma_ion_O2_to_O_cm2( 218) =    2.63
  given_sigma_ion_O2_to_O2_cm2( 219) =     2.34 ; given_sigma_ion_O2_to_O_cm2( 219) =    2.65
  given_sigma_ion_O2_to_O2_cm2( 220) =     2.38 ; given_sigma_ion_O2_to_O_cm2( 220) =    2.67
  given_sigma_ion_O2_to_O2_cm2( 221) =     2.49 ; given_sigma_ion_O2_to_O_cm2( 221) =    2.72
  given_sigma_ion_O2_to_O2_cm2( 222) =     2.56 ; given_sigma_ion_O2_to_O_cm2( 222) =    2.76
  given_sigma_ion_O2_to_O2_cm2( 223) =     2.56 ; given_sigma_ion_O2_to_O_cm2( 223) =    2.76
  given_sigma_ion_O2_to_O2_cm2( 224) =     2.66 ; given_sigma_ion_O2_to_O_cm2( 224) =    2.82
  given_sigma_ion_O2_to_O2_cm2( 225) =     2.69 ; given_sigma_ion_O2_to_O_cm2( 225) =    2.84
  given_sigma_ion_O2_to_O2_cm2( 226) =     2.75 ; given_sigma_ion_O2_to_O_cm2( 226) =    2.87
  given_sigma_ion_O2_to_O2_cm2( 227) =     2.78 ; given_sigma_ion_O2_to_O_cm2( 227) =    2.89
  given_sigma_ion_O2_to_O2_cm2( 228) =     2.91 ; given_sigma_ion_O2_to_O_cm2( 228) =    2.96
  given_sigma_ion_O2_to_O2_cm2( 229) =     2.94 ; given_sigma_ion_O2_to_O_cm2( 229) =    2.98
  given_sigma_ion_O2_to_O2_cm2( 230) =     3.01 ; given_sigma_ion_O2_to_O_cm2( 230) =    3.01
  given_sigma_ion_O2_to_O2_cm2( 231) =     3.03 ; given_sigma_ion_O2_to_O_cm2( 231) =    3.02
  given_sigma_ion_O2_to_O2_cm2( 232) =     3.21 ; given_sigma_ion_O2_to_O_cm2( 232) =    3.09
  given_sigma_ion_O2_to_O2_cm2( 233) =     3.26 ; given_sigma_ion_O2_to_O_cm2( 233) =    3.11
  given_sigma_ion_O2_to_O2_cm2( 234) =     3.27 ; given_sigma_ion_O2_to_O_cm2( 234) =    3.11
  given_sigma_ion_O2_to_O2_cm2( 235) =     3.38 ; given_sigma_ion_O2_to_O_cm2( 235) =    3.15
  given_sigma_ion_O2_to_O2_cm2( 236) =     3.41 ; given_sigma_ion_O2_to_O_cm2( 236) =    3.16
  given_sigma_ion_O2_to_O2_cm2( 237) =     3.43 ; given_sigma_ion_O2_to_O_cm2( 237) =    3.17
  given_sigma_ion_O2_to_O2_cm2( 238) =     3.45 ; given_sigma_ion_O2_to_O_cm2( 238) =    3.18
  given_sigma_ion_O2_to_O2_cm2( 239) =     3.49 ; given_sigma_ion_O2_to_O_cm2( 239) =    3.19
  given_sigma_ion_O2_to_O2_cm2( 240) =     3.50 ; given_sigma_ion_O2_to_O_cm2( 240) =    3.20
  given_sigma_ion_O2_to_O2_cm2( 241) =     3.56 ; given_sigma_ion_O2_to_O_cm2( 241) =    3.20
  given_sigma_ion_O2_to_O2_cm2( 242) =     3.62 ; given_sigma_ion_O2_to_O_cm2( 242) =    3.21
  given_sigma_ion_O2_to_O2_cm2( 243) =     3.67 ; given_sigma_ion_O2_to_O_cm2( 243) =    3.22
  given_sigma_ion_O2_to_O2_cm2( 244) =     3.67 ; given_sigma_ion_O2_to_O_cm2( 244) =    3.22
  given_sigma_ion_O2_to_O2_cm2( 245) =     3.76 ; given_sigma_ion_O2_to_O_cm2( 245) =    3.23
  given_sigma_ion_O2_to_O2_cm2( 246) =     3.78 ; given_sigma_ion_O2_to_O_cm2( 246) =    3.23
  given_sigma_ion_O2_to_O2_cm2( 247) =     3.80 ; given_sigma_ion_O2_to_O_cm2( 247) =    3.23
  given_sigma_ion_O2_to_O2_cm2( 248) =     3.81 ; given_sigma_ion_O2_to_O_cm2( 248) =    3.24
  given_sigma_ion_O2_to_O2_cm2( 249) =     3.90 ; given_sigma_ion_O2_to_O_cm2( 249) =    3.25
  given_sigma_ion_O2_to_O2_cm2( 250) =     3.91 ; given_sigma_ion_O2_to_O_cm2( 250) =    3.25
  given_sigma_ion_O2_to_O2_cm2( 251) =     3.95 ; given_sigma_ion_O2_to_O_cm2( 251) =    3.25
  given_sigma_ion_O2_to_O2_cm2( 252) =     4.03 ; given_sigma_ion_O2_to_O_cm2( 252) =    3.26
  given_sigma_ion_O2_to_O2_cm2( 253) =     4.05 ; given_sigma_ion_O2_to_O_cm2( 253) =    3.26
  given_sigma_ion_O2_to_O2_cm2( 254) =     4.07 ; given_sigma_ion_O2_to_O_cm2( 254) =    3.26
  given_sigma_ion_O2_to_O2_cm2( 255) =     4.09 ; given_sigma_ion_O2_to_O_cm2( 255) =    3.27
  given_sigma_ion_O2_to_O2_cm2( 256) =     4.11 ; given_sigma_ion_O2_to_O_cm2( 256) =    3.28
  given_sigma_ion_O2_to_O2_cm2( 257) =     4.13 ; given_sigma_ion_O2_to_O_cm2( 257) =    3.29
  given_sigma_ion_O2_to_O2_cm2( 258) =     4.18 ; given_sigma_ion_O2_to_O_cm2( 258) =    3.32
  given_sigma_ion_O2_to_O2_cm2( 259) =     4.20 ; given_sigma_ion_O2_to_O_cm2( 259) =    3.33
  given_sigma_ion_O2_to_O2_cm2( 260) =     4.25 ; given_sigma_ion_O2_to_O_cm2( 260) =    3.36
  given_sigma_ion_O2_to_O2_cm2( 261) =     4.28 ; given_sigma_ion_O2_to_O_cm2( 261) =    3.38
  given_sigma_ion_O2_to_O2_cm2( 262) =     4.29 ; given_sigma_ion_O2_to_O_cm2( 262) =    3.38
  given_sigma_ion_O2_to_O2_cm2( 263) =     4.31 ; given_sigma_ion_O2_to_O_cm2( 263) =    3.39
  given_sigma_ion_O2_to_O2_cm2( 264) =     4.32 ; given_sigma_ion_O2_to_O_cm2( 264) =    3.40
  given_sigma_ion_O2_to_O2_cm2( 265) =     4.33 ; given_sigma_ion_O2_to_O_cm2( 265) =    3.41
  given_sigma_ion_O2_to_O2_cm2( 266) =     4.34 ; given_sigma_ion_O2_to_O_cm2( 266) =    3.41
  given_sigma_ion_O2_to_O2_cm2( 267) =     4.42 ; given_sigma_ion_O2_to_O_cm2( 267) =    3.45
  given_sigma_ion_O2_to_O2_cm2( 268) =     4.43 ; given_sigma_ion_O2_to_O_cm2( 268) =    3.46
  given_sigma_ion_O2_to_O2_cm2( 269) =     4.49 ; given_sigma_ion_O2_to_O_cm2( 269) =    3.49
  given_sigma_ion_O2_to_O2_cm2( 270) =     4.51 ; given_sigma_ion_O2_to_O_cm2( 270) =    3.50
  given_sigma_ion_O2_to_O2_cm2( 271) =     4.53 ; given_sigma_ion_O2_to_O_cm2( 271) =    3.51
  given_sigma_ion_O2_to_O2_cm2( 272) =     4.60 ; given_sigma_ion_O2_to_O_cm2( 272) =    3.55
  given_sigma_ion_O2_to_O2_cm2( 273) =     4.63 ; given_sigma_ion_O2_to_O_cm2( 273) =    3.57
  given_sigma_ion_O2_to_O2_cm2( 274) =     4.65 ; given_sigma_ion_O2_to_O_cm2( 274) =    3.57
  given_sigma_ion_O2_to_O2_cm2( 275) =     4.66 ; given_sigma_ion_O2_to_O_cm2( 275) =    3.58
  given_sigma_ion_O2_to_O2_cm2( 276) =     4.72 ; given_sigma_ion_O2_to_O_cm2( 276) =    3.60
  given_sigma_ion_O2_to_O2_cm2( 277) =     4.74 ; given_sigma_ion_O2_to_O_cm2( 277) =    3.61
  given_sigma_ion_O2_to_O2_cm2( 278) =     4.77 ; given_sigma_ion_O2_to_O_cm2( 278) =    3.63
  given_sigma_ion_O2_to_O2_cm2( 279) =     4.84 ; given_sigma_ion_O2_to_O_cm2( 279) =    3.66
  given_sigma_ion_O2_to_O2_cm2( 280) =     4.85 ; given_sigma_ion_O2_to_O_cm2( 280) =    3.66
  given_sigma_ion_O2_to_O2_cm2( 281) =     4.92 ; given_sigma_ion_O2_to_O_cm2( 281) =    3.70
  given_sigma_ion_O2_to_O2_cm2( 282) =     4.92 ; given_sigma_ion_O2_to_O_cm2( 282) =    3.70
  given_sigma_ion_O2_to_O2_cm2( 283) =     4.96 ; given_sigma_ion_O2_to_O_cm2( 283) =    3.71
  given_sigma_ion_O2_to_O2_cm2( 284) =     5.02 ; given_sigma_ion_O2_to_O_cm2( 284) =    3.74
  given_sigma_ion_O2_to_O2_cm2( 285) =     5.09 ; given_sigma_ion_O2_to_O_cm2( 285) =    3.77
  given_sigma_ion_O2_to_O2_cm2( 286) =     5.15 ; given_sigma_ion_O2_to_O_cm2( 286) =    3.79
  given_sigma_ion_O2_to_O2_cm2( 287) =     5.20 ; given_sigma_ion_O2_to_O_cm2( 287) =    3.80
  given_sigma_ion_O2_to_O2_cm2( 288) =     5.24 ; given_sigma_ion_O2_to_O_cm2( 288) =    3.81
  given_sigma_ion_O2_to_O2_cm2( 289) =     5.30 ; given_sigma_ion_O2_to_O_cm2( 289) =    3.83
  given_sigma_ion_O2_to_O2_cm2( 290) =     5.33 ; given_sigma_ion_O2_to_O_cm2( 290) =    3.84
  given_sigma_ion_O2_to_O2_cm2( 291) =     5.36 ; given_sigma_ion_O2_to_O_cm2( 291) =    3.85
  given_sigma_ion_O2_to_O2_cm2( 292) =     5.44 ; given_sigma_ion_O2_to_O_cm2( 292) =    3.87
  given_sigma_ion_O2_to_O2_cm2( 293) =     5.44 ; given_sigma_ion_O2_to_O_cm2( 293) =    3.87
  given_sigma_ion_O2_to_O2_cm2( 294) =     5.46 ; given_sigma_ion_O2_to_O_cm2( 294) =    3.87
  given_sigma_ion_O2_to_O2_cm2( 295) =     5.51 ; given_sigma_ion_O2_to_O_cm2( 295) =    3.88
  given_sigma_ion_O2_to_O2_cm2( 296) =     5.55 ; given_sigma_ion_O2_to_O_cm2( 296) =    3.90
  given_sigma_ion_O2_to_O2_cm2( 297) =     5.63 ; given_sigma_ion_O2_to_O_cm2( 297) =    3.91
  given_sigma_ion_O2_to_O2_cm2( 298) =     5.64 ; given_sigma_ion_O2_to_O_cm2( 298) =    3.92
  given_sigma_ion_O2_to_O2_cm2( 299) =     5.65 ; given_sigma_ion_O2_to_O_cm2( 299) =    3.92
  given_sigma_ion_O2_to_O2_cm2( 300) =     5.71 ; given_sigma_ion_O2_to_O_cm2( 300) =    3.95
  given_sigma_ion_O2_to_O2_cm2( 301) =     5.75 ; given_sigma_ion_O2_to_O_cm2( 301) =    3.97
  given_sigma_ion_O2_to_O2_cm2( 302) =     5.83 ; given_sigma_ion_O2_to_O_cm2( 302) =    4.01
  given_sigma_ion_O2_to_O2_cm2( 303) =     5.88 ; given_sigma_ion_O2_to_O_cm2( 303) =    4.04
  given_sigma_ion_O2_to_O2_cm2( 304) =     5.89 ; given_sigma_ion_O2_to_O_cm2( 304) =    4.04
  given_sigma_ion_O2_to_O2_cm2( 305) =     5.90 ; given_sigma_ion_O2_to_O_cm2( 305) =    4.04
  given_sigma_ion_O2_to_O2_cm2( 306) =     5.99 ; given_sigma_ion_O2_to_O_cm2( 306) =    4.09
  given_sigma_ion_O2_to_O2_cm2( 307) =     5.99 ; given_sigma_ion_O2_to_O_cm2( 307) =    4.09
  given_sigma_ion_O2_to_O2_cm2( 308) =     6.05 ; given_sigma_ion_O2_to_O_cm2( 308) =    4.12
  given_sigma_ion_O2_to_O2_cm2( 309) =     6.10 ; given_sigma_ion_O2_to_O_cm2( 309) =    4.14
  given_sigma_ion_O2_to_O2_cm2( 310) =     6.14 ; given_sigma_ion_O2_to_O_cm2( 310) =    4.16
  given_sigma_ion_O2_to_O2_cm2( 311) =     6.22 ; given_sigma_ion_O2_to_O_cm2( 311) =    4.18
  given_sigma_ion_O2_to_O2_cm2( 312) =     6.23 ; given_sigma_ion_O2_to_O_cm2( 312) =    4.18
  given_sigma_ion_O2_to_O2_cm2( 313) =     6.25 ; given_sigma_ion_O2_to_O_cm2( 313) =    4.19
  given_sigma_ion_O2_to_O2_cm2( 314) =     6.35 ; given_sigma_ion_O2_to_O_cm2( 314) =    4.21
  given_sigma_ion_O2_to_O2_cm2( 315) =     6.38 ; given_sigma_ion_O2_to_O_cm2( 315) =    4.22
  given_sigma_ion_O2_to_O2_cm2( 316) =     6.44 ; given_sigma_ion_O2_to_O_cm2( 316) =    4.24
  given_sigma_ion_O2_to_O2_cm2( 317) =     6.46 ; given_sigma_ion_O2_to_O_cm2( 317) =    4.24
  given_sigma_ion_O2_to_O2_cm2( 318) =     6.46 ; given_sigma_ion_O2_to_O_cm2( 318) =    4.24
  given_sigma_ion_O2_to_O2_cm2( 319) =     6.59 ; given_sigma_ion_O2_to_O_cm2( 319) =    4.27
  given_sigma_ion_O2_to_O2_cm2( 320) =     6.60 ; given_sigma_ion_O2_to_O_cm2( 320) =    4.28
  given_sigma_ion_O2_to_O2_cm2( 321) =     6.62 ; given_sigma_ion_O2_to_O_cm2( 321) =    4.28
  given_sigma_ion_O2_to_O2_cm2( 322) =     6.70 ; given_sigma_ion_O2_to_O_cm2( 322) =    4.30
  given_sigma_ion_O2_to_O2_cm2( 323) =     6.75 ; given_sigma_ion_O2_to_O_cm2( 323) =    4.31
  given_sigma_ion_O2_to_O2_cm2( 324) =     6.78 ; given_sigma_ion_O2_to_O_cm2( 324) =    4.32
  given_sigma_ion_O2_to_O2_cm2( 325) =     6.82 ; given_sigma_ion_O2_to_O_cm2( 325) =    4.34
  given_sigma_ion_O2_to_O2_cm2( 326) =     6.87 ; given_sigma_ion_O2_to_O_cm2( 326) =    4.36
  given_sigma_ion_O2_to_O2_cm2( 327) =     6.93 ; given_sigma_ion_O2_to_O_cm2( 327) =    4.38
  given_sigma_ion_O2_to_O2_cm2( 328) =     7.00 ; given_sigma_ion_O2_to_O_cm2( 328) =    4.41
  given_sigma_ion_O2_to_O2_cm2( 329) =     7.03 ; given_sigma_ion_O2_to_O_cm2( 329) =    4.42
  given_sigma_ion_O2_to_O2_cm2( 330) =     7.06 ; given_sigma_ion_O2_to_O_cm2( 330) =    4.44
  given_sigma_ion_O2_to_O2_cm2( 331) =     7.09 ; given_sigma_ion_O2_to_O_cm2( 331) =    4.45
  given_sigma_ion_O2_to_O2_cm2( 332) =     7.20 ; given_sigma_ion_O2_to_O_cm2( 332) =    4.49
  given_sigma_ion_O2_to_O2_cm2( 333) =     7.26 ; given_sigma_ion_O2_to_O_cm2( 333) =    4.51
  given_sigma_ion_O2_to_O2_cm2( 334) =     7.34 ; given_sigma_ion_O2_to_O_cm2( 334) =    4.55
  given_sigma_ion_O2_to_O2_cm2( 335) =     7.35 ; given_sigma_ion_O2_to_O_cm2( 335) =    4.55
  given_sigma_ion_O2_to_O2_cm2( 336) =     7.40 ; given_sigma_ion_O2_to_O_cm2( 336) =    4.56
  given_sigma_ion_O2_to_O2_cm2( 337) =     7.48 ; given_sigma_ion_O2_to_O_cm2( 337) =    4.58
  given_sigma_ion_O2_to_O2_cm2( 338) =     7.57 ; given_sigma_ion_O2_to_O_cm2( 338) =    4.60
  given_sigma_ion_O2_to_O2_cm2( 339) =     7.64 ; given_sigma_ion_O2_to_O_cm2( 339) =    4.61
  given_sigma_ion_O2_to_O2_cm2( 340) =     7.72 ; given_sigma_ion_O2_to_O_cm2( 340) =    4.63
  given_sigma_ion_O2_to_O2_cm2( 341) =     7.79 ; given_sigma_ion_O2_to_O_cm2( 341) =    4.64
  given_sigma_ion_O2_to_O2_cm2( 342) =     7.81 ; given_sigma_ion_O2_to_O_cm2( 342) =    4.65
  given_sigma_ion_O2_to_O2_cm2( 343) =     7.87 ; given_sigma_ion_O2_to_O_cm2( 343) =    4.66
  given_sigma_ion_O2_to_O2_cm2( 344) =     7.89 ; given_sigma_ion_O2_to_O_cm2( 344) =    4.66
  given_sigma_ion_O2_to_O2_cm2( 345) =     7.95 ; given_sigma_ion_O2_to_O_cm2( 345) =    4.67
  given_sigma_ion_O2_to_O2_cm2( 346) =     8.04 ; given_sigma_ion_O2_to_O_cm2( 346) =    4.69
  given_sigma_ion_O2_to_O2_cm2( 347) =     8.10 ; given_sigma_ion_O2_to_O_cm2( 347) =    4.70
  given_sigma_ion_O2_to_O2_cm2( 348) =     8.16 ; given_sigma_ion_O2_to_O_cm2( 348) =    4.74
  given_sigma_ion_O2_to_O2_cm2( 349) =     8.21 ; given_sigma_ion_O2_to_O_cm2( 349) =    4.76
  given_sigma_ion_O2_to_O2_cm2( 350) =     8.22 ; given_sigma_ion_O2_to_O_cm2( 350) =    4.77
  given_sigma_ion_O2_to_O2_cm2( 351) =     8.32 ; given_sigma_ion_O2_to_O_cm2( 351) =    4.82
  given_sigma_ion_O2_to_O2_cm2( 352) =     8.35 ; given_sigma_ion_O2_to_O_cm2( 352) =    4.85
  given_sigma_ion_O2_to_O2_cm2( 353) =     8.38 ; given_sigma_ion_O2_to_O_cm2( 353) =    4.87
  given_sigma_ion_O2_to_O2_cm2( 354) =     8.46 ; given_sigma_ion_O2_to_O_cm2( 354) =    4.91
  given_sigma_ion_O2_to_O2_cm2( 355) =     8.48 ; given_sigma_ion_O2_to_O_cm2( 355) =    4.92
  given_sigma_ion_O2_to_O2_cm2( 356) =     8.51 ; given_sigma_ion_O2_to_O_cm2( 356) =    4.94
  given_sigma_ion_O2_to_O2_cm2( 357) =     8.52 ; given_sigma_ion_O2_to_O_cm2( 357) =    4.94
  given_sigma_ion_O2_to_O2_cm2( 358) =     8.57 ; given_sigma_ion_O2_to_O_cm2( 358) =    4.97
  given_sigma_ion_O2_to_O2_cm2( 359) =     8.64 ; given_sigma_ion_O2_to_O_cm2( 359) =    5.01
  given_sigma_ion_O2_to_O2_cm2( 360) =     8.67 ; given_sigma_ion_O2_to_O_cm2( 360) =    5.03
  given_sigma_ion_O2_to_O2_cm2( 361) =     8.72 ; given_sigma_ion_O2_to_O_cm2( 361) =    5.09
  given_sigma_ion_O2_to_O2_cm2( 362) =     8.81 ; given_sigma_ion_O2_to_O_cm2( 362) =    5.19
  given_sigma_ion_O2_to_O2_cm2( 363) =     8.86 ; given_sigma_ion_O2_to_O_cm2( 363) =    5.26
  given_sigma_ion_O2_to_O2_cm2( 364) =     8.89 ; given_sigma_ion_O2_to_O_cm2( 364) =    5.29
  given_sigma_ion_O2_to_O2_cm2( 365) =     8.90 ; given_sigma_ion_O2_to_O_cm2( 365) =    5.30
  given_sigma_ion_O2_to_O2_cm2( 366) =     9.09 ; given_sigma_ion_O2_to_O_cm2( 366) =    5.56
  given_sigma_ion_O2_to_O2_cm2( 367) =     9.11 ; given_sigma_ion_O2_to_O_cm2( 367) =    5.59
  given_sigma_ion_O2_to_O2_cm2( 368) =     9.14 ; given_sigma_ion_O2_to_O_cm2( 368) =    5.60
  given_sigma_ion_O2_to_O2_cm2( 369) =     9.22 ; given_sigma_ion_O2_to_O_cm2( 369) =    5.64
  given_sigma_ion_O2_to_O2_cm2( 370) =     9.25 ; given_sigma_ion_O2_to_O_cm2( 370) =    5.66
  given_sigma_ion_O2_to_O2_cm2( 371) =     9.34 ; given_sigma_ion_O2_to_O_cm2( 371) =    5.70
  given_sigma_ion_O2_to_O2_cm2( 372) =     9.38 ; given_sigma_ion_O2_to_O_cm2( 372) =    5.72
  given_sigma_ion_O2_to_O2_cm2( 373) =     9.40 ; given_sigma_ion_O2_to_O_cm2( 373) =    5.73
  given_sigma_ion_O2_to_O2_cm2( 374) =     9.42 ; given_sigma_ion_O2_to_O_cm2( 374) =    5.74
  given_sigma_ion_O2_to_O2_cm2( 375) =     9.44 ; given_sigma_ion_O2_to_O_cm2( 375) =    5.75
  given_sigma_ion_O2_to_O2_cm2( 376) =     9.48 ; given_sigma_ion_O2_to_O_cm2( 376) =    5.77
  given_sigma_ion_O2_to_O2_cm2( 377) =     9.49 ; given_sigma_ion_O2_to_O_cm2( 377) =    5.77
  given_sigma_ion_O2_to_O2_cm2( 378) =     9.50 ; given_sigma_ion_O2_to_O_cm2( 378) =    5.78
  given_sigma_ion_O2_to_O2_cm2( 379) =     9.56 ; given_sigma_ion_O2_to_O_cm2( 379) =    5.81
  given_sigma_ion_O2_to_O2_cm2( 380) =     9.65 ; given_sigma_ion_O2_to_O_cm2( 380) =    5.85
  given_sigma_ion_O2_to_O2_cm2( 381) =     9.70 ; given_sigma_ion_O2_to_O_cm2( 381) =    5.86
  given_sigma_ion_O2_to_O2_cm2( 382) =     9.73 ; given_sigma_ion_O2_to_O_cm2( 382) =    5.87
  given_sigma_ion_O2_to_O2_cm2( 383) =     9.89 ; given_sigma_ion_O2_to_O_cm2( 383) =    5.90
  given_sigma_ion_O2_to_O2_cm2( 384) =     9.94 ; given_sigma_ion_O2_to_O_cm2( 384) =    5.91
  given_sigma_ion_O2_to_O2_cm2( 385) =     9.98 ; given_sigma_ion_O2_to_O_cm2( 385) =    5.92
  given_sigma_ion_O2_to_O2_cm2( 386) =     9.99 ; given_sigma_ion_O2_to_O_cm2( 386) =    5.92
  given_sigma_ion_O2_to_O2_cm2( 387) =    10.14 ; given_sigma_ion_O2_to_O_cm2( 387) =    5.94
  given_sigma_ion_O2_to_O2_cm2( 388) =    10.19 ; given_sigma_ion_O2_to_O_cm2( 388) =    5.95
  given_sigma_ion_O2_to_O2_cm2( 389) =    10.24 ; given_sigma_ion_O2_to_O_cm2( 389) =    5.96
  given_sigma_ion_O2_to_O2_cm2( 390) =    10.27 ; given_sigma_ion_O2_to_O_cm2( 390) =    5.97
  given_sigma_ion_O2_to_O2_cm2( 391) =    10.31 ; given_sigma_ion_O2_to_O_cm2( 391) =    5.97
  given_sigma_ion_O2_to_O2_cm2( 392) =    10.33 ; given_sigma_ion_O2_to_O_cm2( 392) =    5.97
  given_sigma_ion_O2_to_O2_cm2( 393) =    10.36 ; given_sigma_ion_O2_to_O_cm2( 393) =    5.98
  given_sigma_ion_O2_to_O2_cm2( 394) =    10.46 ; given_sigma_ion_O2_to_O_cm2( 394) =    5.99
  given_sigma_ion_O2_to_O2_cm2( 395) =    10.48 ; given_sigma_ion_O2_to_O_cm2( 395) =    5.99
  given_sigma_ion_O2_to_O2_cm2( 396) =    10.51 ; given_sigma_ion_O2_to_O_cm2( 396) =    6.00
  given_sigma_ion_O2_to_O2_cm2( 397) =    10.66 ; given_sigma_ion_O2_to_O_cm2( 397) =    6.02
  given_sigma_ion_O2_to_O2_cm2( 398) =    10.68 ; given_sigma_ion_O2_to_O_cm2( 398) =    6.02
  given_sigma_ion_O2_to_O2_cm2( 399) =    10.72 ; given_sigma_ion_O2_to_O_cm2( 399) =    6.00
  given_sigma_ion_O2_to_O2_cm2( 400) =    10.86 ; given_sigma_ion_O2_to_O_cm2( 400) =    5.94
  given_sigma_ion_O2_to_O2_cm2( 401) =    10.88 ; given_sigma_ion_O2_to_O_cm2( 401) =    5.93
  given_sigma_ion_O2_to_O2_cm2( 402) =    11.22 ; given_sigma_ion_O2_to_O_cm2( 402) =    5.78
  given_sigma_ion_O2_to_O2_cm2( 403) =    11.45 ; given_sigma_ion_O2_to_O_cm2( 403) =    5.65
  given_sigma_ion_O2_to_O2_cm2( 404) =    11.50 ; given_sigma_ion_O2_to_O_cm2( 404) =    5.62
  given_sigma_ion_O2_to_O2_cm2( 405) =    11.51 ; given_sigma_ion_O2_to_O_cm2( 405) =    5.62
  given_sigma_ion_O2_to_O2_cm2( 406) =    11.64 ; given_sigma_ion_O2_to_O_cm2( 406) =    5.54
  given_sigma_ion_O2_to_O2_cm2( 407) =    11.68 ; given_sigma_ion_O2_to_O_cm2( 407) =    5.52
  given_sigma_ion_O2_to_O2_cm2( 408) =    11.69 ; given_sigma_ion_O2_to_O_cm2( 408) =    5.51
  given_sigma_ion_O2_to_O2_cm2( 409) =    11.71 ; given_sigma_ion_O2_to_O_cm2( 409) =    5.50
  given_sigma_ion_O2_to_O2_cm2( 410) =    11.89 ; given_sigma_ion_O2_to_O_cm2( 410) =    5.36
  given_sigma_ion_O2_to_O2_cm2( 411) =    12.10 ; given_sigma_ion_O2_to_O_cm2( 411) =    5.20
  given_sigma_ion_O2_to_O2_cm2( 412) =    12.33 ; given_sigma_ion_O2_to_O_cm2( 412) =    5.07
  given_sigma_ion_O2_to_O2_cm2( 413) =    12.35 ; given_sigma_ion_O2_to_O_cm2( 413) =    5.06
  given_sigma_ion_O2_to_O2_cm2( 414) =    12.56 ; given_sigma_ion_O2_to_O_cm2( 414) =    4.94
  given_sigma_ion_O2_to_O2_cm2( 415) =    12.73 ; given_sigma_ion_O2_to_O_cm2( 415) =    4.92
  given_sigma_ion_O2_to_O2_cm2( 416) =    12.73 ; given_sigma_ion_O2_to_O_cm2( 416) =    4.92
  given_sigma_ion_O2_to_O2_cm2( 417) =    12.75 ; given_sigma_ion_O2_to_O_cm2( 417) =    4.92
  given_sigma_ion_O2_to_O2_cm2( 418) =    12.81 ; given_sigma_ion_O2_to_O_cm2( 418) =    4.92
  given_sigma_ion_O2_to_O2_cm2( 419) =    12.89 ; given_sigma_ion_O2_to_O_cm2( 419) =    4.91
  given_sigma_ion_O2_to_O2_cm2( 420) =    12.89 ; given_sigma_ion_O2_to_O_cm2( 420) =    4.91
  given_sigma_ion_O2_to_O2_cm2( 421) =    13.05 ; given_sigma_ion_O2_to_O_cm2( 421) =    4.83
  given_sigma_ion_O2_to_O2_cm2( 422) =    13.14 ; given_sigma_ion_O2_to_O_cm2( 422) =    4.78
  given_sigma_ion_O2_to_O2_cm2( 423) =    13.30 ; given_sigma_ion_O2_to_O_cm2( 423) =    4.70
  given_sigma_ion_O2_to_O2_cm2( 424) =    13.34 ; given_sigma_ion_O2_to_O_cm2( 424) =    4.69
  given_sigma_ion_O2_to_O2_cm2( 425) =    13.55 ; given_sigma_ion_O2_to_O_cm2( 425) =    4.63
  given_sigma_ion_O2_to_O2_cm2( 426) =    13.57 ; given_sigma_ion_O2_to_O_cm2( 426) =    4.62
  given_sigma_ion_O2_to_O2_cm2( 427) =    13.76 ; given_sigma_ion_O2_to_O_cm2( 427) =    4.56
  given_sigma_ion_O2_to_O2_cm2( 428) =    13.87 ; given_sigma_ion_O2_to_O_cm2( 428) =    4.53
  given_sigma_ion_O2_to_O2_cm2( 429) =    13.97 ; given_sigma_ion_O2_to_O_cm2( 429) =    4.62
  given_sigma_ion_O2_to_O2_cm2( 430) =    14.07 ; given_sigma_ion_O2_to_O_cm2( 430) =    4.73
  given_sigma_ion_O2_to_O2_cm2( 431) =    14.40 ; given_sigma_ion_O2_to_O_cm2( 431) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 432) =    14.80 ; given_sigma_ion_O2_to_O_cm2( 432) =    4.79
  given_sigma_ion_O2_to_O2_cm2( 433) =    14.81 ; given_sigma_ion_O2_to_O_cm2( 433) =    4.79
  given_sigma_ion_O2_to_O2_cm2( 434) =    14.85 ; given_sigma_ion_O2_to_O_cm2( 434) =    4.79
  given_sigma_ion_O2_to_O2_cm2( 435) =    14.87 ; given_sigma_ion_O2_to_O_cm2( 435) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 436) =    14.88 ; given_sigma_ion_O2_to_O_cm2( 436) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 437) =    14.93 ; given_sigma_ion_O2_to_O_cm2( 437) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 438) =    15.00 ; given_sigma_ion_O2_to_O_cm2( 438) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 439) =    15.04 ; given_sigma_ion_O2_to_O_cm2( 439) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 440) =    15.08 ; given_sigma_ion_O2_to_O_cm2( 440) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 441) =    15.12 ; given_sigma_ion_O2_to_O_cm2( 441) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 442) =    15.15 ; given_sigma_ion_O2_to_O_cm2( 442) =    4.81
  given_sigma_ion_O2_to_O2_cm2( 443) =    15.19 ; given_sigma_ion_O2_to_O_cm2( 443) =    4.81
  given_sigma_ion_O2_to_O2_cm2( 444) =    15.23 ; given_sigma_ion_O2_to_O_cm2( 444) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 445) =    15.26 ; given_sigma_ion_O2_to_O_cm2( 445) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 446) =    15.29 ; given_sigma_ion_O2_to_O_cm2( 446) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 447) =    15.32 ; given_sigma_ion_O2_to_O_cm2( 447) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 448) =    15.36 ; given_sigma_ion_O2_to_O_cm2( 448) =    4.79
  given_sigma_ion_O2_to_O2_cm2( 449) =    15.39 ; given_sigma_ion_O2_to_O_cm2( 449) =    4.79
  given_sigma_ion_O2_to_O2_cm2( 450) =    15.42 ; given_sigma_ion_O2_to_O_cm2( 450) =    4.79
  given_sigma_ion_O2_to_O2_cm2( 451) =    15.43 ; given_sigma_ion_O2_to_O_cm2( 451) =    4.79
  given_sigma_ion_O2_to_O2_cm2( 452) =    15.45 ; given_sigma_ion_O2_to_O_cm2( 452) =    4.78
  given_sigma_ion_O2_to_O2_cm2( 453) =    15.46 ; given_sigma_ion_O2_to_O_cm2( 453) =    4.78
  given_sigma_ion_O2_to_O2_cm2( 454) =    15.49 ; given_sigma_ion_O2_to_O_cm2( 454) =    4.78
  given_sigma_ion_O2_to_O2_cm2( 455) =    15.52 ; given_sigma_ion_O2_to_O_cm2( 455) =    4.78
  given_sigma_ion_O2_to_O2_cm2( 456) =    15.55 ; given_sigma_ion_O2_to_O_cm2( 456) =    4.79
  given_sigma_ion_O2_to_O2_cm2( 457) =    15.58 ; given_sigma_ion_O2_to_O_cm2( 457) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 458) =    15.61 ; given_sigma_ion_O2_to_O_cm2( 458) =    4.81
  given_sigma_ion_O2_to_O2_cm2( 459) =    15.64 ; given_sigma_ion_O2_to_O_cm2( 459) =    4.82
  given_sigma_ion_O2_to_O2_cm2( 460) =    15.67 ; given_sigma_ion_O2_to_O_cm2( 460) =    4.83
  given_sigma_ion_O2_to_O2_cm2( 461) =    15.70 ; given_sigma_ion_O2_to_O_cm2( 461) =    4.84
  given_sigma_ion_O2_to_O2_cm2( 462) =    15.73 ; given_sigma_ion_O2_to_O_cm2( 462) =    4.85
  given_sigma_ion_O2_to_O2_cm2( 463) =    15.76 ; given_sigma_ion_O2_to_O_cm2( 463) =    4.86
  given_sigma_ion_O2_to_O2_cm2( 464) =    15.79 ; given_sigma_ion_O2_to_O_cm2( 464) =    4.87
  given_sigma_ion_O2_to_O2_cm2( 465) =    15.82 ; given_sigma_ion_O2_to_O_cm2( 465) =    4.88
  given_sigma_ion_O2_to_O2_cm2( 466) =    15.84 ; given_sigma_ion_O2_to_O_cm2( 466) =    4.87
  given_sigma_ion_O2_to_O2_cm2( 467) =    15.86 ; given_sigma_ion_O2_to_O_cm2( 467) =    4.87
  given_sigma_ion_O2_to_O2_cm2( 468) =    15.90 ; given_sigma_ion_O2_to_O_cm2( 468) =    4.86
  given_sigma_ion_O2_to_O2_cm2( 469) =    15.94 ; given_sigma_ion_O2_to_O_cm2( 469) =    4.85
  given_sigma_ion_O2_to_O2_cm2( 470) =    15.98 ; given_sigma_ion_O2_to_O_cm2( 470) =    4.84
  given_sigma_ion_O2_to_O2_cm2( 471) =    16.01 ; given_sigma_ion_O2_to_O_cm2( 471) =    4.84
  given_sigma_ion_O2_to_O2_cm2( 472) =    16.05 ; given_sigma_ion_O2_to_O_cm2( 472) =    4.83
  given_sigma_ion_O2_to_O2_cm2( 473) =    16.06 ; given_sigma_ion_O2_to_O_cm2( 473) =    4.83
  given_sigma_ion_O2_to_O2_cm2( 474) =    16.08 ; given_sigma_ion_O2_to_O_cm2( 474) =    4.82
  given_sigma_ion_O2_to_O2_cm2( 475) =    16.09 ; given_sigma_ion_O2_to_O_cm2( 475) =    4.82
  given_sigma_ion_O2_to_O2_cm2( 476) =    16.13 ; given_sigma_ion_O2_to_O_cm2( 476) =    4.81
  given_sigma_ion_O2_to_O2_cm2( 477) =    16.17 ; given_sigma_ion_O2_to_O_cm2( 477) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 478) =    16.21 ; given_sigma_ion_O2_to_O_cm2( 478) =    4.79
  given_sigma_ion_O2_to_O2_cm2( 479) =    16.23 ; given_sigma_ion_O2_to_O_cm2( 479) =    4.79
  given_sigma_ion_O2_to_O2_cm2( 480) =    16.25 ; given_sigma_ion_O2_to_O_cm2( 480) =    4.79
  given_sigma_ion_O2_to_O2_cm2( 481) =    16.28 ; given_sigma_ion_O2_to_O_cm2( 481) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 482) =    16.29 ; given_sigma_ion_O2_to_O_cm2( 482) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 483) =    16.31 ; given_sigma_ion_O2_to_O_cm2( 483) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 484) =    16.32 ; given_sigma_ion_O2_to_O_cm2( 484) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 485) =    16.34 ; given_sigma_ion_O2_to_O_cm2( 485) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 486) =    16.36 ; given_sigma_ion_O2_to_O_cm2( 486) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 487) =    16.40 ; given_sigma_ion_O2_to_O_cm2( 487) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 488) =    16.40 ; given_sigma_ion_O2_to_O_cm2( 488) =    4.80
  given_sigma_ion_O2_to_O2_cm2( 489) =    16.43 ; given_sigma_ion_O2_to_O_cm2( 489) =    4.81
  given_sigma_ion_O2_to_O2_cm2( 490) =    16.44 ; given_sigma_ion_O2_to_O_cm2( 490) =    4.81
  given_sigma_ion_O2_to_O2_cm2( 491) =    16.46 ; given_sigma_ion_O2_to_O_cm2( 491) =    4.81
  given_sigma_ion_O2_to_O2_cm2( 492) =    16.47 ; given_sigma_ion_O2_to_O_cm2( 492) =    4.81
  given_sigma_ion_O2_to_O2_cm2( 493) =    16.49 ; given_sigma_ion_O2_to_O_cm2( 493) =    4.81
  given_sigma_ion_O2_to_O2_cm2( 494) =    16.51 ; given_sigma_ion_O2_to_O_cm2( 494) =    4.81
  given_sigma_ion_O2_to_O2_cm2( 495) =    16.53 ; given_sigma_ion_O2_to_O_cm2( 495) =    4.81
  given_sigma_ion_O2_to_O2_cm2( 496) =    16.55 ; given_sigma_ion_O2_to_O_cm2( 496) =    4.81
  given_sigma_ion_O2_to_O2_cm2( 497) =    16.58 ; given_sigma_ion_O2_to_O_cm2( 497) =    4.82
  given_sigma_ion_O2_to_O2_cm2( 498) =    16.59 ; given_sigma_ion_O2_to_O_cm2( 498) =    4.82
  given_sigma_ion_O2_to_O2_cm2( 499) =    16.60 ; given_sigma_ion_O2_to_O_cm2( 499) =    4.83
  given_sigma_ion_O2_to_O2_cm2( 500) =    16.62 ; given_sigma_ion_O2_to_O_cm2( 500) =    4.84
  given_sigma_ion_O2_to_O2_cm2( 501) =    16.63 ; given_sigma_ion_O2_to_O_cm2( 501) =    4.84
  given_sigma_ion_O2_to_O2_cm2( 502) =    16.64 ; given_sigma_ion_O2_to_O_cm2( 502) =    4.85
  given_sigma_ion_O2_to_O2_cm2( 503) =    16.66 ; given_sigma_ion_O2_to_O_cm2( 503) =    4.86
  given_sigma_ion_O2_to_O2_cm2( 504) =    16.67 ; given_sigma_ion_O2_to_O_cm2( 504) =    4.86
  given_sigma_ion_O2_to_O2_cm2( 505) =    16.68 ; given_sigma_ion_O2_to_O_cm2( 505) =    4.87
  given_sigma_ion_O2_to_O2_cm2( 506) =    16.68 ; given_sigma_ion_O2_to_O_cm2( 506) =    4.87
  given_sigma_ion_O2_to_O2_cm2( 507) =    16.69 ; given_sigma_ion_O2_to_O_cm2( 507) =    4.87
  given_sigma_ion_O2_to_O2_cm2( 508) =    16.70 ; given_sigma_ion_O2_to_O_cm2( 508) =    4.88
  given_sigma_ion_O2_to_O2_cm2( 509) =    16.71 ; given_sigma_ion_O2_to_O_cm2( 509) =    4.89
  given_sigma_ion_O2_to_O2_cm2( 510) =    16.72 ; given_sigma_ion_O2_to_O_cm2( 510) =    4.89
  given_sigma_ion_O2_to_O2_cm2( 511) =    16.74 ; given_sigma_ion_O2_to_O_cm2( 511) =    4.90
  given_sigma_ion_O2_to_O2_cm2( 512) =    16.76 ; given_sigma_ion_O2_to_O_cm2( 512) =    4.91
  given_sigma_ion_O2_to_O2_cm2( 513) =    16.76 ; given_sigma_ion_O2_to_O_cm2( 513) =    4.91
  given_sigma_ion_O2_to_O2_cm2( 514) =    16.78 ; given_sigma_ion_O2_to_O_cm2( 514) =    4.92
  given_sigma_ion_O2_to_O2_cm2( 515) =    16.81 ; given_sigma_ion_O2_to_O_cm2( 515) =    4.93
  given_sigma_ion_O2_to_O2_cm2( 516) =    16.85 ; given_sigma_ion_O2_to_O_cm2( 516) =    4.93
  given_sigma_ion_O2_to_O2_cm2( 517) =    16.89 ; given_sigma_ion_O2_to_O_cm2( 517) =    4.93
  given_sigma_ion_O2_to_O2_cm2( 518) =    16.92 ; given_sigma_ion_O2_to_O_cm2( 518) =    4.94
  given_sigma_ion_O2_to_O2_cm2( 519) =    16.96 ; given_sigma_ion_O2_to_O_cm2( 519) =    4.94
  given_sigma_ion_O2_to_O2_cm2( 520) =    16.97 ; given_sigma_ion_O2_to_O_cm2( 520) =    4.94
  given_sigma_ion_O2_to_O2_cm2( 521) =    16.99 ; given_sigma_ion_O2_to_O_cm2( 521) =    4.95
  given_sigma_ion_O2_to_O2_cm2( 522) =    17.03 ; given_sigma_ion_O2_to_O_cm2( 522) =    4.95
  given_sigma_ion_O2_to_O2_cm2( 523) =    17.07 ; given_sigma_ion_O2_to_O_cm2( 523) =    4.95
  given_sigma_ion_O2_to_O2_cm2( 524) =    17.10 ; given_sigma_ion_O2_to_O_cm2( 524) =    4.96
  given_sigma_ion_O2_to_O2_cm2( 525) =    17.13 ; given_sigma_ion_O2_to_O_cm2( 525) =    4.96
  given_sigma_ion_O2_to_O2_cm2( 526) =    17.14 ; given_sigma_ion_O2_to_O_cm2( 526) =    4.96
  given_sigma_ion_O2_to_O2_cm2( 527) =    17.16 ; given_sigma_ion_O2_to_O_cm2( 527) =    4.99
  given_sigma_ion_O2_to_O2_cm2( 528) =    17.18 ; given_sigma_ion_O2_to_O_cm2( 528) =    5.02
  given_sigma_ion_O2_to_O2_cm2( 529) =    17.20 ; given_sigma_ion_O2_to_O_cm2( 529) =    5.05
  given_sigma_ion_O2_to_O2_cm2( 530) =    17.20 ; given_sigma_ion_O2_to_O_cm2( 530) =    5.06
  given_sigma_ion_O2_to_O2_cm2( 531) =    17.21 ; given_sigma_ion_O2_to_O_cm2( 531) =    5.08
  given_sigma_ion_O2_to_O2_cm2( 532) =    17.22 ; given_sigma_ion_O2_to_O_cm2( 532) =    5.08
  given_sigma_ion_O2_to_O2_cm2( 533) =    17.24 ; given_sigma_ion_O2_to_O_cm2( 533) =    5.11
  given_sigma_ion_O2_to_O2_cm2( 534) =    17.25 ; given_sigma_ion_O2_to_O_cm2( 534) =    5.14
  given_sigma_ion_O2_to_O2_cm2( 535) =    17.25 ; given_sigma_ion_O2_to_O_cm2( 535) =    5.15
  given_sigma_ion_O2_to_O2_cm2( 536) =    17.26 ; given_sigma_ion_O2_to_O_cm2( 536) =    5.16
  given_sigma_ion_O2_to_O2_cm2( 537) =    17.27 ; given_sigma_ion_O2_to_O_cm2( 537) =    5.18
  given_sigma_ion_O2_to_O2_cm2( 538) =    17.29 ; given_sigma_ion_O2_to_O_cm2( 538) =    5.21
  given_sigma_ion_O2_to_O2_cm2( 539) =    17.30 ; given_sigma_ion_O2_to_O_cm2( 539) =    5.22
  given_sigma_ion_O2_to_O2_cm2( 540) =    17.31 ; given_sigma_ion_O2_to_O_cm2( 540) =    5.24
  given_sigma_ion_O2_to_O2_cm2( 541) =    17.33 ; given_sigma_ion_O2_to_O_cm2( 541) =    5.27
  given_sigma_ion_O2_to_O2_cm2( 542) =    17.36 ; given_sigma_ion_O2_to_O_cm2( 542) =    5.27
  given_sigma_ion_O2_to_O2_cm2( 543) =    17.37 ; given_sigma_ion_O2_to_O_cm2( 543) =    5.27
  given_sigma_ion_O2_to_O2_cm2( 544) =    17.39 ; given_sigma_ion_O2_to_O_cm2( 544) =    5.27
  given_sigma_ion_O2_to_O2_cm2( 545) =    17.41 ; given_sigma_ion_O2_to_O_cm2( 545) =    5.27
  given_sigma_ion_O2_to_O2_cm2( 546) =    17.41 ; given_sigma_ion_O2_to_O_cm2( 546) =    5.27
  given_sigma_ion_O2_to_O2_cm2( 547) =    17.43 ; given_sigma_ion_O2_to_O_cm2( 547) =    5.28
  given_sigma_ion_O2_to_O2_cm2( 548) =    17.44 ; given_sigma_ion_O2_to_O_cm2( 548) =    5.28
  given_sigma_ion_O2_to_O2_cm2( 549) =    17.45 ; given_sigma_ion_O2_to_O_cm2( 549) =    5.28
  given_sigma_ion_O2_to_O2_cm2( 550) =    17.48 ; given_sigma_ion_O2_to_O_cm2( 550) =    5.28
  given_sigma_ion_O2_to_O2_cm2( 551) =    17.51 ; given_sigma_ion_O2_to_O_cm2( 551) =    5.28
  given_sigma_ion_O2_to_O2_cm2( 552) =    17.52 ; given_sigma_ion_O2_to_O_cm2( 552) =    5.28
  given_sigma_ion_O2_to_O2_cm2( 553) =    17.56 ; given_sigma_ion_O2_to_O_cm2( 553) =    5.28
  given_sigma_ion_O2_to_O2_cm2( 554) =    17.58 ; given_sigma_ion_O2_to_O_cm2( 554) =    5.28
  given_sigma_ion_O2_to_O2_cm2( 555) =    17.59 ; given_sigma_ion_O2_to_O_cm2( 555) =    5.29
  given_sigma_ion_O2_to_O2_cm2( 556) =    17.63 ; given_sigma_ion_O2_to_O_cm2( 556) =    5.29
  given_sigma_ion_O2_to_O2_cm2( 557) =    17.67 ; given_sigma_ion_O2_to_O_cm2( 557) =    5.29
  given_sigma_ion_O2_to_O2_cm2( 558) =    17.69 ; given_sigma_ion_O2_to_O_cm2( 558) =    5.29
  given_sigma_ion_O2_to_O2_cm2( 559) =    17.71 ; given_sigma_ion_O2_to_O_cm2( 559) =    5.29
  given_sigma_ion_O2_to_O2_cm2( 560) =    17.75 ; given_sigma_ion_O2_to_O_cm2( 560) =    5.31
  given_sigma_ion_O2_to_O2_cm2( 561) =    17.77 ; given_sigma_ion_O2_to_O_cm2( 561) =    5.33
  given_sigma_ion_O2_to_O2_cm2( 562) =    17.78 ; given_sigma_ion_O2_to_O_cm2( 562) =    5.34
  given_sigma_ion_O2_to_O2_cm2( 563) =    17.82 ; given_sigma_ion_O2_to_O_cm2( 563) =    5.36
  given_sigma_ion_O2_to_O2_cm2( 564) =    17.86 ; given_sigma_ion_O2_to_O_cm2( 564) =    5.38
  given_sigma_ion_O2_to_O2_cm2( 565) =    17.90 ; given_sigma_ion_O2_to_O_cm2( 565) =    5.40
  given_sigma_ion_O2_to_O2_cm2( 566) =    17.94 ; given_sigma_ion_O2_to_O_cm2( 566) =    5.42
  given_sigma_ion_O2_to_O2_cm2( 567) =    17.98 ; given_sigma_ion_O2_to_O_cm2( 567) =    5.44
  given_sigma_ion_O2_to_O2_cm2( 568) =    18.02 ; given_sigma_ion_O2_to_O_cm2( 568) =    5.46
  given_sigma_ion_O2_to_O2_cm2( 569) =    18.05 ; given_sigma_ion_O2_to_O_cm2( 569) =    5.49
  given_sigma_ion_O2_to_O2_cm2( 570) =    18.06 ; given_sigma_ion_O2_to_O_cm2( 570) =    5.49
  given_sigma_ion_O2_to_O2_cm2( 571) =    18.07 ; given_sigma_ion_O2_to_O_cm2( 571) =    5.49
  given_sigma_ion_O2_to_O2_cm2( 572) =    18.09 ; given_sigma_ion_O2_to_O_cm2( 572) =    5.51
  given_sigma_ion_O2_to_O2_cm2( 573) =    18.13 ; given_sigma_ion_O2_to_O_cm2( 573) =    5.51
  given_sigma_ion_O2_to_O2_cm2( 574) =    18.13 ; given_sigma_ion_O2_to_O_cm2( 574) =    5.51
  given_sigma_ion_O2_to_O2_cm2( 575) =    18.17 ; given_sigma_ion_O2_to_O_cm2( 575) =    5.51
  given_sigma_ion_O2_to_O2_cm2( 576) =    18.21 ; given_sigma_ion_O2_to_O_cm2( 576) =    5.51
  given_sigma_ion_O2_to_O2_cm2( 577) =    18.25 ; given_sigma_ion_O2_to_O_cm2( 577) =    5.51
  given_sigma_ion_O2_to_O2_cm2( 578) =    18.28 ; given_sigma_ion_O2_to_O_cm2( 578) =    5.52
  given_sigma_ion_O2_to_O2_cm2( 579) =    18.40 ; given_sigma_ion_O2_to_O_cm2( 579) =    5.52
  given_sigma_ion_O2_to_O2_cm2( 580) =    18.48 ; given_sigma_ion_O2_to_O_cm2( 580) =    5.52
  given_sigma_ion_O2_to_O2_cm2( 581) =    18.81 ; given_sigma_ion_O2_to_O_cm2( 581) =    5.47
  given_sigma_ion_O2_to_O2_cm2( 582) =    19.08 ; given_sigma_ion_O2_to_O_cm2( 582) =    5.42
  given_sigma_ion_O2_to_O2_cm2( 583) =    19.14 ; given_sigma_ion_O2_to_O_cm2( 583) =    5.39
  given_sigma_ion_O2_to_O2_cm2( 584) =    19.17 ; given_sigma_ion_O2_to_O_cm2( 584) =    5.37
  given_sigma_ion_O2_to_O2_cm2( 585) =    19.60 ; given_sigma_ion_O2_to_O_cm2( 585) =    5.14
  given_sigma_ion_O2_to_O2_cm2( 586) =    19.98 ; given_sigma_ion_O2_to_O_cm2( 586) =    4.92
  given_sigma_ion_O2_to_O2_cm2( 587) =    20.50 ; given_sigma_ion_O2_to_O_cm2( 587) =    4.75
  given_sigma_ion_O2_to_O2_cm2( 588) =    20.72 ; given_sigma_ion_O2_to_O_cm2( 588) =    4.68
  given_sigma_ion_O2_to_O2_cm2( 589) =    20.92 ; given_sigma_ion_O2_to_O_cm2( 589) =    4.36
  given_sigma_ion_O2_to_O2_cm2( 590) =    21.06 ; given_sigma_ion_O2_to_O_cm2( 590) =    4.14
  given_sigma_ion_O2_to_O2_cm2( 591) =    21.34 ; given_sigma_ion_O2_to_O_cm2( 591) =    4.56
  given_sigma_ion_O2_to_O2_cm2( 592) =    21.48 ; given_sigma_ion_O2_to_O_cm2( 592) =    4.42
  given_sigma_ion_O2_to_O2_cm2( 593) =    21.67 ; given_sigma_ion_O2_to_O_cm2( 593) =    4.23
  given_sigma_ion_O2_to_O2_cm2( 594) =    21.60 ; given_sigma_ion_O2_to_O_cm2( 594) =    4.44
  given_sigma_ion_O2_to_O2_cm2( 595) =    21.60 ; given_sigma_ion_O2_to_O_cm2( 595) =    4.45
  given_sigma_ion_O2_to_O2_cm2( 596) =    21.58 ; given_sigma_ion_O2_to_O_cm2( 596) =    4.52
  given_sigma_ion_O2_to_O2_cm2( 597) =    21.58 ; given_sigma_ion_O2_to_O_cm2( 597) =    4.52
  given_sigma_ion_O2_to_O2_cm2( 598) =    21.59 ; given_sigma_ion_O2_to_O_cm2( 598) =    4.51
  given_sigma_ion_O2_to_O2_cm2( 599) =    21.55 ; given_sigma_ion_O2_to_O_cm2( 599) =    4.50
  given_sigma_ion_O2_to_O2_cm2( 600) =    21.52 ; given_sigma_ion_O2_to_O_cm2( 600) =    4.48
  given_sigma_ion_O2_to_O2_cm2( 601) =    21.41 ; given_sigma_ion_O2_to_O_cm2( 601) =    4.45
  given_sigma_ion_O2_to_O2_cm2( 602) =    21.36 ; given_sigma_ion_O2_to_O_cm2( 602) =    4.44
  given_sigma_ion_O2_to_O2_cm2( 603) =    21.25 ; given_sigma_ion_O2_to_O_cm2( 603) =    4.41
  given_sigma_ion_O2_to_O2_cm2( 604) =    21.13 ; given_sigma_ion_O2_to_O_cm2( 604) =    4.37
  given_sigma_ion_O2_to_O2_cm2( 605) =    20.56 ; given_sigma_ion_O2_to_O_cm2( 605) =    4.24
  given_sigma_ion_O2_to_O2_cm2( 606) =    20.40 ; given_sigma_ion_O2_to_O_cm2( 606) =    4.21
  given_sigma_ion_O2_to_O2_cm2( 607) =    18.98 ; given_sigma_ion_O2_to_O_cm2( 607) =    3.90
  given_sigma_ion_O2_to_O2_cm2( 608) =    18.84 ; given_sigma_ion_O2_to_O_cm2( 608) =    3.88
  given_sigma_ion_O2_to_O2_cm2( 609) =    18.58 ; given_sigma_ion_O2_to_O_cm2( 609) =    3.82
  given_sigma_ion_O2_to_O2_cm2( 610) =    18.08 ; given_sigma_ion_O2_to_O_cm2( 610) =    3.72
  given_sigma_ion_O2_to_O2_cm2( 611) =    16.12 ; given_sigma_ion_O2_to_O_cm2( 611) =    3.18
  given_sigma_ion_O2_to_O2_cm2( 612) =    16.12 ; given_sigma_ion_O2_to_O_cm2( 612) =    3.28
  given_sigma_ion_O2_to_O2_cm2( 613) =    17.09 ; given_sigma_ion_O2_to_O_cm2( 613) =    3.72
  given_sigma_ion_O2_to_O2_cm2( 614) =    17.76 ; given_sigma_ion_O2_to_O_cm2( 614) =    4.04
  given_sigma_ion_O2_to_O2_cm2( 615) =    19.98 ; given_sigma_ion_O2_to_O_cm2( 615) =    3.82
  given_sigma_ion_O2_to_O2_cm2( 616) =    24.56 ; given_sigma_ion_O2_to_O_cm2( 616) =    3.14
  given_sigma_ion_O2_to_O2_cm2( 617) =    25.54 ; given_sigma_ion_O2_to_O_cm2( 617) =    2.95
  given_sigma_ion_O2_to_O2_cm2( 618) =    25.67 ; given_sigma_ion_O2_to_O_cm2( 618) =    2.93
  given_sigma_ion_O2_to_O2_cm2( 619) =    26.62 ; given_sigma_ion_O2_to_O_cm2( 619) =    2.57
  given_sigma_ion_O2_to_O2_cm2( 620) =    23.63 ; given_sigma_ion_O2_to_O_cm2( 620) =    1.94
  given_sigma_ion_O2_to_O2_cm2( 621) =    23.32 ; given_sigma_ion_O2_to_O_cm2( 621) =    1.88
  given_sigma_ion_O2_to_O2_cm2( 622) =    23.01 ; given_sigma_ion_O2_to_O_cm2( 622) =    1.81
  given_sigma_ion_O2_to_O2_cm2( 623) =    24.43 ; given_sigma_ion_O2_to_O_cm2( 623) =    1.85
  given_sigma_ion_O2_to_O2_cm2( 624) =    24.97 ; given_sigma_ion_O2_to_O_cm2( 624) =    1.86
  given_sigma_ion_O2_to_O2_cm2( 625) =    26.79 ; given_sigma_ion_O2_to_O_cm2( 625) =    1.92
  given_sigma_ion_O2_to_O2_cm2( 626) =    26.65 ; given_sigma_ion_O2_to_O_cm2( 626) =    1.89
  given_sigma_ion_O2_to_O2_cm2( 627) =    23.81 ; given_sigma_ion_O2_to_O_cm2( 627) =    1.59
  given_sigma_ion_O2_to_O2_cm2( 628) =    22.72 ; given_sigma_ion_O2_to_O_cm2( 628) =    1.49
  given_sigma_ion_O2_to_O2_cm2( 629) =    21.66 ; given_sigma_ion_O2_to_O_cm2( 629) =    1.39
  given_sigma_ion_O2_to_O2_cm2( 630) =    25.29 ; given_sigma_ion_O2_to_O_cm2( 630) =    1.60
  given_sigma_ion_O2_to_O2_cm2( 631) =    26.55 ; given_sigma_ion_O2_to_O_cm2( 631) =    1.58
  given_sigma_ion_O2_to_O2_cm2( 632) =    26.59 ; given_sigma_ion_O2_to_O_cm2( 632) =    1.48
  given_sigma_ion_O2_to_O2_cm2( 633) =    26.60 ; given_sigma_ion_O2_to_O_cm2( 633) =    1.46
  given_sigma_ion_O2_to_O2_cm2( 634) =    23.56 ; given_sigma_ion_O2_to_O_cm2( 634) =    1.23
  given_sigma_ion_O2_to_O2_cm2( 635) =    22.43 ; given_sigma_ion_O2_to_O_cm2( 635) =    1.14
  given_sigma_ion_O2_to_O2_cm2( 636) =    24.47 ; given_sigma_ion_O2_to_O_cm2( 636) =    1.22
  given_sigma_ion_O2_to_O2_cm2( 637) =    26.50 ; given_sigma_ion_O2_to_O_cm2( 637) =    1.29
  given_sigma_ion_O2_to_O2_cm2( 638) =    31.98 ; given_sigma_ion_O2_to_O_cm2( 638) =    1.47
  given_sigma_ion_O2_to_O2_cm2( 639) =    28.47 ; given_sigma_ion_O2_to_O_cm2( 639) =    1.28
  given_sigma_ion_O2_to_O2_cm2( 640) =    26.13 ; given_sigma_ion_O2_to_O_cm2( 640) =    1.15
  given_sigma_ion_O2_to_O2_cm2( 641) =    27.49 ; given_sigma_ion_O2_to_O_cm2( 641) =    1.16
  given_sigma_ion_O2_to_O2_cm2( 642) =    26.78 ; given_sigma_ion_O2_to_O_cm2( 642) =    1.12
  given_sigma_ion_O2_to_O2_cm2( 643) =    23.59 ; given_sigma_ion_O2_to_O_cm2( 643) =    0.97
  given_sigma_ion_O2_to_O2_cm2( 644) =    22.89 ; given_sigma_ion_O2_to_O_cm2( 644) =    0.94
  given_sigma_ion_O2_to_O2_cm2( 645) =    21.13 ; given_sigma_ion_O2_to_O_cm2( 645) =    0.86
  given_sigma_ion_O2_to_O2_cm2( 646) =    24.20 ; given_sigma_ion_O2_to_O_cm2( 646) =    0.98
  given_sigma_ion_O2_to_O2_cm2( 647) =    29.70 ; given_sigma_ion_O2_to_O_cm2( 647) =    1.18
  given_sigma_ion_O2_to_O2_cm2( 648) =    26.71 ; given_sigma_ion_O2_to_O_cm2( 648) =    1.04
  given_sigma_ion_O2_to_O2_cm2( 649) =    24.02 ; given_sigma_ion_O2_to_O_cm2( 649) =    0.92
  given_sigma_ion_O2_to_O2_cm2( 650) =    22.68 ; given_sigma_ion_O2_to_O_cm2( 650) =    0.86
  given_sigma_ion_O2_to_O2_cm2( 651) =    23.85 ; given_sigma_ion_O2_to_O_cm2( 651) =    0.90
  given_sigma_ion_O2_to_O2_cm2( 652) =    25.71 ; given_sigma_ion_O2_to_O_cm2( 652) =    0.95
  given_sigma_ion_O2_to_O2_cm2( 653) =    26.91 ; given_sigma_ion_O2_to_O_cm2( 653) =    0.98
  given_sigma_ion_O2_to_O2_cm2( 654) =    23.90 ; given_sigma_ion_O2_to_O_cm2( 654) =    0.87
  given_sigma_ion_O2_to_O2_cm2( 655) =    25.73 ; given_sigma_ion_O2_to_O_cm2( 655) =    0.91
  given_sigma_ion_O2_to_O2_cm2( 656) =    28.63 ; given_sigma_ion_O2_to_O_cm2( 656) =    1.00
  given_sigma_ion_O2_to_O2_cm2( 657) =    30.37 ; given_sigma_ion_O2_to_O_cm2( 657) =    1.05
  given_sigma_ion_O2_to_O2_cm2( 658) =    30.14 ; given_sigma_ion_O2_to_O_cm2( 658) =    1.04
  given_sigma_ion_O2_to_O2_cm2( 659) =    29.43 ; given_sigma_ion_O2_to_O_cm2( 659) =    1.00
  given_sigma_ion_O2_to_O2_cm2( 660) =    28.73 ; given_sigma_ion_O2_to_O_cm2( 660) =    0.97
  given_sigma_ion_O2_to_O2_cm2( 661) =    27.04 ; given_sigma_ion_O2_to_O_cm2( 661) =    0.94
  given_sigma_ion_O2_to_O2_cm2( 662) =    25.58 ; given_sigma_ion_O2_to_O_cm2( 662) =    0.91
  given_sigma_ion_O2_to_O2_cm2( 663) =    24.61 ; given_sigma_ion_O2_to_O_cm2( 663) =    0.89
  given_sigma_ion_O2_to_O2_cm2( 664) =    22.17 ; given_sigma_ion_O2_to_O_cm2( 664) =    0.84
  given_sigma_ion_O2_to_O2_cm2( 665) =    27.65 ; given_sigma_ion_O2_to_O_cm2( 665) =    1.08
  given_sigma_ion_O2_to_O2_cm2( 666) =    29.86 ; given_sigma_ion_O2_to_O_cm2( 666) =    1.18
  given_sigma_ion_O2_to_O2_cm2( 667) =    28.64 ; given_sigma_ion_O2_to_O_cm2( 667) =    1.15
  given_sigma_ion_O2_to_O2_cm2( 668) =    27.43 ; given_sigma_ion_O2_to_O_cm2( 668) =    1.09
  given_sigma_ion_O2_to_O2_cm2( 669) =    24.20 ; given_sigma_ion_O2_to_O_cm2( 669) =    0.95
  given_sigma_ion_O2_to_O2_cm2( 670) =    25.22 ; given_sigma_ion_O2_to_O_cm2( 670) =    0.98
  given_sigma_ion_O2_to_O2_cm2( 671) =    26.75 ; given_sigma_ion_O2_to_O_cm2( 671) =    1.04
  given_sigma_ion_O2_to_O2_cm2( 672) =    25.09 ; given_sigma_ion_O2_to_O_cm2( 672) =    0.96
  given_sigma_ion_O2_to_O2_cm2( 673) =    24.38 ; given_sigma_ion_O2_to_O_cm2( 673) =    0.92
  given_sigma_ion_O2_to_O2_cm2( 674) =    26.43 ; given_sigma_ion_O2_to_O_cm2( 674) =    1.04
  given_sigma_ion_O2_to_O2_cm2( 675) =    24.98 ; given_sigma_ion_O2_to_O_cm2( 675) =    1.00
  given_sigma_ion_O2_to_O2_cm2( 676) =    24.12 ; given_sigma_ion_O2_to_O_cm2( 676) =    0.98
  given_sigma_ion_O2_to_O2_cm2( 677) =    24.34 ; given_sigma_ion_O2_to_O_cm2( 677) =    1.00
  given_sigma_ion_O2_to_O2_cm2( 678) =    25.08 ; given_sigma_ion_O2_to_O_cm2( 678) =    1.08
  given_sigma_ion_O2_to_O2_cm2( 679) =    25.37 ; given_sigma_ion_O2_to_O_cm2( 679) =    1.11
  given_sigma_ion_O2_to_O2_cm2( 680) =    25.75 ; given_sigma_ion_O2_to_O_cm2( 680) =    1.15
  given_sigma_ion_O2_to_O2_cm2( 681) =    26.37 ; given_sigma_ion_O2_to_O_cm2( 681) =    1.22
  given_sigma_ion_O2_to_O2_cm2( 682) =    26.50 ; given_sigma_ion_O2_to_O_cm2( 682) =    1.24
  given_sigma_ion_O2_to_O2_cm2( 683) =    26.85 ; given_sigma_ion_O2_to_O_cm2( 683) =    1.28
  given_sigma_ion_O2_to_O2_cm2( 684) =    24.65 ; given_sigma_ion_O2_to_O_cm2( 684) =    1.20
  given_sigma_ion_O2_to_O2_cm2( 685) =    20.22 ; given_sigma_ion_O2_to_O_cm2( 685) =    1.02
  given_sigma_ion_O2_to_O2_cm2( 686) =    19.34 ; given_sigma_ion_O2_to_O_cm2( 686) =    0.98
  given_sigma_ion_O2_to_O2_cm2( 687) =    20.39 ; given_sigma_ion_O2_to_O_cm2( 687) =    1.04
  given_sigma_ion_O2_to_O2_cm2( 688) =    23.61 ; given_sigma_ion_O2_to_O_cm2( 688) =    1.19
  given_sigma_ion_O2_to_O2_cm2( 689) =    23.20 ; given_sigma_ion_O2_to_O_cm2( 689) =    1.13
  given_sigma_ion_O2_to_O2_cm2( 690) =    27.41 ; given_sigma_ion_O2_to_O_cm2( 690) =    1.31
  given_sigma_ion_O2_to_O2_cm2( 691) =    28.01 ; given_sigma_ion_O2_to_O_cm2( 691) =    1.33
  given_sigma_ion_O2_to_O2_cm2( 692) =    27.51 ; given_sigma_ion_O2_to_O_cm2( 692) =    1.30
  given_sigma_ion_O2_to_O2_cm2( 693) =    26.66 ; given_sigma_ion_O2_to_O_cm2( 693) =    1.24
  given_sigma_ion_O2_to_O2_cm2( 694) =    25.81 ; given_sigma_ion_O2_to_O_cm2( 694) =    1.18
  given_sigma_ion_O2_to_O2_cm2( 695) =    24.12 ; given_sigma_ion_O2_to_O_cm2( 695) =    1.06
  given_sigma_ion_O2_to_O2_cm2( 696) =    23.43 ; given_sigma_ion_O2_to_O_cm2( 696) =    1.03
  given_sigma_ion_O2_to_O2_cm2( 697) =    24.43 ; given_sigma_ion_O2_to_O_cm2( 697) =    1.08
  given_sigma_ion_O2_to_O2_cm2( 698) =    24.93 ; given_sigma_ion_O2_to_O_cm2( 698) =    1.10
  given_sigma_ion_O2_to_O2_cm2( 699) =    26.08 ; given_sigma_ion_O2_to_O_cm2( 699) =    1.15
  given_sigma_ion_O2_to_O2_cm2( 700) =    27.42 ; given_sigma_ion_O2_to_O_cm2( 700) =    1.18
  given_sigma_ion_O2_to_O2_cm2( 701) =    27.37 ; given_sigma_ion_O2_to_O_cm2( 701) =    1.17
  given_sigma_ion_O2_to_O2_cm2( 702) =    27.09 ; given_sigma_ion_O2_to_O_cm2( 702) =    1.13
  given_sigma_ion_O2_to_O2_cm2( 703) =    26.81 ; given_sigma_ion_O2_to_O_cm2( 703) =    1.08
  given_sigma_ion_O2_to_O2_cm2( 704) =    26.53 ; given_sigma_ion_O2_to_O_cm2( 704) =    1.03
  given_sigma_ion_O2_to_O2_cm2( 705) =    26.25 ; given_sigma_ion_O2_to_O_cm2( 705) =    0.98
  given_sigma_ion_O2_to_O2_cm2( 706) =    25.96 ; given_sigma_ion_O2_to_O_cm2( 706) =    0.94
  given_sigma_ion_O2_to_O2_cm2( 707) =    25.87 ; given_sigma_ion_O2_to_O_cm2( 707) =    0.92
  given_sigma_ion_O2_to_O2_cm2( 708) =    25.66 ; given_sigma_ion_O2_to_O_cm2( 708) =    0.89
  given_sigma_ion_O2_to_O2_cm2( 709) =    25.37 ; given_sigma_ion_O2_to_O_cm2( 709) =    0.85
  given_sigma_ion_O2_to_O2_cm2( 710) =    25.07 ; given_sigma_ion_O2_to_O_cm2( 710) =    0.80
  given_sigma_ion_O2_to_O2_cm2( 711) =    24.98 ; given_sigma_ion_O2_to_O_cm2( 711) =    0.79
  given_sigma_ion_O2_to_O2_cm2( 712) =    24.77 ; given_sigma_ion_O2_to_O_cm2( 712) =    0.76
  given_sigma_ion_O2_to_O2_cm2( 713) =    24.65 ; given_sigma_ion_O2_to_O_cm2( 713) =    0.75
  given_sigma_ion_O2_to_O2_cm2( 714) =    24.49 ; given_sigma_ion_O2_to_O_cm2( 714) =    0.72
  given_sigma_ion_O2_to_O2_cm2( 715) =    23.97 ; given_sigma_ion_O2_to_O_cm2( 715) =    0.67
  given_sigma_ion_O2_to_O2_cm2( 716) =    23.50 ; given_sigma_ion_O2_to_O_cm2( 716) =    0.63
  given_sigma_ion_O2_to_O2_cm2( 717) =    23.22 ; given_sigma_ion_O2_to_O_cm2( 717) =    0.60
  given_sigma_ion_O2_to_O2_cm2( 718) =    23.11 ; given_sigma_ion_O2_to_O_cm2( 718) =    0.59
  given_sigma_ion_O2_to_O2_cm2( 719) =    22.96 ; given_sigma_ion_O2_to_O_cm2( 719) =    0.58
  given_sigma_ion_O2_to_O2_cm2( 720) =    22.77 ; given_sigma_ion_O2_to_O_cm2( 720) =    0.56
  given_sigma_ion_O2_to_O2_cm2( 721) =    22.72 ; given_sigma_ion_O2_to_O_cm2( 721) =    0.55
  given_sigma_ion_O2_to_O2_cm2( 722) =    22.47 ; given_sigma_ion_O2_to_O_cm2( 722) =    0.51
  given_sigma_ion_O2_to_O2_cm2( 723) =    22.40 ; given_sigma_ion_O2_to_O_cm2( 723) =    0.50
  given_sigma_ion_O2_to_O2_cm2( 724) =    22.22 ; given_sigma_ion_O2_to_O_cm2( 724) =    0.48
  given_sigma_ion_O2_to_O2_cm2( 725) =    21.97 ; given_sigma_ion_O2_to_O_cm2( 725) =    0.44
  given_sigma_ion_O2_to_O2_cm2( 726) =    21.82 ; given_sigma_ion_O2_to_O_cm2( 726) =    0.42
  given_sigma_ion_O2_to_O2_cm2( 727) =    20.01 ; given_sigma_ion_O2_to_O_cm2( 727) =    0.38
  given_sigma_ion_O2_to_O2_cm2( 728) =    20.63 ; given_sigma_ion_O2_to_O_cm2( 728) =    0.38
  given_sigma_ion_O2_to_O2_cm2( 729) =    20.44 ; given_sigma_ion_O2_to_O_cm2( 729) =    0.36
  given_sigma_ion_O2_to_O2_cm2( 730) =    20.10 ; given_sigma_ion_O2_to_O_cm2( 730) =    0.34
  given_sigma_ion_O2_to_O2_cm2( 731) =    20.23 ; given_sigma_ion_O2_to_O_cm2( 731) =    0.33
  given_sigma_ion_O2_to_O2_cm2( 732) =    19.32 ; given_sigma_ion_O2_to_O_cm2( 732) =    0.29
  given_sigma_ion_O2_to_O2_cm2( 733) =    18.69 ; given_sigma_ion_O2_to_O_cm2( 733) =    0.27
  given_sigma_ion_O2_to_O2_cm2( 734) =    18.97 ; given_sigma_ion_O2_to_O_cm2( 734) =    0.27
  given_sigma_ion_O2_to_O2_cm2( 735) =    19.25 ; given_sigma_ion_O2_to_O_cm2( 735) =    0.26
  given_sigma_ion_O2_to_O2_cm2( 736) =    19.81 ; given_sigma_ion_O2_to_O_cm2( 736) =    0.26
  given_sigma_ion_O2_to_O2_cm2( 737) =    20.64 ; given_sigma_ion_O2_to_O_cm2( 737) =    0.26
  given_sigma_ion_O2_to_O2_cm2( 738) =    20.19 ; given_sigma_ion_O2_to_O_cm2( 738) =    0.25
  given_sigma_ion_O2_to_O2_cm2( 739) =    19.08 ; given_sigma_ion_O2_to_O_cm2( 739) =    0.22
  given_sigma_ion_O2_to_O2_cm2( 740) =    19.99 ; given_sigma_ion_O2_to_O_cm2( 740) =    0.22
  given_sigma_ion_O2_to_O2_cm2( 741) =    20.59 ; given_sigma_ion_O2_to_O_cm2( 741) =    0.22
  given_sigma_ion_O2_to_O2_cm2( 742) =    19.24 ; given_sigma_ion_O2_to_O_cm2( 742) =    0.20
  given_sigma_ion_O2_to_O2_cm2( 743) =    20.36 ; given_sigma_ion_O2_to_O_cm2( 743) =    0.20
  given_sigma_ion_O2_to_O2_cm2( 744) =    19.16 ; given_sigma_ion_O2_to_O_cm2( 744) =    0.18
  given_sigma_ion_O2_to_O2_cm2( 745) =    18.22 ; given_sigma_ion_O2_to_O_cm2( 745) =    0.16
  given_sigma_ion_O2_to_O2_cm2( 746) =    17.95 ; given_sigma_ion_O2_to_O_cm2( 746) =    0.15
  given_sigma_ion_O2_to_O2_cm2( 747) =    17.59 ; given_sigma_ion_O2_to_O_cm2( 747) =    0.13
  given_sigma_ion_O2_to_O2_cm2( 748) =    17.45 ; given_sigma_ion_O2_to_O_cm2( 748) =    0.13
  given_sigma_ion_O2_to_O2_cm2( 749) =    17.30 ; given_sigma_ion_O2_to_O_cm2( 749) =    0.12
  given_sigma_ion_O2_to_O2_cm2( 750) =    23.11 ; given_sigma_ion_O2_to_O_cm2( 750) =    0.14
  given_sigma_ion_O2_to_O2_cm2( 751) =    21.67 ; given_sigma_ion_O2_to_O_cm2( 751) =    0.13
  given_sigma_ion_O2_to_O2_cm2( 752) =    20.69 ; given_sigma_ion_O2_to_O_cm2( 752) =    0.11
  given_sigma_ion_O2_to_O2_cm2( 753) =    19.20 ; given_sigma_ion_O2_to_O_cm2( 753) =    0.10
  given_sigma_ion_O2_to_O2_cm2( 754) =    20.32 ; given_sigma_ion_O2_to_O_cm2( 754) =    0.10
  given_sigma_ion_O2_to_O2_cm2( 755) =    20.84 ; given_sigma_ion_O2_to_O_cm2( 755) =    0.09
  given_sigma_ion_O2_to_O2_cm2( 756) =    20.93 ; given_sigma_ion_O2_to_O_cm2( 756) =    0.09
  given_sigma_ion_O2_to_O2_cm2( 757) =    21.77 ; given_sigma_ion_O2_to_O_cm2( 757) =    0.09
  given_sigma_ion_O2_to_O2_cm2( 758) =    21.13 ; given_sigma_ion_O2_to_O_cm2( 758) =    0.08
  given_sigma_ion_O2_to_O2_cm2( 759) =    19.25 ; given_sigma_ion_O2_to_O_cm2( 759) =    0.07
  given_sigma_ion_O2_to_O2_cm2( 760) =    17.43 ; given_sigma_ion_O2_to_O_cm2( 760) =    0.05
  given_sigma_ion_O2_to_O2_cm2( 761) =    20.97 ; given_sigma_ion_O2_to_O_cm2( 761) =    0.05
  given_sigma_ion_O2_to_O2_cm2( 762) =    20.21 ; given_sigma_ion_O2_to_O_cm2( 762) =    0.05
  given_sigma_ion_O2_to_O2_cm2( 763) =    19.06 ; given_sigma_ion_O2_to_O_cm2( 763) =    0.04
  given_sigma_ion_O2_to_O2_cm2( 764) =    17.11 ; given_sigma_ion_O2_to_O_cm2( 764) =    0.02
  given_sigma_ion_O2_to_O2_cm2( 765) =    18.01 ; given_sigma_ion_O2_to_O_cm2( 765) =    0.02
  given_sigma_ion_O2_to_O2_cm2( 766) =    21.64 ; given_sigma_ion_O2_to_O_cm2( 766) =    0.01
  given_sigma_ion_O2_to_O2_cm2( 767) =    23.48 ; given_sigma_ion_O2_to_O_cm2( 767) =    0.01
  given_sigma_ion_O2_to_O2_cm2( 768) =    26.27 ; given_sigma_ion_O2_to_O_cm2( 768) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 769) =    27.20 ; given_sigma_ion_O2_to_O_cm2( 769) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 770) =    21.95 ; given_sigma_ion_O2_to_O_cm2( 770) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 771) =    20.11 ; given_sigma_ion_O2_to_O_cm2( 771) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 772) =    17.34 ; given_sigma_ion_O2_to_O_cm2( 772) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 773) =    18.99 ; given_sigma_ion_O2_to_O_cm2( 773) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 774) =    22.85 ; given_sigma_ion_O2_to_O_cm2( 774) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 775) =    22.22 ; given_sigma_ion_O2_to_O_cm2( 775) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 776) =    20.35 ; given_sigma_ion_O2_to_O_cm2( 776) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 777) =    18.53 ; given_sigma_ion_O2_to_O_cm2( 777) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 778) =    15.03 ; given_sigma_ion_O2_to_O_cm2( 778) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 779) =    16.98 ; given_sigma_ion_O2_to_O_cm2( 779) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 780) =    18.10 ; given_sigma_ion_O2_to_O_cm2( 780) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 781) =    15.75 ; given_sigma_ion_O2_to_O_cm2( 781) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 782) =    23.46 ; given_sigma_ion_O2_to_O_cm2( 782) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 783) =    21.08 ; given_sigma_ion_O2_to_O_cm2( 783) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 784) =    19.75 ; given_sigma_ion_O2_to_O_cm2( 784) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 785) =    18.43 ; given_sigma_ion_O2_to_O_cm2( 785) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 786) =    17.12 ; given_sigma_ion_O2_to_O_cm2( 786) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 787) =    17.83 ; given_sigma_ion_O2_to_O_cm2( 787) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 788) =    25.88 ; given_sigma_ion_O2_to_O_cm2( 788) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 789) =    28.52 ; given_sigma_ion_O2_to_O_cm2( 789) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 790) =    31.14 ; given_sigma_ion_O2_to_O_cm2( 790) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 791) =    28.49 ; given_sigma_ion_O2_to_O_cm2( 791) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 792) =    25.82 ; given_sigma_ion_O2_to_O_cm2( 792) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 793) =    18.58 ; given_sigma_ion_O2_to_O_cm2( 793) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 794) =    16.74 ; given_sigma_ion_O2_to_O_cm2( 794) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 795) =    19.49 ; given_sigma_ion_O2_to_O_cm2( 795) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 796) =    24.92 ; given_sigma_ion_O2_to_O_cm2( 796) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 797) =    27.60 ; given_sigma_ion_O2_to_O_cm2( 797) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 798) =    17.27 ; given_sigma_ion_O2_to_O_cm2( 798) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 799) =    20.98 ; given_sigma_ion_O2_to_O_cm2( 799) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 800) =    24.64 ; given_sigma_ion_O2_to_O_cm2( 800) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 801) =    26.80 ; given_sigma_ion_O2_to_O_cm2( 801) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 802) =    27.21 ; given_sigma_ion_O2_to_O_cm2( 802) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 803) =    27.61 ; given_sigma_ion_O2_to_O_cm2( 803) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 804) =    28.80 ; given_sigma_ion_O2_to_O_cm2( 804) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 805) =    29.97 ; given_sigma_ion_O2_to_O_cm2( 805) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 806) =    25.27 ; given_sigma_ion_O2_to_O_cm2( 806) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 807) =    21.42 ; given_sigma_ion_O2_to_O_cm2( 807) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 808) =    19.89 ; given_sigma_ion_O2_to_O_cm2( 808) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 809) =    20.71 ; given_sigma_ion_O2_to_O_cm2( 809) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 810) =    23.18 ; given_sigma_ion_O2_to_O_cm2( 810) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 811) =    25.64 ; given_sigma_ion_O2_to_O_cm2( 811) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 812) =    26.46 ; given_sigma_ion_O2_to_O_cm2( 812) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 813) =    29.75 ; given_sigma_ion_O2_to_O_cm2( 813) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 814) =    26.35 ; given_sigma_ion_O2_to_O_cm2( 814) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 815) =    22.95 ; given_sigma_ion_O2_to_O_cm2( 815) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 816) =    16.15 ; given_sigma_ion_O2_to_O_cm2( 816) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 817) =    17.28 ; given_sigma_ion_O2_to_O_cm2( 817) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 818) =    18.42 ; given_sigma_ion_O2_to_O_cm2( 818) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 819) =    24.08 ; given_sigma_ion_O2_to_O_cm2( 819) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 820) =    23.75 ; given_sigma_ion_O2_to_O_cm2( 820) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 821) =    22.88 ; given_sigma_ion_O2_to_O_cm2( 821) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 822) =    31.20 ; given_sigma_ion_O2_to_O_cm2( 822) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 823) =    42.04 ; given_sigma_ion_O2_to_O_cm2( 823) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 824) =    48.92 ; given_sigma_ion_O2_to_O_cm2( 824) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 825) =    52.46 ; given_sigma_ion_O2_to_O_cm2( 825) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 826) =    38.20 ; given_sigma_ion_O2_to_O_cm2( 826) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 827) =    21.17 ; given_sigma_ion_O2_to_O_cm2( 827) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 828) =    21.64 ; given_sigma_ion_O2_to_O_cm2( 828) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 829) =    22.35 ; given_sigma_ion_O2_to_O_cm2( 829) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 830) =    22.70 ; given_sigma_ion_O2_to_O_cm2( 830) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 831) =    21.85 ; given_sigma_ion_O2_to_O_cm2( 831) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 832) =    22.48 ; given_sigma_ion_O2_to_O_cm2( 832) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 833) =    21.65 ; given_sigma_ion_O2_to_O_cm2( 833) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 834) =    24.68 ; given_sigma_ion_O2_to_O_cm2( 834) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 835) =    32.64 ; given_sigma_ion_O2_to_O_cm2( 835) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 836) =    34.61 ; given_sigma_ion_O2_to_O_cm2( 836) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 837) =    25.60 ; given_sigma_ion_O2_to_O_cm2( 837) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 838) =    28.85 ; given_sigma_ion_O2_to_O_cm2( 838) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 839) =    29.96 ; given_sigma_ion_O2_to_O_cm2( 839) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 840) =    31.08 ; given_sigma_ion_O2_to_O_cm2( 840) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 841) =    29.12 ; given_sigma_ion_O2_to_O_cm2( 841) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 842) =    27.80 ; given_sigma_ion_O2_to_O_cm2( 842) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 843) =    29.40 ; given_sigma_ion_O2_to_O_cm2( 843) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 844) =    28.00 ; given_sigma_ion_O2_to_O_cm2( 844) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 845) =    25.27 ; given_sigma_ion_O2_to_O_cm2( 845) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 846) =    24.82 ; given_sigma_ion_O2_to_O_cm2( 846) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 847) =    23.49 ; given_sigma_ion_O2_to_O_cm2( 847) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 848) =    25.01 ; given_sigma_ion_O2_to_O_cm2( 848) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 849) =    21.43 ; given_sigma_ion_O2_to_O_cm2( 849) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 850) =    20.08 ; given_sigma_ion_O2_to_O_cm2( 850) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 851) =    25.93 ; given_sigma_ion_O2_to_O_cm2( 851) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 852) =    24.75 ; given_sigma_ion_O2_to_O_cm2( 852) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 853) =    23.02 ; given_sigma_ion_O2_to_O_cm2( 853) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 854) =    21.27 ; given_sigma_ion_O2_to_O_cm2( 854) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 855) =    20.59 ; given_sigma_ion_O2_to_O_cm2( 855) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 856) =    19.92 ; given_sigma_ion_O2_to_O_cm2( 856) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 857) =    24.67 ; given_sigma_ion_O2_to_O_cm2( 857) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 858) =    27.90 ; given_sigma_ion_O2_to_O_cm2( 858) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 859) =    24.78 ; given_sigma_ion_O2_to_O_cm2( 859) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 860) =    22.24 ; given_sigma_ion_O2_to_O_cm2( 860) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 861) =    22.47 ; given_sigma_ion_O2_to_O_cm2( 861) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 862) =    23.38 ; given_sigma_ion_O2_to_O_cm2( 862) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 863) =    22.65 ; given_sigma_ion_O2_to_O_cm2( 863) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 864) =    20.90 ; given_sigma_ion_O2_to_O_cm2( 864) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 865) =    19.50 ; given_sigma_ion_O2_to_O_cm2( 865) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 866) =    21.25 ; given_sigma_ion_O2_to_O_cm2( 866) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 867) =    22.41 ; given_sigma_ion_O2_to_O_cm2( 867) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 868) =    23.10 ; given_sigma_ion_O2_to_O_cm2( 868) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 869) =    21.91 ; given_sigma_ion_O2_to_O_cm2( 869) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 870) =    21.71 ; given_sigma_ion_O2_to_O_cm2( 870) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 871) =    28.97 ; given_sigma_ion_O2_to_O_cm2( 871) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 872) =    31.99 ; given_sigma_ion_O2_to_O_cm2( 872) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 873) =    36.56 ; given_sigma_ion_O2_to_O_cm2( 873) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 874) =    42.72 ; given_sigma_ion_O2_to_O_cm2( 874) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 875) =    31.70 ; given_sigma_ion_O2_to_O_cm2( 875) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 876) =    26.19 ; given_sigma_ion_O2_to_O_cm2( 876) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 877) =    26.40 ; given_sigma_ion_O2_to_O_cm2( 877) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 878) =    26.93 ; given_sigma_ion_O2_to_O_cm2( 878) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 879) =    27.25 ; given_sigma_ion_O2_to_O_cm2( 879) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 880) =    26.20 ; given_sigma_ion_O2_to_O_cm2( 880) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 881) =    25.39 ; given_sigma_ion_O2_to_O_cm2( 881) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 882) =    25.03 ; given_sigma_ion_O2_to_O_cm2( 882) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 883) =    26.10 ; given_sigma_ion_O2_to_O_cm2( 883) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 884) =    26.81 ; given_sigma_ion_O2_to_O_cm2( 884) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 885) =    28.87 ; given_sigma_ion_O2_to_O_cm2( 885) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 886) =    29.87 ; given_sigma_ion_O2_to_O_cm2( 886) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 887) =    29.37 ; given_sigma_ion_O2_to_O_cm2( 887) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 888) =    28.20 ; given_sigma_ion_O2_to_O_cm2( 888) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 889) =    33.41 ; given_sigma_ion_O2_to_O_cm2( 889) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 890) =    38.68 ; given_sigma_ion_O2_to_O_cm2( 890) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 891) =    46.68 ; given_sigma_ion_O2_to_O_cm2( 891) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 892) =    36.35 ; given_sigma_ion_O2_to_O_cm2( 892) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 893) =    30.08 ; given_sigma_ion_O2_to_O_cm2( 893) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 894) =    30.22 ; given_sigma_ion_O2_to_O_cm2( 894) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 895) =    30.29 ; given_sigma_ion_O2_to_O_cm2( 895) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 896) =    30.34 ; given_sigma_ion_O2_to_O_cm2( 896) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 897) =    30.36 ; given_sigma_ion_O2_to_O_cm2( 897) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 898) =    31.98 ; given_sigma_ion_O2_to_O_cm2( 898) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 899) =    30.82 ; given_sigma_ion_O2_to_O_cm2( 899) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 900) =    29.86 ; given_sigma_ion_O2_to_O_cm2( 900) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 901) =    29.48 ; given_sigma_ion_O2_to_O_cm2( 901) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 902) =    30.10 ; given_sigma_ion_O2_to_O_cm2( 902) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 903) =    27.73 ; given_sigma_ion_O2_to_O_cm2( 903) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 904) =    25.86 ; given_sigma_ion_O2_to_O_cm2( 904) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 905) =    23.06 ; given_sigma_ion_O2_to_O_cm2( 905) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 906) =    22.13 ; given_sigma_ion_O2_to_O_cm2( 906) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 907) =    19.53 ; given_sigma_ion_O2_to_O_cm2( 907) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 908) =    19.98 ; given_sigma_ion_O2_to_O_cm2( 908) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 909) =    17.29 ; given_sigma_ion_O2_to_O_cm2( 909) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 910) =    16.48 ; given_sigma_ion_O2_to_O_cm2( 910) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 911) =    17.13 ; given_sigma_ion_O2_to_O_cm2( 911) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 912) =    17.28 ; given_sigma_ion_O2_to_O_cm2( 912) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 913) =    16.65 ; given_sigma_ion_O2_to_O_cm2( 913) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 914) =    14.76 ; given_sigma_ion_O2_to_O_cm2( 914) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 915) =    14.14 ; given_sigma_ion_O2_to_O_cm2( 915) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 916) =    15.58 ; given_sigma_ion_O2_to_O_cm2( 916) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 917) =    16.57 ; given_sigma_ion_O2_to_O_cm2( 917) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 918) =    16.30 ; given_sigma_ion_O2_to_O_cm2( 918) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 919) =    15.78 ; given_sigma_ion_O2_to_O_cm2( 919) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 920) =    15.30 ; given_sigma_ion_O2_to_O_cm2( 920) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 921) =    12.83 ; given_sigma_ion_O2_to_O_cm2( 921) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 922) =    11.72 ; given_sigma_ion_O2_to_O_cm2( 922) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 923) =    12.91 ; given_sigma_ion_O2_to_O_cm2( 923) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 924) =    14.15 ; given_sigma_ion_O2_to_O_cm2( 924) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 925) =    12.86 ; given_sigma_ion_O2_to_O_cm2( 925) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 926) =    11.39 ; given_sigma_ion_O2_to_O_cm2( 926) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 927) =    12.74 ; given_sigma_ion_O2_to_O_cm2( 927) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 928) =    14.05 ; given_sigma_ion_O2_to_O_cm2( 928) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 929) =    14.57 ; given_sigma_ion_O2_to_O_cm2( 929) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 930) =    16.58 ; given_sigma_ion_O2_to_O_cm2( 930) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 931) =    13.28 ; given_sigma_ion_O2_to_O_cm2( 931) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 932) =    11.90 ; given_sigma_ion_O2_to_O_cm2( 932) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 933) =    12.74 ; given_sigma_ion_O2_to_O_cm2( 933) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 934) =    11.78 ; given_sigma_ion_O2_to_O_cm2( 934) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 935) =    11.07 ; given_sigma_ion_O2_to_O_cm2( 935) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 936) =     9.68 ; given_sigma_ion_O2_to_O_cm2( 936) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 937) =    12.27 ; given_sigma_ion_O2_to_O_cm2( 937) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 938) =    13.18 ; given_sigma_ion_O2_to_O_cm2( 938) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 939) =    14.51 ; given_sigma_ion_O2_to_O_cm2( 939) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 940) =    13.60 ; given_sigma_ion_O2_to_O_cm2( 940) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 941) =    12.30 ; given_sigma_ion_O2_to_O_cm2( 941) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 942) =    11.05 ; given_sigma_ion_O2_to_O_cm2( 942) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 943) =    10.26 ; given_sigma_ion_O2_to_O_cm2( 943) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 944) =    10.48 ; given_sigma_ion_O2_to_O_cm2( 944) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 945) =    10.40 ; given_sigma_ion_O2_to_O_cm2( 945) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 946) =    10.27 ; given_sigma_ion_O2_to_O_cm2( 946) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 947) =    10.21 ; given_sigma_ion_O2_to_O_cm2( 947) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 948) =    10.13 ; given_sigma_ion_O2_to_O_cm2( 948) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 949) =    10.07 ; given_sigma_ion_O2_to_O_cm2( 949) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 950) =     9.97 ; given_sigma_ion_O2_to_O_cm2( 950) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 951) =    10.74 ; given_sigma_ion_O2_to_O_cm2( 951) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 952) =    11.00 ; given_sigma_ion_O2_to_O_cm2( 952) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 953) =    10.65 ; given_sigma_ion_O2_to_O_cm2( 953) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 954) =    10.33 ; given_sigma_ion_O2_to_O_cm2( 954) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 955) =     9.92 ; given_sigma_ion_O2_to_O_cm2( 955) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 956) =     9.73 ; given_sigma_ion_O2_to_O_cm2( 956) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 957) =     9.62 ; given_sigma_ion_O2_to_O_cm2( 957) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 958) =     9.92 ; given_sigma_ion_O2_to_O_cm2( 958) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 959) =    12.51 ; given_sigma_ion_O2_to_O_cm2( 959) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 960) =    10.86 ; given_sigma_ion_O2_to_O_cm2( 960) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 961) =     9.76 ; given_sigma_ion_O2_to_O_cm2( 961) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 962) =     9.82 ; given_sigma_ion_O2_to_O_cm2( 962) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 963) =    10.19 ; given_sigma_ion_O2_to_O_cm2( 963) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 964) =    10.61 ; given_sigma_ion_O2_to_O_cm2( 964) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 965) =    11.18 ; given_sigma_ion_O2_to_O_cm2( 965) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 966) =    11.39 ; given_sigma_ion_O2_to_O_cm2( 966) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 967) =    12.23 ; given_sigma_ion_O2_to_O_cm2( 967) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 968) =    11.94 ; given_sigma_ion_O2_to_O_cm2( 968) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 969) =    11.49 ; given_sigma_ion_O2_to_O_cm2( 969) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 970) =    10.88 ; given_sigma_ion_O2_to_O_cm2( 970) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 971) =     9.05 ; given_sigma_ion_O2_to_O_cm2( 971) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 972) =    10.55 ; given_sigma_ion_O2_to_O_cm2( 972) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 973) =    13.74 ; given_sigma_ion_O2_to_O_cm2( 973) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 974) =    15.16 ; given_sigma_ion_O2_to_O_cm2( 974) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 975) =    14.40 ; given_sigma_ion_O2_to_O_cm2( 975) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 976) =    13.65 ; given_sigma_ion_O2_to_O_cm2( 976) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 977) =    12.44 ; given_sigma_ion_O2_to_O_cm2( 977) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 978) =    11.97 ; given_sigma_ion_O2_to_O_cm2( 978) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 979) =    12.26 ; given_sigma_ion_O2_to_O_cm2( 979) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 980) =    12.55 ; given_sigma_ion_O2_to_O_cm2( 980) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 981) =    12.17 ; given_sigma_ion_O2_to_O_cm2( 981) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 982) =    11.79 ; given_sigma_ion_O2_to_O_cm2( 982) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 983) =    11.42 ; given_sigma_ion_O2_to_O_cm2( 983) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 984) =    10.68 ; given_sigma_ion_O2_to_O_cm2( 984) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 985) =    10.50 ; given_sigma_ion_O2_to_O_cm2( 985) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 986) =    12.14 ; given_sigma_ion_O2_to_O_cm2( 986) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 987) =    11.61 ; given_sigma_ion_O2_to_O_cm2( 987) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 988) =    10.56 ; given_sigma_ion_O2_to_O_cm2( 988) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 989) =    11.72 ; given_sigma_ion_O2_to_O_cm2( 989) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 990) =     8.89 ; given_sigma_ion_O2_to_O_cm2( 990) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 991) =     7.45 ; given_sigma_ion_O2_to_O_cm2( 991) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 992) =     8.43 ; given_sigma_ion_O2_to_O_cm2( 992) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 993) =    10.44 ; given_sigma_ion_O2_to_O_cm2( 993) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 994) =     9.48 ; given_sigma_ion_O2_to_O_cm2( 994) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 995) =     6.32 ; given_sigma_ion_O2_to_O_cm2( 995) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 996) =    12.24 ; given_sigma_ion_O2_to_O_cm2( 996) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 997) =    18.80 ; given_sigma_ion_O2_to_O_cm2( 997) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 998) =    16.76 ; given_sigma_ion_O2_to_O_cm2( 998) =    0.00
  given_sigma_ion_O2_to_O2_cm2( 999) =    14.90 ; given_sigma_ion_O2_to_O_cm2( 999) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1000) =    13.05 ; given_sigma_ion_O2_to_O_cm2(1000) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1001) =    11.15 ; given_sigma_ion_O2_to_O_cm2(1001) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1002) =     7.90 ; given_sigma_ion_O2_to_O_cm2(1002) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1003) =     8.92 ; given_sigma_ion_O2_to_O_cm2(1003) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1004) =     9.77 ; given_sigma_ion_O2_to_O_cm2(1004) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1005) =    11.11 ; given_sigma_ion_O2_to_O_cm2(1005) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1006) =    13.50 ; given_sigma_ion_O2_to_O_cm2(1006) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1007) =    16.11 ; given_sigma_ion_O2_to_O_cm2(1007) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1008) =    18.92 ; given_sigma_ion_O2_to_O_cm2(1008) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1009) =    21.93 ; given_sigma_ion_O2_to_O_cm2(1009) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1010) =    14.68 ; given_sigma_ion_O2_to_O_cm2(1010) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1011) =    11.19 ; given_sigma_ion_O2_to_O_cm2(1011) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1012) =     8.27 ; given_sigma_ion_O2_to_O_cm2(1012) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1013) =     8.74 ; given_sigma_ion_O2_to_O_cm2(1013) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1014) =     9.02 ; given_sigma_ion_O2_to_O_cm2(1014) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1015) =     9.11 ; given_sigma_ion_O2_to_O_cm2(1015) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1016) =    11.10 ; given_sigma_ion_O2_to_O_cm2(1016) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1017) =    10.15 ; given_sigma_ion_O2_to_O_cm2(1017) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1018) =     7.90 ; given_sigma_ion_O2_to_O_cm2(1018) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1019) =     7.04 ; given_sigma_ion_O2_to_O_cm2(1019) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1020) =     5.80 ; given_sigma_ion_O2_to_O_cm2(1020) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1021) =     6.94 ; given_sigma_ion_O2_to_O_cm2(1021) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1022) =     8.14 ; given_sigma_ion_O2_to_O_cm2(1022) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1023) =     9.15 ; given_sigma_ion_O2_to_O_cm2(1023) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1024) =     8.82 ; given_sigma_ion_O2_to_O_cm2(1024) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1025) =     9.81 ; given_sigma_ion_O2_to_O_cm2(1025) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1026) =    10.86 ; given_sigma_ion_O2_to_O_cm2(1026) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1027) =    11.58 ; given_sigma_ion_O2_to_O_cm2(1027) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1028) =    13.11 ; given_sigma_ion_O2_to_O_cm2(1028) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1029) =    16.18 ; given_sigma_ion_O2_to_O_cm2(1029) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1030) =    16.20 ; given_sigma_ion_O2_to_O_cm2(1030) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1031) =    16.18 ; given_sigma_ion_O2_to_O_cm2(1031) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1032) =    16.08 ; given_sigma_ion_O2_to_O_cm2(1032) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1033) =    11.63 ; given_sigma_ion_O2_to_O_cm2(1033) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1034) =     7.88 ; given_sigma_ion_O2_to_O_cm2(1034) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1035) =     9.31 ; given_sigma_ion_O2_to_O_cm2(1035) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1036) =    10.54 ; given_sigma_ion_O2_to_O_cm2(1036) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1037) =    10.76 ; given_sigma_ion_O2_to_O_cm2(1037) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1038) =    10.65 ; given_sigma_ion_O2_to_O_cm2(1038) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1039) =    12.98 ; given_sigma_ion_O2_to_O_cm2(1039) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1040) =    11.52 ; given_sigma_ion_O2_to_O_cm2(1040) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1041) =    10.88 ; given_sigma_ion_O2_to_O_cm2(1041) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1042) =    10.19 ; given_sigma_ion_O2_to_O_cm2(1042) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1043) =     9.01 ; given_sigma_ion_O2_to_O_cm2(1043) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1044) =     8.32 ; given_sigma_ion_O2_to_O_cm2(1044) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1045) =     7.87 ; given_sigma_ion_O2_to_O_cm2(1045) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1046) =     7.21 ; given_sigma_ion_O2_to_O_cm2(1046) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1047) =     7.24 ; given_sigma_ion_O2_to_O_cm2(1047) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1048) =     7.25 ; given_sigma_ion_O2_to_O_cm2(1048) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1049) =     6.97 ; given_sigma_ion_O2_to_O_cm2(1049) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1050) =     8.47 ; given_sigma_ion_O2_to_O_cm2(1050) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1051) =     9.15 ; given_sigma_ion_O2_to_O_cm2(1051) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1052) =    10.47 ; given_sigma_ion_O2_to_O_cm2(1052) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1053) =    12.81 ; given_sigma_ion_O2_to_O_cm2(1053) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1054) =    16.45 ; given_sigma_ion_O2_to_O_cm2(1054) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1055) =    14.28 ; given_sigma_ion_O2_to_O_cm2(1055) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1056) =    10.84 ; given_sigma_ion_O2_to_O_cm2(1056) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1057) =    10.76 ; given_sigma_ion_O2_to_O_cm2(1057) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1058) =     9.44 ; given_sigma_ion_O2_to_O_cm2(1058) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1059) =    10.91 ; given_sigma_ion_O2_to_O_cm2(1059) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1060) =    11.69 ; given_sigma_ion_O2_to_O_cm2(1060) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1061) =     9.94 ; given_sigma_ion_O2_to_O_cm2(1061) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1062) =    10.72 ; given_sigma_ion_O2_to_O_cm2(1062) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1063) =    11.14 ; given_sigma_ion_O2_to_O_cm2(1063) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1064) =    12.23 ; given_sigma_ion_O2_to_O_cm2(1064) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1065) =    13.00 ; given_sigma_ion_O2_to_O_cm2(1065) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1066) =    13.46 ; given_sigma_ion_O2_to_O_cm2(1066) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1067) =    13.61 ; given_sigma_ion_O2_to_O_cm2(1067) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1068) =    12.96 ; given_sigma_ion_O2_to_O_cm2(1068) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1069) =    11.29 ; given_sigma_ion_O2_to_O_cm2(1069) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1070) =     9.00 ; given_sigma_ion_O2_to_O_cm2(1070) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1071) =     6.50 ; given_sigma_ion_O2_to_O_cm2(1071) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1072) =     9.55 ; given_sigma_ion_O2_to_O_cm2(1072) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1073) =    11.40 ; given_sigma_ion_O2_to_O_cm2(1073) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1074) =    12.23 ; given_sigma_ion_O2_to_O_cm2(1074) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1075) =    13.14 ; given_sigma_ion_O2_to_O_cm2(1075) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1076) =    15.35 ; given_sigma_ion_O2_to_O_cm2(1076) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1077) =    17.28 ; given_sigma_ion_O2_to_O_cm2(1077) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1078) =    17.60 ; given_sigma_ion_O2_to_O_cm2(1078) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1079) =    17.17 ; given_sigma_ion_O2_to_O_cm2(1079) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1080) =    16.22 ; given_sigma_ion_O2_to_O_cm2(1080) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1081) =    15.38 ; given_sigma_ion_O2_to_O_cm2(1081) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1082) =    14.58 ; given_sigma_ion_O2_to_O_cm2(1082) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1083) =    13.91 ; given_sigma_ion_O2_to_O_cm2(1083) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1084) =    11.72 ; given_sigma_ion_O2_to_O_cm2(1084) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1085) =    10.61 ; given_sigma_ion_O2_to_O_cm2(1085) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1086) =    10.19 ; given_sigma_ion_O2_to_O_cm2(1086) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1087) =     9.12 ; given_sigma_ion_O2_to_O_cm2(1087) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1088) =     8.04 ; given_sigma_ion_O2_to_O_cm2(1088) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1089) =     7.32 ; given_sigma_ion_O2_to_O_cm2(1089) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1090) =     6.24 ; given_sigma_ion_O2_to_O_cm2(1090) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1091) =     7.42 ; given_sigma_ion_O2_to_O_cm2(1091) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1092) =    10.75 ; given_sigma_ion_O2_to_O_cm2(1092) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1093) =    14.65 ; given_sigma_ion_O2_to_O_cm2(1093) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1094) =    16.07 ; given_sigma_ion_O2_to_O_cm2(1094) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1095) =    17.73 ; given_sigma_ion_O2_to_O_cm2(1095) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1096) =    19.13 ; given_sigma_ion_O2_to_O_cm2(1096) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1097) =    21.69 ; given_sigma_ion_O2_to_O_cm2(1097) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1098) =    24.06 ; given_sigma_ion_O2_to_O_cm2(1098) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1099) =    25.38 ; given_sigma_ion_O2_to_O_cm2(1099) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1100) =    27.50 ; given_sigma_ion_O2_to_O_cm2(1100) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1101) =    20.85 ; given_sigma_ion_O2_to_O_cm2(1101) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1102) =    13.32 ; given_sigma_ion_O2_to_O_cm2(1102) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1103) =     8.23 ; given_sigma_ion_O2_to_O_cm2(1103) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1104) =     9.16 ; given_sigma_ion_O2_to_O_cm2(1104) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1105) =     9.81 ; given_sigma_ion_O2_to_O_cm2(1105) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1106) =    11.53 ; given_sigma_ion_O2_to_O_cm2(1106) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1107) =    12.25 ; given_sigma_ion_O2_to_O_cm2(1107) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1108) =    11.90 ; given_sigma_ion_O2_to_O_cm2(1108) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1109) =    10.67 ; given_sigma_ion_O2_to_O_cm2(1109) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1110) =     9.11 ; given_sigma_ion_O2_to_O_cm2(1110) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1111) =    12.77 ; given_sigma_ion_O2_to_O_cm2(1111) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1112) =    15.92 ; given_sigma_ion_O2_to_O_cm2(1112) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1113) =    18.63 ; given_sigma_ion_O2_to_O_cm2(1113) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1114) =    20.91 ; given_sigma_ion_O2_to_O_cm2(1114) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1115) =    22.47 ; given_sigma_ion_O2_to_O_cm2(1115) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1116) =    23.76 ; given_sigma_ion_O2_to_O_cm2(1116) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1117) =    20.02 ; given_sigma_ion_O2_to_O_cm2(1117) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1118) =    16.44 ; given_sigma_ion_O2_to_O_cm2(1118) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1119) =    13.51 ; given_sigma_ion_O2_to_O_cm2(1119) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1120) =    10.85 ; given_sigma_ion_O2_to_O_cm2(1120) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1121) =    11.06 ; given_sigma_ion_O2_to_O_cm2(1121) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1122) =    11.27 ; given_sigma_ion_O2_to_O_cm2(1122) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1123) =    11.63 ; given_sigma_ion_O2_to_O_cm2(1123) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1124) =    11.59 ; given_sigma_ion_O2_to_O_cm2(1124) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1125) =    11.38 ; given_sigma_ion_O2_to_O_cm2(1125) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1126) =    10.61 ; given_sigma_ion_O2_to_O_cm2(1126) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1127) =     8.71 ; given_sigma_ion_O2_to_O_cm2(1127) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1128) =     6.82 ; given_sigma_ion_O2_to_O_cm2(1128) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1129) =     5.69 ; given_sigma_ion_O2_to_O_cm2(1129) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1130) =     6.14 ; given_sigma_ion_O2_to_O_cm2(1130) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1131) =     7.26 ; given_sigma_ion_O2_to_O_cm2(1131) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1132) =     9.47 ; given_sigma_ion_O2_to_O_cm2(1132) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1133) =     9.90 ; given_sigma_ion_O2_to_O_cm2(1133) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1134) =     9.01 ; given_sigma_ion_O2_to_O_cm2(1134) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1135) =     9.78 ; given_sigma_ion_O2_to_O_cm2(1135) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1136) =    10.49 ; given_sigma_ion_O2_to_O_cm2(1136) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1137) =     7.64 ; given_sigma_ion_O2_to_O_cm2(1137) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1138) =     8.09 ; given_sigma_ion_O2_to_O_cm2(1138) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1139) =     8.50 ; given_sigma_ion_O2_to_O_cm2(1139) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1140) =     6.87 ; given_sigma_ion_O2_to_O_cm2(1140) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1141) =     5.00 ; given_sigma_ion_O2_to_O_cm2(1141) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1142) =     4.96 ; given_sigma_ion_O2_to_O_cm2(1142) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1143) =     4.86 ; given_sigma_ion_O2_to_O_cm2(1143) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1144) =     4.80 ; given_sigma_ion_O2_to_O_cm2(1144) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1145) =     4.72 ; given_sigma_ion_O2_to_O_cm2(1145) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1146) =     4.59 ; given_sigma_ion_O2_to_O_cm2(1146) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1147) =     7.35 ; given_sigma_ion_O2_to_O_cm2(1147) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1148) =     8.83 ; given_sigma_ion_O2_to_O_cm2(1148) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1149) =     8.59 ; given_sigma_ion_O2_to_O_cm2(1149) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1150) =     8.81 ; given_sigma_ion_O2_to_O_cm2(1150) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1151) =     5.95 ; given_sigma_ion_O2_to_O_cm2(1151) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1152) =     3.87 ; given_sigma_ion_O2_to_O_cm2(1152) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1153) =     8.62 ; given_sigma_ion_O2_to_O_cm2(1153) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1154) =    10.79 ; given_sigma_ion_O2_to_O_cm2(1154) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1155) =    10.84 ; given_sigma_ion_O2_to_O_cm2(1155) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1156) =    10.31 ; given_sigma_ion_O2_to_O_cm2(1156) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1157) =     7.27 ; given_sigma_ion_O2_to_O_cm2(1157) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1158) =     6.23 ; given_sigma_ion_O2_to_O_cm2(1158) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1159) =     4.56 ; given_sigma_ion_O2_to_O_cm2(1159) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1160) =     4.03 ; given_sigma_ion_O2_to_O_cm2(1160) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1161) =     4.10 ; given_sigma_ion_O2_to_O_cm2(1161) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1162) =     3.85 ; given_sigma_ion_O2_to_O_cm2(1162) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1163) =     3.77 ; given_sigma_ion_O2_to_O_cm2(1163) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1164) =     3.70 ; given_sigma_ion_O2_to_O_cm2(1164) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1165) =     5.45 ; given_sigma_ion_O2_to_O_cm2(1165) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1166) =     6.33 ; given_sigma_ion_O2_to_O_cm2(1166) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1167) =     5.31 ; given_sigma_ion_O2_to_O_cm2(1167) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1168) =     4.58 ; given_sigma_ion_O2_to_O_cm2(1168) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1169) =     4.14 ; given_sigma_ion_O2_to_O_cm2(1169) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1170) =     5.37 ; given_sigma_ion_O2_to_O_cm2(1170) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1171) =     9.07 ; given_sigma_ion_O2_to_O_cm2(1171) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1172) =    10.48 ; given_sigma_ion_O2_to_O_cm2(1172) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1173) =     8.33 ; given_sigma_ion_O2_to_O_cm2(1173) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1174) =     6.97 ; given_sigma_ion_O2_to_O_cm2(1174) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1175) =     6.44 ; given_sigma_ion_O2_to_O_cm2(1175) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1176) =     5.67 ; given_sigma_ion_O2_to_O_cm2(1176) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1177) =     4.43 ; given_sigma_ion_O2_to_O_cm2(1177) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1178) =     3.26 ; given_sigma_ion_O2_to_O_cm2(1178) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1179) =     3.33 ; given_sigma_ion_O2_to_O_cm2(1179) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1180) =     3.67 ; given_sigma_ion_O2_to_O_cm2(1180) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1181) =     4.00 ; given_sigma_ion_O2_to_O_cm2(1181) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1182) =     4.18 ; given_sigma_ion_O2_to_O_cm2(1182) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1183) =     4.07 ; given_sigma_ion_O2_to_O_cm2(1183) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1184) =     3.77 ; given_sigma_ion_O2_to_O_cm2(1184) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1185) =     4.33 ; given_sigma_ion_O2_to_O_cm2(1185) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1186) =     4.91 ; given_sigma_ion_O2_to_O_cm2(1186) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1187) =     5.21 ; given_sigma_ion_O2_to_O_cm2(1187) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1188) =     3.55 ; given_sigma_ion_O2_to_O_cm2(1188) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1189) =     2.44 ; given_sigma_ion_O2_to_O_cm2(1189) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1190) =     2.67 ; given_sigma_ion_O2_to_O_cm2(1190) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1191) =     2.96 ; given_sigma_ion_O2_to_O_cm2(1191) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1192) =     3.34 ; given_sigma_ion_O2_to_O_cm2(1192) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1193) =     3.48 ; given_sigma_ion_O2_to_O_cm2(1193) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1194) =     3.61 ; given_sigma_ion_O2_to_O_cm2(1194) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1195) =     3.80 ; given_sigma_ion_O2_to_O_cm2(1195) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1196) =     3.98 ; given_sigma_ion_O2_to_O_cm2(1196) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1197) =     3.87 ; given_sigma_ion_O2_to_O_cm2(1197) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1198) =     3.73 ; given_sigma_ion_O2_to_O_cm2(1198) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1199) =     3.65 ; given_sigma_ion_O2_to_O_cm2(1199) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1200) =     3.91 ; given_sigma_ion_O2_to_O_cm2(1200) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1201) =     4.16 ; given_sigma_ion_O2_to_O_cm2(1201) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1202) =     4.58 ; given_sigma_ion_O2_to_O_cm2(1202) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1203) =     5.26 ; given_sigma_ion_O2_to_O_cm2(1203) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1204) =     5.54 ; given_sigma_ion_O2_to_O_cm2(1204) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1205) =     5.10 ; given_sigma_ion_O2_to_O_cm2(1205) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1206) =     4.38 ; given_sigma_ion_O2_to_O_cm2(1206) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1207) =     3.67 ; given_sigma_ion_O2_to_O_cm2(1207) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1208) =     2.98 ; given_sigma_ion_O2_to_O_cm2(1208) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1209) =     3.26 ; given_sigma_ion_O2_to_O_cm2(1209) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1210) =     3.55 ; given_sigma_ion_O2_to_O_cm2(1210) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1211) =     3.66 ; given_sigma_ion_O2_to_O_cm2(1211) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1212) =     3.83 ; given_sigma_ion_O2_to_O_cm2(1212) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1213) =     4.12 ; given_sigma_ion_O2_to_O_cm2(1213) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1214) =     4.30 ; given_sigma_ion_O2_to_O_cm2(1214) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1215) =     3.94 ; given_sigma_ion_O2_to_O_cm2(1215) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1216) =     3.68 ; given_sigma_ion_O2_to_O_cm2(1216) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1217) =     3.25 ; given_sigma_ion_O2_to_O_cm2(1217) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1218) =     2.84 ; given_sigma_ion_O2_to_O_cm2(1218) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1219) =     2.68 ; given_sigma_ion_O2_to_O_cm2(1219) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1220) =     3.55 ; given_sigma_ion_O2_to_O_cm2(1220) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1221) =     4.01 ; given_sigma_ion_O2_to_O_cm2(1221) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1222) =     4.74 ; given_sigma_ion_O2_to_O_cm2(1222) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1223) =     4.88 ; given_sigma_ion_O2_to_O_cm2(1223) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1224) =     5.03 ; given_sigma_ion_O2_to_O_cm2(1224) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1225) =     5.32 ; given_sigma_ion_O2_to_O_cm2(1225) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1226) =     5.37 ; given_sigma_ion_O2_to_O_cm2(1226) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1227) =     5.60 ; given_sigma_ion_O2_to_O_cm2(1227) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1228) =     5.77 ; given_sigma_ion_O2_to_O_cm2(1228) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1229) =     5.42 ; given_sigma_ion_O2_to_O_cm2(1229) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1230) =     5.26 ; given_sigma_ion_O2_to_O_cm2(1230) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1231) =     5.09 ; given_sigma_ion_O2_to_O_cm2(1231) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1232) =     4.59 ; given_sigma_ion_O2_to_O_cm2(1232) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1233) =     4.18 ; given_sigma_ion_O2_to_O_cm2(1233) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1234) =     3.78 ; given_sigma_ion_O2_to_O_cm2(1234) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1235) =     3.38 ; given_sigma_ion_O2_to_O_cm2(1235) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1236) =     3.56 ; given_sigma_ion_O2_to_O_cm2(1236) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1237) =     3.88 ; given_sigma_ion_O2_to_O_cm2(1237) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1238) =     4.20 ; given_sigma_ion_O2_to_O_cm2(1238) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1239) =     4.52 ; given_sigma_ion_O2_to_O_cm2(1239) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1240) =     4.84 ; given_sigma_ion_O2_to_O_cm2(1240) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1241) =     5.44 ; given_sigma_ion_O2_to_O_cm2(1241) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1242) =     5.60 ; given_sigma_ion_O2_to_O_cm2(1242) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1243) =     5.90 ; given_sigma_ion_O2_to_O_cm2(1243) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1244) =     5.56 ; given_sigma_ion_O2_to_O_cm2(1244) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1245) =     5.27 ; given_sigma_ion_O2_to_O_cm2(1245) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1246) =     4.98 ; given_sigma_ion_O2_to_O_cm2(1246) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1247) =     4.67 ; given_sigma_ion_O2_to_O_cm2(1247) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1248) =     4.36 ; given_sigma_ion_O2_to_O_cm2(1248) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1249) =     4.04 ; given_sigma_ion_O2_to_O_cm2(1249) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1250) =     3.72 ; given_sigma_ion_O2_to_O_cm2(1250) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1251) =     3.58 ; given_sigma_ion_O2_to_O_cm2(1251) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1252) =     4.13 ; given_sigma_ion_O2_to_O_cm2(1252) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1253) =     5.07 ; given_sigma_ion_O2_to_O_cm2(1253) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1254) =     5.45 ; given_sigma_ion_O2_to_O_cm2(1254) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1255) =     6.03 ; given_sigma_ion_O2_to_O_cm2(1255) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1256) =     7.01 ; given_sigma_ion_O2_to_O_cm2(1256) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1257) =     8.01 ; given_sigma_ion_O2_to_O_cm2(1257) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1258) =     8.46 ; given_sigma_ion_O2_to_O_cm2(1258) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1259) =     8.75 ; given_sigma_ion_O2_to_O_cm2(1259) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1260) =     8.23 ; given_sigma_ion_O2_to_O_cm2(1260) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1261) =     7.89 ; given_sigma_ion_O2_to_O_cm2(1261) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1262) =     7.53 ; given_sigma_ion_O2_to_O_cm2(1262) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1263) =     7.15 ; given_sigma_ion_O2_to_O_cm2(1263) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1264) =     6.93 ; given_sigma_ion_O2_to_O_cm2(1264) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1265) =     6.33 ; given_sigma_ion_O2_to_O_cm2(1265) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1266) =     5.74 ; given_sigma_ion_O2_to_O_cm2(1266) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1267) =     5.16 ; given_sigma_ion_O2_to_O_cm2(1267) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1268) =     4.59 ; given_sigma_ion_O2_to_O_cm2(1268) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1269) =     4.02 ; given_sigma_ion_O2_to_O_cm2(1269) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1270) =     3.79 ; given_sigma_ion_O2_to_O_cm2(1270) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1271) =     3.45 ; given_sigma_ion_O2_to_O_cm2(1271) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1272) =     3.12 ; given_sigma_ion_O2_to_O_cm2(1272) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1273) =     3.79 ; given_sigma_ion_O2_to_O_cm2(1273) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1274) =     5.46 ; given_sigma_ion_O2_to_O_cm2(1274) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1275) =     7.10 ; given_sigma_ion_O2_to_O_cm2(1275) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1276) =     8.72 ; given_sigma_ion_O2_to_O_cm2(1276) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1277) =    10.33 ; given_sigma_ion_O2_to_O_cm2(1277) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1278) =    11.28 ; given_sigma_ion_O2_to_O_cm2(1278) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1279) =    10.95 ; given_sigma_ion_O2_to_O_cm2(1279) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1280) =    10.07 ; given_sigma_ion_O2_to_O_cm2(1280) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1281) =     9.29 ; given_sigma_ion_O2_to_O_cm2(1281) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1282) =     8.00 ; given_sigma_ion_O2_to_O_cm2(1282) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1283) =     6.81 ; given_sigma_ion_O2_to_O_cm2(1283) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1284) =     5.52 ; given_sigma_ion_O2_to_O_cm2(1284) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1285) =     4.13 ; given_sigma_ion_O2_to_O_cm2(1285) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1286) =     4.80 ; given_sigma_ion_O2_to_O_cm2(1286) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1287) =     6.00 ; given_sigma_ion_O2_to_O_cm2(1287) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1288) =     7.20 ; given_sigma_ion_O2_to_O_cm2(1288) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1289) =     8.40 ; given_sigma_ion_O2_to_O_cm2(1289) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1290) =     7.81 ; given_sigma_ion_O2_to_O_cm2(1290) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1291) =     6.77 ; given_sigma_ion_O2_to_O_cm2(1291) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1292) =     5.73 ; given_sigma_ion_O2_to_O_cm2(1292) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1293) =     4.69 ; given_sigma_ion_O2_to_O_cm2(1293) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1294) =     6.11 ; given_sigma_ion_O2_to_O_cm2(1294) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1295) =     7.84 ; given_sigma_ion_O2_to_O_cm2(1295) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1296) =     7.85 ; given_sigma_ion_O2_to_O_cm2(1296) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1297) =     5.76 ; given_sigma_ion_O2_to_O_cm2(1297) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1298) =     6.34 ; given_sigma_ion_O2_to_O_cm2(1298) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1299) =     4.29 ; given_sigma_ion_O2_to_O_cm2(1299) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1300) =     3.25 ; given_sigma_ion_O2_to_O_cm2(1300) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1301) =     3.41 ; given_sigma_ion_O2_to_O_cm2(1301) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1302) =     3.83 ; given_sigma_ion_O2_to_O_cm2(1302) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1303) =     4.25 ; given_sigma_ion_O2_to_O_cm2(1303) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1304) =     4.43 ; given_sigma_ion_O2_to_O_cm2(1304) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1305) =     4.29 ; given_sigma_ion_O2_to_O_cm2(1305) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1306) =     3.74 ; given_sigma_ion_O2_to_O_cm2(1306) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1307) =     3.18 ; given_sigma_ion_O2_to_O_cm2(1307) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1308) =     2.95 ; given_sigma_ion_O2_to_O_cm2(1308) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1309) =     3.46 ; given_sigma_ion_O2_to_O_cm2(1309) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1310) =     4.31 ; given_sigma_ion_O2_to_O_cm2(1310) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1311) =     5.17 ; given_sigma_ion_O2_to_O_cm2(1311) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1312) =     5.52 ; given_sigma_ion_O2_to_O_cm2(1312) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1313) =     6.40 ; given_sigma_ion_O2_to_O_cm2(1313) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1314) =     5.83 ; given_sigma_ion_O2_to_O_cm2(1314) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1315) =     7.57 ; given_sigma_ion_O2_to_O_cm2(1315) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1316) =     7.03 ; given_sigma_ion_O2_to_O_cm2(1316) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1317) =     5.72 ; given_sigma_ion_O2_to_O_cm2(1317) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1318) =     5.29 ; given_sigma_ion_O2_to_O_cm2(1318) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1319) =     8.12 ; given_sigma_ion_O2_to_O_cm2(1319) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1320) =     8.64 ; given_sigma_ion_O2_to_O_cm2(1320) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1321) =     8.10 ; given_sigma_ion_O2_to_O_cm2(1321) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1322) =     6.81 ; given_sigma_ion_O2_to_O_cm2(1322) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1323) =     5.60 ; given_sigma_ion_O2_to_O_cm2(1323) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1324) =     4.46 ; given_sigma_ion_O2_to_O_cm2(1324) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1325) =     3.41 ; given_sigma_ion_O2_to_O_cm2(1325) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1326) =     2.62 ; given_sigma_ion_O2_to_O_cm2(1326) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1327) =     4.24 ; given_sigma_ion_O2_to_O_cm2(1327) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1328) =     5.34 ; given_sigma_ion_O2_to_O_cm2(1328) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1329) =     7.02 ; given_sigma_ion_O2_to_O_cm2(1329) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1330) =     8.45 ; given_sigma_ion_O2_to_O_cm2(1330) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1331) =     9.91 ; given_sigma_ion_O2_to_O_cm2(1331) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1332) =    11.40 ; given_sigma_ion_O2_to_O_cm2(1332) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1333) =    11.10 ; given_sigma_ion_O2_to_O_cm2(1333) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1334) =    10.49 ; given_sigma_ion_O2_to_O_cm2(1334) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1335) =    10.80 ; given_sigma_ion_O2_to_O_cm2(1335) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1336) =     9.58 ; given_sigma_ion_O2_to_O_cm2(1336) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1337) =     8.27 ; given_sigma_ion_O2_to_O_cm2(1337) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1338) =     7.57 ; given_sigma_ion_O2_to_O_cm2(1338) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1339) =     6.92 ; given_sigma_ion_O2_to_O_cm2(1339) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1340) =     6.03 ; given_sigma_ion_O2_to_O_cm2(1340) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1341) =     5.53 ; given_sigma_ion_O2_to_O_cm2(1341) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1342) =     4.09 ; given_sigma_ion_O2_to_O_cm2(1342) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1343) =     3.22 ; given_sigma_ion_O2_to_O_cm2(1343) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1344) =     3.56 ; given_sigma_ion_O2_to_O_cm2(1344) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1345) =     4.48 ; given_sigma_ion_O2_to_O_cm2(1345) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1346) =     5.49 ; given_sigma_ion_O2_to_O_cm2(1346) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1347) =     5.92 ; given_sigma_ion_O2_to_O_cm2(1347) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1348) =     4.66 ; given_sigma_ion_O2_to_O_cm2(1348) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1349) =     2.86 ; given_sigma_ion_O2_to_O_cm2(1349) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1350) =     5.71 ; given_sigma_ion_O2_to_O_cm2(1350) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1351) =     9.95 ; given_sigma_ion_O2_to_O_cm2(1351) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1352) =    14.48 ; given_sigma_ion_O2_to_O_cm2(1352) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1353) =    16.38 ; given_sigma_ion_O2_to_O_cm2(1353) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1354) =    15.88 ; given_sigma_ion_O2_to_O_cm2(1354) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1355) =    14.50 ; given_sigma_ion_O2_to_O_cm2(1355) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1356) =    13.12 ; given_sigma_ion_O2_to_O_cm2(1356) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1357) =    11.56 ; given_sigma_ion_O2_to_O_cm2(1357) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1358) =     7.51 ; given_sigma_ion_O2_to_O_cm2(1358) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1359) =     5.07 ; given_sigma_ion_O2_to_O_cm2(1359) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1360) =     2.31 ; given_sigma_ion_O2_to_O_cm2(1360) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1361) =     3.41 ; given_sigma_ion_O2_to_O_cm2(1361) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1362) =     4.49 ; given_sigma_ion_O2_to_O_cm2(1362) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1363) =     5.51 ; given_sigma_ion_O2_to_O_cm2(1363) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1364) =     6.90 ; given_sigma_ion_O2_to_O_cm2(1364) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1365) =     7.74 ; given_sigma_ion_O2_to_O_cm2(1365) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1366) =     8.32 ; given_sigma_ion_O2_to_O_cm2(1366) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1367) =     5.70 ; given_sigma_ion_O2_to_O_cm2(1367) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1368) =    10.62 ; given_sigma_ion_O2_to_O_cm2(1368) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1369) =    13.97 ; given_sigma_ion_O2_to_O_cm2(1369) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1370) =    17.01 ; given_sigma_ion_O2_to_O_cm2(1370) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1371) =    20.59 ; given_sigma_ion_O2_to_O_cm2(1371) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1372) =    15.02 ; given_sigma_ion_O2_to_O_cm2(1372) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1373) =    10.03 ; given_sigma_ion_O2_to_O_cm2(1373) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1374) =     5.60 ; given_sigma_ion_O2_to_O_cm2(1374) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1375) =     4.00 ; given_sigma_ion_O2_to_O_cm2(1375) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1376) =     2.46 ; given_sigma_ion_O2_to_O_cm2(1376) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1377) =     2.67 ; given_sigma_ion_O2_to_O_cm2(1377) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1378) =     2.86 ; given_sigma_ion_O2_to_O_cm2(1378) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1379) =     2.45 ; given_sigma_ion_O2_to_O_cm2(1379) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1380) =     2.59 ; given_sigma_ion_O2_to_O_cm2(1380) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1381) =     2.91 ; given_sigma_ion_O2_to_O_cm2(1381) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1382) =     2.74 ; given_sigma_ion_O2_to_O_cm2(1382) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1383) =    16.69 ; given_sigma_ion_O2_to_O_cm2(1383) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1384) =    17.86 ; given_sigma_ion_O2_to_O_cm2(1384) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1385) =    15.08 ; given_sigma_ion_O2_to_O_cm2(1385) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1386) =     9.41 ; given_sigma_ion_O2_to_O_cm2(1386) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1387) =    15.09 ; given_sigma_ion_O2_to_O_cm2(1387) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1388) =    22.31 ; given_sigma_ion_O2_to_O_cm2(1388) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1389) =    18.70 ; given_sigma_ion_O2_to_O_cm2(1389) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1390) =    16.38 ; given_sigma_ion_O2_to_O_cm2(1390) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1391) =    15.64 ; given_sigma_ion_O2_to_O_cm2(1391) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1392) =    12.54 ; given_sigma_ion_O2_to_O_cm2(1392) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1393) =     9.39 ; given_sigma_ion_O2_to_O_cm2(1393) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1394) =     6.19 ; given_sigma_ion_O2_to_O_cm2(1394) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1395) =     4.25 ; given_sigma_ion_O2_to_O_cm2(1395) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1396) =     2.94 ; given_sigma_ion_O2_to_O_cm2(1396) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1397) =     6.05 ; given_sigma_ion_O2_to_O_cm2(1397) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1398) =    10.81 ; given_sigma_ion_O2_to_O_cm2(1398) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1399) =    15.64 ; given_sigma_ion_O2_to_O_cm2(1399) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1400) =    20.53 ; given_sigma_ion_O2_to_O_cm2(1400) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1401) =    23.50 ; given_sigma_ion_O2_to_O_cm2(1401) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1402) =    24.50 ; given_sigma_ion_O2_to_O_cm2(1402) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1403) =    30.52 ; given_sigma_ion_O2_to_O_cm2(1403) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1404) =    35.61 ; given_sigma_ion_O2_to_O_cm2(1404) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1405) =    38.70 ; given_sigma_ion_O2_to_O_cm2(1405) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1406) =    36.59 ; given_sigma_ion_O2_to_O_cm2(1406) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1407) =    31.54 ; given_sigma_ion_O2_to_O_cm2(1407) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1408) =    26.83 ; given_sigma_ion_O2_to_O_cm2(1408) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1409) =    22.47 ; given_sigma_ion_O2_to_O_cm2(1409) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1410) =    18.44 ; given_sigma_ion_O2_to_O_cm2(1410) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1411) =    14.75 ; given_sigma_ion_O2_to_O_cm2(1411) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1412) =    12.05 ; given_sigma_ion_O2_to_O_cm2(1412) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1413) =     8.40 ; given_sigma_ion_O2_to_O_cm2(1413) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1414) =     6.76 ; given_sigma_ion_O2_to_O_cm2(1414) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1415) =     6.56 ; given_sigma_ion_O2_to_O_cm2(1415) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1416) =     5.28 ; given_sigma_ion_O2_to_O_cm2(1416) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1417) =     2.90 ; given_sigma_ion_O2_to_O_cm2(1417) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1418) =     8.25 ; given_sigma_ion_O2_to_O_cm2(1418) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1419) =    22.38 ; given_sigma_ion_O2_to_O_cm2(1419) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1420) =    35.39 ; given_sigma_ion_O2_to_O_cm2(1420) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1421) =    43.84 ; given_sigma_ion_O2_to_O_cm2(1421) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1422) =    38.92 ; given_sigma_ion_O2_to_O_cm2(1422) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1423) =    30.75 ; given_sigma_ion_O2_to_O_cm2(1423) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1424) =    22.63 ; given_sigma_ion_O2_to_O_cm2(1424) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1425) =    14.55 ; given_sigma_ion_O2_to_O_cm2(1425) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1426) =    10.68 ; given_sigma_ion_O2_to_O_cm2(1426) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1427) =     6.50 ; given_sigma_ion_O2_to_O_cm2(1427) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1428) =     1.69 ; given_sigma_ion_O2_to_O_cm2(1428) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1429) =     3.14 ; given_sigma_ion_O2_to_O_cm2(1429) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1430) =     6.75 ; given_sigma_ion_O2_to_O_cm2(1430) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1431) =    13.97 ; given_sigma_ion_O2_to_O_cm2(1431) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1432) =    21.18 ; given_sigma_ion_O2_to_O_cm2(1432) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1433) =    28.40 ; given_sigma_ion_O2_to_O_cm2(1433) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1434) =    35.62 ; given_sigma_ion_O2_to_O_cm2(1434) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1435) =    42.12 ; given_sigma_ion_O2_to_O_cm2(1435) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1436) =    27.28 ; given_sigma_ion_O2_to_O_cm2(1436) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1437) =    22.33 ; given_sigma_ion_O2_to_O_cm2(1437) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1438) =    29.65 ; given_sigma_ion_O2_to_O_cm2(1438) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1439) =    24.31 ; given_sigma_ion_O2_to_O_cm2(1439) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1440) =    19.25 ; given_sigma_ion_O2_to_O_cm2(1440) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1441) =    17.30 ; given_sigma_ion_O2_to_O_cm2(1441) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1442) =    14.46 ; given_sigma_ion_O2_to_O_cm2(1442) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1443) =    11.72 ; given_sigma_ion_O2_to_O_cm2(1443) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1444) =     9.94 ; given_sigma_ion_O2_to_O_cm2(1444) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1445) =     5.70 ; given_sigma_ion_O2_to_O_cm2(1445) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1446) =     1.73 ; given_sigma_ion_O2_to_O_cm2(1446) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1447) =     4.17 ; given_sigma_ion_O2_to_O_cm2(1447) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1448) =     6.68 ; given_sigma_ion_O2_to_O_cm2(1448) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1449) =     9.27 ; given_sigma_ion_O2_to_O_cm2(1449) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1450) =    11.40 ; given_sigma_ion_O2_to_O_cm2(1450) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1451) =     7.19 ; given_sigma_ion_O2_to_O_cm2(1451) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1452) =     5.03 ; given_sigma_ion_O2_to_O_cm2(1452) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1453) =     7.55 ; given_sigma_ion_O2_to_O_cm2(1453) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1454) =    20.85 ; given_sigma_ion_O2_to_O_cm2(1454) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1455) =    35.24 ; given_sigma_ion_O2_to_O_cm2(1455) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1456) =    42.84 ; given_sigma_ion_O2_to_O_cm2(1456) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1457) =    38.11 ; given_sigma_ion_O2_to_O_cm2(1457) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1458) =    33.50 ; given_sigma_ion_O2_to_O_cm2(1458) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1459) =    29.02 ; given_sigma_ion_O2_to_O_cm2(1459) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1460) =    24.65 ; given_sigma_ion_O2_to_O_cm2(1460) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1461) =    20.41 ; given_sigma_ion_O2_to_O_cm2(1461) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1462) =    12.28 ; given_sigma_ion_O2_to_O_cm2(1462) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1463) =     8.40 ; given_sigma_ion_O2_to_O_cm2(1463) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1464) =     4.65 ; given_sigma_ion_O2_to_O_cm2(1464) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1465) =     1.73 ; given_sigma_ion_O2_to_O_cm2(1465) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1466) =     7.43 ; given_sigma_ion_O2_to_O_cm2(1466) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1467) =    11.81 ; given_sigma_ion_O2_to_O_cm2(1467) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1468) =    15.87 ; given_sigma_ion_O2_to_O_cm2(1468) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1469) =    19.59 ; given_sigma_ion_O2_to_O_cm2(1469) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1470) =    24.31 ; given_sigma_ion_O2_to_O_cm2(1470) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1471) =    14.14 ; given_sigma_ion_O2_to_O_cm2(1471) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1472) =     7.32 ; given_sigma_ion_O2_to_O_cm2(1472) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1473) =     2.02 ; given_sigma_ion_O2_to_O_cm2(1473) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1474) =    10.67 ; given_sigma_ion_O2_to_O_cm2(1474) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1475) =    18.49 ; given_sigma_ion_O2_to_O_cm2(1475) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1476) =    18.27 ; given_sigma_ion_O2_to_O_cm2(1476) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1477) =    17.58 ; given_sigma_ion_O2_to_O_cm2(1477) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1478) =    16.70 ; given_sigma_ion_O2_to_O_cm2(1478) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1479) =    15.62 ; given_sigma_ion_O2_to_O_cm2(1479) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1480) =    13.63 ; given_sigma_ion_O2_to_O_cm2(1480) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1481) =    11.65 ; given_sigma_ion_O2_to_O_cm2(1481) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1482) =     9.68 ; given_sigma_ion_O2_to_O_cm2(1482) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1483) =     7.73 ; given_sigma_ion_O2_to_O_cm2(1483) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1484) =     5.79 ; given_sigma_ion_O2_to_O_cm2(1484) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1485) =     3.86 ; given_sigma_ion_O2_to_O_cm2(1485) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1486) =     1.94 ; given_sigma_ion_O2_to_O_cm2(1486) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1487) =     8.37 ; given_sigma_ion_O2_to_O_cm2(1487) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1488) =    14.97 ; given_sigma_ion_O2_to_O_cm2(1488) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1489) =    21.74 ; given_sigma_ion_O2_to_O_cm2(1489) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1490) =    28.67 ; given_sigma_ion_O2_to_O_cm2(1490) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1491) =    35.77 ; given_sigma_ion_O2_to_O_cm2(1491) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1492) =    40.11 ; given_sigma_ion_O2_to_O_cm2(1492) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1493) =    36.29 ; given_sigma_ion_O2_to_O_cm2(1493) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1494) =    26.79 ; given_sigma_ion_O2_to_O_cm2(1494) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1495) =    17.32 ; given_sigma_ion_O2_to_O_cm2(1495) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1496) =     7.90 ; given_sigma_ion_O2_to_O_cm2(1496) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1497) =     4.14 ; given_sigma_ion_O2_to_O_cm2(1497) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1498) =     5.10 ; given_sigma_ion_O2_to_O_cm2(1498) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1499) =     6.36 ; given_sigma_ion_O2_to_O_cm2(1499) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1500) =     5.72 ; given_sigma_ion_O2_to_O_cm2(1500) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1501) =     4.54 ; given_sigma_ion_O2_to_O_cm2(1501) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1502) =     2.76 ; given_sigma_ion_O2_to_O_cm2(1502) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1503) =     4.03 ; given_sigma_ion_O2_to_O_cm2(1503) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1504) =     2.58 ; given_sigma_ion_O2_to_O_cm2(1504) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1505) =     1.07 ; given_sigma_ion_O2_to_O_cm2(1505) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1506) =     1.91 ; given_sigma_ion_O2_to_O_cm2(1506) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1507) =     2.87 ; given_sigma_ion_O2_to_O_cm2(1507) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1508) =     7.45 ; given_sigma_ion_O2_to_O_cm2(1508) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1509) =     9.74 ; given_sigma_ion_O2_to_O_cm2(1509) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1510) =    12.95 ; given_sigma_ion_O2_to_O_cm2(1510) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1511) =    19.61 ; given_sigma_ion_O2_to_O_cm2(1511) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1512) =    18.79 ; given_sigma_ion_O2_to_O_cm2(1512) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1513) =    20.83 ; given_sigma_ion_O2_to_O_cm2(1513) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1514) =    18.03 ; given_sigma_ion_O2_to_O_cm2(1514) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1515) =    12.57 ; given_sigma_ion_O2_to_O_cm2(1515) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1516) =     7.29 ; given_sigma_ion_O2_to_O_cm2(1516) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1517) =     2.16 ; given_sigma_ion_O2_to_O_cm2(1517) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1518) =     1.16 ; given_sigma_ion_O2_to_O_cm2(1518) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1519) =     1.18 ; given_sigma_ion_O2_to_O_cm2(1519) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1520) =     1.20 ; given_sigma_ion_O2_to_O_cm2(1520) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1521) =     1.22 ; given_sigma_ion_O2_to_O_cm2(1521) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1522) =     2.24 ; given_sigma_ion_O2_to_O_cm2(1522) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1523) =     3.29 ; given_sigma_ion_O2_to_O_cm2(1523) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1524) =     4.37 ; given_sigma_ion_O2_to_O_cm2(1524) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1525) =     5.48 ; given_sigma_ion_O2_to_O_cm2(1525) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1526) =     4.85 ; given_sigma_ion_O2_to_O_cm2(1526) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1527) =     5.44 ; given_sigma_ion_O2_to_O_cm2(1527) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1528) =     4.65 ; given_sigma_ion_O2_to_O_cm2(1528) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1529) =     2.69 ; given_sigma_ion_O2_to_O_cm2(1529) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1530) =     1.16 ; given_sigma_ion_O2_to_O_cm2(1530) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1531) =     1.22 ; given_sigma_ion_O2_to_O_cm2(1531) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1532) =     1.51 ; given_sigma_ion_O2_to_O_cm2(1532) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1533) =     1.18 ; given_sigma_ion_O2_to_O_cm2(1533) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1534) =     1.24 ; given_sigma_ion_O2_to_O_cm2(1534) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1535) =     1.17 ; given_sigma_ion_O2_to_O_cm2(1535) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1536) =     1.15 ; given_sigma_ion_O2_to_O_cm2(1536) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1537) =     1.07 ; given_sigma_ion_O2_to_O_cm2(1537) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1538) =     1.03 ; given_sigma_ion_O2_to_O_cm2(1538) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1539) =     0.91 ; given_sigma_ion_O2_to_O_cm2(1539) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1540) =     0.84 ; given_sigma_ion_O2_to_O_cm2(1540) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1541) =     0.95 ; given_sigma_ion_O2_to_O_cm2(1541) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1542) =     1.03 ; given_sigma_ion_O2_to_O_cm2(1542) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1543) =     0.86 ; given_sigma_ion_O2_to_O_cm2(1543) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1544) =     1.12 ; given_sigma_ion_O2_to_O_cm2(1544) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1545) =     1.30 ; given_sigma_ion_O2_to_O_cm2(1545) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1546) =     0.88 ; given_sigma_ion_O2_to_O_cm2(1546) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1547) =     1.02 ; given_sigma_ion_O2_to_O_cm2(1547) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1548) =     0.79 ; given_sigma_ion_O2_to_O_cm2(1548) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1549) =     1.17 ; given_sigma_ion_O2_to_O_cm2(1549) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1550) =     0.70 ; given_sigma_ion_O2_to_O_cm2(1550) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1551) =     0.86 ; given_sigma_ion_O2_to_O_cm2(1551) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1552) =     1.09 ; given_sigma_ion_O2_to_O_cm2(1552) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1553) =     0.77 ; given_sigma_ion_O2_to_O_cm2(1553) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1554) =     0.85 ; given_sigma_ion_O2_to_O_cm2(1554) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1555) =     1.00 ; given_sigma_ion_O2_to_O_cm2(1555) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1556) =     0.90 ; given_sigma_ion_O2_to_O_cm2(1556) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1557) =     0.83 ; given_sigma_ion_O2_to_O_cm2(1557) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1558) =     1.11 ; given_sigma_ion_O2_to_O_cm2(1558) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1559) =     0.97 ; given_sigma_ion_O2_to_O_cm2(1559) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1560) =     1.12 ; given_sigma_ion_O2_to_O_cm2(1560) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1561) =     1.00 ; given_sigma_ion_O2_to_O_cm2(1561) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1562) =     0.88 ; given_sigma_ion_O2_to_O_cm2(1562) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1563) =     1.09 ; given_sigma_ion_O2_to_O_cm2(1563) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1564) =     1.24 ; given_sigma_ion_O2_to_O_cm2(1564) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1565) =     0.94 ; given_sigma_ion_O2_to_O_cm2(1565) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1566) =     0.65 ; given_sigma_ion_O2_to_O_cm2(1566) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1567) =     0.94 ; given_sigma_ion_O2_to_O_cm2(1567) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1568) =     1.15 ; given_sigma_ion_O2_to_O_cm2(1568) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1569) =     1.05 ; given_sigma_ion_O2_to_O_cm2(1569) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1570) =     0.00 ; given_sigma_ion_O2_to_O_cm2(1570) =    0.00
  given_sigma_ion_O2_to_O2_cm2(1571) =     0.00 ; given_sigma_ion_O2_to_O_cm2(1571) =    0.00







 given_sigma_tot_O2_cm2(   1) =       0.07
 given_sigma_tot_O2_cm2(   2) =       0.13
 given_sigma_tot_O2_cm2(   3) =       0.13
 given_sigma_tot_O2_cm2(   4) =       0.14
 given_sigma_tot_O2_cm2(   5) =       0.14
 given_sigma_tot_O2_cm2(   6) =       0.15
 given_sigma_tot_O2_cm2(   7) =       0.15
 given_sigma_tot_O2_cm2(   8) =       0.16
 given_sigma_tot_O2_cm2(   9) =       0.16
 given_sigma_tot_O2_cm2(  10) =       0.19
 given_sigma_tot_O2_cm2(  11) =       0.21
 given_sigma_tot_O2_cm2(  12) =       0.28
 given_sigma_tot_O2_cm2(  13) =       0.28
 given_sigma_tot_O2_cm2(  14) =       0.31
 given_sigma_tot_O2_cm2(  15) =       0.31
 given_sigma_tot_O2_cm2(  16) =       0.31
 given_sigma_tot_O2_cm2(  17) =       0.32
 given_sigma_tot_O2_cm2(  18) =       0.33
 given_sigma_tot_O2_cm2(  19) =       0.34
 given_sigma_tot_O2_cm2(  20) =       0.35
 given_sigma_tot_O2_cm2(  21) =       0.36
 given_sigma_tot_O2_cm2(  22) =       0.36
 given_sigma_tot_O2_cm2(  23) =       0.38
 given_sigma_tot_O2_cm2(  24) =       0.39
 given_sigma_tot_O2_cm2(  25) =       0.42
 given_sigma_tot_O2_cm2(  26) =       0.43
 given_sigma_tot_O2_cm2(  27) =       0.44
 given_sigma_tot_O2_cm2(  28) =       0.44
 given_sigma_tot_O2_cm2(  29) =       0.45
 given_sigma_tot_O2_cm2(  30) =       0.47
 given_sigma_tot_O2_cm2(  31) =       0.48
 given_sigma_tot_O2_cm2(  32) =       0.50
 given_sigma_tot_O2_cm2(  33) =       0.52
 given_sigma_tot_O2_cm2(  34) =       0.53
 given_sigma_tot_O2_cm2(  35) =       0.54
 given_sigma_tot_O2_cm2(  36) =       0.54
 given_sigma_tot_O2_cm2(  37) =       0.55
 given_sigma_tot_O2_cm2(  38) =       0.56
 given_sigma_tot_O2_cm2(  39) =       0.58
 given_sigma_tot_O2_cm2(  40) =       0.59
 given_sigma_tot_O2_cm2(  41) =       0.60
 given_sigma_tot_O2_cm2(  42) =       0.61
 given_sigma_tot_O2_cm2(  43) =       0.62
 given_sigma_tot_O2_cm2(  44) =       0.62
 given_sigma_tot_O2_cm2(  45) =       0.65
 given_sigma_tot_O2_cm2(  46) =       0.67
 given_sigma_tot_O2_cm2(  47) =       0.68
 given_sigma_tot_O2_cm2(  48) =       0.70
 given_sigma_tot_O2_cm2(  49) =       0.70
 given_sigma_tot_O2_cm2(  50) =       0.72
 given_sigma_tot_O2_cm2(  51) =       0.73
 given_sigma_tot_O2_cm2(  52) =       0.74
 given_sigma_tot_O2_cm2(  53) =       0.75
 given_sigma_tot_O2_cm2(  54) =       0.75
 given_sigma_tot_O2_cm2(  55) =       0.76
 given_sigma_tot_O2_cm2(  56) =       0.76
 given_sigma_tot_O2_cm2(  57) =       0.77
 given_sigma_tot_O2_cm2(  58) =       0.78
 given_sigma_tot_O2_cm2(  59) =       0.79
 given_sigma_tot_O2_cm2(  60) =       0.80
 given_sigma_tot_O2_cm2(  61) =       0.82
 given_sigma_tot_O2_cm2(  62) =       0.84
 given_sigma_tot_O2_cm2(  63) =       0.84
 given_sigma_tot_O2_cm2(  64) =       0.86
 given_sigma_tot_O2_cm2(  65) =       0.89
 given_sigma_tot_O2_cm2(  66) =       0.89
 given_sigma_tot_O2_cm2(  67) =       0.90
 given_sigma_tot_O2_cm2(  68) =       0.93
 given_sigma_tot_O2_cm2(  69) =       0.94
 given_sigma_tot_O2_cm2(  70) =       0.96
 given_sigma_tot_O2_cm2(  71) =       0.98
 given_sigma_tot_O2_cm2(  72) =       0.99
 given_sigma_tot_O2_cm2(  73) =       1.00
 given_sigma_tot_O2_cm2(  74) =       1.00
 given_sigma_tot_O2_cm2(  75) =       1.03
 given_sigma_tot_O2_cm2(  76) =       1.04
 given_sigma_tot_O2_cm2(  77) =       1.04
 given_sigma_tot_O2_cm2(  78) =       1.06
 given_sigma_tot_O2_cm2(  79) =       1.06
 given_sigma_tot_O2_cm2(  80) =       1.09
 given_sigma_tot_O2_cm2(  81) =       1.11
 given_sigma_tot_O2_cm2(  82) =       1.12
 given_sigma_tot_O2_cm2(  83) =       1.14
 given_sigma_tot_O2_cm2(  84) =       1.14
 given_sigma_tot_O2_cm2(  85) =       1.15
 given_sigma_tot_O2_cm2(  86) =       1.16
 given_sigma_tot_O2_cm2(  87) =       1.17
 given_sigma_tot_O2_cm2(  88) =       1.18
 given_sigma_tot_O2_cm2(  89) =       1.19
 given_sigma_tot_O2_cm2(  90) =       1.21
 given_sigma_tot_O2_cm2(  91) =       1.21
 given_sigma_tot_O2_cm2(  92) =       1.22
 given_sigma_tot_O2_cm2(  93) =       1.23
 given_sigma_tot_O2_cm2(  94) =       1.24
 given_sigma_tot_O2_cm2(  95) =       1.27
 given_sigma_tot_O2_cm2(  96) =       1.28
 given_sigma_tot_O2_cm2(  97) =       1.29
 given_sigma_tot_O2_cm2(  98) =       1.30
 given_sigma_tot_O2_cm2(  99) =       1.31
 given_sigma_tot_O2_cm2( 100) =       1.32
 given_sigma_tot_O2_cm2( 101) =       1.33
 given_sigma_tot_O2_cm2( 102) =       1.34
 given_sigma_tot_O2_cm2( 103) =       1.36
 given_sigma_tot_O2_cm2( 104) =       1.36
 given_sigma_tot_O2_cm2( 105) =       1.38
 given_sigma_tot_O2_cm2( 106) =       1.39
 given_sigma_tot_O2_cm2( 107) =       1.40
 given_sigma_tot_O2_cm2( 108) =       1.41
 given_sigma_tot_O2_cm2( 109) =       1.43
 given_sigma_tot_O2_cm2( 110) =       1.45
 given_sigma_tot_O2_cm2( 111) =       1.46
 given_sigma_tot_O2_cm2( 112) =       1.47
 given_sigma_tot_O2_cm2( 113) =       1.49
 given_sigma_tot_O2_cm2( 114) =       1.50
 given_sigma_tot_O2_cm2( 115) =       1.51
 given_sigma_tot_O2_cm2( 116) =       1.53
 given_sigma_tot_O2_cm2( 117) =       1.53
 given_sigma_tot_O2_cm2( 118) =       1.55
 given_sigma_tot_O2_cm2( 119) =       1.56
 given_sigma_tot_O2_cm2( 120) =       1.57
 given_sigma_tot_O2_cm2( 121) =       1.58
 given_sigma_tot_O2_cm2( 122) =       1.60
 given_sigma_tot_O2_cm2( 123) =       1.61
 given_sigma_tot_O2_cm2( 124) =       1.62
 given_sigma_tot_O2_cm2( 125) =       1.64
 given_sigma_tot_O2_cm2( 126) =       1.65
 given_sigma_tot_O2_cm2( 127) =       1.67
 given_sigma_tot_O2_cm2( 128) =       1.69
 given_sigma_tot_O2_cm2( 129) =       1.71
 given_sigma_tot_O2_cm2( 130) =       1.71
 given_sigma_tot_O2_cm2( 131) =       1.73
 given_sigma_tot_O2_cm2( 132) =       1.74
 given_sigma_tot_O2_cm2( 133) =       1.77
 given_sigma_tot_O2_cm2( 134) =       1.79
 given_sigma_tot_O2_cm2( 135) =       1.80
 given_sigma_tot_O2_cm2( 136) =       1.82
 given_sigma_tot_O2_cm2( 137) =       1.83
 given_sigma_tot_O2_cm2( 138) =       1.85
 given_sigma_tot_O2_cm2( 139) =       1.86
 given_sigma_tot_O2_cm2( 140) =       1.87
 given_sigma_tot_O2_cm2( 141) =       1.88
 given_sigma_tot_O2_cm2( 142) =       1.91
 given_sigma_tot_O2_cm2( 143) =       1.92
 given_sigma_tot_O2_cm2( 144) =       1.96
 given_sigma_tot_O2_cm2( 145) =       1.98
 given_sigma_tot_O2_cm2( 146) =       1.99
 given_sigma_tot_O2_cm2( 147) =       2.00
 given_sigma_tot_O2_cm2( 148) =       2.02
 given_sigma_tot_O2_cm2( 149) =       2.04
 given_sigma_tot_O2_cm2( 150) =       2.05
 given_sigma_tot_O2_cm2( 151) =       2.06
 given_sigma_tot_O2_cm2( 152) =       2.08
 given_sigma_tot_O2_cm2( 153) =       2.10
 given_sigma_tot_O2_cm2( 154) =       2.11
 given_sigma_tot_O2_cm2( 155) =       2.13
 given_sigma_tot_O2_cm2( 156) =       2.15
 given_sigma_tot_O2_cm2( 157) =       2.16
 given_sigma_tot_O2_cm2( 158) =       2.18
 given_sigma_tot_O2_cm2( 159) =       2.18
 given_sigma_tot_O2_cm2( 160) =       2.20
 given_sigma_tot_O2_cm2( 161) =       2.21
 given_sigma_tot_O2_cm2( 162) =       2.24
 given_sigma_tot_O2_cm2( 163) =       2.25
 given_sigma_tot_O2_cm2( 164) =       2.27
 given_sigma_tot_O2_cm2( 165) =       2.29
 given_sigma_tot_O2_cm2( 166) =       2.31
 given_sigma_tot_O2_cm2( 167) =       2.34
 given_sigma_tot_O2_cm2( 168) =       2.37
 given_sigma_tot_O2_cm2( 169) =       2.41
 given_sigma_tot_O2_cm2( 170) =       2.42
 given_sigma_tot_O2_cm2( 171) =       2.43
 given_sigma_tot_O2_cm2( 172) =       2.44
 given_sigma_tot_O2_cm2( 173) =       2.46
 given_sigma_tot_O2_cm2( 174) =       2.47
 given_sigma_tot_O2_cm2( 175) =       2.50
 given_sigma_tot_O2_cm2( 176) =       2.52
 given_sigma_tot_O2_cm2( 177) =       2.57
 given_sigma_tot_O2_cm2( 178) =       2.59
 given_sigma_tot_O2_cm2( 179) =       2.60
 given_sigma_tot_O2_cm2( 180) =       2.65
 given_sigma_tot_O2_cm2( 181) =       2.66
 given_sigma_tot_O2_cm2( 182) =       2.68
 given_sigma_tot_O2_cm2( 183) =       2.73
 given_sigma_tot_O2_cm2( 184) =       2.75
 given_sigma_tot_O2_cm2( 185) =       2.78
 given_sigma_tot_O2_cm2( 186) =       2.79
 given_sigma_tot_O2_cm2( 187) =       2.81
 given_sigma_tot_O2_cm2( 188) =       2.89
 given_sigma_tot_O2_cm2( 189) =       2.94
 given_sigma_tot_O2_cm2( 190) =       2.95
 given_sigma_tot_O2_cm2( 191) =       2.96
 given_sigma_tot_O2_cm2( 192) =       3.02
 given_sigma_tot_O2_cm2( 193) =       3.04
 given_sigma_tot_O2_cm2( 194) =       3.08
 given_sigma_tot_O2_cm2( 195) =       3.10
 given_sigma_tot_O2_cm2( 196) =       3.15
 given_sigma_tot_O2_cm2( 197) =       3.24
 given_sigma_tot_O2_cm2( 198) =       3.27
 given_sigma_tot_O2_cm2( 199) =       3.32
 given_sigma_tot_O2_cm2( 200) =       3.37
 given_sigma_tot_O2_cm2( 201) =       3.43
 given_sigma_tot_O2_cm2( 202) =       3.49
 given_sigma_tot_O2_cm2( 203) =       3.52
 given_sigma_tot_O2_cm2( 204) =       3.60
 given_sigma_tot_O2_cm2( 205) =       3.78
 given_sigma_tot_O2_cm2( 206) =       3.94
 given_sigma_tot_O2_cm2( 207) =       3.95
 given_sigma_tot_O2_cm2( 208) =       3.97
 given_sigma_tot_O2_cm2( 209) =       3.99
 given_sigma_tot_O2_cm2( 210) =       4.02
 given_sigma_tot_O2_cm2( 211) =       4.03
 given_sigma_tot_O2_cm2( 212) =       4.30
 given_sigma_tot_O2_cm2( 213) =       4.38
 given_sigma_tot_O2_cm2( 214) =       4.40
 given_sigma_tot_O2_cm2( 215) =       4.50
 given_sigma_tot_O2_cm2( 216) =       4.65
 given_sigma_tot_O2_cm2( 217) =       4.73
 given_sigma_tot_O2_cm2( 218) =       4.93
 given_sigma_tot_O2_cm2( 219) =       4.99
 given_sigma_tot_O2_cm2( 220) =       5.05
 given_sigma_tot_O2_cm2( 221) =       5.21
 given_sigma_tot_O2_cm2( 222) =       5.32
 given_sigma_tot_O2_cm2( 223) =       5.33
 given_sigma_tot_O2_cm2( 224) =       5.47
 given_sigma_tot_O2_cm2( 225) =       5.52
 given_sigma_tot_O2_cm2( 226) =       5.62
 given_sigma_tot_O2_cm2( 227) =       5.68
 given_sigma_tot_O2_cm2( 228) =       5.87
 given_sigma_tot_O2_cm2( 229) =       5.91
 given_sigma_tot_O2_cm2( 230) =       6.03
 given_sigma_tot_O2_cm2( 231) =       6.05
 given_sigma_tot_O2_cm2( 232) =       6.31
 given_sigma_tot_O2_cm2( 233) =       6.37
 given_sigma_tot_O2_cm2( 234) =       6.39
 given_sigma_tot_O2_cm2( 235) =       6.53
 given_sigma_tot_O2_cm2( 236) =       6.58
 given_sigma_tot_O2_cm2( 237) =       6.60
 given_sigma_tot_O2_cm2( 238) =       6.63
 given_sigma_tot_O2_cm2( 239) =       6.68
 given_sigma_tot_O2_cm2( 240) =       6.70
 given_sigma_tot_O2_cm2( 241) =       6.77
 given_sigma_tot_O2_cm2( 242) =       6.83
 given_sigma_tot_O2_cm2( 243) =       6.88
 given_sigma_tot_O2_cm2( 244) =       6.89
 given_sigma_tot_O2_cm2( 245) =       6.99
 given_sigma_tot_O2_cm2( 246) =       7.01
 given_sigma_tot_O2_cm2( 247) =       7.03
 given_sigma_tot_O2_cm2( 248) =       7.04
 given_sigma_tot_O2_cm2( 249) =       7.15
 given_sigma_tot_O2_cm2( 250) =       7.15
 given_sigma_tot_O2_cm2( 251) =       7.21
 given_sigma_tot_O2_cm2( 252) =       7.28
 given_sigma_tot_O2_cm2( 253) =       7.31
 given_sigma_tot_O2_cm2( 254) =       7.33
 given_sigma_tot_O2_cm2( 255) =       7.36
 given_sigma_tot_O2_cm2( 256) =       7.39
 given_sigma_tot_O2_cm2( 257) =       7.42
 given_sigma_tot_O2_cm2( 258) =       7.51
 given_sigma_tot_O2_cm2( 259) =       7.53
 given_sigma_tot_O2_cm2( 260) =       7.61
 given_sigma_tot_O2_cm2( 261) =       7.65
 given_sigma_tot_O2_cm2( 262) =       7.67
 given_sigma_tot_O2_cm2( 263) =       7.70
 given_sigma_tot_O2_cm2( 264) =       7.72
 given_sigma_tot_O2_cm2( 265) =       7.74
 given_sigma_tot_O2_cm2( 266) =       7.76
 given_sigma_tot_O2_cm2( 267) =       7.87
 given_sigma_tot_O2_cm2( 268) =       7.89
 given_sigma_tot_O2_cm2( 269) =       7.98
 given_sigma_tot_O2_cm2( 270) =       8.00
 given_sigma_tot_O2_cm2( 271) =       8.04
 given_sigma_tot_O2_cm2( 272) =       8.15
 given_sigma_tot_O2_cm2( 273) =       8.20
 given_sigma_tot_O2_cm2( 274) =       8.22
 given_sigma_tot_O2_cm2( 275) =       8.24
 given_sigma_tot_O2_cm2( 276) =       8.32
 given_sigma_tot_O2_cm2( 277) =       8.35
 given_sigma_tot_O2_cm2( 278) =       8.40
 given_sigma_tot_O2_cm2( 279) =       8.51
 given_sigma_tot_O2_cm2( 280) =       8.51
 given_sigma_tot_O2_cm2( 281) =       8.61
 given_sigma_tot_O2_cm2( 282) =       8.62
 given_sigma_tot_O2_cm2( 283) =       8.68
 given_sigma_tot_O2_cm2( 284) =       8.76
 given_sigma_tot_O2_cm2( 285) =       8.86
 given_sigma_tot_O2_cm2( 286) =       8.94
 given_sigma_tot_O2_cm2( 287) =       9.01
 given_sigma_tot_O2_cm2( 288) =       9.05
 given_sigma_tot_O2_cm2( 289) =       9.13
 given_sigma_tot_O2_cm2( 290) =       9.16
 given_sigma_tot_O2_cm2( 291) =       9.21
 given_sigma_tot_O2_cm2( 292) =       9.30
 given_sigma_tot_O2_cm2( 293) =       9.31
 given_sigma_tot_O2_cm2( 294) =       9.33
 given_sigma_tot_O2_cm2( 295) =       9.39
 given_sigma_tot_O2_cm2( 296) =       9.45
 given_sigma_tot_O2_cm2( 297) =       9.54
 given_sigma_tot_O2_cm2( 298) =       9.55
 given_sigma_tot_O2_cm2( 299) =       9.57
 given_sigma_tot_O2_cm2( 300) =       9.65
 given_sigma_tot_O2_cm2( 301) =       9.73
 given_sigma_tot_O2_cm2( 302) =       9.85
 given_sigma_tot_O2_cm2( 303) =       9.92
 given_sigma_tot_O2_cm2( 304) =       9.94
 given_sigma_tot_O2_cm2( 305) =       9.95
 given_sigma_tot_O2_cm2( 306) =      10.07
 given_sigma_tot_O2_cm2( 307) =      10.08
 given_sigma_tot_O2_cm2( 308) =      10.17
 given_sigma_tot_O2_cm2( 309) =      10.23
 given_sigma_tot_O2_cm2( 310) =      10.30
 given_sigma_tot_O2_cm2( 311) =      10.40
 given_sigma_tot_O2_cm2( 312) =      10.41
 given_sigma_tot_O2_cm2( 313) =      10.45
 given_sigma_tot_O2_cm2( 314) =      10.56
 given_sigma_tot_O2_cm2( 315) =      10.60
 given_sigma_tot_O2_cm2( 316) =      10.68
 given_sigma_tot_O2_cm2( 317) =      10.70
 given_sigma_tot_O2_cm2( 318) =      10.71
 given_sigma_tot_O2_cm2( 319) =      10.86
 given_sigma_tot_O2_cm2( 320) =      10.88
 given_sigma_tot_O2_cm2( 321) =      10.90
 given_sigma_tot_O2_cm2( 322) =      11.00
 given_sigma_tot_O2_cm2( 323) =      11.07
 given_sigma_tot_O2_cm2( 324) =      11.10
 given_sigma_tot_O2_cm2( 325) =      11.15
 given_sigma_tot_O2_cm2( 326) =      11.22
 given_sigma_tot_O2_cm2( 327) =      11.31
 given_sigma_tot_O2_cm2( 328) =      11.41
 given_sigma_tot_O2_cm2( 329) =      11.45
 given_sigma_tot_O2_cm2( 330) =      11.50
 given_sigma_tot_O2_cm2( 331) =      11.54
 given_sigma_tot_O2_cm2( 332) =      11.69
 given_sigma_tot_O2_cm2( 333) =      11.77
 given_sigma_tot_O2_cm2( 334) =      11.89
 given_sigma_tot_O2_cm2( 335) =      11.90
 given_sigma_tot_O2_cm2( 336) =      11.96
 given_sigma_tot_O2_cm2( 337) =      12.06
 given_sigma_tot_O2_cm2( 338) =      12.17
 given_sigma_tot_O2_cm2( 339) =      12.25
 given_sigma_tot_O2_cm2( 340) =      12.34
 given_sigma_tot_O2_cm2( 341) =      12.43
 given_sigma_tot_O2_cm2( 342) =      12.46
 given_sigma_tot_O2_cm2( 343) =      12.52
 given_sigma_tot_O2_cm2( 344) =      12.55
 given_sigma_tot_O2_cm2( 345) =      12.62
 given_sigma_tot_O2_cm2( 346) =      12.73
 given_sigma_tot_O2_cm2( 347) =      12.80
 given_sigma_tot_O2_cm2( 348) =      12.90
 given_sigma_tot_O2_cm2( 349) =      12.98
 given_sigma_tot_O2_cm2( 350) =      13.00
 given_sigma_tot_O2_cm2( 351) =      13.14
 given_sigma_tot_O2_cm2( 352) =      13.20
 given_sigma_tot_O2_cm2( 353) =      13.25
 given_sigma_tot_O2_cm2( 354) =      13.37
 given_sigma_tot_O2_cm2( 355) =      13.40
 given_sigma_tot_O2_cm2( 356) =      13.44
 given_sigma_tot_O2_cm2( 357) =      13.47
 given_sigma_tot_O2_cm2( 358) =      13.55
 given_sigma_tot_O2_cm2( 359) =      13.65
 given_sigma_tot_O2_cm2( 360) =      13.70
 given_sigma_tot_O2_cm2( 361) =      13.80
 given_sigma_tot_O2_cm2( 362) =      14.00
 given_sigma_tot_O2_cm2( 363) =      14.12
 given_sigma_tot_O2_cm2( 364) =      14.18
 given_sigma_tot_O2_cm2( 365) =      14.20
 given_sigma_tot_O2_cm2( 366) =      14.65
 given_sigma_tot_O2_cm2( 367) =      14.70
 given_sigma_tot_O2_cm2( 368) =      14.74
 given_sigma_tot_O2_cm2( 369) =      14.86
 given_sigma_tot_O2_cm2( 370) =      14.91
 given_sigma_tot_O2_cm2( 371) =      15.04
 given_sigma_tot_O2_cm2( 372) =      15.10
 given_sigma_tot_O2_cm2( 373) =      15.13
 given_sigma_tot_O2_cm2( 374) =      15.15
 given_sigma_tot_O2_cm2( 375) =      15.19
 given_sigma_tot_O2_cm2( 376) =      15.24
 given_sigma_tot_O2_cm2( 377) =      15.26
 given_sigma_tot_O2_cm2( 378) =      15.28
 given_sigma_tot_O2_cm2( 379) =      15.37
 given_sigma_tot_O2_cm2( 380) =      15.50
 given_sigma_tot_O2_cm2( 381) =      15.56
 given_sigma_tot_O2_cm2( 382) =      15.60
 given_sigma_tot_O2_cm2( 383) =      15.79
 given_sigma_tot_O2_cm2( 384) =      15.85
 given_sigma_tot_O2_cm2( 385) =      15.90
 given_sigma_tot_O2_cm2( 386) =      15.91
 given_sigma_tot_O2_cm2( 387) =      16.09
 given_sigma_tot_O2_cm2( 388) =      16.14
 given_sigma_tot_O2_cm2( 389) =      16.20
 given_sigma_tot_O2_cm2( 390) =      16.23
 given_sigma_tot_O2_cm2( 391) =      16.28
 given_sigma_tot_O2_cm2( 392) =      16.30
 given_sigma_tot_O2_cm2( 393) =      16.34
 given_sigma_tot_O2_cm2( 394) =      16.45
 given_sigma_tot_O2_cm2( 395) =      16.48
 given_sigma_tot_O2_cm2( 396) =      16.51
 given_sigma_tot_O2_cm2( 397) =      16.68
 given_sigma_tot_O2_cm2( 398) =      16.70
 given_sigma_tot_O2_cm2( 399) =      16.72
 given_sigma_tot_O2_cm2( 400) =      16.80
 given_sigma_tot_O2_cm2( 401) =      16.81
 given_sigma_tot_O2_cm2( 402) =      17.00
 given_sigma_tot_O2_cm2( 403) =      17.10
 given_sigma_tot_O2_cm2( 404) =      17.12
 given_sigma_tot_O2_cm2( 405) =      17.12
 given_sigma_tot_O2_cm2( 406) =      17.18
 given_sigma_tot_O2_cm2( 407) =      17.20
 given_sigma_tot_O2_cm2( 408) =      17.20
 given_sigma_tot_O2_cm2( 409) =      17.21
 given_sigma_tot_O2_cm2( 410) =      17.25
 given_sigma_tot_O2_cm2( 411) =      17.30
 given_sigma_tot_O2_cm2( 412) =      17.40
 given_sigma_tot_O2_cm2( 413) =      17.41
 given_sigma_tot_O2_cm2( 414) =      17.50
 given_sigma_tot_O2_cm2( 415) =      17.65
 given_sigma_tot_O2_cm2( 416) =      17.65
 given_sigma_tot_O2_cm2( 417) =      17.67
 given_sigma_tot_O2_cm2( 418) =      17.72
 given_sigma_tot_O2_cm2( 419) =      17.80
 given_sigma_tot_O2_cm2( 420) =      17.80
 given_sigma_tot_O2_cm2( 421) =      17.88
 given_sigma_tot_O2_cm2( 422) =      17.92
 given_sigma_tot_O2_cm2( 423) =      18.00
 given_sigma_tot_O2_cm2( 424) =      18.03
 given_sigma_tot_O2_cm2( 425) =      18.18
 given_sigma_tot_O2_cm2( 426) =      18.19
 given_sigma_tot_O2_cm2( 427) =      18.32
 given_sigma_tot_O2_cm2( 428) =      18.40
 given_sigma_tot_O2_cm2( 429) =      18.59
 given_sigma_tot_O2_cm2( 430) =      18.80
 given_sigma_tot_O2_cm2( 431) =      19.20
 given_sigma_tot_O2_cm2( 432) =      19.59
 given_sigma_tot_O2_cm2( 433) =      19.60
 given_sigma_tot_O2_cm2( 434) =      19.65
 given_sigma_tot_O2_cm2( 435) =      19.67
 given_sigma_tot_O2_cm2( 436) =      19.68
 given_sigma_tot_O2_cm2( 437) =      19.73
 given_sigma_tot_O2_cm2( 438) =      19.80
 given_sigma_tot_O2_cm2( 439) =      19.84
 given_sigma_tot_O2_cm2( 440) =      19.88
 given_sigma_tot_O2_cm2( 441) =      19.92
 given_sigma_tot_O2_cm2( 442) =      19.96
 given_sigma_tot_O2_cm2( 443) =      20.00
 given_sigma_tot_O2_cm2( 444) =      20.03
 given_sigma_tot_O2_cm2( 445) =      20.06
 given_sigma_tot_O2_cm2( 446) =      20.09
 given_sigma_tot_O2_cm2( 447) =      20.12
 given_sigma_tot_O2_cm2( 448) =      20.15
 given_sigma_tot_O2_cm2( 449) =      20.18
 given_sigma_tot_O2_cm2( 450) =      20.21
 given_sigma_tot_O2_cm2( 451) =      20.22
 given_sigma_tot_O2_cm2( 452) =      20.23
 given_sigma_tot_O2_cm2( 453) =      20.24
 given_sigma_tot_O2_cm2( 454) =      20.27
 given_sigma_tot_O2_cm2( 455) =      20.30
 given_sigma_tot_O2_cm2( 456) =      20.34
 given_sigma_tot_O2_cm2( 457) =      20.38
 given_sigma_tot_O2_cm2( 458) =      20.42
 given_sigma_tot_O2_cm2( 459) =      20.46
 given_sigma_tot_O2_cm2( 460) =      20.50
 given_sigma_tot_O2_cm2( 461) =      20.54
 given_sigma_tot_O2_cm2( 462) =      20.58
 given_sigma_tot_O2_cm2( 463) =      20.62
 given_sigma_tot_O2_cm2( 464) =      20.66
 given_sigma_tot_O2_cm2( 465) =      20.70
 given_sigma_tot_O2_cm2( 466) =      20.71
 given_sigma_tot_O2_cm2( 467) =      20.73
 given_sigma_tot_O2_cm2( 468) =      20.76
 given_sigma_tot_O2_cm2( 469) =      20.79
 given_sigma_tot_O2_cm2( 470) =      20.82
 given_sigma_tot_O2_cm2( 471) =      20.85
 given_sigma_tot_O2_cm2( 472) =      20.88
 given_sigma_tot_O2_cm2( 473) =      20.88
 given_sigma_tot_O2_cm2( 474) =      20.90
 given_sigma_tot_O2_cm2( 475) =      20.91
 given_sigma_tot_O2_cm2( 476) =      20.94
 given_sigma_tot_O2_cm2( 477) =      20.97
 given_sigma_tot_O2_cm2( 478) =      21.00
 given_sigma_tot_O2_cm2( 479) =      21.02
 given_sigma_tot_O2_cm2( 480) =      21.04
 given_sigma_tot_O2_cm2( 481) =      21.08
 given_sigma_tot_O2_cm2( 482) =      21.09
 given_sigma_tot_O2_cm2( 483) =      21.11
 given_sigma_tot_O2_cm2( 484) =      21.12
 given_sigma_tot_O2_cm2( 485) =      21.14
 given_sigma_tot_O2_cm2( 486) =      21.16
 given_sigma_tot_O2_cm2( 487) =      21.20
 given_sigma_tot_O2_cm2( 488) =      21.20
 given_sigma_tot_O2_cm2( 489) =      21.24
 given_sigma_tot_O2_cm2( 490) =      21.25
 given_sigma_tot_O2_cm2( 491) =      21.27
 given_sigma_tot_O2_cm2( 492) =      21.28
 given_sigma_tot_O2_cm2( 493) =      21.30
 given_sigma_tot_O2_cm2( 494) =      21.32
 given_sigma_tot_O2_cm2( 495) =      21.34
 given_sigma_tot_O2_cm2( 496) =      21.36
 given_sigma_tot_O2_cm2( 497) =      21.40
 given_sigma_tot_O2_cm2( 498) =      21.41
 given_sigma_tot_O2_cm2( 499) =      21.43
 given_sigma_tot_O2_cm2( 500) =      21.46
 given_sigma_tot_O2_cm2( 501) =      21.48
 given_sigma_tot_O2_cm2( 502) =      21.49
 given_sigma_tot_O2_cm2( 503) =      21.52
 given_sigma_tot_O2_cm2( 504) =      21.53
 given_sigma_tot_O2_cm2( 505) =      21.54
 given_sigma_tot_O2_cm2( 506) =      21.55
 given_sigma_tot_O2_cm2( 507) =      21.56
 given_sigma_tot_O2_cm2( 508) =      21.58
 given_sigma_tot_O2_cm2( 509) =      21.60
 given_sigma_tot_O2_cm2( 510) =      21.61
 given_sigma_tot_O2_cm2( 511) =      21.64
 given_sigma_tot_O2_cm2( 512) =      21.67
 given_sigma_tot_O2_cm2( 513) =      21.68
 given_sigma_tot_O2_cm2( 514) =      21.70
 given_sigma_tot_O2_cm2( 515) =      21.74
 given_sigma_tot_O2_cm2( 516) =      21.78
 given_sigma_tot_O2_cm2( 517) =      21.82
 given_sigma_tot_O2_cm2( 518) =      21.86
 given_sigma_tot_O2_cm2( 519) =      21.90
 given_sigma_tot_O2_cm2( 520) =      21.91
 given_sigma_tot_O2_cm2( 521) =      21.94
 given_sigma_tot_O2_cm2( 522) =      21.98
 given_sigma_tot_O2_cm2( 523) =      22.02
 given_sigma_tot_O2_cm2( 524) =      22.06
 given_sigma_tot_O2_cm2( 525) =      22.09
 given_sigma_tot_O2_cm2( 526) =      22.10
 given_sigma_tot_O2_cm2( 527) =      22.15
 given_sigma_tot_O2_cm2( 528) =      22.20
 given_sigma_tot_O2_cm2( 529) =      22.25
 given_sigma_tot_O2_cm2( 530) =      22.26
 given_sigma_tot_O2_cm2( 531) =      22.29
 given_sigma_tot_O2_cm2( 532) =      22.30
 given_sigma_tot_O2_cm2( 533) =      22.35
 given_sigma_tot_O2_cm2( 534) =      22.40
 given_sigma_tot_O2_cm2( 535) =      22.40
 given_sigma_tot_O2_cm2( 536) =      22.42
 given_sigma_tot_O2_cm2( 537) =      22.45
 given_sigma_tot_O2_cm2( 538) =      22.50
 given_sigma_tot_O2_cm2( 539) =      22.52
 given_sigma_tot_O2_cm2( 540) =      22.55
 given_sigma_tot_O2_cm2( 541) =      22.60
 given_sigma_tot_O2_cm2( 542) =      22.63
 given_sigma_tot_O2_cm2( 543) =      22.64
 given_sigma_tot_O2_cm2( 544) =      22.67
 given_sigma_tot_O2_cm2( 545) =      22.68
 given_sigma_tot_O2_cm2( 546) =      22.68
 given_sigma_tot_O2_cm2( 547) =      22.70
 given_sigma_tot_O2_cm2( 548) =      22.72
 given_sigma_tot_O2_cm2( 549) =      22.73
 given_sigma_tot_O2_cm2( 550) =      22.76
 given_sigma_tot_O2_cm2( 551) =      22.79
 given_sigma_tot_O2_cm2( 552) =      22.80
 given_sigma_tot_O2_cm2( 553) =      22.84
 given_sigma_tot_O2_cm2( 554) =      22.86
 given_sigma_tot_O2_cm2( 555) =      22.88
 given_sigma_tot_O2_cm2( 556) =      22.92
 given_sigma_tot_O2_cm2( 557) =      22.96
 given_sigma_tot_O2_cm2( 558) =      22.98
 given_sigma_tot_O2_cm2( 559) =      23.00
 given_sigma_tot_O2_cm2( 560) =      23.06
 given_sigma_tot_O2_cm2( 561) =      23.10
 given_sigma_tot_O2_cm2( 562) =      23.12
 given_sigma_tot_O2_cm2( 563) =      23.18
 given_sigma_tot_O2_cm2( 564) =      23.24
 given_sigma_tot_O2_cm2( 565) =      23.30
 given_sigma_tot_O2_cm2( 566) =      23.36
 given_sigma_tot_O2_cm2( 567) =      23.42
 given_sigma_tot_O2_cm2( 568) =      23.48
 given_sigma_tot_O2_cm2( 569) =      23.54
 given_sigma_tot_O2_cm2( 570) =      23.56
 given_sigma_tot_O2_cm2( 571) =      23.56
 given_sigma_tot_O2_cm2( 572) =      23.60
 given_sigma_tot_O2_cm2( 573) =      23.64
 given_sigma_tot_O2_cm2( 574) =      23.64
 given_sigma_tot_O2_cm2( 575) =      23.68
 given_sigma_tot_O2_cm2( 576) =      23.72
 given_sigma_tot_O2_cm2( 577) =      23.76
 given_sigma_tot_O2_cm2( 578) =      23.80
 given_sigma_tot_O2_cm2( 579) =      23.92
 given_sigma_tot_O2_cm2( 580) =      24.00
 given_sigma_tot_O2_cm2( 581) =      24.28
 given_sigma_tot_O2_cm2( 582) =      24.50
 given_sigma_tot_O2_cm2( 583) =      24.53
 given_sigma_tot_O2_cm2( 584) =      24.54
 given_sigma_tot_O2_cm2( 585) =      24.73
 given_sigma_tot_O2_cm2( 586) =      24.90
 given_sigma_tot_O2_cm2( 587) =      25.25
 given_sigma_tot_O2_cm2( 588) =      25.40
 given_sigma_tot_O2_cm2( 589) =      25.28
 given_sigma_tot_O2_cm2( 590) =      25.20
 given_sigma_tot_O2_cm2( 591) =      25.90
 given_sigma_tot_O2_cm2( 592) =      25.90
 given_sigma_tot_O2_cm2( 593) =      25.90
 given_sigma_tot_O2_cm2( 594) =      26.04
 given_sigma_tot_O2_cm2( 595) =      26.05
 given_sigma_tot_O2_cm2( 596) =      26.10
 given_sigma_tot_O2_cm2( 597) =      26.10
 given_sigma_tot_O2_cm2( 598) =      26.10
 given_sigma_tot_O2_cm2( 599) =      26.04
 given_sigma_tot_O2_cm2( 600) =      26.00
 given_sigma_tot_O2_cm2( 601) =      25.86
 given_sigma_tot_O2_cm2( 602) =      25.80
 given_sigma_tot_O2_cm2( 603) =      25.66
 given_sigma_tot_O2_cm2( 604) =      25.50
 given_sigma_tot_O2_cm2( 605) =      24.80
 given_sigma_tot_O2_cm2( 606) =      24.61
 given_sigma_tot_O2_cm2( 607) =      22.88
 given_sigma_tot_O2_cm2( 608) =      22.72
 given_sigma_tot_O2_cm2( 609) =      22.40
 given_sigma_tot_O2_cm2( 610) =      21.80
 given_sigma_tot_O2_cm2( 611) =      19.30
 given_sigma_tot_O2_cm2( 612) =      19.40
 given_sigma_tot_O2_cm2( 613) =      20.80
 given_sigma_tot_O2_cm2( 614) =      21.80
 given_sigma_tot_O2_cm2( 615) =      23.80
 given_sigma_tot_O2_cm2( 616) =      27.70
 given_sigma_tot_O2_cm2( 617) =      28.49
 given_sigma_tot_O2_cm2( 618) =      28.60
 given_sigma_tot_O2_cm2( 619) =      30.10
 given_sigma_tot_O2_cm2( 620) =      27.12
 given_sigma_tot_O2_cm2( 621) =      26.80
 given_sigma_tot_O2_cm2( 622) =      26.40
 given_sigma_tot_O2_cm2( 623) =      27.96
 given_sigma_tot_O2_cm2( 624) =      28.54
 given_sigma_tot_O2_cm2( 625) =      29.91
 given_sigma_tot_O2_cm2( 626) =      29.73
 given_sigma_tot_O2_cm2( 627) =      25.68
 given_sigma_tot_O2_cm2( 628) =      24.21
 given_sigma_tot_O2_cm2( 629) =      23.10
 given_sigma_tot_O2_cm2( 630) =      27.00
 given_sigma_tot_O2_cm2( 631) =      28.43
 given_sigma_tot_O2_cm2( 632) =      28.57
 given_sigma_tot_O2_cm2( 633) =      28.60
 given_sigma_tot_O2_cm2( 634) =      25.40
 given_sigma_tot_O2_cm2( 635) =      24.20
 given_sigma_tot_O2_cm2( 636) =      26.43
 given_sigma_tot_O2_cm2( 637) =      28.66
 given_sigma_tot_O2_cm2( 638) =      34.60
 given_sigma_tot_O2_cm2( 639) =      30.82
 given_sigma_tot_O2_cm2( 640) =      28.30
 given_sigma_tot_O2_cm2( 641) =      29.80
 given_sigma_tot_O2_cm2( 642) =      29.06
 given_sigma_tot_O2_cm2( 643) =      25.71
 given_sigma_tot_O2_cm2( 644) =      24.96
 given_sigma_tot_O2_cm2( 645) =      23.10
 given_sigma_tot_O2_cm2( 646) =      26.53
 given_sigma_tot_O2_cm2( 647) =      32.70
 given_sigma_tot_O2_cm2( 648) =      29.57
 given_sigma_tot_O2_cm2( 649) =      26.72
 given_sigma_tot_O2_cm2( 650) =      25.30
 given_sigma_tot_O2_cm2( 651) =      26.06
 given_sigma_tot_O2_cm2( 652) =      27.77
 given_sigma_tot_O2_cm2( 653) =      29.00
 given_sigma_tot_O2_cm2( 654) =      25.72
 given_sigma_tot_O2_cm2( 655) =      27.60
 given_sigma_tot_O2_cm2( 656) =      30.60
 given_sigma_tot_O2_cm2( 657) =      32.40
 given_sigma_tot_O2_cm2( 658) =      32.14
 given_sigma_tot_O2_cm2( 659) =      31.34
 given_sigma_tot_O2_cm2( 660) =      30.55
 given_sigma_tot_O2_cm2( 661) =      28.69
 given_sigma_tot_O2_cm2( 662) =      27.11
 given_sigma_tot_O2_cm2( 663) =      26.05
 given_sigma_tot_O2_cm2( 664) =      23.40
 given_sigma_tot_O2_cm2( 665) =      29.11
 given_sigma_tot_O2_cm2( 666) =      31.40
 given_sigma_tot_O2_cm2( 667) =      30.09
 given_sigma_tot_O2_cm2( 668) =      28.79
 given_sigma_tot_O2_cm2( 669) =      25.30
 given_sigma_tot_O2_cm2( 670) =      26.34
 given_sigma_tot_O2_cm2( 671) =      27.90
 given_sigma_tot_O2_cm2( 672) =      26.08
 given_sigma_tot_O2_cm2( 673) =      25.30
 given_sigma_tot_O2_cm2( 674) =      27.58
 given_sigma_tot_O2_cm2( 675) =      26.17
 given_sigma_tot_O2_cm2( 676) =      25.30
 given_sigma_tot_O2_cm2( 677) =      25.59
 given_sigma_tot_O2_cm2( 678) =      26.57
 given_sigma_tot_O2_cm2( 679) =      26.96
 given_sigma_tot_O2_cm2( 680) =      27.47
 given_sigma_tot_O2_cm2( 681) =      28.33
 given_sigma_tot_O2_cm2( 682) =      28.51
 given_sigma_tot_O2_cm2( 683) =      29.00
 given_sigma_tot_O2_cm2( 684) =      26.59
 given_sigma_tot_O2_cm2( 685) =      21.56
 given_sigma_tot_O2_cm2( 686) =      20.80
 given_sigma_tot_O2_cm2( 687) =      21.93
 given_sigma_tot_O2_cm2( 688) =      25.30
 given_sigma_tot_O2_cm2( 689) =      24.90
 given_sigma_tot_O2_cm2( 690) =      29.45
 given_sigma_tot_O2_cm2( 691) =      30.10
 given_sigma_tot_O2_cm2( 692) =      29.57
 given_sigma_tot_O2_cm2( 693) =      28.68
 given_sigma_tot_O2_cm2( 694) =      27.79
 given_sigma_tot_O2_cm2( 695) =      26.01
 given_sigma_tot_O2_cm2( 696) =      25.30
 given_sigma_tot_O2_cm2( 697) =      26.42
 given_sigma_tot_O2_cm2( 698) =      26.99
 given_sigma_tot_O2_cm2( 699) =      28.30
 given_sigma_tot_O2_cm2( 700) =      29.80
 given_sigma_tot_O2_cm2( 701) =      29.71
 given_sigma_tot_O2_cm2( 702) =      29.29
 given_sigma_tot_O2_cm2( 703) =      28.86
 given_sigma_tot_O2_cm2( 704) =      28.43
 given_sigma_tot_O2_cm2( 705) =      28.00
 given_sigma_tot_O2_cm2( 706) =      27.57
 given_sigma_tot_O2_cm2( 707) =      27.44
 given_sigma_tot_O2_cm2( 708) =      27.14
 given_sigma_tot_O2_cm2( 709) =      26.71
 given_sigma_tot_O2_cm2( 710) =      26.29
 given_sigma_tot_O2_cm2( 711) =      26.16
 given_sigma_tot_O2_cm2( 712) =      25.86
 given_sigma_tot_O2_cm2( 713) =      25.69
 given_sigma_tot_O2_cm2( 714) =      25.45
 given_sigma_tot_O2_cm2( 715) =      25.00
 given_sigma_tot_O2_cm2( 716) =      24.57
 given_sigma_tot_O2_cm2( 717) =      24.31
 given_sigma_tot_O2_cm2( 718) =      24.19
 given_sigma_tot_O2_cm2( 719) =      24.01
 given_sigma_tot_O2_cm2( 720) =      23.80
 given_sigma_tot_O2_cm2( 721) =      23.74
 given_sigma_tot_O2_cm2( 722) =      23.45
 given_sigma_tot_O2_cm2( 723) =      23.37
 given_sigma_tot_O2_cm2( 724) =      23.16
 given_sigma_tot_O2_cm2( 725) =      22.87
 given_sigma_tot_O2_cm2( 726) =      22.70
 given_sigma_tot_O2_cm2( 727) =      20.80
 given_sigma_tot_O2_cm2( 728) =      21.44
 given_sigma_tot_O2_cm2( 729) =      21.23
 given_sigma_tot_O2_cm2( 730) =      20.86
 given_sigma_tot_O2_cm2( 731) =      20.56
 given_sigma_tot_O2_cm2( 732) =      19.82
 given_sigma_tot_O2_cm2( 733) =      19.30
 given_sigma_tot_O2_cm2( 734) =      19.63
 given_sigma_tot_O2_cm2( 735) =      19.96
 given_sigma_tot_O2_cm2( 736) =      20.61
 given_sigma_tot_O2_cm2( 737) =      21.60
 given_sigma_tot_O2_cm2( 738) =      21.17
 given_sigma_tot_O2_cm2( 739) =      20.10
 given_sigma_tot_O2_cm2( 740) =      21.00
 given_sigma_tot_O2_cm2( 741) =      21.60
 given_sigma_tot_O2_cm2( 742) =      20.10
 given_sigma_tot_O2_cm2( 743) =      21.20
 given_sigma_tot_O2_cm2( 744) =      20.52
 given_sigma_tot_O2_cm2( 745) =      19.97
 given_sigma_tot_O2_cm2( 746) =      19.43
 given_sigma_tot_O2_cm2( 747) =      18.75
 given_sigma_tot_O2_cm2( 748) =      18.47
 given_sigma_tot_O2_cm2( 749) =      18.20
 given_sigma_tot_O2_cm2( 750) =      23.80
 given_sigma_tot_O2_cm2( 751) =      22.11
 given_sigma_tot_O2_cm2( 752) =      20.99
 given_sigma_tot_O2_cm2( 753) =      19.30
 given_sigma_tot_O2_cm2( 754) =      21.55
 given_sigma_tot_O2_cm2( 755) =      22.03
 given_sigma_tot_O2_cm2( 756) =      21.90
 given_sigma_tot_O2_cm2( 757) =      22.30
 given_sigma_tot_O2_cm2( 758) =      21.77
 given_sigma_tot_O2_cm2( 759) =      20.19
 given_sigma_tot_O2_cm2( 760) =      18.60
 given_sigma_tot_O2_cm2( 761) =      21.90
 given_sigma_tot_O2_cm2( 762) =      21.02
 given_sigma_tot_O2_cm2( 763) =      19.70
 given_sigma_tot_O2_cm2( 764) =      17.50
 given_sigma_tot_O2_cm2( 765) =      18.38
 given_sigma_tot_O2_cm2( 766) =      21.91
 given_sigma_tot_O2_cm2( 767) =      23.67
 given_sigma_tot_O2_cm2( 768) =      26.32
 given_sigma_tot_O2_cm2( 769) =      27.20
 given_sigma_tot_O2_cm2( 770) =      21.95
 given_sigma_tot_O2_cm2( 771) =      20.11
 given_sigma_tot_O2_cm2( 772) =      17.57
 given_sigma_tot_O2_cm2( 773) =      19.54
 given_sigma_tot_O2_cm2( 774) =      23.80
 given_sigma_tot_O2_cm2( 775) =      23.25
 given_sigma_tot_O2_cm2( 776) =      21.62
 given_sigma_tot_O2_cm2( 777) =      19.98
 given_sigma_tot_O2_cm2( 778) =      16.70
 given_sigma_tot_O2_cm2( 779) =      19.26
 given_sigma_tot_O2_cm2( 780) =      20.80
 given_sigma_tot_O2_cm2( 781) =      17.50
 given_sigma_tot_O2_cm2( 782) =      25.04
 given_sigma_tot_O2_cm2( 783) =      22.63
 given_sigma_tot_O2_cm2( 784) =      21.27
 given_sigma_tot_O2_cm2( 785) =      19.92
 given_sigma_tot_O2_cm2( 786) =      18.57
 given_sigma_tot_O2_cm2( 787) =      19.42
 given_sigma_tot_O2_cm2( 788) =      28.53
 given_sigma_tot_O2_cm2( 789) =      31.56
 given_sigma_tot_O2_cm2( 790) =      34.60
 given_sigma_tot_O2_cm2( 791) =      31.53
 given_sigma_tot_O2_cm2( 792) =      28.45
 given_sigma_tot_O2_cm2( 793) =      20.25
 given_sigma_tot_O2_cm2( 794) =      18.20
 given_sigma_tot_O2_cm2( 795) =      21.27
 given_sigma_tot_O2_cm2( 796) =      27.43
 given_sigma_tot_O2_cm2( 797) =      30.50
 given_sigma_tot_O2_cm2( 798) =      19.30
 given_sigma_tot_O2_cm2( 799) =      23.61
 given_sigma_tot_O2_cm2( 800) =      27.92
 given_sigma_tot_O2_cm2( 801) =      30.50
 given_sigma_tot_O2_cm2( 802) =      31.05
 given_sigma_tot_O2_cm2( 803) =      31.60
 given_sigma_tot_O2_cm2( 804) =      33.10
 given_sigma_tot_O2_cm2( 805) =      34.60
 given_sigma_tot_O2_cm2( 806) =      29.43
 given_sigma_tot_O2_cm2( 807) =      25.12
 given_sigma_tot_O2_cm2( 808) =      23.40
 given_sigma_tot_O2_cm2( 809) =      24.37
 given_sigma_tot_O2_cm2( 810) =      25.27
 given_sigma_tot_O2_cm2( 811) =      30.17
 given_sigma_tot_O2_cm2( 812) =      31.13
 given_sigma_tot_O2_cm2( 813) =      35.00
 given_sigma_tot_O2_cm2( 814) =      31.00
 given_sigma_tot_O2_cm2( 815) =      27.00
 given_sigma_tot_O2_cm2( 816) =      19.00
 given_sigma_tot_O2_cm2( 817) =      20.33
 given_sigma_tot_O2_cm2( 818) =      21.67
 given_sigma_tot_O2_cm2( 819) =      28.33
 given_sigma_tot_O2_cm2( 820) =      27.44
 given_sigma_tot_O2_cm2( 821) =      26.00
 given_sigma_tot_O2_cm2( 822) =      34.75
 given_sigma_tot_O2_cm2( 823) =      45.69
 given_sigma_tot_O2_cm2( 824) =      56.64
 given_sigma_tot_O2_cm2( 825) =      63.20
 given_sigma_tot_O2_cm2( 826) =      45.52
 given_sigma_tot_O2_cm2( 827) =      24.90
 given_sigma_tot_O2_cm2( 828) =      25.36
 given_sigma_tot_O2_cm2( 829) =      26.05
 given_sigma_tot_O2_cm2( 830) =      26.40
 given_sigma_tot_O2_cm2( 831) =      25.70
 given_sigma_tot_O2_cm2( 832) =      26.53
 given_sigma_tot_O2_cm2( 833) =      25.70
 given_sigma_tot_O2_cm2( 834) =      29.39
 given_sigma_tot_O2_cm2( 835) =      39.24
 given_sigma_tot_O2_cm2( 836) =      41.70
 given_sigma_tot_O2_cm2( 837) =      32.40
 given_sigma_tot_O2_cm2( 838) =      35.70
 given_sigma_tot_O2_cm2( 839) =      36.80
 given_sigma_tot_O2_cm2( 840) =      37.90
 given_sigma_tot_O2_cm2( 841) =      35.28
 given_sigma_tot_O2_cm2( 842) =      33.54
 given_sigma_tot_O2_cm2( 843) =      35.00
 given_sigma_tot_O2_cm2( 844) =      33.62
 given_sigma_tot_O2_cm2( 845) =      30.85
 given_sigma_tot_O2_cm2( 846) =      30.38
 given_sigma_tot_O2_cm2( 847) =      29.00
 given_sigma_tot_O2_cm2( 848) =      30.50
 given_sigma_tot_O2_cm2( 849) =      28.36
 given_sigma_tot_O2_cm2( 850) =      27.50
 given_sigma_tot_O2_cm2( 851) =      29.80
 given_sigma_tot_O2_cm2( 852) =      28.71
 given_sigma_tot_O2_cm2( 853) =      27.35
 given_sigma_tot_O2_cm2( 854) =      25.99
 given_sigma_tot_O2_cm2( 855) =      25.44
 given_sigma_tot_O2_cm2( 856) =      24.90
 given_sigma_tot_O2_cm2( 857) =      30.48
 given_sigma_tot_O2_cm2( 858) =      34.20
 given_sigma_tot_O2_cm2( 859) =      30.09
 given_sigma_tot_O2_cm2( 860) =      26.80
 given_sigma_tot_O2_cm2( 861) =      26.94
 given_sigma_tot_O2_cm2( 862) =      27.50
 given_sigma_tot_O2_cm2( 863) =      27.05
 given_sigma_tot_O2_cm2( 864) =      25.92
 given_sigma_tot_O2_cm2( 865) =      26.00
 given_sigma_tot_O2_cm2( 866) =      25.30
 given_sigma_tot_O2_cm2( 867) =      26.68
 given_sigma_tot_O2_cm2( 868) =      27.50
 given_sigma_tot_O2_cm2( 869) =      25.85
 given_sigma_tot_O2_cm2( 870) =      25.57
 given_sigma_tot_O2_cm2( 871) =      33.80
 given_sigma_tot_O2_cm2( 872) =      37.20
 given_sigma_tot_O2_cm2( 873) =      42.30
 given_sigma_tot_O2_cm2( 874) =      49.10
 given_sigma_tot_O2_cm2( 875) =      36.43
 given_sigma_tot_O2_cm2( 876) =      30.10
 given_sigma_tot_O2_cm2( 877) =      30.18
 given_sigma_tot_O2_cm2( 878) =      30.38
 given_sigma_tot_O2_cm2( 879) =      30.50
 given_sigma_tot_O2_cm2( 880) =      29.01
 given_sigma_tot_O2_cm2( 881) =      25.90
 given_sigma_tot_O2_cm2( 882) =      29.80
 given_sigma_tot_O2_cm2( 883) =      29.00
 given_sigma_tot_O2_cm2( 884) =      29.97
 given_sigma_tot_O2_cm2( 885) =      32.88
 given_sigma_tot_O2_cm2( 886) =      34.33
 given_sigma_tot_O2_cm2( 887) =      33.19
 given_sigma_tot_O2_cm2( 888) =      31.60
 given_sigma_tot_O2_cm2( 889) =      37.23
 given_sigma_tot_O2_cm2( 890) =      42.86
 given_sigma_tot_O2_cm2( 891) =      51.30
 given_sigma_tot_O2_cm2( 892) =      39.67
 given_sigma_tot_O2_cm2( 893) =      32.70
 given_sigma_tot_O2_cm2( 894) =      33.61
 given_sigma_tot_O2_cm2( 895) =      34.26
 given_sigma_tot_O2_cm2( 896) =      34.91
 given_sigma_tot_O2_cm2( 897) =      35.30
 given_sigma_tot_O2_cm2( 898) =      34.38
 given_sigma_tot_O2_cm2( 899) =      33.47
 given_sigma_tot_O2_cm2( 900) =      32.71
 given_sigma_tot_O2_cm2( 901) =      32.40
 given_sigma_tot_O2_cm2( 902) =      34.20
 given_sigma_tot_O2_cm2( 903) =      31.62
 given_sigma_tot_O2_cm2( 904) =      29.55
 given_sigma_tot_O2_cm2( 905) =      26.45
 given_sigma_tot_O2_cm2( 906) =      25.42
 given_sigma_tot_O2_cm2( 907) =      25.70
 given_sigma_tot_O2_cm2( 908) =      24.77
 given_sigma_tot_O2_cm2( 909) =      21.03
 given_sigma_tot_O2_cm2( 910) =      20.10
 given_sigma_tot_O2_cm2( 911) =      21.30
 given_sigma_tot_O2_cm2( 912) =      21.60
 given_sigma_tot_O2_cm2( 913) =      20.86
 given_sigma_tot_O2_cm2( 914) =      18.64
 given_sigma_tot_O2_cm2( 915) =      17.90
 given_sigma_tot_O2_cm2( 916) =      18.33
 given_sigma_tot_O2_cm2( 917) =      18.62
 given_sigma_tot_O2_cm2( 918) =      19.06
 given_sigma_tot_O2_cm2( 919) =      19.78
 given_sigma_tot_O2_cm2( 920) =      20.36
 given_sigma_tot_O2_cm2( 921) =      19.14
 given_sigma_tot_O2_cm2( 922) =      18.60
 given_sigma_tot_O2_cm2( 923) =      19.90
 given_sigma_tot_O2_cm2( 924) =      21.20
 given_sigma_tot_O2_cm2( 925) =      18.40
 given_sigma_tot_O2_cm2( 926) =      15.60
 given_sigma_tot_O2_cm2( 927) =      17.65
 given_sigma_tot_O2_cm2( 928) =      19.70
 given_sigma_tot_O2_cm2( 929) =      20.52
 given_sigma_tot_O2_cm2( 930) =      23.80
 given_sigma_tot_O2_cm2( 931) =      19.39
 given_sigma_tot_O2_cm2( 932) =      17.50
 given_sigma_tot_O2_cm2( 933) =      19.00
 given_sigma_tot_O2_cm2( 934) =      17.74
 given_sigma_tot_O2_cm2( 935) =      16.79
 given_sigma_tot_O2_cm2( 936) =      14.90
 given_sigma_tot_O2_cm2( 937) =      19.35
 given_sigma_tot_O2_cm2( 938) =      20.97
 given_sigma_tot_O2_cm2( 939) =      23.40
 given_sigma_tot_O2_cm2( 940) =      22.52
 given_sigma_tot_O2_cm2( 941) =      21.20
 given_sigma_tot_O2_cm2( 942) =      19.88
 given_sigma_tot_O2_cm2( 943) =      19.00
 given_sigma_tot_O2_cm2( 944) =      19.30
 given_sigma_tot_O2_cm2( 945) =      19.00
 given_sigma_tot_O2_cm2( 946) =      18.50
 given_sigma_tot_O2_cm2( 947) =      18.30
 given_sigma_tot_O2_cm2( 948) =      18.00
 given_sigma_tot_O2_cm2( 949) =      17.80
 given_sigma_tot_O2_cm2( 950) =      17.50
 given_sigma_tot_O2_cm2( 951) =      18.85
 given_sigma_tot_O2_cm2( 952) =      19.30
 given_sigma_tot_O2_cm2( 953) =      18.68
 given_sigma_tot_O2_cm2( 954) =      17.98
 given_sigma_tot_O2_cm2( 955) =      17.10
 given_sigma_tot_O2_cm2( 956) =      18.99
 given_sigma_tot_O2_cm2( 957) =      19.62
 given_sigma_tot_O2_cm2( 958) =      20.25
 given_sigma_tot_O2_cm2( 959) =      21.20
 given_sigma_tot_O2_cm2( 960) =      20.49
 given_sigma_tot_O2_cm2( 961) =      19.30
 given_sigma_tot_O2_cm2( 962) =      19.63
 given_sigma_tot_O2_cm2( 963) =      19.60
 given_sigma_tot_O2_cm2( 964) =      19.30
 given_sigma_tot_O2_cm2( 965) =      21.16
 given_sigma_tot_O2_cm2( 966) =      21.90
 given_sigma_tot_O2_cm2( 967) =      20.96
 given_sigma_tot_O2_cm2( 968) =      20.33
 given_sigma_tot_O2_cm2( 969) =      19.39
 given_sigma_tot_O2_cm2( 970) =      18.13
 given_sigma_tot_O2_cm2( 971) =      18.10
 given_sigma_tot_O2_cm2( 972) =      19.90
 given_sigma_tot_O2_cm2( 973) =      21.70
 given_sigma_tot_O2_cm2( 974) =      22.30
 given_sigma_tot_O2_cm2( 975) =      21.54
 given_sigma_tot_O2_cm2( 976) =      20.78
 given_sigma_tot_O2_cm2( 977) =      19.51
 given_sigma_tot_O2_cm2( 978) =      19.00
 given_sigma_tot_O2_cm2( 979) =      19.55
 given_sigma_tot_O2_cm2( 980) =      20.10
 given_sigma_tot_O2_cm2( 981) =      19.55
 given_sigma_tot_O2_cm2( 982) =      19.01
 given_sigma_tot_O2_cm2( 983) =      18.46
 given_sigma_tot_O2_cm2( 984) =      17.37
 given_sigma_tot_O2_cm2( 985) =      17.10
 given_sigma_tot_O2_cm2( 986) =      19.88
 given_sigma_tot_O2_cm2( 987) =      19.35
 given_sigma_tot_O2_cm2( 988) =      17.90
 given_sigma_tot_O2_cm2( 989) =      20.38
 given_sigma_tot_O2_cm2( 990) =      17.20
 given_sigma_tot_O2_cm2( 991) =      15.20
 given_sigma_tot_O2_cm2( 992) =      18.73
 given_sigma_tot_O2_cm2( 993) =      19.70
 given_sigma_tot_O2_cm2( 994) =      20.83
 given_sigma_tot_O2_cm2( 995) =      23.64
 given_sigma_tot_O2_cm2( 996) =      24.00
 given_sigma_tot_O2_cm2( 997) =      23.80
 given_sigma_tot_O2_cm2( 998) =      26.37
 given_sigma_tot_O2_cm2( 999) =      24.20
 given_sigma_tot_O2_cm2(1000) =      20.95
 given_sigma_tot_O2_cm2(1001) =      17.70
 given_sigma_tot_O2_cm2(1002) =      14.45
 given_sigma_tot_O2_cm2(1003) =      15.87
 given_sigma_tot_O2_cm2(1004) =      16.90
 given_sigma_tot_O2_cm2(1005) =      18.45
 given_sigma_tot_O2_cm2(1006) =      21.03
 given_sigma_tot_O2_cm2(1007) =      23.62
 given_sigma_tot_O2_cm2(1008) =      26.20
 given_sigma_tot_O2_cm2(1009) =      28.78
 given_sigma_tot_O2_cm2(1010) =      26.56
 given_sigma_tot_O2_cm2(1011) =      25.19
 given_sigma_tot_O2_cm2(1012) =      25.07
 given_sigma_tot_O2_cm2(1013) =      26.48
 given_sigma_tot_O2_cm2(1014) =      27.33
 given_sigma_tot_O2_cm2(1015) =      27.62
 given_sigma_tot_O2_cm2(1016) =      25.22
 given_sigma_tot_O2_cm2(1017) =      23.43
 given_sigma_tot_O2_cm2(1018) =      18.97
 given_sigma_tot_O2_cm2(1019) =      17.18
 given_sigma_tot_O2_cm2(1020) =      14.50
 given_sigma_tot_O2_cm2(1021) =      16.75
 given_sigma_tot_O2_cm2(1022) =      19.00
 given_sigma_tot_O2_cm2(1023) =      20.80
 given_sigma_tot_O2_cm2(1024) =      20.50
 given_sigma_tot_O2_cm2(1025) =      21.80
 given_sigma_tot_O2_cm2(1026) =      23.10
 given_sigma_tot_O2_cm2(1027) =      23.97
 given_sigma_tot_O2_cm2(1028) =      25.70
 given_sigma_tot_O2_cm2(1029) =      23.80
 given_sigma_tot_O2_cm2(1030) =      24.30
 given_sigma_tot_O2_cm2(1031) =      25.55
 given_sigma_tot_O2_cm2(1032) =      26.80
 given_sigma_tot_O2_cm2(1033) =      23.25
 given_sigma_tot_O2_cm2(1034) =      19.70
 given_sigma_tot_O2_cm2(1035) =      22.70
 given_sigma_tot_O2_cm2(1036) =      25.20
 given_sigma_tot_O2_cm2(1037) =      25.32
 given_sigma_tot_O2_cm2(1038) =      24.20
 given_sigma_tot_O2_cm2(1039) =      26.95
 given_sigma_tot_O2_cm2(1040) =      26.80
 given_sigma_tot_O2_cm2(1041) =      27.55
 given_sigma_tot_O2_cm2(1042) =      28.30
 given_sigma_tot_O2_cm2(1043) =      25.57
 given_sigma_tot_O2_cm2(1044) =      23.93
 given_sigma_tot_O2_cm2(1045) =      22.84
 given_sigma_tot_O2_cm2(1046) =      21.20
 given_sigma_tot_O2_cm2(1047) =      21.81
 given_sigma_tot_O2_cm2(1048) =      24.25
 given_sigma_tot_O2_cm2(1049) =      27.90
 given_sigma_tot_O2_cm2(1050) =      24.20
 given_sigma_tot_O2_cm2(1051) =      25.09
 given_sigma_tot_O2_cm2(1052) =      26.71
 given_sigma_tot_O2_cm2(1053) =      29.37
 given_sigma_tot_O2_cm2(1054) =      33.06
 given_sigma_tot_O2_cm2(1055) =      29.04
 given_sigma_tot_O2_cm2(1056) =      23.10
 given_sigma_tot_O2_cm2(1057) =      23.40
 given_sigma_tot_O2_cm2(1058) =      21.90
 given_sigma_tot_O2_cm2(1059) =      27.61
 given_sigma_tot_O2_cm2(1060) =      31.60
 given_sigma_tot_O2_cm2(1061) =      28.00
 given_sigma_tot_O2_cm2(1062) =      35.16
 given_sigma_tot_O2_cm2(1063) =      39.80
 given_sigma_tot_O2_cm2(1064) =      36.80
 given_sigma_tot_O2_cm2(1065) =      33.80
 given_sigma_tot_O2_cm2(1066) =      30.80
 given_sigma_tot_O2_cm2(1067) =      27.80
 given_sigma_tot_O2_cm2(1068) =      30.00
 given_sigma_tot_O2_cm2(1069) =      34.20
 given_sigma_tot_O2_cm2(1070) =      30.56
 given_sigma_tot_O2_cm2(1071) =      26.00
 given_sigma_tot_O2_cm2(1072) =      33.22
 given_sigma_tot_O2_cm2(1073) =      37.07
 given_sigma_tot_O2_cm2(1074) =      38.71
 given_sigma_tot_O2_cm2(1075) =      40.44
 given_sigma_tot_O2_cm2(1076) =      44.41
 given_sigma_tot_O2_cm2(1077) =      47.66
 given_sigma_tot_O2_cm2(1078) =      46.92
 given_sigma_tot_O2_cm2(1079) =      45.31
 given_sigma_tot_O2_cm2(1080) =      41.90
 given_sigma_tot_O2_cm2(1081) =      39.06
 given_sigma_tot_O2_cm2(1082) =      36.48
 given_sigma_tot_O2_cm2(1083) =      34.37
 given_sigma_tot_O2_cm2(1084) =      27.90
 given_sigma_tot_O2_cm2(1085) =      34.98
 given_sigma_tot_O2_cm2(1086) =      36.40
 given_sigma_tot_O2_cm2(1087) =      32.55
 given_sigma_tot_O2_cm2(1088) =      28.71
 given_sigma_tot_O2_cm2(1089) =      26.15
 given_sigma_tot_O2_cm2(1090) =      22.30
 given_sigma_tot_O2_cm2(1091) =      24.92
 given_sigma_tot_O2_cm2(1092) =      31.46
 given_sigma_tot_O2_cm2(1093) =      38.00
 given_sigma_tot_O2_cm2(1094) =      40.15
 given_sigma_tot_O2_cm2(1095) =      42.57
 given_sigma_tot_O2_cm2(1096) =      44.54
 given_sigma_tot_O2_cm2(1097) =      47.94
 given_sigma_tot_O2_cm2(1098) =      50.92
 given_sigma_tot_O2_cm2(1099) =      52.51
 given_sigma_tot_O2_cm2(1100) =      55.00
 given_sigma_tot_O2_cm2(1101) =      47.69
 given_sigma_tot_O2_cm2(1102) =      37.81
 given_sigma_tot_O2_cm2(1103) =      29.40
 given_sigma_tot_O2_cm2(1104) =      30.80
 given_sigma_tot_O2_cm2(1105) =      31.73
 given_sigma_tot_O2_cm2(1106) =      34.07
 given_sigma_tot_O2_cm2(1107) =      35.00
 given_sigma_tot_O2_cm2(1108) =      30.90
 given_sigma_tot_O2_cm2(1109) =      24.07
 given_sigma_tot_O2_cm2(1110) =      18.60
 given_sigma_tot_O2_cm2(1111) =      25.59
 given_sigma_tot_O2_cm2(1112) =      31.42
 given_sigma_tot_O2_cm2(1113) =      36.31
 given_sigma_tot_O2_cm2(1114) =      40.34
 given_sigma_tot_O2_cm2(1115) =      43.07
 given_sigma_tot_O2_cm2(1116) =      45.28
 given_sigma_tot_O2_cm2(1117) =      37.81
 given_sigma_tot_O2_cm2(1118) =      30.80
 given_sigma_tot_O2_cm2(1119) =      25.16
 given_sigma_tot_O2_cm2(1120) =      20.10
 given_sigma_tot_O2_cm2(1121) =      21.11
 given_sigma_tot_O2_cm2(1122) =      22.26
 given_sigma_tot_O2_cm2(1123) =      25.85
 given_sigma_tot_O2_cm2(1124) =      29.44
 given_sigma_tot_O2_cm2(1125) =      31.60
 given_sigma_tot_O2_cm2(1126) =      29.52
 given_sigma_tot_O2_cm2(1127) =      24.32
 given_sigma_tot_O2_cm2(1128) =      19.12
 given_sigma_tot_O2_cm2(1129) =      16.00
 given_sigma_tot_O2_cm2(1130) =      17.29
 given_sigma_tot_O2_cm2(1131) =      20.53
 given_sigma_tot_O2_cm2(1132) =      27.01
 given_sigma_tot_O2_cm2(1133) =      28.30
 given_sigma_tot_O2_cm2(1134) =      17.37
 given_sigma_tot_O2_cm2(1135) =      20.10
 given_sigma_tot_O2_cm2(1136) =      24.20
 given_sigma_tot_O2_cm2(1137) =      20.10
 given_sigma_tot_O2_cm2(1138) =      22.64
 given_sigma_tot_O2_cm2(1139) =      29.00
 given_sigma_tot_O2_cm2(1140) =      18.31
 given_sigma_tot_O2_cm2(1141) =      11.90
 given_sigma_tot_O2_cm2(1142) =      11.98
 given_sigma_tot_O2_cm2(1143) =      12.18
 given_sigma_tot_O2_cm2(1144) =      12.30
 given_sigma_tot_O2_cm2(1145) =      11.86
 given_sigma_tot_O2_cm2(1146) =      11.20
 given_sigma_tot_O2_cm2(1147) =      18.52
 given_sigma_tot_O2_cm2(1148) =      22.70
 given_sigma_tot_O2_cm2(1149) =      22.30
 given_sigma_tot_O2_cm2(1150) =      23.10
 given_sigma_tot_O2_cm2(1151) =      21.05
 given_sigma_tot_O2_cm2(1152) =      10.80
 given_sigma_tot_O2_cm2(1153) =      25.40
 given_sigma_tot_O2_cm2(1154) =      32.70
 given_sigma_tot_O2_cm2(1155) =      28.52
 given_sigma_tot_O2_cm2(1156) =      25.12
 given_sigma_tot_O2_cm2(1157) =      18.56
 given_sigma_tot_O2_cm2(1158) =      15.97
 given_sigma_tot_O2_cm2(1159) =      11.59
 given_sigma_tot_O2_cm2(1160) =      10.50
 given_sigma_tot_O2_cm2(1161) =      10.80
 given_sigma_tot_O2_cm2(1162) =      10.36
 given_sigma_tot_O2_cm2(1163) =      10.18
 given_sigma_tot_O2_cm2(1164) =      10.00
 given_sigma_tot_O2_cm2(1165) =      14.73
 given_sigma_tot_O2_cm2(1166) =      17.10
 given_sigma_tot_O2_cm2(1167) =      14.35
 given_sigma_tot_O2_cm2(1168) =      12.38
 given_sigma_tot_O2_cm2(1169) =      11.20
 given_sigma_tot_O2_cm2(1170) =      14.53
 given_sigma_tot_O2_cm2(1171) =      24.50
 given_sigma_tot_O2_cm2(1172) =      23.40
 given_sigma_tot_O2_cm2(1173) =      19.89
 given_sigma_tot_O2_cm2(1174) =      17.11
 given_sigma_tot_O2_cm2(1175) =      15.99
 given_sigma_tot_O2_cm2(1176) =      14.32
 given_sigma_tot_O2_cm2(1177) =      11.54
 given_sigma_tot_O2_cm2(1178) =       8.76
 given_sigma_tot_O2_cm2(1179) =       9.16
 given_sigma_tot_O2_cm2(1180) =      10.37
 given_sigma_tot_O2_cm2(1181) =      11.58
 given_sigma_tot_O2_cm2(1182) =      12.30
 given_sigma_tot_O2_cm2(1183) =      11.72
 given_sigma_tot_O2_cm2(1184) =      10.29
 given_sigma_tot_O2_cm2(1185) =      12.65
 given_sigma_tot_O2_cm2(1186) =      15.95
 given_sigma_tot_O2_cm2(1187) =      18.60
 given_sigma_tot_O2_cm2(1188) =      11.35
 given_sigma_tot_O2_cm2(1189) =       7.40
 given_sigma_tot_O2_cm2(1190) =       7.58
 given_sigma_tot_O2_cm2(1191) =       7.80
 given_sigma_tot_O2_cm2(1192) =       7.51
 given_sigma_tot_O2_cm2(1193) =       7.40
 given_sigma_tot_O2_cm2(1194) =       7.89
 given_sigma_tot_O2_cm2(1195) =       8.71
 given_sigma_tot_O2_cm2(1196) =       9.70
 given_sigma_tot_O2_cm2(1197) =       9.33
 given_sigma_tot_O2_cm2(1198) =       8.87
 given_sigma_tot_O2_cm2(1199) =       8.60
 given_sigma_tot_O2_cm2(1200) =       9.17
 given_sigma_tot_O2_cm2(1201) =       9.71
 given_sigma_tot_O2_cm2(1202) =      10.60
 given_sigma_tot_O2_cm2(1203) =      12.03
 given_sigma_tot_O2_cm2(1204) =      12.60
 given_sigma_tot_O2_cm2(1205) =      11.68
 given_sigma_tot_O2_cm2(1206) =      10.16
 given_sigma_tot_O2_cm2(1207) =       8.63
 given_sigma_tot_O2_cm2(1208) =       7.10
 given_sigma_tot_O2_cm2(1209) =       7.73
 given_sigma_tot_O2_cm2(1210) =       8.36
 given_sigma_tot_O2_cm2(1211) =       8.61
 given_sigma_tot_O2_cm2(1212) =       8.99
 given_sigma_tot_O2_cm2(1213) =       9.62
 given_sigma_tot_O2_cm2(1214) =      10.00
 given_sigma_tot_O2_cm2(1215) =       9.31
 given_sigma_tot_O2_cm2(1216) =       8.78
 given_sigma_tot_O2_cm2(1217) =       7.92
 given_sigma_tot_O2_cm2(1218) =       7.05
 given_sigma_tot_O2_cm2(1219) =       6.70
 given_sigma_tot_O2_cm2(1220) =       7.01
 given_sigma_tot_O2_cm2(1221) =       7.17
 given_sigma_tot_O2_cm2(1222) =       7.40
 given_sigma_tot_O2_cm2(1223) =       7.66
 given_sigma_tot_O2_cm2(1224) =       7.93
 given_sigma_tot_O2_cm2(1225) =       8.46
 given_sigma_tot_O2_cm2(1226) =       8.56
 given_sigma_tot_O2_cm2(1227) =       8.98
 given_sigma_tot_O2_cm2(1228) =       9.30
 given_sigma_tot_O2_cm2(1229) =       8.81
 given_sigma_tot_O2_cm2(1230) =       8.56
 given_sigma_tot_O2_cm2(1231) =       8.31
 given_sigma_tot_O2_cm2(1232) =       7.57
 given_sigma_tot_O2_cm2(1233) =       6.96
 given_sigma_tot_O2_cm2(1234) =       6.34
 given_sigma_tot_O2_cm2(1235) =       5.72
 given_sigma_tot_O2_cm2(1236) =       6.03
 given_sigma_tot_O2_cm2(1237) =       6.58
 given_sigma_tot_O2_cm2(1238) =       7.12
 given_sigma_tot_O2_cm2(1239) =       7.66
 given_sigma_tot_O2_cm2(1240) =       8.20
 given_sigma_tot_O2_cm2(1241) =       9.23
 given_sigma_tot_O2_cm2(1242) =       9.49
 given_sigma_tot_O2_cm2(1243) =      10.00
 given_sigma_tot_O2_cm2(1244) =       9.31
 given_sigma_tot_O2_cm2(1245) =       8.73
 given_sigma_tot_O2_cm2(1246) =       8.15
 given_sigma_tot_O2_cm2(1247) =       7.57
 given_sigma_tot_O2_cm2(1248) =       6.99
 given_sigma_tot_O2_cm2(1249) =       6.41
 given_sigma_tot_O2_cm2(1250) =       5.83
 given_sigma_tot_O2_cm2(1251) =       5.60
 given_sigma_tot_O2_cm2(1252) =       6.41
 given_sigma_tot_O2_cm2(1253) =       7.75
 given_sigma_tot_O2_cm2(1254) =       8.29
 given_sigma_tot_O2_cm2(1255) =       9.10
 given_sigma_tot_O2_cm2(1256) =      10.44
 given_sigma_tot_O2_cm2(1257) =      11.79
 given_sigma_tot_O2_cm2(1258) =      12.38
 given_sigma_tot_O2_cm2(1259) =      12.74
 given_sigma_tot_O2_cm2(1260) =      12.04
 given_sigma_tot_O2_cm2(1261) =      11.58
 given_sigma_tot_O2_cm2(1262) =      11.08
 given_sigma_tot_O2_cm2(1263) =      10.59
 given_sigma_tot_O2_cm2(1264) =      10.26
 given_sigma_tot_O2_cm2(1265) =       9.43
 given_sigma_tot_O2_cm2(1266) =       8.60
 given_sigma_tot_O2_cm2(1267) =       7.78
 given_sigma_tot_O2_cm2(1268) =       6.95
 given_sigma_tot_O2_cm2(1269) =       6.12
 given_sigma_tot_O2_cm2(1270) =       5.79
 given_sigma_tot_O2_cm2(1271) =       5.30
 given_sigma_tot_O2_cm2(1272) =       4.80
 given_sigma_tot_O2_cm2(1273) =       5.85
 given_sigma_tot_O2_cm2(1274) =       8.47
 given_sigma_tot_O2_cm2(1275) =      11.09
 given_sigma_tot_O2_cm2(1276) =      13.71
 given_sigma_tot_O2_cm2(1277) =      16.33
 given_sigma_tot_O2_cm2(1278) =      17.90
 given_sigma_tot_O2_cm2(1279) =      17.11
 given_sigma_tot_O2_cm2(1280) =      15.12
 given_sigma_tot_O2_cm2(1281) =      13.53
 given_sigma_tot_O2_cm2(1282) =      11.15
 given_sigma_tot_O2_cm2(1283) =       9.17
 given_sigma_tot_O2_cm2(1284) =       7.18
 given_sigma_tot_O2_cm2(1285) =       5.20
 given_sigma_tot_O2_cm2(1286) =       6.00
 given_sigma_tot_O2_cm2(1287) =       7.50
 given_sigma_tot_O2_cm2(1288) =       9.00
 given_sigma_tot_O2_cm2(1289) =      10.50
 given_sigma_tot_O2_cm2(1290) =       9.76
 given_sigma_tot_O2_cm2(1291) =       8.46
 given_sigma_tot_O2_cm2(1292) =       7.16
 given_sigma_tot_O2_cm2(1293) =       5.86
 given_sigma_tot_O2_cm2(1294) =       8.09
 given_sigma_tot_O2_cm2(1295) =      11.20
 given_sigma_tot_O2_cm2(1296) =       9.64
 given_sigma_tot_O2_cm2(1297) =       9.02
 given_sigma_tot_O2_cm2(1298) =       8.09
 given_sigma_tot_O2_cm2(1299) =       6.53
 given_sigma_tot_O2_cm2(1300) =       5.60
 given_sigma_tot_O2_cm2(1301) =       5.84
 given_sigma_tot_O2_cm2(1302) =       6.44
 given_sigma_tot_O2_cm2(1303) =       7.04
 given_sigma_tot_O2_cm2(1304) =       7.28
 given_sigma_tot_O2_cm2(1305) =       7.03
 given_sigma_tot_O2_cm2(1306) =       6.10
 given_sigma_tot_O2_cm2(1307) =       5.17
 given_sigma_tot_O2_cm2(1308) =       4.80
 given_sigma_tot_O2_cm2(1309) =       5.62
 given_sigma_tot_O2_cm2(1310) =       6.99
 given_sigma_tot_O2_cm2(1311) =       8.35
 given_sigma_tot_O2_cm2(1312) =       8.90
 given_sigma_tot_O2_cm2(1313) =      10.40
 given_sigma_tot_O2_cm2(1314) =      12.90
 given_sigma_tot_O2_cm2(1315) =      12.58
 given_sigma_tot_O2_cm2(1316) =      11.56
 given_sigma_tot_O2_cm2(1317) =       9.70
 given_sigma_tot_O2_cm2(1318) =      10.31
 given_sigma_tot_O2_cm2(1319) =      10.62
 given_sigma_tot_O2_cm2(1320) =      10.80
 given_sigma_tot_O2_cm2(1321) =      10.28
 given_sigma_tot_O2_cm2(1322) =       9.00
 given_sigma_tot_O2_cm2(1323) =       7.71
 given_sigma_tot_O2_cm2(1324) =       6.42
 given_sigma_tot_O2_cm2(1325) =       5.13
 given_sigma_tot_O2_cm2(1326) =       4.10
 given_sigma_tot_O2_cm2(1327) =       6.54
 given_sigma_tot_O2_cm2(1328) =       8.16
 given_sigma_tot_O2_cm2(1329) =      10.60
 given_sigma_tot_O2_cm2(1330) =      12.63
 given_sigma_tot_O2_cm2(1331) =      14.66
 given_sigma_tot_O2_cm2(1332) =      16.69
 given_sigma_tot_O2_cm2(1333) =      16.15
 given_sigma_tot_O2_cm2(1334) =      15.20
 given_sigma_tot_O2_cm2(1335) =      17.04
 given_sigma_tot_O2_cm2(1336) =      15.23
 given_sigma_tot_O2_cm2(1337) =      12.96
 given_sigma_tot_O2_cm2(1338) =      11.78
 given_sigma_tot_O2_cm2(1339) =      10.70
 given_sigma_tot_O2_cm2(1340) =       9.24
 given_sigma_tot_O2_cm2(1341) =       8.43
 given_sigma_tot_O2_cm2(1342) =       6.16
 given_sigma_tot_O2_cm2(1343) =       4.80
 given_sigma_tot_O2_cm2(1344) =       5.17
 given_sigma_tot_O2_cm2(1345) =       6.10
 given_sigma_tot_O2_cm2(1346) =       7.03
 given_sigma_tot_O2_cm2(1347) =       7.40
 given_sigma_tot_O2_cm2(1348) =       6.30
 given_sigma_tot_O2_cm2(1349) =       4.47
 given_sigma_tot_O2_cm2(1350) =       8.93
 given_sigma_tot_O2_cm2(1351) =      14.96
 given_sigma_tot_O2_cm2(1352) =      20.99
 given_sigma_tot_O2_cm2(1353) =      23.40
 given_sigma_tot_O2_cm2(1354) =      21.49
 given_sigma_tot_O2_cm2(1355) =      18.30
 given_sigma_tot_O2_cm2(1356) =      15.11
 given_sigma_tot_O2_cm2(1357) =      12.56
 given_sigma_tot_O2_cm2(1358) =       8.74
 given_sigma_tot_O2_cm2(1359) =       6.19
 given_sigma_tot_O2_cm2(1360) =       3.00
 given_sigma_tot_O2_cm2(1361) =       4.39
 given_sigma_tot_O2_cm2(1362) =       5.73
 given_sigma_tot_O2_cm2(1363) =       6.97
 given_sigma_tot_O2_cm2(1364) =       8.21
 given_sigma_tot_O2_cm2(1365) =       8.96
 given_sigma_tot_O2_cm2(1366) =       9.45
 given_sigma_tot_O2_cm2(1367) =       8.90
 given_sigma_tot_O2_cm2(1368) =      11.80
 given_sigma_tot_O2_cm2(1369) =      16.15
 given_sigma_tot_O2_cm2(1370) =      20.50
 given_sigma_tot_O2_cm2(1371) =      23.40
 given_sigma_tot_O2_cm2(1372) =      18.22
 given_sigma_tot_O2_cm2(1373) =      13.03
 given_sigma_tot_O2_cm2(1374) =       7.85
 given_sigma_tot_O2_cm2(1375) =       5.77
 given_sigma_tot_O2_cm2(1376) =       3.70
 given_sigma_tot_O2_cm2(1377) =       3.90
 given_sigma_tot_O2_cm2(1378) =       4.10
 given_sigma_tot_O2_cm2(1379) =       3.46
 given_sigma_tot_O2_cm2(1380) =       3.62
 given_sigma_tot_O2_cm2(1381) =       4.02
 given_sigma_tot_O2_cm2(1382) =       3.70
 given_sigma_tot_O2_cm2(1383) =      22.28
 given_sigma_tot_O2_cm2(1384) =      23.72
 given_sigma_tot_O2_cm2(1385) =      19.91
 given_sigma_tot_O2_cm2(1386) =      12.30
 given_sigma_tot_O2_cm2(1387) =      19.54
 given_sigma_tot_O2_cm2(1388) =      28.60
 given_sigma_tot_O2_cm2(1389) =      23.76
 given_sigma_tot_O2_cm2(1390) =      20.70
 given_sigma_tot_O2_cm2(1391) =      19.73
 given_sigma_tot_O2_cm2(1392) =      15.70
 given_sigma_tot_O2_cm2(1393) =      11.67
 given_sigma_tot_O2_cm2(1394) =       7.64
 given_sigma_tot_O2_cm2(1395) =       5.22
 given_sigma_tot_O2_cm2(1396) =       3.61
 given_sigma_tot_O2_cm2(1397) =       7.36
 given_sigma_tot_O2_cm2(1398) =      13.07
 given_sigma_tot_O2_cm2(1399) =      18.77
 given_sigma_tot_O2_cm2(1400) =      24.47
 given_sigma_tot_O2_cm2(1401) =      27.89
 given_sigma_tot_O2_cm2(1402) =      29.03
 given_sigma_tot_O2_cm2(1403) =      35.88
 given_sigma_tot_O2_cm2(1404) =      41.58
 given_sigma_tot_O2_cm2(1405) =      45.00
 given_sigma_tot_O2_cm2(1406) =      43.40
 given_sigma_tot_O2_cm2(1407) =      39.40
 given_sigma_tot_O2_cm2(1408) =      35.40
 given_sigma_tot_O2_cm2(1409) =      31.40
 given_sigma_tot_O2_cm2(1410) =      27.40
 given_sigma_tot_O2_cm2(1411) =      23.40
 given_sigma_tot_O2_cm2(1412) =      20.20
 given_sigma_tot_O2_cm2(1413) =      15.40
 given_sigma_tot_O2_cm2(1414) =      13.00
 given_sigma_tot_O2_cm2(1415) =      11.40
 given_sigma_tot_O2_cm2(1416) =       7.40
 given_sigma_tot_O2_cm2(1417) =       3.40
 given_sigma_tot_O2_cm2(1418) =       9.52
 given_sigma_tot_O2_cm2(1419) =      26.81
 given_sigma_tot_O2_cm2(1420) =      44.10
 given_sigma_tot_O2_cm2(1421) =      56.20
 given_sigma_tot_O2_cm2(1422) =      49.97
 given_sigma_tot_O2_cm2(1423) =      39.58
 given_sigma_tot_O2_cm2(1424) =      29.20
 given_sigma_tot_O2_cm2(1425) =      18.82
 given_sigma_tot_O2_cm2(1426) =      13.83
 given_sigma_tot_O2_cm2(1427) =       8.43
 given_sigma_tot_O2_cm2(1428) =       2.20
 given_sigma_tot_O2_cm2(1429) =       4.08
 given_sigma_tot_O2_cm2(1430) =       8.76
 given_sigma_tot_O2_cm2(1431) =      18.14
 given_sigma_tot_O2_cm2(1432) =      27.51
 given_sigma_tot_O2_cm2(1433) =      36.89
 given_sigma_tot_O2_cm2(1434) =      46.26
 given_sigma_tot_O2_cm2(1435) =      54.70
 given_sigma_tot_O2_cm2(1436) =      35.43
 given_sigma_tot_O2_cm2(1437) =      29.00
 given_sigma_tot_O2_cm2(1438) =      35.30
 given_sigma_tot_O2_cm2(1439) =      29.83
 given_sigma_tot_O2_cm2(1440) =      24.37
 given_sigma_tot_O2_cm2(1441) =      22.18
 given_sigma_tot_O2_cm2(1442) =      18.90
 given_sigma_tot_O2_cm2(1443) =      15.62
 given_sigma_tot_O2_cm2(1444) =      13.43
 given_sigma_tot_O2_cm2(1445) =       7.97
 given_sigma_tot_O2_cm2(1446) =       2.50
 given_sigma_tot_O2_cm2(1447) =       5.95
 given_sigma_tot_O2_cm2(1448) =       9.39
 given_sigma_tot_O2_cm2(1449) =      12.84
 given_sigma_tot_O2_cm2(1450) =      15.60
 given_sigma_tot_O2_cm2(1451) =       9.67
 given_sigma_tot_O2_cm2(1452) =       6.70
 given_sigma_tot_O2_cm2(1453) =       9.98
 given_sigma_tot_O2_cm2(1454) =      26.39
 given_sigma_tot_O2_cm2(1455) =      42.80
 given_sigma_tot_O2_cm2(1456) =      51.00
 given_sigma_tot_O2_cm2(1457) =      46.04
 given_sigma_tot_O2_cm2(1458) =      41.08
 given_sigma_tot_O2_cm2(1459) =      36.12
 given_sigma_tot_O2_cm2(1460) =      31.16
 given_sigma_tot_O2_cm2(1461) =      26.20
 given_sigma_tot_O2_cm2(1462) =      16.29
 given_sigma_tot_O2_cm2(1463) =      11.33
 given_sigma_tot_O2_cm2(1464) =       6.37
 given_sigma_tot_O2_cm2(1465) =       2.40
 given_sigma_tot_O2_cm2(1466) =      10.74
 given_sigma_tot_O2_cm2(1467) =      17.70
 given_sigma_tot_O2_cm2(1468) =      24.65
 given_sigma_tot_O2_cm2(1469) =      31.60
 given_sigma_tot_O2_cm2(1470) =      42.80
 given_sigma_tot_O2_cm2(1471) =      28.85
 given_sigma_tot_O2_cm2(1472) =      17.23
 given_sigma_tot_O2_cm2(1473) =       5.60
 given_sigma_tot_O2_cm2(1474) =      18.85
 given_sigma_tot_O2_cm2(1475) =      26.80
 given_sigma_tot_O2_cm2(1476) =      25.86
 given_sigma_tot_O2_cm2(1477) =      23.52
 given_sigma_tot_O2_cm2(1478) =      21.17
 given_sigma_tot_O2_cm2(1479) =      18.82
 given_sigma_tot_O2_cm2(1480) =      16.48
 given_sigma_tot_O2_cm2(1481) =      14.13
 given_sigma_tot_O2_cm2(1482) =      11.78
 given_sigma_tot_O2_cm2(1483) =       9.44
 given_sigma_tot_O2_cm2(1484) =       7.09
 given_sigma_tot_O2_cm2(1485) =       4.75
 given_sigma_tot_O2_cm2(1486) =       2.40
 given_sigma_tot_O2_cm2(1487) =      10.20
 given_sigma_tot_O2_cm2(1488) =      18.01
 given_sigma_tot_O2_cm2(1489) =      25.81
 given_sigma_tot_O2_cm2(1490) =      33.61
 given_sigma_tot_O2_cm2(1491) =      41.42
 given_sigma_tot_O2_cm2(1492) =      46.10
 given_sigma_tot_O2_cm2(1493) =      41.75
 given_sigma_tot_O2_cm2(1494) =      30.88
 given_sigma_tot_O2_cm2(1495) =      20.02
 given_sigma_tot_O2_cm2(1496) =       9.15
 given_sigma_tot_O2_cm2(1497) =       4.80
 given_sigma_tot_O2_cm2(1498) =       5.91
 given_sigma_tot_O2_cm2(1499) =       7.40
 given_sigma_tot_O2_cm2(1500) =       6.56
 given_sigma_tot_O2_cm2(1501) =       5.10
 given_sigma_tot_O2_cm2(1502) =       3.00
 given_sigma_tot_O2_cm2(1503) =       4.80
 given_sigma_tot_O2_cm2(1504) =       3.30
 given_sigma_tot_O2_cm2(1505) =       1.50
 given_sigma_tot_O2_cm2(1506) =       2.76
 given_sigma_tot_O2_cm2(1507) =       4.15
 given_sigma_tot_O2_cm2(1508) =      10.79
 given_sigma_tot_O2_cm2(1509) =      14.11
 given_sigma_tot_O2_cm2(1510) =      17.43
 given_sigma_tot_O2_cm2(1511) =      23.40
 given_sigma_tot_O2_cm2(1512) =      21.60
 given_sigma_tot_O2_cm2(1513) =      24.50
 given_sigma_tot_O2_cm2(1514) =      21.39
 given_sigma_tot_O2_cm2(1515) =      15.16
 given_sigma_tot_O2_cm2(1516) =       8.93
 given_sigma_tot_O2_cm2(1517) =       2.70
 given_sigma_tot_O2_cm2(1518) =       1.45
 given_sigma_tot_O2_cm2(1519) =       1.46
 given_sigma_tot_O2_cm2(1520) =       1.48
 given_sigma_tot_O2_cm2(1521) =       1.49
 given_sigma_tot_O2_cm2(1522) =       2.69
 given_sigma_tot_O2_cm2(1523) =       3.90
 given_sigma_tot_O2_cm2(1524) =       5.10
 given_sigma_tot_O2_cm2(1525) =       6.30
 given_sigma_tot_O2_cm2(1526) =       5.60
 given_sigma_tot_O2_cm2(1527) =       6.30
 given_sigma_tot_O2_cm2(1528) =       5.41
 given_sigma_tot_O2_cm2(1529) =       3.17
 given_sigma_tot_O2_cm2(1530) =       1.38
 given_sigma_tot_O2_cm2(1531) =       1.46
 given_sigma_tot_O2_cm2(1532) =       1.82
 given_sigma_tot_O2_cm2(1533) =       1.44
 given_sigma_tot_O2_cm2(1534) =       1.52
 given_sigma_tot_O2_cm2(1535) =       1.45
 given_sigma_tot_O2_cm2(1536) =       1.43
 given_sigma_tot_O2_cm2(1537) =       1.34
 given_sigma_tot_O2_cm2(1538) =       1.30
 given_sigma_tot_O2_cm2(1539) =       1.15
 given_sigma_tot_O2_cm2(1540) =       1.08
 given_sigma_tot_O2_cm2(1541) =       1.23
 given_sigma_tot_O2_cm2(1542) =       1.34
 given_sigma_tot_O2_cm2(1543) =       1.12
 given_sigma_tot_O2_cm2(1544) =       1.48
 given_sigma_tot_O2_cm2(1545) =       1.75
 given_sigma_tot_O2_cm2(1546) =       1.19
 given_sigma_tot_O2_cm2(1547) =       1.38
 given_sigma_tot_O2_cm2(1548) =       1.08
 given_sigma_tot_O2_cm2(1549) =       1.60
 given_sigma_tot_O2_cm2(1550) =       0.97
 given_sigma_tot_O2_cm2(1551) =       1.19
 given_sigma_tot_O2_cm2(1552) =       1.52
 given_sigma_tot_O2_cm2(1553) =       1.08
 given_sigma_tot_O2_cm2(1554) =       1.19
 given_sigma_tot_O2_cm2(1555) =       1.41
 given_sigma_tot_O2_cm2(1556) =       1.28
 given_sigma_tot_O2_cm2(1557) =       1.19
 given_sigma_tot_O2_cm2(1558) =       1.60
 given_sigma_tot_O2_cm2(1559) =       1.40
 given_sigma_tot_O2_cm2(1560) =       1.64
 given_sigma_tot_O2_cm2(1561) =       1.47
 given_sigma_tot_O2_cm2(1562) =       1.30
 given_sigma_tot_O2_cm2(1563) =       1.64
 given_sigma_tot_O2_cm2(1564) =       1.86
 given_sigma_tot_O2_cm2(1565) =       1.43
 given_sigma_tot_O2_cm2(1566) =       1.00
 given_sigma_tot_O2_cm2(1567) =       1.45
 given_sigma_tot_O2_cm2(1568) =       1.79
 given_sigma_tot_O2_cm2(1569) =       1.64
 given_sigma_tot_O2_cm2(1570) =       1.52
 given_sigma_tot_O2_cm2(1571) =       1.11

  DO i = 1, N_solar_bins

     CALL get_the_j_and_a(lambda_EUV_A(i), N_given_xs_points, given_lambda_A, j, aj, ajp1)

     sigma_tot_O_cm2(i)  = aj * given_sigma_tot_O_cm2(j)  + ajp1 * given_sigma_tot_O_cm2(j+1)
     sigma_tot_N2_cm2(i) = aj * given_sigma_tot_N2_cm2(j) + ajp1 * given_sigma_tot_N2_cm2(j+1)
     sigma_tot_O2_cm2(i) = aj * given_sigma_tot_O2_cm2(j) + ajp1 * given_sigma_tot_O2_cm2(j+1)

     sigma_ion_O_4S_cm2(i)   = aj * given_sigma_ion_O_4S_cm2(j)   + ajp1 * given_sigma_ion_O_4S_cm2(j+1)
     sigma_ion_O_2D_cm2(i)   = aj * given_sigma_ion_O_2D_cm2(j)   + ajp1 * given_sigma_ion_O_2D_cm2(j+1)
     sigma_ion_O_2P_cm2(i)   = aj * given_sigma_ion_O_2P_cm2(j)   + ajp1 * given_sigma_ion_O_2P_cm2(j+1)
     sigma_ion_O_4Pst_cm2(i) = aj * given_sigma_ion_O_4Pst_cm2(j) + ajp1 * given_sigma_ion_O_4Pst_cm2(j+1)
     sigma_ion_O_2Pst_cm2(i) = aj * given_sigma_ion_O_2Pst_cm2(j) + ajp1 * given_sigma_ion_O_2Pst_cm2(j+1)

     sigma_ion_N2_to_N2_cm2(i) = aj * given_sigma_ion_N2_to_N2_cm2(j) + ajp1 * given_sigma_ion_N2_to_N2_cm2(j+1)
     sigma_ion_N2_to_N_cm2(i)  = aj * given_sigma_ion_N2_to_N_cm2(j)  + ajp1 * given_sigma_ion_N2_to_N_cm2(j+1)

     sigma_ion_O2_to_O2_cm2(i) = aj * given_sigma_ion_O2_to_O2_cm2(j) + ajp1 * given_sigma_ion_O2_to_O2_cm2(j+1)
     sigma_ion_O2_to_O_cm2(i)  = aj * given_sigma_ion_O2_to_O_cm2(j)  + ajp1 * given_sigma_ion_O2_to_O_cm2(j+1)

  END DO

!###########

  IF (Rank_of_process.EQ.0) THEN

     SELECT CASE (EUV_spectrum_model)
     CASE (f76)
        OPEN (10, FILE = 'solar_f76_fluxes_crossections.dat')
     CASE (f74113)
        OPEN (10, FILE = 'solar_f74113_fluxes_crossections.dat')
     END SELECT

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

     DO i = 1, N_solar_bins

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

     SELECT CASE (EUV_spectrum_model)
     CASE (f76)
        PRINT '("### file solar_f76_fluxes_crossections.dat is ready")'
     CASE (f74113)
        PRINT '("### file solar_f74113_fluxes_crossections.dat is ready")'
     END SELECT

  END IF   !###   IF (Rank_of_process.EQ.0) THEN

! final cleanup

  DEALLOCATE(lambda_EUV_A, STAT = ALLOC_ERR)

  DEALLOCATE(given_lambda_A, STAT = ALLOC_ERR)

  DEALLOCATE(given_sigma_tot_O_cm2, STAT = ALLOC_ERR)            ! total cross section [10^-18 cm^2] for atomic oxygen          
  DEALLOCATE(given_sigma_tot_N2_cm2, STAT = ALLOC_ERR)           ! - " -, molecular nitrogen             
  DEALLOCATE(given_sigma_tot_O2_cm2, STAT = ALLOC_ERR)           ! - " -, molecular oxygen

! used to calculate ionization rates
  DEALLOCATE(given_sigma_ion_N2_to_N2_cm2, STAT = ALLOC_ERR)     ! ionization cross section [10^-18 cm^2] for molecular nitrogen, N2+
  DEALLOCATE(given_sigma_ion_N2_to_N_cm2, STAT = ALLOC_ERR)      ! - " -, molecular nitrogen, N+

  DEALLOCATE(given_sigma_ion_O_4S_cm2, STAT = ALLOC_ERR)         ! - " -, atomic oxygen, O+4S
  DEALLOCATE(given_sigma_ion_O_2D_cm2, STAT = ALLOC_ERR)         ! - " -, atomic oxygen, O+2D
  DEALLOCATE(given_sigma_ion_O_2P_cm2, STAT = ALLOC_ERR)         ! - " -, atomic oxygen, O+2P
  DEALLOCATE(given_sigma_ion_O_4Pst_cm2, STAT = ALLOC_ERR)       ! - " -, atomic oxygen, O+4P*
  DEALLOCATE(given_sigma_ion_O_2Pst_cm2, STAT = ALLOC_ERR)       ! - " -, atomic oxygen, O+2P*

  DEALLOCATE(given_sigma_ion_O2_to_O2_cm2, STAT = ALLOC_ERR)     ! - " -, molecular oxygen, O2+
  DEALLOCATE(given_sigma_ion_O2_to_O_cm2, STAT = ALLOC_ERR)      ! - " -, molecular oxygen, O+

END SUBROUTINE Initiate_EUV_fluxes_crossections

!-----------------------------------
! 
subroutine get_the_j_and_a(x0, Nx, xarr, j, aj, ajp1)

  implicit none

  real,    intent(in) :: x0
  integer, intent(in) :: Nx
  real,    intent(in) :: xarr(1:Nx)

  integer, intent(out) :: j
  real, intent(out) :: aj, ajp1

  integer jleft, jmid, jright

! default values
  j=1
  aj = 0.0
  ajp1 = 0.0

  if (x0.lt.xarr(1)) return
  if (x0.gt.xarr(Nx)) return

  if (x0.eq.xarr(1)) then
     j=1
     aj   = 1.0
     ajp1 = 0.0
     return
  end if

  if (x0.eq.xarr(Nx)) then
     j=Nx-1
     aj   = 0.0
     ajp1 = 1.0
     return
  end if

  jleft=1
  jright=Nx
  jmid = max(jleft,min(jright-1,(jleft+jright)/2))
  do while ((jright-jleft).gt.1)
     if (x0.lt.xarr(jmid)) then
        jright = jmid
     else
        jleft = jmid
     end if
     jmid = max(jleft,min(jright-1,(jleft+jright)/2))
  end do
  j = jleft
  
  if ((x0.lt.xarr(jleft)).or.(x0.gt.xarr(jleft+1))) then
     print '("error in get_the_j_and_a",2(2x,i5),3(2x,f9.3))', jleft, jright, xarr(jleft), xarr(jright), x0
     stop
  end if

  aj = max(0.0, (xarr(j+1) - x0) / (xarr(j+1) - xarr(j)))
  ajp1 = max(0.0, min(1.0, 1.0-aj))

  return

end subroutine get_the_j_and_a
