
!---------------------------------------------------------------------------------------------------- #1
!
! cross_section_check/output/data-01/Cross-section-N2-elastic-02.dat
! A. Schmalzried, A. Luque and N. Lehtinen,  IAA Database on lxcat, www.lxcat.net/IAA, August 2023,
!                  Instituto de Astrofísica de Andalucía.
!
module N2_elastic
!  integer, parameter :: N_td=35  ! number of table data points
  integer, parameter :: N_td=67  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_elastic
!
!-------------------------------------
!
subroutine Prepare_N2_elastic

  use N2_elastic
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n
  td_E_eV(  1)  = 0.000000d+00 ; td_cs_cm2(  1) = 1.080000d-20   !### in m^2 here ###
  td_E_eV(  2)  = 1.000000d-03 ; td_cs_cm2(  2) = 1.168220d-20
  td_E_eV(  3)  = 1.000000d-02 ; td_cs_cm2(  3) = 2.961780d-20
  td_E_eV(  4)  = 1.000000d-01 ; td_cs_cm2(  4) = 4.755480d-20
  td_E_eV(  5)  = 2.000000d-01 ; td_cs_cm2(  5) = 6.341130d-20
  td_E_eV(  6)  = 3.000000d-01 ; td_cs_cm2(  6) = 7.268810d-20
  td_E_eV(  7)  = 3.500000d-01 ; td_cs_cm2(  7) = 7.621590d-20
  td_E_eV(  8)  = 4.000000d-01 ; td_cs_cm2(  8) = 7.939550d-20
  td_E_eV(  9)  = 5.000000d-01 ; td_cs_cm2(  9) = 8.470200d-20
  td_E_eV( 10)  = 5.500000d-01 ; td_cs_cm2( 10) = 8.697660d-20
  td_E_eV( 11)  = 6.000000d-01 ; td_cs_cm2( 11) = 8.918840d-20
  td_E_eV( 12)  = 7.000000d-01 ; td_cs_cm2( 12) = 9.311810d-20
  td_E_eV( 13)  = 8.000000d-01 ; td_cs_cm2( 13) = 9.651760d-20
  td_E_eV( 14)  = 9.000000d-01 ; td_cs_cm2( 14) = 9.952120d-20
  td_E_eV( 15)  = 1.000000d+00 ; td_cs_cm2( 15) = 1.022530d-19
  td_E_eV( 16)  = 1.500000d+00 ; td_cs_cm2( 16) = 1.175130d-19
  td_E_eV( 17)  = 1.900000d+00 ; td_cs_cm2( 17) = 1.704020d-19
  td_E_eV( 18)  = 1.950000d+00 ; td_cs_cm2( 18) = 1.739240d-19
  td_E_eV( 19)  = 2.460000d+00 ; td_cs_cm2( 19) = 2.458630d-19
  td_E_eV( 20)  = 2.610000d+00 ; td_cs_cm2( 20) = 1.531640d-19
  td_E_eV( 21)  = 3.000000d+00 ; td_cs_cm2( 21) = 1.623140d-19
  td_E_eV( 22)  = 4.000000d+00 ; td_cs_cm2( 22) = 1.249160d-19
  td_E_eV( 23)  = 5.000000d+00 ; td_cs_cm2( 23) = 1.118010d-19
  td_E_eV( 24)  = 6.000000d+00 ; td_cs_cm2( 24) = 1.070430d-19
  td_E_eV( 25)  = 7.000000d+00 ; td_cs_cm2( 25) = 1.030850d-19
  td_E_eV( 26)  = 8.000000d+00 ; td_cs_cm2( 26) = 1.024230d-19
  td_E_eV( 27)  = 1.000000d+01 ; td_cs_cm2( 27) = 1.089360d-19
  td_E_eV( 28)  = 1.140000d+01 ; td_cs_cm2( 28) = 1.040610d-19
  td_E_eV( 29)  = 1.500000d+01 ; td_cs_cm2( 29) = 1.125470d-19
  td_E_eV( 30)  = 2.000000d+01 ; td_cs_cm2( 30) = 1.093640d-19
  td_E_eV( 31)  = 3.000000d+01 ; td_cs_cm2( 31) = 1.164420d-19
  td_E_eV( 32)  = 3.500000d+01 ; td_cs_cm2( 32) = 1.084220d-19
  td_E_eV( 33)  = 4.000000d+01 ; td_cs_cm2( 33) = 9.888670d-20
  td_E_eV( 34)  = 4.500000d+01 ; td_cs_cm2( 34) = 8.960070d-20
  td_E_eV( 35)  = 5.000000d+01 ; td_cs_cm2( 35) = 8.236830d-20
  td_E_eV( 36)  = 5.500000d+01 ; td_cs_cm2( 36) = 7.676580d-20
  td_E_eV( 37)  = 6.000000d+01 ; td_cs_cm2( 37) = 7.224120d-20
  td_E_eV( 38)  = 6.500000d+01 ; td_cs_cm2( 38) = 6.847180d-20
  td_E_eV( 39)  = 7.000000d+01 ; td_cs_cm2( 39) = 6.527350d-20
  td_E_eV( 40)  = 7.500000d+01 ; td_cs_cm2( 40) = 6.251960d-20
  td_E_eV( 41)  = 8.000000d+01 ; td_cs_cm2( 41) = 6.013290d-20
  td_E_eV( 42)  = 8.500000d+01 ; td_cs_cm2( 42) = 5.804660d-20
  td_E_eV( 43)  = 9.000000d+01 ; td_cs_cm2( 43) = 5.621150d-20
  td_E_eV( 44)  = 9.500000d+01 ; td_cs_cm2( 44) = 5.457800d-20
  td_E_eV( 45)  = 1.000000d+02 ; td_cs_cm2( 45) = 5.310570d-20
  td_E_eV( 46)  = 1.100000d+02 ; td_cs_cm2( 46) = 5.053860d-20
  td_E_eV( 47)  = 1.200000d+02 ; td_cs_cm2( 47) = 4.831070d-20
  td_E_eV( 48)  = 1.300000d+02 ; td_cs_cm2( 48) = 4.631010d-20
  td_E_eV( 49)  = 1.400000d+02 ; td_cs_cm2( 49) = 4.447190d-20
  td_E_eV( 50)  = 1.500000d+02 ; td_cs_cm2( 50) = 4.276690d-20
  td_E_eV( 51)  = 2.000000d+02 ; td_cs_cm2( 51) = 3.839290d-20
  td_E_eV( 52)  = 2.500000d+02 ; td_cs_cm2( 52) = 3.281180d-20
  td_E_eV( 53)  = 3.000000d+02 ; td_cs_cm2( 53) = 2.896860d-20
  td_E_eV( 54)  = 3.500000d+02 ; td_cs_cm2( 54) = 2.601360d-20
  td_E_eV( 55)  = 4.000000d+02 ; td_cs_cm2( 55) = 2.371470d-20
  td_E_eV( 56)  = 4.500000d+02 ; td_cs_cm2( 56) = 2.179860d-20
  td_E_eV( 57)  = 5.000000d+02 ; td_cs_cm2( 57) = 2.019370d-20
  td_E_eV( 58)  = 5.500000d+02 ; td_cs_cm2( 58) = 1.882640d-20
  td_E_eV( 59)  = 6.000000d+02 ; td_cs_cm2( 59) = 1.764220d-20
  td_E_eV( 60)  = 6.500000d+02 ; td_cs_cm2( 60) = 1.660500d-20
  td_E_eV( 61)  = 7.000000d+02 ; td_cs_cm2( 61) = 1.568630d-20
  td_E_eV( 62)  = 7.500000d+02 ; td_cs_cm2( 62) = 1.486770d-20
  td_E_eV( 63)  = 8.000000d+02 ; td_cs_cm2( 63) = 1.413360d-20
  td_E_eV( 64)  = 8.500000d+02 ; td_cs_cm2( 64) = 1.362870d-20
  td_E_eV( 65)  = 9.000000d+02 ; td_cs_cm2( 65) = 1.317180d-20
  td_E_eV( 66)  = 9.500000d+02 ; td_cs_cm2( 66) = 1.282650d-20
  td_E_eV( 67)  = 1.000000d+03 ; td_cs_cm2( 67) = 1.251150d-20

  td_cs_cm2 = 1.0d4 * td_cs_cm2 ! the original data are in m^2, convert to cm^2

! N2, elastic scattering cross section [NOT the momentum transfer cross section]  
! Table 3 of [Y.Itikawa, J.Phys.Chem.Ref.Data, v.35, p.31-53, 2006]! 
!
!  td_E_eV(1)   =    0.55_8  ; td_cs_cm2(1)   = 8.390d-16 
!  td_E_eV(2)   =    0.70_8  ; td_cs_cm2(2)   = 9.030d-16
!  td_E_eV(3)   =    0.90_8  ; td_cs_cm2(3)   = 9.620d-16
!  td_E_eV(4)   =    1.00_8  ; td_cs_cm2(4)   = 9.830d-16
!  td_E_eV(5)   =    1.50_8  ; td_cs_cm2(5)   = 1.053d-15
!  td_E_eV(6)   =    2.00_8  ; td_cs_cm2(6)   = 1.793d-15
!  td_E_eV(7)   =    2.20_8  ; td_cs_cm2(7)   = 1.950d-15
!  td_E_eV(8)   =    2.35_8  ; td_cs_cm2(8)   = 2.050d-15
!  td_E_eV(9)   =    2.50_8  ; td_cs_cm2(9)   = 2.100d-15
!  td_E_eV(10)  =    2.70_8  ; td_cs_cm2(10)  = 1.750d-15
!  td_E_eV(11)  =    3.00_8  ; td_cs_cm2(11)  = 1.500d-15
!  td_E_eV(12)  =    4.00_8  ; td_cs_cm2(12)  = 1.160d-15
!  td_E_eV(13)  =    5.00_8  ; td_cs_cm2(13)  = 1.075d-15
!  td_E_eV(14)  =    6.00_8  ; td_cs_cm2(14)  = 1.060d-15
!  td_E_eV(15)  =    8.00_8  ; td_cs_cm2(15)  = 1.060d-15
!  td_E_eV(16)  =   10.00_8  ; td_cs_cm2(16)  = 1.140d-15
!  td_E_eV(17)  =   15.00_8  ; td_cs_cm2(17)  = 1.180d-15
!  td_E_eV(18)  =   20.00_8  ; td_cs_cm2(18)  = 1.115d-15
!  td_E_eV(19)  =   25.00_8  ; td_cs_cm2(19)  = 1.025d-15
!  td_E_eV(20)  =   30.00_8  ; td_cs_cm2(20)  = 9.650d-16
!  td_E_eV(21)  =   40.00_8  ; td_cs_cm2(21)  = 8.850d-16
!  td_E_eV(22)  =   50.00_8  ; td_cs_cm2(22)  = 8.200d-16
!  td_E_eV(23)  =   60.00_8  ; td_cs_cm2(23)  = 7.400d-16
!  td_E_eV(24)  =   80.00_8  ; td_cs_cm2(24)  = 6.250d-16
!  td_E_eV(25)  =  100.00_8  ; td_cs_cm2(25)  = 5.600d-16
!  td_E_eV(26)  =  120.00_8  ; td_cs_cm2(26)  = 4.900d-16
!  td_E_eV(27)  =  150.00_8  ; td_cs_cm2(27)  = 4.200d-16
!  td_E_eV(28)  =  200.00_8  ; td_cs_cm2(28)  = 3.500d-16
!  td_E_eV(29)  =  250.00_8  ; td_cs_cm2(29)  = 3.000d-16
!  td_E_eV(30)  =  300.00_8  ; td_cs_cm2(30)  = 2.650d-16
!  td_E_eV(31)  =  400.00_8  ; td_cs_cm2(31)  = 2.150d-16
!  td_E_eV(32)  =  500.00_8  ; td_cs_cm2(32)  = 1.850d-16
!  td_E_eV(33)  =  600.00_8  ; td_cs_cm2(33)  = 1.600d-16
!  td_E_eV(34)  =  800.00_8  ; td_cs_cm2(34)  = 1.250d-16
!  td_E_eV(35)  = 1000.00_8  ; td_cs_cm2(35)  = 1.000d-16
  
  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_elastic
!
!-------------------------------------
! above 1000 eV the cross section is calculated using
! an interplation formula derived upon assuming that in log-log scale the 
! cross section decays with energy as a straight line 
! continuing the tabulated data
! (see Figure 2 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
real(8) function CSV_N2_elastic_m3s(E_eV)

  use N2_elastic
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  real(8) CS_N2_elastic_cm2

  if (E_eV.lt.td_E_eV(1)) then
     CSV_N2_elastic_m3s = 1.0d-4 * td_cs_cm2(1) * sqrt(E_eV) * efactor_eV_to_ms
     return
  end if

  if (E_eV.GT.td_E_eV(N_td)) then
     CS_N2_elastic_cm2 = (1.0d-15 * log(E_eV) + 2.4d-14) / E_eV**0.83
     CSV_N2_elastic_m3s = max(0.0_8, 1.0d-4 * CS_N2_elastic_cm2 * sqrt(E_eV) * efactor_eV_to_ms)
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_elastic_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_elastic_m3s

!---------------------------------------------------------------------------------------------------- #2
! N2, {N2->N2+} ionization cross-section
! Tables 15,16,17 of [Y.Itikawa, J.Phys.Chem.Ref.Data, v.35, p.31-53, 2006]
!
module N2_ionN2
  integer, parameter :: N_td=58  ! number of table data points
  ! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_ionN2
!
!-------------------------------------
!
subroutine Prepare_N2_ionN2

  use N2_ionN2
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

! table 15 of [Itikawa]
  td_E_eV(1)  =   16.0_8; td_cs_cm2(1)  = 2.11d-18
  td_E_eV(2)  =   16.5_8; td_cs_cm2(2)  = 4.66d-18
  td_E_eV(3)  =   17.0_8; td_cs_cm2(3)  = 7.13d-18
  td_E_eV(4)  =   17.5_8; td_cs_cm2(4)  = 9.85d-18
  td_E_eV(5)  =   18.0_8; td_cs_cm2(5)  = 1.29d-17
  td_E_eV(6)  =   18.5_8; td_cs_cm2(6)  = 1.64d-17
  td_E_eV(7)  =   19.0_8; td_cs_cm2(7)  = 1.99d-17
  td_E_eV(8)  =   19.5_8; td_cs_cm2(8)  = 2.30d-17
  td_E_eV(9)  =   20.0_8; td_cs_cm2(9)  = 2.70d-17
  td_E_eV(10) =   20.5_8; td_cs_cm2(10) = 3.08d-17
  td_E_eV(11) =   21.0_8; td_cs_cm2(11) = 3.44d-17
  td_E_eV(12) =   21.5_8; td_cs_cm2(12) = 3.80d-17
  td_E_eV(13) =   22.0_8; td_cs_cm2(13) = 4.18d-17
  td_E_eV(14) =   22.5_8; td_cs_cm2(14) = 4.55d-17
  td_E_eV(15) =   23.0_8; td_cs_cm2(15) = 4.92d-17
  td_E_eV(16) =   23.5_8; td_cs_cm2(16) = 5.28d-17
  td_E_eV(17) =   24.0_8; td_cs_cm2(17) = 5.65d-17
  td_E_eV(18) =   24.5_8; td_cs_cm2(18) = 6.03d-17
  td_E_eV(19) =   25.0_8; td_cs_cm2(19) = 6.40d-17
  td_E_eV(20) =   30.0_8; td_cs_cm2(20) = 9.29d-17
  td_E_eV(21) =   35.0_8; td_cs_cm2(21) = 1.16d-16
  td_E_eV(22) =   40.0_8; td_cs_cm2(22) = 1.37d-16
! table 16 of [Itikawa]
  td_E_eV(23) =   45.0_8; td_cs_cm2(23) = 1.52d-16
  td_E_eV(24) =   50.0_8; td_cs_cm2(24) = 1.60d-16
  td_E_eV(25) =   55.0_8; td_cs_cm2(25) = 1.66d-16
  td_E_eV(26) =   60.0_8; td_cs_cm2(26) = 1.72d-16
  td_E_eV(27) =   65.0_8; td_cs_cm2(27) = 1.74d-16
  td_E_eV(28) =   70.0_8; td_cs_cm2(28) = 1.78d-16
  td_E_eV(29) =   75.0_8; td_cs_cm2(29) = 1.80d-16
  td_E_eV(30) =   80.0_8; td_cs_cm2(30) = 1.81d-16
  td_E_eV(31) =   85.0_8; td_cs_cm2(31) = 1.82d-16
  td_E_eV(32) =   90.0_8; td_cs_cm2(32) = 1.83d-16
  td_E_eV(33) =   95.0_8; td_cs_cm2(33) = 1.85d-16
  td_E_eV(34) =  100.0_8; td_cs_cm2(34) = 1.85d-16
  td_E_eV(35) =  110.0_8; td_cs_cm2(35) = 1.83d-16
  td_E_eV(36) =  120.0_8; td_cs_cm2(36) = 1.81d-16
  td_E_eV(37) =  140.0_8; td_cs_cm2(37) = 1.78d-16
  td_E_eV(38) =  160.0_8; td_cs_cm2(38) = 1.72d-16
  td_E_eV(39) =  180.0_8; td_cs_cm2(39) = 1.67d-16
  td_E_eV(40) =  200.0_8; td_cs_cm2(40) = 1.61d-16
  td_E_eV(41) =  225.0_8; td_cs_cm2(41) = 1.55d-16
  td_E_eV(42) =  250.0_8; td_cs_cm2(42) = 1.48d-16
  td_E_eV(43) =  275.0_8; td_cs_cm2(43) = 1.41d-16
  td_E_eV(44) =  300.0_8; td_cs_cm2(44) = 1.37d-16
! table 17 of [Itikawa]
  td_E_eV(45) =  350.0_8; td_cs_cm2(45) = 1.28d-16
  td_E_eV(46) =  400.0_8; td_cs_cm2(46) = 1.20d-16
  td_E_eV(47) =  450.0_8; td_cs_cm2(47) = 1.11d-16
  td_E_eV(48) =  500.0_8; td_cs_cm2(48) = 1.05d-16
  td_E_eV(49) =  550.0_8; td_cs_cm2(49) = 9.98d-17
  td_E_eV(50) =  600.0_8; td_cs_cm2(50) = 9.43d-17
  td_E_eV(51) =  650.0_8; td_cs_cm2(51) = 8.80d-17
  td_E_eV(52) =  700.0_8; td_cs_cm2(52) = 8.44d-17
  td_E_eV(53) =  750.0_8; td_cs_cm2(53) = 7.96d-17
  td_E_eV(54) =  800.0_8; td_cs_cm2(54) = 7.65d-17
  td_E_eV(55) =  850.0_8; td_cs_cm2(55) = 7.38d-17
  td_E_eV(56) =  900.0_8; td_cs_cm2(56) = 7.19d-17
  td_E_eV(57) =  950.0_8; td_cs_cm2(57) = 6.98d-17
  td_E_eV(58) = 1000.0_8; td_cs_cm2(58) = 6.76d-17

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_ionN2
!
!-------------------------------------
!
real(8) function CSV_N2_ionN2_m3s(E_eV)

  use N2_ionN2
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  real(8) CS_N2_ionN2_cm2  ! cross section, used outside of the table data range
  

  IF (E_eV.LT.td_E_eV(1)) THEN
     CSV_N2_ionN2_m3s = 0.0_8
     RETURN
  END IF

  IF (E_eV.GT.td_E_eV(N_td)) THEN
     CS_N2_ionN2_cm2 = (2.14d-14 * log(E_eV) - 8.02d-14) / E_eV
     CSV_N2_ionN2_m3s = max(0.0_8, 1.0d-4 * CS_N2_ionN2_cm2 * sqrt(E_eV) * efactor_eV_to_ms)
     RETURN
  END IF

  DO n = 1, N_td-1
     IF ((E_eV.GE.td_E_eV(n)).AND.(E_eV.LE.td_E_eV(n+1))) THEN
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_ionN2_m3s = MAX(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        RETURN
     END IF
  END DO

end function CSV_N2_ionN2_m3s

!---------------------------------------------------------------------------------------------------- #3
! N2, {N2->N+} ionization cross-section 
! Tables 15,16,17 of [Y.Itikawa, J.Phys.Chem.Ref.Data, v.35, p.31-53, 2006] 
!
module N2_ionN
  integer, parameter :: N_td=39  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_ionN
!
!-------------------------------------
!
subroutine Prepare_N2_ionN

  use N2_ionN
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

! table 15 of [Itikawa]
  td_E_eV(1)  =   30.0_8; td_cs_cm2( 1) = 3.25d-18
  td_E_eV(2)  =   35.0_8; td_cs_cm2( 2) = 9.04d-18
  td_E_eV(3)  =   40.0_8; td_cs_cm2( 3) = 1.66d-17
! table 16 of [Itikawa]
  td_E_eV(4)  =   45.0_8; td_cs_cm2( 4) = 2.45d-17
  td_E_eV(5)  =   50.0_8; td_cs_cm2( 5) = 3.19d-17
  td_E_eV(6)  =   55.0_8; td_cs_cm2( 6) = 3.90d-17
  td_E_eV(7)  =   60.0_8; td_cs_cm2( 7) = 4.38d-17
  td_E_eV(8)  =   65.0_8; td_cs_cm2( 8) = 4.82d-17
  td_E_eV(9)  =   70.0_8; td_cs_cm2( 9) = 5.23d-17
  td_E_eV(10) =   75.0_8; td_cs_cm2(10) = 5.61d-17
  td_E_eV(11) =   80.0_8; td_cs_cm2(11) = 5.87d-17
  td_E_eV(12) =   85.0_8; td_cs_cm2(12) = 6.05d-17
  td_E_eV(13) =   90.0_8; td_cs_cm2(13) = 6.32d-17
  td_E_eV(14) =   95.0_8; td_cs_cm2(14) = 6.45d-17
  td_E_eV(15) =  100.0_8; td_cs_cm2(15) = 6.56d-17
  td_E_eV(16) =  110.0_8; td_cs_cm2(16) = 6.60d-17
  td_E_eV(17) =  120.0_8; td_cs_cm2(17) = 6.61d-17
  td_E_eV(18) =  140.0_8; td_cs_cm2(18) = 6.52d-17
  td_E_eV(19) =  160.0_8; td_cs_cm2(19) = 6.33d-17
  td_E_eV(20) =  180.0_8; td_cs_cm2(20) = 5.95d-17
  td_E_eV(21) =  200.0_8; td_cs_cm2(21) = 5.66d-17
  td_E_eV(22) =  225.0_8; td_cs_cm2(22) = 5.16d-17
  td_E_eV(23) =  250.0_8; td_cs_cm2(23) = 4.93d-17
  td_E_eV(24) =  275.0_8; td_cs_cm2(24) = 4.58d-17
  td_E_eV(25) =  300.0_8; td_cs_cm2(25) = 4.38d-17
! table 17 of [Itikawa]
  td_E_eV(26) =  350.0_8; td_cs_cm2(26) = 3.93d-17
  td_E_eV(27) =  400.0_8; td_cs_cm2(27) = 3.51d-17
  td_E_eV(28) =  450.0_8; td_cs_cm2(28) = 3.24d-17
  td_E_eV(29) =  500.0_8; td_cs_cm2(29) = 2.99d-17
  td_E_eV(30) =  550.0_8; td_cs_cm2(30) = 2.74d-17
  td_E_eV(31) =  600.0_8; td_cs_cm2(31) = 2.48d-17
  td_E_eV(32) =  650.0_8; td_cs_cm2(32) = 2.34d-17
  td_E_eV(33) =  700.0_8; td_cs_cm2(33) = 2.17d-17
  td_E_eV(34) =  750.0_8; td_cs_cm2(34) = 2.05d-17
  td_E_eV(35) =  800.0_8; td_cs_cm2(35) = 2.00d-17
  td_E_eV(36) =  850.0_8; td_cs_cm2(36) = 1.92d-17
  td_E_eV(37) =  900.0_8; td_cs_cm2(37) = 1.83d-17
  td_E_eV(38) =  950.0_8; td_cs_cm2(38) = 1.76d-17
  td_E_eV(39) = 1000.0_8; td_cs_cm2(39) = 1.67d-17

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_ionN
!
!-------------------------------------
!
real(8) function CSV_N2_ionN_m3s(E_eV)

  use N2_ionN
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  real(8) CS_N2_ionN_cm2  ! cross section, used outside of the table data range
  
  IF (E_eV.LT.td_E_eV(1)) THEN
     CSV_N2_ionN_m3s = 0.0_8
     RETURN
  END IF

  IF (E_eV.GT.td_E_eV(N_td)) THEN
     CS_N2_ionN_cm2 = (2.9d-15 * log(E_eV) - 3.35d-15) / E_eV
     CSV_N2_ionN_m3s = max(0.0_8, 1.0d-4 * CS_N2_ionN_cm2 * sqrt(E_eV) * efactor_eV_to_ms)
     RETURN
  END IF

  DO n = 1, N_td-1
     IF ((E_eV.GE.td_E_eV(n)).AND.(E_eV.LE.td_E_eV(n+1))) THEN
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_ionN_m3s = MAX(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        RETURN
     END IF
  END DO

end function CSV_N2_ionN_m3s







!---------------------------------------------------------------------------------------------------- #4
!#VIBRATIONAL
!#SPECIES: e / N2
!#PROCESS: E + N2 <-> E + N2*(v=0-1), Vibrational
!#PARAM.:  E = 0.288 eV, g1/g0 = 1, complete set
!#COMMENT: Vibrational excitation : N2(1Σg+:v=0) -> N2*(v=1) | Source: resonant CS from [Laporta
!#COMMENT: et al, Plasma Sources Sci. Technol. 23, 065002 (2014)] and outside of resonance: [Sohn et
!#COMMENT: al. J. Phys. B: At. Mol. Phys. 19, 4017 (1986)] corrected, [LinertampZubek, J. Phys. B:
!#COMMENT: At. Mol. Opt. Phys. 42, 085203 (2009)] and [Tanaka et al, J. Phys. B: At. Mol. Phys. 14,
!#COMMENT: 2081 (1981)].
!!#UPDATED: 2023-12-23 21:19:17
!
module N2_vib01
  integer, parameter :: N_td=738  !8   ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_vib01
!
!-------------------------------------
!
subroutine Prepare_N2_vib01

  use N2_vib01
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

!  td_E_eV(1) = 1.3_8; td_cs_cm2(1) = 1.11d-18
!  td_E_eV(2) = 1.6_8; td_cs_cm2(2) = 8.84d-17
!  td_E_eV(3) = 2.0_8; td_cs_cm2(3) = 1.01d-15
!  td_E_eV(4) = 2.4_8; td_cs_cm2(4) = 1.35d-15
!  td_E_eV(5) = 3.0_8; td_cs_cm2(5) = 5.86d-16
!  td_E_eV(6) = 4.0_8; td_cs_cm2(6) = 1.37d-16
!  td_E_eV(7) = 5.0_8; td_cs_cm2(7) = 2.34d-17
!  td_E_eV(8) = 6.0_8; td_cs_cm2(8) = 3.66d-18

  td_E_eV(  1) = 2.880000d-1 ; td_cs_cm2(  1) = 0.000000d+0
  td_E_eV(  2) = 3.000000d-1 ; td_cs_cm2(  2) = 4.622222d-24
  td_E_eV(  3) = 3.200000d-1 ; td_cs_cm2(  3) = 9.244444d-24
  td_E_eV(  4) = 3.400000d-1 ; td_cs_cm2(  4) = 1.386667d-23
  td_E_eV(  5) = 3.600000d-1 ; td_cs_cm2(  5) = 1.848889d-23
  td_E_eV(  6) = 3.800000d-1 ; td_cs_cm2(  6) = 2.311111d-23
  td_E_eV(  7) = 4.000000d-1 ; td_cs_cm2(  7) = 2.773333d-23
  td_E_eV(  8) = 4.200000d-1 ; td_cs_cm2(  8) = 3.235556d-23
  td_E_eV(  9) = 4.400000d-1 ; td_cs_cm2(  9) = 3.697778d-23
  td_E_eV( 10) = 4.600000d-1 ; td_cs_cm2( 10) = 4.160000d-23
  td_E_eV( 11) = 4.800000d-1 ; td_cs_cm2( 11) = 4.622222d-23
  td_E_eV( 12) = 5.000000d-1 ; td_cs_cm2( 12) = 5.084444d-23
  td_E_eV( 13) = 5.200000d-1 ; td_cs_cm2( 13) = 5.546667d-23
  td_E_eV( 14) = 5.400000d-1 ; td_cs_cm2( 14) = 6.008889d-23
  td_E_eV( 15) = 5.600000d-1 ; td_cs_cm2( 15) = 6.379872d-23
  td_E_eV( 16) = 5.800000d-1 ; td_cs_cm2( 16) = 6.661336d-23
  td_E_eV( 17) = 6.000000d-1 ; td_cs_cm2( 17) = 6.945045d-23
  td_E_eV( 18) = 6.200000d-1 ; td_cs_cm2( 18) = 7.230941d-23
  td_E_eV( 19) = 6.400000d-1 ; td_cs_cm2( 19) = 7.518968d-23
  td_E_eV( 20) = 6.600000d-1 ; td_cs_cm2( 20) = 7.809076d-23
  td_E_eV( 21) = 6.800000d-1 ; td_cs_cm2( 21) = 8.101215d-23
  td_E_eV( 22) = 7.000000d-1 ; td_cs_cm2( 22) = 8.395341d-23
  td_E_eV( 23) = 7.200000d-1 ; td_cs_cm2( 23) = 8.691408d-23
  td_E_eV( 24) = 7.400000d-1 ; td_cs_cm2( 24) = 8.989375d-23
  td_E_eV( 25) = 7.600000d-1 ; td_cs_cm2( 25) = 9.289202d-23
  td_E_eV( 26) = 7.800000d-1 ; td_cs_cm2( 26) = 9.590853d-23
  td_E_eV( 27) = 8.000000d-1 ; td_cs_cm2( 27) = 9.894289d-23
  td_E_eV( 28) = 8.200000d-1 ; td_cs_cm2( 28) = 1.019948d-22
  td_E_eV( 29) = 8.400000d-1 ; td_cs_cm2( 29) = 1.050639d-22
  td_E_eV( 30) = 8.600000d-1 ; td_cs_cm2( 30) = 1.081498d-22
  td_E_eV( 31) = 8.800000d-1 ; td_cs_cm2( 31) = 1.112523d-22
  td_E_eV( 32) = 9.000000d-1 ; td_cs_cm2( 32) = 1.143712d-22
  td_E_eV( 33) = 9.200000d-1 ; td_cs_cm2( 33) = 1.175060d-22
  td_E_eV( 34) = 9.400000d-1 ; td_cs_cm2( 34) = 1.206565d-22
  td_E_eV( 35) = 9.600000d-1 ; td_cs_cm2( 35) = 1.238225d-22
  td_E_eV( 36) = 9.800000d-1 ; td_cs_cm2( 36) = 1.270038d-22
  td_E_eV( 37) = 1.000000d+0 ; td_cs_cm2( 37) = 1.302000d-22
  td_E_eV( 38) = 1.020000d+0 ; td_cs_cm2( 38) = 1.420963d-22
  td_E_eV( 39) = 1.040000d+0 ; td_cs_cm2( 39) = 1.548165d-22
  td_E_eV( 40) = 1.060000d+0 ; td_cs_cm2( 40) = 1.684001d-22
  td_E_eV( 41) = 1.080000d+0 ; td_cs_cm2( 41) = 1.828877d-22
  td_E_eV( 42) = 1.100000d+0 ; td_cs_cm2( 42) = 1.983213d-22
  td_E_eV( 43) = 1.120000d+0 ; td_cs_cm2( 43) = 2.147435d-22
  td_E_eV( 44) = 1.140000d+0 ; td_cs_cm2( 44) = 2.321983d-22
  td_E_eV( 45) = 1.160000d+0 ; td_cs_cm2( 45) = 2.507309d-22
  td_E_eV( 46) = 1.180000d+0 ; td_cs_cm2( 46) = 2.703875d-22
  td_E_eV( 47) = 1.200000d+0 ; td_cs_cm2( 47) = 2.912155d-22
  td_E_eV( 48) = 1.220000d+0 ; td_cs_cm2( 48) = 3.132634d-22
  td_E_eV( 49) = 1.240000d+0 ; td_cs_cm2( 49) = 3.365808d-22
  td_E_eV( 50) = 1.260000d+0 ; td_cs_cm2( 50) = 3.612186d-22
  td_E_eV( 51) = 1.280000d+0 ; td_cs_cm2( 51) = 3.872289d-22
  td_E_eV( 52) = 1.300000d+0 ; td_cs_cm2( 52) = 4.146649d-22
  td_E_eV( 53) = 1.320000d+0 ; td_cs_cm2( 53) = 4.435808d-22
  td_E_eV( 54) = 1.340000d+0 ; td_cs_cm2( 54) = 4.740324d-22
  td_E_eV( 55) = 1.360000d+0 ; td_cs_cm2( 55) = 5.060765d-22
  td_E_eV( 56) = 1.380000d+0 ; td_cs_cm2( 56) = 5.397709d-22
  td_E_eV( 57) = 1.400000d+0 ; td_cs_cm2( 57) = 5.751751d-22
  td_E_eV( 58) = 1.420000d+0 ; td_cs_cm2( 58) = 6.123493d-22
  td_E_eV( 59) = 1.440000d+0 ; td_cs_cm2( 59) = 6.513554d-22
  td_E_eV( 60) = 1.460000d+0 ; td_cs_cm2( 60) = 6.922562d-22
  td_E_eV( 61) = 1.480000d+0 ; td_cs_cm2( 61) = 7.351159d-22
  td_E_eV( 62) = 1.500000d+0 ; td_cs_cm2( 62) = 7.800000d-22
  td_E_eV( 63) = 1.520000d+0 ; td_cs_cm2( 63) = 8.370184d-22
  td_E_eV( 64) = 1.540000d+0 ; td_cs_cm2( 64) = 8.973769d-22
  td_E_eV( 65) = 1.560000d+0 ; td_cs_cm2( 65) = 9.612239d-22
  td_E_eV( 66) = 1.580000d+0 ; td_cs_cm2( 66) = 1.028712d-21
  td_E_eV( 67) = 1.600000d+0 ; td_cs_cm2( 67) = 1.100000d-21
  td_E_eV( 68) = 1.620000d+0 ; td_cs_cm2( 68) = 1.353300d-21
  td_E_eV( 69) = 1.640000d+0 ; td_cs_cm2( 69) = 1.673500d-21
  td_E_eV( 70) = 1.660000d+0 ; td_cs_cm2( 70) = 2.077100d-21
  td_E_eV( 71) = 1.680000d+0 ; td_cs_cm2( 71) = 2.589800d-21
  td_E_eV( 72) = 1.700000d+0 ; td_cs_cm2( 72) = 3.245900d-21
  td_E_eV( 73) = 1.720000d+0 ; td_cs_cm2( 73) = 4.092900d-21
  td_E_eV( 74) = 1.740000d+0 ; td_cs_cm2( 74) = 5.196300d-21
  td_E_eV( 75) = 1.760000d+0 ; td_cs_cm2( 75) = 6.647500d-21
  td_E_eV( 76) = 1.780000d+0 ; td_cs_cm2( 76) = 8.574300d-21
  td_E_eV( 77) = 1.800000d+0 ; td_cs_cm2( 77) = 1.115400d-20
  td_E_eV( 78) = 1.820000d+0 ; td_cs_cm2( 78) = 1.462800d-20
  td_E_eV( 79) = 1.840000d+0 ; td_cs_cm2( 79) = 1.930200d-20
  td_E_eV( 80) = 1.860000d+0 ; td_cs_cm2( 80) = 2.551300d-20
  td_E_eV( 81) = 1.880000d+0 ; td_cs_cm2( 81) = 3.347000d-20
  td_E_eV( 82) = 1.900000d+0 ; td_cs_cm2( 82) = 4.287900d-20
  td_E_eV( 83) = 1.920000d+0 ; td_cs_cm2( 83) = 5.234900d-20
  td_E_eV( 84) = 1.940000d+0 ; td_cs_cm2( 84) = 5.918400d-20
  td_E_eV( 85) = 1.960000d+0 ; td_cs_cm2( 85) = 6.067600d-20
  td_E_eV( 86) = 1.980000d+0 ; td_cs_cm2( 86) = 5.634500d-20
  td_E_eV( 87) = 2.000000d+0 ; td_cs_cm2( 87) = 4.823800d-20
  td_E_eV( 88) = 2.020000d+0 ; td_cs_cm2( 88) = 3.898000d-20
  td_E_eV( 89) = 2.040000d+0 ; td_cs_cm2( 89) = 3.030900d-20
  td_E_eV( 90) = 2.060000d+0 ; td_cs_cm2( 90) = 2.302900d-20
  td_E_eV( 91) = 2.080000d+0 ; td_cs_cm2( 91) = 1.751600d-20
  td_E_eV( 92) = 2.100000d+0 ; td_cs_cm2( 92) = 1.416500d-20
  td_E_eV( 93) = 2.120000d+0 ; td_cs_cm2( 93) = 1.362500d-20
  td_E_eV( 94) = 2.140000d+0 ; td_cs_cm2( 94) = 1.672000d-20
  td_E_eV( 95) = 2.160000d+0 ; td_cs_cm2( 95) = 2.383000d-20
  td_E_eV( 96) = 2.180000d+0 ; td_cs_cm2( 96) = 3.389700d-20
  td_E_eV( 97) = 2.200000d+0 ; td_cs_cm2( 97) = 4.430100d-20
  td_E_eV( 98) = 2.220000d+0 ; td_cs_cm2( 98) = 5.245800d-20
  td_E_eV( 99) = 2.240000d+0 ; td_cs_cm2( 99) = 5.735100d-20
  td_E_eV(100) = 2.260000d+0 ; td_cs_cm2(100) = 5.933800d-20
  td_E_eV(101) = 2.280000d+0 ; td_cs_cm2(101) = 5.916500d-20
  td_E_eV(102) = 2.300000d+0 ; td_cs_cm2(102) = 5.734300d-20
  td_E_eV(103) = 2.320000d+0 ; td_cs_cm2(103) = 5.399700d-20
  td_E_eV(104) = 2.340000d+0 ; td_cs_cm2(104) = 4.895600d-20
  td_E_eV(105) = 2.360000d+0 ; td_cs_cm2(105) = 4.205300d-20
  td_E_eV(106) = 2.380000d+0 ; td_cs_cm2(106) = 3.367200d-20
  td_E_eV(107) = 2.400000d+0 ; td_cs_cm2(107) = 2.527900d-20
  td_E_eV(108) = 2.420000d+0 ; td_cs_cm2(108) = 1.909300d-20
  td_E_eV(109) = 2.440000d+0 ; td_cs_cm2(109) = 1.667400d-20
  td_E_eV(110) = 2.460000d+0 ; td_cs_cm2(110) = 1.801900d-20
  td_E_eV(111) = 2.480000d+0 ; td_cs_cm2(111) = 2.209700d-20
  td_E_eV(112) = 2.500000d+0 ; td_cs_cm2(112) = 2.778000d-20
  td_E_eV(113) = 2.520000d+0 ; td_cs_cm2(113) = 3.420500d-20
  td_E_eV(114) = 2.540000d+0 ; td_cs_cm2(114) = 4.065600d-20
  td_E_eV(115) = 2.560000d+0 ; td_cs_cm2(115) = 4.627000d-20
  td_E_eV(116) = 2.580000d+0 ; td_cs_cm2(116) = 4.978800d-20
  td_E_eV(117) = 2.600000d+0 ; td_cs_cm2(117) = 4.966900d-20
  td_E_eV(118) = 2.620000d+0 ; td_cs_cm2(118) = 4.494100d-20
  td_E_eV(119) = 2.640000d+0 ; td_cs_cm2(119) = 3.642000d-20
  td_E_eV(120) = 2.660000d+0 ; td_cs_cm2(120) = 2.672900d-20
  td_E_eV(121) = 2.680000d+0 ; td_cs_cm2(121) = 1.861300d-20
  td_E_eV(122) = 2.700000d+0 ; td_cs_cm2(122) = 1.351900d-20
  td_E_eV(123) = 2.720000d+0 ; td_cs_cm2(123) = 1.163100d-20
  td_E_eV(124) = 2.740000d+0 ; td_cs_cm2(124) = 1.254500d-20
  td_E_eV(125) = 2.760000d+0 ; td_cs_cm2(125) = 1.570700d-20
  td_E_eV(126) = 2.780000d+0 ; td_cs_cm2(126) = 2.046300d-20
  td_E_eV(127) = 2.800000d+0 ; td_cs_cm2(127) = 2.585900d-20
  td_E_eV(128) = 2.820000d+0 ; td_cs_cm2(128) = 3.044300d-20
  td_E_eV(129) = 2.840000d+0 ; td_cs_cm2(129) = 3.252900d-20
  td_E_eV(130) = 2.860000d+0 ; td_cs_cm2(130) = 3.114700d-20
  td_E_eV(131) = 2.880000d+0 ; td_cs_cm2(131) = 2.684600d-20
  td_E_eV(132) = 2.900000d+0 ; td_cs_cm2(132) = 2.126000d-20
  td_E_eV(133) = 2.920000d+0 ; td_cs_cm2(133) = 1.597200d-20
  td_E_eV(134) = 2.940000d+0 ; td_cs_cm2(134) = 1.191400d-20
  td_E_eV(135) = 2.960000d+0 ; td_cs_cm2(135) = 9.453100d-21
  td_E_eV(136) = 2.980000d+0 ; td_cs_cm2(136) = 8.660600d-21
  td_E_eV(137) = 3.000000d+0 ; td_cs_cm2(137) = 9.445100d-21
  td_E_eV(138) = 3.020000d+0 ; td_cs_cm2(138) = 1.150700d-20
  td_E_eV(139) = 3.040000d+0 ; td_cs_cm2(139) = 1.419300d-20
  td_E_eV(140) = 3.060000d+0 ; td_cs_cm2(140) = 1.649500d-20
  td_E_eV(141) = 3.080000d+0 ; td_cs_cm2(141) = 1.747400d-20
  td_E_eV(142) = 3.100000d+0 ; td_cs_cm2(142) = 1.682400d-20
  td_E_eV(143) = 3.120000d+0 ; td_cs_cm2(143) = 1.495100d-20
  td_E_eV(144) = 3.140000d+0 ; td_cs_cm2(144) = 1.254500d-20
  td_E_eV(145) = 3.160000d+0 ; td_cs_cm2(145) = 1.018800d-20
  td_E_eV(146) = 3.180000d+0 ; td_cs_cm2(146) = 8.255700d-21
  td_E_eV(147) = 3.200000d+0 ; td_cs_cm2(147) = 6.966800d-21
  td_E_eV(148) = 3.220000d+0 ; td_cs_cm2(148) = 6.442600d-21
  td_E_eV(149) = 3.240000d+0 ; td_cs_cm2(149) = 6.686300d-21
  td_E_eV(150) = 3.260000d+0 ; td_cs_cm2(150) = 7.502300d-21
  td_E_eV(151) = 3.280000d+0 ; td_cs_cm2(151) = 8.468600d-21
  td_E_eV(152) = 3.300000d+0 ; td_cs_cm2(152) = 9.107200d-21
  td_E_eV(153) = 3.320000d+0 ; td_cs_cm2(153) = 9.154000d-21
  td_E_eV(154) = 3.340000d+0 ; td_cs_cm2(154) = 8.646800d-21
  td_E_eV(155) = 3.360000d+0 ; td_cs_cm2(155) = 7.791200d-21
  td_E_eV(156) = 3.380000d+0 ; td_cs_cm2(156) = 6.809700d-21
  td_E_eV(157) = 3.400000d+0 ; td_cs_cm2(157) = 5.881000d-21
  td_E_eV(158) = 3.420000d+0 ; td_cs_cm2(158) = 5.140600d-21
  td_E_eV(159) = 3.440000d+0 ; td_cs_cm2(159) = 4.688200d-21
  td_E_eV(160) = 3.460000d+0 ; td_cs_cm2(160) = 4.566400d-21
  td_E_eV(161) = 3.480000d+0 ; td_cs_cm2(161) = 4.716300d-21
  td_E_eV(162) = 3.500000d+0 ; td_cs_cm2(162) = 4.973000d-21
  td_E_eV(163) = 3.520000d+0 ; td_cs_cm2(163) = 5.148700d-21
  td_E_eV(164) = 3.540000d+0 ; td_cs_cm2(164) = 5.137700d-21
  td_E_eV(165) = 3.560000d+0 ; td_cs_cm2(165) = 4.938000d-21
  td_E_eV(166) = 3.580000d+0 ; td_cs_cm2(166) = 4.605600d-21
  td_E_eV(167) = 3.600000d+0 ; td_cs_cm2(167) = 4.209700d-21
  td_E_eV(168) = 3.620000d+0 ; td_cs_cm2(168) = 3.814400d-21
  td_E_eV(169) = 3.640000d+0 ; td_cs_cm2(169) = 3.477000d-21
  td_E_eV(170) = 3.660000d+0 ; td_cs_cm2(170) = 3.243600d-21
  td_E_eV(171) = 3.680000d+0 ; td_cs_cm2(171) = 3.133100d-21
  td_E_eV(172) = 3.700000d+0 ; td_cs_cm2(172) = 3.120100d-21
  td_E_eV(173) = 3.720000d+0 ; td_cs_cm2(173) = 3.142700d-21
  td_E_eV(174) = 3.740000d+0 ; td_cs_cm2(174) = 3.140300d-21
  td_E_eV(175) = 3.760000d+0 ; td_cs_cm2(175) = 3.082300d-21
  td_E_eV(176) = 3.780000d+0 ; td_cs_cm2(176) = 2.968300d-21
  td_E_eV(177) = 3.800000d+0 ; td_cs_cm2(177) = 2.813500d-21
  td_E_eV(178) = 3.820000d+0 ; td_cs_cm2(178) = 2.638100d-21
  td_E_eV(179) = 3.840000d+0 ; td_cs_cm2(179) = 2.464200d-21
  td_E_eV(180) = 3.860000d+0 ; td_cs_cm2(180) = 2.313500d-21
  td_E_eV(181) = 3.880000d+0 ; td_cs_cm2(181) = 2.202800d-21
  td_E_eV(182) = 3.900000d+0 ; td_cs_cm2(182) = 2.135000d-21
  td_E_eV(183) = 3.920000d+0 ; td_cs_cm2(183) = 2.096200d-21
  td_E_eV(184) = 3.940000d+0 ; td_cs_cm2(184) = 2.064600d-21
  td_E_eV(185) = 3.960000d+0 ; td_cs_cm2(185) = 2.023100d-21
  td_E_eV(186) = 3.980000d+0 ; td_cs_cm2(186) = 1.965300d-21
  td_E_eV(187) = 4.000000d+0 ; td_cs_cm2(187) = 1.892000d-21
  td_E_eV(188) = 4.020000d+0 ; td_cs_cm2(188) = 1.808300d-21
  td_E_eV(189) = 4.040000d+0 ; td_cs_cm2(189) = 1.721100d-21
  td_E_eV(190) = 4.060000d+0 ; td_cs_cm2(190) = 1.638600d-21
  td_E_eV(191) = 4.080000d+0 ; td_cs_cm2(191) = 1.568300d-21
  td_E_eV(192) = 4.100000d+0 ; td_cs_cm2(192) = 1.513600d-21
  td_E_eV(193) = 4.120000d+0 ; td_cs_cm2(193) = 1.472000d-21
  td_E_eV(194) = 4.140000d+0 ; td_cs_cm2(194) = 1.436500d-21
  td_E_eV(195) = 4.160000d+0 ; td_cs_cm2(195) = 1.401100d-21
  td_E_eV(196) = 4.180000d+0 ; td_cs_cm2(196) = 1.362100d-21
  td_E_eV(197) = 4.200000d+0 ; td_cs_cm2(197) = 1.318900d-21
  td_E_eV(198) = 4.220000d+0 ; td_cs_cm2(198) = 1.272300d-21
  td_E_eV(199) = 4.240000d+0 ; td_cs_cm2(199) = 1.224400d-21
  td_E_eV(200) = 4.260000d+0 ; td_cs_cm2(200) = 1.177900d-21
  td_E_eV(201) = 4.280000d+0 ; td_cs_cm2(201) = 1.135500d-21
  td_E_eV(202) = 4.300000d+0 ; td_cs_cm2(202) = 1.098800d-21
  td_E_eV(203) = 4.320000d+0 ; td_cs_cm2(203) = 1.067300d-21
  td_E_eV(204) = 4.340000d+0 ; td_cs_cm2(204) = 1.038900d-21
  td_E_eV(205) = 4.360000d+0 ; td_cs_cm2(205) = 1.011600d-21
  td_E_eV(206) = 4.380000d+0 ; td_cs_cm2(206) = 9.839300d-22
  td_E_eV(207) = 4.400000d+0 ; td_cs_cm2(207) = 9.554900d-22
  td_E_eV(208) = 4.420000d+0 ; td_cs_cm2(208) = 9.263700d-22
  td_E_eV(209) = 4.440000d+0 ; td_cs_cm2(209) = 8.971600d-22
  td_E_eV(210) = 4.460000d+0 ; td_cs_cm2(210) = 8.688100d-22
  td_E_eV(211) = 4.480000d+0 ; td_cs_cm2(211) = 8.423000d-22
  td_E_eV(212) = 4.500000d+0 ; td_cs_cm2(212) = 8.160000d-22
  td_E_eV(213) = 4.520000d+0 ; td_cs_cm2(213) = 8.135934d-22
  td_E_eV(214) = 4.540000d+0 ; td_cs_cm2(214) = 8.112044d-22
  td_E_eV(215) = 4.560000d+0 ; td_cs_cm2(215) = 8.088329d-22
  td_E_eV(216) = 4.580000d+0 ; td_cs_cm2(216) = 8.064787d-22
  td_E_eV(217) = 4.600000d+0 ; td_cs_cm2(217) = 8.041416d-22
  td_E_eV(218) = 4.620000d+0 ; td_cs_cm2(218) = 8.018213d-22
  td_E_eV(219) = 4.640000d+0 ; td_cs_cm2(219) = 7.995177d-22
  td_E_eV(220) = 4.660000d+0 ; td_cs_cm2(220) = 7.972305d-22
  td_E_eV(221) = 4.680000d+0 ; td_cs_cm2(221) = 7.949597d-22
  td_E_eV(222) = 4.700000d+0 ; td_cs_cm2(222) = 7.927050d-22
  td_E_eV(223) = 4.720000d+0 ; td_cs_cm2(223) = 7.904662d-22
  td_E_eV(224) = 4.740000d+0 ; td_cs_cm2(224) = 7.882432d-22
  td_E_eV(225) = 4.760000d+0 ; td_cs_cm2(225) = 7.860357d-22
  td_E_eV(226) = 4.780000d+0 ; td_cs_cm2(226) = 7.838436d-22
  td_E_eV(227) = 4.800000d+0 ; td_cs_cm2(227) = 7.816668d-22
  td_E_eV(228) = 4.820000d+0 ; td_cs_cm2(228) = 7.795050d-22
  td_E_eV(229) = 4.840000d+0 ; td_cs_cm2(229) = 7.773581d-22
  td_E_eV(230) = 4.860000d+0 ; td_cs_cm2(230) = 7.752259d-22
  td_E_eV(231) = 4.880000d+0 ; td_cs_cm2(231) = 7.731083d-22
  td_E_eV(232) = 4.900000d+0 ; td_cs_cm2(232) = 7.710052d-22
  td_E_eV(233) = 4.920000d+0 ; td_cs_cm2(233) = 7.689162d-22
  td_E_eV(234) = 4.940000d+0 ; td_cs_cm2(234) = 7.668414d-22
  td_E_eV(235) = 4.960000d+0 ; td_cs_cm2(235) = 7.647805d-22
  td_E_eV(236) = 4.980000d+0 ; td_cs_cm2(236) = 7.627334d-22
  td_E_eV(237) = 5.000000d+0 ; td_cs_cm2(237) = 7.607000d-22
  td_E_eV(238) = 5.020000d+0 ; td_cs_cm2(238) = 7.538975d-22
  td_E_eV(239) = 5.040000d+0 ; td_cs_cm2(239) = 7.471825d-22
  td_E_eV(240) = 5.060000d+0 ; td_cs_cm2(240) = 7.405535d-22
  td_E_eV(241) = 5.080000d+0 ; td_cs_cm2(241) = 7.340092d-22
  td_E_eV(242) = 5.100000d+0 ; td_cs_cm2(242) = 7.275480d-22
  td_E_eV(243) = 5.120000d+0 ; td_cs_cm2(243) = 7.211687d-22
  td_E_eV(244) = 5.140000d+0 ; td_cs_cm2(244) = 7.148699d-22
  td_E_eV(245) = 5.160000d+0 ; td_cs_cm2(245) = 7.086502d-22
  td_E_eV(246) = 5.180000d+0 ; td_cs_cm2(246) = 7.025084d-22
  td_E_eV(247) = 5.200000d+0 ; td_cs_cm2(247) = 6.964432d-22
  td_E_eV(248) = 5.220000d+0 ; td_cs_cm2(248) = 6.904534d-22
  td_E_eV(249) = 5.240000d+0 ; td_cs_cm2(249) = 6.845376d-22
  td_E_eV(250) = 5.260000d+0 ; td_cs_cm2(250) = 6.786948d-22
  td_E_eV(251) = 5.280000d+0 ; td_cs_cm2(251) = 6.729238d-22
  td_E_eV(252) = 5.300000d+0 ; td_cs_cm2(252) = 6.672233d-22
  td_E_eV(253) = 5.320000d+0 ; td_cs_cm2(253) = 6.615924d-22
  td_E_eV(254) = 5.340000d+0 ; td_cs_cm2(254) = 6.560298d-22
  td_E_eV(255) = 5.360000d+0 ; td_cs_cm2(255) = 6.505346d-22
  td_E_eV(256) = 5.380000d+0 ; td_cs_cm2(256) = 6.451056d-22
  td_E_eV(257) = 5.400000d+0 ; td_cs_cm2(257) = 6.397418d-22
  td_E_eV(258) = 5.420000d+0 ; td_cs_cm2(258) = 6.344421d-22
  td_E_eV(259) = 5.440000d+0 ; td_cs_cm2(259) = 6.292057d-22
  td_E_eV(260) = 5.460000d+0 ; td_cs_cm2(260) = 6.240314d-22
  td_E_eV(261) = 5.480000d+0 ; td_cs_cm2(261) = 6.189184d-22
  td_E_eV(262) = 5.500000d+0 ; td_cs_cm2(262) = 6.138657d-22
  td_E_eV(263) = 5.520000d+0 ; td_cs_cm2(263) = 6.088723d-22
  td_E_eV(264) = 5.540000d+0 ; td_cs_cm2(264) = 6.039374d-22
  td_E_eV(265) = 5.560000d+0 ; td_cs_cm2(265) = 5.990601d-22
  td_E_eV(266) = 5.580000d+0 ; td_cs_cm2(266) = 5.942394d-22
  td_E_eV(267) = 5.600000d+0 ; td_cs_cm2(267) = 5.894746d-22
  td_E_eV(268) = 5.620000d+0 ; td_cs_cm2(268) = 5.847648d-22
  td_E_eV(269) = 5.640000d+0 ; td_cs_cm2(269) = 5.801091d-22
  td_E_eV(270) = 5.660000d+0 ; td_cs_cm2(270) = 5.755068d-22
  td_E_eV(271) = 5.680000d+0 ; td_cs_cm2(271) = 5.709570d-22
  td_E_eV(272) = 5.700000d+0 ; td_cs_cm2(272) = 5.664591d-22
  td_E_eV(273) = 5.720000d+0 ; td_cs_cm2(273) = 5.620121d-22
  td_E_eV(274) = 5.740000d+0 ; td_cs_cm2(274) = 5.576153d-22
  td_E_eV(275) = 5.760000d+0 ; td_cs_cm2(275) = 5.532681d-22
  td_E_eV(276) = 5.780000d+0 ; td_cs_cm2(276) = 5.489697d-22
  td_E_eV(277) = 5.800000d+0 ; td_cs_cm2(277) = 5.447193d-22
  td_E_eV(278) = 5.820000d+0 ; td_cs_cm2(278) = 5.405163d-22
  td_E_eV(279) = 5.840000d+0 ; td_cs_cm2(279) = 5.363600d-22
  td_E_eV(280) = 5.860000d+0 ; td_cs_cm2(280) = 5.322497d-22
  td_E_eV(281) = 5.880000d+0 ; td_cs_cm2(281) = 5.281847d-22
  td_E_eV(282) = 5.900000d+0 ; td_cs_cm2(282) = 5.241645d-22
  td_E_eV(283) = 5.920000d+0 ; td_cs_cm2(283) = 5.201882d-22
  td_E_eV(284) = 5.940000d+0 ; td_cs_cm2(284) = 5.162554d-22
  td_E_eV(285) = 5.960000d+0 ; td_cs_cm2(285) = 5.123654d-22
  td_E_eV(286) = 5.980000d+0 ; td_cs_cm2(286) = 5.085176d-22
  td_E_eV(287) = 6.000000d+0 ; td_cs_cm2(287) = 5.047114d-22
  td_E_eV(288) = 6.020000d+0 ; td_cs_cm2(288) = 5.009462d-22
  td_E_eV(289) = 6.040000d+0 ; td_cs_cm2(289) = 4.972215d-22
  td_E_eV(290) = 6.060000d+0 ; td_cs_cm2(290) = 4.935366d-22
  td_E_eV(291) = 6.080000d+0 ; td_cs_cm2(291) = 4.898910d-22
  td_E_eV(292) = 6.100000d+0 ; td_cs_cm2(292) = 4.862842d-22
  td_E_eV(293) = 6.120000d+0 ; td_cs_cm2(293) = 4.827157d-22
  td_E_eV(294) = 6.140000d+0 ; td_cs_cm2(294) = 4.791848d-22
  td_E_eV(295) = 6.160000d+0 ; td_cs_cm2(295) = 4.756911d-22
  td_E_eV(296) = 6.180000d+0 ; td_cs_cm2(296) = 4.722341d-22
  td_E_eV(297) = 6.200000d+0 ; td_cs_cm2(297) = 4.688133d-22
  td_E_eV(298) = 6.220000d+0 ; td_cs_cm2(298) = 4.654281d-22
  td_E_eV(299) = 6.240000d+0 ; td_cs_cm2(299) = 4.620781d-22
  td_E_eV(300) = 6.260000d+0 ; td_cs_cm2(300) = 4.587629d-22
  td_E_eV(301) = 6.280000d+0 ; td_cs_cm2(301) = 4.554819d-22
  td_E_eV(302) = 6.300000d+0 ; td_cs_cm2(302) = 4.522347d-22
  td_E_eV(303) = 6.320000d+0 ; td_cs_cm2(303) = 4.490208d-22
  td_E_eV(304) = 6.340000d+0 ; td_cs_cm2(304) = 4.458398d-22
  td_E_eV(305) = 6.360000d+0 ; td_cs_cm2(305) = 4.426913d-22
  td_E_eV(306) = 6.380000d+0 ; td_cs_cm2(306) = 4.395747d-22
  td_E_eV(307) = 6.400000d+0 ; td_cs_cm2(307) = 4.364898d-22
  td_E_eV(308) = 6.420000d+0 ; td_cs_cm2(308) = 4.334360d-22
  td_E_eV(309) = 6.440000d+0 ; td_cs_cm2(309) = 4.304130d-22
  td_E_eV(310) = 6.460000d+0 ; td_cs_cm2(310) = 4.274204d-22
  td_E_eV(311) = 6.480000d+0 ; td_cs_cm2(311) = 4.244577d-22
  td_E_eV(312) = 6.500000d+0 ; td_cs_cm2(312) = 4.215246d-22
  td_E_eV(313) = 6.520000d+0 ; td_cs_cm2(313) = 4.186207d-22
  td_E_eV(314) = 6.540000d+0 ; td_cs_cm2(314) = 4.157455d-22
  td_E_eV(315) = 6.560000d+0 ; td_cs_cm2(315) = 4.128989d-22
  td_E_eV(316) = 6.580000d+0 ; td_cs_cm2(316) = 4.100803d-22
  td_E_eV(317) = 6.600000d+0 ; td_cs_cm2(317) = 4.072894d-22
  td_E_eV(318) = 6.620000d+0 ; td_cs_cm2(318) = 4.045258d-22
  td_E_eV(319) = 6.640000d+0 ; td_cs_cm2(319) = 4.017893d-22
  td_E_eV(320) = 6.660000d+0 ; td_cs_cm2(320) = 3.990794d-22
  td_E_eV(321) = 6.680000d+0 ; td_cs_cm2(321) = 3.963958d-22
  td_E_eV(322) = 6.700000d+0 ; td_cs_cm2(322) = 3.937382d-22
  td_E_eV(323) = 6.720000d+0 ; td_cs_cm2(323) = 3.911063d-22
  td_E_eV(324) = 6.740000d+0 ; td_cs_cm2(324) = 3.884997d-22
  td_E_eV(325) = 6.760000d+0 ; td_cs_cm2(325) = 3.859182d-22
  td_E_eV(326) = 6.780000d+0 ; td_cs_cm2(326) = 3.833613d-22
  td_E_eV(327) = 6.800000d+0 ; td_cs_cm2(327) = 3.808288d-22
  td_E_eV(328) = 6.820000d+0 ; td_cs_cm2(328) = 3.783205d-22
  td_E_eV(329) = 6.840000d+0 ; td_cs_cm2(329) = 3.758359d-22
  td_E_eV(330) = 6.860000d+0 ; td_cs_cm2(330) = 3.733748d-22
  td_E_eV(331) = 6.880000d+0 ; td_cs_cm2(331) = 3.709369d-22
  td_E_eV(332) = 6.900000d+0 ; td_cs_cm2(332) = 3.685220d-22
  td_E_eV(333) = 6.920000d+0 ; td_cs_cm2(333) = 3.661297d-22
  td_E_eV(334) = 6.940000d+0 ; td_cs_cm2(334) = 3.637598d-22
  td_E_eV(335) = 6.960000d+0 ; td_cs_cm2(335) = 3.614119d-22
  td_E_eV(336) = 6.980000d+0 ; td_cs_cm2(336) = 3.590859d-22
  td_E_eV(337) = 7.000000d+0 ; td_cs_cm2(337) = 3.567815d-22
  td_E_eV(338) = 7.020000d+0 ; td_cs_cm2(338) = 3.544983d-22
  td_E_eV(339) = 7.040000d+0 ; td_cs_cm2(339) = 3.522362d-22
  td_E_eV(340) = 7.060000d+0 ; td_cs_cm2(340) = 3.499949d-22
  td_E_eV(341) = 7.080000d+0 ; td_cs_cm2(341) = 3.477741d-22
  td_E_eV(342) = 7.100000d+0 ; td_cs_cm2(342) = 3.455737d-22
  td_E_eV(343) = 7.120000d+0 ; td_cs_cm2(343) = 3.433932d-22
  td_E_eV(344) = 7.140000d+0 ; td_cs_cm2(344) = 3.412326d-22
  td_E_eV(345) = 7.160000d+0 ; td_cs_cm2(345) = 3.390916d-22
  td_E_eV(346) = 7.180000d+0 ; td_cs_cm2(346) = 3.369699d-22
  td_E_eV(347) = 7.200000d+0 ; td_cs_cm2(347) = 3.348674d-22
  td_E_eV(348) = 7.220000d+0 ; td_cs_cm2(348) = 3.327837d-22
  td_E_eV(349) = 7.240000d+0 ; td_cs_cm2(349) = 3.307187d-22
  td_E_eV(350) = 7.260000d+0 ; td_cs_cm2(350) = 3.286722d-22
  td_E_eV(351) = 7.280000d+0 ; td_cs_cm2(351) = 3.266439d-22
  td_E_eV(352) = 7.300000d+0 ; td_cs_cm2(352) = 3.246337d-22
  td_E_eV(353) = 7.320000d+0 ; td_cs_cm2(353) = 3.226413d-22
  td_E_eV(354) = 7.340000d+0 ; td_cs_cm2(354) = 3.206664d-22
  td_E_eV(355) = 7.360000d+0 ; td_cs_cm2(355) = 3.187090d-22
  td_E_eV(356) = 7.380000d+0 ; td_cs_cm2(356) = 3.167688d-22
  td_E_eV(357) = 7.400000d+0 ; td_cs_cm2(357) = 3.148457d-22
  td_E_eV(358) = 7.420000d+0 ; td_cs_cm2(358) = 3.129393d-22
  td_E_eV(359) = 7.440000d+0 ; td_cs_cm2(359) = 3.110496d-22
  td_E_eV(360) = 7.460000d+0 ; td_cs_cm2(360) = 3.091763d-22
  td_E_eV(361) = 7.480000d+0 ; td_cs_cm2(361) = 3.073192d-22
  td_E_eV(362) = 7.500000d+0 ; td_cs_cm2(362) = 3.054783d-22
  td_E_eV(363) = 7.520000d+0 ; td_cs_cm2(363) = 3.036532d-22
  td_E_eV(364) = 7.540000d+0 ; td_cs_cm2(364) = 3.018438d-22
  td_E_eV(365) = 7.560000d+0 ; td_cs_cm2(365) = 3.000500d-22
  td_E_eV(366) = 7.580000d+0 ; td_cs_cm2(366) = 2.982715d-22
  td_E_eV(367) = 7.600000d+0 ; td_cs_cm2(367) = 2.965082d-22
  td_E_eV(368) = 7.620000d+0 ; td_cs_cm2(368) = 2.947599d-22
  td_E_eV(369) = 7.640000d+0 ; td_cs_cm2(369) = 2.930265d-22
  td_E_eV(370) = 7.660000d+0 ; td_cs_cm2(370) = 2.913077d-22
  td_E_eV(371) = 7.680000d+0 ; td_cs_cm2(371) = 2.896035d-22
  td_E_eV(372) = 7.700000d+0 ; td_cs_cm2(372) = 2.879136d-22
  td_E_eV(373) = 7.720000d+0 ; td_cs_cm2(373) = 2.862380d-22
  td_E_eV(374) = 7.740000d+0 ; td_cs_cm2(374) = 2.845764d-22
  td_E_eV(375) = 7.760000d+0 ; td_cs_cm2(375) = 2.829287d-22
  td_E_eV(376) = 7.780000d+0 ; td_cs_cm2(376) = 2.812947d-22
  td_E_eV(377) = 7.800000d+0 ; td_cs_cm2(377) = 2.796743d-22
  td_E_eV(378) = 7.820000d+0 ; td_cs_cm2(378) = 2.780674d-22
  td_E_eV(379) = 7.840000d+0 ; td_cs_cm2(379) = 2.764738d-22
  td_E_eV(380) = 7.860000d+0 ; td_cs_cm2(380) = 2.748933d-22
  td_E_eV(381) = 7.880000d+0 ; td_cs_cm2(381) = 2.733259d-22
  td_E_eV(382) = 7.900000d+0 ; td_cs_cm2(382) = 2.717713d-22
  td_E_eV(383) = 7.920000d+0 ; td_cs_cm2(383) = 2.702295d-22
  td_E_eV(384) = 7.940000d+0 ; td_cs_cm2(384) = 2.687003d-22
  td_E_eV(385) = 7.960000d+0 ; td_cs_cm2(385) = 2.671835d-22
  td_E_eV(386) = 7.980000d+0 ; td_cs_cm2(386) = 2.656791d-22
  td_E_eV(387) = 8.000000d+0 ; td_cs_cm2(387) = 2.641869d-22
  td_E_eV(388) = 8.020000d+0 ; td_cs_cm2(388) = 2.627068d-22
  td_E_eV(389) = 8.040000d+0 ; td_cs_cm2(389) = 2.612386d-22
  td_E_eV(390) = 8.060000d+0 ; td_cs_cm2(390) = 2.597822d-22
  td_E_eV(391) = 8.080000d+0 ; td_cs_cm2(391) = 2.583375d-22
  td_E_eV(392) = 8.100000d+0 ; td_cs_cm2(392) = 2.569044d-22
  td_E_eV(393) = 8.120000d+0 ; td_cs_cm2(393) = 2.554828d-22
  td_E_eV(394) = 8.140000d+0 ; td_cs_cm2(394) = 2.540725d-22
  td_E_eV(395) = 8.160000d+0 ; td_cs_cm2(395) = 2.526734d-22
  td_E_eV(396) = 8.180000d+0 ; td_cs_cm2(396) = 2.512854d-22
  td_E_eV(397) = 8.200000d+0 ; td_cs_cm2(397) = 2.499084d-22
  td_E_eV(398) = 8.220000d+0 ; td_cs_cm2(398) = 2.485423d-22
  td_E_eV(399) = 8.240000d+0 ; td_cs_cm2(399) = 2.471869d-22
  td_E_eV(400) = 8.260000d+0 ; td_cs_cm2(400) = 2.458422d-22
  td_E_eV(401) = 8.280000d+0 ; td_cs_cm2(401) = 2.445080d-22
  td_E_eV(402) = 8.300000d+0 ; td_cs_cm2(402) = 2.431843d-22
  td_E_eV(403) = 8.320000d+0 ; td_cs_cm2(403) = 2.418709d-22
  td_E_eV(404) = 8.340000d+0 ; td_cs_cm2(404) = 2.405677d-22
  td_E_eV(405) = 8.360000d+0 ; td_cs_cm2(405) = 2.392746d-22
  td_E_eV(406) = 8.380000d+0 ; td_cs_cm2(406) = 2.379915d-22
  td_E_eV(407) = 8.400000d+0 ; td_cs_cm2(407) = 2.367184d-22
  td_E_eV(408) = 8.420000d+0 ; td_cs_cm2(408) = 2.354551d-22
  td_E_eV(409) = 8.440000d+0 ; td_cs_cm2(409) = 2.342014d-22
  td_E_eV(410) = 8.460000d+0 ; td_cs_cm2(410) = 2.329574d-22
  td_E_eV(411) = 8.480000d+0 ; td_cs_cm2(411) = 2.317230d-22
  td_E_eV(412) = 8.500000d+0 ; td_cs_cm2(412) = 2.304979d-22
  td_E_eV(413) = 8.520000d+0 ; td_cs_cm2(413) = 2.292822d-22
  td_E_eV(414) = 8.540000d+0 ; td_cs_cm2(414) = 2.280757d-22
  td_E_eV(415) = 8.560000d+0 ; td_cs_cm2(415) = 2.268784d-22
  td_E_eV(416) = 8.580000d+0 ; td_cs_cm2(416) = 2.256901d-22
  td_E_eV(417) = 8.600000d+0 ; td_cs_cm2(417) = 2.245108d-22
  td_E_eV(418) = 8.620000d+0 ; td_cs_cm2(418) = 2.233404d-22
  td_E_eV(419) = 8.640000d+0 ; td_cs_cm2(419) = 2.221788d-22
  td_E_eV(420) = 8.660000d+0 ; td_cs_cm2(420) = 2.210258d-22
  td_E_eV(421) = 8.680000d+0 ; td_cs_cm2(421) = 2.198815d-22
  td_E_eV(422) = 8.700000d+0 ; td_cs_cm2(422) = 2.187458d-22
  td_E_eV(423) = 8.720000d+0 ; td_cs_cm2(423) = 2.176185d-22
  td_E_eV(424) = 8.740000d+0 ; td_cs_cm2(424) = 2.164995d-22
  td_E_eV(425) = 8.760000d+0 ; td_cs_cm2(425) = 2.153889d-22
  td_E_eV(426) = 8.780000d+0 ; td_cs_cm2(426) = 2.142865d-22
  td_E_eV(427) = 8.800000d+0 ; td_cs_cm2(427) = 2.131921d-22
  td_E_eV(428) = 8.820000d+0 ; td_cs_cm2(428) = 2.121059d-22
  td_E_eV(429) = 8.840000d+0 ; td_cs_cm2(429) = 2.110276d-22
  td_E_eV(430) = 8.860000d+0 ; td_cs_cm2(430) = 2.099573d-22
  td_E_eV(431) = 8.880000d+0 ; td_cs_cm2(431) = 2.088947d-22
  td_E_eV(432) = 8.900000d+0 ; td_cs_cm2(432) = 2.078399d-22
  td_E_eV(433) = 8.920000d+0 ; td_cs_cm2(433) = 2.067928d-22
  td_E_eV(434) = 8.940000d+0 ; td_cs_cm2(434) = 2.057533d-22
  td_E_eV(435) = 8.960000d+0 ; td_cs_cm2(435) = 2.047213d-22
  td_E_eV(436) = 8.980000d+0 ; td_cs_cm2(436) = 2.036967d-22
  td_E_eV(437) = 9.000000d+0 ; td_cs_cm2(437) = 2.026796d-22
  td_E_eV(438) = 9.020000d+0 ; td_cs_cm2(438) = 2.016698d-22
  td_E_eV(439) = 9.040000d+0 ; td_cs_cm2(439) = 2.006672d-22
  td_E_eV(440) = 9.060000d+0 ; td_cs_cm2(440) = 1.996718d-22
  td_E_eV(441) = 9.080000d+0 ; td_cs_cm2(441) = 1.986836d-22
  td_E_eV(442) = 9.100000d+0 ; td_cs_cm2(442) = 1.977023d-22
  td_E_eV(443) = 9.120000d+0 ; td_cs_cm2(443) = 1.967281d-22
  td_E_eV(444) = 9.140000d+0 ; td_cs_cm2(444) = 1.957608d-22
  td_E_eV(445) = 9.160000d+0 ; td_cs_cm2(445) = 1.948003d-22
  td_E_eV(446) = 9.180000d+0 ; td_cs_cm2(446) = 1.938466d-22
  td_E_eV(447) = 9.200000d+0 ; td_cs_cm2(447) = 1.928997d-22
  td_E_eV(448) = 9.220000d+0 ; td_cs_cm2(448) = 1.919594d-22
  td_E_eV(449) = 9.240000d+0 ; td_cs_cm2(449) = 1.910258d-22
  td_E_eV(450) = 9.260000d+0 ; td_cs_cm2(450) = 1.900986d-22
  td_E_eV(451) = 9.280000d+0 ; td_cs_cm2(451) = 1.891780d-22
  td_E_eV(452) = 9.300000d+0 ; td_cs_cm2(452) = 1.882638d-22
  td_E_eV(453) = 9.320000d+0 ; td_cs_cm2(453) = 1.873559d-22
  td_E_eV(454) = 9.340000d+0 ; td_cs_cm2(454) = 1.864544d-22
  td_E_eV(455) = 9.360000d+0 ; td_cs_cm2(455) = 1.855591d-22
  td_E_eV(456) = 9.380000d+0 ; td_cs_cm2(456) = 1.846700d-22
  td_E_eV(457) = 9.400000d+0 ; td_cs_cm2(457) = 1.837871d-22
  td_E_eV(458) = 9.420000d+0 ; td_cs_cm2(458) = 1.829102d-22
  td_E_eV(459) = 9.440000d+0 ; td_cs_cm2(459) = 1.820394d-22
  td_E_eV(460) = 9.460000d+0 ; td_cs_cm2(460) = 1.811746d-22
  td_E_eV(461) = 9.480000d+0 ; td_cs_cm2(461) = 1.803156d-22
  td_E_eV(462) = 9.500000d+0 ; td_cs_cm2(462) = 1.794626d-22
  td_E_eV(463) = 9.520000d+0 ; td_cs_cm2(463) = 1.786153d-22
  td_E_eV(464) = 9.540000d+0 ; td_cs_cm2(464) = 1.777738d-22
  td_E_eV(465) = 9.560000d+0 ; td_cs_cm2(465) = 1.769381d-22
  td_E_eV(466) = 9.580000d+0 ; td_cs_cm2(466) = 1.761080d-22
  td_E_eV(467) = 9.600000d+0 ; td_cs_cm2(467) = 1.752835d-22
  td_E_eV(468) = 9.620000d+0 ; td_cs_cm2(468) = 1.744646d-22
  td_E_eV(469) = 9.640000d+0 ; td_cs_cm2(469) = 1.736511d-22
  td_E_eV(470) = 9.660000d+0 ; td_cs_cm2(470) = 1.728432d-22
  td_E_eV(471) = 9.680000d+0 ; td_cs_cm2(471) = 1.720407d-22
  td_E_eV(472) = 9.700000d+0 ; td_cs_cm2(472) = 1.712435d-22
  td_E_eV(473) = 9.720000d+0 ; td_cs_cm2(473) = 1.704517d-22
  td_E_eV(474) = 9.740000d+0 ; td_cs_cm2(474) = 1.696651d-22
  td_E_eV(475) = 9.760000d+0 ; td_cs_cm2(475) = 1.688838d-22
  td_E_eV(476) = 9.780000d+0 ; td_cs_cm2(476) = 1.681077d-22
  td_E_eV(477) = 9.800000d+0 ; td_cs_cm2(477) = 1.673367d-22
  td_E_eV(478) = 9.820000d+0 ; td_cs_cm2(478) = 1.665708d-22
  td_E_eV(479) = 9.840000d+0 ; td_cs_cm2(479) = 1.658100d-22
  td_E_eV(480) = 9.860000d+0 ; td_cs_cm2(480) = 1.650541d-22
  td_E_eV(481) = 9.880000d+0 ; td_cs_cm2(481) = 1.643033d-22
  td_E_eV(482) = 9.900000d+0 ; td_cs_cm2(482) = 1.635573d-22
  td_E_eV(483) = 9.920000d+0 ; td_cs_cm2(483) = 1.628163d-22
  td_E_eV(484) = 9.940000d+0 ; td_cs_cm2(484) = 1.620800d-22
  td_E_eV(485) = 9.960000d+0 ; td_cs_cm2(485) = 1.613486d-22
  td_E_eV(486) = 9.980000d+0 ; td_cs_cm2(486) = 1.606219d-22
  td_E_eV(487) = 1.000000d+1 ; td_cs_cm2(487) = 1.599000d-22
  td_E_eV(488) = 1.002000d+1 ; td_cs_cm2(488) = 1.600850d-22
  td_E_eV(489) = 1.004000d+1 ; td_cs_cm2(489) = 1.602699d-22
  td_E_eV(490) = 1.006000d+1 ; td_cs_cm2(490) = 1.604547d-22
  td_E_eV(491) = 1.008000d+1 ; td_cs_cm2(491) = 1.606392d-22
  td_E_eV(492) = 1.010000d+1 ; td_cs_cm2(492) = 1.608237d-22
  td_E_eV(493) = 1.012000d+1 ; td_cs_cm2(493) = 1.610079d-22
  td_E_eV(494) = 1.014000d+1 ; td_cs_cm2(494) = 1.611920d-22
  td_E_eV(495) = 1.016000d+1 ; td_cs_cm2(495) = 1.613760d-22
  td_E_eV(496) = 1.018000d+1 ; td_cs_cm2(496) = 1.615598d-22
  td_E_eV(497) = 1.020000d+1 ; td_cs_cm2(497) = 1.617435d-22
  td_E_eV(498) = 1.022000d+1 ; td_cs_cm2(498) = 1.619270d-22
  td_E_eV(499) = 1.024000d+1 ; td_cs_cm2(499) = 1.621103d-22
  td_E_eV(500) = 1.026000d+1 ; td_cs_cm2(500) = 1.622935d-22
  td_E_eV(501) = 1.028000d+1 ; td_cs_cm2(501) = 1.624766d-22
  td_E_eV(502) = 1.030000d+1 ; td_cs_cm2(502) = 1.626595d-22
  td_E_eV(503) = 1.032000d+1 ; td_cs_cm2(503) = 1.628423d-22
  td_E_eV(504) = 1.034000d+1 ; td_cs_cm2(504) = 1.630249d-22
  td_E_eV(505) = 1.036000d+1 ; td_cs_cm2(505) = 1.632073d-22
  td_E_eV(506) = 1.038000d+1 ; td_cs_cm2(506) = 1.633896d-22
  td_E_eV(507) = 1.040000d+1 ; td_cs_cm2(507) = 1.635718d-22
  td_E_eV(508) = 1.042000d+1 ; td_cs_cm2(508) = 1.637538d-22
  td_E_eV(509) = 1.044000d+1 ; td_cs_cm2(509) = 1.639357d-22
  td_E_eV(510) = 1.046000d+1 ; td_cs_cm2(510) = 1.641174d-22
  td_E_eV(511) = 1.048000d+1 ; td_cs_cm2(511) = 1.642990d-22
  td_E_eV(512) = 1.050000d+1 ; td_cs_cm2(512) = 1.644804d-22
  td_E_eV(513) = 1.052000d+1 ; td_cs_cm2(513) = 1.646617d-22
  td_E_eV(514) = 1.054000d+1 ; td_cs_cm2(514) = 1.648428d-22
  td_E_eV(515) = 1.056000d+1 ; td_cs_cm2(515) = 1.650238d-22
  td_E_eV(516) = 1.058000d+1 ; td_cs_cm2(516) = 1.652047d-22
  td_E_eV(517) = 1.060000d+1 ; td_cs_cm2(517) = 1.653854d-22
  td_E_eV(518) = 1.062000d+1 ; td_cs_cm2(518) = 1.655659d-22
  td_E_eV(519) = 1.064000d+1 ; td_cs_cm2(519) = 1.657463d-22
  td_E_eV(520) = 1.066000d+1 ; td_cs_cm2(520) = 1.659266d-22
  td_E_eV(521) = 1.068000d+1 ; td_cs_cm2(521) = 1.661068d-22
  td_E_eV(522) = 1.070000d+1 ; td_cs_cm2(522) = 1.662867d-22
  td_E_eV(523) = 1.072000d+1 ; td_cs_cm2(523) = 1.664666d-22
  td_E_eV(524) = 1.074000d+1 ; td_cs_cm2(524) = 1.666463d-22
  td_E_eV(525) = 1.076000d+1 ; td_cs_cm2(525) = 1.668259d-22
  td_E_eV(526) = 1.078000d+1 ; td_cs_cm2(526) = 1.670053d-22
  td_E_eV(527) = 1.080000d+1 ; td_cs_cm2(527) = 1.671846d-22
  td_E_eV(528) = 1.082000d+1 ; td_cs_cm2(528) = 1.673637d-22
  td_E_eV(529) = 1.084000d+1 ; td_cs_cm2(529) = 1.675427d-22
  td_E_eV(530) = 1.086000d+1 ; td_cs_cm2(530) = 1.677216d-22
  td_E_eV(531) = 1.088000d+1 ; td_cs_cm2(531) = 1.679003d-22
  td_E_eV(532) = 1.090000d+1 ; td_cs_cm2(532) = 1.680789d-22
  td_E_eV(533) = 1.092000d+1 ; td_cs_cm2(533) = 1.682574d-22
  td_E_eV(534) = 1.094000d+1 ; td_cs_cm2(534) = 1.684357d-22
  td_E_eV(535) = 1.096000d+1 ; td_cs_cm2(535) = 1.686139d-22
  td_E_eV(536) = 1.098000d+1 ; td_cs_cm2(536) = 1.687919d-22
  td_E_eV(537) = 1.100000d+1 ; td_cs_cm2(537) = 1.689698d-22
  td_E_eV(538) = 1.102000d+1 ; td_cs_cm2(538) = 1.691476d-22
  td_E_eV(539) = 1.104000d+1 ; td_cs_cm2(539) = 1.693252d-22
  td_E_eV(540) = 1.106000d+1 ; td_cs_cm2(540) = 1.695027d-22
  td_E_eV(541) = 1.108000d+1 ; td_cs_cm2(541) = 1.696801d-22
  td_E_eV(542) = 1.110000d+1 ; td_cs_cm2(542) = 1.698573d-22
  td_E_eV(543) = 1.112000d+1 ; td_cs_cm2(543) = 1.700344d-22
  td_E_eV(544) = 1.114000d+1 ; td_cs_cm2(544) = 1.702114d-22
  td_E_eV(545) = 1.116000d+1 ; td_cs_cm2(545) = 1.703882d-22
  td_E_eV(546) = 1.118000d+1 ; td_cs_cm2(546) = 1.705649d-22
  td_E_eV(547) = 1.120000d+1 ; td_cs_cm2(547) = 1.707414d-22
  td_E_eV(548) = 1.122000d+1 ; td_cs_cm2(548) = 1.709179d-22
  td_E_eV(549) = 1.124000d+1 ; td_cs_cm2(549) = 1.710942d-22
  td_E_eV(550) = 1.126000d+1 ; td_cs_cm2(550) = 1.712703d-22
  td_E_eV(551) = 1.128000d+1 ; td_cs_cm2(551) = 1.714464d-22
  td_E_eV(552) = 1.130000d+1 ; td_cs_cm2(552) = 1.716223d-22
  td_E_eV(553) = 1.132000d+1 ; td_cs_cm2(553) = 1.717980d-22
  td_E_eV(554) = 1.134000d+1 ; td_cs_cm2(554) = 1.719737d-22
  td_E_eV(555) = 1.136000d+1 ; td_cs_cm2(555) = 1.721492d-22
  td_E_eV(556) = 1.138000d+1 ; td_cs_cm2(556) = 1.723245d-22
  td_E_eV(557) = 1.140000d+1 ; td_cs_cm2(557) = 1.724998d-22
  td_E_eV(558) = 1.142000d+1 ; td_cs_cm2(558) = 1.726749d-22
  td_E_eV(559) = 1.144000d+1 ; td_cs_cm2(559) = 1.728499d-22
  td_E_eV(560) = 1.146000d+1 ; td_cs_cm2(560) = 1.730248d-22
  td_E_eV(561) = 1.148000d+1 ; td_cs_cm2(561) = 1.731995d-22
  td_E_eV(562) = 1.150000d+1 ; td_cs_cm2(562) = 1.733741d-22
  td_E_eV(563) = 1.152000d+1 ; td_cs_cm2(563) = 1.735486d-22
  td_E_eV(564) = 1.154000d+1 ; td_cs_cm2(564) = 1.737229d-22
  td_E_eV(565) = 1.156000d+1 ; td_cs_cm2(565) = 1.738971d-22
  td_E_eV(566) = 1.158000d+1 ; td_cs_cm2(566) = 1.740712d-22
  td_E_eV(567) = 1.160000d+1 ; td_cs_cm2(567) = 1.742452d-22
  td_E_eV(568) = 1.162000d+1 ; td_cs_cm2(568) = 1.744190d-22
  td_E_eV(569) = 1.164000d+1 ; td_cs_cm2(569) = 1.745927d-22
  td_E_eV(570) = 1.166000d+1 ; td_cs_cm2(570) = 1.747663d-22
  td_E_eV(571) = 1.168000d+1 ; td_cs_cm2(571) = 1.749398d-22
  td_E_eV(572) = 1.170000d+1 ; td_cs_cm2(572) = 1.751131d-22
  td_E_eV(573) = 1.172000d+1 ; td_cs_cm2(573) = 1.752863d-22
  td_E_eV(574) = 1.174000d+1 ; td_cs_cm2(574) = 1.754594d-22
  td_E_eV(575) = 1.176000d+1 ; td_cs_cm2(575) = 1.756324d-22
  td_E_eV(576) = 1.178000d+1 ; td_cs_cm2(576) = 1.758052d-22
  td_E_eV(577) = 1.180000d+1 ; td_cs_cm2(577) = 1.759780d-22
  td_E_eV(578) = 1.182000d+1 ; td_cs_cm2(578) = 1.761506d-22
  td_E_eV(579) = 1.184000d+1 ; td_cs_cm2(579) = 1.763230d-22
  td_E_eV(580) = 1.186000d+1 ; td_cs_cm2(580) = 1.764954d-22
  td_E_eV(581) = 1.188000d+1 ; td_cs_cm2(581) = 1.766676d-22
  td_E_eV(582) = 1.190000d+1 ; td_cs_cm2(582) = 1.768397d-22
  td_E_eV(583) = 1.192000d+1 ; td_cs_cm2(583) = 1.770117d-22
  td_E_eV(584) = 1.194000d+1 ; td_cs_cm2(584) = 1.771836d-22
  td_E_eV(585) = 1.196000d+1 ; td_cs_cm2(585) = 1.773553d-22
  td_E_eV(586) = 1.198000d+1 ; td_cs_cm2(586) = 1.775269d-22
  td_E_eV(587) = 1.200000d+1 ; td_cs_cm2(587) = 1.776984d-22
  td_E_eV(588) = 1.202000d+1 ; td_cs_cm2(588) = 1.778698d-22
  td_E_eV(589) = 1.204000d+1 ; td_cs_cm2(589) = 1.780410d-22
  td_E_eV(590) = 1.206000d+1 ; td_cs_cm2(590) = 1.782122d-22
  td_E_eV(591) = 1.208000d+1 ; td_cs_cm2(591) = 1.783832d-22
  td_E_eV(592) = 1.210000d+1 ; td_cs_cm2(592) = 1.785541d-22
  td_E_eV(593) = 1.212000d+1 ; td_cs_cm2(593) = 1.787249d-22
  td_E_eV(594) = 1.214000d+1 ; td_cs_cm2(594) = 1.788955d-22
  td_E_eV(595) = 1.216000d+1 ; td_cs_cm2(595) = 1.790661d-22
  td_E_eV(596) = 1.218000d+1 ; td_cs_cm2(596) = 1.792365d-22
  td_E_eV(597) = 1.220000d+1 ; td_cs_cm2(597) = 1.794068d-22
  td_E_eV(598) = 1.222000d+1 ; td_cs_cm2(598) = 1.795770d-22
  td_E_eV(599) = 1.224000d+1 ; td_cs_cm2(599) = 1.797471d-22
  td_E_eV(600) = 1.226000d+1 ; td_cs_cm2(600) = 1.799170d-22
  td_E_eV(601) = 1.228000d+1 ; td_cs_cm2(601) = 1.800869d-22
  td_E_eV(602) = 1.230000d+1 ; td_cs_cm2(602) = 1.802566d-22
  td_E_eV(603) = 1.232000d+1 ; td_cs_cm2(603) = 1.804262d-22
  td_E_eV(604) = 1.234000d+1 ; td_cs_cm2(604) = 1.805957d-22
  td_E_eV(605) = 1.236000d+1 ; td_cs_cm2(605) = 1.807651d-22
  td_E_eV(606) = 1.238000d+1 ; td_cs_cm2(606) = 1.809343d-22
  td_E_eV(607) = 1.240000d+1 ; td_cs_cm2(607) = 1.811035d-22
  td_E_eV(608) = 1.242000d+1 ; td_cs_cm2(608) = 1.812725d-22
  td_E_eV(609) = 1.244000d+1 ; td_cs_cm2(609) = 1.814414d-22
  td_E_eV(610) = 1.246000d+1 ; td_cs_cm2(610) = 1.816102d-22
  td_E_eV(611) = 1.248000d+1 ; td_cs_cm2(611) = 1.817789d-22
  td_E_eV(612) = 1.250000d+1 ; td_cs_cm2(612) = 1.819475d-22
  td_E_eV(613) = 1.252000d+1 ; td_cs_cm2(613) = 1.821159d-22
  td_E_eV(614) = 1.254000d+1 ; td_cs_cm2(614) = 1.822843d-22
  td_E_eV(615) = 1.256000d+1 ; td_cs_cm2(615) = 1.824525d-22
  td_E_eV(616) = 1.258000d+1 ; td_cs_cm2(616) = 1.826206d-22
  td_E_eV(617) = 1.260000d+1 ; td_cs_cm2(617) = 1.827887d-22
  td_E_eV(618) = 1.262000d+1 ; td_cs_cm2(618) = 1.829566d-22
  td_E_eV(619) = 1.264000d+1 ; td_cs_cm2(619) = 1.831243d-22
  td_E_eV(620) = 1.266000d+1 ; td_cs_cm2(620) = 1.832920d-22
  td_E_eV(621) = 1.268000d+1 ; td_cs_cm2(621) = 1.834596d-22
  td_E_eV(622) = 1.270000d+1 ; td_cs_cm2(622) = 1.836270d-22
  td_E_eV(623) = 1.272000d+1 ; td_cs_cm2(623) = 1.837944d-22
  td_E_eV(624) = 1.274000d+1 ; td_cs_cm2(624) = 1.839616d-22
  td_E_eV(625) = 1.276000d+1 ; td_cs_cm2(625) = 1.841287d-22
  td_E_eV(626) = 1.278000d+1 ; td_cs_cm2(626) = 1.842957d-22
  td_E_eV(627) = 1.280000d+1 ; td_cs_cm2(627) = 1.844626d-22
  td_E_eV(628) = 1.282000d+1 ; td_cs_cm2(628) = 1.846294d-22
  td_E_eV(629) = 1.284000d+1 ; td_cs_cm2(629) = 1.847961d-22
  td_E_eV(630) = 1.286000d+1 ; td_cs_cm2(630) = 1.849626d-22
  td_E_eV(631) = 1.288000d+1 ; td_cs_cm2(631) = 1.851291d-22
  td_E_eV(632) = 1.290000d+1 ; td_cs_cm2(632) = 1.852954d-22
  td_E_eV(633) = 1.292000d+1 ; td_cs_cm2(633) = 1.854617d-22
  td_E_eV(634) = 1.294000d+1 ; td_cs_cm2(634) = 1.856278d-22
  td_E_eV(635) = 1.296000d+1 ; td_cs_cm2(635) = 1.857938d-22
  td_E_eV(636) = 1.298000d+1 ; td_cs_cm2(636) = 1.859598d-22
  td_E_eV(637) = 1.300000d+1 ; td_cs_cm2(637) = 1.861256d-22
  td_E_eV(638) = 1.302000d+1 ; td_cs_cm2(638) = 1.862913d-22
  td_E_eV(639) = 1.304000d+1 ; td_cs_cm2(639) = 1.864569d-22
  td_E_eV(640) = 1.306000d+1 ; td_cs_cm2(640) = 1.866224d-22
  td_E_eV(641) = 1.308000d+1 ; td_cs_cm2(641) = 1.867877d-22
  td_E_eV(642) = 1.310000d+1 ; td_cs_cm2(642) = 1.869530d-22
  td_E_eV(643) = 1.312000d+1 ; td_cs_cm2(643) = 1.871182d-22
  td_E_eV(644) = 1.314000d+1 ; td_cs_cm2(644) = 1.872832d-22
  td_E_eV(645) = 1.316000d+1 ; td_cs_cm2(645) = 1.874482d-22
  td_E_eV(646) = 1.318000d+1 ; td_cs_cm2(646) = 1.876130d-22
  td_E_eV(647) = 1.320000d+1 ; td_cs_cm2(647) = 1.877778d-22
  td_E_eV(648) = 1.322000d+1 ; td_cs_cm2(648) = 1.879424d-22
  td_E_eV(649) = 1.324000d+1 ; td_cs_cm2(649) = 1.881070d-22
  td_E_eV(650) = 1.326000d+1 ; td_cs_cm2(650) = 1.882714d-22
  td_E_eV(651) = 1.328000d+1 ; td_cs_cm2(651) = 1.884357d-22
  td_E_eV(652) = 1.330000d+1 ; td_cs_cm2(652) = 1.886000d-22
  td_E_eV(653) = 1.332000d+1 ; td_cs_cm2(653) = 1.887641d-22
  td_E_eV(654) = 1.334000d+1 ; td_cs_cm2(654) = 1.889281d-22
  td_E_eV(655) = 1.336000d+1 ; td_cs_cm2(655) = 1.890920d-22
  td_E_eV(656) = 1.338000d+1 ; td_cs_cm2(656) = 1.892558d-22
  td_E_eV(657) = 1.340000d+1 ; td_cs_cm2(657) = 1.894195d-22
  td_E_eV(658) = 1.342000d+1 ; td_cs_cm2(658) = 1.895831d-22
  td_E_eV(659) = 1.344000d+1 ; td_cs_cm2(659) = 1.897466d-22
  td_E_eV(660) = 1.346000d+1 ; td_cs_cm2(660) = 1.899100d-22
  td_E_eV(661) = 1.348000d+1 ; td_cs_cm2(661) = 1.900733d-22
  td_E_eV(662) = 1.350000d+1 ; td_cs_cm2(662) = 1.902365d-22
  td_E_eV(663) = 1.352000d+1 ; td_cs_cm2(663) = 1.903996d-22
  td_E_eV(664) = 1.354000d+1 ; td_cs_cm2(664) = 1.905626d-22
  td_E_eV(665) = 1.356000d+1 ; td_cs_cm2(665) = 1.907255d-22
  td_E_eV(666) = 1.358000d+1 ; td_cs_cm2(666) = 1.908883d-22
  td_E_eV(667) = 1.360000d+1 ; td_cs_cm2(667) = 1.910509d-22
  td_E_eV(668) = 1.362000d+1 ; td_cs_cm2(668) = 1.912135d-22
  td_E_eV(669) = 1.364000d+1 ; td_cs_cm2(669) = 1.913760d-22
  td_E_eV(670) = 1.366000d+1 ; td_cs_cm2(670) = 1.915384d-22
  td_E_eV(671) = 1.368000d+1 ; td_cs_cm2(671) = 1.917007d-22
  td_E_eV(672) = 1.370000d+1 ; td_cs_cm2(672) = 1.918629d-22
  td_E_eV(673) = 1.372000d+1 ; td_cs_cm2(673) = 1.920249d-22
  td_E_eV(674) = 1.374000d+1 ; td_cs_cm2(674) = 1.921869d-22
  td_E_eV(675) = 1.376000d+1 ; td_cs_cm2(675) = 1.923488d-22
  td_E_eV(676) = 1.378000d+1 ; td_cs_cm2(676) = 1.925106d-22
  td_E_eV(677) = 1.380000d+1 ; td_cs_cm2(677) = 1.926723d-22
  td_E_eV(678) = 1.382000d+1 ; td_cs_cm2(678) = 1.928339d-22
  td_E_eV(679) = 1.384000d+1 ; td_cs_cm2(679) = 1.929954d-22
  td_E_eV(680) = 1.386000d+1 ; td_cs_cm2(680) = 1.931568d-22
  td_E_eV(681) = 1.388000d+1 ; td_cs_cm2(681) = 1.933181d-22
  td_E_eV(682) = 1.390000d+1 ; td_cs_cm2(682) = 1.934793d-22
  td_E_eV(683) = 1.392000d+1 ; td_cs_cm2(683) = 1.936404d-22
  td_E_eV(684) = 1.394000d+1 ; td_cs_cm2(684) = 1.938014d-22
  td_E_eV(685) = 1.396000d+1 ; td_cs_cm2(685) = 1.939623d-22
  td_E_eV(686) = 1.398000d+1 ; td_cs_cm2(686) = 1.941231d-22
  td_E_eV(687) = 1.400000d+1 ; td_cs_cm2(687) = 1.942838d-22
  td_E_eV(688) = 1.402000d+1 ; td_cs_cm2(688) = 1.944444d-22
  td_E_eV(689) = 1.404000d+1 ; td_cs_cm2(689) = 1.946049d-22
  td_E_eV(690) = 1.406000d+1 ; td_cs_cm2(690) = 1.947653d-22
  td_E_eV(691) = 1.408000d+1 ; td_cs_cm2(691) = 1.949257d-22
  td_E_eV(692) = 1.410000d+1 ; td_cs_cm2(692) = 1.950859d-22
  td_E_eV(693) = 1.412000d+1 ; td_cs_cm2(693) = 1.952460d-22
  td_E_eV(694) = 1.414000d+1 ; td_cs_cm2(694) = 1.954061d-22
  td_E_eV(695) = 1.416000d+1 ; td_cs_cm2(695) = 1.955660d-22
  td_E_eV(696) = 1.418000d+1 ; td_cs_cm2(696) = 1.957259d-22
  td_E_eV(697) = 1.420000d+1 ; td_cs_cm2(697) = 1.958856d-22
  td_E_eV(698) = 1.422000d+1 ; td_cs_cm2(698) = 1.960453d-22
  td_E_eV(699) = 1.424000d+1 ; td_cs_cm2(699) = 1.962048d-22
  td_E_eV(700) = 1.426000d+1 ; td_cs_cm2(700) = 1.963643d-22
  td_E_eV(701) = 1.428000d+1 ; td_cs_cm2(701) = 1.965237d-22
  td_E_eV(702) = 1.430000d+1 ; td_cs_cm2(702) = 1.966830d-22
  td_E_eV(703) = 1.432000d+1 ; td_cs_cm2(703) = 1.968421d-22
  td_E_eV(704) = 1.434000d+1 ; td_cs_cm2(704) = 1.970012d-22
  td_E_eV(705) = 1.436000d+1 ; td_cs_cm2(705) = 1.971602d-22
  td_E_eV(706) = 1.438000d+1 ; td_cs_cm2(706) = 1.973191d-22
  td_E_eV(707) = 1.440000d+1 ; td_cs_cm2(707) = 1.974780d-22
  td_E_eV(708) = 1.442000d+1 ; td_cs_cm2(708) = 1.976367d-22
  td_E_eV(709) = 1.444000d+1 ; td_cs_cm2(709) = 1.977953d-22
  td_E_eV(710) = 1.446000d+1 ; td_cs_cm2(710) = 1.979538d-22
  td_E_eV(711) = 1.448000d+1 ; td_cs_cm2(711) = 1.981123d-22
  td_E_eV(712) = 1.450000d+1 ; td_cs_cm2(712) = 1.982706d-22
  td_E_eV(713) = 1.452000d+1 ; td_cs_cm2(713) = 1.984289d-22
  td_E_eV(714) = 1.454000d+1 ; td_cs_cm2(714) = 1.985871d-22
  td_E_eV(715) = 1.456000d+1 ; td_cs_cm2(715) = 1.987451d-22
  td_E_eV(716) = 1.458000d+1 ; td_cs_cm2(716) = 1.989031d-22
  td_E_eV(717) = 1.460000d+1 ; td_cs_cm2(717) = 1.990610d-22
  td_E_eV(718) = 1.462000d+1 ; td_cs_cm2(718) = 1.992188d-22
  td_E_eV(719) = 1.464000d+1 ; td_cs_cm2(719) = 1.993765d-22
  td_E_eV(720) = 1.466000d+1 ; td_cs_cm2(720) = 1.995342d-22
  td_E_eV(721) = 1.468000d+1 ; td_cs_cm2(721) = 1.996917d-22
  td_E_eV(722) = 1.470000d+1 ; td_cs_cm2(722) = 1.998491d-22
  td_E_eV(723) = 1.472000d+1 ; td_cs_cm2(723) = 2.000065d-22
  td_E_eV(724) = 1.474000d+1 ; td_cs_cm2(724) = 2.001637d-22
  td_E_eV(725) = 1.476000d+1 ; td_cs_cm2(725) = 2.003209d-22
  td_E_eV(726) = 1.478000d+1 ; td_cs_cm2(726) = 2.004780d-22
  td_E_eV(727) = 1.480000d+1 ; td_cs_cm2(727) = 2.006350d-22
  td_E_eV(728) = 1.482000d+1 ; td_cs_cm2(728) = 2.007919d-22
  td_E_eV(729) = 1.484000d+1 ; td_cs_cm2(729) = 2.009487d-22
  td_E_eV(730) = 1.486000d+1 ; td_cs_cm2(730) = 2.011054d-22
  td_E_eV(731) = 1.488000d+1 ; td_cs_cm2(731) = 2.012620d-22
  td_E_eV(732) = 1.490000d+1 ; td_cs_cm2(732) = 2.014186d-22
  td_E_eV(733) = 1.492000d+1 ; td_cs_cm2(733) = 2.015751d-22
  td_E_eV(734) = 1.494000d+1 ; td_cs_cm2(734) = 2.017314d-22
  td_E_eV(735) = 1.496000d+1 ; td_cs_cm2(735) = 2.018877d-22
  td_E_eV(736) = 1.498000d+1 ; td_cs_cm2(736) = 2.020439d-22
  td_E_eV(737) = 1.500000d+1 ; td_cs_cm2(737) = 2.022000d-22
  td_E_eV(738) = 2.500000d+1 ; td_cs_cm2(738) = 8.220000d-22

  td_cs_cm2 = 1.0d4 * td_cs_cm2 ! the original data are in m^2, convert to cm^2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_vib01
!
!-------------------------------------
!
real(8) function CSV_N2_vib01_m3s(E_eV)

  use N2_vib01
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_vib01_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_vib01_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_vib01_m3s










!---------------------------------------------------------------------------------------------------- #5
!#VIBRATIONAL
!#SPECIES: e / N2
!#PROCESS: E + N2 <-> E + N2*(v=0-2), Vibrational
!#PARAM.:  E = 0.573 eV, g1/g0 = 1, complete set
!#COMMENT: Vibrational excitation : N2(1Σg+:v=0) -> N2*(v=2) | Source: [Laporta et al, Plasma
!#COMMENT: Sources Sci. Technol. 23, 065002 (2014)].
!#UPDATED: 2023-12-23 21:19:17!
!
module N2_vib02
  integer, parameter :: N_td=722  !8   ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_vib02
!
!-------------------------------------
!
subroutine Prepare_N2_vib02

  use N2_vib02
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

!  td_E_eV(1) = 1.3_8; td_cs_cm2(1) = 1.11d-18
!  td_E_eV(2) = 1.6_8; td_cs_cm2(2) = 8.84d-17
!  td_E_eV(3) = 2.0_8; td_cs_cm2(3) = 1.01d-15
!  td_E_eV(4) = 2.4_8; td_cs_cm2(4) = 1.35d-15
!  td_E_eV(5) = 3.0_8; td_cs_cm2(5) = 5.86d-16
!  td_E_eV(6) = 4.0_8; td_cs_cm2(6) = 1.37d-16
!  td_E_eV(7) = 5.0_8; td_cs_cm2(7) = 2.34d-17
!  td_E_eV(8) = 6.0_8; td_cs_cm2(8) = 3.66d-18

  td_E_eV(  1) = 5.730000d-1 ; td_cs_cm2(  1) = 0.000000d+0
  td_E_eV(  2) = 5.900000d-1 ; td_cs_cm2(  2) = 3.582400d-26
  td_E_eV(  3) = 6.100000d-1 ; td_cs_cm2(  3) = 5.397000d-26
  td_E_eV(  4) = 6.300000d-1 ; td_cs_cm2(  4) = 6.810700d-26
  td_E_eV(  5) = 6.500000d-1 ; td_cs_cm2(  5) = 8.064300d-26
  td_E_eV(  6) = 6.700000d-1 ; td_cs_cm2(  6) = 9.250900d-26
  td_E_eV(  7) = 6.900000d-1 ; td_cs_cm2(  7) = 1.042300d-25
  td_E_eV(  8) = 7.100000d-1 ; td_cs_cm2(  8) = 1.161900d-25
  td_E_eV(  9) = 7.300000d-1 ; td_cs_cm2(  9) = 1.286800d-25
  td_E_eV( 10) = 7.500000d-1 ; td_cs_cm2( 10) = 1.420100d-25
  td_E_eV( 11) = 7.700000d-1 ; td_cs_cm2( 11) = 1.564500d-25
  td_E_eV( 12) = 7.900000d-1 ; td_cs_cm2( 12) = 1.723200d-25
  td_E_eV( 13) = 8.100000d-1 ; td_cs_cm2( 13) = 1.899600d-25
  td_E_eV( 14) = 8.300000d-1 ; td_cs_cm2( 14) = 2.097300d-25
  td_E_eV( 15) = 8.500000d-1 ; td_cs_cm2( 15) = 2.321000d-25
  td_E_eV( 16) = 8.700000d-1 ; td_cs_cm2( 16) = 2.575500d-25
  td_E_eV( 17) = 8.900000d-1 ; td_cs_cm2( 17) = 2.866900d-25
  td_E_eV( 18) = 9.100000d-1 ; td_cs_cm2( 18) = 3.202200d-25
  td_E_eV( 19) = 9.300000d-1 ; td_cs_cm2( 19) = 3.589700d-25
  td_E_eV( 20) = 9.500000d-1 ; td_cs_cm2( 20) = 4.039300d-25
  td_E_eV( 21) = 9.700000d-1 ; td_cs_cm2( 21) = 4.562800d-25
  td_E_eV( 22) = 9.900000d-1 ; td_cs_cm2( 22) = 5.174200d-25
  td_E_eV( 23) = 1.010000d+0 ; td_cs_cm2( 23) = 5.890600d-25
  td_E_eV( 24) = 1.030000d+0 ; td_cs_cm2( 24) = 6.732100d-25
  td_E_eV( 25) = 1.050000d+0 ; td_cs_cm2( 25) = 7.723400d-25
  td_E_eV( 26) = 1.070000d+0 ; td_cs_cm2( 26) = 8.893900d-25
  td_E_eV( 27) = 1.090000d+0 ; td_cs_cm2( 27) = 1.027900d-24
  td_E_eV( 28) = 1.110000d+0 ; td_cs_cm2( 28) = 1.192300d-24
  td_E_eV( 29) = 1.130000d+0 ; td_cs_cm2( 29) = 1.387900d-24
  td_E_eV( 30) = 1.150000d+0 ; td_cs_cm2( 30) = 1.621000d-24
  td_E_eV( 31) = 1.170000d+0 ; td_cs_cm2( 31) = 1.899500d-24
  td_E_eV( 32) = 1.190000d+0 ; td_cs_cm2( 32) = 2.233100d-24
  td_E_eV( 33) = 1.210000d+0 ; td_cs_cm2( 33) = 2.633600d-24
  td_E_eV( 34) = 1.230000d+0 ; td_cs_cm2( 34) = 3.115700d-24
  td_E_eV( 35) = 1.250000d+0 ; td_cs_cm2( 35) = 3.697500d-24
  td_E_eV( 36) = 1.270000d+0 ; td_cs_cm2( 36) = 4.401300d-24
  td_E_eV( 37) = 1.290000d+0 ; td_cs_cm2( 37) = 5.255100d-24
  td_E_eV( 38) = 1.310000d+0 ; td_cs_cm2( 38) = 6.293800d-24
  td_E_eV( 39) = 1.330000d+0 ; td_cs_cm2( 39) = 7.561300d-24
  td_E_eV( 40) = 1.350000d+0 ; td_cs_cm2( 40) = 9.112600d-24
  td_E_eV( 41) = 1.370000d+0 ; td_cs_cm2( 41) = 1.101800d-23
  td_E_eV( 42) = 1.390000d+0 ; td_cs_cm2( 42) = 1.336500d-23
  td_E_eV( 43) = 1.410000d+0 ; td_cs_cm2( 43) = 1.626900d-23
  td_E_eV( 44) = 1.430000d+0 ; td_cs_cm2( 44) = 1.987400d-23
  td_E_eV( 45) = 1.450000d+0 ; td_cs_cm2( 45) = 2.436900d-23
  td_E_eV( 46) = 1.470000d+0 ; td_cs_cm2( 46) = 2.999800d-23
  td_E_eV( 47) = 1.490000d+0 ; td_cs_cm2( 47) = 3.708100d-23
  td_E_eV( 48) = 1.510000d+0 ; td_cs_cm2( 48) = 4.603900d-23
  td_E_eV( 49) = 1.530000d+0 ; td_cs_cm2( 49) = 5.742900d-23
  td_E_eV( 50) = 1.550000d+0 ; td_cs_cm2( 50) = 7.199900d-23
  td_E_eV( 51) = 1.570000d+0 ; td_cs_cm2( 51) = 9.075400d-23
  td_E_eV( 52) = 1.590000d+0 ; td_cs_cm2( 52) = 1.150600d-22
  td_E_eV( 53) = 1.610000d+0 ; td_cs_cm2( 53) = 1.468100d-22
  td_E_eV( 54) = 1.630000d+0 ; td_cs_cm2( 54) = 1.886000d-22
  td_E_eV( 55) = 1.650000d+0 ; td_cs_cm2( 55) = 2.441000d-22
  td_E_eV( 56) = 1.670000d+0 ; td_cs_cm2( 56) = 3.185200d-22
  td_E_eV( 57) = 1.690000d+0 ; td_cs_cm2( 57) = 4.193500d-22
  td_E_eV( 58) = 1.710000d+0 ; td_cs_cm2( 58) = 5.575100d-22
  td_E_eV( 59) = 1.730000d+0 ; td_cs_cm2( 59) = 7.491500d-22
  td_E_eV( 60) = 1.750000d+0 ; td_cs_cm2( 60) = 1.018400d-21
  td_E_eV( 61) = 1.770000d+0 ; td_cs_cm2( 61) = 1.402000d-21
  td_E_eV( 62) = 1.790000d+0 ; td_cs_cm2( 62) = 1.955900d-21
  td_E_eV( 63) = 1.810000d+0 ; td_cs_cm2( 63) = 2.766300d-21
  td_E_eV( 64) = 1.830000d+0 ; td_cs_cm2( 64) = 3.964100d-21
  td_E_eV( 65) = 1.850000d+0 ; td_cs_cm2( 65) = 5.741200d-21
  td_E_eV( 66) = 1.870000d+0 ; td_cs_cm2( 66) = 8.355200d-21
  td_E_eV( 67) = 1.890000d+0 ; td_cs_cm2( 67) = 1.207800d-20
  td_E_eV( 68) = 1.910000d+0 ; td_cs_cm2( 68) = 1.701200d-20
  td_E_eV( 69) = 1.930000d+0 ; td_cs_cm2( 69) = 2.274400d-20
  td_E_eV( 70) = 1.950000d+0 ; td_cs_cm2( 70) = 2.815100d-20
  td_E_eV( 71) = 1.970000d+0 ; td_cs_cm2( 71) = 3.194300d-20
  td_E_eV( 72) = 1.990000d+0 ; td_cs_cm2( 72) = 3.362900d-20
  td_E_eV( 73) = 2.010000d+0 ; td_cs_cm2( 73) = 3.368600d-20
  td_E_eV( 74) = 2.030000d+0 ; td_cs_cm2( 74) = 3.290100d-20
  td_E_eV( 75) = 2.050000d+0 ; td_cs_cm2( 75) = 3.183800d-20
  td_E_eV( 76) = 2.070000d+0 ; td_cs_cm2( 76) = 3.071800d-20
  td_E_eV( 77) = 2.090000d+0 ; td_cs_cm2( 77) = 2.948100d-20
  td_E_eV( 78) = 2.110000d+0 ; td_cs_cm2( 78) = 2.782900d-20
  td_E_eV( 79) = 2.130000d+0 ; td_cs_cm2( 79) = 2.530000d-20
  td_E_eV( 80) = 2.150000d+0 ; td_cs_cm2( 80) = 2.146800d-20
  td_E_eV( 81) = 2.170000d+0 ; td_cs_cm2( 81) = 1.639500d-20
  td_E_eV( 82) = 2.190000d+0 ; td_cs_cm2( 82) = 1.099600d-20
  td_E_eV( 83) = 2.210000d+0 ; td_cs_cm2( 83) = 6.647600d-21
  td_E_eV( 84) = 2.230000d+0 ; td_cs_cm2( 84) = 4.265000d-21
  td_E_eV( 85) = 2.250000d+0 ; td_cs_cm2( 85) = 3.974400d-21
  td_E_eV( 86) = 2.270000d+0 ; td_cs_cm2( 86) = 5.495500d-21
  td_E_eV( 87) = 2.290000d+0 ; td_cs_cm2( 87) = 8.548500d-21
  td_E_eV( 88) = 2.310000d+0 ; td_cs_cm2( 88) = 1.298300d-20
  td_E_eV( 89) = 2.330000d+0 ; td_cs_cm2( 89) = 1.868500d-20
  td_E_eV( 90) = 2.350000d+0 ; td_cs_cm2( 90) = 2.532200d-20
  td_E_eV( 91) = 2.370000d+0 ; td_cs_cm2( 91) = 3.197600d-20
  td_E_eV( 92) = 2.390000d+0 ; td_cs_cm2( 92) = 3.696800d-20
  td_E_eV( 93) = 2.410000d+0 ; td_cs_cm2( 93) = 3.849600d-20
  td_E_eV( 94) = 2.430000d+0 ; td_cs_cm2( 94) = 3.597800d-20
  td_E_eV( 95) = 2.450000d+0 ; td_cs_cm2( 95) = 3.055000d-20
  td_E_eV( 96) = 2.470000d+0 ; td_cs_cm2( 96) = 2.402300d-20
  td_E_eV( 97) = 2.490000d+0 ; td_cs_cm2( 97) = 1.774900d-20
  td_E_eV( 98) = 2.510000d+0 ; td_cs_cm2( 98) = 1.238900d-20
  td_E_eV( 99) = 2.530000d+0 ; td_cs_cm2( 99) = 8.188400d-21
  td_E_eV(100) = 2.550000d+0 ; td_cs_cm2(100) = 5.287400d-21
  td_E_eV(101) = 2.570000d+0 ; td_cs_cm2(101) = 3.882500d-21
  td_E_eV(102) = 2.590000d+0 ; td_cs_cm2(102) = 4.213000d-21
  td_E_eV(103) = 2.610000d+0 ; td_cs_cm2(103) = 6.316100d-21
  td_E_eV(104) = 2.630000d+0 ; td_cs_cm2(104) = 9.675200d-21
  td_E_eV(105) = 2.650000d+0 ; td_cs_cm2(105) = 1.322900d-20
  td_E_eV(106) = 2.670000d+0 ; td_cs_cm2(106) = 1.596500d-20
  td_E_eV(107) = 2.690000d+0 ; td_cs_cm2(107) = 1.745000d-20
  td_E_eV(108) = 2.710000d+0 ; td_cs_cm2(108) = 1.777100d-20
  td_E_eV(109) = 2.730000d+0 ; td_cs_cm2(109) = 1.718800d-20
  td_E_eV(110) = 2.750000d+0 ; td_cs_cm2(110) = 1.591500d-20
  td_E_eV(111) = 2.770000d+0 ; td_cs_cm2(111) = 1.406800d-20
  td_E_eV(112) = 2.790000d+0 ; td_cs_cm2(112) = 1.171800d-20
  td_E_eV(113) = 2.810000d+0 ; td_cs_cm2(113) = 9.019100d-21
  td_E_eV(114) = 2.830000d+0 ; td_cs_cm2(114) = 6.335300d-21
  td_E_eV(115) = 2.850000d+0 ; td_cs_cm2(115) = 4.222000d-21
  td_E_eV(116) = 2.870000d+0 ; td_cs_cm2(116) = 3.141100d-21
  td_E_eV(117) = 2.890000d+0 ; td_cs_cm2(117) = 3.159800d-21
  td_E_eV(118) = 2.910000d+0 ; td_cs_cm2(118) = 3.988100d-21
  td_E_eV(119) = 2.930000d+0 ; td_cs_cm2(119) = 5.242800d-21
  td_E_eV(120) = 2.950000d+0 ; td_cs_cm2(120) = 6.618000d-21
  td_E_eV(121) = 2.970000d+0 ; td_cs_cm2(121) = 7.895800d-21
  td_E_eV(122) = 2.990000d+0 ; td_cs_cm2(122) = 8.885900d-21
  td_E_eV(123) = 3.010000d+0 ; td_cs_cm2(123) = 9.374700d-21
  td_E_eV(124) = 3.030000d+0 ; td_cs_cm2(124) = 9.141600d-21
  td_E_eV(125) = 3.050000d+0 ; td_cs_cm2(125) = 8.084000d-21
  td_E_eV(126) = 3.070000d+0 ; td_cs_cm2(126) = 6.390100d-21
  td_E_eV(127) = 3.090000d+0 ; td_cs_cm2(127) = 4.534800d-21
  td_E_eV(128) = 3.110000d+0 ; td_cs_cm2(128) = 3.011700d-21
  td_E_eV(129) = 3.130000d+0 ; td_cs_cm2(129) = 2.082900d-21
  td_E_eV(130) = 3.150000d+0 ; td_cs_cm2(130) = 1.766800d-21
  td_E_eV(131) = 3.170000d+0 ; td_cs_cm2(131) = 1.958600d-21
  td_E_eV(132) = 3.190000d+0 ; td_cs_cm2(132) = 2.519500d-21
  td_E_eV(133) = 3.210000d+0 ; td_cs_cm2(133) = 3.295100d-21
  td_E_eV(134) = 3.230000d+0 ; td_cs_cm2(134) = 4.094100d-21
  td_E_eV(135) = 3.250000d+0 ; td_cs_cm2(135) = 4.675500d-21
  td_E_eV(136) = 3.270000d+0 ; td_cs_cm2(136) = 4.810200d-21
  td_E_eV(137) = 3.290000d+0 ; td_cs_cm2(137) = 4.417600d-21
  td_E_eV(138) = 3.310000d+0 ; td_cs_cm2(138) = 3.642000d-21
  td_E_eV(139) = 3.330000d+0 ; td_cs_cm2(139) = 2.751600d-21
  td_E_eV(140) = 3.350000d+0 ; td_cs_cm2(140) = 1.975400d-21
  td_E_eV(141) = 3.370000d+0 ; td_cs_cm2(141) = 1.434700d-21
  td_E_eV(142) = 3.390000d+0 ; td_cs_cm2(142) = 1.165300d-21
  td_E_eV(143) = 3.410000d+0 ; td_cs_cm2(143) = 1.156700d-21
  td_E_eV(144) = 3.430000d+0 ; td_cs_cm2(144) = 1.366200d-21
  td_E_eV(145) = 3.450000d+0 ; td_cs_cm2(145) = 1.710200d-21
  td_E_eV(146) = 3.470000d+0 ; td_cs_cm2(146) = 2.055800d-21
  td_E_eV(147) = 3.490000d+0 ; td_cs_cm2(147) = 2.256600d-21
  td_E_eV(148) = 3.510000d+0 ; td_cs_cm2(148) = 2.230200d-21
  td_E_eV(149) = 3.530000d+0 ; td_cs_cm2(149) = 2.002400d-21
  td_E_eV(150) = 3.550000d+0 ; td_cs_cm2(150) = 1.668100d-21
  td_E_eV(151) = 3.570000d+0 ; td_cs_cm2(151) = 1.323300d-21
  td_E_eV(152) = 3.590000d+0 ; td_cs_cm2(152) = 1.034400d-21
  td_E_eV(153) = 3.610000d+0 ; td_cs_cm2(153) = 8.381800d-22
  td_E_eV(154) = 3.630000d+0 ; td_cs_cm2(154) = 7.516800d-22
  td_E_eV(155) = 3.650000d+0 ; td_cs_cm2(155) = 7.727200d-22
  td_E_eV(156) = 3.670000d+0 ; td_cs_cm2(156) = 8.727200d-22
  td_E_eV(157) = 3.690000d+0 ; td_cs_cm2(157) = 9.948100d-22
  td_E_eV(158) = 3.710000d+0 ; td_cs_cm2(158) = 1.075100d-21
  td_E_eV(159) = 3.730000d+0 ; td_cs_cm2(159) = 1.077500d-21
  td_E_eV(160) = 3.750000d+0 ; td_cs_cm2(160) = 1.006100d-21
  td_E_eV(161) = 3.770000d+0 ; td_cs_cm2(161) = 8.892300d-22
  td_E_eV(162) = 3.790000d+0 ; td_cs_cm2(162) = 7.581200d-22
  td_E_eV(163) = 3.810000d+0 ; td_cs_cm2(163) = 6.378000d-22
  td_E_eV(164) = 3.830000d+0 ; td_cs_cm2(164) = 5.464100d-22
  td_E_eV(165) = 3.850000d+0 ; td_cs_cm2(165) = 4.956200d-22
  td_E_eV(166) = 3.870000d+0 ; td_cs_cm2(166) = 4.877300d-22
  td_E_eV(167) = 3.890000d+0 ; td_cs_cm2(167) = 5.112200d-22
  td_E_eV(168) = 3.910000d+0 ; td_cs_cm2(168) = 5.427400d-22
  td_E_eV(169) = 3.930000d+0 ; td_cs_cm2(169) = 5.593700d-22
  td_E_eV(170) = 3.950000d+0 ; td_cs_cm2(170) = 5.505000d-22
  td_E_eV(171) = 3.970000d+0 ; td_cs_cm2(171) = 5.184200d-22
  td_E_eV(172) = 3.990000d+0 ; td_cs_cm2(172) = 4.718400d-22
  td_E_eV(173) = 4.010000d+0 ; td_cs_cm2(173) = 4.204900d-22
  td_E_eV(174) = 4.030000d+0 ; td_cs_cm2(174) = 3.731700d-22
  td_E_eV(175) = 4.050000d+0 ; td_cs_cm2(175) = 3.371100d-22
  td_E_eV(176) = 4.070000d+0 ; td_cs_cm2(176) = 3.167500d-22
  td_E_eV(177) = 4.090000d+0 ; td_cs_cm2(177) = 3.113300d-22
  td_E_eV(178) = 4.110000d+0 ; td_cs_cm2(178) = 3.143300d-22
  td_E_eV(179) = 4.130000d+0 ; td_cs_cm2(179) = 3.170200d-22
  td_E_eV(180) = 4.150000d+0 ; td_cs_cm2(180) = 3.134300d-22
  td_E_eV(181) = 4.170000d+0 ; td_cs_cm2(181) = 3.020400d-22
  td_E_eV(182) = 4.190000d+0 ; td_cs_cm2(182) = 2.844300d-22
  td_E_eV(183) = 4.210000d+0 ; td_cs_cm2(183) = 2.634000d-22
  td_E_eV(184) = 4.230000d+0 ; td_cs_cm2(184) = 2.420800d-22
  td_E_eV(185) = 4.250000d+0 ; td_cs_cm2(185) = 2.235000d-22
  td_E_eV(186) = 4.270000d+0 ; td_cs_cm2(186) = 2.099900d-22
  td_E_eV(187) = 4.290000d+0 ; td_cs_cm2(187) = 2.021000d-22
  td_E_eV(188) = 4.310000d+0 ; td_cs_cm2(188) = 1.981900d-22
  td_E_eV(189) = 4.330000d+0 ; td_cs_cm2(189) = 1.953700d-22
  td_E_eV(190) = 4.350000d+0 ; td_cs_cm2(190) = 1.913300d-22
  td_E_eV(191) = 4.370000d+0 ; td_cs_cm2(191) = 1.850800d-22
  td_E_eV(192) = 4.390000d+0 ; td_cs_cm2(192) = 1.767300d-22
  td_E_eV(193) = 4.410000d+0 ; td_cs_cm2(193) = 1.670100d-22
  td_E_eV(194) = 4.430000d+0 ; td_cs_cm2(194) = 1.569600d-22
  td_E_eV(195) = 4.450000d+0 ; td_cs_cm2(195) = 1.477200d-22
  td_E_eV(196) = 4.470000d+0 ; td_cs_cm2(196) = 1.402400d-22
  td_E_eV(197) = 4.490000d+0 ; td_cs_cm2(197) = 1.348000d-22
  td_E_eV(198) = 4.510000d+0 ; td_cs_cm2(198) = 1.308800d-22
  td_E_eV(199) = 4.530000d+0 ; td_cs_cm2(199) = 1.275100d-22
  td_E_eV(200) = 4.550000d+0 ; td_cs_cm2(200) = 1.239300d-22
  td_E_eV(201) = 4.570000d+0 ; td_cs_cm2(201) = 1.197400d-22
  td_E_eV(202) = 4.590000d+0 ; td_cs_cm2(202) = 1.149400d-22
  td_E_eV(203) = 4.610000d+0 ; td_cs_cm2(203) = 1.097500d-22
  td_E_eV(204) = 4.630000d+0 ; td_cs_cm2(204) = 1.045200d-22
  td_E_eV(205) = 4.650000d+0 ; td_cs_cm2(205) = 9.966700d-23
  td_E_eV(206) = 4.670000d+0 ; td_cs_cm2(206) = 9.552200d-23
  td_E_eV(207) = 4.690000d+0 ; td_cs_cm2(207) = 9.213600d-23
  td_E_eV(208) = 4.710000d+0 ; td_cs_cm2(208) = 8.928300d-23
  td_E_eV(209) = 4.730000d+0 ; td_cs_cm2(209) = 8.663600d-23
  td_E_eV(210) = 4.750000d+0 ; td_cs_cm2(210) = 8.394800d-23
  td_E_eV(211) = 4.770000d+0 ; td_cs_cm2(211) = 8.110900d-23
  td_E_eV(212) = 4.790000d+0 ; td_cs_cm2(212) = 7.812000d-23
  td_E_eV(213) = 4.810000d+0 ; td_cs_cm2(213) = 7.506500d-23
  td_E_eV(214) = 4.830000d+0 ; td_cs_cm2(214) = 7.208500d-23
  td_E_eV(215) = 4.850000d+0 ; td_cs_cm2(215) = 6.933400d-23
  td_E_eV(216) = 4.870000d+0 ; td_cs_cm2(216) = 6.689900d-23
  td_E_eV(217) = 4.890000d+0 ; td_cs_cm2(217) = 6.475400d-23
  td_E_eV(218) = 4.910000d+0 ; td_cs_cm2(218) = 6.279400d-23
  td_E_eV(219) = 4.930000d+0 ; td_cs_cm2(219) = 6.090700d-23
  td_E_eV(220) = 4.950000d+0 ; td_cs_cm2(220) = 5.902100d-23
  td_E_eV(221) = 4.970000d+0 ; td_cs_cm2(221) = 5.710800d-23
  td_E_eV(222) = 4.990000d+0 ; td_cs_cm2(222) = 5.518100d-23
  td_E_eV(223) = 5.010000d+0 ; td_cs_cm2(223) = 5.328000d-23
  td_E_eV(224) = 5.030000d+0 ; td_cs_cm2(224) = 5.146700d-23
  td_E_eV(225) = 5.050000d+0 ; td_cs_cm2(225) = 4.979000d-23
  td_E_eV(226) = 5.070000d+0 ; td_cs_cm2(226) = 4.825700d-23
  td_E_eV(227) = 5.090000d+0 ; td_cs_cm2(227) = 4.683900d-23
  td_E_eV(228) = 5.110000d+0 ; td_cs_cm2(228) = 4.549000d-23
  td_E_eV(229) = 5.130000d+0 ; td_cs_cm2(229) = 4.417400d-23
  td_E_eV(230) = 5.150000d+0 ; td_cs_cm2(230) = 4.287100d-23
  td_E_eV(231) = 5.170000d+0 ; td_cs_cm2(231) = 4.157700d-23
  td_E_eV(232) = 5.190000d+0 ; td_cs_cm2(232) = 4.030400d-23
  td_E_eV(233) = 5.210000d+0 ; td_cs_cm2(233) = 3.907500d-23
  td_E_eV(234) = 5.230000d+0 ; td_cs_cm2(234) = 3.791100d-23
  td_E_eV(235) = 5.250000d+0 ; td_cs_cm2(235) = 3.682200d-23
  td_E_eV(236) = 5.270000d+0 ; td_cs_cm2(236) = 3.580000d-23
  td_E_eV(237) = 5.290000d+0 ; td_cs_cm2(237) = 3.482600d-23
  td_E_eV(238) = 5.310000d+0 ; td_cs_cm2(238) = 3.388400d-23
  td_E_eV(239) = 5.330000d+0 ; td_cs_cm2(239) = 3.296300d-23
  td_E_eV(240) = 5.350000d+0 ; td_cs_cm2(240) = 3.205700d-23
  td_E_eV(241) = 5.370000d+0 ; td_cs_cm2(241) = 3.117000d-23
  td_E_eV(242) = 5.390000d+0 ; td_cs_cm2(242) = 3.030900d-23
  td_E_eV(243) = 5.410000d+0 ; td_cs_cm2(243) = 2.948400d-23
  td_E_eV(244) = 5.430000d+0 ; td_cs_cm2(244) = 2.870100d-23
  td_E_eV(245) = 5.450000d+0 ; td_cs_cm2(245) = 2.795600d-23
  td_E_eV(246) = 5.470000d+0 ; td_cs_cm2(246) = 2.724400d-23
  td_E_eV(247) = 5.490000d+0 ; td_cs_cm2(247) = 2.655700d-23
  td_E_eV(248) = 5.510000d+0 ; td_cs_cm2(248) = 2.588800d-23
  td_E_eV(249) = 5.530000d+0 ; td_cs_cm2(249) = 2.523400d-23
  td_E_eV(250) = 5.550000d+0 ; td_cs_cm2(250) = 2.459600d-23
  td_E_eV(251) = 5.570000d+0 ; td_cs_cm2(251) = 2.397500d-23
  td_E_eV(252) = 5.590000d+0 ; td_cs_cm2(252) = 2.337600d-23
  td_E_eV(253) = 5.610000d+0 ; td_cs_cm2(253) = 2.280200d-23
  td_E_eV(254) = 5.630000d+0 ; td_cs_cm2(254) = 2.225200d-23
  td_E_eV(255) = 5.650000d+0 ; td_cs_cm2(255) = 2.172400d-23
  td_E_eV(256) = 5.670000d+0 ; td_cs_cm2(256) = 2.121300d-23
  td_E_eV(257) = 5.690000d+0 ; td_cs_cm2(257) = 2.071700d-23
  td_E_eV(258) = 5.710000d+0 ; td_cs_cm2(258) = 2.023300d-23
  td_E_eV(259) = 5.730000d+0 ; td_cs_cm2(259) = 1.976100d-23
  td_E_eV(260) = 5.750000d+0 ; td_cs_cm2(260) = 1.930200d-23
  td_E_eV(261) = 5.770000d+0 ; td_cs_cm2(261) = 1.885700d-23
  td_E_eV(262) = 5.790000d+0 ; td_cs_cm2(262) = 1.842900d-23
  td_E_eV(263) = 5.810000d+0 ; td_cs_cm2(263) = 1.801600d-23
  td_E_eV(264) = 5.830000d+0 ; td_cs_cm2(264) = 1.761700d-23
  td_E_eV(265) = 5.850000d+0 ; td_cs_cm2(265) = 1.723100d-23
  td_E_eV(266) = 5.870000d+0 ; td_cs_cm2(266) = 1.685500d-23
  td_E_eV(267) = 5.890000d+0 ; td_cs_cm2(267) = 1.649000d-23
  td_E_eV(268) = 5.910000d+0 ; td_cs_cm2(268) = 1.613300d-23
  td_E_eV(269) = 5.930000d+0 ; td_cs_cm2(269) = 1.578600d-23
  td_E_eV(270) = 5.950000d+0 ; td_cs_cm2(270) = 1.545000d-23
  td_E_eV(271) = 5.970000d+0 ; td_cs_cm2(271) = 1.512300d-23
  td_E_eV(272) = 5.990000d+0 ; td_cs_cm2(272) = 1.480800d-23
  td_E_eV(273) = 6.010000d+0 ; td_cs_cm2(273) = 1.450200d-23
  td_E_eV(274) = 6.030000d+0 ; td_cs_cm2(274) = 1.420500d-23
  td_E_eV(275) = 6.050000d+0 ; td_cs_cm2(275) = 1.391700d-23
  td_E_eV(276) = 6.070000d+0 ; td_cs_cm2(276) = 1.363500d-23
  td_E_eV(277) = 6.090000d+0 ; td_cs_cm2(277) = 1.336100d-23
  td_E_eV(278) = 6.110000d+0 ; td_cs_cm2(278) = 1.309300d-23
  td_E_eV(279) = 6.130000d+0 ; td_cs_cm2(279) = 1.283300d-23
  td_E_eV(280) = 6.150000d+0 ; td_cs_cm2(280) = 1.258100d-23
  td_E_eV(281) = 6.170000d+0 ; td_cs_cm2(281) = 1.233600d-23
  td_E_eV(282) = 6.190000d+0 ; td_cs_cm2(282) = 1.209800d-23
  td_E_eV(283) = 6.210000d+0 ; td_cs_cm2(283) = 1.186600d-23
  td_E_eV(284) = 6.230000d+0 ; td_cs_cm2(284) = 1.164100d-23
  td_E_eV(285) = 6.250000d+0 ; td_cs_cm2(285) = 1.142000d-23
  td_E_eV(286) = 6.270000d+0 ; td_cs_cm2(286) = 1.120600d-23
  td_E_eV(287) = 6.290000d+0 ; td_cs_cm2(287) = 1.099600d-23
  td_E_eV(288) = 6.310000d+0 ; td_cs_cm2(288) = 1.079200d-23
  td_E_eV(289) = 6.330000d+0 ; td_cs_cm2(289) = 1.059400d-23
  td_E_eV(290) = 6.350000d+0 ; td_cs_cm2(290) = 1.040100d-23
  td_E_eV(291) = 6.370000d+0 ; td_cs_cm2(291) = 1.021200d-23
  td_E_eV(292) = 6.390000d+0 ; td_cs_cm2(292) = 1.002900d-23
  td_E_eV(293) = 6.410000d+0 ; td_cs_cm2(293) = 9.850100d-24
  td_E_eV(294) = 6.430000d+0 ; td_cs_cm2(294) = 9.675400d-24
  td_E_eV(295) = 6.450000d+0 ; td_cs_cm2(295) = 9.504800d-24
  td_E_eV(296) = 6.470000d+0 ; td_cs_cm2(296) = 9.338300d-24
  td_E_eV(297) = 6.490000d+0 ; td_cs_cm2(297) = 9.175900d-24
  td_E_eV(298) = 6.510000d+0 ; td_cs_cm2(298) = 9.017500d-24
  td_E_eV(299) = 6.530000d+0 ; td_cs_cm2(299) = 8.863000d-24
  td_E_eV(300) = 6.550000d+0 ; td_cs_cm2(300) = 8.712300d-24
  td_E_eV(301) = 6.570000d+0 ; td_cs_cm2(301) = 8.565200d-24
  td_E_eV(302) = 6.590000d+0 ; td_cs_cm2(302) = 8.421500d-24
  td_E_eV(303) = 6.610000d+0 ; td_cs_cm2(303) = 8.281000d-24
  td_E_eV(304) = 6.630000d+0 ; td_cs_cm2(304) = 8.143700d-24
  td_E_eV(305) = 6.650000d+0 ; td_cs_cm2(305) = 8.009500d-24
  td_E_eV(306) = 6.670000d+0 ; td_cs_cm2(306) = 7.878400d-24
  td_E_eV(307) = 6.690000d+0 ; td_cs_cm2(307) = 7.750400d-24
  td_E_eV(308) = 6.710000d+0 ; td_cs_cm2(308) = 7.625400d-24
  td_E_eV(309) = 6.730000d+0 ; td_cs_cm2(309) = 7.503200d-24
  td_E_eV(310) = 6.750000d+0 ; td_cs_cm2(310) = 7.383700d-24
  td_E_eV(311) = 6.770000d+0 ; td_cs_cm2(311) = 7.266800d-24
  td_E_eV(312) = 6.790000d+0 ; td_cs_cm2(312) = 7.152500d-24
  td_E_eV(313) = 6.810000d+0 ; td_cs_cm2(313) = 7.040600d-24
  td_E_eV(314) = 6.830000d+0 ; td_cs_cm2(314) = 6.931200d-24
  td_E_eV(315) = 6.850000d+0 ; td_cs_cm2(315) = 6.824200d-24
  td_E_eV(316) = 6.870000d+0 ; td_cs_cm2(316) = 6.719500d-24
  td_E_eV(317) = 6.890000d+0 ; td_cs_cm2(317) = 6.617100d-24
  td_E_eV(318) = 6.910000d+0 ; td_cs_cm2(318) = 6.516900d-24
  td_E_eV(319) = 6.930000d+0 ; td_cs_cm2(319) = 6.418800d-24
  td_E_eV(320) = 6.950000d+0 ; td_cs_cm2(320) = 6.322800d-24
  td_E_eV(321) = 6.970000d+0 ; td_cs_cm2(321) = 6.228700d-24
  td_E_eV(322) = 6.990000d+0 ; td_cs_cm2(322) = 6.136600d-24
  td_E_eV(323) = 7.010000d+0 ; td_cs_cm2(323) = 6.046400d-24
  td_E_eV(324) = 7.030000d+0 ; td_cs_cm2(324) = 5.958000d-24
  td_E_eV(325) = 7.050000d+0 ; td_cs_cm2(325) = 5.871500d-24
  td_E_eV(326) = 7.070000d+0 ; td_cs_cm2(326) = 5.786800d-24
  td_E_eV(327) = 7.090000d+0 ; td_cs_cm2(327) = 5.703800d-24
  td_E_eV(328) = 7.110000d+0 ; td_cs_cm2(328) = 5.622400d-24
  td_E_eV(329) = 7.130000d+0 ; td_cs_cm2(329) = 5.542700d-24
  td_E_eV(330) = 7.150000d+0 ; td_cs_cm2(330) = 5.464500d-24
  td_E_eV(331) = 7.170000d+0 ; td_cs_cm2(331) = 5.387900d-24
  td_E_eV(332) = 7.190000d+0 ; td_cs_cm2(332) = 5.312800d-24
  td_E_eV(333) = 7.210000d+0 ; td_cs_cm2(333) = 5.239200d-24
  td_E_eV(334) = 7.230000d+0 ; td_cs_cm2(334) = 5.167000d-24
  td_E_eV(335) = 7.250000d+0 ; td_cs_cm2(335) = 5.096200d-24
  td_E_eV(336) = 7.270000d+0 ; td_cs_cm2(336) = 5.026800d-24
  td_E_eV(337) = 7.290000d+0 ; td_cs_cm2(337) = 4.958700d-24
  td_E_eV(338) = 7.310000d+0 ; td_cs_cm2(338) = 4.891900d-24
  td_E_eV(339) = 7.330000d+0 ; td_cs_cm2(339) = 4.826400d-24
  td_E_eV(340) = 7.350000d+0 ; td_cs_cm2(340) = 4.762000d-24
  td_E_eV(341) = 7.370000d+0 ; td_cs_cm2(341) = 4.698900d-24
  td_E_eV(342) = 7.390000d+0 ; td_cs_cm2(342) = 4.637000d-24
  td_E_eV(343) = 7.410000d+0 ; td_cs_cm2(343) = 4.576300d-24
  td_E_eV(344) = 7.430000d+0 ; td_cs_cm2(344) = 4.516700d-24
  td_E_eV(345) = 7.450000d+0 ; td_cs_cm2(345) = 4.458100d-24
  td_E_eV(346) = 7.470000d+0 ; td_cs_cm2(346) = 4.400700d-24
  td_E_eV(347) = 7.490000d+0 ; td_cs_cm2(347) = 4.344200d-24
  td_E_eV(348) = 7.510000d+0 ; td_cs_cm2(348) = 4.288800d-24
  td_E_eV(349) = 7.530000d+0 ; td_cs_cm2(349) = 4.234400d-24
  td_E_eV(350) = 7.550000d+0 ; td_cs_cm2(350) = 4.180900d-24
  td_E_eV(351) = 7.570000d+0 ; td_cs_cm2(351) = 4.128500d-24
  td_E_eV(352) = 7.590000d+0 ; td_cs_cm2(352) = 4.076900d-24
  td_E_eV(353) = 7.610000d+0 ; td_cs_cm2(353) = 4.026300d-24
  td_E_eV(354) = 7.630000d+0 ; td_cs_cm2(354) = 3.976500d-24
  td_E_eV(355) = 7.650000d+0 ; td_cs_cm2(355) = 3.927600d-24
  td_E_eV(356) = 7.670000d+0 ; td_cs_cm2(356) = 3.879600d-24
  td_E_eV(357) = 7.690000d+0 ; td_cs_cm2(357) = 3.832400d-24
  td_E_eV(358) = 7.710000d+0 ; td_cs_cm2(358) = 3.786000d-24
  td_E_eV(359) = 7.730000d+0 ; td_cs_cm2(359) = 3.740400d-24
  td_E_eV(360) = 7.750000d+0 ; td_cs_cm2(360) = 3.695600d-24
  td_E_eV(361) = 7.770000d+0 ; td_cs_cm2(361) = 3.651500d-24
  td_E_eV(362) = 7.790000d+0 ; td_cs_cm2(362) = 3.608200d-24
  td_E_eV(363) = 7.810000d+0 ; td_cs_cm2(363) = 3.565600d-24
  td_E_eV(364) = 7.830000d+0 ; td_cs_cm2(364) = 3.523700d-24
  td_E_eV(365) = 7.850000d+0 ; td_cs_cm2(365) = 3.482600d-24
  td_E_eV(366) = 7.870000d+0 ; td_cs_cm2(366) = 3.442100d-24
  td_E_eV(367) = 7.890000d+0 ; td_cs_cm2(367) = 3.402300d-24
  td_E_eV(368) = 7.910000d+0 ; td_cs_cm2(368) = 3.363100d-24
  td_E_eV(369) = 7.930000d+0 ; td_cs_cm2(369) = 3.324600d-24
  td_E_eV(370) = 7.950000d+0 ; td_cs_cm2(370) = 3.286600d-24
  td_E_eV(371) = 7.970000d+0 ; td_cs_cm2(371) = 3.249300d-24
  td_E_eV(372) = 7.990000d+0 ; td_cs_cm2(372) = 3.212600d-24
  td_E_eV(373) = 8.010000d+0 ; td_cs_cm2(373) = 3.176500d-24
  td_E_eV(374) = 8.030000d+0 ; td_cs_cm2(374) = 3.141000d-24
  td_E_eV(375) = 8.050000d+0 ; td_cs_cm2(375) = 3.106000d-24
  td_E_eV(376) = 8.070000d+0 ; td_cs_cm2(376) = 3.071600d-24
  td_E_eV(377) = 8.090000d+0 ; td_cs_cm2(377) = 3.037700d-24
  td_E_eV(378) = 8.110000d+0 ; td_cs_cm2(378) = 3.004400d-24
  td_E_eV(379) = 8.130000d+0 ; td_cs_cm2(379) = 2.971600d-24
  td_E_eV(380) = 8.150000d+0 ; td_cs_cm2(380) = 2.939300d-24
  td_E_eV(381) = 8.170000d+0 ; td_cs_cm2(381) = 2.907500d-24
  td_E_eV(382) = 8.190000d+0 ; td_cs_cm2(382) = 2.876200d-24
  td_E_eV(383) = 8.210000d+0 ; td_cs_cm2(383) = 2.845300d-24
  td_E_eV(384) = 8.230000d+0 ; td_cs_cm2(384) = 2.814900d-24
  td_E_eV(385) = 8.250000d+0 ; td_cs_cm2(385) = 2.785000d-24
  td_E_eV(386) = 8.270000d+0 ; td_cs_cm2(386) = 2.755600d-24
  td_E_eV(387) = 8.290000d+0 ; td_cs_cm2(387) = 2.726600d-24
  td_E_eV(388) = 8.310000d+0 ; td_cs_cm2(388) = 2.698000d-24
  td_E_eV(389) = 8.330000d+0 ; td_cs_cm2(389) = 2.669800d-24
  td_E_eV(390) = 8.350000d+0 ; td_cs_cm2(390) = 2.642100d-24
  td_E_eV(391) = 8.370000d+0 ; td_cs_cm2(391) = 2.614800d-24
  td_E_eV(392) = 8.390000d+0 ; td_cs_cm2(392) = 2.587800d-24
  td_E_eV(393) = 8.410000d+0 ; td_cs_cm2(393) = 2.561300d-24
  td_E_eV(394) = 8.430000d+0 ; td_cs_cm2(394) = 2.535200d-24
  td_E_eV(395) = 8.450000d+0 ; td_cs_cm2(395) = 2.509400d-24
  td_E_eV(396) = 8.470000d+0 ; td_cs_cm2(396) = 2.484100d-24
  td_E_eV(397) = 8.490000d+0 ; td_cs_cm2(397) = 2.459000d-24
  td_E_eV(398) = 8.510000d+0 ; td_cs_cm2(398) = 2.434400d-24
  td_E_eV(399) = 8.530000d+0 ; td_cs_cm2(399) = 2.410100d-24
  td_E_eV(400) = 8.550000d+0 ; td_cs_cm2(400) = 2.386100d-24
  td_E_eV(401) = 8.570000d+0 ; td_cs_cm2(401) = 2.362500d-24
  td_E_eV(402) = 8.590000d+0 ; td_cs_cm2(402) = 2.339200d-24
  td_E_eV(403) = 8.610000d+0 ; td_cs_cm2(403) = 2.316300d-24
  td_E_eV(404) = 8.630000d+0 ; td_cs_cm2(404) = 2.293600d-24
  td_E_eV(405) = 8.650000d+0 ; td_cs_cm2(405) = 2.271300d-24
  td_E_eV(406) = 8.670000d+0 ; td_cs_cm2(406) = 2.249300d-24
  td_E_eV(407) = 8.690000d+0 ; td_cs_cm2(407) = 2.227600d-24
  td_E_eV(408) = 8.710000d+0 ; td_cs_cm2(408) = 2.206200d-24
  td_E_eV(409) = 8.730000d+0 ; td_cs_cm2(409) = 2.185100d-24
  td_E_eV(410) = 8.750000d+0 ; td_cs_cm2(410) = 2.164300d-24
  td_E_eV(411) = 8.770000d+0 ; td_cs_cm2(411) = 2.143700d-24
  td_E_eV(412) = 8.790000d+0 ; td_cs_cm2(412) = 2.123500d-24
  td_E_eV(413) = 8.810000d+0 ; td_cs_cm2(413) = 2.103500d-24
  td_E_eV(414) = 8.830000d+0 ; td_cs_cm2(414) = 2.083800d-24
  td_E_eV(415) = 8.850000d+0 ; td_cs_cm2(415) = 2.064300d-24
  td_E_eV(416) = 8.870000d+0 ; td_cs_cm2(416) = 2.045100d-24
  td_E_eV(417) = 8.890000d+0 ; td_cs_cm2(417) = 2.026200d-24
  td_E_eV(418) = 8.910000d+0 ; td_cs_cm2(418) = 2.007500d-24
  td_E_eV(419) = 8.930000d+0 ; td_cs_cm2(419) = 1.989100d-24
  td_E_eV(420) = 8.950000d+0 ; td_cs_cm2(420) = 1.970900d-24
  td_E_eV(421) = 8.970000d+0 ; td_cs_cm2(421) = 1.952900d-24
  td_E_eV(422) = 8.990000d+0 ; td_cs_cm2(422) = 1.935200d-24
  td_E_eV(423) = 9.010000d+0 ; td_cs_cm2(423) = 1.917700d-24
  td_E_eV(424) = 9.030000d+0 ; td_cs_cm2(424) = 1.900500d-24
  td_E_eV(425) = 9.050000d+0 ; td_cs_cm2(425) = 1.883400d-24
  td_E_eV(426) = 9.070000d+0 ; td_cs_cm2(426) = 1.866600d-24
  td_E_eV(427) = 9.090000d+0 ; td_cs_cm2(427) = 1.850000d-24
  td_E_eV(428) = 9.110000d+0 ; td_cs_cm2(428) = 1.833600d-24
  td_E_eV(429) = 9.130000d+0 ; td_cs_cm2(429) = 1.817400d-24
  td_E_eV(430) = 9.150000d+0 ; td_cs_cm2(430) = 1.801400d-24
  td_E_eV(431) = 9.170000d+0 ; td_cs_cm2(431) = 1.785600d-24
  td_E_eV(432) = 9.190000d+0 ; td_cs_cm2(432) = 1.770100d-24
  td_E_eV(433) = 9.210000d+0 ; td_cs_cm2(433) = 1.754700d-24
  td_E_eV(434) = 9.230000d+0 ; td_cs_cm2(434) = 1.739500d-24
  td_E_eV(435) = 9.250000d+0 ; td_cs_cm2(435) = 1.724500d-24
  td_E_eV(436) = 9.270000d+0 ; td_cs_cm2(436) = 1.709700d-24
  td_E_eV(437) = 9.290000d+0 ; td_cs_cm2(437) = 1.695000d-24
  td_E_eV(438) = 9.310000d+0 ; td_cs_cm2(438) = 1.680600d-24
  td_E_eV(439) = 9.330000d+0 ; td_cs_cm2(439) = 1.666300d-24
  td_E_eV(440) = 9.350000d+0 ; td_cs_cm2(440) = 1.652200d-24
  td_E_eV(441) = 9.370000d+0 ; td_cs_cm2(441) = 1.638300d-24
  td_E_eV(442) = 9.390000d+0 ; td_cs_cm2(442) = 1.624500d-24
  td_E_eV(443) = 9.410000d+0 ; td_cs_cm2(443) = 1.610900d-24
  td_E_eV(444) = 9.430000d+0 ; td_cs_cm2(444) = 1.597500d-24
  td_E_eV(445) = 9.450000d+0 ; td_cs_cm2(445) = 1.584200d-24
  td_E_eV(446) = 9.470000d+0 ; td_cs_cm2(446) = 1.571100d-24
  td_E_eV(447) = 9.490000d+0 ; td_cs_cm2(447) = 1.558100d-24
  td_E_eV(448) = 9.510000d+0 ; td_cs_cm2(448) = 1.545300d-24
  td_E_eV(449) = 9.530000d+0 ; td_cs_cm2(449) = 1.532700d-24
  td_E_eV(450) = 9.550000d+0 ; td_cs_cm2(450) = 1.520100d-24
  td_E_eV(451) = 9.570000d+0 ; td_cs_cm2(451) = 1.507800d-24
  td_E_eV(452) = 9.590000d+0 ; td_cs_cm2(452) = 1.495600d-24
  td_E_eV(453) = 9.610000d+0 ; td_cs_cm2(453) = 1.483500d-24
  td_E_eV(454) = 9.630000d+0 ; td_cs_cm2(454) = 1.471600d-24
  td_E_eV(455) = 9.650000d+0 ; td_cs_cm2(455) = 1.459800d-24
  td_E_eV(456) = 9.670000d+0 ; td_cs_cm2(456) = 1.448100d-24
  td_E_eV(457) = 9.690000d+0 ; td_cs_cm2(457) = 1.436600d-24
  td_E_eV(458) = 9.710000d+0 ; td_cs_cm2(458) = 1.425200d-24
  td_E_eV(459) = 9.730000d+0 ; td_cs_cm2(459) = 1.414000d-24
  td_E_eV(460) = 9.750000d+0 ; td_cs_cm2(460) = 1.402800d-24
  td_E_eV(461) = 9.770000d+0 ; td_cs_cm2(461) = 1.391800d-24
  td_E_eV(462) = 9.790000d+0 ; td_cs_cm2(462) = 1.381000d-24
  td_E_eV(463) = 9.810000d+0 ; td_cs_cm2(463) = 1.370200d-24
  td_E_eV(464) = 9.830000d+0 ; td_cs_cm2(464) = 1.359600d-24
  td_E_eV(465) = 9.850000d+0 ; td_cs_cm2(465) = 1.349000d-24
  td_E_eV(466) = 9.870000d+0 ; td_cs_cm2(466) = 1.338600d-24
  td_E_eV(467) = 9.890000d+0 ; td_cs_cm2(467) = 1.328400d-24
  td_E_eV(468) = 9.910000d+0 ; td_cs_cm2(468) = 1.318200d-24
  td_E_eV(469) = 9.930000d+0 ; td_cs_cm2(469) = 1.308100d-24
  td_E_eV(470) = 9.950000d+0 ; td_cs_cm2(470) = 1.298200d-24
  td_E_eV(471) = 9.970000d+0 ; td_cs_cm2(471) = 1.288300d-24
  td_E_eV(472) = 9.990000d+0 ; td_cs_cm2(472) = 1.278600d-24
  td_E_eV(473) = 1.001000d+1 ; td_cs_cm2(473) = 1.269000d-24
  td_E_eV(474) = 1.003000d+1 ; td_cs_cm2(474) = 1.259500d-24
  td_E_eV(475) = 1.005000d+1 ; td_cs_cm2(475) = 1.250100d-24
  td_E_eV(476) = 1.007000d+1 ; td_cs_cm2(476) = 1.240700d-24
  td_E_eV(477) = 1.009000d+1 ; td_cs_cm2(477) = 1.231500d-24
  td_E_eV(478) = 1.011000d+1 ; td_cs_cm2(478) = 1.222400d-24
  td_E_eV(479) = 1.013000d+1 ; td_cs_cm2(479) = 1.213400d-24
  td_E_eV(480) = 1.015000d+1 ; td_cs_cm2(480) = 1.204500d-24
  td_E_eV(481) = 1.017000d+1 ; td_cs_cm2(481) = 1.195700d-24
  td_E_eV(482) = 1.019000d+1 ; td_cs_cm2(482) = 1.186900d-24
  td_E_eV(483) = 1.021000d+1 ; td_cs_cm2(483) = 1.178300d-24
  td_E_eV(484) = 1.023000d+1 ; td_cs_cm2(484) = 1.169700d-24
  td_E_eV(485) = 1.025000d+1 ; td_cs_cm2(485) = 1.161300d-24
  td_E_eV(486) = 1.027000d+1 ; td_cs_cm2(486) = 1.152900d-24
  td_E_eV(487) = 1.029000d+1 ; td_cs_cm2(487) = 1.144600d-24
  td_E_eV(488) = 1.031000d+1 ; td_cs_cm2(488) = 1.136400d-24
  td_E_eV(489) = 1.033000d+1 ; td_cs_cm2(489) = 1.128300d-24
  td_E_eV(490) = 1.035000d+1 ; td_cs_cm2(490) = 1.120300d-24
  td_E_eV(491) = 1.037000d+1 ; td_cs_cm2(491) = 1.112300d-24
  td_E_eV(492) = 1.039000d+1 ; td_cs_cm2(492) = 1.104500d-24
  td_E_eV(493) = 1.041000d+1 ; td_cs_cm2(493) = 1.096700d-24
  td_E_eV(494) = 1.043000d+1 ; td_cs_cm2(494) = 1.089000d-24
  td_E_eV(495) = 1.045000d+1 ; td_cs_cm2(495) = 1.081300d-24
  td_E_eV(496) = 1.047000d+1 ; td_cs_cm2(496) = 1.073800d-24
  td_E_eV(497) = 1.049000d+1 ; td_cs_cm2(497) = 1.066300d-24
  td_E_eV(498) = 1.051000d+1 ; td_cs_cm2(498) = 1.058900d-24
  td_E_eV(499) = 1.053000d+1 ; td_cs_cm2(499) = 1.051600d-24
  td_E_eV(500) = 1.055000d+1 ; td_cs_cm2(500) = 1.044300d-24
  td_E_eV(501) = 1.057000d+1 ; td_cs_cm2(501) = 1.037200d-24
  td_E_eV(502) = 1.059000d+1 ; td_cs_cm2(502) = 1.030100d-24
  td_E_eV(503) = 1.061000d+1 ; td_cs_cm2(503) = 1.023000d-24
  td_E_eV(504) = 1.063000d+1 ; td_cs_cm2(504) = 1.016100d-24
  td_E_eV(505) = 1.065000d+1 ; td_cs_cm2(505) = 1.009200d-24
  td_E_eV(506) = 1.067000d+1 ; td_cs_cm2(506) = 1.002300d-24
  td_E_eV(507) = 1.069000d+1 ; td_cs_cm2(507) = 9.955700d-25
  td_E_eV(508) = 1.071000d+1 ; td_cs_cm2(508) = 9.888700d-25
  td_E_eV(509) = 1.073000d+1 ; td_cs_cm2(509) = 9.822400d-25
  td_E_eV(510) = 1.075000d+1 ; td_cs_cm2(510) = 9.756700d-25
  td_E_eV(511) = 1.077000d+1 ; td_cs_cm2(511) = 9.691700d-25
  td_E_eV(512) = 1.079000d+1 ; td_cs_cm2(512) = 9.627300d-25
  td_E_eV(513) = 1.081000d+1 ; td_cs_cm2(513) = 9.563500d-25
  td_E_eV(514) = 1.083000d+1 ; td_cs_cm2(514) = 9.500300d-25
  td_E_eV(515) = 1.085000d+1 ; td_cs_cm2(515) = 9.437800d-25
  td_E_eV(516) = 1.087000d+1 ; td_cs_cm2(516) = 9.375800d-25
  td_E_eV(517) = 1.089000d+1 ; td_cs_cm2(517) = 9.314400d-25
  td_E_eV(518) = 1.091000d+1 ; td_cs_cm2(518) = 9.253600d-25
  td_E_eV(519) = 1.093000d+1 ; td_cs_cm2(519) = 9.193400d-25
  td_E_eV(520) = 1.095000d+1 ; td_cs_cm2(520) = 9.133800d-25
  td_E_eV(521) = 1.097000d+1 ; td_cs_cm2(521) = 9.074700d-25
  td_E_eV(522) = 1.099000d+1 ; td_cs_cm2(522) = 9.016200d-25
  td_E_eV(523) = 1.101000d+1 ; td_cs_cm2(523) = 8.958200d-25
  td_E_eV(524) = 1.103000d+1 ; td_cs_cm2(524) = 8.900700d-25
  td_E_eV(525) = 1.105000d+1 ; td_cs_cm2(525) = 8.843800d-25
  td_E_eV(526) = 1.107000d+1 ; td_cs_cm2(526) = 8.787500d-25
  td_E_eV(527) = 1.109000d+1 ; td_cs_cm2(527) = 8.731600d-25
  td_E_eV(528) = 1.111000d+1 ; td_cs_cm2(528) = 8.676300d-25
  td_E_eV(529) = 1.113000d+1 ; td_cs_cm2(529) = 8.621400d-25
  td_E_eV(530) = 1.115000d+1 ; td_cs_cm2(530) = 8.567100d-25
  td_E_eV(531) = 1.117000d+1 ; td_cs_cm2(531) = 8.513300d-25
  td_E_eV(532) = 1.119000d+1 ; td_cs_cm2(532) = 8.460000d-25
  td_E_eV(533) = 1.121000d+1 ; td_cs_cm2(533) = 8.407100d-25
  td_E_eV(534) = 1.123000d+1 ; td_cs_cm2(534) = 8.354700d-25
  td_E_eV(535) = 1.125000d+1 ; td_cs_cm2(535) = 8.302800d-25
  td_E_eV(536) = 1.127000d+1 ; td_cs_cm2(536) = 8.251400d-25
  td_E_eV(537) = 1.129000d+1 ; td_cs_cm2(537) = 8.200400d-25
  td_E_eV(538) = 1.131000d+1 ; td_cs_cm2(538) = 8.149900d-25
  td_E_eV(539) = 1.133000d+1 ; td_cs_cm2(539) = 8.099900d-25
  td_E_eV(540) = 1.135000d+1 ; td_cs_cm2(540) = 8.050300d-25
  td_E_eV(541) = 1.137000d+1 ; td_cs_cm2(541) = 8.001100d-25
  td_E_eV(542) = 1.139000d+1 ; td_cs_cm2(542) = 7.952400d-25
  td_E_eV(543) = 1.141000d+1 ; td_cs_cm2(543) = 7.904100d-25
  td_E_eV(544) = 1.143000d+1 ; td_cs_cm2(544) = 7.856200d-25
  td_E_eV(545) = 1.145000d+1 ; td_cs_cm2(545) = 7.808800d-25
  td_E_eV(546) = 1.147000d+1 ; td_cs_cm2(546) = 7.761700d-25
  td_E_eV(547) = 1.149000d+1 ; td_cs_cm2(547) = 7.715100d-25
  td_E_eV(548) = 1.151000d+1 ; td_cs_cm2(548) = 7.668900d-25
  td_E_eV(549) = 1.153000d+1 ; td_cs_cm2(549) = 7.623100d-25
  td_E_eV(550) = 1.155000d+1 ; td_cs_cm2(550) = 7.577700d-25
  td_E_eV(551) = 1.157000d+1 ; td_cs_cm2(551) = 7.532700d-25
  td_E_eV(552) = 1.159000d+1 ; td_cs_cm2(552) = 7.488000d-25
  td_E_eV(553) = 1.161000d+1 ; td_cs_cm2(553) = 7.443800d-25
  td_E_eV(554) = 1.163000d+1 ; td_cs_cm2(554) = 7.399900d-25
  td_E_eV(555) = 1.165000d+1 ; td_cs_cm2(555) = 7.356500d-25
  td_E_eV(556) = 1.167000d+1 ; td_cs_cm2(556) = 7.313300d-25
  td_E_eV(557) = 1.169000d+1 ; td_cs_cm2(557) = 7.270600d-25
  td_E_eV(558) = 1.171000d+1 ; td_cs_cm2(558) = 7.228200d-25
  td_E_eV(559) = 1.173000d+1 ; td_cs_cm2(559) = 7.186200d-25
  td_E_eV(560) = 1.175000d+1 ; td_cs_cm2(560) = 7.144500d-25
  td_E_eV(561) = 1.177000d+1 ; td_cs_cm2(561) = 7.103200d-25
  td_E_eV(562) = 1.179000d+1 ; td_cs_cm2(562) = 7.062300d-25
  td_E_eV(563) = 1.181000d+1 ; td_cs_cm2(563) = 7.021600d-25
  td_E_eV(564) = 1.183000d+1 ; td_cs_cm2(564) = 6.981400d-25
  td_E_eV(565) = 1.185000d+1 ; td_cs_cm2(565) = 6.941400d-25
  td_E_eV(566) = 1.187000d+1 ; td_cs_cm2(566) = 6.901800d-25
  td_E_eV(567) = 1.189000d+1 ; td_cs_cm2(567) = 6.862500d-25
  td_E_eV(568) = 1.191000d+1 ; td_cs_cm2(568) = 6.823600d-25
  td_E_eV(569) = 1.193000d+1 ; td_cs_cm2(569) = 6.784900d-25
  td_E_eV(570) = 1.195000d+1 ; td_cs_cm2(570) = 6.746600d-25
  td_E_eV(571) = 1.197000d+1 ; td_cs_cm2(571) = 6.708600d-25
  td_E_eV(572) = 1.199000d+1 ; td_cs_cm2(572) = 6.670900d-25
  td_E_eV(573) = 1.201000d+1 ; td_cs_cm2(573) = 6.633600d-25
  td_E_eV(574) = 1.203000d+1 ; td_cs_cm2(574) = 6.596500d-25
  td_E_eV(575) = 1.205000d+1 ; td_cs_cm2(575) = 6.559700d-25
  td_E_eV(576) = 1.207000d+1 ; td_cs_cm2(576) = 6.523300d-25
  td_E_eV(577) = 1.209000d+1 ; td_cs_cm2(577) = 6.487100d-25
  td_E_eV(578) = 1.211000d+1 ; td_cs_cm2(578) = 6.451200d-25
  td_E_eV(579) = 1.213000d+1 ; td_cs_cm2(579) = 6.415600d-25
  td_E_eV(580) = 1.215000d+1 ; td_cs_cm2(580) = 6.380300d-25
  td_E_eV(581) = 1.217000d+1 ; td_cs_cm2(581) = 6.345300d-25
  td_E_eV(582) = 1.219000d+1 ; td_cs_cm2(582) = 6.310600d-25
  td_E_eV(583) = 1.221000d+1 ; td_cs_cm2(583) = 6.276100d-25
  td_E_eV(584) = 1.223000d+1 ; td_cs_cm2(584) = 6.241900d-25
  td_E_eV(585) = 1.225000d+1 ; td_cs_cm2(585) = 6.208000d-25
  td_E_eV(586) = 1.227000d+1 ; td_cs_cm2(586) = 6.174300d-25
  td_E_eV(587) = 1.229000d+1 ; td_cs_cm2(587) = 6.141000d-25
  td_E_eV(588) = 1.231000d+1 ; td_cs_cm2(588) = 6.107800d-25
  td_E_eV(589) = 1.233000d+1 ; td_cs_cm2(589) = 6.075000d-25
  td_E_eV(590) = 1.235000d+1 ; td_cs_cm2(590) = 6.042400d-25
  td_E_eV(591) = 1.237000d+1 ; td_cs_cm2(591) = 6.010000d-25
  td_E_eV(592) = 1.239000d+1 ; td_cs_cm2(592) = 5.977900d-25
  td_E_eV(593) = 1.241000d+1 ; td_cs_cm2(593) = 5.946100d-25
  td_E_eV(594) = 1.243000d+1 ; td_cs_cm2(594) = 5.914500d-25
  td_E_eV(595) = 1.245000d+1 ; td_cs_cm2(595) = 5.883200d-25
  td_E_eV(596) = 1.247000d+1 ; td_cs_cm2(596) = 5.852100d-25
  td_E_eV(597) = 1.249000d+1 ; td_cs_cm2(597) = 5.821200d-25
  td_E_eV(598) = 1.251000d+1 ; td_cs_cm2(598) = 5.790600d-25
  td_E_eV(599) = 1.253000d+1 ; td_cs_cm2(599) = 5.760200d-25
  td_E_eV(600) = 1.255000d+1 ; td_cs_cm2(600) = 5.730000d-25
  td_E_eV(601) = 1.257000d+1 ; td_cs_cm2(601) = 5.700100d-25
  td_E_eV(602) = 1.259000d+1 ; td_cs_cm2(602) = 5.670400d-25
  td_E_eV(603) = 1.261000d+1 ; td_cs_cm2(603) = 5.640900d-25
  td_E_eV(604) = 1.263000d+1 ; td_cs_cm2(604) = 5.611700d-25
  td_E_eV(605) = 1.265000d+1 ; td_cs_cm2(605) = 5.582700d-25
  td_E_eV(606) = 1.267000d+1 ; td_cs_cm2(606) = 5.553800d-25
  td_E_eV(607) = 1.269000d+1 ; td_cs_cm2(607) = 5.525300d-25
  td_E_eV(608) = 1.271000d+1 ; td_cs_cm2(608) = 5.496900d-25
  td_E_eV(609) = 1.273000d+1 ; td_cs_cm2(609) = 5.468700d-25
  td_E_eV(610) = 1.275000d+1 ; td_cs_cm2(610) = 5.440800d-25
  td_E_eV(611) = 1.277000d+1 ; td_cs_cm2(611) = 5.413000d-25
  td_E_eV(612) = 1.279000d+1 ; td_cs_cm2(612) = 5.385500d-25
  td_E_eV(613) = 1.281000d+1 ; td_cs_cm2(613) = 5.358200d-25
  td_E_eV(614) = 1.283000d+1 ; td_cs_cm2(614) = 5.331000d-25
  td_E_eV(615) = 1.285000d+1 ; td_cs_cm2(615) = 5.304100d-25
  td_E_eV(616) = 1.287000d+1 ; td_cs_cm2(616) = 5.277400d-25
  td_E_eV(617) = 1.289000d+1 ; td_cs_cm2(617) = 5.250900d-25
  td_E_eV(618) = 1.291000d+1 ; td_cs_cm2(618) = 5.224500d-25
  td_E_eV(619) = 1.293000d+1 ; td_cs_cm2(619) = 5.198400d-25
  td_E_eV(620) = 1.295000d+1 ; td_cs_cm2(620) = 5.172400d-25
  td_E_eV(621) = 1.297000d+1 ; td_cs_cm2(621) = 5.146700d-25
  td_E_eV(622) = 1.299000d+1 ; td_cs_cm2(622) = 5.121100d-25
  td_E_eV(623) = 1.301000d+1 ; td_cs_cm2(623) = 5.095700d-25
  td_E_eV(624) = 1.303000d+1 ; td_cs_cm2(624) = 5.070500d-25
  td_E_eV(625) = 1.305000d+1 ; td_cs_cm2(625) = 5.045500d-25
  td_E_eV(626) = 1.307000d+1 ; td_cs_cm2(626) = 5.020700d-25
  td_E_eV(627) = 1.309000d+1 ; td_cs_cm2(627) = 4.996000d-25
  td_E_eV(628) = 1.311000d+1 ; td_cs_cm2(628) = 4.971500d-25
  td_E_eV(629) = 1.313000d+1 ; td_cs_cm2(629) = 4.947200d-25
  td_E_eV(630) = 1.315000d+1 ; td_cs_cm2(630) = 4.923100d-25
  td_E_eV(631) = 1.317000d+1 ; td_cs_cm2(631) = 4.899100d-25
  td_E_eV(632) = 1.319000d+1 ; td_cs_cm2(632) = 4.875300d-25
  td_E_eV(633) = 1.321000d+1 ; td_cs_cm2(633) = 4.851700d-25
  td_E_eV(634) = 1.323000d+1 ; td_cs_cm2(634) = 4.828300d-25
  td_E_eV(635) = 1.325000d+1 ; td_cs_cm2(635) = 4.805000d-25
  td_E_eV(636) = 1.327000d+1 ; td_cs_cm2(636) = 4.781800d-25
  td_E_eV(637) = 1.329000d+1 ; td_cs_cm2(637) = 4.758900d-25
  td_E_eV(638) = 1.331000d+1 ; td_cs_cm2(638) = 4.736100d-25
  td_E_eV(639) = 1.333000d+1 ; td_cs_cm2(639) = 4.713500d-25
  td_E_eV(640) = 1.335000d+1 ; td_cs_cm2(640) = 4.691000d-25
  td_E_eV(641) = 1.337000d+1 ; td_cs_cm2(641) = 4.668700d-25
  td_E_eV(642) = 1.339000d+1 ; td_cs_cm2(642) = 4.646500d-25
  td_E_eV(643) = 1.341000d+1 ; td_cs_cm2(643) = 4.624500d-25
  td_E_eV(644) = 1.343000d+1 ; td_cs_cm2(644) = 4.602600d-25
  td_E_eV(645) = 1.345000d+1 ; td_cs_cm2(645) = 4.580900d-25
  td_E_eV(646) = 1.347000d+1 ; td_cs_cm2(646) = 4.559400d-25
  td_E_eV(647) = 1.349000d+1 ; td_cs_cm2(647) = 4.537900d-25
  td_E_eV(648) = 1.351000d+1 ; td_cs_cm2(648) = 4.516700d-25
  td_E_eV(649) = 1.353000d+1 ; td_cs_cm2(649) = 4.495600d-25
  td_E_eV(650) = 1.355000d+1 ; td_cs_cm2(650) = 4.474600d-25
  td_E_eV(651) = 1.357000d+1 ; td_cs_cm2(651) = 4.453800d-25
  td_E_eV(652) = 1.359000d+1 ; td_cs_cm2(652) = 4.433100d-25
  td_E_eV(653) = 1.361000d+1 ; td_cs_cm2(653) = 4.412500d-25
  td_E_eV(654) = 1.363000d+1 ; td_cs_cm2(654) = 4.392100d-25
  td_E_eV(655) = 1.365000d+1 ; td_cs_cm2(655) = 4.371900d-25
  td_E_eV(656) = 1.367000d+1 ; td_cs_cm2(656) = 4.351700d-25
  td_E_eV(657) = 1.369000d+1 ; td_cs_cm2(657) = 4.331700d-25
  td_E_eV(658) = 1.371000d+1 ; td_cs_cm2(658) = 4.311900d-25
  td_E_eV(659) = 1.373000d+1 ; td_cs_cm2(659) = 4.292200d-25
  td_E_eV(660) = 1.375000d+1 ; td_cs_cm2(660) = 4.272600d-25
  td_E_eV(661) = 1.377000d+1 ; td_cs_cm2(661) = 4.253100d-25
  td_E_eV(662) = 1.379000d+1 ; td_cs_cm2(662) = 4.233800d-25
  td_E_eV(663) = 1.381000d+1 ; td_cs_cm2(663) = 4.214600d-25
  td_E_eV(664) = 1.383000d+1 ; td_cs_cm2(664) = 4.195500d-25
  td_E_eV(665) = 1.385000d+1 ; td_cs_cm2(665) = 4.176600d-25
  td_E_eV(666) = 1.387000d+1 ; td_cs_cm2(666) = 4.157700d-25
  td_E_eV(667) = 1.389000d+1 ; td_cs_cm2(667) = 4.139000d-25
  td_E_eV(668) = 1.391000d+1 ; td_cs_cm2(668) = 4.120500d-25
  td_E_eV(669) = 1.393000d+1 ; td_cs_cm2(669) = 4.102000d-25
  td_E_eV(670) = 1.395000d+1 ; td_cs_cm2(670) = 4.083700d-25
  td_E_eV(671) = 1.397000d+1 ; td_cs_cm2(671) = 4.065500d-25
  td_E_eV(672) = 1.399000d+1 ; td_cs_cm2(672) = 4.047400d-25
  td_E_eV(673) = 1.401000d+1 ; td_cs_cm2(673) = 4.029400d-25
  td_E_eV(674) = 1.403000d+1 ; td_cs_cm2(674) = 4.011600d-25
  td_E_eV(675) = 1.405000d+1 ; td_cs_cm2(675) = 3.993800d-25
  td_E_eV(676) = 1.407000d+1 ; td_cs_cm2(676) = 3.976200d-25
  td_E_eV(677) = 1.409000d+1 ; td_cs_cm2(677) = 3.958700d-25
  td_E_eV(678) = 1.411000d+1 ; td_cs_cm2(678) = 3.941300d-25
  td_E_eV(679) = 1.413000d+1 ; td_cs_cm2(679) = 3.924000d-25
  td_E_eV(680) = 1.415000d+1 ; td_cs_cm2(680) = 3.906800d-25
  td_E_eV(681) = 1.417000d+1 ; td_cs_cm2(681) = 3.889800d-25
  td_E_eV(682) = 1.419000d+1 ; td_cs_cm2(682) = 3.872800d-25
  td_E_eV(683) = 1.421000d+1 ; td_cs_cm2(683) = 3.856000d-25
  td_E_eV(684) = 1.423000d+1 ; td_cs_cm2(684) = 3.839200d-25
  td_E_eV(685) = 1.425000d+1 ; td_cs_cm2(685) = 3.822600d-25
  td_E_eV(686) = 1.427000d+1 ; td_cs_cm2(686) = 3.806100d-25
  td_E_eV(687) = 1.429000d+1 ; td_cs_cm2(687) = 3.789600d-25
  td_E_eV(688) = 1.431000d+1 ; td_cs_cm2(688) = 3.773300d-25
  td_E_eV(689) = 1.433000d+1 ; td_cs_cm2(689) = 3.757100d-25
  td_E_eV(690) = 1.435000d+1 ; td_cs_cm2(690) = 3.741000d-25
  td_E_eV(691) = 1.437000d+1 ; td_cs_cm2(691) = 3.725000d-25
  td_E_eV(692) = 1.439000d+1 ; td_cs_cm2(692) = 3.709100d-25
  td_E_eV(693) = 1.441000d+1 ; td_cs_cm2(693) = 3.693300d-25
  td_E_eV(694) = 1.443000d+1 ; td_cs_cm2(694) = 3.677600d-25
  td_E_eV(695) = 1.445000d+1 ; td_cs_cm2(695) = 3.662000d-25
  td_E_eV(696) = 1.447000d+1 ; td_cs_cm2(696) = 3.646500d-25
  td_E_eV(697) = 1.449000d+1 ; td_cs_cm2(697) = 3.631000d-25
  td_E_eV(698) = 1.451000d+1 ; td_cs_cm2(698) = 3.615700d-25
  td_E_eV(699) = 1.453000d+1 ; td_cs_cm2(699) = 3.600500d-25
  td_E_eV(700) = 1.455000d+1 ; td_cs_cm2(700) = 3.585400d-25
  td_E_eV(701) = 1.457000d+1 ; td_cs_cm2(701) = 3.570300d-25
  td_E_eV(702) = 1.459000d+1 ; td_cs_cm2(702) = 3.555400d-25
  td_E_eV(703) = 1.461000d+1 ; td_cs_cm2(703) = 3.540500d-25
  td_E_eV(704) = 1.463000d+1 ; td_cs_cm2(704) = 3.525800d-25
  td_E_eV(705) = 1.465000d+1 ; td_cs_cm2(705) = 3.511100d-25
  td_E_eV(706) = 1.467000d+1 ; td_cs_cm2(706) = 3.496500d-25
  td_E_eV(707) = 1.469000d+1 ; td_cs_cm2(707) = 3.482000d-25
  td_E_eV(708) = 1.471000d+1 ; td_cs_cm2(708) = 3.467600d-25
  td_E_eV(709) = 1.473000d+1 ; td_cs_cm2(709) = 3.453300d-25
  td_E_eV(710) = 1.475000d+1 ; td_cs_cm2(710) = 3.439100d-25
  td_E_eV(711) = 1.477000d+1 ; td_cs_cm2(711) = 3.424900d-25
  td_E_eV(712) = 1.479000d+1 ; td_cs_cm2(712) = 3.410900d-25
  td_E_eV(713) = 1.481000d+1 ; td_cs_cm2(713) = 3.396900d-25
  td_E_eV(714) = 1.483000d+1 ; td_cs_cm2(714) = 3.383000d-25
  td_E_eV(715) = 1.485000d+1 ; td_cs_cm2(715) = 3.369200d-25
  td_E_eV(716) = 1.487000d+1 ; td_cs_cm2(716) = 3.355500d-25
  td_E_eV(717) = 1.489000d+1 ; td_cs_cm2(717) = 3.341800d-25
  td_E_eV(718) = 1.491000d+1 ; td_cs_cm2(718) = 3.328300d-25
  td_E_eV(719) = 1.493000d+1 ; td_cs_cm2(719) = 3.314800d-25
  td_E_eV(720) = 1.495000d+1 ; td_cs_cm2(720) = 3.301400d-25
  td_E_eV(721) = 1.497000d+1 ; td_cs_cm2(721) = 3.288100d-25
  td_E_eV(722) = 1.499000d+1 ; td_cs_cm2(722) = 3.274800d-25

  td_cs_cm2 = 1.0d4 * td_cs_cm2 ! the original data are in m^2, convert to cm^2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_vib02
!
!-------------------------------------
!
real(8) function CSV_N2_vib02_m3s(E_eV)

  use N2_vib02
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_vib02_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_vib02_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_vib02_m3s










!---------------------------------------------------------------------------------------------------- #6
!
module N2_vib03
  integer, parameter :: N_td=690  !8   ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_vib03
!
!-------------------------------------
!
subroutine Prepare_N2_vib03

  use N2_vib03
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

!  td_E_eV(1) = 1.3_8; td_cs_cm2(1) = 1.11d-18
!  td_E_eV(2) = 1.6_8; td_cs_cm2(2) = 8.84d-17
!  td_E_eV(3) = 2.0_8; td_cs_cm2(3) = 1.01d-15
!  td_E_eV(4) = 2.4_8; td_cs_cm2(4) = 1.35d-15
!  td_E_eV(5) = 3.0_8; td_cs_cm2(5) = 5.86d-16
!  td_E_eV(6) = 4.0_8; td_cs_cm2(6) = 1.37d-16
!  td_E_eV(7) = 5.0_8; td_cs_cm2(7) = 2.34d-17
!  td_E_eV(8) = 6.0_8; td_cs_cm2(8) = 3.66d-18

  td_E_eV(  1) = 1.210000d+0 ; td_cs_cm2(  1) = 0.000000d+0
  td_E_eV(  2) = 1.230000d+0 ; td_cs_cm2(  2) = 1.383800d-26
  td_E_eV(  3) = 1.250000d+0 ; td_cs_cm2(  3) = 2.044900d-26
  td_E_eV(  4) = 1.270000d+0 ; td_cs_cm2(  4) = 2.982600d-26
  td_E_eV(  5) = 1.290000d+0 ; td_cs_cm2(  5) = 4.298500d-26
  td_E_eV(  6) = 1.310000d+0 ; td_cs_cm2(  6) = 6.130200d-26
  td_E_eV(  7) = 1.330000d+0 ; td_cs_cm2(  7) = 8.663700d-26
  td_E_eV(  8) = 1.350000d+0 ; td_cs_cm2(  8) = 1.215200d-25
  td_E_eV(  9) = 1.370000d+0 ; td_cs_cm2(  9) = 1.694000d-25
  td_E_eV( 10) = 1.390000d+0 ; td_cs_cm2( 10) = 2.349600d-25
  td_E_eV( 11) = 1.410000d+0 ; td_cs_cm2( 11) = 3.246400d-25
  td_E_eV( 12) = 1.430000d+0 ; td_cs_cm2( 12) = 4.472900d-25
  td_E_eV( 13) = 1.450000d+0 ; td_cs_cm2( 13) = 6.151100d-25
  td_E_eV( 14) = 1.470000d+0 ; td_cs_cm2( 14) = 8.450600d-25
  td_E_eV( 15) = 1.490000d+0 ; td_cs_cm2( 15) = 1.160800d-24
  td_E_eV( 16) = 1.510000d+0 ; td_cs_cm2( 16) = 1.595500d-24
  td_E_eV( 17) = 1.530000d+0 ; td_cs_cm2( 17) = 2.196200d-24
  td_E_eV( 18) = 1.550000d+0 ; td_cs_cm2( 18) = 3.029600d-24
  td_E_eV( 19) = 1.570000d+0 ; td_cs_cm2( 19) = 4.191800d-24
  td_E_eV( 20) = 1.590000d+0 ; td_cs_cm2( 20) = 5.821900d-24
  td_E_eV( 21) = 1.610000d+0 ; td_cs_cm2( 21) = 8.123100d-24
  td_E_eV( 22) = 1.630000d+0 ; td_cs_cm2( 22) = 1.139600d-23
  td_E_eV( 23) = 1.650000d+0 ; td_cs_cm2( 23) = 1.609000d-23
  td_E_eV( 24) = 1.670000d+0 ; td_cs_cm2( 24) = 2.288400d-23
  td_E_eV( 25) = 1.690000d+0 ; td_cs_cm2( 25) = 3.282200d-23
  td_E_eV( 26) = 1.710000d+0 ; td_cs_cm2( 26) = 4.752400d-23
  td_E_eV( 27) = 1.730000d+0 ; td_cs_cm2( 27) = 6.955300d-23
  td_E_eV( 28) = 1.750000d+0 ; td_cs_cm2( 28) = 1.030100d-22
  td_E_eV( 29) = 1.770000d+0 ; td_cs_cm2( 29) = 1.545800d-22
  td_E_eV( 30) = 1.790000d+0 ; td_cs_cm2( 30) = 2.352800d-22
  td_E_eV( 31) = 1.810000d+0 ; td_cs_cm2( 31) = 3.634700d-22
  td_E_eV( 32) = 1.830000d+0 ; td_cs_cm2( 32) = 5.697300d-22
  td_E_eV( 33) = 1.850000d+0 ; td_cs_cm2( 33) = 9.042600d-22
  td_E_eV( 34) = 1.870000d+0 ; td_cs_cm2( 34) = 1.445300d-21
  td_E_eV( 35) = 1.890000d+0 ; td_cs_cm2( 35) = 2.300900d-21
  td_E_eV( 36) = 1.910000d+0 ; td_cs_cm2( 36) = 3.580200d-21
  td_E_eV( 37) = 1.930000d+0 ; td_cs_cm2( 37) = 5.307500d-21
  td_E_eV( 38) = 1.950000d+0 ; td_cs_cm2( 38) = 7.316600d-21
  td_E_eV( 39) = 1.970000d+0 ; td_cs_cm2( 39) = 9.295000d-21
  td_E_eV( 40) = 1.990000d+0 ; td_cs_cm2( 40) = 1.102400d-20
  td_E_eV( 41) = 2.010000d+0 ; td_cs_cm2( 41) = 1.253400d-20
  td_E_eV( 42) = 2.030000d+0 ; td_cs_cm2( 42) = 1.402000d-20
  td_E_eV( 43) = 2.050000d+0 ; td_cs_cm2( 43) = 1.570800d-20
  td_E_eV( 44) = 2.070000d+0 ; td_cs_cm2( 44) = 1.778400d-20
  td_E_eV( 45) = 2.090000d+0 ; td_cs_cm2( 45) = 2.036500d-20
  td_E_eV( 46) = 2.110000d+0 ; td_cs_cm2( 46) = 2.342600d-20
  td_E_eV( 47) = 2.130000d+0 ; td_cs_cm2( 47) = 2.666200d-20
  td_E_eV( 48) = 2.150000d+0 ; td_cs_cm2( 48) = 2.931800d-20
  td_E_eV( 49) = 2.170000d+0 ; td_cs_cm2( 49) = 3.029000d-20
  td_E_eV( 50) = 2.190000d+0 ; td_cs_cm2( 50) = 2.880100d-20
  td_E_eV( 51) = 2.210000d+0 ; td_cs_cm2( 51) = 2.513400d-20
  td_E_eV( 52) = 2.230000d+0 ; td_cs_cm2( 52) = 2.041300d-20
  td_E_eV( 53) = 2.250000d+0 ; td_cs_cm2( 53) = 1.571100d-20
  td_E_eV( 54) = 2.270000d+0 ; td_cs_cm2( 54) = 1.157500d-20
  td_E_eV( 55) = 2.290000d+0 ; td_cs_cm2( 55) = 8.137700d-21
  td_E_eV( 56) = 2.310000d+0 ; td_cs_cm2( 56) = 5.361300d-21
  td_E_eV( 57) = 2.330000d+0 ; td_cs_cm2( 57) = 3.217600d-21
  td_E_eV( 58) = 2.350000d+0 ; td_cs_cm2( 58) = 1.786200d-21
  td_E_eV( 59) = 2.370000d+0 ; td_cs_cm2( 59) = 1.277900d-21
  td_E_eV( 60) = 2.390000d+0 ; td_cs_cm2( 60) = 1.920200d-21
  td_E_eV( 61) = 2.410000d+0 ; td_cs_cm2( 61) = 3.692800d-21
  td_E_eV( 62) = 2.430000d+0 ; td_cs_cm2( 62) = 6.168400d-21
  td_E_eV( 63) = 2.450000d+0 ; td_cs_cm2( 63) = 8.750200d-21
  td_E_eV( 64) = 2.470000d+0 ; td_cs_cm2( 64) = 1.104200d-20
  td_E_eV( 65) = 2.490000d+0 ; td_cs_cm2( 65) = 1.293700d-20
  td_E_eV( 66) = 2.510000d+0 ; td_cs_cm2( 66) = 1.447500d-20
  td_E_eV( 67) = 2.530000d+0 ; td_cs_cm2( 67) = 1.569300d-20
  td_E_eV( 68) = 2.550000d+0 ; td_cs_cm2( 68) = 1.652600d-20
  td_E_eV( 69) = 2.570000d+0 ; td_cs_cm2( 69) = 1.676500d-20
  td_E_eV( 70) = 2.590000d+0 ; td_cs_cm2( 70) = 1.607000d-20
  td_E_eV( 71) = 2.610000d+0 ; td_cs_cm2( 71) = 1.415600d-20
  td_E_eV( 72) = 2.630000d+0 ; td_cs_cm2( 72) = 1.112900d-20
  td_E_eV( 73) = 2.650000d+0 ; td_cs_cm2( 73) = 7.640600d-21
  td_E_eV( 74) = 2.670000d+0 ; td_cs_cm2( 74) = 4.535800d-21
  td_E_eV( 75) = 2.690000d+0 ; td_cs_cm2( 75) = 2.345700d-21
  td_E_eV( 76) = 2.710000d+0 ; td_cs_cm2( 76) = 1.179900d-21
  td_E_eV( 77) = 2.730000d+0 ; td_cs_cm2( 77) = 9.271300d-22
  td_E_eV( 78) = 2.750000d+0 ; td_cs_cm2( 78) = 1.440400d-21
  td_E_eV( 79) = 2.770000d+0 ; td_cs_cm2( 79) = 2.598600d-21
  td_E_eV( 80) = 2.790000d+0 ; td_cs_cm2( 80) = 4.268200d-21
  td_E_eV( 81) = 2.810000d+0 ; td_cs_cm2( 81) = 6.205500d-21
  td_E_eV( 82) = 2.830000d+0 ; td_cs_cm2( 82) = 7.978000d-21
  td_E_eV( 83) = 2.850000d+0 ; td_cs_cm2( 83) = 9.056400d-21
  td_E_eV( 84) = 2.870000d+0 ; td_cs_cm2( 84) = 9.119300d-21
  td_E_eV( 85) = 2.890000d+0 ; td_cs_cm2( 85) = 8.268300d-21
  td_E_eV( 86) = 2.910000d+0 ; td_cs_cm2( 86) = 6.883700d-21
  td_E_eV( 87) = 2.930000d+0 ; td_cs_cm2( 87) = 5.340900d-21
  td_E_eV( 88) = 2.950000d+0 ; td_cs_cm2( 88) = 3.878300d-21
  td_E_eV( 89) = 2.970000d+0 ; td_cs_cm2( 89) = 2.618100d-21
  td_E_eV( 90) = 2.990000d+0 ; td_cs_cm2( 90) = 1.631200d-21
  td_E_eV( 91) = 3.010000d+0 ; td_cs_cm2( 91) = 9.846000d-22
  td_E_eV( 92) = 3.030000d+0 ; td_cs_cm2( 92) = 7.468700d-22
  td_E_eV( 93) = 3.050000d+0 ; td_cs_cm2( 93) = 9.406000d-22
  td_E_eV( 94) = 3.070000d+0 ; td_cs_cm2( 94) = 1.471400d-21
  td_E_eV( 95) = 3.090000d+0 ; td_cs_cm2( 95) = 2.129200d-21
  td_E_eV( 96) = 3.110000d+0 ; td_cs_cm2( 96) = 2.698000d-21
  td_E_eV( 97) = 3.130000d+0 ; td_cs_cm2( 97) = 3.058700d-21
  td_E_eV( 98) = 3.150000d+0 ; td_cs_cm2( 98) = 3.191600d-21
  td_E_eV( 99) = 3.170000d+0 ; td_cs_cm2( 99) = 3.123400d-21
  td_E_eV(100) = 3.190000d+0 ; td_cs_cm2(100) = 2.885500d-21
  td_E_eV(101) = 3.210000d+0 ; td_cs_cm2(101) = 2.502800d-21
  td_E_eV(102) = 3.230000d+0 ; td_cs_cm2(102) = 2.005600d-21
  td_E_eV(103) = 3.250000d+0 ; td_cs_cm2(103) = 1.453500d-21
  td_E_eV(104) = 3.270000d+0 ; td_cs_cm2(104) = 9.474000d-22
  td_E_eV(105) = 3.290000d+0 ; td_cs_cm2(105) = 5.981600d-22
  td_E_eV(106) = 3.310000d+0 ; td_cs_cm2(106) = 4.626000d-22
  td_E_eV(107) = 3.330000d+0 ; td_cs_cm2(107) = 5.167800d-22
  td_E_eV(108) = 3.350000d+0 ; td_cs_cm2(108) = 6.908400d-22
  td_E_eV(109) = 3.370000d+0 ; td_cs_cm2(109) = 9.141500d-22
  td_E_eV(110) = 3.390000d+0 ; td_cs_cm2(110) = 1.132600d-21
  td_E_eV(111) = 3.410000d+0 ; td_cs_cm2(111) = 1.303900d-21
  td_E_eV(112) = 3.430000d+0 ; td_cs_cm2(112) = 1.389600d-21
  td_E_eV(113) = 3.450000d+0 ; td_cs_cm2(113) = 1.355500d-21
  td_E_eV(114) = 3.470000d+0 ; td_cs_cm2(114) = 1.190300d-21
  td_E_eV(115) = 3.490000d+0 ; td_cs_cm2(115) = 9.284900d-22
  td_E_eV(116) = 3.510000d+0 ; td_cs_cm2(116) = 6.459100d-22
  td_E_eV(117) = 3.530000d+0 ; td_cs_cm2(117) = 4.184700d-22
  td_E_eV(118) = 3.550000d+0 ; td_cs_cm2(118) = 2.847300d-22
  td_E_eV(119) = 3.570000d+0 ; td_cs_cm2(119) = 2.451400d-22
  td_E_eV(120) = 3.590000d+0 ; td_cs_cm2(120) = 2.805100d-22
  td_E_eV(121) = 3.610000d+0 ; td_cs_cm2(121) = 3.657700d-22
  td_E_eV(122) = 3.630000d+0 ; td_cs_cm2(122) = 4.731000d-22
  td_E_eV(123) = 3.650000d+0 ; td_cs_cm2(123) = 5.702900d-22
  td_E_eV(124) = 3.670000d+0 ; td_cs_cm2(124) = 6.229900d-22
  td_E_eV(125) = 3.690000d+0 ; td_cs_cm2(125) = 6.075600d-22
  td_E_eV(126) = 3.710000d+0 ; td_cs_cm2(126) = 5.272900d-22
  td_E_eV(127) = 3.730000d+0 ; td_cs_cm2(127) = 4.119500d-22
  td_E_eV(128) = 3.750000d+0 ; td_cs_cm2(128) = 2.977600d-22
  td_E_eV(129) = 3.770000d+0 ; td_cs_cm2(129) = 2.095400d-22
  td_E_eV(130) = 3.790000d+0 ; td_cs_cm2(130) = 1.575200d-22
  td_E_eV(131) = 3.810000d+0 ; td_cs_cm2(131) = 1.421600d-22
  td_E_eV(132) = 3.830000d+0 ; td_cs_cm2(132) = 1.582700d-22
  td_E_eV(133) = 3.850000d+0 ; td_cs_cm2(133) = 1.956600d-22
  td_E_eV(134) = 3.870000d+0 ; td_cs_cm2(134) = 2.384600d-22
  td_E_eV(135) = 3.890000d+0 ; td_cs_cm2(135) = 2.678500d-22
  td_E_eV(136) = 3.910000d+0 ; td_cs_cm2(136) = 2.704500d-22
  td_E_eV(137) = 3.930000d+0 ; td_cs_cm2(137) = 2.459800d-22
  td_E_eV(138) = 3.950000d+0 ; td_cs_cm2(138) = 2.051100d-22
  td_E_eV(139) = 3.970000d+0 ; td_cs_cm2(139) = 1.610100d-22
  td_E_eV(140) = 3.990000d+0 ; td_cs_cm2(140) = 1.235000d-22
  td_E_eV(141) = 4.010000d+0 ; td_cs_cm2(141) = 9.810600d-23
  td_E_eV(142) = 4.030000d+0 ; td_cs_cm2(142) = 8.704400d-23
  td_E_eV(143) = 4.050000d+0 ; td_cs_cm2(143) = 8.959500d-23
  td_E_eV(144) = 4.070000d+0 ; td_cs_cm2(144) = 1.015700d-22
  td_E_eV(145) = 4.090000d+0 ; td_cs_cm2(145) = 1.155400d-22
  td_E_eV(146) = 4.110000d+0 ; td_cs_cm2(146) = 1.236800d-22
  td_E_eV(147) = 4.130000d+0 ; td_cs_cm2(147) = 1.220000d-22
  td_E_eV(148) = 4.150000d+0 ; td_cs_cm2(148) = 1.115600d-22
  td_E_eV(149) = 4.170000d+0 ; td_cs_cm2(149) = 9.624000d-23
  td_E_eV(150) = 4.190000d+0 ; td_cs_cm2(150) = 8.013300d-23
  td_E_eV(151) = 4.210000d+0 ; td_cs_cm2(151) = 6.638200d-23
  td_E_eV(152) = 4.230000d+0 ; td_cs_cm2(152) = 5.706700d-23
  td_E_eV(153) = 4.250000d+0 ; td_cs_cm2(153) = 5.314100d-23
  td_E_eV(154) = 4.270000d+0 ; td_cs_cm2(154) = 5.405600d-23
  td_E_eV(155) = 4.290000d+0 ; td_cs_cm2(155) = 5.756900d-23
  td_E_eV(156) = 4.310000d+0 ; td_cs_cm2(156) = 6.061700d-23
  td_E_eV(157) = 4.330000d+0 ; td_cs_cm2(157) = 6.098100d-23
  td_E_eV(158) = 4.350000d+0 ; td_cs_cm2(158) = 5.818100d-23
  td_E_eV(159) = 4.370000d+0 ; td_cs_cm2(159) = 5.304900d-23
  td_E_eV(160) = 4.390000d+0 ; td_cs_cm2(160) = 4.687700d-23
  td_E_eV(161) = 4.410000d+0 ; td_cs_cm2(161) = 4.091100d-23
  td_E_eV(162) = 4.430000d+0 ; td_cs_cm2(162) = 3.617200d-23
  td_E_eV(163) = 4.450000d+0 ; td_cs_cm2(163) = 3.331500d-23
  td_E_eV(164) = 4.470000d+0 ; td_cs_cm2(164) = 3.238800d-23
  td_E_eV(165) = 4.490000d+0 ; td_cs_cm2(165) = 3.270200d-23
  td_E_eV(166) = 4.510000d+0 ; td_cs_cm2(166) = 3.316200d-23
  td_E_eV(167) = 4.530000d+0 ; td_cs_cm2(167) = 3.290700d-23
  td_E_eV(168) = 4.550000d+0 ; td_cs_cm2(168) = 3.164500d-23
  td_E_eV(169) = 4.570000d+0 ; td_cs_cm2(169) = 2.955100d-23
  td_E_eV(170) = 4.590000d+0 ; td_cs_cm2(170) = 2.701300d-23
  td_E_eV(171) = 4.610000d+0 ; td_cs_cm2(171) = 2.446900d-23
  td_E_eV(172) = 4.630000d+0 ; td_cs_cm2(172) = 2.232500d-23
  td_E_eV(173) = 4.650000d+0 ; td_cs_cm2(173) = 2.086100d-23
  td_E_eV(174) = 4.670000d+0 ; td_cs_cm2(174) = 2.010400d-23
  td_E_eV(175) = 4.690000d+0 ; td_cs_cm2(175) = 1.979500d-23
  td_E_eV(176) = 4.710000d+0 ; td_cs_cm2(176) = 1.954800d-23
  td_E_eV(177) = 4.730000d+0 ; td_cs_cm2(177) = 1.907700d-23
  td_E_eV(178) = 4.750000d+0 ; td_cs_cm2(178) = 1.829000d-23
  td_E_eV(179) = 4.770000d+0 ; td_cs_cm2(179) = 1.723700d-23
  td_E_eV(180) = 4.790000d+0 ; td_cs_cm2(180) = 1.605000d-23
  td_E_eV(181) = 4.810000d+0 ; td_cs_cm2(181) = 1.488800d-23
  td_E_eV(182) = 4.830000d+0 ; td_cs_cm2(182) = 1.390400d-23
  td_E_eV(183) = 4.850000d+0 ; td_cs_cm2(183) = 1.319100d-23
  td_E_eV(184) = 4.870000d+0 ; td_cs_cm2(184) = 1.273000d-23
  td_E_eV(185) = 4.890000d+0 ; td_cs_cm2(185) = 1.240300d-23
  td_E_eV(186) = 4.910000d+0 ; td_cs_cm2(186) = 1.207600d-23
  td_E_eV(187) = 4.930000d+0 ; td_cs_cm2(187) = 1.166600d-23
  td_E_eV(188) = 4.950000d+0 ; td_cs_cm2(188) = 1.115400d-23
  td_E_eV(189) = 4.970000d+0 ; td_cs_cm2(189) = 1.056600d-23
  td_E_eV(190) = 4.990000d+0 ; td_cs_cm2(190) = 9.953700d-24
  td_E_eV(191) = 5.010000d+0 ; td_cs_cm2(191) = 9.380700d-24
  td_E_eV(192) = 5.030000d+0 ; td_cs_cm2(192) = 8.901400d-24
  td_E_eV(193) = 5.050000d+0 ; td_cs_cm2(193) = 8.534000d-24
  td_E_eV(194) = 5.070000d+0 ; td_cs_cm2(194) = 8.250400d-24
  td_E_eV(195) = 5.090000d+0 ; td_cs_cm2(195) = 7.999000d-24
  td_E_eV(196) = 5.110000d+0 ; td_cs_cm2(196) = 7.736600d-24
  td_E_eV(197) = 5.130000d+0 ; td_cs_cm2(197) = 7.442900d-24
  td_E_eV(198) = 5.150000d+0 ; td_cs_cm2(198) = 7.117800d-24
  td_E_eV(199) = 5.170000d+0 ; td_cs_cm2(199) = 6.776100d-24
  td_E_eV(200) = 5.190000d+0 ; td_cs_cm2(200) = 6.442000d-24
  td_E_eV(201) = 5.210000d+0 ; td_cs_cm2(201) = 6.141200d-24
  td_E_eV(202) = 5.230000d+0 ; td_cs_cm2(202) = 5.888200d-24
  td_E_eV(203) = 5.250000d+0 ; td_cs_cm2(203) = 5.678400d-24
  td_E_eV(204) = 5.270000d+0 ; td_cs_cm2(204) = 5.493400d-24
  td_E_eV(205) = 5.290000d+0 ; td_cs_cm2(205) = 5.314100d-24
  td_E_eV(206) = 5.310000d+0 ; td_cs_cm2(206) = 5.128100d-24
  td_E_eV(207) = 5.330000d+0 ; td_cs_cm2(207) = 4.931700d-24
  td_E_eV(208) = 5.350000d+0 ; td_cs_cm2(208) = 4.728000d-24
  td_E_eV(209) = 5.370000d+0 ; td_cs_cm2(209) = 4.525700d-24
  td_E_eV(210) = 5.390000d+0 ; td_cs_cm2(210) = 4.336100d-24
  td_E_eV(211) = 5.410000d+0 ; td_cs_cm2(211) = 4.167500d-24
  td_E_eV(212) = 5.430000d+0 ; td_cs_cm2(212) = 4.020400d-24
  td_E_eV(213) = 5.450000d+0 ; td_cs_cm2(213) = 3.888500d-24
  td_E_eV(214) = 5.470000d+0 ; td_cs_cm2(214) = 3.763700d-24
  td_E_eV(215) = 5.490000d+0 ; td_cs_cm2(215) = 3.639600d-24
  td_E_eV(216) = 5.510000d+0 ; td_cs_cm2(216) = 3.513200d-24
  td_E_eV(217) = 5.530000d+0 ; td_cs_cm2(217) = 3.384600d-24
  td_E_eV(218) = 5.550000d+0 ; td_cs_cm2(218) = 3.256600d-24
  td_E_eV(219) = 5.570000d+0 ; td_cs_cm2(219) = 3.134400d-24
  td_E_eV(220) = 5.590000d+0 ; td_cs_cm2(220) = 3.021900d-24
  td_E_eV(221) = 5.610000d+0 ; td_cs_cm2(221) = 2.920200d-24
  td_E_eV(222) = 5.630000d+0 ; td_cs_cm2(222) = 2.827200d-24
  td_E_eV(223) = 5.650000d+0 ; td_cs_cm2(223) = 2.739500d-24
  td_E_eV(224) = 5.670000d+0 ; td_cs_cm2(224) = 2.654000d-24
  td_E_eV(225) = 5.690000d+0 ; td_cs_cm2(225) = 2.569000d-24
  td_E_eV(226) = 5.710000d+0 ; td_cs_cm2(226) = 2.483700d-24
  td_E_eV(227) = 5.730000d+0 ; td_cs_cm2(227) = 2.399400d-24
  td_E_eV(228) = 5.750000d+0 ; td_cs_cm2(228) = 2.318100d-24
  td_E_eV(229) = 5.770000d+0 ; td_cs_cm2(229) = 2.241900d-24
  td_E_eV(230) = 5.790000d+0 ; td_cs_cm2(230) = 2.171400d-24
  td_E_eV(231) = 5.810000d+0 ; td_cs_cm2(231) = 2.105700d-24
  td_E_eV(232) = 5.830000d+0 ; td_cs_cm2(232) = 2.043500d-24
  td_E_eV(233) = 5.850000d+0 ; td_cs_cm2(233) = 1.983400d-24
  td_E_eV(234) = 5.870000d+0 ; td_cs_cm2(234) = 1.924200d-24
  td_E_eV(235) = 5.890000d+0 ; td_cs_cm2(235) = 1.865500d-24
  td_E_eV(236) = 5.910000d+0 ; td_cs_cm2(236) = 1.807800d-24
  td_E_eV(237) = 5.930000d+0 ; td_cs_cm2(237) = 1.752200d-24
  td_E_eV(238) = 5.950000d+0 ; td_cs_cm2(238) = 1.699400d-24
  td_E_eV(239) = 5.970000d+0 ; td_cs_cm2(239) = 1.649800d-24
  td_E_eV(240) = 5.990000d+0 ; td_cs_cm2(240) = 1.603000d-24
  td_E_eV(241) = 6.010000d+0 ; td_cs_cm2(241) = 1.558400d-24
  td_E_eV(242) = 6.030000d+0 ; td_cs_cm2(242) = 1.515200d-24
  td_E_eV(243) = 6.050000d+0 ; td_cs_cm2(243) = 1.473000d-24
  td_E_eV(244) = 6.070000d+0 ; td_cs_cm2(244) = 1.431500d-24
  td_E_eV(245) = 6.090000d+0 ; td_cs_cm2(245) = 1.390800d-24
  td_E_eV(246) = 6.110000d+0 ; td_cs_cm2(246) = 1.351600d-24
  td_E_eV(247) = 6.130000d+0 ; td_cs_cm2(247) = 1.314200d-24
  td_E_eV(248) = 6.150000d+0 ; td_cs_cm2(248) = 1.278600d-24
  td_E_eV(249) = 6.170000d+0 ; td_cs_cm2(249) = 1.244700d-24
  td_E_eV(250) = 6.190000d+0 ; td_cs_cm2(250) = 1.212200d-24
  td_E_eV(251) = 6.210000d+0 ; td_cs_cm2(251) = 1.180700d-24
  td_E_eV(252) = 6.230000d+0 ; td_cs_cm2(252) = 1.150000d-24
  td_E_eV(253) = 6.250000d+0 ; td_cs_cm2(253) = 1.119800d-24
  td_E_eV(254) = 6.270000d+0 ; td_cs_cm2(254) = 1.090500d-24
  td_E_eV(255) = 6.290000d+0 ; td_cs_cm2(255) = 1.062200d-24
  td_E_eV(256) = 6.310000d+0 ; td_cs_cm2(256) = 1.035000d-24
  td_E_eV(257) = 6.330000d+0 ; td_cs_cm2(257) = 1.008900d-24
  td_E_eV(258) = 6.350000d+0 ; td_cs_cm2(258) = 9.840000d-25
  td_E_eV(259) = 6.370000d+0 ; td_cs_cm2(259) = 9.599400d-25
  td_E_eV(260) = 6.390000d+0 ; td_cs_cm2(260) = 9.365600d-25
  td_E_eV(261) = 6.410000d+0 ; td_cs_cm2(261) = 9.137200d-25
  td_E_eV(262) = 6.430000d+0 ; td_cs_cm2(262) = 8.914300d-25
  td_E_eV(263) = 6.450000d+0 ; td_cs_cm2(263) = 8.697700d-25
  td_E_eV(264) = 6.470000d+0 ; td_cs_cm2(264) = 8.488300d-25
  td_E_eV(265) = 6.490000d+0 ; td_cs_cm2(265) = 8.286700d-25
  td_E_eV(266) = 6.510000d+0 ; td_cs_cm2(266) = 8.092700d-25
  td_E_eV(267) = 6.530000d+0 ; td_cs_cm2(267) = 7.905600d-25
  td_E_eV(268) = 6.550000d+0 ; td_cs_cm2(268) = 7.724400d-25
  td_E_eV(269) = 6.570000d+0 ; td_cs_cm2(269) = 7.547800d-25
  td_E_eV(270) = 6.590000d+0 ; td_cs_cm2(270) = 7.375500d-25
  td_E_eV(271) = 6.610000d+0 ; td_cs_cm2(271) = 7.207600d-25
  td_E_eV(272) = 6.630000d+0 ; td_cs_cm2(272) = 7.044400d-25
  td_E_eV(273) = 6.650000d+0 ; td_cs_cm2(273) = 6.886600d-25
  td_E_eV(274) = 6.670000d+0 ; td_cs_cm2(274) = 6.734100d-25
  td_E_eV(275) = 6.690000d+0 ; td_cs_cm2(275) = 6.586900d-25
  td_E_eV(276) = 6.710000d+0 ; td_cs_cm2(276) = 6.444200d-25
  td_E_eV(277) = 6.730000d+0 ; td_cs_cm2(277) = 6.305400d-25
  td_E_eV(278) = 6.750000d+0 ; td_cs_cm2(278) = 6.170000d-25
  td_E_eV(279) = 6.770000d+0 ; td_cs_cm2(279) = 6.037900d-25
  td_E_eV(280) = 6.790000d+0 ; td_cs_cm2(280) = 5.909200d-25
  td_E_eV(281) = 6.810000d+0 ; td_cs_cm2(281) = 5.784200d-25
  td_E_eV(282) = 6.830000d+0 ; td_cs_cm2(282) = 5.663100d-25
  td_E_eV(283) = 6.850000d+0 ; td_cs_cm2(283) = 5.545800d-25
  td_E_eV(284) = 6.870000d+0 ; td_cs_cm2(284) = 5.432100d-25
  td_E_eV(285) = 6.890000d+0 ; td_cs_cm2(285) = 5.321500d-25
  td_E_eV(286) = 6.910000d+0 ; td_cs_cm2(286) = 5.213600d-25
  td_E_eV(287) = 6.930000d+0 ; td_cs_cm2(287) = 5.108200d-25
  td_E_eV(288) = 6.950000d+0 ; td_cs_cm2(288) = 5.005400d-25
  td_E_eV(289) = 6.970000d+0 ; td_cs_cm2(289) = 4.905300d-25
  td_E_eV(290) = 6.990000d+0 ; td_cs_cm2(290) = 4.808100d-25
  td_E_eV(291) = 7.010000d+0 ; td_cs_cm2(291) = 4.713800d-25
  td_E_eV(292) = 7.030000d+0 ; td_cs_cm2(292) = 4.622100d-25
  td_E_eV(293) = 7.050000d+0 ; td_cs_cm2(293) = 4.532900d-25
  td_E_eV(294) = 7.070000d+0 ; td_cs_cm2(294) = 4.445800d-25
  td_E_eV(295) = 7.090000d+0 ; td_cs_cm2(295) = 4.360600d-25
  td_E_eV(296) = 7.110000d+0 ; td_cs_cm2(296) = 4.277500d-25
  td_E_eV(297) = 7.130000d+0 ; td_cs_cm2(297) = 4.196500d-25
  td_E_eV(298) = 7.150000d+0 ; td_cs_cm2(298) = 4.117700d-25
  td_E_eV(299) = 7.170000d+0 ; td_cs_cm2(299) = 4.041000d-25
  td_E_eV(300) = 7.190000d+0 ; td_cs_cm2(300) = 3.966400d-25
  td_E_eV(301) = 7.210000d+0 ; td_cs_cm2(301) = 3.893700d-25
  td_E_eV(302) = 7.230000d+0 ; td_cs_cm2(302) = 3.822600d-25
  td_E_eV(303) = 7.250000d+0 ; td_cs_cm2(303) = 3.753100d-25
  td_E_eV(304) = 7.270000d+0 ; td_cs_cm2(304) = 3.685200d-25
  td_E_eV(305) = 7.290000d+0 ; td_cs_cm2(305) = 3.618900d-25
  td_E_eV(306) = 7.310000d+0 ; td_cs_cm2(306) = 3.554400d-25
  td_E_eV(307) = 7.330000d+0 ; td_cs_cm2(307) = 3.491500d-25
  td_E_eV(308) = 7.350000d+0 ; td_cs_cm2(308) = 3.430100d-25
  td_E_eV(309) = 7.370000d+0 ; td_cs_cm2(309) = 3.370200d-25
  td_E_eV(310) = 7.390000d+0 ; td_cs_cm2(310) = 3.311700d-25
  td_E_eV(311) = 7.410000d+0 ; td_cs_cm2(311) = 3.254400d-25
  td_E_eV(312) = 7.430000d+0 ; td_cs_cm2(312) = 3.198400d-25
  td_E_eV(313) = 7.450000d+0 ; td_cs_cm2(313) = 3.143700d-25
  td_E_eV(314) = 7.470000d+0 ; td_cs_cm2(314) = 3.090400d-25
  td_E_eV(315) = 7.490000d+0 ; td_cs_cm2(315) = 3.038300d-25
  td_E_eV(316) = 7.510000d+0 ; td_cs_cm2(316) = 2.987400d-25
  td_E_eV(317) = 7.530000d+0 ; td_cs_cm2(317) = 2.937600d-25
  td_E_eV(318) = 7.550000d+0 ; td_cs_cm2(318) = 2.888900d-25
  td_E_eV(319) = 7.570000d+0 ; td_cs_cm2(319) = 2.841300d-25
  td_E_eV(320) = 7.590000d+0 ; td_cs_cm2(320) = 2.794700d-25
  td_E_eV(321) = 7.610000d+0 ; td_cs_cm2(321) = 2.749200d-25
  td_E_eV(322) = 7.630000d+0 ; td_cs_cm2(322) = 2.704700d-25
  td_E_eV(323) = 7.650000d+0 ; td_cs_cm2(323) = 2.661200d-25
  td_E_eV(324) = 7.670000d+0 ; td_cs_cm2(324) = 2.618600d-25
  td_E_eV(325) = 7.690000d+0 ; td_cs_cm2(325) = 2.577000d-25
  td_E_eV(326) = 7.710000d+0 ; td_cs_cm2(326) = 2.536200d-25
  td_E_eV(327) = 7.730000d+0 ; td_cs_cm2(327) = 2.496300d-25
  td_E_eV(328) = 7.750000d+0 ; td_cs_cm2(328) = 2.457200d-25
  td_E_eV(329) = 7.770000d+0 ; td_cs_cm2(329) = 2.419000d-25
  td_E_eV(330) = 7.790000d+0 ; td_cs_cm2(330) = 2.381600d-25
  td_E_eV(331) = 7.810000d+0 ; td_cs_cm2(331) = 2.344900d-25
  td_E_eV(332) = 7.830000d+0 ; td_cs_cm2(332) = 2.309100d-25
  td_E_eV(333) = 7.850000d+0 ; td_cs_cm2(333) = 2.274000d-25
  td_E_eV(334) = 7.870000d+0 ; td_cs_cm2(334) = 2.239600d-25
  td_E_eV(335) = 7.890000d+0 ; td_cs_cm2(335) = 2.205900d-25
  td_E_eV(336) = 7.910000d+0 ; td_cs_cm2(336) = 2.172800d-25
  td_E_eV(337) = 7.930000d+0 ; td_cs_cm2(337) = 2.140500d-25
  td_E_eV(338) = 7.950000d+0 ; td_cs_cm2(338) = 2.108800d-25
  td_E_eV(339) = 7.970000d+0 ; td_cs_cm2(339) = 2.077700d-25
  td_E_eV(340) = 7.990000d+0 ; td_cs_cm2(340) = 2.047300d-25
  td_E_eV(341) = 8.010000d+0 ; td_cs_cm2(341) = 2.017500d-25
  td_E_eV(342) = 8.030000d+0 ; td_cs_cm2(342) = 1.988300d-25
  td_E_eV(343) = 8.050000d+0 ; td_cs_cm2(343) = 1.959600d-25
  td_E_eV(344) = 8.070000d+0 ; td_cs_cm2(344) = 1.931500d-25
  td_E_eV(345) = 8.090000d+0 ; td_cs_cm2(345) = 1.904000d-25
  td_E_eV(346) = 8.110000d+0 ; td_cs_cm2(346) = 1.877000d-25
  td_E_eV(347) = 8.130000d+0 ; td_cs_cm2(347) = 1.850500d-25
  td_E_eV(348) = 8.150000d+0 ; td_cs_cm2(348) = 1.824500d-25
  td_E_eV(349) = 8.170000d+0 ; td_cs_cm2(349) = 1.799100d-25
  td_E_eV(350) = 8.190000d+0 ; td_cs_cm2(350) = 1.774100d-25
  td_E_eV(351) = 8.210000d+0 ; td_cs_cm2(351) = 1.749500d-25
  td_E_eV(352) = 8.230000d+0 ; td_cs_cm2(352) = 1.725400d-25
  td_E_eV(353) = 8.250000d+0 ; td_cs_cm2(353) = 1.701800d-25
  td_E_eV(354) = 8.270000d+0 ; td_cs_cm2(354) = 1.678700d-25
  td_E_eV(355) = 8.290000d+0 ; td_cs_cm2(355) = 1.656000d-25
  td_E_eV(356) = 8.310000d+0 ; td_cs_cm2(356) = 1.633700d-25
  td_E_eV(357) = 8.330000d+0 ; td_cs_cm2(357) = 1.611700d-25
  td_E_eV(358) = 8.350000d+0 ; td_cs_cm2(358) = 1.590200d-25
  td_E_eV(359) = 8.370000d+0 ; td_cs_cm2(359) = 1.569100d-25
  td_E_eV(360) = 8.390000d+0 ; td_cs_cm2(360) = 1.548400d-25
  td_E_eV(361) = 8.410000d+0 ; td_cs_cm2(361) = 1.528000d-25
  td_E_eV(362) = 8.430000d+0 ; td_cs_cm2(362) = 1.508100d-25
  td_E_eV(363) = 8.450000d+0 ; td_cs_cm2(363) = 1.488400d-25
  td_E_eV(364) = 8.470000d+0 ; td_cs_cm2(364) = 1.469100d-25
  td_E_eV(365) = 8.490000d+0 ; td_cs_cm2(365) = 1.450200d-25
  td_E_eV(366) = 8.510000d+0 ; td_cs_cm2(366) = 1.431600d-25
  td_E_eV(367) = 8.530000d+0 ; td_cs_cm2(367) = 1.413300d-25
  td_E_eV(368) = 8.550000d+0 ; td_cs_cm2(368) = 1.395400d-25
  td_E_eV(369) = 8.570000d+0 ; td_cs_cm2(369) = 1.377700d-25
  td_E_eV(370) = 8.590000d+0 ; td_cs_cm2(370) = 1.360400d-25
  td_E_eV(371) = 8.610000d+0 ; td_cs_cm2(371) = 1.343300d-25
  td_E_eV(372) = 8.630000d+0 ; td_cs_cm2(372) = 1.326600d-25
  td_E_eV(373) = 8.650000d+0 ; td_cs_cm2(373) = 1.310100d-25
  td_E_eV(374) = 8.670000d+0 ; td_cs_cm2(374) = 1.294000d-25
  td_E_eV(375) = 8.690000d+0 ; td_cs_cm2(375) = 1.278100d-25
  td_E_eV(376) = 8.710000d+0 ; td_cs_cm2(376) = 1.262400d-25
  td_E_eV(377) = 8.730000d+0 ; td_cs_cm2(377) = 1.247000d-25
  td_E_eV(378) = 8.750000d+0 ; td_cs_cm2(378) = 1.231900d-25
  td_E_eV(379) = 8.770000d+0 ; td_cs_cm2(379) = 1.217100d-25
  td_E_eV(380) = 8.790000d+0 ; td_cs_cm2(380) = 1.202500d-25
  td_E_eV(381) = 8.810000d+0 ; td_cs_cm2(381) = 1.188100d-25
  td_E_eV(382) = 8.830000d+0 ; td_cs_cm2(382) = 1.173900d-25
  td_E_eV(383) = 8.850000d+0 ; td_cs_cm2(383) = 1.160000d-25
  td_E_eV(384) = 8.870000d+0 ; td_cs_cm2(384) = 1.146400d-25
  td_E_eV(385) = 8.890000d+0 ; td_cs_cm2(385) = 1.132900d-25
  td_E_eV(386) = 8.910000d+0 ; td_cs_cm2(386) = 1.119700d-25
  td_E_eV(387) = 8.930000d+0 ; td_cs_cm2(387) = 1.106600d-25
  td_E_eV(388) = 8.950000d+0 ; td_cs_cm2(388) = 1.093800d-25
  td_E_eV(389) = 8.970000d+0 ; td_cs_cm2(389) = 1.081200d-25
  td_E_eV(390) = 8.990000d+0 ; td_cs_cm2(390) = 1.068800d-25
  td_E_eV(391) = 9.010000d+0 ; td_cs_cm2(391) = 1.056600d-25
  td_E_eV(392) = 9.030000d+0 ; td_cs_cm2(392) = 1.044600d-25
  td_E_eV(393) = 9.050000d+0 ; td_cs_cm2(393) = 1.032700d-25
  td_E_eV(394) = 9.070000d+0 ; td_cs_cm2(394) = 1.021100d-25
  td_E_eV(395) = 9.090000d+0 ; td_cs_cm2(395) = 1.009600d-25
  td_E_eV(396) = 9.110000d+0 ; td_cs_cm2(396) = 9.983500d-26
  td_E_eV(397) = 9.130000d+0 ; td_cs_cm2(397) = 9.872200d-26
  td_E_eV(398) = 9.150000d+0 ; td_cs_cm2(398) = 9.762700d-26
  td_E_eV(399) = 9.170000d+0 ; td_cs_cm2(399) = 9.655100d-26
  td_E_eV(400) = 9.190000d+0 ; td_cs_cm2(400) = 9.549100d-26
  td_E_eV(401) = 9.210000d+0 ; td_cs_cm2(401) = 9.444600d-26
  td_E_eV(402) = 9.230000d+0 ; td_cs_cm2(402) = 9.341500d-26
  td_E_eV(403) = 9.250000d+0 ; td_cs_cm2(403) = 9.240200d-26
  td_E_eV(404) = 9.270000d+0 ; td_cs_cm2(404) = 9.140500d-26
  td_E_eV(405) = 9.290000d+0 ; td_cs_cm2(405) = 9.042200d-26
  td_E_eV(406) = 9.310000d+0 ; td_cs_cm2(406) = 8.945200d-26
  td_E_eV(407) = 9.330000d+0 ; td_cs_cm2(407) = 8.849700d-26
  td_E_eV(408) = 9.350000d+0 ; td_cs_cm2(408) = 8.755800d-26
  td_E_eV(409) = 9.370000d+0 ; td_cs_cm2(409) = 8.663200d-26
  td_E_eV(410) = 9.390000d+0 ; td_cs_cm2(410) = 8.571800d-26
  td_E_eV(411) = 9.410000d+0 ; td_cs_cm2(411) = 8.481800d-26
  td_E_eV(412) = 9.430000d+0 ; td_cs_cm2(412) = 8.393300d-26
  td_E_eV(413) = 9.450000d+0 ; td_cs_cm2(413) = 8.306000d-26
  td_E_eV(414) = 9.470000d+0 ; td_cs_cm2(414) = 8.219800d-26
  td_E_eV(415) = 9.490000d+0 ; td_cs_cm2(415) = 8.135000d-26
  td_E_eV(416) = 9.510000d+0 ; td_cs_cm2(416) = 8.051500d-26
  td_E_eV(417) = 9.530000d+0 ; td_cs_cm2(417) = 7.968900d-26
  td_E_eV(418) = 9.550000d+0 ; td_cs_cm2(418) = 7.887700d-26
  td_E_eV(419) = 9.570000d+0 ; td_cs_cm2(419) = 7.807700d-26
  td_E_eV(420) = 9.590000d+0 ; td_cs_cm2(420) = 7.728500d-26
  td_E_eV(421) = 9.610000d+0 ; td_cs_cm2(421) = 7.650700d-26
  td_E_eV(422) = 9.630000d+0 ; td_cs_cm2(422) = 7.574000d-26
  td_E_eV(423) = 9.650000d+0 ; td_cs_cm2(423) = 7.498200d-26
  td_E_eV(424) = 9.670000d+0 ; td_cs_cm2(424) = 7.423600d-26
  td_E_eV(425) = 9.690000d+0 ; td_cs_cm2(425) = 7.350000d-26
  td_E_eV(426) = 9.710000d+0 ; td_cs_cm2(426) = 7.277300d-26
  td_E_eV(427) = 9.730000d+0 ; td_cs_cm2(427) = 7.205900d-26
  td_E_eV(428) = 9.750000d+0 ; td_cs_cm2(428) = 7.135100d-26
  td_E_eV(429) = 9.770000d+0 ; td_cs_cm2(429) = 7.065600d-26
  td_E_eV(430) = 9.790000d+0 ; td_cs_cm2(430) = 6.996800d-26
  td_E_eV(431) = 9.810000d+0 ; td_cs_cm2(431) = 6.929100d-26
  td_E_eV(432) = 9.830000d+0 ; td_cs_cm2(432) = 6.862100d-26
  td_E_eV(433) = 9.850000d+0 ; td_cs_cm2(433) = 6.796300d-26
  td_E_eV(434) = 9.870000d+0 ; td_cs_cm2(434) = 6.731100d-26
  td_E_eV(435) = 9.890000d+0 ; td_cs_cm2(435) = 6.666900d-26
  td_E_eV(436) = 9.910000d+0 ; td_cs_cm2(436) = 6.603600d-26
  td_E_eV(437) = 9.930000d+0 ; td_cs_cm2(437) = 6.541100d-26
  td_E_eV(438) = 9.950000d+0 ; td_cs_cm2(438) = 6.479300d-26
  td_E_eV(439) = 9.970000d+0 ; td_cs_cm2(439) = 6.418500d-26
  td_E_eV(440) = 9.990000d+0 ; td_cs_cm2(440) = 6.358400d-26
  td_E_eV(441) = 1.001000d+1 ; td_cs_cm2(441) = 6.299200d-26
  td_E_eV(442) = 1.003000d+1 ; td_cs_cm2(442) = 6.240600d-26
  td_E_eV(443) = 1.005000d+1 ; td_cs_cm2(443) = 6.182900d-26
  td_E_eV(444) = 1.007000d+1 ; td_cs_cm2(444) = 6.125900d-26
  td_E_eV(445) = 1.009000d+1 ; td_cs_cm2(445) = 6.069600d-26
  td_E_eV(446) = 1.011000d+1 ; td_cs_cm2(446) = 6.014000d-26
  td_E_eV(447) = 1.013000d+1 ; td_cs_cm2(447) = 5.959200d-26
  td_E_eV(448) = 1.015000d+1 ; td_cs_cm2(448) = 5.905000d-26
  td_E_eV(449) = 1.017000d+1 ; td_cs_cm2(449) = 5.851600d-26
  td_E_eV(450) = 1.019000d+1 ; td_cs_cm2(450) = 5.798800d-26
  td_E_eV(451) = 1.021000d+1 ; td_cs_cm2(451) = 5.746700d-26
  td_E_eV(452) = 1.023000d+1 ; td_cs_cm2(452) = 5.695200d-26
  td_E_eV(453) = 1.025000d+1 ; td_cs_cm2(453) = 5.644400d-26
  td_E_eV(454) = 1.027000d+1 ; td_cs_cm2(454) = 5.594300d-26
  td_E_eV(455) = 1.029000d+1 ; td_cs_cm2(455) = 5.544700d-26
  td_E_eV(456) = 1.031000d+1 ; td_cs_cm2(456) = 5.495800d-26
  td_E_eV(457) = 1.033000d+1 ; td_cs_cm2(457) = 5.447500d-26
  td_E_eV(458) = 1.035000d+1 ; td_cs_cm2(458) = 5.399700d-26
  td_E_eV(459) = 1.037000d+1 ; td_cs_cm2(459) = 5.352600d-26
  td_E_eV(460) = 1.039000d+1 ; td_cs_cm2(460) = 5.306000d-26
  td_E_eV(461) = 1.041000d+1 ; td_cs_cm2(461) = 5.260000d-26
  td_E_eV(462) = 1.043000d+1 ; td_cs_cm2(462) = 5.214600d-26
  td_E_eV(463) = 1.045000d+1 ; td_cs_cm2(463) = 5.169700d-26
  td_E_eV(464) = 1.047000d+1 ; td_cs_cm2(464) = 5.125300d-26
  td_E_eV(465) = 1.049000d+1 ; td_cs_cm2(465) = 5.081500d-26
  td_E_eV(466) = 1.051000d+1 ; td_cs_cm2(466) = 5.038200d-26
  td_E_eV(467) = 1.053000d+1 ; td_cs_cm2(467) = 4.995400d-26
  td_E_eV(468) = 1.055000d+1 ; td_cs_cm2(468) = 4.953200d-26
  td_E_eV(469) = 1.057000d+1 ; td_cs_cm2(469) = 4.911400d-26
  td_E_eV(470) = 1.059000d+1 ; td_cs_cm2(470) = 4.870200d-26
  td_E_eV(471) = 1.061000d+1 ; td_cs_cm2(471) = 4.829400d-26
  td_E_eV(472) = 1.063000d+1 ; td_cs_cm2(472) = 4.789100d-26
  td_E_eV(473) = 1.065000d+1 ; td_cs_cm2(473) = 4.749300d-26
  td_E_eV(474) = 1.067000d+1 ; td_cs_cm2(474) = 4.709900d-26
  td_E_eV(475) = 1.069000d+1 ; td_cs_cm2(475) = 4.671000d-26
  td_E_eV(476) = 1.071000d+1 ; td_cs_cm2(476) = 4.632600d-26
  td_E_eV(477) = 1.073000d+1 ; td_cs_cm2(477) = 4.594600d-26
  td_E_eV(478) = 1.075000d+1 ; td_cs_cm2(478) = 4.557000d-26
  td_E_eV(479) = 1.077000d+1 ; td_cs_cm2(479) = 4.519900d-26
  td_E_eV(480) = 1.079000d+1 ; td_cs_cm2(480) = 4.483200d-26
  td_E_eV(481) = 1.081000d+1 ; td_cs_cm2(481) = 4.446900d-26
  td_E_eV(482) = 1.083000d+1 ; td_cs_cm2(482) = 4.411000d-26
  td_E_eV(483) = 1.085000d+1 ; td_cs_cm2(483) = 4.375500d-26
  td_E_eV(484) = 1.087000d+1 ; td_cs_cm2(484) = 4.340500d-26
  td_E_eV(485) = 1.089000d+1 ; td_cs_cm2(485) = 4.305800d-26
  td_E_eV(486) = 1.091000d+1 ; td_cs_cm2(486) = 4.271500d-26
  td_E_eV(487) = 1.093000d+1 ; td_cs_cm2(487) = 4.237600d-26
  td_E_eV(488) = 1.095000d+1 ; td_cs_cm2(488) = 4.204100d-26
  td_E_eV(489) = 1.097000d+1 ; td_cs_cm2(489) = 4.171000d-26
  td_E_eV(490) = 1.099000d+1 ; td_cs_cm2(490) = 4.138200d-26
  td_E_eV(491) = 1.101000d+1 ; td_cs_cm2(491) = 4.105800d-26
  td_E_eV(492) = 1.103000d+1 ; td_cs_cm2(492) = 4.073700d-26
  td_E_eV(493) = 1.105000d+1 ; td_cs_cm2(493) = 4.042000d-26
  td_E_eV(494) = 1.107000d+1 ; td_cs_cm2(494) = 4.010700d-26
  td_E_eV(495) = 1.109000d+1 ; td_cs_cm2(495) = 3.979700d-26
  td_E_eV(496) = 1.111000d+1 ; td_cs_cm2(496) = 3.949000d-26
  td_E_eV(497) = 1.113000d+1 ; td_cs_cm2(497) = 3.918700d-26
  td_E_eV(498) = 1.115000d+1 ; td_cs_cm2(498) = 3.888700d-26
  td_E_eV(499) = 1.117000d+1 ; td_cs_cm2(499) = 3.859000d-26
  td_E_eV(500) = 1.119000d+1 ; td_cs_cm2(500) = 3.829600d-26
  td_E_eV(501) = 1.121000d+1 ; td_cs_cm2(501) = 3.800600d-26
  td_E_eV(502) = 1.123000d+1 ; td_cs_cm2(502) = 3.771900d-26
  td_E_eV(503) = 1.125000d+1 ; td_cs_cm2(503) = 3.743400d-26
  td_E_eV(504) = 1.127000d+1 ; td_cs_cm2(504) = 3.715300d-26
  td_E_eV(505) = 1.129000d+1 ; td_cs_cm2(505) = 3.687500d-26
  td_E_eV(506) = 1.131000d+1 ; td_cs_cm2(506) = 3.660000d-26
  td_E_eV(507) = 1.133000d+1 ; td_cs_cm2(507) = 3.632700d-26
  td_E_eV(508) = 1.135000d+1 ; td_cs_cm2(508) = 3.605800d-26
  td_E_eV(509) = 1.137000d+1 ; td_cs_cm2(509) = 3.579100d-26
  td_E_eV(510) = 1.139000d+1 ; td_cs_cm2(510) = 3.552700d-26
  td_E_eV(511) = 1.141000d+1 ; td_cs_cm2(511) = 3.526600d-26
  td_E_eV(512) = 1.143000d+1 ; td_cs_cm2(512) = 3.500700d-26
  td_E_eV(513) = 1.145000d+1 ; td_cs_cm2(513) = 3.475200d-26
  td_E_eV(514) = 1.147000d+1 ; td_cs_cm2(514) = 3.449900d-26
  td_E_eV(515) = 1.149000d+1 ; td_cs_cm2(515) = 3.424800d-26
  td_E_eV(516) = 1.151000d+1 ; td_cs_cm2(516) = 3.400000d-26
  td_E_eV(517) = 1.153000d+1 ; td_cs_cm2(517) = 3.375500d-26
  td_E_eV(518) = 1.155000d+1 ; td_cs_cm2(518) = 3.351200d-26
  td_E_eV(519) = 1.157000d+1 ; td_cs_cm2(519) = 3.327100d-26
  td_E_eV(520) = 1.159000d+1 ; td_cs_cm2(520) = 3.303400d-26
  td_E_eV(521) = 1.161000d+1 ; td_cs_cm2(521) = 3.279800d-26
  td_E_eV(522) = 1.163000d+1 ; td_cs_cm2(522) = 3.256500d-26
  td_E_eV(523) = 1.165000d+1 ; td_cs_cm2(523) = 3.233400d-26
  td_E_eV(524) = 1.167000d+1 ; td_cs_cm2(524) = 3.210600d-26
  td_E_eV(525) = 1.169000d+1 ; td_cs_cm2(525) = 3.187900d-26
  td_E_eV(526) = 1.171000d+1 ; td_cs_cm2(526) = 3.165500d-26
  td_E_eV(527) = 1.173000d+1 ; td_cs_cm2(527) = 3.143400d-26
  td_E_eV(528) = 1.175000d+1 ; td_cs_cm2(528) = 3.121400d-26
  td_E_eV(529) = 1.177000d+1 ; td_cs_cm2(529) = 3.099700d-26
  td_E_eV(530) = 1.179000d+1 ; td_cs_cm2(530) = 3.078200d-26
  td_E_eV(531) = 1.181000d+1 ; td_cs_cm2(531) = 3.056900d-26
  td_E_eV(532) = 1.183000d+1 ; td_cs_cm2(532) = 3.035800d-26
  td_E_eV(533) = 1.185000d+1 ; td_cs_cm2(533) = 3.014900d-26
  td_E_eV(534) = 1.187000d+1 ; td_cs_cm2(534) = 2.994200d-26
  td_E_eV(535) = 1.189000d+1 ; td_cs_cm2(535) = 2.973700d-26
  td_E_eV(536) = 1.191000d+1 ; td_cs_cm2(536) = 2.953400d-26
  td_E_eV(537) = 1.193000d+1 ; td_cs_cm2(537) = 2.933300d-26
  td_E_eV(538) = 1.195000d+1 ; td_cs_cm2(538) = 2.913400d-26
  td_E_eV(539) = 1.197000d+1 ; td_cs_cm2(539) = 2.893700d-26
  td_E_eV(540) = 1.199000d+1 ; td_cs_cm2(540) = 2.874200d-26
  td_E_eV(541) = 1.201000d+1 ; td_cs_cm2(541) = 2.854900d-26
  td_E_eV(542) = 1.203000d+1 ; td_cs_cm2(542) = 2.835700d-26
  td_E_eV(543) = 1.205000d+1 ; td_cs_cm2(543) = 2.816800d-26
  td_E_eV(544) = 1.207000d+1 ; td_cs_cm2(544) = 2.798000d-26
  td_E_eV(545) = 1.209000d+1 ; td_cs_cm2(545) = 2.779400d-26
  td_E_eV(546) = 1.211000d+1 ; td_cs_cm2(546) = 2.761000d-26
  td_E_eV(547) = 1.213000d+1 ; td_cs_cm2(547) = 2.742700d-26
  td_E_eV(548) = 1.215000d+1 ; td_cs_cm2(548) = 2.724600d-26
  td_E_eV(549) = 1.217000d+1 ; td_cs_cm2(549) = 2.706700d-26
  td_E_eV(550) = 1.219000d+1 ; td_cs_cm2(550) = 2.689000d-26
  td_E_eV(551) = 1.221000d+1 ; td_cs_cm2(551) = 2.671400d-26
  td_E_eV(552) = 1.223000d+1 ; td_cs_cm2(552) = 2.654000d-26
  td_E_eV(553) = 1.225000d+1 ; td_cs_cm2(553) = 2.636800d-26
  td_E_eV(554) = 1.227000d+1 ; td_cs_cm2(554) = 2.619700d-26
  td_E_eV(555) = 1.229000d+1 ; td_cs_cm2(555) = 2.602700d-26
  td_E_eV(556) = 1.231000d+1 ; td_cs_cm2(556) = 2.585900d-26
  td_E_eV(557) = 1.233000d+1 ; td_cs_cm2(557) = 2.569300d-26
  td_E_eV(558) = 1.235000d+1 ; td_cs_cm2(558) = 2.552900d-26
  td_E_eV(559) = 1.237000d+1 ; td_cs_cm2(559) = 2.536500d-26
  td_E_eV(560) = 1.239000d+1 ; td_cs_cm2(560) = 2.520400d-26
  td_E_eV(561) = 1.241000d+1 ; td_cs_cm2(561) = 2.504300d-26
  td_E_eV(562) = 1.243000d+1 ; td_cs_cm2(562) = 2.488400d-26
  td_E_eV(563) = 1.245000d+1 ; td_cs_cm2(563) = 2.472700d-26
  td_E_eV(564) = 1.247000d+1 ; td_cs_cm2(564) = 2.457100d-26
  td_E_eV(565) = 1.249000d+1 ; td_cs_cm2(565) = 2.441700d-26
  td_E_eV(566) = 1.251000d+1 ; td_cs_cm2(566) = 2.426300d-26
  td_E_eV(567) = 1.253000d+1 ; td_cs_cm2(567) = 2.411200d-26
  td_E_eV(568) = 1.255000d+1 ; td_cs_cm2(568) = 2.396100d-26
  td_E_eV(569) = 1.257000d+1 ; td_cs_cm2(569) = 2.381200d-26
  td_E_eV(570) = 1.259000d+1 ; td_cs_cm2(570) = 2.366400d-26
  td_E_eV(571) = 1.261000d+1 ; td_cs_cm2(571) = 2.351800d-26
  td_E_eV(572) = 1.263000d+1 ; td_cs_cm2(572) = 2.337300d-26
  td_E_eV(573) = 1.265000d+1 ; td_cs_cm2(573) = 2.322900d-26
  td_E_eV(574) = 1.267000d+1 ; td_cs_cm2(574) = 2.308600d-26
  td_E_eV(575) = 1.269000d+1 ; td_cs_cm2(575) = 2.294500d-26
  td_E_eV(576) = 1.271000d+1 ; td_cs_cm2(576) = 2.280500d-26
  td_E_eV(577) = 1.273000d+1 ; td_cs_cm2(577) = 2.266600d-26
  td_E_eV(578) = 1.275000d+1 ; td_cs_cm2(578) = 2.252800d-26
  td_E_eV(579) = 1.277000d+1 ; td_cs_cm2(579) = 2.239200d-26
  td_E_eV(580) = 1.279000d+1 ; td_cs_cm2(580) = 2.225600d-26
  td_E_eV(581) = 1.281000d+1 ; td_cs_cm2(581) = 2.212200d-26
  td_E_eV(582) = 1.283000d+1 ; td_cs_cm2(582) = 2.198900d-26
  td_E_eV(583) = 1.285000d+1 ; td_cs_cm2(583) = 2.185700d-26
  td_E_eV(584) = 1.287000d+1 ; td_cs_cm2(584) = 2.172600d-26
  td_E_eV(585) = 1.289000d+1 ; td_cs_cm2(585) = 2.159700d-26
  td_E_eV(586) = 1.291000d+1 ; td_cs_cm2(586) = 2.146800d-26
  td_E_eV(587) = 1.293000d+1 ; td_cs_cm2(587) = 2.134100d-26
  td_E_eV(588) = 1.295000d+1 ; td_cs_cm2(588) = 2.121500d-26
  td_E_eV(589) = 1.297000d+1 ; td_cs_cm2(589) = 2.108900d-26
  td_E_eV(590) = 1.299000d+1 ; td_cs_cm2(590) = 2.096500d-26
  td_E_eV(591) = 1.301000d+1 ; td_cs_cm2(591) = 2.084200d-26
  td_E_eV(592) = 1.303000d+1 ; td_cs_cm2(592) = 2.072000d-26
  td_E_eV(593) = 1.305000d+1 ; td_cs_cm2(593) = 2.059900d-26
  td_E_eV(594) = 1.307000d+1 ; td_cs_cm2(594) = 2.047900d-26
  td_E_eV(595) = 1.309000d+1 ; td_cs_cm2(595) = 2.036000d-26
  td_E_eV(596) = 1.311000d+1 ; td_cs_cm2(596) = 2.024100d-26
  td_E_eV(597) = 1.313000d+1 ; td_cs_cm2(597) = 2.012400d-26
  td_E_eV(598) = 1.315000d+1 ; td_cs_cm2(598) = 2.000800d-26
  td_E_eV(599) = 1.317000d+1 ; td_cs_cm2(599) = 1.989300d-26
  td_E_eV(600) = 1.319000d+1 ; td_cs_cm2(600) = 1.977900d-26
  td_E_eV(601) = 1.321000d+1 ; td_cs_cm2(601) = 1.966500d-26
  td_E_eV(602) = 1.323000d+1 ; td_cs_cm2(602) = 1.955300d-26
  td_E_eV(603) = 1.325000d+1 ; td_cs_cm2(603) = 1.944200d-26
  td_E_eV(604) = 1.327000d+1 ; td_cs_cm2(604) = 1.933100d-26
  td_E_eV(605) = 1.329000d+1 ; td_cs_cm2(605) = 1.922200d-26
  td_E_eV(606) = 1.331000d+1 ; td_cs_cm2(606) = 1.911300d-26
  td_E_eV(607) = 1.333000d+1 ; td_cs_cm2(607) = 1.900500d-26
  td_E_eV(608) = 1.335000d+1 ; td_cs_cm2(608) = 1.889800d-26
  td_E_eV(609) = 1.337000d+1 ; td_cs_cm2(609) = 1.879200d-26
  td_E_eV(610) = 1.339000d+1 ; td_cs_cm2(610) = 1.868700d-26
  td_E_eV(611) = 1.341000d+1 ; td_cs_cm2(611) = 1.858200d-26
  td_E_eV(612) = 1.343000d+1 ; td_cs_cm2(612) = 1.847800d-26
  td_E_eV(613) = 1.345000d+1 ; td_cs_cm2(613) = 1.837600d-26
  td_E_eV(614) = 1.347000d+1 ; td_cs_cm2(614) = 1.827400d-26
  td_E_eV(615) = 1.349000d+1 ; td_cs_cm2(615) = 1.817300d-26
  td_E_eV(616) = 1.351000d+1 ; td_cs_cm2(616) = 1.807200d-26
  td_E_eV(617) = 1.353000d+1 ; td_cs_cm2(617) = 1.797300d-26
  td_E_eV(618) = 1.355000d+1 ; td_cs_cm2(618) = 1.787400d-26
  td_E_eV(619) = 1.357000d+1 ; td_cs_cm2(619) = 1.777600d-26
  td_E_eV(620) = 1.359000d+1 ; td_cs_cm2(620) = 1.767900d-26
  td_E_eV(621) = 1.361000d+1 ; td_cs_cm2(621) = 1.758200d-26
  td_E_eV(622) = 1.363000d+1 ; td_cs_cm2(622) = 1.748700d-26
  td_E_eV(623) = 1.365000d+1 ; td_cs_cm2(623) = 1.739200d-26
  td_E_eV(624) = 1.367000d+1 ; td_cs_cm2(624) = 1.729800d-26
  td_E_eV(625) = 1.369000d+1 ; td_cs_cm2(625) = 1.720400d-26
  td_E_eV(626) = 1.371000d+1 ; td_cs_cm2(626) = 1.711100d-26
  td_E_eV(627) = 1.373000d+1 ; td_cs_cm2(627) = 1.701900d-26
  td_E_eV(628) = 1.375000d+1 ; td_cs_cm2(628) = 1.692800d-26
  td_E_eV(629) = 1.377000d+1 ; td_cs_cm2(629) = 1.683800d-26
  td_E_eV(630) = 1.379000d+1 ; td_cs_cm2(630) = 1.674800d-26
  td_E_eV(631) = 1.381000d+1 ; td_cs_cm2(631) = 1.665800d-26
  td_E_eV(632) = 1.383000d+1 ; td_cs_cm2(632) = 1.657000d-26
  td_E_eV(633) = 1.385000d+1 ; td_cs_cm2(633) = 1.648200d-26
  td_E_eV(634) = 1.387000d+1 ; td_cs_cm2(634) = 1.639500d-26
  td_E_eV(635) = 1.389000d+1 ; td_cs_cm2(635) = 1.630800d-26
  td_E_eV(636) = 1.391000d+1 ; td_cs_cm2(636) = 1.622200d-26
  td_E_eV(637) = 1.393000d+1 ; td_cs_cm2(637) = 1.613700d-26
  td_E_eV(638) = 1.395000d+1 ; td_cs_cm2(638) = 1.605300d-26
  td_E_eV(639) = 1.397000d+1 ; td_cs_cm2(639) = 1.596900d-26
  td_E_eV(640) = 1.399000d+1 ; td_cs_cm2(640) = 1.588500d-26
  td_E_eV(641) = 1.401000d+1 ; td_cs_cm2(641) = 1.580300d-26
  td_E_eV(642) = 1.403000d+1 ; td_cs_cm2(642) = 1.572100d-26
  td_E_eV(643) = 1.405000d+1 ; td_cs_cm2(643) = 1.563900d-26
  td_E_eV(644) = 1.407000d+1 ; td_cs_cm2(644) = 1.555800d-26
  td_E_eV(645) = 1.409000d+1 ; td_cs_cm2(645) = 1.547800d-26
  td_E_eV(646) = 1.411000d+1 ; td_cs_cm2(646) = 1.539900d-26
  td_E_eV(647) = 1.413000d+1 ; td_cs_cm2(647) = 1.532000d-26
  td_E_eV(648) = 1.415000d+1 ; td_cs_cm2(648) = 1.524100d-26
  td_E_eV(649) = 1.417000d+1 ; td_cs_cm2(649) = 1.516300d-26
  td_E_eV(650) = 1.419000d+1 ; td_cs_cm2(650) = 1.508600d-26
  td_E_eV(651) = 1.421000d+1 ; td_cs_cm2(651) = 1.500900d-26
  td_E_eV(652) = 1.423000d+1 ; td_cs_cm2(652) = 1.493300d-26
  td_E_eV(653) = 1.425000d+1 ; td_cs_cm2(653) = 1.485700d-26
  td_E_eV(654) = 1.427000d+1 ; td_cs_cm2(654) = 1.478200d-26
  td_E_eV(655) = 1.429000d+1 ; td_cs_cm2(655) = 1.470800d-26
  td_E_eV(656) = 1.431000d+1 ; td_cs_cm2(656) = 1.463400d-26
  td_E_eV(657) = 1.433000d+1 ; td_cs_cm2(657) = 1.456000d-26
  td_E_eV(658) = 1.435000d+1 ; td_cs_cm2(658) = 1.448700d-26
  td_E_eV(659) = 1.437000d+1 ; td_cs_cm2(659) = 1.441500d-26
  td_E_eV(660) = 1.439000d+1 ; td_cs_cm2(660) = 1.434300d-26
  td_E_eV(661) = 1.441000d+1 ; td_cs_cm2(661) = 1.427200d-26
  td_E_eV(662) = 1.443000d+1 ; td_cs_cm2(662) = 1.420100d-26
  td_E_eV(663) = 1.445000d+1 ; td_cs_cm2(663) = 1.413000d-26
  td_E_eV(664) = 1.447000d+1 ; td_cs_cm2(664) = 1.406100d-26
  td_E_eV(665) = 1.449000d+1 ; td_cs_cm2(665) = 1.399100d-26
  td_E_eV(666) = 1.451000d+1 ; td_cs_cm2(666) = 1.392200d-26
  td_E_eV(667) = 1.453000d+1 ; td_cs_cm2(667) = 1.385400d-26
  td_E_eV(668) = 1.455000d+1 ; td_cs_cm2(668) = 1.378600d-26
  td_E_eV(669) = 1.457000d+1 ; td_cs_cm2(669) = 1.371900d-26
  td_E_eV(670) = 1.459000d+1 ; td_cs_cm2(670) = 1.365200d-26
  td_E_eV(671) = 1.461000d+1 ; td_cs_cm2(671) = 1.358500d-26
  td_E_eV(672) = 1.463000d+1 ; td_cs_cm2(672) = 1.351900d-26
  td_E_eV(673) = 1.465000d+1 ; td_cs_cm2(673) = 1.345400d-26
  td_E_eV(674) = 1.467000d+1 ; td_cs_cm2(674) = 1.338900d-26
  td_E_eV(675) = 1.469000d+1 ; td_cs_cm2(675) = 1.332400d-26
  td_E_eV(676) = 1.471000d+1 ; td_cs_cm2(676) = 1.326000d-26
  td_E_eV(677) = 1.473000d+1 ; td_cs_cm2(677) = 1.319600d-26
  td_E_eV(678) = 1.475000d+1 ; td_cs_cm2(678) = 1.313300d-26
  td_E_eV(679) = 1.477000d+1 ; td_cs_cm2(679) = 1.307000d-26
  td_E_eV(680) = 1.479000d+1 ; td_cs_cm2(680) = 1.300800d-26
  td_E_eV(681) = 1.481000d+1 ; td_cs_cm2(681) = 1.294600d-26
  td_E_eV(682) = 1.483000d+1 ; td_cs_cm2(682) = 1.288400d-26
  td_E_eV(683) = 1.485000d+1 ; td_cs_cm2(683) = 1.282300d-26
  td_E_eV(684) = 1.487000d+1 ; td_cs_cm2(684) = 1.276300d-26
  td_E_eV(685) = 1.489000d+1 ; td_cs_cm2(685) = 1.270200d-26
  td_E_eV(686) = 1.491000d+1 ; td_cs_cm2(686) = 1.264200d-26
  td_E_eV(687) = 1.493000d+1 ; td_cs_cm2(687) = 1.258300d-26
  td_E_eV(688) = 1.495000d+1 ; td_cs_cm2(688) = 1.252400d-26
  td_E_eV(689) = 1.497000d+1 ; td_cs_cm2(689) = 1.246500d-26
  td_E_eV(690) = 1.499000d+1 ; td_cs_cm2(690) = 1.240700d-26

  td_cs_cm2 = 1.0d4 * td_cs_cm2 ! the original data are in m^2, convert to cm^2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_vib03
!
!-------------------------------------
!
real(8) function CSV_N2_vib03_m3s(E_eV)

  use N2_vib03
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_vib03_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_vib03_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_vib03_m3s











!---------------------------------------------------------------------------------------------------- #7
!#VIBRATIONAL
!#SPECIES: e / N2
!#PROCESS: E + N2 <-> E + N2*(v=0-4), Vibrational
!#PARAM.:  E = 1.133 eV, g1/g0 = 1, complete set
!#COMMENT: Vibrational excitation : N2(1Σg+:v=0) -> N2*(v=4) | Source: [Laporta et al, Plasma
!#COMMENT: Sources Sci. Technol. 23, 065002 (2014)].
!#UPDATED: 2023-12-23 21:19:17
!
module N2_vib04
  integer, parameter :: N_td=748  !8   ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_vib04
!
!-------------------------------------
!
subroutine Prepare_N2_vib04

  use N2_vib04
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

!  td_E_eV(1) = 1.3_8; td_cs_cm2(1) = 1.11d-18
!  td_E_eV(2) = 1.6_8; td_cs_cm2(2) = 8.84d-17
!  td_E_eV(3) = 2.0_8; td_cs_cm2(3) = 1.01d-15
!  td_E_eV(4) = 2.4_8; td_cs_cm2(4) = 1.35d-15
!  td_E_eV(5) = 3.0_8; td_cs_cm2(5) = 5.86d-16
!  td_E_eV(6) = 4.0_8; td_cs_cm2(6) = 1.37d-16
!  td_E_eV(7) = 5.0_8; td_cs_cm2(7) = 2.34d-17
!  td_E_eV(8) = 6.0_8; td_cs_cm2(8) = 3.66d-18

  td_E_eV(  1) = 1.420000d+0 ; td_cs_cm2(  1) = 0.000000d+0
  td_E_eV(  2) = 1.430000d+0 ; td_cs_cm2(  2) = 1.138800d-26
  td_E_eV(  3) = 1.440000d+0 ; td_cs_cm2(  3) = 1.362700d-26
  td_E_eV(  4) = 1.450000d+0 ; td_cs_cm2(  4) = 1.632100d-26
  td_E_eV(  5) = 1.460000d+0 ; td_cs_cm2(  5) = 1.956600d-26
  td_E_eV(  6) = 1.470000d+0 ; td_cs_cm2(  6) = 2.347800d-26
  td_E_eV(  7) = 1.480000d+0 ; td_cs_cm2(  7) = 2.820000d-26
  td_E_eV(  8) = 1.490000d+0 ; td_cs_cm2(  8) = 3.390400d-26
  td_E_eV(  9) = 1.500000d+0 ; td_cs_cm2(  9) = 4.080300d-26
  td_E_eV( 10) = 1.510000d+0 ; td_cs_cm2( 10) = 4.915600d-26
  td_E_eV( 11) = 1.520000d+0 ; td_cs_cm2( 11) = 5.928300d-26
  td_E_eV( 12) = 1.530000d+0 ; td_cs_cm2( 12) = 7.157400d-26
  td_E_eV( 13) = 1.540000d+0 ; td_cs_cm2( 13) = 8.651300d-26
  td_E_eV( 14) = 1.550000d+0 ; td_cs_cm2( 14) = 1.047000d-25
  td_E_eV( 15) = 1.560000d+0 ; td_cs_cm2( 15) = 1.268600d-25
  td_E_eV( 16) = 1.570000d+0 ; td_cs_cm2( 16) = 1.539200d-25
  td_E_eV( 17) = 1.580000d+0 ; td_cs_cm2( 17) = 1.870100d-25
  td_E_eV( 18) = 1.590000d+0 ; td_cs_cm2( 18) = 2.275400d-25
  td_E_eV( 19) = 1.600000d+0 ; td_cs_cm2( 19) = 2.772900d-25
  td_E_eV( 20) = 1.610000d+0 ; td_cs_cm2( 20) = 3.384500d-25
  td_E_eV( 21) = 1.620000d+0 ; td_cs_cm2( 21) = 4.138200d-25
  td_E_eV( 22) = 1.630000d+0 ; td_cs_cm2( 22) = 5.068700d-25
  td_E_eV( 23) = 1.640000d+0 ; td_cs_cm2( 23) = 6.220200d-25
  td_E_eV( 24) = 1.650000d+0 ; td_cs_cm2( 24) = 7.648700d-25
  td_E_eV( 25) = 1.660000d+0 ; td_cs_cm2( 25) = 9.425100d-25
  td_E_eV( 26) = 1.670000d+0 ; td_cs_cm2( 26) = 1.164000d-24
  td_E_eV( 27) = 1.680000d+0 ; td_cs_cm2( 27) = 1.440900d-24
  td_E_eV( 28) = 1.690000d+0 ; td_cs_cm2( 28) = 1.788200d-24
  td_E_eV( 29) = 1.700000d+0 ; td_cs_cm2( 29) = 2.225000d-24
  td_E_eV( 30) = 1.710000d+0 ; td_cs_cm2( 30) = 2.776100d-24
  td_E_eV( 31) = 1.720000d+0 ; td_cs_cm2( 31) = 3.473800d-24
  td_E_eV( 32) = 1.730000d+0 ; td_cs_cm2( 32) = 4.360300d-24
  td_E_eV( 33) = 1.740000d+0 ; td_cs_cm2( 33) = 5.490700d-24
  td_E_eV( 34) = 1.750000d+0 ; td_cs_cm2( 34) = 6.937500d-24
  td_E_eV( 35) = 1.760000d+0 ; td_cs_cm2( 35) = 8.796500d-24
  td_E_eV( 36) = 1.770000d+0 ; td_cs_cm2( 36) = 1.119500d-23
  td_E_eV( 37) = 1.780000d+0 ; td_cs_cm2( 37) = 1.430100d-23
  td_E_eV( 38) = 1.790000d+0 ; td_cs_cm2( 38) = 1.834200d-23
  td_E_eV( 39) = 1.800000d+0 ; td_cs_cm2( 39) = 2.361800d-23
  td_E_eV( 40) = 1.810000d+0 ; td_cs_cm2( 40) = 3.053400d-23
  td_E_eV( 41) = 1.820000d+0 ; td_cs_cm2( 41) = 3.963200d-23
  td_E_eV( 42) = 1.830000d+0 ; td_cs_cm2( 42) = 5.163800d-23
  td_E_eV( 43) = 1.840000d+0 ; td_cs_cm2( 43) = 6.751500d-23
  td_E_eV( 44) = 1.850000d+0 ; td_cs_cm2( 44) = 8.853700d-23
  td_E_eV( 45) = 1.860000d+0 ; td_cs_cm2( 45) = 1.163600d-22
  td_E_eV( 46) = 1.870000d+0 ; td_cs_cm2( 46) = 1.530900d-22
  td_E_eV( 47) = 1.880000d+0 ; td_cs_cm2( 47) = 2.013100d-22
  td_E_eV( 48) = 1.890000d+0 ; td_cs_cm2( 48) = 2.640400d-22
  td_E_eV( 49) = 1.900000d+0 ; td_cs_cm2( 49) = 3.445500d-22
  td_E_eV( 50) = 1.910000d+0 ; td_cs_cm2( 50) = 4.459000d-22
  td_E_eV( 51) = 1.920000d+0 ; td_cs_cm2( 51) = 5.703900d-22
  td_E_eV( 52) = 1.930000d+0 ; td_cs_cm2( 52) = 7.187700d-22
  td_E_eV( 53) = 1.940000d+0 ; td_cs_cm2( 53) = 8.896900d-22
  td_E_eV( 54) = 1.950000d+0 ; td_cs_cm2( 54) = 1.079700d-21
  td_E_eV( 55) = 1.960000d+0 ; td_cs_cm2( 55) = 1.284100d-21
  td_E_eV( 56) = 1.970000d+0 ; td_cs_cm2( 56) = 1.498200d-21
  td_E_eV( 57) = 1.980000d+0 ; td_cs_cm2( 57) = 1.719100d-21
  td_E_eV( 58) = 1.990000d+0 ; td_cs_cm2( 58) = 1.946100d-21
  td_E_eV( 59) = 2.000000d+0 ; td_cs_cm2( 59) = 2.181800d-21
  td_E_eV( 60) = 2.010000d+0 ; td_cs_cm2( 60) = 2.430800d-21
  td_E_eV( 61) = 2.020000d+0 ; td_cs_cm2( 61) = 2.699900d-21
  td_E_eV( 62) = 2.030000d+0 ; td_cs_cm2( 62) = 2.997400d-21
  td_E_eV( 63) = 2.040000d+0 ; td_cs_cm2( 63) = 3.332900d-21
  td_E_eV( 64) = 2.050000d+0 ; td_cs_cm2( 64) = 3.716900d-21
  td_E_eV( 65) = 2.060000d+0 ; td_cs_cm2( 65) = 4.161400d-21
  td_E_eV( 66) = 2.070000d+0 ; td_cs_cm2( 66) = 4.679200d-21
  td_E_eV( 67) = 2.080000d+0 ; td_cs_cm2( 67) = 5.284200d-21
  td_E_eV( 68) = 2.090000d+0 ; td_cs_cm2( 68) = 5.990200d-21
  td_E_eV( 69) = 2.100000d+0 ; td_cs_cm2( 69) = 6.809800d-21
  td_E_eV( 70) = 2.110000d+0 ; td_cs_cm2( 70) = 7.752100d-21
  td_E_eV( 71) = 2.120000d+0 ; td_cs_cm2( 71) = 8.818900d-21
  td_E_eV( 72) = 2.130000d+0 ; td_cs_cm2( 72) = 9.999900d-21
  td_E_eV( 73) = 2.140000d+0 ; td_cs_cm2( 73) = 1.126800d-20
  td_E_eV( 74) = 2.150000d+0 ; td_cs_cm2( 74) = 1.257400d-20
  td_E_eV( 75) = 2.160000d+0 ; td_cs_cm2( 75) = 1.385100d-20
  td_E_eV( 76) = 2.170000d+0 ; td_cs_cm2( 76) = 1.501800d-20
  td_E_eV( 77) = 2.180000d+0 ; td_cs_cm2( 77) = 1.599600d-20
  td_E_eV( 78) = 2.190000d+0 ; td_cs_cm2( 78) = 1.672700d-20
  td_E_eV( 79) = 2.200000d+0 ; td_cs_cm2( 79) = 1.718700d-20
  td_E_eV( 80) = 2.210000d+0 ; td_cs_cm2( 80) = 1.738600d-20
  td_E_eV( 81) = 2.220000d+0 ; td_cs_cm2( 81) = 1.736600d-20
  td_E_eV( 82) = 2.230000d+0 ; td_cs_cm2( 82) = 1.718000d-20
  td_E_eV( 83) = 2.240000d+0 ; td_cs_cm2( 83) = 1.688800d-20
  td_E_eV( 84) = 2.250000d+0 ; td_cs_cm2( 84) = 1.653900d-20
  td_E_eV( 85) = 2.260000d+0 ; td_cs_cm2( 85) = 1.617300d-20
  td_E_eV( 86) = 2.270000d+0 ; td_cs_cm2( 86) = 1.581800d-20
  td_E_eV( 87) = 2.280000d+0 ; td_cs_cm2( 87) = 1.549000d-20
  td_E_eV( 88) = 2.290000d+0 ; td_cs_cm2( 88) = 1.519400d-20
  td_E_eV( 89) = 2.300000d+0 ; td_cs_cm2( 89) = 1.492800d-20
  td_E_eV( 90) = 2.310000d+0 ; td_cs_cm2( 90) = 1.468300d-20
  td_E_eV( 91) = 2.320000d+0 ; td_cs_cm2( 91) = 1.444300d-20
  td_E_eV( 92) = 2.330000d+0 ; td_cs_cm2( 92) = 1.418500d-20
  td_E_eV( 93) = 2.340000d+0 ; td_cs_cm2( 93) = 1.387700d-20
  td_E_eV( 94) = 2.350000d+0 ; td_cs_cm2( 94) = 1.348600d-20
  td_E_eV( 95) = 2.360000d+0 ; td_cs_cm2( 95) = 1.297300d-20
  td_E_eV( 96) = 2.370000d+0 ; td_cs_cm2( 96) = 1.230300d-20
  td_E_eV( 97) = 2.380000d+0 ; td_cs_cm2( 97) = 1.145400d-20
  td_E_eV( 98) = 2.390000d+0 ; td_cs_cm2( 98) = 1.042200d-20
  td_E_eV( 99) = 2.400000d+0 ; td_cs_cm2( 99) = 9.229600d-21
  td_E_eV(100) = 2.410000d+0 ; td_cs_cm2(100) = 7.928400d-21
  td_E_eV(101) = 2.420000d+0 ; td_cs_cm2(101) = 6.587700d-21
  td_E_eV(102) = 2.430000d+0 ; td_cs_cm2(102) = 5.282500d-21
  td_E_eV(103) = 2.440000d+0 ; td_cs_cm2(103) = 4.078200d-21
  td_E_eV(104) = 2.450000d+0 ; td_cs_cm2(104) = 3.021900d-21
  td_E_eV(105) = 2.460000d+0 ; td_cs_cm2(105) = 2.139100d-21
  td_E_eV(106) = 2.470000d+0 ; td_cs_cm2(106) = 1.437500d-21
  td_E_eV(107) = 2.480000d+0 ; td_cs_cm2(107) = 9.125800d-22
  td_E_eV(108) = 2.490000d+0 ; td_cs_cm2(108) = 5.534700d-22
  td_E_eV(109) = 2.500000d+0 ; td_cs_cm2(109) = 3.478500d-22
  td_E_eV(110) = 2.510000d+0 ; td_cs_cm2(110) = 2.848000d-22
  td_E_eV(111) = 2.520000d+0 ; td_cs_cm2(111) = 3.564700d-22
  td_E_eV(112) = 2.530000d+0 ; td_cs_cm2(112) = 5.585300d-22
  td_E_eV(113) = 2.540000d+0 ; td_cs_cm2(113) = 8.897200d-22
  td_E_eV(114) = 2.550000d+0 ; td_cs_cm2(114) = 1.350500d-21
  td_E_eV(115) = 2.560000d+0 ; td_cs_cm2(115) = 1.941000d-21
  td_E_eV(116) = 2.570000d+0 ; td_cs_cm2(116) = 2.657700d-21
  td_E_eV(117) = 2.580000d+0 ; td_cs_cm2(117) = 3.489400d-21
  td_E_eV(118) = 2.590000d+0 ; td_cs_cm2(118) = 4.413100d-21
  td_E_eV(119) = 2.600000d+0 ; td_cs_cm2(119) = 5.390900d-21
  td_E_eV(120) = 2.610000d+0 ; td_cs_cm2(120) = 6.369900d-21
  td_E_eV(121) = 2.620000d+0 ; td_cs_cm2(121) = 7.287300d-21
  td_E_eV(122) = 2.630000d+0 ; td_cs_cm2(122) = 8.080100d-21
  td_E_eV(123) = 2.640000d+0 ; td_cs_cm2(123) = 8.697100d-21
  td_E_eV(124) = 2.650000d+0 ; td_cs_cm2(124) = 9.108600d-21
  td_E_eV(125) = 2.660000d+0 ; td_cs_cm2(125) = 9.310100d-21
  td_E_eV(126) = 2.670000d+0 ; td_cs_cm2(126) = 9.318700d-21
  td_E_eV(127) = 2.680000d+0 ; td_cs_cm2(127) = 9.165100d-21
  td_E_eV(128) = 2.690000d+0 ; td_cs_cm2(128) = 8.885200d-21
  td_E_eV(129) = 2.700000d+0 ; td_cs_cm2(129) = 8.513100d-21
  td_E_eV(130) = 2.710000d+0 ; td_cs_cm2(130) = 8.077300d-21
  td_E_eV(131) = 2.720000d+0 ; td_cs_cm2(131) = 7.599300d-21
  td_E_eV(132) = 2.730000d+0 ; td_cs_cm2(132) = 7.093600d-21
  td_E_eV(133) = 2.740000d+0 ; td_cs_cm2(133) = 6.568700d-21
  td_E_eV(134) = 2.750000d+0 ; td_cs_cm2(134) = 6.028800d-21
  td_E_eV(135) = 2.760000d+0 ; td_cs_cm2(135) = 5.475200d-21
  td_E_eV(136) = 2.770000d+0 ; td_cs_cm2(136) = 4.907800d-21
  td_E_eV(137) = 2.780000d+0 ; td_cs_cm2(137) = 4.326900d-21
  td_E_eV(138) = 2.790000d+0 ; td_cs_cm2(138) = 3.735000d-21
  td_E_eV(139) = 2.800000d+0 ; td_cs_cm2(139) = 3.138300d-21
  td_E_eV(140) = 2.810000d+0 ; td_cs_cm2(140) = 2.548800d-21
  td_E_eV(141) = 2.820000d+0 ; td_cs_cm2(141) = 1.984300d-21
  td_E_eV(142) = 2.830000d+0 ; td_cs_cm2(142) = 1.467500d-21
  td_E_eV(143) = 2.840000d+0 ; td_cs_cm2(143) = 1.023200d-21
  td_E_eV(144) = 2.850000d+0 ; td_cs_cm2(144) = 6.736200d-22
  td_E_eV(145) = 2.860000d+0 ; td_cs_cm2(145) = 4.337500d-22
  td_E_eV(146) = 2.870000d+0 ; td_cs_cm2(146) = 3.083400d-22
  td_E_eV(147) = 2.880000d+0 ; td_cs_cm2(147) = 2.918600d-22
  td_E_eV(148) = 2.890000d+0 ; td_cs_cm2(148) = 3.708300d-22
  td_E_eV(149) = 2.900000d+0 ; td_cs_cm2(149) = 5.275700d-22
  td_E_eV(150) = 2.910000d+0 ; td_cs_cm2(150) = 7.436900d-22
  td_E_eV(151) = 2.920000d+0 ; td_cs_cm2(151) = 1.002500d-21
  td_E_eV(152) = 2.930000d+0 ; td_cs_cm2(152) = 1.290300d-21
  td_E_eV(153) = 2.940000d+0 ; td_cs_cm2(153) = 1.595900d-21
  td_E_eV(154) = 2.950000d+0 ; td_cs_cm2(154) = 1.911000d-21
  td_E_eV(155) = 2.960000d+0 ; td_cs_cm2(155) = 2.228400d-21
  td_E_eV(156) = 2.970000d+0 ; td_cs_cm2(156) = 2.541300d-21
  td_E_eV(157) = 2.980000d+0 ; td_cs_cm2(157) = 2.842300d-21
  td_E_eV(158) = 2.990000d+0 ; td_cs_cm2(158) = 3.122600d-21
  td_E_eV(159) = 3.000000d+0 ; td_cs_cm2(159) = 3.371200d-21
  td_E_eV(160) = 3.010000d+0 ; td_cs_cm2(160) = 3.574600d-21
  td_E_eV(161) = 3.020000d+0 ; td_cs_cm2(161) = 3.717600d-21
  td_E_eV(162) = 3.030000d+0 ; td_cs_cm2(162) = 3.785100d-21
  td_E_eV(163) = 3.040000d+0 ; td_cs_cm2(163) = 3.764100d-21
  td_E_eV(164) = 3.050000d+0 ; td_cs_cm2(164) = 3.647700d-21
  td_E_eV(165) = 3.060000d+0 ; td_cs_cm2(165) = 3.437800d-21
  td_E_eV(166) = 3.070000d+0 ; td_cs_cm2(166) = 3.146300d-21
  td_E_eV(167) = 3.080000d+0 ; td_cs_cm2(167) = 2.794000d-21
  td_E_eV(168) = 3.090000d+0 ; td_cs_cm2(168) = 2.406500d-21
  td_E_eV(169) = 3.100000d+0 ; td_cs_cm2(169) = 2.010500d-21
  td_E_eV(170) = 3.110000d+0 ; td_cs_cm2(170) = 1.629000d-21
  td_E_eV(171) = 3.120000d+0 ; td_cs_cm2(171) = 1.279400d-21
  td_E_eV(172) = 3.130000d+0 ; td_cs_cm2(172) = 9.727900d-22
  td_E_eV(173) = 3.140000d+0 ; td_cs_cm2(173) = 7.152400d-22
  td_E_eV(174) = 3.150000d+0 ; td_cs_cm2(174) = 5.087300d-22
  td_E_eV(175) = 3.160000d+0 ; td_cs_cm2(175) = 3.528200d-22
  td_E_eV(176) = 3.170000d+0 ; td_cs_cm2(176) = 2.457900d-22
  td_E_eV(177) = 3.180000d+0 ; td_cs_cm2(177) = 1.854000d-22
  td_E_eV(178) = 3.190000d+0 ; td_cs_cm2(178) = 1.692800d-22
  td_E_eV(179) = 3.200000d+0 ; td_cs_cm2(179) = 1.949300d-22
  td_E_eV(180) = 3.210000d+0 ; td_cs_cm2(180) = 2.594800d-22
  td_E_eV(181) = 3.220000d+0 ; td_cs_cm2(181) = 3.591100d-22
  td_E_eV(182) = 3.230000d+0 ; td_cs_cm2(182) = 4.884300d-22
  td_E_eV(183) = 3.240000d+0 ; td_cs_cm2(183) = 6.398100d-22
  td_E_eV(184) = 3.250000d+0 ; td_cs_cm2(184) = 8.031000d-22
  td_E_eV(185) = 3.260000d+0 ; td_cs_cm2(185) = 9.660100d-22
  td_E_eV(186) = 3.270000d+0 ; td_cs_cm2(186) = 1.115300d-21
  td_E_eV(187) = 3.280000d+0 ; td_cs_cm2(187) = 1.238700d-21
  td_E_eV(188) = 3.290000d+0 ; td_cs_cm2(188) = 1.326800d-21
  td_E_eV(189) = 3.300000d+0 ; td_cs_cm2(189) = 1.374700d-21
  td_E_eV(190) = 3.310000d+0 ; td_cs_cm2(190) = 1.382100d-21
  td_E_eV(191) = 3.320000d+0 ; td_cs_cm2(191) = 1.352700d-21
  td_E_eV(192) = 3.330000d+0 ; td_cs_cm2(192) = 1.292600d-21
  td_E_eV(193) = 3.340000d+0 ; td_cs_cm2(193) = 1.209000d-21
  td_E_eV(194) = 3.350000d+0 ; td_cs_cm2(194) = 1.109000d-21
  td_E_eV(195) = 3.360000d+0 ; td_cs_cm2(195) = 9.986500d-22
  td_E_eV(196) = 3.370000d+0 ; td_cs_cm2(196) = 8.831700d-22
  td_E_eV(197) = 3.380000d+0 ; td_cs_cm2(197) = 7.666100d-22
  td_E_eV(198) = 3.390000d+0 ; td_cs_cm2(198) = 6.522400d-22
  td_E_eV(199) = 3.400000d+0 ; td_cs_cm2(199) = 5.428000d-22
  td_E_eV(200) = 3.410000d+0 ; td_cs_cm2(200) = 4.408100d-22
  td_E_eV(201) = 3.420000d+0 ; td_cs_cm2(201) = 3.487400d-22
  td_E_eV(202) = 3.430000d+0 ; td_cs_cm2(202) = 2.691300d-22
  td_E_eV(203) = 3.440000d+0 ; td_cs_cm2(203) = 2.045300d-22
  td_E_eV(204) = 3.450000d+0 ; td_cs_cm2(204) = 1.572800d-22
  td_E_eV(205) = 3.460000d+0 ; td_cs_cm2(205) = 1.290700d-22
  td_E_eV(206) = 3.470000d+0 ; td_cs_cm2(206) = 1.204500d-22
  td_E_eV(207) = 3.480000d+0 ; td_cs_cm2(207) = 1.304200d-22
  td_E_eV(208) = 3.490000d+0 ; td_cs_cm2(208) = 1.562700d-22
  td_E_eV(209) = 3.500000d+0 ; td_cs_cm2(209) = 1.939100d-22
  td_E_eV(210) = 3.510000d+0 ; td_cs_cm2(210) = 2.385200d-22
  td_E_eV(211) = 3.520000d+0 ; td_cs_cm2(211) = 2.853000d-22
  td_E_eV(212) = 3.530000d+0 ; td_cs_cm2(212) = 3.301400d-22
  td_E_eV(213) = 3.540000d+0 ; td_cs_cm2(213) = 3.699400d-22
  td_E_eV(214) = 3.550000d+0 ; td_cs_cm2(214) = 4.027200d-22
  td_E_eV(215) = 3.560000d+0 ; td_cs_cm2(215) = 4.274000d-22
  td_E_eV(216) = 3.570000d+0 ; td_cs_cm2(216) = 4.435800d-22
  td_E_eV(217) = 3.580000d+0 ; td_cs_cm2(217) = 4.512700d-22
  td_E_eV(218) = 3.590000d+0 ; td_cs_cm2(218) = 4.506900d-22
  td_E_eV(219) = 3.600000d+0 ; td_cs_cm2(219) = 4.421800d-22
  td_E_eV(220) = 3.610000d+0 ; td_cs_cm2(220) = 4.261000d-22
  td_E_eV(221) = 3.620000d+0 ; td_cs_cm2(221) = 4.029100d-22
  td_E_eV(222) = 3.630000d+0 ; td_cs_cm2(222) = 3.731900d-22
  td_E_eV(223) = 3.640000d+0 ; td_cs_cm2(223) = 3.378300d-22
  td_E_eV(224) = 3.650000d+0 ; td_cs_cm2(224) = 2.980800d-22
  td_E_eV(225) = 3.660000d+0 ; td_cs_cm2(225) = 2.556800d-22
  td_E_eV(226) = 3.670000d+0 ; td_cs_cm2(226) = 2.128500d-22
  td_E_eV(227) = 3.680000d+0 ; td_cs_cm2(227) = 1.720700d-22
  td_E_eV(228) = 3.690000d+0 ; td_cs_cm2(228) = 1.358100d-22
  td_E_eV(229) = 3.700000d+0 ; td_cs_cm2(229) = 1.061500d-22
  td_E_eV(230) = 3.710000d+0 ; td_cs_cm2(230) = 8.440300d-23
  td_E_eV(231) = 3.720000d+0 ; td_cs_cm2(231) = 7.102100d-23
  td_E_eV(232) = 3.730000d+0 ; td_cs_cm2(232) = 6.560700d-23
  td_E_eV(233) = 3.740000d+0 ; td_cs_cm2(233) = 6.714800d-23
  td_E_eV(234) = 3.750000d+0 ; td_cs_cm2(234) = 7.427900d-23
  td_E_eV(235) = 3.760000d+0 ; td_cs_cm2(235) = 8.553300d-23
  td_E_eV(236) = 3.770000d+0 ; td_cs_cm2(236) = 9.950000d-23
  td_E_eV(237) = 3.780000d+0 ; td_cs_cm2(237) = 1.149200d-22
  td_E_eV(238) = 3.790000d+0 ; td_cs_cm2(238) = 1.306700d-22
  td_E_eV(239) = 3.800000d+0 ; td_cs_cm2(239) = 1.457900d-22
  td_E_eV(240) = 3.810000d+0 ; td_cs_cm2(240) = 1.593900d-22
  td_E_eV(241) = 3.820000d+0 ; td_cs_cm2(241) = 1.706600d-22
  td_E_eV(242) = 3.830000d+0 ; td_cs_cm2(242) = 1.788400d-22
  td_E_eV(243) = 3.840000d+0 ; td_cs_cm2(243) = 1.832100d-22
  td_E_eV(244) = 3.850000d+0 ; td_cs_cm2(244) = 1.832300d-22
  td_E_eV(245) = 3.860000d+0 ; td_cs_cm2(245) = 1.785500d-22
  td_E_eV(246) = 3.870000d+0 ; td_cs_cm2(246) = 1.692100d-22
  td_E_eV(247) = 3.880000d+0 ; td_cs_cm2(247) = 1.556800d-22
  td_E_eV(248) = 3.890000d+0 ; td_cs_cm2(248) = 1.388900d-22
  td_E_eV(249) = 3.900000d+0 ; td_cs_cm2(249) = 1.201600d-22
  td_E_eV(250) = 3.910000d+0 ; td_cs_cm2(250) = 1.009900d-22
  td_E_eV(251) = 3.920000d+0 ; td_cs_cm2(251) = 8.279600d-23
  td_E_eV(252) = 3.930000d+0 ; td_cs_cm2(252) = 6.675900d-23
  td_E_eV(253) = 3.940000d+0 ; td_cs_cm2(253) = 5.366400d-23
  td_E_eV(254) = 3.950000d+0 ; td_cs_cm2(254) = 4.390400d-23
  td_E_eV(255) = 3.960000d+0 ; td_cs_cm2(255) = 3.753500d-23
  td_E_eV(256) = 3.970000d+0 ; td_cs_cm2(256) = 3.436800d-23
  td_E_eV(257) = 3.980000d+0 ; td_cs_cm2(257) = 3.406100d-23
  td_E_eV(258) = 3.990000d+0 ; td_cs_cm2(258) = 3.618700d-23
  td_E_eV(259) = 4.000000d+0 ; td_cs_cm2(259) = 4.026900d-23
  td_E_eV(260) = 4.010000d+0 ; td_cs_cm2(260) = 4.580500d-23
  td_E_eV(261) = 4.020000d+0 ; td_cs_cm2(261) = 5.226200d-23
  td_E_eV(262) = 4.030000d+0 ; td_cs_cm2(262) = 5.907600d-23
  td_E_eV(263) = 4.040000d+0 ; td_cs_cm2(263) = 6.565000d-23
  td_E_eV(264) = 4.050000d+0 ; td_cs_cm2(264) = 7.137300d-23
  td_E_eV(265) = 4.060000d+0 ; td_cs_cm2(265) = 7.565700d-23
  td_E_eV(266) = 4.070000d+0 ; td_cs_cm2(266) = 7.801100d-23
  td_E_eV(267) = 4.080000d+0 ; td_cs_cm2(267) = 7.811700d-23
  td_E_eV(268) = 4.090000d+0 ; td_cs_cm2(268) = 7.590200d-23
  td_E_eV(269) = 4.100000d+0 ; td_cs_cm2(269) = 7.156400d-23
  td_E_eV(270) = 4.110000d+0 ; td_cs_cm2(270) = 6.554100d-23
  td_E_eV(271) = 4.120000d+0 ; td_cs_cm2(271) = 5.842400d-23
  td_E_eV(272) = 4.130000d+0 ; td_cs_cm2(272) = 5.085300d-23
  td_E_eV(273) = 4.140000d+0 ; td_cs_cm2(273) = 4.341300d-23
  td_E_eV(274) = 4.150000d+0 ; td_cs_cm2(274) = 3.657900d-23
  td_E_eV(275) = 4.160000d+0 ; td_cs_cm2(275) = 3.068600d-23
  td_E_eV(276) = 4.170000d+0 ; td_cs_cm2(276) = 2.593800d-23
  td_E_eV(277) = 4.180000d+0 ; td_cs_cm2(277) = 2.243300d-23
  td_E_eV(278) = 4.190000d+0 ; td_cs_cm2(278) = 2.018300d-23
  td_E_eV(279) = 4.200000d+0 ; td_cs_cm2(279) = 1.913700d-23
  td_E_eV(280) = 4.210000d+0 ; td_cs_cm2(280) = 1.919400d-23
  td_E_eV(281) = 4.220000d+0 ; td_cs_cm2(281) = 2.020900d-23
  td_E_eV(282) = 4.230000d+0 ; td_cs_cm2(282) = 2.199200d-23
  td_E_eV(283) = 4.240000d+0 ; td_cs_cm2(283) = 2.430700d-23
  td_E_eV(284) = 4.250000d+0 ; td_cs_cm2(284) = 2.687500d-23
  td_E_eV(285) = 4.260000d+0 ; td_cs_cm2(285) = 2.938800d-23
  td_E_eV(286) = 4.270000d+0 ; td_cs_cm2(286) = 3.153500d-23
  td_E_eV(287) = 4.280000d+0 ; td_cs_cm2(287) = 3.304600d-23
  td_E_eV(288) = 4.290000d+0 ; td_cs_cm2(288) = 3.372700d-23
  td_E_eV(289) = 4.300000d+0 ; td_cs_cm2(289) = 3.349800d-23
  td_E_eV(290) = 4.310000d+0 ; td_cs_cm2(290) = 3.239400d-23
  td_E_eV(291) = 4.320000d+0 ; td_cs_cm2(291) = 3.054900d-23
  td_E_eV(292) = 4.330000d+0 ; td_cs_cm2(292) = 2.816400d-23
  td_E_eV(293) = 4.340000d+0 ; td_cs_cm2(293) = 2.546300d-23
  td_E_eV(294) = 4.350000d+0 ; td_cs_cm2(294) = 2.266100d-23
  td_E_eV(295) = 4.360000d+0 ; td_cs_cm2(295) = 1.994400d-23
  td_E_eV(296) = 4.370000d+0 ; td_cs_cm2(296) = 1.745900d-23
  td_E_eV(297) = 4.380000d+0 ; td_cs_cm2(297) = 1.531500d-23
  td_E_eV(298) = 4.390000d+0 ; td_cs_cm2(298) = 1.358600d-23
  td_E_eV(299) = 4.400000d+0 ; td_cs_cm2(299) = 1.231600d-23
  td_E_eV(300) = 4.410000d+0 ; td_cs_cm2(300) = 1.152000d-23
  td_E_eV(301) = 4.420000d+0 ; td_cs_cm2(301) = 1.118900d-23
  td_E_eV(302) = 4.430000d+0 ; td_cs_cm2(302) = 1.128200d-23
  td_E_eV(303) = 4.440000d+0 ; td_cs_cm2(303) = 1.173100d-23
  td_E_eV(304) = 4.450000d+0 ; td_cs_cm2(304) = 1.243600d-23
  td_E_eV(305) = 4.460000d+0 ; td_cs_cm2(305) = 1.327300d-23
  td_E_eV(306) = 4.470000d+0 ; td_cs_cm2(306) = 1.410600d-23
  td_E_eV(307) = 4.480000d+0 ; td_cs_cm2(307) = 1.480400d-23
  td_E_eV(308) = 4.490000d+0 ; td_cs_cm2(308) = 1.526000d-23
  td_E_eV(309) = 4.500000d+0 ; td_cs_cm2(309) = 1.540900d-23
  td_E_eV(310) = 4.510000d+0 ; td_cs_cm2(310) = 1.523200d-23
  td_E_eV(311) = 4.520000d+0 ; td_cs_cm2(311) = 1.475200d-23
  td_E_eV(312) = 4.530000d+0 ; td_cs_cm2(312) = 1.402100d-23
  td_E_eV(313) = 4.540000d+0 ; td_cs_cm2(313) = 1.311100d-23
  td_E_eV(314) = 4.550000d+0 ; td_cs_cm2(314) = 1.209400d-23
  td_E_eV(315) = 4.560000d+0 ; td_cs_cm2(315) = 1.104400d-23
  td_E_eV(316) = 4.570000d+0 ; td_cs_cm2(316) = 1.002100d-23
  td_E_eV(317) = 4.580000d+0 ; td_cs_cm2(317) = 9.078700d-24
  td_E_eV(318) = 4.590000d+0 ; td_cs_cm2(318) = 8.258600d-24
  td_E_eV(319) = 4.600000d+0 ; td_cs_cm2(319) = 7.592000d-24
  td_E_eV(320) = 4.610000d+0 ; td_cs_cm2(320) = 7.099700d-24
  td_E_eV(321) = 4.620000d+0 ; td_cs_cm2(321) = 6.790000d-24
  td_E_eV(322) = 4.630000d+0 ; td_cs_cm2(322) = 6.657600d-24
  td_E_eV(323) = 4.640000d+0 ; td_cs_cm2(323) = 6.681800d-24
  td_E_eV(324) = 4.650000d+0 ; td_cs_cm2(324) = 6.826700d-24
  td_E_eV(325) = 4.660000d+0 ; td_cs_cm2(325) = 7.044500d-24
  td_E_eV(326) = 4.670000d+0 ; td_cs_cm2(326) = 7.280900d-24
  td_E_eV(327) = 4.680000d+0 ; td_cs_cm2(327) = 7.484000d-24
  td_E_eV(328) = 4.690000d+0 ; td_cs_cm2(328) = 7.611800d-24
  td_E_eV(329) = 4.700000d+0 ; td_cs_cm2(329) = 7.637500d-24
  td_E_eV(330) = 4.710000d+0 ; td_cs_cm2(330) = 7.551200d-24
  td_E_eV(331) = 4.720000d+0 ; td_cs_cm2(331) = 7.358000d-24
  td_E_eV(332) = 4.730000d+0 ; td_cs_cm2(332) = 7.073500d-24
  td_E_eV(333) = 4.740000d+0 ; td_cs_cm2(333) = 6.719900d-24
  td_E_eV(334) = 4.750000d+0 ; td_cs_cm2(334) = 6.321800d-24
  td_E_eV(335) = 4.760000d+0 ; td_cs_cm2(335) = 5.903900d-24
  td_E_eV(336) = 4.770000d+0 ; td_cs_cm2(336) = 5.489500d-24
  td_E_eV(337) = 4.780000d+0 ; td_cs_cm2(337) = 5.099300d-24
  td_E_eV(338) = 4.790000d+0 ; td_cs_cm2(338) = 4.751400d-24
  td_E_eV(339) = 4.800000d+0 ; td_cs_cm2(339) = 4.460300d-24
  td_E_eV(340) = 4.810000d+0 ; td_cs_cm2(340) = 4.235900d-24
  td_E_eV(341) = 4.820000d+0 ; td_cs_cm2(341) = 4.082700d-24
  td_E_eV(342) = 4.830000d+0 ; td_cs_cm2(342) = 3.998500d-24
  td_E_eV(343) = 4.840000d+0 ; td_cs_cm2(343) = 3.973900d-24
  td_E_eV(344) = 4.850000d+0 ; td_cs_cm2(344) = 3.993200d-24
  td_E_eV(345) = 4.860000d+0 ; td_cs_cm2(345) = 4.036500d-24
  td_E_eV(346) = 4.870000d+0 ; td_cs_cm2(346) = 4.082600d-24
  td_E_eV(347) = 4.880000d+0 ; td_cs_cm2(347) = 4.113100d-24
  td_E_eV(348) = 4.890000d+0 ; td_cs_cm2(348) = 4.114100d-24
  td_E_eV(349) = 4.900000d+0 ; td_cs_cm2(349) = 4.078000d-24
  td_E_eV(350) = 4.910000d+0 ; td_cs_cm2(350) = 4.002900d-24
  td_E_eV(351) = 4.920000d+0 ; td_cs_cm2(351) = 3.891500d-24
  td_E_eV(352) = 4.930000d+0 ; td_cs_cm2(352) = 3.749900d-24
  td_E_eV(353) = 4.940000d+0 ; td_cs_cm2(353) = 3.586000d-24
  td_E_eV(354) = 4.950000d+0 ; td_cs_cm2(354) = 3.408400d-24
  td_E_eV(355) = 4.960000d+0 ; td_cs_cm2(355) = 3.226200d-24
  td_E_eV(356) = 4.970000d+0 ; td_cs_cm2(356) = 3.048100d-24
  td_E_eV(357) = 4.980000d+0 ; td_cs_cm2(357) = 2.882200d-24
  td_E_eV(358) = 4.990000d+0 ; td_cs_cm2(358) = 2.735700d-24
  td_E_eV(359) = 5.000000d+0 ; td_cs_cm2(359) = 2.613900d-24
  td_E_eV(360) = 5.010000d+0 ; td_cs_cm2(360) = 2.520200d-24
  td_E_eV(361) = 5.020000d+0 ; td_cs_cm2(361) = 2.454900d-24
  td_E_eV(362) = 5.030000d+0 ; td_cs_cm2(362) = 2.415300d-24
  td_E_eV(363) = 5.040000d+0 ; td_cs_cm2(363) = 2.395600d-24
  td_E_eV(364) = 5.050000d+0 ; td_cs_cm2(364) = 2.388300d-24
  td_E_eV(365) = 5.060000d+0 ; td_cs_cm2(365) = 2.385200d-24
  td_E_eV(366) = 5.070000d+0 ; td_cs_cm2(366) = 2.378700d-24
  td_E_eV(367) = 5.080000d+0 ; td_cs_cm2(367) = 2.363100d-24
  td_E_eV(368) = 5.090000d+0 ; td_cs_cm2(368) = 2.334800d-24
  td_E_eV(369) = 5.100000d+0 ; td_cs_cm2(369) = 2.292500d-24
  td_E_eV(370) = 5.110000d+0 ; td_cs_cm2(370) = 2.236500d-24
  td_E_eV(371) = 5.120000d+0 ; td_cs_cm2(371) = 2.168600d-24
  td_E_eV(372) = 5.130000d+0 ; td_cs_cm2(372) = 2.091200d-24
  td_E_eV(373) = 5.140000d+0 ; td_cs_cm2(373) = 2.007700d-24
  td_E_eV(374) = 5.150000d+0 ; td_cs_cm2(374) = 1.921600d-24
  td_E_eV(375) = 5.160000d+0 ; td_cs_cm2(375) = 1.836300d-24
  td_E_eV(376) = 5.170000d+0 ; td_cs_cm2(376) = 1.755500d-24
  td_E_eV(377) = 5.180000d+0 ; td_cs_cm2(377) = 1.682400d-24
  td_E_eV(378) = 5.190000d+0 ; td_cs_cm2(378) = 1.619500d-24
  td_E_eV(379) = 5.200000d+0 ; td_cs_cm2(379) = 1.568300d-24
  td_E_eV(380) = 5.210000d+0 ; td_cs_cm2(380) = 1.529000d-24
  td_E_eV(381) = 5.220000d+0 ; td_cs_cm2(381) = 1.500500d-24
  td_E_eV(382) = 5.230000d+0 ; td_cs_cm2(382) = 1.480200d-24
  td_E_eV(383) = 5.240000d+0 ; td_cs_cm2(383) = 1.465200d-24
  td_E_eV(384) = 5.250000d+0 ; td_cs_cm2(384) = 1.452200d-24
  td_E_eV(385) = 5.260000d+0 ; td_cs_cm2(385) = 1.438300d-24
  td_E_eV(386) = 5.270000d+0 ; td_cs_cm2(386) = 1.421200d-24
  td_E_eV(387) = 5.280000d+0 ; td_cs_cm2(387) = 1.399500d-24
  td_E_eV(388) = 5.290000d+0 ; td_cs_cm2(388) = 1.372700d-24
  td_E_eV(389) = 5.300000d+0 ; td_cs_cm2(389) = 1.340700d-24
  td_E_eV(390) = 5.310000d+0 ; td_cs_cm2(390) = 1.304000d-24
  td_E_eV(391) = 5.320000d+0 ; td_cs_cm2(391) = 1.263700d-24
  td_E_eV(392) = 5.330000d+0 ; td_cs_cm2(392) = 1.221000d-24
  td_E_eV(393) = 5.340000d+0 ; td_cs_cm2(393) = 1.177200d-24
  td_E_eV(394) = 5.350000d+0 ; td_cs_cm2(394) = 1.134000d-24
  td_E_eV(395) = 5.360000d+0 ; td_cs_cm2(395) = 1.092800d-24
  td_E_eV(396) = 5.370000d+0 ; td_cs_cm2(396) = 1.055100d-24
  td_E_eV(397) = 5.380000d+0 ; td_cs_cm2(397) = 1.022000d-24
  td_E_eV(398) = 5.390000d+0 ; td_cs_cm2(398) = 9.940600d-25
  td_E_eV(399) = 5.400000d+0 ; td_cs_cm2(399) = 9.711900d-25
  td_E_eV(400) = 5.410000d+0 ; td_cs_cm2(400) = 9.527800d-25
  td_E_eV(401) = 5.420000d+0 ; td_cs_cm2(401) = 9.377400d-25
  td_E_eV(402) = 5.430000d+0 ; td_cs_cm2(402) = 9.248000d-25
  td_E_eV(403) = 5.440000d+0 ; td_cs_cm2(403) = 9.126500d-25
  td_E_eV(404) = 5.450000d+0 ; td_cs_cm2(404) = 9.001800d-25
  td_E_eV(405) = 5.460000d+0 ; td_cs_cm2(405) = 8.865500d-25
  td_E_eV(406) = 5.470000d+0 ; td_cs_cm2(406) = 8.712100d-25
  td_E_eV(407) = 5.480000d+0 ; td_cs_cm2(407) = 8.539200d-25
  td_E_eV(408) = 5.490000d+0 ; td_cs_cm2(408) = 8.346500d-25
  td_E_eV(409) = 5.500000d+0 ; td_cs_cm2(409) = 8.136000d-25
  td_E_eV(410) = 5.510000d+0 ; td_cs_cm2(410) = 7.911400d-25
  td_E_eV(411) = 5.520000d+0 ; td_cs_cm2(411) = 7.677900d-25
  td_E_eV(412) = 5.530000d+0 ; td_cs_cm2(412) = 7.442000d-25
  td_E_eV(413) = 5.540000d+0 ; td_cs_cm2(413) = 7.210800d-25
  td_E_eV(414) = 5.550000d+0 ; td_cs_cm2(414) = 6.991400d-25
  td_E_eV(415) = 5.560000d+0 ; td_cs_cm2(415) = 6.789800d-25
  td_E_eV(416) = 5.570000d+0 ; td_cs_cm2(416) = 6.610000d-25
  td_E_eV(417) = 5.580000d+0 ; td_cs_cm2(417) = 6.453400d-25
  td_E_eV(418) = 5.590000d+0 ; td_cs_cm2(418) = 6.318600d-25
  td_E_eV(419) = 5.600000d+0 ; td_cs_cm2(419) = 6.202200d-25
  td_E_eV(420) = 5.610000d+0 ; td_cs_cm2(420) = 6.099100d-25
  td_E_eV(421) = 5.620000d+0 ; td_cs_cm2(421) = 6.004100d-25
  td_E_eV(422) = 5.630000d+0 ; td_cs_cm2(422) = 5.911900d-25
  td_E_eV(423) = 5.640000d+0 ; td_cs_cm2(423) = 5.818400d-25
  td_E_eV(424) = 5.650000d+0 ; td_cs_cm2(424) = 5.720500d-25
  td_E_eV(425) = 5.660000d+0 ; td_cs_cm2(425) = 5.615800d-25
  td_E_eV(426) = 5.670000d+0 ; td_cs_cm2(426) = 5.503500d-25
  td_E_eV(427) = 5.680000d+0 ; td_cs_cm2(427) = 5.383400d-25
  td_E_eV(428) = 5.690000d+0 ; td_cs_cm2(428) = 5.256400d-25
  td_E_eV(429) = 5.700000d+0 ; td_cs_cm2(429) = 5.124500d-25
  td_E_eV(430) = 5.710000d+0 ; td_cs_cm2(430) = 4.990200d-25
  td_E_eV(431) = 5.720000d+0 ; td_cs_cm2(431) = 4.856900d-25
  td_E_eV(432) = 5.730000d+0 ; td_cs_cm2(432) = 4.727900d-25
  td_E_eV(433) = 5.740000d+0 ; td_cs_cm2(433) = 4.606300d-25
  td_E_eV(434) = 5.750000d+0 ; td_cs_cm2(434) = 4.494500d-25
  td_E_eV(435) = 5.760000d+0 ; td_cs_cm2(435) = 4.393300d-25
  td_E_eV(436) = 5.770000d+0 ; td_cs_cm2(436) = 4.302600d-25
  td_E_eV(437) = 5.780000d+0 ; td_cs_cm2(437) = 4.221200d-25
  td_E_eV(438) = 5.790000d+0 ; td_cs_cm2(438) = 4.147200d-25
  td_E_eV(439) = 5.800000d+0 ; td_cs_cm2(439) = 4.078200d-25
  td_E_eV(440) = 5.810000d+0 ; td_cs_cm2(440) = 4.012100d-25
  td_E_eV(441) = 5.820000d+0 ; td_cs_cm2(441) = 3.946900d-25
  td_E_eV(442) = 5.830000d+0 ; td_cs_cm2(442) = 3.881000d-25
  td_E_eV(443) = 5.840000d+0 ; td_cs_cm2(443) = 3.812900d-25
  td_E_eV(444) = 5.850000d+0 ; td_cs_cm2(444) = 3.742000d-25
  td_E_eV(445) = 5.860000d+0 ; td_cs_cm2(445) = 3.667900d-25
  td_E_eV(446) = 5.870000d+0 ; td_cs_cm2(446) = 3.590800d-25
  td_E_eV(447) = 5.880000d+0 ; td_cs_cm2(447) = 3.511300d-25
  td_E_eV(448) = 5.890000d+0 ; td_cs_cm2(448) = 3.430800d-25
  td_E_eV(449) = 5.900000d+0 ; td_cs_cm2(449) = 3.350700d-25
  td_E_eV(450) = 5.910000d+0 ; td_cs_cm2(450) = 3.272800d-25
  td_E_eV(451) = 5.920000d+0 ; td_cs_cm2(451) = 3.198400d-25
  td_E_eV(452) = 5.930000d+0 ; td_cs_cm2(452) = 3.128700d-25
  td_E_eV(453) = 5.940000d+0 ; td_cs_cm2(453) = 3.064200d-25
  td_E_eV(454) = 5.950000d+0 ; td_cs_cm2(454) = 3.004700d-25
  td_E_eV(455) = 5.960000d+0 ; td_cs_cm2(455) = 2.949800d-25
  td_E_eV(456) = 5.970000d+0 ; td_cs_cm2(456) = 2.898700d-25
  td_E_eV(457) = 5.980000d+0 ; td_cs_cm2(457) = 2.850500d-25
  td_E_eV(458) = 5.990000d+0 ; td_cs_cm2(458) = 2.804100d-25
  td_E_eV(459) = 6.000000d+0 ; td_cs_cm2(459) = 2.758800d-25
  td_E_eV(460) = 6.010000d+0 ; td_cs_cm2(460) = 2.713500d-25
  td_E_eV(461) = 6.020000d+0 ; td_cs_cm2(461) = 2.667600d-25
  td_E_eV(462) = 6.030000d+0 ; td_cs_cm2(462) = 2.620800d-25
  td_E_eV(463) = 6.040000d+0 ; td_cs_cm2(463) = 2.572600d-25
  td_E_eV(464) = 6.050000d+0 ; td_cs_cm2(464) = 2.523100d-25
  td_E_eV(465) = 6.060000d+0 ; td_cs_cm2(465) = 2.472800d-25
  td_E_eV(466) = 6.070000d+0 ; td_cs_cm2(466) = 2.422300d-25
  td_E_eV(467) = 6.080000d+0 ; td_cs_cm2(467) = 2.372200d-25
  td_E_eV(468) = 6.090000d+0 ; td_cs_cm2(468) = 2.323500d-25
  td_E_eV(469) = 6.100000d+0 ; td_cs_cm2(469) = 2.276700d-25
  td_E_eV(470) = 6.110000d+0 ; td_cs_cm2(470) = 2.232300d-25
  td_E_eV(471) = 6.120000d+0 ; td_cs_cm2(471) = 2.190400d-25
  td_E_eV(472) = 6.130000d+0 ; td_cs_cm2(472) = 2.151100d-25
  td_E_eV(473) = 6.140000d+0 ; td_cs_cm2(473) = 2.114200d-25
  td_E_eV(474) = 6.150000d+0 ; td_cs_cm2(474) = 2.079200d-25
  td_E_eV(475) = 6.160000d+0 ; td_cs_cm2(475) = 2.045800d-25
  td_E_eV(476) = 6.170000d+0 ; td_cs_cm2(476) = 2.013500d-25
  td_E_eV(477) = 6.180000d+0 ; td_cs_cm2(477) = 1.981900d-25
  td_E_eV(478) = 6.190000d+0 ; td_cs_cm2(478) = 1.950500d-25
  td_E_eV(479) = 6.200000d+0 ; td_cs_cm2(479) = 1.918900d-25
  td_E_eV(480) = 6.210000d+0 ; td_cs_cm2(480) = 1.886900d-25
  td_E_eV(481) = 6.220000d+0 ; td_cs_cm2(481) = 1.854400d-25
  td_E_eV(482) = 6.230000d+0 ; td_cs_cm2(482) = 1.821500d-25
  td_E_eV(483) = 6.240000d+0 ; td_cs_cm2(483) = 1.788300d-25
  td_E_eV(484) = 6.250000d+0 ; td_cs_cm2(484) = 1.755400d-25
  td_E_eV(485) = 6.260000d+0 ; td_cs_cm2(485) = 1.722900d-25
  td_E_eV(486) = 6.270000d+0 ; td_cs_cm2(486) = 1.691300d-25
  td_E_eV(487) = 6.280000d+0 ; td_cs_cm2(487) = 1.660900d-25
  td_E_eV(488) = 6.290000d+0 ; td_cs_cm2(488) = 1.631800d-25
  td_E_eV(489) = 6.300000d+0 ; td_cs_cm2(489) = 1.604100d-25
  td_E_eV(490) = 6.310000d+0 ; td_cs_cm2(490) = 1.577800d-25
  td_E_eV(491) = 6.320000d+0 ; td_cs_cm2(491) = 1.552600d-25
  td_E_eV(492) = 6.330000d+0 ; td_cs_cm2(492) = 1.528500d-25
  td_E_eV(493) = 6.340000d+0 ; td_cs_cm2(493) = 1.505300d-25
  td_E_eV(494) = 6.350000d+0 ; td_cs_cm2(494) = 1.482600d-25
  td_E_eV(495) = 6.360000d+0 ; td_cs_cm2(495) = 1.460300d-25
  td_E_eV(496) = 6.370000d+0 ; td_cs_cm2(496) = 1.438100d-25
  td_E_eV(497) = 6.380000d+0 ; td_cs_cm2(497) = 1.415900d-25
  td_E_eV(498) = 6.390000d+0 ; td_cs_cm2(498) = 1.393500d-25
  td_E_eV(499) = 6.400000d+0 ; td_cs_cm2(499) = 1.371000d-25
  td_E_eV(500) = 6.410000d+0 ; td_cs_cm2(500) = 1.348400d-25
  td_E_eV(501) = 6.420000d+0 ; td_cs_cm2(501) = 1.325900d-25
  td_E_eV(502) = 6.430000d+0 ; td_cs_cm2(502) = 1.303700d-25
  td_E_eV(503) = 6.440000d+0 ; td_cs_cm2(503) = 1.282000d-25
  td_E_eV(504) = 6.450000d+0 ; td_cs_cm2(504) = 1.260900d-25
  td_E_eV(505) = 6.460000d+0 ; td_cs_cm2(505) = 1.240500d-25
  td_E_eV(506) = 6.470000d+0 ; td_cs_cm2(506) = 1.220900d-25
  td_E_eV(507) = 6.480000d+0 ; td_cs_cm2(507) = 1.202000d-25
  td_E_eV(508) = 6.490000d+0 ; td_cs_cm2(508) = 1.183900d-25
  td_E_eV(509) = 6.500000d+0 ; td_cs_cm2(509) = 1.166400d-25
  td_E_eV(510) = 6.510000d+0 ; td_cs_cm2(510) = 1.149600d-25
  td_E_eV(511) = 6.520000d+0 ; td_cs_cm2(511) = 1.133100d-25
  td_E_eV(512) = 6.530000d+0 ; td_cs_cm2(512) = 1.117000d-25
  td_E_eV(513) = 6.540000d+0 ; td_cs_cm2(513) = 1.101000d-25
  td_E_eV(514) = 6.550000d+0 ; td_cs_cm2(514) = 1.085100d-25
  td_E_eV(515) = 6.560000d+0 ; td_cs_cm2(515) = 1.069100d-25
  td_E_eV(516) = 6.570000d+0 ; td_cs_cm2(516) = 1.053200d-25
  td_E_eV(517) = 6.580000d+0 ; td_cs_cm2(517) = 1.037200d-25
  td_E_eV(518) = 6.590000d+0 ; td_cs_cm2(518) = 1.021400d-25
  td_E_eV(519) = 6.600000d+0 ; td_cs_cm2(519) = 1.005800d-25
  td_E_eV(520) = 6.610000d+0 ; td_cs_cm2(520) = 9.904200d-26
  td_E_eV(521) = 6.620000d+0 ; td_cs_cm2(521) = 9.754300d-26
  td_E_eV(522) = 6.630000d+0 ; td_cs_cm2(522) = 9.608500d-26
  td_E_eV(523) = 6.640000d+0 ; td_cs_cm2(523) = 9.467400d-26
  td_E_eV(524) = 6.650000d+0 ; td_cs_cm2(524) = 9.331000d-26
  td_E_eV(525) = 6.660000d+0 ; td_cs_cm2(525) = 9.199300d-26
  td_E_eV(526) = 6.670000d+0 ; td_cs_cm2(526) = 9.072000d-26
  td_E_eV(527) = 6.680000d+0 ; td_cs_cm2(527) = 8.948300d-26
  td_E_eV(528) = 6.690000d+0 ; td_cs_cm2(528) = 8.827500d-26
  td_E_eV(529) = 6.700000d+0 ; td_cs_cm2(529) = 8.708700d-26
  td_E_eV(530) = 6.710000d+0 ; td_cs_cm2(530) = 8.591100d-26
  td_E_eV(531) = 6.720000d+0 ; td_cs_cm2(531) = 8.474200d-26
  td_E_eV(532) = 6.730000d+0 ; td_cs_cm2(532) = 8.357600d-26
  td_E_eV(533) = 6.740000d+0 ; td_cs_cm2(533) = 8.241400d-26
  td_E_eV(534) = 6.750000d+0 ; td_cs_cm2(534) = 8.125700d-26
  td_E_eV(535) = 6.760000d+0 ; td_cs_cm2(535) = 8.011000d-26
  td_E_eV(536) = 6.770000d+0 ; td_cs_cm2(536) = 7.897700d-26
  td_E_eV(537) = 6.780000d+0 ; td_cs_cm2(537) = 7.786400d-26
  td_E_eV(538) = 6.790000d+0 ; td_cs_cm2(538) = 7.677400d-26
  td_E_eV(539) = 6.800000d+0 ; td_cs_cm2(539) = 7.571300d-26
  td_E_eV(540) = 6.810000d+0 ; td_cs_cm2(540) = 7.468100d-26
  td_E_eV(541) = 6.820000d+0 ; td_cs_cm2(541) = 7.368200d-26
  td_E_eV(542) = 6.830000d+0 ; td_cs_cm2(542) = 7.271300d-26
  td_E_eV(543) = 6.840000d+0 ; td_cs_cm2(543) = 7.177200d-26
  td_E_eV(544) = 6.850000d+0 ; td_cs_cm2(544) = 7.085300d-26
  td_E_eV(545) = 6.860000d+0 ; td_cs_cm2(545) = 6.995200d-26
  td_E_eV(546) = 6.870000d+0 ; td_cs_cm2(546) = 6.906300d-26
  td_E_eV(547) = 6.880000d+0 ; td_cs_cm2(547) = 6.818200d-26
  td_E_eV(548) = 6.890000d+0 ; td_cs_cm2(548) = 6.730600d-26
  td_E_eV(549) = 6.900000d+0 ; td_cs_cm2(549) = 6.643500d-26
  td_E_eV(550) = 6.910000d+0 ; td_cs_cm2(550) = 6.556700d-26
  td_E_eV(551) = 6.920000d+0 ; td_cs_cm2(551) = 6.470600d-26
  td_E_eV(552) = 6.930000d+0 ; td_cs_cm2(552) = 6.385400d-26
  td_E_eV(553) = 6.940000d+0 ; td_cs_cm2(553) = 6.301400d-26
  td_E_eV(554) = 6.950000d+0 ; td_cs_cm2(554) = 6.218900d-26
  td_E_eV(555) = 6.960000d+0 ; td_cs_cm2(555) = 6.138200d-26
  td_E_eV(556) = 6.970000d+0 ; td_cs_cm2(556) = 6.059600d-26
  td_E_eV(557) = 6.980000d+0 ; td_cs_cm2(557) = 5.983100d-26
  td_E_eV(558) = 6.990000d+0 ; td_cs_cm2(558) = 5.908800d-26
  td_E_eV(559) = 7.000000d+0 ; td_cs_cm2(559) = 5.836400d-26
  td_E_eV(560) = 7.010000d+0 ; td_cs_cm2(560) = 5.765700d-26
  td_E_eV(561) = 7.020000d+0 ; td_cs_cm2(561) = 5.696300d-26
  td_E_eV(562) = 7.030000d+0 ; td_cs_cm2(562) = 5.627900d-26
  td_E_eV(563) = 7.040000d+0 ; td_cs_cm2(563) = 5.560100d-26
  td_E_eV(564) = 7.050000d+0 ; td_cs_cm2(564) = 5.492800d-26
  td_E_eV(565) = 7.060000d+0 ; td_cs_cm2(565) = 5.425900d-26
  td_E_eV(566) = 7.070000d+0 ; td_cs_cm2(566) = 5.359400d-26
  td_E_eV(567) = 7.080000d+0 ; td_cs_cm2(567) = 5.293400d-26
  td_E_eV(568) = 7.090000d+0 ; td_cs_cm2(568) = 5.228100d-26
  td_E_eV(569) = 7.100000d+0 ; td_cs_cm2(569) = 5.163600d-26
  td_E_eV(570) = 7.110000d+0 ; td_cs_cm2(570) = 5.100300d-26
  td_E_eV(571) = 7.120000d+0 ; td_cs_cm2(571) = 5.038200d-26
  td_E_eV(572) = 7.130000d+0 ; td_cs_cm2(572) = 4.977600d-26
  td_E_eV(573) = 7.140000d+0 ; td_cs_cm2(573) = 4.918600d-26
  td_E_eV(574) = 7.150000d+0 ; td_cs_cm2(574) = 4.861000d-26
  td_E_eV(575) = 7.160000d+0 ; td_cs_cm2(575) = 4.804800d-26
  td_E_eV(576) = 7.170000d+0 ; td_cs_cm2(576) = 4.749600d-26
  td_E_eV(577) = 7.180000d+0 ; td_cs_cm2(577) = 4.695400d-26
  td_E_eV(578) = 7.190000d+0 ; td_cs_cm2(578) = 4.641800d-26
  td_E_eV(579) = 7.200000d+0 ; td_cs_cm2(579) = 4.588800d-26
  td_E_eV(580) = 7.210000d+0 ; td_cs_cm2(580) = 4.536200d-26
  td_E_eV(581) = 7.220000d+0 ; td_cs_cm2(581) = 4.483900d-26
  td_E_eV(582) = 7.230000d+0 ; td_cs_cm2(582) = 4.432000d-26
  td_E_eV(583) = 7.240000d+0 ; td_cs_cm2(583) = 4.380500d-26
  td_E_eV(584) = 7.250000d+0 ; td_cs_cm2(584) = 4.329700d-26
  td_E_eV(585) = 7.260000d+0 ; td_cs_cm2(585) = 4.279500d-26
  td_E_eV(586) = 7.270000d+0 ; td_cs_cm2(586) = 4.230200d-26
  td_E_eV(587) = 7.280000d+0 ; td_cs_cm2(587) = 4.181900d-26
  td_E_eV(588) = 7.290000d+0 ; td_cs_cm2(588) = 4.134700d-26
  td_E_eV(589) = 7.300000d+0 ; td_cs_cm2(589) = 4.088600d-26
  td_E_eV(590) = 7.310000d+0 ; td_cs_cm2(590) = 4.043400d-26
  td_E_eV(591) = 7.320000d+0 ; td_cs_cm2(591) = 3.999100d-26
  td_E_eV(592) = 7.330000d+0 ; td_cs_cm2(592) = 3.955600d-26
  td_E_eV(593) = 7.340000d+0 ; td_cs_cm2(593) = 3.912600d-26
  td_E_eV(594) = 7.350000d+0 ; td_cs_cm2(594) = 3.870100d-26
  td_E_eV(595) = 7.360000d+0 ; td_cs_cm2(595) = 3.828000d-26
  td_E_eV(596) = 7.370000d+0 ; td_cs_cm2(596) = 3.786300d-26
  td_E_eV(597) = 7.380000d+0 ; td_cs_cm2(597) = 3.744800d-26
  td_E_eV(598) = 7.390000d+0 ; td_cs_cm2(598) = 3.703700d-26
  td_E_eV(599) = 7.400000d+0 ; td_cs_cm2(599) = 3.663100d-26
  td_E_eV(600) = 7.410000d+0 ; td_cs_cm2(600) = 3.623000d-26
  td_E_eV(601) = 7.420000d+0 ; td_cs_cm2(601) = 3.583500d-26
  td_E_eV(602) = 7.430000d+0 ; td_cs_cm2(602) = 3.544700d-26
  td_E_eV(603) = 7.440000d+0 ; td_cs_cm2(603) = 3.506700d-26
  td_E_eV(604) = 7.450000d+0 ; td_cs_cm2(604) = 3.469400d-26
  td_E_eV(605) = 7.460000d+0 ; td_cs_cm2(605) = 3.432800d-26
  td_E_eV(606) = 7.470000d+0 ; td_cs_cm2(606) = 3.396900d-26
  td_E_eV(607) = 7.480000d+0 ; td_cs_cm2(607) = 3.361600d-26
  td_E_eV(608) = 7.490000d+0 ; td_cs_cm2(608) = 3.326700d-26
  td_E_eV(609) = 7.500000d+0 ; td_cs_cm2(609) = 3.292300d-26
  td_E_eV(610) = 7.510000d+0 ; td_cs_cm2(610) = 3.258200d-26
  td_E_eV(611) = 7.520000d+0 ; td_cs_cm2(611) = 3.224400d-26
  td_E_eV(612) = 7.530000d+0 ; td_cs_cm2(612) = 3.190900d-26
  td_E_eV(613) = 7.540000d+0 ; td_cs_cm2(613) = 3.157700d-26
  td_E_eV(614) = 7.550000d+0 ; td_cs_cm2(614) = 3.124800d-26
  td_E_eV(615) = 7.560000d+0 ; td_cs_cm2(615) = 3.092400d-26
  td_E_eV(616) = 7.570000d+0 ; td_cs_cm2(616) = 3.060400d-26
  td_E_eV(617) = 7.580000d+0 ; td_cs_cm2(617) = 3.029000d-26
  td_E_eV(618) = 7.590000d+0 ; td_cs_cm2(618) = 2.998100d-26
  td_E_eV(619) = 7.600000d+0 ; td_cs_cm2(619) = 2.967800d-26
  td_E_eV(620) = 7.610000d+0 ; td_cs_cm2(620) = 2.937900d-26
  td_E_eV(621) = 7.620000d+0 ; td_cs_cm2(621) = 2.908600d-26
  td_E_eV(622) = 7.630000d+0 ; td_cs_cm2(622) = 2.879700d-26
  td_E_eV(623) = 7.640000d+0 ; td_cs_cm2(623) = 2.851100d-26
  td_E_eV(624) = 7.650000d+0 ; td_cs_cm2(624) = 2.822900d-26
  td_E_eV(625) = 7.660000d+0 ; td_cs_cm2(625) = 2.794900d-26
  td_E_eV(626) = 7.670000d+0 ; td_cs_cm2(626) = 2.767300d-26
  td_E_eV(627) = 7.680000d+0 ; td_cs_cm2(627) = 2.739900d-26
  td_E_eV(628) = 7.690000d+0 ; td_cs_cm2(628) = 2.712700d-26
  td_E_eV(629) = 7.700000d+0 ; td_cs_cm2(629) = 2.685900d-26
  td_E_eV(630) = 7.710000d+0 ; td_cs_cm2(630) = 2.659400d-26
  td_E_eV(631) = 7.720000d+0 ; td_cs_cm2(631) = 2.633300d-26
  td_E_eV(632) = 7.730000d+0 ; td_cs_cm2(632) = 2.607600d-26
  td_E_eV(633) = 7.740000d+0 ; td_cs_cm2(633) = 2.582300d-26
  td_E_eV(634) = 7.750000d+0 ; td_cs_cm2(634) = 2.557400d-26
  td_E_eV(635) = 7.760000d+0 ; td_cs_cm2(635) = 2.532900d-26
  td_E_eV(636) = 7.770000d+0 ; td_cs_cm2(636) = 2.508700d-26
  td_E_eV(637) = 7.780000d+0 ; td_cs_cm2(637) = 2.484800d-26
  td_E_eV(638) = 7.790000d+0 ; td_cs_cm2(638) = 2.461200d-26
  td_E_eV(639) = 7.800000d+0 ; td_cs_cm2(639) = 2.437900d-26
  td_E_eV(640) = 7.810000d+0 ; td_cs_cm2(640) = 2.414800d-26
  td_E_eV(641) = 7.820000d+0 ; td_cs_cm2(641) = 2.392000d-26
  td_E_eV(642) = 7.830000d+0 ; td_cs_cm2(642) = 2.369400d-26
  td_E_eV(643) = 7.840000d+0 ; td_cs_cm2(643) = 2.347000d-26
  td_E_eV(644) = 7.850000d+0 ; td_cs_cm2(644) = 2.324900d-26
  td_E_eV(645) = 7.860000d+0 ; td_cs_cm2(645) = 2.303100d-26
  td_E_eV(646) = 7.870000d+0 ; td_cs_cm2(646) = 2.281600d-26
  td_E_eV(647) = 7.880000d+0 ; td_cs_cm2(647) = 2.260400d-26
  td_E_eV(648) = 7.890000d+0 ; td_cs_cm2(648) = 2.239400d-26
  td_E_eV(649) = 7.900000d+0 ; td_cs_cm2(649) = 2.218800d-26
  td_E_eV(650) = 7.910000d+0 ; td_cs_cm2(650) = 2.198400d-26
  td_E_eV(651) = 7.920000d+0 ; td_cs_cm2(651) = 2.178300d-26
  td_E_eV(652) = 7.930000d+0 ; td_cs_cm2(652) = 2.158500d-26
  td_E_eV(653) = 7.940000d+0 ; td_cs_cm2(653) = 2.138800d-26
  td_E_eV(654) = 7.950000d+0 ; td_cs_cm2(654) = 2.119400d-26
  td_E_eV(655) = 7.960000d+0 ; td_cs_cm2(655) = 2.100200d-26
  td_E_eV(656) = 7.970000d+0 ; td_cs_cm2(656) = 2.081200d-26
  td_E_eV(657) = 7.980000d+0 ; td_cs_cm2(657) = 2.062400d-26
  td_E_eV(658) = 7.990000d+0 ; td_cs_cm2(658) = 2.043900d-26
  td_E_eV(659) = 8.000000d+0 ; td_cs_cm2(659) = 2.025500d-26
  td_E_eV(660) = 8.010000d+0 ; td_cs_cm2(660) = 2.007400d-26
  td_E_eV(661) = 8.020000d+0 ; td_cs_cm2(661) = 1.989500d-26
  td_E_eV(662) = 8.030000d+0 ; td_cs_cm2(662) = 1.971800d-26
  td_E_eV(663) = 8.040000d+0 ; td_cs_cm2(663) = 1.954300d-26
  td_E_eV(664) = 8.050000d+0 ; td_cs_cm2(664) = 1.937100d-26
  td_E_eV(665) = 8.060000d+0 ; td_cs_cm2(665) = 1.920100d-26
  td_E_eV(666) = 8.070000d+0 ; td_cs_cm2(666) = 1.903200d-26
  td_E_eV(667) = 8.080000d+0 ; td_cs_cm2(667) = 1.886600d-26
  td_E_eV(668) = 8.090000d+0 ; td_cs_cm2(668) = 1.870200d-26
  td_E_eV(669) = 8.100000d+0 ; td_cs_cm2(669) = 1.853900d-26
  td_E_eV(670) = 8.110000d+0 ; td_cs_cm2(670) = 1.837800d-26
  td_E_eV(671) = 8.120000d+0 ; td_cs_cm2(671) = 1.821900d-26
  td_E_eV(672) = 8.130000d+0 ; td_cs_cm2(672) = 1.806200d-26
  td_E_eV(673) = 8.140000d+0 ; td_cs_cm2(673) = 1.790700d-26
  td_E_eV(674) = 8.150000d+0 ; td_cs_cm2(674) = 1.775300d-26
  td_E_eV(675) = 8.160000d+0 ; td_cs_cm2(675) = 1.760100d-26
  td_E_eV(676) = 8.170000d+0 ; td_cs_cm2(676) = 1.745100d-26
  td_E_eV(677) = 8.180000d+0 ; td_cs_cm2(677) = 1.730200d-26
  td_E_eV(678) = 8.190000d+0 ; td_cs_cm2(678) = 1.715500d-26
  td_E_eV(679) = 8.200000d+0 ; td_cs_cm2(679) = 1.701000d-26
  td_E_eV(680) = 8.210000d+0 ; td_cs_cm2(680) = 1.686700d-26
  td_E_eV(681) = 8.220000d+0 ; td_cs_cm2(681) = 1.672500d-26
  td_E_eV(682) = 8.230000d+0 ; td_cs_cm2(682) = 1.658500d-26
  td_E_eV(683) = 8.240000d+0 ; td_cs_cm2(683) = 1.644700d-26
  td_E_eV(684) = 8.250000d+0 ; td_cs_cm2(684) = 1.631000d-26
  td_E_eV(685) = 8.260000d+0 ; td_cs_cm2(685) = 1.617400d-26
  td_E_eV(686) = 8.270000d+0 ; td_cs_cm2(686) = 1.604000d-26
  td_E_eV(687) = 8.280000d+0 ; td_cs_cm2(687) = 1.590800d-26
  td_E_eV(688) = 8.290000d+0 ; td_cs_cm2(688) = 1.577600d-26
  td_E_eV(689) = 8.300000d+0 ; td_cs_cm2(689) = 1.564600d-26
  td_E_eV(690) = 8.310000d+0 ; td_cs_cm2(690) = 1.551800d-26
  td_E_eV(691) = 8.320000d+0 ; td_cs_cm2(691) = 1.539100d-26
  td_E_eV(692) = 8.330000d+0 ; td_cs_cm2(692) = 1.526500d-26
  td_E_eV(693) = 8.340000d+0 ; td_cs_cm2(693) = 1.514100d-26
  td_E_eV(694) = 8.350000d+0 ; td_cs_cm2(694) = 1.501800d-26
  td_E_eV(695) = 8.360000d+0 ; td_cs_cm2(695) = 1.489700d-26
  td_E_eV(696) = 8.370000d+0 ; td_cs_cm2(696) = 1.477700d-26
  td_E_eV(697) = 8.380000d+0 ; td_cs_cm2(697) = 1.465900d-26
  td_E_eV(698) = 8.390000d+0 ; td_cs_cm2(698) = 1.454100d-26
  td_E_eV(699) = 8.400000d+0 ; td_cs_cm2(699) = 1.442500d-26
  td_E_eV(700) = 8.410000d+0 ; td_cs_cm2(700) = 1.431000d-26
  td_E_eV(701) = 8.420000d+0 ; td_cs_cm2(701) = 1.419600d-26
  td_E_eV(702) = 8.430000d+0 ; td_cs_cm2(702) = 1.408300d-26
  td_E_eV(703) = 8.440000d+0 ; td_cs_cm2(703) = 1.397200d-26
  td_E_eV(704) = 8.450000d+0 ; td_cs_cm2(704) = 1.386100d-26
  td_E_eV(705) = 8.460000d+0 ; td_cs_cm2(705) = 1.375200d-26
  td_E_eV(706) = 8.470000d+0 ; td_cs_cm2(706) = 1.364400d-26
  td_E_eV(707) = 8.480000d+0 ; td_cs_cm2(707) = 1.353700d-26
  td_E_eV(708) = 8.490000d+0 ; td_cs_cm2(708) = 1.343200d-26
  td_E_eV(709) = 8.500000d+0 ; td_cs_cm2(709) = 1.332800d-26
  td_E_eV(710) = 8.510000d+0 ; td_cs_cm2(710) = 1.322500d-26
  td_E_eV(711) = 8.520000d+0 ; td_cs_cm2(711) = 1.312200d-26
  td_E_eV(712) = 8.530000d+0 ; td_cs_cm2(712) = 1.302100d-26
  td_E_eV(713) = 8.540000d+0 ; td_cs_cm2(713) = 1.292100d-26
  td_E_eV(714) = 8.550000d+0 ; td_cs_cm2(714) = 1.282100d-26
  td_E_eV(715) = 8.560000d+0 ; td_cs_cm2(715) = 1.272300d-26
  td_E_eV(716) = 8.570000d+0 ; td_cs_cm2(716) = 1.262500d-26
  td_E_eV(717) = 8.580000d+0 ; td_cs_cm2(717) = 1.252900d-26
  td_E_eV(718) = 8.590000d+0 ; td_cs_cm2(718) = 1.243300d-26
  td_E_eV(719) = 8.600000d+0 ; td_cs_cm2(719) = 1.233900d-26
  td_E_eV(720) = 8.610000d+0 ; td_cs_cm2(720) = 1.224600d-26
  td_E_eV(721) = 8.620000d+0 ; td_cs_cm2(721) = 1.215400d-26
  td_E_eV(722) = 8.630000d+0 ; td_cs_cm2(722) = 1.206300d-26
  td_E_eV(723) = 8.640000d+0 ; td_cs_cm2(723) = 1.197300d-26
  td_E_eV(724) = 8.650000d+0 ; td_cs_cm2(724) = 1.188400d-26
  td_E_eV(725) = 8.660000d+0 ; td_cs_cm2(725) = 1.179500d-26
  td_E_eV(726) = 8.670000d+0 ; td_cs_cm2(726) = 1.170700d-26
  td_E_eV(727) = 8.680000d+0 ; td_cs_cm2(727) = 1.161900d-26
  td_E_eV(728) = 8.690000d+0 ; td_cs_cm2(728) = 1.153300d-26
  td_E_eV(729) = 8.700000d+0 ; td_cs_cm2(729) = 1.144800d-26
  td_E_eV(730) = 8.710000d+0 ; td_cs_cm2(730) = 1.136300d-26
  td_E_eV(731) = 8.720000d+0 ; td_cs_cm2(731) = 1.128000d-26
  td_E_eV(732) = 8.730000d+0 ; td_cs_cm2(732) = 1.119800d-26
  td_E_eV(733) = 8.740000d+0 ; td_cs_cm2(733) = 1.111600d-26
  td_E_eV(734) = 8.750000d+0 ; td_cs_cm2(734) = 1.103500d-26
  td_E_eV(735) = 8.760000d+0 ; td_cs_cm2(735) = 1.095500d-26
  td_E_eV(736) = 8.770000d+0 ; td_cs_cm2(736) = 1.087600d-26
  td_E_eV(737) = 8.780000d+0 ; td_cs_cm2(737) = 1.079700d-26
  td_E_eV(738) = 8.790000d+0 ; td_cs_cm2(738) = 1.071900d-26
  td_E_eV(739) = 8.800000d+0 ; td_cs_cm2(739) = 1.064100d-26
  td_E_eV(740) = 8.810000d+0 ; td_cs_cm2(740) = 1.056500d-26
  td_E_eV(741) = 8.820000d+0 ; td_cs_cm2(741) = 1.048900d-26
  td_E_eV(742) = 8.830000d+0 ; td_cs_cm2(742) = 1.041400d-26
  td_E_eV(743) = 8.840000d+0 ; td_cs_cm2(743) = 1.034000d-26
  td_E_eV(744) = 8.850000d+0 ; td_cs_cm2(744) = 1.026700d-26
  td_E_eV(745) = 8.860000d+0 ; td_cs_cm2(745) = 1.019400d-26
  td_E_eV(746) = 8.870000d+0 ; td_cs_cm2(746) = 1.012300d-26
  td_E_eV(747) = 8.880000d+0 ; td_cs_cm2(747) = 1.005100d-26
  td_E_eV(748) = 8.890000d+0 ; td_cs_cm2(748) = 0.000000d+0

  td_cs_cm2 = 1.0d4 * td_cs_cm2 ! the original data are in m^2, convert to cm^2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_vib04
!
!-------------------------------------
!
real(8) function CSV_N2_vib04_m3s(E_eV)

  use N2_vib04
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_vib04_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_vib04_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_vib04_m3s












!---------------------------------------------------------------------------------------------------- #8
!#VIBRATIONAL
!#SPECIES: e / N2
!#PROCESS: E + N2 <-> E + N2*(v=0-5), Vibrational
!#PARAM.:  E = 1.408 eV, g1/g0 = 1, complete set
!#COMMENT: Vibrational excitation : N2(1Σg+:v=0) -> N2*(v=5) | Source: [Laporta et al, Plasma
!#COMMENT: Sources Sci. Technol. 23, 065002 (2014)].
!#UPDATED: 2023-12-23 21:19:17
!
module N2_vib05
  integer, parameter :: N_td=555  !8   ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_vib05
!
!-------------------------------------
!
subroutine Prepare_N2_vib05

  use N2_vib05
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

!  td_E_eV(1) = 1.3_8; td_cs_cm2(1) = 1.11d-18
!  td_E_eV(2) = 1.6_8; td_cs_cm2(2) = 8.84d-17
!  td_E_eV(3) = 2.0_8; td_cs_cm2(3) = 1.01d-15
!  td_E_eV(4) = 2.4_8; td_cs_cm2(4) = 1.35d-15
!  td_E_eV(5) = 3.0_8; td_cs_cm2(5) = 5.86d-16
!  td_E_eV(6) = 4.0_8; td_cs_cm2(6) = 1.37d-16
!  td_E_eV(7) = 5.0_8; td_cs_cm2(7) = 2.34d-17
!  td_E_eV(8) = 6.0_8; td_cs_cm2(8) = 3.66d-18

  td_E_eV(  1) = 1.630000d+0 ; td_cs_cm2(  1) = 0.000000d+0
  td_E_eV(  2) = 1.640000d+0 ; td_cs_cm2(  2) = 1.056800d-26
  td_E_eV(  3) = 1.650000d+0 ; td_cs_cm2(  3) = 1.383300d-26
  td_E_eV(  4) = 1.660000d+0 ; td_cs_cm2(  4) = 1.811400d-26
  td_E_eV(  5) = 1.670000d+0 ; td_cs_cm2(  5) = 2.373400d-26
  td_E_eV(  6) = 1.680000d+0 ; td_cs_cm2(  6) = 3.112500d-26
  td_E_eV(  7) = 1.690000d+0 ; td_cs_cm2(  7) = 4.086500d-26
  td_E_eV(  8) = 1.700000d+0 ; td_cs_cm2(  8) = 5.372600d-26
  td_E_eV(  9) = 1.710000d+0 ; td_cs_cm2(  9) = 7.075200d-26
  td_E_eV( 10) = 1.720000d+0 ; td_cs_cm2( 10) = 9.334600d-26
  td_E_eV( 11) = 1.730000d+0 ; td_cs_cm2( 11) = 1.234200d-25
  td_E_eV( 12) = 1.740000d+0 ; td_cs_cm2( 12) = 1.635500d-25
  td_E_eV( 13) = 1.750000d+0 ; td_cs_cm2( 13) = 2.173000d-25
  td_E_eV( 14) = 1.760000d+0 ; td_cs_cm2( 14) = 2.895200d-25
  td_E_eV( 15) = 1.770000d+0 ; td_cs_cm2( 15) = 3.868900d-25
  td_E_eV( 16) = 1.780000d+0 ; td_cs_cm2( 16) = 5.186600d-25
  td_E_eV( 17) = 1.790000d+0 ; td_cs_cm2( 17) = 6.976300d-25
  td_E_eV( 18) = 1.800000d+0 ; td_cs_cm2( 18) = 9.416200d-25
  td_E_eV( 19) = 1.810000d+0 ; td_cs_cm2( 19) = 1.275400d-24
  td_E_eV( 20) = 1.820000d+0 ; td_cs_cm2( 20) = 1.733700d-24
  td_E_eV( 21) = 1.830000d+0 ; td_cs_cm2( 21) = 2.364700d-24
  td_E_eV( 22) = 1.840000d+0 ; td_cs_cm2( 22) = 3.235400d-24
  td_E_eV( 23) = 1.850000d+0 ; td_cs_cm2( 23) = 4.438500d-24
  td_E_eV( 24) = 1.860000d+0 ; td_cs_cm2( 24) = 6.100600d-24
  td_E_eV( 25) = 1.870000d+0 ; td_cs_cm2( 25) = 8.391800d-24
  td_E_eV( 26) = 1.880000d+0 ; td_cs_cm2( 26) = 1.153500d-23
  td_E_eV( 27) = 1.890000d+0 ; td_cs_cm2( 27) = 1.581300d-23
  td_E_eV( 28) = 1.900000d+0 ; td_cs_cm2( 28) = 2.156200d-23
  td_E_eV( 29) = 1.910000d+0 ; td_cs_cm2( 29) = 2.915500d-23
  td_E_eV( 30) = 1.920000d+0 ; td_cs_cm2( 30) = 3.896300d-23
  td_E_eV( 31) = 1.930000d+0 ; td_cs_cm2( 31) = 5.129100d-23
  td_E_eV( 32) = 1.940000d+0 ; td_cs_cm2( 32) = 6.632000d-23
  td_E_eV( 33) = 1.950000d+0 ; td_cs_cm2( 33) = 8.407400d-23
  td_E_eV( 34) = 1.960000d+0 ; td_cs_cm2( 34) = 1.044500d-22
  td_E_eV( 35) = 1.970000d+0 ; td_cs_cm2( 35) = 1.273100d-22
  td_E_eV( 36) = 1.980000d+0 ; td_cs_cm2( 36) = 1.526100d-22
  td_E_eV( 37) = 1.990000d+0 ; td_cs_cm2( 37) = 1.805200d-22
  td_E_eV( 38) = 2.000000d+0 ; td_cs_cm2( 38) = 2.114800d-22
  td_E_eV( 39) = 2.010000d+0 ; td_cs_cm2( 39) = 2.462500d-22
  td_E_eV( 40) = 2.020000d+0 ; td_cs_cm2( 40) = 2.859100d-22
  td_E_eV( 41) = 2.030000d+0 ; td_cs_cm2( 41) = 3.318800d-22
  td_E_eV( 42) = 2.040000d+0 ; td_cs_cm2( 42) = 3.859300d-22
  td_E_eV( 43) = 2.050000d+0 ; td_cs_cm2( 43) = 4.502600d-22
  td_E_eV( 44) = 2.060000d+0 ; td_cs_cm2( 44) = 5.275300d-22
  td_E_eV( 45) = 2.070000d+0 ; td_cs_cm2( 45) = 6.209600d-22
  td_E_eV( 46) = 2.080000d+0 ; td_cs_cm2( 46) = 7.343800d-22
  td_E_eV( 47) = 2.090000d+0 ; td_cs_cm2( 47) = 8.722200d-22
  td_E_eV( 48) = 2.100000d+0 ; td_cs_cm2( 48) = 1.039400d-21
  td_E_eV( 49) = 2.110000d+0 ; td_cs_cm2( 49) = 1.240900d-21
  td_E_eV( 50) = 2.120000d+0 ; td_cs_cm2( 50) = 1.481300d-21
  td_E_eV( 51) = 2.130000d+0 ; td_cs_cm2( 51) = 1.763700d-21
  td_E_eV( 52) = 2.140000d+0 ; td_cs_cm2( 52) = 2.088100d-21
  td_E_eV( 53) = 2.150000d+0 ; td_cs_cm2( 53) = 2.450200d-21
  td_E_eV( 54) = 2.160000d+0 ; td_cs_cm2( 54) = 2.840300d-21
  td_E_eV( 55) = 2.170000d+0 ; td_cs_cm2( 55) = 3.243400d-21
  td_E_eV( 56) = 2.180000d+0 ; td_cs_cm2( 56) = 3.642000d-21
  td_E_eV( 57) = 2.190000d+0 ; td_cs_cm2( 57) = 4.019200d-21
  td_E_eV( 58) = 2.200000d+0 ; td_cs_cm2( 58) = 4.363000d-21
  td_E_eV( 59) = 2.210000d+0 ; td_cs_cm2( 59) = 4.668500d-21
  td_E_eV( 60) = 2.220000d+0 ; td_cs_cm2( 60) = 4.938700d-21
  td_E_eV( 61) = 2.230000d+0 ; td_cs_cm2( 61) = 5.182400d-21
  td_E_eV( 62) = 2.240000d+0 ; td_cs_cm2( 62) = 5.411700d-21
  td_E_eV( 63) = 2.250000d+0 ; td_cs_cm2( 63) = 5.639900d-21
  td_E_eV( 64) = 2.260000d+0 ; td_cs_cm2( 64) = 5.880200d-21
  td_E_eV( 65) = 2.270000d+0 ; td_cs_cm2( 65) = 6.144300d-21
  td_E_eV( 66) = 2.280000d+0 ; td_cs_cm2( 66) = 6.442500d-21
  td_E_eV( 67) = 2.290000d+0 ; td_cs_cm2( 67) = 6.783800d-21
  td_E_eV( 68) = 2.300000d+0 ; td_cs_cm2( 68) = 7.175100d-21
  td_E_eV( 69) = 2.310000d+0 ; td_cs_cm2( 69) = 7.621100d-21
  td_E_eV( 70) = 2.320000d+0 ; td_cs_cm2( 70) = 8.123300d-21
  td_E_eV( 71) = 2.330000d+0 ; td_cs_cm2( 71) = 8.678600d-21
  td_E_eV( 72) = 2.340000d+0 ; td_cs_cm2( 72) = 9.277200d-21
  td_E_eV( 73) = 2.350000d+0 ; td_cs_cm2( 73) = 9.900300d-21
  td_E_eV( 74) = 2.360000d+0 ; td_cs_cm2( 74) = 1.051800d-20
  td_E_eV( 75) = 2.370000d+0 ; td_cs_cm2( 75) = 1.109000d-20
  td_E_eV( 76) = 2.380000d+0 ; td_cs_cm2( 76) = 1.156500d-20
  td_E_eV( 77) = 2.390000d+0 ; td_cs_cm2( 77) = 1.189300d-20
  td_E_eV( 78) = 2.400000d+0 ; td_cs_cm2( 78) = 1.203000d-20
  td_E_eV( 79) = 2.410000d+0 ; td_cs_cm2( 79) = 1.195100d-20
  td_E_eV( 80) = 2.420000d+0 ; td_cs_cm2( 80) = 1.165800d-20
  td_E_eV( 81) = 2.430000d+0 ; td_cs_cm2( 81) = 1.118000d-20
  td_E_eV( 82) = 2.440000d+0 ; td_cs_cm2( 82) = 1.055900d-20
  td_E_eV( 83) = 2.450000d+0 ; td_cs_cm2( 83) = 9.848400d-21
  td_E_eV( 84) = 2.460000d+0 ; td_cs_cm2( 84) = 9.094800d-21
  td_E_eV( 85) = 2.470000d+0 ; td_cs_cm2( 85) = 8.336400d-21
  td_E_eV( 86) = 2.480000d+0 ; td_cs_cm2( 86) = 7.599400d-21
  td_E_eV( 87) = 2.490000d+0 ; td_cs_cm2( 87) = 6.899100d-21
  td_E_eV( 88) = 2.500000d+0 ; td_cs_cm2( 88) = 6.242000d-21
  td_E_eV( 89) = 2.510000d+0 ; td_cs_cm2( 89) = 5.628300d-21
  td_E_eV( 90) = 2.520000d+0 ; td_cs_cm2( 90) = 5.054400d-21
  td_E_eV( 91) = 2.530000d+0 ; td_cs_cm2( 91) = 4.514200d-21
  td_E_eV( 92) = 2.540000d+0 ; td_cs_cm2( 92) = 4.000600d-21
  td_E_eV( 93) = 2.550000d+0 ; td_cs_cm2( 93) = 3.507000d-21
  td_E_eV( 94) = 2.560000d+0 ; td_cs_cm2( 94) = 3.027800d-21
  td_E_eV( 95) = 2.570000d+0 ; td_cs_cm2( 95) = 2.559800d-21
  td_E_eV( 96) = 2.580000d+0 ; td_cs_cm2( 96) = 2.103600d-21
  td_E_eV( 97) = 2.590000d+0 ; td_cs_cm2( 97) = 1.664000d-21
  td_E_eV( 98) = 2.600000d+0 ; td_cs_cm2( 98) = 1.250900d-21
  td_E_eV( 99) = 2.610000d+0 ; td_cs_cm2( 99) = 8.785500d-22
  td_E_eV(100) = 2.620000d+0 ; td_cs_cm2(100) = 5.634600d-22
  td_E_eV(101) = 2.630000d+0 ; td_cs_cm2(101) = 3.209500d-22
  td_E_eV(102) = 2.640000d+0 ; td_cs_cm2(102) = 1.616500d-22
  td_E_eV(103) = 2.650000d+0 ; td_cs_cm2(103) = 8.912800d-23
  td_E_eV(104) = 2.660000d+0 ; td_cs_cm2(104) = 9.965200d-23
  td_E_eV(105) = 2.670000d+0 ; td_cs_cm2(105) = 1.839300d-22
  td_E_eV(106) = 2.680000d+0 ; td_cs_cm2(106) = 3.298800d-22
  td_E_eV(107) = 2.690000d+0 ; td_cs_cm2(107) = 5.252800d-22
  td_E_eV(108) = 2.700000d+0 ; td_cs_cm2(108) = 7.595700d-22
  td_E_eV(109) = 2.710000d+0 ; td_cs_cm2(109) = 1.024600d-21
  td_E_eV(110) = 2.720000d+0 ; td_cs_cm2(110) = 1.314800d-21
  td_E_eV(111) = 2.730000d+0 ; td_cs_cm2(111) = 1.626500d-21
  td_E_eV(112) = 2.740000d+0 ; td_cs_cm2(112) = 1.957400d-21
  td_E_eV(113) = 2.750000d+0 ; td_cs_cm2(113) = 2.305400d-21
  td_E_eV(114) = 2.760000d+0 ; td_cs_cm2(114) = 2.667900d-21
  td_E_eV(115) = 2.770000d+0 ; td_cs_cm2(115) = 3.040700d-21
  td_E_eV(116) = 2.780000d+0 ; td_cs_cm2(116) = 3.416500d-21
  td_E_eV(117) = 2.790000d+0 ; td_cs_cm2(117) = 3.784200d-21
  td_E_eV(118) = 2.800000d+0 ; td_cs_cm2(118) = 4.128100d-21
  td_E_eV(119) = 2.810000d+0 ; td_cs_cm2(119) = 4.428000d-21
  td_E_eV(120) = 2.820000d+0 ; td_cs_cm2(120) = 4.661100d-21
  td_E_eV(121) = 2.830000d+0 ; td_cs_cm2(121) = 4.804900d-21
  td_E_eV(122) = 2.840000d+0 ; td_cs_cm2(122) = 4.841900d-21
  td_E_eV(123) = 2.850000d+0 ; td_cs_cm2(123) = 4.764500d-21
  td_E_eV(124) = 2.860000d+0 ; td_cs_cm2(124) = 4.576700d-21
  td_E_eV(125) = 2.870000d+0 ; td_cs_cm2(125) = 4.294600d-21
  td_E_eV(126) = 2.880000d+0 ; td_cs_cm2(126) = 3.941800d-21
  td_E_eV(127) = 2.890000d+0 ; td_cs_cm2(127) = 3.545000d-21
  td_E_eV(128) = 2.900000d+0 ; td_cs_cm2(128) = 3.129100d-21
  td_E_eV(129) = 2.910000d+0 ; td_cs_cm2(129) = 2.714400d-21
  td_E_eV(130) = 2.920000d+0 ; td_cs_cm2(130) = 2.315600d-21
  td_E_eV(131) = 2.930000d+0 ; td_cs_cm2(131) = 1.941900d-21
  td_E_eV(132) = 2.940000d+0 ; td_cs_cm2(132) = 1.598500d-21
  td_E_eV(133) = 2.950000d+0 ; td_cs_cm2(133) = 1.287800d-21
  td_E_eV(134) = 2.960000d+0 ; td_cs_cm2(134) = 1.010500d-21
  td_E_eV(135) = 2.970000d+0 ; td_cs_cm2(135) = 7.668900d-22
  td_E_eV(136) = 2.980000d+0 ; td_cs_cm2(136) = 5.571000d-22
  td_E_eV(137) = 2.990000d+0 ; td_cs_cm2(137) = 3.820100d-22
  td_E_eV(138) = 3.000000d+0 ; td_cs_cm2(138) = 2.432600d-22
  td_E_eV(139) = 3.010000d+0 ; td_cs_cm2(139) = 1.431500d-22
  td_E_eV(140) = 3.020000d+0 ; td_cs_cm2(140) = 8.424700d-23
  td_E_eV(141) = 3.030000d+0 ; td_cs_cm2(141) = 6.854800d-23
  td_E_eV(142) = 3.040000d+0 ; td_cs_cm2(142) = 9.642000d-23
  td_E_eV(143) = 3.050000d+0 ; td_cs_cm2(143) = 1.655100d-22
  td_E_eV(144) = 3.060000d+0 ; td_cs_cm2(144) = 2.701300d-22
  td_E_eV(145) = 3.070000d+0 ; td_cs_cm2(145) = 4.013800d-22
  td_E_eV(146) = 3.080000d+0 ; td_cs_cm2(146) = 5.483500d-22
  td_E_eV(147) = 3.090000d+0 ; td_cs_cm2(147) = 6.997700d-22
  td_E_eV(148) = 3.100000d+0 ; td_cs_cm2(148) = 8.457600d-22
  td_E_eV(149) = 3.110000d+0 ; td_cs_cm2(149) = 9.789100d-22
  td_E_eV(150) = 3.120000d+0 ; td_cs_cm2(150) = 1.094600d-21
  td_E_eV(151) = 3.130000d+0 ; td_cs_cm2(151) = 1.190700d-21
  td_E_eV(152) = 3.140000d+0 ; td_cs_cm2(152) = 1.266800d-21
  td_E_eV(153) = 3.150000d+0 ; td_cs_cm2(153) = 1.323600d-21
  td_E_eV(154) = 3.160000d+0 ; td_cs_cm2(154) = 1.362100d-21
  td_E_eV(155) = 3.170000d+0 ; td_cs_cm2(155) = 1.383100d-21
  td_E_eV(156) = 3.180000d+0 ; td_cs_cm2(156) = 1.387100d-21
  td_E_eV(157) = 3.190000d+0 ; td_cs_cm2(157) = 1.374300d-21
  td_E_eV(158) = 3.200000d+0 ; td_cs_cm2(158) = 1.344200d-21
  td_E_eV(159) = 3.210000d+0 ; td_cs_cm2(159) = 1.295800d-21
  td_E_eV(160) = 3.220000d+0 ; td_cs_cm2(160) = 1.228500d-21
  td_E_eV(161) = 3.230000d+0 ; td_cs_cm2(161) = 1.141900d-21
  td_E_eV(162) = 3.240000d+0 ; td_cs_cm2(162) = 1.036700d-21
  td_E_eV(163) = 3.250000d+0 ; td_cs_cm2(163) = 9.153100d-22
  td_E_eV(164) = 3.260000d+0 ; td_cs_cm2(164) = 7.820200d-22
  td_E_eV(165) = 3.270000d+0 ; td_cs_cm2(165) = 6.431900d-22
  td_E_eV(166) = 3.280000d+0 ; td_cs_cm2(166) = 5.064100d-22
  td_E_eV(167) = 3.290000d+0 ; td_cs_cm2(167) = 3.794100d-22
  td_E_eV(168) = 3.300000d+0 ; td_cs_cm2(168) = 2.687900d-22
  td_E_eV(169) = 3.310000d+0 ; td_cs_cm2(169) = 1.790900d-22
  td_E_eV(170) = 3.320000d+0 ; td_cs_cm2(170) = 1.124800d-22
  td_E_eV(171) = 3.330000d+0 ; td_cs_cm2(171) = 6.894300d-23
  td_E_eV(172) = 3.340000d+0 ; td_cs_cm2(172) = 4.693700d-23
  td_E_eV(173) = 3.350000d+0 ; td_cs_cm2(173) = 4.401400d-23
  td_E_eV(174) = 3.360000d+0 ; td_cs_cm2(174) = 5.740100d-23
  td_E_eV(175) = 3.370000d+0 ; td_cs_cm2(175) = 8.435300d-23
  td_E_eV(176) = 3.380000d+0 ; td_cs_cm2(176) = 1.223200d-22
  td_E_eV(177) = 3.390000d+0 ; td_cs_cm2(177) = 1.689600d-22
  td_E_eV(178) = 3.400000d+0 ; td_cs_cm2(178) = 2.220300d-22
  td_E_eV(179) = 3.410000d+0 ; td_cs_cm2(179) = 2.792400d-22
  td_E_eV(180) = 3.420000d+0 ; td_cs_cm2(180) = 3.380900d-22
  td_E_eV(181) = 3.430000d+0 ; td_cs_cm2(181) = 3.957000d-22
  td_E_eV(182) = 3.440000d+0 ; td_cs_cm2(182) = 4.487600d-22
  td_E_eV(183) = 3.450000d+0 ; td_cs_cm2(183) = 4.936700d-22
  td_E_eV(184) = 3.460000d+0 ; td_cs_cm2(184) = 5.268300d-22
  td_E_eV(185) = 3.470000d+0 ; td_cs_cm2(185) = 5.452000d-22
  td_E_eV(186) = 3.480000d+0 ; td_cs_cm2(186) = 5.469000d-22
  td_E_eV(187) = 3.490000d+0 ; td_cs_cm2(187) = 5.316900d-22
  td_E_eV(188) = 3.500000d+0 ; td_cs_cm2(188) = 5.011400d-22
  td_E_eV(189) = 3.510000d+0 ; td_cs_cm2(189) = 4.582600d-22
  td_E_eV(190) = 3.520000d+0 ; td_cs_cm2(190) = 4.069800d-22
  td_E_eV(191) = 3.530000d+0 ; td_cs_cm2(191) = 3.513500d-22
  td_E_eV(192) = 3.540000d+0 ; td_cs_cm2(192) = 2.950200d-22
  td_E_eV(193) = 3.550000d+0 ; td_cs_cm2(193) = 2.409000d-22
  td_E_eV(194) = 3.560000d+0 ; td_cs_cm2(194) = 1.910800d-22
  td_E_eV(195) = 3.570000d+0 ; td_cs_cm2(195) = 1.469400d-22
  td_E_eV(196) = 3.580000d+0 ; td_cs_cm2(196) = 1.092800d-22
  td_E_eV(197) = 3.590000d+0 ; td_cs_cm2(197) = 7.856000d-23
  td_E_eV(198) = 3.600000d+0 ; td_cs_cm2(198) = 5.497400d-23
  td_E_eV(199) = 3.610000d+0 ; td_cs_cm2(199) = 3.859400d-23
  td_E_eV(200) = 3.620000d+0 ; td_cs_cm2(200) = 2.938800d-23
  td_E_eV(201) = 3.630000d+0 ; td_cs_cm2(201) = 2.720200d-23
  td_E_eV(202) = 3.640000d+0 ; td_cs_cm2(202) = 3.169900d-23
  td_E_eV(203) = 3.650000d+0 ; td_cs_cm2(203) = 4.225900d-23
  td_E_eV(204) = 3.660000d+0 ; td_cs_cm2(204) = 5.788300d-23
  td_E_eV(205) = 3.670000d+0 ; td_cs_cm2(205) = 7.715200d-23
  td_E_eV(206) = 3.680000d+0 ; td_cs_cm2(206) = 9.827500d-23
  td_E_eV(207) = 3.690000d+0 ; td_cs_cm2(207) = 1.192700d-22
  td_E_eV(208) = 3.700000d+0 ; td_cs_cm2(208) = 1.382200d-22
  td_E_eV(209) = 3.710000d+0 ; td_cs_cm2(209) = 1.535800d-22
  td_E_eV(210) = 3.720000d+0 ; td_cs_cm2(210) = 1.643400d-22
  td_E_eV(211) = 3.730000d+0 ; td_cs_cm2(211) = 1.701100d-22
  td_E_eV(212) = 3.740000d+0 ; td_cs_cm2(212) = 1.710200d-22
  td_E_eV(213) = 3.750000d+0 ; td_cs_cm2(213) = 1.675700d-22
  td_E_eV(214) = 3.760000d+0 ; td_cs_cm2(214) = 1.604800d-22
  td_E_eV(215) = 3.770000d+0 ; td_cs_cm2(215) = 1.505000d-22
  td_E_eV(216) = 3.780000d+0 ; td_cs_cm2(216) = 1.383600d-22
  td_E_eV(217) = 3.790000d+0 ; td_cs_cm2(217) = 1.247100d-22
  td_E_eV(218) = 3.800000d+0 ; td_cs_cm2(218) = 1.101400d-22
  td_E_eV(219) = 3.810000d+0 ; td_cs_cm2(219) = 9.514300d-23
  td_E_eV(220) = 3.820000d+0 ; td_cs_cm2(220) = 8.021300d-23
  td_E_eV(221) = 3.830000d+0 ; td_cs_cm2(221) = 6.582300d-23
  td_E_eV(222) = 3.840000d+0 ; td_cs_cm2(222) = 5.246600d-23
  td_E_eV(223) = 3.850000d+0 ; td_cs_cm2(223) = 4.065000d-23
  td_E_eV(224) = 3.860000d+0 ; td_cs_cm2(224) = 3.087000d-23
  td_E_eV(225) = 3.870000d+0 ; td_cs_cm2(225) = 2.355200d-23
  td_E_eV(226) = 3.880000d+0 ; td_cs_cm2(226) = 1.897900d-23
  td_E_eV(227) = 3.890000d+0 ; td_cs_cm2(227) = 1.721400d-23
  td_E_eV(228) = 3.900000d+0 ; td_cs_cm2(228) = 1.806200d-23
  td_E_eV(229) = 3.910000d+0 ; td_cs_cm2(229) = 2.109300d-23
  td_E_eV(230) = 3.920000d+0 ; td_cs_cm2(230) = 2.570800d-23
  td_E_eV(231) = 3.930000d+0 ; td_cs_cm2(231) = 3.124600d-23
  td_E_eV(232) = 3.940000d+0 ; td_cs_cm2(232) = 3.708000d-23
  td_E_eV(233) = 3.950000d+0 ; td_cs_cm2(233) = 4.268800d-23
  td_E_eV(234) = 3.960000d+0 ; td_cs_cm2(234) = 4.767500d-23
  td_E_eV(235) = 3.970000d+0 ; td_cs_cm2(235) = 5.177000d-23
  td_E_eV(236) = 3.980000d+0 ; td_cs_cm2(236) = 5.480800d-23
  td_E_eV(237) = 3.990000d+0 ; td_cs_cm2(237) = 5.669300d-23
  td_E_eV(238) = 4.000000d+0 ; td_cs_cm2(238) = 5.738200d-23
  td_E_eV(239) = 4.010000d+0 ; td_cs_cm2(239) = 5.686500d-23
  td_E_eV(240) = 4.020000d+0 ; td_cs_cm2(240) = 5.515900d-23
  td_E_eV(241) = 4.030000d+0 ; td_cs_cm2(241) = 5.231600d-23
  td_E_eV(242) = 4.040000d+0 ; td_cs_cm2(242) = 4.843200d-23
  td_E_eV(243) = 4.050000d+0 ; td_cs_cm2(243) = 4.366300d-23
  td_E_eV(244) = 4.060000d+0 ; td_cs_cm2(244) = 3.824300d-23
  td_E_eV(245) = 4.070000d+0 ; td_cs_cm2(245) = 3.248000d-23
  td_E_eV(246) = 4.080000d+0 ; td_cs_cm2(246) = 2.673800d-23
  td_E_eV(247) = 4.090000d+0 ; td_cs_cm2(247) = 2.140000d-23
  td_E_eV(248) = 4.100000d+0 ; td_cs_cm2(248) = 1.680800d-23
  td_E_eV(249) = 4.110000d+0 ; td_cs_cm2(249) = 1.321300d-23
  td_E_eV(250) = 4.120000d+0 ; td_cs_cm2(250) = 1.074200d-23
  td_E_eV(251) = 4.130000d+0 ; td_cs_cm2(251) = 9.400400d-24
  td_E_eV(252) = 4.140000d+0 ; td_cs_cm2(252) = 9.087100d-24
  td_E_eV(253) = 4.150000d+0 ; td_cs_cm2(253) = 9.636500d-24
  td_E_eV(254) = 4.160000d+0 ; td_cs_cm2(254) = 1.085100d-23
  td_E_eV(255) = 4.170000d+0 ; td_cs_cm2(255) = 1.252800d-23
  td_E_eV(256) = 4.180000d+0 ; td_cs_cm2(256) = 1.447600d-23
  td_E_eV(257) = 4.190000d+0 ; td_cs_cm2(257) = 1.652300d-23
  td_E_eV(258) = 4.200000d+0 ; td_cs_cm2(258) = 1.851200d-23
  td_E_eV(259) = 4.210000d+0 ; td_cs_cm2(259) = 2.030200d-23
  td_E_eV(260) = 4.220000d+0 ; td_cs_cm2(260) = 2.176200d-23
  td_E_eV(261) = 4.230000d+0 ; td_cs_cm2(261) = 2.277600d-23
  td_E_eV(262) = 4.240000d+0 ; td_cs_cm2(262) = 2.324800d-23
  td_E_eV(263) = 4.250000d+0 ; td_cs_cm2(263) = 2.310800d-23
  td_E_eV(264) = 4.260000d+0 ; td_cs_cm2(264) = 2.233300d-23
  td_E_eV(265) = 4.270000d+0 ; td_cs_cm2(265) = 2.095400d-23
  td_E_eV(266) = 4.280000d+0 ; td_cs_cm2(266) = 1.906700d-23
  td_E_eV(267) = 4.290000d+0 ; td_cs_cm2(267) = 1.682400d-23
  td_E_eV(268) = 4.300000d+0 ; td_cs_cm2(268) = 1.441400d-23
  td_E_eV(269) = 4.310000d+0 ; td_cs_cm2(269) = 1.203100d-23
  td_E_eV(270) = 4.320000d+0 ; td_cs_cm2(270) = 9.849200d-24
  td_E_eV(271) = 4.330000d+0 ; td_cs_cm2(271) = 7.998300d-24
  td_E_eV(272) = 4.340000d+0 ; td_cs_cm2(272) = 6.556800d-24
  td_E_eV(273) = 4.350000d+0 ; td_cs_cm2(273) = 5.555100d-24
  td_E_eV(274) = 4.360000d+0 ; td_cs_cm2(274) = 4.985000d-24
  td_E_eV(275) = 4.370000d+0 ; td_cs_cm2(275) = 4.810800d-24
  td_E_eV(276) = 4.380000d+0 ; td_cs_cm2(276) = 4.979400d-24
  td_E_eV(277) = 4.390000d+0 ; td_cs_cm2(277) = 5.427300d-24
  td_E_eV(278) = 4.400000d+0 ; td_cs_cm2(278) = 6.083800d-24
  td_E_eV(279) = 4.410000d+0 ; td_cs_cm2(279) = 6.872900d-24
  td_E_eV(280) = 4.420000d+0 ; td_cs_cm2(280) = 7.713700d-24
  td_E_eV(281) = 4.430000d+0 ; td_cs_cm2(281) = 8.522100d-24
  td_E_eV(282) = 4.440000d+0 ; td_cs_cm2(282) = 9.214000d-24
  td_E_eV(283) = 4.450000d+0 ; td_cs_cm2(283) = 9.711800d-24
  td_E_eV(284) = 4.460000d+0 ; td_cs_cm2(284) = 9.953900d-24
  td_E_eV(285) = 4.470000d+0 ; td_cs_cm2(285) = 9.904800d-24
  td_E_eV(286) = 4.480000d+0 ; td_cs_cm2(286) = 9.562500d-24
  td_E_eV(287) = 4.490000d+0 ; td_cs_cm2(287) = 8.960100d-24
  td_E_eV(288) = 4.500000d+0 ; td_cs_cm2(288) = 8.159800d-24
  td_E_eV(289) = 4.510000d+0 ; td_cs_cm2(289) = 7.240500d-24
  td_E_eV(290) = 4.520000d+0 ; td_cs_cm2(290) = 6.284300d-24
  td_E_eV(291) = 4.530000d+0 ; td_cs_cm2(291) = 5.365000d-24
  td_E_eV(292) = 4.540000d+0 ; td_cs_cm2(292) = 4.540600d-24
  td_E_eV(293) = 4.550000d+0 ; td_cs_cm2(293) = 3.851200d-24
  td_E_eV(294) = 4.560000d+0 ; td_cs_cm2(294) = 3.320300d-24
  td_E_eV(295) = 4.570000d+0 ; td_cs_cm2(295) = 2.957000d-24
  td_E_eV(296) = 4.580000d+0 ; td_cs_cm2(296) = 2.759500d-24
  td_E_eV(297) = 4.590000d+0 ; td_cs_cm2(297) = 2.716600d-24
  td_E_eV(298) = 4.600000d+0 ; td_cs_cm2(298) = 2.809800d-24
  td_E_eV(299) = 4.610000d+0 ; td_cs_cm2(299) = 3.013100d-24
  td_E_eV(300) = 4.620000d+0 ; td_cs_cm2(300) = 3.294000d-24
  td_E_eV(301) = 4.630000d+0 ; td_cs_cm2(301) = 3.613700d-24
  td_E_eV(302) = 4.640000d+0 ; td_cs_cm2(302) = 3.929700d-24
  td_E_eV(303) = 4.650000d+0 ; td_cs_cm2(303) = 4.199600d-24
  td_E_eV(304) = 4.660000d+0 ; td_cs_cm2(304) = 4.386100d-24
  td_E_eV(305) = 4.670000d+0 ; td_cs_cm2(305) = 4.463300d-24
  td_E_eV(306) = 4.680000d+0 ; td_cs_cm2(306) = 4.420000d-24
  td_E_eV(307) = 4.690000d+0 ; td_cs_cm2(307) = 4.261200d-24
  td_E_eV(308) = 4.700000d+0 ; td_cs_cm2(308) = 4.005600d-24
  td_E_eV(309) = 4.710000d+0 ; td_cs_cm2(309) = 3.680600d-24
  td_E_eV(310) = 4.720000d+0 ; td_cs_cm2(310) = 3.317200d-24
  td_E_eV(311) = 4.730000d+0 ; td_cs_cm2(311) = 2.945000d-24
  td_E_eV(312) = 4.740000d+0 ; td_cs_cm2(312) = 2.589700d-24
  td_E_eV(313) = 4.750000d+0 ; td_cs_cm2(313) = 2.271500d-24
  td_E_eV(314) = 4.760000d+0 ; td_cs_cm2(314) = 2.005100d-24
  td_E_eV(315) = 4.770000d+0 ; td_cs_cm2(315) = 1.800100d-24
  td_E_eV(316) = 4.780000d+0 ; td_cs_cm2(316) = 1.661300d-24
  td_E_eV(317) = 4.790000d+0 ; td_cs_cm2(317) = 1.589100d-24
  td_E_eV(318) = 4.800000d+0 ; td_cs_cm2(318) = 1.579700d-24
  td_E_eV(319) = 4.810000d+0 ; td_cs_cm2(319) = 1.624500d-24
  td_E_eV(320) = 4.820000d+0 ; td_cs_cm2(320) = 1.710500d-24
  td_E_eV(321) = 4.830000d+0 ; td_cs_cm2(321) = 1.820800d-24
  td_E_eV(322) = 4.840000d+0 ; td_cs_cm2(322) = 1.936100d-24
  td_E_eV(323) = 4.850000d+0 ; td_cs_cm2(323) = 2.036800d-24
  td_E_eV(324) = 4.860000d+0 ; td_cs_cm2(324) = 2.106600d-24
  td_E_eV(325) = 4.870000d+0 ; td_cs_cm2(325) = 2.134100d-24
  td_E_eV(326) = 4.880000d+0 ; td_cs_cm2(326) = 2.114700d-24
  td_E_eV(327) = 4.890000d+0 ; td_cs_cm2(327) = 2.050400d-24
  td_E_eV(328) = 4.900000d+0 ; td_cs_cm2(328) = 1.948100d-24
  td_E_eV(329) = 4.910000d+0 ; td_cs_cm2(329) = 1.818200d-24
  td_E_eV(330) = 4.920000d+0 ; td_cs_cm2(330) = 1.671700d-24
  td_E_eV(331) = 4.930000d+0 ; td_cs_cm2(331) = 1.519900d-24
  td_E_eV(332) = 4.940000d+0 ; td_cs_cm2(332) = 1.372600d-24
  td_E_eV(333) = 4.950000d+0 ; td_cs_cm2(333) = 1.238200d-24
  td_E_eV(334) = 4.960000d+0 ; td_cs_cm2(334) = 1.123100d-24
  td_E_eV(335) = 4.970000d+0 ; td_cs_cm2(335) = 1.032300d-24
  td_E_eV(336) = 4.980000d+0 ; td_cs_cm2(336) = 9.686100d-25
  td_E_eV(337) = 4.990000d+0 ; td_cs_cm2(337) = 9.330100d-25
  td_E_eV(338) = 5.000000d+0 ; td_cs_cm2(338) = 9.239900d-25
  td_E_eV(339) = 5.010000d+0 ; td_cs_cm2(339) = 9.375900d-25
  td_E_eV(340) = 5.020000d+0 ; td_cs_cm2(340) = 9.675300d-25
  td_E_eV(341) = 5.030000d+0 ; td_cs_cm2(341) = 1.005800d-24
  td_E_eV(342) = 5.040000d+0 ; td_cs_cm2(342) = 1.043900d-24
  td_E_eV(343) = 5.050000d+0 ; td_cs_cm2(343) = 1.073700d-24
  td_E_eV(344) = 5.060000d+0 ; td_cs_cm2(344) = 1.089400d-24
  td_E_eV(345) = 5.070000d+0 ; td_cs_cm2(345) = 1.087500d-24
  td_E_eV(346) = 5.080000d+0 ; td_cs_cm2(346) = 1.067400d-24
  td_E_eV(347) = 5.090000d+0 ; td_cs_cm2(347) = 1.030700d-24
  td_E_eV(348) = 5.100000d+0 ; td_cs_cm2(348) = 9.804800d-25
  td_E_eV(349) = 5.110000d+0 ; td_cs_cm2(349) = 9.207900d-25
  td_E_eV(350) = 5.120000d+0 ; td_cs_cm2(350) = 8.559200d-25
  td_E_eV(351) = 5.130000d+0 ; td_cs_cm2(351) = 7.900400d-25
  td_E_eV(352) = 5.140000d+0 ; td_cs_cm2(352) = 7.270000d-25
  td_E_eV(353) = 5.150000d+0 ; td_cs_cm2(353) = 6.701700d-25
  td_E_eV(354) = 5.160000d+0 ; td_cs_cm2(354) = 6.223400d-25
  td_E_eV(355) = 5.170000d+0 ; td_cs_cm2(355) = 5.855200d-25
  td_E_eV(356) = 5.180000d+0 ; td_cs_cm2(356) = 5.607800d-25
  td_E_eV(357) = 5.190000d+0 ; td_cs_cm2(357) = 5.480700d-25
  td_E_eV(358) = 5.200000d+0 ; td_cs_cm2(358) = 5.460300d-25
  td_E_eV(359) = 5.210000d+0 ; td_cs_cm2(359) = 5.521700d-25
  td_E_eV(360) = 5.220000d+0 ; td_cs_cm2(360) = 5.630900d-25
  td_E_eV(361) = 5.230000d+0 ; td_cs_cm2(361) = 5.750900d-25
  td_E_eV(362) = 5.240000d+0 ; td_cs_cm2(362) = 5.847200d-25
  td_E_eV(363) = 5.250000d+0 ; td_cs_cm2(363) = 5.893100d-25
  td_E_eV(364) = 5.260000d+0 ; td_cs_cm2(364) = 5.872400d-25
  td_E_eV(365) = 5.270000d+0 ; td_cs_cm2(365) = 5.779700d-25
  td_E_eV(366) = 5.280000d+0 ; td_cs_cm2(366) = 5.618400d-25
  td_E_eV(367) = 5.290000d+0 ; td_cs_cm2(367) = 5.398600d-25
  td_E_eV(368) = 5.300000d+0 ; td_cs_cm2(368) = 5.134400d-25
  td_E_eV(369) = 5.310000d+0 ; td_cs_cm2(369) = 4.842100d-25
  td_E_eV(370) = 5.320000d+0 ; td_cs_cm2(370) = 4.538800d-25
  td_E_eV(371) = 5.330000d+0 ; td_cs_cm2(371) = 4.241400d-25
  td_E_eV(372) = 5.340000d+0 ; td_cs_cm2(372) = 3.965800d-25
  td_E_eV(373) = 5.350000d+0 ; td_cs_cm2(373) = 3.726100d-25
  td_E_eV(374) = 5.360000d+0 ; td_cs_cm2(374) = 3.533100d-25
  td_E_eV(375) = 5.370000d+0 ; td_cs_cm2(375) = 3.393200d-25
  td_E_eV(376) = 5.380000d+0 ; td_cs_cm2(376) = 3.306900d-25
  td_E_eV(377) = 5.390000d+0 ; td_cs_cm2(377) = 3.268600d-25
  td_E_eV(378) = 5.400000d+0 ; td_cs_cm2(378) = 3.267200d-25
  td_E_eV(379) = 5.410000d+0 ; td_cs_cm2(379) = 3.287900d-25
  td_E_eV(380) = 5.420000d+0 ; td_cs_cm2(380) = 3.314900d-25
  td_E_eV(381) = 5.430000d+0 ; td_cs_cm2(381) = 3.333900d-25
  td_E_eV(382) = 5.440000d+0 ; td_cs_cm2(382) = 3.333800d-25
  td_E_eV(383) = 5.450000d+0 ; td_cs_cm2(383) = 3.307900d-25
  td_E_eV(384) = 5.460000d+0 ; td_cs_cm2(384) = 3.253400d-25
  td_E_eV(385) = 5.470000d+0 ; td_cs_cm2(385) = 3.170800d-25
  td_E_eV(386) = 5.480000d+0 ; td_cs_cm2(386) = 3.063600d-25
  td_E_eV(387) = 5.490000d+0 ; td_cs_cm2(387) = 2.936700d-25
  td_E_eV(388) = 5.500000d+0 ; td_cs_cm2(388) = 2.796600d-25
  td_E_eV(389) = 5.510000d+0 ; td_cs_cm2(389) = 2.650500d-25
  td_E_eV(390) = 5.520000d+0 ; td_cs_cm2(390) = 2.506200d-25
  td_E_eV(391) = 5.530000d+0 ; td_cs_cm2(391) = 2.371100d-25
  td_E_eV(392) = 5.540000d+0 ; td_cs_cm2(392) = 2.252100d-25
  td_E_eV(393) = 5.550000d+0 ; td_cs_cm2(393) = 2.154400d-25
  td_E_eV(394) = 5.560000d+0 ; td_cs_cm2(394) = 2.080800d-25
  td_E_eV(395) = 5.570000d+0 ; td_cs_cm2(395) = 2.031100d-25
  td_E_eV(396) = 5.580000d+0 ; td_cs_cm2(396) = 2.002200d-25
  td_E_eV(397) = 5.590000d+0 ; td_cs_cm2(397) = 1.988700d-25
  td_E_eV(398) = 5.600000d+0 ; td_cs_cm2(398) = 1.984200d-25
  td_E_eV(399) = 5.610000d+0 ; td_cs_cm2(399) = 1.981800d-25
  td_E_eV(400) = 5.620000d+0 ; td_cs_cm2(400) = 1.976100d-25
  td_E_eV(401) = 5.630000d+0 ; td_cs_cm2(401) = 1.962500d-25
  td_E_eV(402) = 5.640000d+0 ; td_cs_cm2(402) = 1.938600d-25
  td_E_eV(403) = 5.650000d+0 ; td_cs_cm2(403) = 1.902900d-25
  td_E_eV(404) = 5.660000d+0 ; td_cs_cm2(404) = 1.855600d-25
  td_E_eV(405) = 5.670000d+0 ; td_cs_cm2(405) = 1.797800d-25
  td_E_eV(406) = 5.680000d+0 ; td_cs_cm2(406) = 1.731600d-25
  td_E_eV(407) = 5.690000d+0 ; td_cs_cm2(407) = 1.659700d-25
  td_E_eV(408) = 5.700000d+0 ; td_cs_cm2(408) = 1.585700d-25
  td_E_eV(409) = 5.710000d+0 ; td_cs_cm2(409) = 1.513300d-25
  td_E_eV(410) = 5.720000d+0 ; td_cs_cm2(410) = 1.446100d-25
  td_E_eV(411) = 5.730000d+0 ; td_cs_cm2(411) = 1.387400d-25
  td_E_eV(412) = 5.740000d+0 ; td_cs_cm2(412) = 1.339200d-25
  td_E_eV(413) = 5.750000d+0 ; td_cs_cm2(413) = 1.302200d-25
  td_E_eV(414) = 5.760000d+0 ; td_cs_cm2(414) = 1.275500d-25
  td_E_eV(415) = 5.770000d+0 ; td_cs_cm2(415) = 1.257400d-25
  td_E_eV(416) = 5.780000d+0 ; td_cs_cm2(416) = 1.245000d-25
  td_E_eV(417) = 5.790000d+0 ; td_cs_cm2(417) = 1.235600d-25
  td_E_eV(418) = 5.800000d+0 ; td_cs_cm2(418) = 1.226500d-25
  td_E_eV(419) = 5.810000d+0 ; td_cs_cm2(419) = 1.215400d-25
  td_E_eV(420) = 5.820000d+0 ; td_cs_cm2(420) = 1.200800d-25
  td_E_eV(421) = 5.830000d+0 ; td_cs_cm2(421) = 1.181500d-25
  td_E_eV(422) = 5.840000d+0 ; td_cs_cm2(422) = 1.157000d-25
  td_E_eV(423) = 5.850000d+0 ; td_cs_cm2(423) = 1.127500d-25
  td_E_eV(424) = 5.860000d+0 ; td_cs_cm2(424) = 1.093600d-25
  td_E_eV(425) = 5.870000d+0 ; td_cs_cm2(425) = 1.056300d-25
  td_E_eV(426) = 5.880000d+0 ; td_cs_cm2(426) = 1.017200d-25
  td_E_eV(427) = 5.890000d+0 ; td_cs_cm2(427) = 9.782100d-26
  td_E_eV(428) = 5.900000d+0 ; td_cs_cm2(428) = 9.412100d-26
  td_E_eV(429) = 5.910000d+0 ; td_cs_cm2(429) = 9.078300d-26
  td_E_eV(430) = 5.920000d+0 ; td_cs_cm2(430) = 8.791900d-26
  td_E_eV(431) = 5.930000d+0 ; td_cs_cm2(431) = 8.557300d-26
  td_E_eV(432) = 5.940000d+0 ; td_cs_cm2(432) = 8.372100d-26
  td_E_eV(433) = 5.950000d+0 ; td_cs_cm2(433) = 8.228800d-26
  td_E_eV(434) = 5.960000d+0 ; td_cs_cm2(434) = 8.116500d-26
  td_E_eV(435) = 5.970000d+0 ; td_cs_cm2(435) = 8.022900d-26
  td_E_eV(436) = 5.980000d+0 ; td_cs_cm2(436) = 7.936100d-26
  td_E_eV(437) = 5.990000d+0 ; td_cs_cm2(437) = 7.845600d-26
  td_E_eV(438) = 6.000000d+0 ; td_cs_cm2(438) = 7.742600d-26
  td_E_eV(439) = 6.010000d+0 ; td_cs_cm2(439) = 7.620800d-26
  td_E_eV(440) = 6.020000d+0 ; td_cs_cm2(440) = 7.476100d-26
  td_E_eV(441) = 6.030000d+0 ; td_cs_cm2(441) = 7.307400d-26
  td_E_eV(442) = 6.040000d+0 ; td_cs_cm2(442) = 7.116800d-26
  td_E_eV(443) = 6.050000d+0 ; td_cs_cm2(443) = 6.909500d-26
  td_E_eV(444) = 6.060000d+0 ; td_cs_cm2(444) = 6.693000d-26
  td_E_eV(445) = 6.070000d+0 ; td_cs_cm2(445) = 6.476900d-26
  td_E_eV(446) = 6.080000d+0 ; td_cs_cm2(446) = 6.270300d-26
  td_E_eV(447) = 6.090000d+0 ; td_cs_cm2(447) = 6.081200d-26
  td_E_eV(448) = 6.100000d+0 ; td_cs_cm2(448) = 5.914500d-26
  td_E_eV(449) = 6.110000d+0 ; td_cs_cm2(449) = 5.772100d-26
  td_E_eV(450) = 6.120000d+0 ; td_cs_cm2(450) = 5.653200d-26
  td_E_eV(451) = 6.130000d+0 ; td_cs_cm2(451) = 5.554400d-26
  td_E_eV(452) = 6.140000d+0 ; td_cs_cm2(452) = 5.471100d-26
  td_E_eV(453) = 6.150000d+0 ; td_cs_cm2(453) = 5.398000d-26
  td_E_eV(454) = 6.160000d+0 ; td_cs_cm2(454) = 5.329700d-26
  td_E_eV(455) = 6.170000d+0 ; td_cs_cm2(455) = 5.261000d-26
  td_E_eV(456) = 6.180000d+0 ; td_cs_cm2(456) = 5.187300d-26
  td_E_eV(457) = 6.190000d+0 ; td_cs_cm2(457) = 5.104900d-26
  td_E_eV(458) = 6.200000d+0 ; td_cs_cm2(458) = 5.011500d-26
  td_E_eV(459) = 6.210000d+0 ; td_cs_cm2(459) = 4.906400d-26
  td_E_eV(460) = 6.220000d+0 ; td_cs_cm2(460) = 4.790900d-26
  td_E_eV(461) = 6.230000d+0 ; td_cs_cm2(461) = 4.667900d-26
  td_E_eV(462) = 6.240000d+0 ; td_cs_cm2(462) = 4.541800d-26
  td_E_eV(463) = 6.250000d+0 ; td_cs_cm2(463) = 4.417300d-26
  td_E_eV(464) = 6.260000d+0 ; td_cs_cm2(464) = 4.298800d-26
  td_E_eV(465) = 6.270000d+0 ; td_cs_cm2(465) = 4.189800d-26
  td_E_eV(466) = 6.280000d+0 ; td_cs_cm2(466) = 4.092100d-26
  td_E_eV(467) = 6.290000d+0 ; td_cs_cm2(467) = 4.006400d-26
  td_E_eV(468) = 6.300000d+0 ; td_cs_cm2(468) = 3.931900d-26
  td_E_eV(469) = 6.310000d+0 ; td_cs_cm2(469) = 3.867200d-26
  td_E_eV(470) = 6.320000d+0 ; td_cs_cm2(470) = 3.810100d-26
  td_E_eV(471) = 6.330000d+0 ; td_cs_cm2(471) = 3.757900d-26
  td_E_eV(472) = 6.340000d+0 ; td_cs_cm2(472) = 3.707900d-26
  td_E_eV(473) = 6.350000d+0 ; td_cs_cm2(473) = 3.657300d-26
  td_E_eV(474) = 6.360000d+0 ; td_cs_cm2(474) = 3.603700d-26
  td_E_eV(475) = 6.370000d+0 ; td_cs_cm2(475) = 3.545300d-26
  td_E_eV(476) = 6.380000d+0 ; td_cs_cm2(476) = 3.481100d-26
  td_E_eV(477) = 6.390000d+0 ; td_cs_cm2(477) = 3.411400d-26
  td_E_eV(478) = 6.400000d+0 ; td_cs_cm2(478) = 3.337300d-26
  td_E_eV(479) = 6.410000d+0 ; td_cs_cm2(479) = 3.260900d-26
  td_E_eV(480) = 6.420000d+0 ; td_cs_cm2(480) = 3.184500d-26
  td_E_eV(481) = 6.430000d+0 ; td_cs_cm2(481) = 3.110300d-26
  td_E_eV(482) = 6.440000d+0 ; td_cs_cm2(482) = 3.040300d-26
  td_E_eV(483) = 6.450000d+0 ; td_cs_cm2(483) = 2.975600d-26
  td_E_eV(484) = 6.460000d+0 ; td_cs_cm2(484) = 2.916900d-26
  td_E_eV(485) = 6.470000d+0 ; td_cs_cm2(485) = 2.864300d-26
  td_E_eV(486) = 6.480000d+0 ; td_cs_cm2(486) = 2.817300d-26
  td_E_eV(487) = 6.490000d+0 ; td_cs_cm2(487) = 2.775000d-26
  td_E_eV(488) = 6.500000d+0 ; td_cs_cm2(488) = 2.736300d-26
  td_E_eV(489) = 6.510000d+0 ; td_cs_cm2(489) = 2.699600d-26
  td_E_eV(490) = 6.520000d+0 ; td_cs_cm2(490) = 2.663300d-26
  td_E_eV(491) = 6.530000d+0 ; td_cs_cm2(491) = 2.626100d-26
  td_E_eV(492) = 6.540000d+0 ; td_cs_cm2(492) = 2.586700d-26
  td_E_eV(493) = 6.550000d+0 ; td_cs_cm2(493) = 2.544400d-26
  td_E_eV(494) = 6.560000d+0 ; td_cs_cm2(494) = 2.499300d-26
  td_E_eV(495) = 6.570000d+0 ; td_cs_cm2(495) = 2.451900d-26
  td_E_eV(496) = 6.580000d+0 ; td_cs_cm2(496) = 2.403100d-26
  td_E_eV(497) = 6.590000d+0 ; td_cs_cm2(497) = 2.354200d-26
  td_E_eV(498) = 6.600000d+0 ; td_cs_cm2(498) = 2.306300d-26
  td_E_eV(499) = 6.610000d+0 ; td_cs_cm2(499) = 2.260300d-26
  td_E_eV(500) = 6.620000d+0 ; td_cs_cm2(500) = 2.217100d-26
  td_E_eV(501) = 6.630000d+0 ; td_cs_cm2(501) = 2.177200d-26
  td_E_eV(502) = 6.640000d+0 ; td_cs_cm2(502) = 2.140600d-26
  td_E_eV(503) = 6.650000d+0 ; td_cs_cm2(503) = 2.107400d-26
  td_E_eV(504) = 6.660000d+0 ; td_cs_cm2(504) = 2.077100d-26
  td_E_eV(505) = 6.670000d+0 ; td_cs_cm2(505) = 2.049100d-26
  td_E_eV(506) = 6.680000d+0 ; td_cs_cm2(506) = 2.022400d-26
  td_E_eV(507) = 6.690000d+0 ; td_cs_cm2(507) = 1.996100d-26
  td_E_eV(508) = 6.700000d+0 ; td_cs_cm2(508) = 1.969300d-26
  td_E_eV(509) = 6.710000d+0 ; td_cs_cm2(509) = 1.941300d-26
  td_E_eV(510) = 6.720000d+0 ; td_cs_cm2(510) = 1.911700d-26
  td_E_eV(511) = 6.730000d+0 ; td_cs_cm2(511) = 1.880700d-26
  td_E_eV(512) = 6.740000d+0 ; td_cs_cm2(512) = 1.848400d-26
  td_E_eV(513) = 6.750000d+0 ; td_cs_cm2(513) = 1.815500d-26
  td_E_eV(514) = 6.760000d+0 ; td_cs_cm2(514) = 1.782600d-26
  td_E_eV(515) = 6.770000d+0 ; td_cs_cm2(515) = 1.750400d-26
  td_E_eV(516) = 6.780000d+0 ; td_cs_cm2(516) = 1.719300d-26
  td_E_eV(517) = 6.790000d+0 ; td_cs_cm2(517) = 1.690000d-26
  td_E_eV(518) = 6.800000d+0 ; td_cs_cm2(518) = 1.662600d-26
  td_E_eV(519) = 6.810000d+0 ; td_cs_cm2(519) = 1.637400d-26
  td_E_eV(520) = 6.820000d+0 ; td_cs_cm2(520) = 1.614200d-26
  td_E_eV(521) = 6.830000d+0 ; td_cs_cm2(521) = 1.592800d-26
  td_E_eV(522) = 6.840000d+0 ; td_cs_cm2(522) = 1.572600d-26
  td_E_eV(523) = 6.850000d+0 ; td_cs_cm2(523) = 1.553100d-26
  td_E_eV(524) = 6.860000d+0 ; td_cs_cm2(524) = 1.533600d-26
  td_E_eV(525) = 6.870000d+0 ; td_cs_cm2(525) = 1.513700d-26
  td_E_eV(526) = 6.880000d+0 ; td_cs_cm2(526) = 1.492900d-26
  td_E_eV(527) = 6.890000d+0 ; td_cs_cm2(527) = 1.471300d-26
  td_E_eV(528) = 6.900000d+0 ; td_cs_cm2(528) = 1.448800d-26
  td_E_eV(529) = 6.910000d+0 ; td_cs_cm2(529) = 1.425800d-26
  td_E_eV(530) = 6.920000d+0 ; td_cs_cm2(530) = 1.402600d-26
  td_E_eV(531) = 6.930000d+0 ; td_cs_cm2(531) = 1.379600d-26
  td_E_eV(532) = 6.940000d+0 ; td_cs_cm2(532) = 1.357100d-26
  td_E_eV(533) = 6.950000d+0 ; td_cs_cm2(533) = 1.335600d-26
  td_E_eV(534) = 6.960000d+0 ; td_cs_cm2(534) = 1.315400d-26
  td_E_eV(535) = 6.970000d+0 ; td_cs_cm2(535) = 1.296600d-26
  td_E_eV(536) = 6.980000d+0 ; td_cs_cm2(536) = 1.279100d-26
  td_E_eV(537) = 6.990000d+0 ; td_cs_cm2(537) = 1.262900d-26
  td_E_eV(538) = 7.000000d+0 ; td_cs_cm2(538) = 1.247600d-26
  td_E_eV(539) = 7.010000d+0 ; td_cs_cm2(539) = 1.232800d-26
  td_E_eV(540) = 7.020000d+0 ; td_cs_cm2(540) = 1.218100d-26
  td_E_eV(541) = 7.030000d+0 ; td_cs_cm2(541) = 1.203200d-26
  td_E_eV(542) = 7.040000d+0 ; td_cs_cm2(542) = 1.187800d-26
  td_E_eV(543) = 7.050000d+0 ; td_cs_cm2(543) = 1.171800d-26
  td_E_eV(544) = 7.060000d+0 ; td_cs_cm2(544) = 1.155300d-26
  td_E_eV(545) = 7.070000d+0 ; td_cs_cm2(545) = 1.138500d-26
  td_E_eV(546) = 7.080000d+0 ; td_cs_cm2(546) = 1.121600d-26
  td_E_eV(547) = 7.090000d+0 ; td_cs_cm2(547) = 1.104700d-26
  td_E_eV(548) = 7.100000d+0 ; td_cs_cm2(548) = 1.088300d-26
  td_E_eV(549) = 7.110000d+0 ; td_cs_cm2(549) = 1.072500d-26
  td_E_eV(550) = 7.120000d+0 ; td_cs_cm2(550) = 1.057600d-26
  td_E_eV(551) = 7.130000d+0 ; td_cs_cm2(551) = 1.043700d-26
  td_E_eV(552) = 7.140000d+0 ; td_cs_cm2(552) = 1.030600d-26
  td_E_eV(553) = 7.150000d+0 ; td_cs_cm2(553) = 1.018400d-26
  td_E_eV(554) = 7.160000d+0 ; td_cs_cm2(554) = 1.006700d-26
  td_E_eV(555) = 7.170000d+0 ; td_cs_cm2(555) = 0.000000d+0

  td_cs_cm2 = 1.0d4 * td_cs_cm2 ! the original data are in m^2, convert to cm^2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_vib05
!
!-------------------------------------
!
real(8) function CSV_N2_vib05_m3s(E_eV)

  use N2_vib05
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_vib05_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_vib05_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_vib05_m3s










!!!---------------------------------------------------------------------------------------------------- #4
!!! N2, total vibrational cross-section
!!! Table 4 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!! threhold is 0.9 eV, see table 1 of [Majeed and Strickland]
!!
!module N2_vibtot
!  integer, parameter :: N_td=8   ! number of table data points
!! table data
!  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
!  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
!  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
!end module N2_vibtot
!!
!!-------------------------------------
!!
!subroutine Prepare_N2_vibtot
!
!  use N2_vibtot
!  use photoelectrons, only : efactor_eV_to_ms
!
!  implicit none
!
!  integer n
!
!  td_E_eV(1) = 1.3_8; td_cs_cm2(1) = 1.11d-18
!  td_E_eV(2) = 1.6_8; td_cs_cm2(2) = 8.84d-17
!  td_E_eV(3) = 2.0_8; td_cs_cm2(3) = 1.01d-15
!  td_E_eV(4) = 2.4_8; td_cs_cm2(4) = 1.35d-15
!  td_E_eV(5) = 3.0_8; td_cs_cm2(5) = 5.86d-16
!  td_E_eV(6) = 4.0_8; td_cs_cm2(6) = 1.37d-16
!  td_E_eV(7) = 5.0_8; td_cs_cm2(7) = 2.34d-17
!  td_E_eV(8) = 6.0_8; td_cs_cm2(8) = 3.66d-18
!
!  do n = 1, N_td
!     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
!  end do
!
!end subroutine Prepare_N2_vibtot
!!
!!-------------------------------------
!!
!real(8) function CSV_N2_vibtot_m3s(E_eV)
!
!  use N2_vibtot
!  use photoelectrons, only : efactor_eV_to_ms
!
!  implicit none
!
!  real(8) E_eV  ! energy of the electron
!  integer n
!  real(8) left
!
!  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
!     CSV_N2_vibtot_m3s = 0.0_8
!     return
!  end if
!
!  do n = 1, N_td-1
!     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
!        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
!        CSV_N2_vibtot_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
!        return
!     end if
!  end do
!
!end function CSV_N2_vibtot_m3s

!---------------------------------------------------------------------------------------------------- #9  !5
! N2, triplet A^3\Sigma_u^+ cross-section 
! Table 5 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_triplet_A3Zup
  integer, parameter :: N_td=18  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_triplet_A3Zup
!
!-------------------------------------
!
subroutine Prepare_N2_triplet_A3Zup

  use N2_triplet_A3Zup
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   6.5_8; td_cs_cm2(1)  = 1.00d-18 
  td_E_eV(2)  =   7.0_8; td_cs_cm2(2)  = 4.00d-18
  td_E_eV(3)  =   8.0_8; td_cs_cm2(3)  = 7.00d-18
  td_E_eV(4)  =   9.0_8; td_cs_cm2(4)  = 1.00d-17
  td_E_eV(5)  =  10.0_8; td_cs_cm2(5)  = 1.23d-17
  td_E_eV(6)  =  12.0_8; td_cs_cm2(6)  = 1.65d-17
  td_E_eV(7)  =  14.0_8; td_cs_cm2(7)  = 2.00d-17
  td_E_eV(8)  =  16.0_8; td_cs_cm2(8)  = 2.13d-17
  td_E_eV(9)  =  18.0_8; td_cs_cm2(9)  = 2.10d-17
  td_E_eV(10) =  20.0_8; td_cs_cm2(10) = 1.90d-17
  td_E_eV(11) =  24.0_8; td_cs_cm2(11) = 1.40d-17
  td_E_eV(12) =  30.0_8; td_cs_cm2(12) = 9.19d-18
  td_E_eV(13) =  40.0_8; td_cs_cm2(13) = 5.00d-18
  td_E_eV(14) =  50.0_8; td_cs_cm2(14) = 2.62d-18
  td_E_eV(15) =  70.0_8; td_cs_cm2(15) = 9.60d-19
  td_E_eV(16) = 100.0_8; td_cs_cm2(16) = 3.20d-19
  td_E_eV(17) = 150.0_8; td_cs_cm2(17) = 9.00d-20
  td_E_eV(18) = 200.0_8; td_cs_cm2(18) = 4.00d-20

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_triplet_A3Zup
!
!-------------------------------------
!
real(8) function CSV_N2_triplet_A3Zup_m3s(E_eV)

  use N2_triplet_A3Zup
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_triplet_A3Zup_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_triplet_A3Zup_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_triplet_A3Zup_m3s

!---------------------------------------------------------------------------------------------------- #10  !6
! N2, triplet B^3\Pi_g cross-section 
! Table 5 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_triplet_B3Pg
  integer, parameter :: N_td=17  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_triplet_B3Pg
!
!-------------------------------------
!
subroutine Prepare_N2_triplet_B3Pg

  use N2_triplet_B3Pg
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   7.6_8; td_cs_cm2(1)  = 5.30d-19
  td_E_eV(2)  =   8.0_8; td_cs_cm2(2)  = 3.77d-18
  td_E_eV(3)  =   9.0_8; td_cs_cm2(3)  = 1.33d-17
  td_E_eV(4)  =  10.0_8; td_cs_cm2(4)  = 2.19d-17
  td_E_eV(5)  =  12.0_8; td_cs_cm2(5)  = 2.93d-17
  td_E_eV(6)  =  14.0_8; td_cs_cm2(6)  = 2.70d-17
  td_E_eV(7)  =  16.0_8; td_cs_cm2(7)  = 2.16d-17
  td_E_eV(8)  =  18.0_8; td_cs_cm2(8)  = 1.84d-17
  td_E_eV(9)  =  20.0_8; td_cs_cm2(9)  = 1.60d-17
  td_E_eV(10) =  24.0_8; td_cs_cm2(10) = 1.31d-17
  td_E_eV(11) =  30.0_8; td_cs_cm2(11) = 9.73d-18
  td_E_eV(12) =  40.0_8; td_cs_cm2(12) = 5.92d-18
  td_E_eV(13) =  50.0_8; td_cs_cm2(13) = 3.04d-18
  td_E_eV(14) =  70.0_8; td_cs_cm2(14) = 1.25d-18
  td_E_eV(15) = 100.0_8; td_cs_cm2(15) = 4.20d-19
  td_E_eV(16) = 150.0_8; td_cs_cm2(16) = 1.20d-19
  td_E_eV(17) = 200.0_8; td_cs_cm2(17) = 3.00d-20

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_triplet_B3Pg
!
!-------------------------------------
!
real(8) function CSV_N2_triplet_B3Pg_m3s(E_eV)

  use N2_triplet_B3Pg
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_triplet_B3Pg_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_triplet_B3Pg_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_triplet_B3Pg_m3s

!---------------------------------------------------------------------------------------------------- #11  !7
! N2, triplet C^3\Pi_g cross-section 
! Table 5 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_triplet_C3Pg
  integer, parameter :: N_td=14  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_triplet_C3Pg
!
!-------------------------------------
!
subroutine Prepare_N2_triplet_C3Pg

  use N2_triplet_C3Pg
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =  12.0_8; td_cs_cm2(1)  = 5.80d-18
  td_E_eV(2)  =  13.0_8; td_cs_cm2(2)  = 2.53d-17
  td_E_eV(3)  =  14.0_8; td_cs_cm2(3)  = 4.23d-17
  td_E_eV(4)  =  16.0_8; td_cs_cm2(4)  = 2.70d-17
  td_E_eV(5)  =  18.0_8; td_cs_cm2(5)  = 2.07d-17
  td_E_eV(6)  =  20.0_8; td_cs_cm2(6)  = 1.73d-17
  td_E_eV(7)  =  24.0_8; td_cs_cm2(7)  = 1.21d-17
  td_E_eV(8)  =  30.0_8; td_cs_cm2(8)  = 7.70d-18
  td_E_eV(9)  =  40.0_8; td_cs_cm2(9)  = 4.00d-18
  td_E_eV(10) =  50.0_8; td_cs_cm2(10) = 2.50d-18
  td_E_eV(11) =  70.0_8; td_cs_cm2(11) = 1.20d-18
  td_E_eV(12) = 100.0_8; td_cs_cm2(12) = 6.90d-19
  td_E_eV(13) = 150.0_8; td_cs_cm2(13) = 3.10d-19   ! in the table it is -18 but in the corresponding fig.3 of M&S97 it is more like -19
  td_E_eV(14) = 200.0_8; td_cs_cm2(14) = 2.00d-20

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_triplet_C3Pg
!
!-------------------------------------
!
real(8) function CSV_N2_triplet_C3Pg_m3s(E_eV)

  use N2_triplet_C3Pg
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_triplet_C3Pg_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_triplet_C3Pg_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_triplet_C3Pg_m3s

!---------------------------------------------------------------------------------------------------- #12  !8
! N2, triplet W^3\Delta_u cross-section 
! Table 5 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_triplet_W3Du
  integer, parameter :: N_td=16  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_triplet_W3Du
!
!-------------------------------------
!
subroutine Prepare_N2_triplet_W3Du

  use N2_triplet_W3Du
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   8.0_8; td_cs_cm2(1)  = 2.00d-18
  td_E_eV(2)  =   9.0_8; td_cs_cm2(2)  = 7.40d-18
  td_E_eV(3)  =  10.0_8; td_cs_cm2(3)  = 1.20d-17
  td_E_eV(4)  =  12.0_8; td_cs_cm2(4)  = 2.10d-17
  td_E_eV(5)  =  14.0_8; td_cs_cm2(5)  = 3.06d-17
  td_E_eV(6)  =  16.0_8; td_cs_cm2(6)  = 3.73d-17
  td_E_eV(7)  =  18.0_8; td_cs_cm2(7)  = 3.50d-17
  td_E_eV(8)  =  20.0_8; td_cs_cm2(8)  = 2.57d-17
  td_E_eV(9)  =  24.0_8; td_cs_cm2(9)  = 1.57d-17
  td_E_eV(10) =  30.0_8; td_cs_cm2(10) = 9.72d-18
  td_E_eV(11) =  40.0_8; td_cs_cm2(11) = 5.00d-18
  td_E_eV(12) =  50.0_8; td_cs_cm2(12) = 2.62d-18
  td_E_eV(13) =  70.0_8; td_cs_cm2(13) = 9.60d-19
  td_E_eV(14) = 100.0_8; td_cs_cm2(14) = 3.20d-19
  td_E_eV(15) = 150.0_8; td_cs_cm2(15) = 1.00d-19
  td_E_eV(16) = 200.0_8; td_cs_cm2(16) = 4.00d-20

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_triplet_W3Du
!
!-------------------------------------
!
real(8) function CSV_N2_triplet_W3Du_m3s(E_eV)

  use N2_triplet_W3Du
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_triplet_W3Du_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_triplet_W3Du_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_triplet_W3Du_m3s

!---------------------------------------------------------------------------------------------------- #13  !9
! N2, triplet B^prime^3\Sigma_u^- cross-section 
! Table 5 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]! 
!
module N2_triplet_Bp3Zum
  integer, parameter :: N_td=15  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_triplet_Bp3Zum
!
!-------------------------------------
!
subroutine Prepare_N2_triplet_Bp3Zum

  use N2_triplet_Bp3Zum
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   9.0_8; td_cs_cm2(1)  = 1.60d-18
  td_E_eV(2)  =  10.0_8; td_cs_cm2(2)  = 3.50d-18
  td_E_eV(3)  =  12.0_8; td_cs_cm2(3)  = 7.40d-18
  td_E_eV(4)  =  14.0_8; td_cs_cm2(4)  = 1.13d-17
  td_E_eV(5)  =  16.0_8; td_cs_cm2(5)  = 1.14d-17
  td_E_eV(6)  =  18.0_8; td_cs_cm2(6)  = 7.30d-18
  td_E_eV(7)  =  20.0_8; td_cs_cm2(7)  = 5.40d-18
  td_E_eV(8)  =  24.0_8; td_cs_cm2(8)  = 4.30d-18
  td_E_eV(9)  =  30.0_8; td_cs_cm2(9)  = 3.37d-18
  td_E_eV(10) =  40.0_8; td_cs_cm2(10) = 2.45d-18
  td_E_eV(11) =  50.0_8; td_cs_cm2(11) = 1.90d-18
  td_E_eV(12) =  70.0_8; td_cs_cm2(12) = 1.14d-18
  td_E_eV(13) = 100.0_8; td_cs_cm2(13) = 5.30d-19
  td_E_eV(14) = 150.0_8; td_cs_cm2(14) = 1.50d-19
  td_E_eV(15) = 200.0_8; td_cs_cm2(15) = 4.00d-20

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_triplet_Bp3Zum
!
!-------------------------------------
!
real(8) function CSV_N2_triplet_Bp3Zum_m3s(E_eV)

  use N2_triplet_Bp3Zum
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_triplet_Bp3Zum_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_triplet_Bp3Zum_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_triplet_Bp3Zum_m3s

!---------------------------------------------------------------------------------------------------- #14  !10
! N2, singlet a^1\Pi_g cross-section 
! Table 6 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_singlet_a1Pg
  integer, parameter :: N_td=9   ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_singlet_a1Pg
!
!-------------------------------------
!
subroutine Prepare_N2_singlet_a1Pg

  use N2_singlet_a1Pg
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1) =   10.0_8; td_cs_cm2(1) = 2.20d-18
  td_E_eV(2) =   12.0_8; td_cs_cm2(2) = 1.15d-17
  td_E_eV(3) =   18.0_8; td_cs_cm2(3) = 2.69d-17
  td_E_eV(4) =   30.0_8; td_cs_cm2(4) = 1.86d-17
  td_E_eV(5) =   50.0_8; td_cs_cm2(5) = 1.12d-17
  td_E_eV(6) =  100.0_8; td_cs_cm2(6) = 5.59d-18
  td_E_eV(7) =  200.0_8; td_cs_cm2(7) = 2.80d-18
  td_E_eV(8) =  500.0_8; td_cs_cm2(8) = 1.11d-18
  td_E_eV(9) = 1000.0_8; td_cs_cm2(9) = 5.60d-19

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_singlet_a1Pg
!
!-------------------------------------
!
real(8) function CSV_N2_singlet_a1Pg_m3s(E_eV)

  use N2_singlet_a1Pg
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_singlet_a1Pg_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_singlet_a1Pg_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_singlet_a1Pg_m3s

!---------------------------------------------------------------------------------------------------- #15  !11
! N2, singlet b^1\Pi_u cross-section 
! Table 6 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_singlet_b1Pu
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_singlet_b1Pu
!
!-------------------------------------
!
subroutine Prepare_N2_singlet_b1Pu

  use N2_singlet_b1Pu
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   13.0_8; td_cs_cm2(1)  = 2.00d-20
  td_E_eV(2)  =   14.0_8; td_cs_cm2(2)  = 1.57d-18
  td_E_eV(3)  =   16.0_8; td_cs_cm2(3)  = 5.26d-18
  td_E_eV(4)  =   20.0_8; td_cs_cm2(4)  = 1.21d-17
  td_E_eV(5)  =   30.0_8; td_cs_cm2(5)  = 2.08d-17
  td_E_eV(6)  =   50.0_8; td_cs_cm2(6)  = 2.04d-17
  td_E_eV(7)  =  100.0_8; td_cs_cm2(7)  = 1.63d-17
  td_E_eV(8)  =  200.0_8; td_cs_cm2(8)  = 1.16d-17
  td_E_eV(9)  =  500.0_8; td_cs_cm2(9)  = 6.76d-18
  td_E_eV(10) = 1000.0_8; td_cs_cm2(10) = 4.40d-18

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_singlet_b1Pu
!
!-------------------------------------
!
real(8) function CSV_N2_singlet_b1Pu_m3s(E_eV)

  use N2_singlet_b1Pu
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_singlet_b1Pu_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_singlet_b1Pu_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function 

!---------------------------------------------------------------------------------------------------- #16  !12
! N2, singlet b^prime^1\Sigma_u^+ cross-section 
! Table 6 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_singlet_bp1Zup
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_singlet_bp1Zup
!
!-------------------------------------
!
subroutine Prepare_N2_singlet_bp1Zup

  use N2_singlet_bp1Zup
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   15.0_8; td_cs_cm2(1)  = 1.30d-19
  td_E_eV(2)  =   16.0_8; td_cs_cm2(2)  = 8.90d-19
  td_E_eV(3)  =   18.0_8; td_cs_cm2(3)  = 2.56d-18
  td_E_eV(4)  =   20.0_8; td_cs_cm2(4)  = 4.36d-18
  td_E_eV(5)  =   30.0_8; td_cs_cm2(5)  = 1.01d-17
  td_E_eV(6)  =   50.0_8; td_cs_cm2(6)  = 1.32d-17
  td_E_eV(7)  =  100.0_8; td_cs_cm2(7)  = 1.26d-17
  td_E_eV(8)  =  200.0_8; td_cs_cm2(8)  = 9.78d-18
  td_E_eV(9)  =  500.0_8; td_cs_cm2(9)  = 6.23d-18
  td_E_eV(10) = 1000.0_8; td_cs_cm2(10) = 4.09d-18

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_singlet_bp1Zup
!
!-------------------------------------
!
real(8) function CSV_N2_singlet_bp1Zup_m3s(E_eV)

  use N2_singlet_bp1Zup
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_singlet_bp1Zup_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_singlet_bp1Zup_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_singlet_bp1Zup_m3s

!---------------------------------------------------------------------------------------------------- #17  !13
! N2, singlet c_4^prime^1\Sigma_u cross-section 
! Table 6 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_singlet_c4p1Zu
  integer, parameter :: N_td=9   ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_singlet_c4p1Zu
!
!-------------------------------------
!
subroutine Prepare_N2_singlet_c4p1Zu

  use N2_singlet_c4p1Zu
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   14.0_8; td_cs_cm2(1) = 1.20d-19
  td_E_eV(2)  =   16.0_8; td_cs_cm2(2) = 8.10d-19
  td_E_eV(3)  =   20.0_8; td_cs_cm2(3) = 2.98d-18
  td_E_eV(4)  =   30.0_8; td_cs_cm2(4) = 7.19d-18
  td_E_eV(5)  =   50.0_8; td_cs_cm2(5) = 1.10d-17
  td_E_eV(6)  =  100.0_8; td_cs_cm2(6) = 1.20d-17
  td_E_eV(7)  =  200.0_8; td_cs_cm2(7) = 9.95d-18
  td_E_eV(8)  =  500.0_8; td_cs_cm2(8) = 6.29d-18
  td_E_eV(9)  = 1000.0_8; td_cs_cm2(9) = 3.94d-18

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_singlet_c4p1Zu
!
!-------------------------------------
!
real(8) function CSV_N2_singlet_c4p1Zu_m3s(E_eV)

  use N2_singlet_c4p1Zu
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_singlet_c4p1Zu_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_singlet_c4p1Zu_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_singlet_c4p1Zu_m3s

!---------------------------------------------------------------------------------------------------- #18  !14
! N2, singlet w^1\Delta_u cross-section 
! Table 6 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_singlet_w1Du
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_singlet_w1Du
!
!-------------------------------------
!
subroutine Prepare_N2_singlet_w1Du

  use N2_singlet_w1Du
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   9.0_8; td_cs_cm2(1)  = 1.80d-19
  td_E_eV(2)  =  10.0_8; td_cs_cm2(2)  = 3.62d-18
  td_E_eV(3)  =  12.0_8; td_cs_cm2(3)  = 9.81d-18
  td_E_eV(4)  =  14.0_8; td_cs_cm2(4)  = 1.15d-17
  td_E_eV(5)  =  16.0_8; td_cs_cm2(5)  = 8.06d-18
  td_E_eV(6)  =  20.0_8; td_cs_cm2(6)  = 4.30d-18
  td_E_eV(7)  =  30.0_8; td_cs_cm2(7)  = 2.31d-18
  td_E_eV(8)  =  50.0_8; td_cs_cm2(8)  = 7.10d-19
  td_E_eV(9)  = 100.0_8; td_cs_cm2(9)  = 1.30d-19
  td_E_eV(10) = 200.0_8; td_cs_cm2(10) = 1.00d-20

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_singlet_w1Du
!
!-------------------------------------
!
real(8) function CSV_N2_singlet_w1Du_m3s(E_eV)

  use N2_singlet_w1Du
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_singlet_w1Du_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_singlet_w1Du_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_singlet_w1Du_m3s

!---------------------------------------------------------------------------------------------------- #19  !15
! N2, singlet c^1\Pi_u cross-section 
! Table 7 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_singlet_c1Pu
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_singlet_c1Pu
!
!-------------------------------------
!
subroutine Prepare_N2_singlet_c1Pu

  use N2_singlet_c1Pu
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   13.0_8; td_cs_cm2(1)  = 2.00d-20
  td_E_eV(2)  =   14.0_8; td_cs_cm2(2)  = 1.63d-18
  td_E_eV(3)  =   16.0_8; td_cs_cm2(3)  = 5.49d-18
  td_E_eV(4)  =   20.0_8; td_cs_cm2(4)  = 1.26d-17
  td_E_eV(5)  =   30.0_8; td_cs_cm2(5)  = 2.17d-17
  td_E_eV(6)  =   50.0_8; td_cs_cm2(6)  = 2.13d-17
  td_E_eV(7)  =  100.0_8; td_cs_cm2(7)  = 1.70d-17
  td_E_eV(8)  =  200.0_8; td_cs_cm2(8)  = 1.21d-17
  td_E_eV(9)  =  500.0_8; td_cs_cm2(9)  = 7.05d-18
  td_E_eV(10) = 1000.0_8; td_cs_cm2(10) = 4.58d-18

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_singlet_c1Pu
!
!-------------------------------------
!
real(8) function CSV_N2_singlet_c1Pu_m3s(E_eV)

  use N2_singlet_c1Pu
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_singlet_c1Pu_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_singlet_c1Pu_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_singlet_c1Pu_m3s

!---------------------------------------------------------------------------------------------------- #20  !16
! N2, singlet a^prime^1\Sigma_u^- cross-section 
! Table 7 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_singlet_ap1Zum
  integer, parameter :: N_td=11  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_singlet_ap1Zum
!
!-------------------------------------
!
subroutine Prepare_N2_singlet_ap1Zum

  use N2_singlet_ap1Zum
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   11.0_8; td_cs_cm2(1)  = 4.80d-18
  td_E_eV(2)  =   12.0_8; td_cs_cm2(2)  = 6.52d-18
  td_E_eV(3)  =   14.0_8; td_cs_cm2(3)  = 9.50d-18
  td_E_eV(4)  =   16.0_8; td_cs_cm2(4)  = 8.50d-18
  td_E_eV(5)  =   20.0_8; td_cs_cm2(5)  = 4.80d-18
  td_E_eV(6)  =   30.0_8; td_cs_cm2(6)  = 2.30d-18
  td_E_eV(7)  =   50.0_8; td_cs_cm2(7)  = 1.93d-18
  td_E_eV(8)  =  100.0_8; td_cs_cm2(8)  = 9.60d-19
  td_E_eV(9)  =  200.0_8; td_cs_cm2(9)  = 4.50d-19
  td_E_eV(10) =  500.0_8; td_cs_cm2(10) = 1.60d-19
  td_E_eV(11) = 1000.0_8; td_cs_cm2(11) = 8.00d-20

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_singlet_ap1Zum
!
!-------------------------------------
!
real(8) function CSV_N2_singlet_ap1Zum_m3s(E_eV)

  use N2_singlet_ap1Zum
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_singlet_ap1Zum_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_singlet_ap1Zum_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_singlet_ap1Zum_m3s

!---------------------------------------------------------------------------------------------------- #21  !17
! N2, singlet a^prime^prime^1\Sigma_g^+ cross-section 
! Table 7 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_singlet_app1Zgp
  integer, parameter :: N_td=11  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_singlet_app1Zgp
!
!-------------------------------------
!
subroutine Prepare_N2_singlet_app1Zgp

  use N2_singlet_app1Zgp
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   13.0_8; td_cs_cm2(1)  = 3.50d-19
  td_E_eV(2)  =   14.0_8; td_cs_cm2(2)  = 1.82d-18
  td_E_eV(3)  =   16.0_8; td_cs_cm2(3)  = 3.68d-18
  td_E_eV(4)  =   20.0_8; td_cs_cm2(4)  = 5.51d-18
  td_E_eV(5)  =   25.0_8; td_cs_cm2(5)  = 4.26d-18
  td_E_eV(6)  =   30.0_8; td_cs_cm2(6)  = 3.04d-18
  td_E_eV(7)  =   50.0_8; td_cs_cm2(7)  = 1.52d-18
  td_E_eV(8)  =  100.0_8; td_cs_cm2(8)  = 6.90d-19
  td_E_eV(9)  =  200.0_8; td_cs_cm2(9)  = 3.30d-19
  td_E_eV(10) =  500.0_8; td_cs_cm2(10) = 1.30d-19
  td_E_eV(11) = 1000.0_8; td_cs_cm2(11) = 7.00d-20

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_singlet_app1Zgp
!
!-------------------------------------
!
real(8) function CSV_N2_singlet_app1Zgp_m3s(E_eV)

  use N2_singlet_app1Zgp
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_singlet_app1Zgp_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_singlet_app1Zgp_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_singlet_app1Zgp_m3s

!---------------------------------------------------------------------------------------------------- #22  !18
! N2, singlet other ^1\Pi_u cross-section 
! Table 7 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_singlet_other1Pu
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_singlet_other1Pu
!
!-------------------------------------
!
subroutine Prepare_N2_singlet_other1Pu

  use N2_singlet_other1Pu
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   13.0_8; td_cs_cm2(1)  = 2.00d-20
  td_E_eV(2)  =   14.0_8; td_cs_cm2(2)  = 2.16d-18
  td_E_eV(3)  =   16.0_8; td_cs_cm2(3)  = 7.24d-18
  td_E_eV(4)  =   20.0_8; td_cs_cm2(4)  = 1.66d-17
  td_E_eV(5)  =   30.0_8; td_cs_cm2(5)  = 2.86d-17
  td_E_eV(6)  =   50.0_8; td_cs_cm2(6)  = 2.81d-17
  td_E_eV(7)  =  100.0_8; td_cs_cm2(7)  = 2.25d-17
  td_E_eV(8)  =  200.0_8; td_cs_cm2(8)  = 1.59d-17
  td_E_eV(9)  =  500.0_8; td_cs_cm2(9)  = 9.30d-18
  td_E_eV(10) = 1000.0_8; td_cs_cm2(10) = 6.05d-18

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_singlet_other1Pu
!
!-------------------------------------
!
real(8) function CSV_N2_singlet_other1Pu_m3s(E_eV)

  use N2_singlet_other1Pu
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_singlet_other1Pu_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_singlet_other1Pu_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_singlet_other1Pu_m3s

!---------------------------------------------------------------------------------------------------- #23  !19
! N2, high lying state 15.8 eV peak cross-section 
! Table 8 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_hls_15p8eV
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_hls_15p8eV
!
!-------------------------------------
!
subroutine Prepare_N2_hls_15p8eV

  use N2_hls_15p8eV
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   17.0_8; td_cs_cm2(1)  = 2.25d-18
  td_E_eV(2)  =   18.0_8; td_cs_cm2(2)  = 4.53d-18
  td_E_eV(3)  =   20.0_8; td_cs_cm2(3)  = 8.52d-18
  td_E_eV(4)  =   25.0_8; td_cs_cm2(4)  = 1.74d-17
  td_E_eV(5)  =   30.0_8; td_cs_cm2(5)  = 2.28d-17
  td_E_eV(6)  =   50.0_8; td_cs_cm2(6)  = 2.39d-17
  td_E_eV(7)  =  100.0_8; td_cs_cm2(7)  = 1.92d-17
  td_E_eV(8)  =  200.0_8; td_cs_cm2(8)  = 1.40d-17
  td_E_eV(9)  =  500.0_8; td_cs_cm2(9)  = 8.00d-18
  td_E_eV(10) = 1000.0_8; td_cs_cm2(10) = 5.10d-18

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_hls_15p8eV
!
!-------------------------------------
!
real(8) function CSV_N2_hls_15p8eV_m3s(E_eV)

  use N2_hls_15p8eV
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_hls_15p8eV_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_hls_15p8eV_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_hls_15p8eV_m3s

!---------------------------------------------------------------------------------------------------- #24  !20
! N2, high lying state VUV cross-section 
! Table 8 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_hls_VUV
  integer, parameter :: N_td=9  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_hls_VUV
!
!-------------------------------------
!
subroutine Prepare_N2_hls_VUV

  use N2_hls_VUV
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   25.0_8; td_cs_cm2(1)  = 1.61d-18
  td_E_eV(2)  =   30.0_8; td_cs_cm2(2)  = 3.28d-18
  td_E_eV(3)  =   40.0_8; td_cs_cm2(3)  = 5.29d-18
  td_E_eV(4)  =   50.0_8; td_cs_cm2(4)  = 8.97d-18
  td_E_eV(5)  =   70.0_8; td_cs_cm2(5)  = 1.51d-17
  td_E_eV(6)  =  100.0_8; td_cs_cm2(6)  = 1.60d-17
  td_E_eV(7)  =  200.0_8; td_cs_cm2(7)  = 1.24d-17
  td_E_eV(8)  =  500.0_8; td_cs_cm2(8)  = 7.12d-18
  td_E_eV(9)  = 1000.0_8; td_cs_cm2(9)  = 4.56d-18

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_hls_VUV
!
!-------------------------------------
!
real(8) function CSV_N2_hls_VUV_m3s(E_eV)

  use N2_hls_VUV
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_hls_VUV_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_hls_VUV_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_hls_VUV_m3s

!---------------------------------------------------------------------------------------------------- #25  !21
! N2, high lying state 17.3 eV peak cross-section 
! Table 8 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_hls_17p3eV
  integer, parameter :: N_td=9  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_hls_17p3eV
!
!-------------------------------------
!
subroutine Prepare_N2_hls_17p3eV

  use N2_hls_17p3eV
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1) =   18.0_8; td_cs_cm2(1) = 1.20d-18
  td_E_eV(2) =   20.0_8; td_cs_cm2(2) = 2.70d-18
  td_E_eV(3) =   25.0_8; td_cs_cm2(3) = 7.20d-18
  td_E_eV(4) =   30.0_8; td_cs_cm2(4) = 9.70d-18
  td_E_eV(5) =   50.0_8; td_cs_cm2(5) = 1.02d-17
  td_E_eV(6) =  100.0_8; td_cs_cm2(6) = 8.50d-18
  td_E_eV(7) =  200.0_8; td_cs_cm2(7) = 6.26d-18
  td_E_eV(8) =  500.0_8; td_cs_cm2(8) = 3.77d-18
  td_E_eV(9) = 1000.0_8; td_cs_cm2(9) = 2.32d-18

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_hls_17p3eV
!
!-------------------------------------
!
real(8) function CSV_N2_hls_17p3eV_m3s(E_eV)

  use N2_hls_17p3eV
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_hls_17p3eV_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_hls_17p3eV_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_hls_17p3eV_m3s

!---------------------------------------------------------------------------------------------------- #26  !22
! N2, high lying state Rydberg cross-section 
! Table 8 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_hls_Rydberg
  integer, parameter :: N_td=6   ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_hls_Rydberg
!
!-------------------------------------
!
subroutine Prepare_N2_hls_Rydberg

  use N2_hls_Rydberg
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1) =   50.0_8; td_cs_cm2(1) = 1.86d-18
  td_E_eV(2) =   70.0_8; td_cs_cm2(2) = 3.27d-18
  td_E_eV(3) =  100.0_8; td_cs_cm2(3) = 3.48d-18
  td_E_eV(4) =  200.0_8; td_cs_cm2(4) = 2.75d-18
  td_E_eV(5) =  500.0_8; td_cs_cm2(5) = 1.44d-18
  td_E_eV(6) = 1000.0_8; td_cs_cm2(6) = 6.10d-19

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_hls_Rydberg
!
!-------------------------------------
!
real(8) function CSV_N2_hls_Rydberg_m3s(E_eV)

  use N2_hls_Rydberg
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_hls_Rydberg_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_hls_Rydberg_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function 

!---------------------------------------------------------------------------------------------------- #27  !23
! N2, high lying state triplet manifold cross-section 
! Table 8 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module N2_hls_tripman
  integer, parameter :: N_td=12  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_hls_tripman
!
!-------------------------------------
!
subroutine Prepare_N2_hls_tripman

  use N2_hls_tripman
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   12.0_8; td_cs_cm2(1)  = 1.23d-18
  td_E_eV(2)  =   13.0_8; td_cs_cm2(2)  = 5.19d-18
  td_E_eV(3)  =   14.0_8; td_cs_cm2(3)  = 1.03d-17
  td_E_eV(4)  =   16.0_8; td_cs_cm2(4)  = 1.31d-17
  td_E_eV(5)  =   18.0_8; td_cs_cm2(5)  = 9.86d-18
  td_E_eV(6)  =   20.0_8; td_cs_cm2(6)  = 6.45d-18
  td_E_eV(7)  =   25.0_8; td_cs_cm2(7)  = 3.13d-18
  td_E_eV(8)  =   30.0_8; td_cs_cm2(8)  = 1.92d-18
  td_E_eV(9)  =   40.0_8; td_cs_cm2(9)  = 1.03d-18
  td_E_eV(10) =   50.0_8; td_cs_cm2(10) = 6.60d-19
  td_E_eV(11) =   70.0_8; td_cs_cm2(11) = 2.70d-19
  td_E_eV(12) = 1000.0_8; td_cs_cm2(12) = 7.00d-20

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_hls_tripman
!
!-------------------------------------
!
real(8) function CSV_N2_hls_tripman_m3s(E_eV)

  use N2_hls_tripman
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_hls_tripman_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_hls_tripman_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_hls_tripman_m3s

!---------------------------------------------------------------------------------------------------- #28  !24
!
! from file cross_sections_check/output/data-24/cross-section-O2-elastic-01.dat
!
! A. Schmalzried, A. Luque and N. Lehtinen,  IAA Database on lxcat, www.lxcat.net/IAA, August 2023,
!                  Instituto de Astrofísica de Andalucía.
!
module O2_elastic
!  integer, parameter :: N_td=30  ! number of table data points
  integer, parameter :: N_td=41  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O2_elastic
!
!-------------------------------------
!
subroutine Prepare_O2_elastic

  use O2_elastic
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV( 1) = 1.000000d-03 ;   td_cs_cm2( 1) = 5.671630d-21   !###   here in m^2
  td_E_eV( 2) = 1.000000d-02 ;   td_cs_cm2( 2) = 1.179840d-20
  td_E_eV( 3) = 1.000000d-01 ;   td_cs_cm2( 3) = 2.172730d-20
  td_E_eV( 4) = 1.000000d+00 ;   td_cs_cm2( 4) = 6.114150d-20
  td_E_eV( 5) = 2.000000d+00 ;   td_cs_cm2( 5) = 6.882130d-20
  td_E_eV( 6) = 3.000000d+00 ;   td_cs_cm2( 6) = 6.849890d-20
  td_E_eV( 7) = 4.000000d+00 ;   td_cs_cm2( 7) = 7.020770d-20
  td_E_eV( 8) = 5.000000d+00 ;   td_cs_cm2( 8) = 6.972390d-20
  td_E_eV( 9) = 7.000000d+00 ;   td_cs_cm2( 9) = 7.504820d-20
  td_E_eV(10) = 8.000000d+00 ;   td_cs_cm2(10) = 7.922820d-20
  td_E_eV(11) = 9.000000d+00 ;   td_cs_cm2(11) = 7.998070d-20
  td_E_eV(12) = 1.000000d+01 ;   td_cs_cm2(12) = 8.835780d-20
  td_E_eV(13) = 1.500000d+01 ;   td_cs_cm2(13) = 9.230460d-20
  td_E_eV(14) = 2.000000d+01 ;   td_cs_cm2(14) = 9.153310d-20
  td_E_eV(15) = 3.000000d+01 ;   td_cs_cm2(15) = 9.720510d-20
  td_E_eV(16) = 3.500000d+01 ;   td_cs_cm2(16) = 9.038340d-20
  td_E_eV(17) = 4.000000d+01 ;   td_cs_cm2(17) = 8.498410d-20
  td_E_eV(18) = 4.500000d+01 ;   td_cs_cm2(18) = 8.038850d-20
  td_E_eV(19) = 5.000000d+01 ;   td_cs_cm2(19) = 7.629050d-20
  td_E_eV(20) = 5.500000d+01 ;   td_cs_cm2(20) = 7.266870d-20
  td_E_eV(21) = 6.000000d+01 ;   td_cs_cm2(21) = 6.953450d-20
  td_E_eV(22) = 6.500000d+01 ;   td_cs_cm2(22) = 6.685070d-20
  td_E_eV(23) = 7.000000d+01 ;   td_cs_cm2(23) = 6.456950d-20
  td_E_eV(24) = 7.500000d+01 ;   td_cs_cm2(24) = 6.261970d-20
  td_E_eV(25) = 8.000000d+01 ;   td_cs_cm2(25) = 6.092450d-20
  td_E_eV(26) = 8.500000d+01 ;   td_cs_cm2(26) = 5.938580d-20
  td_E_eV(27) = 9.000000d+01 ;   td_cs_cm2(27) = 5.794050d-20
  td_E_eV(28) = 9.500000d+01 ;   td_cs_cm2(28) = 5.653660d-20
  td_E_eV(29) = 1.000000d+02 ;   td_cs_cm2(29) = 5.514560d-20
  td_E_eV(30) = 1.250000d+02 ;   td_cs_cm2(30) = 4.859740d-20
  td_E_eV(31) = 1.500000d+02 ;   td_cs_cm2(31) = 4.360340d-20
  td_E_eV(32) = 2.000000d+02 ;   td_cs_cm2(32) = 3.897740d-20
  td_E_eV(33) = 2.500000d+02 ;   td_cs_cm2(33) = 3.525910d-20
  td_E_eV(34) = 3.000000d+02 ;   td_cs_cm2(34) = 3.189290d-20
  td_E_eV(35) = 4.000000d+02 ;   td_cs_cm2(35) = 2.636590d-20
  td_E_eV(36) = 5.000000d+02 ;   td_cs_cm2(36) = 2.268070d-20
  td_E_eV(37) = 6.000000d+02 ;   td_cs_cm2(37) = 2.019160d-20
  td_E_eV(38) = 7.000000d+02 ;   td_cs_cm2(38) = 1.788780d-20
  td_E_eV(39) = 8.000000d+02 ;   td_cs_cm2(39) = 1.621740d-20
  td_E_eV(40) = 9.000000d+02 ;   td_cs_cm2(40) = 1.484520d-20
  td_E_eV(41) = 1.000000d+03 ;   td_cs_cm2(41) = 1.369500d-20

  td_cs_cm2 = 1.0d4 * td_cs_cm2  ! the original data are in m^2, convert to cm^2

! O2, elastic scattering cross-section 
! Table 2 of [Y.Itikawa, J.Phys.Chem.Ref.Data, v.38, p.1-20, 2009] 
!  td_E_eV(1)   =    1.0_8  ; td_cs_cm2(1)   = 5.97d-16
!  td_E_eV(2)   =    2.0_8  ; td_cs_cm2(2)   = 6.45d-16
!  td_E_eV(3)   =    3.0_8  ; td_cs_cm2(3)   = 6.74d-16
!  td_E_eV(4)   =    4.0_8  ; td_cs_cm2(4)   = 6.93d-16
!  td_E_eV(5)   =    5.0_8  ; td_cs_cm2(5)   = 7.20d-16
!  td_E_eV(6)   =    6.0_8  ; td_cs_cm2(6)   = 7.52d-16
!  td_E_eV(7)   =    7.0_8  ; td_cs_cm2(7)   = 7.86d-16
!  td_E_eV(8)   =    8.0_8  ; td_cs_cm2(8)   = 8.21d-16
!  td_E_eV(9)   =    9.0_8  ; td_cs_cm2(9)   = 8.49d-16
!  td_E_eV(10)  =   10.0_8  ; td_cs_cm2(10)  = 8.80d-16
!  td_E_eV(11)  =   12.0_8  ; td_cs_cm2(11)  = 9.00d-16
!  td_E_eV(12)  =   15.0_8  ; td_cs_cm2(12)  = 8.89d-16
!  td_E_eV(13)  =   20.0_8  ; td_cs_cm2(13)  = 8.60d-16
!  td_E_eV(14)  =   30.0_8  ; td_cs_cm2(14)  = 8.09d-16
!  td_E_eV(15)  =   40.0_8  ; td_cs_cm2(15)  = 7.30d-16
!  td_E_eV(16)  =   50.0_8  ; td_cs_cm2(16)  = 6.59d-16
!  td_E_eV(17)  =   60.0_8  ; td_cs_cm2(17)  = 6.08d-16
!  td_E_eV(18)  =   70.0_8  ; td_cs_cm2(18)  = 5.63d-16
!  td_E_eV(19)  =   80.0_8  ; td_cs_cm2(19)  = 5.29d-16
!  td_E_eV(20)  =   90.0_8  ; td_cs_cm2(20)  = 5.01d-16
!  td_E_eV(21)  =  100.0_8  ; td_cs_cm2(21)  = 4.78d-16
!  td_E_eV(22)  =  200.0_8  ; td_cs_cm2(22)  = 3.15d-16
!  td_E_eV(23)  =  300.0_8  ; td_cs_cm2(23)  = 2.40d-16
!  td_E_eV(24)  =  400.0_8  ; td_cs_cm2(24)  = 2.00d-16
!  td_E_eV(25)  =  500.0_8  ; td_cs_cm2(25)  = 1.72d-16
!  td_E_eV(26)  =  600.0_8  ; td_cs_cm2(26)  = 1.53d-16
!  td_E_eV(27)  =  700.0_8  ; td_cs_cm2(27)  = 1.37d-16
!  td_E_eV(28)  =  800.0_8  ; td_cs_cm2(28)  = 1.27d-16
!  td_E_eV(29)  =  900.0_8  ; td_cs_cm2(29)  = 1.18d-16
!  td_E_eV(30)  = 1000.0_8  ; td_cs_cm2(30)  = 1.10d-16
  
  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O2_elastic
!
!-------------------------------------
! above 1000 eV the cross section is calculated using
! an interplation formula derived upon assuming that in log-log scale the 
! cross section decays with energy as a straight line 
! continuing the tabulated data
! (see Figure 9 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
real(8) function CSV_O2_elastic_m3s(E_eV)

  use O2_elastic
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  real(8) CS_O2_elastic_cm2

  if (E_eV.lt.td_E_eV(1)) then
     CSV_O2_elastic_m3s = 1.0d-4 * td_cs_cm2(1) * sqrt(E_eV) * efactor_eV_to_ms
     return
  end if

  if (E_eV.GT.td_E_eV(N_td)) then
     CS_O2_elastic_cm2 = (1.0d-15 * log(E_eV) + 1.889d-14) / E_eV**0.79
     CSV_O2_elastic_m3s = max(0.0_8, 1.0d-4 * CS_O2_elastic_cm2 * sqrt(E_eV) * efactor_eV_to_ms)
     RETURN
  END IF

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O2_elastic_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O2_elastic_m3s

!---------------------------------------------------------------------------------------------------- #29  !25
! O2, {O2->O2+} ionization cross-section 
! Table 11 of [Y.Itikawa, J.Phys.Chem.Ref.Data, v.38, p.1-20, 2009] 
!
module O2_ionO2
  integer, parameter :: N_td=43  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O2_ionO2
!
!-------------------------------------
!
subroutine Prepare_O2_ionO2

  use O2_ionO2
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =  13.0_8; td_cs_cm2(1)  = 1.17d-18
  td_E_eV(2)  =  15.5_8; td_cs_cm2(2)  = 7.30d-18
  td_E_eV(3)  =  18.0_8; td_cs_cm2(3)  = 1.64d-17
  td_E_eV(4)  =  23.0_8; td_cs_cm2(4)  = 3.66d-17
  td_E_eV(5)  =  28.0_8; td_cs_cm2(5)  = 5.63d-17
  td_E_eV(6)  =  33.0_8; td_cs_cm2(6)  = 7.58d-17
  td_E_eV(7)  =  38.0_8; td_cs_cm2(7)  = 9.29d-17
  td_E_eV(8)  =  43.0_8; td_cs_cm2(8)  = 1.08d-16
  td_E_eV(9)  =  48.0_8; td_cs_cm2(9)  = 1.19d-16
  td_E_eV(10) =  53.0_8; td_cs_cm2(10) = 1.29d-16
  td_E_eV(11) =  58.0_8; td_cs_cm2(11) = 1.36d-16
  td_E_eV(12) =  63.5_8; td_cs_cm2(12) = 1.42d-16
  td_E_eV(13) =  68.0_8; td_cs_cm2(13) = 1.47d-16
  td_E_eV(14) =  73.0_8; td_cs_cm2(14) = 1.50d-16
  td_E_eV(15) =  78.0_8; td_cs_cm2(15) = 1.51d-16
  td_E_eV(16) =  83.0_8; td_cs_cm2(16) = 1.53d-16
  td_E_eV(17) =  88.0_8; td_cs_cm2(17) = 1.55d-16
  td_E_eV(18) =  93.0_8; td_cs_cm2(18) = 1.56d-16
  td_E_eV(19) =  98.0_8; td_cs_cm2(19) = 1.56d-16
  td_E_eV(20) = 108.0_8; td_cs_cm2(20) = 1.54d-16
  td_E_eV(21) = 118.0_8; td_cs_cm2(21) = 1.53d-16
  td_E_eV(22) = 138.0_8; td_cs_cm2(22) = 1.50d-16
  td_E_eV(23) = 158.0_8; td_cs_cm2(23) = 1.48d-16
  td_E_eV(24) = 178.0_8; td_cs_cm2(24) = 1.43d-16
  td_E_eV(25) = 198.0_8; td_cs_cm2(25) = 1.39d-16
  td_E_eV(26) = 223.0_8; td_cs_cm2(26) = 1.34d-16
  td_E_eV(27) = 248.0_8; td_cs_cm2(27) = 1.31d-16
  td_E_eV(28) = 273.0_8; td_cs_cm2(28) = 1.24d-16
  td_E_eV(29) = 298.0_8; td_cs_cm2(29) = 1.20d-16
  td_E_eV(30) = 348.0_8; td_cs_cm2(30) = 1.13d-16
  td_E_eV(31) = 398.0_8; td_cs_cm2(31) = 1.05d-16
  td_E_eV(32) = 448.0_8; td_cs_cm2(32) = 9.83d-17
  td_E_eV(33) = 498.0_8; td_cs_cm2(33) = 9.23d-17
  td_E_eV(34) = 548.0_8; td_cs_cm2(34) = 8.82d-17
  td_E_eV(35) = 598.0_8; td_cs_cm2(35) = 8.27d-17
  td_E_eV(36) = 648.0_8; td_cs_cm2(36) = 8.00d-17
  td_E_eV(37) = 698.0_8; td_cs_cm2(37) = 7.61d-17
  td_E_eV(38) = 748.0_8; td_cs_cm2(38) = 7.20d-17
  td_E_eV(39) = 798.0_8; td_cs_cm2(39) = 6.86d-17
  td_E_eV(40) = 848.0_8; td_cs_cm2(40) = 6.71d-17
  td_E_eV(41) = 898.0_8; td_cs_cm2(41) = 6.43d-17
  td_E_eV(42) = 948.0_8; td_cs_cm2(42) = 6.17d-17
  td_E_eV(43) = 998.0_8; td_cs_cm2(43) = 5.97d-17

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O2_ionO2
!
!-------------------------------------
!
real(8) function CSV_O2_ionO2_m3s(E_eV)

  use O2_ionO2
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  real(8) CS_O2_ionO2_cm2

  if (E_eV.lt.td_E_eV(1)) then
     CSV_O2_ionO2_m3s = 0.0_8
     return
  end if

  if (E_eV.gt.td_E_eV(N_td)) then
     CS_O2_ionO2_cm2 = (1.94d-14 * log(E_eV) - 7.41d-14) / E_eV
     CSV_O2_ionO2_m3s = max(0.0_8, 1.0d-4 * CS_O2_ionO2_cm2 * sqrt(E_eV) * efactor_eV_to_ms)
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O2_ionO2_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O2_ionO2_m3s

!---------------------------------------------------------------------------------------------------- #30  !26
! O2, {O2->O+} ionization cross-section 
! Table 11 of [Y.Itikawa, J.Phys.Chem.Ref.Data, v.38, p.1-20, 2009] 
!
module O2_ionO
  integer, parameter :: N_td=40  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O2_ionO
!
!-------------------------------------
!
subroutine Prepare_O2_ionO

  use O2_ionO
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =  23.0_8; td_cs_cm2(1)  = 1.67d-18
  td_E_eV(2)  =  28.0_8; td_cs_cm2(2)  = 7.81d-18
  td_E_eV(3)  =  33.0_8; td_cs_cm2(3)  = 1.69d-17
  td_E_eV(4)  =  38.0_8; td_cs_cm2(4)  = 2.58d-17
  td_E_eV(5)  =  43.0_8; td_cs_cm2(5)  = 3.33d-17
  td_E_eV(6)  =  48.0_8; td_cs_cm2(6)  = 4.19d-17
  td_E_eV(7)  =  53.0_8; td_cs_cm2(7)  = 4.90d-17
  td_E_eV(8)  =  58.0_8; td_cs_cm2(8)  = 5.53d-17
  td_E_eV(9)  =  63.5_8; td_cs_cm2(9)  = 6.21d-17
  td_E_eV(10) =  68.0_8; td_cs_cm2(10) = 6.79d-17
  td_E_eV(11) =  73.0_8; td_cs_cm2(11) = 7.17d-17
  td_E_eV(12) =  78.0_8; td_cs_cm2(12) = 7.51d-17
  td_E_eV(13) =  83.0_8; td_cs_cm2(13) = 8.01d-17
  td_E_eV(14) =  88.0_8; td_cs_cm2(14) = 8.27d-17
  td_E_eV(15) =  93.0_8; td_cs_cm2(15) = 8.55d-17
  td_E_eV(16) =  98.0_8; td_cs_cm2(16) = 8.71d-17
  td_E_eV(17) = 108.0_8; td_cs_cm2(17) = 9.00d-17
  td_E_eV(18) = 118.0_8; td_cs_cm2(18) = 9.10d-17
  td_E_eV(19) = 138.0_8; td_cs_cm2(19) = 9.13d-17
  td_E_eV(20) = 158.0_8; td_cs_cm2(20) = 9.05d-17
  td_E_eV(21) = 178.0_8; td_cs_cm2(21) = 8.91d-17
  td_E_eV(22) = 198.0_8; td_cs_cm2(22) = 8.64d-17
  td_E_eV(23) = 223.0_8; td_cs_cm2(23) = 8.30d-17
  td_E_eV(24) = 248.0_8; td_cs_cm2(24) = 7.94d-17
  td_E_eV(25) = 273.0_8; td_cs_cm2(25) = 7.55d-17
  td_E_eV(26) = 298.0_8; td_cs_cm2(26) = 7.21d-17
  td_E_eV(27) = 348.0_8; td_cs_cm2(27) = 6.59d-17
  td_E_eV(28) = 398.0_8; td_cs_cm2(28) = 6.11d-17
  td_E_eV(29) = 448.0_8; td_cs_cm2(29) = 5.62d-17
  td_E_eV(30) = 498.0_8; td_cs_cm2(30) = 5.26d-17
  td_E_eV(31) = 548.0_8; td_cs_cm2(31) = 4.87d-17
  td_E_eV(32) = 598.0_8; td_cs_cm2(32) = 4.57d-17
  td_E_eV(33) = 648.0_8; td_cs_cm2(33) = 4.32d-17
  td_E_eV(34) = 698.0_8; td_cs_cm2(34) = 4.15d-17
  td_E_eV(35) = 748.0_8; td_cs_cm2(35) = 3.88d-17
  td_E_eV(36) = 798.0_8; td_cs_cm2(36) = 3.69d-17
  td_E_eV(37) = 848.0_8; td_cs_cm2(37) = 3.55d-17
  td_E_eV(38) = 898.0_8; td_cs_cm2(38) = 3.36d-17
  td_E_eV(39) = 948.0_8; td_cs_cm2(39) = 3.26d-17
  td_E_eV(40) = 998.0_8; td_cs_cm2(40) = 3.17d-17

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O2_ionO
!
!-------------------------------------
!
real(8) function CSV_O2_ionO_m3s(E_eV)

  use O2_ionO
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  real(8) CS_O2_ionO_cm2

  if (E_eV.lt.td_E_eV(1)) then
     CSV_O2_ionO_m3s = 0.0_8
     return
  end if

  if (E_eV.GT.td_E_eV(N_td)) then
     CS_O2_ionO_cm2 = (7.96d-15 * log(E_eV) - 2.33d-14) / E_eV
     CSV_O2_ionO_m3s = max(0.0_8, 1.0d-4 * CS_O2_ionO_cm2 * sqrt(E_eV) * efactor_eV_to_ms)
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O2_ionO_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O2_ionO_m3s

!---------------------------------------------------------------------------------------------------- #31  !27
! O2, total vibrational cross-section 
! Table 9 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
! threshold is 0.5 eV ~ average of thresholds 0.3(v=1) 0.4(v=2) 0.6(v=3) 0.8(v=4)
!
module O2_vibtot
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O2_vibtot
!
!-------------------------------------
!
subroutine Prepare_O2_vibtot

  use O2_vibtot
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =  1.0_8; td_cs_cm2(1)  = 1.45d-18
  td_E_eV(2)  =  2.0_8; td_cs_cm2(2)  = 3.81d-18
  td_E_eV(3)  =  3.0_8; td_cs_cm2(3)  = 7.20d-18
  td_E_eV(4)  =  4.0_8; td_cs_cm2(4)  = 1.17d-17
  td_E_eV(5)  =  6.0_8; td_cs_cm2(5)  = 3.24d-17
  td_E_eV(6)  =  9.0_8; td_cs_cm2(6)  = 5.98d-17
  td_E_eV(7)  = 12.0_8; td_cs_cm2(7)  = 2.59d-17
  td_E_eV(8)  = 15.0_8; td_cs_cm2(8)  = 7.55d-18
  td_E_eV(9)  = 20.0_8; td_cs_cm2(9)  = 1.91d-18
  td_E_eV(10) = 30.0_8; td_cs_cm2(10) = 3.00d-19

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O2_vibtot
!
!-------------------------------------
!
real(8) function CSV_O2_vibtot_m3s(E_eV)

  use O2_vibtot
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O2_vibtot_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O2_vibtot_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O2_vibtot_m3s

!---------------------------------------------------------------------------------------------------- #32  !28
! O2, Rydberg cross-section 
! Table 9 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module O2_Rydberg
  integer, parameter :: N_td=8  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O2_Rydberg
!
!-------------------------------------
!
subroutine Prepare_O2_Rydberg

  use O2_Rydberg
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1) =   17.0_8; td_cs_cm2(1) = 2.810d-18
  td_E_eV(2) =   20.0_8; td_cs_cm2(2) = 4.180d-17
  td_E_eV(3) =   30.0_8; td_cs_cm2(3) = 1.409d-16
  td_E_eV(4) =   50.0_8; td_cs_cm2(4) = 1.236d-16
  td_E_eV(5) =  100.0_8; td_cs_cm2(5) = 7.590d-17
  td_E_eV(6) =  200.0_8; td_cs_cm2(6) = 4.340d-17
  td_E_eV(7) =  500.0_8; td_cs_cm2(7) = 2.180d-17
  td_E_eV(8) = 1000.0_8; td_cs_cm2(8) = 1.260d-17

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O2_Rydberg
!
!-------------------------------------
!
real(8) function CSV_O2_Rydberg_m3s(E_eV)

  use O2_Rydberg
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O2_Rydberg_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O2_Rydberg_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O2_Rydberg_m3s

!---------------------------------------------------------------------------------------------------- #33  !29
! O2, electronic state A+A'+c cross-section 
! Table 10 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module O2_AApc
  integer, parameter :: N_td=11  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O2_AApc
!
!-------------------------------------
!
subroutine Prepare_O2_AApc

  use O2_AApc
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   6.0_8; td_cs_cm2(1)  = 2.11d-18
  td_E_eV(2)  =   7.0_8; td_cs_cm2(2)  = 8.46d-18
  td_E_eV(3)  =   8.0_8; td_cs_cm2(3)  = 1.41d-17
  td_E_eV(4)  =  10.0_8; td_cs_cm2(4)  = 1.73d-17
  td_E_eV(5)  =  15.0_8; td_cs_cm2(5)  = 1.33d-17
  td_E_eV(6)  =  20.0_8; td_cs_cm2(6)  = 9.74d-18
  td_E_eV(7)  =  30.0_8; td_cs_cm2(7)  = 5.80d-18
  td_E_eV(8)  =  50.0_8; td_cs_cm2(8)  = 2.97d-18
  td_E_eV(9)  = 100.0_8; td_cs_cm2(9)  = 1.31d-18
  td_E_eV(10) = 150.0_8; td_cs_cm2(10) = 8.80d-19
  td_E_eV(11) = 200.0_8; td_cs_cm2(11) = 6.80d-19

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O2_AApc
!
!-------------------------------------
!
real(8) function CSV_O2_AApc_m3s(E_eV)

  use O2_AApc
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O2_AApc_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O2_AApc_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O2_AApc_m3s

!---------------------------------------------------------------------------------------------------- #34  !30
! O2, electronic state a^1\Delta_u cross-section 
! Table 10 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module O2_a1Du
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O2_a1Du
!
!-------------------------------------
!
subroutine Prepare_O2_a1Du

  use O2_a1Du
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   2.0_8; td_cs_cm2(1)  = 9.80d-19
  td_E_eV(2)  =   3.0_8; td_cs_cm2(2)  = 2.83d-18
  td_E_eV(3)  =   4.0_8; td_cs_cm2(3)  = 5.48d-18
  td_E_eV(4)  =   6.0_8; td_cs_cm2(4)  = 8.70d-18
  td_E_eV(5)  =  10.0_8; td_cs_cm2(5)  = 6.61d-18
  td_E_eV(6)  =  20.0_8; td_cs_cm2(6)  = 3.16d-18
  td_E_eV(7)  =  50.0_8; td_cs_cm2(7)  = 1.13d-18
  td_E_eV(8)  = 100.0_8; td_cs_cm2(8)  = 5.50d-19
  td_E_eV(9)  = 150.0_8; td_cs_cm2(9)  = 3.50d-19
  td_E_eV(10) = 200.0_8; td_cs_cm2(10) = 2.50d-19

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O2_a1Du
!
!-------------------------------------
!
real(8) function CSV_O2_a1Du_m3s(E_eV)

  use O2_a1Du
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O2_a1Du_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O2_a1Du_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O2_a1Du_m3s

!---------------------------------------------------------------------------------------------------- #35  !31
! O2, electronic state b^1\Sigma_g^+ cross-section 
! Table 10 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module O2_b1Zgp
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O2_b1Zgp
!
!-------------------------------------
!
subroutine Prepare_O2_b1Zgp

  use O2_b1Zgp
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   2.0_8; td_cs_cm2(1)  = 2.10d-19
  td_E_eV(2)  =   3.0_8; td_cs_cm2(2)  = 6.40d-19
  td_E_eV(3)  =   5.0_8; td_cs_cm2(3)  = 2.14d-18
  td_E_eV(4)  =   7.0_8; td_cs_cm2(4)  = 3.24d-18
  td_E_eV(5)  =   9.0_8; td_cs_cm2(5)  = 2.47d-18
  td_E_eV(6)  =  13.0_8; td_cs_cm2(6)  = 1.17d-18
  td_E_eV(7)  =  20.0_8; td_cs_cm2(7)  = 5.20d-19
  td_E_eV(8)  =  30.0_8; td_cs_cm2(8)  = 2.40d-19
  td_E_eV(9)  =  50.0_8; td_cs_cm2(9)  = 9.00d-20
  td_E_eV(10) = 100.0_8; td_cs_cm2(10) = 2.00d-20

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O2_b1Zgp
!
!-------------------------------------
!
real(8) function CSV_O2_b1Zgp_m3s(E_eV)

  use O2_b1Zgp
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O2_b1Zgp_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O2_b1Zgp_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O2_b1Zgp_m3s

!---------------------------------------------------------------------------------------------------- #36  !32
! O2, longest band cross-section 
! Table 11 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module O2_longband
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O2_longband
!
!-------------------------------------
!
subroutine Prepare_O2_longband

  use O2_longband
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   11.0_8; td_cs_cm2(1)  = 8.30d-19
  td_E_eV(2)  =   12.0_8; td_cs_cm2(2)  = 2.12d-18
  td_E_eV(3)  =   15.0_8; td_cs_cm2(3)  = 5.19d-18
  td_E_eV(4)  =   20.0_8; td_cs_cm2(4)  = 6.70d-18
  td_E_eV(5)  =   30.0_8; td_cs_cm2(5)  = 5.10d-18
  td_E_eV(6)  =   50.0_8; td_cs_cm2(6)  = 2.64d-18
  td_E_eV(7)  =  100.0_8; td_cs_cm2(7)  = 1.00d-18
  td_E_eV(8)  =  200.0_8; td_cs_cm2(8)  = 3.80d-19
  td_E_eV(9)  =  500.0_8; td_cs_cm2(9)  = 1.10d-19
  td_E_eV(10) = 1000.0_8; td_cs_cm2(10) = 2.00d-20

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O2_longband
!
!-------------------------------------
!
real(8) function CSV_O2_longband_m3s(E_eV)

  use O2_longband
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O2_longband_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O2_longband_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O2_longband_m3s

!---------------------------------------------------------------------------------------------------- #37  !33
! O2, second band cross-section 
! Table 11 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module O2_secondband
  integer, parameter :: N_td=8  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O2_secondband
!
!-------------------------------------
!
subroutine Prepare_O2_secondband

  use O2_secondband
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1) =  12.0_8; td_cs_cm2(1) = 4.60d-19
  td_E_eV(2) =  15.0_8; td_cs_cm2(2) = 7.10d-19
  td_E_eV(3) =  20.0_8; td_cs_cm2(3) = 9.30d-19
  td_E_eV(4) =  30.0_8; td_cs_cm2(4) = 8.60d-19
  td_E_eV(5) =  50.0_8; td_cs_cm2(5) = 4.70d-19
  td_E_eV(6) = 100.0_8; td_cs_cm2(6) = 2.00d-19
  td_E_eV(7) = 200.0_8; td_cs_cm2(7) = 8.00d-20
  td_E_eV(8) = 500.0_8; td_cs_cm2(8) = 1.00d-20

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O2_secondband
!
!-------------------------------------
!
real(8) function CSV_O2_secondband_m3s(E_eV)

  use O2_secondband
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O2_secondband_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O2_secondband_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O2_secondband_m3s

!---------------------------------------------------------------------------------------------------- #38  !34
! O2, 1 ^3\Pi_g cross-section 
! Table 11 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module O2_13Pg
  integer, parameter :: N_td=11  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O2_13Pg
!
!-------------------------------------
!
subroutine Prepare_O2_13Pg

  use O2_13Pg
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   8.0_8; td_cs_cm2(1)  = 1.70d-19
  td_E_eV(2)  =   9.0_8; td_cs_cm2(2)  = 5.90d-19
  td_E_eV(3)  =  10.0_8; td_cs_cm2(3)  = 1.37d-18
  td_E_eV(4)  =  11.0_8; td_cs_cm2(4)  = 2.10d-18
  td_E_eV(5)  =  13.0_8; td_cs_cm2(5)  = 4.31d-18
  td_E_eV(6)  =  18.0_8; td_cs_cm2(6)  = 7.00d-18
  td_E_eV(7)  =  24.0_8; td_cs_cm2(7)  = 5.02d-18
  td_E_eV(8)  =  30.0_8; td_cs_cm2(8)  = 3.07d-18
  td_E_eV(9)  =  50.0_8; td_cs_cm2(9)  = 9.50d-19
  td_E_eV(10) = 100.0_8; td_cs_cm2(10) = 2.10d-19
  td_E_eV(11) = 200.0_8; td_cs_cm2(11) = 3.00d-20

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O2_13Pg
!
!-------------------------------------
!
real(8) function CSV_O2_13Pg_m3s(E_eV)

  use O2_13Pg
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O2_13Pg_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O2_13Pg_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O2_13Pg_m3s

!---------------------------------------------------------------------------------------------------- #39  !35
! O2, 8.9 eV peak cross-section 
! Table 11 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module O2_8p9eV
  integer, parameter :: N_td=12  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O2_8p9eV
!
!-------------------------------------
!
subroutine Prepare_O2_8p9eV

  use O2_8p9eV
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =    9.0_8; td_cs_cm2(1)  = 5.80d-19
  td_E_eV(2)  =   10.0_8; td_cs_cm2(2)  = 1.62d-18
  td_E_eV(3)  =   12.0_8; td_cs_cm2(3)  = 4.13d-18
  td_E_eV(4)  =   15.0_8; td_cs_cm2(4)  = 9.18d-18
  td_E_eV(5)  =   20.0_8; td_cs_cm2(5)  = 1.42d-17
  td_E_eV(6)  =   25.0_8; td_cs_cm2(6)  = 1.49d-17
  td_E_eV(7)  =   30.0_8; td_cs_cm2(7)  = 1.27d-17
  td_E_eV(8)  =   50.0_8; td_cs_cm2(8)  = 6.71d-18
  td_E_eV(9)  =  100.0_8; td_cs_cm2(9)  = 2.71d-18
  td_E_eV(10) =  200.0_8; td_cs_cm2(10) = 1.09d-18
  td_E_eV(11) =  500.0_8; td_cs_cm2(11) = 3.40d-19
  td_E_eV(12) = 1000.0_8; td_cs_cm2(12) = 1.40d-19

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O2_8p9eV
!
!-------------------------------------
!
real(8) function CSV_O2_8p9eV_m3s(E_eV)

  use O2_8p9eV
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O2_8p9eV_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O2_8p9eV_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O2_8p9eV_m3s

!---------------------------------------------------------------------------------------------------- #40  !36
! O2, B ^3\Sigma_u cross-section 
! Table 11 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
!
module O2_B3Zu
  integer, parameter :: N_td=11  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O2_B3Zu
!
!-------------------------------------
!
subroutine Prepare_O2_B3Zu

  use O2_B3Zu
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =    9.0_8; td_cs_cm2(1)  = 7.50d-19
  td_E_eV(2)  =   10.0_8; td_cs_cm2(2)  = 6.39d-18
  td_E_eV(3)  =   12.0_8; td_cs_cm2(3)  = 2.98d-17
  td_E_eV(4)  =   15.0_8; td_cs_cm2(4)  = 5.29d-17
  td_E_eV(5)  =   20.0_8; td_cs_cm2(5)  = 5.58d-17
  td_E_eV(6)  =   30.0_8; td_cs_cm2(6)  = 3.91d-17
  td_E_eV(7)  =   50.0_8; td_cs_cm2(7)  = 1.99d-17
  td_E_eV(8)  =  100.0_8; td_cs_cm2(8)  = 7.48d-18
  td_E_eV(9)  =  200.0_8; td_cs_cm2(9)  = 2.85d-18
  td_E_eV(10) =  500.0_8; td_cs_cm2(10) = 8.20d-19
  td_E_eV(11) = 1000.0_8; td_cs_cm2(11) = 3.00d-19

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O2_B3Zu
!
!-------------------------------------
!
real(8) function CSV_O2_B3Zu_m3s(E_eV)

  use O2_B3Zu
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O2_B3Zu_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O2_B3Zu_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O2_B3Zu_m3s

!---------------------------------------------------------------------------------------------------- #41  !37
!
! A. Schmalzried, A. Luque and N. Lehtinen,  IAA Database on lxcat, www.lxcat.net/IAA, August 2023,
!                  Instituto de Astrofísica de Andalucía.
!
module O_elastic
!  integer, parameter :: N_td=28  ! number of table data points
  integer, parameter :: N_td=61  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_elastic
!
!-------------------------------------
!
subroutine Prepare_O_elastic

  use O_elastic
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(  1) = 0.000000d+00 ; td_cs_cm2(  1) = 1.400000d-20   !### original data in m^2
  td_E_eV(  2) = 1.000000d-03 ; td_cs_cm2(  2) = 1.610000d-20
  td_E_eV(  3) = 1.000000d-02 ; td_cs_cm2(  3) = 2.060000d-20
  td_E_eV(  4) = 5.400000d-01 ; td_cs_cm2(  4) = 2.864490d-20
  td_E_eV(  5) = 2.180000d+00 ; td_cs_cm2(  5) = 5.774350d-20
  td_E_eV(  6) = 3.400000d+00 ; td_cs_cm2(  6) = 7.063500d-20
  td_E_eV(  7) = 4.900000d+00 ; td_cs_cm2(  7) = 7.498520d-20
  td_E_eV(  8) = 8.700000d+00 ; td_cs_cm2(  8) = 8.400360d-20
  td_E_eV(  9) = 1.000000d+01 ; td_cs_cm2(  9) = 7.350530d-20
  td_E_eV( 10) = 1.100000d+01 ; td_cs_cm2( 10) = 7.295590d-20
  td_E_eV( 11) = 1.200000d+01 ; td_cs_cm2( 11) = 7.217360d-20
  td_E_eV( 12) = 1.300000d+01 ; td_cs_cm2( 12) = 7.123100d-20
  td_E_eV( 13) = 1.400000d+01 ; td_cs_cm2( 13) = 7.017720d-20
  td_E_eV( 14) = 1.500000d+01 ; td_cs_cm2( 14) = 6.905010d-20
  td_E_eV( 15) = 1.600000d+01 ; td_cs_cm2( 15) = 6.787900d-20
  td_E_eV( 16) = 1.700000d+01 ; td_cs_cm2( 16) = 6.668690d-20
  td_E_eV( 17) = 1.800000d+01 ; td_cs_cm2( 17) = 6.548620d-20
  td_E_eV( 18) = 1.900000d+01 ; td_cs_cm2( 18) = 6.429090d-20
  td_E_eV( 19) = 2.000000d+01 ; td_cs_cm2( 19) = 6.310860d-20
  td_E_eV( 20) = 2.500000d+01 ; td_cs_cm2( 20) = 5.758450d-20
  td_E_eV( 21) = 3.000000d+01 ; td_cs_cm2( 21) = 5.285330d-20
  td_E_eV( 22) = 3.500000d+01 ; td_cs_cm2( 22) = 4.887070d-20
  td_E_eV( 23) = 4.000000d+01 ; td_cs_cm2( 23) = 4.548370d-20
  td_E_eV( 24) = 4.500000d+01 ; td_cs_cm2( 24) = 4.256410d-20
  td_E_eV( 25) = 5.000000d+01 ; td_cs_cm2( 25) = 4.002110d-20
  td_E_eV( 26) = 5.500000d+01 ; td_cs_cm2( 26) = 3.779270d-20
  td_E_eV( 27) = 6.000000d+01 ; td_cs_cm2( 27) = 3.583240d-20
  td_E_eV( 28) = 6.500000d+01 ; td_cs_cm2( 28) = 3.409700d-20
  td_E_eV( 29) = 7.000000d+01 ; td_cs_cm2( 29) = 3.255030d-20
  td_E_eV( 30) = 7.500000d+01 ; td_cs_cm2( 30) = 3.116150d-20
  td_E_eV( 31) = 8.000000d+01 ; td_cs_cm2( 31) = 2.990720d-20
  td_E_eV( 32) = 8.500000d+01 ; td_cs_cm2( 32) = 2.876730d-20
  td_E_eV( 33) = 9.000000d+01 ; td_cs_cm2( 33) = 2.772530d-20
  td_E_eV( 34) = 9.500000d+01 ; td_cs_cm2( 34) = 2.677080d-20
  td_E_eV( 35) = 1.000000d+02 ; td_cs_cm2( 35) = 2.589340d-20
  td_E_eV( 36) = 1.100000d+02 ; td_cs_cm2( 36) = 2.433430d-20
  td_E_eV( 37) = 1.200000d+02 ; td_cs_cm2( 37) = 2.299300d-20
  td_E_eV( 38) = 1.300000d+02 ; td_cs_cm2( 38) = 2.182630d-20
  td_E_eV( 39) = 1.400000d+02 ; td_cs_cm2( 39) = 2.080030d-20
  td_E_eV( 40) = 1.500000d+02 ; td_cs_cm2( 40) = 1.988950d-20
  td_E_eV( 41) = 1.600000d+02 ; td_cs_cm2( 41) = 1.907450d-20
  td_E_eV( 42) = 1.700000d+02 ; td_cs_cm2( 42) = 1.834010d-20
  td_E_eV( 43) = 1.800000d+02 ; td_cs_cm2( 43) = 1.767460d-20
  td_E_eV( 44) = 1.900000d+02 ; td_cs_cm2( 44) = 1.706760d-20
  td_E_eV( 45) = 2.000000d+02 ; td_cs_cm2( 45) = 1.651170d-20
  td_E_eV( 46) = 2.500000d+02 ; td_cs_cm2( 46) = 1.429610d-20
  td_E_eV( 47) = 3.000000d+02 ; td_cs_cm2( 47) = 1.270460d-20
  td_E_eV( 48) = 3.500000d+02 ; td_cs_cm2( 48) = 1.149850d-20
  td_E_eV( 49) = 4.000000d+02 ; td_cs_cm2( 49) = 1.054770d-20
  td_E_eV( 50) = 4.500000d+02 ; td_cs_cm2( 50) = 9.768040d-21
  td_E_eV( 51) = 5.000000d+02 ; td_cs_cm2( 51) = 9.110740d-21
  td_E_eV( 52) = 5.500000d+02 ; td_cs_cm2( 52) = 8.546710d-21
  td_E_eV( 53) = 6.000000d+02 ; td_cs_cm2( 53) = 8.055050d-21
  td_E_eV( 54) = 6.500000d+02 ; td_cs_cm2( 54) = 7.621540d-21
  td_E_eV( 55) = 7.000000d+02 ; td_cs_cm2( 55) = 7.235750d-21
  td_E_eV( 56) = 7.500000d+02 ; td_cs_cm2( 56) = 6.889220d-21
  td_E_eV( 57) = 8.000000d+02 ; td_cs_cm2( 57) = 6.576020d-21
  td_E_eV( 58) = 8.500000d+02 ; td_cs_cm2( 58) = 6.291600d-21
  td_E_eV( 59) = 9.000000d+02 ; td_cs_cm2( 59) = 6.031970d-21
  td_E_eV( 60) = 9.500000d+02 ; td_cs_cm2( 60) = 5.793590d-21
  td_E_eV( 61) = 1.000000d+03 ; td_cs_cm2( 61) = 5.573770d-21

  td_cs_cm2 = td_cs_cm2 * 1.0d4   ! convert m^2 to cm^2

! I could not find table data however I found the following paper
! Y. Itikawa and A. Ichimura, "Cross Sections for Collisions and Photons with Atomic Oxygen"
! J. Phys. Chem. Ref. Data, Vol.19, No.3, pp.637-651, 1990
! where Figures 4.1 and 4.2 show curves of the elastic cross sections
! so the data below are retrieved manually from these figures.
!  td_E_eV(1)   =    1.0_8  ; td_cs_cm2(1)   = 4.2e-16
!  td_E_eV(2)   =    2.0_8  ; td_cs_cm2(2)   = 6.0e-16
!  td_E_eV(3)   =    3.0_8  ; td_cs_cm2(3)   = 6.7e-16
!  td_E_eV(4)   =    4.0_8  ; td_cs_cm2(4)   = 7.0e-16
!  td_E_eV(5)   =    5.0_8  ; td_cs_cm2(5)   = 7.3e-16
!  td_E_eV(6)   =    6.0_8  ; td_cs_cm2(6)   = 7.3e-16
!  td_E_eV(7)   =    7.0_8  ; td_cs_cm2(7)   = 7.3e-16
!  td_E_eV(8)   =    8.0_8  ; td_cs_cm2(8)   = 7.3e-16
!  td_E_eV(9)   =    9.0_8  ; td_cs_cm2(9)   = 7.2e-16
!  td_E_eV(10)  =   10.0_8  ; td_cs_cm2(10)  = 7.1e-16
!  td_E_eV(11)  =   20.0_8  ; td_cs_cm2(11)  = 6.0e-16
!  td_E_eV(12)  =   30.0_8  ; td_cs_cm2(12)  = 5.3e-16
!  td_E_eV(13)  =   40.0_8  ; td_cs_cm2(13)  = 4.6e-16
!  td_E_eV(14)  =   50.0_8  ; td_cs_cm2(14)  = 4.1e-16
!  td_E_eV(15)  =   60.0_8  ; td_cs_cm2(15)  = 3.7e-16
!  td_E_eV(16)  =   70.0_8  ; td_cs_cm2(16)  = 3.4e-16
!  td_E_eV(17)  =   80.0_8  ; td_cs_cm2(17)  = 3.2e-16
!  td_E_eV(18)  =   90.0_8  ; td_cs_cm2(18)  = 2.9e-16
!  td_E_eV(19)  =  100.0_8  ; td_cs_cm2(19)  = 2.7e-16
!  td_E_eV(20)  =  200.0_8  ; td_cs_cm2(20)  = 1.7e-16
!  td_E_eV(21)  =  300.0_8  ; td_cs_cm2(21)  = 1.3e-16
!  td_E_eV(22)  =  400.0_8  ; td_cs_cm2(22)  = 1.05e-16
!  td_E_eV(23)  =  500.0_8  ; td_cs_cm2(23)  = 9.0e-17
!  td_E_eV(24)  =  600.0_8  ; td_cs_cm2(24)  = 7.8e-17
!  td_E_eV(25)  =  700.0_8  ; td_cs_cm2(25)  = 7.0e-17
!  td_E_eV(26)  =  800.0_8  ; td_cs_cm2(26)  = 6.2e-17
!  td_E_eV(27)  =  900.0_8  ; td_cs_cm2(27)  = 5.8e-17
!  td_E_eV(28)  = 1000.0_8  ; td_cs_cm2(28)  = 5.3e-17
  
  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_elastic
!
!-------------------------------------
!
real(8) function CSV_O_elastic_m3s(E_eV)

  use O_elastic
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  real(8) CS_O_elastic_cm2

  if (E_eV.lt.td_E_eV(1)) then
     CSV_O_elastic_m3s = 1.0d-4 * td_cs_cm2(1) * sqrt(E_eV) * efactor_eV_to_ms
     return
  end if

  if (E_eV.GT.td_E_eV(N_td)) then
     CS_O_elastic_cm2 = (1.0d-15 * log(E_eV) + 1.966d-14) / E_eV**0.9
     CSV_O_elastic_m3s = max(0.0_8, 1.0d-4 * CS_O_elastic_cm2 * sqrt(E_eV) * efactor_eV_to_ms)
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_elastic_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_elastic_m3s















!---------------------------------------------------------------------------------------------------- #42  !38
! O -> O+(4S0) ionization cross-section
! Table 23 of "Updated excitation and ionization cross sections for electron impact on atomic oxygen" by
! R. R. Laher and F. R. Gilmore, technical report AD-A206 987, 1988.
!
module O_ion4S0
 integer, parameter :: N_td=24  ! number of table data points
! integer, parameter :: N_td=10  ! number of table data points

! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_ion4S0
!
!-------------------------------------
!
subroutine Prepare_O_ion4S0

  use O_ion4S0
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV( 1) =   14.0_8 ; td_cs_cm2( 1) = 2.60d-18
  td_E_eV( 2) =   16.0_8 ; td_cs_cm2( 2) = 1.63d-17
  td_E_eV( 3) =   18.0_8 ; td_cs_cm2( 3) = 2.30d-17
  td_E_eV( 4) =   20.0_8 ; td_cs_cm2( 4) = 2.76d-17
  td_E_eV( 5) =   25.0_8 ; td_cs_cm2( 5) = 3.61d-17
  td_E_eV( 6) =   30.0_8 ; td_cs_cm2( 6) = 4.48d-17
  td_E_eV( 7) =   35.0_8 ; td_cs_cm2( 7) = 5.27d-17
  td_E_eV( 8) =   40.0_8 ; td_cs_cm2( 8) = 6.08d-17
  td_E_eV( 9) =   50.0_8 ; td_cs_cm2( 9) = 6.85d-17
  td_E_eV(10) =   60.0_8 ; td_cs_cm2(10) = 7.01d-17
  td_E_eV(11) =   70.0_8 ; td_cs_cm2(11) = 6.80d-17
  td_E_eV(12) =   80.0_8 ; td_cs_cm2(12) = 6.61d-17
  td_E_eV(13) =   90.0_8 ; td_cs_cm2(13) = 6.30d-17
  td_E_eV(14) =  100.0_8 ; td_cs_cm2(14) = 6.16d-17
  td_E_eV(15) =  150.0_8 ; td_cs_cm2(15) = 5.50d-17
  td_E_eV(16) =  200.0_8 ; td_cs_cm2(16) = 4.84d-17
  td_E_eV(17) =  300.0_8 ; td_cs_cm2(17) = 3.84d-17
  td_E_eV(18) =  400.0_8 ; td_cs_cm2(18) = 3.23d-17
  td_E_eV(19) =  500.0_8 ; td_cs_cm2(19) = 2.82d-17
  td_E_eV(20) =  600.0_8 ; td_cs_cm2(20) = 2.47d-17
  td_E_eV(21) =  700.0_8 ; td_cs_cm2(21) = 2.21d-17
  td_E_eV(22) =  800.0_8 ; td_cs_cm2(22) = 2.01d-17
  td_E_eV(23) =  900.0_8 ; td_cs_cm2(23) = 1.84d-17
  td_E_eV(24) = 1000.0_8 ; td_cs_cm2(24) = 1.71d-17

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_ion4S0
!
!-------------------------------------
!
real(8) function CSV_O_ion4S0_m3s(E_eV)

  use O_ion4S0
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  real(8) CS_O_ion4S0_cm2

  if (E_eV.lt.td_E_eV(1)) then
     CSV_O_ion4S0_m3s = 0.0_8
     return
  end if

  if (E_eV.gt.td_E_eV(N_td)) then
!     CS_O_ion4S0_cm2 = (1.49d-14 * log(E_eV) - 5.27d-14) / E_eV
!     CSV_O_ion4S0_m3s = max(0.0_8, 1.0d-4 * CS_O_ion4S0_cm2 * sqrt(E_eV) * efactor_eV_to_ms)
     CSV_O_ion4S0_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_ion4S0_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_ion4S0_m3s



!---------------------------------------------------------------------------------------------------- #43  !38
! O -> O+(2D0) ionization cross-section
! Table 23 of "Updated excitation and ionization cross sections for electron impact on atomic oxygen" by
! R. R. Laher and F. R. Gilmore, technical report AD-A206 987, 1988.
!
module O_ion2D0
 integer, parameter :: N_td=22  ! number of table data points
! integer, parameter :: N_td=10  ! number of table data points

! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_ion2D0
!
!-------------------------------------
!
subroutine Prepare_O_ion2D0

  use O_ion2D0
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV( 1) =   18.0_8 ; td_cs_cm2( 1) = 2.20d-18
  td_E_eV( 2) =   20.0_8 ; td_cs_cm2( 2) = 7.00d-18
  td_E_eV( 3) =   25.0_8 ; td_cs_cm2( 3) = 1.28d-17
  td_E_eV( 4) =   30.0_8 ; td_cs_cm2( 4) = 1.63d-17
  td_E_eV( 5) =   35.0_8 ; td_cs_cm2( 5) = 1.85d-17
  td_E_eV( 6) =   40.0_8 ; td_cs_cm2( 6) = 2.06d-17
  td_E_eV( 7) =   50.0_8 ; td_cs_cm2( 7) = 2.38d-17
  td_E_eV( 8) =   60.0_8 ; td_cs_cm2( 8) = 2.58d-17
  td_E_eV( 9) =   70.0_8 ; td_cs_cm2( 9) = 2.82d-17
  td_E_eV(10) =   80.0_8 ; td_cs_cm2(10) = 3.06d-17
  td_E_eV(11) =   90.0_8 ; td_cs_cm2(11) = 3.23d-17
  td_E_eV(12) =  100.0_8 ; td_cs_cm2(12) = 3.40d-17
  td_E_eV(13) =  150.0_8 ; td_cs_cm2(13) = 3.62d-17
  td_E_eV(14) =  200.0_8 ; td_cs_cm2(14) = 3.60d-17
  td_E_eV(15) =  300.0_8 ; td_cs_cm2(15) = 3.03d-17
  td_E_eV(16) =  400.0_8 ; td_cs_cm2(16) = 2.60d-17
  td_E_eV(17) =  500.0_8 ; td_cs_cm2(17) = 2.30d-17
  td_E_eV(18) =  600.0_8 ; td_cs_cm2(18) = 2.03d-17
  td_E_eV(19) =  700.0_8 ; td_cs_cm2(19) = 1.82d-17
  td_E_eV(20) =  800.0_8 ; td_cs_cm2(20) = 1.66d-17
  td_E_eV(21) =  900.0_8 ; td_cs_cm2(21) = 1.53d-17
  td_E_eV(22) = 1000.0_8 ; td_cs_cm2(22) = 1.42d-17

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_ion2D0
!
!-------------------------------------
!
real(8) function CSV_O_ion2D0_m3s(E_eV)

  use O_ion2D0
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  real(8) CS_O_ion2D0_cm2

  if (E_eV.lt.td_E_eV(1)) then
     CSV_O_ion2D0_m3s = 0.0_8
     return
  end if

  if (E_eV.gt.td_E_eV(N_td)) then
!     CS_O_ion2D0_cm2 = (1.49d-14 * log(E_eV) - 5.27d-14) / E_eV
!     CSV_O_ion2D0_m3s = max(0.0_8, 1.0d-4 * CS_O_ion2D0_cm2 * sqrt(E_eV) * efactor_eV_to_ms)
     CSV_O_ion2D0_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_ion2D0_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_ion2D0_m3s



!---------------------------------------------------------------------------------------------------- #44  !38
! O -> O+(2P0) ionization cross-section
! Table 23 of "Updated excitation and ionization cross sections for electron impact on atomic oxygen" by
! R. R. Laher and F. R. Gilmore, technical report AD-A206 987, 1988.
!
module O_ion2P0
 integer, parameter :: N_td=21  ! number of table data points
! integer, parameter :: N_td=10  ! number of table data points

! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_ion2P0
!
!-------------------------------------
!
subroutine Prepare_O_ion2P0

  use O_ion2P0
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV( 1) =   20.0_8 ; td_cs_cm2( 1) = 1.40d-18
  td_E_eV( 2) =   25.0_8 ; td_cs_cm2( 2) = 4.20d-18
  td_E_eV( 3) =   30.0_8 ; td_cs_cm2( 3) = 6.60d-18
  td_E_eV( 4) =   35.0_8 ; td_cs_cm2( 4) = 8.00d-18
  td_E_eV( 5) =   40.0_8 ; td_cs_cm2( 5) = 9.40d-18
  td_E_eV( 6) =   50.0_8 ; td_cs_cm2( 6) = 1.15d-17
  td_E_eV( 7) =   60.0_8 ; td_cs_cm2( 7) = 1.33d-17
  td_E_eV( 8) =   70.0_8 ; td_cs_cm2( 8) = 1.50d-17
  td_E_eV( 9) =   80.0_8 ; td_cs_cm2( 9) = 1.62d-17
  td_E_eV(10) =   90.0_8 ; td_cs_cm2(10) = 1.79d-17
  td_E_eV(11) =  100.0_8 ; td_cs_cm2(11) = 1.89d-17
  td_E_eV(12) =  150.0_8 ; td_cs_cm2(12) = 2.04d-17
  td_E_eV(13) =  200.0_8 ; td_cs_cm2(13) = 2.04d-17
  td_E_eV(14) =  300.0_8 ; td_cs_cm2(14) = 1.85d-17
  td_E_eV(15) =  400.0_8 ; td_cs_cm2(15) = 1.58d-17
  td_E_eV(16) =  500.0_8 ; td_cs_cm2(16) = 1.39d-17
  td_E_eV(17) =  600.0_8 ; td_cs_cm2(17) = 1.21d-17
  td_E_eV(18) =  700.0_8 ; td_cs_cm2(18) = 1.08d-17
  td_E_eV(19) =  800.0_8 ; td_cs_cm2(19) = 9.70d-18
  td_E_eV(20) =  900.0_8 ; td_cs_cm2(20) = 8.90d-18
  td_E_eV(21) = 1000.0_8 ; td_cs_cm2(21) = 8.20d-18

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_ion2P0
!
!-------------------------------------
!
real(8) function CSV_O_ion2P0_m3s(E_eV)

  use O_ion2P0
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  real(8) CS_O_ion2P0_cm2

  if (E_eV.lt.td_E_eV(1)) then
     CSV_O_ion2P0_m3s = 0.0_8
     return
  end if

  if (E_eV.gt.td_E_eV(N_td)) then
!     CS_O_ion2P0_cm2 = (1.49d-14 * log(E_eV) - 5.27d-14) / E_eV
!     CSV_O_ion2P0_m3s = max(0.0_8, 1.0d-4 * CS_O_ion2P0_cm2 * sqrt(E_eV) * efactor_eV_to_ms)
     CSV_O_ion2P0_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_ion2P0_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_ion2P0_m3s




!---------------------------------------------------------------------------------------------------- #45  !38
! O -> O+(4P) ionization cross-section
! Table 23 of "Updated excitation and ionization cross sections for electron impact on atomic oxygen" by
! R. R. Laher and F. R. Gilmore, technical report AD-A206 987, 1988.
!
module O_ion4P
 integer, parameter :: N_td=19  ! number of table data points
! integer, parameter :: N_td=10  ! number of table data points

! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_ion4P
!
!-------------------------------------
!
subroutine Prepare_O_ion4P

  use O_ion4P
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV( 1) =   30.0_8 ; td_cs_cm2( 1) = 1.00d-18
  td_E_eV( 2) =   35.0_8 ; td_cs_cm2( 2) = 2.50d-18
  td_E_eV( 3) =   40.0_8 ; td_cs_cm2( 3) = 5.70d-18
  td_E_eV( 4) =   50.0_8 ; td_cs_cm2( 4) = 1.15d-17
  td_E_eV( 5) =   60.0_8 ; td_cs_cm2( 5) = 1.51d-17
  td_E_eV( 6) =   70.0_8 ; td_cs_cm2( 6) = 1.70d-17
  td_E_eV( 7) =   80.0_8 ; td_cs_cm2( 7) = 1.81d-17
  td_E_eV( 8) =   90.0_8 ; td_cs_cm2( 8) = 1.88d-17
  td_E_eV( 9) =  100.0_8 ; td_cs_cm2( 9) = 1.91d-17
  td_E_eV(10) =  150.0_8 ; td_cs_cm2(10) = 1.95d-17
  td_E_eV(11) =  200.0_8 ; td_cs_cm2(11) = 1.92d-17
  td_E_eV(12) =  300.0_8 ; td_cs_cm2(12) = 1.70d-17
  td_E_eV(13) =  400.0_8 ; td_cs_cm2(13) = 1.53d-17
  td_E_eV(14) =  500.0_8 ; td_cs_cm2(14) = 1.40d-17
  td_E_eV(15) =  600.0_8 ; td_cs_cm2(15) = 1.24d-17
  td_E_eV(16) =  700.0_8 ; td_cs_cm2(16) = 1.11d-17
  td_E_eV(17) =  800.0_8 ; td_cs_cm2(17) = 1.01d-17
  td_E_eV(18) =  900.0_8 ; td_cs_cm2(18) = 9.30d-18
  td_E_eV(19) = 1000.0_8 ; td_cs_cm2(19) = 8.70d-18

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_ion4P
!
!-------------------------------------
!
real(8) function CSV_O_ion4P_m3s(E_eV)

  use O_ion4P
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  real(8) CS_O_ion4P_cm2

  if (E_eV.lt.td_E_eV(1)) then
     CSV_O_ion4P_m3s = 0.0_8
     return
  end if

  if (E_eV.gt.td_E_eV(N_td)) then
!     CS_O_ion4P_cm2 = (1.49d-14 * log(E_eV) - 5.27d-14) / E_eV
!     CSV_O_ion4P_m3s = max(0.0_8, 1.0d-4 * CS_O_ion4P_cm2 * sqrt(E_eV) * efactor_eV_to_ms)
     CSV_O_ion4P_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_ion4P_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_ion4P_m3s






!!---------------------------------------------------------------------------------------------------- #42  !38
!! O, total {O->(O+,O++(?))} ionization cross-section
!! Table 12 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!!
!module O_iontot
! integer, parameter :: N_td=200  ! number of table data points
!! integer, parameter :: N_td=10  ! number of table data points
!
!! table data
!  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
!  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
!  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
!end module O_iontot
!!
!!-------------------------------------
!!
!subroutine Prepare_O_iontot
!
!  use O_iontot
!  use photoelectrons, only : efactor_eV_to_ms
!
!  implicit none
!
!  integer n
!
! td_E_eV(  1) = 1.360000e+1 ; td_cs_cm2(  1) = 0.000000e+0
! td_E_eV(  2) = 1.389691e+1 ; td_cs_cm2(  2) = 1.083592e-22
! td_E_eV(  3) = 1.420029e+1 ; td_cs_cm2(  3) = 2.185236e-22
! td_E_eV(  4) = 1.451030e+1 ; td_cs_cm2(  4) = 3.266507e-22
! td_E_eV(  5) = 1.482708e+1 ; td_cs_cm2(  5) = 4.327528e-22
! td_E_eV(  6) = 1.515078e+1 ; td_cs_cm2(  6) = 5.368437e-22
! td_E_eV(  7) = 1.548154e+1 ; td_cs_cm2(  7) = 6.389386e-22
! td_E_eV(  8) = 1.581952e+1 ; td_cs_cm2(  8) = 7.390539e-22
! td_E_eV(  9) = 1.616488e+1 ; td_cs_cm2(  9) = 8.372069e-22
! td_E_eV( 10) = 1.651778e+1 ; td_cs_cm2( 10) = 9.334160e-22
! td_E_eV( 11) = 1.687839e+1 ; td_cs_cm2( 11) = 1.027700e-21
! td_E_eV( 12) = 1.724686e+1 ; td_cs_cm2( 12) = 1.125028e-21
! td_E_eV( 13) = 1.762338e+1 ; td_cs_cm2( 13) = 1.228676e-21
! td_E_eV( 14) = 1.800813e+1 ; td_cs_cm2( 14) = 1.337865e-21
! td_E_eV( 15) = 1.840127e+1 ; td_cs_cm2( 15) = 1.452085e-21
! td_E_eV( 16) = 1.880299e+1 ; td_cs_cm2( 16) = 1.576383e-21
! td_E_eV( 17) = 1.921348e+1 ; td_cs_cm2( 17) = 1.714696e-21
! td_E_eV( 18) = 1.963294e+1 ; td_cs_cm2( 18) = 1.860014e-21
! td_E_eV( 19) = 2.006155e+1 ; td_cs_cm2( 19) = 2.011685e-21
! td_E_eV( 20) = 2.049952e+1 ; td_cs_cm2( 20) = 2.169091e-21
! td_E_eV( 21) = 2.094705e+1 ; td_cs_cm2( 21) = 2.331650e-21
! td_E_eV( 22) = 2.140435e+1 ; td_cs_cm2( 22) = 2.498812e-21
! td_E_eV( 23) = 2.187164e+1 ; td_cs_cm2( 23) = 2.670061e-21
! td_E_eV( 24) = 2.234913e+1 ; td_cs_cm2( 24) = 2.844906e-21
! td_E_eV( 25) = 2.283704e+1 ; td_cs_cm2( 25) = 3.022887e-21
! td_E_eV( 26) = 2.333560e+1 ; td_cs_cm2( 26) = 3.203568e-21
! td_E_eV( 27) = 2.384504e+1 ; td_cs_cm2( 27) = 3.386538e-21
! td_E_eV( 28) = 2.436561e+1 ; td_cs_cm2( 28) = 3.571410e-21
! td_E_eV( 29) = 2.489755e+1 ; td_cs_cm2( 29) = 3.757818e-21
! td_E_eV( 30) = 2.544109e+1 ; td_cs_cm2( 30) = 3.945415e-21
! td_E_eV( 31) = 2.599650e+1 ; td_cs_cm2( 31) = 4.133876e-21
! td_E_eV( 32) = 2.656404e+1 ; td_cs_cm2( 32) = 4.322892e-21
! td_E_eV( 33) = 2.714397e+1 ; td_cs_cm2( 33) = 4.512173e-21
! td_E_eV( 34) = 2.773656e+1 ; td_cs_cm2( 34) = 4.701444e-21
! td_E_eV( 35) = 2.834208e+1 ; td_cs_cm2( 35) = 4.890446e-21
! td_E_eV( 36) = 2.896083e+1 ; td_cs_cm2( 36) = 5.078934e-21
! td_E_eV( 37) = 2.959308e+1 ; td_cs_cm2( 37) = 5.266676e-21
! td_E_eV( 38) = 3.023914e+1 ; td_cs_cm2( 38) = 5.453454e-21
! td_E_eV( 39) = 3.089929e+1 ; td_cs_cm2( 39) = 5.639062e-21
! td_E_eV( 40) = 3.157387e+1 ; td_cs_cm2( 40) = 5.823306e-21
! td_E_eV( 41) = 3.226317e+1 ; td_cs_cm2( 41) = 6.006000e-21
! td_E_eV( 42) = 3.296751e+1 ; td_cs_cm2( 42) = 6.186973e-21
! td_E_eV( 43) = 3.368724e+1 ; td_cs_cm2( 43) = 6.366059e-21
! td_E_eV( 44) = 3.442267e+1 ; td_cs_cm2( 44) = 6.544246e-21
! td_E_eV( 45) = 3.517416e+1 ; td_cs_cm2( 45) = 6.723017e-21
! td_E_eV( 46) = 3.594206e+1 ; td_cs_cm2( 46) = 6.901962e-21
! td_E_eV( 47) = 3.672672e+1 ; td_cs_cm2( 47) = 7.080762e-21
! td_E_eV( 48) = 3.752851e+1 ; td_cs_cm2( 48) = 7.259119e-21
! td_E_eV( 49) = 3.834781e+1 ; td_cs_cm2( 49) = 7.436749e-21
! td_E_eV( 50) = 3.918499e+1 ; td_cs_cm2( 50) = 7.613388e-21
! td_E_eV( 51) = 4.004045e+1 ; td_cs_cm2( 51) = 7.788784e-21
! td_E_eV( 52) = 4.091459e+1 ; td_cs_cm2( 52) = 7.962702e-21
! td_E_eV( 53) = 4.180781e+1 ; td_cs_cm2( 53) = 8.134919e-21
! td_E_eV( 54) = 4.272053e+1 ; td_cs_cm2( 54) = 8.305227e-21
! td_E_eV( 55) = 4.365317e+1 ; td_cs_cm2( 55) = 8.473429e-21
! td_E_eV( 56) = 4.460618e+1 ; td_cs_cm2( 56) = 8.639341e-21
! td_E_eV( 57) = 4.557999e+1 ; td_cs_cm2( 57) = 8.802791e-21
! td_E_eV( 58) = 4.657506e+1 ; td_cs_cm2( 58) = 8.963616e-21
! td_E_eV( 59) = 4.759185e+1 ; td_cs_cm2( 59) = 9.121665e-21
! td_E_eV( 60) = 4.863084e+1 ; td_cs_cm2( 60) = 9.276796e-21
! td_E_eV( 61) = 4.969252e+1 ; td_cs_cm2( 61) = 9.428878e-21
! td_E_eV( 62) = 5.077737e+1 ; td_cs_cm2( 62) = 9.577788e-21
! td_E_eV( 63) = 5.188591e+1 ; td_cs_cm2( 63) = 9.723411e-21
! td_E_eV( 64) = 5.301864e+1 ; td_cs_cm2( 64) = 9.865641e-21
! td_E_eV( 65) = 5.417611e+1 ; td_cs_cm2( 65) = 1.000438e-20
! td_E_eV( 66) = 5.535884e+1 ; td_cs_cm2( 66) = 1.013954e-20
! td_E_eV( 67) = 5.656740e+1 ; td_cs_cm2( 67) = 1.027103e-20
! td_E_eV( 68) = 5.780234e+1 ; td_cs_cm2( 68) = 1.039879e-20
! td_E_eV( 69) = 5.906424e+1 ; td_cs_cm2( 69) = 1.052273e-20
! td_E_eV( 70) = 6.035369e+1 ; td_cs_cm2( 70) = 1.064281e-20
! td_E_eV( 71) = 6.167129e+1 ; td_cs_cm2( 71) = 1.075895e-20
! td_E_eV( 72) = 6.301765e+1 ; td_cs_cm2( 72) = 1.087112e-20
! td_E_eV( 73) = 6.439341e+1 ; td_cs_cm2( 73) = 1.097926e-20
! td_E_eV( 74) = 6.579920e+1 ; td_cs_cm2( 74) = 1.108335e-20
! td_E_eV( 75) = 6.723568e+1 ; td_cs_cm2( 75) = 1.118334e-20
! td_E_eV( 76) = 6.870352e+1 ; td_cs_cm2( 76) = 1.127920e-20
! td_E_eV( 77) = 7.020341e+1 ; td_cs_cm2( 77) = 1.137092e-20
! td_E_eV( 78) = 7.173604e+1 ; td_cs_cm2( 78) = 1.145848e-20
! td_E_eV( 79) = 7.330214e+1 ; td_cs_cm2( 79) = 1.154185e-20
! td_E_eV( 80) = 7.490242e+1 ; td_cs_cm2( 80) = 1.162104e-20
! td_E_eV( 81) = 7.653763e+1 ; td_cs_cm2( 81) = 1.169604e-20
! td_E_eV( 82) = 7.820855e+1 ; td_cs_cm2( 82) = 1.176684e-20
! td_E_eV( 83) = 7.991594e+1 ; td_cs_cm2( 83) = 1.183346e-20
! td_E_eV( 84) = 8.166061e+1 ; td_cs_cm2( 84) = 1.189590e-20
! td_E_eV( 85) = 8.344337e+1 ; td_cs_cm2( 85) = 1.195416e-20
! td_E_eV( 86) = 8.526505e+1 ; td_cs_cm2( 86) = 1.200827e-20
! td_E_eV( 87) = 8.712649e+1 ; td_cs_cm2( 87) = 1.205825e-20
! td_E_eV( 88) = 8.902858e+1 ; td_cs_cm2( 88) = 1.210411e-20
! td_E_eV( 89) = 9.097219e+1 ; td_cs_cm2( 89) = 1.214588e-20
! td_E_eV( 90) = 9.295823e+1 ; td_cs_cm2( 90) = 1.218359e-20
! td_E_eV( 91) = 9.498763e+1 ; td_cs_cm2( 91) = 1.221728e-20
! td_E_eV( 92) = 9.706133e+1 ; td_cs_cm2( 92) = 1.224696e-20
! td_E_eV( 93) = 9.918031e+1 ; td_cs_cm2( 93) = 1.227270e-20
! td_E_eV( 94) = 1.013455e+2 ; td_cs_cm2( 94) = 1.229451e-20
! td_E_eV( 95) = 1.035580e+2 ; td_cs_cm2( 95) = 1.231245e-20
! td_E_eV( 96) = 1.058189e+2 ; td_cs_cm2( 96) = 1.232656e-20
! td_E_eV( 97) = 1.081290e+2 ; td_cs_cm2( 97) = 1.233689e-20
! td_E_eV( 98) = 1.104896e+2 ; td_cs_cm2( 98) = 1.234348e-20
! td_E_eV( 99) = 1.129017e+2 ; td_cs_cm2( 99) = 1.234640e-20
! td_E_eV(100) = 1.153665e+2 ; td_cs_cm2(100) = 1.234568e-20
! td_E_eV(101) = 1.178851e+2 ; td_cs_cm2(101) = 1.234140e-20
! td_E_eV(102) = 1.204587e+2 ; td_cs_cm2(102) = 1.233360e-20
! td_E_eV(103) = 1.230885e+2 ; td_cs_cm2(103) = 1.232235e-20
! td_E_eV(104) = 1.257757e+2 ; td_cs_cm2(104) = 1.230770e-20
! td_E_eV(105) = 1.285215e+2 ; td_cs_cm2(105) = 1.228972e-20
! td_E_eV(106) = 1.313273e+2 ; td_cs_cm2(106) = 1.226847e-20
! td_E_eV(107) = 1.341944e+2 ; td_cs_cm2(107) = 1.224401e-20
! td_E_eV(108) = 1.371240e+2 ; td_cs_cm2(108) = 1.221643e-20
! td_E_eV(109) = 1.401176e+2 ; td_cs_cm2(109) = 1.218577e-20
! td_E_eV(110) = 1.431765e+2 ; td_cs_cm2(110) = 1.215210e-20
! td_E_eV(111) = 1.463023e+2 ; td_cs_cm2(111) = 1.211551e-20
! td_E_eV(112) = 1.494962e+2 ; td_cs_cm2(112) = 1.207606e-20
! td_E_eV(113) = 1.527599e+2 ; td_cs_cm2(113) = 1.203381e-20
! td_E_eV(114) = 1.560949e+2 ; td_cs_cm2(114) = 1.198884e-20
! td_E_eV(115) = 1.595026e+2 ; td_cs_cm2(115) = 1.194122e-20
! td_E_eV(116) = 1.629848e+2 ; td_cs_cm2(116) = 1.189103e-20
! td_E_eV(117) = 1.665430e+2 ; td_cs_cm2(117) = 1.183833e-20
! td_E_eV(118) = 1.701788e+2 ; td_cs_cm2(118) = 1.178321e-20
! td_E_eV(119) = 1.738940e+2 ; td_cs_cm2(119) = 1.172572e-20
! td_E_eV(120) = 1.776904e+2 ; td_cs_cm2(120) = 1.166596e-20
! td_E_eV(121) = 1.815696e+2 ; td_cs_cm2(121) = 1.160398e-20
! td_E_eV(122) = 1.855335e+2 ; td_cs_cm2(122) = 1.153987e-20
! td_E_eV(123) = 1.895839e+2 ; td_cs_cm2(123) = 1.147369e-20
! td_E_eV(124) = 1.937228e+2 ; td_cs_cm2(124) = 1.140552e-20
! td_E_eV(125) = 1.979520e+2 ; td_cs_cm2(125) = 1.133543e-20
! td_E_eV(126) = 2.022735e+2 ; td_cs_cm2(126) = 1.126350e-20
! td_E_eV(127) = 2.066894e+2 ; td_cs_cm2(127) = 1.118980e-20
! td_E_eV(128) = 2.112017e+2 ; td_cs_cm2(128) = 1.111439e-20
! td_E_eV(129) = 2.158125e+2 ; td_cs_cm2(129) = 1.103735e-20
! td_E_eV(130) = 2.205240e+2 ; td_cs_cm2(130) = 1.095876e-20
! td_E_eV(131) = 2.253383e+2 ; td_cs_cm2(131) = 1.087867e-20
! td_E_eV(132) = 2.302578e+2 ; td_cs_cm2(132) = 1.079717e-20
! td_E_eV(133) = 2.352846e+2 ; td_cs_cm2(133) = 1.071431e-20
! td_E_eV(134) = 2.404212e+2 ; td_cs_cm2(134) = 1.063017e-20
! td_E_eV(135) = 2.456699e+2 ; td_cs_cm2(135) = 1.054482e-20
! td_E_eV(136) = 2.510332e+2 ; td_cs_cm2(136) = 1.045831e-20
! td_E_eV(137) = 2.565135e+2 ; td_cs_cm2(137) = 1.037072e-20
! td_E_eV(138) = 2.621136e+2 ; td_cs_cm2(138) = 1.028211e-20
! td_E_eV(139) = 2.678359e+2 ; td_cs_cm2(139) = 1.019255e-20
! td_E_eV(140) = 2.736831e+2 ; td_cs_cm2(140) = 1.010209e-20
! td_E_eV(141) = 2.796579e+2 ; td_cs_cm2(141) = 1.001079e-20
! td_E_eV(142) = 2.857632e+2 ; td_cs_cm2(142) = 9.918728e-21
! td_E_eV(143) = 2.920018e+2 ; td_cs_cm2(143) = 9.825949e-21
! td_E_eV(144) = 2.983766e+2 ; td_cs_cm2(144) = 9.732516e-21
! td_E_eV(145) = 3.048905e+2 ; td_cs_cm2(145) = 9.638485e-21
! td_E_eV(146) = 3.115467e+2 ; td_cs_cm2(146) = 9.543912e-21
! td_E_eV(147) = 3.183481e+2 ; td_cs_cm2(147) = 9.448852e-21
! td_E_eV(148) = 3.252981e+2 ; td_cs_cm2(148) = 9.353359e-21
! td_E_eV(149) = 3.323998e+2 ; td_cs_cm2(149) = 9.257484e-21
! td_E_eV(150) = 3.396565e+2 ; td_cs_cm2(150) = 9.161280e-21
! td_E_eV(151) = 3.470716e+2 ; td_cs_cm2(151) = 9.064796e-21
! td_E_eV(152) = 3.546487e+2 ; td_cs_cm2(152) = 8.968080e-21
! td_E_eV(153) = 3.623911e+2 ; td_cs_cm2(153) = 8.871181e-21
! td_E_eV(154) = 3.703026e+2 ; td_cs_cm2(154) = 8.774144e-21
! td_E_eV(155) = 3.783868e+2 ; td_cs_cm2(155) = 8.677014e-21
! td_E_eV(156) = 3.866474e+2 ; td_cs_cm2(156) = 8.579836e-21
! td_E_eV(157) = 3.950885e+2 ; td_cs_cm2(157) = 8.482651e-21
! td_E_eV(158) = 4.037137e+2 ; td_cs_cm2(158) = 8.385502e-21
! td_E_eV(159) = 4.125273e+2 ; td_cs_cm2(159) = 8.288427e-21
! td_E_eV(160) = 4.215334e+2 ; td_cs_cm2(160) = 8.191466e-21
! td_E_eV(161) = 4.307360e+2 ; td_cs_cm2(161) = 8.094657e-21
! td_E_eV(162) = 4.401395e+2 ; td_cs_cm2(162) = 7.998035e-21
! td_E_eV(163) = 4.497483e+2 ; td_cs_cm2(163) = 7.901637e-21
! td_E_eV(164) = 4.595669e+2 ; td_cs_cm2(164) = 7.805495e-21
! td_E_eV(165) = 4.695998e+2 ; td_cs_cm2(165) = 7.709643e-21
! td_E_eV(166) = 4.798518e+2 ; td_cs_cm2(166) = 7.614112e-21
! td_E_eV(167) = 4.903276e+2 ; td_cs_cm2(167) = 7.518933e-21
! td_E_eV(168) = 5.010321e+2 ; td_cs_cm2(168) = 7.424135e-21
! td_E_eV(169) = 5.119703e+2 ; td_cs_cm2(169) = 7.329746e-21
! td_E_eV(170) = 5.231473e+2 ; td_cs_cm2(170) = 7.235794e-21
! td_E_eV(171) = 5.345682e+2 ; td_cs_cm2(171) = 7.142303e-21
! td_E_eV(172) = 5.462386e+2 ; td_cs_cm2(172) = 7.049300e-21
! td_E_eV(173) = 5.581637e+2 ; td_cs_cm2(173) = 6.956807e-21
! td_E_eV(174) = 5.703491e+2 ; td_cs_cm2(174) = 6.864834e-21
! td_E_eV(175) = 5.828006e+2 ; td_cs_cm2(175) = 6.773420e-21
! td_E_eV(176) = 5.955239e+2 ; td_cs_cm2(176) = 6.682596e-21
! td_E_eV(177) = 6.085249e+2 ; td_cs_cm2(177) = 6.592380e-21
! td_E_eV(178) = 6.218098e+2 ; td_cs_cm2(178) = 6.502789e-21
! td_E_eV(179) = 6.353847e+2 ; td_cs_cm2(179) = 6.413842e-21
! td_E_eV(180) = 6.492560e+2 ; td_cs_cm2(180) = 6.325553e-21
! td_E_eV(181) = 6.634301e+2 ; td_cs_cm2(181) = 6.237938e-21
! td_E_eV(182) = 6.779136e+2 ; td_cs_cm2(182) = 6.151011e-21
! td_E_eV(183) = 6.927134e+2 ; td_cs_cm2(183) = 6.064786e-21
! td_E_eV(184) = 7.078362e+2 ; td_cs_cm2(184) = 5.979274e-21
! td_E_eV(185) = 7.232892e+2 ; td_cs_cm2(185) = 5.894488e-21
! td_E_eV(186) = 7.390795e+2 ; td_cs_cm2(186) = 5.810439e-21
! td_E_eV(187) = 7.552146e+2 ; td_cs_cm2(187) = 5.727137e-21
! td_E_eV(188) = 7.717019e+2 ; td_cs_cm2(188) = 5.644591e-21
! td_E_eV(189) = 7.885492e+2 ; td_cs_cm2(189) = 5.562811e-21
! td_E_eV(190) = 8.057642e+2 ; td_cs_cm2(190) = 5.481804e-21
! td_E_eV(191) = 8.233551e+2 ; td_cs_cm2(191) = 5.401578e-21
! td_E_eV(192) = 8.413300e+2 ; td_cs_cm2(192) = 5.322140e-21
! td_E_eV(193) = 8.596974e+2 ; td_cs_cm2(193) = 5.243495e-21
! td_E_eV(194) = 8.784657e+2 ; td_cs_cm2(194) = 5.165650e-21
! td_E_eV(195) = 8.976437e+2 ; td_cs_cm2(195) = 5.088609e-21
! td_E_eV(196) = 9.172404e+2 ; td_cs_cm2(196) = 5.012376e-21
! td_E_eV(197) = 9.372650e+2 ; td_cs_cm2(197) = 4.936955e-21
! td_E_eV(198) = 9.577267e+2 ; td_cs_cm2(198) = 4.862350e-21
! td_E_eV(199) = 9.786351e+2 ; td_cs_cm2(199) = 4.788563e-21
! td_E_eV(200) = 1.000000e+3 ; td_cs_cm2(200) = 4.715596e-21
!
! td_cs_cm2 = td_cs_cm2 * 1.0d4 !* 0.5_8  ! convert m^2 to cm^2
!
!!  td_E_eV( 1) =   13.618_8 ; td_cs_cm2( 1) = 0.0
!!  td_E_eV( 2) =   14.0_8   ; td_cs_cm2( 2) = 2.300d-18
!!  td_E_eV( 3) =   16.0_8   ; td_cs_cm2( 3) = 1.720d-17
!!  td_E_eV( 4) =   20.0_8   ; td_cs_cm2( 4) = 3.300d-17
!!  td_E_eV( 5) =   30.0_8   ; td_cs_cm2( 5) = 7.060d-17
!!  td_E_eV( 6) =   50.0_8   ; td_cs_cm2( 6) = 1.106d-16
!!  td_E_eV( 7) =  100.0_8   ; td_cs_cm2( 7) = 1.380d-16
!!  td_E_eV( 8) =  200.0_8   ; td_cs_cm2( 8) = 1.222d-16
!!  td_E_eV( 9) =  500.0_8   ; td_cs_cm2( 9) = 7.920d-17
!!  td_E_eV(10) = 1000.0_8   ; td_cs_cm2(10) = 4.990d-17
!
!  do n = 1, N_td
!     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
!  end do
!
!end subroutine Prepare_O_iontot
!!
!!-------------------------------------
!!
!real(8) function CSV_O_iontot_m3s(E_eV)
!
!  use O_iontot
!  use photoelectrons, only : efactor_eV_to_ms
!
!  implicit none
!
!  real(8) E_eV  ! energy of the electron
!  integer n
!  real(8) left
!
!  real(8) CS_O_iontot_cm2
!
!  if (E_eV.lt.td_E_eV(1)) then
!     CSV_O_iontot_m3s = 0.0_8
!     return
!  end if
!
!  if (E_eV.gt.td_E_eV(N_td)) then
!     CS_O_iontot_cm2 = (1.49d-14 * log(E_eV) - 5.27d-14) / E_eV
!     CSV_O_iontot_m3s = max(0.0_8, 1.0d-4 * CS_O_iontot_cm2 * sqrt(E_eV) * efactor_eV_to_ms)
!     return
!  end if
!
!  do n = 1, N_td-1
!     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
!        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
!        CSV_O_iontot_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
!        return
!     end if
!  end do
!
!end function CSV_O_iontot_m3s





!---------------------------------------------------------------------------------------------------- #46  !39
! O, Rydberg cross-section 
! Table 12 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
module O_Rydberg
  integer, parameter :: N_td=11  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_Rydberg
!
!-------------------------------------
!
subroutine Prepare_O_Rydberg

  use O_Rydberg
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   15.0_8; td_cs_cm2(1)  = 5.56d-18 
  td_E_eV(2)  =   16.0_8; td_cs_cm2(2)  = 1.53d-17
  td_E_eV(3)  =   18.0_8; td_cs_cm2(3)  = 2.34d-17
  td_E_eV(4)  =   20.0_8; td_cs_cm2(4)  = 2.85d-17
  td_E_eV(5)  =   25.0_8; td_cs_cm2(5)  = 3.19d-17
  td_E_eV(6)  =   30.0_8; td_cs_cm2(6)  = 3.14d-17
  td_E_eV(7)  =   50.0_8; td_cs_cm2(7)  = 2.58d-17
  td_E_eV(8)  =  100.0_8; td_cs_cm2(8)  = 1.90d-17
  td_E_eV(9)  =  200.0_8; td_cs_cm2(9)  = 1.28d-17
  td_E_eV(10) =  500.0_8; td_cs_cm2(10) = 7.12d-18
  td_E_eV(11) = 1000.0_8; td_cs_cm2(11) = 4.36d-18

!td_cs_cm2 = 0.5_8 * td_cs_cm2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_Rydberg
!
!-------------------------------------
!
real(8) function CSV_O_Rydberg_m3s(E_eV)

  use O_Rydberg
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O_Rydberg_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_Rydberg_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_Rydberg_m3s

!---------------------------------------------------------------------------------------------------- #47  !40
! O, ^1D cross-section
! Table 12 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
module O_1D
  integer, parameter :: N_td=11  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_1D
!
!-------------------------------------
!
subroutine Prepare_O_1D

  use O_1D
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   2.3_8; td_cs_cm2(1)  = 2.74d-18 
  td_E_eV(2)  =   3.0_8; td_cs_cm2(2)  = 1.55d-17
  td_E_eV(3)  =   4.0_8; td_cs_cm2(3)  = 3.95d-17
  td_E_eV(4)  =   5.5_8; td_cs_cm2(4)  = 5.41d-17
  td_E_eV(5)  =   8.0_8; td_cs_cm2(5)  = 4.59d-17
  td_E_eV(6)  =  10.0_8; td_cs_cm2(6)  = 3.61d-17
  td_E_eV(7)  =  15.0_8; td_cs_cm2(7)  = 2.10d-17
  td_E_eV(8)  =  20.0_8; td_cs_cm2(8)  = 1.33d-17
  td_E_eV(9)  =  30.0_8; td_cs_cm2(9)  = 6.81d-18
  td_E_eV(10) =  50.0_8; td_cs_cm2(10) = 2.90d-18
  td_E_eV(11) = 100.0_8; td_cs_cm2(11) = 8.10d-19

!td_cs_cm2 = 0.5_8 * td_cs_cm2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_1D
!
!-------------------------------------
!
real(8) function CSV_O_1D_m3s(E_eV)

  use O_1D
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O_1D_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_1D_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_1D_m3s

!---------------------------------------------------------------------------------------------------- #48  !41
! O, ^1S cross-section
! Table 12 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
module O_1S
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_1S
!
!-------------------------------------
!
subroutine Prepare_O_1S

  use O_1S
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   4.5_8; td_cs_cm2(1)  = 3.20d-19 
  td_E_eV(2)  =   5.0_8; td_cs_cm2(2)  = 9.00d-19
  td_E_eV(3)  =   6.0_8; td_cs_cm2(3)  = 1.90d-18
  td_E_eV(4)  =   8.0_8; td_cs_cm2(4)  = 3.02d-18
  td_E_eV(5)  =  11.0_8; td_cs_cm2(5)  = 3.36d-18
  td_E_eV(6)  =  15.0_8; td_cs_cm2(6)  = 2.93d-18
  td_E_eV(7)  =  20.0_8; td_cs_cm2(7)  = 2.23d-18
  td_E_eV(8)  =  30.0_8; td_cs_cm2(8)  = 1.31d-18
  td_E_eV(9)  =  50.0_8; td_cs_cm2(9)  = 5.30d-19
  td_E_eV(10) = 100.0_8; td_cs_cm2(10) = 1.40d-19

!td_cs_cm2 = 0.5_8 * td_cs_cm2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_1S
!
!-------------------------------------
!
real(8) function CSV_O_1S_m3s(E_eV)

  use O_1S
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O_1S_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_1S_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_1S_m3s

!---------------------------------------------------------------------------------------------------- #49  !42
! O, ^5P cross-section
! Table 12 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
module O_5P
  integer, parameter :: N_td=9  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_5P
!
!-------------------------------------
!
subroutine Prepare_O_5P

  use O_5P
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1) =  11.0_8; td_cs_cm2(1) = 4.00d-19 
  td_E_eV(2) =  12.0_8; td_cs_cm2(2) = 1.34d-18
  td_E_eV(3) =  14.0_8; td_cs_cm2(3) = 2.23d-18
  td_E_eV(4) =  16.0_8; td_cs_cm2(4) = 2.46d-18
  td_E_eV(5) =  20.0_8; td_cs_cm2(5) = 1.94d-18
  td_E_eV(6) =  25.0_8; td_cs_cm2(6) = 1.15d-18
  td_E_eV(7) =  30.0_8; td_cs_cm2(7) = 7.00d-19
  td_E_eV(8) =  50.0_8; td_cs_cm2(8) = 1.90d-19
  td_E_eV(9) = 100.0_8; td_cs_cm2(9) = 2.00d-20

!td_cs_cm2 = 0.5_8 * td_cs_cm2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_5P
!
!-------------------------------------
!
real(8) function CSV_O_5P_m3s(E_eV)

  use O_5P
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O_5P_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_5P_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_5P_m3s

!---------------------------------------------------------------------------------------------------- #50  !43
! O, ^5S cross-section
! Table 12 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
module O_5S
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_5S
!
!-------------------------------------
!
subroutine Prepare_O_5S

  use O_5S
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =  9.5_8; td_cs_cm2(1)  = 7.70d-19 
  td_E_eV(2)  = 10.0_8; td_cs_cm2(2)  = 1.40d-18
  td_E_eV(3)  = 11.0_8; td_cs_cm2(3)  = 2.05d-18
  td_E_eV(4)  = 12.0_8; td_cs_cm2(4)  = 2.61d-18
  td_E_eV(5)  = 14.0_8; td_cs_cm2(5)  = 3.07d-18
  td_E_eV(6)  = 17.0_8; td_cs_cm2(6)  = 2.10d-18
  td_E_eV(7)  = 20.0_8; td_cs_cm2(7)  = 1.19d-18
  td_E_eV(8)  = 25.0_8; td_cs_cm2(8)  = 5.80d-19
  td_E_eV(9)  = 30.0_8; td_cs_cm2(9)  = 2.90d-19
  td_E_eV(10) = 50.0_8; td_cs_cm2(10) = 4.00d-20

!td_cs_cm2 = 0.5_8 * td_cs_cm2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_5S
!
!-------------------------------------
!
real(8) function CSV_O_5S_m3s(E_eV)

  use O_5S
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O_5S_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_5S_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_5S_m3s

!---------------------------------------------------------------------------------------------------- #51  !44
! O, 3s^3S cross-section
! Table 13 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
module O_3s3S
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_3s3S
!
!-------------------------------------
!
subroutine Prepare_O_3s3S

  use O_3s3S
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   10.0_8; td_cs_cm2(1)  = 3.02d-18 
  td_E_eV(2)  =   12.0_8; td_cs_cm2(2)  = 7.32d-18
  td_E_eV(3)  =   15.0_8; td_cs_cm2(3)  = 9.98d-18
  td_E_eV(4)  =   20.0_8; td_cs_cm2(4)  = 1.11d-17
  td_E_eV(5)  =   30.0_8; td_cs_cm2(5)  = 1.05d-17
  td_E_eV(6)  =   50.0_8; td_cs_cm2(6)  = 8.36d-18
  td_E_eV(7)  =  100.0_8; td_cs_cm2(7)  = 5.57d-18
  td_E_eV(8)  =  200.0_8; td_cs_cm2(8)  = 3.69d-18
  td_E_eV(9)  =  500.0_8; td_cs_cm2(9)  = 2.06d-18
  td_E_eV(10) = 1000.0_8; td_cs_cm2(10) = 1.31d-18

!td_cs_cm2 = 0.5_8 * td_cs_cm2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_3s3S
!
!-------------------------------------
!
real(8) function CSV_O_3s3S_m3s(E_eV)

  use O_3s3S
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O_3s3S_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_3s3S_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_3s3S_m3s

!---------------------------------------------------------------------------------------------------- #52  !45
! O, 3d^3D cross-section
! Table 13 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
module O_3d3D
  integer, parameter :: N_td=11  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_3d3D
!
!-------------------------------------
!
subroutine Prepare_O_3d3D

  use O_3d3D
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   12.0_8; td_cs_cm2(1)  = 1.00d-19 
  td_E_eV(2)  =   13.0_8; td_cs_cm2(2)  = 3.70d-19
  td_E_eV(3)  =   14.0_8; td_cs_cm2(3)  = 5.70d-19
  td_E_eV(4)  =   16.0_8; td_cs_cm2(4)  = 9.10d-19
  td_E_eV(5)  =   20.0_8; td_cs_cm2(5)  = 1.70d-18
  td_E_eV(6)  =   30.0_8; td_cs_cm2(6)  = 2.93d-18
  td_E_eV(7)  =   50.0_8; td_cs_cm2(7)  = 3.47d-18
  td_E_eV(8)  =  100.0_8; td_cs_cm2(8)  = 2.32d-18
  td_E_eV(9)  =  200.0_8; td_cs_cm2(9)  = 1.25d-18
  td_E_eV(10) =  500.0_8; td_cs_cm2(10) = 5.50d-19
  td_E_eV(11) = 1000.0_8; td_cs_cm2(11) = 3.00d-19

!td_cs_cm2 = 0.5_8 * td_cs_cm2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_3d3D
!
!-------------------------------------
!
real(8) function CSV_O_3d3D_m3s(E_eV)

  use O_3d3D
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O_3d3D_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_3d3D_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_3d3D_m3s

!---------------------------------------------------------------------------------------------------- #53  !46
! O, 3s'^3D cross-section
! Table 13 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
module O_3sp3D
  integer, parameter :: N_td=12  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_3sp3D
!
!-------------------------------------
!
subroutine Prepare_O_3sp3D

  use O_3sp3D
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   13.0_8; td_cs_cm2(1)  = 1.40d-19 
  td_E_eV(2)  =   14.0_8; td_cs_cm2(2)  = 2.80d-19
  td_E_eV(3)  =   16.0_8; td_cs_cm2(3)  = 5.70d-19
  td_E_eV(4)  =   18.0_8; td_cs_cm2(4)  = 1.20d-18
  td_E_eV(5)  =   20.0_8; td_cs_cm2(5)  = 1.87d-18
  td_E_eV(6)  =   25.0_8; td_cs_cm2(6)  = 4.15d-18
  td_E_eV(7)  =   30.0_8; td_cs_cm2(7)  = 5.32d-18
  td_E_eV(8)  =   50.0_8; td_cs_cm2(8)  = 6.42d-18
  td_E_eV(9)  =  100.0_8; td_cs_cm2(9)  = 4.70d-18
  td_E_eV(10) =  200.0_8; td_cs_cm2(10) = 2.88d-18
  td_E_eV(11) =  500.0_8; td_cs_cm2(11) = 1.33d-18
  td_E_eV(12) = 1000.0_8; td_cs_cm2(12) = 7.00d-19

!td_cs_cm2 = 0.5_8 * td_cs_cm2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_3sp3D
!
!-------------------------------------
!
real(8) function CSV_O_3sp3D_m3s(E_eV)

  use O_3sp3D
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O_3sp3D_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_3sp3D_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_3sp3D_m3s

!---------------------------------------------------------------------------------------------------- #54  !47
! O, 3p^3P cross-section
! Table 13 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
module O_3p3P
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_3p3P
!
!-------------------------------------
!
subroutine Prepare_O_3p3P

  use O_3p3P
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =  12.0_8; td_cs_cm2(1)  = 1.20d-18 
  td_E_eV(2)  =  14.0_8; td_cs_cm2(2)  = 3.44d-18
  td_E_eV(3)  =  16.0_8; td_cs_cm2(3)  = 5.55d-18
  td_E_eV(4)  =  20.0_8; td_cs_cm2(4)  = 7.55d-18
  td_E_eV(5)  =  25.0_8; td_cs_cm2(5)  = 6.41d-18
  td_E_eV(6)  =  30.0_8; td_cs_cm2(6)  = 5.17d-18
  td_E_eV(7)  =  50.0_8; td_cs_cm2(7)  = 2.75d-18
  td_E_eV(8)  = 100.0_8; td_cs_cm2(8)  = 1.14d-18
  td_E_eV(9)  = 200.0_8; td_cs_cm2(9)  = 5.00d-19
  td_E_eV(10) = 500.0_8; td_cs_cm2(10) = 1.40d-19

!td_cs_cm2 = 0.5_8 * td_cs_cm2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_3p3P
!
!-------------------------------------
!
real(8) function CSV_O_3p3P_m3s(E_eV)

  use O_3p3P
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O_3p3P_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_3p3P_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_3p3P_m3s

!---------------------------------------------------------------------------------------------------- #55  !48
! O, 5d^3D cross-section
! Table 13 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
module O_5d3D
  integer, parameter :: N_td=9  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_5d3D
!
!-------------------------------------
!
subroutine Prepare_O_5d3D

  use O_5d3D
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1) =   14.0_8; td_cs_cm2(1) = 2.40d-19 
  td_E_eV(2) =   16.0_8; td_cs_cm2(2) = 4.30d-19
  td_E_eV(3) =   20.0_8; td_cs_cm2(3) = 6.50d-19
  td_E_eV(4) =   30.0_8; td_cs_cm2(4) = 1.00d-18
  td_E_eV(5) =   50.0_8; td_cs_cm2(5) = 1.13d-18
  td_E_eV(6) =  100.0_8; td_cs_cm2(6) = 7.30d-19
  td_E_eV(7) =  200.0_8; td_cs_cm2(7) = 4.00d-19
  td_E_eV(8) =  500.0_8; td_cs_cm2(8) = 2.00d-19
  td_E_eV(9) = 1000.0_8; td_cs_cm2(9) = 1.00d-19

!td_cs_cm2 = 0.5_8 * td_cs_cm2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_5d3D
!
!-------------------------------------
!
real(8) function CSV_O_5d3D_m3s(E_eV)

  use O_5d3D
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O_5d3D_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_5d3D_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_5d3D_m3s

!---------------------------------------------------------------------------------------------------- #56  !49
! O, 4d^3D cross-section
! Table 13 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
module O_4d3D
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_4d3D
!
!-------------------------------------
!
subroutine Prepare_O_4d3D

  use O_4d3D
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   13.0_8; td_cs_cm2(1)  = 1.50d-19 
  td_E_eV(2)  =   14.0_8; td_cs_cm2(2)  = 3.30d-19
  td_E_eV(3)  =   16.0_8; td_cs_cm2(3)  = 5.70d-19
  td_E_eV(4)  =   20.0_8; td_cs_cm2(4)  = 1.00d-18
  td_E_eV(5)  =   30.0_8; td_cs_cm2(5)  = 1.68d-18
  td_E_eV(6)  =   50.0_8; td_cs_cm2(6)  = 1.99d-18
  td_E_eV(7)  =  100.0_8; td_cs_cm2(7)  = 1.28d-18
  td_E_eV(8)  =  200.0_8; td_cs_cm2(8)  = 7.10d-19
  td_E_eV(9)  =  500.0_8; td_cs_cm2(9)  = 3.20d-19
  td_E_eV(10) = 1000.0_8; td_cs_cm2(10) = 1.80d-19

!td_cs_cm2 = 0.5_8 * td_cs_cm2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_4d3D
!
!-------------------------------------
!
real(8) function CSV_O_4d3D_m3s(E_eV)

  use O_4d3D
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O_4d3D_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_4d3D_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_4d3D_m3s

!---------------------------------------------------------------------------------------------------- #57  !50
! O, 2p^5 ^3P cross-section
! Table 14 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
module O_2p53P
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_2p53P
!
!-------------------------------------
!
subroutine Prepare_O_2p53P

  use O_2p53P
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   16.0_8; td_cs_cm2(1)  = 2.00d-19 
  td_E_eV(2)  =   18.0_8; td_cs_cm2(2)  = 6.90d-19
  td_E_eV(3)  =   20.0_8; td_cs_cm2(3)  = 1.23d-18
  td_E_eV(4)  =   25.0_8; td_cs_cm2(4)  = 3.71d-18
  td_E_eV(5)  =   30.0_8; td_cs_cm2(5)  = 6.48d-18
  td_E_eV(6)  =   50.0_8; td_cs_cm2(6)  = 1.15d-17
  td_E_eV(7)  =  100.0_8; td_cs_cm2(7)  = 8.28d-18
  td_E_eV(8)  =  200.0_8; td_cs_cm2(8)  = 4.20d-18
  td_E_eV(9)  =  500.0_8; td_cs_cm2(9)  = 1.20d-18
  td_E_eV(10) = 1000.0_8; td_cs_cm2(10) = 4.10d-19

!td_cs_cm2 = 0.5_8 * td_cs_cm2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_2p53P
!
!-------------------------------------
!
real(8) function CSV_O_2p53P_m3s(E_eV)

  use O_2p53P
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O_2p53P_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_2p53P_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_2p53P_m3s

!---------------------------------------------------------------------------------------------------- #58  !51
! O, 3s'' ^3P cross-section
! Table 14 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
module O_3spp3P
  integer, parameter :: N_td=10  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_3spp3P
!
!-------------------------------------
!
subroutine Prepare_O_3spp3P

  use O_3spp3P
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   16.0_8; td_cs_cm2(1)  = 3.50d-19 
  td_E_eV(2)  =   18.0_8; td_cs_cm2(2)  = 8.40d-19
  td_E_eV(3)  =   20.0_8; td_cs_cm2(3)  = 1.37d-18
  td_E_eV(4)  =   25.0_8; td_cs_cm2(4)  = 3.60d-18
  td_E_eV(5)  =   30.0_8; td_cs_cm2(5)  = 6.32d-18
  td_E_eV(6)  =   50.0_8; td_cs_cm2(6)  = 1.10d-17
  td_E_eV(7)  =  100.0_8; td_cs_cm2(7)  = 7.94d-18
  td_E_eV(8)  =  200.0_8; td_cs_cm2(8)  = 4.00d-18
  td_E_eV(9)  =  500.0_8; td_cs_cm2(9)  = 1.14d-18
  td_E_eV(10) = 1000.0_8; td_cs_cm2(10) = 3.90d-19

!td_cs_cm2 = 0.5_8 * td_cs_cm2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_3spp3P
!
!-------------------------------------
!
real(8) function CSV_O_3spp3P_m3s(E_eV)

  use O_3spp3P
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O_3spp3P_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_3spp3P_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_3spp3P_m3s

!---------------------------------------------------------------------------------------------------- #59  !52
! O, 4d' ^3P cross-section
! Table 14 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
module O_4dp3P
  integer, parameter :: N_td=9  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_4dp3P
!
!-------------------------------------
!
subroutine Prepare_O_4dp3P

  use O_4dp3P
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1) =  18.0_8; td_cs_cm2(1) = 1.90d-19 
  td_E_eV(2) =  20.0_8; td_cs_cm2(2) = 3.80d-19
  td_E_eV(3) =  25.0_8; td_cs_cm2(3) = 8.50d-19
  td_E_eV(4) =  30.0_8; td_cs_cm2(4) = 1.40d-18
  td_E_eV(5) =  50.0_8; td_cs_cm2(5) = 2.64d-18
  td_E_eV(6) = 100.0_8; td_cs_cm2(6) = 1.78d-18
  td_E_eV(7) = 200.0_8; td_cs_cm2(7) = 8.10d-19
  td_E_eV(8) = 500.0_8; td_cs_cm2(8) = 2.10d-19
  td_E_eV(9) = 800.0_8; td_cs_cm2(9) = 1.00d-19

!td_cs_cm2 = 0.5_8 * td_cs_cm2

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_4dp3P
!
!-------------------------------------
!
real(8) function CSV_O_4dp3P_m3s(E_eV)

  use O_4dp3P
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_O_4dp3P_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_4dp3P_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_4dp3P_m3s
