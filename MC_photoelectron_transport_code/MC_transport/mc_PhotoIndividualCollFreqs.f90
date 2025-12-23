
!---------------------------------------------------------------------------------------------------- #1
! N2, elastic scattering cross section [NOT the momentum transfer cross section]  
! Table 3 of [Y.Itikawa, J.Phys.Chem.Ref.Data, v.35, p.31-53, 2006]! 
!
module N2_elastic
  integer, parameter :: N_td=35  ! number of table data points
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

  td_E_eV(1)   =    0.55_8  ; td_cs_cm2(1)   = 8.390d-16 
  td_E_eV(2)   =    0.70_8  ; td_cs_cm2(2)   = 9.030d-16
  td_E_eV(3)   =    0.90_8  ; td_cs_cm2(3)   = 9.620d-16
  td_E_eV(4)   =    1.00_8  ; td_cs_cm2(4)   = 9.830d-16
  td_E_eV(5)   =    1.50_8  ; td_cs_cm2(5)   = 1.053d-15
  td_E_eV(6)   =    2.00_8  ; td_cs_cm2(6)   = 1.793d-15
  td_E_eV(7)   =    2.20_8  ; td_cs_cm2(7)   = 1.950d-15
  td_E_eV(8)   =    2.35_8  ; td_cs_cm2(8)   = 2.050d-15
  td_E_eV(9)   =    2.50_8  ; td_cs_cm2(9)   = 2.100d-15
  td_E_eV(10)  =    2.70_8  ; td_cs_cm2(10)  = 1.750d-15
  td_E_eV(11)  =    3.00_8  ; td_cs_cm2(11)  = 1.500d-15
  td_E_eV(12)  =    4.00_8  ; td_cs_cm2(12)  = 1.160d-15
  td_E_eV(13)  =    5.00_8  ; td_cs_cm2(13)  = 1.075d-15
  td_E_eV(14)  =    6.00_8  ; td_cs_cm2(14)  = 1.060d-15
  td_E_eV(15)  =    8.00_8  ; td_cs_cm2(15)  = 1.060d-15
  td_E_eV(16)  =   10.00_8  ; td_cs_cm2(16)  = 1.140d-15
  td_E_eV(17)  =   15.00_8  ; td_cs_cm2(17)  = 1.180d-15
  td_E_eV(18)  =   20.00_8  ; td_cs_cm2(18)  = 1.115d-15
  td_E_eV(19)  =   25.00_8  ; td_cs_cm2(19)  = 1.025d-15
  td_E_eV(20)  =   30.00_8  ; td_cs_cm2(20)  = 9.650d-16
  td_E_eV(21)  =   40.00_8  ; td_cs_cm2(21)  = 8.850d-16
  td_E_eV(22)  =   50.00_8  ; td_cs_cm2(22)  = 8.200d-16
  td_E_eV(23)  =   60.00_8  ; td_cs_cm2(23)  = 7.400d-16
  td_E_eV(24)  =   80.00_8  ; td_cs_cm2(24)  = 6.250d-16
  td_E_eV(25)  =  100.00_8  ; td_cs_cm2(25)  = 5.600d-16
  td_E_eV(26)  =  120.00_8  ; td_cs_cm2(26)  = 4.900d-16
  td_E_eV(27)  =  150.00_8  ; td_cs_cm2(27)  = 4.200d-16
  td_E_eV(28)  =  200.00_8  ; td_cs_cm2(28)  = 3.500d-16
  td_E_eV(29)  =  250.00_8  ; td_cs_cm2(29)  = 3.000d-16
  td_E_eV(30)  =  300.00_8  ; td_cs_cm2(30)  = 2.650d-16
  td_E_eV(31)  =  400.00_8  ; td_cs_cm2(31)  = 2.150d-16
  td_E_eV(32)  =  500.00_8  ; td_cs_cm2(32)  = 1.850d-16
  td_E_eV(33)  =  600.00_8  ; td_cs_cm2(33)  = 1.600d-16
  td_E_eV(34)  =  800.00_8  ; td_cs_cm2(34)  = 1.250d-16
  td_E_eV(35)  = 1000.00_8  ; td_cs_cm2(35)  = 1.000d-16
  
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
! N2, total vibrational cross-section 
! Table 4 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997] 
! threhold is 0.9 eV, see table 1 of [Majeed and Strickland] 
!
module N2_vibtot
  integer, parameter :: N_td=8   ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module N2_vibtot
!
!-------------------------------------
!
subroutine Prepare_N2_vibtot

  use N2_vibtot
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1) = 1.3_8; td_cs_cm2(1) = 1.11d-18
  td_E_eV(2) = 1.6_8; td_cs_cm2(2) = 8.84d-17
  td_E_eV(3) = 2.0_8; td_cs_cm2(3) = 1.01d-15
  td_E_eV(4) = 2.4_8; td_cs_cm2(4) = 1.35d-15
  td_E_eV(5) = 3.0_8; td_cs_cm2(5) = 5.86d-16
  td_E_eV(6) = 4.0_8; td_cs_cm2(6) = 1.37d-16
  td_E_eV(7) = 5.0_8; td_cs_cm2(7) = 2.34d-17
  td_E_eV(8) = 6.0_8; td_cs_cm2(8) = 3.66d-18

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_N2_vibtot
!
!-------------------------------------
!
real(8) function CSV_N2_vibtot_m3s(E_eV)

  use N2_vibtot
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_td))) then
     CSV_N2_vibtot_m3s = 0.0_8
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_N2_vibtot_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_N2_vibtot_m3s

!---------------------------------------------------------------------------------------------------- #5
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

!---------------------------------------------------------------------------------------------------- #6
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

!---------------------------------------------------------------------------------------------------- #7
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
  td_E_eV(13) = 150.0_8; td_cs_cm2(13) = 3.10d-18
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

!---------------------------------------------------------------------------------------------------- #8
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

!---------------------------------------------------------------------------------------------------- #9
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

!---------------------------------------------------------------------------------------------------- #10
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

!---------------------------------------------------------------------------------------------------- #11
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

!---------------------------------------------------------------------------------------------------- #12
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

!---------------------------------------------------------------------------------------------------- #13
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

!---------------------------------------------------------------------------------------------------- #14
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

!---------------------------------------------------------------------------------------------------- #15
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

!---------------------------------------------------------------------------------------------------- #16
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

!---------------------------------------------------------------------------------------------------- #17
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

!---------------------------------------------------------------------------------------------------- #18
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

!---------------------------------------------------------------------------------------------------- #19
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

!---------------------------------------------------------------------------------------------------- #20
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

!---------------------------------------------------------------------------------------------------- #21
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

!---------------------------------------------------------------------------------------------------- #22
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

!---------------------------------------------------------------------------------------------------- #23
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

!---------------------------------------------------------------------------------------------------- #24
! O2, elastic scattering cross-section 
! Table 2 of [Y.Itikawa, J.Phys.Chem.Ref.Data, v.38, p.1-20, 2009] 
!
module O2_elastic
  integer, parameter :: N_td=30  ! number of table data points
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

  td_E_eV(1)   =    1.0_8  ; td_cs_cm2(1)   = 5.97d-16
  td_E_eV(2)   =    2.0_8  ; td_cs_cm2(2)   = 6.45d-16
  td_E_eV(3)   =    3.0_8  ; td_cs_cm2(3)   = 6.74d-16
  td_E_eV(4)   =    4.0_8  ; td_cs_cm2(4)   = 6.93d-16
  td_E_eV(5)   =    5.0_8  ; td_cs_cm2(5)   = 7.20d-16
  td_E_eV(6)   =    6.0_8  ; td_cs_cm2(6)   = 7.52d-16
  td_E_eV(7)   =    7.0_8  ; td_cs_cm2(7)   = 7.86d-16
  td_E_eV(8)   =    8.0_8  ; td_cs_cm2(8)   = 8.21d-16
  td_E_eV(9)   =    9.0_8  ; td_cs_cm2(9)   = 8.49d-16
  td_E_eV(10)  =   10.0_8  ; td_cs_cm2(10)  = 8.80d-16
  td_E_eV(11)  =   12.0_8  ; td_cs_cm2(11)  = 9.00d-16
  td_E_eV(12)  =   15.0_8  ; td_cs_cm2(12)  = 8.89d-16
  td_E_eV(13)  =   20.0_8  ; td_cs_cm2(13)  = 8.60d-16
  td_E_eV(14)  =   30.0_8  ; td_cs_cm2(14)  = 8.09d-16
  td_E_eV(15)  =   40.0_8  ; td_cs_cm2(15)  = 7.30d-16
  td_E_eV(16)  =   50.0_8  ; td_cs_cm2(16)  = 6.59d-16
  td_E_eV(17)  =   60.0_8  ; td_cs_cm2(17)  = 6.08d-16
  td_E_eV(18)  =   70.0_8  ; td_cs_cm2(18)  = 5.63d-16
  td_E_eV(19)  =   80.0_8  ; td_cs_cm2(19)  = 5.29d-16
  td_E_eV(20)  =   90.0_8  ; td_cs_cm2(20)  = 5.01d-16
  td_E_eV(21)  =  100.0_8  ; td_cs_cm2(21)  = 4.78d-16
  td_E_eV(22)  =  200.0_8  ; td_cs_cm2(22)  = 3.15d-16
  td_E_eV(23)  =  300.0_8  ; td_cs_cm2(23)  = 2.40d-16
  td_E_eV(24)  =  400.0_8  ; td_cs_cm2(24)  = 2.00d-16
  td_E_eV(25)  =  500.0_8  ; td_cs_cm2(25)  = 1.72d-16
  td_E_eV(26)  =  600.0_8  ; td_cs_cm2(26)  = 1.53d-16
  td_E_eV(27)  =  700.0_8  ; td_cs_cm2(27)  = 1.37d-16
  td_E_eV(28)  =  800.0_8  ; td_cs_cm2(28)  = 1.27d-16
  td_E_eV(29)  =  900.0_8  ; td_cs_cm2(29)  = 1.18d-16
  td_E_eV(30)  = 1000.0_8  ; td_cs_cm2(30)  = 1.10d-16
  
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

!---------------------------------------------------------------------------------------------------- #25
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

!---------------------------------------------------------------------------------------------------- #26
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

!---------------------------------------------------------------------------------------------------- #27
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

!---------------------------------------------------------------------------------------------------- #28
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

!---------------------------------------------------------------------------------------------------- #29
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

!---------------------------------------------------------------------------------------------------- #30
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

!---------------------------------------------------------------------------------------------------- #31
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

!---------------------------------------------------------------------------------------------------- #32
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

!---------------------------------------------------------------------------------------------------- #33
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

!---------------------------------------------------------------------------------------------------- #34
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

!---------------------------------------------------------------------------------------------------- #35
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

!---------------------------------------------------------------------------------------------------- #36
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

!---------------------------------------------------------------------------------------------------- #37
! I could not find table data however I found the following paper
! Y. Itikawa and A. Ichimura, "Cross Sections for Collisions and Photons with Atomic Oxygen"
! J. Phys. Chem. Ref. Data, Vol.19, No.3, pp.637-651, 1990
! where Figures 4.1 and 4.2 show curves of the elastic cross sections
! so the data below are retrieved manually from these figures.
!
module O_elastic
  integer, parameter :: N_td=28  ! number of table data points
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

  td_E_eV(1)   =    1.0_8  ; td_cs_cm2(1)   = 4.2e-16
  td_E_eV(2)   =    2.0_8  ; td_cs_cm2(2)   = 6.0e-16
  td_E_eV(3)   =    3.0_8  ; td_cs_cm2(3)   = 6.7e-16
  td_E_eV(4)   =    4.0_8  ; td_cs_cm2(4)   = 7.0e-16
  td_E_eV(5)   =    5.0_8  ; td_cs_cm2(5)   = 7.3e-16
  td_E_eV(6)   =    6.0_8  ; td_cs_cm2(6)   = 7.3e-16
  td_E_eV(7)   =    7.0_8  ; td_cs_cm2(7)   = 7.3e-16
  td_E_eV(8)   =    8.0_8  ; td_cs_cm2(8)   = 7.3e-16
  td_E_eV(9)   =    9.0_8  ; td_cs_cm2(9)   = 7.2e-16
  td_E_eV(10)  =   10.0_8  ; td_cs_cm2(10)  = 7.1e-16
  td_E_eV(11)  =   20.0_8  ; td_cs_cm2(11)  = 6.0e-16
  td_E_eV(12)  =   30.0_8  ; td_cs_cm2(12)  = 5.3e-16
  td_E_eV(13)  =   40.0_8  ; td_cs_cm2(13)  = 4.6e-16
  td_E_eV(14)  =   50.0_8  ; td_cs_cm2(14)  = 4.1e-16
  td_E_eV(15)  =   60.0_8  ; td_cs_cm2(15)  = 3.7e-16
  td_E_eV(16)  =   70.0_8  ; td_cs_cm2(16)  = 3.4e-16
  td_E_eV(17)  =   80.0_8  ; td_cs_cm2(17)  = 3.2e-16
  td_E_eV(18)  =   90.0_8  ; td_cs_cm2(18)  = 2.9e-16
  td_E_eV(19)  =  100.0_8  ; td_cs_cm2(19)  = 2.7e-16
  td_E_eV(20)  =  200.0_8  ; td_cs_cm2(20)  = 1.7e-16
  td_E_eV(21)  =  300.0_8  ; td_cs_cm2(21)  = 1.3e-16
  td_E_eV(22)  =  400.0_8  ; td_cs_cm2(22)  = 1.05e-16
  td_E_eV(23)  =  500.0_8  ; td_cs_cm2(23)  = 9.0e-17
  td_E_eV(24)  =  600.0_8  ; td_cs_cm2(24)  = 7.8e-17
  td_E_eV(25)  =  700.0_8  ; td_cs_cm2(25)  = 7.0e-17
  td_E_eV(26)  =  800.0_8  ; td_cs_cm2(26)  = 6.2e-17
  td_E_eV(27)  =  900.0_8  ; td_cs_cm2(27)  = 5.8e-17
  td_E_eV(28)  = 1000.0_8  ; td_cs_cm2(28)  = 5.3e-17
  
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

!---------------------------------------------------------------------------------------------------- #38
! O, total {O->(O+,O++(?))} ionization cross-section 
! Table 12 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
module O_iontot
  integer, parameter :: N_td=9  ! number of table data points
! table data
  real(8) td_E_eV(1:N_td)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_td)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_td)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module O_iontot
!
!-------------------------------------
!
subroutine Prepare_O_iontot

  use O_iontot
  use photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1) =   14.0_8; td_cs_cm2(1) = 2.300d-18 
  td_E_eV(2) =   16.0_8; td_cs_cm2(2) = 1.720d-17
  td_E_eV(3) =   20.0_8; td_cs_cm2(3) = 3.300d-17
  td_E_eV(4) =   30.0_8; td_cs_cm2(4) = 7.060d-17
  td_E_eV(5) =   50.0_8; td_cs_cm2(5) = 1.106d-16
  td_E_eV(6) =  100.0_8; td_cs_cm2(6) = 1.380d-16
  td_E_eV(7) =  200.0_8; td_cs_cm2(7) = 1.222d-16
  td_E_eV(8) =  500.0_8; td_cs_cm2(8) = 7.920d-17
  td_E_eV(9) = 1000.0_8; td_cs_cm2(9) = 4.990d-17

  do n = 1, N_td
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_O_iontot
!
!-------------------------------------
!
real(8) function CSV_O_iontot_m3s(E_eV)

  use O_iontot
  use photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  real(8) CS_O_iontot_cm2

  if (E_eV.lt.td_E_eV(1)) then
     CSV_O_iontot_m3s = 0.0_8
     return
  end if

  if (E_eV.gt.td_E_eV(N_td)) then
     CS_O_iontot_cm2 = (1.49d-14 * log(E_eV) - 5.27d-14) / E_eV
     CSV_O_iontot_m3s = max(0.0_8, 1.0d-4 * CS_O_iontot_cm2 * sqrt(E_eV) * efactor_eV_to_ms)
     return
  end if

  do n = 1, N_td-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_O_iontot_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_O_iontot_m3s

!---------------------------------------------------------------------------------------------------- #39
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

!---------------------------------------------------------------------------------------------------- #40
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

!---------------------------------------------------------------------------------------------------- #41
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

!---------------------------------------------------------------------------------------------------- #42
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

!---------------------------------------------------------------------------------------------------- #43
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

!---------------------------------------------------------------------------------------------------- #44
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

!---------------------------------------------------------------------------------------------------- #45
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

!---------------------------------------------------------------------------------------------------- #46
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

!---------------------------------------------------------------------------------------------------- #47
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

!---------------------------------------------------------------------------------------------------- #48
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

!---------------------------------------------------------------------------------------------------- #49
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

!---------------------------------------------------------------------------------------------------- #50
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

!---------------------------------------------------------------------------------------------------- #51
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

!---------------------------------------------------------------------------------------------------- #52
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

!---------------------------------------------------------------------------------------------------- #53
! He, elastic scattering cross section (not the momentum transfer)
!
! 0.5-18 eV is Table II of R. E. Kennerly and R. A. Bonham, Phys.Rev.A, vol.17, pp.1844-1854, 1978
! full energy range in this table is 0.5-50 eV
! note: Table II shows "electron-He total cross section", NOT total elastic cross-section as table X of Nesbet
! where it is compared with (and is very close to) total elastic cross-section calculated by Nesbet
! [R. K. Nesbet, "Variational calculations of accurate e-He cross sections below 19 eV", v.20, pp.58-70, 1979)
!
! 20 eV - 1000 eV is from table III of D. F. Register, S. Trajmar, and S. K. Srivastava, Phys.Rev.A, v.21, pp.1134-1151, 1980
!
module He_elastic
  integer, parameter :: N_tdp=37  ! number of table data points
! table data
  real(8) td_E_eV(1:N_tdp)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_tdp)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_tdp)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module He_elastic
!
!-------------------------------------
!
subroutine Prepare_He_elastic

  use He_elastic
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =    0.5_8 ; td_cs_cm2(1)  =   6.24d-16
  td_E_eV(2)  =    1.0_8 ; td_cs_cm2(2)  =   6.23d-16
  td_E_eV(3)  =    1.5_8 ; td_cs_cm2(3)  =   6.18d-16
  td_E_eV(4)  =    2.0_8 ; td_cs_cm2(4)  =   6.06d-16
  td_E_eV(5)  =    2.5_8 ; td_cs_cm2(5)  =   5.92d-16
  td_E_eV(6)  =    3.0_8 ; td_cs_cm2(6)  =   5.78d-16
  td_E_eV(7)  =    4.0_8 ; td_cs_cm2(7)  =   5.50d-16
  td_E_eV(8)  =    5.0_8 ; td_cs_cm2(8)  =   5.25d-16
  td_E_eV(9)  =    6.0_8 ; td_cs_cm2(9)  =   5.04d-16
  td_E_eV(10) =    7.0_8 ; td_cs_cm2(10) =   4.83d-16
  td_E_eV(11) =    8.0_8 ; td_cs_cm2(11) =   4.64d-16
  td_E_eV(12) =    9.0_8 ; td_cs_cm2(12) =   4.46d-16
  td_E_eV(13) =   10.0_8 ; td_cs_cm2(13) =   4.30d-16
  td_E_eV(14) =   12.0_8 ; td_cs_cm2(14) =   3.96d-16
  td_E_eV(15) =   14.0_8 ; td_cs_cm2(15) =   3.69d-16
  td_E_eV(16) =   16.0_8 ; td_cs_cm2(16) =   3.43d-16
  td_E_eV(17) =   18.0_8 ; td_cs_cm2(17) =   3.22d-16

  td_E_eV(18) =   20.0_8 ; td_cs_cm2(18) =   3.00d-16
  td_E_eV(19) =   25.0_8 ; td_cs_cm2(19) =   2.51d-16
  td_E_eV(20) =   30.0_8 ; td_cs_cm2(20) =   2.11d-16
  td_E_eV(21) =   40.0_8 ; td_cs_cm2(21) =   1.58d-16
  td_E_eV(22) =   50.0_8 ; td_cs_cm2(22) =   1.26d-16
  td_E_eV(23) =   60.0_8 ; td_cs_cm2(23) =   1.02d-16
  td_E_eV(24) =   70.0_8 ; td_cs_cm2(24) =   0.85d-16
  td_E_eV(25) =   75.0_8 ; td_cs_cm2(25) =   0.77d-16
  td_E_eV(26) =   80.0_8 ; td_cs_cm2(26) =   0.71d-16
  td_E_eV(27) =   90.0_8 ; td_cs_cm2(27) =   0.62d-16
  td_E_eV(28) =  100.0_8 ; td_cs_cm2(28) =   0.56d-16
  td_E_eV(29) =  150.0_8 ; td_cs_cm2(29) =   0.35d-16
  td_E_eV(30) =  200.0_8 ; td_cs_cm2(30) =   0.25d-16
  td_E_eV(31) =  300.0_8 ; td_cs_cm2(31) =   0.155d-16
  td_E_eV(32) =  400.0_8 ; td_cs_cm2(32) =   0.115d-16
  td_E_eV(33) =  500.0_8 ; td_cs_cm2(33) =   0.084d-16
  td_E_eV(34) =  600.0_8 ; td_cs_cm2(34) =   0.076d-16
  td_E_eV(35) =  700.0_8 ; td_cs_cm2(35) =   0.065d-16
  td_E_eV(36) =  800.0_8 ; td_cs_cm2(36) =   0.060d-16
  td_E_eV(37) = 1000.0_8 ; td_cs_cm2(37) =   0.041d-16

  do n = 1, N_tdp
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_He_elastic
!
!-------------------------------------
! above 1000 eV the cross section is calculated using
! an interplation formula derived upon assuming that in log-log scale the 
! cross section decays with energy as a straight line 
! continuing the tabulated data
! (see Figure 2 of [T.Majeed and D.Strickland, J.Phys.Chem.Ref.Data, v.26, p.335-349, 1997]
!
real(8) function CSV_He_elastic_m3s(E_eV)

  use He_elastic
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  real(8) CS_cm2

  if (E_eV.lt.td_E_eV(1)) then
     CSV_He_elastic_m3s = 1.0d-4 * td_cs_cm2(1) * sqrt(E_eV) * efactor_eV_to_ms
     return
  end if

  if (E_eV.GT.td_E_eV(N_tdp)) then
     CS_cm2 = (1.0d-15 * log(E_eV) + 2.4d-14) / E_eV**0.83
     CSV_He_elastic_m3s = max(0.0_8, 1.0d-4 * CS_cm2 * sqrt(E_eV) * efactor_eV_to_ms)
     return
  end if

  do n = 1, N_tdp-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_He_elastic_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_He_elastic_m3s

!---------------------------------------------------------------------------------------------------- #54
! He -> He+
!
! From Table V of D. Rapp and P. Englander-Golden, "Total Cross Sections for Ionization and Attachment 
! in Gases by Electron Impact. I. Positive Ionization", Journal of Chemical Physics, vol.43, pp.1464-1479, 1965.
!
module He_ionHe
  integer, parameter :: N_tdp=63  ! number of table data points
! table data
  real(8) td_E_eV(1:N_tdp)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_tdp)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_tdp)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module He_ionHe
!
!-------------------------------------
!
subroutine Prepare_He_ionHe

  use He_ionHe
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)  =   24.6_8 ; td_cs_cm2(1)  = 2.10d-19
  td_E_eV(2)  =   25.0_8 ; td_cs_cm2(2)  = 5.19d-19
  td_E_eV(3)  =   25.5_8 ; td_cs_cm2(3)  = 1.13d-18
  td_E_eV(4)  =   26.0_8 ; td_cs_cm2(4)  = 1.75d-18
  td_E_eV(5)  =   26.5_8 ; td_cs_cm2(5)  = 2.36d-18
  td_E_eV(6)  =   27.0_8 ; td_cs_cm2(6)  = 3.03d-18
  td_E_eV(7)  =   27.5_8 ; td_cs_cm2(7)  = 3.63d-18
  td_E_eV(8)  =   28.0_8 ; td_cs_cm2(8)  = 4.25d-18
  td_E_eV(9)  =   28.5_8 ; td_cs_cm2(9)  = 4.86d-18
  td_E_eV(10) =   29.0_8 ; td_cs_cm2(10) = 5.52d-18
  td_E_eV(11) =   29.5_8 ; td_cs_cm2(11) = 6.14d-18
  td_E_eV(12) =   30.0_8 ; td_cs_cm2(12) = 6.73d-18
  td_E_eV(13) =   30.5_8 ; td_cs_cm2(13) = 7.37d-18
  td_E_eV(14) =   31.0_8 ; td_cs_cm2(14) = 8.02d-18
  td_E_eV(15) =   31.5_8 ; td_cs_cm2(15) = 8.64d-18
  td_E_eV(16) =   32.0_8 ; td_cs_cm2(16) = 9.24d-18
  td_E_eV(17) =   32.5_8 ; td_cs_cm2(17) = 9.85d-18
  td_E_eV(18) =   33.0_8 ; td_cs_cm2(18) = 1.04d-17
  td_E_eV(19) =   33.5_8 ; td_cs_cm2(19) = 1.09d-17
  td_E_eV(20) =   34.0_8 ; td_cs_cm2(20) = 1.14d-17
  td_E_eV(21) =   36.0_8 ; td_cs_cm2(21) = 1.35d-17
  td_E_eV(22) =   38.0_8 ; td_cs_cm2(22) = 1.55d-17
  td_E_eV(23) =   40.0_8 ; td_cs_cm2(23) = 1.72d-17
  td_E_eV(24) =   45.0_8 ; td_cs_cm2(24) = 2.10d-17
  td_E_eV(25) =   50.0_8 ; td_cs_cm2(25) = 2.43d-17
  td_E_eV(26) =   55.0_8 ; td_cs_cm2(26) = 2.71d-17
  td_E_eV(27) =   60.0_8 ; td_cs_cm2(27) = 2.90d-17
  td_E_eV(28) =   65.0_8 ; td_cs_cm2(28) = 3.08d-17
  td_E_eV(29) =   70.0_8 ; td_cs_cm2(29) = 3.21d-17
  td_E_eV(30) =   75.0_8 ; td_cs_cm2(30) = 3.34d-17
  td_E_eV(31) =   80.0_8 ; td_cs_cm2(31) = 3.44d-17
  td_E_eV(32) =   85.0_8 ; td_cs_cm2(32) = 3.51d-17
  td_E_eV(33) =   90.0_8 ; td_cs_cm2(33) = 3.57d-17
  td_E_eV(34) =   95.0_8 ; td_cs_cm2(34) = 3.62d-17
  td_E_eV(35) =  100.0_8 ; td_cs_cm2(35) = 3.66d-17
  td_E_eV(36) =  105.0_8 ; td_cs_cm2(36) = 3.69d-17
  td_E_eV(37) =  110.0_8 ; td_cs_cm2(37) = 3.71d-17
  td_E_eV(38) =  115.0_8 ; td_cs_cm2(38) = 3.72d-17
  td_E_eV(39) =  120.0_8 ; td_cs_cm2(39) = 3.73d-17
  td_E_eV(40) =  125.0_8 ; td_cs_cm2(40) = 3.74d-17
  td_E_eV(41) =  130.0_8 ; td_cs_cm2(41) = 3.74d-17
  td_E_eV(42) =  135.0_8 ; td_cs_cm2(42) = 3.73d-17
  td_E_eV(43) =  140.0_8 ; td_cs_cm2(43) = 3.72d-17
  td_E_eV(44) =  145.0_8 ; td_cs_cm2(44) = 3.70d-17
  td_E_eV(45) =  150.0_8 ; td_cs_cm2(45) = 3.69d-17
  td_E_eV(46) =  175.0_8 ; td_cs_cm2(46) = 3.59d-17
  td_E_eV(47) =  200.0_8 ; td_cs_cm2(47) = 3.47d-17
  td_E_eV(48) =  250.0_8 ; td_cs_cm2(48) = 3.21d-17
  td_E_eV(49) =  300.0_8 ; td_cs_cm2(49) = 2.96d-17
  td_E_eV(50) =  350.0_8 ; td_cs_cm2(50) = 2.75d-17
  td_E_eV(51) =  400.0_8 ; td_cs_cm2(51) = 2.57d-17
  td_E_eV(52) =  450.0_8 ; td_cs_cm2(52) = 2.39d-17
  td_E_eV(53) =  500.0_8 ; td_cs_cm2(53) = 2.24d-17
  td_E_eV(54) =  550.0_8 ; td_cs_cm2(54) = 2.11d-17
  td_E_eV(55) =  600.0_8 ; td_cs_cm2(55) = 2.00d-17
  td_E_eV(56) =  650.0_8 ; td_cs_cm2(56) = 1.90d-17
  td_E_eV(57) =  700.0_8 ; td_cs_cm2(57) = 1.80d-17
  td_E_eV(58) =  750.0_8 ; td_cs_cm2(58) = 1.71d-17
  td_E_eV(59) =  800.0_8 ; td_cs_cm2(59) = 1.65d-17
  td_E_eV(60) =  850.0_8 ; td_cs_cm2(60) = 1.57d-17
  td_E_eV(61) =  900.0_8 ; td_cs_cm2(61) = 1.50d-17
  td_E_eV(62) =  950.0_8 ; td_cs_cm2(62) = 1.45d-17
  td_E_eV(63) = 1000.0_8 ; td_cs_cm2(63) = 1.41d-17

  do n = 1, N_tdp
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_He_ionHe
!
!-------------------------------------
!
real(8) function CSV_He_ionHe_m3s(E_eV)

  use He_ionHe
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  real(8) CS_cm2  ! cross section, used outside of the table data range

  IF (E_eV.LT.td_E_eV(1)) THEN
     CSV_He_ionHe_m3s = 0.0_8
     RETURN
  END IF

  IF (E_eV.GT.td_E_eV(N_tdp)) THEN
     CS_cm2 = (2.14d-14 * log(E_eV) - 8.02d-14) / E_eV
     CSV_He_ionHe_m3s = max(0.0_8, 1.0d-4 * CS_cm2 * sqrt(E_eV) * efactor_eV_to_ms)
     RETURN
  END IF

  DO n = 1, N_tdp-1
     IF ((E_eV.GE.td_E_eV(n)).AND.(E_eV.LE.td_E_eV(n+1))) THEN
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_He_ionHe_m3s = MAX(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        RETURN
     END IF
  END DO

end function CSV_He_ionHe_m3s

!---------------------------------------------------------------------------------------------------- #55
! He+ excitation 1^1S --> 2^1P dipole allowed
! 
! Eqs (1), (2), and Table 1 of Ralchenko et al., Atomic Data and Nuclear Data Tables, 94 (2008) 603-622
! threshold is from Biagi database, www.lxcat.net
!
module He_21P
  integer, parameter :: N_tdp=58   ! number of table data points
! table data
  real(8) td_E_eV(1:N_tdp)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_tdp)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_tdp)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module He_21P

!-------------------------------------
!
subroutine Prepare_He_21P

  use He_21P
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)   =    21.220_8 ; td_cs_cm2(1)   = 0.132d-18
  td_E_eV(2)   =    21.320_8 ; td_cs_cm2(2)   = 0.216d-18
  td_E_eV(3)   =    21.433_8 ; td_cs_cm2(3)   = 0.307d-18
  td_E_eV(4)   =    21.561_8 ; td_cs_cm2(4)   = 0.404d-18
  td_E_eV(5)   =    21.707_8 ; td_cs_cm2(5)   = 0.508d-18
  td_E_eV(6)   =    21.871_8 ; td_cs_cm2(6)   = 0.619d-18
  td_E_eV(7)   =    22.057_8 ; td_cs_cm2(7)   = 0.735d-18
  td_E_eV(8)   =    22.268_8 ; td_cs_cm2(8)   = 0.856d-18
  td_E_eV(9)   =    22.507_8 ; td_cs_cm2(9)   = 0.982d-18
  td_E_eV(10)  =    22.777_8 ; td_cs_cm2(10)  = 0.111d-17
  td_E_eV(11)  =    23.083_8 ; td_cs_cm2(11)  = 0.125d-17
  td_E_eV(12)  =    23.430_8 ; td_cs_cm2(12)  = 0.138d-17
  td_E_eV(13)  =    23.822_8 ; td_cs_cm2(13)  = 0.153d-17
  td_E_eV(14)  =    24.266_8 ; td_cs_cm2(14)  = 0.167d-17
  td_E_eV(15)  =    24.770_8 ; td_cs_cm2(15)  = 0.182d-17
  td_E_eV(16)  =    25.339_8 ; td_cs_cm2(16)  = 0.198d-17
  td_E_eV(17)  =    25.985_8 ; td_cs_cm2(17)  = 0.216d-17
  td_E_eV(18)  =    26.715_8 ; td_cs_cm2(18)  = 0.235d-17
  td_E_eV(19)  =    27.543_8 ; td_cs_cm2(19)  = 0.257d-17
  td_E_eV(20)  =    28.479_8 ; td_cs_cm2(20)  = 0.282d-17
  td_E_eV(21)  =    29.540_8 ; td_cs_cm2(21)  = 0.311d-17
  td_E_eV(22)  =    30.741_8 ; td_cs_cm2(22)  = 0.345d-17
  td_E_eV(23)  =    32.101_8 ; td_cs_cm2(23)  = 0.385d-17
  td_E_eV(24)  =    33.641_8 ; td_cs_cm2(24)  = 0.431d-17
  td_E_eV(25)  =    35.385_8 ; td_cs_cm2(25)  = 0.482d-17
  td_E_eV(26)  =    37.360_8 ; td_cs_cm2(26)  = 0.539d-17
  td_E_eV(27)  =    39.596_8 ; td_cs_cm2(27)  = 0.601d-17
  td_E_eV(28)  =    42.128_8 ; td_cs_cm2(28)  = 0.665d-17
  td_E_eV(29)  =    44.995_8 ; td_cs_cm2(29)  = 0.732d-17
  td_E_eV(30)  =    48.242_8 ; td_cs_cm2(30)  = 0.797d-17
  td_E_eV(31)  =    51.918_8 ; td_cs_cm2(31)  = 0.859d-17
  td_E_eV(32)  =    56.081_8 ; td_cs_cm2(32)  = 0.917d-17
  td_E_eV(33)  =    60.794_8 ; td_cs_cm2(33)  = 0.967d-17
  td_E_eV(34)  =    66.132_8 ; td_cs_cm2(34)  = 0.101d-16
  td_E_eV(35)  =    72.176_8 ; td_cs_cm2(35)  = 0.104d-16
  td_E_eV(36)  =    79.020_8 ; td_cs_cm2(36)  = 0.106d-16
  td_E_eV(37)  =    86.769_8 ; td_cs_cm2(37)  = 0.108d-16
  td_E_eV(38)  =    95.544_8 ; td_cs_cm2(38)  = 0.108d-16
  td_E_eV(39)  =   105.481_8 ; td_cs_cm2(39)  = 0.107d-16
  td_E_eV(40)  =   116.733_8 ; td_cs_cm2(40)  = 0.105d-16
  td_E_eV(41)  =   129.474_8 ; td_cs_cm2(41)  = 0.103d-16
  td_E_eV(42)  =   143.900_8 ; td_cs_cm2(42)  = 0.993d-17
  td_E_eV(43)  =   160.237_8 ; td_cs_cm2(43)  = 0.957d-17
  td_E_eV(44)  =   178.735_8 ; td_cs_cm2(44)  = 0.916d-17
  td_E_eV(45)  =   199.682_8 ; td_cs_cm2(45)  = 0.872d-17
  td_E_eV(46)  =   223.401_8 ; td_cs_cm2(46)  = 0.826d-17
  td_E_eV(47)  =   250.259_8 ; td_cs_cm2(47)  = 0.779d-17
  td_E_eV(48)  =   280.671_8 ; td_cs_cm2(48)  = 0.732d-17
  td_E_eV(49)  =   315.109_8 ; td_cs_cm2(49)  = 0.685d-17
  td_E_eV(50)  =   354.104_8 ; td_cs_cm2(50)  = 0.639d-17
  td_E_eV(51)  =   398.260_8 ; td_cs_cm2(51)  = 0.594d-17
  td_E_eV(52)  =   448.260_8 ; td_cs_cm2(52)  = 0.551d-17
  td_E_eV(53)  =   504.877_8 ; td_cs_cm2(53)  = 0.510d-17
  td_E_eV(54)  =   568.988_8 ; td_cs_cm2(54)  = 0.471d-17
  td_E_eV(55)  =   641.583_8 ; td_cs_cm2(55)  = 0.434d-17
  td_E_eV(56)  =   723.786_8 ; td_cs_cm2(56)  = 0.399d-17
  td_E_eV(57)  =   816.868_8 ; td_cs_cm2(57)  = 0.366d-17
  td_E_eV(58)  =   922.269_8 ; td_cs_cm2(58)  = 0.335d-17

  do n = 1, N_tdp
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_He_21P
!
!-------------------------------------
!
real(8) function CSV_He_21P_m3s(E_eV)

  use He_21P
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none
  
  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_tdp))) then
     CSV_He_21P_m3s = 0.0_8
     return
  end if

  do n = 1, N_tdp-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_He_21P_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_He_21P_m3s

!---------------------------------------------------------------------------------------------------- #56
! He+ excitation 1^1S --> 3^1P dipole allowed
!
! Eqs (1), (2), and Table 1 of Ralchenko et al., Atomic Data and Nuclear Data Tables, 94 (2008) 603-622
! threshold is from Biagi database, www.lxcat.net
!
module He_31P
  integer, parameter :: N_tdp=58   ! number of table data points
! table data
  real(8) td_E_eV(1:N_tdp)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_tdp)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_tdp)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module He_31P

!-------------------------------------
!
subroutine Prepare_He_31P

  use He_31P
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)   =   23.087_8 ; td_cs_cm2(1)   = 0.335d-19
  td_E_eV(2)   =   23.187_8 ; td_cs_cm2(2)   = 0.502d-19
  td_E_eV(3)   =   23.300_8 ; td_cs_cm2(3)   = 0.685d-19
  td_E_eV(4)   =   23.428_8 ; td_cs_cm2(4)   = 0.885d-19
  td_E_eV(5)   =   23.574_8 ; td_cs_cm2(5)   = 0.110d-18
  td_E_eV(6)   =   23.738_8 ; td_cs_cm2(6)   = 0.134d-18
  td_E_eV(7)   =   23.924_8 ; td_cs_cm2(7)   = 0.159d-18
  td_E_eV(8)   =   24.135_8 ; td_cs_cm2(8)   = 0.186d-18
  td_E_eV(9)   =   24.374_8 ; td_cs_cm2(9)   = 0.216d-18
  td_E_eV(10)  =   24.644_8 ; td_cs_cm2(10)  = 0.247d-18
  td_E_eV(11)  =   24.950_8 ; td_cs_cm2(11)  = 0.280d-18
  td_E_eV(12)  =   25.297_8 ; td_cs_cm2(12)  = 0.315d-18
  td_E_eV(13)  =   25.689_8 ; td_cs_cm2(13)  = 0.352d-18
  td_E_eV(14)  =   26.133_8 ; td_cs_cm2(14)  = 0.391d-18
  td_E_eV(15)  =   26.637_8 ; td_cs_cm2(15)  = 0.433d-18
  td_E_eV(16)  =   27.206_8 ; td_cs_cm2(16)  = 0.478d-18
  td_E_eV(17)  =   27.852_8 ; td_cs_cm2(17)  = 0.526d-18
  td_E_eV(18)  =   28.582_8 ; td_cs_cm2(18)  = 0.578d-18
  td_E_eV(19)  =   29.410_8 ; td_cs_cm2(19)  = 0.635d-18
  td_E_eV(20)  =   30.346_8 ; td_cs_cm2(20)  = 0.699d-18
  td_E_eV(21)  =   31.407_8 ; td_cs_cm2(21)  = 0.770d-18
  td_E_eV(22)  =   32.608_8 ; td_cs_cm2(22)  = 0.850d-18
  td_E_eV(23)  =   33.968_8 ; td_cs_cm2(23)  = 0.939d-18
  td_E_eV(24)  =   35.508_8 ; td_cs_cm2(24)  = 0.104d-17
  td_E_eV(25)  =   37.252_8 ; td_cs_cm2(25)  = 0.115d-17
  td_E_eV(26)  =   39.227_8 ; td_cs_cm2(26)  = 0.127d-17
  td_E_eV(27)  =   41.463_8 ; td_cs_cm2(27)  = 0.140d-17
  td_E_eV(28)  =   43.995_8 ; td_cs_cm2(28)  = 0.154d-17
  td_E_eV(29)  =   46.862_8 ; td_cs_cm2(29)  = 0.168d-17
  td_E_eV(30)  =   50.109_8 ; td_cs_cm2(30)  = 0.183d-17
  td_E_eV(31)  =   53.785_8 ; td_cs_cm2(31)  = 0.197d-17
  td_E_eV(32)  =   57.948_8 ; td_cs_cm2(32)  = 0.210d-17
  td_E_eV(33)  =   62.661_8 ; td_cs_cm2(33)  = 0.222d-17
  td_E_eV(34)  =   67.999_8 ; td_cs_cm2(34)  = 0.232d-17
  td_E_eV(35)  =   74.043_8 ; td_cs_cm2(35)  = 0.241d-17
  td_E_eV(36)  =   80.887_8 ; td_cs_cm2(36)  = 0.247d-17
  td_E_eV(37)  =   88.636_8 ; td_cs_cm2(37)  = 0.251d-17
  td_E_eV(38)  =   97.411_8 ; td_cs_cm2(38)  = 0.253d-17
  td_E_eV(39)  =  107.348_8 ; td_cs_cm2(39)  = 0.253d-17
  td_E_eV(40)  =  118.600_8 ; td_cs_cm2(40)  = 0.250d-17
  td_E_eV(41)  =  131.341_8 ; td_cs_cm2(41)  = 0.246d-17
  td_E_eV(42)  =  145.768_8 ; td_cs_cm2(42)  = 0.240d-17
  td_E_eV(43)  =  162.104_8 ; td_cs_cm2(43)  = 0.232d-17
  td_E_eV(44)  =  180.602_8 ; td_cs_cm2(44)  = 0.223d-17
  td_E_eV(45)  =  201.549_8 ; td_cs_cm2(45)  = 0.214d-17
  td_E_eV(46)  =  225.268_8 ; td_cs_cm2(46)  = 0.203d-17
  td_E_eV(47)  =  252.126_8 ; td_cs_cm2(47)  = 0.193d-17
  td_E_eV(48)  =  282.538_8 ; td_cs_cm2(48)  = 0.182d-17
  td_E_eV(49)  =  316.976_8 ; td_cs_cm2(49)  = 0.170d-17
  td_E_eV(50)  =  355.971_8 ; td_cs_cm2(50)  = 0.159d-17
  td_E_eV(51)  =  400.127_8 ; td_cs_cm2(51)  = 0.149d-17
  td_E_eV(52)  =  450.127_8 ; td_cs_cm2(52)  = 0.138d-17
  td_E_eV(53)  =  506.744_8 ; td_cs_cm2(53)  = 0.128d-17
  td_E_eV(54)  =  570.855_8 ; td_cs_cm2(54)  = 0.118d-17
  td_E_eV(55)  =  643.450_8 ; td_cs_cm2(55)  = 0.109d-17
  td_E_eV(56)  =  725.653_8 ; td_cs_cm2(56)  = 0.100d-17
  td_E_eV(57)  =  818.735_8 ; td_cs_cm2(57)  = 0.923d-18
  td_E_eV(58)  =  924.136_8 ; td_cs_cm2(58)  = 0.846d-18

  do n = 1, N_tdp
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_He_31P
!
!-------------------------------------
!
real(8) function CSV_He_31P_m3s(E_eV)

  use He_31P
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none

  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_tdp))) then
     CSV_He_31P_m3s = 0.0_8
     return
  end if

  do n = 1, N_tdp-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_He_31P_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_He_31P_m3s

!---------------------------------------------------------------------------------------------------- #57
! He+ excitation 1^1S --> 4^1P dipole allowed
!
! Eqs (1), (2), and Table 1 of Ralchenko et al., Atomic Data and Nuclear Data Tables, 94 (2008) 603-622
! threshold is from Biagi database, www.lxcat.net
!
module He_41P
  integer, parameter :: N_tdp=58   ! number of table data points
! table data
  real(8) td_E_eV(1:N_tdp)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_tdp)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_tdp)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module He_41P

!-------------------------------------
!
subroutine Prepare_He_41P

  use He_41P
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)   =   23.724_8 ; td_cs_cm2(1)   = 0.900d-20
  td_E_eV(2)   =   23.824_8 ; td_cs_cm2(2)   = 0.188d-19
  td_E_eV(3)   =   23.937_8 ; td_cs_cm2(3)   = 0.294d-19
  td_E_eV(4)   =   24.065_8 ; td_cs_cm2(4)   = 0.408d-19
  td_E_eV(5)   =   24.211_8 ; td_cs_cm2(5)   = 0.530d-19
  td_E_eV(6)   =   24.375_8 ; td_cs_cm2(6)   = 0.659d-19
  td_E_eV(7)   =   24.561_8 ; td_cs_cm2(7)   = 0.796d-19
  td_E_eV(8)   =   24.772_8 ; td_cs_cm2(8)   = 0.940d-19
  td_E_eV(9)   =   25.011_8 ; td_cs_cm2(9)   = 0.109d-18
  td_E_eV(10)  =   25.281_8 ; td_cs_cm2(10)  = 0.124d-18
  td_E_eV(11)  =   25.587_8 ; td_cs_cm2(11)  = 0.140d-18
  td_E_eV(12)  =   25.934_8 ; td_cs_cm2(12)  = 0.156d-18
  td_E_eV(13)  =   26.326_8 ; td_cs_cm2(13)  = 0.173d-18
  td_E_eV(14)  =   26.770_8 ; td_cs_cm2(14)  = 0.189d-18
  td_E_eV(15)  =   27.274_8 ; td_cs_cm2(15)  = 0.206d-18
  td_E_eV(16)  =   27.843_8 ; td_cs_cm2(16)  = 0.224d-18
  td_E_eV(17)  =   28.489_8 ; td_cs_cm2(17)  = 0.242d-18
  td_E_eV(18)  =   29.219_8 ; td_cs_cm2(18)  = 0.262d-18
  td_E_eV(19)  =   30.047_8 ; td_cs_cm2(19)  = 0.283d-18
  td_E_eV(20)  =   30.983_8 ; td_cs_cm2(20)  = 0.307d-18
  td_E_eV(21)  =   32.044_8 ; td_cs_cm2(21)  = 0.335d-18
  td_E_eV(22)  =   33.245_8 ; td_cs_cm2(22)  = 0.367d-18
  td_E_eV(23)  =   34.605_8 ; td_cs_cm2(23)  = 0.404d-18
  td_E_eV(24)  =   36.145_8 ; td_cs_cm2(24)  = 0.447d-18
  td_E_eV(25)  =   37.889_8 ; td_cs_cm2(25)  = 0.496d-18
  td_E_eV(26)  =   39.864_8 ; td_cs_cm2(26)  = 0.550d-18
  td_E_eV(27)  =   42.100_8 ; td_cs_cm2(27)  = 0.610d-18
  td_E_eV(28)  =   44.632_8 ; td_cs_cm2(28)  = 0.673d-18
  td_E_eV(29)  =   47.499_8 ; td_cs_cm2(29)  = 0.739d-18
  td_E_eV(30)  =   50.746_8 ; td_cs_cm2(30)  = 0.805d-18
  td_E_eV(31)  =   54.422_8 ; td_cs_cm2(31)  = 0.868d-18
  td_E_eV(32)  =   58.585_8 ; td_cs_cm2(32)  = 0.927d-18
  td_E_eV(33)  =   63.298_8 ; td_cs_cm2(33)  = 0.980d-18
  td_E_eV(34)  =   68.636_8 ; td_cs_cm2(34)  = 0.102d-17
  td_E_eV(35)  =   74.680_8 ; td_cs_cm2(35)  = 0.106d-17
  td_E_eV(36)  =   81.524_8 ; td_cs_cm2(36)  = 0.108d-17
  td_E_eV(37)  =   89.273_8 ; td_cs_cm2(37)  = 0.110d-17
  td_E_eV(38)  =   98.048_8 ; td_cs_cm2(38)  = 0.110d-17
  td_E_eV(39)  =  107.985_8 ; td_cs_cm2(39)  = 0.109d-17
  td_E_eV(40)  =  119.237_8 ; td_cs_cm2(40)  = 0.107d-17
  td_E_eV(41)  =  131.978_8 ; td_cs_cm2(41)  = 0.105d-17
  td_E_eV(42)  =  146.404_8 ; td_cs_cm2(42)  = 0.102d-17
  td_E_eV(43)  =  162.741_8 ; td_cs_cm2(43)  = 0.979d-18
  td_E_eV(44)  =  181.239_8 ; td_cs_cm2(44)  = 0.937d-18
  td_E_eV(45)  =  202.186_8 ; td_cs_cm2(45)  = 0.891d-18
  td_E_eV(46)  =  225.905_8 ; td_cs_cm2(46)  = 0.844d-18
  td_E_eV(47)  =  252.763_8 ; td_cs_cm2(47)  = 0.795d-18
  td_E_eV(48)  =  283.175_8 ; td_cs_cm2(48)  = 0.746d-18
  td_E_eV(49)  =  317.613_8 ; td_cs_cm2(49)  = 0.698d-18
  td_E_eV(50)  =  356.608_8 ; td_cs_cm2(50)  = 0.650d-18
  td_E_eV(51)  =  400.764_8 ; td_cs_cm2(51)  = 0.604d-18
  td_E_eV(52)  =  450.764_8 ; td_cs_cm2(52)  = 0.560d-18
  td_E_eV(53)  =  507.381_8 ; td_cs_cm2(53)  = 0.517d-18
  td_E_eV(54)  =  571.492_8 ; td_cs_cm2(54)  = 0.477d-18
  td_E_eV(55)  =  644.087_8 ; td_cs_cm2(55)  = 0.439d-18
  td_E_eV(56)  =  726.290_8 ; td_cs_cm2(56)  = 0.403d-18
  td_E_eV(57)  =  819.372_8 ; td_cs_cm2(57)  = 0.370d-18
  td_E_eV(58)  =  924.773_8 ; td_cs_cm2(58)  = 0.339d-18

  do n = 1, N_tdp
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_He_41P
!
!-------------------------------------
!
real(8) function CSV_He_41P_m3s(E_eV)

  use He_41P
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none

  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_tdp))) then
     CSV_He_41P_m3s = 0.0_8
     return
  end if

  do n = 1, N_tdp-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_He_41P_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_He_41P_m3s


!---------------------------------------------------------------------------------------------------- #58
! He+ excitation 1^1S --> 2^1S dipole forbidden
!
! Eqs (1), (3), and Table 2 of Ralchenko et al., Atomic Data and Nuclear Data Tables, 94 (2008) 603-622
! threshold is from Biagi database, www.lxcat.net
!
module He_21S
  integer, parameter :: N_tdp=58   ! number of table data points
! table data
  real(8) td_E_eV(1:N_tdp)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_tdp)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_tdp)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module He_21S

!-------------------------------------
!
subroutine Prepare_He_21S

  use He_21S
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)   =   20.620_8 ; td_cs_cm2(1)   = 0.212E-17
  td_E_eV(2)   =   20.720_8 ; td_cs_cm2(2)   = 0.213E-17
  td_E_eV(3)   =   20.833_8 ; td_cs_cm2(3)   = 0.214E-17
  td_E_eV(4)   =   20.961_8 ; td_cs_cm2(4)   = 0.215E-17
  td_E_eV(5)   =   21.107_8 ; td_cs_cm2(5)   = 0.216E-17
  td_E_eV(6)   =   21.271_8 ; td_cs_cm2(6)   = 0.217E-17
  td_E_eV(7)   =   21.457_8 ; td_cs_cm2(7)   = 0.218E-17
  td_E_eV(8)   =   21.668_8 ; td_cs_cm2(8)   = 0.220E-17
  td_E_eV(9)   =   21.907_8 ; td_cs_cm2(9)   = 0.221E-17
  td_E_eV(10)  =   22.177_8 ; td_cs_cm2(10)  = 0.223E-17
  td_E_eV(11)  =   22.483_8 ; td_cs_cm2(11)  = 0.224E-17
  td_E_eV(12)  =   22.830_8 ; td_cs_cm2(12)  = 0.226E-17
  td_E_eV(13)  =   23.222_8 ; td_cs_cm2(13)  = 0.227E-17
  td_E_eV(14)  =   23.666_8 ; td_cs_cm2(14)  = 0.228E-17
  td_E_eV(15)  =   24.170_8 ; td_cs_cm2(15)  = 0.229E-17
  td_E_eV(16)  =   24.739_8 ; td_cs_cm2(16)  = 0.230E-17
  td_E_eV(17)  =   25.385_8 ; td_cs_cm2(17)  = 0.231E-17
  td_E_eV(18)  =   26.115_8 ; td_cs_cm2(18)  = 0.231E-17
  td_E_eV(19)  =   26.943_8 ; td_cs_cm2(19)  = 0.231E-17
  td_E_eV(20)  =   27.879_8 ; td_cs_cm2(20)  = 0.230E-17
  td_E_eV(21)  =   28.940_8 ; td_cs_cm2(21)  = 0.229E-17
  td_E_eV(22)  =   30.141_8 ; td_cs_cm2(22)  = 0.226E-17
  td_E_eV(23)  =   31.501_8 ; td_cs_cm2(23)  = 0.224E-17
  td_E_eV(24)  =   33.041_8 ; td_cs_cm2(24)  = 0.220E-17
  td_E_eV(25)  =   34.785_8 ; td_cs_cm2(25)  = 0.215E-17
  td_E_eV(26)  =   36.760_8 ; td_cs_cm2(26)  = 0.210E-17
  td_E_eV(27)  =   38.996_8 ; td_cs_cm2(27)  = 0.203E-17
  td_E_eV(28)  =   41.528_8 ; td_cs_cm2(28)  = 0.196E-17
  td_E_eV(29)  =   44.395_8 ; td_cs_cm2(29)  = 0.189E-17
  td_E_eV(30)  =   47.642_8 ; td_cs_cm2(30)  = 0.181E-17
  td_E_eV(31)  =   51.318_8 ; td_cs_cm2(31)  = 0.172E-17
  td_E_eV(32)  =   55.481_8 ; td_cs_cm2(32)  = 0.164E-17
  td_E_eV(33)  =   60.194_8 ; td_cs_cm2(33)  = 0.155E-17
  td_E_eV(34)  =   65.532_8 ; td_cs_cm2(34)  = 0.147E-17
  td_E_eV(35)  =   71.576_8 ; td_cs_cm2(35)  = 0.139E-17
  td_E_eV(36)  =   78.420_8 ; td_cs_cm2(36)  = 0.131E-17
  td_E_eV(37)  =   86.169_8 ; td_cs_cm2(37)  = 0.124E-17
  td_E_eV(38)  =   94.944_8 ; td_cs_cm2(38)  = 0.117E-17
  td_E_eV(39)  =  104.881_8 ; td_cs_cm2(39)  = 0.111E-17
  td_E_eV(40)  =  116.133_8 ; td_cs_cm2(40)  = 0.105E-17
  td_E_eV(41)  =  128.874_8 ; td_cs_cm2(41)  = 0.991E-18
  td_E_eV(42)  =  143.301_8 ; td_cs_cm2(42)  = 0.935E-18
  td_E_eV(43)  =  159.637_8 ; td_cs_cm2(43)  = 0.881E-18
  td_E_eV(44)  =  178.135_8 ; td_cs_cm2(44)  = 0.827E-18
  td_E_eV(45)  =  199.082_8 ; td_cs_cm2(45)  = 0.774E-18
  td_E_eV(46)  =  222.801_8 ; td_cs_cm2(46)  = 0.721E-18
  td_E_eV(47)  =  249.659_8 ; td_cs_cm2(47)  = 0.668E-18
  td_E_eV(48)  =  280.071_8 ; td_cs_cm2(48)  = 0.617E-18
  td_E_eV(49)  =  314.509_8 ; td_cs_cm2(49)  = 0.567E-18
  td_E_eV(50)  =  353.504_8 ; td_cs_cm2(50)  = 0.519E-18
  td_E_eV(51)  =  397.660_8 ; td_cs_cm2(51)  = 0.473E-18
  td_E_eV(52)  =  447.660_8 ; td_cs_cm2(52)  = 0.430E-18
  td_E_eV(53)  =  504.277_8 ; td_cs_cm2(53)  = 0.389E-18
  td_E_eV(54)  =  568.388_8 ; td_cs_cm2(54)  = 0.351E-18
  td_E_eV(55)  =  640.983_8 ; td_cs_cm2(55)  = 0.316E-18
  td_E_eV(56)  =  723.186_8 ; td_cs_cm2(56)  = 0.284E-18
  td_E_eV(57)  =  816.268_8 ; td_cs_cm2(57)  = 0.255E-18
  td_E_eV(58)  =  921.669_8 ; td_cs_cm2(58)  = 0.228E-18

  do n = 1, N_tdp
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_He_21S
!
!-------------------------------------
!
real(8) function CSV_He_21S_m3s(E_eV)

  use He_21S
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none

  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_tdp))) then
     CSV_He_21S_m3s = 0.0_8
     return
  end if

  do n = 1, N_tdp-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_He_21S_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_He_21S_m3s

!---------------------------------------------------------------------------------------------------- #59
! He+ excitation 1^1S --> 2^3S spin forbidden
!
! Eqs (1), (4), and Table 3 of Ralchenko et al., Atomic Data and Nuclear Data Tables, 94 (2008) 603-622
! threshold is from Biagi database, www.lxcat.net
!
module He_23S
  integer, parameter :: N_tdp=58   ! number of table data points
! table data
  real(8) td_E_eV(1:N_tdp)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_tdp)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_tdp)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module He_23S

!-------------------------------------
!
subroutine Prepare_He_23S

  use He_23S
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)   =   19.820_8 ; td_cs_cm2(1)   = 0.388d-17
  td_E_eV(2)   =   19.920_8 ; td_cs_cm2(2)   = 0.386d-17
  td_E_eV(3)   =   20.033_8 ; td_cs_cm2(3)   = 0.383d-17
  td_E_eV(4)   =   20.161_8 ; td_cs_cm2(4)   = 0.381d-17
  td_E_eV(5)   =   20.307_8 ; td_cs_cm2(5)   = 0.378d-17
  td_E_eV(6)   =   20.471_8 ; td_cs_cm2(6)   = 0.374d-17
  td_E_eV(7)   =   20.657_8 ; td_cs_cm2(7)   = 0.370d-17
  td_E_eV(8)   =   20.868_8 ; td_cs_cm2(8)   = 0.365d-17
  td_E_eV(9)   =   21.107_8 ; td_cs_cm2(9)   = 0.360d-17
  td_E_eV(10)  =   21.377_8 ; td_cs_cm2(10)  = 0.353d-17
  td_E_eV(11)  =   21.683_8 ; td_cs_cm2(11)  = 0.347d-17
  td_E_eV(12)  =   22.030_8 ; td_cs_cm2(12)  = 0.339d-17
  td_E_eV(13)  =   22.422_8 ; td_cs_cm2(13)  = 0.330d-17
  td_E_eV(14)  =   22.866_8 ; td_cs_cm2(14)  = 0.320d-17
  td_E_eV(15)  =   23.370_8 ; td_cs_cm2(15)  = 0.309d-17
  td_E_eV(16)  =   23.939_8 ; td_cs_cm2(16)  = 0.297d-17
  td_E_eV(17)  =   24.585_8 ; td_cs_cm2(17)  = 0.284d-17
  td_E_eV(18)  =   25.315_8 ; td_cs_cm2(18)  = 0.270d-17
  td_E_eV(19)  =   26.143_8 ; td_cs_cm2(19)  = 0.255d-17
  td_E_eV(20)  =   27.079_8 ; td_cs_cm2(20)  = 0.239d-17
  td_E_eV(21)  =   28.140_8 ; td_cs_cm2(21)  = 0.223d-17
  td_E_eV(22)  =   29.341_8 ; td_cs_cm2(22)  = 0.205d-17
  td_E_eV(23)  =   30.701_8 ; td_cs_cm2(23)  = 0.188d-17
  td_E_eV(24)  =   32.241_8 ; td_cs_cm2(24)  = 0.170d-17
  td_E_eV(25)  =   33.985_8 ; td_cs_cm2(25)  = 0.153d-17
  td_E_eV(26)  =   35.960_8 ; td_cs_cm2(26)  = 0.136d-17
  td_E_eV(27)  =   38.196_8 ; td_cs_cm2(27)  = 0.120d-17
  td_E_eV(28)  =   40.728_8 ; td_cs_cm2(28)  = 0.105d-17
  td_E_eV(29)  =   43.595_8 ; td_cs_cm2(29)  = 0.910d-18
  td_E_eV(30)  =   46.842_8 ; td_cs_cm2(30)  = 0.781d-18
  td_E_eV(31)  =   50.518_8 ; td_cs_cm2(31)  = 0.665d-18
  td_E_eV(32)  =   54.681_8 ; td_cs_cm2(32)  = 0.562d-18
  td_E_eV(33)  =   59.394_8 ; td_cs_cm2(33)  = 0.471d-18
  td_E_eV(34)  =   64.732_8 ; td_cs_cm2(34)  = 0.392d-18
  td_E_eV(35)  =   70.776_8 ; td_cs_cm2(35)  = 0.325d-18
  td_E_eV(36)  =   77.620_8 ; td_cs_cm2(36)  = 0.267d-18
  td_E_eV(37)  =   85.369_8 ; td_cs_cm2(37)  = 0.217d-18
  td_E_eV(38)  =   94.144_8 ; td_cs_cm2(38)  = 0.176d-18
  td_E_eV(39)  =  104.081_8 ; td_cs_cm2(39)  = 0.142d-18
  td_E_eV(40)  =  115.333_8 ; td_cs_cm2(40)  = 0.113d-18
  td_E_eV(41)  =  128.074_8 ; td_cs_cm2(41)  = 0.891d-19
  td_E_eV(42)  =  142.501_8 ; td_cs_cm2(42)  = 0.696d-19
  td_E_eV(43)  =  158.837_8 ; td_cs_cm2(43)  = 0.539d-19
  td_E_eV(44)  =  177.335_8 ; td_cs_cm2(44)  = 0.412d-19
  td_E_eV(45)  =  198.282_8 ; td_cs_cm2(45)  = 0.312d-19
  td_E_eV(46)  =  222.001_8 ; td_cs_cm2(46)  = 0.234d-19
  td_E_eV(47)  =  248.859_8 ; td_cs_cm2(47)  = 0.173d-19
  td_E_eV(48)  =  279.271_8 ; td_cs_cm2(48)  = 0.127d-19
  td_E_eV(49)  =  313.709_8 ; td_cs_cm2(49)  = 0.925d-20
  td_E_eV(50)  =  352.704_8 ; td_cs_cm2(50)  = 0.668d-20
  td_E_eV(51)  =  396.860_8 ; td_cs_cm2(51)  = 0.479d-20
  td_E_eV(52)  =  446.860_8 ; td_cs_cm2(52)  = 0.341d-20
  td_E_eV(53)  =  503.477_8 ; td_cs_cm2(53)  = 0.242d-20
  td_E_eV(54)  =  567.588_8 ; td_cs_cm2(54)  = 0.171d-20
  td_E_eV(55)  =  640.183_8 ; td_cs_cm2(55)  = 0.120d-20
  td_E_eV(56)  =  722.386_8 ; td_cs_cm2(56)  = 0.841d-21
  td_E_eV(57)  =  815.468_8 ; td_cs_cm2(57)  = 0.588d-21
  td_E_eV(58)  =  920.869_8 ; td_cs_cm2(58)  = 0.410d-21

  do n = 1, N_tdp
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_He_23S
!
!-------------------------------------
!
real(8) function CSV_He_23S_m3s(E_eV)

  use He_23S
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none

  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_tdp))) then
     CSV_He_23S_m3s = 0.0_8
     return
  end if

  do n = 1, N_tdp-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_He_23S_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_He_23S_m3s

!---------------------------------------------------------------------------------------------------- #60
! He+ excitation 1^1S --> 2^3P spin forbidden
!
! Eqs (1), (4), and Table 3 of Ralchenko et al., Atomic Data and Nuclear Data Tables, 94 (2008) 603-622
! threshold is from Biagi database, www.lxcat.net
!
module He_23P
  integer, parameter :: N_tdp=58   ! number of table data points
! table data
  real(8) td_E_eV(1:N_tdp)        ! electron energy [eV]
  real(8) td_cs_cm2(1:N_tdp)      ! cross section [cm^2]
  real(8) td_csv_m3s(1:N_tdp)     ! cross section [cm^2] times velocity [m/s] times 1.0e-4 [m^2/cm^2]
end module He_23P

!-------------------------------------
!
subroutine Prepare_He_23P

  use He_23P
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none

  integer n

  td_E_eV(1)   =   20.960_8 ; td_cs_cm2(1)   = 0.513d-18
  td_E_eV(2)   =   21.060_8 ; td_cs_cm2(2)   = 0.600d-18
  td_E_eV(3)   =   21.173_8 ; td_cs_cm2(3)   = 0.695d-18
  td_E_eV(4)   =   21.301_8 ; td_cs_cm2(4)   = 0.797d-18
  td_E_eV(5)   =   21.447_8 ; td_cs_cm2(5)   = 0.907d-18
  td_E_eV(6)   =   21.611_8 ; td_cs_cm2(6)   = 0.102d-17
  td_E_eV(7)   =   21.797_8 ; td_cs_cm2(7)   = 0.115d-17
  td_E_eV(8)   =   22.008_8 ; td_cs_cm2(8)   = 0.128d-17
  td_E_eV(9)   =   22.247_8 ; td_cs_cm2(9)   = 0.141d-17
  td_E_eV(10)  =   22.517_8 ; td_cs_cm2(10)  = 0.155d-17
  td_E_eV(11)  =   22.823_8 ; td_cs_cm2(11)  = 0.169d-17
  td_E_eV(12)  =   23.170_8 ; td_cs_cm2(12)  = 0.183d-17
  td_E_eV(13)  =   23.562_8 ; td_cs_cm2(13)  = 0.197d-17
  td_E_eV(14)  =   24.006_8 ; td_cs_cm2(14)  = 0.209d-17
  td_E_eV(15)  =   24.510_8 ; td_cs_cm2(15)  = 0.221d-17
  td_E_eV(16)  =   25.079_8 ; td_cs_cm2(16)  = 0.231d-17
  td_E_eV(17)  =   25.725_8 ; td_cs_cm2(17)  = 0.239d-17
  td_E_eV(18)  =   26.455_8 ; td_cs_cm2(18)  = 0.245d-17
  td_E_eV(19)  =   27.283_8 ; td_cs_cm2(19)  = 0.248d-17
  td_E_eV(20)  =   28.219_8 ; td_cs_cm2(20)  = 0.248d-17
  td_E_eV(21)  =   29.280_8 ; td_cs_cm2(21)  = 0.246d-17
  td_E_eV(22)  =   30.481_8 ; td_cs_cm2(22)  = 0.240d-17
  td_E_eV(23)  =   31.841_8 ; td_cs_cm2(23)  = 0.231d-17
  td_E_eV(24)  =   33.381_8 ; td_cs_cm2(24)  = 0.220d-17
  td_E_eV(25)  =   35.125_8 ; td_cs_cm2(25)  = 0.206d-17
  td_E_eV(26)  =   37.100_8 ; td_cs_cm2(26)  = 0.191d-17
  td_E_eV(27)  =   39.336_8 ; td_cs_cm2(27)  = 0.174d-17
  td_E_eV(28)  =   41.868_8 ; td_cs_cm2(28)  = 0.156d-17
  td_E_eV(29)  =   44.735_8 ; td_cs_cm2(29)  = 0.138d-17
  td_E_eV(30)  =   47.982_8 ; td_cs_cm2(30)  = 0.121d-17
  td_E_eV(31)  =   51.658_8 ; td_cs_cm2(31)  = 0.104d-17
  td_E_eV(32)  =   55.821_8 ; td_cs_cm2(32)  = 0.878d-18
  td_E_eV(33)  =   60.534_8 ; td_cs_cm2(33)  = 0.732d-18
  td_E_eV(34)  =   65.872_8 ; td_cs_cm2(34)  = 0.602d-18
  td_E_eV(35)  =   71.916_8 ; td_cs_cm2(35)  = 0.488d-18
  td_E_eV(36)  =   78.760_8 ; td_cs_cm2(36)  = 0.390d-18
  td_E_eV(37)  =   86.509_8 ; td_cs_cm2(37)  = 0.307d-18
  td_E_eV(38)  =   95.284_8 ; td_cs_cm2(38)  = 0.238d-18
  td_E_eV(39)  =  105.221_8 ; td_cs_cm2(39)  = 0.182d-18
  td_E_eV(40)  =  116.473_8 ; td_cs_cm2(40)  = 0.138d-18
  td_E_eV(41)  =  129.214_8 ; td_cs_cm2(41)  = 0.102d-18
  td_E_eV(42)  =  143.641_8 ; td_cs_cm2(42)  = 0.753d-19
  td_E_eV(43)  =  159.977_8 ; td_cs_cm2(43)  = 0.547d-19
  td_E_eV(44)  =  178.475_8 ; td_cs_cm2(44)  = 0.392d-19
  td_E_eV(45)  =  199.422_8 ; td_cs_cm2(45)  = 0.279d-19
  td_E_eV(46)  =  223.141_8 ; td_cs_cm2(46)  = 0.197d-19
  td_E_eV(47)  =  249.999_8 ; td_cs_cm2(47)  = 0.137d-19
  td_E_eV(48)  =  280.411_8 ; td_cs_cm2(48)  = 0.956d-20
  td_E_eV(49)  =  314.849_8 ; td_cs_cm2(49)  = 0.661d-20
  td_E_eV(50)  =  353.844_8 ; td_cs_cm2(50)  = 0.456d-20
  td_E_eV(51)  =  398.000_8 ; td_cs_cm2(51)  = 0.313d-20
  td_E_eV(52)  =  448.000_8 ; td_cs_cm2(52)  = 0.215d-20
  td_E_eV(53)  =  504.617_8 ; td_cs_cm2(53)  = 0.147d-20
  td_E_eV(54)  =  568.728_8 ; td_cs_cm2(54)  = 0.100d-20
  td_E_eV(55)  =  641.323_8 ; td_cs_cm2(55)  = 0.687d-21
  td_E_eV(56)  =  723.526_8 ; td_cs_cm2(56)  = 0.469d-21
  td_E_eV(57)  =  816.608_8 ; td_cs_cm2(57)  = 0.321d-21
  td_E_eV(58)  =  922.009_8 ; td_cs_cm2(58)  = 0.219d-21

  do n = 1, N_tdp
     td_csv_m3s(n) = 1.0d-4 * td_cs_cm2(n) * sqrt(td_E_eV(n)) * efactor_eV_to_ms
  end do

end subroutine Prepare_He_23P
!
!-------------------------------------
!
real(8) function CSV_He_23P_m3s(E_eV)

  use He_23P
  use Photoelectrons, only : efactor_eV_to_ms

  implicit none

  real(8) E_eV  ! energy of the electron
  integer n
  real(8) left

  if ((E_eV.lt.td_E_eV(1)).or.(E_eV.gt.td_E_eV(N_tdp))) then
     CSV_He_23P_m3s = 0.0_8
     return
  end if

  do n = 1, N_tdp-1
     if ((E_eV.ge.td_E_eV(n)).and.(E_eV.le.td_E_eV(n+1))) then
        left  = (td_E_eV(n+1)-E_eV) / (td_E_eV(n+1) - td_E_eV(n))
        CSV_He_23P_m3s = max(0.0_8, left * td_csv_m3s(n) + (1.0_8 - left) * td_csv_m3s(n+1))
        return
     end if
  end do

end function CSV_He_23P_m3s
