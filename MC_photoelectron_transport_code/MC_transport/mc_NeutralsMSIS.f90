
!-------------------------------------
!
subroutine read_apf107dat(year, month, day_of_month, ihr, msisAp, f107, f107_pd, f107_81)

  implicit none

  integer, intent(in) :: year, month, day_of_month, ihr
  real, intent(out) :: msisAp(7), f107, f107_pd, f107_81

  logical exists
  integer ierr

  integer myyear, mymonth, myday, myAp(-23:8), myApd, idummy
  real myf107(1:2), myf10781, rdummy

  integer indx_ihr, i

  msisAp=0.0
  msisAp(2) = -5.0
  f107_pd=90.0
  f107_81=90.0

  inquire(file = 'apf107.dat', exist=exists)

  if (.not.exists) then
     print '("file apf107.dat not found, return defaults")'
     return
  end if

  open (13, file = 'apf107.dat')

  myAp=0.0
  myf107 = 0.0
  myf10781 = 0.0

  do 
     read (13, '(3(i3),9(i3),i3,3(f5.1))', iostat=ierr) myyear, mymonth, myday, myAp(1:8), myApd, idummy, myf107(2), myf10781, rdummy
     if (ierr.ne.0) exit

     if (myyear.lt.58) then
        myyear = myyear + 2000
     else
        myyear = myyear + 1900
     end if

     if ((myyear.eq.year).and.(mymonth.eq.month).and.(myday.eq.day_of_month)) then

! found the day, prepare the output

        f107    = myf107(2)
        f107_pd = myf107(1)
        f107_81 = myf10781
        msisAp(1) = real(myApd)
        indx_ihr = min(max(1,1+ihr/3),8)

        msisAp(2) = myAp(indx_ihr) 
        msisAp(3) = myAp(indx_ihr-1) 
        msisAp(4) = myAp(indx_ihr-2) 
        msisAp(5) = myAp(indx_ihr-3)

        msisAp(6) = 0.0
        msisAp(7) = 0.0
        do i = 1, 8
           msisAp(6) = msisAp(6) + myAp(indx_ihr-3-i)
           msisAp(7) = msisAp(7) + myAp(indx_ihr-11-i)
        end do
        msisAp(6) = msisAp(6) * 0.125
        msisAp(7) = msisAp(7) * 0.125
           
        close (13, status = 'keep')
        
        return
     end if      

! day not found yet, keep searching

! shift the data
     myAp(-23:-16) = myAp(-15:-8)
     myAp(-15:-8) = myAp(-7:0)
     myAp(-7:0) = myAp(1:8)
     myAp(1:8) = 0.0
     myf107(1) = myf107(2)
          
  end do
  
  close (13, status = 'keep')
! if we are here, the day has not been found

  print '("read_apf107dat :: date yy/mm/dd ",i2,"/",i2,"/",i2," not found in apf107.dat, return defaults")', year, month, day_of_month
  
end subroutine read_apf107dat

!----------------------------------------------------
!
! modified to use NRL MSIS 2.1
!
subroutine msis_geo( yearday, time_ut_s, alt_km, lat_geo_deg, long_geo_deg, f10p7_81, f10p7_pd, Ap, Nnm3, TnK )

!  use msis_init, only : msisinit

  implicit none

! input
  INTEGER yearday             ! year and day as YYDDD
  REAL    time_ut_s           ! universal time [s]

  REAL alt_km                 ! altitude [km]
  REAL lat_geo_deg            ! geographic latitude [deg]
  REAL long_geo_deg           ! geographic longtude [deg]

  REAL f10p7_81               ! F10.7 flux averaged over 81 days CENTERED on the requested date [SFU, or solar flux units]. 1 SFU = 10^-22 W/m^2 Hz
  REAL f10p7_pd               ! F10.7 flux for the day before the requested day [SFU]
  REAL Ap(1:7)                ! (1) DAILY AP
                              ! (2) 3 HR AP INDEX FOR CURRENT TIME
                              ! (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
                              ! (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
                              ! (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
                              ! (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR TO CURRENT TIME
                              ! (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR TO CURRENT TIME

! output
  real Nnm3(7)                ! number densities of different neutral species [m^-3]
  real TnK                    ! temperature of neutrals [K]

! parameters calculated below for the NRL-MSISE-00 procedure GTD7

  REAL    time_local_h           ! local time [h] 

  REAL    D(1:10)                ! D(1) - HE NUMBER DENSITY  [cm^-3]
                                 ! D(2) - O NUMBER DENSITY   [cm^-3]  !###
                                 ! D(3) - N2 NUMBER DENSITY  [cm^-3]  !###
                                 ! D(4) - O2 NUMBER DENSITY  [cm^-3]  !###
                                 ! D(5) - AR NUMBER DENSITY  [cm^-3]
                                 ! D(6) - TOTAL MASS DENSITY [g/cm^3]
                                 ! D(7) - H NUMBER DENSITY   [cm^-3]  !###
                                 ! D(8) - N NUMBER DENSITY   [cm^-3]  !###
                                 ! D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)   NEW!!!
                                 ! D(10) - NO number density (cm^{-3}) !!!  new, introduced in msis 2.1

  REAL    T(1:2)                 ! T(1) - EXOSPHERIC TEMPERATURE [K]
                                 ! T(2) - TEMPERATURE AT ALT     [K]  !###

!  real SW(25)          ! MSIS switches
!!!  data SW/8*1.,-1.,16*1./
!  SW = 1.0
!!!  SW(9) = -1.0
!!!  CALL TSELEC(SW)
!  call msisinit(switch_legacy=SW)

! local time in hours, decimal
  time_local_h = time_ut_s / 3600.0 + long_geo_deg / 15.0
  if (time_local_h.LT.0.0) time_local_h = time_local_h + 24.0
  if (time_local_h.GE.24.0)  time_local_h = time_local_h - 24.0

!########## AP is an array of size 7 #######!!!!!!!!!!!!

!  CALL GTD7(yearday, time_ut_s, alt_km, lat_geo_deg, long_geo_deg, time_local_h, f10p7_81, f10p7_pd, Ap, 48, D, T)
  CALL gtd8d(yearday, time_ut_s, alt_km, lat_geo_deg, long_geo_deg, time_local_h, f10p7_81, f10p7_pd, Ap, 48, D, T)

! initiate the ionospheric data structure:
! neutral densities [m^-3]
  Nnm3(1) = 1.0e6 * D(7)    ! hydrogen
  Nnm3(2) = 1.0e6 * D(1)    ! helium
  Nnm3(3) = 1.0e6 * D(8)    ! nitrogen
  Nnm3(4) = 1.0e6 * D(2)    ! oxygen
  Nnm3(5) = 1.0e6 * D(3)    ! nitrogen molecular
  Nnm3(6) = 1.0e6 * D(10)   ! NO
  Nnm3(7) = 1.0e6 * D(4)    ! oxygen molecular
! neutral temperature [Kelvins]
  TnK = T(2)
! calculate the nitric oxide density [m^-3] as described in [Mitra, J.Atmos.Terr.Phys., 30, 1065 (1968)]
! ##### NOTE, alternative formulas are given in page 45 of STEP handbook
!  Nnm3(6) = 0.4 * EXP(-3700.0 / TnK) * Nnm3(7) + 5.0E-7 * Nnm3(4)

  return

end subroutine msis_geo
