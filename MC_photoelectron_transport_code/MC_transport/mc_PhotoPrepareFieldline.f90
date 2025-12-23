!-----------------------------------------------
! there is only one field line here
! 
! the field line has to fit below 2000 km - then it is considered as a whole
! field lines reaching higer than 2000 km are terminated at thsi altitude and considered open
! thus the model should be be applied either at low latitudes or at high latitudes where the field lines can be considered open
!
! The velocity distribution is calculated at the spacecraft location. 
! In view of the altitude limit mentioned above, the spacecraft location has to be below 2000 km.
!
! all processes share the same fieldline to improve statistics and reduce noise
! 
! process with rank zero also works on the field line
! intialization of the field line is performed by process with rank 0
! values are set using IRI and MSIS models
! then the copy of the field line is distributed among other processes
!
SUBROUTINE PREPARE_FIELDLINE(i_op)

  USE CoordinateTransformations
  USE field_line
  USE PhysicalConstants
  USE SpacecraftValues
  USE TimeValues
  USE GlobalIndices

  use mpi

  use, intrinsic :: ieee_arithmetic

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i_op

  INTEGER ierr

  INTEGER, PARAMETER :: max_N_fine_fieldline = 20000
  REAL(8), PARAMETER :: step_km = 1.0_8                ! step to create IGRF field line
  REAL(8), PARAMETER :: R_bottom_limit_km = 6470.0_8   ! bottom of ionosphere (alt = 99 km)

  REAL(8), PARAMETER :: max_IRI_alt_km = 2000.0_8      ! maximal altitude of the apex of a field line
                                                       ! the limit corresponds to the altitude limit of the IRI
                                                       ! field lines reaching above this altitude will be terminated here
                                                       ! and treated as open field lines

  REAL, PARAMETER :: alt_thermosphere_km         = 1000.0  ! neutral densities are zero above this altitude (no e-n collisions)

  REAL, PARAMETER :: minimal_alt_of_EUVionization_km = 800.0  ! no primary photoelectrons produced above this altitude
                                                              ! required: minimal_alt_of_EUVionization_km < alt_thermosphere_km

  REAL(8), PARAMETER :: ds_km = 2.0_8   ! step along the direction towards the Sun for calculation of columnar content

  REAL(8), ALLOCATABLE :: x_geo_km(:)
  REAL(8), ALLOCATABLE :: y_geo_km(:)
  REAL(8), ALLOCATABLE :: z_geo_km(:)
  REAL(8), ALLOCATABLE :: B_nT(:)
  INTEGER ALLOC_ERR

  REAL(8) r_km, theta_geo, phi_geo
  REAL(8) colat_geo_deg, lon_geo_deg

  INTEGER count_plus, count_minus

  REAL(8) Bnorth_nT, Beast_nT, Bvert_nT, intensity_nT
  REAL(8) Br_nT, Btheta_nT, Bphi_nT
  REAL(8) Bx_nT, By_nT, Bz_nT
  REAL(8) dir_x, dir_y, dir_z
  REAL(8) rvec(3)

  INTEGER n_skip
  INTEGER count_minus_use, count_plus_use
  INTEGER ipfl, i

  REAL(8) pfl_step_m

  integer yearday
  REAL time_ut_s

  REAL alt_km

  real Nnm3(1:7)  ! densities of 7 neutral species (H, He, N, O, N2, NO, O2), [m^-3], msis_sm
  real TnK        ! neutral temperature, [K]

  logical jf(50)
  
  data jf / .true., .true., .true., .false.,.false., &
       &    .false.,.true., .true., .true., .true., &
       &    .true., .true., .true., .true., .true., &
       &    .true., .true., .true., .true., .true., &
       &    .true., .true., .false.,.true., .true., &
       &    .true., .true., .true., .true., .false., &
       &    .true., .true., .false.,.true., .false., &
       &    .true., .true., .true., .false.,.false., &
       &    .true., .true., .true., .true., .true., &
       &    .true., .false.,.true., .true., .true. /

  real oarr(100)         ! output arrays used by IRI_SUB 2020
  real outf(20,1000)     !

  REAL lat_geo_deg_4, long_geo_deg_4

  REAL(8), ALLOCATABLE :: rbufer(:)
  REAL(8), ALLOCATABLE :: rbufer2(:)

  INTEGER count_bad

  INTEGER istart, iend

  REAL vector_to_sun(1:3)   ! in GEO

  REAL columnar_content_He_m2, columnar_content_O_m2, columnar_content_N2_m2, columnar_content_O2_m2

  REAL msis_NnHe_m3, msis_NnO_m3, msis_NnN2_m3, msis_NnO2_m3

  IF ((orbit_point(i_op)%r_km - R_Earth_km).GE.max_IRI_alt_km) THEN
     PRINT '("orbit point ",i4," has altitude ",f8.2," km")', i_op, orbit_point(i_op)%r_km - R_Earth_km
     PRINT '("this version is limited to altitudes below the maximal IRI altitude of ",f8.2," km")', max_IRI_alt_km
     PRINT '("please, use orbit points with altitudes below this limit")'
     PRINT '("the program will be terminated, sorry")'
     CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
  END IF

! arrays to store field line points calculated with fine spatial resolution
  ALLOCATE( x_geo_km(-max_N_fine_fieldline:max_N_fine_fieldline), STAT = ALLOC_ERR)
  ALLOCATE( y_geo_km(-max_N_fine_fieldline:max_N_fine_fieldline), STAT = ALLOC_ERR)
  ALLOCATE( z_geo_km(-max_N_fine_fieldline:max_N_fine_fieldline), STAT = ALLOC_ERR)
  ALLOCATE(     B_nT(-max_N_fine_fieldline:max_N_fine_fieldline), STAT = ALLOC_ERR)

! set coordinates of points
! start in the spacecraft location
  x_geo_km(0) = orbit_point(i_op)%x_geo_km
  y_geo_km(0) = orbit_point(i_op)%y_geo_km
  z_geo_km(0) = orbit_point(i_op)%z_geo_km

  r_km      = orbit_point(i_op)%r_km
  theta_geo = orbit_point(i_op)%theta_geo
  phi_geo   = orbit_point(i_op)%phi_geo

  colat_geo_deg = orbit_point(i_op)%theta_geo * rad_to_deg   ! IGRF needs this
  lon_geo_deg   = orbit_point(i_op)%phi_geo   * rad_to_deg   ! 

  count_plus = 0

! default assumption
  pfl_closed = .TRUE.

! walk ALONG the direction of the geomagnetic field

  DO
! 0 means main-field values are required
! 2 means geocentric CS (the Earth is a sphere)
     CALL igrf13syn(0, DBLE(year), 2, r_km, colat_geo_deg, lon_geo_deg, Bnorth_nT, Beast_nT, Bvert_nT, intensity_nT)

     Br_nT     = -Bvert_nT
     Btheta_nT = -Bnorth_nT
     Bphi_nT   =  Beast_nT

     CALL convert_spherical_vector_to_cart(theta_geo, phi_geo, Br_nT, Btheta_nT, Bphi_nT, Bx_nT, By_nT, Bz_nT)  ! gives B{xyz}_nT in GEO

     B_nT(count_plus) = sqrt(Bx_nT**2 + By_nT**2 + Bz_nT**2)
     dir_x = Bx_nT / B_nT(count_plus)   !
     dir_y = By_nT / B_nT(count_plus)   ! direction along the field line towards the northern hemisphere 
     dir_z = Bz_nT / B_nT(count_plus)   !

     x_geo_km(count_plus+1) = x_geo_km(count_plus) + step_km * dir_x
     y_geo_km(count_plus+1) = y_geo_km(count_plus) + step_km * dir_y   ! walk ALONG the direction of the geomagnetic field
     z_geo_km(count_plus+1) = z_geo_km(count_plus) + step_km * dir_z

     rvec(1) = x_geo_km(count_plus+1)
     rvec(2) = y_geo_km(count_plus+1)
     rvec(3) = z_geo_km(count_plus+1)
     CALL CONVERT_CART_TO_SPHERICAL(rvec, r_km, theta_geo, phi_geo) 

!print *, rvec, r_km, theta_geo, phi_geo

     alt_km =  r_km-R_Earth_km

     IF (r_km.LT.R_bottom_limit_km) EXIT

     IF (alt_km.GT.max_IRI_alt_km) THEN
        pfl_closed = .FALSE.
        PRINT '("maximal altitude reached while moving northward... must be an open field line in the southern hemisphere")'
        EXIT
     END IF

     count_plus = count_plus + 1

     colat_geo_deg = theta_geo * rad_to_deg
     lon_geo_deg   = phi_geo   * rad_to_deg
  END DO

!print *, count_plus

! return to the spacecraft location
  r_km      = orbit_point(i_op)%r_km
  theta_geo = orbit_point(i_op)%theta_geo
  phi_geo   = orbit_point(i_op)%phi_geo

  colat_geo_deg = orbit_point(i_op)%theta_geo * rad_to_deg   ! IGRF needs this
  lon_geo_deg   = orbit_point(i_op)%phi_geo   * rad_to_deg   ! 

  count_minus = 0

! walk in the OPPOSITE direction

  DO
! 0 means main-field values are required
! 2 means geocentric CS (the Earth is a sphere)
     CALL igrf13syn(0, DBLE(year), 2, r_km, colat_geo_deg, lon_geo_deg, Bnorth_nT, Beast_nT, Bvert_nT, intensity_nT)

     Br_nT     = -Bvert_nT
     Btheta_nT = -Bnorth_nT
     Bphi_nT   =  Beast_nT

     CALL convert_spherical_vector_to_cart(theta_geo, phi_geo, Br_nT, Btheta_nT, Bphi_nT, Bx_nT, By_nT, Bz_nT)  ! gives B{xyz}_nT in GEO

     B_nT(count_minus) = sqrt(Bx_nT**2 + By_nT**2 + Bz_nT**2)
     dir_x = Bx_nT / B_nT(count_minus)   !
     dir_y = By_nT / B_nT(count_minus)   ! direction along the field line towards the northern hemisphere 
     dir_z = Bz_nT / B_nT(count_minus)   !

     x_geo_km(count_minus-1) = x_geo_km(count_minus) - step_km * dir_x
     y_geo_km(count_minus-1) = y_geo_km(count_minus) - step_km * dir_y   ! walk in the direction OPPOSITE to the direction of the geomagnetic field
     z_geo_km(count_minus-1) = z_geo_km(count_minus) - step_km * dir_z

     rvec(1) = x_geo_km(count_minus-1)
     rvec(2) = y_geo_km(count_minus-1)
     rvec(3) = z_geo_km(count_minus-1)
     CALL CONVERT_CART_TO_SPHERICAL(rvec, r_km, theta_geo, phi_geo) 

!print *, rvec, r_km, theta_geo, phi_geo

     alt_km =  r_km - R_Earth_km

     IF (r_km.LT.R_bottom_limit_km) EXIT

     IF (alt_km.GT.max_IRI_alt_km) THEN
        pfl_closed = .FALSE.
        PRINT '("maximal altitude reached while moving southward... must be an open field line in northern hemisphere")'
        EXIT
     END IF

     count_minus = count_minus - 1

     colat_geo_deg = theta_geo * rad_to_deg
     lon_geo_deg   = phi_geo   * rad_to_deg
  END DO

!print *, count_minus

! transfer necessary points to the pfl object, note that we want to keep the point i=0 where the spacecraft is

  n_skip = MAX(1, INT(pfl_desired_dL_km / step_km))

  count_minus_use = INT(-count_minus/n_skip) * n_skip  ! >0
  count_plus_use  = INT(  count_plus/n_skip) * n_skip
  pfl_N_of_points = 1 + count_plus_use/n_skip + count_minus_use/n_skip
  pfl_i_spacecraft = count_minus_use/n_skip + 1               ! index of spacecraft location point in pfl arrays

!print *, n_skip, count_minus, count_plus, count_minus_use, count_plus_use

  ALLOCATE(pfl_point_coords(1:pfl_N_of_points), STAT = ALLOC_ERR)
  ALLOCATE(shared_pfl_point(1:pfl_N_of_points), STAT = ALLOC_ERR)

! set coordinates and magnetic field of field line points
  ipfl = 0
  DO i = -count_minus_use, count_plus_use, n_skip

     ipfl = ipfl + 1

     rvec(1) = x_geo_km(i) / R_Earth_km
     rvec(2) = y_geo_km(i) / R_Earth_km
     rvec(3) = z_geo_km(i) / R_Earth_km

     CALL CONVERT_CART_TO_SPHERICAL(rvec, pfl_point_coords(ipfl)%r, pfl_point_coords(ipfl)%theta_geo, pfl_point_coords(ipfl)%phi_geo)   ! r in units of RE, theta and phi in MAG CS

     pfl_point_coords(ipfl)%x_geo_km = x_geo_km(i)
     pfl_point_coords(ipfl)%y_geo_km = y_geo_km(i)
     pfl_point_coords(ipfl)%z_geo_km = z_geo_km(i)

     shared_pfl_point(ipfl)%geoB_T = B_nT(i) * 1.0d-9

  END DO

! set distances between nodes so that the distance in the first node is zero
  shared_pfl_point(1)%L_m = 0.0_8
!print *, shared_pfl_point(1)%L_m
  pfl_step_m = n_skip * step_km * 1000.0_8
  DO i = 2, pfl_N_of_points
     shared_pfl_point(i)%L_m = shared_pfl_point(i-1)%L_m + pfl_step_m
!print *, i, shared_pfl_point(i)%L_m
  END DO

! partial cleanup
  DEALLOCATE( x_geo_km, STAT = ALLOC_ERR)
  DEALLOCATE( y_geo_km, STAT = ALLOC_ERR)
  DEALLOCATE( z_geo_km, STAT = ALLOC_ERR)
  DEALLOCATE(     B_nT, STAT = ALLOC_ERR)

! calculate cross sections of flux tubes

! assume that h2h3*B=const=1, so that h2h3=1/B
  DO i = 1, pfl_N_of_points
     shared_pfl_point(i)%fts = 1.0_8 / shared_pfl_point(i)%geoB_T
  END DO

! set default zero values of neutral densities and -1 columnar content
  DO i = 1, pfl_N_of_points
     shared_pfl_point(i)%Nn_He_m3 = 0.0_8
     shared_pfl_point(i)%Nn_O_m3 = 0.0_8
     shared_pfl_point(i)%Nn_N2_m3 = 0.0_8
     shared_pfl_point(i)%Nn_O2_m3 = 0.0_8
     shared_pfl_point(i)%columnar_content_He_m2 = -1.0_8
     shared_pfl_point(i)%columnar_content_O_m2 = -1.0_8
     shared_pfl_point(i)%columnar_content_N2_m2 = -1.0_8
     shared_pfl_point(i)%columnar_content_O2_m2 = -1.0_8
  END DO

! set electron density, electron temperature, and neutral densities along the field line -----------------------------------------------------------------------

! prepare yearday in YYDDD format for MSIS
  IF (year.lt.2000) THEN
     yearday = 1000 * (year-1900) + day_of_year
  ELSE
     yearday = 1000 * (year-2000) + day_of_year
  END IF

! universal time in seconds for MSIS
  time_ut_s = time_ut_h * 3600.0

! required by IRI before use of IRI_SUB
!#######     call read_ig_rz
!#######     call readapf107
! required by MSIS
!#######     CALL read_apf107dat(year, month, day_of_month, INT(time_ut_h), msis_Ap, msis_f10p7_pd, msis_f10p7_81)
!#######     call msisinit(switch_legacy=SW)
!####### called outside in the main program now

  print '("calculating plasma and neutral densities along the field line ...")'

  DO i = 1, pfl_N_of_points

     alt_km = (REAL(pfl_point_coords(i)%r)-1.0) * REAL(R_Earth_km)

     lat_geo_deg_4 = MAX(-90.0, 90.0 - REAL(rad_to_deg * pfl_point_coords(i)%theta_geo))
     long_geo_deg_4 = REAL(pfl_point_coords(i)%phi_geo * rad_to_deg)
           
     IF (long_geo_deg_4.LT.0.0)   long_geo_deg_4 = long_geo_deg_4 + 360.0
     IF (long_geo_deg_4.GE.360.0) long_geo_deg_4 = long_geo_deg_4 - 360.0

! note, most input variables are REAL !!!
     IF (alt_km.LE.alt_thermosphere_km) THEN
       CALL msis_geo( yearday, time_ut_s, alt_km, lat_geo_deg_4, long_geo_deg_4, f10p7_81, f10p7_pd, Ap, Nnm3, TnK )

        shared_pfl_point(i)%Nn_He_m3 = Nnm3(2)
        shared_pfl_point(i)%Nn_O_m3  = Nnm3(4)
        shared_pfl_point(i)%Nn_N2_m3 = Nnm3(5)
        shared_pfl_point(i)%Nn_O2_m3 = Nnm3(7)
     END IF

! 0 in the second argument means "use geographic coordinates"
! "alt_km, alt_km, 1.0" means start/end altitude and step
     call IRI_SUB(jf, 0, lat_geo_deg_4, long_geo_deg_4, year, -day_of_year, time_ut_h+25.0, alt_km, alt_km, 1.0, outf, oarr)

! IRI electron density is outf(1,1), in 1/m^3
! IRI electron temperature is outf(4,1), in K

     shared_pfl_point(i)%Ne_m3 = outf(1,1)
     shared_pfl_point(i)%Te_K  = outf(4,1)

!print *, i, alt_km, lat_geo_deg_4, long_geo_deg_4, shared_pfl_point(i)%Ne_m3, shared_pfl_point(i)%Te_K

  END DO   !###   DO i = 1, pfl_N_of_points

! clean PLASMA DENSITY from possible NaNs ##########
  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
  IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT=ALLOC_ERR)
  ALLOCATE(rbufer(1:pfl_N_of_points), STAT = ALLOC_ERR)
  ALLOCATE(rbufer2(1:pfl_N_of_points), STAT = ALLOC_ERR)
  DO i = 1, pfl_N_of_points
     rbufer(i) = shared_pfl_point(i)%Ne_m3
     rbufer2(i) = shared_pfl_point(i)%L_m   ! need the parallel coordinate for interpolation
  END DO
  CALL Clean_array_from_NaNs_log(1, pfl_N_of_points, rbufer, rbufer2, count_bad)   ! uses linear interpolation of logarithms to fill gaps
  IF (count_bad.GT.0) THEN 
     PRINT '("Clean_array_from_NaNs_log-1 found and cleaned ",i6," bad Ne values")', count_bad
     DO i = 1, pfl_N_of_points
        shared_pfl_point(i)%Ne_m3 = rbufer(i)
     END DO
  END IF

! clean ELECTRON TEMPERATURE from possible NaNs ##########
  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
  IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT=ALLOC_ERR)
  ALLOCATE(rbufer(1:pfl_N_of_points), STAT = ALLOC_ERR)
  ALLOCATE(rbufer2(1:pfl_N_of_points), STAT = ALLOC_ERR)
  DO i = 1, pfl_N_of_points
     rbufer(i) = shared_pfl_point(i)%Te_K
     rbufer2(i) = shared_pfl_point(i)%L_m   ! need the parallel coordinate for interpolation
  END DO
  CALL Clean_array_from_NaNs(1, pfl_N_of_points, rbufer, rbufer2, count_bad)   ! uses linear interpolation of logarithms to fill gaps
  IF (count_bad.GT.0) THEN 
     PRINT '("Clean_array_from_NaNs found and cleaned ",i6," bad Te values")', count_bad
     DO i = 1, pfl_N_of_points
        shared_pfl_point(i)%Te_K = rbufer(i)
     END DO
  END IF

! ensure continuity of NEUTRAL DENSITIES for altitudes below threshold

  IF (pfl_closed) THEN
! the segment near the beginning of the field line
     istart = 1
     iend = pfl_N_of_points
     DO i = 1, pfl_N_of_points-1
        alt_km = (REAL(pfl_point_coords(i+1)%r) - 1.0) * REAL(R_Earth_km)
        IF (alt_km.GT.alt_thermosphere_km) THEN
           iend = i
           EXIT
        END IF
     END DO
! even if the whole field line is below the threshold altitude, still istart < iend 
     IF (ALLOCATED(rbufer))  DEALLOCATE(rbufer,  STAT = ALLOC_ERR)
     IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT = ALLOC_ERR)
     ALLOCATE( rbufer(istart:iend), STAT = ALLOC_ERR)
     ALLOCATE(rbufer2(istart:iend), STAT = ALLOC_ERR)
! He
     DO i = istart, iend
        rbufer(i)  = shared_pfl_point(i)%Nn_He_m3
        rbufer2(i) = shared_pfl_point(i)%L_m   ! needed for interpolation
     END DO
     CALL Clean_array_from_NaNs_log(istart, iend, rbufer, rbufer2, count_bad)
     IF (count_bad.GT.0) THEN
        PRINT '("Clean_array_from_NaNs_log-3 found and cleaned ",i6," bad Nn_He_m3 values")', count_bad
        DO i = istart, iend
           shared_pfl_point(i)%Nn_He_m3 = rbufer(i)
        END DO
     END IF
! O
     DO i = istart, iend
        rbufer(i)  = shared_pfl_point(i)%Nn_O_m3
        rbufer2(i) = shared_pfl_point(i)%L_m   ! needed for interpolation
     END DO
     CALL Clean_array_from_NaNs_log(istart, iend, rbufer, rbufer2, count_bad)
     IF (count_bad.GT.0) THEN
        PRINT '("Clean_array_from_NaNs_log-4 found and cleaned ",i6," bad Nn_O_m3 values")', count_bad
        DO i = istart, iend
           shared_pfl_point(i)%Nn_O_m3 = rbufer(i)
        END DO
     END IF
! N2
     DO i = istart, iend
        rbufer(i)  = shared_pfl_point(i)%Nn_N2_m3
        rbufer2(i) = shared_pfl_point(i)%L_m   ! needed for interpolation
     END DO
     CALL Clean_array_from_NaNs_log(istart, iend, rbufer, rbufer2, count_bad)
     IF (count_bad.GT.0) THEN
        PRINT '("Clean_array_from_NaNs_log-5 found and cleaned ",i6," bad Nn_N2_m3 values")', count_bad
        DO i = istart, iend
           shared_pfl_point(i)%Nn_N2_m3 = rbufer(i)
        END DO
     END IF
! O2
     DO i = istart, iend
        rbufer(i)  = shared_pfl_point(i)%Nn_O2_m3
        rbufer2(i) = shared_pfl_point(i)%L_m   ! needed for interpolation
     END DO
     CALL Clean_array_from_NaNs_log(istart, iend, rbufer, rbufer2, count_bad)
     IF (count_bad.GT.0) THEN 
        PRINT '("Clean_array_from_NaNs_log-6 found and cleaned ",i6," bad Nn_O2_m3 values")', count_bad
        DO i = istart, iend
           shared_pfl_point(i)%Nn_O2_m3 = rbufer(i)
        END DO
     END IF

! the segment near the end of the field line only
     istart = pfl_N_of_points+1
     iend = pfl_N_of_points
     DO i = pfl_N_of_points, 2, -1
        alt_km = (REAL(pfl_point_coords(i-1)%r) - 1.0) * REAL(R_Earth_km)
        IF (alt_km.GT.alt_thermosphere_km) THEN
           istart = i
           EXIT
        END IF
     END DO
! note, the case when the whole field line is below the threshold altitude has been already processed above
! repeating this step is avoided here since if the whole field line is below the threshold altitude, then istart > iend
     IF (iend.GE.istart) THEN 
        IF (ALLOCATED(rbufer))  DEALLOCATE(rbufer,  STAT = ALLOC_ERR)
        IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT = ALLOC_ERR)
        ALLOCATE( rbufer(istart:iend), STAT = ALLOC_ERR)
        ALLOCATE(rbufer2(istart:iend), STAT = ALLOC_ERR)
! He
        DO i = istart, iend
           rbufer(i)  = shared_pfl_point(i)%Nn_He_m3
           rbufer2(i) = shared_pfl_point(i)%L_m   ! needed for interpolation
        END DO
        CALL Clean_array_from_NaNs_log(istart, iend, rbufer, rbufer2, count_bad)
        IF (count_bad.GT.0) THEN
           PRINT '("Clean_array_from_NaNs_log-7 found and cleaned ",i6," bad Nn_He_m3 values")', count_bad
           DO i = istart, iend
              shared_pfl_point(i)%Nn_He_m3 = rbufer(i)
           END DO
        END IF
! O
        DO i = istart, iend
           rbufer(i)  = shared_pfl_point(i)%Nn_O_m3
           rbufer2(i) = shared_pfl_point(i)%L_m   ! needed for interpolation
        END DO
        CALL Clean_array_from_NaNs_log(istart, iend, rbufer, rbufer2, count_bad)
        IF (count_bad.GT.0) THEN
           PRINT '("Clean_array_from_NaNs_log-8 found and cleaned ",i6," bad Nn_O_m3 values")', count_bad
           DO i = istart, iend
              shared_pfl_point(i)%Nn_O_m3 = rbufer(i)
           END DO
        END IF
! N2
        DO i = istart, iend
           rbufer(i)  = shared_pfl_point(i)%Nn_N2_m3
           rbufer2(i) = shared_pfl_point(i)%L_m   ! needed for interpolation
        END DO
        CALL Clean_array_from_NaNs_log(istart, iend, rbufer, rbufer2, count_bad)
        IF (count_bad.GT.0) THEN
           PRINT '("Clean_array_from_NaNs_log-9 found and cleaned ",i6," bad Nn_N2_m3 values")', count_bad
           DO i = istart, iend
              shared_pfl_point(i)%Nn_N2_m3 = rbufer(i)
           END DO
        END IF
! O2
        DO i = istart, iend
           rbufer(i)  = shared_pfl_point(i)%Nn_O2_m3
           rbufer2(i) = shared_pfl_point(i)%L_m   ! needed for interpolation
        END DO
        CALL Clean_array_from_NaNs_log(istart, iend, rbufer, rbufer2, count_bad)
        IF (count_bad.GT.0) THEN
           PRINT '("Clean_array_from_NaNs_log-10 found and cleaned ",i6," bad Nn_O2_m3 values")', count_bad
           DO i = istart, iend
              shared_pfl_point(i)%Nn_O2_m3 = rbufer(i)
           END DO
        END IF
     END IF   !###   IF (iend.GE.istart) THEN
     
  ELSE   !###   IF (pfl_closed) THEN
! open field line

     istart = pfl_N_of_points+1
     iend = 0
! scan forward to find the first element
     DO i = 1, pfl_N_of_points
        alt_km = (REAL(pfl_point_coords(i)%r) - 1.0) * REAL(R_Earth_km)
        IF (alt_km.LE.alt_thermosphere_km) THEN
           istart = i
           EXIT
        END IF
     END DO
! scan backward to find the last element
     DO i = pfl_N_of_points, 1, -1
        alt_km = (REAL(pfl_point_coords(i)%r) - 1.0) * REAL(R_Earth_km)
        IF (alt_km.LE.alt_thermosphere_km) THEN
           iend = i
           EXIT
        END IF
     END DO
     IF (istart.GT.iend) THEN
! just in case, me being paranoidal
        PRINT '("Error-2")'
        CALL MPI_ABORT(MPI_COMM_WORLD, 22, ierr)
     END IF
     IF (ALLOCATED(rbufer))  DEALLOCATE(rbufer,  STAT = ALLOC_ERR)
     IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT = ALLOC_ERR)
     ALLOCATE( rbufer(istart:iend), STAT = ALLOC_ERR)
     ALLOCATE(rbufer2(istart:iend), STAT = ALLOC_ERR)
! He
     DO i = istart, iend
        rbufer(i)  = shared_pfl_point(i)%Nn_He_m3
        rbufer2(i) = shared_pfl_point(i)%L_m   ! needed for interpolation
     END DO
     CALL Clean_array_from_NaNs_log(istart, iend, rbufer, rbufer2, count_bad)
     IF (count_bad.GT.0) THEN
        PRINT '("Clean_array_from_NaNs_log-11 found and cleaned ",i6," bad Nn_He_m3 values")', count_bad
        DO i = istart, iend
           shared_pfl_point(i)%Nn_He_m3 = rbufer(i)
        END DO
     END IF
! O
     DO i = istart, iend
        rbufer(i)  = shared_pfl_point(i)%Nn_O_m3
        rbufer2(i) = shared_pfl_point(i)%L_m   ! needed for interpolation
     END DO
     CALL Clean_array_from_NaNs_log(istart, iend, rbufer, rbufer2, count_bad)
     IF (count_bad.GT.0) THEN
        PRINT '("Clean_array_from_NaNs_log-12 found and cleaned ",i6," bad Nn_O_m3 values")', count_bad
        DO i = istart, iend
           shared_pfl_point(i)%Nn_O_m3 = rbufer(i)
        END DO
     END IF
! N2
     DO i = istart, iend
        rbufer(i)  = shared_pfl_point(i)%Nn_N2_m3
        rbufer2(i) = shared_pfl_point(i)%L_m   ! needed for interpolation
     END DO
     CALL Clean_array_from_NaNs_log(istart, iend, rbufer, rbufer2, count_bad)
     IF (count_bad.GT.0) THEN
        PRINT '("Clean_array_from_NaNs_log-13 found and cleaned ",i6," bad Nn_N2_m3 values")', count_bad
        DO i = istart, iend
           shared_pfl_point(i)%Nn_N2_m3 = rbufer(i)
        END DO
     END IF
! O2
     DO i = istart, iend
        rbufer(i)  = shared_pfl_point(i)%Nn_O2_m3
        rbufer2(i) = shared_pfl_point(i)%L_m   ! needed for interpolation
     END DO
     CALL Clean_array_from_NaNs_log(istart, iend, rbufer, rbufer2, count_bad)
     IF (count_bad.GT.0) THEN 
        PRINT '("Clean_array_from_NaNs_log-14 found and cleaned ",i6," bad Nn_O2_m3 values")', count_bad
        DO i = istart, iend
           shared_pfl_point(i)%Nn_O2_m3 = rbufer(i)
        END DO
     END IF

  END IF   !###   IF (pfl_closed) THEN

! cleanup
  IF (ALLOCATED(rbufer))  DEALLOCATE(rbufer,  STAT = ALLOC_ERR)
  IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT = ALLOC_ERR)

! calculate the columnar content ----------------------------------------------------------------------------------------------------
! note that yearday and UT in seconds are already calculated above

! vector to the Sun in GEO is the first column of the GSE -> GEO transformation 

  vector_to_sun(1) = gse_2_geo_11
  vector_to_sun(2) = gse_2_geo_21
  vector_to_sun(3) = gse_2_geo_31

  print '("calculating columnar content along the field line ...")'

  DO i = 1, pfl_N_of_points

     alt_km = (REAL(pfl_point_coords(i)%r)-1.0) * REAL(R_Earth_km)
     IF (alt_km.GT.minimal_alt_of_EUVionization_km) CYCLE

! note, if a point is above the EUV ionization threshold altitude, 
! the columnar content keeps its default value of -1 which prevents photoelectron emission from this point

     columnar_content_He_m2 = shared_pfl_point(i)%Nn_He_m3
     columnar_content_O_m2  = shared_pfl_point(i)%Nn_O_m3
     columnar_content_N2_m2 = shared_pfl_point(i)%Nn_N2_m3
     columnar_content_O2_m2 = shared_pfl_point(i)%Nn_O2_m3

     rvec(1) = pfl_point_coords(i)%x_geo_km
     rvec(2) = pfl_point_coords(i)%y_geo_km
     rvec(3) = pfl_point_coords(i)%z_geo_km

! walk along the ray
     DO
        rvec = rvec + ds_km * vector_to_sun   ! GEO
           
        CALL CONVERT_CART_TO_SPHERICAL(rvec, r_km, theta_geo, phi_geo)   !#### also see convert_cart_xyz_to_spherical above, maybe use just one??????
        alt_km = REAL(r_km - R_Earth_km)
! ray to Sun exits thermosphere, columnar content will not change now
        IF (alt_km.GT.alt_thermosphere_km) EXIT

! ray to Sun passes through dense atmosphere or the Earth, no direct EUV effects in the origin point
        IF (r_km.LT.R_bottom_limit_km) THEN
           columnar_content_He_m2 = -1.0
           columnar_content_O_m2  = -1.0
           columnar_content_N2_m2 = -1.0
           columnar_content_O2_m2 = -1.0
           EXIT
        END IF

        msis_NnHe_m3 = 0.0
        msis_NnO_m3  = 0.0
        msis_NnN2_m3 = 0.0
        msis_NnO2_m3 = 0.0

        lat_geo_deg_4 = MAX(-90.0, 90.0 - REAL(rad_to_deg * theta_geo))
        long_geo_deg_4 = REAL(phi_geo * rad_to_deg)
           
        IF (long_geo_deg_4.LT.0.0)   long_geo_deg_4 = long_geo_deg_4 + 360.0
        IF (long_geo_deg_4.GE.360.0) long_geo_deg_4 = long_geo_deg_4 - 360.0

        CALL msis_geo( yearday, time_ut_s, alt_km, lat_geo_deg_4, long_geo_deg_4, f10p7_81, f10p7_pd, Ap, Nnm3, TnK )   ! note that there is also msis_sm ##########

        IF (ieee_is_finite(Nnm3(2))) msis_NnHe_m3 = Nnm3(2)
        IF (ieee_is_finite(Nnm3(4))) msis_NnO_m3  = Nnm3(4)
        IF (ieee_is_finite(Nnm3(5))) msis_NnN2_m3 = Nnm3(5)
        IF (ieee_is_finite(Nnm3(7))) msis_NnO2_m3 = Nnm3(7)

        columnar_content_He_m2 = columnar_content_He_m2 + msis_NnHe_m3 !Nnm3(2)
        columnar_content_O_m2  = columnar_content_O_m2  + msis_NnO_m3  !Nnm3(4)
        columnar_content_N2_m2 = columnar_content_N2_m2 + msis_NnN2_m3 !Nnm3(5)
        columnar_content_O2_m2 = columnar_content_O2_m2 + msis_NnO2_m3 !Nnm3(7)

     END DO   ! end of cycle advancing one ray

! account for previously omitted factors
     IF (columnar_content_He_m2.GE.0.0) columnar_content_He_m2 = columnar_content_He_m2 * ds_km * 1.0e3
     IF (columnar_content_O_m2.GE.0.0)  columnar_content_O_m2  = columnar_content_O_m2  * ds_km * 1.0e3
     IF (columnar_content_N2_m2.GE.0.0) columnar_content_N2_m2 = columnar_content_N2_m2 * ds_km * 1.0e3
     IF (columnar_content_O2_m2.GE.0.0) columnar_content_O2_m2 = columnar_content_O2_m2 * ds_km * 1.0e3

     shared_pfl_point(i)%columnar_content_He_m2 = columnar_content_He_m2
     shared_pfl_point(i)%columnar_content_O_m2  = columnar_content_O_m2
     shared_pfl_point(i)%columnar_content_N2_m2 = columnar_content_N2_m2
     shared_pfl_point(i)%columnar_content_O2_m2 = columnar_content_O2_m2

  END DO   !###   DO i = 1, pfl_N_of_points

END SUBROUTINE PREPARE_FIELDLINE

!-----------------------------------------------------
!
SUBROUTINE Clean_array_from_NaNs_log(i_start, i_end, arr, x, count)

  USE, INTRINSIC :: ieee_arithmetic

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i_start, i_end
  REAL(8), INTENT(INOUT) :: arr(i_start:i_end)
  REAL(8), INTENT(IN) :: x(i_start:i_end)
  INTEGER, INTENT(OUT) :: count

  INTEGER icheck
  INTEGER ibegin_correct
  INTEGER iend_correct

  INTEGER i
  REAL(8) log_left, log_right, slope, log_of_arr

  icheck = i_start-1
  count = 0

  DO WHILE (icheck.LT.i_end)

     icheck = icheck + 1

     IF (arr(icheck).GT.0.0_8) CYCLE

! found a bad value (negative or NaN)
! find continuous interval of bad values
     ibegin_correct = icheck
     iend_correct = icheck
     DO WHILE (icheck.LT.i_end)
        icheck = icheck+1
        IF (arr(icheck).GT.0.0_8) EXIT
        iend_correct = icheck
     END DO

! correct bad values within the interval
     IF (iend_correct.EQ.i_end) THEN
! interval adjacent to the end of the field line
        DO i = ibegin_correct, iend_correct
           arr(i) = arr(ibegin_correct-1)
        END DO
     ELSE IF (ibegin_correct.EQ.i_start) THEN
! interval adjacent to the beginning of the field line
        DO i = ibegin_correct, iend_correct
           arr(i) = arr(iend_correct+1)
        END DO
     ELSE
! interval begins and ends inside the field line (excluding the end points)
        log_left = LOG(arr(ibegin_correct-1))
        log_right = LOG(arr(iend_correct+1))
        slope = (log_right - log_left) / (x(iend_correct+1) - x(ibegin_correct-1))
        DO i = ibegin_correct, iend_correct
           log_of_arr = log_left + slope * (x(i) - x(ibegin_correct-1))
           arr(i) = EXP(log_of_arr)
        END DO
     END IF
     count = count + (iend_correct - ibegin_correct + 1)
  END DO

END SUBROUTINE Clean_array_from_NaNs_log

!-----------------------------------------------------
!
SUBROUTINE Clean_array_from_NaNs(i_start, i_end, arr, x, count)

  USE, INTRINSIC :: ieee_arithmetic

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i_start, i_end
  REAL(8), INTENT(INOUT) :: arr(i_start:i_end)
  REAL(8), INTENT(IN) :: x(i_start:i_end)
  INTEGER, INTENT(OUT) :: count

  INTEGER icheck
  INTEGER ibegin_correct
  INTEGER iend_correct

  INTEGER i
  REAL(8) arr_left, arr_right, slope

  icheck = i_start-1
  count = 0

  DO WHILE (icheck.LT.i_end)

     icheck = icheck + 1

     IF (arr(icheck).GT.0.0_8) CYCLE

! found a bad value (negative or NaN)
! find continuous interval of bad values
     ibegin_correct = icheck
     iend_correct = icheck
     DO WHILE (icheck.LT.i_end)
        icheck = icheck+1
        IF (arr(icheck).GT.0.0_8) EXIT
        iend_correct = icheck
     END DO

! correct bad values within the interval
     IF (iend_correct.EQ.i_end) THEN
! interval adjacent to the end of the field line
        DO i = ibegin_correct, iend_correct
           arr(i) = arr(ibegin_correct-1)
        END DO
     ELSE IF (ibegin_correct.EQ.i_start) THEN
! interval adjacent to the beginning of the field line
        DO i = ibegin_correct, iend_correct
           arr(i) = arr(iend_correct+1)
        END DO
     ELSE
! interval begins and ends inside the field line (excluding the end points)
        arr_left = arr(ibegin_correct-1)
        arr_right = arr(iend_correct+1)
        slope = (arr_right - arr_left) / (x(iend_correct+1) - x(ibegin_correct-1))
        DO i = ibegin_correct, iend_correct
           arr(i) = arr_left + slope * (x(i) - x(ibegin_correct-1))
        END DO
     END IF
     count = count + (iend_correct - ibegin_correct + 1)
  END DO

END SUBROUTINE Clean_array_from_NaNs

