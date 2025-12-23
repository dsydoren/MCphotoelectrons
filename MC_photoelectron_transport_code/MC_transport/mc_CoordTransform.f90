!-------------------------------------------------
!
SUBROUTINE CONVERT_CART_TO_SPHERICAL(rvec, r, theta, phi)

  USE PhysicalConstants

  IMPLICIT NONE

  REAL(8), INTENT(IN) :: rvec(1:3)
  REAL(8), INTENT(OUT) :: r, theta, phi

  REAL(8) rxy

  r = SQRT(rvec(1)**2 + rvec(2)**2 + rvec(3)**2)

  theta = 0.0_8
  phi = 0.0_8

  IF (r.EQ.0.0_8) RETURN

  theta = ACOS(MIN(1.0_8,MAX(rvec(3)/r,-1.0_8)))              ! [0:pi]

  rxy = SQRT(rvec(1)**2 + rvec(2)**2)

  IF (rxy.EQ.0.0_8) RETURN

  IF (rvec(2).GE.0.0_8) THEN
     phi = ACOS(MIN(1.0_8,MAX(rvec(1)/rxy,-1.0_8)))        ! [0:pi]
  ELSE
     IF (rvec(1).GT.0.0_8) THEN
        phi = -ACOS(MIN(1.0_8,MAX(rvec(1)/rxy,-1.0_8)))    ! ]-pi/2:0[
     ELSE 
        phi = pi + ACOS(MIN(1.0_8,MAX(ABS(rvec(1))/rxy,-1.0_8)))  ! ]pi:3pi/2]
     END IF
  END IF

  RETURN

END SUBROUTINE CONVERT_CART_TO_SPHERICAL

!--------------------------------------------
!
subroutine convert_spherical_vector_to_cart(theta, phi, ar, atheta, aphi, ax, ay, az)

  implicit none

  real(8), intent(in) :: theta, phi, ar, atheta, aphi
  real(8), intent(out) :: ax, ay, az

  real(8) sin_theta, cos_theta, sin_phi, cos_phi
  real(8) atr11, atr12, atr13
  real(8) atr21, atr22, atr23
  real(8) atr31, atr32, atr33

  sin_theta = sin(theta)
  cos_theta = cos(theta)
  sin_phi = sin(phi)
  cos_phi = cos(phi)

! if 
!
! r^                  x^
! theta^ = matrix_A * y^
! phi^                z^
!
! where ^ denotes a unitary basis vector in the corresponding direction
! 
! then
!
! (ax,ay,az) = (ar,atheta,aphi) * matrix_A
!
! or
!
! ax                         ar
! ay = matrix_A_transposed * atheta
! az                         aphi
!
! below, atrnm is the element in row n column m of the transposed matrix

  atr11 =  sin_theta * cos_phi
  atr12 =  cos_theta * cos_phi
  atr13 = -sin_phi

  atr21 = sin_theta * sin_phi
  atr22 = cos_theta * sin_phi
  atr23 = cos_phi

  atr31 =  cos_theta
  atr32 = -sin_theta
  atr33 =  0.0_8

  ax = atr11 * ar + atr12 * atheta + atr13 * aphi
  ay = atr21 * ar + atr22 * atheta + atr23 * aphi
  az = atr31 * ar + atr32 * atheta + atr33 * aphi

  return

end subroutine convert_spherical_vector_to_cart

!-----------------------------------------
!
! as described in 
! "Space Physics Coordinate Transfromations: a User Guide" by M. A. Hapgood
! Planet. Space. Sci., Vol.40, pp.711-717, 1992
!
SUBROUTINE Prepare_GSE_to_GEO_transformation

  USE TimeValues
  USE CoordinateTransformations

  IMPLICIT NONE

  REAL, PARAMETER :: deg_to_rad = 0.01745329252   ! factor for conversion of degrees into radians

  REAL theta_deg

! T1^{-1}
  REAL T1inv_11, T1inv_12, T1inv_13
  REAL T1inv_21, T1inv_22, T1inv_23
  REAL T1inv_31, T1inv_32, T1inv_33

  REAL M_deg
  REAL lambda_deg
 
  REAL Z_11, Z_12, Z_13
  REAL Z_21, Z_22, Z_23
  REAL Z_31, Z_32, Z_33

  REAL epsilon_deg
 
  REAL X_11, X_12, X_13
  REAL X_21, X_22, X_23
  REAL X_31, X_32, X_33

! T2
  REAL T2_11, T2_12, T2_13
  REAL T2_21, T2_22, T2_23
  REAL T2_31, T2_32, T2_33

!----------------------------------
! below we need :: 
!
! T0_jcent
! time_ut_h
! MJD_days
!----------------------------------

!###### T1^{-1} ########## GEO -> GEI

! calculate the Greenwich sidereal time
  theta_deg = 100.461 + 36000.770 * T0_jcent + 15.04107 * time_ut_h

!T1=<theta,Z>, T1^{-1}=T1^transposed
  
  T1inv_11 = COS(theta_deg * deg_to_rad) ; T1inv_12 =-SIN(theta_deg * deg_to_rad) ; T1inv_13 = 0.0
  T1inv_21 = SIN(theta_deg * deg_to_rad) ; T1inv_22 = COS(theta_deg * deg_to_rad) ; T1inv_23 = 0.0
  T1inv_31 = 0.0                         ; T1inv_32 = 0.0                         ; T1inv_33 = 1.0

!###### T2 ########## GEI -> GSE

! calculate the Sun's mean anomaly
  M_deg = 357.528 + 35999.050 * T0_jcent + 0.04107 * time_ut_h

! calculate the Sun's ecliptic longitude
  lambda_deg = 280.460 + 36000.772 * T0_jcent + 0.04107 * time_ut_h + &
             & (1.915 - 0.0048 * T0_jcent) * SIN(M_deg * deg_to_rad) + &
             & 0.020 * SIN(2.0 * M_deg * deg_to_rad)

  Z_11 = COS(lambda_deg * deg_to_rad) ; Z_12 = SIN(lambda_deg * deg_to_rad) ; Z_13 = 0.0
  Z_21 =-SIN(lambda_deg * deg_to_rad) ; Z_22 = COS(lambda_deg * deg_to_rad) ; Z_23 = 0.0
  Z_31 = 0.0                          ; Z_32 = 0.0                          ; Z_33 = 1.0

! calculate the obliquity of the ecliptic
  epsilon_deg = 23.439 - 0.013 * T0_jcent

  X_11 = 1.0 ; X_12 = 0.0                           ; X_13 = 0.0
  X_21 = 0.0 ; X_22 = COS(epsilon_deg * deg_to_rad) ; X_23 = SIN(epsilon_deg * deg_to_rad)
  X_31 = 0.0 ; X_32 =-SIN(epsilon_deg * deg_to_rad) ; X_33 = COS(epsilon_deg * deg_to_rad)

! T2=<lambda,Z>*<epsilon,X>=Z*X

  T2_11 = Z_11 * X_11 + Z_12 * X_21 + Z_13 * X_31
  T2_12 = Z_11 * X_12 + Z_12 * X_22 + Z_13 * X_32
  T2_13 = Z_11 * X_13 + Z_12 * X_23 + Z_13 * X_33

  T2_21 = Z_21 * X_11 + Z_22 * X_21 + Z_23 * X_31
  T2_22 = Z_21 * X_12 + Z_22 * X_22 + Z_23 * X_32
  T2_23 = Z_21 * X_13 + Z_22 * X_23 + Z_23 * X_33

  T2_31 = Z_31 * X_11 + Z_32 * X_21 + Z_33 * X_31
  T2_32 = Z_31 * X_12 + Z_32 * X_22 + Z_33 * X_32
  T2_33 = Z_31 * X_13 + Z_32 * X_23 + Z_33 * X_33

!##### GSE -> GEO #### T1 * T2^{-1}) #### GSE -> GEI then GEI -> GEO
  gse_2_geo_11 = T1inv_11 * T2_11 + T1inv_21 * T2_12 + T1inv_31 * T2_13
  gse_2_geo_12 = T1inv_11 * T2_21 + T1inv_21 * T2_22 + T1inv_31 * T2_23
  gse_2_geo_13 = T1inv_11 * T2_31 + T1inv_21 * T2_32 + T1inv_31 * T2_33

  gse_2_geo_21 = T1inv_12 * T2_11 + T1inv_22 * T2_12 + T1inv_32 * T2_13
  gse_2_geo_22 = T1inv_12 * T2_21 + T1inv_22 * T2_22 + T1inv_32 * T2_23
  gse_2_geo_23 = T1inv_12 * T2_31 + T1inv_22 * T2_32 + T1inv_32 * T2_33

  gse_2_geo_31 = T1inv_13 * T2_11 + T1inv_23 * T2_12 + T1inv_33 * T2_13
  gse_2_geo_32 = T1inv_13 * T2_21 + T1inv_23 * T2_22 + T1inv_33 * T2_23
  gse_2_geo_33 = T1inv_13 * T2_31 + T1inv_23 * T2_32 + T1inv_33 * T2_33

END SUBROUTINE Prepare_GSE_To_GEO_transformation

!---------------------------------------------------------
!
subroutine Convert_GSE8_to_GEO8(xgse, ygse, zgse, xgeo, ygeo, zgeo)

  use CoordinateTransformations

  implicit none

  real(8), intent(in)  :: xgse, ygse, zgse
  real(8), intent(out) :: xgeo, ygeo, zgeo

  xgeo = gse_2_geo_11 * xgse + gse_2_geo_12 * ygse + gse_2_geo_13 * zgse
  ygeo = gse_2_geo_21 * xgse + gse_2_geo_22 * ygse + gse_2_geo_23 * zgse
  zgeo = gse_2_geo_31 * xgse + gse_2_geo_32 * ygse + gse_2_geo_33 * zgse

  return

end subroutine Convert_GSE8_to_GEO8
