!-----------------------------------------------
! there is only one field line here
! both eedf and evdf are calculated
! all processes share the same fieldline to improv statistics and reduce noise
! work sharing is off
! process with rank zero also works on the field line
! intialization of the field line is performed by process with rank 0
! values are set using IRI, GCPM, MSIS models
! then the copy of the field line is distributed among other processes
!
SUBROUTINE CALCULATE_PHOTOELECTRON_ENERGY_SPECTRA

  USE field_line
  USE photoelectrons
  USE ParallelOperationValues

  USE MPI

  use, intrinsic :: ieee_arithmetic

  IMPLICIT NONE

  INTEGER ierr
  INTEGER stattus(MPI_STATUS_SIZE)
  INTEGER request

  INTEGER ibufer4(4)

  REAL(8), ALLOCATABLE :: rbufer(:)
  REAL(8), ALLOCATABLE :: rbufer2(:)
  INTEGER ALLOC_ERR

  INTEGER buflen
  INTEGER pos, i, j, nr

  INTEGER particle_count
  INTEGER(8) step_count
  INTEGER mirror_count

  INTEGER jsol, jchannel

  REAL(8) temp

  IF (Rank_of_process.EQ.0) THEN
! distribute field line data to other processes

     ibufer4(1) = pfl_N_of_points
     IF (pfl_closed) THEN
        ibufer4(2) = 1
     ELSE
        ibufer4(2) = 0
     END IF
     ibufer4(3) = pfl_i_spacecraft
     ibufer4(4) = N_ranges

     CALL MPI_BCAST(ibufer4, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

     buflen = 13*pfl_N_of_points + N_ranges*4
     ALLOCATE(rbufer(1:buflen), STAT = ALLOC_ERR)
     pos=1
     DO i = 1, pfl_N_of_points
        rbufer(pos)    = shared_pfl_point(i)%L_m
        rbufer(pos+1)  = shared_pfl_point(i)%fts
        rbufer(pos+2)  = shared_pfl_point(i)%columnar_content_He_m2
        rbufer(pos+3)  = shared_pfl_point(i)%columnar_content_O_m2
        rbufer(pos+4)  = shared_pfl_point(i)%columnar_content_N2_m2
        rbufer(pos+5)  = shared_pfl_point(i)%columnar_content_O2_m2
        rbufer(pos+6)  = shared_pfl_point(i)%Nn_He_m3
        rbufer(pos+7)  = shared_pfl_point(i)%Nn_O_m3
        rbufer(pos+8)  = shared_pfl_point(i)%Nn_N2_m3
        rbufer(pos+9)  = shared_pfl_point(i)%Nn_O2_m3
        rbufer(pos+10) = shared_pfl_point(i)%geoB_T
        rbufer(pos+11) = shared_pfl_point(i)%Te_K
        rbufer(pos+12) = shared_pfl_point(i)%Ne_m3
        pos = pos+13
     END DO
     DO nr = 1, N_ranges
        rbufer(pos)    = min_energy_eV_range(nr)
        rbufer(pos+1)  = max_energy_eV_range(nr)
        rbufer(pos+2)  = min_pitch_angle_range(nr)
        rbufer(pos+3)  = max_pitch_angle_range(nr)
        pos = pos+4
     END DO
     CALL MPI_BCAST(rbufer, buflen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  ELSE     !###   IF (Rank_of_process.EQ.0) THEN
! receive the field line

     CALL MPI_BCAST(ibufer4, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

     pfl_N_of_points = ibufer4(1)
     pfl_closed = .TRUE.
     IF (ibufer4(2).EQ.0) pfl_closed = .FALSE.
     pfl_i_spacecraft = ibufer4(3)
     N_ranges = ibufer4(4)

! allocate field line
     ALLOCATE(shared_pfl_point(1:pfl_N_of_points), STAT = ALLOC_ERR)
     buflen = 13*pfl_N_of_points + N_ranges*4
     ALLOCATE(rbufer(1:buflen), STAT = ALLOC_ERR)
        
     IF (ALLOCATED(min_energy_eV_range)) DEALLOCATE(min_energy_eV_range)
     IF (ALLOCATED(max_energy_eV_range)) DEALLOCATE(max_energy_eV_range)
     IF (ALLOCATED(min_pitch_angle_range)) DEALLOCATE(min_pitch_angle_range)
     IF (ALLOCATED(max_pitch_angle_range)) DEALLOCATE(max_pitch_angle_range)

     IF (N_ranges.GE.1) THEN
        ALLOCATE(min_energy_eV_range(1:N_ranges), STAT = ALLOC_ERR)
        ALLOCATE(max_energy_eV_range(1:N_ranges), STAT = ALLOC_ERR)
        ALLOCATE(min_pitch_angle_range(1:N_ranges), STAT = ALLOC_ERR)
        ALLOCATE(max_pitch_angle_range(1:N_ranges), STAT = ALLOC_ERR)
     END IF

     CALL MPI_BCAST(rbufer, buflen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
! translate the message
     pos=1
     DO i = 1, pfl_N_of_points
        shared_pfl_point(i)%L_m                    = rbufer(pos)
        shared_pfl_point(i)%fts                    = rbufer(pos+1)
        shared_pfl_point(i)%columnar_content_He_m2 = rbufer(pos+2)
        shared_pfl_point(i)%columnar_content_O_m2  = rbufer(pos+3)
        shared_pfl_point(i)%columnar_content_N2_m2 = rbufer(pos+4)
        shared_pfl_point(i)%columnar_content_O2_m2 = rbufer(pos+5)
        shared_pfl_point(i)%Nn_He_m3               = rbufer(pos+6)
        shared_pfl_point(i)%Nn_O_m3                = rbufer(pos+7)
        shared_pfl_point(i)%Nn_N2_m3               = rbufer(pos+8)
        shared_pfl_point(i)%Nn_O2_m3               = rbufer(pos+9)
        shared_pfl_point(i)%geoB_T                 = rbufer(pos+10)
        shared_pfl_point(i)%Te_K                   = rbufer(pos+11)
        shared_pfl_point(i)%Ne_m3                  = rbufer(pos+12)
        pos = pos+13
     END DO
     DO nr = 1, N_ranges
        min_energy_eV_range(nr) = rbufer(pos) 
        max_energy_eV_range(nr) = rbufer(pos+1)
        min_pitch_angle_range(nr) = rbufer(pos+2)
        max_pitch_angle_range(nr) = rbufer(pos+3)
        pos = pos+4
     END DO

  END IF   !###   IF (Rank_of_process.EQ.0) THEN

! cleanup
  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT = ALLOC_ERR)

  IF (N_ranges.GE.1) THEN
! arrays for fluxes, if necessary
     DO i = 1, pfl_N_of_points
        ALLOCATE(shared_pfl_point(i)%Ge_m2s_range(1:N_ranges), STAT = ALLOC_ERR)
        ALLOCATE(shared_pfl_point(i)%Ge_par_m2s_range(1:N_ranges), STAT = ALLOC_ERR)

        shared_pfl_point(i)%Ge_m2s_range = 0.0_8
        shared_pfl_point(i)%Ge_par_m2s_range = 0.0_8

        ALLOCATE(shared_pfl_point(i)%source_EUV_m3s_range(1:N_ranges), STAT = ALLOC_ERR)
        ALLOCATE(shared_pfl_point(i)%source_photoe_m3s_range(1:N_ranges), STAT = ALLOC_ERR)
        ALLOCATE(shared_pfl_point(i)%source_neutrals_m3s_range(1:N_ranges), STAT = ALLOC_ERR)
        ALLOCATE(shared_pfl_point(i)%source_plasma_m3s_range(1:N_ranges), STAT = ALLOC_ERR)
        ALLOCATE(shared_pfl_point(i)%sink_neutrals_m3s_range(1:N_ranges), STAT = ALLOC_ERR)
        ALLOCATE(shared_pfl_point(i)%sink_plasma_m3s_range(1:N_ranges), STAT = ALLOC_ERR)

        shared_pfl_point(i)%source_EUV_m3s_range = 0.0_8
        shared_pfl_point(i)%source_photoe_m3s_range = 0.0_8
        shared_pfl_point(i)%source_neutrals_m3s_range = 0.0_8
        shared_pfl_point(i)%source_plasma_m3s_range = 0.0_8
        shared_pfl_point(i)%sink_neutrals_m3s_range = 0.0_8
        shared_pfl_point(i)%sink_plasma_m3s_range = 0.0_8

     END DO
  END IF

  ALLOCATE(evdf(-max_evdf_N_vpar:max_evdf_N_vpar, 0:max_evdf_N_vperp), STAT = ALLOC_ERR)

  evdf = 0.0_8

  particle_count = 0
  step_count = 0_8
  mirror_count = 0

! clear rates 
  DO i = 1, pfl_N_of_points
     shared_pfl_point(i)%rate_iHe_production_m3s = 0.0_8
     shared_pfl_point(i)%rate_iN_production_m3s = 0.0_8
     shared_pfl_point(i)%rate_iO_production_m3s = 0.0_8
     shared_pfl_point(i)%rate_iN2_production_m3s = 0.0_8
     shared_pfl_point(i)%rate_iO2_production_m3s = 0.0_8
     shared_pfl_point(i)%kinetic_ne_m3 = 0.0_8
     shared_pfl_point(i)%Ge_plus_m2s = 0.0_8
     shared_pfl_point(i)%Ge_minus_m2s = 0.0_8
  END DO

  DO i = 1, pfl_N_of_points

     if (Rank_of_process.eq.1) print '(2x,i4,2x,i6,2x,i6)', Rank_of_process, i, pfl_N_of_points

! special case:
! either the observation point is in the Earth's shadow
! or it is above the EUV ionization threshold altitude
     IF ( (shared_pfl_point(i)%columnar_content_He_m2.LT.0.0_8).OR. &
        & (shared_pfl_point(i)%columnar_content_O_m2.LT.0.0_8).OR. &
        & (shared_pfl_point(i)%columnar_content_N2_m2.LT.0.0_8).OR. &
        & (shared_pfl_point(i)%columnar_content_O2_m2.LT.0.0_8) ) CYCLE

     DO jsol = 1, N_solar_bins
        DO jchannel = 1, N_photoe_producing_channels
           CALL Emit_Photoelectrons_From_This_PFL_Point( i, &
                                                       & jsol, &
                                                       & jchannel, &
                                                       & particle_count, &
                                                       & step_count, &
                                                       & mirror_count )
        END DO
     END DO
  END DO

  print '("proc ",i4," did its line with ",i7," particles and ",i12," steps and ",i9," mirrors")', Rank_of_process, particle_count, step_count, mirror_count

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !##########################################################################################################

! assemble all calculated values in process with rank zero

  temp = 1.0_8 / DBLE(N_of_processes)    ! averaging factor

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
  buflen = (8 + N_ranges*8) * pfl_N_of_points
  ALLOCATE(rbufer(buflen), STAT = ALLOC_ERR)

  IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT=ALLOC_ERR)
  ALLOCATE(rbufer2(buflen), STAT = ALLOC_ERR)
  rbufer2 = 0.0_8

 ! save its own values 
  pos=1
  DO i = 1, pfl_N_of_points
     rbufer(pos)   = shared_pfl_point(i)%rate_iHe_production_m3s
     rbufer(pos+1) = shared_pfl_point(i)%rate_iN_production_m3s
     rbufer(pos+2) = shared_pfl_point(i)%rate_iO_production_m3s
     rbufer(pos+3) = shared_pfl_point(i)%rate_iN2_production_m3s
     rbufer(pos+4) = shared_pfl_point(i)%rate_iO2_production_m3s
     rbufer(pos+5) = shared_pfl_point(i)%kinetic_ne_m3
     rbufer(pos+6) = shared_pfl_point(i)%Ge_plus_m2s
     rbufer(pos+7) = shared_pfl_point(i)%Ge_minus_m2s

     pos = pos+7
     DO nr = 1, N_ranges
        rbufer(pos+nr) = shared_pfl_point(i)%Ge_m2s_range(nr)
     END DO

     pos = pos+N_ranges
     DO nr = 1, N_ranges
        rbufer(pos+nr) = shared_pfl_point(i)%Ge_par_m2s_range(nr)
     END DO

     pos = pos+N_ranges
     DO nr = 1, N_ranges
        rbufer(pos+nr) = shared_pfl_point(i)%source_EUV_m3s_range(nr)
     END DO

     pos = pos+N_ranges
     DO nr = 1, N_ranges
        rbufer(pos+nr) = shared_pfl_point(i)%source_photoe_m3s_range(nr)
     END DO

     pos = pos+N_ranges
     DO nr = 1, N_ranges
        rbufer(pos+nr) = shared_pfl_point(i)%source_neutrals_m3s_range(nr)
     END DO

     pos = pos+N_ranges
     DO nr = 1, N_ranges
        rbufer(pos+nr) = shared_pfl_point(i)%source_plasma_m3s_range(nr)
     END DO

     pos = pos+N_ranges
     DO nr = 1, N_ranges
        rbufer(pos+nr) = shared_pfl_point(i)%sink_neutrals_m3s_range(nr)
     END DO

     pos = pos+N_ranges
     DO nr = 1, N_ranges
        rbufer(pos+nr) = shared_pfl_point(i)%sink_plasma_m3s_range(nr)
     END DO

     pos = pos+1+N_ranges
  END DO

  CALL MPI_REDUCE(rbufer, rbufer2, buflen, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  IF (Rank_of_process.EQ.0) THEN
! translate the message
     pos=1
     DO i = 1, pfl_N_of_points
        shared_pfl_point(i)%rate_iHe_production_m3s = rbufer2(pos)   * temp
        shared_pfl_point(i)%rate_iN_production_m3s  = rbufer2(pos+1) * temp
        shared_pfl_point(i)%rate_iO_production_m3s  = rbufer2(pos+2) * temp
        shared_pfl_point(i)%rate_iN2_production_m3s = rbufer2(pos+3) * temp
        shared_pfl_point(i)%rate_iO2_production_m3s = rbufer2(pos+4) * temp
        shared_pfl_point(i)%kinetic_ne_m3           = rbufer2(pos+5) * temp
        shared_pfl_point(i)%Ge_plus_m2s             = rbufer2(pos+6) * temp
        shared_pfl_point(i)%Ge_minus_m2s            = rbufer2(pos+7) * temp

        pos = pos+7
        DO nr = 1, N_ranges
           shared_pfl_point(i)%Ge_m2s_range(nr) = rbufer2(pos+nr) * temp
        END DO

        pos = pos+N_ranges
        DO nr = 1, N_ranges
           shared_pfl_point(i)%Ge_par_m2s_range(nr) = rbufer2(pos+nr) * temp
        END DO

        pos = pos+N_ranges
        DO nr = 1, N_ranges
           shared_pfl_point(i)%source_EUV_m3s_range(nr) = rbufer2(pos+nr) * temp
        END DO

        pos = pos+N_ranges
        DO nr = 1, N_ranges
           shared_pfl_point(i)%source_photoe_m3s_range(nr) = rbufer2(pos+nr) * temp
        END DO

        pos = pos+N_ranges
        DO nr = 1, N_ranges
           shared_pfl_point(i)%source_neutrals_m3s_range(nr) = rbufer2(pos+nr) * temp
        END DO

        pos = pos+N_ranges
        DO nr = 1, N_ranges
           shared_pfl_point(i)%source_plasma_m3s_range(nr) = rbufer2(pos+nr) * temp
        END DO

        pos = pos+N_ranges
        DO nr = 1, N_ranges
           shared_pfl_point(i)%sink_neutrals_m3s_range(nr) = rbufer2(pos+nr) * temp
        END DO

        pos = pos+N_ranges
        DO nr = 1, N_ranges
           shared_pfl_point(i)%sink_plasma_m3s_range(nr) = rbufer2(pos+nr) * temp
        END DO

        pos = pos+1+N_ranges
     END DO
  END IF   !###   IF (Rank_of_process.EQ.0) THEN

! assemble evdf 

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
  buflen = (2 * max_evdf_N_vpar + 1) * (max_evdf_N_vperp + 1)
  ALLOCATE(rbufer(buflen), STAT = ALLOC_ERR)

  IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT=ALLOC_ERR)
  ALLOCATE(rbufer2(buflen), STAT = ALLOC_ERR)
  rbufer2 = 0.0_8

 ! save its own values
  pos=1
  DO j = 0, max_evdf_N_vperp
     DO i = -max_evdf_N_vpar, max_evdf_N_vpar
        rbufer(pos)   = evdf(i,j)
        pos=pos+1
     END DO
  END DO
  
  CALL MPI_REDUCE(rbufer, rbufer2, buflen, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  IF (Rank_of_process.EQ.0) THEN
! translate the message
     pos=1
     DO j = 0, max_evdf_N_vperp
        DO i = -max_evdf_N_vpar, max_evdf_N_vpar
           evdf(i,j) = rbufer2(pos) * temp
           pos = pos+1
        END DO
     END DO
  END IF
 
  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
  IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT=ALLOC_ERR)

  RETURN

END SUBROUTINE CALCULATE_PHOTOELECTRON_ENERGY_SPECTRA

!---------------------------------------------
! save data
!
SUBROUTINE SAVE_PHOTOELECTRON_DATA(i_orbit_point)

  USE photoelectrons
  USE PhysicalConstants
  USE field_line
  USE SpacecraftValues

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i_orbit_point

!                                  ! ----x----I----
  CHARACTER(14) pfl_filename       ! op_NNN_pfl.dat
                                   ! op_NNN_osp.dat   ! "osp" stands for "only spacecraft"

  REAL(8) temp_euv_rates_m3s(1:10) ! 1 = He+; 2-6 = O+(4S/2D/2P/4P/2P*) from O; 7,8 = N+,N2+ from N2; 9,10 = O+,O2+ from O2

  REAL(8) vec(2:1233)              ! don't expect more than 99 energy ranges but rhere are 12 parameters for each range to save in each point
  REAL(8) avg_factor_eV_sr         !
  INTEGER nr, i, shift, j
  INTEGER, PARAMETER :: Ntheta = 200
  REAL(8) dtheta, theta
  REAL(8) loc_min_angle, loc_max_angle

!                                  ! ----x----I----x
  CHARACTER(15) op_evdf_filename   ! op_NNN_evdf.dat

!                                  ! ----x----I----x---
  CHARACTER(18) pfl_vtk_filename   ! op_NNN_pfl_GEO.vtk
                                   ! op_NNN_ter_GEO.vtk    terminator

  INTEGER, PARAMETER :: N_ter_points = 361
  REAL(8), PARAMETER :: dphi_ter = 2.0_8 * pi / (N_ter_points-1)
  REAL(8), PARAMETER :: R_ter_km = R_Earth_km + 100.0_8

  REAL(8) y_ter_gse_km, z_ter_gse_km  ! x_ter_gse_km = 0.0_8
  REAL(8) x_ter_geo_km, y_ter_geo_km, z_ter_geo_km

  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string

     FUNCTION Get_EUV_production_rates_m3s( columnar_content_He_m2, &
                                          & columnar_content_O_m2, &
                                          & columnar_content_N2_m2, &
                                          & columnar_content_O2_m2, &
                                          & Nn_He_m3, Nn_O_m3, Nn_N2_m3, Nn_O2_m3)
       REAL(8) Get_EUV_production_rates_m3s(1:10)
       REAL(8) columnar_content_He_m2
       REAL(8) columnar_content_O_m2
       REAL(8) columnar_content_N2_m2
       REAL(8) columnar_content_O2_m2
       REAL(8) Nn_He_m3, Nn_O_m3, Nn_N2_m3, Nn_O2_m3
     END FUNCTION Get_EUV_production_rates_m3s
  end interface

! general field line characteristics ---------------------------------------------------------------------------------------------------------------------------------
! general field line characteristics ---------------------------------------------------------------------------------------------------------------------------------
! general field line characteristics ---------------------------------------------------------------------------------------------------------------------------------
! general field line characteristics ---------------------------------------------------------------------------------------------------------------------------------
! general field line characteristics ---------------------------------------------------------------------------------------------------------------------------------

!               ! ----x----I----
  pfl_filename = 'op_NNN_pfl.dat'
  pfl_filename(4:6)   = convert_int_to_txt_string(i_orbit_point, 3)

  open(19, file = pfl_filename)

! header
  write (19, '("# col  1 is point number")')

  write (19, '("# col  2 is distance along the field line (km)")')
  write (19, '("# col  3 is altitude (km)")')
  write (19, '("# col  4 is GEO latitude (deg)")')
  write (19, '("# col  5 is GEO longitude (deg)")')

  write (19, '("# col  6 is columnar content of He (m^2)")')
  write (19, '("# col  7 is columnar content of O (m^2)")')
  write (19, '("# col  8 is columnar content of N2 (m^2)")')
  write (19, '("# col  9 is columnar content of O2 (m^2)")')

  write (19, '("# col 10 is density of He (m^{-3})")')
  write (19, '("# col 11 is density of O (m^{-3})")')
  write (19, '("# col 12 is density of N2 (m^{-3})")')
  write (19, '("# col 13 is density of O2 (m^{-3})")')

  write (19, '("# col 14 is geomagnetic field magnitude (T)")')
  write (19, '("# col 15 is electron temperature (K)")')
  write (19, '("# col 16 is electron density (m^{-3})")')

  write (19, '("# col 17 is EUV ionization rate, He+ from He,  (m^{-3}s^{-1})")')
  write (19, '("# col 18 is EUV ionization rate, O+4S from O,  (m^{-3}s^{-1})")')
  write (19, '("# col 19 is EUV ionization rate, O+2D from O,  (m^{-3}s^{-1})")')
  write (19, '("# col 20 is EUV ionization rate, O+2P from O,  (m^{-3}s^{-1})")')
  write (19, '("# col 21 is EUV ionization rate, O+4P from O,  (m^{-3}s^{-1})")')
  write (19, '("# col 22 is EUV ionization rate, O+2P* from O, (m^{-3}s^{-1})")')
  write (19, '("# col 23 is EUV ionization rate, N+,           (m^{-3}s^{-1})")')
  write (19, '("# col 24 is EUV ionization rate, N2+,          (m^{-3}s^{-1})")')
  write (19, '("# col 25 is EUV ionization rate, O+ from O2,   (m^{-3}s^{-1})")')
  write (19, '("# col 26 is EUV ionization rate, O2+,          (m^{-3}s^{-1})")')

  write (19, '("# col 27 is photoelectron He+ ionization rate  (m^{-3}s^{-1})")')
  write (19, '("# col 28 is photoelectron N+  ionization rate  (m^{-3}s^{-1})")')
  write (19, '("# col 29 is photoelectron O+  ionization rate  (m^{-3}s^{-1})")')
  write (19, '("# col 30 is photoelectron N2+ ionization rate  (m^{-3}s^{-1})")')
  write (19, '("# col 31 is photoelectron O2+ ionization rate  (m^{-3}s^{-1})")')

  write (19, '("# col 32 is photoelectron number density (m^{-3})")')
  write (19, '("# col 33 is photoelectron particle flux in the positive direction, northward, (m^{-2}s^{-1})")')
  write (19, '("# col 34 is photoelectron particle flux in the negative direction, southward, (m^{-2}s^{-1})")')

  shift = 34
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is diff.photoe.part.flux (1/m2/s/eV/sr) averaged over range ",i2," :: energy ",f10.3," : ",f10.3," eV, pitch angle (in spacecraft location) ",f7.3," : ",f7.3," deg")') &
          & shift+nr, nr, min_energy_eV_range(nr), max_energy_eV_range(nr), min_pitch_angle_range(nr) * rad_to_deg, max_pitch_angle_range(nr) * rad_to_deg
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is the flux averaging factor (denominator) (eV sr) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is local minimal pitch angle (deg) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is local maximal pitch angle (deg) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is source due to ionization by EUV (1/m^3/s) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is source due to ionization by photoelectrons (1/m^3/s) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is source due to scattering by neutrals (1/m^3/s) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is source due to interaction with plasma (1/m^3/s) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is sink due to scattering by neutrals (1/m^3/s) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is sink due to interaction with plasma (1/m^3/s) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is parallel particle flux (1/m^2/s) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is divergence of the parallel particle flux (1/m^3/s) for range ",i2)') shift+nr, nr
  end do

  do i = 1, pfl_N_of_points

     temp_euv_rates_m3s = Get_EUV_production_rates_m3s( shared_pfl_point(i)%columnar_content_He_m2, &
                                                      & shared_pfl_point(i)%columnar_content_O_m2, &
                                                      & shared_pfl_point(i)%columnar_content_N2_m2, &
                                                      & shared_pfl_point(i)%columnar_content_O2_m2, &
                                                      & shared_pfl_point(i)%Nn_He_m3, &
                                                      & shared_pfl_point(i)%Nn_O_m3, &
                                                      & shared_pfl_point(i)%Nn_N2_m3, &
                                                      & shared_pfl_point(i)%Nn_O2_m3 )

     vec(2) = shared_pfl_point(i)%L_m * 0.001_8                     ! 2
     vec(3) = (pfl_point_coords(i)%r-1.0_8) * R_Earth_km            ! 3
     vec(4) = 90.0_8 - pfl_point_coords(i)%theta_geo * rad_to_deg   ! 4
     vec(5) = pfl_point_coords(i)%phi_geo * rad_to_deg              ! 5
!
     vec(6) = shared_pfl_point(i)%columnar_content_He_m2           ! 6
     vec(7) = shared_pfl_point(i)%columnar_content_O_m2            ! 7
     vec(8) = shared_pfl_point(i)%columnar_content_N2_m2           ! 8
     vec(9) = shared_pfl_point(i)%columnar_content_O2_m2           ! 9
!
     vec(10) = shared_pfl_point(i)%Nn_He_m3      ! 10
     vec(11) = shared_pfl_point(i)%Nn_O_m3       ! 11
     vec(12) = shared_pfl_point(i)%Nn_N2_m3      ! 12
     vec(13) = shared_pfl_point(i)%Nn_O2_m3      ! 13
!
     vec(14) = shared_pfl_point(i)%geoB_T        ! 14
     vec(15) = shared_pfl_point(i)%Te_K          ! 15
     vec(16) = shared_pfl_point(i)%Ne_m3         ! 16
!
     vec(17) = temp_euv_rates_m3s(1)         ! He+                  ! 17
     vec(18) = temp_euv_rates_m3s(2)         ! O+4S  [from O]       ! 18
     vec(19) = temp_euv_rates_m3s(3)         ! O+2D  [from O]       ! 19
     vec(20) = temp_euv_rates_m3s(4)         ! O+2P  [from O]       ! 20
     vec(21) = temp_euv_rates_m3s(5)         ! O+4P  [from O]       ! 21
     vec(22) = temp_euv_rates_m3s(6)         ! O+2P* [from O]       ! 22
     vec(23) = temp_euv_rates_m3s(7)         ! N+                   ! 23
     vec(24) = temp_euv_rates_m3s(8)         ! N2+                  ! 24
     vec(25) = temp_euv_rates_m3s(9)         ! O+ [from O2]         ! 25
     vec(26) = temp_euv_rates_m3s(10)        ! O2+                  ! 26
!
     vec(27) = shared_pfl_point(i)%rate_iHe_production_m3s      ! 27
     vec(28) = shared_pfl_point(i)%rate_iN_production_m3s       ! 28
     vec(29) = shared_pfl_point(i)%rate_iO_production_m3s       ! 29
     vec(30) = shared_pfl_point(i)%rate_iN2_production_m3s      ! 30
     vec(31) = shared_pfl_point(i)%rate_iO2_production_m3s      ! 31
!
     vec(32) = shared_pfl_point(i)%kinetic_ne_m3            ! 32
     vec(33) = shared_pfl_point(i)%Ge_plus_m2s              ! 33
     vec(34) = shared_pfl_point(i)%Ge_minus_m2s             ! 34

     shift=34
     DO nr = 1, N_ranges
        call get_local_pitch_angle_range( loc_min_angle, &
                                        & loc_max_angle, &
                                        & shared_pfl_point(i)%geoB_T, &
                                        & min_pitch_angle_range(nr), &
                                        & max_pitch_angle_range(nr), &
                                        & shared_pfl_point(pfl_i_spacecraft)%geoB_T )

        avg_factor_eV_sr = 0.0_8
        dtheta = (loc_max_angle - loc_min_angle) / Ntheta
        do j = 1, Ntheta
           theta = loc_min_angle + (dble(j) - 0.5_8) * dtheta
           avg_factor_eV_sr = avg_factor_eV_sr + sin(theta)
        end do
        avg_factor_eV_sr = 2.0_8 * pi * (max_energy_eV_range(nr) - min_energy_eV_range(nr)) * avg_factor_eV_sr * dtheta

        vec(shift+           nr) = shared_pfl_point(i)%Ge_m2s_range(nr) / avg_factor_eV_sr
        vec(shift+  N_ranges+nr) = avg_factor_eV_sr
        vec(shift+2*N_ranges+nr) = loc_min_angle * rad_to_deg
        vec(shift+3*N_ranges+nr) = loc_max_angle * rad_to_deg
     END DO

     shift = shift+4*N_ranges
     do nr = 1, N_ranges
        vec(shift+nr) = shared_pfl_point(i)%source_EUV_m3s_range(nr)
     end do

     shift = shift+N_ranges
     do nr = 1, N_ranges
        vec(shift+nr) = shared_pfl_point(i)%source_photoe_m3s_range(nr)
     end do

     shift = shift+N_ranges
     do nr = 1, N_ranges
        vec(shift+nr) = shared_pfl_point(i)%source_neutrals_m3s_range(nr)
     end do

     shift = shift+N_ranges
     do nr = 1, N_ranges
        vec(shift+nr) = shared_pfl_point(i)%source_plasma_m3s_range(nr)
     end do

     shift = shift+N_ranges
     do nr = 1, N_ranges
        vec(shift+nr) = shared_pfl_point(i)%sink_neutrals_m3s_range(nr)
     end do

     shift = shift+N_ranges
     do nr = 1, N_ranges
        vec(shift+nr) = shared_pfl_point(i)%sink_plasma_m3s_range(nr)
     end do

     shift = shift+N_ranges
     do nr = 1, N_ranges
        vec(shift+nr) = shared_pfl_point(i)%Ge_par_m2s_range(nr)
     end do

     shift = shift+N_ranges
     do nr = 1, N_ranges
        if (i.eq.1) then
           vec(shift+nr) = 0.0_8
        else if (i.eq.pfl_N_of_points) then
           vec(shift+nr) = 0.0_8
        else
           vec(shift+nr) = ( shared_pfl_point(i+1)%Ge_par_m2s_range(nr) * shared_pfl_point(i+1)%fts - &
                         &   shared_pfl_point(i-1)%Ge_par_m2s_range(nr) * shared_pfl_point(i-1)%fts ) &
                         & / &
                         & ((shared_pfl_point(i+1)%L_m - shared_pfl_point(i-1)%L_m) * shared_pfl_point(i)%fts)
        end if
     end do

     write (19, '(2x,i6,4(2x,f12.3),29(2x,e14.7),4x,1200(2x,e14.7))') i, vec(2:(shift+N_ranges)) ! vec(2:(34+12*N_ranges))

  end do   !###   do i = 1, pfl_N_of_points
  
  close (19, status = 'keep')
  print '("file ",A14," is ready")', pfl_filename

!------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
!
! save a file with the same data s above in the location of the spacecraft only (single point)
!
!               ! ----x----I----
  pfl_filename = 'op_NNN_osp.dat'
  pfl_filename(4:6)   = convert_int_to_txt_string(i_orbit_point, 3)

  open(19, file = pfl_filename)

  i = pfl_i_spacecraft

! header
  write (19, '("# col  1 is point number (*** this is the point where the spacecraft is ***)")')

  write (19, '("# col  2 is distance along the field line (km)")')
  write (19, '("# col  3 is altitude (km)")')
  write (19, '("# col  4 is GEO latitude (deg)")')
  write (19, '("# col  5 is GEO longitude (deg)")')

  write (19, '("# col  6 is columnar content of He (m^2)")')
  write (19, '("# col  7 is columnar content of O (m^2)")')
  write (19, '("# col  8 is columnar content of N2 (m^2)")')
  write (19, '("# col  9 is columnar content of O2 (m^2)")')

  write (19, '("# col 10 is density of He (m^{-3})")')
  write (19, '("# col 11 is density of O (m^{-3})")')
  write (19, '("# col 12 is density of N2 (m^{-3})")')
  write (19, '("# col 13 is density of O2 (m^{-3})")')

  write (19, '("# col 14 is geomagnetic field magnitude (T)")')
  write (19, '("# col 15 is electron temperature (K)")')
  write (19, '("# col 16 is electron density (m^{-3})")')

  write (19, '("# col 17 is EUV ionization rate, He+ from He,  (m^{-3}s^{-1})")')
  write (19, '("# col 18 is EUV ionization rate, O+4S from O,  (m^{-3}s^{-1})")')
  write (19, '("# col 19 is EUV ionization rate, O+2D from O,  (m^{-3}s^{-1})")')
  write (19, '("# col 20 is EUV ionization rate, O+2P from O,  (m^{-3}s^{-1})")')
  write (19, '("# col 21 is EUV ionization rate, O+4P from O,  (m^{-3}s^{-1})")')
  write (19, '("# col 22 is EUV ionization rate, O+2P* from O, (m^{-3}s^{-1})")')
  write (19, '("# col 23 is EUV ionization rate, N+,           (m^{-3}s^{-1})")')
  write (19, '("# col 24 is EUV ionization rate, N2+,          (m^{-3}s^{-1})")')
  write (19, '("# col 25 is EUV ionization rate, O+ from O2,   (m^{-3}s^{-1})")')
  write (19, '("# col 26 is EUV ionization rate, O2+,          (m^{-3}s^{-1})")')

  write (19, '("# col 27 is photoelectron He+ ionization rate  (m^{-3}s^{-1})")')
  write (19, '("# col 28 is photoelectron N+  ionization rate  (m^{-3}s^{-1})")')
  write (19, '("# col 29 is photoelectron O+  ionization rate  (m^{-3}s^{-1})")')
  write (19, '("# col 30 is photoelectron N2+ ionization rate  (m^{-3}s^{-1})")')
  write (19, '("# col 31 is photoelectron O2+ ionization rate  (m^{-3}s^{-1})")')

  write (19, '("# col 32 is photoelectron number density (m^{-3})")')
  write (19, '("# col 33 is photoelectron particle flux in the positive direction, northward, (m^{-2}s^{-1})")')
  write (19, '("# col 34 is photoelectron particle flux in the negative direction, southward, (m^{-2}s^{-1})")')

  shift = 34
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is diff.photoe.part.flux (1/m2/s/eV/sr) averaged over range ",i2," :: energy ",f10.3," : ",f10.3," eV, pitch angle (in spacecraft location) ",f7.3," : ",f7.3," deg")') &
          & shift+nr, nr, min_energy_eV_range(nr), max_energy_eV_range(nr), min_pitch_angle_range(nr) * rad_to_deg, max_pitch_angle_range(nr) * rad_to_deg
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is the flux averaging factor (denominator) (eV sr) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is local minimal pitch angle (deg) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is local maximal pitch angle (deg) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is source due to ionization by EUV (1/m^3/s) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is source due to ionization by photoelectrons (1/m^3/s) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is source due to scattering by neutrals (1/m^3/s) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is source due to interaction with plasma (1/m^3/s) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is sink due to scattering by neutrals (1/m^3/s) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is sink due to interaction with plasma (1/m^3/s) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is parallel particle flux (1/m^2/s) for range ",i2)') shift+nr, nr
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     write (19, '("# col ",i4," is divergence of the parallel particle flux (1/m^3/s) for range ",i2)') shift+nr, nr
  end do

  temp_euv_rates_m3s = Get_EUV_production_rates_m3s( shared_pfl_point(i)%columnar_content_He_m2, &
                                                   & shared_pfl_point(i)%columnar_content_O_m2, &
                                                   & shared_pfl_point(i)%columnar_content_N2_m2, &
                                                   & shared_pfl_point(i)%columnar_content_O2_m2, &
                                                   & shared_pfl_point(i)%Nn_He_m3, &
                                                   & shared_pfl_point(i)%Nn_O_m3, &
                                                   & shared_pfl_point(i)%Nn_N2_m3, &
                                                   & shared_pfl_point(i)%Nn_O2_m3 )

  vec(2) = shared_pfl_point(i)%L_m * 0.001_8                     ! 2
  vec(3) = (pfl_point_coords(i)%r-1.0_8) * R_Earth_km            ! 3
  vec(4) = 90.0_8 - pfl_point_coords(i)%theta_geo * rad_to_deg   ! 4
  vec(5) = pfl_point_coords(i)%phi_geo * rad_to_deg              ! 5
!
  vec(6) = shared_pfl_point(i)%columnar_content_He_m2           ! 6
  vec(7) = shared_pfl_point(i)%columnar_content_O_m2            ! 7
  vec(8) = shared_pfl_point(i)%columnar_content_N2_m2           ! 8
  vec(9) = shared_pfl_point(i)%columnar_content_O2_m2           ! 9
!
  vec(10) = shared_pfl_point(i)%Nn_He_m3      ! 10
  vec(11) = shared_pfl_point(i)%Nn_O_m3       ! 11
  vec(12) = shared_pfl_point(i)%Nn_N2_m3      ! 12
  vec(13) = shared_pfl_point(i)%Nn_O2_m3      ! 13
!
  vec(14) = shared_pfl_point(i)%geoB_T        ! 14
  vec(15) = shared_pfl_point(i)%Te_K          ! 15
  vec(16) = shared_pfl_point(i)%Ne_m3         ! 16
!
  vec(17) = temp_euv_rates_m3s(1)         ! He+                  ! 17
  vec(18) = temp_euv_rates_m3s(2)         ! O+4S  [from O]       ! 18
  vec(19) = temp_euv_rates_m3s(3)         ! O+2D  [from O]       ! 19
  vec(20) = temp_euv_rates_m3s(4)         ! O+2P  [from O]       ! 20
  vec(21) = temp_euv_rates_m3s(5)         ! O+4P  [from O]       ! 21
  vec(22) = temp_euv_rates_m3s(6)         ! O+2P* [from O]       ! 22
  vec(23) = temp_euv_rates_m3s(7)         ! N+                   ! 23
  vec(24) = temp_euv_rates_m3s(8)         ! N2+                  ! 24
  vec(25) = temp_euv_rates_m3s(9)         ! O+ [from O2]         ! 25
  vec(26) = temp_euv_rates_m3s(10)        ! O2+                  ! 26
!
  vec(27) = shared_pfl_point(i)%rate_iHe_production_m3s      ! 27
  vec(28) = shared_pfl_point(i)%rate_iN_production_m3s       ! 28
  vec(29) = shared_pfl_point(i)%rate_iO_production_m3s       ! 29
  vec(30) = shared_pfl_point(i)%rate_iN2_production_m3s      ! 30
  vec(31) = shared_pfl_point(i)%rate_iO2_production_m3s      ! 31
!
  vec(32) = shared_pfl_point(i)%kinetic_ne_m3            ! 32
  vec(33) = shared_pfl_point(i)%Ge_plus_m2s              ! 33
  vec(34) = shared_pfl_point(i)%Ge_minus_m2s             ! 34
  
  shift=34
  DO nr = 1, N_ranges
     call get_local_pitch_angle_range( loc_min_angle, &
                                     & loc_max_angle, &
                                     & shared_pfl_point(i)%geoB_T, &
                                     & min_pitch_angle_range(nr), &
                                     & max_pitch_angle_range(nr), &
                                     & shared_pfl_point(pfl_i_spacecraft)%geoB_T )

     avg_factor_eV_sr = 0.0_8
     dtheta = (loc_max_angle - loc_min_angle) / Ntheta
     do j = 1, Ntheta
        theta = loc_min_angle + (dble(j) - 0.5_8) * dtheta
        avg_factor_eV_sr = avg_factor_eV_sr + sin(theta)
     end do
     avg_factor_eV_sr = 2.0_8 * pi * (max_energy_eV_range(nr) - min_energy_eV_range(nr)) * avg_factor_eV_sr * dtheta

     vec(shift+           nr) = shared_pfl_point(i)%Ge_m2s_range(nr) / avg_factor_eV_sr
     vec(shift+  N_ranges+nr) = avg_factor_eV_sr
     vec(shift+2*N_ranges+nr) = loc_min_angle * rad_to_deg
     vec(shift+3*N_ranges+nr) = loc_max_angle * rad_to_deg
  END DO

  shift = shift+4*N_ranges
  do nr = 1, N_ranges
     vec(shift+nr) = shared_pfl_point(i)%source_EUV_m3s_range(nr)
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     vec(shift+nr) = shared_pfl_point(i)%source_photoe_m3s_range(nr)
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     vec(shift+nr) = shared_pfl_point(i)%source_neutrals_m3s_range(nr)
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     vec(shift+nr) = shared_pfl_point(i)%source_plasma_m3s_range(nr)
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     vec(shift+nr) = shared_pfl_point(i)%sink_neutrals_m3s_range(nr)
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     vec(shift+nr) = shared_pfl_point(i)%sink_plasma_m3s_range(nr)
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     vec(shift+nr) = shared_pfl_point(i)%Ge_par_m2s_range(nr)
  end do

  shift = shift+N_ranges
  do nr = 1, N_ranges
     if (i.eq.1) then
        vec(shift+nr) = 0.0_8
     else if (i.eq.pfl_N_of_points) then
        vec(shift+nr) = 0.0_8
     else
        vec(shift+nr) = ( shared_pfl_point(i+1)%Ge_par_m2s_range(nr) * shared_pfl_point(i+1)%fts - &
                      &   shared_pfl_point(i-1)%Ge_par_m2s_range(nr) * shared_pfl_point(i-1)%fts ) &
                      & / &
                      & ((shared_pfl_point(i+1)%L_m - shared_pfl_point(i-1)%L_m) * shared_pfl_point(i)%fts)
     end if
  end do

  write (19, '(2x,i6,4(2x,f12.3),29(2x,e14.7),4x,1200(2x,e14.7))') i, vec(2:(shift+N_ranges)) ! vec(2:(34+12*N_ranges))

  close (19, status = 'keep')
  print '("file ",A14," is ready")', pfl_filename

! evdf ----------------------------------------------------------------------------------------------------------------------------------------------------
! evdf ----------------------------------------------------------------------------------------------------------------------------------------------------
! evdf ----------------------------------------------------------------------------------------------------------------------------------------------------
! evdf ----------------------------------------------------------------------------------------------------------------------------------------------------
! evdf ----------------------------------------------------------------------------------------------------------------------------------------------------


!                   ! ----x----I----x
  op_evdf_filename = 'op_NNN_evdf.dat'
  op_evdf_filename(4:6)   = convert_int_to_txt_string(i_orbit_point, 3)

  open(19, file = op_evdf_filename)

  write (19, '("# max_evdf_N_vperp = ",i4)') max_evdf_N_vperp
  write (19, '("# max_evdf_N_vpar = ",i4)') max_evdf_N_vpar
  write (19, '("# dvpar_ms = ",f10.2)') dvpar_ms
  write (19, '("# dvperp_ms = ",f10.2)') dvperp_ms  

  write (19, '("# col 1 is the parallel (relative to the magnetic field) velocity bin number")')
  write (19, '("# col 2 is the transverse velocity bin number")')
  write (19, '("# col 3 is the middle of the parallel velocity bin (m/s)")')
  write (19, '("# col 4 is the middle of the transverse velocity bin (m/s)")')
  write (19, '("# col 5 is the photoelectron velocity distribution function ((1/m^3)/(m/s)^2)")')

  do j = 0, max_evdf_N_vperp
     do i = -max_evdf_N_vpar, max_evdf_N_vpar
        write (19, '(2x,i4,2x,i4,2x,f12.1,2x,f12.1,2x,e12.5)') &
             & i, &
             & j, &
             & i * dvpar_ms, &
             & (real(j) + 0.5) * dvperp_ms, &
             & real(evdf(i,j) / (dvpar_ms * dvperp_ms))
     end do
     write (19, '(" ")')
  end do

  close (19, status = 'keep')
  print '("file ",A15," is ready")', op_evdf_filename

! field line for paraview ---------------------------------------------------------------------------------------------------------------------------------
! field line for paraview ---------------------------------------------------------------------------------------------------------------------------------
! field line for paraview ---------------------------------------------------------------------------------------------------------------------------------
! field line for paraview ---------------------------------------------------------------------------------------------------------------------------------
! field line for paraview ---------------------------------------------------------------------------------------------------------------------------------

!                   ! ----x----I----x---
  pfl_vtk_filename = 'op_NNN_pfl_GEO.vtk'
  pfl_vtk_filename(4:6) = convert_int_to_txt_string(i_orbit_point, 3)

  open (19, file = pfl_vtk_filename)

  write (19, '("# vtk DataFile Version 4.2")')
  write (19, '(A14)') pfl_vtk_filename(1:14)
  write (19, '("ASCII")')
  write (19, '("DATASET POLYDATA")')
  write (19, '("POINTS ",i6," float")') pfl_N_of_points

  do i = 1, pfl_N_of_points
     write (19, '(3(2x,f10.2))') pfl_point_coords(i)%x_geo_km, pfl_point_coords(i)%y_geo_km, pfl_point_coords(i)%z_geo_km
  end do

  write (19, '("LINES ",i4,2x,i6)') 1, pfl_N_of_points+1
  write (19, '(2x,i6)') pfl_N_of_points
  do i = 0, pfl_N_of_points-1
     write (19, '(2x,i6)') i
  end do

  write (19, '("POINT_DATA ",i6)') pfl_N_of_points
  write (19, '("SCALARS totEUVm3s float")')
  write (19, '("LOOKUP_TABLE default")')

  do i = 1, pfl_N_of_points
     temp_euv_rates_m3s = Get_EUV_production_rates_m3s( shared_pfl_point(i)%columnar_content_He_m2, &
                                                      & shared_pfl_point(i)%columnar_content_O_m2, &
                                                      & shared_pfl_point(i)%columnar_content_N2_m2, &
                                                      & shared_pfl_point(i)%columnar_content_O2_m2, &
                                                      & shared_pfl_point(i)%Nn_He_m3, &
                                                      & shared_pfl_point(i)%Nn_O_m3, &
                                                      & shared_pfl_point(i)%Nn_N2_m3, &
                                                      & shared_pfl_point(i)%Nn_O2_m3 )
     write (19, '(2x,e12.5)') real( temp_euv_rates_m3s(1) + &
                                  & temp_euv_rates_m3s(2) + &
                                  & temp_euv_rates_m3s(3) + &
                                  & temp_euv_rates_m3s(4) + &
                                  & temp_euv_rates_m3s(5) + &
                                  & temp_euv_rates_m3s(6) + &
                                  & temp_euv_rates_m3s(7) + &
                                  & temp_euv_rates_m3s(8) + &
                                  & temp_euv_rates_m3s(9) + &
                                  & temp_euv_rates_m3s(10) )
  end do
  
  close (19, status = 'keep')
  print '("file ",A18," is ready")', pfl_vtk_filename
  
! terminator for paraview ---------------------------------------------------------------------------------------------------------------------------------
! terminator for paraview ---------------------------------------------------------------------------------------------------------------------------------
! terminator for paraview ---------------------------------------------------------------------------------------------------------------------------------
! terminator for paraview ---------------------------------------------------------------------------------------------------------------------------------
! terminator for paraview ---------------------------------------------------------------------------------------------------------------------------------

!                   ! ----x----I----x---
  pfl_vtk_filename = 'op_NNN_ter_GEO.vtk'
  pfl_vtk_filename(4:6) = convert_int_to_txt_string(i_orbit_point, 3)

  open (19, file = pfl_vtk_filename)

  write (19, '("# vtk DataFile Version 4.2")')
  write (19, '(A14)') pfl_vtk_filename(1:14)
  write (19, '("ASCII")')
  write (19, '("DATASET POLYDATA")')
  write (19, '("POINTS ",i6," float")') N_ter_points

!  x_ter_gse_km = 0.0_8
  do i = 1, N_ter_points

     y_ter_gse_km = R_ter_km * COS((i-1) * dphi_ter)
     z_ter_gse_km = R_ter_km * SIN((i-1) * dphi_ter)

     CALL Convert_GSE8_to_GEO8(0.0_8, y_ter_gse_km, z_ter_gse_km, x_ter_geo_km, y_ter_geo_km, z_ter_geo_km)

     write (19, '(3(2x,f10.2))') x_ter_geo_km, y_ter_geo_km, z_ter_geo_km
  end do

  write (19, '("LINES ",i4,2x,i6)') 1, N_ter_points+1
  write (19, '(2x,i6)') N_ter_points
  do i = 0, N_ter_points-1
     write (19, '(2x,i6)') i
  end do

  close (19, status = 'keep')
  print '("file ",A18," is ready")', pfl_vtk_filename
  
END SUBROUTINE SAVE_PHOTOELECTRON_DATA
