!-----------------------------------------------
! there is only one field line here
! only the evdf in spacecraft location is calculated
! all processes share the same fieldline to improve statistics and reduce noise
!
! intialization of the field line is performed by process with rank 0
! values are set using IRI and MSIS models
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

  INTEGER ibufer5(5)

  REAL(8), ALLOCATABLE :: rbufer(:)
  REAL(8), ALLOCATABLE :: rbufer2(:)
  INTEGER ALLOC_ERR

  INTEGER buflen
  INTEGER pos, i, j, nr, kk

  INTEGER particle_count
  INTEGER(8) step_count
  INTEGER mirror_count

  INTEGER jsol, jchannel

  REAL(8) temp

  IF (Rank_of_process.EQ.0) THEN
! distribute field line data to other processes

     ibufer5(1) = pfl_N_of_points
     IF (pfl_closed) THEN
        ibufer5(2) = 1
     ELSE
        ibufer5(2) = 0
     END IF
     ibufer5(3) = pfl_i_spacecraft
     ibufer5(4) = N_ranges
     ibufer5(5) = N_coll_ranges

     CALL MPI_BCAST(ibufer5, 5, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

     buflen = 11*pfl_N_of_points + N_ranges*4 + N_coll_ranges*2
     ALLOCATE(rbufer(1:buflen), STAT = ALLOC_ERR)
     pos=1
     DO i = 1, pfl_N_of_points
        rbufer(pos)    = shared_pfl_point(i)%L_m
        rbufer(pos+1)  = shared_pfl_point(i)%fts
        rbufer(pos+2)  = shared_pfl_point(i)%columnar_content_O_m2
        rbufer(pos+3)  = shared_pfl_point(i)%columnar_content_N2_m2
        rbufer(pos+4)  = shared_pfl_point(i)%columnar_content_O2_m2
        rbufer(pos+5)  = shared_pfl_point(i)%Nn_O_m3
        rbufer(pos+6)  = shared_pfl_point(i)%Nn_N2_m3
        rbufer(pos+7)  = shared_pfl_point(i)%Nn_O2_m3
        rbufer(pos+8)  = shared_pfl_point(i)%geoB_T
        rbufer(pos+9)  = shared_pfl_point(i)%Te_K
        rbufer(pos+10) = shared_pfl_point(i)%Ne_m3
        pos = pos+11
     END DO
     DO nr = 1, N_ranges
        rbufer(pos)    = min_energy_eV_range(nr)
        rbufer(pos+1)  = max_energy_eV_range(nr)
        rbufer(pos+2)  = min_pitch_angle_range(nr)
        rbufer(pos+3)  = max_pitch_angle_range(nr)
        pos = pos+4
     END DO
     DO nr = 1, N_coll_ranges
        rbufer(pos)    = min_w_eV_colrange(nr)
        rbufer(pos+1)  = max_w_eV_colrange(nr)
        pos = pos+2
     END DO
     CALL MPI_BCAST(rbufer, buflen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  ELSE     !###   IF (Rank_of_process.EQ.0) THEN
! receive the field line

     CALL MPI_BCAST(ibufer5, 5, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

     pfl_N_of_points = ibufer5(1)
     pfl_closed = .TRUE.
     IF (ibufer5(2).EQ.0) pfl_closed = .FALSE.
     pfl_i_spacecraft = ibufer5(3)
     N_ranges = ibufer5(4)
     N_coll_ranges = ibufer5(5)

! allocate field line
     ALLOCATE(shared_pfl_point(1:pfl_N_of_points), STAT = ALLOC_ERR)
     buflen = 11*pfl_N_of_points + N_ranges*4 + N_coll_ranges*2
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
        
     IF (ALLOCATED(min_w_eV_colrange)) DEALLOCATE(min_w_eV_colrange)
     IF (ALLOCATED(max_w_eV_colrange)) DEALLOCATE(max_w_eV_colrange)

     if (N_coll_ranges.ge.1) then
        ALLOCATE(min_w_eV_colrange(1:N_coll_ranges), STAT = ALLOC_ERR)
        ALLOCATE(max_w_eV_colrange(1:N_coll_ranges), STAT = ALLOC_ERR)
     end if

     CALL MPI_BCAST(rbufer, buflen, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
! translate the message
     pos=1
     DO i = 1, pfl_N_of_points
        shared_pfl_point(i)%L_m                    = rbufer(pos)
        shared_pfl_point(i)%fts                    = rbufer(pos+1)
        shared_pfl_point(i)%columnar_content_O_m2  = rbufer(pos+2)
        shared_pfl_point(i)%columnar_content_N2_m2 = rbufer(pos+3)
        shared_pfl_point(i)%columnar_content_O2_m2 = rbufer(pos+4)
        shared_pfl_point(i)%Nn_O_m3                = rbufer(pos+5)
        shared_pfl_point(i)%Nn_N2_m3               = rbufer(pos+6)
        shared_pfl_point(i)%Nn_O2_m3               = rbufer(pos+7)
        shared_pfl_point(i)%geoB_T                 = rbufer(pos+8)
        shared_pfl_point(i)%Te_K                   = rbufer(pos+9)
        shared_pfl_point(i)%Ne_m3                  = rbufer(pos+10)
        pos = pos+11
     END DO
     DO nr = 1, N_ranges
        min_energy_eV_range(nr) = rbufer(pos) 
        max_energy_eV_range(nr) = rbufer(pos+1)
        min_pitch_angle_range(nr) = rbufer(pos+2)
        max_pitch_angle_range(nr) = rbufer(pos+3)
        pos = pos+4
     END DO
     DO nr = 1, N_coll_ranges
        min_w_eV_colrange(nr) = rbufer(pos)
        max_w_eV_colrange(nr) = rbufer(pos+1)
        pos = pos+2
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

        ALLOCATE(shared_pfl_point(i)%source_EUV_m3s_range_channel(   1:N_ranges,1:N_photoe_producing_channels), STAT = ALLOC_ERR)
        ALLOCATE(shared_pfl_point(i)%source_photoe_m3s_range_channel(1:N_ranges,1:5), STAT = ALLOC_ERR)

        ALLOCATE(shared_pfl_point(i)%source_neutrals_m3s_range_coll(1:N_ranges,1:N_colkind), STAT = ALLOC_ERR)
        ALLOCATE(shared_pfl_point(i)%sink_neutrals_m3s_range_coll(  1:N_ranges,1:N_colkind), STAT = ALLOC_ERR)

        ALLOCATE(shared_pfl_point(i)%source_plasma_m3s_range(1:N_ranges), STAT = ALLOC_ERR)
        ALLOCATE(shared_pfl_point(i)%sink_plasma_m3s_range(  1:N_ranges), STAT = ALLOC_ERR)

        ALLOCATE(shared_pfl_point(i)%source_Coulomb_m3s_range(1:N_ranges), STAT = ALLOC_ERR)
        ALLOCATE(shared_pfl_point(i)%sink_Coulomb_m3s_range(  1:N_ranges), STAT = ALLOC_ERR)

        shared_pfl_point(i)%source_EUV_m3s_range_channel    = 0.0_8
        shared_pfl_point(i)%source_photoe_m3s_range_channel = 0.0_8

        shared_pfl_point(i)%source_neutrals_m3s_range_coll = 0.0_8
        shared_pfl_point(i)%sink_neutrals_m3s_range_coll   = 0.0_8

        shared_pfl_point(i)%source_plasma_m3s_range = 0.0_8
        shared_pfl_point(i)%sink_plasma_m3s_range   = 0.0_8

        shared_pfl_point(i)%source_Coulomb_m3s_range = 0.0_8
        shared_pfl_point(i)%sink_Coulomb_m3s_range   = 0.0_8

     END DO
  END IF

  if (N_coll_ranges.gt.0) then
     allocate(coll_freq_s1_L_colkind_colrange(1:pfl_N_of_points, 1:N_colkind, 1:N_coll_ranges), stat = alloc_err)
     allocate(          density_m3_L_colrange(1:pfl_N_of_points,              1:N_coll_ranges), stat = alloc_err)

     coll_freq_s1_L_colkind_colrange = 0.0_8
     density_m3_L_colrange = 0.0_8
  end if

  ALLOCATE(evdf(-max_evdf_N_vpar:max_evdf_N_vpar, 0:max_evdf_N_vperp), STAT = ALLOC_ERR)

  evdf = 0.0_8

  particle_count = 0
  step_count = 0_8
  mirror_count = 0

! clear rates 
  DO i = 1, pfl_N_of_points
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
     IF ( (shared_pfl_point(i)%columnar_content_O_m2.LT.0.0_8).OR. &
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

  print '("proc ",i4," did its line with ",i9," particles and ",i12," steps and ",i9," mirrors")', Rank_of_process, particle_count, step_count, mirror_count

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr) !##########################################################################################################

! assemble all calculated values in process with rank zero

  temp = 1.0_8 / DBLE(N_of_processes)    ! averaging factor

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)

  buflen = (7 + N_ranges * (6 + N_photoe_producing_channels + 5 + N_colkind + N_colkind)) * pfl_N_of_points

  ALLOCATE(rbufer(buflen), STAT = ALLOC_ERR)

  IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT=ALLOC_ERR)
  ALLOCATE(rbufer2(buflen), STAT = ALLOC_ERR)
  rbufer2 = 0.0_8

 ! save its own values 
  pos=1
  DO i = 1, pfl_N_of_points

     rbufer(pos)   = shared_pfl_point(i)%rate_iN_production_m3s
     rbufer(pos+1) = shared_pfl_point(i)%rate_iO_production_m3s
     rbufer(pos+2) = shared_pfl_point(i)%rate_iN2_production_m3s
     rbufer(pos+3) = shared_pfl_point(i)%rate_iO2_production_m3s
     rbufer(pos+4) = shared_pfl_point(i)%kinetic_ne_m3
     rbufer(pos+5) = shared_pfl_point(i)%Ge_plus_m2s
     rbufer(pos+6) = shared_pfl_point(i)%Ge_minus_m2s

     pos = pos+6
     DO nr = 1, N_ranges
        rbufer(pos+nr) = shared_pfl_point(i)%Ge_m2s_range(nr)
     END DO

     pos = pos+N_ranges
     DO nr = 1, N_ranges
        rbufer(pos+nr) = shared_pfl_point(i)%Ge_par_m2s_range(nr)
     END DO

     pos = pos+N_ranges
     DO nr = 1, N_ranges
        rbufer(pos+nr) = shared_pfl_point(i)%source_plasma_m3s_range(nr)
     END DO

     pos = pos+N_ranges
     DO nr = 1, N_ranges
        rbufer(pos+nr) = shared_pfl_point(i)%sink_plasma_m3s_range(nr)
     END DO

     pos = pos+N_ranges
     DO nr = 1, N_ranges
        rbufer(pos+nr) = shared_pfl_point(i)%source_Coulomb_m3s_range(nr)
     END DO

     pos = pos+N_ranges
     DO nr = 1, N_ranges
        rbufer(pos+nr) = shared_pfl_point(i)%sink_Coulomb_m3s_range(nr)
     END DO

     pos = pos+N_ranges

     DO nr = 1, N_ranges
        do j = 1, N_photoe_producing_channels
           pos = pos+1
           rbufer(pos) = shared_pfl_point(i)%source_EUV_m3s_range_channel(nr,j)
        end do
     END DO

     DO nr = 1, N_ranges
        do j = 1, 5
           pos = pos+1
           rbufer(pos) = shared_pfl_point(i)%source_photoe_m3s_range_channel(nr,j)
        end do
     END DO

     DO nr = 1, N_ranges
        do j = 1, N_colkind
           pos = pos+1
           rbufer(pos) = shared_pfl_point(i)%source_neutrals_m3s_range_coll(nr,j)
        end do
     END DO

     DO nr = 1, N_ranges
        do j = 1, N_colkind
           pos = pos+1
           rbufer(pos) = shared_pfl_point(i)%sink_neutrals_m3s_range_coll(nr,j)
        end do
     END DO

     pos = pos+1

  END DO

  CALL MPI_REDUCE(rbufer, rbufer2, buflen, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  IF (Rank_of_process.EQ.0) THEN
! translate the message
     pos=1
     DO i = 1, pfl_N_of_points

        shared_pfl_point(i)%rate_iN_production_m3s  = rbufer2(pos) * temp
        shared_pfl_point(i)%rate_iO_production_m3s  = rbufer2(pos+1) * temp
        shared_pfl_point(i)%rate_iN2_production_m3s = rbufer2(pos+2) * temp
        shared_pfl_point(i)%rate_iO2_production_m3s = rbufer2(pos+3) * temp
        shared_pfl_point(i)%kinetic_ne_m3           = rbufer2(pos+4) * temp
        shared_pfl_point(i)%Ge_plus_m2s             = rbufer2(pos+5) * temp
        shared_pfl_point(i)%Ge_minus_m2s            = rbufer2(pos+6) * temp

        pos = pos+6
        DO nr = 1, N_ranges
           shared_pfl_point(i)%Ge_m2s_range(nr) = rbufer2(pos+nr) * temp
        END DO

        pos = pos+N_ranges
        DO nr = 1, N_ranges
           shared_pfl_point(i)%Ge_par_m2s_range(nr) = rbufer2(pos+nr) * temp
        END DO

        pos = pos+N_ranges
        DO nr = 1, N_ranges
           shared_pfl_point(i)%source_plasma_m3s_range(nr) = rbufer2(pos+nr) * temp
        END DO

        pos = pos+N_ranges
        DO nr = 1, N_ranges
           shared_pfl_point(i)%sink_plasma_m3s_range(nr) = rbufer2(pos+nr) * temp
        END DO

        pos = pos+N_ranges
        DO nr = 1, N_ranges
           shared_pfl_point(i)%source_Coulomb_m3s_range(nr) = rbufer2(pos+nr) * temp
        END DO

        pos = pos+N_ranges
        DO nr = 1, N_ranges
           shared_pfl_point(i)%sink_Coulomb_m3s_range(nr) = rbufer2(pos+nr) * temp
        END DO

        pos = pos+N_ranges

        DO nr = 1, N_ranges
           do j = 1, N_photoe_producing_channels
              pos = pos+1
              shared_pfl_point(i)%source_EUV_m3s_range_channel(nr,j) = rbufer2(pos) * temp
           end do
        END DO

        DO nr = 1, N_ranges
           do j = 1, 5
              pos = pos+1
              shared_pfl_point(i)%source_photoe_m3s_range_channel(nr,j) = rbufer2(pos) * temp
           end do
        END DO

        DO nr = 1, N_ranges
           do j = 1, N_colkind
              pos = pos+1
              shared_pfl_point(i)%source_neutrals_m3s_range_coll(nr,j) = rbufer2(pos) * temp
           end do
        END DO

        DO nr = 1, N_ranges
           do j = 1, N_colkind
              pos = pos+1
              shared_pfl_point(i)%sink_neutrals_m3s_range_coll(nr,j) = rbufer2(pos) * temp
           end do
        END DO

        pos = pos+1
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

! assemble collision frequencies and corresponding partial densities

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
  IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT=ALLOC_ERR)

  buflen = pfl_N_of_points * N_coll_ranges * (N_colkind + 1)  ! +1 includes partial densities for each range

  IF (buflen.EQ.0) RETURN

  ALLOCATE(rbufer(buflen), STAT = ALLOC_ERR)
  ALLOCATE(rbufer2(buflen), STAT = ALLOC_ERR)
  rbufer2 = 0.0_8

 ! save its own values 
  pos=1
  do nr = 1, N_coll_ranges
     do kk = 1, N_colkind
        do i = 1, pfl_N_of_points
           rbufer(pos) = coll_freq_s1_L_colkind_colrange(i, kk, nr)
           pos = pos+1
        end do
     end do
  end do
  do nr = 1, N_coll_ranges
     do i = 1, pfl_N_of_points
        rbufer(pos) = density_m3_L_colrange(i, nr)
        pos = pos+1
     end do
  end do

  CALL MPI_REDUCE(rbufer, rbufer2, buflen, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  IF (Rank_of_process.EQ.0) THEN
! translate the message
     pos=1
     do nr = 1, N_coll_ranges
        do kk = 1, N_colkind
           do i = 1, pfl_N_of_points
              coll_freq_s1_L_colkind_colrange(i, kk, nr) = rbufer2(pos)
              pos = pos+1
           end do
        end do
     end do
     do nr = 1, N_coll_ranges
        do i = 1, pfl_N_of_points
           density_m3_L_colrange(i, nr) = rbufer2(pos)
           pos = pos+1
        end do
     end do

     do nr = 1, N_coll_ranges
        do kk = 1, N_colkind
           do i = 1, pfl_N_of_points
              if (density_m3_L_colrange(i, nr).gt.0.0_8) then
                 coll_freq_s1_L_colkind_colrange(i, kk, nr) = coll_freq_s1_L_colkind_colrange(i, kk, nr) / density_m3_L_colrange(i, nr)
              else
                 coll_freq_s1_L_colkind_colrange(i, kk, nr) = 0.0_8
              end if
           end do
        end do
     end do
     do nr = 1, N_coll_ranges
        do i = 1, pfl_N_of_points
           density_m3_L_colrange(i, nr) = density_m3_L_colrange(i, nr) * temp
        end do
     end do
  end IF
 
! final cleanup

  IF (ALLOCATED(rbufer)) DEALLOCATE(rbufer, STAT=ALLOC_ERR)
  IF (ALLOCATED(rbufer2)) DEALLOCATE(rbufer2, STAT=ALLOC_ERR)

  RETURN

END SUBROUTINE CALCULATE_PHOTOELECTRON_ENERGY_SPECTRA

!---------------------------------------------
! save data
!
SUBROUTINE SAVE_PHOTOELECTRON_DATA(i_orbit_point)

  USE Photoelectrons
  USE PhysicalConstants
  USE field_line
  USE SpacecraftValues

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i_orbit_point

!                                  ! ----x----I----
  CHARACTER(14) pfl_filename       ! op_NNN_pfl.dat
                                   ! op_NNN_osp.dat   ! "osp" stands for "only spacecraft"

!                                  ! ----x----I----x----I---
  CHARACTER(23) range_filename     ! op_NNN_range_NN_pfl.dat
                                   ! op_NNN_range_NN_osp.dat   ! "osp" stands for "only spacecraft"

  REAL(8) temp_euv_rates_m3s(1:9) ! 1-5 = O+(4S/2D/2P/4P/2P*) from O; 6,7 = N+,N2+ from N2; 8,9 = O+,O2+ from O2

  REAL(8) vec(2:137)               ! for 52 e-n collisions, 5 e-n ionization and 9 EUV photoionization processes
  real(8) vec_osp(2:137)

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

! collision frequencies

  real(8) w_eV
  character(38) collfreq_filename   ! op_NNN_collfreq_s1_avgw_NNN_eV_pfl.dat
!                                   ! ----x----I----x----I----x----I----x---
  real(8) dw_eV
  integer, parameter :: n_w_points_colrange = 20
  real(8) Nn_O_m3, Nn_N2_m3, Nn_O2_m3
  real(8) avgcoll_s1(1:N_colkind)
  integer n, kk
  real(8) sumcoll_s1

  real(8) avgcoll_s1_osp(1:N_colkind)

  character(47) collfreqmeas_filename    ! op_NNN_collfreq_s1_avgw_NNN_eV_measured_pfl.dat
!                                        ! ----x----I----x----I----x----I----x----I----x--

! functions
  real(8) fl_Frequency_of_Collision_s1
  real(8) fl_Max_frequency_of_en_collisions_s1

  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string

     FUNCTION Get_EUV_production_rates_m3s( columnar_content_O_m2, &
                                          & columnar_content_N2_m2, &
                                          & columnar_content_O2_m2, &
                                          & Nn_O_m3, Nn_N2_m3, Nn_O2_m3)
       REAL(8) Get_EUV_production_rates_m3s(1:9)
       REAL(8) columnar_content_O_m2
       REAL(8) columnar_content_N2_m2
       REAL(8) columnar_content_O2_m2
       REAL(8) Nn_O_m3, Nn_N2_m3, Nn_O2_m3
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

  write (19, '("# col  6 is columnar content of O (m^2)")')
  write (19, '("# col  7 is columnar content of N2 (m^2)")')
  write (19, '("# col  8 is columnar content of O2 (m^2)")')

  write (19, '("# col  9 is density of O (m^{-3})")')
  write (19, '("# col 10 is density of N2 (m^{-3})")')
  write (19, '("# col 11 is density of O2 (m^{-3})")')

  write (19, '("# col 12 is geomagnetic field magnitude (T)")')
  write (19, '("# col 13 is electron temperature (K)")')
  write (19, '("# col 14 is electron density (m^{-3})")')

  write (19, '("# col 15 is EUV ionization rate, O+4S from O,  (m^{-3}s^{-1})")')
  write (19, '("# col 16 is EUV ionization rate, O+2D from O,  (m^{-3}s^{-1})")')
  write (19, '("# col 17 is EUV ionization rate, O+2P from O,  (m^{-3}s^{-1})")')
  write (19, '("# col 18 is EUV ionization rate, O+4P* from O, (m^{-3}s^{-1})")')
  write (19, '("# col 19 is EUV ionization rate, O+2P* from O, (m^{-3}s^{-1})")')
  write (19, '("# col 20 is EUV ionization rate, N+,           (m^{-3}s^{-1})")')
  write (19, '("# col 21 is EUV ionization rate, N2+,          (m^{-3}s^{-1})")')
  write (19, '("# col 22 is EUV ionization rate, O+ from O2,   (m^{-3}s^{-1})")')
  write (19, '("# col 23 is EUV ionization rate, O2+,          (m^{-3}s^{-1})")')

  write (19, '("# col 24 is photoelectron N+  ionization rate  (m^{-3}s^{-1})")')
  write (19, '("# col 25 is photoelectron O+  ionization rate  (m^{-3}s^{-1})")')
  write (19, '("# col 26 is photoelectron N2+ ionization rate  (m^{-3}s^{-1})")')
  write (19, '("# col 27 is photoelectron O2+ ionization rate  (m^{-3}s^{-1})")')

  write (19, '("# col 28 is photoelectron number density (m^{-3})")')
  write (19, '("# col 29 is photoelectron particle flux in the positive direction, northward, (m^{-2}s^{-1})")')
  write (19, '("# col 30 is photoelectron particle flux in the negative direction, southward, (m^{-2}s^{-1})")')

  do i = 1, pfl_N_of_points

     vec = 0.0_8

     temp_euv_rates_m3s = Get_EUV_production_rates_m3s( shared_pfl_point(i)%columnar_content_O_m2, &
                                                      & shared_pfl_point(i)%columnar_content_N2_m2, &
                                                      & shared_pfl_point(i)%columnar_content_O2_m2, &
                                                      & shared_pfl_point(i)%Nn_O_m3, &
                                                      & shared_pfl_point(i)%Nn_N2_m3, &
                                                      & shared_pfl_point(i)%Nn_O2_m3 )

     vec(2) = shared_pfl_point(i)%L_m * 0.001_8                     ! 2
     vec(3) = (pfl_point_coords(i)%r-1.0_8) * R_Earth_km            ! 3
     vec(4) = 90.0_8 - pfl_point_coords(i)%theta_geo * rad_to_deg   ! 4
     vec(5) = pfl_point_coords(i)%phi_geo * rad_to_deg              ! 5
!
     vec(6) = shared_pfl_point(i)%columnar_content_O_m2            ! 6
     vec(7) = shared_pfl_point(i)%columnar_content_N2_m2           ! 7
     vec(8) = shared_pfl_point(i)%columnar_content_O2_m2           ! 8
!
     vec( 9) = shared_pfl_point(i)%Nn_O_m3       !  9
     vec(10) = shared_pfl_point(i)%Nn_N2_m3      ! 10
     vec(11) = shared_pfl_point(i)%Nn_O2_m3      ! 11
!
     vec(12) = shared_pfl_point(i)%geoB_T        ! 12
     vec(13) = shared_pfl_point(i)%Te_K          ! 13
     vec(14) = shared_pfl_point(i)%Ne_m3         ! 14
!
     vec(15) = temp_euv_rates_m3s(1)         ! O+4S  [from O]       ! 15
     vec(16) = temp_euv_rates_m3s(2)         ! O+2D  [from O]       ! 16
     vec(17) = temp_euv_rates_m3s(3)         ! O+2P  [from O]       ! 17
     vec(18) = temp_euv_rates_m3s(4)         ! O+4P  [from O]       ! 18
     vec(19) = temp_euv_rates_m3s(5)         ! O+2P* [from O]       ! 19
     vec(20) = temp_euv_rates_m3s(6)         ! N+                   ! 20
     vec(21) = temp_euv_rates_m3s(7)         ! N2+                  ! 21
     vec(22) = temp_euv_rates_m3s(8)         ! O+ [from O2]         ! 22
     vec(23) = temp_euv_rates_m3s(9)         ! O2+                  ! 23
!
     vec(24) = shared_pfl_point(i)%rate_iN_production_m3s       ! 24
     vec(25) = shared_pfl_point(i)%rate_iO_production_m3s       ! 25
     vec(26) = shared_pfl_point(i)%rate_iN2_production_m3s      ! 26
     vec(27) = shared_pfl_point(i)%rate_iO2_production_m3s      ! 27
!
     vec(28) = shared_pfl_point(i)%kinetic_ne_m3            ! 28
     vec(29) = shared_pfl_point(i)%Ge_plus_m2s              ! 29
     vec(30) = shared_pfl_point(i)%Ge_minus_m2s             ! 30

     write (19, '(2x,i6,4(2x,f12.3),25(2x,e14.7))') i, vec(2:30)

     if (i.eq.pfl_i_spacecraft) then
        vec_osp = vec
     end if

  end do   !###   do i = 1, pfl_N_of_points
  
  close (19, status = 'keep')
  print '("file ",A14," is ready")', pfl_filename

!------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------
!
! save a file with the same data as above in the location of the spacecraft only (single point)
!
!               ! ----x----I----
  pfl_filename = 'op_NNN_osp.dat'
  pfl_filename(4:6)   = convert_int_to_txt_string(i_orbit_point, 3)

  open(19, file = pfl_filename)

! skip the header

  write (19, '(2x,i6,4(2x,f12.3),25(2x,e14.7))') pfl_i_spacecraft, vec_osp(2:30)

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
     temp_euv_rates_m3s = Get_EUV_production_rates_m3s( shared_pfl_point(i)%columnar_content_O_m2, &
                                                      & shared_pfl_point(i)%columnar_content_N2_m2, &
                                                      & shared_pfl_point(i)%columnar_content_O2_m2, &
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
                                  & temp_euv_rates_m3s(9) )
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

! ranges ---------------------------------------------------------------------------------------------------------------------------------------------------------
! ranges ---------------------------------------------------------------------------------------------------------------------------------------------------------
! ranges ---------------------------------------------------------------------------------------------------------------------------------------------------------
! ranges ---------------------------------------------------------------------------------------------------------------------------------------------------------
! ranges ---------------------------------------------------------------------------------------------------------------------------------------------------------
! ranges ---------------------------------------------------------------------------------------------------------------------------------------------------------

  do nr = 1, N_ranges

!                    ! ----x----I----x----I---
     range_filename = 'op_NNN_range_NN_pfl.dat'
     range_filename(4:6)   = convert_int_to_txt_string(i_orbit_point, 3)
     range_filename(14:15) = convert_int_to_txt_string(nr, 2)

     open(19, file = range_filename)

! header

     write (19, '("# range ",i2," :: energy ",f10.3," : ",f10.3," eV, pitch angle (in spacecraft location) ",f7.3," : ",f7.3," deg")') &
          & nr, min_energy_eV_range(nr), max_energy_eV_range(nr), min_pitch_angle_range(nr) * rad_to_deg, max_pitch_angle_range(nr) * rad_to_deg

     write (19, '("# col   1 is point number")')
     
     write (19, '("# col   2 is distance along the field line (km)")')
     write (19, '("# col   3 is altitude (km)")')
     write (19, '("# col   4 is GEO latitude (deg)")')
     write (19, '("# col   5 is GEO longitude (deg)")')

     write (19, '("# col   6 is differential photoelectron particle flux (1/m2/s/eV/sr) averaged over the range")') 
     write (19, '("# col   7 is flux averaging factor (denominator) (eV sr)")')
     write (19, '("# col   8 is local minimal pitch angle (deg)")')
     write (19, '("# col   9 is local maximal pitch angle (deg)")')

     write (19, '("# col  10 is the parallel photoelectron particle flux (1/m^2/s)")')

     write (19, '("# col  11 is divergence of the parallel particle flux (1/m^3/s)")')

     write (19, '("#")')
     write (19, '("# all sources and sinks below are in (1/m^3/s)")')
     write (19, '("#")')

     write (19, '("# col  12 is total source due to ionization by EUV for all ionization channels")')
     write (19, '("# col  13 is total source due to ionization by photoelectrons for all ionization channels")')

     write (19, '("# col  14 is total source due to scattering by neutrals for all collisions")')
     write (19, '("# col  15 is source due to interaction with plasma")')
     write (19, '("# col  16 is source due to Coulomb scattering")')
     write (19, '("#")')

     write (19, '("# col  17 is total sink due to scattering by neutrals for all collisions")')
     write (19, '("# col  18 is sink due to interaction with plasma")')
     write (19, '("# col  19 is sink due to Coulomb scattering")')

     write (19, '("#")')
     write (19, '("# SOURCES due to EUV ionization :")')
     write (19, '("#")')

     write (19, '("# col  20 is O+(4S)  from O")')
     write (19, '("# col  21 is O+(2D)  from O")')
     write (19, '("# col  22 is O+(2P)  from O")')
     write (19, '("# col  23 is O+(4P*) from O")')
     write (19, '("# col  24 is O+(2P*) from O")')
     write (19, '("# col  25 is N2+(combined) from N2")')
     write (19, '("# col  26 is N+(combined)  from N2")')
     write (19, '("# col  27 is O2+(combined) from O2")')
     write (19, '("# col  28 is O+(combined)  from O2")')
     write (19, '("#")')

     write (19, '("#")')
     write (19, '("# SOURCES due to ionization by photoelectrons :")')
     write (19, '("#")')

     write (19, '("# col  29 is N2+ from N2s")')
     write (19, '("# col  30 is N+  from N2")')
     write (19, '("# col  31 is O2+ from O2")')
     write (19, '("# col  32 is O+  from O2")')
     write (19, '("# col  33 is O+  from O")')

     write (19, '("#")')
     write (19, '("# SOURCES due to electron-neutral N2 collisions :")')
     write (19, '("#")')

     write (19, '("# col  34 is collision #  1, elastic")')
     write (19, '("# col  35 is collision #  2, ionization N2+")')
     write (19, '("# col  36 is collision #  3, ionization N+")')
     write (19, '("# col  37 is collisions # 4-8, vibrational total")')
     write (19, '("# col  38 is collision #  9, A3Zup")')
     write (19, '("# col  39 is collision # 10, B3Pg")')
     write (19, '("# col  40 is collision # 11, C3Pg")')
     write (19, '("# col  41 is collision # 12, W3Du")')
     write (19, '("# col  42 is collision # 13, Bp3Zum")')
     write (19, '("# col  43 is collision # 14, a1Pg")')
     write (19, '("# col  44 is collision # 15, b1Pu")')
     write (19, '("# col  45 is collision # 16, bp1Zup")')
     write (19, '("# col  46 is collision # 17, c4p1Zu")')
     write (19, '("# col  47 is collision # 18, w1Du")')
     write (19, '("# col  48 is collision # 19, c1Pu")')
     write (19, '("# col  49 is collision # 20, ap1Zum")')
     write (19, '("# col  50 is collision # 21, app1Zgp")')
     write (19, '("# col  51 is collision # 22, other1Pu")')
     write (19, '("# col  52 is collision # 23, hls_15p8eV")')
     write (19, '("# col  53 is collision # 24, hls_VUV")')
     write (19, '("# col  54 is collision # 25, hls_17p3eV")')
     write (19, '("# col  55 is collision # 26, hls_Rydberg")')
     write (19, '("# col  56 is collision # 27, hls_tripman")')

     write (19, '("#")')
     write (19, '("# SOURCES due to electron-neutral O2 collisions :")')
     write (19, '("#")')

     write (19, '("# col  57 is collision # 28, elastic")')
     write (19, '("# col  58 is collision # 29, ionization O2+")')
     write (19, '("# col  59 is collision # 30, ionization O+")')
     write (19, '("# col  60 is collision # 31, vibrational")')
     write (19, '("# col  61 is collision # 32, Rydberg")')
     write (19, '("# col  62 is collision # 33, AApc")')
     write (19, '("# col  63 is collision # 34, a1Du")')
     write (19, '("# col  64 is collision # 35, b1Zgp")')
     write (19, '("# col  65 is collision # 36, longband")')
     write (19, '("# col  66 is collision # 37, secondband")')
     write (19, '("# col  67 is collision # 38, 13Pg")')
     write (19, '("# col  68 is collision # 39, 8p9eV")')
     write (19, '("# col  69 is collision # 40, B3Zu")')

     write (19, '("#")')
     write (19, '("# SOURCES due to electron-neutral O collisions :")')
     write (19, '("#")')

     write (19, '("# col  70 is collision # 41, elastic")')
     write (19, '("# col  71 is collision # 42-45, ionization 4S0/2D0/2P0/4P, total")')
     write (19, '("# col  72 is collision # 46, Rydberg")')
     write (19, '("# col  73 is collision # 47, 1D")')
     write (19, '("# col  74 is collision # 48, 1S")')
     write (19, '("# col  75 is collision # 49, 5P")')
     write (19, '("# col  76 is collision # 50, 5S")')
     write (19, '("# col  77 is collision # 51, 3s3S")')
     write (19, '("# col  78 is collision # 52, 3d3D")')
     write (19, '("# col  79 is collision # 53, 3sp3D")')
     write (19, '("# col  80 is collision # 54, 3p3P")')
     write (19, '("# col  81 is collision # 55, 5d3D")')
     write (19, '("# col  82 is collision # 56, 4d3D")')
     write (19, '("# col  83 is collision # 57, 2p53P")')
     write (19, '("# col  84 is collision # 58, 3spp3P")')
     write (19, '("# col  85 is collision # 59, 4dp3P")')

     write (19, '("#")')
     write (19, '("# LOSSES due to electron-neutral N2 collisions :")')
     write (19, '("#")')

     write (19, '("# col  86 is collision #  1, elastic")')
     write (19, '("# col  87 is collision #  2, ionization N2+")')
     write (19, '("# col  88 is collision #  3, ionization N+")')
     write (19, '("# col  89 is collisions # 4-8, vibrational total")')
     write (19, '("# col  90 is collision #  9, A3Zup")')
     write (19, '("# col  91 is collision # 10, B3Pg")')
     write (19, '("# col  92 is collision # 11, C3Pg")')
     write (19, '("# col  93 is collision # 12, W3Du")')
     write (19, '("# col  94 is collision # 13, Bp3Zum")')
     write (19, '("# col  95 is collision # 14, a1Pg")')
     write (19, '("# col  96 is collision # 15, b1Pu")')
     write (19, '("# col  97 is collision # 16, bp1Zup")')
     write (19, '("# col  98 is collision # 17, c4p1Zu")')
     write (19, '("# col  99 is collision # 18, w1Du")')
     write (19, '("# col 100 is collision # 19, c1Pu")')
     write (19, '("# col 101 is collision # 20, ap1Zum")')
     write (19, '("# col 102 is collision # 21, app1Zgp")')
     write (19, '("# col 103 is collision # 22, other1Pu")')
     write (19, '("# col 104 is collision # 23, hls_15p8eV")')
     write (19, '("# col 105 is collision # 24, hls_VUV")')
     write (19, '("# col 106 is collision # 25, hls_17p3eV")')
     write (19, '("# col 107 is collision # 26, hls_Rydberg")')
     write (19, '("# col 108 is collision # 27, hls_tripman")')

     write (19, '("#")')
     write (19, '("# LOSSES due to electron-neutral O2 collisions :")')
     write (19, '("#")')

     write (19, '("# col 109 is collision # 28, elastic")')
     write (19, '("# col 110 is collision # 29, ionization O2+")')
     write (19, '("# col 111 is collision # 30, ionization O+")')
     write (19, '("# col 112 is collision # 31, vibrational")')
     write (19, '("# col 113 is collision # 32, Rydberg")')
     write (19, '("# col 114 is collision # 33, AApc")')
     write (19, '("# col 115 is collision # 34, a1Du")')
     write (19, '("# col 116 is collision # 35, b1Zgp")')
     write (19, '("# col 117 is collision # 36, longband")')
     write (19, '("# col 118 is collision # 37, secondband")')
     write (19, '("# col 119 is collision # 38, 13Pg")')
     write (19, '("# col 120 is collision # 39, 8p9eV")')
     write (19, '("# col 121 is collision # 40, B3Zu")')

     write (19, '("#")')
     write (19, '("# LOSSES due to electron-neutral O collisions :")')
     write (19, '("#")')

     write (19, '("# col 122 is collision # 41, elastic")')
     write (19, '("# col 123 is collisions # 42-45, ionization 4S0/2D0/2P0/4P, total")')
     write (19, '("# col 124 is collision # 46, Rydberg")')
     write (19, '("# col 125 is collision # 47, 1D")')
     write (19, '("# col 126 is collision # 48, 1S")')
     write (19, '("# col 127 is collision # 49, 5P")')
     write (19, '("# col 128 is collision # 50, 5S")')
     write (19, '("# col 129 is collision # 51, 3s3S")')
     write (19, '("# col 130 is collision # 52, 3d3D")')
     write (19, '("# col 131 is collision # 53, 3sp3D")')
     write (19, '("# col 132 is collision # 54, 3p3P")')
     write (19, '("# col 133 is collision # 55, 5d3D")')
     write (19, '("# col 134 is collision # 56, 4d3D")')
     write (19, '("# col 135 is collision # 57, 2p53P")')
     write (19, '("# col 136 is collision # 58, 3spp3P")')
     write (19, '("# col 137 is collision # 59, 4dp3P")')

     do i = 1, pfl_N_of_points

        vec(2) = shared_pfl_point(i)%L_m * 0.001_8                     ! 2
        vec(3) = (pfl_point_coords(i)%r-1.0_8) * R_Earth_km            ! 3
        vec(4) = 90.0_8 - pfl_point_coords(i)%theta_geo * rad_to_deg   ! 4
        vec(5) = pfl_point_coords(i)%phi_geo * rad_to_deg              ! 5

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

        vec(6) = shared_pfl_point(i)%Ge_m2s_range(nr) / avg_factor_eV_sr
        vec(7) = avg_factor_eV_sr
        vec(8) = loc_min_angle * rad_to_deg
        vec(9) = loc_max_angle * rad_to_deg

        vec(10) = shared_pfl_point(i)%Ge_par_m2s_range(nr)

        if (i.eq.1) then
           vec(11) = 0.0_8
        else if (i.eq.pfl_N_of_points) then
           vec(11) = 0.0_8
        else
           vec(11) = ( shared_pfl_point(i+1)%Ge_par_m2s_range(nr) * shared_pfl_point(i+1)%fts - &
                   &   shared_pfl_point(i-1)%Ge_par_m2s_range(nr) * shared_pfl_point(i-1)%fts ) &
                   & / &
                   & ((shared_pfl_point(i+1)%L_m - shared_pfl_point(i-1)%L_m) * shared_pfl_point(i)%fts)
        end if

        vec(12) = 0.0_8
        do j = 1, N_photoe_producing_channels
           vec(12) = vec(12) + shared_pfl_point(i)%source_EUV_m3s_range_channel(nr,j)
        end do

        vec(13) = 0.0_8
        do j = 1, 5
           vec(13) = vec(13) + shared_pfl_point(i)%source_photoe_m3s_range_channel(nr,j)
        end do

        vec(14) = 0.0_8
        do j = 1, N_colkind
           vec(14) = vec(14) + shared_pfl_point(i)%source_neutrals_m3s_range_coll(nr,j)
        end do

        vec(15) = shared_pfl_point(i)%source_plasma_m3s_range(nr)
        vec(16) = shared_pfl_point(i)%source_Coulomb_m3s_range(nr)

        vec(17) = 0.0_8
        do j = 1, N_colkind
           vec(17) = vec(17) + shared_pfl_point(i)%sink_neutrals_m3s_range_coll(nr,j)
        end do

        vec(18) = shared_pfl_point(i)%sink_plasma_m3s_range(nr)
        vec(19) = shared_pfl_point(i)%sink_Coulomb_m3s_range(nr)

        do j = 1, N_photoe_producing_channels
           vec(19+j) = shared_pfl_point(i)%source_EUV_m3s_range_channel(nr,j)     ! 20 to 28
        end do

        do j = 1, 5
           vec(28+j) = shared_pfl_point(i)%source_photoe_m3s_range_channel(nr,j)  ! 29 to 33
        end do

        do j = 1, 3
           vec(33+j) = shared_pfl_point(i)%source_neutrals_m3s_range_coll(nr,j)   ! 34 to 36
        end do
        vec(37) = 0.0_8
        do j = 4, 8
           vec(37) = vec(37) + shared_pfl_point(i)%source_neutrals_m3s_range_coll(nr,j)   ! 37, total vibrational N2
        end do
        do j = 9, 41
           vec(33+j-4) = shared_pfl_point(i)%source_neutrals_m3s_range_coll(nr,j)   ! 38 to 70
        end do
        vec(71) = 0.0_8
        do j = 42, 45
           vec(71) = vec(71) + shared_pfl_point(i)%source_neutrals_m3s_range_coll(nr,j)   ! 71, total ionization O+
        end do
        do j = 46, N_colkind
           vec(33+j-4-3) = shared_pfl_point(i)%source_neutrals_m3s_range_coll(nr,j)   ! 72 to 85
        end do


        do j = 1, 3
           vec(85+j) = shared_pfl_point(i)%sink_neutrals_m3s_range_coll(nr,j)   ! 86 to 88
        end do
        vec(89) = 0.0_8
        do j = 4, 8
           vec(89) = vec(89) + shared_pfl_point(i)%sink_neutrals_m3s_range_coll(nr,j)   ! 89, total vibrational N2
        end do
        do j = 9, 41 !N_colkind
           vec(85+j-4) = shared_pfl_point(i)%sink_neutrals_m3s_range_coll(nr,j)   ! 90 to 122
        end do
        vec(123) = 0.0_8
        do j = 42, 45
           vec(123) = vec(123) + shared_pfl_point(i)%sink_neutrals_m3s_range_coll(nr,j)  ! 123, total ionization O+
        end do
        do j = 46, N_colkind
           vec(85+j-4-3) = shared_pfl_point(i)%sink_neutrals_m3s_range_coll(nr,j)   ! 124 to 137
        end do

        write (19, '(2x,i6,4(2x,f12.3),132(2x,e16.9))') i, vec(2:137)

        if (i.eq.pfl_i_spacecraft) then
           vec_osp = vec
        end if

     end do   !###   do i = 1, pfl_N_of_points

     close (19, status = 'keep')
     print '("file ",A23," is ready")', range_filename

! values in spacecraft location only ----------------------------------------------------
! values in spacecraft location only ----------------------------------------------------
! values in spacecraft location only ----------------------------------------------------
! values in spacecraft location only ----------------------------------------------------
! values in spacecraft location only ----------------------------------------------------
     
!                    ! ----x----I----x----I---
     range_filename = 'op_NNN_range_NN_osp.dat'
     range_filename(4:6)   = convert_int_to_txt_string(i_orbit_point, 3)
     range_filename(14:15) = convert_int_to_txt_string(nr, 2)

     open(19, file = range_filename)

! skip the header

     write (19, '(2x,i6,4(2x,f12.3),132(2x,e16.9))') pfl_i_spacecraft, vec_osp(2:137)

     close (19, status = 'keep')
     print '("file ",A23," is ready")', range_filename

  end do   !###   do nr = 1, N_ranges

! theoretical collision frequencies -------------------------------------------------------------------------------------------------------------------------------------------
! theoretical collision frequencies -------------------------------------------------------------------------------------------------------------------------------------------
! theoretical collision frequencies -------------------------------------------------------------------------------------------------------------------------------------------
! theoretical collision frequencies -------------------------------------------------------------------------------------------------------------------------------------------
! theoretical collision frequencies -------------------------------------------------------------------------------------------------------------------------------------------
! theoretical collision frequencies -------------------------------------------------------------------------------------------------------------------------------------------

  do nr = 1, N_coll_ranges

     w_eV = 0.5_8 * (min_w_eV_colrange(nr) + max_w_eV_colrange(nr))

!                       ! ----x----I----x----I----x----I----x---
     collfreq_filename = 'op_NNN_collfreq_s1_avgw_NNN_eV_pfl.dat'
     collfreq_filename(4:6)   = convert_int_to_txt_string(i_orbit_point, 3)
     collfreq_filename(25:27) = convert_int_to_txt_string(int(w_eV), 3)

     open(19, file = collfreq_filename)

     write (19, '("# theoretical collision frequencies (s^-1) averaged over energy interval ",i2," from ",f8.3," eV to ",f8.3," eV using ",i3," uniformly distributed energy values")') nr, min_w_eV_colrange(nr), max_w_eV_colrange(nr), n_w_points_colrange

     write (19, '("# electron-neutral N2 collisions --------------------------------------------- :")')

     write (19, '("# col  1, elastic")')
     write (19, '("# col  2, ionization N2+")')
     write (19, '("# col  3, ionization N+")')
     write (19, '("# col  4, vibrational 0-1")')
     write (19, '("# col  5, vibrational 0-2")')
     write (19, '("# col  6, vibrational 0-3")')
     write (19, '("# col  7, vibrational 0-4")')
     write (19, '("# col  8, vibrational 0-5")')
     write (19, '("# col  9, A3Zup")')
     write (19, '("# col 10, B3Pg")')
     write (19, '("# col 11, C3Pg")')
     write (19, '("# col 12, W3Du")')
     write (19, '("# col 13, Bp3Zum")')
     write (19, '("# col 14, a1Pg")')
     write (19, '("# col 15, b1Pu")')
     write (19, '("# col 16, bp1Zup")')
     write (19, '("# col 17, c4p1Zu")')
     write (19, '("# col 18, w1Du")')
     write (19, '("# col 19, c1Pu")')
     write (19, '("# col 20, ap1Zum")')
     write (19, '("# col 21, app1Zgp")')
     write (19, '("# col 22, other1Pu")')
     write (19, '("# col 23, hls_15p8eV")')
     write (19, '("# col 24, hls_VUV")')
     write (19, '("# col 25, hls_17p3eV")')
     write (19, '("# col 26, hls_Rydberg")')
     write (19, '("# col 27, hls_tripman")')

     write (19, '("# electron-neutral O2 collisions --------------------------------------------- :")')

     write (19, '("# col 28, elastic")')
     write (19, '("# col 29, ionization O2+")')
     write (19, '("# col 30, ionization O+")')
     write (19, '("# col 31, vibrational")')
     write (19, '("# col 32, Rydberg")')
     write (19, '("# col 33, AApc")')
     write (19, '("# col 34, a1Du")')
     write (19, '("# col 35, b1Zgp")')
     write (19, '("# col 36, longband")')
     write (19, '("# col 37, secondband")')
     write (19, '("# col 38, 13Pg")')
     write (19, '("# col 39, 8p9eV")')
     write (19, '("# col 40, B3Zu")')

     write (19, '("# electron-neutral O collisions ---------------------------------------------- :")')

     write (19, '("# col 41, elastic")')
     write (19, '("# col 42, ionization 4S0")')
     write (19, '("# col 43, ionization 2D0")')
     write (19, '("# col 44, ionization 2P0")')
     write (19, '("# col 45, ionization 4P")')
     write (19, '("# col 46, Rydberg")')
     write (19, '("# col 47, 1D")')
     write (19, '("# col 48, 1S")')
     write (19, '("# col 49, 5P")')
     write (19, '("# col 50, 5S")')
     write (19, '("# col 51, 3s3S")')
     write (19, '("# col 52, 3d3D")')
     write (19, '("# col 53, 3sp3D")')
     write (19, '("# col 54, 3p3P")')
     write (19, '("# col 55, 5d3D")')
     write (19, '("# col 56, 4d3D")')
     write (19, '("# col 57, 2p53P")')
     write (19, '("# col 58, 3spp3P")')
     write (19, '("# col 59, 4dp3P")')

     write (19, '("# ---------------------------------------------- ")')

     write (19, '("# col 60 is point number")')
     write (19, '("# col 61 is distance along the field line (km)")')
     write (19, '("# col 62 is altitude (km)")')
     write (19, '("# col 63 is the middle of the energy range (eV)")')
     write (19, '("# col 64 is the SUM of columns 1-59 (s^-1)")')
     write (19, '("# col 65 is the maximal frequency of e-n collisions (s^-1)")')

     dw_eV = (max_w_eV_colrange(nr) - min_w_eV_colrange(nr)) / n_w_points_colrange

     do i = 1, pfl_N_of_points
        Nn_O_m3  = shared_pfl_point(i)%Nn_O_m3
        Nn_N2_m3 = shared_pfl_point(i)%Nn_N2_m3
        Nn_O2_m3 = shared_pfl_point(i)%Nn_O2_m3

        avgcoll_s1 = 0.0_8
        do n = 1, n_w_points_colrange
           w_eV = min_w_eV_colrange(nr) + (dble(n) - 0.5_8) * dw_eV
           do kk = 1, N_colkind
              avgcoll_s1(kk) = avgcoll_s1(kk) + fl_Frequency_of_Collision_s1(kk, w_eV, Nn_O_m3, Nn_N2_m3, Nn_O2_m3)
           end do
        end do
        avgcoll_s1 = avgcoll_s1 / n_w_points_colrange

        if (i.eq.pfl_i_spacecraft) then
           avgcoll_s1_osp = avgcoll_s1
        end if

        sumcoll_s1 = 0.0_8
        do kk = 1, N_colkind
           sumcoll_s1 = sumcoll_s1 + avgcoll_s1(kk)
        end do

        w_eV = 0.5_8 * (min_w_eV_colrange(nr) + max_w_eV_colrange(nr))

        write (19, '(59(2x,e14.7),2x,i4,2(2x,f12.3),2x,f6.1,2x,2(2x,e14.7))') &
             & avgcoll_s1(1:N_colkind), &   !1-59
             & i, &                                   ! 60
             & shared_pfl_point(i)%L_m * 0.001_8, &           ! 61
             & (pfl_point_coords(i)%r-1.0_8) * R_Earth_km, &  ! 62
             & w_eV, &       ! 63
             & sumcoll_s1, &                                                                     ! 64
             & fl_Max_frequency_of_en_collisions_s1(w_eV, Nn_O_m3, Nn_N2_m3, Nn_O2_m3) ! 65

     end do   !###   do i = 1, pfl_N_of_points

     close (19, status = 'keep')
     print '("file ",A38," is ready")', collfreq_filename   

! values in spacecraft location only ----------------------------------------------------
! values in spacecraft location only ----------------------------------------------------
! values in spacecraft location only ----------------------------------------------------

     w_eV = 0.5_8 * (min_w_eV_colrange(nr) + max_w_eV_colrange(nr))

!                       ! ----x----I----x----I----x----I----x---
     collfreq_filename = 'op_NNN_collfreq_s1_avgw_NNN_eV_osp.dat'
     collfreq_filename(4:6)   = convert_int_to_txt_string(i_orbit_point, 3)
     collfreq_filename(25:27) = convert_int_to_txt_string(int(w_eV), 3)

     i = pfl_i_spacecraft

     Nn_O_m3  = shared_pfl_point(i)%Nn_O_m3
     Nn_N2_m3 = shared_pfl_point(i)%Nn_N2_m3
     Nn_O2_m3 = shared_pfl_point(i)%Nn_O2_m3

     sumcoll_s1 = 0.0_8
     do kk = 1, N_colkind
        sumcoll_s1 = sumcoll_s1 + avgcoll_s1_osp(kk)
     end do

     open(19, file = collfreq_filename)
     write (19, '("# theoretical collision frequencies averaged over energy interval ",i2," from ",f8.3," eV to ",f8.3," eV using ",i3," uniformly distributed energy values")') nr, min_w_eV_colrange(nr), max_w_eV_colrange(nr), n_w_points_colrange
     write (19, '(59(2x,e14.7),2x,i4,2(2x,f12.3),2x,f6.1,2x,2(2x,e14.7))') &
          & avgcoll_s1_osp(1:N_colkind), &   !1-59
          & i, &                                   ! 60
          & shared_pfl_point(i)%L_m * 0.001_8, &           ! 61
          & (pfl_point_coords(i)%r-1.0_8) * R_Earth_km, &  ! 62
          & w_eV, &       ! 63
          & sumcoll_s1, &                                                                     ! 64
          & fl_Max_frequency_of_en_collisions_s1(w_eV, Nn_O_m3, Nn_N2_m3, Nn_O2_m3) ! 65

     close (19, status = 'keep')
     print '("file ",A38," is ready")', collfreq_filename   

  end do   !###  do nr = 1, N_coll_ranges 

! calculated average collision frequencies ---------------------------------------------------------------------------------------------------------------------------------
! calculated average collision frequencies ---------------------------------------------------------------------------------------------------------------------------------
! calculated average collision frequencies ---------------------------------------------------------------------------------------------------------------------------------
! calculated average collision frequencies ---------------------------------------------------------------------------------------------------------------------------------
! calculated average collision frequencies ---------------------------------------------------------------------------------------------------------------------------------
! calculated average collision frequencies ---------------------------------------------------------------------------------------------------------------------------------

  do nr = 1, N_coll_ranges

     w_eV = 0.5_8 * (min_w_eV_colrange(nr) + max_w_eV_colrange(nr))

!                           ! ----x----I----x----I----x----I----x----I----x--
     collfreqmeas_filename = 'op_NNN_collfreq_s1_avgw_NNN_eV_measured_pfl.dat'
     collfreqmeas_filename(4:6)   = convert_int_to_txt_string(i_orbit_point, 3)
     collfreqmeas_filename(25:27) = convert_int_to_txt_string(int(w_eV), 3)

     open(19, file = collfreqmeas_filename)

     write (19, '("# collision frequencies averaged over energy interval ",i2," from ",f8.3," eV to ",f8.3," eV, calculated in simulation")') nr, min_w_eV_colrange(nr), max_w_eV_colrange(nr)     

     write (19, '("# electron-neutral N2 collisions --------------------------------------------- :")')

     write (19, '("# col  1, elastic")')
     write (19, '("# col  2, ionization N2+")')
     write (19, '("# col  3, ionization N+")')
     write (19, '("# col  4, vibrational 0-1")')
     write (19, '("# col  5, vibrational 0-2")')
     write (19, '("# col  6, vibrational 0-3")')
     write (19, '("# col  7, vibrational 0-4")')
     write (19, '("# col  8, vibrational 0-5")')
     write (19, '("# col  9, A3Zup")')
     write (19, '("# col 10, B3Pg")')
     write (19, '("# col 11, C3Pg")')
     write (19, '("# col 12, W3Du")')
     write (19, '("# col 13, Bp3Zum")')
     write (19, '("# col 14, a1Pg")')
     write (19, '("# col 15, b1Pu")')
     write (19, '("# col 16, bp1Zup")')
     write (19, '("# col 17, c4p1Zu")')
     write (19, '("# col 18, w1Du")')
     write (19, '("# col 19, c1Pu")')
     write (19, '("# col 20, ap1Zum")')
     write (19, '("# col 21, app1Zgp")')
     write (19, '("# col 22, other1Pu")')
     write (19, '("# col 23, hls_15p8eV")')
     write (19, '("# col 24, hls_VUV")')
     write (19, '("# col 25, hls_17p3eV")')
     write (19, '("# col 26, hls_Rydberg")')
     write (19, '("# col 27, hls_tripman")')

     write (19, '("# electron-neutral O2 collisions --------------------------------------------- :")')

     write (19, '("# col 28, elastic")')
     write (19, '("# col 29, ionization O2+")')
     write (19, '("# col 30, ionization O+")')
     write (19, '("# col 31, vibrational")')
     write (19, '("# col 32, Rydberg")')
     write (19, '("# col 33, AApc")')
     write (19, '("# col 34, a1Du")')
     write (19, '("# col 35, b1Zgp")')
     write (19, '("# col 36, longband")')
     write (19, '("# col 37, secondband")')
     write (19, '("# col 38, 13Pg")')
     write (19, '("# col 39, 8p9eV")')
     write (19, '("# col 40, B3Zu")')

     write (19, '("# electron-neutral O collisions ---------------------------------------------- :")')

     write (19, '("# col 41, elastic")')
     write (19, '("# col 42, ionization 4S0")')
     write (19, '("# col 43, ionization 2D0")')
     write (19, '("# col 44, ionization 2P0")')
     write (19, '("# col 45, ionization 4P")')
     write (19, '("# col 46, Rydberg")')
     write (19, '("# col 47, 1D")')
     write (19, '("# col 48, 1S")')
     write (19, '("# col 49, 5P")')
     write (19, '("# col 50, 5S")')
     write (19, '("# col 51, 3s3S")')
     write (19, '("# col 52, 3d3D")')
     write (19, '("# col 53, 3sp3D")')
     write (19, '("# col 54, 3p3P")')
     write (19, '("# col 55, 5d3D")')
     write (19, '("# col 56, 4d3D")')
     write (19, '("# col 57, 2p53P")')
     write (19, '("# col 58, 3spp3P")')
     write (19, '("# col 59, 4dp3P")')

     write (19, '("# ---------------------------------------------- ")')

     write (19, '("# col 60 is point number")')
     write (19, '("# col 61 is distance along the field line (km)")')
     write (19, '("# col 62 is altitude (km)")')
     write (19, '("# col 63 is the middle of the energy range (eV)")')
     write (19, '("# col 64 is the SUM of columns 1-59 (s^-1)")')
     write (19, '("# col 65 is the density of photoelectrons with energy within the aforementioned energy arnge (m^-3)")')

     do i = 1, pfl_N_of_points

        sumcoll_s1 = 0.0_8
        do kk = 1, N_colkind
           sumcoll_s1 = sumcoll_s1 + coll_freq_s1_L_colkind_colrange(i, kk, nr)
        end do

        write (19, '(59(2x,e14.7),2x,i4,2(2x,f12.3),2x,f6.1,2x,2(2x,e14.7))') &
             & coll_freq_s1_L_colkind_colrange(i, 1:N_colkind, nr), &   !1-59
             & i, &                                   ! 60
             & shared_pfl_point(i)%L_m * 0.001_8, &           ! 61
             & (pfl_point_coords(i)%r-1.0_8) * R_Earth_km, &  ! 62
             & w_eV, &       ! 63
             & sumcoll_s1, &    ! 64
             & density_m3_L_colrange(i,nr)  ! 65
        
     end do   !###   do i = 1, pfl_N_of_points

     close (19, status = 'keep')
     print '("file ",A47," is ready")', collfreqmeas_filename   

! values in spacecraft location only ----------------------------------------------------
! values in spacecraft location only ----------------------------------------------------
! values in spacecraft location only ----------------------------------------------------

!                           ! ----x----I----x----I----x----I----x----I----x--
     collfreqmeas_filename = 'op_NNN_collfreq_s1_avgw_NNN_eV_measured_osp.dat'
     collfreqmeas_filename(4:6)   = convert_int_to_txt_string(i_orbit_point, 3)
     collfreqmeas_filename(25:27) = convert_int_to_txt_string(int(w_eV), 3)

     i = pfl_i_spacecraft

     sumcoll_s1 = 0.0_8
     do kk = 1, N_colkind
        sumcoll_s1 = sumcoll_s1 + coll_freq_s1_L_colkind_colrange(i, kk, nr)
     end do

     open(19, file = collfreqmeas_filename)
     write (19, '("# collision frequencies averaged over energy interval ",i2," from ",f8.3," eV to ",f8.3," eV, calculated in simulation")') nr, min_w_eV_colrange(nr), max_w_eV_colrange(nr)
     write (19, '(59(2x,e14.7),2x,i4,2(2x,f12.3),2x,f6.1,2x,2(2x,e14.7))') &
          & coll_freq_s1_L_colkind_colrange(i, 1:N_colkind, nr), &   !1-59
          & i, &                                   ! 60
          & shared_pfl_point(i)%L_m * 0.001_8, &           ! 61
          & (pfl_point_coords(i)%r-1.0_8) * R_Earth_km, &  ! 62
          & w_eV, &       ! 63
          & sumcoll_s1, &    ! 64
          & density_m3_L_colrange(i,kk)  ! 65
        
     close (19, status = 'keep')
     print '("file ",A47," is ready")', collfreqmeas_filename   

  end do

END SUBROUTINE SAVE_PHOTOELECTRON_DATA
