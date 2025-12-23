!===========================================
PROGRAM MAIN

  USE ParallelOperationValues
  USE TimeValues
  USE GlobalIndices
  USE SpacecraftValues
  USE field_line
  USE Photoelectrons, ONLY : evdf

  USE MPI

  use msis_init, only : msisinit

  IMPLICIT NONE

  INTEGER ierr

  real SW(25)          ! MSIS switches

  INTEGER i, n
  INTEGER ALLOC_ERR

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, Rank_of_process, ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, N_of_processes, ierr)

  CALL PREPARE_PHOTOELECTRONS   ! reads input file
                                ! prepares e-n collision frequencies
                                ! sets random numbers (involves MPI communications)

  IF (Rank_of_process.eq.0) THEN
     CALL PREPARE_SPACECRAFT_DATA
     CALL UPDATE_TIME(orbit_point(1)%year, orbit_point(1)%month, orbit_point(1)%day_of_month, orbit_point(1)%ut_h)
     call read_apf107dat(year, month, day_of_month, INT(time_ut_h), Ap, f10p7, f10p7_pd, f10p7_81)   ! Ap, f10p7_pd, f10p7_81 used by gt8d (MSIS), f10p7 is for EUVAC
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  CALL MPI_BCAST(N_of_orbit_points, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(f10p7,    1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_BCAST(f10p7_81, 1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

  CALL Initiate_EUVAC_flux_crossections   ! uses f10p7 and f10p7_81

  DO i = 1, N_of_orbit_points

     if (Rank_of_process.eq.0) then
        CALL UPDATE_TIME(orbit_point(i)%year, orbit_point(i)%month, orbit_point(i)%day_of_month, orbit_point(i)%ut_h)
        CALL Prepare_GSE_to_GEO_transformation
        call read_ig_rz   ! from irifun.for
        call readapf107   ! from irifun.for
        SW = 1.0
        call msisinit(switch_legacy=SW)   ! required before first call of MSIS procedure
        call read_apf107dat(year, month, day_of_month, INT(time_ut_h), Ap, f10p7, f10p7_pd, f10p7_81)   ! Ap, f10p7_pd, f10p7_81 used by gtd8d (MSIS); f10p7 and f10p7_81 used by EUVAC
        CALL PREPARE_ENERGY_PITCH_ANGLE_RANGES(i)
        CALL PREPARE_FIELDLINE(i)   !####### transfer global indices here !########
     end if

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     CALL CALCULATE_PHOTOELECTRON_ENERGY_SPECTRA

!     if (Rank_of_process.eq.0) 
!     print *, Rank_of_process, " done CALCULATE_PHOTOELECTRON_ENERGY_SPECTRA"

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

     if (Rank_of_process.eq.0) then
        CALL SAVE_PHOTOELECTRON_DATA(i)
     end if

     IF (ALLOCATED(pfl_point_coords)) DEALLOCATE(pfl_point_coords, STAT = ALLOC_ERR)
     IF (ALLOCATED(shared_pfl_point)) THEN
        IF (N_ranges.GT.0) THEN
           DO n = 1, pfl_N_of_points
              DEALLOCATE(shared_pfl_point(n)%Ge_m2s_range, STAT = ALLOC_ERR)
              DEALLOCATE(shared_pfl_point(n)%Ge_par_m2s_range, STAT = ALLOC_ERR)

              DEALLOCATE(shared_pfl_point(n)%source_EUV_m3s_range, STAT = ALLOC_ERR)
              DEALLOCATE(shared_pfl_point(n)%source_photoe_m3s_range, STAT = ALLOC_ERR)
              DEALLOCATE(shared_pfl_point(n)%source_neutrals_m3s_range, STAT = ALLOC_ERR)
              DEALLOCATE(shared_pfl_point(n)%source_plasma_m3s_range, STAT = ALLOC_ERR)
              DEALLOCATE(shared_pfl_point(n)%sink_neutrals_m3s_range, STAT = ALLOC_ERR)
              DEALLOCATE(shared_pfl_point(n)%sink_plasma_m3s_range, STAT = ALLOC_ERR)
           END DO
        END IF
        DEALLOCATE(shared_pfl_point, STAT = ALLOC_ERR)
     END IF
     IF (ALLOCATED(evdf)) DEALLOCATE(evdf, STAT = ALLOC_ERR)

     CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

  END DO

  if (Rank_of_process.eq.0) then
     CALL SAVE_ORBIT
  end if

  CALL MPI_FINALIZE(ierr)

END PROGRAM MAIN

!-------------------------------
!
SUBROUTINE PREPARE_SPACECRAFT_DATA

  USE PhysicalConstants
  USE ParallelOperationValues
  USE SpacecraftValues

  USE MPI

  IMPLICIT NONE

  INTEGER ierr

  LOGICAL exists
  CHARACTER(1) buf
  INTEGER ALLOC_ERR
  INTEGER n

  INTEGER hours, minutes, seconds
  REAL(8) sat_h_km, sat_lat_deg, sat_lon_deg
  REAL loctime, sza, orbnum

  INQUIRE (FILE = 'input_spacecraft_orbit.dat', EXIST = exists)

  IF (.NOT.exists) THEN

     PRINT '(2x,"Process ",i4," : ERROR : input_spacecraft_orbit.dat not found. Program terminated")', Rank_of_process
     CALL MPI_ABORT(MPI_COMM_WORLD, 1, ierr)

  END IF

  IF (Rank_of_process.EQ.0) THEN
     PRINT '(2x,"Process ",i3," : input_spacecraft_orbit.dat is found. Reading the data file...")', Rank_of_process
  END IF

  OPEN (9, FILE = 'input_spacecraft_orbit.dat')

  READ (9, '(A1)') buf ! number of points along the spacecraft trajectory
  READ (9, *) N_of_orbit_points, orbit_number

  ALLOCATE(orbit_point(1:N_of_orbit_points), STAT = ALLOC_ERR)
 
  READ (9, '(A1)') buf ! year month day-of-month UT-hours/minutes/second altitutude(km) GEO-latitude/longitude(deg) LT SZA Num
 
  DO n = 1, N_of_orbit_points

     READ(9, *) orbit_point(n)%year, &
              & orbit_point(n)%month, &
              & orbit_point(n)%day_of_month, &
              & hours, &
              & minutes, &
              & seconds, &
              & sat_h_km, &
              & sat_lat_deg, &
              & sat_lon_deg, &
              & loctime, &     ! not used now
              & sza, &         ! not used now
              & orbnum

     IF (orbit_point(n)%year.LT.40) THEN
        orbit_point(n)%year = orbit_point(n)%year + 2000
     ELSE IF (orbit_point(n)%year.LT.99) THEN
        orbit_point(n)%year = orbit_point(n)%year + 1900
     END IF

     orbit_point(n)%ut_h = REAL(hours) + REAL(minutes * 60 + seconds) / 3600.0

     orbit_point(n)%r_km      = R_Earth_km + sat_h_km
     orbit_point(n)%theta_geo = MAX(0.0_8, MIN(pi, 0.5_8 * pi - sat_lat_deg * deg_to_rad))
     orbit_point(n)%phi_geo   = sat_lon_deg * deg_to_rad

     orbit_point(n)%x_geo_km = orbit_point(n)%r_km * SIN(orbit_point(n)%theta_geo) * COS(orbit_point(n)%phi_geo)
     orbit_point(n)%y_geo_km = orbit_point(n)%r_km * SIN(orbit_point(n)%theta_geo) * SIN(orbit_point(n)%phi_geo)
     orbit_point(n)%z_geo_km = orbit_point(n)%r_km * COS(orbit_point(n)%theta_geo)

  END DO

  orbit_number = INT(orbnum)

  CLOSE (9, STATUS = 'KEEP')

END SUBROUTINE PREPARE_SPACECRAFT_DATA

!-------------------------------
!
SUBROUTINE PREPARE_ENERGY_PITCH_ANGLE_RANGES(op_num)

  USE ParallelOperationValues
  USE field_line
  USE PhysicalConstants, ONLY : deg_to_rad

  IMPLICIT NONE

  INTEGER ierr

  INTEGER, INTENT(IN) :: op_num

  LOGICAL exists
  CHARACTER(1) buf
  INTEGER ALLOC_ERR

  CHARACTER(33) filename   ! input_energy_pa_ranges_op_NNN.dat
                           ! ----x----I----x----I----x----I---

  INTEGER i, ios

  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string
  end interface

! default values, expect no ranges

  IF (ALLOCATED(min_energy_eV_range)) DEALLOCATE(min_energy_eV_range)
  IF (ALLOCATED(max_energy_eV_range)) DEALLOCATE(max_energy_eV_range)
  IF (ALLOCATED(min_pitch_angle_range)) DEALLOCATE(min_pitch_angle_range)
  IF (ALLOCATED(max_pitch_angle_range)) DEALLOCATE(max_pitch_angle_range)

  N_ranges = 0

!             ----x----I----x----I----x----I---
  filename = 'input_energy_pa_ranges_op_NNN.dat'
  filename(27:29) = convert_int_to_txt_string(op_num, 3)

  INQUIRE (FILE = filename, EXIST = exists)
  IF (.NOT.exists) THEN
     
     PRINT '(2x,"##### file ",A33," not found #########")', filename
     RETURN

  END IF

  PRINT '(2x,"##### found file ",A33," , reading the data file...")', filename

  OPEN (9, FILE = filename)

  READ (9, '(A1)') buf ! # below is the number of energy ranges
  READ (9, *) N_ranges

  IF (N_ranges.lt.1) THEN
     N_ranges = MAX(0, N_ranges)  ! to avoid thinking what a negative N_ranges can do in the future
     RETURN
  END IF

  ALLOCATE(min_energy_eV_range(1:N_ranges), STAT = ALLOC_ERR)
  ALLOCATE(max_energy_eV_range(1:N_ranges), STAT = ALLOC_ERR)
  ALLOCATE(min_pitch_angle_range(1:N_ranges), STAT = ALLOC_ERR)
  ALLOCATE(max_pitch_angle_range(1:N_ranges), STAT = ALLOC_ERR)

  READ (9, '(A1)') buf ! # below, for each range provide in one line minimal and maximal energy (eV) and minimal and maximal pitch angles (deg)

  DO i = 1, N_ranges
     READ (9, *, iostat=ios) min_energy_eV_range(i), max_energy_eV_range(i), min_pitch_angle_range(i), max_pitch_angle_range(i)
     IF (ios.ne.0) THEN
        N_ranges = i-1
        EXIT
     END IF
  END DO

! convert degrees ot radians
  DO i = 1, N_ranges
     min_pitch_angle_range(i) = min_pitch_angle_range(i) * deg_to_rad
     max_pitch_angle_range(i) = max_pitch_angle_range(i) * deg_to_rad
  END DO

END SUBROUTINE PREPARE_ENERGY_PITCH_ANGLE_RANGES

!------------------------------------------
!
SUBROUTINE SAVE_ORBIT

  USE SpacecraftValues

  IMPLICIT NONE

!                                  ! ----x----I----x----I
  CHARACTER(20) orb_vtk_filename   ! orbit_NNNNNN_GEO.vtk

  INTEGER i

  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
       integer length_of_string
     end function convert_int_to_txt_string
  end interface

!                   ! ----x----I----x----I
  orb_vtk_filename = 'orbit_NNNNNN_GEO.vtk'
  orb_vtk_filename(7:12) = convert_int_to_txt_string(orbit_number, 6)

  open (19, file = orb_vtk_filename)

  write (19, '("# vtk DataFile Version 4.2")')
  write (19, '(A14)') orb_vtk_filename(1:16)
  write (19, '("ASCII")')
  write (19, '("DATASET POLYDATA")')
  write (19, '("POINTS ",i6," float")') N_of_orbit_points

  do i = 1, N_of_orbit_points
     write (19, '(3(2x,f10.2))') orbit_point(i)%x_geo_km, orbit_point(i)%y_geo_km, orbit_point(i)%z_geo_km
  end do

  write (19, '("LINES ",i4,2x,i6)') 1, N_of_orbit_points+1
  write (19, '(2x,i6)') N_of_orbit_points
  do i = 0, N_of_orbit_points-1
     write (19, '(2x,i6)') i
  end do

  write (19, '("POINT_DATA ",i6)') N_of_orbit_points
  write (19, '("SCALARS totEUVm3s float")')
  write (19, '("LOOKUP_TABLE default")')

  do i = 1, N_of_orbit_points
     write (19, '(2x,e12.5)') orbit_point(i)%total_euv_rate_m3s
  end do
  
  close (19, status = 'keep')
  print '("file ",A20," is ready")', orb_vtk_filename
  
END SUBROUTINE SAVE_ORBIT
