!-------------------------------
!
program get_photoelectron_fluxes

  implicit none

  real, parameter :: R_Earth_km  = 6371.0

  REAL(8), PARAMETER :: e_Cl     = 1.602189d-19                          ! Elementary charge [Cl]
  REAL(8), PARAMETER :: m_e_kg   = 9.109534d-31                          ! Electron mass [kg]

  real pi

  real rad_to_deg, deg_to_rad
  real factor_eedf

  logical exists
  character(1) buf

  real instr_fov_angle
  real fov_axis_length_km

  integer max_eedf_N_w
  real dw_eV
  integer max_eedf_N_theta

  integer save_edfepa_flag
  logical save_edfepa

  real, allocatable :: w_eV(:)
  integer alloc_err
  integer m

  real, allocatable :: theta(:)
  real dtheta
  integer n

  integer N_ranges

  integer, parameter :: max_N_ranges = 20

  integer m_start(1:max_N_ranges)
  integer m_end(1:max_N_ranges)
  integer i
  real energy_min_eV, energy_max_eV

  character(38) flux_filename   ! photoelectron_fluxes_range_NN_igrf.dat
!                               ! ----x----I----x----I----x----I----x---

  integer N_of_orbit_points

  real, allocatable :: theta_instr(:)  ! angle between the axis of the field of view of photoelectron instrument and the goemagnetic field

  real, allocatable :: x_geo_km(:)
  real, allocatable :: y_geo_km(:)
  real, allocatable :: z_geo_km(:)

  real, allocatable :: bx(:)
  real, allocatable :: by(:)
  real, allocatable :: bz(:)

  real, allocatable :: cont_time_h(:)  ! continuous time, useful when the UT changes from 23:59:59 to 00:00:00 during the orbit

  integer k
  integer year, month, day_of_month, hours, minutes, seconds
  real sat_h_km, sat_lat_deg, sat_lon_deg

  real r_km, theta_geo, phi_geo

  REAL(8) Bnorth_nT, Beast_nT, Bvert_nT, intensity_nT

  REAL Br_nT, Btheta_nT, Bphi_nT
  real temp
  REAL Bx_nT, By_nT, Bz_nT

  real vx, vy, vz
  real zsp_x, zsp_y, zsp_z   ! z-ort of spacecraft CS in GEO
  real xsp_x, xsp_y, xsp_z   ! x-ort of spacecraft CS in GEO
  real ysp_x, ysp_y, ysp_z   ! y-ort of spacecraft CS in GEO

  character(18) viewray_filename

  real, allocatable :: eedf(:,:)

  real ut_h

  character(15) op_evdf_filename

  integer max_evdf_N_vpar
  integer max_evdf_N_vperp
  real dvpar_ms
  real dvperp_ms

  real, allocatable :: vpar_ms(:)
  real, allocatable :: vperp_ms(:)
  real, allocatable :: evdf(:,:)

  integer npar, nperp

  integer idummy
  real rdummy

  real myvpar_ms
  real myvperp_ms

  character(17) edfepa_filename   ! op_NNN_edfepa.dat'
!                                 ! ----x----I----x--

  real theta_min, theta_max

  integer n_start
  integer n_end

  character(21) edfeavgpa_filename      ! op_NNN_edfe_avgpa.dat
!                                       ! ----x----I----x----I-

  real temp1, temp2, temp3

  real avg_factor_sr

  real particle_flux_m2s
  real energy_flux_eVm2s
  real avg_diff_particle_flux_m2seVsr
  real avg_diff_energy_flux_eVm2seVsr

! functions

  real get_evdf_value

  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       integer length_of_string
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
     end function convert_int_to_txt_string
  end interface

  pi = acos(-1.0)
  rad_to_deg = 180.0 / pi
  deg_to_rad = pi / 180.0
  factor_eedf = sqrt(2.0) * (e_cl / m_e_kg)**1.5

! read the control file

  inquire (file = 'input_photoelectron_fluxes_igrf.dat', exist = exists)
  if (.not.exists) then
     print '("error : file input_photoelectron_fluxes_igrf.dat not found")'
     stop
  end if

  open (9, file = 'input_photoelectron_fluxes_igrf.dat')

  read (9, '(A1)') buf   ! # spacecraft frame: z is vertical
  read (9, '(A1)') buf   ! # z is in the vertical direction
  read (9, '(A1)') buf   ! # x is in vertical plane containing spacecraft velocity, in the direction of spacecraft motion
  read (9, '(A1)') buf   ! # y is vector product of z and x
  read (9, '(A1)') buf   ! # it is assumed that photoelectron sensor is pointed along the y-axis in the spacecraft frame, as sensor #1 on AE-E
  read (9, '(A1)') buf   ! # full angle of the field of view of the spacecraft photoelectron sensor (deg)
  read (9, *) instr_fov_angle

  instr_fov_angle = 0.5 * instr_fov_angle * deg_to_rad  ! we need half of the full fov in radians

  read (9, '(A1)') buf   ! # for plotting, length of line representing sensor's field of view axis (km)
  read (9, *) fov_axis_length_km

  read (9, '(A1)') buf   ! # for electron distribution over energy and pitch angle define below
  read (9, '(A1)') buf   ! # number of energy bins
  read (9, *) max_eedf_N_w

  read (9, '(A1)') buf   ! # size of energy bin (eV)
  read (9, *) dw_eV

  read (9, '(A1)') buf   ! # number of pitch angle bins (pitch angle varies between 0 and pi)
  read (9, *) max_eedf_N_theta

  read (9, '(A1)') buf   ! # save the electron distribution over energy and pitch angle into a file for each orbit point (1/0 = Yes/No)
  read (9, *) save_edfepa_flag

  save_edfepa = .true.
  if (save_edfepa_flag.eq.0) save_edfepa = .false.

! define energy bins
  allocate(w_eV(0:max_eedf_N_w), stat = alloc_err)
  do m = 0, max_eedf_N_w
     w_eV(m) = (REAL(m) + 0.5) * dw_eV
  end do

! define pitch angle bins
  allocate(theta(0:max_eedf_N_theta), stat = alloc_err)
  dtheta = pi/(max_eedf_N_theta+1)
  do n = 0, max_eedf_N_theta
     theta(n) = (REAL(n) + 0.5) * dtheta
  end do

  read (9, '(A1)') buf   ! # number of energy-pitch angle ranges where particle fluxes will be calculated (<=20)
  read (9, *) N_ranges

  N_ranges = min(N_ranges, max_N_ranges)

  read (9, '(A1)') buf   ! # for each range give (in one line) minimal/ maximal energy (eV) 

  m_start = -1
  m_end = -1

  do i = 1, N_ranges
     
     read (9, *) energy_min_eV, energy_max_eV

! calculate range of energy bins
     m_start(i) = max_eedf_N_w
     do m = 0, max_eedf_N_w
        if (w_eV(m).ge.energy_min_eV) then
           m_start(i) = m
           exit
        end if
     end do

     m_end(i) = 0
     do m = max_eedf_N_w, 0, -1
        if (w_eV(m).le.energy_max_eV) then
           m_end(i) = m
           exit
        end if
     end do

! consistency check
     if ( (m_start(i).lt.0).or. &
        & (m_end(i).gt.max_eedf_N_w).or. &
        & (m_start(i).gt.m_end(i)) ) then
        print '("error : inconsistent energy limits for range ",i2," :: ",2(2x,f10.3),2x,2(2x,i4))', &
             & i, energy_min_eV, energy_max_eV, m_start(i), m_end(i)
        stop
     end if
 
     flux_filename = 'photoelectron_fluxes_range_NN_igrf.dat'
!                   ! ----x----I----x----I----x----I----x---
     flux_filename(28:29) = convert_int_to_txt_string(i, 2)

     open (10, file = flux_filename, status = 'replace')
! header
     write (10, '("# range ",i2," parameters ::")') i
     write (10, '("# energy_min = ",f8.2," eV, m_start = ",i4," w_eV(m_start) = ",f8.2," eV")') energy_min_eV, m_start(i), w_eV(m_start(i))
     write (10, '("# energy_max = ",f8.2," eV, m_end   = ",i4," w_eV(m_end)   = ",f8.2," eV")') energy_max_eV, m_end(i),   w_eV(m_end(i))
     write (10, '("# col 1 is orbit point number (int)")') 
     write (10, '("# col 2 is UT (h)")')
     write (10, '("# col 3 is continuous UT (h)")')
     write (10, '("# col 4 is integrated particle flux (1/m^2/s)")')
     write (10, '("# col 5 is integrated particle energy flux (eV/m^2/s)")')
     write (10, '("# col 6 is differential particle flux averaged over the energy and solid angle (1/m^2/s/eV/sr)")')
     write (10, '("# col 7 is differential particle energy flux averaged over the energy and solid angle (eV/m^2/s/eV/sr)")')
     write (10, '("# col 8 is angle between the axis of the photoelectron instrument sensor and the positive direction of geomagnetic field (deg)")')
     write (10, '("# col 9 is minimal desired pitch angle (deg)")')
     write (10, '("# col 10 is maximal desired pitch angle (deg)")')
     write (10, '("# col 11 is minimal pitch angle bin (deg)")')
     write (10, '("# col 12 is maximal pitch angle bin (deg)")')
     write (10, '("# col 13 is minimal pitch angle bin number (int)")')
     write (10, '("# col 14 is maximal pitch angle bin number (int)")')
     close (10, status = 'keep')

  end do   !###   do i = 1, N_ranges

  close (9, status = 'keep')   !###   open (9, file = 'input_photoelectron_fluxes.dat')

  print '("finished reading input_photoelectron_fluxes.dat")'

  INQUIRE (FILE = 'input_spacecraft_orbit.dat', EXIST = exists)

  IF (.NOT.exists) THEN
     PRINT '(2x,"error : input_spacecraft_orbit.dat not found")'
     STOP
  END IF

  PRINT '(2x,"input_spacecraft_orbit.dat is found. Reading the data file...")'

! prepare data of spacecraft instrument orientation relative to the geomagnetic field

! first read the orbit data file to get coordiates of all orbit points, necessary to claculate velocity vector
! also calculate vector of direction of geomagnetic field

  OPEN (9, FILE = 'input_spacecraft_orbit.dat')

  READ (9, '(A1)') buf ! number of points along the spacecraft trajectory
  READ (9, *) N_of_orbit_points
  READ (9, '(A1)') buf ! year month day-of-month UT-hours/minutes/second altitutude(km) GEO-latitude/longitude(deg)

  ALLOCATE(theta_instr(1:N_of_orbit_points), STAT = ALLOC_ERR)

  ALLOCATE(x_geo_km(1:N_of_orbit_points), STAT = ALLOC_ERR)
  ALLOCATE(y_geo_km(1:N_of_orbit_points), STAT = ALLOC_ERR)
  ALLOCATE(z_geo_km(1:N_of_orbit_points), STAT = ALLOC_ERR)

  ALLOCATE(bx(1:N_of_orbit_points), STAT = ALLOC_ERR)
  ALLOCATE(by(1:N_of_orbit_points), STAT = ALLOC_ERR)
  ALLOCATE(bz(1:N_of_orbit_points), STAT = ALLOC_ERR)

  ALLOCATE(cont_time_h(1:N_of_orbit_points), STAT = ALLOC_ERR)

  DO k = 1, N_of_orbit_points

! read orbit point ----------------------------------------------------------------------------------------------
     READ(9, *) year, &
              & month, &
              & day_of_month, &
              & hours, &
              & minutes, &
              & seconds, &
              & sat_h_km, &
              & sat_lat_deg, &
              & sat_lon_deg

     cont_time_h(k) = REAL(hours) + REAL(minutes * 60 + seconds) / 3600.0

     IF (year.LT.40) THEN
        year = year + 2000
     ELSE IF (year.LT.99) THEN
        year = year + 1900
     END IF

     r_km      = R_Earth_km + sat_h_km
     theta_geo = real(MAX(0.0, MIN(pi, 0.5 * pi - sat_lat_deg * deg_to_rad)))
     phi_geo   = sat_lon_deg * deg_to_rad

     x_geo_km(k) = r_km * SIN(theta_geo) * COS(phi_geo)
     y_geo_km(k) = r_km * SIN(theta_geo) * SIN(phi_geo)
     z_geo_km(k) = r_km * COS(theta_geo)

! 0 means main-field values are required
! 2 means geocentric CS (the Earth is a sphere)
     CALL igrf13syn(0, dble(year), 2, dble(r_km), dble(90.0 - sat_lat_deg), dble(sat_lon_deg), Bnorth_nT, Beast_nT, Bvert_nT, intensity_nT)

     Br_nT     = -Bvert_nT
     Btheta_nT = -Bnorth_nT
     Bphi_nT   =  Beast_nT

     CALL convert_spherical_vector_to_cart_4(theta_geo, phi_geo, Br_nT, Btheta_nT, Bphi_nT, Bx_nT, By_nT, Bz_nT)  ! gives B{xyz}_nT in GEO

     temp = sqrt(Bx_nT**2 + By_nT**2 + Bz_nT**2)
     bx(k) = Bx_nT / temp   !
     by(k) = By_nT / temp   ! direction along the field line towards the northern hemisphere 
     bz(k) = Bz_nT / temp   !

  END DO

  CLOSE (9, STATUS = 'KEEP')

! calculate angle between the axis of the instrument and the geomagnetic field
  
  do k = 1, N_of_orbit_points
     
     vx = x_geo_km(min(k+1,N_of_orbit_points)) - x_geo_km(max(k-1,1))
     vy = y_geo_km(min(k+1,N_of_orbit_points)) - y_geo_km(max(k-1,1))
     vz = z_geo_km(min(k+1,N_of_orbit_points)) - z_geo_km(max(k-1,1))

     temp = sqrt(vx**2 + vy**2 + vz**2)
     vx = vx / temp
     vy = vy / temp
     vz = vz / temp
         
! unitary vector in the radial direction, z-ort in spacecraft CS
     temp = sqrt(x_geo_km(k)**2 + y_geo_km(k)**2 + z_geo_km(k)**2)
     zsp_x = x_geo_km(k) / temp
     zsp_y = y_geo_km(k) / temp
     zsp_z = z_geo_km(k) / temp

! x-ort in spacecraft CS
     temp = zsp_x*vx + zsp_y*vy + zsp_z*vz
     xsp_x = vx - zsp_x * temp
     xsp_y = vy - zsp_y * temp
     xsp_z = vz - zsp_z * temp
     temp = sqrt(xsp_x**2 + xsp_y**2 + xsp_z**2)
     xsp_x = xsp_x / temp
     xsp_y = xsp_y / temp
     xsp_z = xsp_z / temp

! y-ort in spacecraft CS
! this is the direction of the axis of the PES sensor 1 on AE-E
     ysp_x = zsp_y * xsp_z - zsp_z * xsp_y
     ysp_y = zsp_z * xsp_x - zsp_x * xsp_z
     ysp_z = zsp_x * xsp_y - zsp_y * xsp_x
     temp = sqrt(ysp_x**2 + ysp_y**2 + ysp_z**2)
     ysp_x = ysp_x / temp
     ysp_y = ysp_y / temp
     ysp_z = ysp_z / temp

! angle between the geomagnetic field and the sensor axis
     theta_instr(k) = acos(max(-1.0, min(1.0, bx(k) * ysp_x + by(k) * ysp_y + bz(k) * ysp_z)))

! save instrument's field-of-view axis
     viewray_filename = 'op_NNN_viewray.vtk'
!                      ! ----x----I----x---
     viewray_filename(4:6) = convert_int_to_txt_string(k, 3)

     open (10, file = viewray_filename)
     write (10, '("# vtk DataFile Version 4.2")')
     write (10, '(A14)') viewray_filename(1:14)
     write (10, '("ASCII")')
     write (10, '("DATASET POLYDATA")')
     write (10, '("POINTS ",i6," float")') 2 !total_N_of_ray_points

     write (10, '(3(2x,f10.2))') x_geo_km(k), y_geo_km(k), z_geo_km(k)
     write (10, '(3(2x,f10.2))') x_geo_km(k) + fov_axis_length_km * ysp_x, &
                               & y_geo_km(k) + fov_axis_length_km * ysp_y, &
                               & z_geo_km(k) + fov_axis_length_km * ysp_z

     write (10, '("LINES ",i4,2x,i6)') 1, 2+1

!     pos=0   ! yes, zero, because numbering of points in paraview starts with zero
     write (10, '(2x,i6)') 2
     write (10, '(2x,i6)') 0
     write (10, '(2x,i6)') 1
!     do nn = 1, 2
!        write (10, '(2x,i6)') pos
!        pos=pos+1
!     end do

     close (10, STATUS = 'KEEP')

     print '("### file ",A18," is ready")', viewray_filename

  end do

! below we only need theta_instr
! cleanup
  deallocate(x_geo_km, stat = alloc_err)
  deallocate(y_geo_km, stat = alloc_err)
  deallocate(z_geo_km, stat = alloc_err)

  deallocate(bx, stat = alloc_err)
  deallocate(by, stat = alloc_err)
  deallocate(bz, stat = alloc_err)

! prepare the continuous time
  if (cont_time_h(1).gt.cont_time_h(N_of_orbit_points)) then
! the orbit segment begins before UT 00:00 and ends after UT 00:00 
     do k = 1, N_of_orbit_points-1
        if (cont_time_h(k).gt.cont_time_h(N_of_orbit_points)) cont_time_h(k) = cont_time_h(k)-24.0
     end do
  end if

! read the orbit data file the second time
! now we can calculate the fluxes

  OPEN (9, FILE = 'input_spacecraft_orbit.dat')

  READ (9, '(A1)') buf ! number of points along the spacecraft trajectory
  READ (9, *) N_of_orbit_points

  READ (9, '(A1)') buf ! year month day-of-month UT-hours/minutes/second altitutude(km) GEO-latitude/longitude(deg)

  allocate(eedf(0:max_eedf_N_w, 0:max_eedf_N_theta), stat = alloc_err)

! prepare output files for photoelectron fluxes

  DO k = 1, N_of_orbit_points

! read orbit point ----------------------------------------------------------------------------------------------
     READ(9, *) year, &
              & month, &
              & day_of_month, &
              & hours, &
              & minutes, &
              & seconds, &
              & sat_h_km, &
              & sat_lat_deg, &
              & sat_lon_deg

     ut_h = REAL(hours) + REAL(minutes * 60 + seconds) / 3600.0

! read and process EVDF for each orbit point --------------------------------------------------------------------

! read evdf
     op_evdf_filename = 'op_NNN_evdf.dat'
!                      ! ----x----I----x
     op_evdf_filename(4:6) = convert_int_to_txt_string(k, 3)

     inquire(file = op_evdf_filename, exist = exists)
     if (.not.exists) then
        print '("file ",A15," not found, skip")', op_evdf_filename
        cycle
     end if

     OPEN (20, FILE = op_evdf_filename)
        
!     read (20, '("# max_evdf_N_vperp = ",i4)') max_evdf_N_vperp
!                  ----x----I----x----I-
!     read (20, '("# max_evdf_N_vpar = ",i4)') max_evdf_N_vpar
!                  ----x----I----x----I
!     read (20, '("# dvpar_ms = ",f10.2)') dvpar_ms
!                  ----x----I---
!     read (20, '("# dvperp_ms = ",f10.2)') dvperp_ms  
!                  ----x----I----

     read (20, '(21x,i4)') max_evdf_N_vperp
     read (20, '(20x,i4)') max_evdf_N_vpar
     read (20, '(13x,f10.2)') dvpar_ms
     read (20, '(14x,f10.2)') dvperp_ms  
     read (20, '(A1)') buf  ! skip the rest of the header with the description of columns
     read (20, '(A1)') buf
     read (20, '(A1)') buf
     read (20, '(A1)') buf
     read (20, '(A1)') buf

     allocate(vpar_ms(-max_evdf_N_vpar:max_evdf_N_vpar), stat = alloc_err)
     allocate(vperp_ms(0:max_evdf_N_vperp), stat = alloc_err)
     allocate(evdf(-max_evdf_N_vpar:max_evdf_N_vpar,0:max_evdf_N_vperp), stat = alloc_err)

     do npar = -max_evdf_N_vpar, max_evdf_N_vpar
        vpar_ms(npar) = npar * dvpar_ms
     end do

     do nperp = 0, max_evdf_N_vperp
        vperp_ms(nperp) = (REAL(nperp) + 0.5) * dvperp_ms
     end do

     do nperp = 0, max_evdf_N_vperp
        do npar = -max_evdf_N_vpar, max_evdf_N_vpar

           read (20, '(2x,i4,2x,i4,2x,f12.1,2x,f12.1,2x,e12.5)') &
                & idummy, &
                & idummy, &
                & rdummy, &
                & rdummy, &
                & evdf(npar,nperp)
        end do
        read (20, '(A1)') buf
     end do

     close (20, status = 'keep')
     print '("reading file ",A15," is done")', op_evdf_filename

! convert into evdf in cylindrical CS in the velocity space
  
!     do nperp = 0, max_evdf_N_vperp
!        do npar = -max_evdf_N_vpar, max_evdf_N_vpar
!           evdf(npar,nperp) = evdf(npar,nperp) / (2.0 * pi * vperp_ms(nperp))
!        end do
!     end do

! get distribution over energy and pitch angle

     do n = 0, max_eedf_N_theta
        do m = 0, max_eedf_N_w
           myvpar_ms  = real(sqrt(2.0_8 * e_Cl * dble(w_eV(m)) / m_e_kg)) * cos(theta(n))
           myvperp_ms = real(sqrt(2.0_8 * e_Cl * dble(w_eV(m)) / m_e_kg)) * sin(theta(n))

!           eedf(m,n) = get_evdf_value(myvpar_ms, myvperp_ms, max_evdf_N_vpar, max_evdf_N_vperp, vpar_ms, vperp_ms, evdf) * sqrt(w_eV(m)) * factor_eedf

           eedf(m,n) = get_evdf_value(myvpar_ms, myvperp_ms, max_evdf_N_vpar, max_evdf_N_vperp, vpar_ms, vperp_ms, evdf) * (e_Cl / (2.0_8 * pi * m_e_kg * sin(theta(n))))

        end do
     end do

! if requested, save electron distribution function over energy and pitch angle (edfepa)

     if (save_edfepa) then
        edfepa_filename = 'op_NNN_edfepa.dat'
!                        ! ----x----I----x--
        edfepa_filename(4:6) = convert_int_to_txt_string(k, 3)
        open (21, file = edfepa_filename, status = 'replace')
! header
        write (21, '("# col 1 is energy bin number (int)")')
        write (21, '("# col 2 is pitch angle bin number (int)")')
        write (21, '("# col 3 is energy bin middle energy (eV)")')
        write (21, '("# col 4 is pitch angle bin middle angle (rad)")')
        write (21, '("# col 5 is electron distribution function over energy and pitch angle (1/m^3/eV/sr)")')
        write (21, '("# col 6 is electron distribution function over energy and pitch angle multiplied by speed (1/m^2/s/eV/sr)")')
        write (21, '("# col 7 is electron distribution function over energy and pitch angle multiplied by speed and energy (eV/m^2/s/eV/sr)")')
        do n = 0, max_eedf_N_theta
           do m = 0, max_eedf_N_w
              write (21, '(2x,i4,2x,i4,2x,f8.2,2x,f8.6,3(2x,e12.5))') &
                   & m, &
                   & n, &
                   & w_eV(m), &
                   & theta(n), &
                   & eedf(m,n), &                                                    ! m^-3/eV/sr
                   & eedf(m,n) * sqrt(2.0_8 * e_Cl * w_eV(m) / m_e_kg), &            ! (m/s)m^-3/eV/sr = m^-2/eV/s/sr
                   & eedf(m,n) * w_eV(m) * sqrt(2.0_8 * e_Cl * w_eV(m) / m_e_kg)     ! eV(m/s)m^-3/eV/sr = m^-2/s/sr
           end do
           write (21, '(" ")')
        end do
        close (21, status = 'keep')
        print '("file ",A20," is ready")', edfepa_filename
     end if

! pitch angle limits
     theta_min = max(0.0, min(pi, pi - theta_instr(k) - instr_fov_angle))
     theta_max = max(0.0, min(pi, pi - theta_instr(k) + instr_fov_angle))

! calculate range of theta bins
     n_start = max_eedf_N_theta
     do n = 0, max_eedf_N_theta
        if (theta(n).ge.theta_min) then
           n_start = n
           exit
        end if
     end do

     n_end = 0
     do n = max_eedf_N_theta, 0, -1
        if (theta(n).le.theta_max) then
           n_end = n
           exit
        end if
     end do

! consistency check
     if ( (n_start.lt.0).or. &
        & (n_start.gt.max_eedf_N_theta).or. &
        & (n_end.lt.0).or. &
        & (n_end.gt.max_eedf_N_theta).or. &
        & (n_start.gt.n_end) ) then
        print '("error : inconsistent limits for pitch angle range ",i2," :: ",2(2x,f10.3),2x,2(2x,i4))', &
             & i, theta_min, theta_max, n_start, n_end
        stop
     end if
 
! prepare averaging factor for averaging over the solid angle
     avg_factor_sr = 0.0
     do n = n_start, n_end
        avg_factor_sr = avg_factor_sr + sin(theta(n))
     end do
     avg_factor_sr = avg_factor_sr * 2.0 * pi * dtheta

! if requested, save electron distribution function over energy averaged over pitch angle (edfe_avgpa)

     if (save_edfepa) then

! averaging over the Solid Angle

        edfeavgpa_filename = 'op_NNN_edfe_avgSA.dat'
!                           ! ----x----I----x----I-
        edfeavgpa_filename(4:6) = convert_int_to_txt_string(k, 3)
        open (21, file = edfeavgpa_filename, status = 'replace')
! header
        write (21, '("# averaging over the SOLID ANGLE")')
        write (21, '("# col 1 is energy bin number (int)")')
        write (21, '("# col 2 is energy bin middle energy (eV)")')
        write (21, '("# values in columns 3-5 are averaged over pitch angle bins ",i4," to ",i4,", corresponding angles ",f7.3," deg to ",f7.3," deg")') n_start, n_end, theta(n_start) * rad_to_deg, theta(n_end) * rad_to_deg
        write (21, '("# col 3 is electron distribution function over energy and pitch angle (1/m^3/eV/sr)")')
        write (21, '("# col 4 is electron distribution function over energy and pitch angle multiplied by speed (1/m^2/s/eV/sr)")')
        write (21, '("# col 5 is electron distribution function over energy and pitch angle multiplied by speed and energy (eV/m^2/s/eV/sr)")')
        do m = 0, max_eedf_N_w
           temp1 = 0.0
           temp2 = 0.0
           temp3 = 0.0
           do n =  n_start, n_end  !0, max_eedf_N_theta
              temp1 = temp1 + eedf(m,n) * sin(theta(n))                                                    ! m^-3/eV/sr
              temp2 = temp2 + eedf(m,n) * sin(theta(n)) * sqrt(2.0_8 * e_Cl * w_eV(m) / m_e_kg)            ! (m/s)m^-3/eV/sr = m^-2/eV/s/sr
              temp3 = temp3 + eedf(m,n) * sin(theta(n)) * w_eV(m) * sqrt(2.0_8 * e_Cl * w_eV(m) / m_e_kg)  ! eV(m/s)m^-3/eV/sr = m^-2/s/sr
           end do
           temp1 = temp1 * 2.0 * pi * dtheta / avg_factor_sr
           temp2 = temp2 * 2.0 * pi * dtheta / avg_factor_sr
           temp3 = temp3 * 2.0 * pi * dtheta / avg_factor_sr
           write (21, '(2x,i4,2x,f8.2,3(2x,e12.5))') &
                   & m, &
                   & w_eV(m), &
                   & temp1, temp2, temp3
        end do
        close (21, status = 'keep')
        print '("file ",A21," is ready")', edfeavgpa_filename

     end if

! calculate fluxes for each range

     do i = 1, N_ranges

        particle_flux_m2s = 0.0
        energy_flux_eVm2s = 0.0

        do m = m_start(i), m_end(i)
           do n = n_start, n_end
              particle_flux_m2s = particle_flux_m2s + eedf(m,n)           * sqrt(w_eV(m)) * sin(theta(n))
              energy_flux_eVm2s = energy_flux_eVm2s + eedf(m,n) * w_eV(m) * sqrt(w_eV(m)) * sin(theta(n))

           end do
        end do
     
        particle_flux_m2s = particle_flux_m2s * 2.0 * pi * dtheta * dw_eV * sqrt(2.0_8 * e_Cl / m_e_kg)
        energy_flux_eVm2s = energy_flux_eVm2s * 2.0 * pi * dtheta * dw_eV * sqrt(2.0_8 * e_Cl / m_e_kg)

        avg_diff_particle_flux_m2seVsr = particle_flux_m2s / (avg_factor_sr * dw_eV * (m_end(i) - m_start(i) + 1))
        avg_diff_energy_flux_eVm2seVsr = energy_flux_eVm2s / (avg_factor_sr * dw_eV * (m_end(i) - m_start(i) + 1))

        flux_filename = 'photoelectron_fluxes_range_NN_igrf.dat'
!                      ! ----x----I----x----I----x----I----x---
        flux_filename(28:29) = convert_int_to_txt_string(i, 2)

        open (10, file = flux_filename, position = 'append')
        write (10, '(2x,i4,2(2x,f6.3),4(2x,e14.7),5(2x,f6.2),2(2x,i4))') &
             & k, &                       ! 1
             & ut_h, &                          ! 2
             & cont_time_h(k), &                   ! 3
             & particle_flux_m2s, &                   ! 4
             & energy_flux_eVm2s, &                   ! 5
             & avg_diff_particle_flux_m2seVsr, &           ! 6
             & avg_diff_energy_flux_eVm2seVsr, &      ! 7
             & theta_instr(k) * rad_to_deg, &   ! 8
             & theta_min * rad_to_deg, &              ! 9
             & theta_max * rad_to_deg, &              ! 10
             & theta(n_start) * rad_to_deg, &   ! 11
             & theta(n_end) * rad_to_deg, &     ! 12
             & n_start, &                             ! 13
             & n_end                                  ! 14
        close (10, status = 'keep')

     end do

! cleanup
     if (allocated(vpar_ms)) deallocate(vpar_ms, stat = alloc_err)
     if (allocated(vperp_ms)) deallocate(vperp_ms, stat = alloc_err)
     if (allocated(evdf)) deallocate(evdf, stat = alloc_err)

  END DO   !###     DO k = 1, N_of_orbit_points

  close (9, status = 'keep')   ! OPEN (9, FILE = 'input_spacecraft_orbit.dat')

! final cleanup
  if (allocated(w_eV)) deallocate(w_eV, stat = alloc_err)
  if (allocated(theta)) deallocate(theta, stat = alloc_err)
  if (allocated(cont_time_h)) deallocate(cont_time_h, stat = alloc_err)
  if (allocated(eedf)) deallocate(eedf, stat = alloc_err)

end program get_photoelectron_fluxes

!---------------------------------
!
real function get_evdf_value(vpar_ms, vperp_ms, max_evdf_N_vpar, max_evdf_N_vperp, arr_vpar_ms, arr_vperp_ms, arr_evdf)

  implicit none

  real vpar_ms, vperp_ms
  integer max_evdf_N_vpar, max_evdf_N_vperp
  real arr_vpar_ms(-max_evdf_N_vpar:max_evdf_N_vpar)
  real arr_vperp_ms(0:max_evdf_N_vperp)
  real arr_evdf(-max_evdf_N_vpar:max_evdf_N_vpar, 0:max_evdf_N_vperp)

  integer i, npar, nperp
  real anpar, anparp1
  real anperp, anperpp1

  get_evdf_value = 0.0

  if (vpar_ms.le.arr_vpar_ms(-max_evdf_N_vpar)) return

  if (vpar_ms.ge.arr_vpar_ms(max_evdf_N_vpar)) return

  if (vperp_ms.ge.arr_vperp_ms(max_evdf_N_vperp)) return

  do i = -max_evdf_N_vpar, max_evdf_N_vpar-1
     if ((vpar_ms.ge.arr_vpar_ms(i)).and.(vpar_ms.le.arr_vpar_ms(i+1))) then
        npar = i

        anpar = (arr_vpar_ms(npar+1) - vpar_ms) / (arr_vpar_ms(npar+1) - arr_vpar_ms(npar))
        anparp1 = (vpar_ms - arr_vpar_ms(npar)) / (arr_vpar_ms(npar+1) - arr_vpar_ms(npar))

        exit
     end if
  end do

  if (vperp_ms.le.arr_vperp_ms(0)) then
     nperp=0
     anperp=1.0
     anperpp1 = 0.0
  else
     do i = 0, max_evdf_N_vperp-1
        if ((vperp_ms.ge.arr_vperp_ms(i)).and.(vperp_ms.le.arr_vperp_ms(i+1))) then
           nperp = i

           anperp = (arr_vperp_ms(nperp+1) - vperp_ms) / (arr_vperp_ms(nperp+1) - arr_vperp_ms(nperp))
           anperpp1 = (vperp_ms - arr_vperp_ms(nperp)) / (arr_vperp_ms(nperp+1) - arr_vperp_ms(nperp))

           exit
        end if
     end do
  end if

  get_evdf_value = anpar   * anperp   * arr_evdf(npar,  nperp) + &
                 & anparp1 * anperp   * arr_evdf(npar+1,nperp) + &
                 & anpar   * anperpp1 * arr_evdf(npar,  nperp+1) + &
                 & anparp1 * anperpp1 * arr_evdf(npar+1,nperp+1)

  return

end function get_evdf_value

!-----------------------------------
! creates a string of length "length_of_string" out of an integer number "int_number"
!
function convert_int_to_txt_string(int_number, length_of_string)

  implicit none

  integer length_of_string
  character*(length_of_string) convert_int_to_txt_string
  integer int_number
  character(5) format_string
  character(2) length2_txt
  character(1) length1_txt

  character*(length_of_string) number_txt
  
  integer blanks_number
  integer i

! produce format string
  if ((length_of_string.gt.0).and.(length_of_string.lt.10)) then
     write (length1_txt, '(i1)') length_of_string
     format_string = '(iN) '
     format_string(3:3) = length1_txt
  else if ((length_of_string.ge.10).and.(length_of_string.lt.100)) then
     write (length2_txt, '(i2)') length_of_string
     format_string = '(iNN)'
     format_string(3:4) = length2_txt
  else
     print *, "ERROR in CONVERT_INT_TO_TXT_STRING:"
     print *, "incorrect string length requested: ", length_of_string
     stop
  end if

  WRITE (number_txt, format_string) int_number
  number_txt = ADJUSTL(TRIM(number_txt))
  blanks_number = length_of_string - LEN_TRIM(number_txt)
  number_txt = ADJUSTR(number_txt)
  do i = 1, blanks_number
     number_txt(i:i) = '0'
  end do

  convert_int_to_txt_string = number_txt

end function convert_int_to_txt_string

!--------------------------------------------
!
subroutine convert_spherical_vector_to_cart_4(theta, phi, ar, atheta, aphi, ax, ay, az)

  implicit none

  real, intent(in) :: theta, phi, ar, atheta, aphi
  real, intent(out) :: ax, ay, az

  real sin_theta, cos_theta, sin_phi, cos_phi
  real atr11, atr12, atr13
  real atr21, atr22, atr23
  real atr31, atr32, atr33

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
  atr33 =  0.0

  ax = atr11 * ar + atr12 * atheta + atr13 * aphi
  ay = atr21 * ar + atr22 * atheta + atr23 * aphi
  az = atr31 * ar + atr32 * atheta + atr33 * aphi

  return

end subroutine convert_spherical_vector_to_cart_4
