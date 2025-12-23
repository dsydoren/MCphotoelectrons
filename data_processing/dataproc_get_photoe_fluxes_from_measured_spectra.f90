!-------------------------------
!
program get_measured_photoe_fluxes_vst

  implicit none

  integer N_orbit_points

  logical exists
  character(1) buf

  real rdummy
  integer idummy

  integer N_ranges

  integer, parameter :: max_N_ranges = 20
  real energy_min_eV(max_N_ranges)
  real energy_max_eV(max_N_ranges)
  integer n

  integer i

  integer hh, mm, ss
  real ut_h

  character(13) spectr_filename

  integer N_spectr_points
  integer ios

  real, allocatable :: w_eV(:)
  real, allocatable :: flux(:)
  integer j

  real flux_avg(max_N_ranges)
  integer count
  real width_eV, step_eV

! functions

  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       integer length_of_string
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
     end function convert_int_to_txt_string
  end interface

  N_orbit_points = 15
  print '("N of orbit points is ",i2)', N_orbit_points

! get energy ranges from file

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
  read (9, *) rdummy !instr_fov_angle
  read (9, '(A1)') buf   ! # for plotting, length of line representing sensor's field of view axis (km)
  read (9, *) rdummy !fov_axis_length_km
  read (9, '(A1)') buf   ! # for electron distribution over energy and pitch angle define below
  read (9, '(A1)') buf   ! # number of energy bins
  read (9, *) idummy !max_eedf_N_w
  read (9, '(A1)') buf   ! # size of energy bin (eV)
  read (9, *) rdummy !dw_eV
  read (9, '(A1)') buf   ! # number of pitch angle bins (pitch angle varies between 0 and pi)
  read (9, *) idummy !max_eedf_N_theta
  read (9, '(A1)') buf   ! # save the electron distribution over energy and pitch angle into a file for each orbit point (1/0 = Yes/No)
  read (9, *) idummy !save_edfepa_flag

  read (9, '(A1)') buf   ! # number of energy-pitch angle ranges where particle fluxes will be calculated (<=20)
  read (9, *) N_ranges
  
  energy_min_eV =  1.0e6
  energy_max_eV = -1.0e6
  read (9, '(A1)') buf   ! # for each range give (in one line) minimal/ maximal energy (eV) 
  do n = 1, N_ranges    
     read (9, *) energy_min_eV(n), energy_max_eV(n)
  end do

  close (9, status = 'keep')

  print '("finished reading input_photoelectron_fluxes_igrf.dat")'

! get orbit points / spectra times from file

  inquire (file = 'spectr_times.dat', exist = exists)
  if (.not.exists) then
     print '("error : file spectr_times.dat not found")'
     stop
  end if

  open (11, file = 'measured_avg_fluxes_vs_t.dat') ! output file with average vs time

  open (9, file = 'spectr_times.dat')   ! spectra times

  do i = 1, N_orbit_points

! read time
     read (9, *) hh, mm, ss
     ut_h = real(hh) + real(mm*60 + ss)/3600.0

! prepare spectrum filename
     spectr_filename = 'spectr_NN.dat'
!                       ----x----I---
     spectr_filename(8:9) = convert_int_to_txt_string(i, 2)

     print '("processing file ",A13)', spectr_filename

     inquire (file = spectr_filename, exist=exists)
     if (.not.exists) then
        print '("file ",A13," not found, skip")', spectr_filename
        cycle
     end if

! find number of points in the spectrum file
     N_spectr_points = 0
     open (10, file = spectr_filename)
     do
        read (10,*,iostat=ios) rdummy, rdummy
        if (ios.ne.0) exit
        N_spectr_points = N_spectr_points+1
     end do
     close (10, status = 'keep')

! allocate arrays
     if (allocated(w_eV)) deallocate(w_eV)
     if (allocated(flux)) deallocate(flux)
     allocate(w_eV(1:N_spectr_points))
     allocate(flux(1:N_spectr_points))
     w_eV = 0.0
     flux = 0.0

! re-read the file
     open (10, file = spectr_filename)
     do j = 1, N_spectr_points
        read (10, *) w_eV(j), flux(j)    ! photoelectron flux in units of 1/cm^2/s/eV/sr
     end do
     close (10, status = 'keep')

! find average fluxes for each range
     flux_avg = 0.0
     do n = 1, N_ranges
        count = 0
        width_eV = 0.0
        do j = 1, N_spectr_points
           if (w_eV(j).lt.energy_min_eV(n)) cycle
           if (w_eV(j).gt.energy_max_eV(n)) cycle
           count = count+1
           if (j.eq.1) then
              step_eV = w_eV(j+1) - w_eV(j)
           else if (j.eq.N_spectr_points) then
              step_eV = w_eV(j) - w_eV(j-1)
           else
              step_eV = 0.5*(w_eV(j+1) - w_eV(j-1))
           end if
           width_eV = width_eV + step_eV
           flux_avg(n) = flux_avg(n) + step_eV * flux(j)
        end do
        if (count.gt.0) then
           flux_avg(n) = flux_avg(n) / width_eV
        end if
     end do

     print '("done")'

! save
     write (11, '(2x,f7.4,20(2x,e12.5))') ut_h, flux_avg(1:N_ranges)

  end do

  close (9, status = 'keep')  ! file with spectra times

  close (11, status = 'keep') ! file with average fluxes vs time

  print '("file measured_avg_fluxes_vs_t.dat is ready")'

! cleanup
  if (allocated(w_eV)) deallocate(w_eV)
  if (allocated(flux)) deallocate(flux)

end program get_measured_photoe_fluxes_vst


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

