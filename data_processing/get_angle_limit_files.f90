program getlimitfile

  implicit none

  integer op

  real thetamin, thetamax, dtheta

  real, parameter :: deg2rad = 3.1415926536/180.0
  real, parameter :: rad2deg = 180.0/3.1415926536

  integer, parameter :: Ntheta=20

  character(37) filename2   ! flux_limits_vperp_vpar_rNN_op_PPP.dat
!                             ----x----I----x----I----x----I----x--

  character(34) filename3   ! flux_limits_w_theta_rNN_op_PPP.dat
!                             ----x----I----x----I----x----I----
 
  REAL(8), PARAMETER :: e_Cl     = 1.602189d-19                          ! Elementary charge [Cl]
  REAL(8), PARAMETER :: m_e_kg   = 9.109534d-31                          ! Electron mass [kg]

  integer r
  real wmineV, wmaxeV, vmin, vmax, theta
  integer i

  interface
     function convert_int_to_txt_string(int_number, length_of_string)
       integer length_of_string
       character*(length_of_string) convert_int_to_txt_string
       integer int_number
     end function convert_int_to_txt_string
  end interface

  do op = 1, 15

     select case (op)
        case (1)
           thetamin = 146.22
           thetamax = 166.22
        case (8)
           thetamin = 151.87
           thetamax = 171.87
        case default
           cycle
     end select

! convert degrees to radians
     thetamin = thetamin * deg2rad
     thetamax = thetamax * deg2rad
     dtheta=(thetamax-thetamin)/(Ntheta-1)

     do r = 1, 10

        wmineV = 5.0+8.0*(r-1)
        wmaxeV = 5.0+8.0*r

        filename2 = 'flux_limits_vperp_vpar_rNN_op_PPP.dat'
!                    ----x----I----x----I----x----I----x--
        filename2(25:26) = convert_int_to_txt_string(r, 2)
        filename2(31:33) = convert_int_to_txt_string(op, 3)

        vmin=sqrt(2.0*e_Cl*wmineV/m_e_kg)
        vmax=sqrt(2.0*e_Cl*wmaxeV/m_e_kg)

        open (9, file = filename2)
        do i = 0, Ntheta-1
           theta=thetamin+i*dtheta
           write (9, '(2x,e12.5,2x,e12.5)') vmin*sin(theta), vmin*cos(theta)
        end do
        do i=Ntheta-1, 0, -1
           theta=thetamin+i*dtheta
           write (9, '(2x,e12.5,2x,e12.5)') vmax*sin(theta), vmax*cos(theta)
        end do
        write (9, '(2x,e12.5,2x,e12.5)') vmin*sin(thetamin), vmin*cos(thetamin)
        close (9, status = 'keep')
        
        print '("file ",A37," is ready")', filename2

        filename3 = 'flux_limits_w_theta_rNN_op_PPP.dat'
!                    ----x----I----x----I----x----I----
        filename3(22:23) = convert_int_to_txt_string(r, 2)
        filename3(28:30) = convert_int_to_txt_string(op, 3)

        open (9, file = filename3)

        write (9, '(2x,f7.2,2x,f7.2)') wmineV, thetamin * rad2deg
        write (9, '(2x,f7.2,2x,f7.2)') wmineV, thetamax * rad2deg
        write (9, '(2x,f7.2,2x,f7.2)') wmaxeV, thetamax * rad2deg
        write (9, '(2x,f7.2,2x,f7.2)') wmaxeV, thetamin * rad2deg
        write (9, '(2x,f7.2,2x,f7.2)') wmineV, thetamin * rad2deg

        print '("file ",A34," is ready")', filename3

     end do   !###   do r = 1, 10

  end do   !###   do op = 1, 15

end program getlimitfile

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
