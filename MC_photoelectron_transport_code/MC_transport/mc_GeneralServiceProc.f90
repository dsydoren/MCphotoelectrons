
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

!============================================================
!
subroutine filter(m_apply, n1, n2, x)

  implicit none

  integer, intent(in) :: m_apply      ! how many times the filter should be applied
  integer, intent(in) :: n1, n2       ! array boundaries
  
  real(8), intent(inout) :: x(n1:n2)  ! data to be filtered

  integer n, m
  real(8) temp

! cycle over filter application
  do m = 1, m_apply

! save the leftmost value
     temp = x(n1)

! sweep right
     do n = n1, n2-1
        x(n) = x(n) + x(n+1)
     end do

! sweep left
     do n = n2-1, n1+1, -1
        x(n) = (x(n) + x(n-1)) * 0.25_8
     end do

! restore the leftmost value
     x(n1) = temp

  end do

end subroutine filter
