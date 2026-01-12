module char_manip

   use precision_mod

   implicit none

   private

   public :: capitalize, length_to_blank, dble2str, write_long_string

contains

   integer function length_to_blank(string)
      !
      !   determines number of characters before first blank in string
      !

      !-- Input variable
      character(len=*), intent(in) :: string

      !-- Local variable
      integer :: i

      length_to_blank=0
      do i=1,len(string)
         if( string(i:i) == ' ' ) then
            length_to_blank=i-1
            exit
         end if
      end do

   end function length_to_blank
!------------------------------------------------------------------------------
   subroutine capitalize(string)
      !
      !   Convert lower-case letters into capital letters
      !

      !-- Input variables
      character(len=*), intent(inout) :: string

      !-- Local variables
      character(len=26), parameter :: LOWER_CASE='abcdefghijklmnopqrstuvwxyz'
      character(len=26), parameter :: UPPER_CASE='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

      integer :: i, n

      do i = 1, len(string)
      ! -- Find location of letter in lower case constant string
         n = index( LOWER_CASE, string( i:i ) )
      ! -- If current substring is a lower case letter,
      !    make it upper case
         if ( n /= 0 ) string( i:i ) = UPPER_CASE( n:n )
      end do

   end subroutine capitalize
!------------------------------------------------------------------------------
   subroutine dble2str(num, str)
      !
      !  converts a dble number num into a character str
      !

      !-- Input variable
      real(cp), intent(in) :: num

      !-- Output variable
      character(len=*), intent(out) :: str

      !-- Local variables
      character(len=72) :: work
      integer :: i

      write(work, '(F20.12)') num
      write(str, '(A8)') trim(adjustl(work))
      i = index(str,'.')
      str(i:i)='_'

   end subroutine dble2str
!------------------------------------------------------------------------------
   subroutine write_long_string(prefix, long_string, out_unit)
      !
      ! This subroutine is used to split a long string (with len(str)>85) into
      ! a multi-lines string for a cleaner printout.
      !
   
      !-- Input variables
      character(len=*), intent(in) :: prefix
      character(len=*), intent(in) :: long_string
      integer,          intent(in) :: out_unit

      !-- Local variables
      integer :: len_long_str, len_prefix, istart, iend, len_per_line
      integer, parameter :: line_width=85
      character(len=8) :: fmtstr
   
      len_long_str = len(long_string)
      len_prefix = len(prefix)
      len_per_line = line_width-len_prefix
      write(fmtstr,'(A,i2)') 'A', len_prefix-2
      
      if ( len_prefix+len_long_str <= 85 ) then
         write(out_unit, '(A,A)') prefix, long_string
      else
         istart = 1; iend = len_per_line
         write(out_unit, '(A,A)') prefix, long_string(istart:iend)
         istart = iend+1
         do while (iend < len_long_str-len_per_line)
            iend = istart+len_per_line
            write(out_unit, '(A,'//trim(fmtstr)//',A)') ' !','',long_string(istart:iend)
            istart = iend+1
         end do
         write(out_unit, '(A,'//trim(fmtstr)//',A)') ' !','',long_string(istart:)
      end if

   end subroutine write_long_string
!------------------------------------------------------------------------------
end module char_manip
