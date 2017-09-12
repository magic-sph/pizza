module char_manip

   use precision_mod

   implicit none

   private

   public :: capitalize, length_to_blank, dble2str

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
end module char_manip
