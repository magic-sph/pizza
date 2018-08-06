module select_time_scheme

   use dirk_schemes, only: type_dirk
   use multistep_schemes, only: type_multistep
   use time_schemes, only: type_tscheme
   use char_manip, only: capitalize
   implicit none

   private

   public :: select_tscheme

contains

   subroutine select_tscheme(scheme_name, tscheme)

      class(type_tscheme), pointer :: tscheme
      character(len=72), intent(inout) :: scheme_name

      call capitalize(scheme_name)

      if ( (index(scheme_name, 'ARS222') /= 0) .or. &
      &    (index(scheme_name, 'ARS443') /= 0) .or. &
      &    (index(scheme_name, 'BPR353') /= 0) .or. &
      &    (index(scheme_name, 'PC2') /= 0)    .or. &
      &    (index(scheme_name, 'CNRKW3') /= 0) .or. &
      &    (index(scheme_name, 'LZ453') /= 0)  .or. &
      &    (index(scheme_name, 'LZ674') /= 0)  .or. &
      &    (index(scheme_name, 'LZ453A') /= 0) .or. &
      &    (index(scheme_name, 'LZ453B') /= 0) .or. &
      &    (index(scheme_name, 'LZ453C') /= 0) .or. &
      &    (index(scheme_name, 'LZ574') /= 0)  .or. &
      &    (index(scheme_name, 'LRR232') /= 0) .or. &
      &    (index(scheme_name, 'FW353') /= 0)  .or. &
      &    (index(scheme_name, 'LZ232') /= 0) ) then
         allocate ( type_dirk :: tscheme )
      else
         allocate ( type_multistep :: tscheme )
      end if

   end subroutine select_tscheme

end module select_time_scheme
