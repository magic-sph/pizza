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
      &    (index(scheme_name, 'ARS232') /= 0) .or. &
      &    (index(scheme_name, 'ARS233') /= 0) .or. &
      &    (index(scheme_name, 'ARS443') /= 0) .or. &
      &    (index(scheme_name, 'ARS343') /= 0) .or. &
      &    (index(scheme_name, 'MARS343') /= 0).or. &
      &    (index(scheme_name, 'KC343') /= 0)  .or. &
      &    (index(scheme_name, 'KC564') /= 0)  .or. &
      &    (index(scheme_name, 'KC785') /= 0)  .or. &
      &    (index(scheme_name, 'DBM453') /= 0) .or. &
      &    (index(scheme_name, 'BPR353') /= 0) .or. &
      &    (index(scheme_name, 'PC2') /= 0)    .or. &
      &    (index(scheme_name, 'LZ453') /= 0)  .or. &
      &    (index(scheme_name, 'CK232') /= 0)  .or. &
      &    (index(scheme_name, 'LZ232') /= 0) ) then
         allocate ( type_dirk :: tscheme )
      else
         allocate ( type_multistep :: tscheme )
      end if

   end subroutine select_tscheme

end module select_time_scheme
