module pre_calculations

   use truncation, only: n_r_max, write_truncation_info
   use truncation_3D, only: write_truncation_3D_info
   use constants, only: one, four, third, pi, vol_oc, surf, two, vol_otc
   use outputs, only: n_log_file
   use radial_functions, only: r, rscheme, radial, radial_3D
   use namelists, only: l_newmap, tag, dtMax, dtMin, pr, l_non_rot, ek, &
       &                r_cmb, r_icb, l_3D
   use horizontal, only: mfunctions, spherical_functions
   use z_functions
   use precision_mod
   use parallel_mod

   implicit none

   private
   
   public :: preCalc

contains

   subroutine preCalc(zinterp)

      !-- Input variable
      type(zfunc_type), intent(inout) :: zinterp

      !-- Local variables
      integer :: file_handle, n_r
      character(len=100) :: file_name

      dtMin = dtMax/1e6_cp

      !-- Set r-dependent arrays
      call radial()

      !-- Set r-dependent arrays for 3-D calculations
      if ( l_3D ) call radial_3D()

      if ( ( l_newmap ) .and. (rank == 0) ) then
         file_name='rNM.'//tag
         open(newunit=file_handle, file=file_name, status='unknown')
         do n_r=1,n_r_max
            write(file_handle,'(I5,3ES16.8)') n_r, r(n_r),     &
            &                                rscheme%drx(n_r), &
            &                                rscheme%ddrx(n_r)
         end do
         close(file_handle)
      end if

      !-- Set m-dependent arrays
      call mfunctions()

      !-- theta and l functions
      if ( l_3D ) then
         call spherical_functions()
         call zinterp%fill_mat()
      end if

      !-- Compute some constants
      vol_oc =four*third*pi*(r_cmb**3-r_icb**3)
      vol_otc=four*third*pi*(r_cmb**2-r_icb**2)**(1.5_cp)
      surf   =pi*(r_cmb**2-r_icb**2)

      !-- Write some informations
      if ( rank == 0 ) then
         call write_info(6)
         call write_info(n_log_file)
      end if

   end subroutine preCalc
!---------------------------------------------------------------------------
   subroutine write_info(n_out)

      integer, intent(in) :: n_out

      write(n_out,*) ''
      write(n_out, '('' ! Spherical shell volume  :'',es14.6)') vol_oc
      write(n_out, '('' ! Volume outside tan. cyl.:'',es14.6)') vol_otc
      write(n_out, '('' ! Annulus surface         :'',es14.6)') surf

      write(n_out,*) ''
      write(n_out,*) '! MPI ranks:'
      write(n_out,'(''  n_procs      ='',i6)') n_procs

      call write_truncation_info(n_out)

      if ( l_3D ) call write_truncation_3D_info(n_out)

      write(n_out,*) ''


   end subroutine write_info
!---------------------------------------------------------------------------

end module pre_calculations
