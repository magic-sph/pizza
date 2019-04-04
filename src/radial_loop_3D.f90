module rloop_3D

   use precision_mod
   use parallel_mod
   use mem_alloc, only: bytes_allocated
   use nonlinear_lm_mod, only: nonlinear_lm_t
   use grid_space_arrays_mod, only: grid_space_arrays_t
   use blocking, only: nRstart3D, nRstop3D
   use truncation_3D, only: lm_max, lmP_max
   use fields_3D, only: ur_3D_Rloc, ut_3D_Rloc, up_3D_Rloc

   implicit none

   private

   type, abstract, public :: rloop_3D_type
      integer :: n_r
      type(grid_space_arrays_t) :: gsa
      type(nonlinear_lm_t) :: nl_lm
   contains
      procedure(empty_if), deferred :: initialize
      procedure :: initialize => initialize_radial_loop_3D
      procedure :: finalize => finalize_radial_loop_3D
      procedure :: radial_loop_3D => radial_loop_3D
      procedure :: transform_to_grid_space_shtns => transform_to_grid_space_shtns
      procedure :: transform_to_lm_space_shtns => transform_to_lm_space_shtns
   end type rloop_3D_type

   public :: radial_loop_3D, initialize_radial_loop_3D, finalize_radial_loop_3D

   interface

      subroutine empty_if(this)
         import
         class(rloop_3D_type) :: this
      end subroutine empty_if
   end interface

contains
!------------------------------------------------------------------------------
   subroutine initialize_radial_loop_3D

      allocate( rloop_3D_type :: this )
      !call this%initialize()

      call this%gsa%initialize()
      call this%nl_lm%initialize(lmP_max)

   end subroutine initialize_radial_loop_3D
!------------------------------------------------------------------------------
   subroutine finalize_radial_loop_3D(this)

      call this%gsa%finalize()
      call this%nl_lm%finalize()

      !call this%finalize()
      deallocate( this )

   end subroutine finalize_radial_loop_3D
!------------------------------------------------------------------------------
   subroutine radial_loop_3D( temp_Rloc, dtempdt_Rloc, dVrT_Rloc)

      !-- Input variables
      complex(cp), intent(in) :: temp_Rloc(lm_max, nRstart3D:nRstop3D)

      !-- Output variables
      complex(cp), intent(out) :: dtempdt_Rloc(lm_max, nRstart3D:nRstop3D)
      complex(cp), intent(out) :: dVrT_Rloc(lm_max, nRstart3D:nRstop3D)

      !-- Local variables
      complex(cp) :: us_fluct
      integer :: n_r

      do n_r=nRstart3D,nRstop3D

         call transform_to_grid_space_shtns(this%gsa, time)

         call gsa%get_nl(ur_3D_Rloc(:,:,n_r), ut_3D_Rloc(:,:,n_r), &
              &          up_3D_Rloc(:,:,n_r), this%n_r)

         call transform_to_lm_space_shtns(this%gsa, this%nl_lm)

         !-- Partial calculation of time derivatives (horizontal parts):
         !   input flm...  is in (l,m) space at radial grid points n_r !
         !   get_td finally calculates the d*dt terms needed for the
         !   time step performed in 's_LMLoop.f' . This should be distributed
         !   over the different models that 's_LMLoop.f' parallelizes over.
         call get_td(this%n_r, dVTrLM, dtempdt_Rloc)

      end do

   end subroutine radial_loop_3D
!-------------------------------------------------------------------------------
   subroutine transform_to_grid_space_shtns(this, temp, gsa)

      class(rloop_3D_type) :: this
      complex(cp), intent(in) :: temp(lm_max)
      type(grid_space_arrays_t) :: gsa

      call scal_to_spat(temp, gsa%Tc)

   end subroutine transform_to_grid_space_shtns
!-------------------------------------------------------------------------------
   subroutine transform_to_lm_space_shtns(this, gsa, nl_lm)

      class(rloop_3D_type) :: this
      type(grid_space_arrays_t) :: gsa
      type(nonlinear_lm_t) :: nl_lm

      call shtns_load_cfg(1)

      call spat_to_SH(gsa%VTr, nl_lm%VTrLM)
      call spat_to_SH(gsa%VTt, nl_lm%VTtLM)
      call spat_to_SH(gsa%VTp, nl_lm%VTpLM)

      call shtns_load_cfg(0)

   end subroutine transform_to_lm_space_shtns
!------------------------------------------------------------------------------
end module rloop_3D
