module courant_mod
 
   use parallel_mod
   use precision_mod
   use outputs, only: n_log_file
   use namelists, only:  dt_fac, tag
   use truncation, only: n_phi_max, n_r_max
   use blocking, only: nRstart, nRstop
   use radial_functions, only: delxr2, delxh2, beta
   use useful, only: logWrite
   use constants, only: half, one, two

   implicit none

   private

   integer :: file_handle

   public :: courant, dt_courant, initialize_courant, finalize_courant

contains

   subroutine initialize_courant(time, dt)

      !-- Input variables
      real(cp), intent(in) :: time ! time
      real(cp), intent(in) :: dt   ! timestep

      if ( rank == 0 ) then
         open(newunit=file_handle, file='timestep.'//tag, status='new')
         write(file_handle, '(1p, es20.12, es16.8)')  time, dt
      end if

   end subroutine initialize_courant
!------------------------------------------------------------------------------
   subroutine finalize_courant()

      if ( rank == 0 ) close(file_handle)

   end subroutine finalize_courant
!------------------------------------------------------------------------------
   subroutine courant(n_r,dtrkc,dthkc,ur,up,courfac)
      !
      !  courant condition check: calculates Courant                      
      !  advection lengths in radial direction dtrkc                      
      !  and in horizontal direction dthkc                                
      !  on the local radial level n_r                                   
      !
    
      !-- Input variable:
      integer,  intent(in) :: n_r             ! radial level
      real(cp), intent(in) :: ur(n_phi_max)   ! radial velocity
      real(cp), intent(in) :: up(n_phi_max)   ! azimuthal velocity
      real(cp), intent(in) :: courfac         ! Courant factor

    
      !-- Output:
      real(cp), intent(out) :: dtrkc    ! Courant step (based on radial advection)
                                        ! for the range of points covered
      real(cp), intent(out) :: dthkc    ! Courant step based on horizontal advection
    
      integer :: n_phi
      real(cp) :: vr2max,vh2max
      real(cp) :: vflr2,vflh2
      real(cp) :: cf2
    
      vr2max=0.0_cp
      vh2max=0.0_cp
      dtrkc =1.0e10_cp
      dthkc =1.0e10_cp
      cf2=courfac*courfac
    
      do n_phi=1,n_phi_max
    
         vflr2 =ur(n_phi)*ur(n_phi)
         vr2max=max(vr2max,cf2*vflr2)
    
         vflh2 =up(n_phi)*up(n_phi) 
         vh2max=max(vh2max,cf2*vflh2)
    
      end do
    
      if ( vr2max /= 0.0_cp ) dtrkc=min(dtrkc,sqrt(delxr2(n_r)/vr2max))
      if ( vh2max /= 0.0_cp ) dthkc=min(dthkc,sqrt(delxh2(n_r)/vh2max))
    
   end subroutine courant
!------------------------------------------------------------------------------
   subroutine dt_courant(dt_r,dt_h,l_new_dt,dt,dt_new,dtMax,dtrkc,dthkc,time)
      !
      ! Check if Courant criterion based on combined
      ! fluid velocity is satisfied
      ! Returns new value of time step dtnew
      !
      ! dtr,dth: (output) radial/horizontal Courant time step
      ! n_time_step: (input) time step number
      ! l_new_dt: (output) flag indicating that time step is changed (=1) or not (=0)
      ! dt: (input) old time step
      ! dtnew: (output) new time step
      ! dtMin: (input) lower limit for time step (termination if dtnew < dtMin)
      ! dtMax: (input) upper limit for time step
      ! dtrkc: (input) radial Courant time step as function of radial level
      ! dthkc: (input) horizontal Courant time step as function of radial level
      !

      !-- Input variables:
      real(cp), intent(in) :: dt
      real(cp), intent(in) :: dtMax
      real(cp), intent(in) :: dtrkc(nRstart:nRstop)
      real(cp), intent(in) :: dthkc(nRstart:nRstop)
      real(cp), intent(in) :: time
    
      !-- Output variables:
      logical,  intent(out) :: l_new_dt
      real(cp), intent(out) :: dt_new
      real(cp), intent(out) :: dt_r,dt_h
    
      !-- Local:
      integer :: n_r
      real(cp) :: dt_rh,dt_2
    
      character(len=200) :: message
    
      dt_r  =1000.0_cp*dtMax
      dt_h  =dt_r
      do n_r=nRstart,nRstop
         dt_r=min(dtrkc(n_r),dt_r)
         dt_h=min(dthkc(n_r),dt_h)
      end do

      call MPI_Allreduce(MPI_IN_PLACE,dt_r,1,MPI_DEF_REAL, &
           &             MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(MPI_IN_PLACE,dt_h,1,MPI_DEF_REAL, &
           &             MPI_MIN,MPI_COMM_WORLD,ierr)

      dt_rh=min(dt_r,dt_h)
      dt_2 =min(half*(one/dt_fac+one)*dt_rh,dtMax)

      if ( dt > dtMax ) then ! Timestep larger than dtMax from Namelist
    
         l_new_dt=.true.
         dt_new=dtMax
         write(message,'(1P," ! COURANT: dt=dtMax =",ES12.4,A)') dtMax,&
              &" ! Think about changing dtMax !"
         call logWrite(message,n_log_file)
    
      else if ( dt > dt_rh ) then ! Timestep decrease
    
         l_new_dt=.true.
         dt_new  =dt_2
         write(message,'(1P," ! COURANT: dt=",ES11.4," > dt_r=",ES12.4, &
              &       " and dt_h=",ES12.4)') dt,dt_r,dt_h
         call logWrite(message,n_log_file)
         if ( rank == 0 ) then
            write(file_handle, '(1p, es20.12, es16.8)')  time, dt_new
         end if
    
      else if ( dt_fac*dt < dt_rh .and. dt < dtMax ) then ! Timestep increase
    
         l_new_dt=.true.
         dt_new=dt_2
         write(message,'(" ! COURANT: ",F4.1,1P,"*dt=",ES11.4, &
              &     " < dt_r=",ES12.4," and dt_h=",ES12.4)') &
              &     dt_fac,dt_fac*dt,dt_r,dt_h
         call logWrite(message,n_log_file)
         if ( rank == 0 ) then
            write(file_handle, '(1p, es20.12, es16.8)')  time, dt_new
         end if

      else 

         l_new_dt=.false.
         dt_new=dt
    
      end if

   end subroutine dt_courant
!-----------------------------------------------------------------------
end module courant_mod
