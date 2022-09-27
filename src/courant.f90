module courant_mod
 
   use parallel_mod
   use precision_mod
   use outputs, only: n_log_file
   use namelists, only:  dt_fac, tag, l_3D, l_mag_3D, l_mag_LF, l_mag_alpha, &
       &                 l_cour_alf_damp, DyMagFac, BdiffFac!=1/pm=eta/nu, 
   use truncation, only: n_phi_max
   use truncation_3D, only: n_theta_max, n_phi_max_3D
   use blocking, only: nRstart3D, nRstop3D, nRstart, nRstop
   use radial_functions, only: delxr2, delxh2, delxr2_3D, delxh2_3D, &
       &                       beta, r_3D, or2_3D!, or1_3D
   !use horizontal, only: cost, sint, osint1
   use useful, only: logWrite
   use constants, only: half, one, two

   implicit none

   private

   integer :: file_handle

   public :: courant, courant_3D, dt_courant, initialize_courant, finalize_courant

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

      !-- Local  variables:
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
   subroutine courant_3D(n_r,dtrkc_3D,dthkc_3D,vr,vt,vp,br,bt,bp,ac,courfac,alffac)
      !
      !  courant condition check: calculates Courant  for Alfven waves                     
      !  advection lengths in radial direction dtrkc_3D                      
      !  and in horizontal direction dthkc_3D                                
      !  on the local radial level n_r                                   
      !
      !  for the effective velocity, the abs. sum of Alfven velocity is taken
      !
    
      !-- Input variable:
      integer,  intent(in) :: n_r             ! radial level

      real(cp), intent(in) :: vr(n_phi_max_3D,n_theta_max)  ! radial 3D velocity
      real(cp), intent(in) :: vt(n_phi_max_3D,n_theta_max)  ! longitudinal 3D velocity
      real(cp), intent(in) :: vp(n_phi_max_3D,n_theta_max)  ! azimuthal 3D velocity
      real(cp), intent(in) :: br(n_phi_max_3D,n_theta_max)  ! radial 3D magnetic field
      real(cp), intent(in) :: bt(n_phi_max_3D,n_theta_max)  ! longitudinal 3D magnetic field
      real(cp), intent(in) :: bp(n_phi_max_3D,n_theta_max)  ! azimuthal 3D magnetic field
      real(cp), intent(in) :: ac(n_phi_max_3D,n_theta_max)  ! Alpha-effect function (if any)
      real(cp), intent(in) :: courfac         ! Courant factor
      real(cp), intent(in) :: alffac          ! Alfven factor

      !-- Output:
      real(cp), intent(out) :: dtrkc_3D    ! Courant step (based on radial 3D advection)
                                             ! for the range of points covered
      real(cp), intent(out) :: dthkc_3D    ! Courant step based on horizontal 3D advection

      !-- Local  variables:
      integer :: n_theta
      integer :: n_phi
      real(cp) :: valr2max,valh2max!, valr2max2,valh2max2
      real(cp) :: valr,valr2
      real(cp) :: valri2,valhi2,valh2,valh2m
      real(cp) :: vflr2,vflh2
      real(cp) :: vapr2,vaph2
      real(cp) :: r2, o_r_2_r, o_r_4_r, r4
      real(cp) :: af2, cf2

      vapr2=0.0_cp
      vaph2=0.0_cp
      valr2max=0.0_cp
      valh2max=0.0_cp
      !valr2max2=0.0_cp
      !valh2max2=0.0_cp
      dtrkc_3D=1.0e10_cp
      dthkc_3D=1.0e10_cp
      af2=alffac*alffac
      cf2=courfac*courfac

      if ( l_cour_alf_damp ) then
         valri2=(half*(one+BdiffFac))**2./delxr2_3D(n_r)
         valhi2=(half*(one+BdiffFac))**2./delxh2_3D(n_r)
      else
         valri2=0.0_cp
         valhi2=0.0_cp
      end if

      r2 = r_3D(n_r)*r_3D(n_r)
      r4 = r2*r2
      o_r_2_r=or2_3D(n_r)
      o_r_4_r=o_r_2_r*o_r_2_r

      if ( l_mag_3D .and. l_mag_LF ) then

         !=Reminder:: Pizza variables == X * Magic variables:
         !-- ._r = r^2 * ._r ; ._t = r*sint * ._t ; ._p = r*sint * ._p
         do n_theta=1,n_theta_max
            do n_phi=1,n_phi_max_3D
               valr =br(n_phi,n_theta)*br(n_phi,n_theta) * &
               &     DyMagFac*r4
               valr2=valr*valr/(valr+valri2)
               vflr2=vr(n_phi,n_theta)*vr(n_phi,n_theta)*r4
               if ( l_mag_alpha ) vapr2=ac(n_phi,n_theta)*ac(n_phi,n_theta)*r4
               valr2max=max(valr2max,o_r_4_r*(af2*valr2 + cf2*vflr2 + cf2*vapr2))!))!
               !valr2max2=max(valr2max2,o_r_4_r*(cf2*vflr2))

               valh2= ( bt(n_phi,n_theta)*bt(n_phi,n_theta) +  &
               &        bp(n_phi,n_theta)*bp(n_phi,n_theta) )* &
               &        DyMagFac*r2
               if ( l_mag_alpha ) valh2=valh2 + &
               &  one*ac(n_phi,n_theta)*ac(n_phi,n_theta) *   &
               &    ( bt(n_phi,n_theta)*bt(n_phi,n_theta) +   &
               &      bp(n_phi,n_theta)*bp(n_phi,n_theta)  )* &
               &  DyMagFac*r2
               valh2m=valh2*valh2/(valh2+valhi2)
               vflh2= ( vt(n_phi,n_theta)*vt(n_phi,n_theta) +  &
               &        vp(n_phi,n_theta)*vp(n_phi,n_theta) )* &
               &        r2
               if ( l_mag_alpha ) vaph2= two*ac(n_phi,n_theta)*ac(n_phi,n_theta)*r2

               valh2max=max(valh2max,o_r_2_r*(af2*valh2m + cf2*vflh2 + cf2*vapr2))!))!
               !valh2max2=max(valh2max2,o_r_2_r*(cf2*vflh2))
            end do
         end do

         !print*, 'l_cour_alf_damp?:: ', l_cour_alf_damp
         !print*, 'n_r; r_3D(n_r)', n_r, '; ', r_3D(n_r)
         !print*, 'DyMagFac; valri2, valhi2', DyMagFac, '; ', valri2, '; ', valhi2
         !print*, 'br; btheta; bphi', r2*br(2,2), '; ', bt(2,2)/(or1_3D(n_r)*osint1(2)), '; ', bp(2,2)/(or1_3D(n_r)*osint1(2))
         !print*, 'valr2max', valr2 
         !print*, 'valh2max', valh2
         !print*, ''
         !print*, 'vr; vtheta; vphi', r2*vr(2,2), '; ', vt(2,2)/(or1_3D(n_r)*osint1(2)), '; ', vp(2,2)/(or1_3D(n_r)*osint1(2))
         !print*, 'vaflr2', vflr2
         !print*, 'vaflh2', vflh2
         !print*, ''
         !print*, "COUR/ALF FAC==", 'valr2max, valh2max', valr2max, valh2max
         !print*, "COUR=", 'vflr2max, vflh2max', valr2max2, valh2max2
         !print*, ''
         !if (n_r == 3) stop

         if ( valr2max /= 0.0_cp ) dtrkc_3D = min(dtrkc_3D,sqrt(delxr2_3D(n_r)/valr2max))
         if ( valh2max /= 0.0_cp ) dthkc_3D = min(dthkc_3D,sqrt(delxh2_3D(n_r)/valh2max))

      else

         dtrkc_3D=min(dtrkc_3D,sqrt(delxr2_3D(n_r)/cf2))
         dthkc_3D=min(dthkc_3D,sqrt(delxr2_3D(n_r)/cf2))

      end if   ! Magnetic field ?

      !print*, 'dtrkc_R', dtrkc_3D
      !print*, 'dthkc_H', dthkc_3D
      !print*, 'COURant LAYER OVER-----------------------------------------------------------------------------'
      !if (n_r == 3) stop

   end subroutine courant_3D
!------------------------------------------------------------------------------
   subroutine dt_courant(dt_r,dt_h,l_new_dt,dt,dt_new,dtMax,dtrkc,dthkc, &
              &          dtrkc_3D,dthkc_3D,time)
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
      real(cp), intent(in) :: dtrkc_3D(nRstart3D:nRstop3D)
      real(cp), intent(in) :: dthkc_3D(nRstart3D:nRstop3D)
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

      if ( l_3D ) then
         do n_r=nRstart3D,nRstop3D
            dt_r=min(dtrkc_3D(n_r),dt_r)
            dt_h=min(dthkc_3D(n_r),dt_h)
         end do
      end if

      call MPI_Allreduce(MPI_IN_PLACE,dt_r,1,MPI_DEF_REAL, &
           &             MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_Allreduce(MPI_IN_PLACE,dt_h,1,MPI_DEF_REAL, &
           &             MPI_MIN,MPI_COMM_WORLD,ierr)

      dt_rh=min(dt_r,dt_h)
      dt_2 =min(half*(one/dt_fac+one)*dt_rh,dtMax)

      !print*, 'dt_R', 'dt_H', dt_r, ';' , dt_h
      !print*, 'dt_RH', 'dt_2', dt_rh, ';' , dt_2

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

      !stop

   end subroutine dt_courant
!-----------------------------------------------------------------------
end module courant_mod
