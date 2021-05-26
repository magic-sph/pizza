module mloop_mod

   use time_schemes, only: type_tscheme
   use time_array, only: type_tarray
   use precision_mod
   use truncation, only: n_r_max
   use vp_balance, only: vp_bal_type
   use vort_balance, only: vort_bal_type
   use blocking, only: nMstart, nMstop
   use namelists, only: l_direct_solve, l_cheb_coll, l_chem, l_heat
   use update_temp_coll, only: update_temp_co, finish_exp_temp_coll, assemble_temp_coll
   use update_xi_coll, only: update_xi_co, finish_exp_xi_coll, assemble_xi_coll
   use update_temp_integ, only: update_temp_int, finish_exp_temp_int, &
       &                        assemble_temp_int
   use update_xi_integ, only: update_xi_int, finish_exp_xi_int, assemble_xi_int
   use update_psi_integ_smat, only: update_psi_int_smat, finish_exp_psi_int_smat, &
       &                            assemble_psi_int_smat
   use update_psi_integ_dmat, only: update_psi_int_dmat, finish_exp_psi_int_dmat
   use update_psi_coll_dmat, only: update_om_coll_dmat, finish_exp_psi_coll_dmat, &
       &                           assemble_psi_coll_dmat
   use update_psi_coll_smat, only: update_om_coll_smat, finish_exp_psi_coll_smat, &
       &                           assemble_psi_coll_smat
   use timers_mod, only: timers_type
   use useful, only: abortRun

   implicit none

   private

   public :: mloop, finish_explicit_assembly, assemble_stage

contains 

   subroutine mloop(temp_hat_Mloc, xi_hat_Mloc, temp_Mloc, dtemp_Mloc, xi_Mloc,&
              &     dxi_Mloc, psi_hat_Mloc, psi_Mloc, om_Mloc,  dom_Mloc,      &
              &     us_Mloc, up_Mloc, dTdt, dxidt, dpsidt, vp_bal, vort_bal,   &
              &     tscheme, lMat, l_log_next, timers)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      logical,             intent(in) :: lMat
      logical,             intent(in) :: l_log_next

      !-- Output variables
      complex(cp),         intent(out) :: temp_hat_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: xi_hat_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: temp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: dtemp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: xi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: dxi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: psi_hat_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: psi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: dom_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(inout) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(inout) :: up_Mloc(nMstart:nMstop,n_r_max)
      type(vp_bal_type),   intent(inout) :: vp_bal
      type(vort_bal_type), intent(inout) :: vort_bal
      type(type_tarray),   intent(inout) :: dpsidt
      type(type_tarray),   intent(inout) :: dTdt
      type(type_tarray),   intent(inout) :: dxidt
      type(timers_type),   intent(inout) :: timers

     if ( l_cheb_coll ) then
         if ( l_heat ) call update_temp_co(temp_Mloc, dtemp_Mloc, dTdt, tscheme, &
                            &              lMat, l_log_next)
         if ( l_chem ) call update_xi_co(xi_Mloc, dxi_Mloc, dxidt, tscheme, lMat, &
                            &            l_log_next)
         if ( l_direct_solve ) then
            call update_om_coll_smat(psi_Mloc, om_Mloc, dom_Mloc, us_Mloc,    &
                 &                   up_Mloc, temp_Mloc, xi_Mloc, dpsidt,     &
                 &                   vp_bal, vort_bal, tscheme, lMat, timers)
         else
            call update_om_coll_dmat(psi_Mloc, om_Mloc, dom_Mloc, us_Mloc,    &
                 &                   up_Mloc, temp_Mloc, xi_Mloc, dpsidt,     &
                 &                   vp_bal, vort_bal, tscheme, lMat, timers)
         end if
      else
         if ( l_heat ) call update_temp_int(temp_hat_Mloc, temp_Mloc, dtemp_Mloc, dTdt, &
                            &               tscheme, lMat, l_log_next)
         if ( l_chem ) call update_xi_int(xi_hat_Mloc, xi_Mloc, dxi_Mloc, dxidt,   &
                            &             tscheme, lMat, l_log_next)
         if ( l_direct_solve ) then
            call update_psi_int_smat(psi_hat_Mloc, psi_Mloc, om_Mloc, us_Mloc,    &
                 &                   up_Mloc, temp_Mloc, xi_Mloc, dpsidt, vp_bal, &
                 &                   vort_bal, tscheme, lMat, timers)
         else
            call update_psi_int_dmat(psi_Mloc, om_Mloc, us_Mloc, up_Mloc,     &
                 &                   temp_Mloc, xi_Mloc, dpsidt, vp_bal,      &
                 &                   vort_bal, tscheme, lMat, timers)
         end if
      end if

   end subroutine mloop
!------------------------------------------------------------------------------
   subroutine assemble_stage(temp_hat_Mloc, xi_hat_Mloc, temp_Mloc, dtemp_Mloc, &
              &              xi_Mloc, dxi_Mloc, psi_hat_Mloc, psi_Mloc,         &
              &              us_Mloc, up_Mloc, om_Mloc, dTdt, dxidt, dpsidt,    &
              &              tscheme, vp_bal, vort_bal, l_log_next, timers)

      !-- Input variables
      type(type_tarray),   intent(inout) :: dpsidt
      type(type_tarray),   intent(inout) :: dTdt
      type(type_tarray),   intent(inout) :: dxidt
      logical,             intent(in) :: l_log_next
      class(type_tscheme), intent(in) :: tscheme

      !-- Output variables
      complex(cp),         intent(out) :: temp_hat_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: xi_hat_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: psi_hat_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: temp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: dtemp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: xi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: dxi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(out) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(inout) :: psi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(inout) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(inout) :: up_Mloc(nMstart:nMstop,n_r_max)
      type(vp_bal_type),   intent(inout) :: vp_bal
      type(vort_bal_type), intent(inout) :: vort_bal
      type(timers_type),   intent(inout) :: timers

      if ( l_cheb_coll ) then
         if ( l_heat ) call assemble_temp_coll(temp_Mloc, dtemp_Mloc, dTdt, &
                            &                  tscheme, l_log_next)
         if ( l_chem ) call assemble_xi_coll(xi_Mloc, dxi_Mloc, dxidt, tscheme, &
                            &                l_log_next)

         if ( l_direct_solve ) then
            call assemble_psi_coll_smat(psi_Mloc, us_Mloc, up_Mloc, om_Mloc, temp_Mloc, &
                 &                      xi_Mloc, dpsidt, tscheme, vp_bal, vort_bal)
         else
            call assemble_psi_coll_dmat(psi_Mloc, us_Mloc, up_Mloc, om_Mloc, temp_Mloc, &
                 &                      xi_Mloc, dpsidt, tscheme, vp_bal, vort_bal)
         end if
      else
         if ( l_heat ) call assemble_temp_int(temp_hat_Mloc, temp_Mloc, dtemp_Mloc, &
                            &                 dTdt, tscheme, l_log_next)
         if ( l_chem) call assemble_xi_int(xi_hat_Mloc, xi_Mloc, dxi_Mloc, dxidt, &
                           &               tscheme, l_log_next)

         if ( l_direct_solve ) then
            call assemble_psi_int_smat(psi_hat_Mloc, psi_Mloc, om_Mloc, us_Mloc,    &
              &                        up_Mloc, temp_Mloc, xi_Mloc, dpsidt, vp_bal, &
              &                        vort_bal, tscheme, timers)
         else
            call abortRun('Assembly stage with integration method not implemented yet!')
         end if
      end if

   end subroutine assemble_stage
!------------------------------------------------------------------------------
   subroutine finish_explicit_assembly(temp_Mloc, xi_Mloc, psi_Mloc, us_Mloc,    &
              &                        up_Mloc, om_Mloc, dVsT_Mloc, dVsXi_Mloc,  &
              &                        dVsOm_Mloc, dTdt, dxidt, dpsidt, tscheme, &
              &                        vp_bal, vort_bal)

      !-- Input variables
      class(type_tscheme), intent(in) :: tscheme
      complex(cp),         intent(in) :: temp_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: xi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: psi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: us_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: up_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(in) :: om_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(inout) :: dVsT_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(inout) :: dVsXi_Mloc(nMstart:nMstop,n_r_max)
      complex(cp),         intent(inout) :: dVsOm_Mloc(nMstart:nMstop,n_r_max)

      !-- Output variables
      type(type_tarray),   intent(inout) :: dTdt
      type(type_tarray),   intent(inout) :: dxidt
      type(type_tarray),   intent(inout) :: dpsidt
      type(vp_bal_type),   intent(inout) :: vp_bal
      type(vort_bal_type), intent(inout) :: vort_bal


      if ( l_cheb_coll ) then
         if ( l_heat ) call finish_exp_temp_coll(us_Mloc, dVsT_Mloc, &
                            &                    dTdt%expl(:,:,tscheme%istage))
         if ( l_chem ) call finish_exp_xi_coll(us_Mloc, dVsXi_Mloc, &
                            &                  dxidt%expl(:,:,tscheme%istage))
         if ( l_direct_solve ) then
            call finish_exp_psi_coll_smat(us_Mloc, dVsOm_Mloc, temp_Mloc, xi_Mloc, &
                 &                        dpsidt%expl(:,:,tscheme%istage),&
                 &                        vort_bal)
         else
            call finish_exp_psi_coll_dmat(us_Mloc, up_Mloc, om_Mloc, dVsOm_Mloc, &
                 &                        temp_Mloc, xi_Mloc,                    &
                 &                        dpsidt%expl(:,:,tscheme%istage),       &
                 &                        vort_bal)
         end if
      else
         if ( l_heat ) call finish_exp_temp_int(psi_Mloc, dVsT_Mloc,   &
                            &                   dTdt%expl(:,:,tscheme%istage))
         if ( l_chem ) call finish_exp_xi_int(psi_Mloc, dVsXi_Mloc,  &
                            &                 dxidt%expl(:,:,tscheme%istage))
         if ( l_direct_solve ) then
            call finish_exp_psi_int_smat(psi_Mloc, us_Mloc, up_Mloc, om_Mloc, &
                 &                       dVsOm_Mloc, temp_Mloc, xi_Mloc,      &
                 &                       dpsidt%expl(:,:,tscheme%istage),     &
                 &                       vp_bal, vort_bal)
         else
            call finish_exp_psi_int_dmat(psi_Mloc, us_Mloc, up_Mloc, om_Mloc, &
                 &                       dVsOm_Mloc, temp_Mloc, xi_Mloc,      &
                 &                       dpsidt%expl(:,:,tscheme%istage),     &
                 &                       vp_bal, vort_bal)
         end if
      end if

   end subroutine finish_explicit_assembly
!------------------------------------------------------------------------------
end module mloop_mod
