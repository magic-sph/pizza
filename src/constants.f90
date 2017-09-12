module constants
   !
   ! module containing constants and parameters
   ! used in the code.
   !

   use precision_mod

   real(cp), parameter :: pi=4.0_cp*atan(1.0_cp) ! pi=3.1415..
   real(cp), parameter :: one  =1.0_cp ! 1
   real(cp), parameter :: two  =2.0_cp ! 2
   real(cp), parameter :: three=3.0_cp ! 3
   real(cp), parameter :: four =4.0_cp ! 4
   real(cp), parameter :: half =0.5_cp ! 0.5
   real(cp), parameter :: third=one/three ! 1/3
   complex(cp), parameter :: zero=(0.0_cp,0.0_cp) ! cmplx(0.0, 0.0)
   complex(cp), parameter :: ci=(0.0_cp,1.0_cp) ! cmplx(0.0, 1.0)

   real(cp), parameter :: tiny_number=10.0_cp*epsilon(1.0_cp)
   real(cp) :: vol_oc  ! Volume of the outer core
   real(cp) :: vol_otc ! Volume outside the tangent cylinder
   real(cp) :: surf    ! Disk surface

end module constants
