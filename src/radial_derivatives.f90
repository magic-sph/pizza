module radial_der
   !
   ! Radial derivatives functions
   !

   use parallel_mod
   use constants, only: zero, one, three
   use precision_mod
   use mem_alloc
   use radial_scheme, only: type_rscheme
   use useful, only: abortRun

   implicit none

   private

   interface get_dcheb
      module procedure get_dcheb_real_1d
      module procedure get_dcheb_complex_2d
   end interface get_dcheb

   interface get_ddcheb
      module procedure get_ddcheb_real_1d
      module procedure get_ddcheb_complex_2d
   end interface get_ddcheb

   interface get_dr
      module procedure get_dr_real_1d
      module procedure get_dr_complex_2d
   end interface get_dr

   interface get_ddr
      module procedure get_ddr_real_1d
      module procedure get_ddr_complex_2d
   end interface get_ddr

   public :: get_ddr, get_dcheb, get_dr, initialize_der_arrays, exch_ghosts, &
   &         finalize_der_arrays, get_ddr_ghost, get_dr_Rloc, bulk_to_ghost, &
   &         get_ddddr_ghost

   complex(cp), allocatable :: work(:,:)
   real(cp), allocatable :: work_1d_real(:)
   real(cp) :: thr

contains

!------------------------------------------------------------------------------
   subroutine initialize_der_arrays(n_r_max, nMstart, nMstop, l_rerror, rerror_fac)
      !
      ! Allocate work arrays to compute derivatives
      !

      integer,  intent(in) :: n_r_max
      integer,  intent(in) :: nMstart
      integer,  intent(in) :: nMstop
      logical,  intent(in) :: l_rerror
      real(cp), intent(in) :: rerror_fac

      if ( l_rerror ) then
         thr = rerror_fac * epsilon(1.0_cp)
      else
         thr = 0.0_cp
      end if

      allocate( work_1d_real(n_r_max) )
      allocate( work(nMstart:nMstop,n_r_max) )
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL+&
      &                 n_r_max*(nMstop-nMstart+1)*SIZEOF_DEF_COMPLEX

   end subroutine initialize_der_arrays
!------------------------------------------------------------------------------
   subroutine finalize_der_arrays
      !
      ! Deallocate work arrays
      !
      deallocate( work_1d_real, work )

   end subroutine finalize_der_arrays
!------------------------------------------------------------------------------
   subroutine get_dcheb_complex_2d(f,df,nMstart,nMstop,n_r_max,n_cheb_max)
      !
      !  Returns chebychev coeffitiens of first derivative df and second  
      !  derivative ddf for a function whose cheb-coeff. are given as     
      !  columns in array f(nMstart:nMstop,n_r_max).                             
      !

      !-- Input variables:
      integer,     intent(in) :: nMstart  ! No of function to start with
      integer,     intent(in) :: nMstop   ! No of function to stop with
      integer,     intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,     intent(in) :: n_cheb_max ! Number of cheb modes
      complex(cp), intent(in) :: f(nMstart:nMstop,n_r_max)

      !-- Output variables:
      complex(cp), intent(out) ::  df(nMstart:nMstop,n_r_max)

      !-- Local variables:
      integer :: n_m,n_cheb
      real(cp) :: fac_cheb
      real(cp) :: threshold(nMstart:nMstop)

      do n_m=nMstart, nMstop
         threshold(n_m) = maxval(abs(f(n_m,:)))*thr
      end do

      !-- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         do n_m=nMstart,nMstop
            df(n_m,n_cheb)=zero
         end do
      end do
      n_cheb  =n_cheb_max-1
      if ( n_r_max == n_cheb_max ) then
         fac_cheb=real(n_cheb,kind=cp)
      else
         fac_cheb=real(2*n_cheb,kind=cp)
      end if
      do n_m=nMstart,nMstop
         if ( abs(df(n_m,n_cheb)) >= threshold(n_m) ) then
            df(n_m,n_cheb)=fac_cheb*f(n_m,n_cheb+1)
         else
            df(n_m,n_cheb)=zero
         end if
      end do

      !----- Recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=real(2*n_cheb,kind=cp)
         do n_m=nMstart,nMstop
            if ( abs(f(n_m,n_cheb+1)) >= threshold(n_m) ) then
               df(n_m,n_cheb)=df(n_m,n_cheb+2) + fac_cheb*f(n_m,n_cheb+1)
            else
               df(n_m,n_cheb)=df(n_m,n_cheb+2)
            end if
         end do
      end do

   end subroutine get_dcheb_complex_2d
!------------------------------------------------------------------------------
   subroutine get_dcheb_real_1d(f,df,n_r_max,n_cheb_max)

      !-- Input variables:
      integer,  intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,  intent(in) :: n_cheb_max ! Number of cheb modes
      real(cp), intent(in) :: f(n_r_max)

      !-- Output variables:
      real(cp), intent(out) ::  df(n_r_max)

      !-- Local variables:
      integer :: n_cheb
      real(cp) :: fac_cheb, threshold

      threshold = maxval(abs(f))*thr

      !-- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         df(n_cheb)=0.0_cp
      end do
      n_cheb  =n_cheb_max-1
      if ( n_r_max == n_cheb_max ) then
         fac_cheb=real(n_cheb,kind=cp)
      else
         fac_cheb=real(2*n_cheb,kind=cp)
      end if
      if ( abs(f(n_cheb+1)) >= threshold ) then
         df(n_cheb)=fac_cheb*f(n_cheb+1)
      else
         df(n_cheb)=0.0_cp
      end if

      !----- Recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=real(2*n_cheb,kind=cp)
         if ( abs(f(n_cheb+1)) >= threshold ) then
            df(n_cheb)=df(n_cheb+2) + fac_cheb*f(n_cheb+1)
         else
            df(n_cheb)=0.0_cp
         end if
      end do

   end subroutine get_dcheb_real_1d
!------------------------------------------------------------------------------
   subroutine get_ddcheb_complex_2d(f,df,ddf,nMstart,nMstop,n_r_max,n_cheb_max)
      !
      !  Returns chebychev coefficents of first derivative df and second  
      !  derivative ddf for a function whose cheb-coeff. are given as     
      !  columns in array f(n_c_tot,n_r_max).                             
      !
    
      !-- Input variables:
      integer,     intent(in) :: nMstart  ! No of column to start with
      integer,     intent(in) :: nMstop   ! No of column to stop with
      integer,     intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,     intent(in) :: n_cheb_max ! Number of cheb modes
      complex(cp), intent(in) :: f(nMstart:nMstop,n_r_max)
    
      !-- Output variables:
      complex(cp), intent(out) ::  df(nMstart:nMstop,n_r_max)
      complex(cp), intent(out) ::  ddf(nMstart:nMstop,n_r_max)
    
      !-- local variables:
      integer :: n_m,n_cheb
      real(cp) :: fac_cheb
      real(cp) :: threshold(nMstart:nMstop)

      do n_m=nMstart, nMstop
         threshold(n_m) = maxval(abs(f(n_m,:)))*thr
      end do

      !----- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         do n_m=nMstart,nMstop
            df(n_m,n_cheb) =zero
            ddf(n_m,n_cheb)=zero
         end do
      end do
      n_cheb=n_cheb_max-1
      if ( n_cheb_max == n_r_max ) then
         fac_cheb=real(n_cheb,kind=cp)
      else
         fac_cheb=real(2*n_cheb,kind=cp)
      end if
      do n_m=nMstart,nMstop
         if ( abs(df(n_m,n_cheb)) >= threshold(n_m) ) then
            df(n_m,n_cheb)=fac_cheb*f(n_m,n_cheb+1)
         else
            df(n_m,n_cheb)=zero
         end if
         ddf(n_m,n_cheb)=zero
      end do
    
      !----- recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=real(2*n_cheb,kind=cp)
         do n_m=nMstart,nMstop
            if ( abs(f(n_m,n_cheb+1)) >= threshold(n_m) ) then
               df(n_m,n_cheb) = df(n_m,n_cheb+2) + fac_cheb* f(n_m,n_cheb+1)
            else
               df(n_m,n_cheb) = df(n_m,n_cheb+2)
            end if
            if ( abs(df(n_m,n_cheb+1)) >= threshold(n_m) ) then
               ddf(n_m,n_cheb)=ddf(n_m,n_cheb+2) + fac_cheb*df(n_m,n_cheb+1)
            else
               ddf(n_m,n_cheb)=ddf(n_m,n_cheb+2)
            end if
         end do
      end do

   end subroutine get_ddcheb_complex_2d
!------------------------------------------------------------------------------
   subroutine get_ddcheb_real_1d(f,df,ddf,n_r_max,n_cheb_max)
      !
      !  Returns chebychev coefficents of first derivative df and second  
      !  derivative ddf for a function whose cheb-coeff. are given as     
      !  columns in array f(n_c_tot,n_r_max).                             
      !
    
      !-- Input variables:
      integer,  intent(in) :: n_r_max    ! second dimension of f,df,ddf
      integer,  intent(in) :: n_cheb_max ! Number of cheb modes
      real(cp), intent(in) :: f(n_r_max)
    
      !-- Output variables:
      real(cp), intent(out) ::  df(n_r_max)
      real(cp), intent(out) ::  ddf(n_r_max)
    
      !-- local variables:
      integer :: n_cheb
      real(cp) :: fac_cheb, threshold

      threshold = maxval(abs(f))*thr
    
      !----- initialize derivatives:
      do n_cheb=n_cheb_max,n_r_max
         df(n_cheb) =0.0_cp
         ddf(n_cheb)=0.0_cp
      end do
      n_cheb=n_cheb_max-1
      if ( n_cheb_max == n_r_max ) then
         fac_cheb=real(n_cheb,kind=cp)
      else
         fac_cheb=real(2*n_cheb,kind=cp)
      end if
      if ( abs(f(n_cheb+1)) >= threshold ) then
         df(n_cheb) =fac_cheb*f(n_cheb+1)
      else
         df(n_cheb) =0.0_cp
      end if
      ddf(n_cheb)=0.0_cp
    
      !----- recursion
      do n_cheb=n_cheb_max-2,1,-1
         fac_cheb=real(2*n_cheb,kind=cp)
         if ( abs(f(n_cheb+1)) >= threshold ) then
            df(n_cheb) = df(n_cheb+2) + fac_cheb* f(n_cheb+1)
         else
            df(n_cheb) = df(n_cheb+2)
         end if

         if ( abs(df(n_cheb+1)) >= threshold ) then
            ddf(n_cheb)=ddf(n_cheb+2) + fac_cheb*df(n_cheb+1)
         else
            ddf(n_cheb)=ddf(n_cheb+2)
         end if
      end do

   end subroutine get_ddcheb_real_1d
!------------------------------------------------------------------------------
   subroutine get_dr_complex_2d(f,df,nMstart,nMstop,n_r_max,r_scheme,nocopy,  &
              &                 l_dct_in,l_dct_out)
      !
      !  Returns first radial derivative df of the input function f.      
      !  Array f(nMstart:nMstop,*) may contain several functions numbered by     
      !  the first index. The subroutine calculates the derivaties of     
      !  the functions f(nMstart,*) to f(nMstop) by transforming      
      !  to a Chebychev representation using n_r_max radial grid points . 
      !
    
      !-- Input variables:
      integer,             intent(in) :: n_r_max  ! number of radial grid points
      integer,             intent(in) :: nMstart  ! first function to be treated
      integer,             intent(in) :: nMstop   ! last function to be treated
      complex(cp),         intent(inout) :: f(nMstart:nMstop,n_r_max)
      class(type_rscheme), intent(in) :: r_scheme
      logical, optional,   intent(in) :: nocopy
      logical, optional,   intent(in) :: l_dct_in
      logical, optional,   intent(in) :: l_dct_out
    
      !-- Output variables:
      complex(cp), intent(out) :: df(nMstart:nMstop,n_r_max) ! first derivative of f
    
      !-- Local:
      integer :: n_r,n_f,od
      logical :: copy_array, l_dct_in_loc, l_dct_out_loc

      if ( present(l_dct_in) ) then
         l_dct_in_loc=l_dct_in
      else
         l_dct_in_loc=.true.
      end if

      if ( present(l_dct_out) ) then
         l_dct_out_loc=l_dct_out
      else
         l_dct_out_loc=.true.
      end if
    
      if ( r_scheme%version == 'cheb' ) then

         if ( present(nocopy) ) then
            copy_array=.false.
         else
            copy_array=.true.
         end if
    
         if ( copy_array )  then
            do n_r=1,n_r_max
               do n_f=nMstart,nMstop
                  work(n_f,n_r)=f(n_f,n_r)
               end do
            end do
       
            !-- Transform f to cheb space:
            if ( l_dct_in_loc ) call r_scheme%costf1(work,nMstart,nMstop,n_r_max)
          
            !-- Get derivatives:
            !call get_dcheb(work,df,nMstart,nMstop,n_r_max,r_scheme%n_max)
            call get_dcheb(work,df,nMstart,nMstop,n_r_max,n_r_max)
          
            !-- Transform back:
            if ( l_dct_out_loc ) call r_scheme%costf1(df,nMstart,nMstop,n_r_max)

         else

            !-- Transform f to cheb space:
            if ( l_dct_in_loc ) call r_scheme%costf1(f,nMstart,nMstop,n_r_max)
          
            !-- Get derivatives:
            !call get_dcheb(f,df,nMstart,nMstop,n_r_max,r_scheme%n_max)
            call get_dcheb(f,df,nMstart,nMstop,n_r_max,n_r_max)
          
            !-- Transform back:
            if ( l_dct_in_loc ) call r_scheme%costf1(f,nMstart,nMstop,n_r_max)
            if ( l_dct_out_loc ) call r_scheme%costf1(df,nMstart,nMstop,n_r_max)

         end if
       
         !-- New map:
         do n_r=1,n_r_max
            do n_f=nMstart,nMstop
               df(n_f,n_r)=r_scheme%drx(n_r)*df(n_f,n_r)
            end do
         end do

      else

         !-- Initialise to zero:
         do n_r=1,n_r_max
            do n_f=nMstart,nMstop
               df(n_f,n_r) =zero
            end do
         end do

         !-- Bulk points for 1st derivative
         do od=0,r_scheme%order
            do n_r=1+r_scheme%order/2,n_r_max-r_scheme%order/2
               do n_f=nMstart,nMstop
                  df(n_f,n_r)=df(n_f,n_r)+r_scheme%dr(n_r,od)*f(n_f,n_r-r_scheme%order/2+od)
               end do
            end do
         end do

         !-- Boundary points for 1st derivative
         do od=0,r_scheme%order_boundary
            do n_r=1,r_scheme%order/2
               do n_f=nMstart,nMstop
                  df(n_f,n_r) = df(n_f,n_r)+r_scheme%dr_top(n_r,od) * f(n_f,od+1)
               end do
            end do
            do n_r=1,r_scheme%order/2
               do n_f=nMstart,nMstop
                  df(n_f,n_r_max-n_r+1) = df(n_f,n_r_max-n_r+1)+               &
                  &                       r_scheme%dr_bot(n_r,od)*f(n_f,n_r_max-od)
               end do
            end do
         end do

      end if

   end subroutine get_dr_complex_2d
!------------------------------------------------------------------------------
   subroutine get_dr_real_1d(f,df,n_r_max,r_scheme,l_dct)
    
      !-- Input variables:
      integer,             intent(in) :: n_r_max  ! number of radial grid points
      real(cp),            intent(inout) :: f(n_r_max)
      class(type_rscheme), intent(in) :: r_scheme
      logical, optional,   intent(in) :: l_dct
    
      !-- Output variables:
      real(cp), intent(out) :: df(n_r_max) ! first derivative of f
    
      !-- Local:
      integer :: n_r,od
      logical :: l_dct_loc

      if ( present(l_dct) ) then
         l_dct_loc=l_dct
      else
         l_dct_loc=.true.
      end if
    
      if ( r_scheme%version == 'cheb' ) then

         do n_r=1,n_r_max
            work_1d_real(n_r)=f(n_r)
         end do
       
         !-- Transform f to cheb space:
         if ( l_dct_loc ) call r_scheme%costf1(work_1d_real,n_r_max)
          
         !-- Get derivatives:
         !call get_dcheb(work_1d_real,df,n_r_max,r_scheme%n_max)
         call get_dcheb(work_1d_real,df,n_r_max,n_r_max)
          
         !-- Transform back:
         call r_scheme%costf1(df,n_r_max)

         !-- New map:
         do n_r=1,n_r_max
            df(n_r)=r_scheme%drx(n_r)*df(n_r)
         end do

      else

         !-- Initialise to zero:
         do n_r=1,n_r_max
            df(n_r)=0.0_cp
         end do

         !-- Bulk points for 1st derivative
         do od=0,r_scheme%order
            do n_r=1+r_scheme%order/2,n_r_max-r_scheme%order/2
               df(n_r)=df(n_r)+r_scheme%dr(n_r,od)*f(n_r-r_scheme%order/2+od)
            end do
         end do

         !-- Boundary points for 1st derivative
         do od=0,r_scheme%order_boundary
            do n_r=1,r_scheme%order/2
               df(n_r) = df(n_r)+r_scheme%dr_top(n_r,od) * f(od+1)
            end do
            do n_r=1,r_scheme%order/2
               df(n_r_max-n_r+1) = df(n_r_max-n_r+1)+               &
               &                   r_scheme%dr_bot(n_r,od)*f(n_r_max-od)
            end do
         end do

      end if

   end subroutine get_dr_real_1d
!------------------------------------------------------------------------------
   subroutine get_ddr_complex_2d(f,df,ddf,nMstart,nMstop,n_r_max,r_scheme,l_dct)
      !
      !  Returns first radial derivative df and second radial             
      !  derivative ddf of the input function f.                          
      !  Array f(nMstart:nMstop,*) may contain several functions numbered by     
      !  the first index. The subroutine calculates the derivatives of    
      !  the functions f(nMstart,*) to f(nMstop) by transforming      
      !  to a Chebychev representation using n_r_max radial grid points.  
      !
    
      !-- Input variables:
      integer,             intent(in) :: n_r_max  ! number of radial grid points
      integer,             intent(in) :: nMstart  ! first function to be treated
      integer,             intent(in) :: nMstop   ! last function to be treated
      complex(cp),         intent(in) :: f(nMstart:nMstop,n_r_max)
      class(type_rscheme), intent(in) :: r_scheme
      logical, optional,   intent(in) :: l_dct
    
      !-- Output variables:
      complex(cp), intent(out) :: df(nMstart:nMstop,n_r_max) ! first derivative of f
      complex(cp), intent(out) :: ddf(nMstart:nMstop,n_r_max)! second derivative of f
    
      !-- Local variables:
      logical :: l_dct_loc
      integer :: n_r,n_f,od

      if ( present(l_dct) ) then
         l_dct_loc=l_dct
      else
         l_dct_loc=.true.
      end if

      if ( r_scheme%version == 'cheb' ) then
    
         !-- Copy input functions:
         do n_r=1,n_r_max
            do n_f=nMstart,nMstop
               work(n_f,n_r)=f(n_f,n_r)
            end do
         end do
    
         !-- Transform f to cheb space:
         if ( l_dct_loc ) call r_scheme%costf1(work,nMstart,nMstop,n_r_max)
    
         !-- Get derivatives:
         !call get_ddcheb(work,df,ddf,nMstart,nMstop,n_r_max,r_scheme%n_max)
         call get_ddcheb(work,df,ddf,nMstart,nMstop,n_r_max,n_r_max)
    
         !-- Transform back:
         call r_scheme%costf1(df,nMstart,nMstop,n_r_max)
         call r_scheme%costf1(ddf,nMstart,nMstop,n_r_max)
    
         !-- New map:
         do n_r=1,n_r_max
            do n_f=nMstart,nMstop
               ddf(n_f,n_r)=r_scheme%ddrx(n_r)*df(n_f,n_r)+&
               &            r_scheme%drx(n_r)*r_scheme%drx(n_r)*ddf(n_f,n_r)
               df(n_f,n_r) =r_scheme%drx(n_r)*df(n_f,n_r)
            end do
         end do

      else

         !-- Initialise to zero:
         do n_r=1,n_r_max
            do n_f=nMstart,nMstop
               df(n_f,n_r) =zero
               ddf(n_f,n_r)=zero
            end do
         end do

         !-- Bulk points for 1st and 2nd derivatives
         do od=0,r_scheme%order
            do n_r=1+r_scheme%order/2,n_r_max-r_scheme%order/2
               do n_f=nMstart,nMstop
                  df(n_f,n_r)  = df(n_f,n_r) + r_scheme%dr(n_r,od) * f(n_f,n_r-r_scheme%order/2+od)
                  ddf(n_f,n_r) = ddf(n_f,n_r)+r_scheme%ddr(n_r,od) * f(n_f,n_r-r_scheme%order/2+od)
               end do
            end do
         end do

         !-- Boundary points for 1st derivative
         do od=0,r_scheme%order_boundary
            do n_r=1,r_scheme%order/2
               do n_f=nMstart,nMstop
                  df(n_f,n_r) = df(n_f,n_r)+r_scheme%dr_top(n_r,od) * f(n_f,od+1)
               end do
            end do
            do n_r=1,r_scheme%order/2
               do n_f=nMstart,nMstop
                  df(n_f,n_r_max-n_r+1) = df(n_f,n_r_max-n_r+1)+               &
                  &                       r_scheme%dr_bot(n_r,od)*f(n_f,n_r_max-od)
               end do
            end do
         end do

         !-- Boundary points for 2nd derivative
         do od=0,r_scheme%order_boundary+1
            do n_r=1,r_scheme%order/2
               do n_f=nMstart,nMstop
                  ddf(n_f,n_r) = ddf(n_f,n_r)+r_scheme%ddr_top(n_r,od) * f(n_f,od+1)
               end do
            end do
            do n_r=1,r_scheme%order/2
               do n_f=nMstart,nMstop
                  ddf(n_f,n_r_max-n_r+1) = ddf(n_f,n_r_max-n_r+1)+               &
                  &                       r_scheme%ddr_bot(n_r,od)*f(n_f,n_r_max-od)
               end do
            end do
         end do

      end if

   end subroutine get_ddr_complex_2d
!------------------------------------------------------------------------------
   subroutine get_ddr_real_1d(f,df,ddf,n_r_max,r_scheme)
    
      !-- Input variables:
      integer,             intent(in) :: n_r_max  ! number of radial grid points
      real(cp),            intent(in) :: f(n_r_max)
      class(type_rscheme), intent(in) :: r_scheme
    
      !-- Output variables:
      real(cp), intent(out) :: df(n_r_max) ! first derivative of f
      real(cp), intent(out) :: ddf(n_r_max)! second derivative of f
    
      !-- Local variables:
      integer :: n_r,od

      if ( r_scheme%version == 'cheb' ) then
    
         !-- Copy input functions:
         do n_r=1,n_r_max
            work_1d_real(n_r)=f(n_r)
         end do
    
         !-- Transform f to cheb space:
         call r_scheme%costf1(work_1d_real,n_r_max)
    
         !-- Get derivatives:
         !call get_ddcheb(work_1d_real,df,ddf,n_r_max,r_scheme%n_max)
         call get_ddcheb(work_1d_real,df,ddf,n_r_max,n_r_max)
    
         !-- Transform back:
         call r_scheme%costf1(df,n_r_max)
         call r_scheme%costf1(ddf,n_r_max)
    
         !-- New map:
         do n_r=1,n_r_max
            ddf(n_r)=r_scheme%ddrx(n_r)*df(n_r)+&
            &         r_scheme%drx(n_r)*r_scheme%drx(n_r)*ddf(n_r)
            df(n_r) =r_scheme%drx(n_r)*df(n_r)
         end do

      else

         !-- Initialise to zero:
         do n_r=1,n_r_max
            df(n_r) =0.0_cp
            ddf(n_r)=0.0_cp
         end do

         !-- Bulk points for 1st and 2nd derivatives
         do od=0,r_scheme%order
            do n_r=1+r_scheme%order/2,n_r_max-r_scheme%order/2
               df(n_r)  = df(n_r) + r_scheme%dr(n_r,od) * f(n_r-r_scheme%order/2+od)
               ddf(n_r) = ddf(n_r)+r_scheme%ddr(n_r,od) * f(n_r-r_scheme%order/2+od)
            end do
         end do

         !-- Boundary points for 1st derivative
         do od=0,r_scheme%order_boundary
            do n_r=1,r_scheme%order/2
               df(n_r) = df(n_r)+r_scheme%dr_top(n_r,od) * f(od+1)
            end do
            do n_r=1,r_scheme%order/2
               df(n_r_max-n_r+1) = df(n_r_max-n_r+1)+               &
               &                   r_scheme%dr_bot(n_r,od)*f(n_r_max-od)
            end do
         end do

         !-- Boundary points for 2nd derivative
         do od=0,r_scheme%order_boundary+1
            do n_r=1,r_scheme%order/2
               ddf(n_r) = ddf(n_r)+r_scheme%ddr_top(n_r,od) * f(od+1)
            end do
            do n_r=1,r_scheme%order/2
               ddf(n_r_max-n_r+1) = ddf(n_r_max-n_r+1)+               &
               &                    r_scheme%ddr_bot(n_r,od)*f(n_r_max-od)
            end do
         end do

      end if

   end subroutine get_ddr_real_1d
!------------------------------------------------------------------------------
   subroutine get_dr_Rloc(f_Rloc, df_Rloc, n_m_max, nRstart, nRstop, n_r_max, r_scheme)
      !
      ! Purpose of this subroutine is to take the first radial derivative of an input
      ! complex array distributed over radius. This can only be used with
      ! finite differences.
      !

      !-- Input variables
      integer,             intent(in) :: n_m_max, nRstart, nRstop, n_r_max
      class(type_rscheme), intent(in) :: r_scheme
      complex(cp),         intent(in) :: f_Rloc(n_m_max,nRstart:nRstop)

      !-- Output variable
      complex(cp), intent(out) ::  df_Rloc(n_m_max,nRstart:nRstop)

      !-- Local variables:
      complex(cp) :: work_ghost(n_m_max,nRstart-r_scheme%order/2:nRstop+r_scheme%order/2)
      complex(cp) :: ftop(n_m_max,r_scheme%order_boundary+1)
      complex(cp) :: fbot(n_m_max,n_r_max-r_scheme%order_boundary:n_r_max)
      integer :: n_r, od, start_m, stop_m, n_m

      if ( (r_scheme%order>2 .or. r_scheme%order_boundary>2) .and. &
      &    (nRstop-nRstart+1)<r_scheme%order ) then
         call abortRun('Distributed r-der not implemented in this case yet!')
      end if

      !$omp parallel default(shared) private(start_m,stop_m,n_m)
      start_m=1; stop_m=n_m_max
      call get_openmp_blocks(start_m,stop_m)

      !-- Copy input array
      work_ghost(start_m:stop_m,nRstart:nRstop)=f_Rloc(start_m:stop_m,:)
      do n_r=1,r_scheme%order_boundary+1
         if (n_r >= nRstart .and. n_r <= nRstop) then
            ftop(start_m:stop_m,n_r)=f_Rloc(start_m:stop_m,n_r)
         end if
      end do

      do n_r=n_r_max-r_scheme%order_boundary,n_r_max
         if (n_r >= nRstart .and. n_r <= nRstop) then
            fbot(start_m:stop_m,n_r)=f_Rloc(start_m:stop_m,n_r)
         end if
      end do

     !-- Exchange the ghost zones
      !$omp barrier
      !$omp master
      call exch_ghosts(work_ghost, n_m_max, nRstart, nRstop, r_scheme%order/2)
      !$omp end master
      !$omp barrier

      !-- Bulk points for 1st derivative
      do n_r=nRstart,nRstop
         do n_m=start_m,stop_m
            df_Rloc(n_m,n_r)=r_scheme%dr(n_r,0)*work_ghost(n_m,n_r-1)+ &
            &                r_scheme%dr(n_r,1)*work_ghost(n_m,n_r)+   &
            &                r_scheme%dr(n_r,2)*work_ghost(n_m,n_r+1)
         end do
      end do

      !-- Exchange boundary values
      !$omp barrier
      !$omp master
      call get_bound_vals(fbot, ftop, n_m_max, nRstart, nRstop, n_r_max, &
           &              r_scheme%order_boundary+1)
      !$omp end master
      !$omp barrier

      !-- Boundary points for 1st derivative
      if ( rank == 0 ) then
         do n_m=start_m,stop_m
            df_Rloc(n_m,1)=zero
            do od=0,r_scheme%order_boundary
               df_Rloc(n_m,1)=df_Rloc(n_m,1) + r_scheme%dr_top(1,od)*ftop(n_m,od+1)
            end do
         end do
      end if

      if ( rank == n_procs -1 ) then
         do n_m=start_m,stop_m
            df_Rloc(n_m,n_r_max)=zero
            do od=0,r_scheme%order_boundary
               df_Rloc(n_m,n_r_max)=df_Rloc(n_m,n_r_max)+r_scheme%dr_bot(1,od)* &
               &                    fbot(n_m,n_r_max-od)
            end do
         end do
      end if
      !$omp end parallel

   end subroutine get_dr_Rloc
!------------------------------------------------------------------------------
   subroutine get_ddr_ghost(f_Rloc, df_Rloc, ddf_Rloc, n_m_max, start_m, stop_m, &
              &             nRstart, nRstop, r_scheme)
      !
      ! Purpose of this subroutine is to take the first and second
      ! radial derivatives of an input complex array distributed over radius that
      ! has the ghost zones properly filled.
      !

      !-- Input variables
      integer,             intent(in) :: n_m_max, nRstart, nRstop, start_m, stop_m
      class(type_rscheme), intent(in) :: r_scheme
      complex(cp),         intent(in) :: f_Rloc(n_m_max,nRstart-1:nRstop+1)

      !-- Output variable
      complex(cp), intent(out) ::  df_Rloc(n_m_max,nRstart:nRstop)
      complex(cp), intent(out) ::  ddf_Rloc(n_m_max,nRstart:nRstop)

      !-- Local variables:
      integer :: n_r, n_m

      if ( (r_scheme%order>2 .or. r_scheme%order_boundary>2) .and. &
      &    (nRstop-nRstart+1)<r_scheme%order ) then
         call abortRun('Distributed r-der not implemented in this case yet!')
      end if

      !-- Bulk points for 1st and 2nd derivatives
      do n_r=nRstart,nRstop
         do n_m=start_m,stop_m
            df_Rloc(n_m,n_r)=r_scheme%dr(n_r,0)*f_Rloc(n_m,n_r-1) + &
            &                r_scheme%dr(n_r,1)*f_Rloc(n_m,n_r)   + &
            &                r_scheme%dr(n_r,2)*f_Rloc(n_m,n_r+1)
            ddf_Rloc(n_m,n_r)=r_scheme%ddr(n_r,0)*f_Rloc(n_m,n_r-1) + &
            &                 r_scheme%ddr(n_r,1)*f_Rloc(n_m,n_r)   + &
            &                 r_scheme%ddr(n_r,2)*f_Rloc(n_m,n_r+1)
         end do
      end do

   end subroutine get_ddr_ghost
!------------------------------------------------------------------------------
   subroutine get_ddddr_ghost(f_Rloc, df_Rloc, ddf_Rloc, dddf_Rloc, ddddf_Rloc, &
              &               n_m_max, start_m, stop_m, nRstart, nRstop, r_scheme)
      !
      ! Purpose of this subroutine is to take the first and second
      ! radial derivatives of an input complex array distributed over radius that
      ! has the ghost zones properly filled.
      !

      !-- Input variables
      integer,             intent(in) :: n_m_max, nRstart, nRstop, start_m, stop_m
      class(type_rscheme), intent(in) :: r_scheme
      complex(cp),         intent(in) :: f_Rloc(n_m_max,nRstart-2:nRstop+2)

      !-- Output variable
      complex(cp), intent(out) ::  df_Rloc(n_m_max,nRstart:nRstop)
      complex(cp), intent(out) ::  ddf_Rloc(n_m_max,nRstart:nRstop)
      complex(cp), intent(out) ::  dddf_Rloc(n_m_max,nRstart:nRstop)
      complex(cp), intent(out) ::  ddddf_Rloc(n_m_max,nRstart:nRstop)

      !-- Local variables:
      integer :: n_r, lm

      if ( (r_scheme%order>2 .or. r_scheme%order_boundary>2) .and. &
      &    (nRstop-nRstart+1)<r_scheme%order ) then
         call abortRun('Distributed r-der not implemented in this case yet!')
      end if

      !-- 1st and 2nd derivatives
      do n_r=nRstart,nRstop
         do lm=start_m,stop_m
            df_Rloc(lm,n_r)=r_scheme%dr(n_r,0)*f_Rloc(lm,n_r-1) + &
            &               r_scheme%dr(n_r,1)*f_Rloc(lm,n_r)   + &
            &               r_scheme%dr(n_r,2)*f_Rloc(lm,n_r+1)
            ddf_Rloc(lm,n_r)=r_scheme%ddr(n_r,0)*f_Rloc(lm,n_r-1) + &
            &                r_scheme%ddr(n_r,1)*f_Rloc(lm,n_r)   + &
            &                r_scheme%ddr(n_r,2)*f_Rloc(lm,n_r+1)
            dddf_Rloc(lm,n_r)=r_scheme%dddr(n_r,0)*f_Rloc(lm,n_r-2) + &
            &                 r_scheme%dddr(n_r,1)*f_Rloc(lm,n_r-1) + &
            &                 r_scheme%dddr(n_r,2)*f_Rloc(lm,n_r)   + &
            &                 r_scheme%dddr(n_r,3)*f_Rloc(lm,n_r+1) + &
            &                 r_scheme%dddr(n_r,4)*f_Rloc(lm,n_r+2)
            ddddf_Rloc(lm,n_r)=r_scheme%ddddr(n_r,0)*f_Rloc(lm,n_r-2) + &
            &                  r_scheme%ddddr(n_r,1)*f_Rloc(lm,n_r-1) + &
            &                  r_scheme%ddddr(n_r,2)*f_Rloc(lm,n_r)   + &
            &                  r_scheme%ddddr(n_r,3)*f_Rloc(lm,n_r+1) + &
            &                  r_scheme%ddddr(n_r,4)*f_Rloc(lm,n_r+2)
         end do
      end do

   end subroutine get_ddddr_ghost
!------------------------------------------------------------------------------
   subroutine exch_ghosts(f, n_m_max, nRstart, nRstop, nghosts)

      integer, intent(in) :: n_m_max, nRstart, nRstop, nghosts

      complex(cp), intent(inout) :: f(n_m_max, nRstart-nghosts:nRstop+nghosts)

      integer :: n_counts, rightProc, leftProc, st(MPI_STATUS_SIZE)

      leftProc=rank-1
      if ( leftProc < 0 ) leftProc = MPI_PROC_NULL
      rightProc=rank+1
      if ( rightProc >= n_procs ) rightProc = MPI_PROC_NULL
      n_counts = n_m_max*nghosts

      !-- Halo exchange in forward direction
      call MPI_SendRecv(f(:,nRstop-nghosts+1:nRstop), n_counts, MPI_DEF_COMPLEX,    &
           &            rightProc, 12345, f(:,nRstart-nghosts:nRstart-1), n_counts, &
           &            MPI_DEF_COMPLEX, leftProc, 12345, MPI_COMM_WORLD, st, ierr)

      !-- Halo exchange in backward direction
      call MPI_SendRecv(f(:,nRstart:nRstart+nghosts-1), n_counts, MPI_DEF_COMPLEX, &
           &            leftProc, 12345, f(:,nRstop+1:nRstop+nghosts), n_counts,   &
           &            MPI_DEF_COMPLEX, rightProc, 12345, MPI_COMM_WORLD, st, ierr)

   end subroutine exch_ghosts
!------------------------------------------------------------------------------
   subroutine get_bound_vals(fbot, ftop, n_m_max, nRstart, nRstop, n_r_max, nbounds)

      !-- Input variables
      integer, intent(in) :: n_m_max, nRstart, nRstop
      integer, intent(in) :: nbounds, n_r_max

      !-- Output boundary values
      complex(cp), intent(inout) :: ftop(n_m_max,nbounds)
      complex(cp), intent(inout) :: fbot(n_m_max,n_r_max-nbounds+1:n_r_max)

      !-- Local variables
      integer :: nreq, nR, tag, req(2*nbounds)

      if ( nRstart > nbounds .and. nRstop < n_r_max-nbounds ) return

      req(:)=MPI_REQUEST_NULL
      nreq=0
      do nR=1,nbounds
         tag = 754432+nR
         if ( rank /= 0 .and. nRstart<=nR .and. nRstop>=nR) then
            call MPI_Send(ftop(:,nR), n_m_max, MPI_DEF_COMPLEX, 0, tag, &
                 &        MPI_COMM_WORLD, ierr)
         else if ( rank == 0 .and. nRstop < nR ) then
            nreq=nreq+1
            call MPI_IRecv(ftop(:,nR), n_m_max, MPI_DEF_COMPLEX, MPI_ANY_SOURCE, tag, &
                 &         MPI_COMM_WORLD, req(nreq), ierr)
         end if
      end do
      call MPI_Waitall(nreq, req(1:nreq), MPI_STATUSES_IGNORE, ierr)

      req(:)=MPI_REQUEST_NULL
      nreq=0
      do nR=n_r_max-nbounds+1,n_r_max
         tag = 92113+nR
         if ( rank /= n_procs-1 .and. nRstart<=nR .and. nRstop>=nR) then
            call MPI_Send(fbot(:,nR), n_m_max, MPI_DEF_COMPLEX, n_procs-1, tag, &
                 &        MPI_COMM_WORLD, ierr)
         else if ( rank == n_procs-1 .and. nR < nRstart ) then
            nreq=nreq+1
            call MPI_IRecv(fbot(:,nR), n_m_max, MPI_DEF_COMPLEX, MPI_ANY_SOURCE, tag, &
                 &         MPI_COMM_WORLD, req(nreq), ierr)
         end if
      end do
      call MPI_Waitall(nreq, req(1:nreq), MPI_STATUSES_IGNORE, ierr)

   end subroutine get_bound_vals
!------------------------------------------------------------------------------
   subroutine bulk_to_ghost(x, x_g, ng, nRstart, nRstop, n_m_max, start_m, stop_m)
      !
      ! This subroutine is used to copy an array that is defined from nRstart to
      ! nRstop to an array that is defined from nRstart-1 to nRstop+1
      !

      !-- Input variables
      integer,     intent(in) :: start_m, stop_m, nRstart, nRstop
      integer,     intent(in) :: n_m_max
      integer,     intent(in) :: ng ! Number of ghost zones
      complex(cp), intent(in) :: x(n_m_max,nRstart:nRstop)

      !-- Output variable
      complex(cp), intent(out) :: x_g(n_m_max,nRstart-ng:nRstop+ng)

      !-- Local variables
      integer :: n_r, n_m

      do n_r=nRstart,nRstop
         do n_m=start_m,stop_m
            x_g(n_m,n_r)=x(n_m,n_r)
         end do
      end do

   end subroutine bulk_to_ghost
!------------------------------------------------------------------------------
end module radial_der
