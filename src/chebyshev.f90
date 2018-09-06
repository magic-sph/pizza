module chebyshev

   use precision_mod
   use mem_alloc, only: bytes_allocated
   use constants, only: half, one, two, three, four, pi, third
   use blocking, only: nMstart,nMstop
   use radial_scheme, only: type_rscheme
   use dct_fftw, only: costf_t
   use namelists, only: map_function

   implicit none

   private

   type, public, extends(type_rscheme) :: type_cheb

      real(cp) :: alpha1 !Input parameter for non-linear map to define degree of spacing (0.0:2.0)
      real(cp) :: alpha2 !Input parameter for non-linear map to define central point of different spacing (-1.0:1.0)
      logical :: l_map
      type(costf_t) :: chebt
      real(cp), allocatable :: r_cheb(:)

   contains

      procedure :: initialize
      procedure :: finalize
      procedure :: get_der_mat
      procedure :: get_grid => initialize_mapping
      procedure :: costf1_complex_2d
      procedure :: costf1_real_1d

   end type type_cheb

contains

   subroutine initialize(this, n_r_max, order, order_boundary, l_cheb_coll, &
              &          no_work_array)
      !
      !  Purpose of this subroutine is to calculate and store several     
      !  values that will be needed for a fast cosine transform of the    
      !  first kind. The actual transform is performed by the             
      !  subroutine costf1.                                               
      !

      class(type_cheb) :: this
      
      !-- Input variables
      integer, intent(in) :: n_r_max
      integer, intent(in) :: order ! This is going to be n_cheb_max
      integer, intent(in) :: order_boundary ! this is used to determine whether mappings are used
      logical, intent(in) :: l_cheb_coll
      logical, optional, intent(in) :: no_work_array

      !-- Local variable
      logical :: l_work_array
            
      this%rnorm = sqrt(two/real(n_r_max-1,kind=cp))
      this%n_max = order  ! n_cheb_max
      this%boundary_fac = half
      this%version = 'cheb'
      this%nRmax = n_r_max
      this%order_boundary=order_boundary

      if ( order_boundary == 1 ) then
         this%l_map=.true.
      else
         this%l_map=.false.
      end if

      if ( l_cheb_coll ) then
         allocate( this%rMat(n_r_max,n_r_max) )
         allocate( this%drMat(n_r_max,n_r_max) )
         allocate( this%d2rMat(n_r_max,n_r_max) )
         bytes_allocated=bytes_allocated+3*n_r_max*n_r_max*SIZEOF_DEF_REAL
      else
         allocate( this%rMat(2,n_r_max) )
         allocate( this%drMat(2,n_r_max) )
         allocate( this%d2rMat(2,n_r_max) )
         bytes_allocated=bytes_allocated+6*n_r_max*SIZEOF_DEF_REAL
      end if

      if ( present(no_work_array) ) then
         l_work_array=no_work_array
      else
         l_work_array=.false.
      end if

      call this%chebt%initialize(nMstart, nMstop, n_r_max, l_work_array)

      allocate( this%r_cheb(n_r_max), this%drx(n_r_max), this%ddrx(n_r_max) )
      bytes_allocated=bytes_allocated+3*n_r_max*SIZEOF_DEF_REAL

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine initialize_mapping(this, n_r_max, ricb, rcmb, ratio1, ratio2, r)

      class(type_cheb) :: this

      !-- Input variables:
      integer,  intent(in) :: n_r_max
      real(cp), intent(in) :: ricb
      real(cp), intent(in) :: rcmb
      real(cp), intent(inout) :: ratio1
      real(cp), intent(in) :: ratio2

      !-- Output variable:
      real(cp), intent(out) :: r(n_r_max)

      !-- Local variables:
      integer :: n_r
      real(cp) :: lambd,paraK,paraX0 !parameters of the nonlinear mapping

      !--
      !-- There's possibly an issue when the Chebyshev mapping was used in
      !-- the old grid and a different mapping is used on the new one
      !--

      if ( this%l_map ) then
         this%alpha1=ratio1
         this%alpha2=ratio2
         paraK=atan(this%alpha1*(one+this%alpha2))/atan(this%alpha1*(one-this%alpha2))
         paraX0=(paraK-one)/(paraK+one)
         lambd=atan(this%alpha1*(one-this%alpha2))/(one-paraX0)
      else
         this%alpha1=0.0_cp
         this%alpha2=0.0_cp
      end if

      call cheb_grid(ricb,rcmb,n_r_max-1,r,this%r_cheb,this%alpha1,this%alpha2, &
           &         paraX0,lambd,this%l_map)

      if ( this%l_map ) then

         !-- Tangent mapping (see Bayliss et al. 1992)
         if ( index(map_function, 'TAN') /= 0 .or.      &
         &    index(map_function, 'BAY') /= 0 ) then

            do n_r=1,n_r_max
               this%drx(n_r)  =                         (two*this%alpha1) /        &
               &    ((one+this%alpha1**2*(two*r(n_r)-ricb-rcmb-this%alpha2)**2)*   &
               &    lambd)
               this%ddrx(n_r) =-(8.0_cp*this%alpha1**3*(two*r(n_r)-ricb-rcmb-      &
               &               this%alpha2)) / ((one+this%alpha1**2*(-two*r(n_r)+  &
               &               ricb+rcmb+this%alpha2)**2)**2*lambd)
            end do

         !-- Arcsin mapping (see Kosloff and Tal-Ezer, 1993)
         else if ( index(map_function, 'ARCSIN') /= 0 .or. &
         &         index(map_function, 'KTL') /= 0 ) then

            do n_r=1,n_r_max
               this%drx(n_r)  =two*asin(this%alpha1)/this%alpha1*sqrt(one-     &
               &               this%alpha1**2*this%r_cheb(n_r)**2)
               this%ddrx(n_r) =-four*asin(this%alpha1)**2*this%r_cheb(n_r)
            end do

         end if

      else !-- No mapping is used: this is the regular Gauss-Lobatto grid

         do n_r=1,n_r_max
            this%drx(n_r)  =two/(rcmb-ricb)
            this%ddrx(n_r) =0.0_cp
         end do

      end if

   end subroutine initialize_mapping
!------------------------------------------------------------------------------
   subroutine finalize(this, no_work_array)

      class(type_cheb) :: this

      !-- Input variable
      logical, optional, intent(in) :: no_work_array

      !-- Local variable
      logical :: l_work_array


      deallocate( this%rMat, this%drMat, this%d2rMat )
      deallocate( this%r_cheb, this%drx, this%ddrx)

      if ( present(no_work_array) ) then
         l_work_array=no_work_array
      else
         l_work_array=.false.
      end if

      call this%chebt%finalize(l_work_array)

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine cheb_grid(a,b,n,x,y,a1,a2,x0,lbd,l_map)
      !
      !   Given the interval [a,b] the routine returns the
      !   n+1 points that should be used to support a
      !   Chebychev expansion. These are the n+1 extrema y(i) of
      !   the Chebychev polynomial of degree n in the
      !   interval [-1,1].
      !   The respective points mapped into the interval of
      !   question [a,b] are the x(i).
      !
      !   .. note:: x(i) and y(i) are stored in the reversed order:
      !             x(1)=b, x(n+1)=a, y(1)=1, y(n+1)=-1
      !
       
      !-- Input variables
      real(cp), intent(in) :: a,b   ! interval boundaries
      integer,  intent(in) :: n ! degree of Cheb polynomial to be represented by the grid points
      real(cp), intent(in) :: a1,a2,x0,lbd
      logical,  intent(in) :: l_map ! Chebyshev mapping

      !-- Output variables
      real(cp), intent(out) :: x(*) ! grid points in interval [a,b]
      real(cp), intent(out) :: y(*) ! grid points in interval [-1,1]
       
      !-- Local variables:
      real(cp) :: bpa,bma
      integer :: k
       
      bma=half*(b-a)
      bpa=half*(a+b)

      do k=1,n+1
         y(k)=cos( pi*real(k-1,cp)/real(n,cp) )
         if ( l_map ) then
            if ( index(map_function, 'TAN') /= 0 .or.    &
            &   index(map_function, 'BAY') /= 0 ) then
               x(k)=half*(a2+tan(lbd*(y(k)-x0))/a1) + bpa
            else if ( index(map_function, 'ARCSIN') /= 0 .or. &
            &         index(map_function, 'KTL') /= 0 ) then
               x(k)=half*asin(a1*y(k))/asin(a1)+bpa
            end if
         else
            x(k)=bma * y(k) + bpa
         end if
      end do
        
   end subroutine cheb_grid
!------------------------------------------------------------------------------
   subroutine get_der_mat(this, n_r_max, l_cheb_coll)
      !
      !  Construct Chebychev polynomials and their first, second,
      !  and third derivative up to degree n_r at n_r points x
      !  in the interval [a,b]. Since the Chebs are only defined
      !  in [-1,1] we have to use a map, mapping the points x
      !  points y in the interval [-1,1]. This map is executed
      !  by the subroutine cheb_grid and has to be done
      !  before calling this program.
      !

      class(type_cheb) :: this

      !-- Input variables:
      integer, intent(in) :: n_r_max
      logical, intent(in) :: l_cheb_coll

      !-- Local variables:
      real(cp) :: dn2, monen
      integer :: n,k   ! counter

      if ( l_cheb_coll ) then
         !-- construction of chebs and derivatives with recursion:
         do k=1,n_r_max  ! do loop over the n_r grid points !

            !----- set first two chebs:
            this%rMat(1,k)  =one
            this%rMat(2,k)  =this%r_cheb(k)
            this%drMat(1,k) =0.0_cp
            this%drMat(2,k) =this%drx(k)
            this%d2rMat(1,k)=0.0_cp
            this%d2rMat(2,k)=this%ddrx(k)

            !----- now construct the rest with a recursion:
            do n=3,n_r_max ! do loop over the (n-1) order of the chebs

               this%rMat(n,k)  =    two*this%r_cheb(k)*this%rMat(n-1,k) - &
               &                                       this%rMat(n-2,k)
               this%drMat(n,k) =       two*this%drx(k)*this%rMat(n-1,k) + &
               &                   two*this%r_cheb(k)*this%drMat(n-1,k) - &
               &                                      this%drMat(n-2,k)
               this%d2rMat(n,k)=      two*this%ddrx(k)*this%rMat(n-1,k) + &
               &                     four*this%drx(k)*this%drMat(n-1,k) + &
               &                  two*this%r_cheb(k)*this%d2rMat(n-1,k) - &
               &                                     this%d2rMat(n-2,k)
            end do

         end do

         !-- This transposition is needed to bring those matrices in alignement
         !-- with the fortran column-major storage (see update routines)
         this%rMat  =transpose(this%rMat)
         this%drMat =transpose(this%drMat)
         this%d2rMat=transpose(this%d2rMat)

      else ! When Collocation is not used one just needs to store the tau-lines
           ! For each derivative the first row corresponds to the outer boundary
           ! the second row to the inner boundary

         do n=1,n_r_max
            dn2 = real(n-1,cp)*real(n-1,cp)
            this%rMat(1,n)  = one
            if ( mod(n,2)==1 ) then
               monen = one
            else
               monen = -one
            end if
            this%rMat(2,n)  = monen
            this%drMat(1,n) = two*dn2 ! Replace by drx maybe
            this%drMat(2,n) = two*dn2*(-one)*monen ! drx
            this%d2rMat(1,n)= four*third*dn2*(dn2-one) ! drx
            this%d2rMat(2,n)= four*third*dn2*(dn2-one)*monen
         end do

      end if

   end subroutine get_der_mat
!------------------------------------------------------------------------------
   subroutine costf1_complex_2d(this,f,nMstart,nMstop,n_r_max)

      class(type_cheb) :: this

      !-- Input variables:
      integer,  intent(in) :: nMstart
      integer,  intent(in) :: nMstop 
      integer,  intent(in) :: n_r_max            ! number of columns in f,f2
    
      !-- Output variables:
      complex(cp), intent(inout) :: f(nMstart:nMstop,n_r_max) ! data/coeff input
    
      call this%chebt%costf(f,nMstart,nMstop,n_r_max)
    
   end subroutine costf1_complex_2d
!------------------------------------------------------------------------------
   subroutine costf1_real_1d(this,f,n_r_max)

      class(type_cheb) :: this

      !-- Input variables:
      integer, intent(in) :: n_r_max

      !-- Output variables:
      real(cp), intent(inout) :: f(n_r_max)   ! data/coeff input

      call this%chebt%costf(f, n_r_max)

   end subroutine costf1_real_1d
!------------------------------------------------------------------------------
end module chebyshev
