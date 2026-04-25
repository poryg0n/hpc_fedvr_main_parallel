      module math_util
      use iso_fortran_env, only : dp => real64
      implicit none
      private
      
      public :: differentiate
      public :: differentiate4
      public :: apply_momentum_operator
      public :: composite_simpson_uniform
      public :: composite_simpson_18
      public :: composite_simpson_18c
      public :: integr_over_range
      public :: varkap
      public :: sgn
      
      interface differentiate
         module procedure differentiate_real
         module procedure differentiate_complex
      end interface
      
      interface differentiate4
         module procedure differentiate4_real
         module procedure differentiate4_complex
      end interface

      contains

      subroutine differentiate_real(x, f, df)
        implicit none
      
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: f(:)
        real(dp), intent(out) :: df(:)
        
        integer :: n, i
        real(dp) :: dx1, dx2

        n = size(x)

      
        if (n == 1) then
           df(1) = (0.0d0,0.0d0)
           return
        end if

      
        ! --- forward difference
        dx1 = x(2)-x(1)
        df(1) = (f(2)-f(1)) / dx1
      
        ! --- central differences
        do i=2,n-1
           dx1 = x(i)-x(i-1)
           dx2 = x(i+1)-x(i)
           df(i) = 0.5d0 * ( (f(i+1)-f(i))/dx2 + (f(i)-f(i-1))/dx1 )
        end do
      
        ! --- backward difference
        dx1 = x(n)-x(n-1)
        df(n) = (f(n)-f(n-1)) / dx1
      
      end subroutine differentiate_real


      subroutine differentiate_complex(x, f, df)
        implicit none
      
        real(dp), intent(in) :: x(:)
        complex(dp), intent(in) :: f(:)
        complex(dp), intent(out) :: df(:)
        
        integer :: n, i
        real(dp) :: dx1, dx2

        n = size(x)
      
      
        if (n == 1) then
           df(1) = (0.0d0,0.0d0)
           return
        end if
      
        ! --- forward difference
        dx1 = x(2)-x(1)
        df(1) = (f(2)-f(1)) / dx1
      
        ! --- central differences
        do i=2,n-1
           dx1 = x(i)-x(i-1)
           dx2 = x(i+1)-x(i)
           df(i) = 0.5d0 * ( (f(i+1)-f(i))/dx2 + (f(i)-f(i-1))/dx1 )
        end do
      
        ! --- backward difference
        dx1 = x(n)-x(n-1)
        df(n) = (f(n)-f(n-1)) / dx1
      
      end subroutine differentiate_complex


      subroutine differentiate4_real(x, f, df)
      implicit none
      
      real(dp), intent(in) :: x(:)
      real(dp), intent(in) :: f(:)
      real(dp), intent(out) :: df(:)
      
      integer :: n, i
      real(dp) :: h

      n = size(x)
      
      if (n < 5) then
         print *, 'differentiate4: need at least 5 points'
         stop
      end if
      
      h = x(2) - x(1)
      
      ! ---- forward edges ----
      df(1) = (-25.d0*f(1) + 48.d0*f(2) - 36.d0*f(3) +                &
           16.d0*f(4) - 3.d0*f(5)) / (12.d0*h)
    
      df(2) = (-3.d0*f(1) - 10.d0*f(2) + 18.d0*f(3) -                 & 
           6.d0*f(4) + f(5)) / (12.d0*h)
      
      ! ---- interior stencil ----
      do i = 3, n-2
         df(i) = (-f(i+2) + 8.d0*f(i+1) - 8.d0*f(i-1) + f(i-2))       &
                       / (12.d0*h)
      end do
      
      ! ---- backward edges ----
      df(n-1) = (-f(n-4) + 6.d0*f(n-3) - 18.d0*f(n-2) +               &
           10.d0*f(n-1) + 3.d0*f(n)) / (12.d0*h)
      
      df(n) = (3.d0*f(n-4) - 16.d0*f(n-3) + 36.d0*f(n-2) -            &
           48.d0*f(n-1) + 25.d0*f(n)) / (12.d0*h)
      
      end subroutine differentiate4_real


      subroutine differentiate4_complex(x, f, df)
      implicit none
      
      real(dp), intent(in) :: x(:)
      complex(dp), intent(in) :: f(:)
      complex(dp), intent(out) :: df(:)
      
      integer :: n, i
      real(dp) :: h
      
      n = size(x)

      if (n < 5) then
         print *, 'differentiate4: need at least 5 points'
         stop
      end if
      
      h = x(2) - x(1)
      
      ! ---- forward edges ----
      df(1) = (-25.d0*f(1) + 48.d0*f(2) - 36.d0*f(3) +                &
           16.d0*f(4) - 3.d0*f(5)) / (12.d0*h)
    
      df(2) = (-3.d0*f(1) - 10.d0*f(2) + 18.d0*f(3) -                 & 
           6.d0*f(4) + f(5)) / (12.d0*h)
      
      ! ---- interior stencil ----
      do i = 3, n-2
         df(i) = (-f(i+2) + 8.d0*f(i+1) - 8.d0*f(i-1) + f(i-2))       &
                       / (12.d0*h)
      end do
      
      ! ---- backward edges ----
      df(n-1) = (-f(n-4) + 6.d0*f(n-3) - 18.d0*f(n-2) +               &
           10.d0*f(n-1) + 3.d0*f(n)) / (12.d0*h)
      
      df(n) = (3.d0*f(n-4) - 16.d0*f(n-3) + 36.d0*f(n-2) -            &
           48.d0*f(n-1) + 25.d0*f(n)) / (12.d0*h)
      
      end subroutine differentiate4_complex


      subroutine apply_momentum_operator(nmax, eigvec, xx, wx,        &
                          jac, psi_in, psi_out, method)
      
      implicit none
      
      integer, intent(in) :: nmax
      integer, intent(in) :: method
      
      complex*16, intent(in)  :: psi_in(nmax)
      complex*16, intent(out) :: psi_out(nmax)
      
      real*8, intent(in) :: eigvec(nmax,nmax)
      real*8, intent(in) :: xx(nmax)
      real*8, intent(in) :: wx(nmax)
      real*8, intent(in) :: jac
      
      complex*16 :: psi_x(nmax)
      complex*16 :: dpsi_x(nmax)
      
      integer :: n
      complex*16, parameter :: ci=(0d0,1d0)
      
      select case(method)
      
      !-----------------------------------------
      ! METHOD 1
      ! rotate -> differentiate -> rotate back
      !-----------------------------------------
      case(0)
      
         psi_x = matmul(eigvec,psi_in)
         psi_x = psi_x / wx / dsqrt(jac)
      
         call differentiate(xx,psi_x,dpsi_x)
      
         dpsi_x = -ci*dpsi_x
      
         psi_out = matmul(transpose(eigvec),dpsi_x)
      
      !-----------------------------------------
      ! METHOD 2
      ! spectral HO momentum operator
      !-----------------------------------------
      case(1)
      
         psi_out = 0d0
      
         do n=1,nmax
      
            if(n>1) then
               psi_out(n) = psi_out(n)                                &
                       + ci/sqrt(2d0) * sqrt(dble(n-1)) * psi_in(n-1)
            endif
      
            if(n<nmax) then
               psi_out(n) = psi_out(n)                                &
                       - ci/sqrt(2d0) * sqrt(dble(n)) * psi_in(n+1)
            endif
      
         enddo
      
      end select
      
      end subroutine apply_momentum_operator


      subroutine composite_simpson_uniform(n, h, f, integral)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: h
        real(8), intent(in) :: f(n)
        real(8), intent(out) :: integral
      
        integer :: i, n_use
      
        ! --- enforce odd number of points ---
        if (mod(n,2) == 0) then
!          print*, "Warning: Simpson requires odd n, \\dropping last point"
           n_use = n - 1
        else
           n_use = n
        end if
      
        if (n_use < 3) then
           integral = 0.d0
           return
        end if
      
        integral = f(1) + f(n_use)
      
        do i = 2, n_use-1, 2
           integral = integral + 4.d0 * f(i)
        end do
      
        do i = 3, n_use-2, 2
           integral = integral + 2.d0 * f(i)
        end do
      
        integral = integral * h / 3.d0
      
      end subroutine



      function sgn(x) result(sgnx)
        implicit none
        real(8), intent(in) :: x
        real(8) :: sgnx
        
        if (x.lt.0.d0) then
           sgnx = -1.d0
        else 
           sgnx = +1.d0
        end if  

        return
      end function


      function varkap(kap,omega) result(vk)
        implicit none
        real(8), intent(in) :: omega
        real(8), intent(in) :: kap
        real(8) :: vk
        
        vk = dsqrt(kap**2 + 2.d0*omega)     
        return
      end function




      subroutine composite_simpson_18c(ndim, xx, df, res, ff)
      implicit none
      integer :: i, ndim, nn, nm, even
      real*8 :: hh
      complex*16 :: res, res1
      real*8, dimension(*) :: xx
      complex*16, dimension(*) :: df

      real*8, dimension(4) :: xxl
      complex*16, dimension(4) :: dfl
      complex*16, optional, dimension(*) :: ff

      nn=ndim-1
      if (mod(nn,2).eq.0) then
         even = 1 
!        write(*,*) "N is even", nn
      else
!        write(*,*) "Number of subintervals is not even"
         even = 0
         nn=nn-3
      end if

      hh = (xx(2)-xx(1))
!     write(*,*) hh, nn
      
      res = 0.d0
      res1 = 0.d0
      do i=1,nn/2
         res = res +                                                  &
                hh/3*( df(2*i-1) + 4.d0*df(2*i) + df(2*i+1) )
!        write(*,*) nn/2, i, hh, res
      enddo
!      if (even.eq.0) then
!!        write(*,*) "Not even even"
!         do i=1,4
!            xxl(i)=xx(ndim-4+i)
!            dfl(i)=df(ndim-4+i)
!         enddo
!         call composite_simpson_38c(4, xxl, dfl, res1) 
!      end if
      res = res + res1
!     write(*,*) res

      if(present(ff)) then
         ff(1)=df(1)
         do i=1,ndim
            ff(i+1) = ff(i) + .5d0*(df(i+1)+df(i))*(xx(i+1)-xx(i))
         enddo
      end if
      end


      subroutine composite_simpson_18(ndim, xx, df, res, ff)
      implicit none
      integer :: i, ndim, nn, nm, even
      real*8 :: hh, res, res1
      real*8, dimension(*) :: xx, df
      real*8, dimension(4) :: xxl, dfl
      real*8, optional, dimension(*) :: ff


      nn=ndim-1
      if (mod(nn,2).eq.0) then
         even = 1 
!        write(*,*) "N is even", nn
      else
!        write(*,*) "Number of subintervals is not even"
         even = 0
         nn=nn-3
      end if

      hh = (xx(2)-xx(1))
!     write(*,*) hh, nn
      
      res = 0.d0
      res1 = 0.d0
      do i=1,nn/2
         res = res +                                                  &
                hh/3*( df(2*i-1) + 4.d0*df(2*i) + df(2*i+1) )
!        write(*,*) nn/2, i, hh, res
      enddo
!      if (even.eq.0) then
!!        write(*,*) "Not even even"
!         do i=1,4
!            xxl(i)=xx(ndim-4+i)
!            dfl(i)=df(ndim-4+i)
!         enddo
!         call composite_simpson_38(4, xxl, dfl, res1) 
!      end if
      res = res + res1
!     write(*,*) res

      if(present(ff)) then
         ff(1)=df(1)
         do i=1,ndim
            ff(i+1) = ff(i) + .5d0*(df(i+1)+df(i))*(xx(i+1)-xx(i))
         enddo
      end if
      end


      subroutine integr_over_range(krange, kk, vec_k, res)
         implicit none
         integer, intent(in) :: krange
         real(8), intent(in) :: kk(krange)
         complex(8), intent(in) :: vec_k(krange)
         complex(8), intent(out) :: res

         integer :: n_cont, j
         real(8) :: res1, res2
         real(8) :: auxr(krange/2)
         real(8) :: auxc(krange/2)

         n_cont = krange/2

        do j=1,n_cont
           auxr(j) = kk(j)
           auxc(j) = abs(vec_k(j))**2
        enddo
        call composite_simpson_18(n_cont, auxr, auxc, res1)

        do j=1,n_cont
           auxr(j) = kk(j+n_cont)
           auxc(j) = abs(vec_k(j+n_cont))**2
        enddo
        call composite_simpson_18(n_cont, auxr, auxc, res2)

        res = res1 + res2

      end subroutine

      end module
