      module conv_exact_tests
        use constants, only : omg_ => omega_qho
        use propagation
        use conv_tests
        implicit none
        private
        public :: exact_closed
        public :: exact_incremental

        integer, parameter :: wp = kind(1.0d0)
        complex(wp), parameter :: ci = (0.0_wp,1.0_wp)
      
      contains


!  ***   s = svec with svec = (1, 1, 1, ..., 1, 1)

!     subroutine exact_closed(nmax, eigval, omega, svec,      &
!                                          t, t0, psi0, phi0, psi, phi)
!     
!     implicit none
!     integer, intent(in) :: nmax
!     real(8), intent(in)  :: eigval(nmax)
!     complex(8), intent(in)  :: psi0(nmax)
!     complex(8), intent(in)  :: phi0(nmax)
!     complex(8), intent(in)  :: svec(nmax)
!     real(8),    intent(in)  :: omega
!     real(8),    intent(in)  :: t0, t
!     complex(8), intent(out) :: psi(nmax)
!     complex(8), intent(out) :: phi(nmax)
!     
!     integer :: k
!     complex(8) :: expA, expWt, expWt0, denom
!     real(8) :: dt
!     
!     dt = t - t0
!     
!     expWt  = exp(  (0d0,1d0) * omega * t )
!     expWt0 = exp(  (0d0,1d0) * omega * t0 )

!     do k = 1, nmax
!        expA = exp( -(0d0,1d0) * eigval(k) * dt )
!     
!        denom = eigval(k) + omega
!     
!        psi(k) = expA * psi0(k)    
!        phi(k) = expA * phi0(k)                                      &
!               - svec(k)/denom * ( expWt - expA * expWt0 )
!     
!     end do
!     
!     end subroutine



!     subroutine exact_incremental(nmax, eigval, svec, omega,  &
!                         dt, tn, psi_in, phi_in, psi_out, phi_out)
!     
!     implicit none
!     integer, intent(in) :: nmax
!     real(8), intent(in)    :: eigval(nmax)
!     real(8),    intent(in)    :: omega
!     real(8),    intent(in)    :: dt
!     real(8),    intent(in)    :: tn
!     complex(8), intent(in)    :: svec(nmax)
!     complex(8), intent(in) :: psi_in(nmax)
!     complex(8), intent(in) :: phi_in(nmax)
!     complex(8), intent(out) :: psi_out(nmax)
!     complex(8), intent(out) :: phi_out(nmax)
!     
!     integer :: n, k
!     real(8) :: tnp1
!     complex(8) :: phaseH, expWn, expWnp1, denom
!     
!     
!     tnp1 = tn + dt
!     expWn   = exp(  (0d0,1d0) * omega * tn )
!     expWnp1 = exp(  (0d0,1d0) * omega * tnp1 )
!     
!     do k = 1, nmax
!     
!        phaseH  = exp( -(0d0,1d0) * eigval(k) * dt )
!        denom = eigval(k) + omega
!     
!        psi_out(k) = phaseH * psi_in(k)
!        phi_out(k) = phaseH * phi_in(k)                           &
!               - svec(k)/denom * ( expWnp1 - phaseH * expWn )
!     
!     end do
!     
!     end subroutine


!  ***    s(t) = cos(Omega_1 t)e_1 + sin(Omega_2 t )e_2

      subroutine exact_closed(nmax, eigval, omega, svec,      &
                                           t, t0, psi0, phi0, psi, phi)
      
      implicit none
      integer, intent(in) :: nmax
      real(8), intent(in)  :: eigval(nmax)
      complex(8), intent(in)  :: psi0(nmax)
      complex(8), intent(in)  :: phi0(nmax)
      complex(8), intent(in)  :: svec(nmax)
      real(8),    intent(in)  :: omega
      real(8),    intent(in)  :: t0, t
      complex(8), intent(out) :: psi(nmax)
      complex(8), intent(out) :: phi(nmax)
      
      integer :: k
      real(8) :: alphap, alpham, betap, betam
      real(8) :: omega1, omega2
      real(8) :: denomp, denomm
      real(8) :: dt

      complex(8) :: ci
      complex(8) :: expA
      
      dt = t - t0
      ci = (0.d0, 1.d0)
      
      omega1 = omega 
      omega2 = 2.d0 * omega 

      do k = 1, nmax
         expA = exp( -ci * eigval(k) * dt )
      
         psi(k) = expA * psi0(k)    
         phi(k) = expA * phi0(k)                
      
      end do

      alphap = omega + omega1
      alpham = omega - omega1
      betap  = omega + omega2
      betam  = omega - omega2

      ! mode 1 (cos forcing)
      denomp = eigval(1) + alphap
      denomm = eigval(1) + alpham
      phi(1) = phi(1) - 1.d0/2 * svec(1) * (                           &
       ( exp(ci*alphap*t) - exp(ci*alphap*t0) * exp(-ci*eigval(1)*dt)) &
                  /   denomp  +                                        &
       ( exp(ci*alpham*t) - exp(ci*alpham*t0) * exp(-ci*eigval(1)*dt)) &
                  /   denomm                                           &
       )
      
      ! mode 2 (sin forcing)
      denomp = eigval(2) + betap
      denomm = eigval(2) + betam
      phi(2) = phi(2) + ci/2 * svec(2) * (                             &
       ( exp(ci*betap*t) - exp(ci*betap*t0) * exp(-ci*eigval(2)*dt))   &
                  /   denomp  -                                        &
       ( exp(ci*betam*t) - exp(ci*betam*t0) * exp(-ci*eigval(2)*dt))   &
                  /   denomm                                           &
       )

      end subroutine exact_closed

      subroutine exact_increment(nmax, eigval, omega, svec,      &
                              tn, dt, psi_in, phi_in, psi_out, phi_out)
      
      implicit none
      integer, intent(in) :: nmax
      real(8), intent(in) :: eigval(nmax)
      complex(8), intent(in) :: psi_in(nmax), phi_in(nmax)
      complex(8), intent(in) :: svec(nmax)
      real(8), intent(in) :: omega, tn, dt
      complex(8), intent(out) :: psi_out(nmax), phi_out(nmax)
      
      integer :: k
      real(8) :: omega1, omega2
      real(8) :: alphap, alpham, betap, betam
      real(8) :: denomp, denomm
      complex(8) :: ci, expA
      real(8) :: tnp1
      
      ci = (0.d0,1.d0)
      tnp1 = tn + dt
      
      omega1 = omega
      omega2 = 2.d0*omega
      
      ! propagate homogeneous part
      do k = 1, nmax
         expA = exp(-ci*eigval(k)*dt)
         psi_out(k) = expA * psi_in(k)
         phi_out(k) = expA * phi_in(k)
      end do
      
      alphap = omega + omega1
      alpham = omega - omega1
      betap  = omega + omega2
      betam  = omega - omega2
      
      ! mode 1 (cos forcing)
      denomp = eigval(1) + alphap
      denomm = eigval(1) + alpham
      
      phi_out(1) = phi_out(1) - 0.5d0 * svec(1) * (                  &
         ( exp(ci*alphap*tnp1) -                                     &
               exp(ci*alphap*tn) * exp(-ci*eigval(1)*dt) ) / denomp  &
       + ( exp(ci*alpham*tnp1) -                                     &
               exp(ci*alpham*tn) * exp(-ci*eigval(1)*dt) ) / denomm  &
         )
      
      ! mode 2 (sin forcing)
      denomp = eigval(2) + betap
      denomm = eigval(2) + betam
      
      phi_out(2) = phi_out(2) + ci*0.5d0 * svec(2) * (               &
         ( exp(ci*betap*tnp1) -                                      &
               exp(ci*betap*tn) * exp(-ci*eigval(2)*dt) ) / denomp   &
       - ( exp(ci*betam*tnp1) -                                      &
               exp(ci*betam*tn) * exp(-ci*eigval(2)*dt) ) / denomm   &
         )
      
      end subroutine exact_increment




      end module conv_exact_tests
