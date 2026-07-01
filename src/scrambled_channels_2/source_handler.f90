      module source_handler
        use iso_fortran_env, only : dp => real64
        implicit none
      contains
        subroutine build_source(nmax, t, omega, tc, tau, phi, src)
          implicit none
          integer, intent(in) :: nmax
          real(dp), intent(in) :: t, omega, tc, tau
          complex(dp), intent(in) :: phi(nmax)
          complex(dp), intent(out) :: src(nmax)
        
          real(dp) :: g
        
          g = sin(omega*t) * exp( -((t-tc)/tau)**2 )
          src = g * phi
        end subroutine

      end module source_handler

