      module constants
        use iso_fortran_env, only : dp => real64
        implicit none
        real(dp), parameter :: hbar = 1.d0           ! Planck’s constant in a.u.
        real(dp), parameter :: mass = 1.d0           ! Mass (set to 1 for QHO a.u.)
        real(dp), parameter :: omega_qho = 1.d0      ! Oscillator frequency

        real(dp), parameter :: kapp = 1.d0           ! kappa

! ***   Maths constant
        real(dp), parameter :: ppi = 4.d0*datan(1.d0)
        complex(dp), parameter :: c0 = (0.d0, 0.d0)
        complex(dp), parameter :: c1 = (1.d0, 0.d0)
        complex(dp), parameter :: ci = (0.d0, 1.d0)
      end module

