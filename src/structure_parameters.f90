      module structure_parameters
      implicit none
      private
      
      !-----------------------------
      ! Public variables
      !-----------------------------
      
      integer, public :: snbr      ! number of elements
      integer, public :: nnbr      ! lobatto points per element
      integer, public :: nmax      ! lobatto points per element
      
      real(8), public :: xmin
      real(8), public :: xmax
      real(8), public :: xc
      real(8), public :: dx
      real(8), public :: xrange
      real(8), public :: jac
      real(8), public :: inv_jac

      logical, public :: qho
      
      
      public :: set_grid_params
      public :: print_structure_parameters
      
      
      contains
      
      
        subroutine set_grid_params(np, ns, xc_in, xmin_in, xmax_in)
        
        integer,intent(in) :: ns, np
        real(8),intent(in) :: xmin_in, xmax_in, xc_in
        
        nnbr = np
        snbr = ns
        nmax = snbr*(nnbr-1)-1
        
        xmin = xmin_in
        xmax = xmax_in
        xc =  xc_in
        
        xrange = xmax - xmin
        dx     = xrange / snbr
        jac    = 0.5d0 * dx        ! = xmax / snbr
        inv_jac = 1.d0 / jac
        
        end subroutine
      
      
      

        subroutine print_structure_parameters()
        
        print*, "----- Static parameters -----"
        print*, "elements  :", snbr
        print*, "nodes     :", nnbr
        print*, "nbr. pts  :", nmax
        print*, "xmin      :", xmin
        print*, "xmax      :", xmax
        
        end subroutine
      
      
      end module
