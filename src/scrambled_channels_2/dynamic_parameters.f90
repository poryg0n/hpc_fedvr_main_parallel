      module dynamic_parameters
        use constants, only : ppi
        implicit none
        public :: force_t, set_force_params
      
        real(8), public :: f0        ! field amplitude
        real(8), public :: omega0    ! field frequency
        real(8), public :: pfai      ! field cep phase
        real(8), public :: period    ! field period
        real(8), public :: tau       ! pulse active duration
  
        real(8), public :: t_end
        real(8), public :: t_ini
        real(8), public :: trange
        real(8), public :: dt0
        real(8), public :: t
  
        integer, public :: nt        ! number of time steps
        integer, public :: noc
        integer, public :: ntau

        integer, public :: src_type
        integer, public :: run
        integer, public :: need
        integer, public :: nch
        real(8), public :: omg_max
        real(8), public :: omg_min
        real(8), public :: wstep
        real(8), public :: wrange

        integer, public :: order

        logical, public :: do_time_obs
        integer, public :: obs_stride

      
      contains

  
        subroutine init_time_grid(noc_, ntau_, nsteps)
  
           integer, intent(in) :: nsteps, noc_, ntau_
 
           noc   = noc_
           ntau  = ntau_
           nt    = nsteps

           tau    = 2.d0*ppi*noc/omega0
           trange = ntau * tau
           t_end  =  trange/2.d0
           t_ini  = -t_end

           dt0   = trange/nt
  
        end subroutine


        subroutine init_src(src_type_, nch_,                          &
                                omg_0, omg_1, wstep_, run_)
  
           integer, intent(in) :: src_type_
           integer, intent(in) :: run_
           integer, intent(in) :: nch_
           real(8), intent(in) :: omg_1, omg_0
           real(8), intent(in) :: wstep_
           integer :: need 

  
           src_type =  src_type_
           nch = nch_
           omg_min = omg_0
           omg_max = omg_1
           wstep = wstep_
           run = run_

           wrange = (omg_max-omg_min)
  
           need = ceiling(wrange/(wstep*(nch-1)))

        end subroutine


        subroutine set_resolution_order(order_)
  
           integer, intent(in) :: order_
     
           order =  order_
  
        end subroutine

        subroutine set_other_dyn_params(do_time_obs_, obs_stride_)
  
           logical, intent(in) :: do_time_obs_
           integer, intent(in) :: obs_stride_
     
           do_time_obs  =  do_time_obs_
           obs_stride   =  obs_stride_
  
        end subroutine


        logical function advance_time(t_out)
          ! Optional argument: user can ask for current time
          real(8), intent(out), optional :: t_out
      
          t = t + dt0
          advance_time = (t <= t_end)
      
          if (present(t_out)) t_out = t
        end function
      
        function current_time() result(t_now)
          ! Just in case you want to peek without advancing
          real(8) :: t_now
          t_now = t
        end function
      

      
        subroutine set_force_params(f0_in, omega_in, pfai_in)
          real(8), intent(in) :: f0_in, omega_in, pfai_in
          f0    = f0_in
          omega0 = omega_in
          pfai  = pfai_in
        end subroutine

      
        real(8) function force_t(tt)
          real(8), intent(in) :: tt

          force_t = f0 * exp(-(2*tt/tau)**2) * cos(omega0 * tt)
        end function



        subroutine plot_force(unit, t1, t0, dt)
          implicit none
          real(8), intent(in) :: t0, t1, dt
          integer, intent(in) :: unit
          real(8) :: tt
        
          tt = t0
          do while (tt <= t1)
             write(unit,'(2e20.10)') tt, force_t(tt)
             tt = tt + dt
          end do
        end subroutine


        subroutine print_dynamic_parameters()
  
           print*, "----- Dynamic parameters -----"
           print*, "field f0    :", f0
           print*, "omega0      :", omega0
           print*, "pfai_       :", pfai
           print*, "final time  :", t_end
           print*, "init. time  :", t_ini
           print*, "nsteps      :", nt
           print*, "dt0         :", dt0
           print*, "noc         :", noc
           print*, "ntau        :", ntau
           print*, "             "
           print*, "source_type :", src_type

           print*, "nchan         = ", nch
           print*, "omg_max       = ", omg_max
           print*, "omg_min       = ", omg_min
           print*, "order       = ", order
  
        end subroutine


      end module dynamic_parameters
