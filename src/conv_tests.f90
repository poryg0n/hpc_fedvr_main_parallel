      module conv_tests
        use constants, only : omg_ => omega_qho
        use constants, only : ppi
        use dynamic_parameters
        use propagation
!       use force_field
!       use time_grid
!       use source_handler
        use fedvr_derivative_ops
        implicit none
        private
        public :: exact_closed_duhamel
        public :: exact_increment_duhamel
        public :: test_against_exact_solution
        public :: test_richardson_inhomogeneous

        integer, parameter :: wp = kind(1.0d0)
        complex(wp), parameter :: ci = (0.0_wp,1.0_wp)
      
      contains


!         ***  Richardson tests ******


      subroutine test_richardson_inhomogeneous( log_unit,            &
                                           nmax, lnbr, nnbr,   &
                                           jac,                   &
                                           xs, xx, wx,                &
                                           map, Dref,                 &
                                           nsteps, omega0_,            &
                                           eigval, eigvec,            &
                                           psi0, phi0,               &
                                           src_type, omega, order)
      
         implicit none
         integer, intent(in) :: nmax, lnbr, nnbr, order
         integer, intent(in) :: src_type, log_unit
         integer, intent(in) :: nsteps
         real(8), intent(in) :: omega0_
         real(8), intent(in) :: omega, jac
         real(8), intent(in) :: xs(:), xx(:), wx(:)
         integer, intent(in) :: map(:,:)
         real(8), intent(in) :: Dref(:,:), eigval(:), eigvec(:,:)
         complex(8), intent(in) :: psi0(:), phi0(:)
      
         integer, parameter :: nrun = 4
         real(8) :: dt(nrun)
         real(8) :: t0, t1, dt0, trange, tau
         complex(8), allocatable :: psi(:,:), phi(:,:)
         real(8) :: errL2_psi(nrun), errL2_phi(nrun)
         real(8) :: errMax_psi(nrun), errMax_phi(nrun)
         real(8) :: norm_psi(nrun), norm_phi(nrun)
         real(8) :: rel_psi(nrun), rel_phi(nrun)
         real(8) :: order_psi, order_phi

         real(8) :: err_psi(nrun-1), err_phi(nrun-1)
         integer :: k
      
         allocate(psi(nmax,nrun), phi(nmax,nrun))

         tau    = 5 * 2.d0*ppi/omega0
         trange = 3 * tau
         t1     =  trange/2.d0
         t0     = -t1

         dt0    = trange/nsteps

      
         ! --- time steps ---
         do k = 1, nrun
            dt(k) = dt0 / 2.0d0**(k-1)
         end do
      
         ! --- run propagations ---
         do k = 1, nrun
            call run_inhomogeneous_duhamel_driver( log_unit,       &
                                     nmax, lnbr, nnbr,             &
                                     jac,                          &
                                     xs, xx, wx,                   &
                                     map, Dref,                    &
                                     dt(k), t1, t0,                &
                                     eigval, eigvec,               &
                                     psi0, phi0,                   &
                                     psi(:,k), phi(:,k),           &
                                     src_type, omega, order)
         end do


         ! --- errors ---
         do k = 1, nrun-1
            err_psi(k) = sqrt(sum(abs(psi(:,k) - psi(:,k+1))**2))
            err_phi(k) = sqrt(sum(abs(phi(:,k) - phi(:,k+1))**2))
         end do
      
         ! --- observed order (last two refinements) ---
         order_psi = log(err_psi(nrun-2)/err_psi(nrun-1)) / log(2.0d0)
         order_phi = log(err_phi(nrun-2)/err_phi(nrun-1)) / log(2.0d0)
      
         write(log_unit,'(a)') '--- Richardson convergence test ---'

         if (order.eq.2) then
            write(log_unit,*) "order 2 duhamel test"
         else
            write(log_unit,*) "order 4 duhamel test"
         end if

         do k = 1, nrun-1
            write(log_unit,'(a,i1,2e16.8)') 'dt level ',               &
                                          k, err_psi(k), err_phi(k)
         end do

         write(log_unit,'(a,2f8.4)') 'Observed order (psi, phi): ',    &
                                               order_psi, order_phi
      
         deallocate(psi, phi)

      
      end subroutine test_richardson_inhomogeneous



!     ***     Absolute convergence tests     ***

      subroutine test_against_exact_solution( log_unit,                &
                                           nmax, lnbr, nnbr,        &
                                           jac,                   &
                                           xs, xx, wx,                &
                                           map, Dref,                 &
                                             dt0, tf, t0,        &
                                             eigval, eigvec,          &
                                             psi0, phi0,              &
                                            src_type, omega, order)
      
         implicit none
      
         integer, intent(in) :: nmax, lnbr, nnbr
         integer, intent(in) :: order, src_type, log_unit
         real(8), intent(in) :: t0, tf, dt0, omega, jac
         real(8), intent(in) :: xs(:), xx(:), wx(:)
         integer, intent(in) :: map(:,:)
         real(8), intent(in) :: Dref(:,:), eigval(:), eigvec(:,:)
         complex(8), intent(in) :: psi0(:), phi0(:)
      
         integer, parameter :: nrun = 4
         real(8) :: dt(nrun)
         complex(8), allocatable :: psi(:,:), phi(:,:)
         complex(8), allocatable :: psi_exact(:), phi_exact(:), svec(:)
         real(8) :: errL2_psi(nrun), errL2_phi(nrun)
         real(8) :: errMax_psi(nrun), errMax_phi(nrun)
         real(8) :: norm_psi(nrun), norm_phi(nrun)
         real(8) :: rel_psi(nrun), rel_phi(nrun)
         real(8) :: order_psi, order_phi

         integer :: k
         real(8) :: err_psi(nrun), err_phi(nrun)
      
         allocate(psi(nmax,nrun), phi(nmax,nrun))
         allocate(psi_exact(nmax), phi_exact(nmax))
         allocate(svec(nmax))
      
         !--------------------------------------------
         ! Time steps
         !--------------------------------------------
         do k = 1, nrun
            dt(k) = dt0 / 2.0d0**(k-1)
         end do


         !--------------------------------------------
         ! Compute exact solution once
         !--------------------------------------------

         call exact_closed_duhamel(nmax, tf, t0,                       &
                eigval, psi0, phi0, psi_exact, phi_exact, src_type, omega )
      
         !--------------------------------------------
         ! Run numerical propagators
         !--------------------------------------------
            do k = 1, nrun
               call run_inhomogeneous_duhamel_driver( log_unit,       &
                                           nmax, lnbr, nnbr,   &
                                           jac,                   &
                                           xs, xx, wx,                &
                                           map, Dref,                 &
                                           dt(k), tf, t0,             &
                                           eigval, eigvec,            &
                                           psi0, phi0,                &
                                           psi(:,k), phi(:,k),        &
                                           src_type, omega, order )


            end do
      
         !--------------------------------------------
         ! Compute error vs exact
         !--------------------------------------------
         do k = 1, nrun
            err_psi(k) = sqrt(sum(abs(psi(:,k) - psi_exact(:))**2))
            err_phi(k) = sqrt(sum(abs(phi(:,k) - phi_exact(:))**2))
         end do
      
         !--------------------------------------------
         ! Observed order (last two refinements)
         !--------------------------------------------
         order_psi = log(err_psi(nrun-1)/err_psi(nrun)) / log(2.0d0)
         order_phi = log(err_phi(nrun-1)/err_phi(nrun)) / log(2.0d0)
      
         write(log_unit,'(a)') '--- Convergence vs Exact Solution ---'
         do k = 1, nrun
            write(log_unit,'(a,i1,2e16.8)') 'dt level ',               &
                                             k, err_psi(k), err_phi(k)
         end do
         write(log_unit,'(a,2f8.4)') 'Observed order (psi, phi): ',    &
                                           order_psi, order_phi


         deallocate(psi, phi, psi_exact, phi_exact)
      
      end subroutine



      

!  *** One loop full run drivers 

      subroutine run_split_operator_driver( nmax, xx, dt, tf, t0,     &
                               eigval, eigvec,                        & 
                               psi0, phi0,                            &
                               psi, phi, order )

         integer, intent(in)    :: nmax, order
         real(8), intent(in)    :: t0, tf, dt
         real(8), intent(in)    :: xx(nmax), eigval(nmax)
         real(8), intent(in)    :: eigvec(nmax,nmax)

         complex(8), intent(in) :: psi0(nmax), phi0(nmax)
         complex(8), intent(out):: psi(nmax), phi(nmax)

         real(8)    :: t, tmid, fmid
         real(8)    :: vext_mid(nmax)
         complex(8) :: psi_in(nmax)
         complex(8) :: phi_in(nmax)
         integer    :: nt, i
      
         psi_in = psi0
         phi_in = phi0
         t      = t0
         nt     = int((tf - t0)/dt)

      
         write(*,'(a)') '# t  Re[psi(1)] Im[psi(1)] ||psi||          '
         do i = 1, nt
      
            call split_operator( nmax, dt, t, xx, eigval, eigvec,      &
                                      psi_in, psi, order )
            call split_operator( nmax, dt, t, xx, eigval, eigvec,      &
                                      phi_in, phi, order )
      
            psi_in = psi
            phi_in = phi
            t      = t + dt

            ! === diagnostics ===
            if (mod(i,10) == 0) then
               write(*,'(f8.4,3x,6e16.8)') t,                       &
                    real(psi(1)), aimag(psi(1)),                    &
                                         sqrt(sum(abs(psi)**2)),    &
                    real(phi(1)), aimag(phi(1)),                    &
                                         sqrt(sum(abs(phi)**2))
            end if


         end do

      end subroutine run_split_operator_driver





      subroutine run_inhomogeneous_duhamel_driver( log_unit,           &
                                 nmax, lnbr, nnbr,                     &
                                 jac,                                  &
                                 xs, xx, wx,                           &
                                 map, Dref,                            &
                                 dt, tf, t0,                           &
                                 eigval, eigvec,                       &
                                 psi0, phi0,                           &
                                 psi, phi,                             &
                                 src_type, omega, order)
      
         implicit none
         integer, intent(in) :: nmax, lnbr, nnbr
         integer, intent(in) :: order, src_type, log_unit
         real(8), intent(in) :: t0, tf, dt, omega, jac
         integer, intent(in) :: map(:,:)
         real(8), intent(in) :: xs(:), xx(:), wx(:)
         real(8), intent(in) :: Dref(:,:), eigval(:), eigvec(:,:)
         complex(8), intent(in) :: psi0(:), phi0(:)
      
         complex(8), intent(out) :: psi(:), phi(:)

         real(8) :: t, tmid, nt
         real(8) :: fmid
         real(8) :: vext_mid(nmax)
         complex(8), dimension(nmax) :: psi_in, phi_in, src
         complex(8) :: svec(nmax,3)
         integer :: i, j, k
      

! --- initial condition ---
       t = t0  ! start time
       psi_in = psi0           ! initial wavefunction
       phi_in = phi0

       nt = int((tf - t0)/dt)

       j=1
       if (src_type.eq.3) j=2

!      write(*,'(a)') '         Inside of convergence test         '
       write(log_unit,*) 'nt = ', nt
       write(log_unit,'(a)') '# t  Re[psi(1)] Im[psi(1)] ||psi||  '// & 
                      '                   Re[phi(2)] Im[phi(2)] ||phi||'
             write(log_unit,'(f10.4,3x,6e16.8)') t,               &
                  real(psi_in(1)), aimag(psi_in(1)),              &
                                     sqrt(sum(abs(psi_in)**2)),   &
                  real(phi_in(j)), aimag(phi_in(j)),              &
                                       sqrt(sum(abs(phi_in)**2))


 
       ! ------------------------------------------------------------
       ! Main time loop
       ! ------------------------------------------------------------
       do i = 1, nt

          ! === coupled inhomogeneous step ===

         call process_src_ingredients ( nmax, lnbr, nnbr,              &
                                 jac,                   &
                                 xs, xx, wx,                &
                                 map, Dref,                 &
                                 dt, t,                             &
                                 eigval, eigvec, psi_in, svec,         &
                                 src_type, order)

         call build_source_quadrature ( nmax, lnbr, nnbr,             &
                                               xs, xx, map, Dref,     &
                                               dt, t,                 &
                                               eigval, eigvec,        &
                                               svec,                  &
                                               src, omega,           &
                                               order )

         call split_operator(nmax, dt, t, xx, eigval, eigvec,         &
                                                   psi_in, psi, order)
         call split_operator(nmax, dt, t, xx, eigval, eigvec,         &
                                                   phi_in, phi, order)

!        !--------------------------------------------
         ! Add source contribution
         !--------------------------------------------
         phi = phi - ci * src

         ! === state update ===
         psi_in = psi
         phi_in = phi

         t     = t + dt

         ! === diagnostics ===
         if (mod(i,10) == 0) then
            write(log_unit,'(f10.4,3x,6e16.8)') t,               &
                 real(psi(1)), aimag(psi(1)),                    &
                                      sqrt(sum(abs(psi)**2)),    &
                 real(phi(j)), aimag(phi(j)),                    &
                                      sqrt(sum(abs(phi)**2))
         end if


      end do

      end subroutine run_inhomogeneous_duhamel_driver



! **********  Exact test subroutines for simple cases ***********


!!  ***    S(t) = exp(i omega t) F(x,t)

       subroutine exact_closed_duhamel(nmax, t,t0,eigval,    &
                                         psi0, phi0,              &
                                         psi_out, phi_out, src_type, omega)
       
       implicit none
       
       integer, intent(in) :: nmax, src_type
       real(8), intent(in) :: t,t0,omega
       real(8), intent(in) :: eigval(nmax)
       
       complex(8), intent(in)  :: psi0(nmax),phi0(nmax)
       complex(8), intent(out) :: psi_out(nmax),phi_out(nmax)
       
       integer :: i, n
       complex(8) :: evol, delta
       real(8) :: dt, En, Em, phase
       
       dt = t - t0
       
       !-----------------------------------------
       ! homogeneous propagation
       !-----------------------------------------
       do n = 1,nmax
          En = eigval(n)
          psi_out(n) = exp(-ci*En*dt)* psi0(n)
          phi_out(n) = exp(-ci*En*dt)* phi0(n)
       enddo
       
       !-----------------------------------------
       ! duhamel term
       !-----------------------------------------
       select case(src_type)
       
       !-----------------------------------------
       ! 1 : constant forcing  F(t) = psi0
       !-----------------------------------------
       case(1)
       
          do n = 1,nmax
       
             En = eigval(n)
             evol = exp(-ci*En*dt)
       
             delta = (exp(ci*omega*t) - exp(ci*omega*t0)*evol) &
                     / (En + omega)
       
             phi_out(n) = phi_out(n) - delta * psi0(n)
       
          enddo
       
       
       !-----------------------------------------
       ! 2 : stationary forcing  F(t) = psi(t)
       !-----------------------------------------
       case(2)
       
          delta = (exp(ci*omega*t) - exp(ci*omega*t0)) / omega
       
          do n = 1,nmax
             phi_out(n) = phi_out(n) - delta * psi_out(n)
          enddo


       !-----------------------------------------
       ! 3 : momentum forcing  F(t) = -i exp(i ω t) ∂x ψ(t)
       !-----------------------------------------
       case(3)

       do n=1,nmax
          En = eigval(n)

          if(n > 1) then
       
             Em = eigval(n-1)
       
             phase = omega + En - Em
       
             delta = exp(-ci*En*t) * &
                     ( exp(ci*phase*t) - exp(ci*phase*t0) ) &
                     /(ci*phase)
       
             phi_out(n) = phi_out(n) &
                + ( sqrt(dble(n)) * psi0(n-1) * delta ) / sqrt(2.d0)
       
          endif

         !-------------------
         ! coupling n+1
         !-------------------
      
         if(n < nmax) then
      
            Em = eigval(n+1)
      
            phase = omega + En - Em
      
            delta = exp(-ci*En*t) * &
                    ( exp(ci*phase*t) - exp(ci*phase*t0) ) &
                    /(ci*phase)
      
            phi_out(n) = phi_out(n) &
               - ( sqrt(dble(n+1)) * psi0(n+1) * delta ) / sqrt(2.d0)
      
         endif
       enddo

           
       end select
       
       end subroutine exact_closed_duhamel


       subroutine exact_increment_duhamel(nmax, dt, tn,        &
                                  eigval, psi_in, phi_in,             &
                                  psi_out, phi_out, src_type, omega)
       
       implicit none
       
       integer, intent(in) :: nmax, src_type
       real(8), intent(in) :: omega, dt, tn
       real(8), intent(in) :: eigval(nmax)
       
       complex(8), intent(in)  :: psi_in(nmax), phi_in(nmax)
       complex(8), intent(out) :: psi_out(nmax), phi_out(nmax)
       
       integer :: n
       real(8) :: En, Em, phase
       complex(8) :: evol, delta
       real(8) :: tnp1
       
       tnp1 = tn + dt
       
       !-----------------------------------
       ! homogeneous propagation
       !-----------------------------------
       do n = 1,nmax
          En = eigval(n)
          evol = exp(-ci*En*dt)
       
          psi_out(n) = evol * psi_in(n)
          phi_out(n) = evol * phi_in(n)
       enddo
       
       !-----------------------------------
       ! duhamel term
       !-----------------------------------
       select case(src_type)
       
       !-----------------------------------
       ! 1 : constant forcing
       !-----------------------------------
       case(1)
       
       do n=1,nmax
       
          En = eigval(n)
       
          delta = exp(-ci*En*tnp1) * &
                  ( exp(ci*(En+omega)*tnp1) - exp(ci*(En+omega)*tn) ) &
                  /(ci*(En+omega))
       
          phi_out(n) = phi_out(n) - delta * psi_in(n)
       
       enddo
       
       
       !-----------------------------------
       ! 2 : F(t)=psi(t)
       !-----------------------------------
       case(2)
       
       delta = (exp(ci*omega*tnp1) - exp(ci*omega*tn))/omega
       
       do n=1,nmax
          phi_out(n) = phi_out(n) - delta * psi_out(n)
       enddo
       
       
       !-----------------------------------
       ! 3 : momentum forcing
       !-----------------------------------
       case(3)
       
       do n=1,nmax
       
          En = eigval(n)
       
          if(n > 1) then
       
             Em = eigval(n-1)
             phase = omega + En - Em
       
             delta = exp(-ci*En*tnp1) * &
                     (exp(ci*phase*tnp1) - exp(ci*phase*tn)) /(ci*phase)
       
             phi_out(n) = phi_out(n) &
                + (sqrt(dble(n))*psi_in(n-1)*delta)/sqrt(2.d0)
       
          endif
       
       
          if(n < nmax) then
       
             Em = eigval(n+1)
             phase = omega + En - Em
       
             delta = exp(-ci*En*tnp1) * &
                     (exp(ci*phase*tnp1) - exp(ci*phase*tn))/(ci*phase)
       
             phi_out(n) = phi_out(n) &
                - (sqrt(dble(n+1))*psi_in(n+1)*delta)/sqrt(2.d0)
       
          endif
       
       enddo
       
       end select
       
       end subroutine exact_increment_duhamel


      end module conv_tests
