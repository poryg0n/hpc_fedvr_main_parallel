      module propagation
        use fedvr_derivative_ops, only : eval_dpsi_fedvr_reduced
!       use force_field, only : force_t
        use dynamic_parameters
        use math_util
        implicit none
        private
        public :: split_operator
        public :: process_src_ingredients
        public :: build_source_quadrature


        integer, parameter :: dp = kind(1.0d0)
        complex(dp), parameter :: ci = (0.0_dp, 1.0_dp)
      
      contains



      subroutine step_strang(nmax, dt, t, xx, eigval, eigvec,    &
                                           vec_in, vec_out )
          implicit none
          integer, intent(in) :: nmax
          real(dp), intent(in) :: dt, t
          real(dp), intent(in) :: eigval(nmax)
          real(dp), intent(in) :: xx(nmax)
          real(dp), intent(in), contiguous :: eigvec(:,:)
          complex(dp), intent(in) :: vec_in(nmax)
          complex(dp), intent(out) :: vec_out(nmax)

          real(dp) :: dt2, tmid
          real(dp) :: vext(nmax)
          complex(dp) :: tmp(nmax)
          complex(dp) :: work(nmax)
      
          complex(dp), parameter :: ci = (0.0_dp, 1.0_dp)

          dt2 = 0.5_dp*dt
          tmid = t+dt2
          vext = force_t(tmid)*xx
!         vext = 0.d0
      
          ! ----- Half-step in energy basis -----
          tmp = vec_in * exp(-ci * eigval * dt2)
      
          ! ----- Rotate to coordinate basis -----
          work = matmul(eigvec, tmp)
      
          ! ----- Full-step in coordinate basis -----
          tmp = work * exp(-ci * vext * dt)
      
          ! ----- Rotate back to energy basis -----
          work = matmul(transpose(eigvec), tmp)
      
          ! ----- Half-step in energy basis -----
          vec_out = work * exp(-ci * eigval * dt2)

      end subroutine step_strang


      subroutine step_yoshida4(                                    &
          nmax, dt, t, xx,                                         &
          eigval, eigvec,                                          &
          vec_in, vec_out )

        implicit none
        integer, intent(in) :: nmax
        real(dp), intent(in) :: dt, t
        real(dp), intent(in) :: xx(nmax)
        real(dp), intent(in) :: eigval(nmax)
        real(dp), intent(in) :: eigvec(nmax,nmax)
        complex(dp), intent(in)  :: vec_in(nmax)
        complex(dp), intent(out) :: vec_out(nmax)

        real(dp), parameter :: a1 = 1.0_dp / (2.0_dp - 2.0_dp**(1.0_dp/3.0_dp))
        real(dp), parameter :: a2 = -2.0_dp**(1.0_dp/3.0_dp) * a1

        complex(dp) :: v1(nmax), v2(nmax)
        real(dp) :: t1, t2, t3
        real(dp) :: vext(nmax)

        ! --- stage times (midpoint freezing per substep) ---

        t1 = t
        t2 = t + a1*dt
        t3 = t + (a1 + a2)*dt
        ! --- γ step ---
        call step_strang(nmax, a1*dt, t1, xx, eigval, eigvec, vec_in, v1)

        ! --- δ step ---
        call step_strang(nmax, a2*dt, t2, xx, eigval, eigvec, v1, v2)

        ! --- γ step ---
        call step_strang(nmax, a1*dt, t3, xx, eigval, eigvec, v2, vec_out)
      end subroutine step_yoshida4



      subroutine split_operator(                                      &
              nmax, dt, t, xx,                                       &
              eigval, eigvec,                                        &
              vec_in, vec_out, order )

        implicit none

        integer, intent(in) :: nmax, order
        real(dp), intent(in) :: dt, t
        real(dp), intent(in) :: xx(nmax)
        real(dp), intent(in) :: eigval(nmax)
        real(dp), intent(in) :: eigvec(nmax,nmax)

        complex(dp), intent(in)  :: vec_in(nmax)
        complex(dp), intent(out) :: vec_out(nmax)

        select case(order)

        case(2)

           call step_strang(                                           &
                nmax, dt, t, xx,                                      &
                eigval, eigvec,                                       &
                vec_in, vec_out )

        case(4)

           call step_yoshida4(                                         &
                nmax, dt, t, xx,                                      &
                eigval, eigvec,                                       &
                vec_in, vec_out )

        case default

           print *, 'split_operator: unsupported order = ', order
           stop

        end select

      end subroutine split_operator


      subroutine build_source_quadrature( nmax, lnbr, nnbr,           &
                                                  xs, xx, map, Dref,  &
                                                  dt, t,              &
                                                  eigval, eigvec,     &
                                                  svec,               &
                                                  src, omega,         &
                                                  order )
         implicit none
         integer, intent(in) :: nmax, lnbr, nnbr, order
         real(8), intent(in) :: dt, t, omega
         real(8), intent(in) :: xs(:), xx(:)
         integer, intent(in) :: map(:,:)
         real(8), intent(in) :: Dref(:,:)
         real(8), intent(in) :: eigval(:), eigvec(:,:)
         complex(8), intent(in)  :: svec(nmax,3)
         complex(8), intent(out) :: src(:)
      

         if (order.eq.2) then
            call midpoint_quadrature( nmax, lnbr, nnbr,               &
                                                  xs, xx, map, Dref,  &
                                                  dt, t,              &
                                                  eigval, eigvec,     &
                                                  svec,               &
                                                  src, omega, order )
         else 

            call simpson_quadrature( nmax, lnbr, nnbr,                &
                                                  xs, xx, map, Dref,  &
                                                  dt, t,              &
                                                  eigval, eigvec,     &
                                                  svec,               &
                                                  src, omega, order )
         end if


!        src = src * exp(-ci * eigval * (t+dt) )

      end subroutine build_source_quadrature




      subroutine midpoint_quadrature( nmax, lnbr, nnbr,               &
                                                  xs, xx, map, Dref,  &
                                                  dt, t,              &
                                                  eigval, eigvec,     &
                                                  svec,            &
                                                  src, omega, order )
      
         implicit none
         integer, intent(in) :: nmax, lnbr, nnbr, order
         real(8), intent(in) :: dt, t, omega
         real(8), intent(in) :: xs(:), xx(:)
         integer, intent(in) :: map(:,:)
         real(8), intent(in) :: Dref(:,:), eigval(:), eigvec(:,:)
         complex(8), intent(in)  :: svec(nmax,3)
         complex(8), intent(out) :: src(:)
      
         complex(8) :: src_mid(nmax)
         complex(8) :: auxm(nmax)
         real(8) :: dt2, tmid
      
         dt2  = 0.5d0 * dt
         tmid = t + dt2

         src_mid = exp(ci * ( omega ) * tmid) * svec(:,2)
         !--------------------------------------------
         ! 3) midpoint quadrature
         !--------------------------------------------
         src = dt * src_mid

      end subroutine midpoint_quadrature



      subroutine simpson_quadrature( nmax, lnbr, nnbr,   &
                                                  xs, xx, map, Dref,  &
                                                  dt, t,              &
                                                  eigval, eigvec,     &
                                                  svec,               &
                                                  src, omega, order )
      
         implicit none
         integer, intent(in) :: nmax, lnbr, nnbr, order
         real(8), intent(in) :: dt, t, omega
         real(8), intent(in) :: xs(:), xx(:)
         integer, intent(in) :: map(:,:)
         real(8), intent(in) :: Dref(:,:), eigval(:), eigvec(:,:)
         complex(8), intent(in)  :: svec(nmax,3)
         complex(8), intent(out)  :: src(:)
      
         real(8) :: vext(nmax)
         complex(8) :: aux(nmax,3)
         complex(8) :: srck(nmax,3)
         real(8) :: tau, tt

         integer :: k
      

         do k=1,3
            tau = 0.5d0*(k-1)*dt
            tt= t + tau

            srck(:,k) = exp(ci * ( omega ) * tt)  * svec(:,k)

         enddo

         src = (dt/6.d0) * ( srck(:,1) + 4.d0 * srck(:,2) + srck(:,3) )
      
      end subroutine simpson_quadrature



      subroutine process_src_ingredients( nmax, lnbr, nnbr,           &
                                         jac,          &
                                         xs, xx, wx,           &
                                         map, Dref,                    &
                                         dt, t,                  &
                                        eigval, eigvec,         &
                                        psi_in, svec,           &
                                        src_type, order )

      implicit none
      integer, intent(in) :: nmax, lnbr, nnbr, order, src_type
      real(8), intent(in) :: dt, t, jac
      real(8), intent(in) :: xs(:), xx(:), wx(:)
      integer, intent(in) :: map(:,:)
      real(8), intent(in) :: Dref(:,:)
      real(8), intent(in) :: eigval(:), eigvec(:,:)
      complex(8), intent(in)  :: psi_in(nmax)
      complex(8), intent(out) :: svec(nmax,3)


      integer :: k
      real(8) :: dt2, tau, delta
      complex(8) :: svec0(nmax)
      complex(8) :: aux0(nmax)
      complex(8) :: auxc(nmax)
      complex(8) :: psi_inx(nmax)
      complex(8) :: dpsi_x(nmax)

      dt2 = 0.5d0*dt
     
      select case(src_type)
         case(1)
            auxc = 0.d0
            auxc(1) = (1.d0, 0.d0)

         case default
            auxc = psi_in
      end select


     ! --- apply whatever function of time to the argument if there is ---
     call apply_stuff_to_arg(nmax, xx, dt, t,      &
              jac, wx, eigval, eigvec, auxc, svec0, src_type, order)

      !-----------------------------------------
      ! build Simpson nodes F(t), F(t+dt/2), F(t+dt)
      !-----------------------------------------

      call build_source_vector(nmax, xx, dt, t,                 &
                               jac, wx, eigval, eigvec, &
                               svec0, svec, src_type, order)
        
         do k = 1,3
            tau = t+0.5d0 * (k-1)*dt
            delta = t+dt-tau
            aux0 = svec(:,k)

            ! 1 -  at t = t,         therefore -> 1.0dt step 
            ! 2 -  at t = t+0.5dt,   therefore -> 0.5dt step 
            ! 3 -  at t = t+dt,      therefore -> 0.0dt step 
            !  Transport
            call split_operator(nmax, delta, tau, xx,                &
                           eigval, eigvec, aux0, auxc, order)

            svec(:,k) = auxc
         enddo
      

      end subroutine process_src_ingredients


      subroutine build_source_vector(nmax, xx, dt, t,       &
                       jac, wx, eigval, eigvec, svec0, svec,   &
                       src_type, order)
        implicit none
        integer, intent(in) :: nmax, src_type, order
        real(8), intent(in) :: t, dt, jac
        real(8), intent(in) :: xx(nmax), wx(nmax), eigval(nmax)
        real(8), intent(in) :: eigvec(nmax,nmax)
        complex(8), intent(in) :: svec0(nmax)
        complex(8), intent(out) :: svec(nmax,3)

        integer :: k
        real(8) :: tau, tt, dt2

        complex(8) :: dpsi_x(nmax), psi_inx(nmax)
        complex(8) :: aux_in(nmax)
        complex(8) :: aux_out(nmax)
      
        do k = 1,3
           tau = 0.5d0 * (k-1)*dt

            svec(:,k) = svec0

            if (src_type.gt.1) then  
              !  Transport  (if the function is not constant, i.e  src_type/=1 ) 
              ! 1 -  transport from t -> t,       therefore tau -> 0.0dt step 
              ! 2 -  transport from t -> 0.5dt,   therefore tau -> 0.5dt step 
              ! 3 -  transport from t -> t+0.5dt, therefore tau -> 1.0dt step 
              call split_operator(nmax, tau, t, xx,                   &
                             eigval, eigvec, svec0, aux_out, order)


               if (src_type.eq.3) then  
                  ! -i \partial_x \psi(x,t)
                  aux_in = aux_out
                  call apply_momentum_operator(nmax, eigvec, xx, wx,  &
                             jac, aux_in, aux_out, 0)
               end if

               svec(:,k) = aux_out
            end if
        enddo

      
      end subroutine build_source_vector


      subroutine apply_stuff_to_arg(nmax, xx, dt, t, jac, wx,          &
                    eigval, eigvec, svec0, svec, src_type, order)
        implicit none
        integer, intent(in) :: nmax, order, src_type
!       real(8), intent(in) :: t, dt, omega, jac
        real(8), intent(in) :: t, dt, jac
        real(8), intent(in) :: xx(nmax), wx(nmax), eigval(nmax)
        real(8), intent(in) :: eigvec(nmax,nmax)
        complex(8), intent(in) :: svec0(nmax)
        complex(8), intent(out) :: svec(nmax)

        integer :: k
        real(8) :: tau, tt, dt2

        svec = 1.d0 * svec0

        ! *** might replace that by a function g(t) later
!       ! svec = g(t) * psi(t)
!                 call apply_momentum_operator(nmax, eigvec, xx, wx,  &
!                            jac, svec0, svec, 0)
          
      
      end subroutine apply_stuff_to_arg




      end module propagation
