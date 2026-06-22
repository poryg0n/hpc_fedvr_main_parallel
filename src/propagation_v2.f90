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



      subroutine step_strang(nmax, nchannel, dt, t, xx,               &
                                           eigval, eigvec,            &
                                           vec_in, vec_out )
          implicit none
          integer, intent(in) :: nmax, nchannel
          real(dp), intent(in) :: dt, t
          real(dp), intent(in) :: eigval(nmax)
          real(dp), intent(in) :: xx(nmax)
          real(dp), intent(in), contiguous :: eigvec(:,:)
          complex(dp), intent(in) :: vec_in(nmax, nchannel+1)
          complex(dp), intent(out) :: vec_out(nmax, nchannel+1)

          real(dp) :: dt2, tmid
          real(dp) :: vext(nmax)
          complex(dp) :: tmp(nmax, nchannel+1)
          complex(dp) :: work(nmax, nchannel+1)
      
          complex(dp), parameter :: c0 = (0.0_dp, 0.0_dp)
          complex(dp), parameter :: c1 = (1.0_dp, 0.0_dp)
          complex(dp), parameter :: ci = (0.0_dp, 1.0_dp)

          complex(dp) :: phase
          integer :: i, k

          dt2 = 0.5_dp*dt
          tmid = t+dt2
          vext = force_t(tmid)*xx
!         vext = 0.d0
      
          ! ----- Half-step in energy basis -----
!         tmp = vec_in * exp(-ci * eigval * dt2)
          do i=1,nmax
             phase = exp(-ci*eigval(i)*dt2)
             tmp(i,:) = phase * vec_in(i,:)  
          enddo
      
          ! ----- Rotate to coordinate basis -----
          do k=1,nchannel+1
             work(:,k) = matmul(eigvec, tmp(:,k))
          enddo
!         call zgemm('n','n', nmax, nchannel+2, nmax,                  &
!                   c1, eigvec, nmax, tmp, nmax, c0, work, nmax)

      
          ! ----- Full-step in coordinate basis -----
!         tmp = work * exp(-ci * vext * dt)
          do i=1,nmax
             phase = exp(-ci*vext(i)* dt)
             tmp(i,:) = phase * work(i,:)  
          enddo
      
          ! ----- Rotate back to energy basis -----
          do k=1,nchannel+1
             work(:,k) = matmul(transpose(eigvec), tmp(:,k))
          enddo
!         call zgemm('t','n', nmax, nchannel+2, nmax,                 &
!                   c1, eigvec, nmax, tmp, nmax, c0, work, nmax)
      
          ! ----- Half-step in energy basis -----
!         vec_out = work * exp(-ci * eigval * dt2)
          do i=1,nmax
             phase = exp(-ci*eigval(i)*dt2)
             vec_out(i,:) = phase * work(i,:)  
          enddo

      end subroutine step_strang


      subroutine step_yoshida4(                                    &
          nmax, nchannel,                                          &
          dt, t, xx,                                               &
          eigval, eigvec,                                          &
          vec_in, vec_out )

        implicit none
        integer, intent(in) :: nmax, nchannel
        real(dp), intent(in) :: dt, t
        real(dp), intent(in) :: xx(nmax)
        real(dp), intent(in) :: eigval(nmax)
        real(dp), intent(in) :: eigvec(nmax,nmax)
        complex(dp), intent(in)  :: vec_in(nmax, nchannel+1)
        complex(dp), intent(out) :: vec_out(nmax, nchannel+1)

        real(dp), parameter :: a1 = 1.0_dp / (2.0_dp - 2.0_dp**(1.0_dp/3.0_dp))
        real(dp), parameter :: a2 = -2.0_dp**(1.0_dp/3.0_dp) * a1

        complex(dp) :: v1(nmax, nchannel+1), v2(nmax, nchannel+1)
        real(dp) :: t1, t2, t3
        real(dp) :: vext(nmax)

        ! --- stage times (midpoint freezing per substep) ---

        t1 = t
        t2 = t + a1*dt
        t3 = t + (a1 + a2)*dt
        ! --- γ step ---
        call step_strang(nmax, nchannel, a1*dt, t1,                    &
                                      xx, eigval, eigvec, vec_in, v1)

        ! --- δ step ---
        call step_strang(nmax, nchannel, a2*dt, t2,                    &
                                      xx, eigval, eigvec, v1, v2)

        ! --- γ step ---
        call step_strang(nmax, nchannel, a1*dt, t3,                    &
                                      xx, eigval, eigvec, v2, vec_out)

      end subroutine step_yoshida4



      subroutine split_operator(                                      &
              nmax, nchan,                                         &
              dt, t, xx,                                              &
              eigval, eigvec,                                        &
              vec_in, vec_out, order )

        implicit none

        integer, intent(in) :: nmax, nchan, order
        real(dp), intent(in) :: dt, t
        real(dp), intent(in) :: xx(nmax)
        real(dp), intent(in) :: eigval(nmax)
        real(dp), intent(in) :: eigvec(nmax,nmax)

        complex(dp), intent(in)  :: vec_in(nmax, nchan+1)
        complex(dp), intent(out) :: vec_out(nmax, nchan+1)

        select case(order)

        case(2)

           call step_strang(                                           &
                nmax, nchan, dt, t, xx,                            &
                eigval, eigvec,                                       &
                vec_in, vec_out )

        case(4)

           call step_yoshida4(                                         &
                nmax, nchan, dt, t, xx,                            &
                eigval, eigvec,                                       &
                vec_in, vec_out )

        case default

           print *, 'split_operator: unsupported order = ', order
           stop

        end select

      end subroutine split_operator


      subroutine build_source_quadrature( nchan, nmax, lnbr, nnbr,    &
                                                  xs, xx, map, Dref,  &
                                                  dt, t, omega,       &
                                                  eigval, eigvec,     &
                                                  svec,               &
                                                  src,                &
                                                  order )
         implicit none
         integer, intent(in) :: nmax, lnbr, nnbr
         integer, intent(in) :: nchan, order
         real(8), intent(in) :: dt, t
         real(8), intent(in) :: xs(:), xx(:)
         real(8), intent(in) :: omega(nchan+1)
         integer, intent(in) :: map(:,:)
         real(8), intent(in) :: Dref(:,:)
         real(8), intent(in) :: eigval(:), eigvec(:,:)
         complex(8), intent(in)  :: svec(nmax,3)
         complex(8), intent(out) :: src(nmax,nchan)
      

         if (order.eq.2) then
            call midpoint_quadrature( nchan, nmax, lnbr, nnbr,        &
                                                  xs, xx, map, Dref,  &
                                                  dt, t, omega,       &
                                                  eigval, eigvec,     &
                                                  svec,               &
                                                  src, order )
         else 

            call simpson_quadrature( nchan, nmax, lnbr, nnbr,         &
                                                  xs, xx, map, Dref,  &
                                                  dt, t, omega,       &
                                                  eigval, eigvec,     &
                                                  svec,               &
                                                  src, order )
         end if


!        src = src * exp(-ci * eigval * (t+dt) )

      end subroutine build_source_quadrature







      subroutine midpoint_quadrature( nchan, nmax, lnbr, nnbr,        &
                                                  xs, xx, map, Dref,  &
                                                  dt, t, omega,       &
                                                  eigval, eigvec,     &
                                                  svec,            &
                                                  src, order )
      
         implicit none
         integer, intent(in) :: nmax, lnbr, nnbr
         integer, intent(in) :: nchan, order
         real(8), intent(in) :: dt, t
         real(8), intent(in) :: xs(:), xx(:)
         real(8), intent(in) :: omega(nchan+1)
         integer, intent(in) :: map(:,:)
         real(8), intent(in) :: Dref(:,:), eigval(:), eigvec(:,:)
         complex(8), intent(in)  :: svec(nmax,3)
         complex(8), intent(out) :: src(nmax, nchan)
      
         complex(8) :: src_mid(nmax, nchan)
         complex(8) :: auxm(nmax)
         real(8) :: dt2, tmid

         integer :: j
      
         dt2  = 0.5d0 * dt
         tmid = t + dt2

         do j=1,nchan
!           src_mid = exp(ci * ( omega + eigval ) * tmid) * svec(:,2)
            src_mid(:,j) = exp(ci * ( omega(j+1) ) * tmid) * svec(:,2)
         enddo
         !--------------------------------------------
         ! 3) midpoint quadrature
         !--------------------------------------------
         src = dt * src_mid

      end subroutine midpoint_quadrature



      subroutine simpson_quadrature( nchan, nmax, lnbr, nnbr,   &
                                                  xs, xx, map, Dref,  &
                                                  dt, t, omega,       &
                                                  eigval, eigvec,     &
                                                  svec,               &
                                                  src, order )
      
         implicit none
         integer, intent(in) :: nmax, lnbr, nnbr
         integer, intent(in) :: nchan, order
         real(8), intent(in) :: dt, t
         real(8), intent(in) :: xs(:), xx(:)
         real(8), intent(in) :: omega(nchan+1)
         integer, intent(in) :: map(:,:)
         real(8), intent(in) :: Dref(:,:), eigval(:), eigvec(:,:)
         complex(8), intent(in)  :: svec(nmax,3)
         complex(8), intent(out) :: src(nmax,nchan)
      
         real(8) :: vext(nmax)
         complex(8) :: aux(nmax,3)
         complex(8) :: srck(nmax,nchan,3)
         real(8) :: tau, tt

         integer :: j, k
      

         do k=1,3
            do j=1,nchan
               tau = 0.5d0*(k-1)*dt
               tt= t + tau

!              srck(:,k) = exp(ci * ( omega + eigval) * tt)  * svec(:,k)
!              srck(:,k) = exp(ci * ( omega ) * tt)  * svec(:,k)
               srck(:,j,k) = exp(ci * ( omega(j+1) ) * tt)  * svec(:,k)

            enddo
         enddo

         src(:,:) = (dt/6.d0) *                                        &
                    ( srck(:,:,1) + 4.d0 * srck(:,:,2) + srck(:,:,3) )
      
      end subroutine simpson_quadrature



      subroutine process_src_ingredients( nchan,                     &
                                         nmax, lnbr, nnbr,           &
                                         jac,                        &
                                         xs, xx, wx,                 &
                                         map, Dref,                  &
                                         dt, t,                      &
                                         eigval, eigvec,         &
                                         omega,               &
                                         psi_in, svec,           &
                                         src_type, order )

      implicit none
      integer, intent(in) :: nmax, lnbr, nnbr
      integer, intent(in) :: nchan, order, src_type
      real(8), intent(in) :: dt, t, jac
      real(8), intent(in) :: omega(:)
      real(8), intent(in) :: xs(:), xx(:), wx(:)
      integer, intent(in) :: map(:,:)
      real(8), intent(in) :: Dref(:,:)
      real(8), intent(in) :: eigval(:), eigvec(:,:)
      complex(8), intent(in) :: psi_in(nmax)
      complex(8), intent(out), target :: svec(nmax,3)


      integer :: k
      real(8) :: dt2, tau, delta
      complex(8) :: svec0(nmax), aux(nmax)
      complex(8), target :: aux0(nmax,3)
      complex(8) :: psi_inx(nmax)
      complex(8) :: dpsi_x(nmax)
!     complex(8), pointer :: psi_mat(nmax,1)
      complex(8), pointer :: psi_mat(:,:)

      complex(8), pointer :: vec_in(:,:), vec_out(:,:)
!     complex(8), target :: svec0_2d(nmax,1)

      dt2 = 0.5d0*dt
     
      select case(src_type)
         case(1)
            aux = 0.d0
            aux(1) = (1.d0, 0.d0)

         case default
            aux = psi_in
      end select


     ! --- apply whatever function of time to the argument if there is ---
     call apply_stuff_to_arg(nmax, nchan, xx, dt, t,                  &
          jac, wx, eigval, eigvec, omega, aux, svec0, src_type, order)

      !-----------------------------------------
      ! build Simpson nodes F(t), F(t+dt/2), F(t+dt)
      !-----------------------------------------

      call build_source_vector(nmax, nchan, xx, dt, t,                 &
                               jac, wx, eigval, eigvec,                &
                               omega, svec0, svec, src_type, order)
        
         do k = 1,3
            tau   = t + 0.5d0*(k-1)*dt
            delta = t + dt - tau
            aux0(:,k) = svec(:,k)
            ! alias directly into svec
!           psi_mat(1:nmax,1:1) => svec(:,k:k)
            vec_in => aux0(:,k:k)
            vec_out => svec(:,k:k)
      

            !  Transport
!           call split_operator(nmax, 0, delta, tau, xx,         &
!                       eigval, eigvec, psi_mat(:,1), svec(:,k), order)
!
         
!           call split_operator(nmax, 0, delta, tau, xx, &
!                               eigval, eigvec, psi_mat, psi_mat, order)
            call split_operator(nmax, 0, delta, tau, xx,              &
                            eigval, eigvec, vec_in, vec_out, order)
         
         enddo
      

      end subroutine process_src_ingredients


      subroutine build_source_vector(nmax, nchan, xx, dt, t,          &
                       jac, wx, eigval, eigvec, omega, svec0, svec,   &
                       src_type, order)
        implicit none
        integer, intent(in) :: nmax
        integer, intent(in) :: nchan, order, src_type
        real(8), intent(in) :: t, dt, jac
        real(8), intent(in) :: omega(nchan+1)
        real(8), intent(in) :: xx(nmax), wx(nmax), eigval(nmax)
        real(8), intent(in) :: eigvec(nmax,nmax)

        complex(8), intent(in), target :: svec0(nmax)
        complex(8), intent(out), target :: svec(nmax,3)

        integer :: k
        real(8) :: tau, tt, dt2

        complex(8) :: dpsi_x(nmax), psi_inx(nmax)
        complex(8) :: tmp(nmax)

        complex(8), pointer :: vec_in(:,:), vec_out(:,:)
        complex(8), target :: svec0_2d(nmax,1)


      
        do k = 1,3
           tau = 0.5d0 * (k-1)*dt

!          svec(:,k) = svec0
           svec0_2d(:,1) = svec0

           ! alias input
           vec_in  => svec0_2d

           ! alias output
           vec_out => svec(:,k:k)


           ! default copy: svec0 --> svec(:,k) 
!          vec_out(1:nmax,1:1) => vec_in(1:nmax,1:1)

           if (src_type.gt.1) then  
             !  Transport  (if the function is not constant, i.e
             !  src_type/=1 ) 
!            call split_operator(nmax, nchan, tau, t, xx,             &
!                           eigval, eigvec, svec0, svec(:,k), order)
             call split_operator(nmax, 0, tau, t, xx,             &
                            eigval, eigvec, vec_in, vec_out, order)


              if (src_type.eq.3) then  
                 ! -i \partial_x \psi(x,t)
                 tmp = vec_out(:,1)
                 call apply_momentum_operator(nmax, eigvec, xx, wx,   &
                                           jac, tmp, vec_out(:,1), 0)
          
              end if
           end if
        enddo

      
      end subroutine build_source_vector


      subroutine apply_stuff_to_arg(nmax, nchan, xx, dt, t, jac,       &
                wx, eigval, eigvec, omega, aux, svec0, src_type, order)
        implicit none
        integer, intent(in) :: nmax, order, nchan, src_type
        real(8), intent(in) :: t, dt, jac
        real(8), intent(in) :: omega(nchan+1)
        real(8), intent(in) :: xx(nmax), wx(nmax), eigval(nmax)
        real(8), intent(in) :: eigvec(nmax,nmax)
        complex(8), intent(in) :: aux(nmax)
        complex(8), intent(out) :: svec0(nmax)

        integer :: k
        real(8) :: tau, tt, dt2

        svec0 = 1.d0 * aux

        ! *** might replace that by a function g(t) later
!       ! svec = g(t)\psi(t)
!                 call apply_momentum_operator(nmax, eigvec, xx, wx,  &
!                            jac, aux, svec0, 0)
          
      
      end subroutine apply_stuff_to_arg




      end module propagation
