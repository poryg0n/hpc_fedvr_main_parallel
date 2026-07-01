      program fedvr_build
      use omp_lib
      use constants, only : omg_ => omega_qho
      use lobatto
      use fedvr_topology
      use fedvr_derivative_ops
      use global_assembly
      use fedvr_conf_struct
      use time_grid
      use force_field
      use propagation
      use conv_tests
      implicit none
      integer :: lwork, info, store_val,                              &
                  i, j, k, ij, p, q, r, s, m, n,                      &
                  jn, ijk,                                            &
                  nmax, nmax_, krange, ios, nt

      real(8) :: t_ini, t_end, emax, tau,                              &
                    step, eft, t, duration,                           &
                    start, finish, lap, aux, aux1, aux2,              &
                    abstol, pemd, prob, p_0, p_ion, dk,               &
                    t0, dt0, tt, tp, tc, tmid,                        &
                    f0, fmid,                                         &
                    err1, err2, order,                                &
                    jac, inv_jac, xx1, xx2, x_min, rowsum

      complex(8) :: cnum, a0, cprobk, auxcc, auxcc1, auxcc2


      real(8), allocatable, dimension(:) :: xa, wa, xs, xx, wx, xg,   &
                                   vec_matup, eigval,                 &
                                   kk, xt,                            & 
                                   auxr1, auxr2, auxr3,               &
                                   time_, omg,                        &
                                   fvec,                              &
                                   p_k, dx, dt,                       &
                                   qvc1, qvc2, qvc3, qvc4,            &
                                   matint_v, vext_mid, vext_,         &
                                   ipiv, work, ifail, dnrm2
  
      complex(8), allocatable, dimension(:) :: wf0, wfc0,            &
                                               dwf0, dwfc0,           &
                                               wf, wfc, wfc_,         &
                                               psi_exact,             &
                                               psi0, psi, psi_x,      &
                                               psi_in, psi_out,       &
                                               phi_in, phi_out,       &
                                               psi_inx, psi_outx,     &
                                               phi_inx, phi_outx,     &
                                               src_out, src_mid,      &
                                               src, src_x,            &
                                               wwfc, dwwfc,           &
                                               dwfc, dwfc2, dwfc3,    & 
                                               fftpsi,                &
                                               init1, init2,          &
                                               phi0, phic0,           &
                                               psi_test, phi_test,    &
                                               phi, phic,             &
                                               auxc, auxck1, auxck2,  &
                                               rr_, ak, cprob2,       &
                                               auxc1, auxc2,          &
                                               auxc3, auxc4,          &
                                               pk0, b0w_1, b0w_2,     &
                                               d_t, dd_t, d_w, dd_w,  &
                                               x_t,                   &
                                               p_t, pp_t, p_w, pp_w,  &
                                               dp_t, dx_t, ddx_t,     &
                                               dp_t_, vextvc, dvextv, &
                                               svec,                  &
                                               ss, t1, t2, ft, intft

      integer, allocatable, dimension(:,:) :: map

      real(8), allocatable, dimension(:,:) :: lu, id, inv,            &
                                       pkk, bkw_1, bkw_2,             &
                                       Dref, Dglobal,                 &
                                       eigvec, basis

      real(8), allocatable, dimension(:,:,:) :: Dloc_all          

      complex(8), allocatable, dimension(:,:) :: in_states,           &
                                                 out_states,          &
                                                 workspace,           &
                                                 cwf, xmat

      character(255) :: workdir


      include 'param'


      call cpu_time(start)

      krange = 2.d0*nmax

      allocate(xa(nnbr),wa(nnbr),xs(1:lnbr+1))
      call lobatto_compute(nnbr,xa,wa) 

!     do i=1,nnbr
!        write(*,*) i, xa(i), wa(i)
!     enddo


      allocate(Dref(nnbr,nnbr)) 
      call build_Dref_lobatto(nnbr, xa, Dref)

!     x_min = -0.5d0 * xrange 
      call build_elements_bounds(lnbr, xrange, xc, xs)

      allocate(Dloc_all(lnbr,nnbr,nnbr))
      do s = 1, lnbr
        call build_Dloc_all(nnbr, xs, s, Dref, Dloc_all(s,:,:))
      end do
      write(*,*) 

      allocate(map(lnbr, nnbr))
      call build_fedvr_map( lnbr, nnbr, map, nmax_ )
      ! dirichlet bc applied map dimension
      if ( maxval(map) /= nmax_ ) stop "Topology mismatch"

      nmax = nmax_ - 2 

      allocate(Dglobal(nmax_,nmax_))
      call assemble_derivative_fedvr(lnbr, nnbr, map, wa,             &
!                                 Dloc_all, Dglobal, FEDVR_SYMMETRIC )
                                   Dloc_all, Dglobal, FEDVR_MUTUAL )
!                                  Dloc_all, Dglobal, FEDVR_LEFT )

      allocate(xg(nmax_))
      call build_x_global_fedvr(lnbr, nnbr, xa, xs, map, xg)

      do i = 2, nmax_
         if (xg(i) <= xg(i-1)) then
            write(*,*) "Grid ordering error at", i
            stop
         end if
      end do

      allocate(psi(nmax_), psi_x(nmax_))

      inv_jac = 0.5d0*xrange/lnbr
      allocate(eigval(nmax), ss(nmax))
      allocate(xx(nmax),wx(nmax))

!      do i=1,nmax
!         write(*,*) i, xx(i), eigval(i)
!      enddo

      ss = 1.d0
      psi = xg
      psi_x = matmul(Dglobal, psi)

      ! ** deriv with operator D
      write(*,*) "psi_x with D operator, blunt"
      do i=1,nmax_
         write(*,*) xg(i), psi(i), psi_x(i)
      enddo
      write(*,*) 

!      call eval_dpsi_fedvr( lnbr, nnbr, xg, xs, map, Dref,    &
!                                             psi, psi_x, FEDVR_LEFT)


      ! ** deriv with eval_psi_full
      write(*,*) "psi_x with eval_dpsi_fedvr  (full)"
      do i=1,nmax_
         write(*,*) xg(i), psi(i), psi_x(i)
      enddo
      write(*,*) 



! *** Construct H = T + V and diagonalize with external legacy procedures
!      allocate(init1(nmax), init2(nmax))
      allocate( auxc(nmax), svec(nmax) )
      allocate(eigvec(nmax,nmax),xmat(nmax,nmax))
      allocate( wf0(nmax), wfc0(nmax), wf(nmax),wfc(nmax), wfc_(nmax) )
      call fedvr_hamilton_conf(inv_jac,xa,wa,xx,wx,eigval,             &
                         eigvec)
       wfc0 = eigvec(:,1)/wx/dsqrt(inv_jac)
       wf0 = wfc0 * wx * dsqrt(inv_jac)
       auxc=wf0
       wf0 = matmul(transpose(eigvec),auxc)

       write(*,*) 
       do i=1,nmax
          write(*,*) i, real(wfc0(i)), real(wf0(i))
       enddo
       write(*,*) 
       write(*,*) 


       allocate( workspace(nmax+2,2) )

!      deallocate(psi, psi_x)
!      allocate(psi(nmax), psi_x(nmax))
!      psi=xx 
!      write(*,*) "psi_x with eval_dpsi_fedvr_reduced)"
!      call eval_dpsi_fedvr_reduced( lnbr, nnbr, xs, map, Dref,        &
!                                  psi, psi_x, workspace )


!      write(*,*)
!      do i=1,nmax
!         write(*,*) i, xx(i), real(psi(i)), real(psi_x(i))
!      enddo
!      write(*,*)
!      pause


! --- time parameters ---
       t_ini =  0.d0
       t_end = 1.d0
       dt0   =  .1d0

       f0    = 1.d0

       xmat = (0.d0, 0.d0)

       call set_force_params(f0, omg_)
       call init_time_grid(t_ini, t_end, dt0)

 
       duration = int((t_end - t_ini))
       nt = int((t_end - t_ini)/dt0) 
       allocate(fvec(nt))
       t0 = t_ini
       tc = 0.d0

 
       allocate(psi0(nmax), phi0(nmax))
       allocate(psi_in(nmax), psi_out(nmax))
       allocate(psi_inx(nmax), psi_outx(nmax))
       allocate(phi_in(nmax), phi_out(nmax))
       allocate(phi_inx(nmax), phi_outx(nmax))
       allocate(psi_exact(nmax))

       allocate(src_mid(nmax), src(nmax))
       allocate(src_x(nmax))

! --- initial condition ---
       tt = t0                ! start time
       psi_in = wf0           ! initial wavefunction

       psi_inx = matmul(eigvec,psi_in)
       psi_inx = psi_inx/wx/dsqrt(inv_jac)


       write(*,*) 
       do i=1,nmax
          write(*,*) i, real(wfc0(i)), real(psi_inx(i))
       enddo
       write(*,*) 

!      call eval_dpsi_fedvr_reduced( lnbr, nnbr, xs, map, Dref,    &
!                                  psi_inx, phi_inx, workspace )

!      phi_inx = -ci * phi_inx
!      phi_in = matmul(transpose(eigvec),phi_inx)
       phi_in = psi_in
       phi0 = phi_in
       psi0 = psi_in

       write(*,'(a)') '# i  E_n'
       do i=1,10
          write(*,*) i, eigval(i)
       enddo

       svec = (0.0d0, 0.d0)
       svec(1) = (1.0d0, 0.d0)
       svec(2) = (1.0d0, 0.d0)

       write(*,'(a)') '# t  Re[psi(1)] Im[psi(1)] ||psi||         '//  &
                      '              Re[phi(1)] Im[phi(1)] ||phi||'

             write(*,'(f8.4,3x,6e16.8)') tt, &
                  real(psi_in(1)), aimag(psi_in(1)),                &
                                       sqrt(sum(abs(psi_in)**2)),    &
                  real(phi_in(1)), aimag(phi_in(1)),                &
                                       sqrt(sum(abs(phi_in)**2))

      write(*,*) 
      write(*,*) 
      write(*,*) "Convergence test"
      call test_against_exact_solution( nmax, lnbr, nnbr,             &
                                          xs, xx, map, Dref,          &
                                          dt0, t_end, t_ini, omega,   &
                                          eigval, eigvec,             &
                                          psi0, phi0,                 &
                                          2 )

      write(*,*)
      write(*,*)

      call test_richardson_inhomogeneous( nmax, lnbr, nnbr,   &
                                           xs, xx, map, Dref,         &
                                           dt0, t_end, t_ini, omega,  &
                                           eigval, eigvec,            &
                                           psi0, phi0,                &
                                           2 )




      pause
      call exact_closed(nmax, eigval, omega, svec,            &
                      t_end, t_ini, psi_in, phi_in, psi_out, phi_out)

      write(*,*) 
      write(*,*) "closed-form" 
      do i=1,nmax
         write(*,*) nt, i, psi_out(i), phi_out(i)
      enddo

      do i=1,nt
!        call exact_increment(nmax, eigval, omega, svec,     &
!                         dt0, tt, psi_in, phi_in, psi_out, phi_out)

         call step_midpoint_duhamel( nmax, lnbr, nnbr,            &
                                              xs, xx, map, Dref,     &
                                              dt0, tt, omega,        &
                                              eigval, eigvec,        &
                                              psi_in, phi_in,        &
                                              psi_out, phi_out )


         psi_in = psi_out
         phi_in = phi_out
         tt = tt + dt0
      enddo

      write(*,*) 
!     write(*,*) "Increment exact" 
      write(*,*) "Midpoint quadrature" 
!     write(*,*) "Simpson quadrature" 
      do i=1,nmax
         write(*,*) nt, i, psi_out(i), phi_out(i)
      enddo

      end program
