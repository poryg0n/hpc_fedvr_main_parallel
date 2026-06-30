      program fedvr_structure_build
      use omp_lib
      use constants, only : ci, omg_ => omega_qho
      use structure_parameters
      use math_util
      use util
      use lobatto
      use fedvr_topology
      use fedvr_derivative_ops
      use global_assembly
      use fedvr_conf_struct
      use propagation
!     use conv_tests
      implicit none
      integer :: lwork, info, store_val,                              &
                  i, j, k, ij, p, q, r, s, m, n,                      &
                  jn, ijk,                                            &
                  krange, ios, order, flag,                           &
                  nmax_
!                 nmax, nmax_, krange, ios, nt, order, flag

      real(8) ::    emax, tau,                              &
                    step, eft, duration,                           &
                    start, finish, lap, aux, aux1, aux2,              &
                    abstol, pemd, prob, p_0, p_ion, dk,               &
                    t0, dt, tt, tp, tc, tmid,                        &
                    fmid,                                             &
                    err1, err2,                                       &
                    xc_in, rowsum
!                   jac, inv_jac, xx1, xx2, x_min, rowsum

      complex(8) :: cnum, a0, cprobk, auxcc, auxcc1, auxcc2


      real(8), allocatable, dimension(:) :: xa, wa, xs, xx, wx, xg,   &
                                   vec_matup, eigval,                 &
                                   kk, xt,                            & 
                                   auxr1, auxr2, auxr3,               &
                                   time_, omg,                        &
                                   fvec,                              &
                                   p_k,                               &
!                                  dx, dt,                            &
!                                  p_k, dx, dt,                       &
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
                                               init1, init2,          &
                                               phi0, phic0,           &
                                               psi_test, phi_test,    &
                                               psi_ex, phi_ex,        &
                                               phi, phic,             &
                                               auxc, auxck1, auxck2,  &
                                               rr_, ak, cprob2,       &
                                               auxc1, auxc2,          &
                                               auxc3, auxc4,          &
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
                                                 svec,                &
                                                 cwf, xmat

      character(255) :: workdir


      include 'param'


      call cpu_time(start)
      call set_grid_params(np, ns, xc_in, xmin_, xmax_)

      call print_structure_parameters()
!     pause

!     krange = 2.d0*nmax

      allocate(xa(nnbr), wa(nnbr), xs(1:snbr+1))
      call lobatto_compute(nnbr, xa, wa) 

!     do i=1,nnbr
!        write(*,*) i, xa(i), wa(i)
!     enddo


      allocate(Dref(nnbr,nnbr)) 
      call build_Dref_lobatto(nnbr, xa, Dref)

!     x_min = -0.5d0 * xrange 
      call build_elements_bounds(snbr, xrange, xc, xs)

      allocate(Dloc_all(nnbr, nnbr, snbr))
      do s = 1, snbr
        call build_Dloc_all(nnbr, xs, s, Dref, Dloc_all(:,:,s))
      end do
      write(*,*) 

      allocate(map(snbr, nnbr))
      call build_fedvr_map( snbr, nnbr, map, nmax_ )
      ! dirichlet bc applied map dimension
      if ( maxval(map) /= nmax_ ) stop "Topology mismatch"

      write(*,*) "nmax", nmax
      write(*,*) "nmax_", nmax_

      if (nmax.ne.(nmax_-2)) stop

      allocate(Dglobal(nmax_,nmax_))
      call assemble_derivative_fedvr(snbr, nnbr, map, wa,             &
!                                 Dloc_all, Dglobal, FEDVR_SYMMETRIC )
                                   Dloc_all, Dglobal, FEDVR_MUTUAL )
!                                  Dloc_all, Dglobal, FEDVR_LEFT )

      allocate(xg(nmax_))
      call build_x_global_fedvr(snbr, nnbr, xa, xs, map, xg)

      do i = 2, nmax_
         if (xg(i) <= xg(i-1)) then
            write(*,*) "Grid ordering error at", i
            stop
         end if
      end do

      allocate(psi(nmax_), psi_x(nmax_))

      allocate(eigval(nmax), ss(nmax))
      allocate(xx(nmax),wx(nmax))


      ss = 1.d0
      psi = xg
      psi_x = matmul(Dglobal, psi)

       ! ** deriv with operator D
       write(*,*) "psi_x with D operator, blunt"
       do i=1,nmax_
          write(*,*) xg(i), psi(i), psi_x(i)
       enddo
       write(*,*) 
 
       call eval_dpsi_fedvr( snbr, nnbr, xg, xs, map, Dref,    &
                                              psi, psi_x, FEDVR_LEFT)


      ! ** deriv with eval_psi_full
      write(*,*) "psi_x with eval_dpsi_fedvr  (full)"
      do i=1,nmax_
         write(*,*) xg(i), psi(i), psi_x(i)
      enddo
      write(*,*)
!     pause



! *** Construct H = T + V and diagonalize with external legacy procedures
!      allocate(init1(nmax), init2(nmax))
      allocate( auxc(nmax), svec(nmax,3) )
      allocate(eigvec(nmax,nmax),xmat(nmax,nmax))
      allocate( wf0(nmax), wfc0(nmax), wf(nmax),wfc(nmax), wfc_(nmax) )
      call fedvr_hamilton_conf(jac,xa,wa,xs,xx,wx,eigval,          &
                         eigvec)

      write(*,*) "Eigenvalues of the time independent problem"
      do i=1,nmax
         write(*,*) i, eigval(i)
      enddo

       wfc0 = eigvec(:,1)/wx/dsqrt(jac)
       wf0 = wfc0 * wx * dsqrt(jac)
       auxc=wf0
       wf0 = matmul(transpose(eigvec),auxc)

       write(*,*) 
       do i=1,nmax
          write(*,*) i, real(wfc0(i)), real(wf0(i))
       enddo
       write(*,*) 
       write(*,*) sum( abs(wfc0)**2 * wx**2 * jac ) 
       write(*,*)     
!      pause


       allocate( workspace(nmax+2,2) )

       deallocate(psi, psi_x)
       allocate(psi(nmax), psi_x(nmax))
       psi=xx 
       write(*,*) "psi_x with eval_dpsi_fedvr_reduced)"
       call eval_dpsi_fedvr_reduced( snbr, nnbr, xs, map, Dref,        &
                                   psi, psi_x, workspace )


       write(*,*)
       do i=1,nmax
          write(*,*) i, xx(i), real(psi(i)), real(psi_x(i))
       enddo
       write(*,*)


! --- time parameters ---
       t_ini = -1.d0
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
       allocate(psi_ex(nmax), phi_ex(nmax))

       allocate(src_mid(nmax), src(nmax))
       allocate(src_x(nmax))

! --- initial condition ---
       tt = t0                ! start time
       psi_in = wf0           ! initial wavefunction


       psi_inx = matmul(eigvec,psi_in)
       psi_inx = psi_inx/wx/dsqrt(inv_jac)


!      write(*,*) 
!      do i=1,nmax
!         write(*,*) i, real(wfc0(i)), real(psi_inx(i))
!      enddo
!      write(*,*) 

       order = 4
       flag = 3
       j=1
       if (flag.eq.3) j=2

       phi_in = psi_in

       if (flag.eq.3) then
          call apply_momentum_operator(nmax, eigvec, xx, wx,  &
                  inv_jac, psi_in, phi_in, 2)
       end if

       phi0 = phi_in
       psi0 = psi_in
       phi_inx = phi_in
       psi_inx = psi_in

!      write(*,'(a)') '# i  E_n'
!      do i=1,10
!         write(*,*) i, eigval(i)
!      enddo


       write(*,*) "Check the initial condition"
       write(*,'(a)') '# t  Re[psi(1)] Im[psi(1)] ||psi||         '//  &
                      '         Re[phi(',j,')] Im[phi(',j,')] ||phi||'

             write(*,'(f8.4,3x,6e16.8)') tt, &
                  real(psi_in(j)), aimag(psi_in(j)),                &
                                       sqrt(sum(abs(psi_in)**2)),    &
                  real(phi_in(j)), aimag(phi_in(j)),                &
                                       sqrt(sum(abs(phi_in)**2))

      write(*,*) 
      write(*,*) 
      write(*,*) "Convergence test"

      call test_against_exact_solution( nmax, snbr, nnbr,             &
                                           inv_jac,                   &
                                           xs, xx, wx,                &
                                           map, Dref,                 &
                                           dt0, t_end, t_ini, omega,   &
                                           eigval, eigvec,             &
                                           psi0, phi0,                 &
                                           flag, order)

       write(*,*)
 
       call test_richardson_inhomogeneous( nmax, snbr, nnbr,   &
                                            inv_jac,                   &
                                            xs, xx, wx,                &
                                            map, Dref,                 &
                                            dt0, t_end, t_ini, omega,  &
                                            eigval, eigvec,            &
                                            psi0, phi0,                &
                                            flag, order )
 
 
      write(*,*) " "
      write(*,*) "exact vs numerical (during propag)"
      write(*,'(a)') '# t  Re[phi(1)] Im[phi(1)] ||phi||         '//  &
                     '           Re[phi_ex(1)] Im[phi_ex(1)] ||phi_ex||'

            write(*,'(f8.4,3x,6e16.8)') tt, &
                 real(psi_in(j)), aimag(psi_in(j)),                &
                                      sqrt(sum(abs(psi_in)**2)),    &
                 real(phi_in(j)), aimag(phi_in(j)),                &
                                      sqrt(sum(abs(phi_in)**2))

!           write(*,'(f8.4,3x,6e16.8)') tt, &
!                real(psi_in(j)), aimag(psi_in(j)),                &
!                                     sqrt(sum(abs(psi_in)**2)),    &
!                real(phi_in(j)), aimag(phi_ex(j)),                &
!                                     sqrt(sum(abs(phi_ex)**2))

      i =  0 
      write(*,*) nt, i, tt, phi_in(j), phi_inx(j)
      do i = 1, nt
!        write(*,*) nt, i, tt, phi_out(j), phi_ex(j)

         call process_src_ingredients ( nmax, snbr, nnbr,              &
                                   inv_jac,                            &
                                   xs, xx, wx, map, Dref,     &
                                   dt0, tt, omega,                     &
                                   eigval, eigvec, psi_in, svec,       &
                                   flag, order)

         call build_source_quadrature ( nmax, snbr, nnbr,             &
                                               xs, xx, map, Dref,     &
                                               dt0, tt, omega,        &
                                               eigval, eigvec,        &
                                               svec,                  &
                                               src,                   &
                                               order )
         
         call split_operator(nmax, dt0, tt, xx, eigval, eigvec,       &
                                            psi_in, psi_out, order)
         call split_operator(nmax, dt0, tt, xx, eigval, eigvec,       &
                                            phi_in, phi_out, order)


         !--------------------------------------------
         ! Add source contribution
         !--------------------------------------------
         phi_out = phi_out - ci * src


         tt = tt + dt0
         call exact_closed_duhamel(nmax, omega,                        &
                tt, t_ini, eigval, psi0, phi0, psi_ex, phi_ex, flag)


!!        write(*,*) "vs exact (during propag)" 
!!        write(*,*) nt, i, tt, psi_out(1), psi_ex(1)
          write(*,*) nt, i, tt, phi_out(j), phi_ex(j)
!         write(*,*) i, tt, psi_out(i),  phi_out(i)
!!        write(*,*) nt, i, psi_out(i), psi_ex(i), phi_out(i), phi_ex(i)

         psi_in = psi_out
         phi_in = phi_out
         psi_inx = psi_ex
         phi_inx = phi_ex
      enddo

!      write(*,*) 
!      write(*,*) "vs exact" 
!!     write(*,*) "Midpoint quadrature" 
!      write(*,*) "Simpson quadrature" 
!       do i=1,nmax
!!         write(*,*) nt, i, psi_out(i), psi_ex(i)
!          write(*,*) nt, i, psi_out(i), phi_out(i)
!!!        write(*,*) nt, i, psi_out(i), psi_ex(i), phi_out(i), phi_ex(i)
!       enddo
!
!!     write(*,*) 
!!     write(*,*) "check derivative p\psi(x,t)" 
!!     call apply_momentum_operator(nmax, eigvec, xx, wx,        &
!!                         inv_jac, psi0, psi_out, 2) 
!!     do i=1,nmax
!!        write(*,*) nt, i, psi0(i), psi_out(i)
!!     enddo

!      call apply_momentum_operator(nmax, eigvec, xx, wx,  &
!                 inv_jac, psi0, psi_out, 2)



          wfc = matmul(eigvec,psi_out)
          wfc = wfc / wx / dsqrt(jac)

          phic = matmul(eigvec,phi_out)
          phic = phic / wx / dsqrt(jac)
       
          call dump_wavefunction("final_psi.dat",nmax, eigvec, tt,    &
                        jac, wx, xx, psi_out)
          call dump_wavefunction("final_phi.dat",nmax, eigvec, tt,    &
                        jac, wx, xx, phi_out)
          call dump_wavefunction("final_psi_ex.dat",nmax, eigvec, tt,  &
                        jac, wx, xx, psi_ex)
          call dump_wavefunction("final_phi_ex.dat",nmax, eigvec, tt,  &
                        jac, wx, xx, phi_ex)
     

      end program
