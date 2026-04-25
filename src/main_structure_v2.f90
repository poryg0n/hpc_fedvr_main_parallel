      program fedvr_build
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
      use space_time_ops
      use io_module
      implicit none
      integer :: lwork, info, store_val,                              &
                  i, j, k, ij, p, q, r, s, m, n,                      &
                  jn, ijk,                                            &
                  ios, order, flag,                                   &
                  nmax_
!                 nmax, nmax_, krange, ios, nt, order, flag

      real(8) ::    start, finish, lap, aux, aux1, aux2,              &
                    abstol,                                           &
                    err1, err2,                                       &
                    xc_in, rowsum

      complex(8) :: cnum, a0, cprobk, auxcc, auxcc1, auxcc2


      real(8), allocatable, dimension(:) :: xa, wa, xs, xx, wx, xg,   &
                                   vec_matup, eigval,                 &
                                   kk, xt,                            & 
                                   auxr1, auxr2, auxr3,               &
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
                                               wwfc, dwwfc,           &
                                               dwfc, dwfc2, dwfc3,    & 
                                               init1, init2,          &
                                               phi0, phic0,           &
                                               psi_test, phi_test,    &
                                               psi_ex, phi_ex,        &
                                               phi, phic,             &
                                               auxc,                  &
                                               auxc1, auxc2,          &
                                               auxc3, auxc4,          &
                                               x_t,                   &
                                               p_t, pp_t, p_w, pp_w,  &
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

      character(255) :: workdir, struct_tag


      include 'param_structure'

      integer :: log_unit
      integer :: obs_unit
      integer :: eigval_unit
      integer :: eigvec_unit
      integer :: fundamental_unit



      struct_tag = "struct/"
      call cpu_time(start)
      call set_grid_params(np, ns, xc_in, xmin_, xmax_)
!     call init_run(workdir, "E0.5_dt0.01")
      call init_run(workdir, struct_tag)


      allocate(xa(nnbr), wa(nnbr), xs(1:snbr+1))
      call lobatto_compute(nnbr, xa, wa) 

      allocate(Dref(nnbr,nnbr)) 
      call build_Dref_lobatto(nnbr, xa, Dref)

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
!     write(*,*) "nmax_", nmax_

!     if (nmax.ne.(nmax_-2)) stop

      allocate(Dglobal(nmax_,nmax_))
      call assemble_derivative_fedvr(snbr, nnbr, map, wa,             &
!                                 Dloc_all, Dglobal, FEDVR_SYMMETRIC )
                                   Dloc_all, Dglobal, FEDVR_MUTUAL )
!                                  Dloc_all, Dglobal, FEDVR_LEFT )

!     allocate(xg(nmax_))
!     call build_x_global_fedvr(snbr, nnbr, xa, xs, map, xg)

!     do i = 2, nmax_
!        if (xg(i) <= xg(i-1)) then
!           write(*,*) "Grid ordering error at", i
!           stop
!        end if
!     end do

!     allocate(psi(nmax_), psi_x(nmax_))

      allocate(eigval(nmax), ss(nmax))
      allocate(xx(nmax),wx(nmax))


!     ss = 1.d0
!     psi = xg
!     psi_x = matmul(Dglobal, psi)

!      ! ** deriv with operator D
!      write(*,*) "psi_x with D operator, blunt"
!      do i=1,nmax_
!         write(*,*) xg(i), psi(i), psi_x(i)
!      enddo
!      write(*,*) 
 
!      call eval_dpsi_fedvr( snbr, nnbr, xg, xs, map, Dref,    &
!                                             psi, psi_x, FEDVR_LEFT)


!      ! ** deriv with eval_psi_full
!      write(*,*) "psi_x with eval_dpsi_fedvr  (full)"
!      do i=1,nmax_
!         write(*,*) xg(i), psi(i), psi_x(i)
!      enddo
!      write(*,*)
!!     pause



      call print_structure_parameters()

! *** Construct H = T + V and diagonalize with external legacy procedures
!     allocate(init1(nmax), init2(nmax))
      allocate(eigvec(nmax,nmax))
      call fedvr_hamilton_conf(jac, xa, wa, xs, xx, wx,                &
                                                  eigval, eigvec)

      call write_eigval_bin(trim(workdir)//"eigval.bin", nmax, eigval)
      call write_eigvec_bin(trim(workdir)//"eigvec.bin", nmax, eigvec)
      call write_problem_bin(trim(workdir)//"problem.bin",       &
                                    workdir, nmax, snbr, nnbr,   &
                                    xmin, xmax, jac, xx, wx)

      call write_problem_input(trim(workdir)//"struct_params.dat",   &
                                    workdir, nmax, snbr, nnbr,       &
                                    xmin, xmax, jac)




      write(*,*) 
      write(*,*) 
      write(*,*) workdir

      open(newunit=fundamental_unit,                              &
                           file=trim(workdir)//"fundamental.dat", &
                           status='replace')

      allocate( wf0(nmax), wfc0(nmax), wf(nmax),wfc(nmax), wfc_(nmax) )
      wf0 = (0.0d0,0.d0)
      wf0(1) = (1.0d0,0.d0)

!     call eigen_to_dvr(nmax, jac, wx, eigvec, wf0, wfc0)
      wfc0 = matmul(eigvec, wf0)
      wfc0 = wfc0 /wx / dsqrt(jac)


      write(*,*) "Norm of the fundamental", sum( abs(wfc0)**2 * wx**2 * jac ) 
      write(*,*) "writing fundamental to file"
      do i=1,nmax
         write(fundamental_unit,*) xx(i),                              &
                  real(wfc0(i)), eigvec(i,1)/wx(i)/dsqrt(jac), wf0(i)
      enddo
      close(fundamental_unit)
      write(*,*) "done"


      open(newunit=eigval_unit, file=trim(workdir)//"eigval.dat",            &
                                                  status='replace')

      write(*,*) "writing eigenvalues to file"
      do i=1,nmax
         write(eigval_unit,*) i, eigval(i)
      enddo
      close(eigval_unit)
      write(*,*) "done"


!      allocate( workspace(nmax+2,2) )

!      deallocate(psi, psi_x)
!      allocate(psi(nmax), psi_x(nmax))
!      psi=xx 
!      write(*,*) "psi_x with eval_dpsi_fedvr_reduced)"
!      call eval_dpsi_fedvr_reduced( snbr, nnbr, xs, map, Dref,        &
!                                  psi, psi_x, workspace )


!      write(*,*)
!      do i=1,nmax
!         write(*,*) i, xx(i), real(psi(i)), real(psi_x(i))
!      enddo
!      write(*,*)


      end program
