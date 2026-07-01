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
                  i, j, k, ij, p, r, s, m, n,                      &
                  jn, ijk,                                            &
                  ios, order, flag,                                   &
                  nmax_

      real(8) ::    start, finish, lap, aux, aux1, aux2,              &
                    abstol,                                           &
                    err1, err2,                                       &
                    xc_in, rowsum



      real(8), allocatable, dimension(:) :: xa, wa, xs, xx, wx, xg,   &
                                   vec_matup, eigval,                 &
                                   auxr1, auxr2, auxr3,               &
                                   matint_v, vext_mid, vext_,         &
                                   ipiv, work, ifail, dnrm2
  
      complex(8), allocatable, dimension(:) :: wf0, wfc0,            &
                                               dwf0, dwfc0,           &
                                               wf, wfc, wfc_,         &
                                               psi_exact,             &
                                               psi0, psi, psi_x,      &
                                               psi_in, psi_out,       &
                                               init1, init2,          &
                                               psi_ex, phi_ex,        &
                                               auxc

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



      struct_tag = "str/"
      call cpu_time(start)
      call set_grid_params(np, ns, xc_in, xmin_, xmax_, q)
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

      allocate(eigval(nmax))
      allocate(xx(nmax),wx(nmax))


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
                                    xmin, xmax, jac, q, xx, wx)

      call write_problem_input(trim(workdir)//"struct_params.dat",   &
                                    workdir, nmax, snbr, nnbr,       &
                                    xmin, xmax, jac, q)


      write(*,*) 
      write(*,*) 
      write(*,*) workdir

      open(newunit=fundamental_unit,                              &
                           file=trim(workdir)//"fundamental.dat", &
                           status='replace')

      allocate( wf0(nmax), wfc0(nmax), wf(nmax), wfc(nmax), wfc_(nmax) )
      wf0 = (0.0d0,0.d0)
      wf0(1) = (1.0d0,0.d0)

      call eigen_to_dvr(nmax, jac, wx, eigvec, wf0, wfc0)

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


      end program
