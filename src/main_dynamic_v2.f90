      program fedvr_dynamic_build
      use omp_lib
      use constants, only : ppi, ci, omg_ => omega_qho
      use dynamic_parameters
      use math_util
      use util
      use propagation
      use observables
!     use conv_tests
      use io_module
      implicit none
      integer :: np, ns, nmax_, nobs, kobs, qho,                      &
                  lwork, info, store_val, ndim,                       &
                  i, j, k, ij, p, q, r, s, m, n,                      &
                  jn, ijk, ios, nch__, nch0, noc__, chan,             &
                  nsteps_, ntau__, src_type__, run__, order__

      real(8) ::    f0__, dt__, pfai__, omega__,                   &
                    step, eft, duration,                           &
                    start, finish, lap,                            &
                    aux, aux1, aux2,                               &
                    abstol,  qq,                                      &
                    tt, tp, tc, tmid,                     &
                    t1, t0, dt_,                          &
                    fmid,                                             &
                    p0, pexc, pion,                             &
                    err1, err2,                                       &
                    rowsum,                                    &
                    kappa_w,                               &
                    src_time, split_time,                             &
                    jacc, xx1, xx2, eps,                              &
                    omg_0, omg_1, dw__,                                &
                    norm_1, omega_k, omega_k0

      complex(8) :: cnum, a0, nrg_, xt_, pt_


      real(8), allocatable, dimension(:) :: xs, xx, wx,    &
                                   vec_matup, eigval,                 &
                                   kk, time_, time_t,                 &
                                   norm_,                &
                                   p0_t, pexc_t, pion_t,       &
                                   norm_ref, norm_refc, norm_saved,   &
                                   omega,     &
                                   omega_h

  
      complex(8), allocatable, dimension(:) :: wf1, wfc1,             &
                                               dwfc1,                 &
                                               wfc0_1, wf1_0,         &
                                               psi, psi0,             &
                                               psi_in,                &
                                               rr_, ak, cprob2,       &
                                               auxc1, auxc2,          &
                                               auxc3, auxc4,          &
                                               pk0, b0w_1, b0w_2,     &
                                               d_t, dd_t, d_w, dd_w,  &
                                               nrg_t, x_t, p_t

      integer, allocatable, dimension(:,:) :: map

      real(8), allocatable, dimension(:,:) :: lu, id, inv,            &
                                       pkk, bkw_1, bkw_2,             &
                                       Dref, Dglobal,                 &
                                       norm_t,                        &
                                       eigvec, basis

      real(8), allocatable, dimension(:,:,:) :: Dloc_all          

      complex(8), allocatable, dimension(:,:) :: svec, wfc,            &
                                                 wf_in, wf,    &
                                                 wf0, wfc0,            &
                                                 wf_saved,                &
                                                 src

      character(255) :: workdir, struct_dir, struct_dir_,             &
                        dyn_dir, dyn_dir_,                    &
                        multichain_tag


      logical :: ready_to_warp_up = .false.
      integer :: tid



!     include 'param_dynamic'
!     include 'param_exploit'
      include 'param_quantum_dynamic'

      integer :: log_unit
      integer :: obs_unit
      integer :: init_unit
      integer :: force_unit
      integer :: unit_pipe

      multichain_tag = "qdyn/"
      log_unit = 20
!     obs_unit = 40

      call cpu_time(start)
      qho = 0 
      eps = 1.e-6
      call get_command_argument(1, struct_dir)
      call get_command_argument(2, dyn_dir)
      call read_struct_bin(trim(struct_dir)//"struct.bin",        &
                             struct_dir_, nmax_, ns, np,             &
                             xx1, xx2, qq, jacc, xx, wx)

!     if (trim(struct_dir) /= trim(struct_dir_)) then
!        stop "Structure mismatch between dynamic and structure run"
!     end if

      call read_eigval_bin(trim(struct_dir)//"/eigval.bin", i, eigval)
      call read_eigvec_bin(trim(struct_dir)//"/eigvec.bin", j, eigvec)
      if (i.ne.j) stop
      if (i.ne.nmax_) stop

!     do i=1,10
!        write(*,*) i, eigval(i)
!     enddo
!     write(*,*)
!     do i=1,5
!        write(*,*) i, eigvec(:,1)
!     enddo



      write(log_unit,*) "Reading the dynamic parameters"
      call read_dynamic_bin(trim(dyn_dir)//"dynamic.bin",         &
                             dyn_dir_, struct_dir_,                    &
                             f0__, omega__, pfai__,                    &
                             t1, t0, nsteps_, dt__,                    &
                             noc__, ntau__, src_type__,                &
                             nch__, omg_max, omg_min,                  &
                             dw__, omega_h, run__,                    &
                             order__)

       call init_src(src_type__, nch_,                          &
                                omg_max, omg_min, dw_, run_)

!      if ((run-run__).ne.1) stop


       call set_force_params(f0__, omega__, pfai__)  
       call init_time_grid(noc__, ntau__, nsteps_)

       call set_resolution_order(order__)
!      call set_other_dyn_params(do_time_obs_, obs_stride_)



 
       allocate(wf_in(nmax_, nch), wf(nmax_, nch))

       allocate(wfc0(nmax_, nch), wf0(nmax_, nch))
       allocate(wfc(nmax_, nch))
       allocate(src(nmax_, nch))
       allocate(svec(nmax_,3))
       allocate(wf_saved(nmax_,nch__))
       allocate(omega(nch))

! --- initial condition ---


       ! --- φ channels ---
       if (src_type.eq.3) then
          !$omp parallel do default(shared) &
          !$omp private(k,i,omega_k,kappa_w)
          do k = 1, nch

             j = 1
             omega_k = omg_min + (k-1)* dw_

             omega(k) = omega_k
             kappa_w = varkap(kapp, omega_k)
          
             do i = 1, nmax_
                 wfc(i,k) = (-ci*kapp**(1.5d0)/omega_k)              &
                 * sgn(xx(i)) *                                       &
                 ( exp(-kapp*abs(xx(i))) - exp(-kappa_w*abs(xx(i))) )
             enddo
          enddo
          !$omp end parallel do
       end if

!      call dvr_to_eigen(nmax_, nch, jacc, wx, eigvec, wfc, wf_in)
       do k=1,nch
          call dvr_to_eigen(nmax_, jacc, wx, eigvec, wfc(:,k), wf(:,k))
          wf_in(:,k) = wf(:,k)
       enddo

       tt = t_ini                ! start time
       wf0 = wf_in          ! initial condition


!     do i=nmax_/2-5, nmax_/2+5
!        write(*,'(E20.10,*(1X,ES20.10))') xx(i), eigvec(i,1), wf1_0(i)
!     enddo
!      write(*,*)


       allocate(norm_ref(nch))
       allocate(norm_(nch), norm_saved(nch__))

       do k=1, nch
          norm_refc(k) = sqrt(sum(abs(wfc(:,k) * wx * dsqrt(jacc))**2))
          norm_(k) = sqrt(sum(abs(wf0(:,k))**2))
!         write(*,*) omega(k)
       enddo


       call init_run(workdir, multichain_tag, extract_name(dyn_dir))
!      call init_run(workdir, multichain_tag)


       open(newunit=init_unit, file=trim(workdir)//"initial_conds.dat", &
                                                status='replace')
       do i=1,nmax_
          write(init_unit,'(E20.10,*(1X,ES20.10))') xx(i), imag(wfc(i,:))
       enddo

       open(newunit=log_unit, file=trim(workdir)//"log.txt",           &
                                                status='replace')

!      do i=1,nmax_
!         write(init_unit,*) i, xx(i), real(wf1_0(i))
!      enddo


       write(*,*) "Saving the dynamic parameters"
       call write_dynamic_bin(trim(workdir)//"dynamic.bin",       &
                        workdir, struct_dir_,                     &
                        f0__, omega__, pfai__,                    &
                        t1, t0, nsteps_, dt__,                    &
                        noc__, ntau__, src_type__,                &
                        nch, omg_max, omg_min,                    &
                        dw_, omega, run,                        &
                        order)

 
        call write_struct_input(trim(workdir)//"param_structure.txt",   &
                                     struct_dir_, nmax_, ns, np,          &
                                     xx1, xx2, qq, jacc)
  
        call write_dynamic_input(trim(workdir)//"param_dynamic.txt",   &
                         workdir, struct_dir_,                         &
                         f0__, omega__, pfai__,                        &
                         t1, t0, nsteps_, dt__,                    &
                         noc__, ntau__, src_type__,                &
                         nch, omg_max, omg_min,                &
                         dw_, run,                           &
                         order)
  
  
  
        write(*,*) "Saving the initial conditions"
        call write_wavefun_bin(trim(workdir)//'initial_state.bin', 0,  &
                                       nch, nmax_, t_ini, omega, wf0)
        write(*,*) "initial conditions - saved"
        deallocate(wf0)
  
  
       write(*,*) "============================================"
       write(*,*) "Structure parameters"
       write(*,*) "============================================"
       write(*,*) "Np          = ", np 
       write(*,*) "Ns          = ", ns
       write(*,*) "N           = ", nmax_
       write(*,*) "xmax        = ", xx2 
       write(*,*) "xmin        = ", xx1
       write(*,*) "q           = ", qq
       write(*,*) "xrange      = ", xx2-xx1
       write(*,*) 
 
 
       call print_dynamic_parameters()
 
       write(*,*) "struct_dir  = ", trim(struct_dir)
       write(*,*) "dyn_dir     = ", trim(workdir)
 
 
       write(*,*) "Starting propagation"
       write(*,*) "nt =", nt
 
       open(newunit=force_unit, file=trim(workdir)//"force.dat",       &
                                                   status='replace')
       call plot_force(force_unit, t_end, t_ini, dt0)
       close(force_unit)
 
!      open(newunit=obs_unit, file=trim(workdir)//"norm.dat",           &
!                                                      status='replace')
!      nobs = nt / obs_stride
!      if (mod(nt, obs_stride) /= 0) nobs = nobs + 1
 
       allocate(time_(nobs))

       call read_wavefun_bin(trim(dyn_dir)//"initial_state.bin", &
                chan, ndim, tt, omega_h, wf_saved)

!      write(obs_unit,'(E20.10,*(1X,ES20.10))') tt,                   &
!                      norm_ref(:), norm_saved(:)

!      write(*,'(E20.10,*(1X,ES20.10))') tt,                   &
!                      norm_ref(:), norm_saved(:)
!      allocate(norm_t(nobs,nch))
!      allocate(p0_t(nobs), pexc_t(nobs), pion_t(nobs))
!      allocate(nrg_t(nobs), x_t(nobs), p_t(nobs))
 
       kobs = 0
       write(*,*) "nobs = ", nobs
 
 
       src_time = 0.d0
       split_time = 0.d0
 
!!     omega(1)=0.d0

       write(*,*) "dt0 = ", dt0

       open(newunit=unit_pipe,                                      &
            file=trim(dyn_dir)//"wf_psi_pipe.bin",                  &
            form="unformatted",                                     &
            status="old",                                           &
            action="read")
 
 


          !$omp parallel private(k) default(shared) 
           do i=1, nt
           ! --- propagation ---

             !$omp single
 
              read(unit_pipe) chan
              read(unit_pipe) ndim
              read(unit_pipe) tt
              read(unit_pipe) omega_k0
              read(unit_pipe) wf_saved

             do j=1,nch__
                norm_saved(j) = sqrt(sum(abs(wf_saved(:,j))**2))
             enddo
  
             write(obs_unit,'(E20.10,*(1X,ES20.10))') tt,               &
                                             norm_(:), norm_saved(:)
    
              
              call process_src_ingredients ( nmax_, ns, np,         &
                                    jacc,                           &
                                    xs, xx, wx, map, Dref,     &
                                    dt0, tt,                        &
                                    eigval, eigvec,                 &
                                    wf_saved(:,1), svec,   &
                                    src_type, order)

             !$omp end single
          

             !$omp do
              do k=1, nch

                 call build_source_quadrature (   nmax_, ns, np,    &
                                             xs, xx, map, Dref,     &
                                             dt0, tt,               &
                                             eigval, eigvec,        &
                                             svec,                  &
                                             src(:,k), omega(k),    &
                                             order )


                 call split_operator(nmax_, dt0, tt, xx, eigval, eigvec,   &
                                            wf_in(:,k), wf(:,k), order)
  
  
  

                 ! Add source contribution
                 wf(:,k) = wf(:,k) - ci * src(:,k)


                 norm_(k) = sqrt(sum(abs(wf(:,k))**2))
                 wf_in(:,k) = wf(:,k)
              enddo
            !$omp end do

       end do
      !$omp end parallel
 
 
 
      call read_observables_bin(trim(dyn_dir)//"/dyn_obs.bin",     &
                 nch0, nobs, time_t, norm_t, p0_t, pexc_t, pion_t,   &
                 nrg_t, x_t, p_t)

       call write_observables_bin(trim(workdir)//"dyn_obs.bin", &
                  nch__, nobs, time_, norm_t,                 &
                  p0_t, pexc_t, pion_t,                          &
                  nrg_t, x_t, p_t)
 
       call write_wavefun_bin(trim(workdir)//'wavfun.bin', 0,   &
                                 nch, nmax_, tt, omega, wf)
 
       close(obs_unit)
       write(log_unit,*) workdir


      end program
