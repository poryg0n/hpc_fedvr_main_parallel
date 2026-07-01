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
                  lwork, info, store_val,                        &
                  i, j, k, ij, p, q, r, s, m, n,                      &
                  jn, ijk, ios, nch__, run_h
!                 nmax, nmax_, krange, ios, nt, order, src_type

      real(8) ::    step, eft, duration,                           &
                    start, finish, lap, aux, aux1, aux2,              &
                    abstol,  qq,                                      &
                    tt, tp, tc, tmid,                     &
                    fmid,                                             &
                    p0, pexc, pion,                             &
                    err1, err2,                                       &
                    rowsum,                                    &
                    jacc, xx1, xx2, eps,                              &
                    norm_1

      complex(8) :: cnum, a0, nrg_, xt_, pt_


      real(8), allocatable, dimension(:) :: xs, xx, wx,    &
                                   vec_matup, eigval,                 &
                                   kk, time_, time_t,                 &
                                   norm_,                &
                                   p0_t, pexc_t, pion_t,       &
                                   norm_ref, norm_refc, omega

  
      complex(8), allocatable, dimension(:) :: wf1, wfc1,             &
                                               dwfc1,                 &
                                               wfc0_1, wf1_0,         &
                                               wfc_in, wf_in,         &
                                               wf_out,                &
                                               auxc1, auxc2,          &
                                               auxc3, auxc4,          &
                                               d_t, dd_t, d_w, dd_w,  &
                                               nrg_t, x_t, p_t

      integer, allocatable, dimension(:,:) :: map

      real(8), allocatable, dimension(:,:) :: lu, id, inv,            &
                                       pkk, bkw_1, bkw_2,             &
                                       Dref, Dglobal,                 &
                                       norm_t,                      &
                                       eigvec, basis

      real(8), allocatable, dimension(:,:,:) :: Dloc_all          

      complex(8), allocatable, dimension(:,:) :: in_states,           &
                                                 out_states,          &
                                                 svec,                &
                                                 wf, wfc,             &
                                                 wf0, wfc0,           &
                                                 cwf

      character(255) :: workdir, struct_dir, struct_dir_,             &
                        dyn_tag

      integer :: tid


      include 'param_dynamic'

      integer :: log_unit
      integer :: obs_unit
      integer :: init_unit
      integer :: force_unit

      dyn_tag = "homo_dyn/"
      log_unit = 20
      obs_unit = 40
      nch__=1
      run_h=0

      call cpu_time(start)
      qho = 0 
      eps = 1.e-6
      call get_command_argument(1, struct_dir)
      call read_structure_bin(trim(struct_dir)//"/structure.bin",    &
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


       call set_force_params(f0_, omega_, pfai_)  
       call init_time_grid(noc_, ntau_, nsteps)
       call init_src(src_type_, nch_,                          &
                                omg_0, omg_1, wstep_, run_h)

       call set_resolution_order(order_)
       call set_other_dyn_params(do_time_obs_, obs_stride_)


       duration = int((t_end - t_ini))
!      allocate(fvec(nt))
       tc = 0.d0
 
       allocate(wfc0_1(nmax_), wf1_0(nmax_))
       allocate(wf_in(nmax_), wf_out(nmax_))

       allocate(wfc0(nmax_,nch__), wf0(nmax_,nch__))
       allocate(wfc(nmax_,nch__), wf(nmax_,nch__))

! --- initial condition ---

       do i=1, nmax_
          wf1_0(i) = kapp**(1.d0/2) * exp(-kapp*abs(xx(i)))
       enddo

!      call dvr_to_eigen(nmax_, 1, jacc, wx, eigvec, wf1_0, wf_in)
       do k=1,nch
          call dvr_to_eigen(nmax_, jacc, wx, eigvec, wf1_0, wf_in)
       enddo

       tt = t_ini                ! start time
       wf0(:,1) = wf_in          ! initial condition


!     do i=nmax_/2-5, nmax_/2+5
!        write(*,'(E20.10,*(1X,ES20.10))') xx(i), eigvec(i,1), wf1_0(i)
!     enddo
!     write(*,*)


      allocate(norm_refc(nch__))
      allocate(norm_ref(nch__))
      allocate(omega(nch__))

      do k=1, nch__
         norm_refc(k) = sqrt(sum(abs(wf1_0    * wx * dsqrt(jacc))**2))
         norm_ref(k) = sqrt(sum(abs(wf_in * wx * dsqrt(jacc))**2))
      enddo

!     norm_ref1(1) = sqrt(sum(abs(wf_in)**2))
!     norm_ref1(2) = sqrt(sum(abs(phi_in)**2))


      call init_run(workdir, dyn_tag, extract_name(struct_dir))

      open(newunit=log_unit, file=trim(workdir)//"log.txt",            &
                                               status='replace')

      open(newunit=init_unit, file=trim(workdir)//"initial_conds.dat", &
                                               status='replace')
      do i=1,nmax_
         write(init_unit,*) i, xx(i), real(wf1_0(i))
      enddo


      write(*,*) "Write structure"
      call write_structure_bin(trim(workdir)//"structure.bin",       &
                                    workdir, nmax_, ns, np,   &
                                    xx1, xx2, qq, jacc, xx, wx)


      write(*,*) "Saving the dynamic parameters"
      call write_dynamic_bin(trim(workdir)//"dynamic.bin",             &
                       workdir, struct_dir_,                           &
                       f0, omega0, pfai,                               &
                       t_end, t_ini, nt, dt0,                          &
                       noc, ntau, src_type,                            &
                       nch, omg_max, omg_min,                          &
                       wstep, omega, run_h,                            &
                       order)


      call write_structure_input(trim(workdir)//"param_structure.txt", &
                                   struct_dir, nmax_, ns, np,          &
                                   xx1, xx2, qq, jacc)

      call write_dynamic_input(trim(workdir)//"param_dynamic.txt",     &
                       workdir, struct_dir_,                           &
                       f0, omega0, pfai,                               &
                       t_end, t_ini, nt, dt0,                          &
                       noc, ntau, src_type,                            &
                       nch, omg_max, omg_min,                  &
                       wstep, run_h,                           &
                       order)




      write(*,*) "Saving the initial conditions"
      call write_wavefun_bin(trim(workdir)//'initial_state.bin', 0,    &
                                     1, nmax_, t_ini, omega, wf0)
      write(*,*) "initial conditions - saved"


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


!     if (do_conv_test) then
!        write(*,*) "performing the richardson test"
!        call test_richardson_inhomogeneous(log_unit, nmax_, ns, np,   &
!                                             jacc,                    &
!                                             xs, xx, wx,              &
!                                             map, Dref,               &
!                                             nt/4, omega0,            &
!                                             eigval, eigvec,          &
!                                             psi0, phi0,              &
!                                             src_type, omeg, order)

!        write(*,*) "richardson test done"
!     end if



      write(*,*) "Starting propagation"
      write(*,*) "nt =", nt

      open(newunit=force_unit, file=trim(workdir)//"force.dat",       &
                                                  status='replace')
      call plot_force(force_unit, t_end, t_ini, dt0)
      close(force_unit)

      open(newunit=obs_unit, file=trim(workdir)//"norm.dat",           &
                                                      status='replace')
      nobs = nt / obs_stride
      if (mod(nt, obs_stride) /= 0) nobs = nobs + 1

      allocate(time_(nobs), norm_(nch__))
      allocate(norm_t(nobs,nch__))
      allocate(p0_t(nobs), pexc_t(nobs), pion_t(nobs))
      allocate(nrg_t(nobs), x_t(nobs), p_t(nobs))

      kobs = 0
      write(*,*) "nobs = ", nobs


      omega(1)=0.d0


      do i=1,nt

         ! --- propagation ---
            call split_operator(nmax_, dt0, tt, xx, eigval, eigvec,   &
                                            wf_in, wf_out, order)


         norm_1 = sqrt(sum(abs(wf_out)**2))

         write(obs_unit,'(4ES20.10)') tt, norm_1

         tt = tt + dt0
         wf_in = wf_out
         wf(:,1) = wf_out
         call write_wavefun_bin(trim(workdir)//'wf_psi_pipe.bin', 1,   &
                                 1, nmax_, tt, omega, wf)


         ! --- observables ---
         if (do_time_obs) then

            if (mod(i, obs_stride) == 0) then

               kobs = kobs + 1
            
               call compute_dyn_observables(nmax_,                  &
                                            xx, wx, jacc,              &
                                            eigval, eigvec,            &
                                            norm_,                     &
                                            p0, pexc, pion,      &
                                            nrg_, xt_, pt_,      &
                                            wf)
            
               call append_dyn_obs_bin(trim(workdir)//"dyn_back.bin",  &
                                   nch__, tt, norm_,                   &
                                   p0, pexc, pion,                     &
                                   nrg_, xt_, pt_ )

               write(log_unit,*) kobs, nobs, tt,                       &
                                        norm_(nch__),                  &
                                        p0, pexc, pion,                &
                                        nrg_, xt_, pt_

               ! --- STORE ---
               time_(kobs)         = tt
               norm_t(kobs,nch__)  = norm_(nch__)
               p0_t(kobs)          = p0
               pexc_t(kobs)        = pexc
               pion_t(kobs)        = pion
               nrg_t(kobs)         = real(nrg_)
               x_t(kobs)           = real(xt_)
               p_t(kobs)           = real(pt_)
            
            end if


         end if


         if (mod(i,100).eq.0) then
            call write_wavefun_bin(trim(workdir)//'wavfun_back.bin', 0, &
                                  1, nmax_, tt,  omega, wf)
         end if


      end do



      call write_observables_bin(trim(workdir)//"dyn_obs.bin", &
                nch__, nobs, time_, norm_t,                 &
                p0_t, pexc_t, pion_t,                          &
                nrg_t, x_t, p_t)

      call write_wavefun_bin(trim(workdir)//'wavfun_psi.bin', 0,   &
                                1, nmax_, tt, omega, wf)

      close(obs_unit)
      write(log_unit,*) workdir
      close(log_unit)


      end program
