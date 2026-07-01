      program fedvr_dynamic_build
      use omp_lib
      use constants, only : ppi, ci, omg_ => omega_qho
      use dynamic_parameters
      use math_util
      use util
      use propagation
      use observables
      use conv_tests
      use io_module
      implicit none
      integer :: np, ns, nmax_, nobs, kobs, qho,                      &
                  lwork, info, store_val,                        &
                  i, j, k, ij, p, q, r, s, m, n,                      &
                  jn, ijk, ios
!                 nmax, nmax_, krange, ios, nt, order, src_type

      real(8) ::    step, eft, duration,                           &
                    start, finish, lap, aux, aux1, aux2,              &
                    abstol,  qq,                                      &
                    tt, tp, tc, tmid,                     &
                    fmid,                                             &
                    norm_1, norm_2,                                   &
                    p0, pexc, pion,                             &
                    err1, err2,                                       &
                    rowsum,                                    &
                    kappa_w,                               &
                    src_time, split_time,                             &
                    jacc, xx1, xx2, eps, norm_ref

      complex(8) :: cnum, a0, nrg_, xt_, pt_


      real(8), allocatable, dimension(:) :: xs, xx, wx,    &
                                   vec_matup, eigval,                 &
                                   kk, time_, time_t,                 &
                                   norm_t1, norm_t2,                  &
                                   norm_ref1, norm_ref2,       &
                                   p0_t, pexc_t, pion_t

  
      complex(8), allocatable, dimension(:) :: wf0, wfc0,            &
                                               dwf0, dwfc0,           &
                                               wf, wfc, wfc_,         &
                                               psi_exact,             &
                                               psi0, psi, psi_x,      &
                                               psi_in, psi_out,       &
                                               phi_in, phi_out,       &
                                               phi_inc,       &
                                               psi_inx, psi_outx,     &
                                               phi_inx, phi_outx,     &
                                               src_out, src_mid,      &
                                               src, src_x,            &
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
                                               pk0, b0w_1, b0w_2,     &
                                               d_t, dd_t, d_w, dd_w,  &
                                               nrg_t, x_t, p_t

      integer, allocatable, dimension(:,:) :: map

      real(8), allocatable, dimension(:,:) :: lu, id, inv,            &
                                       pkk, bkw_1, bkw_2,             &
                                       Dref, Dglobal,                 &
                                       eigvec, basis

      real(8), allocatable, dimension(:,:,:) :: Dloc_all          

      complex(8), allocatable, dimension(:,:) :: in_states,           &
                                                 out_states,          &
                                                 svec,                &
                                                 cwf

      character(255) :: workdir, struct_dir, struct_dir_,             &
                        dyn_tag


      logical :: ready_to_warp_up = .false.
      integer :: tid



      include 'param_dynamic'

      integer :: log_unit
      integer :: obs_unit
      integer :: init_unit
      integer :: force_unit

      dyn_tag = "dyn/"
      log_unit = 20
      obs_unit = 40

      call cpu_time(start)
      qho = 0 
      eps = 1.e-6
      call get_command_argument(1, struct_dir)
      call read_problem_bin(trim(struct_dir)//"/problem.bin",        &
                             struct_dir_, nmax_, ns, np,             &
                             xx1, xx2, jacc, qq, xx, wx)

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


       call set_force_params(f0_, omega_0, pfai_)  
       call init_time_grid(noc_, ntau_, nsteps)
       call init_src(src_type_, omega_)
       call set_resolution_order(order_)
       call set_other_dyn_params(do_time_obs_, obs_stride_)


       duration = int((t_end - t_ini))
!      allocate(fvec(nt))
       tc = 0.d0
 

       allocate( auxc(nmax_), svec(nmax_,3) )
 
       allocate(psi0(nmax_), phi0(nmax_))
       allocate(wfc0(nmax_), wf0(nmax_))
       allocate(psi_in(nmax_), psi_out(nmax_))
       allocate(psi_inx(nmax_), psi_outx(nmax_))
       allocate(phi_in(nmax_), phi_out(nmax_))
       allocate(phi_inx(nmax_), phi_outx(nmax_))
       allocate(phi_inc(nmax_))
       allocate(psi_exact(nmax_))
       allocate(psi_ex(nmax_), phi_ex(nmax_))
 
       allocate(src_mid(nmax_), src(nmax_))
       allocate(src_x(nmax_))

! --- initial condition ---
!      wfc0 = eigvec(:,1)
!      wfc0 = eigvec(:,1)/wx/dsqrt(jacc)
!      psi0 = matmul(transpose(eigvec),wfc0)

!      wfc0 =  kapp**(.5d0) * exp(-kapp*abs(xx))
!      wf0 = wfc0 * wx * dsqrt(jacc)
!      wf0 = matmul(transpose(eigvec),wf0)

       
       
       j=1
       psi0    = (0.d0,0.d0)
       psi0(1) = (1.d0,0.d0)
       wf0 = psi0


       if (src_type.eq.3) then
          j = 2
!         call apply_momentum_operator(nmax_, eigvec, xx, wx,  &
!                 jacc, psi_in, phi_in, qho)

          kappa_w = varkap(kapp, omeg)
          do i=1, nmax_
             phi_inc(i) = (-ci*kapp**(3.d0/2)/omeg) * sgn(xx(i)) *    &
              (  exp(-kapp*abs(xx(i)))  -  exp(-kappa_w*abs(xx(i))) ) 
          enddo

       end if


       do i=1, nmax_
          wfc0(i) = kapp**(1.d0/2) * exp(-kapp*abs(xx(i)))
       enddo

       call dvr_to_eigen(nmax_, jacc, wx, eigvec, wfc0, psi0)
       call dvr_to_eigen(nmax_, jacc, wx, eigvec, phi_inc, phi0)

       tt = t_ini                ! start time
       psi_in = psi0                ! initial wavefunction
       phi_in = phi0                ! initial wavefunction


!     do i=nmax_/2-5, nmax_/2+5
!        write(*,*) xx(i), eigvec(i,1), wfc0(i)
!     enddo
!     write(*,*)
      do i=1, 10
         write(*,'(2E20.10)') psi0(i), phi0(i)
      enddo

!     write(*,*) "src_type is ", src_type
!     pause


      psi_inx = psi_in
      phi_inx = phi_in

      allocate(norm_ref1(2),norm_ref2(2))

      norm_ref2(1) = sqrt(sum(abs(wfc0    * wx * dsqrt(jacc))**2))
      norm_ref2(2) = sqrt(sum(abs(phi_inc * wx * dsqrt(jacc))**2))

      norm_ref1(1) = sqrt(sum(abs(psi_in)**2))
      norm_ref1(2) = sqrt(sum(abs(phi_in)**2))



      write(*,*)
      write(*,*) norm_ref2(1), norm_ref1(1)
      write(*,*) norm_ref2(2), norm_ref1(2)


      call init_run(workdir, dyn_tag, extract_name(struct_dir))

      open(newunit=log_unit, file=trim(workdir)//"log.txt",            &
                                               status='replace')

      open(newunit=init_unit, file=trim(workdir)//"initial_conds.dat", &
                                               status='replace')
      do i=1,nmax_
         write(init_unit,'(E20.10,*(1X,(ES20.10)))') xx(i),            &
                         real(wfc0(i)), imag(phi_inc(i)), omeg
      enddo

!    pause



      write(*,*) "Saving the dynamic parameters"
      call write_dynamic_bin(trim(workdir)//"dynamic.bin",             &
                       workdir, struct_dir_,                           &
                       f0, omega0, pfai,                                &
                       t_end, t_ini, nt, dt0,                          &
                       noc, ntau, src_type,                            &
                       omeg, order)


      call write_problem_input(trim(workdir)//"param_structure.txt",   &
                                   struct_dir, nmax_, ns, np,          &
                                   xx1, xx2, jacc, qq)

      call write_dynamic_input(trim(workdir)//"param_dynamic.txt",     &
                       workdir, struct_dir_,                           &
                       f0, omega0, pfai,                               &
                       t_end, t_ini, nt, dt0,                          &
                       noc, ntau, src_type,                            &
                       omeg, order)



      write(*,*) "Saving the initial conditions"
      call write_wavefunction_bin(trim(workdir)//'initial_state.bin',  &
                                   nmax_, t_ini, omeg ,psi_in, phi_in)
      write(*,*) "initial conditions - saved"


      write(*,*) "============================================"
      write(*,*) "Structure parameters"
      write(*,*) "============================================"
      write(*,*) "Np          = ", np 
      write(*,*) "Ns          = ", ns
      write(*,*) "N           = ", nmax_
      write(*,*) "xmax        = ", xx2 
      write(*,*) "xmin        = ", xx1

      call print_dynamic_parameters()

      write(*,*) "struct_dir  = ", trim(struct_dir)
      write(*,*) "dyn_dir     = ", trim(workdir)


      if (do_conv_test) then
         write(*,*) "performing the richardson test"
         call test_richardson_inhomogeneous(log_unit, nmax_, ns, np,   &
                                              jacc,                    &
                                              xs, xx, wx,              &
                                              map, Dref,               &
                                              nt/4, omega0,            &
                                              eigval, eigvec,          &
                                              psi0, phi0,              &
                                              src_type, omeg, order)

         write(*,*) "richardson test done"
      end if



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

      allocate(time_(nobs), norm_t1(nobs), norm_t2(nobs))
      allocate(p0_t(nobs), pexc_t(nobs), pion_t(nobs))
      allocate(nrg_t(nobs), x_t(nobs), p_t(nobs))

      kobs = 0
      write(*,*) "nobs = ", nobs
!     pause


      i=1

      src_time = 0.d0
      split_time = 0.d0

      ready_to_warp_up = .false.



      !$omp parallel private(tid) default(shared) num_threads(3)
      do i=1,nt

         tid = omp_get_thread_num()
!        write(*,*) i, tid
!        write(*,*) 

         if (tid.eq.0) then
            call process_src_ingredients ( nmax_, ns, np,         &
                               jacc,                              &
                               xs, xx, wx, map, Dref,             &
                               dt0, tt,                           &
                               eigval, eigvec, psi_in, svec,      &
                               src_type, order)

            call build_source_quadrature (   nmax_, ns, np,       &
                                     xs, xx, map, Dref,     &
                                         dt0, tt,               &
                                         eigval, eigvec,        &
                                         svec,                  &
                                         src, omeg,    &
                                         order )

         end if
         !$omp barrier

         if (tid.eq.1) then
            call split_operator(nmax_, dt0, tt, xx, eigval, eigvec,   &
                                         psi_in, psi_out, order)

            norm_1 = sqrt(sum(abs(psi_out)**2))
            psi_in = psi_out

         else if (tid.eq.2) then
            call split_operator(nmax_, dt0, tt, xx, eigval, eigvec,    &
                                         phi_in, phi_out, order)

            phi_out = phi_out - ci * src
            norm_2 = sqrt(sum(abs(phi_out)**2))
            phi_in = phi_out
         end if

         !$omp barrier


!        !$omp single
         if (tid.eq.0) then

            write(obs_unit,'(4ES20.10)') tt, norm_1, norm_2

            tt = tt + dt0


            ! --- observables ---
            if (do_time_obs) then

               if (mod(i, obs_stride) == 0) then

                  kobs = kobs + 1
               
                  call compute_dyn_observables(nmax_,                     &
                                               psi_out, phi_out,          &
                                               xx, wx, jacc,              &
                                               eigval, eigvec,            &
                                               norm_1, norm_2,            &
                                               p0, pexc, pion,      &
                                               nrg_, xt_, pt_)
               
                  call append_dyn_obs_bin(trim(workdir)//"dyn_back.bin",  &
                                      tt, norm_1, norm_2,           &
                                      p0, pexc, pion,                     &
                                      nrg_, xt_, pt_ )

                  write(log_unit,*) kobs, nobs, tt,                       &
                                           norm_1, norm_2,                &
                                           p0, pexc, pion,                &
                                           nrg_, xt_, pt_

                  ! --- STORE ---
                  time_(kobs)   = tt
                  norm_t1(kobs) = norm_1
                  norm_t2(kobs) = norm_2
                  p0_t(kobs)    = p0
                  pexc_t(kobs)  = pexc
                  pion_t(kobs)  = pion
                  nrg_t(kobs)   = real(nrg_)
                  x_t(kobs)     = real(xt_)
                  p_t(kobs)     = real(pt_)
               
               end if


            end if

!           call exact_closed_duhamel(nmax_, tt, t_ini,                   & 
!                   eigval, psi0, phi0, psi_ex, phi_ex, src_type, omeg)


!           write(*,*) nt, i, tt, phi_out(j), phi_ex(j)
            if (mod(i,100).eq.0) then
               call write_wavefunction_bin(trim(workdir)//'wavfun.bin',   &
                                        nmax_, tt, omeg, psi_out, phi_out)
            end if


         end if
!        !$omp end single
      end do
      !$omp end parallel 




      call write_observables_bin(trim(workdir)//"dyn_obs.bin", &
                nobs, time_, norm_t1, norm_t2,                 &
                p0_t, pexc_t, pion_t,                          &
                nrg_t, x_t, p_t)

      call write_wavefunction_bin(trim(workdir)//'wavfun.bin',   &
                                nmax_, tt, omeg, psi_out, phi_out)

      close(obs_unit)
      write(log_unit,*) workdir
      close(log_unit)


      end program
