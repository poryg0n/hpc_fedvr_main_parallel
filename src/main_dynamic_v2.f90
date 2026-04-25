      program fedvr_dynamic_build
      use omp_lib
      use constants, only : ppi, ci, omg_ => omega_qho
      use dynamic_parameters
      use math_util
      use space_time_ops
      use util
      use propagation
      use conv_tests
      use observables
      use io_module
      implicit none
      integer :: np, ns, nmax_, nobs, kobs, qho,                      &
                  lwork, info, store_val,                        &
                  i, j, k, ij, p, q, r, s, m, n, w,                   &
                  jn, ijk, ios
!                 nmax, nmax_, krange, ios, nt, order, src_type

      real(8) ::    step, eft, duration,                           &
                    start, finish, lap, aux, aux1, aux2,              &
                    abstol,                &
                    tt, tp, tc, tmid,                     &
                    fmid,                                             &
                    norm_1, norm_2,                                   &
                    p0, pexc, pion,                             &
                    err1, err2,                                       &
                    rowsum,                                    &
                    kappa_w, omega_k,                               &
                    jacc, xx1, xx2

      complex(8) :: c1, c0, cnum, a0, nrg_, xt_, pt_


      real(8), allocatable, dimension(:) :: xs, xx, wx,    &
                                   vec_matup, eigval,                 &
                                   kk, omega, time_, norm_, norm_x

      real(8), allocatable, dimension(:) :: time_t, norm_t1, norm_t2,  &
                                            p0_t, pexc_t, pion_t

  
      complex(8), allocatable, dimension(:) :: wf0_, wfc0_,            &
                                               dwf0_, dwfc0_,          &
!                                              wf, wfc, wfc_,         &
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
                                               auxc1, auxc2,          &
                                               auxc3, auxc4,          &
                                               d_t, dd_t, d_w, dd_w,  &
                                               nrg_t, x_t, p_t

      integer, allocatable, dimension(:,:) :: map

      real(8), allocatable, dimension(:,:) :: lu, id, inv,            &
                                       Dref, Dglobal,                 &
                                       eigvec, basis, norm_t

      real(8), allocatable, dimension(:,:,:) :: Dloc_all          

      complex(8), allocatable, dimension(:,:) :: in_states,           &
                                                 out_states,          &
                                                 wf0, wfc0,           &
                                                 wf_in, wfc_in,       &
                                                 wf_out, wfc_out,     &
                                                 svec

      character(255) :: workdir, struct_dir, struct_dir_,             &
                        dyn_tag



      include 'param_dynamic'

      integer :: log_unit
      integer :: obs_unit
      integer :: init_unit
      integer :: force_unit

      c0 = (0.d0, 0.d0)
      c1 = (1.d0, 0.d0)
      
      dyn_tag = "dyn/"
      log_unit = 20
      obs_unit = 40

      call cpu_time(start)
      qho = 0 
      call get_command_argument(1, struct_dir)
!     write(*,*) struct_dir
      call read_problem_bin(trim(struct_dir)//"/problem.bin",        &
                             struct_dir_, nmax_, ns, np,             &
                             xx1, xx2, jacc, xx, wx)

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
       call init_src(src_type_, nchannel, omg_1, omg_0)
       call set_resolution_order(order_)
       call set_other_dyn_params(do_time_obs_, obs_stride_)


       duration = int((t_end - t_ini))
!      allocate(fvec(nt))
       tc = 0.d0
 

       allocate( auxc(nmax_), svec(nmax_,3) )
 
!      allocate(psi0(nmax_), phi0(nmax_))
!      allocate(psi_in(nmax_), psi_out(nmax_))
!      allocate(psi_inx(nmax_), psi_outx(nmax_))
!      allocate(phi_in(nmax_), phi_out(nmax_))
!      allocate(phi_inx(nmax_), phi_outx(nmax_))
!      allocate(phi_inc(nmax_))
!      allocate(psi_exact(nmax_))
!      allocate(psi_ex(nmax_), phi_ex(nmax_))
 
       allocate(src_mid(nmax_), src(nmax_))
       allocate(src_x(nmax_))

! --- initial condition ---
       wfc0_ = eigvec(:,1)
       wf0_ = matmul(transpose(eigvec),wfc0_)
       wfc0_ = eigvec(:,1)/wx/dsqrt(jacc)

       allocate(wfc0(nmax_, nchannel+2), wf0(nmax_, nchannel+2))
       allocate(wf_in(nmax_, nchannel+2), wfc_in(nmax_, nchannel+2))
       allocate(wf_out(nmax_, nchannel+2), wfc_out(nmax_, nchannel+2))
       allocate(omega(nchannel+2))
       allocate(norm_(nchannel+2))
       allocate(norm_x(nchannel+2))
       
       ! --- ψ ---
       do k=1,nchannel+2
          wfc0(:,k) = (kapp)**(1.d0/2) *  exp(-kapp*abs(xx))
!         wfc0(:,k) = eigvec(:,1)/wx/dsqrt(jacc)
       enddo
       j=1

!     do i=nmax_/2-5, nmax_/2+5
!        write(*,*) xx(i), eigvec(i,1), wfc0(i,1)
!     enddo
!     pause
      
       omega(1) = 0.d0 
       ! --- φ channels ---
       if (src_type.eq.3) then
          do k = 2, nchannel+2

             j = 2
              omega_k = omg_0 + (k-2)* dw

              omega(k) = omega_k
              kappa_w = varkap(kapp, omega_k)
          
              do i = 1, nmax_
                  wfc0(i,k) = (-ci*kapp**(1.5d0)/omega_k)            &
                    * sgn(xx(i)) *                                   &
                    (exp(-kapp*abs(xx(i))) - exp(-kappa_w*abs(xx(i))))
              enddo
          
          enddo
       end if

       write(*,*) omega
       pause
!     do i=nmax_/2-5, nmax_/2+5
!        write(*,*) xx(i), eigvec(i,1), wfc0(i,1)
!     enddo
!     pause

!     do k=1, nchannel+2
!        wfc0(:,k) = wfc0(:,k) * wx * dsqrt(jacc)
!     enddo

      ! --- rotate to eigenbasis ---
!     call zgemm('c','n', nmax_, nchannel+2, nmax_, c1, eigvec, nmax_, &
!                  wfc0, nmax_, c0, wf0, nmax_)
!     wf0 =  matmul(transpose(eigvec),wfc0)

      call dvr_to_eigen(nmax_, nchannel, jacc,                   &
                                           wx, eigvec, wfc0, wf0)





!     do i=nmax_/2-5, nmax_/2+5
!        write(*,*) xx(i), eigvec(i,1), wfc0(i,1)
!     enddo
!     write(*,*)
!     do i=1, 10
!        write(*,*) i, wf0_(i), wf0(i,1)
!     enddo

      write(*,*) "src_type is ", src_type
!     write(*,*) "norm", sum(conjg(wfc0_)* wfc0_ * wx*wx*jacc),      &
!                    sum((conjg(wfc0(:,1)) * wfc0(:,1)*wx*wx*jacc))
!     write(*,*) "kapp", kapp
!     pause


       tt = t_ini          ! start time
       wf_in = wf0         ! initial wavefunction

      call init_run(workdir, dyn_tag, extract_name(struct_dir))

      open(newunit=log_unit, file=trim(workdir)//"log.txt",            &
                                               status='replace')




!     call test_richardson_inhomogeneous(log_unit, nchannel,           &
!                                          nmax_, ns, np,      &
!                                          jacc,                       &
!                                          xs, xx, wx,                 &
!                                          map, Dref,                  &
!                                          dt0, t_end, t_ini,          &
!                                          eigval, eigvec,             &
!                                          psi0, phi0,                 &
!                                          src_type, omeg, order)

!    pause



       write(*,*) "Saving the dynamic parameters"
       call write_dynamic_bin(trim(workdir)//"dynamic.bin",            &
                        workdir, struct_dir_,                          &
                        f0, omega0, pfai,                              &
                        t_end, t_ini, nt, dt0,                         &
                        noc, ntau, src_type,                           &
                        nchannel, omg_1, omg_0, dw,                    &
                        order_)
 
       call write_dynamic_input(trim(workdir)//"param_dynamic.txt",    &
                        workdir, struct_dir_,                          &
                        f0, omega0, pfai,                              &
                        t_end, t_ini, nt, dt0,                         &
                        noc, ntau, src_type,                           &
                        nchannel, omg_1, omg_0, dw,                    &
                        order_)

 
       call write_problem_input(trim(workdir)//"param_structure.txt",  &
                                    struct_dir, nmax_, ns, np,         &
                                    xx1, xx2, jacc)
 

       write(*,*) "Saving the initial conditions"
       call write_wavefunction_bin(trim(workdir)//'init_state.bin',    &
                                    nmax_, nchannel, t_ini, omega, wf0)
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


      open(newunit=init_unit, file=trim(workdir)//"initial_conds.dat", &
                                               status='replace')
      do i=1,nmax_
         write(init_unit,*) i, xx(i),                                  &
                         real(wfc0(i,1)), imag(wfc0(i,2:))
      enddo
 
 
       open(newunit=force_unit, file=trim(workdir)//"force.dat",       &
                                                   status='replace')
 
       call plot_force(force_unit, t_end, t_ini, dt0)
       close(force_unit)
 
       open(newunit=obs_unit, file=trim(workdir)//"norm.dat",           &
                                                       status='replace')
       nobs = nt / obs_stride
       if (mod(nt, obs_stride) /= 0) nobs = nobs + 1
 
       allocate(time_(nobs), norm_t1(nobs), norm_t2(nobs))
       allocate(norm_t(nchannel+2,nobs))
       allocate(p0_t(nobs), pexc_t(nobs), pion_t(nobs))
       allocate(nrg_t(nobs), x_t(nobs), p_t(nobs))


       kobs = 0
       write(*,*) "Starting propagation"
       write(*,*) "nt          = ", nt
       write(*,*) "nobs        = ", nobs
 
!      pause
 
       do i = 1, nt
 
          ! --- propagation ---
 
!         call process_src_ingredients ( nmax_, ns, np,              &
!                                   jacc,                            &
!                                   xs, xx, wx, map, Dref,     &
!                                   dt0, tt, omega,                     &
!                                   eigval, eigvec, psi_in, svec,       &
!                                   src_type, order)
!
!         call build_source_quadrature ( nmax_, ns, np,             &
!                                               xs, xx, map, Dref,     &
!                                               dt0, tt, omega,        &
!                                               eigval, eigvec,        &
!                                               svec,                  &
!                                               src,                   &
!                                               order )
         
!        call split_operator(nmax_, dt0, tt, xx, eigval, eigvec,       &
!                                           psi_in, psi_out, order)

         call split_operator(nmax_, nchannel,                          &
                              dt0, tt, xx, eigval, eigvec,             &
                              wf_in, wf_out, order)


!        !--------------------------------------------
!        ! Add source contribution
!        !--------------------------------------------
!        phi_out = phi_out - ci * src
!
!        norm_1 = sqrt(sum(abs(psi_out)**2))
!        norm_2 = sqrt(sum(abs(phi_out)**2))


!        write(obs_unit,*) tt, norm_1, norm_2
!
         do k=1,nchannel+2
!           norm_(k)  = sum(abs(wf_out(:,k))**2 * (wx**2) * jac)
            norm_(k)  = sqrt(sum(abs(wf_out(:,k))**2))
         enddo
!        norm_(:) = sum( abs(vec_x(:,:))**2 *                         &
!           spread(wx**2 * jac, dim=2, ncopies=nchannel+2), dim=1 )

          tt = tt + dt0



         ! --- observables ---
         if (do_time_obs) then

            if (mod(i, obs_stride) == 0) then

               kobs = kobs + 1
            
               call compute_dyn_observables(nmax_, nchannel,           &
                                            wf_out,                    &
                                            xx, wx, jacc,              &
                                            eigval, eigvec,            &
!                                           norm_1, norm_2,            &
                                            norm_x,                    &
                                            p0, pexc, pion,      &
                                            nrg_, xt_, pt_)
            
               call append_dyn_obs_bin(trim(workdir)//"dyn_back.bin",  &
!                                  tt, norm_1, norm_2,           &
                                   tt, nchannel, norm_x,               &
                                   p0, pexc, pion,                     &
                                   nrg_, xt_, pt_ )

!              write(log_unit,'(f12.6,1x,7e20.10)') kobs, nobs,        &
               write(log_unit,*) kobs, nobs, tt,                       &
                                        norm_1, norm_2,                &
                                        p0, pexc, pion,                &
                                        nrg_, xt_, pt_

               ! --- STORE ---
               time_(kobs)    = tt
!              norm_t1(kobs)  = norm_1
!              norm_t2(kobs)  = norm_2
               norm_t(:,kobs) = norm_x
               p0_t(kobs)     = p0
               pexc_t(kobs)   = pexc
               pion_t(kobs)   = pion
               nrg_t(kobs)    = real(nrg_)
               x_t(kobs)      = real(xt_)
               p_t(kobs)      = real(pt_)
            
            end if


         end if
!
!         call exact_closed_duhamel(nmax_, omega,                       &
!                tt, t_ini, eigval, psi0, phi0, psi_ex, phi_ex, src_type)
!
!
!        write(*,*) nt, i, tt, phi_out(j), phi_ex(j)
         if (mod(i,100).eq.0) then
!           call write_wavefunction_bin(trim(workdir)//'wavfun.bin',   &
!                                    nmax_, psi_out, phi_out, omg0, tt)
            call write_wavefunction_bin(trim(workdir)//'wavfun.bin',   &
                                    nmax_, nchannel, tt, omega, wf_out)
         end if
 
 
          wf_in = wf_out

       enddo
 
       call write_observables_bin(trim(workdir)//"dyn_obs.bin",        &
                 nchannel, nobs, time_,                                &
!                norm_t1, norm_t2,                                     &
                 norm_t,                                               &
                 p0_t, pexc_t, pion_t,                                 &
                 nrg_t, x_t, p_t)
 
      call write_wavefunction_bin(trim(workdir)//'wavfun.bin',         &
                                  nmax_, nchannel, tt, omega, wf_out)
 
       close(obs_unit)
       write(log_unit,*) workdir
       close(log_unit)


      end program
