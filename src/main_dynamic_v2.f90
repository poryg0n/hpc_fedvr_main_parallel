      program fedvr_dynamic_build
      use omp_lib
      use constants, only : ppi, ci, omg_ => omega_qho
      use dynamic_parameters
      use math_util
      use space_time_ops
      use util
      use propagation
!     use conv_tests
      use observables
      use io_module
      implicit none
      integer :: np, ns, nmax_, nobs, kobs, qho,                      &
                  lwork, info, store_val,                        &
                  i, j, k, ij, p, q, r, s, m, n, w,                   &
                  jn, ijk, ios
!                 nmax, nmax_, krange, ios, nt, order, src_type

      real(8) :: step, eft, duration,                           &
                 start, finish, lap, aux, aux1, aux2,              &
                 abstol, qq,                 &
                 tt, tp, tc, tmid,                     &
                 p0_, pexc_, pion_,                                &
                 fmid,                                             &
                 rowsum,                                    &
                 kappa_w, omega_k, norm_ref,                       &
                 jacc, xx1, xx2, src_time, split_time, eps

      complex(8) :: c1, c0, cnum, a0, nrg_, xt_, pt_


      real(8), allocatable, dimension(:) :: xs, xx, wx,    &
                                   vec_matup, eigval,                 &
                                   kk, omega, time_,                  & 
                                   norm_,                             &
                                   p0_t, pexc_t, pion_t,              &
                                   norm_ref1, norm_ref2

  
      complex(8), allocatable, dimension(:) :: wf0_, wfc0_,            &
                                               dwf0_, dwfc0_,          &
!                                              wf, wfc, wfc_,         &
                                               psi_exact,             &
                                               psi0, psi, psi_x,      &
                                               psi_in, psi_out,       &
                                               src_out, src_mid,      &
!                                              src, src_x,            &
                                               dwfc, dwfc2, dwfc3,    & 
                                               init1, init2,          &
                                               phi0, phic0,           &
                                               psi_test, phi_test,    &
                                               psi_ex, phi_ex,        &
                                               auxc, auxck1, auxck2,  &
                                               auxc1, auxc2,          &
                                               auxc3, auxc4,          &
                                               x_t, p_t, nrg_t,       &
                                               d_t, dd_t, d_w, dd_w

      integer, allocatable, dimension(:,:) :: map

      real(8), allocatable, dimension(:,:) :: lu, id, inv,            &
                                       Dref, Dglobal,                 &
                                       eigvec, basis,                 &
                                       norm_t

      real(8), allocatable, dimension(:,:,:) :: Dloc_all          

      complex(8), allocatable, dimension(:,:) :: in_states,           &
                                                 out_states,          &
                                                 wf0, wfc0,           &
                                                 wf_in, wfc_in,       &
                                                 wf_out, wfc_out,     &
                                                 svec, src

      character(255) :: workdir, struct_dir, struct_dir_,             &
                        dyn_tag

      logical :: done



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
      eps = 1.e-6
      call get_command_argument(1, struct_dir)
      write(log_unit,*) "Reading the structure parameters"
      call read_problem_bin(trim(struct_dir)//"/problem.bin",        &
                             struct_dir_, nmax_, ns, np,             &
                             xx1, xx2, qq, jacc,                     &
                             xx, wx)

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
       call init_src(src_type_, nchan, omg_1, omg_0)
       call set_resolution_order(order_)
       call set_other_dyn_params(do_time_obs_, obs_stride_)


       duration = int((t_end - t_ini))
!      allocate(fvec(nt))
       tc = 0.d0
 

       allocate( auxc(nmax_), svec(nmax_,3) )
 
!      allocate(psi0(nmax_), phi0(nmax_))
       allocate(psi_in(nmax_), psi_out(nmax_))
!      allocate(psi_inx(nmax_), psi_outx(nmax_))
!      allocate(phi_in(nmax_), phi_out(nmax_))
!      allocate(phi_inx(nmax_), phi_outx(nmax_))
!      allocate(phi_inc(nmax_))
!      allocate(psi_exact(nmax_))
!      allocate(psi_ex(nmax_), phi_ex(nmax_))
 

! --- initial condition ---
       wfc0_ = eigvec(:,1)
       wf0_ = matmul(transpose(eigvec),wfc0_)
       wfc0_ = eigvec(:,1)/wx/dsqrt(jacc)


       allocate(wfc0(nmax_, nchan+1), wf0(nmax_, nchan+1))
       allocate(wf_in(nmax_, nchan+1), wfc_in(nmax_, nchan+1))
       allocate(wf_out(nmax_, nchan+1), wfc_out(nmax_, nchan+1))
       allocate(omega(nchan+1))
       allocate(norm_(nchan+1))
       allocate(norm_ref1(nchan+1))
       allocate(norm_ref2(nchan+1))
       
       ! --- ψ ---
!       do k=1,nchan+2
!          wfc0(:,k) = (kapp)**(1.d0/2) *  exp(-kapp*abs(xx))
!!         wfc0(:,k) = eigvec(:,1)/wx/dsqrt(jacc)
!       enddo


!      do i=1, nmax_
!         wfc0(i,1) = kapp**(.5d0) * exp(-kapp*abs(xx(i)))
!      enddo

       j=1
       ! --- φ channels ---
       if (src_type.eq.3) then
          !$omp parallel do default(shared) &
          !$omp private(k,i,omega_k,kappa_w)
          do k = 2, nchan+1

             j = 2
             omega_k = omega_min + (k-2)* dw

             omega(k) = omega_k
             kappa_w = varkap(kapp, omega_k)
          
             do i = 1, nmax_
                 wfc0(i,k) = (-ci*kapp**(1.5d0)/omega_k)              &
                 * sgn(xx(i)) *                                       &
                 ( exp(-kapp*abs(xx(i))) - exp(-kappa_w*abs(xx(i))) )
             enddo
!            write(*,*) "in here?" , wfc0(1300,k), size(wfc0,1)
          enddo
          !$omp end parallel do
       end if


       ! --- ψ ---
       omega(1) = 0.d0 
!      wfc0(:,1) = eigvec(:,1)/wx/dsqrt(jacc)
       wfc0(:,1) = kapp**(.5d0) * exp(-kapp*abs(xx))


!      write(*,*) omega
!      do i=nmax_/2-5, nmax_/2+5
!         write(*,'(5E20.10)') xx(i), wfc0(i,1), wfc0(i,2)
!      enddo


       ! --- rotate to eigenbasis ---
!      call zgemm('c','n', nmax_, nchan+2, nmax_, c1, eigvec, nmax_, &
!                   wfc0, nmax_, c0, wf0, nmax_)
!      wf0 =  matmul(transpose(eigvec),wfc0)

       call dvr_to_eigen(nmax_, nchan, jacc,                   &
                                  wx, eigvec, wfc0, wf0)


!     do i=nmax_/2-5, nmax_/2+5
!        write(*,'(3E20.10)') xx(i), eigvec(i,1)/wx(i)/dsqrt(jacc),    &
!                                                    real(wfc0(i,1))
!     enddo
!     write(*,*)
!     do i=1, 10
!        write(*,'(I8,1x,2E20.10)') i, real(wf0_(i)), real(wf0(i,1))
!     enddo

!     write(*,*) "src_type is ", src_type
!     write(*,*) "norm", sum(conjg(wfc0_)* wfc0_ * wx*wx*jacc),      &
!                    sum((conjg(wfc0(:,1)) * wfc0(:,1)*wx*wx*jacc))
!     write(*,*) "kapp", kapp


       write(*,*)
       do k=1,nchan+1
          norm_ref2(k) = dsqrt(sum(abs((wfc0(:,k)*wx*dsqrt(jacc))**2 )))
          norm_ref1(k) = dsqrt(sum(abs(wf0(:,k)**2)))
          write(*,'(2E20.10)') norm_ref1(k), norm_ref2(k)
       enddo
       
!      write(*,*)
!      do k=1,nchan+1
!         write(*,'(2E20.10)') norm_ref2(k)-norm_ref1(k)
!      enddo
!      pause

       tt = t_ini          ! start time
       wf_in = wf0         ! initial wavefunction

       call init_run(workdir, dyn_tag, extract_name(struct_dir))

       open(newunit=log_unit, file=trim(workdir)//"log.txt",           &
                                                status='replace')


       write(*,*) "Saving the dynamic parameters"
       call write_dynamic_bin(trim(workdir)//"dynamic.bin",            &
                        workdir, struct_dir_,                          &
                        f0, omega0, pfai,                              &
                        t_end, t_ini, nt, dt0,                         &
                        noc, ntau, src_type,                           &
                        nchan, omg_1, omg_0, dw,                    &
                        order_)
 
       call write_dynamic_input(trim(workdir)//"param_dynamic.txt",    &
                        workdir, struct_dir_,                          &
                        f0, omega0, pfai,                              &
                        t_end, t_ini, nt, dt0,                         &
                        noc, ntau, src_type,                           &
                        nchan, omg_1, omg_0, dw,                    &
                        order_)

 
       call write_problem_input(trim(workdir)//"param_structure.txt",  &
                                    struct_dir, nmax_, ns, np,         &
                                    xx1, xx2, qq, jacc)
 

       write(*,*) "Saving the initial conditions"
       call write_wavefunction_bin(trim(workdir)//'init_state.bin',    &
                                    nmax_, nchan, t_ini, omega, wf0)

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


! *** TO FIX LATER
!      if (do_conv_test) then
!        call test_richardson_inhomogeneous(log_unit, nchan,        &
!                                             nmax_, ns, np,           &
!                                             jacc,                    &
!                                             xs, xx, wx,              &
!                                             map, Dref,               &
!                                             dt0, t_end, t_ini,       &
!                                             eigval, eigvec,          &
!                                             wf0,                     &
!                                             src_type, omeg, order)

!        write(*,*) "richardson test done"
!     end if 
   

      open(newunit=init_unit, file=trim(workdir)//"initial_conds.dat", &
                                               status='replace')
      do i=1,nmax_
!        write(init_unit,'()') i, xx(i),                               &
         write(init_unit,'(I8,1X,ES24.15,*(1X,ES24.15,1X,ES24.15))') &
                         i, xx(i), real(wfc0(i,1)), imag(wfc0(i,2:))
      enddo
 
 
       open(newunit=force_unit, file=trim(workdir)//"force.dat",       &
                                                   status='replace')
       call plot_force(force_unit, t_end, t_ini, dt0)
       close(force_unit)
 
       open(newunit=obs_unit, file=trim(workdir)//"norm.dat",          &
                                                  status='replace')

       nobs = nt / obs_stride
       if (mod(nt, obs_stride) /= 0) nobs = nobs + 1
 
       allocate(time_(nobs))
       allocate(norm_t(nchan+1,nobs))
       allocate(p0_t(nobs), pexc_t(nobs), pion_t(nobs))
       allocate(nrg_t(nobs), x_t(nobs), p_t(nobs))
       allocate(src(nmax_, nchan))


       kobs = 0
       write(*,*) "Starting propagation"
       write(*,*) "nt          = ", nt
       write(*,*) "nobs        = ", nobs
 
!      pause
 
       write(obs_unit,'(2ES20.10,*(1X,ES24.15,1X,ES24.15))') 0.d0, &
                  omega(:)
       write(*,'(2ES20.10,*(1X,ES24.15,1X,ES24.15))') 0.d0, &
                  omega(:)




       src_time = 0.d0
       split_time = 0.d0


       do i = 1, nt
  
          ! --- propagation ---

!         call cpu_time(start)

            psi_in = wf_in(:,1)

            call process_src_ingredients ( nmax_, ns, np,              &
                                      jacc,                            &
                                      xs, xx, wx, map, Dref,     &
                                      dt0, tt,                         &
                                      eigval, eigvec, psi_in, svec,    &
                                      src_type, order)

!           write(*,*) svec(:,1)



         call extend_source_quadr_build ( nchan,                     &
                                            nmax_, ns, np,         &
                                            xs, xx, map, Dref,     &
                                            dt0, tt, omega,        &
                                            eigval, eigvec,        &
                                            svec,                  &
                                            src,                   &
                                            order )

         call extend_split_operator(nmax_, nchan,                      &
                              dt0, tt, xx, eigval, eigvec,             &
                              wf_in, wf_out, order)


!         call cpu_time(finish)


       !--------------------------------------------
       ! Add source contribution
       !--------------------------------------------

       do j=2, nchan+1
          wf_out(:,j) = wf_out(:,j) - ci * src(:,j-1)
!         write(*,*) src(:,j-1)
       enddo


!        norm_1 = sqrt(sum(abs(psi_out)**2))
!        norm_2 = sqrt(sum(abs(phi_out)**2))

         do k=1,nchan+1
            norm_(k) = sqrt(sum(abs(wf_out(:,k))**2))
         enddo
         norm_ref = norm_(2)-norm_ref1(2)


!        norm_(:) = sum( abs(vec_x(:,:))**2 *                         &
!           spread(wx**2 * jac, dim=2, ncopies=nchan+2), dim=1 )
         write(obs_unit,'(2ES20.10,*(1X,ES24.15,1X,ES24.15))') tt, &
                  norm_(:)
         write(*,'(2ES20.10,*(1X,ES20.10))') tt, &
                  norm_ref, norm_(:)


         tt = tt + dt0


         ! --- observables ---
         if (do_time_obs) then

            if (mod(i, obs_stride) == 0) then

               kobs = kobs + 1
            
               call compute_dyn_observables(nmax_, nchan,           &
                                            wf_out,                    &
                                            xx, wx, jacc,              &
                                            eigval, eigvec,            &
                                            norm_,                    &
                                            p0_, pexc_, pion_,      &
                                            xt_, pt_, nrg_)
            
               call append_dyn_obs_bin(trim(workdir)//"dyn_back.bin",  &
                                   tt, nchan, norm_,               &
                                   p0_, pexc_, pion_,               &
                                   xt_, pt_, nrg_ )

               write(log_unit,'(f12.6,1x,*(e20.10))') kobs, nobs, tt,  &
                                        norm_(1:4),              &
                                        p0_, pexc_, pion_,          &
                                        xt_, pt_, nrg_

               ! --- STORE ---
               time_(kobs)    = tt
               norm_t(:,kobs) = norm_
               p0_t(kobs)     = p0_
               pexc_t(kobs)   = pexc_
               pion_t(kobs)   = pion_
               nrg_t(kobs)    = real(nrg_)
               x_t(kobs)      = real(xt_)
               p_t(kobs)      = real(pt_)
            
            end if


         end if
!
!         call exact_closed_duhamel(nmax_, omega,                       &
!                tt, t_ini, eigval, psi0, phi0, psi_ex, phi_ex, src_type)
!
         if (mod(i,100).eq.0) then
            call write_wavefunction_bin(trim(workdir)//'wf_back.bin',  &
                                    nmax_, nchan, tt, omega, wf_out)
         end if
 
 
           wf_in = wf_out
!          write(*,*) nt, i, tt

           src_time   = src_time   + (lap-start)
           split_time = split_time + (finish-lap)



       enddo


 
       call write_observables_bin(trim(workdir)//"dyn_obs.bin",        &
                 nchan, nobs, time_,                                &
                 norm_t,                                               &
                 p0_t, pexc_t, pion_t,                                 &
                 x_t, p_t, nrg_t)
 
       call write_wavefunction_bin(trim(workdir)//'wavfun.bin',        &
                                  nmax_, nchan, tt, omega, wf_out)

       write(*,*) split_time
       write(*,*) src_time
 
       close(obs_unit)
       write(log_unit,*) workdir
       close(log_unit)


      end program
