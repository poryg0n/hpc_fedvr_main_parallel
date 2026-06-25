      module observables
      use constants, only : ci
      use constants, only : ppi
      use constants, only : kapp
      use exploit_parameters
      use math_util
      use space_time_ops
      implicit none
      integer :: mode_k = 0                ! 1 for symmetric
      contains

!     subroutine compute_observables(...)
!        call spectral_obs(...)
!        call spatial_obs(...)
!     end subroutine

      subroutine compute_dens_probab(n, nch, jacc, wx, eigvec, wf, rho)

        implicit none
        integer, intent(in) :: n, nch
        real(8), intent(in) :: jacc
        real(8), intent(in) :: wx(n)
        real(8), intent(in) :: eigvec(n,n)
        complex(8), intent(in) :: wf(n,nch+1)
        real(8), allocatable, intent(out) :: rho(:,:)
!       real(8), intent(out) :: norm

        complex(8) :: wfc(n,nch+1)
        integer :: i, k

        allocate(rho(n,nch+1))
        call eigen_to_dvr(n, nch, jacc, wx, eigvec, wf, wfc)
!       wfc = matmul(eigvec, wf)
!       wfc = wfc / wx / dsqrt(jacc)

        do k=1,nch+1
           rho(:,k) = conjg(wfc(:,k))* wfc(:,k) * wx*wx*jacc
        enddo

      end subroutine


      subroutine compute_wf(n, nch,                                   &
                                       jacc, omega, wx, eigvec,       &
                                       wf, re_wf, im_wf, rho, arg)
        implicit none
      
        integer, intent(in) :: n, nch
        real(8), intent(in) :: jacc
        real(8), intent(in) :: wx(n), omega(nch+1)
        real(8), intent(in) :: eigvec(n,n)
        complex(8), intent(in) :: wf(n, nch+1)   ! multi-channel
      
        real(8), intent(out) :: re_wf(n, nch+1)
        real(8), intent(out) :: im_wf(n, nch+1)
        real(8), intent(out) :: rho(n,nch+1), arg(n,nch+1)
      
        complex(8) :: wfc(n, nch+1)
        integer :: i, w
      
      

        call eigen_to_dvr(n, nch, jacc, wx, eigvec, wf, wfc)

      
        ! --- observables
        re_wf = real(wfc)
        im_wf = aimag(wfc)
        rho   = re_wf**2 + im_wf**2
        arg   = atan2(im_wf, re_wf)

!       do w=1,nch
!          do i = 1, n
!             re_wf(i,w) = real(wfc(i,w))
!             im_wf(i,w) = aimag(wfc(i,w))
!             rho(i,w)   = re_wf(i,w)**2 + im_wf(i,w)**2
!             arg(i,w)   = atan2(im_wf(i,w), re_wf(i,w))
!          enddo
!       enddo
      
      end subroutine

      

      subroutine compute_dyn_observables(nmax, nch, vec,      &
                                        xx, wx, jac,               &
                                        eigval, eigvec,            &
                                        norm_t,                     &
                                        p0_t, pexc_t, pion_t,       &
                                        dipole_t, momentum_t, energy_t)
      
        implicit none
        integer, intent(in) :: nmax, nch
        real(8), intent(in) :: jac
        real(8), intent(in) :: xx(nmax), wx(nmax)
        complex(8), intent(in) :: vec(nmax,nch+1)
        real(8), intent(in) :: eigval(nmax)
        real(8), intent(in) :: eigvec(nmax,nmax)
      
        real(8), intent(out) :: p0_t, pexc_t, pion_t
        real(8), intent(out) :: norm_t(nch+1)
        complex(8), intent(out) :: dipole_t, momentum_t, energy_t
      
        complex(8) :: vec_x(nmax, nch+1)
        complex(8) :: p_psi(nmax)


        complex :: c0, c1
        integer :: i, k

        c0 = (0.d0, 0.d0)
        c1 = (1.d0, 0.d0)
      

        ! --- transform to DVR ---
!       psi_x = matmul(eigvec, psi)
!       phi_x = matmul(eigvec, phi)
!       psi_x = psi_x/wx/dsqrt(jac)
!       phi_x = phi_x/wx/dsqrt(jac)
!       call zgemm('n','n', nmax, nch+2, nmax, c1, eigvec,        &
!                  nmax, vec, nmax, c0, vec_x, nmax)
!       vec_x = matmul(eigvec, vec)
        call eigen_to_dvr(nmax, nch, jac, wx, eigvec, vec, vec_x)

        
        ! --- norm ---
        do k=1,nch+1
           norm_t(k)  = sum(abs(vec_x(:,k))**2 * (wx**2) * jac)
        enddo
        
        ! --- populations ---
        p0_t   = abs(vec(1,1))**2
        pexc_t = sum(abs(vec(2:nmax,1))**2)
        pion_t = 1.d0 - norm_t(1)
        
        ! --- dipole ---
        dipole_t = sum(conjg(vec_x(:,1)) * xx * vec_x(:,1) *           &
                                                   (wx**2) * jac)
        
        ! --- energy (field-free) ---
        energy_t = sum(abs(vec(:,1))**2 * eigval)
        
        ! --- momentum (CLEAN) ---
        call apply_momentum_operator(nmax, eigvec, xx, wx, jac,        &
                                                vec(:,1), p_psi, 0)
        momentum_t = sum(conjg(vec(:,1)) * p_psi)
      
      end subroutine




      subroutine compute_pemd_zrp(nmax, nch, krange, t_end,            &
                           xx, wx, jacc,                               &
                           eigvec, eigval,                             &
                           wf0_0, wf,                        &
                           k_max, kk, p_ion, p0,             &  
                           a0, ak,                           &
                           b0wT, bkwT)              


      
        implicit none
        integer, intent(in) :: nmax, krange, nch
        real(8), intent(in) :: jacc
        real(8), intent(in) :: xx(nmax), wx(nmax)
        real(8), intent(in) :: eigval(nmax)
        real(8), intent(in) :: eigvec(nmax,nmax)
        real(8), intent(in) :: t_end, k_max
        complex(8), intent(in) :: wf0_0(nmax)
        complex(8), intent(in) :: wf(nmax,nch+1)

        real(8), intent(out) :: p_ion, p0
        complex(8), intent(out) :: a0
        real(8), allocatable, intent(out) :: kk(:)
        complex(8), allocatable, intent(out) :: ak(:), b0wT(:)
        complex(8), allocatable, intent(out) :: bkwT(:,:)

      
        ! locals
        integer :: j, k, l, ij
        integer :: p, n_cont, ksteps_
        real(8) :: aux1, aux2, dk
        real(8) :: E0
        complex(8) :: auxc_
        complex(8) :: wfc_k(nmax)
        complex(8) :: wfc0_0(nmax)
        complex(8) :: wfc(nmax,nch+1)

        real(8), allocatable :: Ek(:)

        complex(8) :: auxcw(nmax,nch+1)
        complex(8), parameter :: ci = ( 0.d0, 1.d0 )
!       real(8), parameter :: ppi = 3.141592653589793d0
        real(8), parameter :: ppi = 4.d0*datan(1.d0)



        call eigen_to_dvr(nmax, nch, jacc, wx, eigvec, wf, wfc)

        ksteps_ = krange-1

        allocate(kk(krange))
        allocate(ak(krange))
        allocate(b0wT(nch))
        allocate(bkwT(krange,nch))


        ak = (0.d0, 0.d0)
!       dk = k_max/(krange/2-1)

        kk = 0.d0
        do j=1,krange

           if (j.le.(ksteps_/2)) then
              ij = krange-(j-1)
              kk(j)  = -k_max + (j-1)*dk0
              kk(ij) = -kk(j)

           end if

!          call build_wfc_k(xx, kk(j), kapp, mode_k, wfc_k)

           wfc_k = exp(ci*kk(j)*xx) +                                 &
                   (ci*kapp/(-abs(kk(j)) - ci*kapp)) *                &
                   exp(-ci*abs(kk(j)*xx))

           auxcw(:,1)  = conjg(wfc_k) * wfc(:,1) * wx*wx*jacc
           ak(j) = sum(auxcw(:,1))
           
           do l=1,nch
              auxcw(:,l+1)  = conjg(wfc_k) * wfc(:,l+1) * wx*wx*jacc
              bkwT(j,l) = sum(auxcw(:,l+1))
           enddo

        end do


        Ek  = 0.5d0 * kk**2
        E0  = - 0.5d0 * kapp**2

        ak   = exp(ci*Ek*t_end) * ak
        do k=1,nch
           bkwT(:,k) = exp(ci*Ek*t_end) * bkwT(:,k)
        enddo
      
        a0 = (0.d0, 0.d0)
        wfc0_0 = matmul(eigvec, wf0_0)
        wfc0_0 = wfc0_0 / wx /dsqrt(jacc)
        auxcw(:,1) = conjg(wfc0_0) * wfc(:,1) * wx*wx*jacc
        a0 = sum(auxcw(:,1))
        a0 = exp(ci*E0*t_end)*a0


        b0wT = (0.d0, 0.d0)
        do k=1,nch
           auxcw(:,k+1) = conjg(wfc0_0) * wfc(:,k+1) * wx*wx*jacc
           b0wT(k) = sum(auxcw(:,k+1))
           b0wT(k) = exp(ci*E0*t_end)*b0wT(k)
        enddo

        ! --- probabilities ---
        p_ion = 0.d0
        auxcw(:,1) = abs(ak)**2
        call integr_over_krange(ksteps_, kk, auxcw(:,1), auxc_)
        p_ion = real(auxc_) / (2.d0*ppi)

!       do k=1, krange
!          auxcw(:,k)= abs(b0wT(k))**2
!          call integr_over_krange(ksteps_, kk, auxcw(:,k), auxc_)
!       enddo

        p0 = abs(a0)**2

      
      end subroutine




       subroutine compute_phi_elems(workdir,                          &
                            nmax, krange, nch, t_end,            &
                            xx, wx, jacc,                             &
                            eigvec, eigval,                           &
                            wf0_0, wf,                                &
                            omega, k_max, kk,                         &
                            a0, b0wT,                                 &
                            ak, bkwT,                                 &
                            b0w, bkw)
       
         implicit none
         integer, intent(in) :: nmax, krange, nch
         real(8), intent(in) :: jacc
         real(8), intent(in) :: xx(nmax), wx(nmax)
         real(8), intent(in) :: eigval(nmax)
         real(8), intent(in) :: eigvec(nmax,nmax)
         real(8), intent(in) :: t_end, k_max
         real(8), intent(in) :: omega(nch+1)
         complex(8), intent(in) :: a0
         complex(8), intent(in) :: b0wT(nch)

         real(8), intent(in) :: kk(krange)
         complex(8), intent(in) :: wf0_0(nmax)
         complex(8), intent(in) :: wf(nmax,nch+1)
         complex(8), intent(in) :: ak(krange)
         complex(8), intent(in) :: bkwT(krange,nch)

 
         character(255) :: workdir
 
 
         complex(8), allocatable, intent(out) :: b0w(:)
         complex(8), allocatable, intent(out) :: bkw(:, :)

        ! locals
        integer :: j, k, l, ij, w, o
        integer :: p, n_cont, ksteps_
        real(8) :: E0, Ek_, Ekp_
        real(8) :: Ek(krange)
        real(8) :: aux1, aux2, dk, delta_kk
        complex(8) :: factor, denom
        complex(8) :: wfc_k(nmax), wfc_k_(nmax)
        complex(8) :: dwfc_0(nmax), dwfc_k(nmax)
        complex(8) :: pwfc_0(nmax), pwfc_k(nmax)
        complex(8) :: wfc0_0(nmax)
        complex(8) :: wfc(nmax,nch+1)
        complex(8) :: auxc(nmax)
        complex(8) :: vec_1(krange,nch)
        complex(8) :: vec_k(krange,nch)
        complex(8) :: vec_0(nch)

        complex(8), allocatable :: vec_2(:,:)
        complex(8) :: vec_sum

        complex(8), parameter :: ci = ( 0.d0, 1.d0 )
!       real(8), parameter :: ppi = 3.141592653589793d0
        real(8), parameter :: ppi = 4.d0*datan(1.d0)

        complex(8) :: dk0_, dkk_
        complex(8) :: pk0_, p0k_, pkk_
        complex(8) :: auxc_1, auxc_2, auxc_3
        complex(8) :: pk0(krange), p0k(krange), pkk(krange)
        complex(8) :: pkkk(krange), dkk(krange)

       
        real(8) :: eta = 1.d-6


        integer :: unit_pk0, unit_pkk, unit_pkl, unit_vec, unit_b0kw

        open(newunit=unit_pk0, file=trim(workdir)//"/pk0.dat",         &
                                                      status="replace")
        open(newunit=unit_pkk, file=trim(workdir)//"/pkk.dat",         &
                                                      status="replace")
        open(newunit=unit_pkl, file=trim(workdir)//"/pkl.dat",         &
                                                    status="replace")
        open(newunit=unit_vec, file=trim(workdir)//"/vec_01k.dat",     &
                                                    status="replace")
        open(newunit=unit_b0kw, file=trim(workdir)//"/b0wk.dat",       &
                                                    status="replace")

!     
!
        ksteps_=krange-1

        call apply_momentum_operator(nmax, eigvec, xx, wx,        &
                          jacc, wf0_0, auxc, 0)

!       call eigen_to_dvr(nmax, jacc, wx, eigvec, auxc, pwfc_0)
        pwfc_0 = matmul(eigvec, auxc)
        pwfc_0 = pwfc_0 / wx / dsqrt(jacc)

!       call eigen_to_dvr(nmax, jacc, wx, eigvec, wf0, wfc0)
        wfc0_0 = matmul(eigvec, wf0_0)
        wfc0_0 = wfc0_0 / wx / dsqrt(jacc)


        E0  = - 0.5d0 * kapp**2

        
        allocate(vec_2(krange,nch))

        do j=1,krange
!          call build_wfc_k(xx, kk(j), kapp, mode_k, wfc_k)

           wfc_k = exp(ci*kk(j)*xx) +                                 &
                   (ci*kapp/(-abs(kk(j)) - ci*kapp)) *                &
                   exp(-ci*abs(kk(j)*xx))

           Ek_  = 0.5d0 * kk(j)**2

           auxc_1 = 0.d0
           auxc_2 = 0.d0
           auxc_3 = 0.d0

           do l=1,krange
!             if (j == l) cycle
!             call build_wfc_k(xx, kk(l), kapp, mode_k, wfc_k_)

              Ekp_ = 0.5d0 * kk(l)**2

              wfc_k_ = exp(ci*kk(l)*xx) +                              &
                      ( ci*kapp/(-abs(kk(l)) - ci*kapp) ) *            &
                      exp(-ci*abs(kk(l)*xx))

              call differentiate(xx, wfc_k_, dwfc_k)
              pwfc_k = -ci*dwfc_k


              auxc = conjg(wfc_k) * pwfc_k * wx*wx*jacc
              pkk(l) = sum(auxc)

              auxc = conjg(wfc_k) * xx * wfc_k_ * wx*wx*jacc
              dkk_ = sum(auxc)
              dkk(l) = -ci* ( Ek_ - Ekp_ ) * dkk_

              ! *** analytical formula
              if(j.eq.l) then
                 delta_kk = 1.d0/dk0
              else
                 delta_kk = 0.d0
              end if

              denom  = ( kk(j)**2 - kk(l)**2 + ci*eta )
              factor = ( kk(l)*abs(kk(j)) / (abs(kk(j)) -ci*kapp)      &
                      -  kk(j)*abs(kk(l)) / (abs(kk(l)) +ci*kapp) )


              pkkk(l) = 2.d0*ppi*kk(j) * delta_kk                      &
                                   - 2.d0 * kapp * factor/denom

              do w=1,nch       
                 vec_2(l,w) = pkk(l)*ak(l)/(Ek_+omega(w+1)-Ekp_+ci*eta)
                 vec_2(l,w) = exp(ci*(Ek_+omega(w+1)-Ekp_)*t_end)*vec_2(l,w)
              enddo


           enddo

           write(unit_pkl,'(7E20.10)') kk(j),                         &
                               real(auxc_1), imag(auxc_1),             &
                               real(auxc_2), imag(auxc_2),             &
                               real(auxc_3), imag(auxc_3)


           !$omp parallel do default(shared) &
           !$omp private(w)
           do w=1,nch
              call integr_over_krange(ksteps_, kk, vec_2(:,w), vec_k(j,w))
              vec_k(j,w) = vec_k(j,w) / (2.d0*ppi)
           enddo
           !$omp end parallel do
           !$omp barrier


           auxc = conjg(wfc_k) * pwfc_0 * wx*wx*jacc
           pk0(j) = sum(auxc)

           call differentiate(xx, wfc_k, dwfc_k)
           pwfc_k = -ci*dwfc_k


           auxc = conjg(wfc0_0) * pwfc_k * wx*wx*jacc
           p0k(j) = sum(auxc)


           auxc = conjg(wfc_k) * xx * wfc0_0 * wx*wx*jacc
           dk0_ = sum(auxc)
           dk0_ = -ci* ( Ek_ - E0 ) * dk0_

           denom = ( kapp**2 + kk(j)**2 )

           pk0_ = 2.d0 * kk(j) * kapp**(3.d0/2) / denom
           p0k_ = pk0_


           write(unit_pk0,'(5E20.10)') kk(j),                         &
                                     real(pk0(j)), real(p0k(j)),      &
                                     real(dk0_), real(pk0_)

           do w=1,nch
              vec_1(j,w) = p0k(j) * ak(j) / ( E0+omega(w+1)-Ek_+ci*eta )
              vec_1(j,w) = exp(ci*( E0+omega(w+1)-Ek_)*t_end )*vec_1(j,w)
           enddo

           write(*,*) j
        enddo

        deallocate(vec_2)
        allocate(bkw(krange,nch))
        allocate(b0w(nch))


        Ek = 0.5d0 * kk**2

        !$omp parallel do default(shared) &
        !$omp private(w)
        do w=1,nch
           call integr_over_krange(ksteps_, kk, vec_1(:,w), vec_0(w))

           vec_0(w) = vec_0(w) / (2.d0*ppi)

           vec_1(:,w) = pk0 * a0 / ( Ek+omega(w+1)-E0 )
           vec_1(:,w) = exp(ci*( Ek+omega(w+1)-E0 )*t_end) * vec_1(:,w)


           b0w(w) = b0wT(w) + vec_0(w)
           bkw(:,w) = bkwT(:,w) + vec_1(:,w) +  vec_k(:,w)



!          vec_sum = sum(vec_1)
!          write(*,*) "The sum for vec_1 is ", vec_sum

        enddo
        !$omp end parallel do
        !$omp barrier


        o=1

        do j=1,krange
           write(unit_vec,'(7E20.10)') kk(j), real(vec_0(o)),       &
                                              imag(vec_0(o)),       &
                                              real(vec_1(j,o)),     &
                                              imag(vec_1(j,o)),     &
                                              real(vec_k(j,o)),     &
                                              imag(vec_k(j,o))
           write(unit_b0kw,'(5E20.10)') kk(j), real(b0w(o)),        &
                                               imag(b0w(o)),        &
                                               real(bkw(j,o)),      &
                                               imag(bkw(j,o))
        enddo

        close(unit_pk0)
        close(unit_pkk)
        close(unit_pkl)
        close(unit_vec)
        close(unit_b0kw)




       end subroutine




      subroutine compute_Qwc(nt, time, x_t, wsteps, omg,          &
                                               hhg_1, hhg_2)
      
        implicit none
        integer, intent(in) :: nt
        integer, intent(in) :: wsteps
        real(8), intent(in) :: time(nt)
        complex(8), intent(in) :: x_t(nt)
        real(8), allocatable, intent(out) :: omg(:)
        real(8), allocatable, intent(out) :: hhg_1(:)
        real(8), allocatable, intent(out) :: hhg_2(:)
      
        complex(8) :: auxc(nt), integral, amp
        complex(8), parameter :: ci = (0.d0, 1.d0)
        integer :: nw, i, it
        real :: w

        nw = wsteps+1
        allocate(omg(nw), hhg_1(nw), hhg_2(nw))

        do i=1,nw
           omg(i) = wmin + (i-1)*dw0
!          write(*,*) i, omg(i), wmin
        enddo

      
        do i = 1, nw
      
           ! --- integral term ---
           auxc = exp(ci*omg(i)*time) * x_t
           call composite_simpson_18c(nt, time, auxc, amp)
      
!          hhg_2(i) = omg(i)**2 * abs(amp)**2
           ! --- integration by parts ---
           amp = exp(ci*omg(i)*time(nt)) * x_t(nt)   &
               - exp(ci*omg(i)*time(1))  * x_t(1)    &
               - ci*omg(i) * amp
      
           ! --- spectrum ---
           hhg_1(i) = abs(amp)**2



           do it=1,nt
              w = sin(ppi*(time(it)-time(1))/(time(nt)-time(1)))**2
              auxc(it) = exp(ci*omg(i)*time(it))* x_t(it) * w
           enddo
           call composite_simpson_18c(nt, time, auxc, amp)
      
!          hhg_2(i) = abs(amp)**2
           hhg_2(i) = omg(i)**2 * abs(amp)**2
        end do
      
      end subroutine


      subroutine compute_Qw(nch, krange, kk, bkw, b0w, Qw)
        implicit none
        integer, intent(in) :: krange, nch
        real(8), intent(in) :: kk(krange)
        complex(8), intent(in) :: bkw(krange, nch)
        complex(8), intent(in) :: b0w(nch)
        complex(8), allocatable, intent(out) :: Qw(:)

        complex(8) :: sum_kw(nch)
        complex(8) :: auxc(ksteps+1)

        integer :: w

        allocate(Qw(nch))

        do w=1,nch
           auxc = abs(bkw(:,w))**2
           call integr_over_krange(krange, kk, auxc, sum_kw(w))
           sum_kw(w) = sum_kw(w) / (2.d0*ppi)
        enddo

        Qw = abs(b0w)**2 + sum_kw




      end subroutine




      
      end module
