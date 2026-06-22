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


!     subroutine compute_wf_observables(n, nch, ich,                   &
      subroutine compute_wf_observables(n, nch,                        &
                                       jacc, omega, wx, eigvec, wf,    &
                                       re_wf, im_wf, rho, arg)
        implicit none
      
!       integer, intent(in) :: n, nch, ich
        integer, intent(in) :: n, nch
        real(8), intent(in) :: jacc
        real(8), intent(in) :: wx(n), omega(nch+1)
        real(8), intent(in) :: eigvec(n,n)
        complex(8), intent(in) :: wf(n, nch+1)   ! multi-channel
      
        real(8), allocatable, intent(out) :: re_wf(:,:), im_wf(:,:),   &
                                             rho(:,:), arg(:,:)
      
        complex(8) :: wfc(n, nch+1)
        integer :: i, w
      
        allocate(re_wf(n,nch+1), im_wf(n,nch+1))
        allocate(rho(n,nch+1),arg(n,nch+1))
      
!       ! --- select channel
!       wfc = matmul(eigvec, wf)
!     
!       ! --- remove quadrature weights
!       do i = 1, n
!          wfc(i) = wfc(i) / (wx(i) * dsqrt(jacc))
!       enddo

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

      

      subroutine compute_dyn_observables(nmax, nchan, vec,      &
                                        xx, wx, jac,               &
                                        eigval, eigvec,            &
!                                       norm_1, norm_2,            &
                                        norm_x,                     &
                                        p0, pexc, pion,            &
                                        dipole, momentum, energy)
      
        implicit none
        integer, intent(in) :: nmax, nchan
        complex(8), intent(in) :: vec(nmax,nchan+1)
        real(8), intent(in) :: eigval(nmax), eigvec(nmax,nmax)
        real(8), intent(in) :: xx(nmax), wx(nmax), jac
      
!       real(8), intent(out) :: norm_1, norm_2 
        real(8), intent(out) :: p0, pexc, pion
        real(8), intent(out) :: norm_x(nchan+1)
        complex(8), intent(out) :: energy, dipole, momentum
      
!       complex(8) :: psi_x(nmax), p_psi(nmax)
!       complex(8) :: phi_x(nmax), p_phi(nmax)
        complex(8) :: vec_x(nmax, nchan+1)
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
!       call zgemm('n','n', nmax, nchannel+2, nmax, c1, eigvec,        &
!                  nmax, vec, nmax, c0, vec_x, nmax)
!       vec_x = matmul(eigvec, vec)
        call eigen_to_dvr(nmax, nchan, jac, wx, eigvec, vec, vec_x)


!       do k=1,nchannel+1
!          vec_x(:,k) = vec_x(:,k)/wx/dsqrt(jac)
!       enddo

        
        ! --- norm ---
!       norm_1 = sum(abs(psi_x)**2 * (wx**2) * jac)
!       norm_2 = sum(abs(phi_x)**2 * (wx**2) * jac)
        do k=1,nchan+1
           norm_x(k)  = sum(abs(vec_x(:,k))**2 * (wx**2) * jac)
        enddo
        
        ! --- populations ---
        p0   = abs(vec(1,1))**2
        pexc = sum(abs(vec(2:nmax,1))**2)
        pion = 1.d0 - norm_x(1)
        
        ! --- dipole ---
        dipole = sum(conjg(vec_x(:,1)) * xx * vec_x(:,1) *             &
                                                   (wx**2) * jac)
        
        ! --- energy (field-free) ---
        energy = sum(abs(vec(:,1))**2 * eigval)
        
        ! --- momentum (CLEAN) ---
        call apply_momentum_operator(nmax, eigvec, xx, wx, jac,        &
                                                vec(:,1), p_psi, 0)
        momentum = sum(conjg(vec(:,1)) * p_psi)
      
      end subroutine




      subroutine compute_pemd_zrp(nmax, nch, krange, t_end,            &
                           xx, wx, jacc,                               &
                           eigvec, eigval,                             &
                           wf0_0, wf,                        &
                           k_max, kk, p_ion, p0, a0, ak,               &
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
        allocate(bkwT(krange,nch+1))


        ak = (0.d0, 0.d0)
        dk = k_max/(krange/2-1)


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
           
           do l=2,nch+1
              auxcw(:,l)  = conjg(wfc_k) * wfc(:,l) * wx*wx*jacc
              bkwT(j,l) = sum(auxcw(:,l))
           enddo

        end do


      
        ! --- probabilities ---
        p_ion = 0.d0
        call integr_over_range(krange, kk, ak, auxc_)
        p_ion = real(auxc_) / (2.d0*ppi)
        write(*,*) p_ion


        a0 = (0.d0, 0.d0)
!       call eigen_to_dvr(nmax, jacc, wx, eigvec, wf0, wfc0)
        wfc0_0 = matmul(eigvec, wf0_0)
        wfc0_0 = wfc0_0 / wx /dsqrt(jacc)
        auxcw(:,1) = conjg(wfc0_0) * wfc(:,1) * wx*wx*jacc
  
        a0 = sum(auxcw(:,1))
        a0 = exp(ci*eigval(1)*t_end)*a0
      
        p0 = abs(a0)**2
      
      end subroutine



      subroutine compute_transit_mat(workdir,                        &
                           nmax, krange, t_end,            &
                           xx, wx, jacc,                             &
                           eigvec, eigval,                           &
                           wf0,                                      &
                           k_max, kk, a0, ak,                        &
                           pk0, p0k, pkq_aq)
      
        implicit none
        integer, intent(in) :: nmax, krange
        real(8), intent(in) :: jacc
        real(8), intent(in) :: xx(nmax), wx(nmax)
        real(8), intent(in) :: eigval(nmax)
        real(8), intent(in) :: eigvec(nmax,nmax)
        real(8), intent(in) :: t_end, k_max
        complex(8), intent(in) :: a0

        character(255) :: workdir

        real(8), intent(in) :: kk(krange)
        complex(8), intent(in) :: wf0(nmax)
        complex(8), intent(in) :: ak(krange)

        complex(8), allocatable, intent(out) :: pk0(:), p0k(:),        &
                                                   pkq_aq(:)


        ! locals
        integer :: j, k, l, p
        real(8) :: E0, Ek, Ek_
        real(8) :: aux1, aux2, dk, delta_kk
        real(8) :: auxr1(krange/2), auxr2(krange/2)
        complex(8) :: factor, denom
        complex(8) :: wfc_k(nmax), wfc_k_(nmax)
        complex(8) :: dwfc_k(nmax), pwfc_k(nmax)
        complex(8) :: dwfc_0(nmax), pwfc_0(nmax)

        complex(8) :: wfc0(nmax)
        complex(8) :: auxc(nmax)

        complex(8), parameter :: ci = (0.d0,1.d0)
        real(8), parameter :: ppi = 4.d0*datan(1.d0)

        real(8) :: eta = 1.d-6
        complex(8) :: dk0_, dkk_
        complex(8) :: pk0_, p0k_, pkk, pkk_
        complex(8) :: pkqaq, auxc_2, auxc_3
        complex(8) :: dkk(krange)

        integer :: unit_pk0, unit_pkqaq, unit_pkl_
        open(newunit=unit_pk0, file=trim(workdir)//"/pk0.dat",         &
                                                   status="replace")
        open(newunit=unit_pkqaq, file=trim(workdir)//"/pkq_aq.dat",    &
                                                    status="replace")
      
        allocate(pk0(krange), p0k(krange), pkq_aq(krange))

        p=0
        do p=1,nmax
           if (eigval(p).gt.0.0d0) then
              exit
           end if
        enddo
!       write(*,*) p

        dk = k_max/(krange/2-p+1)

!       call eigen_to_dvr(nmax, jacc, wx, eigvec, wf0, wfc0)
        wfc0 = matmul(eigvec, wf0(:))
        wfc0 = wfc0 / wx / dsqrt(jacc) 

        call differentiate(xx, wfc0, dwfc_0)
        pwfc_0 = -ci * dwfc_0

!       Ek = 0.5d0 * kk**2 

        do j=1,krange
      
           wfc_k = exp(ci*kk(j)*xx) +                                  &
                   (ci*kapp/(-abs(kk(j)) - ci*kapp)) *                 &
                   exp(-ci*abs(kk(j)*xx))
       
           Ek  = 0.5d0 * kk(j)**2 

           pkqaq = 0.d0
           auxc_2 = 0.d0
           auxc_3 = 0.d0
 
           do l=1,krange

!             if (j == l) cycle
      
              wfc_k_ = exp(ci*kk(l)*xx) +                              &
                      (ci*kapp/(-abs(kk(l)) - ci*kapp)) *              &
                      exp(-ci*abs(kk(l)*xx))
     
              Ek_ = 0.5d0 * kk(l)**2 

              call differentiate(xx, wfc_k_, dwfc_k)
              pwfc_k = -ci * dwfc_k

              auxc = conjg(wfc_k) * pwfc_k * wx*wx*jacc
              pkk = sum(auxc)

              auxc = conjg(wfc_k) * xx * wfc_k_ * wx*wx*jacc
              dkk_ = sum(auxc)
              dkk_ = -ci* ( Ek - Ek_ ) * dkk_

              ! *** analytical formula
              if(j.eq.l) then
                 delta_kk = 1.d0/dk
              else
                 delta_kk = 0.d0
              end if

!             delta_kk = 0.d0
              denom  = ( kk(j)**2 - kk(l)**2 + ci*eta )
              factor = ( kk(l)*abs(kk(j)) / (abs(kk(j)) -ci*kapp)      &
                      -  kk(j)*abs(kk(l)) / (abs(kk(l)) +ci*kapp) )


              pkk_ = 2.d0*ppi*kk(j) * delta_kk                         &
                                   - 2.d0 * kapp * factor/denom
     

!             pkqaq = pkqaq + pkk * ak(l) * dk
!             auxc_2 = auxc_2 + dkk_* ak(l) * dk
!             auxc_3 = auxc_3 + pkk_* ak(l) * dk

              pkqaq = pkqaq   + pkk * ak(l) 
              auxc_2 = auxc_2 + dkk_* ak(l) 
              auxc_3 = auxc_3 + pkk_* ak(l)





!             write(unit_pkk, *) kk(j), kk(l),                        &
!                        real(pkk(l)), imag(pkk(l)),                  &
!                        real(dkk_), imag(dkk_),                      &
!                        real(pkk_), imag(pkk_)

!             vec_2(l) = pkk(l) * ak(l) / ( Ek_ + omega - Ekp + ci*eta )
!             vec_2(l) = exp(ci*( Ek_+omega-Ekp ) * t_end ) * vec_2(l)

!             if (j.ne.l) then  
!                b_pv(j) = p(j,l) * ak(l) / (kk(j)**2 - kkp(l)**2) * dk
!             end if

           enddo

           pkq_aq(k) = pkqaq

           E0 = -kapp**2 / 2
!          E0  = - 0.5d0 * kapp**2
!          E0  = eigval(1)

           call differentiate(xx, wfc0, dwfc_0)
           pwfc_0 = -ci * dwfc_0
           auxc = conjg(wfc_k) * pwfc_0 * wx*wx*jacc
           pk0(j) = sum(auxc)


           call differentiate(xx, wfc_k, dwfc_k)
           pwfc_k = -ci * dwfc_k
           auxc = conjg(wfc0) * pwfc_k * wx*wx*jacc
           p0k(j) = sum(auxc)


           auxc = conjg(wfc_k) * xx * wfc0 * wx*wx*jacc
           dk0_ = sum(auxc)
           dk0_ = -ci* ( Ek - E0 ) * dk0_

           denom = ( kapp**2 + kk(j)**2 )

           pk0_ = 2.d0 * kk(j) * kapp**(1.5) / denom


!          write(unit_pk0, *) kk(j), real(pk0(j)), real(p0k(j)),      &
!                                             real(dk0_), real(pk0_)

!          write(unit_pkqaq, *) kk(j),                                 &
!                       abs(real(pkqaq)), abs(imag(pkqaq)),            &
!                       abs(real(auxc_2)), abs(imag(auxc_2)),          &
!                       abs(real(auxc_3)), abs(imag(auxc_3))

           write(unit_pk0, *) kk(j), abs(pk0(j))**2 * dk,              &
                                     abs(p0k(j))**2 * dk,              &
                                     abs(dk0_)**2   * dk,              &
                                     abs(pk0_)

           write(unit_pkqaq, *) kk(j),                                 &
                   abs((pkqaq)) **2 *dk,                               &
                   abs((auxc_2))**2 *dk,                               &
                   abs((auxc_3))**2 *dk
        enddo


        close(unit_pk0)
        close(unit_pkqaq)
      
      end subroutine


       subroutine compute_phi_elems(workdir,                          &
                            nmax, krange, nchan, t_end,            &
                            xx, wx, jacc,                             &
                            eigvec, eigval,                           &
                            wf, wf0, k_max, kk, a0, ak,               &
                            p0k, pk0, pkq_aq,                         &
                            omega, b0w, bkw)
       
         implicit none
         integer, intent(in) :: nmax, krange, nchan
         real(8), intent(in) :: jacc
         real(8), intent(in) :: xx(nmax), wx(nmax)
         real(8), intent(in) :: eigval(nmax)
         real(8), intent(in) :: eigvec(nmax,nmax)
         real(8), intent(in) :: t_end, k_max
         real(8), intent(in) :: omega(nchan+2)
         complex(8), intent(in) :: a0
         complex(8), intent(in) :: ak(krange)
         complex(8), intent(in) :: pk0(krange), p0k(krange)
         complex(8), intent(in) :: pkq_aq(krange)
 
         character(255) :: workdir
 
         real(8), intent(in) :: kk(krange)
         complex(8), intent(in) :: wf0(nmax, nchan+2)
         complex(8), intent(in) :: wf(nmax, nchan+2)
 
         complex(8), allocatable, intent(out) :: b0w(:)
         complex(8), allocatable, intent(out) :: bkw(:, :)
       
!        ! locals
         integer :: j, k, l, w
!        integer :: p
         real(8) :: E0, Ek_, Ekp
         real(8) :: omega_w
         real(8) :: Ek(krange)
!        real(8) :: dk, delta_kk
!        real(8) :: auxr1(krange/2), auxr2(krange/2)
!        complex(8) :: factor, denom
         complex(8) :: wfc_k(nmax), wfc_k_(nmax)
         complex(8) :: dwfc_k(nmax), pwfc_k(nmax)
         complex(8) :: dwfc_0(nmax), pwfc_0(nmax)
!
         complex(8) :: wfc0(nmax, nchan+2)
         complex(8) :: wfc(nmax, nchan+2)
!
         complex(8) :: auxc(nmax)
         complex(8) :: vec_1(krange)
         complex(8) :: vec_2(krange)
         complex(8) :: vec_k(krange)
         complex(8) :: vec_0
         complex(8) :: b0wT
         complex(8) :: bkwT(krange, nchan+1)
         complex(8), parameter :: ci = (0.d0,1.d0)
         real(8), parameter :: ppi = 4.d0*datan(1.d0)
 
         real(8) :: eta = 1.d-6

         integer :: unit_b0w, unit_bkw
         open(newunit=unit_b0w, file=trim(workdir)//"/b0w.dat",         &
                                                     status="replace")
         open(newunit=unit_bkw, file=trim(workdir)//"/bkw.dat",         &
                                                     status="replace")

        allocate(bkw(krange, nchan+1))
        allocate(b0w(nchan+1))

!      
!        p=0
!        do p=1,nmax
!           if (eigval(p).gt.0.0d0) then
!              exit
!           end if
!        enddo
!!       write(*,*) p
!
!        dk = k_max/(krange/2-p+1)
!
         call eigen_to_dvr(nmax, nchan, jacc, wx, eigvec, wf0, wfc0)

         call eigen_to_dvr(nmax, nchan, jacc, wx, eigvec, wf, wfc)
!     
!
         Ek = 0.5d0 * kk**2 
         E0 = 0.5d0 * kapp**2


         write(*,*)  "debugging"
         do w=1,nchan+1
            omega_w = omega(w+1)
 
            do k=1,krange
               Ek_ = Ek(k)


               wfc_k = exp(ci*kk(j)*xx) +                              &
                   (ci*kapp/(-abs(kk(j)) - ci*kapp)) *                 &
                   exp(-ci*abs(kk(j)*xx))

!              Ekp = Ek(l)

               vec_2 = pkq_aq / ( Ek_+omega_w-Ek + ci*eta )
               vec_2 = exp(ci*( Ek_+omega_w-Ek ) * t_end ) * vec_2

!              write(*,*)  "inner loop k' passed"
               call integr_over_range(krange, kk, vec_2, vec_k(k))

 
               auxc = conjg(wfc_k) * wfc(:,w+1) * wx*wx*jacc
               bkwT(k,w) = sum(auxc)
               bkwT(k,w) = exp(ci*Ek_*t_end) * bkwT(k,w)

            enddo
            write(*,*)  "counting", w

            vec_1 = p0k * ak  / ( E0+omega_w-Ek + ci*eta )
            vec_1 = exp(ci * ( E0+omega_w-Ek ) * t_end )* vec_1

            call integr_over_range(krange, kk, vec_1, vec_0)

            auxc = conjg(wfc0(:,1)) * wfc(:,w+1) * wx*wx*jacc
            b0wT = sum(auxc)
!           b0wT = exp(ci*eigval(1)*t_end) * b0wT
            b0wT = exp(ci*E0*t_end) * b0wT

            vec_1 = pk0 * a0 / ( Ek+omega_w-E0 )
            vec_1 = exp(ci * ( Ek+omega_w-E0 ) * t_end ) * vec_1
 
            b0w(w) = b0wT
            bkw(:,w) = bkwT(:,w) 
!           b0w(w) = b0wT + vec_0
!           bkw(:,w) = bkwT(:,w) + vec_1(:) +  vec_k(:)
!           write(*,*)  "main loop nchannel  passed"
         enddo
 


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


      subroutine compute_Qw(nchannel, krange, kk, bkw, b0w, Qw)
        implicit none
        integer, intent(in) :: krange, nchannel
        real(8), intent(in) :: kk(krange)
        complex(8), intent(in) :: bkw(krange, nchannel+1)
        complex(8), intent(in) :: b0w(nchannel+1)
        complex(8), allocatable, intent(out) :: Qw(:)

        complex(8) :: auxc(nchannel+1)
        integer :: w

        allocate(Qw(nchannel+1))

        do w=1,nchannel+1
           call integr_over_range(krange, kk, bkw(:,w), auxc(w))
        enddo
        Qw = abs(b0w)**2 + auxc

      end subroutine




      
      end module
