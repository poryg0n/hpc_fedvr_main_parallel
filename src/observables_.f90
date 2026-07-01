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


      subroutine compute_dens_probab(n, jacc, wx, eigvec, wf, rho)

        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: jacc
        real(8), intent(in) :: wx(n)
        real(8), intent(in) :: eigvec(n,n)
        complex(8), intent(in) :: wf(n)
        complex(8), allocatable, intent(out) :: rho(:)
!       real(8), intent(out) :: norm

        complex(8) :: wfc(n)
        integer :: i

        allocate(rho(n))
        call eigen_to_dvr(n, jacc, wx, eigvec, wf, wfc)
        rho = conjg(wfc)* wfc * wx*wx*jacc

      end subroutine

      

      subroutine compute_dyn_observables(nmax,                     &
                                        xx, wx, jac,               &
                                        eigval, eigvec,            &
                                        norm_1,                    &
                                        p0, pexc, pion,            &
                                        dipole, momentum, energy,  &
                                        psi)
      
        implicit none
        integer, intent(in) :: nmax
        complex(8), intent(in) :: psi(nmax)
        real(8), intent(in) :: eigval(nmax), eigvec(nmax,nmax)
        real(8), intent(in) :: xx(nmax), wx(nmax), jac
      
        real(8), intent(out) :: norm_1
        real(8), intent(out) :: p0, pexc, pion
        complex(8), intent(out) :: energy, dipole, momentum
      
        complex(8) :: psi_x(nmax), p_psi(nmax)
!       complex(8) :: phi_x(nmax), p_phi(nmax)
        integer :: i
      

        ! --- transform to DVR ---
!       psi_x = matmul(eigvec, psi)
!       phi_x = matmul(eigvec, phi)
!       psi_x = psi_x/wx/dsqrt(jac)
!       phi_x = phi_x/wx/dsqrt(jac)
        call eigen_to_dvr(nmax, jac, wx, eigvec, psi, psi_x)
!       call eigen_to_dvr(nmax, jac, wx, eigvec, phi, phi_x)

        
        ! --- norm ---
        norm_1 = sum(abs(psi_x)**2 * (wx**2) * jac)
!       norm_2 = sum(abs(phi_x)**2 * (wx**2) * jac)
        
        ! --- populations ---
        p0   = abs(psi(1))**2
        pexc = sum(abs(psi(2:nmax))**2)
        pion = 1.d0 - norm_1
        
        ! --- dipole ---
        dipole = sum(conjg(psi_x) * xx * psi_x * (wx**2) * jac)
        
        ! --- energy (field-free) ---
        energy = sum(abs(psi)**2 * eigval)
        
        ! --- momentum (CLEAN) ---
        call apply_momentum_operator(nmax, eigvec, xx, wx, jac,        &
                                                     psi, p_psi, 0)
        momentum = sum(conjg(psi) * p_psi)
      
      end subroutine


      subroutine compute_pemd_zrp(nmax, krange, t_end,                 &
                           xx, wx, jacc,                               &
                           eigvec, eigval,                             &
                           psi0, phi0, psi, phi,                       &
                           k_max, kk, p_ion, p0,                       &
                           a0, ak,                                     &
                           b0wT, bkwT)
      
        implicit none
        integer, intent(in) :: nmax, krange
        real(8), intent(in) :: jacc
        real(8), intent(in) :: xx(nmax), wx(nmax)
        real(8), intent(in) :: eigval(nmax)
        real(8), intent(in) :: eigvec(nmax,nmax)
        real(8), intent(in) :: t_end, k_max
        complex(8), intent(in) :: psi0(nmax), phi0(nmax)
        complex(8), intent(in) :: psi(nmax), phi(nmax)

        real(8), intent(out) :: p_ion, p0
        complex(8), intent(out) :: a0, b0wT
        real(8), allocatable, intent(out) :: kk(:)
        complex(8), allocatable, intent(out) :: ak(:), bkwT(:)
      
        ! locals
        integer :: j, k, ij
        integer :: p, n_cont, ksteps_
        real(8) :: aux1, aux2, dk
        real(8) :: E0
        complex(8) :: auxc_
        real(8), allocatable :: Ek(:)
        complex(8) :: auxc(nmax), wfc_k(nmax)
        complex(8) :: psic0(nmax), phic0(nmax)
        complex(8) :: psic(nmax), phic(nmax)

        complex(8), parameter :: ci = (0.d0,1.d0)
!       real(8), parameter :: ppi = 3.141592653589793d0
        real(8), parameter :: ppi = 4.d0*datan(1.d0)

      

        call eigen_to_dvr(nmax, jacc, wx, eigvec, psi, psic)
        call eigen_to_dvr(nmax, jacc, wx, eigvec, phi, phic)

        ksteps_  = krange -1

        allocate(kk(krange))
        allocate(ak(krange))
        allocate(bkwT(krange))

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
      
           auxc  = conjg(wfc_k) * psic * wx*wx*jacc
           ak(j) = sum(auxc)

           auxc  = conjg(wfc_k) * phic * wx*wx*jacc
           bkwT(j) = sum(auxc)
     
        end do

        Ek  = 0.5d0 * kk**2
        E0  = - 0.5d0 * kapp**2

        ak   = exp(ci*Ek*t_end) * ak
        bkwT = exp(ci*Ek*t_end) * bkwT

        call eigen_to_dvr(nmax, jacc, wx, eigvec, psi0, psic0)

        a0 = (0.d0, 0.d0)
        auxc = conjg(psic0) * psic * wx*wx*jacc
        a0 = sum(auxc)
        a0 = exp(ci*E0*t_end)*a0

        b0wT = (0.d0, 0.d0)
        auxc = conjg(psic0) * phic * wx*wx*jacc
        b0wT = sum(auxc)
        b0wT = exp(ci*E0*t_end)*b0wT
  
        ! --- probabilities ---
        p_ion = 0.d0
        auxc= abs(ak)**2
        call integr_over_krange(ksteps_, kk, auxc, auxc_)
        p_ion = real(auxc_) / (2.d0*ppi)
      
        p0 = abs(a0)**2
      
      end subroutine


      subroutine compute_phi_elems(workdir,                  &
                           nmax, krange, t_end,              &
                           xx, wx, jacc,                             &
                           eigvec, eigval,                   &
                           wf0_1, wf0_2, wf_1, wf_2,                 &
                           omega, k_max, kk,                          &
                           a0, b0wT,                                  &
                           ak, bkwT,                                  &
                           b0w, bkw)
      
        implicit none
        integer, intent(in) :: nmax, krange
        real(8), intent(in) :: jacc, omega
        real(8), intent(in) :: xx(nmax), wx(nmax)
        real(8), intent(in) :: eigval(nmax)
        real(8), intent(in) :: eigvec(nmax,nmax)
        real(8), intent(in) :: t_end, k_max
        complex(8), intent(in) :: a0
        complex(8), intent(in) :: b0wT

        real(8), intent(in) :: kk(krange)
        complex(8), intent(in) :: wf0_1(nmax)
        complex(8), intent(in) :: wf0_2(nmax)
        complex(8), intent(in) :: wf_1(nmax)
        complex(8), intent(in) :: wf_2(nmax)
        complex(8), intent(in) :: ak(krange)
        complex(8), intent(in) :: bkwT(krange)

        complex(8), intent(out) :: b0w
        complex(8), allocatable, intent(out) :: bkw(:)
      
        character(255), intent(in) :: workdir

        ! locals
        integer :: j, k, l, ij
        integer :: p, n_cont, ksteps_
        real(8) :: E0, Ek_, Ekp_
        real(8) :: Ek(krange)
        real(8) :: aux1, aux2, dk, delta_kk
        complex(8) :: factor, denom
        complex(8) :: wfc_k(nmax), wfc_k_(nmax)
        complex(8) :: dwfc_0(nmax), dwfc_k(nmax)
        complex(8) :: pwfc_0(nmax), pwfc_k(nmax)
        complex(8) :: wfc0_1(nmax), wfc0_2(nmax)
        complex(8) :: wfc_1(nmax), wfc_2(nmax)
        complex(8) :: auxc(nmax)
        complex(8) :: vec_1(krange)
        complex(8) :: vec_2(krange)
        complex(8) :: vec_k(krange)
        complex(8) :: vec_0
                          
        complex(8) :: vec_sum
                                  
        complex(8), parameter :: ci = ( 0.d0, 1.d0 )
!       real(8), parameter :: ppi = 3.141592653589793d0
        real(8), parameter :: ppi = 4.d0*datan(1.d0)

        real(8) :: eta = 1.d-6
        complex(8) :: dk0_, dkk_
        complex(8) :: pk0_, p0k_, pkk_
        complex(8) :: auxc_1, auxc_2, auxc_3
        complex(8) :: pk0(krange), p0k(krange), pkk(krange)
        complex(8) :: pkkk(krange), dkk(krange)

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
      

!       call eigen_to_dvr(nmax, jacc, wx, eigvec, wf_1, wfc_1)
!       call eigen_to_dvr(nmax, jacc, wx, eigvec, wf_2, wfc_2)


        ksteps_=krange-1

        call apply_momentum_operator(nmax, eigvec, xx, wx,        &
                          jacc, wf0_1, auxc, 0)
        call eigen_to_dvr(nmax, jacc, wx, eigvec, auxc, pwfc_0)
     
        call eigen_to_dvr(nmax, jacc, wx, eigvec, wf0_1, wfc0_1)

        E0  = - 0.5d0 * kapp**2

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
     

              auxc_1 = auxc_1 + pkk(l) * ak(l) * dk0
              auxc_2 = auxc_2 + dkk(l) * ak(l) * dk0
              auxc_3 = auxc_3 + pkkk(l) * ak(l) * dk0

!             write(unit_pkk, '(8E20.10)') kk(j), kk(l),               &
!                                  real(pkk(l)), imag(pkk(l)),         &
!                                  real(dkk(l)), imag(dkk(l)),         &
!                                  real(pkkk(l)), imag(pkkk(l))


              vec_2(l) = pkk(l) * ak(l) / ( Ek_+omega-Ekp_ + ci*eta )
              vec_2(l) = exp(ci*( Ek_+omega-Ekp_ ) * t_end ) * vec_2(l)


           enddo

           write(*,*) j
           write(unit_pkl,'(7E20.10)') kk(j),                         &
                               real(auxc_1), imag(auxc_1),             &
                               real(auxc_2), imag(auxc_2),             &
                               real(auxc_3), imag(auxc_3)


           call integr_over_krange(ksteps_, kk, vec_2, vec_k(j))
           vec_k(j) = vec_k(j) / (2.d0*ppi)


           auxc = conjg(wfc_k) * pwfc_0 * wx*wx*jacc
           pk0(j) = sum(auxc)

           call differentiate(xx, wfc_k, dwfc_k)
           pwfc_k = -ci*dwfc_k


           auxc = conjg(wfc0_1) * pwfc_k * wx*wx*jacc
           p0k(j) = sum(auxc)


           auxc = conjg(wfc_k) * xx * wfc0_1 * wx*wx*jacc
           dk0_ = sum(auxc)
           dk0_ = -ci* ( Ek_ - E0 ) * dk0_

           denom = ( kapp**2 + kk(j)**2 )

           pk0_ = 2.d0 * kk(j) * kapp**(3.d0/2) / denom
           p0k_ = pk0_

     
           write(unit_pk0,'(5E20.10)') kk(j),                         &
                                     real(pk0(j)), real(p0k(j)),      & 
                                     real(dk0_), real(pk0_)

           vec_1(j) = p0k(j) * ak(j) / ( E0+omega-Ek_ + ci*eta )
           vec_1(j) = exp( ci * ( E0+omega-Ek_ ) * t_end ) * vec_1(j)

        enddo

        call integr_over_krange(ksteps_, kk, vec_1, vec_0)
        vec_0 = vec_0 / (2.d0*ppi)


        Ek = 0.5d0 * kk**2 

      
        vec_1 = pk0 * a0 / ( Ek+omega-E0 )
        vec_1 = exp(ci * ( Ek+omega-E0 ) * t_end ) * vec_1

        allocate(bkw(krange))

        b0w = b0wT + vec_0
        bkw = bkwT + vec_1 +  vec_k

        do j=1,krange
           write(unit_vec,'(7E20.10)') kk(j), real(vec_0), imag(vec_0),&
                                real(vec_1(j)), imag(vec_1(j)),        &
                                real(vec_k(j)), imag(vec_k(j))
           write(unit_b0kw,'(5E20.10)') kk(j), real(b0w), imag(b0w),   &
                                real(bkw(j)), imag(bkw(j))
        enddo

        vec_sum = sum(vec_1)
        write(*,*) "The sum for vec_1 is ", vec_sum


        close(unit_pk0)
        close(unit_pkk)
        close(unit_pkl)
        close(unit_vec)
        close(unit_b0kw)
      
      end subroutine




      subroutine compute_hhg_ip(nt, time, x_t, wsteps, omg,          &
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


      subroutine compute_Qw(ksteps, kk, bkw, b0w, Qw)
        implicit none
        integer, intent(in) :: ksteps
        real(8), intent(in) :: kk(ksteps+1)
        complex(8), intent(in) :: b0w
        complex(8), intent(in) :: bkw(ksteps+1)
        complex(8), intent(out) :: Qw

        complex(8) :: sum_kw
        complex(8) :: auxc(ksteps+1)


        auxc = abs(bkw)**2
        call integr_over_krange(ksteps, kk, auxc, sum_kw)
        sum_kw = sum_kw / (2.d0*ppi)
        Qw = abs(b0w)**2 + sum_kw

      end subroutine



      end module
