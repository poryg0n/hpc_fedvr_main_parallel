      module io_module
        implicit none
      contains


      !=========================================
      ! Write problem (human-readable)
      !=========================================
      subroutine write_problem_input(filename, struct_dir,       &
                                       nmax, snbr, nnbr,         &
                                       xmin, xmax, jac)
        implicit none
        character(*), intent(in) :: filename, struct_dir
        integer, intent(in) :: nmax, snbr, nnbr
        real(8), intent(in) :: xmin, xmax, jac
      
        integer :: unit
      
        open(newunit=unit, file=filename, status='replace')
      
        write(unit,*) "nmax       =", nmax
        write(unit,*) "snbr       =", snbr
        write(unit,*) "nnbr       =", nnbr
        write(unit,*) "xmin       =", xmin
        write(unit,*) "xmax       =", xmax
        write(unit,*) "jac        =", jac
        write(unit,*) "struct_dir =", struct_dir
      
        close(unit)
      end subroutine


      
      !=========================================
      ! Write problem (binary)
      !=========================================
      subroutine write_problem_bin(filename, workdir,            &
                                  nmax, snbr, nnbr,         &
                                  xmin, xmax, jac, xx, wx)
        implicit none
        character(*), intent(in) :: filename, workdir
        integer, intent(in) :: nmax, snbr, nnbr
        real(8), intent(in) :: xmin, xmax, jac
        real(8), intent(in) :: xx(:), wx(:)
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                   status='replace')
      
        write(unit) 1              ! version
        write(unit) nmax, snbr, nnbr
        write(unit) xmin, xmax, jac
        write(unit) xx
        write(unit) wx
        write(unit) workdir
      
        close(unit)
      end subroutine


      
      !=========================================
      ! Read problem (binary)
      !=========================================
      subroutine read_problem_bin(filename, struct_dir,         &
                                     nmax, snbr, nnbr,          &
                                     xmin, xmax, jac, xx, wx)
        implicit none
        character(*), intent(in) :: filename
        character(*), intent(out) :: struct_dir
        integer, intent(out) :: nmax, snbr, nnbr
        real(8), intent(out) :: xmin, xmax, jac
        real(8), allocatable, intent(out) :: xx(:), wx(:)
      
        integer :: unit, version
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                 status='old')
      
        read(unit) version
        read(unit) nmax, snbr, nnbr
        read(unit) xmin, xmax, jac
      
        allocate(xx(nmax), wx(nmax))
      
        read(unit) xx
        read(unit) wx
        read(unit) struct_dir
      
        close(unit)
      end subroutine
      
      
      

      subroutine write_dynamic_input(filename, dyn_dir, struct_dir,    &
                                 f0, omega0, pfai,                     &
                                 t_end, t_ini, nsteps, dt,             &
                                 noc, ntau, src_type,                  &
                                 nchannel, omg_1, omg_0, dw,           &
                                 order)

        implicit none
        character(*), intent(in) :: filename, struct_dir, dyn_dir
        integer, intent(in) ::  noc, ntau, nsteps, nchannel
        integer, intent(in) :: order, src_type
        real(8), intent(in) :: f0, omega0, pfai
        real(8), intent(in) :: t_end, t_ini, dt
        real(8), intent(in) :: omg_1, omg_0, dw
      
        integer :: unit
      
        open(newunit=unit, file=filename, status='replace')
      
        write(unit,*) "# Dynamic input parameters"
        write(unit,*) "f0        =", f0
        write(unit,*) "omega0    =", omega0
        write(unit,*) "pfai      =", pfai
        write(unit,*) "t_end     =", t_end
        write(unit,*) "t_ini     =", t_ini
        write(unit,*) "nsteps    =", nsteps
        write(unit,*) "dt        =", dt
        write(unit,*) "noc       =", noc
        write(unit,*) "ntau      =", ntau
        write(unit,*) "src_type  =", src_type
        write(unit,*) 
        write(unit,*) "nchannel-2=", nchannel
        write(unit,*) "omg_max   =", omg_1
        write(unit,*) "omg_min   =", omg_0
        write(unit,*) "dw        =", dw
        write(unit,*) "order     =", order
        write(unit,*) 
        write(unit,*) "struct    =", trim(struct_dir)
        write(unit,*) "dyn_dir   =", trim(dyn_dir)
      
        close(unit)
      end subroutine



      subroutine write_dynamic_bin(filename, workdir, struct_dir,      &
                                 f0, omega0, pfai,                     &
                                 t_end, t_ini, nsteps, dt0,            &
                                 noc, ntau, src_type,                  &
                                 nchannel, omg_1, omg_0, dw,           &
                                 order)

        implicit none
        character(*), intent(in) :: filename, workdir, struct_dir
        integer, intent(in) ::  noc, ntau, nsteps, nchannel
        integer, intent(in) :: order, src_type
        real(8), intent(in) :: f0, omega0, pfai
        real(8), intent(in) :: t_end, t_ini, dt0
        real(8), intent(in) :: omg_1, omg_0, dw
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                   status='replace')
      
        write(unit) 2              ! version
!       write(unit) filename
        write(unit) f0, omega0, pfai
        write(unit) t_end, t_ini
        write(unit) noc, ntau
        write(unit) nsteps
        write(unit) dt0
        write(unit) src_type

        write(unit) nchannel
        write(unit) omg_1
        write(unit) omg_0
        write(unit) dw

        write(unit) order
        write(unit) struct_dir
        write(unit) workdir
      
        close(unit)
      end subroutine



      subroutine read_dynamic_bin(filename, dyn_dir, struct_dir,      &
                                 f0, omega0, pfai,                     &
                                 t_end, t_ini, nsteps, dt0,           &
                                 noc, ntau, src_type,                 &
                                 nchannel, omg_1, omg_0, dw,          &
                                 order)
      
        implicit none
      
        character(*), intent(in) :: filename

        character(*), intent(out) :: struct_dir
        character(*), intent(out) :: dyn_dir
        integer, intent(out) :: nsteps, noc, ntau, nchannel
        integer, intent(out) :: order, src_type
        real(8), intent(out) :: f0, omega0, pfai
        real(8), intent(out) :: omg_1, omg_0, dw
        real(8), intent(out) :: t_end, t_ini, dt0
      
        integer :: unit, version
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                          status='old')
      
        read(unit) version
      
        if (version /= 2) then
           write(*,*) "Unsupported dynamic.bin version:", version
           stop
        end if
      
        read(unit) f0, omega0, pfai
        read(unit) t_end, t_ini
        read(unit) noc, ntau
        read(unit) nsteps
        read(unit) dt0
        read(unit) src_type

        read(unit) nchannel
        read(unit) omg_1
        read(unit) omg_0
        read(unit) dw

        read(unit) order
        read(unit) struct_dir
        read(unit) dyn_dir

!       select case(version)
!       case(2)
!          read(unit) f0, omega, pfai
!          read(unit) t_end, t_ini
!          read(unit) nsteps
!          read(unit) src_type
!          read(unit) order
!          read(unit) struct_dir
!       case default
!          write(*,*) "Unsupported dynamic.bin version:", version
!          stop
!       end select
      
        close(unit)

      
      end subroutine





      !=========================================
      ! Write eigenvalues
      !=========================================
      subroutine write_eigval_bin(filename, nmax, eigval)
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: nmax
        real(8), intent(in) :: eigval(:)
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                  status='replace')
      
        write(unit) nmax
        write(unit) eigval
      
        close(unit)
      end subroutine
      
      
      !=========================================
      ! Write eigenvectors
      !=========================================
      subroutine write_eigvec_bin(filename, nmax, eigvec)
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: nmax
        real(8), intent(in) :: eigvec(:,:)
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                   status='replace')
      
        write(unit) nmax
        write(unit) eigvec
      
        close(unit)
      end subroutine





      !=========================================
      ! Read eigenvalues
      !=========================================
      subroutine read_eigval_bin(filename, nmax, eigval)
        implicit none
        character(*), intent(in) :: filename
        integer, intent(out) :: nmax
        real(8), allocatable, intent(out) :: eigval(:)
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                   status='old')
      
        read(unit) nmax
        allocate(eigval(nmax))
        read(unit) eigval
      
        close(unit)
      end subroutine
      
      
      
      
      !=========================================
      ! Read eigenvectors
      !=========================================
      subroutine read_eigvec_bin(filename, nmax, eigvec)
        implicit none
        character(*), intent(in) :: filename
        integer, intent(out) :: nmax
        real(8), allocatable, intent(out) :: eigvec(:,:)
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                   status='old')
      
        read(unit) nmax
        allocate(eigvec(nmax,nmax))
        read(unit) eigvec
      
        close(unit)
      end subroutine




      subroutine write_wavefunction_bin(filename, nmax, nchan, t,   &
                                                omega, wf)
        implicit none
      
        character(*), intent(in) :: filename
        integer, intent(in) :: nmax, nchan
        real(8), intent(in) :: t
        real(8), intent(in) :: omega(nchan+2)
        complex(8), intent(in) :: wf(nmax,nchan+2)
      
        integer :: unit
        integer :: w
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                     status='replace')
      
        write(unit) nmax
        write(unit) nchan
        write(unit) t
        write(unit) omega
        write(unit) wf

        close(unit)
      
      end subroutine


      subroutine read_wavefunction_bin(filename, nmax, nchan, t,    &
                                                        omega, wf)
        implicit none
      
        character(*), intent(in) :: filename
        integer, intent(out) :: nmax, nchan
        real(8), intent(out) :: t
        real(8), allocatable, intent(out) :: omega(:)
        complex(8), allocatable, intent(out) :: wf(:,:)
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                          status='old')

        read(unit) nmax
        read(unit) nchan
        read(unit) t
      
        allocate(wf(nmax,nchan+2), omega(nchan+2))

        read(unit) omega
        read(unit) wf
      
        close(unit)
      
      end subroutine


      subroutine append_dyn_obs_bin(filename, t,                       &
!                            norm_1, norm_2,                           &
                             nchannel,                                 &
                             norm_x,                                   &
                             p0, pexc, pion,                           &
                             energy, dipole, momentum)
      
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: nchannel
        real(8), intent(in) :: t
        real(8), intent(in) :: norm_x(nchannel+2)
        real(8), intent(in) :: p0, pexc, pion
        complex(8), intent(in) :: energy, dipole, momentum
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted', &
             status='unknown', position='append')
      
        write(unit) t, nchannel, norm_x,                               &
                       p0, pexc, pion, energy, dipole, momentum
      
        close(unit)
      
      end subroutine


      subroutine write_observables_bin(filename,                       &
                                      nchan, nobs,                     &
                                      time, norm_t,                    &
                                      p0, pexc, pion,                  &
                                      nrg, dip, mom)
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: nobs, nchan
        real(8), intent(in) :: time(nobs)
        real(8), intent(in) :: norm_t(nchan+2,nobs)
        real(8), intent(in) :: p0(nobs), pexc(nobs), pion(nobs)
        complex(8), intent(in) :: nrg(nobs), dip(nobs), mom(nobs)
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                     status='replace')
      
        write(unit) nobs
        write(unit) nchan
        write(unit) time
!       write(unit) norm_1
!       write(unit) norm_2
        write(unit) norm_t
        write(unit) p0
        write(unit) pexc
        write(unit) pion
        write(unit) nrg
        write(unit) dip
        write(unit) mom
      
        close(unit)
      end subroutine



      subroutine read_observables_bin(filename,                        &
                                     nchan, nobs,                      &
                                     time,                             &
                                     norm_t,                           &
                                     p0, pexc, pion,                   &
                                     nrg, dip, mom)
        implicit none
        character(*), intent(in) :: filename
        integer, intent(out) :: nobs, nchan
        real(8), allocatable, intent(out) :: time(:)
        real(8), allocatable, intent(out) :: norm_t(:,:)

        real(8), allocatable, intent(out) :: p0(:), pexc(:), pion(:)
        complex(8), allocatable, intent(out) :: nrg(:), dip(:), mom(:)
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                  status='old')
      
        read(unit) nobs
        read(unit) nchan
      
        allocate(time(nobs))
        allocate(p0(nobs), pexc(nobs), pion(nobs))
        allocate(nrg(nobs), dip(nobs), mom(nobs))
        allocate(norm_t(nchan+2, nobs))

      
        read(unit) time
!       read(unit) norm_1
!       read(unit) norm_2
        read(unit) norm_t
        read(unit) p0
        read(unit) pexc
        read(unit) pion
        read(unit) nrg
        read(unit) dip
        read(unit) mom
      
        close(unit)

      end subroutine



      subroutine write_observables(filename, nobs,                    &
                                     time, energy, dipole, momentum)
      
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: nobs
        real(8), intent(in) :: time(nobs)
        complex(8), intent(in) :: dipole(nobs)
        complex(8), intent(in) :: momentum(nobs)
        complex(8), intent(in) :: energy(nobs)
      
        integer :: i, unit
      
        open(newunit=unit, file=filename, status='replace')
      
        do i = 1, nobs
           write(unit,'(7E20.10)') time(i),                     &
                                   real(energy(i)),             &
                                   real(dipole(i)),             &
                                   real(momentum(i)),           &
                                   aimag(dipole(i)),            &
                                   aimag(momentum(i)),          &
                                   aimag(energy(i))
        end do
      
        close(unit)
      
      end subroutine

      subroutine write_density_prob(filename, n, jacc,                &
                                              xx, wx, rho1, rho2)
      
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: n
        real(8), intent(in) :: jacc
        real(8), intent(in) :: xx(n), wx(n)
        complex(8), intent(in) :: rho1(n), rho2(n)
      
        integer :: i, unit
      
        open(newunit=unit, file=filename, status='replace')

        do i=1,n
           write(unit,*) xx(i), real(rho1(i)), real(rho2(i))
        enddo
      
        close(unit)
      
      end subroutine


      subroutine write_wf_observables(filename, nch, n, x, t,         &
                                     omega, ich,                      &
                                     re_wf, im_wf,                    &
                                     rho, arg)
        implicit none
      
        character(*), intent(in) :: filename
        integer, intent(in) :: n, nch, ich
        real(8), intent(in) :: t
        real(8), intent(in) :: x(n), omega(nch+2)
        real(8), intent(in) :: re_wf(n,nch+2), im_wf(n,nch+2)
        real(8), intent(in) :: rho(n,nch+2), arg(n, nch+2)
      
        integer :: i, unit
      
        open(newunit=unit, file=filename, status='replace')
      
        ! --- header
        write(unit,*) "# t = ", t
        write(unit,*) "# omega = ", omega(ich)
        write(unit,*) "# channel = ", ich
        write(unit,*) "# x  Re  Im  |psi|^2  arg"
      
!       write(unit,*) omega
        do i = 1, n
           write(unit,*) x(i), re_wf(i, ich), im_wf(i,ich),         &
                                              rho(i,ich), arg(i,ich)
        enddo
      
        close(unit)
      
      end subroutine


      subroutine write_pemd(filename, n, kk, ak, logscale)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: kk(n)
        complex(8), intent(in) :: ak(n)
        logical, intent(in) :: logscale
        character(len=*), intent(in) :: filename
      
        integer :: i, unit
        real(8) :: pk
      
        open(newunit=unit, file=filename, status='replace')
      
        do i = 1, n
           pk = abs(ak(i))**2
      
           if (logscale) then
              if (pk > 1d-20) then
                 write(unit,*) kk(i), log10(pk)
              else
                 write(unit,*) kk(i), -20.d0
              end if
           else
              write(unit,*) kk(i), pk
           end if
        end do
      
        close(unit)
      end subroutine


      subroutine write_hhg(filename, n_omg, omg,                       &
                                    qvc1, qvc2, logscale)
        implicit none
        integer, intent(in) :: n_omg
        real(8), intent(in) :: omg(n_omg)
        real(8), intent(in) :: qvc1(n_omg), qvc2(n_omg)
!       real(8), intent(in) :: qvc3(n_omg), qvc4(n_omg)
        logical, intent(in) :: logscale
        character(len=*), intent(in) :: filename
      
        integer :: i, unit
!       real(8) :: v1, v2, v3
        real(8) :: v1, v2
      
        open(newunit=unit, file=filename, status='replace')
      
        do i = 1, n_omg
      
           v1 = qvc1(i)
           v2 = qvc2(i)
!          v3 = qvc3(i)
      
           if (logscale) then
              write(unit,*) omg(i), log10(max(v1,1d-20)),             &
                                     log10(max(v2,1d-20))
           else
              write(unit,*) omg(i), v1, v2
           end if
      
        end do
      
        close(unit)
      end subroutine

      subroutine write_bkw_b0w(filename, nchannel, omega, krange, kk, &
                                                      b0w, bkw)
      
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: krange, nchannel
        real(8), intent(in) :: kk(krange)
        real(8), intent(in) :: omega(nchannel+2)
        complex(8), intent(in) :: bkw(krange, nchannel+1)
        complex(8), intent(in) :: b0w(nchannel+1)
!       complex(8), intent(in) :: Qw(nchannel+1)
      
        integer :: i, unit
      
        open(newunit=unit, file=filename, status='replace')

!       *** To rewrite ****      
!       write(unit,*) '# k, Re[b_k], Im[b_k], |b_k|^2'
!    
!       do i = 1, krange
!          write(unit,*) kk(i), real(bkw(i,j)), aimag(bkw(i,j)),       &
!                                                     abs(bkw(i,j))**2
!       enddo
!     
!       write(unit,*)
!       write(unit,*) '# b0(w):'
!       write(unit,*) real(b0w(j)), aimag(b0w(j)), abs(b0w(j))**2
!     
!       write(unit,*)
!       write(unit,*) '# Q(w):'
!       write(unit,*) Qw(j)
      
        close(unit)
      
      end subroutine


      subroutine write_Qw(filename, nchan, omega, Qw)
      
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: nchan
        real(8), intent(in) :: omega(nchan+2)
        complex(8), intent(in) :: Qw(nchan+1)
      
        integer :: j, unit
      
        open(newunit=unit, file=filename, status='replace')

!       write(unit,*) '# k, Re[b_k], Im[b_k], |b_k|^2'
        do j=1, nchan+1
           write(unit,*) j, omega(j+1), real(Qw(j)), Qw(j)
        enddo
      
        close(unit)
      
      end subroutine


      end module io_module
