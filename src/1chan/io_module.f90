      module io_module
        implicit none
      contains


      !=========================================
      ! Write problem (human-readable)
      !=========================================
      subroutine write_problem_input(filename, struct_dir,       &
                                       nmax, snbr, nnbr,         &
                                       xmin, xmax, jac, qq)
        implicit none
        character(*), intent(in) :: filename, struct_dir
        integer, intent(in) :: nmax, snbr, nnbr
        real(8), intent(in) :: xmin, xmax, jac, qq
      
        integer :: unit
      
        open(newunit=unit, file=filename, status='replace')
      
        write(unit,*) "nmax       =", nmax
        write(unit,*) "snbr       =", snbr
        write(unit,*) "nnbr       =", nnbr
        write(unit,*) "xmin       =", xmin
        write(unit,*) "xmax       =", xmax
        write(unit,*) "q          =", qq
        write(unit,*) 
        write(unit,*) "jac        =", jac
        write(unit,*) 
        write(unit,*) "struct_dir =", struct_dir
      
        close(unit)
      end subroutine


      
      !=========================================
      ! Write problem (binary)
      !=========================================
      subroutine write_problem_bin(filename, workdir,                 &
                                  nmax, snbr, nnbr,                   &
                                  xmin, xmax, qq, jac,                &
                                  xx, wx)

        implicit none
        character(*), intent(in) :: filename, workdir
        integer, intent(in) :: nmax, snbr, nnbr
        real(8), intent(in) :: xmin, xmax, jac, qq
        real(8), intent(in) :: xx(:), wx(:)
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                   status='replace')
      
        write(unit) 1              ! version
        write(unit) nmax, snbr, nnbr
        write(unit) xmin, xmax, qq, jac 
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
                                     xmin, xmax, qq, jac,       &
                                     xx, wx)
        implicit none
        character(*), intent(in) :: filename
        character(*), intent(out) :: struct_dir
        integer, intent(out) :: nmax, snbr, nnbr
        real(8), intent(out) :: xmin, xmax, jac, qq
        real(8), allocatable, intent(out) :: xx(:), wx(:)
      
        integer :: unit, version
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                 status='old')
      
        read(unit) version
        read(unit) nmax, snbr, nnbr
        read(unit) xmin, xmax, qq, jac
      
        allocate(xx(nmax), wx(nmax))
      
        read(unit) xx
        read(unit) wx
        read(unit) struct_dir
      
        close(unit)
      end subroutine
      
      
      

      subroutine write_dynamic_input(filename, dyn_dir, struct_dir,    &
                                 f0, omega, pfai, t_end, t_ini,        &
                                 nsteps, dt, noc, ntau,                &
                                 src_type, omg, order)
        implicit none
        character(*), intent(in) :: filename, struct_dir, dyn_dir
        integer, intent(in) ::  noc, ntau, nsteps
        integer, intent(in) :: order, src_type
        real(8), intent(in) :: f0, omega, pfai, omg
        real(8), intent(in) :: t_end, t_ini, dt
      
        integer :: unit
      
        open(newunit=unit, file=filename, status='replace')
      
        write(unit,*) "# Dynamic input parameters"
        write(unit,*) "f0        =", f0
        write(unit,*) "omega0    =", omega
        write(unit,*) "pfai      =", pfai
        write(unit,*) "t_end     =", t_end
        write(unit,*) "t_ini     =", t_ini
        write(unit,*) "nsteps    =", nsteps
        write(unit,*) "dt        =", dt
        write(unit,*) "noc       =", noc
        write(unit,*) "ntau      =", ntau
        write(unit,*) "src_type  =", src_type
        write(unit,*) 
        write(unit,*) "omg       =", omg
        write(unit,*) 
        write(unit,*) "order     =", order
        write(unit,*) 
        write(unit,*) "struct    =", trim(struct_dir)
        write(unit,*) "dyn_dir   =", trim(dyn_dir)
      
        close(unit)
      end subroutine



      subroutine write_dynamic_bin(filename, workdir, struct_dir,      &
                                 f0, omega, pfai,                      &
                                 t_end, t_ini, nsteps, dt0,            &
                                 noc, ntau, src_type,                  &
                                 omg, order)

        implicit none
        character(*), intent(in) :: filename, workdir, struct_dir
        integer, intent(in) ::  noc, ntau, nsteps
        integer, intent(in) :: order, src_type
        real(8), intent(in) :: f0, omega, pfai
        real(8), intent(in) :: dt0
        real(8), intent(in) :: omg
        real(8), intent(in) :: t_end, t_ini
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                   status='replace')
      
        write(unit) 2              ! version
!       write(unit) filename
        write(unit) f0, omega, pfai
        write(unit) t_end, t_ini
        write(unit) noc, ntau
        write(unit) nsteps
        write(unit) dt0
        write(unit) src_type
        write(unit) omg
        write(unit) order
        write(unit) struct_dir
        write(unit) workdir
      
        close(unit)
      end subroutine



      subroutine read_dynamic_bin(filename, dyn_dir, struct_dir,      &
                                 f0, omega, pfai,                     &
                                 t_end, t_ini, nsteps, dt0,           &
                                 noc, ntau, src_type,                 &
                                 omg, order)
      
        implicit none
      
        character(*), intent(in) :: filename

        character(*), intent(out) :: struct_dir
        character(*), intent(out) :: dyn_dir
        integer, intent(out) :: nsteps, noc, ntau
        integer, intent(out) :: order, src_type
        real(8), intent(out) :: f0, omega, pfai, omg
        real(8), intent(out) :: t_end, t_ini, dt0
      
        integer :: unit, version
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                          status='old')
      
        read(unit) version
      
        if (version /= 2) then
           write(*,*) "Unsupported dynamic.bin version:", version
           stop
        end if
      
        read(unit) f0, omega, pfai
        read(unit) t_end, t_ini
        read(unit) noc, ntau
        read(unit) nsteps
        read(unit) dt0
        read(unit) src_type

        read(unit) omg

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


      subroutine write_wavefunction_input(filename, t, omg, jac,       &
                                              wx, xx, eigvec, wf)
        implicit none
      
        character(len=*), intent(in) :: filename
        real(8), intent(in) :: t, omg, jac
        real(8), intent(in) :: wx(:), xx(:)
        real(8), intent(in) :: eigvec(:,:)
        complex(8), intent(in) :: wf(:)
      
        integer :: i, n, unit
        complex(8), allocatable :: wfc(:)
      
        open(newunit=unit, file=filename, status='replace')
      
        n = size(wf,1)
        allocate(wfc(n))
      
        ! --- header
        write(unit,*) "# t   = ", t
        write(unit,*) "# omg = ", omg
      
        ! --- eigenbasis → configuration space
        wfc = matmul(eigvec, wf)
      
        do i = 1, n
           wfc(i) = wfc(i) / ( wx(i) * dsqrt(jac) )
!          write(unit,*) x(i), real(wfc(i)), aimag(wfc(i))
           write(unit,'(3E20.10)') xx(i), real(wfc(i)), aimag(wfc(i))
        end do
      
!       write(unit,*) ""  ! separator between snapshots
      
        deallocate(wfc)
        close(unit)
      
      end subroutine


      subroutine write_wavefunction_bin(filename, nmax,               &
                                                t, omg, psi, phi)
        implicit none
      
        character(*), intent(in) :: filename
        integer, intent(in) :: nmax
        real(8), intent(in) :: t
        real(8), intent(in) :: omg
        complex(8), intent(in) :: psi(:), phi(:)
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                     status='replace')
      
        write(unit) nmax
        write(unit) t
        write(unit) omg
        write(unit) psi
        write(unit) phi
      
        close(unit)
      
      end subroutine


      subroutine read_wavefunction_bin(filename, nmax, t, omg,        &
                                                      psi, phi)
        implicit none
      
        character(*), intent(in) :: filename
        integer, intent(out) :: nmax
        complex(8), allocatable, intent(out) :: psi(:), phi(:)
        real(8), intent(out) :: omg
        real(8), intent(out) :: t
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                          status='old')
      
        read(unit) nmax
        read(unit) t
        read(unit) omg
      
        allocate(psi(nmax), phi(nmax))
      
        read(unit) psi
        read(unit) phi
      
        close(unit)
      
      end subroutine


      subroutine append_dyn_obs_bin(filename, t,                       &
                             norm_1, norm_2,                           &
                             p0, pexc, pion,                           &
                             energy, dipole, momentum)
      
        implicit none
        character(*), intent(in) :: filename
        real(8), intent(in) :: t, norm_1, norm_2
        real(8), intent(in) :: p0, pexc, pion
        complex(8), intent(in) :: energy, dipole, momentum
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted', &
             status='unknown', position='append')
      
        write(unit) t, norm_1, norm_2,                      &
                       p0, pexc, pion, energy, dipole, momentum
      
        close(unit)
      
      end subroutine


      subroutine write_observables_bin(filename, nobs, time,      & 
                                      norm_1, norm_2,          &
                                      p0, pexc, pion,          &
                                      nrg, dip, mom)
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: nobs
        real(8), intent(in) :: time(:), norm_1(:), norm_2(:)
        real(8), intent(in) :: p0(:), pexc(:), pion(:)
        complex(8), intent(in) :: nrg(:), dip(:), mom(:)
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                     status='replace')
      
        write(unit) nobs
        write(unit) time
        write(unit) norm_1
        write(unit) norm_2
        write(unit) p0
        write(unit) pexc
        write(unit) pion
        write(unit) nrg
        write(unit) dip
        write(unit) mom
      
        close(unit)
      end subroutine

      subroutine write_observables(filename, nt,                      &
                                     time, energy, dipole, momentum)
      
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: nt
        real(8), intent(in) :: time(nt)
        complex(8), intent(in) :: energy(nt), dipole(nt), momentum(nt)
      
        integer :: i, unit
      
        open(newunit=unit, file=filename, status='replace')
      
        do i = 1, nt
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


      subroutine read_observables_bin(filename, nobs, time,            &
                                     norm_1, norm_2,                   &
                                     p0, pexc, pion,                   &
                                     nrg, dip, mom)
        implicit none
        character(*), intent(in) :: filename
        integer, intent(out) :: nobs
        real(8), allocatable, intent(out) :: time(:),                  &
                                             norm_1(:), norm_2(:)

        real(8), allocatable, intent(out) :: p0(:), pexc(:), pion(:)
        complex(8), allocatable, intent(out) :: nrg(:), dip(:), mom(:)
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                  status='old')
      
        read(unit) nobs
      
        allocate(time(nobs), norm_1(nobs), norm_2(nobs))
        allocate(p0(nobs), pexc(nobs), pion(nobs))
        allocate(nrg(nobs), dip(nobs), mom(nobs))

      
        read(unit) time
        read(unit) norm_1
        read(unit) norm_2
        read(unit) p0
        read(unit) pexc
        read(unit) pion
        read(unit) nrg
        read(unit) dip
        read(unit) mom
      
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
           write(unit,'(3E20.10)') xx(i), real(rho1(i)), real(rho2(i))
        enddo
      
        close(unit)
      
      end subroutine


      subroutine write_pemd(filename, n, kk, ak, bkwT, logscale)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: kk(n)
        complex(8), intent(in) :: ak(n), bkwT(n)
        logical, intent(in) :: logscale
        character(len=*), intent(in) :: filename
      
        integer :: i, unit
        real(8) :: pk, pkwT
      
        open(newunit=unit, file=filename, status='replace')
      
        do i = 1, n
           pk = abs(ak(i))**2
           pkwT = abs(bkwT(i))**2
      
           if (logscale) then
              if (pk > 1d-20) then
                 write(unit,'(3E20.10)') kk(i), log10(pk), log10(pkwT)
              else
                 write(unit,'(3E20.10)') kk(i), -20.d0, log10(pkwT)
              end if
           else
        !     write(unit,'(7E20.10)') kk(i), pk, pkwT,                &
        !                                  real(ak), imag(ak),        &
        !                                  real(bkwT), imag(bkwT)
              write(unit,'(3E20.10)') kk(i), pk, pkwT
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
!             write(unit,'(3E20.10)') omg(i), log10(max(v1,1d-20)),   &
              write(unit,'(3E20.10)') omg(i), log10(max(v1,1d-20)),   &
                                     log10(max(v2,1d-20))
           else
              write(unit,'(3E20.10)') omg(i), v1, v2
           end if
      
        end do
      
        close(unit)
      end subroutine

      subroutine write_Qw(filename, krange, kk, bkw, b0w, Qw)
      
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: krange
        real(8), intent(in) :: kk(krange)
        complex(8), intent(in) :: bkw(krange)
        complex(8), intent(in) :: b0w
        complex(8), intent(in) :: Qw
      
        integer :: i, unit
      
        open(newunit=unit, file=filename, status='replace')
      
        write(unit,*) '# k, Re[b_k], Im[b_k], |b_k|^2'
      
        do i = 1, krange
!          write(unit,*) kk(i), real(bkw(i)), aimag(bkw(i)), abs(bkw(i))**2
           write(unit,'(4E20.10)') kk(i),                              &
                         real(bkw(i)), aimag(bkw(i)), abs(bkw(i))**2
        enddo
      
        write(unit,*)
        write(unit,*) '# b0(w):'
        write(unit,'(3E20.10)') real(b0w), aimag(b0w), abs(b0w)**2
      
        write(unit,*)
        write(unit,*) '# Q(w):'
        write(unit,*) Qw
      
        close(unit)
      
      end subroutine


      end module io_module
