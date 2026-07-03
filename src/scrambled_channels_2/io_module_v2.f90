      module io_module
        use space_time_ops
        implicit none
      contains


      !=========================================
      ! Write problem (human-readable)
      !=========================================
      subroutine write_structure_input(filename, struct_dir,       &
                                       nmax, snbr, nnbr,         &
                                       xmin, xmax, qq, jac)
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
      subroutine write_structure_bin(filename, workdir,            &
                                  nmax, snbr, nnbr,         &
                                  xmin, xmax, qq, jac,      &
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
        write(unit) xmin, xmax, jac, qq
        write(unit) xx
        write(unit) wx
        write(unit) workdir
      
        close(unit)
      end subroutine


      
      !=========================================
      ! Read problem (binary)
      !=========================================
      subroutine read_structure_bin(filename, struct_dir,         &
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
        read(unit) xmin, xmax, jac, qq
      
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
                                 nch, omg_max, omg_min,                &
                                 wstep, run,                           &
                                 order)


        implicit none
        character(*), intent(in) :: filename, struct_dir, dyn_dir
        integer, intent(in) :: noc, ntau, nsteps, nch, run
        integer, intent(in) :: order, src_type
        real(8), intent(in) :: f0, omega0, pfai
        real(8), intent(in) :: t_end, t_ini, dt
        real(8), intent(in) :: omg_max, omg_min, wstep

      
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
        write(unit,*) "nchan     =", nch
        write(unit,*) "omg_max   =", omg_max
        write(unit,*) "omg_min   =", omg_min
        write(unit,*) "wstep     =", wstep
        write(unit,*) 
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
                                 nch, omg_max, omg_min,                &
                                 wstep, omega, run,                    &
                                 order)

        implicit none
        character(*), intent(in) :: filename, workdir, struct_dir
        integer, intent(in) :: noc, ntau, nsteps, nch, run
        integer, intent(in) :: order, src_type
        real(8), intent(in) :: f0, omega0, pfai
        real(8), intent(in) :: omg_max, omg_min, wstep
        real(8), intent(in) :: t_end, t_ini, dt0
        real(8), intent(in) :: omega(nch)
      
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

        write(unit) nch, run
        write(unit) omg_max, omg_min, wstep
        write(unit) omega

        write(unit) order
        write(unit) struct_dir
        write(unit) workdir
      
        close(unit)
      end subroutine



      subroutine read_dynamic_bin(filename, dyn_dir, struct_dir,      &
                                 f0, omega0, pfai,                     &
                                 t_end, t_ini, nsteps, dt0,           &
                                 noc, ntau, src_type,                 &
                                 nch, omg_max, omg_min,                &
                                 wstep, omega, run,                    &
                                 order)
      
        implicit none
      
        character(*), intent(in) :: filename

        character(*), intent(out) :: struct_dir
        character(*), intent(out) :: dyn_dir
        integer, intent(out) :: nsteps, noc, ntau, nch, run
        integer, intent(out) :: order, src_type
        real(8), intent(out) :: f0, omega0, pfai
        real(8), intent(out) :: omg_max, omg_min, wstep
        real(8), allocatable, intent(out) :: omega(:)
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

        read(unit) nch, run

        read(unit) omg_max, omg_min, wstep
        allocate(omega(nch))
        read(unit) omega


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


      subroutine write_wavefun_input(filename, nch, k, t, jac,  &
                                            omega, wx, xx, eigvec, wf)
        implicit none

        character(len=*), intent(in) :: filename
        integer, intent(in) :: k, nch
        real(8), intent(in) :: t, jac
        real(8), intent(in) :: wx(:), xx(:), omega(:)
        real(8), intent(in) :: eigvec(:,:)
        complex(8), intent(in) :: wf(:,:)

        integer :: i, j, n, unit
        complex(8), allocatable :: wfc(:,:)
        complex(8), allocatable :: auxc(:,:)

        open(newunit=unit, file=filename, status='replace')

        n = size(wf,1)
        allocate(wfc(n,nch))
        allocate(auxc(n,nch))

        ! --- header
        write(unit,*) "# t   = ", t
        write(unit,*) "# omg = ", omega(k)

        ! --- eigenbasis → configuration space
!       call eigen_to_dvr(n, nch, jac,                   &
!                                wx, eigvec, wf, wfc)

        do j=1,nch
           call eigen_to_dvr(n, jac,                   &
                                 wx, eigvec, wf(:,j), auxc(:,j))
           wfc(:,j) = auxc(:,j)
        enddo

        do i = 1, n
!          write(unit,*) x(i), real(wfc(i)), aimag(wfc(i))
           write(unit,'(3E20.10)') xx(i), real(wfc(i,k)), aimag(wfc(i,k))
        end do

!       write(unit,*) ""  ! separator between snapshots

        deallocate(wfc)
        close(unit)

      end subroutine





      subroutine write_wavefun_bin(filename, mode, nch, nmax, t,   &
                                                omega, wf)
        implicit none
      
        character(*), intent(in) :: filename
        integer, intent(in) :: nmax, nch, mode
        real(8), intent(in) :: t
        real(8), intent(in) :: omega(nch)
        complex(8), intent(in) :: wf(nmax,nch)
      
        integer :: unit
        integer :: w


        if (mode.eq.1) then
        open(newunit=unit, file=filename, form='unformatted', &
             status='unknown', position='append')
        else
           open(newunit=unit, file=filename, form='unformatted',       &
                                              status='replace')
        end if
      
      
        write(unit) nch
        write(unit) nmax
        write(unit) t
        write(unit) omega
        write(unit) wf

        close(unit)
      
      end subroutine


      subroutine read_wavefun_bin(filename, nch, nmax, t,     &
                                                        omega, wf)
        implicit none
      
        character(*), intent(in) :: filename
        integer, intent(out) :: nmax, nch
        real(8), intent(out) :: t
        real(8), allocatable, intent(out) :: omega(:)
        complex(8), allocatable, intent(out) :: wf(:,:)
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                          status='old')

        read(unit) nch
        read(unit) nmax
        read(unit) t
      
        allocate(wf(nmax,nch), omega(nch))

        read(unit) omega
        read(unit) wf
      
        close(unit)
      
      end subroutine


      subroutine append_dyn_obs_bin(filename,                        &
                             nch,                    &
                             t,                                &
                             norm_,                                   &
                             p0_, pexc_, pion_,                     &
                             energy_, dipole_, momentum_)
      
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: nch
        real(8), intent(in) :: t
        real(8), intent(in) :: norm_(nch)
        real(8), intent(in) :: p0_, pexc_, pion_
        complex(8), intent(in) :: energy_, dipole_, momentum_
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted', &
             status='unknown', position='append')
      
        write(unit) t, nch, norm_,                                   &
                       p0_, pexc_, pion_, dipole_, momentum_, energy_
      
        close(unit)
      
      end subroutine


      subroutine write_observables_bin(filename,                       &
                                      nch, nobs,                     &
                                      time, norm_t,                    &
                                      p0, pexc, pion,                  &
                                      dip, mom, nrg)
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: nobs, nch
        real(8), intent(in) :: time(nobs)
        real(8), intent(in) :: norm_t(nobs,nch)
        real(8), intent(in) :: p0(nobs), pexc(nobs), pion(nobs)
        complex(8), intent(in) :: nrg(nobs), dip(nobs), mom(nobs)
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                     status='replace')
      
        write(unit) nch
        write(unit) nobs
        write(unit) time
        write(unit) norm_t
        write(unit) p0
        write(unit) pexc
        write(unit) pion
        write(unit) dip
        write(unit) mom
        write(unit) nrg
      
        close(unit)
      end subroutine



      subroutine read_observables_bin(filename,                        &
                                     nch, nobs,                      &
                                     time_t,                           &
                                     norm_t,                           &
                                     p0, pexc, pion,                   &
                                     dip, mom, nrg)
        implicit none
        character(*), intent(in) :: filename
        integer, intent(out) :: nobs, nch
        real(8), allocatable, intent(out) :: time_t(:)
        real(8), allocatable, intent(out) :: norm_t(:,:)

        real(8), allocatable, intent(out) :: p0(:), pexc(:), pion(:)
        complex(8), allocatable, intent(out) :: nrg(:), dip(:), mom(:)
      
        integer :: unit
      
        open(newunit=unit, file=filename, form='unformatted',          &
                                                  status='old')
      
        read(unit) nch
        read(unit) nobs
      
        allocate(time_t(nobs))
        allocate(p0(nobs), pexc(nobs), pion(nobs))
        allocate(nrg(nobs), dip(nobs), mom(nobs))
        allocate(norm_t(nch, nobs))

      
        read(unit) time_t
        read(unit) norm_t
        read(unit) p0
        read(unit) pexc
        read(unit) pion
        read(unit) dip
        read(unit) mom
        read(unit) nrg
      
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

      subroutine write_density_prob(filename, n, nch, jacc,          &
                                              xx, wx, rho)
      
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: n, nch
        real(8), intent(in) :: jacc
        real(8), intent(in) :: xx(n), wx(n)
        complex(8), intent(in) :: rho(n, nch)
      
        integer :: i, unit
      
        open(newunit=unit, file=filename, status='replace')

        do i=1,n
           write(unit,'(E20.10,*(1X,ES20.10))') xx(i),                 &
                                                rho(i,1), rho(i,2)
        enddo
      
        close(unit)
      
      end subroutine


      subroutine write_wf(filename, nch, n, x, t,         &
                                     omega, ich,                      &
                                     re_wf, im_wf,                    &
                                     rho, arg)
        implicit none
      
        character(*), intent(in) :: filename
        integer, intent(in) :: n, nch, ich
        real(8), intent(in) :: t
        real(8), intent(in) :: x(n), omega(nch)
        real(8), intent(in) :: re_wf(n,nch), im_wf(n,nch)
        real(8), intent(in) :: rho(n, nch), arg(n, nch)
      
        integer :: i, unit
      
        open(newunit=unit, file=filename, status='replace')
      
        ! --- header
        write(unit,*) "# t = ", t
        write(unit,*) "# omega = ", omega(ich)
        write(unit,*) "# channel = ", ich
        write(unit,*) "# x  Re  Im  |psi|^2  arg"
      
!       write(unit,*) omega
        do i = 1, n
           write(unit,'(E20.10,*(1X,ES20.10))') x(i),                  &
                                      re_wf(i,ich), im_wf(i,ich),      &
                                      rho(i,ich),   arg(i,ich)
        enddo
      
        close(unit)
      
      end subroutine


      subroutine write_pemd(filename, nch, n, kk, ak, bkwT, logscale)
        implicit none
        integer, intent(in) :: n, nch
        real(8), intent(in) :: kk(n)
        complex(8), intent(in) :: ak(n)
        complex(8), intent(in) :: bkwT(n,nch)
        logical, intent(in) :: logscale
        character(len=*), intent(in) :: filename
      
        integer :: i, unit
        real(8) :: pk
        real(8) :: pkwT(nch)
      
        open(newunit=unit, file=filename, status='replace')
      
        do i = 1, n
           pk = abs(ak(i))**2
           pkwT(:) = abs(bkwT(i,:))**2
      
           if (logscale) then
              if (pk > 1d-20) then
                 write(unit,*) kk(i), log10(pk)
              else
                 write(unit,*) kk(i), -20.d0
              end if
           else
              write(unit,'(E20.10,*(1X,ES20.10))') kk(i), pk, pkwT(:)
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
              write(unit,'(E20.10,1X,*(ES20.10))') omg(i), v1, v2
           end if
      
        end do
      
        close(unit)
      end subroutine

      subroutine write_bkw_b0w(filename, nch, omega, krange, kk, &
                                                      b0w, bkw)
      
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: krange, nch
        real(8), intent(in) :: kk(krange)
        real(8), intent(in) :: omega(nch)
        complex(8), intent(in) :: bkw(krange, nch)
        complex(8), intent(in) :: b0w(nch)
!       complex(8), intent(in) :: Qw(nch)
      
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


      subroutine write_Qw(filename, nch, omega, Qw)
      
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: nch
        real(8), intent(in) :: omega(nch)
        complex(8), intent(in) :: Qw(nch)
      
        integer :: j, unit
      
        open(newunit=unit, file=filename, status='replace')

        do j=1, nch
           write(unit,'(I8,1x,*(ES20.10))') j, omega(j), Qw(j)
        enddo
      
        close(unit)
      
      end subroutine


      end module io_module
