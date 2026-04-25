      module util
      implicit none
      contains
       subroutine dump_wavefunction(fname, nmax, eigvec, t, inv_jac,  &
                           wx, xx, psi)
       implicit none
       integer :: nmax, i
       real*8 :: t, inv_jac
       real*8 :: xx(nmax), wx(nmax)
       real*8 :: eigvec(nmax,nmax)
       complex*16 :: psi(nmax)
       complex*16 :: psi_x(nmax)
       character(*) :: fname
       integer unit
       
       open(newunit=unit,file=fname,status='replace')
      
       psi_x = matmul(eigvec,psi)
       psi_x = psi_x/wx/dsqrt(inv_jac)
       do i=1,nmax
          write(unit,'(5E20.12)') t, xx(i),                           &
               real(psi_x(i)), aimag(psi_x(i)), abs(psi_x(i))**2
       enddo
       
       close(unit)
       
       end subroutine dump_wavefunction

       function extract_name(path) result(name)
         implicit none
         character(*), intent(in) :: path
         character(len=256) :: name
         integer :: i
       
         name = trim(path)
       
         do i = len_trim(path), 1, -1
            if (path(i:i) == '/') then
               name = path(i+1:)
               exit
            end if
         end do
       
       end function



       subroutine init_run(workdir, mode, tag)
         implicit none
         character(len=*), intent(out) :: workdir
         character(len=*), intent(in)  :: mode
         character(len=*), intent(in), optional :: tag
       
         character(len=256) :: data_dir, cmd, fulltag
         character(len=32)  :: timestamp
         integer :: values(8), pid
       
         data_dir = "data/"
       
         call date_and_time(values=values)
       
         write(timestamp,'(i4.4,i2.2,i2.2,"_",i2.2,i2.2,i2.2)') &
              values(1), values(2), values(3), &
              values(5), values(6), values(7)
       
         pid = getpid()
       
         if (present(tag)) then
            fulltag = trim(tag)//"_"
         else
            fulltag = ""
         end if
       
         write(workdir,'(a,a,"_",a,a,"_",i0,"/")') &
              trim(data_dir), trim(mode), trim(fulltag),               &
                                             trim(timestamp), pid
       
         cmd = "mkdir -p " // trim(workdir)
         call execute_command_line(cmd)
       
         write(*,*) "Run directory:", trim(workdir)
       
       end subroutine



!      subroutine init_run(workdir, tag)
!        implicit none
!        character(len=*), intent(out) :: workdir
!        character(len=*), intent(in), optional :: tag
!      
!        character(len=256) :: data_dir, cmd, fulltag
!        character(len=32) :: timestamp
!        integer :: values(8), pid
!      
!        data_dir = "data/"
!      
!        call date_and_time(values=values)
!      
!        write(timestamp,'(i4.4,i2.2,i2.2,"_",i2.2,i2.2,i2.2)') &
!             values(1), values(2), values(3), &
!             values(5), values(6), values(7)
!      
!        pid = getpid()
!      
!        if (present(tag)) then
!           fulltag = trim(tag)//"_"
!        else
!           fulltag = ""
!        end if
!      
!        write(workdir,'(a,"dyn_",a,a,"_",i0,"/")') &
!             trim(data_dir), trim(fulltag), timestamp, pid
!      
!        cmd = "mkdir -p " // trim(workdir)
!        call execute_command_line(cmd)
!      
!        write(*,*) "Run directory:", trim(workdir)
!      
!      end subroutine


      end module util
