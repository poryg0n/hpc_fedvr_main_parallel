      module space_time_ops
       implicit none
      contains
         subroutine eigen_to_dvr(nmax, nchan, jac, wx, eigvec, wf, wfc)
           implicit none
           integer, intent(in) :: nmax, nchan
           real(8), intent(in) :: jac
           real(8), intent(in) :: wx(nmax)
           real(8), intent(in) :: eigvec(nmax,nmax)
           complex(8), intent(in) :: wf(nmax,nchan+1)
           complex(8), intent(out) :: wfc(nmax,nchan+1)

           complex(8):: tmp(nmax,nchan+1)
           integer :: k
         
!          wfc = matmul(eigvec, wf)

           do k=1,nchan+1
              tmp(:,k) = matmul(eigvec, wf(:,k))
           enddo

           do k=1,nchan+1
              wfc(:,k) = tmp(:,k) / wx / dsqrt(jac)
           enddo
         
         end subroutine


         subroutine dvr_to_eigen(nmax, nchan, jac, wx, eigvec, wfc, wf)
           implicit none
           integer, intent(in) :: nmax, nchan
           real(8), intent(in) :: jac
           real(8), intent(in) :: wx(nmax)
           real(8), intent(in) :: eigvec(nmax,nmax)
           complex(8), intent(in) :: wfc(nmax,nchan+1)
           complex(8), intent(out) :: wf(nmax,nchan+1)
         
           complex(8):: tmp(nmax,nchan+1)
           integer :: k


           do k=1,nchan+1
              tmp(:,k) = wfc(:,k) * wx * dsqrt(jac)
           enddo

           do k=1,nchan+1
              wf(:,k) = matmul(transpose(eigvec), tmp(:,k))
           enddo

!          pause

!          call zgemm('t', 'n', nmax, nchan+2, nmax,                   &
!         (1.d0,0.d0), eigvec, nmax, wfc, nmax, (0.d0,0.d0), wf, nmax)

         
         end subroutine



         subroutine eigen_to_dvr_1d(nmax, jac, wx, eigvec, psi, psi_dvr)
           implicit none
           integer, intent(in) :: nmax
           real(8), intent(in) :: jac
           real(8), intent(in) :: wx(nmax)
           real(8), intent(in) :: eigvec(nmax,nmax)
           complex(8), intent(in) :: psi(nmax)
           complex(8), intent(out) :: psi_dvr(nmax)

           psi_dvr = matmul(eigvec, psi)
           psi_dvr = psi_dvr/ wx / dsqrt(jac)

         end subroutine

         subroutine dvr_to_eigen_1d(nmax, jac, wx, eigvec, psi_dvr, psi)
           implicit none
           integer, intent(in) :: nmax
           real(8), intent(in) :: jac
           real(8), intent(in) :: wx(nmax)
           real(8), intent(in) :: eigvec(nmax,nmax)
           complex(8), intent(in) :: psi_dvr(nmax)
           complex(8), intent(out) :: psi(nmax)

           complex(8) :: tmp(nmax)

           tmp = psi_dvr * wx * dsqrt(jac)
           psi = matmul(transpose(eigvec), tmp)

         end subroutine


      end module



