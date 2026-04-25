      module space_time_ops
       implicit none
      contains
         subroutine eigen_to_dvr(nmax, nchan, jac, wx, eigvec, wf, wfc)
           implicit none
           integer, intent(in) :: nmax, nchan
           real(8), intent(in) :: jac
           real(8), intent(in) :: wx(nmax)
           real(8), intent(in) :: eigvec(nmax,nmax)
           complex(8), intent(in) :: wf(nmax,nchan+2)
           complex(8), intent(out) :: wfc(nmax,nchan+2)

           integer :: k
         
           wfc = matmul(eigvec, wf)
           do k=1,nchan+2
              wfc(:,k) = wfc(:,k) / wx / dsqrt(jac)
           enddo
         
         end subroutine


         subroutine dvr_to_eigen(nmax, nchan, jac, wx, eigvec, wfc, wf)
           implicit none
           integer, intent(in) :: nmax, nchan
           real(8), intent(in) :: jac
           real(8), intent(in) :: wx(nmax)
           real(8), intent(in) :: eigvec(nmax,nmax)
           complex(8), intent(in) :: wfc(nmax,nchan+2)
           complex(8), intent(out) :: wf(nmax,nchan+2)
         
           integer :: k

           wf = matmul(transpose(eigvec), wfc)
!          call zgemm('t', 'n', nmax, nchan+2, nmax,                   &
!         (1.d0,0.d0), eigvec, nmax, wfc, nmax, (0.d0,0.d0), wf, nmax)

           do k=1,nchan+2
              wf(:,k) = wf(:,k) * wx * dsqrt(jac)
           enddo
         
         end subroutine

      end module



