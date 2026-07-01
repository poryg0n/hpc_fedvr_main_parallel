      module space_time_ops
      implicit none
      contains
         subroutine eigen_to_dvr(nmax, jac, wx, eigvec, psi, psi_dvr)
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

         subroutine dvr_to_eigen(nmax, jac, wx, eigvec, psi_dvr, psi)
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
