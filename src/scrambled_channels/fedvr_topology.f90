      module fedvr_topology
        implicit none
      contains

        subroutine build_fedvr_map(nelem, ng, map, nmax)
          implicit none
          integer, intent(in)  :: nelem, ng
          integer, intent(out) :: map(nelem, ng)
          integer, intent(out) :: nmax
     
          integer :: s, i
     
          ! global basis size
          nmax = nelem * (ng - 1) + 1
     
          do s = 1, nelem
            do i = 1, ng
              map(s,i) = (s - 1) * (ng - 1) + i
            end do
          end do
     
        end subroutine build_fedvr_map


        subroutine build_elements_bounds(Ne, xrange, xc, xs)
          implicit none
          integer, intent(in) :: Ne
          real(8), intent(in) :: xc
          real(8), intent(in) :: xrange
          real(8), intent(out) :: xs(1:Ne+1)
     
          real(8) :: dx, x_min
          integer :: s

          x_min   = xc - 0.5d0*xrange
          dx      = xrange / Ne
          do s = 1, Ne+1
              xs(s) = x_min + (s-1)*dx
          end do

        end subroutine


        subroutine build_x_global_fedvr(lnbr, nnbr, xa, xs, map, xg)
          implicit none
          integer, intent(in) :: lnbr, nnbr
          real(8), intent(in) :: xa(nnbr)
          real(8), intent(in) :: xs(lnbr+1)
          integer, intent(in) :: map(lnbr, nnbr)
          real(8), intent(out) :: xg(:)
        
          integer :: s, i, k
          real(8) :: xL, xR
        
          xg = 0.d0
        
          do s = 1, lnbr
             xL = xs(s)
             xR = xs(s+1)
        
             do i = 1, nnbr
                k = map(s,i)
                if (k > 0) then
                   xg(k) = 0.5d0 * ( (1.d0 - xa(i))*xL + (1.d0 + xa(i))*xR )

                end if
             end do
          end do

!         see also
!         xc = 0.5d0*(xL + xR)
!         h  = 0.5d0*(xR - xL)
!         xg(k) = xc + h * xa(i)

!         or
!         xg(k) = xL + 0.5d0*(1.d0 + xa(i))*(xR - xL)


        end subroutine



        subroutine build_Dloc_all(ng, xs, elem, Dref, Dloc)
          implicit none
          integer, intent(in) :: ng, elem
          real(8), intent(in) :: xs(:), Dref(ng,ng)
          real(8), intent(out):: Dloc(ng,ng)
     
          real(8) :: inv_jac 
         
          inv_jac = 2.d0 / (xs(elem+1) - xs(elem))  
          Dloc = inv_jac * Dref
        end subroutine


      
      end module fedvr_topology

