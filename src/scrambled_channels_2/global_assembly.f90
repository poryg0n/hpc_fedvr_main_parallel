      module global_assembly
        implicit none
        integer, parameter :: FEDVR_SYMMETRIC = 1
        integer, parameter :: FEDVR_LEFT      = 2
        integer, parameter :: FEDVR_MUTUAL    = 3

      contains
      
!       subroutine assemble_derivative(nelem, ng, map, Dloc_all, Dglobal)
!         implicit none
!         integer, intent(in) :: nelem, ng
!         integer, intent(in) :: map(nelem,ng)
!         real(8), intent(in) :: Dloc_all(nelem,ng,ng)
!         real(8), intent(out) :: Dglobal(:,:)
!     
!         integer :: elem, i, j, gI, gJ
!     
!         do elem = 1, nelem
!            do i = 1, ng
!               gI = map(elem,i)
!               do j = 1, ng
!                 gJ = map(elem,j)
!                 Dglobal(gI,gJ) = Dglobal(gI,gJ) + Dloc_all(elem,i,j)
!               end do
!            end do
!         end do
!       end subroutine

        subroutine assemble_interface_row( s, nelem, ng, map, w, &
                                   Dloc_all, Dglobal, mode )
           implicit none
           integer, intent(in) :: s, nelem, ng, mode
           integer, intent(in) :: map(nelem,ng)
           real(8), intent(in) :: w(ng)
           real(8), intent(in) :: Dloc_all(nelem,ng,ng)
           real(8), intent(inout) :: Dglobal(:,:)
         
           integer :: j, gI, gJ
           real(8) :: Omega, facL, facR
         
           ! normalization
           Omega = w(ng) + w(1)
           facL  = sqrt(w(ng)) / sqrt(Omega)
           facR  = sqrt(w(1))  / sqrt(Omega)
         
           select case (mode)
           case (FEDVR_SYMMETRIC)
             facL = 0.5d0 * facL
             facR = 0.5d0 * facR
           case (FEDVR_LEFT)
             facR = 0.d0
           case (FEDVR_MUTUAL)
             ! unchanged
           end select
         
           gI = map(s,ng)
         
           ! overwrite interface row (owned by left element s)
           do j = 1, ng
             gJ = map(s,j)
             Dglobal(gI,gJ) = facL * Dloc_all(s,ng,j)
           end do
         
           if (facR /= 0.d0 .and. s < nelem) then
             do j = 1, ng
               gJ = map(s+1,j)
               Dglobal(gI,gJ) = Dglobal(gI,gJ) + facR * Dloc_all(s+1,1,j)
             end do
           end if
         
         end subroutine assemble_interface_row


         subroutine assemble_derivative_fedvr( nelem, ng, map, w, &
                                               Dloc_all, Dglobal, mode )
           implicit none
           integer, intent(in) :: nelem, ng, mode
           integer, intent(in) :: map(nelem,ng)
           real(8), intent(in) :: w(ng)
           real(8), intent(in) :: Dloc_all(nelem,ng,ng)
           real(8), intent(out) :: Dglobal(:,:)
         
           integer :: s, i, j, gI, gJ
         
           Dglobal = 0.d0
         
           do s = 1, nelem
             do i = 1, ng
         
               gI = map(s,i)
         
               if (i /= ng) then
                 ! interior DOF (owned by this element)
                 do j = 1, ng
                   gJ = map(s,j)
                   Dglobal(gI,gJ) = Dglobal(gI,gJ) + Dloc_all(s,i,j)
                 end do
         
               else
                 ! interface DOF: assemble ONCE
                 if (s < nelem) then
                    call assemble_interface_row( s, nelem, ng, map, w, &
                                              Dloc_all, Dglobal, mode )
                 end if
               end if
         
             end do
           end do
         
         end subroutine

      
      end module

