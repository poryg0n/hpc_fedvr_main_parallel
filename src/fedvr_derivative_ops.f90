      module fedvr_derivative_ops
        use constants, only : wp => dp
        implicit none
        private
        integer, public, parameter :: &
             FEDVR_LOCAL_DIRICHLET = 0, & 
             FEDVR_LEFT_OWNS       = 1, &
             FEDVR_RIGHT_OWNS      = 2, &
             FEDVR_AVERAGED        = 3

        public :: build_Dref_lobatto
        public :: eval_dpsi_fedvr
        public :: eval_dpsi_fedvr_reduced

      
      contains

        subroutine build_barycentric_weights(Np, x, lambda)
          implicit none
          integer, intent(in) :: Np
          real(wp), intent(in) :: x(Np)
          real(wp), intent(out) :: lambda(Np)
          integer :: i, k
        
          do i = 1, Np
             lambda(i) = 1.d0
             do k = 1, Np
                if (k /= i) lambda(i) = lambda(i) / (x(i) - x(k))
             end do
          end do
        end subroutine

        subroutine build_Dref_lobatto(Np, x, D)
          implicit none
          integer, intent(in) :: Np
          real(wp), intent(in) :: x(Np)
          real(wp), intent(out) :: D(Np,Np)
        
          real(wp) :: lambda(Np)
          integer :: i, j
        
          D = 0.d0
          call build_barycentric_weights(Np, x, lambda)
        
          do i = 1, Np
            do j = 1, Np
              if (i /= j) then
                D(i,j) = lambda(j)/lambda(i) / (x(i) - x(j))
              end if
            end do
            D(i,i) = -sum(D(i,:))
          end do
        end subroutine


        subroutine eval_dpsi_fedvr( lnbr, nnbr, xg, xs, map, Dref,    &
                                               psi, dpsi, iface_mode )
           implicit none
        
           integer, intent(in) :: lnbr, nnbr
           real(wp), intent(in) :: xg(:)          ! global grid
           real(wp), intent(in) :: xs(:)          ! element boundaries
           integer, intent(in) :: map(lnbr, nnbr)
           complex(wp), intent(in) :: psi(:)
           complex(wp), intent(out) :: dpsi(:)
           integer, intent(in) :: iface_mode     ! 0=A, 1=B, 2=AVG, 3=ZERO
        
           integer :: e, i, j
           integer :: ig, jg
           real(wp) :: invJ(lnbr)
           real(wp) :: Dref(nnbr, nnbr)
           complex(wp), allocatable :: dpsi_L(:), dpsi_R(:)
        
!          call build_Dref_lobatto(nnbr, xa, Dref)
        
           dpsi = (0.0d0, 0.0d0)
        
           !-------------------------------------------------
           ! Raw elementwise differentiation
           !-------------------------------------------------
           do e = 1, lnbr
              invJ (e) = 2.0d0 / ( xs(e+1) - xs(e) )
        
              do i = 1, nnbr
                 ig = map(e,i)
                 do j = 1, nnbr
                    jg = map(e,j)
                    dpsi(ig) = dpsi(ig) + invJ(e) * Dref(i,j) * psi(jg)
                 end do
              end do
           end do
        
           !-------------------------------------------------
           ! Interface correction (interior Lobatto nodes)
           !-------------------------------------------------
           do e = 1, lnbr-1
              ig = map(e,nnbr)   ! shared node (right of e, left of e+1)
        
              select case (iface_mode)

              case (0)   ! local Dirichlet
                 dpsi(ig) = (0.0d0, 0.0d0)
        
              case (1)   ! left owns
                 dpsi(ig) = invJ(e) * sum( Dref(nnbr,:) * psi(map(e,:)) )
        
              case (2)   ! right owns
                 dpsi(ig) = invJ(e+1) * sum( Dref(1,:) * psi(map(e+1,:)) )
        
              case (3)   ! average
                 dpsi(ig) = 0.5d0 * ( &
                    invJ(e)   * sum(Dref(nnbr,:) * psi(map(e,:))) + &
                    invJ(e+1) * sum(Dref(1,:)    * psi(map(e+1,:))) )
        
        
              end select
           end do
        
        end subroutine eval_dpsi_fedvr


        subroutine eval_dpsi_fedvr_reduced( lnbr, nnbr, xs, map, Dref, &
                                   psi_red, dpsi_red, workspace )
           implicit none
        
           ! ----------------------------
           ! Inputs
           ! ----------------------------
           integer, intent(in) :: lnbr, nnbr
           real(wp), intent(in) :: xs(lnbr+1)
           integer, intent(in) :: map(lnbr, nnbr)
           real(wp), intent(in) :: Dref(nnbr, nnbr)
           complex(wp), intent(in) :: psi_red(:)
        
           ! Workspace: (nfull, 2)
           complex(wp), intent(inout), target :: workspace(:,:)
        
           ! ----------------------------
           ! Output
           ! ----------------------------
           complex(wp), intent(out) :: dpsi_red(size(psi_red))
        
           ! ----------------------------
           ! Locals
           ! ----------------------------
           integer :: nfull, e, i, j, ig, jg, k
           real(wp) :: invJ
           complex(wp), pointer :: psi_full(:), dpsi_full(:)
        
           ! ----------------------------
           ! Sizes & safety
           ! ----------------------------
           nfull = lnbr*(nnbr-1) + 1
        
           if (size(workspace,1) < nfull .or. size(workspace,2) < 2) then
              stop "eval_dpsi_fedvr_reduced: workspace too small"
           end if
        
           psi_full  => workspace(1:nfull,1)
           dpsi_full => workspace(1:nfull,2)
        
           psi_full  = (0.d0,0.d0)
           dpsi_full = (0.d0,0.d0)
        
           ! ----------------------------
           ! Reduced → full injection
           ! ----------------------------
           k = 0
           do e = 1, lnbr
              do i = 2, nnbr-1
                 k = k + 1
                 psi_full(map(e,i)) = psi_red(k)
              end do
           end do
        
           ! ----------------------------
           ! Element-wise derivative
           ! ----------------------------
           do e = 1, lnbr
              invJ = 2.d0 / (xs(e+1) - xs(e))
        
              do i = 1, nnbr
                 ig = map(e,i)
                 do j = 1, nnbr
                    jg = map(e,j)
                    dpsi_full(ig) = dpsi_full(ig) + invJ * Dref(i,j) * psi_full(jg)
                 end do
              end do
           end do
        
           ! ----------------------------
           ! Full → reduced extraction
           ! ----------------------------
           k = 0
           do e = 1, lnbr
              do i = 2, nnbr-1
                 k = k + 1
                 dpsi_red(k) = dpsi_full(map(e,i))
              end do
           end do
        
        end subroutine eval_dpsi_fedvr_reduced



!       subroutine eval_dpsi_fedvr_reduced( lnbr, nnbr, xs, map, Dref, &
!                                          psi_red, dpsi_red )
!          implicit none
!       
!          ! ----------------------------
!          ! Inputs
!          ! ----------------------------
!          integer, intent(in) :: lnbr, nnbr
!          real(wp), intent(in) :: xs(lnbr+1)
!          integer, intent(in) :: map(lnbr, nnbr)
!          real(wp), intent(in) :: Dref(nnbr, nnbr)
!          complex(wp), intent(in) :: psi_red(:)
!       
!          ! ----------------------------
!          ! Output
!          ! ----------------------------
!          complex(wp), intent(out) :: dpsi_red(size(psi_red))
!       
!          ! ----------------------------
!          ! Locals
!          ! ----------------------------
!          integer :: nfull, e, i, j, ig, jg
!          real(wp) :: invJ
!          complex(wp), allocatable :: psi_full(:), dpsi_full(:)
!          integer :: k
!       
!          ! ----------------------------
!          ! Sizes
!          ! ----------------------------
!          nfull = lnbr*(nnbr-1) + 1
!       
!          allocate(psi_full(nfull))
!          allocate(dpsi_full(nfull))
!       
!          psi_full  = (0.d0,0.d0)
!          dpsi_full = (0.d0,0.d0)
!       
!          ! ----------------------------
!          ! Inject reduced → full
!          ! (endpoints = 0 Dirichlet)
!          ! ----------------------------
!          k = 0
!          do e = 1, lnbr
!             do i = 2, nnbr-1
!                k = k + 1
!                psi_full(map(e,i)) = psi_red(k)
!             end do
!          end do
!       
!          ! ----------------------------
!          ! Element-wise derivative
!          ! ----------------------------
!          do e = 1, lnbr
!       
!             invJ = 2.d0 / (xs(e+1) - xs(e))
!       
!             do i = 1, nnbr
!                ig = map(e,i)
!       
!                do j = 1, nnbr
!                   jg = map(e,j)
!                   dpsi_full(ig) = dpsi_full(ig) + invJ * Dref(i,j) * psi_full(jg)
!                end do
!       
!             end do
!          end do
!       
!          ! ----------------------------
!          ! Extract full → reduced
!          ! (interior surfaces already averaged)
!          ! ----------------------------
!          k = 0
!          do e = 1, lnbr
!             do i = 2, nnbr-1
!                k = k + 1
!                dpsi_red(k) = dpsi_full(map(e,i))
!             end do
!          end do
!       
!          deallocate(psi_full, dpsi_full)
!       
!       end subroutine eval_dpsi_fedvr_reduced


!       subroutine eval_dpsi_fedvr_generic(lnbr, nnbr, xg, xs, map,    &
!                                   Dref, psi, dpsi, workspace, bc_mode)
!           implicit none
!       
!           ! === Inputs ===
!           integer, intent(in) :: lnbr, nnbr
!           real(wp), intent(in) :: xg(:)       ! Global FEDVR points
!           real(wp), intent(in) :: xs(:)       ! Element boundaries
!           integer, intent(in) :: map(lnbr, nnbr)
!           real(wp), intent(in) :: Dref(nnbr, nnbr)
!           complex(wp), intent(in), target :: psi(:)
!           complex(wp), intent(out) :: dpsi(:)
!           complex(wp), intent(inout), target :: workspace(:) ! preallocated workspace
!           integer, intent(in) :: bc_mode       ! Boundary mode switch
!       
!           ! === Locals ===
!           integer :: e, i, j, ig, jg
!           integer :: nvec, nmax_
!           real(wp) :: invJ
!           logical :: is_reduced
!           complex(wp), pointer :: psi_full(:)
!       
!           ! -----------------------------
!           ! Determine if vector is reduced
!           nmax_ = size(xg)
!           nvec = size(psi)
!           is_reduced = (nvec < nmax_)
!       
!           ! -----------------------------
!           ! Map reduced vector to full if needed
!           if (is_reduced) then
!               ! Use workspace as psi_full to avoid allocation
!               if (size(workspace) /= nmax_) stop                     &
!                                    "Workspace size mismatch"
!               psi_full => workspace
!               psi_full = (0.d0,0.d0)
!       
!               ! Interior points: insert reduced psi
!               ! Global indexing: map(e,i)
!               do e = 1, lnbr
!                   do i = 2, nnbr-1
!                       ig = map(e,i)
!                       ! reduced vector index
!                       psi_full(ig) = psi(ig - e)  ! careful indexing
!                   end do
!               end do
!           else
!               psi_full => psi
!           end if
!       
!           ! -----------------------------
!           ! Zero derivative output
!           dpsi = (0.d0,0.d0)
!       
!           ! -----------------------------
!           ! Loop over elements
!           do e = 1, lnbr
!               invJ = 2.d0 / (xs(e+1) - xs(e))   ! dx/dξ
!       
!               do i = 1, nnbr
!                   ig = map(e,i)
!                   do j = 1, nnbr
!                       jg = map(e,j)
!                       dpsi(ig) = dpsi(ig) +                          &
!                                invJ * Dref(i,j) * psi_full(jg)
!                   end do
!               end do
!           end do
!       
!           ! -----------------------------
!           ! Apply boundary / surface closures
!           select case(bc_mode)
!           case(0)  ! DIRICHLET: zero derivative at endpoints
!               dpsi(1) = (0.d0,0.d0)
!               dpsi(nmax_) = (0.d0,0.d0)
!           case(1)  ! LEFT_OWNS: left-owned derivative at interior surfaces
!               ! Example: average left-owned at shared points
!               ! TODO: implement for multi-element surface handling
!           case(2)  ! AVERAGED: symmetrize interior surfaces
!               ! TODO: implement averaging
!           end select
!       
!           ! -----------------------------
!           ! Map back to reduced vector if needed
!           if (is_reduced) then
!               do e = 1, lnbr
!                   do i = 2, nnbr-1
!                       ig = map(e,i)
!                       psi_full(ig) = dpsi(ig)
!                   end do
!               end do
!               ! Copy interior-only derivative to output
!               dpsi = psi_full(2:nmax_-1)
!           end if
!       
!       end subroutine eval_dpsi_fedvr_generic


!       subroutine eval_dpsi_fedvr_reduced(lnbr, nnbr, map, Dref,      &
!                                                    psi_in, dpsi_out)
!           implicit none
!       
!           ! === Inputs ===
!           integer, intent(in) :: lnbr        ! number of elements
!           integer, intent(in) :: nnbr        ! number of points per element (Lobatto)
!           integer, intent(in) :: map(lnbr, nnbr)  ! FEDVR global mapping
!           real(wp), intent(in) :: Dref(nnbr, nnbr) ! local derivative matrix
!           complex(wp), intent(in) :: psi_in(:)    ! interior-only vector
!       
!           ! === Output ===
!           complex(wp), intent(out) :: dpsi_out(size(psi_in))  ! derivative
!       
!           ! === Locals ===
!           integer :: e, i, j, ig, jg
!           integer :: idx_left, idx_right
!           integer :: nint      ! total interior points
!           integer :: offset
!           complex(wp) :: contrib
!       
!           ! === Initialize ===
!           nint = size(psi_in)
!           dpsi_out = (0.0d0, 0.0d0)
!       
!           ! --- Loop over elements ---
!           offset = 0   ! running offset in psi_in
!           do e = 1, lnbr
!       
!               ! Loop over interior points of this element (exclude endpoints)
!               do i = 2, nnbr-1
!                   ig = map(e,i)
!       
!                   contrib = (0.0d0,0.0d0)
!                   do j = 2, nnbr-1
!                       jg = map(e,j)
!                       ! offset indexing: interior points are contiguous in psi_in
!                       contrib = contrib + Dref(i,j) * psi_in(offset + j-1)
!                   end do
!       
!                   ! accumulate to output (global interior index)
!                   dpsi_out(offset + i-1) = dpsi_out(offset + i-1) + contrib
!               end do
!       
!               ! Update offset: each element contributes (nnbr-2) interior points
!               offset = offset + (nnbr-2)
!           end do
!       
!           ! --- Simple surface handling ---
!           ! Average contributions at shared interior points between elements
!           do e = 1, lnbr-1
!               ! last interior point of element e = first interior point of element e+1
!               ! global interior indices
!               idx_left  = (e-1)*(nnbr-2) + (nnbr-2)  ! last interior point of element e
!               idx_right = e*(nnbr-2) + 1               ! first interior point of element e+1
!               ! average derivative at shared point
!               dpsi_out(idx_left) = 0.5d0*(dpsi_out(idx_left) + dpsi_out(idx_right))
!               dpsi_out(idx_right) = dpsi_out(idx_left)
!           end do
!       
!       end subroutine eval_dpsi_fedvr_reduced


        ! **** Below to correct

!       subroutine eval_dpsi_raw( lnbr, nnbr, xa, xs, map, psi, dpsi )
!          implicit none
!          integer, intent(in) :: lnbr, nnbr
!          real(wp), intent(in) :: xa(:), xs(:)
!          integer, intent(in) :: map(:)
!          complex(wp), intent(in) :: psi(:)
!          complex(wp), intent(out) :: dpsi(:)
!       
!          integer :: i, e, il, ig
!          real(wp) :: jac_inv
!       
!          dpsi = (0.0d0, 0.0d0)
!       
!          do e = 1, lnbr
!             jac_inv = 2.0d0 / ( xa(e+1) - xa(e) )
!       
!             do il = 1, nnbr
!                ig = map((e-1)*nnbr + il)
!       
!                dpsi(ig) = dpsi(ig) + jac_inv * sum(                  &
!                             Dref(il,1:nnbr) *                        &
!                            psi( map((e-1)*nnbr+1:(e-1)*nnbr+nnbr) ) )
!             end do
!          end do
!       
!       end subroutine eval_dpsi_raw

!       subroutine apply_boundary_closure( psi, dpsi, wa, mode )
!          implicit none
!          complex(wp), intent(in) :: psi(:)
!          complex(wp), intent(inout) :: dpsi(:)
!          real(wp), intent(in) :: wa(:)
!          integer, intent(in) :: mode
!       
!          integer :: n
!          complex(wp) :: corr
!          real(wp) :: norm
!       
!          n = size(psi)
!       
!          select case (mode)
!       
!          case (FEDVR_FULL_A)
!             dpsi(1) = -dpsi(n)
!             dpsi(n) =  dpsi(n)
!       
!          case (FEDVR_FULL_B)
!             corr = psi(n) - psi(1)
!             norm = sqrt(wa(1)) + sqrt(wa(n))
!             dpsi(1) = -corr / norm
!             dpsi(n) =  corr / norm
!       
!          case (FEDVR_FULL_C, FEDVR_FULL_DIRICHLET)
!             dpsi(1) = (0.0d0, 0.0d0)
!             dpsi(n) = (0.0d0, 0.0d0)
!       
!          case default
!             stop "Unknown FEDVR boundary mode"
!       
!          end select
!       
!       end subroutine apply_boundary_closure




!       subroutine eval_dpsi_fedvr_2( lnbr, nnbr, xa, xs, map, wa, &
!                                   psi, dpsi, mode )
!          implicit none
!          integer, intent(in) :: lnbr, nnbr
!          real(wp), intent(in) :: xa(:), xs(:), wa(:)
!          integer, intent(in) :: map(:)
!          complex(wp), intent(in) :: psi(:)
!          complex(wp), intent(out) :: dpsi(:)
!          integer, intent(in) :: mode
!       
!          integer :: n_full, n_red
!          complex(wp), allocatable :: psi_full(:), dpsi_full(:)
!       
!          n_full = size(xa)
!          n_red  = n_full - 2
!       
!          select case (mode)
!       
!          !===========================================
!          ! Dirichlet-reduced vector (INTERIOR ONLY)
!          !===========================================
!          case (FEDVR_REDUCED)
!       
!             if (size(psi) /= n_red) stop                             &
!                                     "FEDVR_REDUCED: wrong psi size"
!             if (size(dpsi) /= n_red) stop                            &
!                                     "FEDVR_REDUCED: wrong dpsi size"
!       
!             allocate(psi_full(n_full), dpsi_full(n_full))
!             psi_full = (0.0d0, 0.0d0)
!       
!             psi_full(2:n_full-1) = psi
!       
!             call eval_dpsi_raw( lnbr, nnbr, xa, xs, map,             &
!                                               psi_full, dpsi_full )
!       
!             dpsi = dpsi_full(2:n_full-1)
!       
!             deallocate(psi_full, dpsi_full)
!       
!          !===========================================
!          ! Full-grid vector
!          !===========================================
!          case default
!       
!             if (size(psi) /= n_full) stop                            &
!                                      "FEDVR_FULL: wrong psi size"
!             if (size(dpsi) /= n_full) stop                           &
!                                       "FEDVR_FULL: wrong dpsi size"
!       
!             call eval_dpsi_raw( lnbr, nnbr, xa, xs, map, psi, dpsi )
!       
!             call apply_boundary_closure( psi, dpsi, wa, mode )
!       
!          end select
!       
!       end subroutine eval_dpsi_fedvr_2


        ! **** Above to correct




        subroutine build_derivative_dvr_lambda(Np, x, w, D)
          implicit none
          integer, intent(in) :: Np
          real(wp), intent(in) :: x(Np), w(Np)
          real(wp), intent(out) :: D(Np, Np)

          real(wp) :: lambda(Np)
        
          integer :: i, j
        
          D = 0.d0
          call build_barycentric_weights(Np, x, lambda)
        
          ! Off-diagonal terms
          do j = 1, Np
             do i = 1, Np
                if (i /= j) then
                   D(j,i) = (lambda(i)/lambda(j)) &
                            / sqrt(w(i)*w(j))    &
                            / (x(j) - x(i))
                end if
             end do
          end do
        
          ! Diagonal terms (boundary condition encoded here)
          do i = 1, Np
             D(i,i) = -sum(D(i,1:Np)) + D(i,i)
          end do
        end subroutine
      
      
      end module fedvr_derivative_ops







