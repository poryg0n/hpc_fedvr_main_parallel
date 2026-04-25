      module fedvr_conf_struct
      use fedvr
      use structure_parameters
!     use fedvr_basis
      implicit none
      contains
      subroutine fedvr_hamilton_conf(jac_,xa,wa,xs,xx,wx,eigval,     &
                         eigenvec)
       integer :: ntot, nout,                             &
                  mdim2, mdim2s, mfnd, lwork2, liwork2, info,        &
                  vl, vu, il, iu, &
                  i, j, k

       real*8 :: jac_, abstol, start, lap, xrange_
       real*8, dimension(:) :: xa, wa, xx, xs, wx, eigval
       real*8, dimension(:,:) :: eigenvec

       real*8, allocatable, dimension(:) :: vec_matup,               &
                                            work2, iwork2, ifail2
  
       real*8, allocatable, dimension(:) :: isuppz
       real*8, allocatable, dimension(:,:) :: ham0, id

!      include 'param_structure'
 
 
       nout = 200000
!      nnbr = size(xa)
!      snbr = size(xs)-1
       ntot = snbr *(nnbr-1)-1

       write(*,*) "nnbr", nnbr
       write(*,*) "snbr", snbr
       write(*,*) "xrange", xrange
       if (ntot.ne.(size(xx))) stop

       mdim2 = nmax

       lwork2 = 26*mdim2
       liwork2 = 10*mdim2
       mdim2s = mdim2*(mdim2+1)/2

       if (lwork2.lt.0.or.liwork2.lt.0) then
          write(*,*) "Problem with the integer lwork2 and liwork2"
          write(*,*) lwork2, liwork2
          call abort
       endif


       write(*,*) "Construction of kinetic energy matrix" 
       call cpu_time(start)
       call kinetic_matrix_lobatto(eigenvec,snbr ,nnbr,xa,wa)
!      call basis_1d(basis,snbr,nnbr,xa,wa)
       call cpu_time(lap)
       write(*,*) "Time required for kinetic energy matrix",         &
                    lap-start 

!         write(*,*) "xa", xa
!         write(*,*) "wa", wa
          write(*,*) "xrange", xrange_
          write(*,*) "jac = 0.5d0*xrange/nelem", jac_
      
          write(111,*) 
           write(*,*) "before jac"
           do i=1,ntot
              write(111,*) eigenvec(i,:)
           enddo
     
       write(*,*) "kinetic matrix built"
       eigenvec=eigenvec/jac_**2

       do i=1,snbr 
          do j=1,nnbr-1
             if((i.ne.snbr).or.(j.ne.(nnbr-1))) then
                k = (nnbr-1)*(i-1)+j
                xx(k) = xa(j+1) + (i-.5d0-.5d0*snbr)*2
                xx(k) = xx(k)*jac_
                if(j.eq.(nnbr-1)) then
                   wx(k) = dsqrt(wa(nnbr)+wa(1))
                else
                   wx(k) = dsqrt(wa(j+1))
                end if 
                eigenvec(k,k) = eigenvec(k,k) +                       &
!                      potx(xx(k),k,wx(ntot/2+1),ntot,jac_)
                       potx(xx(k),k,wx(k),ntot,jac_)
!               write(909,*) k,potx(xx(k),k,wx(ntot/2+1),ntot,jac_)
!               write(909,*) k,potx(xx(k),k,wx(k),ntot,jac_)
             end if 
          enddo
       enddo


       abstol = 1.0d-30
       if (nout.gt.ntot) nout=ntot


       write(*,*) "diagonalization", ntot 
       call cpu_time(start)

      
!      allocate(vec_matup(mdim2s))
       allocate(work2(lwork2),iwork2(liwork2),ifail2(mdim2))
       allocate(isuppz(2*mdim2))

!      k=0
!      do i=1,ntot
!         do j=1,i
!            k=k+1
!            vec_matup(k)=eigenvec(j,i)
!         enddo
!      enddo
       allocate(ham0(ntot,ntot))
       ham0=eigenvec

!      allocate(id(ntot,ntot))
!      id=0.0d0
!      do i=1,ntot
!         id(i,i)=1.0d0
!      enddo

!      call dspevx('v', 'i', 'u', ntot, vec_matup, vl, vu, 1, nout,   &
!               abstol, mfnd, eigval, eigenvec, ntot, &
!               work2, iwork2, ifail2, info) 

!      call dsygv(1,'v','u',ntot,eigenvec,ntot,id,ntot,eigval,        &
!               work2,lwork2,info)

!      call dsyev('v','u',ntot,eigenvec,ntot,eigval,work2,lwork2,info)

!      call dsyevd('v','u',ntot,eigenvec,ntot,eigval,                 &
!                work2,lwork2,iwork2,liwork2,info)

       call dsyevr('v','a','u',ntot,ham0,ntot,vl,vu,il,iu,            &
          abstol,mfnd,eigval,eigenvec,ntot,                           &
          isuppz,work2,lwork2,iwork2,liwork2,                         &
          info)        

!      deallocate(id)
       deallocate(ham0)
      

       write(*,*) info
       if(info /= 0) then
           write(404,*) "dsyevr info", info
           call flush(404)
           call abort
       end if
       call cpu_time(lap)
       write(*,*) "Time required for diagonalization (dsyevr) ",    &
               lap-start 


       return
       end subroutine
      end 

