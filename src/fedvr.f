      MODULE FEDVR
      CONTAINS


      FUNCTION POTX(XX,IX,WX,NX,hscale)
      IMPLICIT REAL*8(A,B,D-H,O-Z)
      IMPLICIT COMPLEX*16(C)

      include 'param_structure'
      REAL*8 :: A
C======================================
C     Potential
C======================================
      
!     The free particle
!     POTX=0.d0

!      1/r Coulomb potential
c      POTX=-1.d0/Dsqrt(1.d0+XX**2)
c      POTX=-1.d0/Dsqrt(1.41d0+XX**2)

!     The zero-range potential
      IF (IX.eq.NX/2+1) THEN
        POTX=-kapp/WX**2/hscale
      ELSE
        POTX=0.d0
      ENDIF


!!     The harmonic oscillator
!      POTX=XX**2
!      POTX=POTX*0.5d0


      RETURN
      END


      SUBROUTINE kinetic_matrix_lobatto(Kinetic,Nelm,Nbas,XA,WA)
      IMPLICIT  NONE
!     EXTERNAL delta, integral
      INTEGER:: I, Ip, m, mp, j, k
      INTEGER NELM,NBAS

!     REAL*8 temp, delta, integral
      REAL*8 kinetic (1:(Nbas-1)*Nelm-1, 1:(Nbas-1)*Nelm-1)
      REAL*8 XA(*),WA(*)

      kinetic = 0.D0

      DO I =1, Nelm
        Ip=I
        DO m=1, Nbas-2
           DO mp=1, Nbas-2
              kinetic((Nbas-1)*(I-1)+m,(Nbas-1)*(Ip-1)+mp)
     &             =integral(m+1,mp+1,Nbas,XA,WA)
           ENDDO
        ENDDO
        m=Nbas-1
        DO mp=1, Nbas-2
           IF ((Nbas-1)*(I-1)+m<(Nbas-1)*Nelm
     &          .AND. (Nbas-1)*(Ip-1)+mp<(Nbas)*Nelm) THEN
              kinetic((Nbas-1)*(I-1)+m,(Nbas-1)*(Ip-1)+mp)
     &             =DSQRT(WA(Nbas)/(WA(1)+WA(Nbas)))
     &             *integral(m+1,mp+1,Nbas,XA,WA)
           ENDIF
        ENDDO

        mp=Nbas-1
        DO m=1, Nbas-2
           IF ((Nbas-1)*(I-1)+m<(Nbas-1)*Nelm
     &          .AND. (Nbas-1)*(Ip-1)+mp<(Nbas-1)*Nelm) THEN
              kinetic((Nbas-1)*(I-1)+m,(Nbas-1)*(Ip-1)+mp)
     &             =DSQRT(WA(Nbas)/(WA(1)+WA(Nbas)))
     &             *integral(m+1,mp+1,Nbas,XA,WA)
           ENDIF
        ENDDO
        m=Nbas-1
        mp=Nbas-1
        IF ((Nbas-1)*(I-1)+m<(Nbas-1)*Nelm
     &       .AND. (Nbas-1)*(Ip-1)+mp<(Nbas)*Nelm) THEN
             kinetic((Nbas-1)*(I-1)+m,(Nbas-1)*(Ip-1)+mp)
     &          =(WA(Nbas)/(WA(1)+WA(Nbas)))
     &          *integral(Nbas,Nbas,Nbas,XA,WA)
     &          +(WA(1)/(WA(1)+WA(Nbas)))
     &          *integral(1,1,Nbas,XA,WA)
          ENDIF
      ENDDO


      DO I =2, Nelm
         IP=I-1
            mp=Nbas-1
            DO m=1, Nbas-2
               IF ((Nbas-1)*(I-1)+m<(Nbas-1)*Nelm
     &              .AND. (Nbas-1)*(Ip-1)+mp<(Nbas)*Nelm) THEN
                  kinetic((Nbas-1)*(I-1)+m,(Nbas-1)*(Ip-1)+mp)
     &                 =DSQRT(WA(1)/(WA(1)+WA(Nbas)))
     &                 *integral(m+1,1,Nbas,XA,WA)
               ENDIF
            ENDDO
            m=Nbas-1
            IF ((Nbas-1)*(I-1)+m<(Nbas-1)*Nelm
     &           .AND. (Nbas-1)*(Ip-1)+mp<(Nbas)*Nelm) THEN
               kinetic((Nbas-1)*(I-1)+m,(Nbas-1)*(Ip-1)+mp)
     &              =DSQRT(WA(1)*WA(Nbas))
     &              /(WA(1)+WA(Nbas))
     &              *integral(1,Nbas,Nbas,XA,WA)
            ENDIF
      ENDDO


      DO I =1, Nelm-1
         Ip=I+1
            m=Nbas-1
            DO mp=1, Nbas-2
              IF ((Nbas-1)*(I-1)+m<(Nbas-1)*Nelm
     &        .AND. (Nbas-1)*(Ip-1)+mp<(Nbas)*Nelm) THEN

                kinetic((Nbas-1)*(I-1)+m,(Nbas-1)*(Ip-1)+mp)
     &          =DSQRT(WA(1)/(WA(1)+WA(Nbas)))
     &          *integral(1,mp+1,Nbas,XA,WA)
              ENDIF
            ENDDO
            mp=Nbas-1
            IF ((Nbas-1)*(I-1)+m<(Nbas-1)*Nelm
     &        .AND. (Nbas-1)*(Ip-1)+mp<(Nbas-1)*Nelm) THEN

              kinetic((Nbas-1)*(I-1)+m,(Nbas-1)*(Ip-1)+mp)
     &        =DSQRT(WA(Nbas)*WA(1))
     &        /(WA(1)+WA(Nbas))
     &        *integral(1,Nbas,Nbas,XA,WA)
             ENDIF

      ENDDO
      RETURN
      END



      FUNCTION product(i, mp, XA,Nbas)
      IMPLICIT NONE

      INTEGER i, mp, k
      INTEGER NBAS

      REAL*8  product, prod1
      REAL*8  XA(*)

      prod1 = 1.d0

      DO  k = 1, Nbas
          IF (k/=i .AND. k/=mp) THEN
            prod1 = prod1 * ((XA(i)-XA(k))
     &              /(XA(mp)-XA(k)))
          ENDIF
      ENDDO
      product = prod1/(XA(mp)-XA(i))
      RETURN
      END


      FUNCTION delta(i, j)

      INTEGER ::  i, j
      REAL*8  :: delta

      IF (i==j) THEN
         delta = 1.d0
      ELSE
         delta = 0.d0
      ENDIF
      RETURN
      END


      FUNCTION integral(m, mp,NBAS,XA,WA)
      IMPLICIT NONE
!     external delta, product
      INTEGER :: m, mp, j
      INTEGER Nbas
      REAL*8 :: integral
!     REAL*8 :: integral, delta, product
      REAL*8 WA(*),XA(*)

      integral =0.D0

      IF (m.eq.mp) THEN
         DO  j = 1, Nbas
            IF (j.ne.m) THEN
               integral = integral  +
     &              0.5d0 *WA(j) *product(j,m,XA,Nbas)**2
     &              / WA(m)
            ENDIF
         ENDDO
         integral = integral +
     &                        (delta(m,Nbas)-delta(m,1))**2
     &                        /(8.D0*WA(m)*WA(m))
      ELSE
         DO j=1, Nbas
            IF (j/=m .AND. j/=mp) THEN
               integral = integral +
     &              0.5d0*WA(j)*product(j,m,XA,Nbas)
     &              *product(j,mp,XA,Nbas)
     &              /DSQRT(WA(m)*WA(mp))

            ENDIF
         ENDDO
         integral = integral+0.25d0*(delta(m, Nbas)
     &      -delta(m,1))*product(m,mp,XA,Nbas)/dsqrt(WA(m)*
     &        WA(mp))+0.25d0*(delta(mp,Nbas)-delta(mp,1))*
     &       product(mp,m,XA,Nbas)/ DSQRT(WA(m)*WA(mp))
      ENDIF
      RETURN
      END FUNCTION integral

      function fact(n)
      implicit none
      integer:: j,n
      real(8) :: fact

      fact = 1
      if (n.gt.1) then
         do j=2,n
            fact=fact*j
         enddo
      end if

      return
      end function


!      SUBROUTINE lobatto_compute ( n, x, w )
!
!!*****************************************************************************80
!!
!!! LOBATTO_COMPUTE computes a Lobatto quadrature rule.
!!
!!  Discussion:
!!
!!    The integration interval is [ -1, 1 ].
!!
!!    The weight function is w(x) = 1.0.
!!
!!    The integral to approximate:
!!
!!      Integral ( -1 <= X <= 1 ) F(X) dX
!!
!!    The quadrature rule:
!!
!!      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) )
!!
!!    The quadrature rule will integrate exactly all polynomials up to
!!    X**(2*N-3).
!!
!!    The Lobatto rule is distinguished by the fact that both endpoints
!!    (-1 and 1) are always abscissas of the rule.
!!
!!  Licensing:
!!
!!    This code is distributed under the GNU LGPL license. 
!!
!!  Modified:
!!
!!    04 February 2007
!!
!!  Author:
!!
!!    Original MATLAB code by Greg von Winckel
!!    FORTRAN90 version by John Burkardt
!!
!!  Reference:
!!
!!    Milton Abramowitz, Irene Stegun,
!!    Handbook of Mathematical Functions,
!!    National Bureau of Standards, 1964,
!!    ISBN: 0-486-61272-4,
!!    LC: QA47.A34.
!!
!!    Claudio Canuto, Yousuff Hussaini, Alfio Quarteroni, Thomas Zang,
!!    Spectral Methods in Fluid Dynamics,
!!    Springer, 1993,
!!    ISNB13: 978-3540522058,
!!    LC: QA377.S676.
!!
!!    Arthur Stroud, Don Secrest,
!!    Gaussian Quadrature Formulas,
!!    Prentice Hall, 1966,
!!    LC: QA299.4G3S7.
!!
!!    Daniel Zwillinger, editor,
!!    CRC Standard Mathematical Tables and Formulae,
!!    30th Edition,
!!    CRC Press, 1996,
!!    ISBN: 0-8493-2479-3,
!!    LC: QA47.M315.
!!
!!  Parameters:
!!
!!    Input, INTEGER ( kind = 4 ) N, the order of the rule.  N must be 
!!    at least 2.
!!
!!    Output, real ( kind = 8 )  X(N), W(N), the abscissas and weights.
!!
!      implicit none
!
!      INTEGER ( kind = 4 ) n
!
!      INTEGER ( kind = 4 ) i
!      INTEGER ( kind = 4 ) j
!      real    ( kind = 8 ) p(n,n)
!      real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
!      real    ( kind = 8 ) tolerance
!      real    ( kind = 8 ) w(n)
!      real    ( kind = 8 ) x(n)
!      real    ( kind = 8 ) xold(n)
!
!      if ( n < 2 ) then
!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'LOBATTO_COMPUTE - Fatal error!'
!        write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', n
!        write ( *, '(a)' ) ' ORDER must be at least 2.'
!        stop
!      END if
!
!      tolerance = 100.0D+00 * epsilon ( tolerance )
!!
!!  Initial estimate for the abscissas is the Chebyshev-Gauss-Lobatto
!!  nodes.
!!
!      DO  i = 1, n
!        x(i) = cos ( pi * real ( i - 1, kind = 8 ) 
!     &       / real ( n - 1, kind = 8 ) )
!      END do
!
!      xold(1:n) = 2.0D+00
!
!      DO  while ( tolerance < maxval ( abs ( x(1:n) - xold(1:n) ) ) )
!
!        xold(1:n) = x(1:n)
!
!        p(1:n,1) = 1.0D+00
!        p(1:n,2) = x(1:n)
!        DO  j = 2, n-1
!          p(1:n,j+1) = ( real ( 2 * j - 1, kind = 8 ) 
!     &  * x(1:n) * p(1:n,j)     
!     & + real (   - j + 1, kind = 8 ) *  p(1:n,j-1) ) 
!     & / real (     j,     kind = 8 )
!        END do
!    
!        x(1:n) = xold(1:n) - ( x(1:n) * p(1:n,n) - p(1:n,n-1) ) 
!     &  / ( real ( n, kind = 8 ) * p(1:n,n) )
!     
!      END do
!
!      x(1:n) = x(n:1:-1)
!      w(1:n) = 2.0D+00 / ( real ( ( n - 1 ) * n, kind = 8 )
!     &    * p(1:n,n)**2 )
!
!      return
!      
!      end


      END




