      module lobatto
      implicit none
      contains
      SUBROUTINE lobatto_compute ( n, x, w )

!*****************************************************************************80
!
!! LOBATTO_COMPUTE computes a Lobatto quadrature rule.
!
!  Discussion:
!
!    The integration interval is [ -1, 1 ].
!
!    The weight function is w(x) = 1.0.
!
!    The integral to approximate:
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    The quadrature rule:
!
!      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) )
!
!    The quadrature rule will integrate exactly all polynomials up to
!    X**(2*N-3).
!
!    The Lobatto rule is distinguished by the fact that both endpoints
!    (-1 and 1) are always abscissas of the rule.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 February 2007
!
!  Author:
!
!    Original MATLAB code by Greg von Winckel
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Claudio Canuto, Yousuff Hussaini, Alfio Quarteroni, Thomas Zang,
!    Spectral Methods in Fluid Dynamics,
!    Springer, 1993,
!    ISNB13: 978-3540522058,
!    LC: QA377.S676.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, INTEGER ( kind = 4 ) N, the order of the rule.  N must be 
!    at least 2.
!
!    Output, real ( kind = 8 )  X(N), W(N), the abscissas and weights.
!
      implicit none

      INTEGER ( kind = 4 ) n

      INTEGER ( kind = 4 ) i
      INTEGER ( kind = 4 ) j
      real    ( kind = 8 ) p(n,n)
      real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
      real    ( kind = 8 ) tolerance
      real    ( kind = 8 ) w(n)
      real    ( kind = 8 ) x(n)
      real    ( kind = 8 ) xold(n)

      if ( n < 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LOBATTO_COMPUTE - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', n
        write ( *, '(a)' ) ' ORDER must be at least 2.'
        stop
      END if

      tolerance = 1.d0
      tolerance = 100.0D+00 * epsilon ( tolerance )
!
!  Initial estimate for the abscissas is the Chebyshev-Gauss-Lobatto
!  nodes.
!
      DO  i = 1, n
        x(i) = cos ( pi * real ( i - 1, kind = 8 ) 
     &       / real ( n - 1, kind = 8 ) )
      END do

      xold(1:n) = 2.0D+00

      DO  while ( tolerance < maxval ( abs ( x(1:n) - xold(1:n) ) ) )

        xold(1:n) = x(1:n)

        p(1:n,1) = 1.0D+00
        p(1:n,2) = x(1:n)
        DO  j = 2, n-1
          p(1:n,j+1) = ( real ( 2 * j - 1, kind = 8 ) 
     &  * x(1:n) * p(1:n,j)     
     & + real (   - j + 1, kind = 8 ) *  p(1:n,j-1) ) 
     & / real (     j,     kind = 8 )
        END do
    
        x(1:n) = xold(1:n) - ( x(1:n) * p(1:n,n) - p(1:n,n-1) ) 
     &  / ( real ( n, kind = 8 ) * p(1:n,n) )
     
      END do

      x(1:n) = x(n:1:-1)
      w(1:n) = 2.0D+00 / ( real ( ( n - 1 ) * n, kind = 8 )
     &    * p(1:n,n)**2 )

      return
      
      end
      end 

