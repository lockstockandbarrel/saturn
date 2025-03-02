!*==besy.f90 processed by SPAG 8.01RF 00:34  2 Mar 2025
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
!     ..................................................................
!
!        SUBROUTINE BESY
!
!        PURPOSE
!           COMPUTE THE Y BESSEL FUNCTION FOR A GIVEN ARGUMENT AND ORDER
!
!        USAGE
!           CALL BESY(X,N,BY,IER)
!
!        DESCRIPTION OF PARAMETERS
!           X  -THE ARGUMENT OF THE Y BESSEL FUNCTION DESIRED
!           N  -THE ORDER OF THE Y BESSEL FUNCTION DESIRED
!           BY -THE RESULTANT Y BESSEL FUNCTION
!           IER-RESULTANT ERROR CODE WHERE
!              IER=0  NO ERROR
!              IER=1  N IS NEGATIVE
!              IER=2  X IS NEGATIVE OR ZERO
!              IER=3  BY HAS EXCEEDED MAGNITUDE OF BIG; 
!                     where BIG=HUGE(0.0)
!
!        REMARKS
!           VERY SMALL VALUES OF X MAY CAUSE THE RANGE OF THE LIBRARY
!           FUNCTION ALOG TO BE EXCEEDED
!           X MUST BE GREATER THAN ZERO
!           N MUST BE GREATER THAN OR EQUAL TO ZERO
!
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
!           NONE
!
!        METHOD
!           RECURRENCE RELATION AND POLYNOMIAL APPROXIMATION TECHNIQUE
!           AS DESCRIBED BY A.J.M.HITCHCOCK,'POLYNOMIAL APPROXIMATIONS
!           TO BESSEL FUNCTIONS OF ORDER ZERO AND ONE AND TO RELATED
!           FUNCTIONS', M.T.A.C., V.11,1957,PP.86-88, AND G.N. WATSON,
!           'A TREATISE ON THE THEORY OF BESSEL FUNCTIONS', CAMBRIDGE
!           UNIVERSITY PRESS, 1958, P. 62
!
!     ..................................................................
!
SUBROUTINE besy(X,N,By,Ier)
IMPLICIT NONE

REAL a , b , By , c , fl , fl1 , p0 , p1 , pi2 , q0 , q1 , sum , t , t1 , t2 , term , ts , X , x2 , xx
REAL y0 , y1 , ya , yb , yc
INTEGER Ier , k , l , N
real,parameter :: big = huge(0.0) ! was 1.0E70; probably needs calculated
!
!     CHECK FOR ERRORS IN N AND X
!
   IF ( N<0 ) THEN
      Ier = 1
      RETURN
   ELSE
      Ier = 0
      IF ( X<=0 ) THEN
         Ier = 2
         RETURN
      ELSE
!
!     BRANCH IF X LESS THAN OR EQUAL 4
!
         IF ( X<=4.0 ) THEN
!
!       COMPUTE Y0 AND Y1 FOR X LESS THAN OR EQUAL TO 4
!
            xx = X/2.
            x2 = xx*xx
            t = alog(xx) + .5772157
            sum = 0.
            term = t
            y0 = t
            DO l = 1 , 15
               IF ( l/=1 ) sum = sum + 1./float(l-1)
               fl = l
               ts = t - sum
               term = (term*(-x2)/fl**2)*(1.-1./(fl*ts))
               y0 = y0 + term
            ENDDO
            term = xx*(t-.5)
            sum = 0.
            y1 = term
            DO l = 2 , 16
               sum = sum + 1./float(l-1)
               fl = l
               fl1 = fl - 1.
               ts = t - sum
               term = (term*(-x2)/(fl1*fl))*((ts-.5/fl)/(ts+.5/fl1))
               y1 = y1 + term
            ENDDO
            pi2 = .6366198
            y0 = pi2*y0
            y1 = -pi2/X + pi2*y1
         ELSE
!
!       COMPUTE Y0 AND Y1 FOR X GREATER THAN 4
!
            t1 = 4.0/X
            t2 = t1*t1
            p0 = ((((-.0000037043*t2+.0000173565)*t2-.0000487613)*t2+.00017343)*t2-.001753062)*t2 + .3989423
            q0 = ((((.0000032312*t2-.0000142078)*t2+.0000342468)*t2-.0000869791)*t2+.0004564324)*t2 - .01246694
            p1 = ((((.0000042414*t2-.0000200920)*t2+.0000580759)*t2-.000223203)*t2+.002921826)*t2 + .3989423
            q1 = ((((-.0000036594*t2+.00001622)*t2-.0000398708)*t2+.0001064741)*t2-.0006390400)*t2 + .03740084
            a = 2.0/sqrt(X)
            b = a*t1
            c = X - .7853982
            y0 = a*p0*sin(c) + b*q0*cos(c)
            y1 = -a*p1*cos(c) + b*q1*sin(c)
         ENDIF
!
!     CHECK IF ONLY Y0 OR Y1 IS DESIRED
!
         IF ( N>1 ) THEN
!
!    PERFORM RECURRENCE OPERATIONS TO FIND YN(X)
!
            ya = y0
            yb = y1
            k = 1
            SPAG_Loop_1_1: DO
               t = float(2*k)/X
               yc = t*yb - ya
               IF ( abs(yc)<=BIG ) THEN
                  k = k + 1
                  IF ( k/=N ) THEN
                     ya = yb
                     yb = yc
                  ELSE
                     By = yc
                     EXIT SPAG_Loop_1_1
                  ENDIF
               ELSE
                  Ier = 3
                  RETURN
               ENDIF
            ENDDO SPAG_Loop_1_1
!
!     RETURN EITHER Y0 OR Y1 AS REQUIRED
!
         ELSEIF ( N/=0 ) THEN
            By = y1
         ELSE
            By = y0
         ENDIF
      ENDIF
   ENDIF
   RETURN
END SUBROUTINE besy
