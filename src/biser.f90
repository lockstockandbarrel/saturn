!*==biser.f90 processed by SPAG 8.01RF 00:35  2 Mar 2025

!
!     ..................................................................
!
!        SUBROUTINE BISER
!
!        PURPOSE
!           TO COMPUTE THE BISERIAL CORRELATION COEFFICIENT BETWEEN TWO
!           CONTINUOUS VARIABLES WHEN ONE OF THEM HAS BEEN ARTIFICIALLY
!           DICHOTOMIZED.
!
!        USAGE
!           CALL BISER (N,A,B,HI,ANS,IER)
!
!        DESCRIPTION OF PARAMETERS
!           N   - NUMBER OF OBSERVATIONS
!           A   - INPUT VECTOR OF LENGTH N CONTAINING THE CONTINUOUS
!                 VARIABLE
!           B   - INPUT VECTOR OF LENGTH N CONTAINING THE DICHOTOMIZED
!                 VARIABLE
!           HI  - INPUT - NUMERICAL CODE TO INDICATE THE HIGHER CATEGORY
!                 OF THE DICHOTOMIZED VARIABLE.  ANY VALUE IN VECTOR B
!                 EQUAL TO OR GREATER THAN HI WILL BE CLASSIFIED INTO
!                 THE HIGHER CATEGORY.
!           ANS - OUTPUT VECTOR OF LENGTH 8 CONTAINING THE FOLLOWING
!                 ANS(1) - MEAN OF VARIABLE A
!                 ANS(2) - STANDARD DEVIATION OF VARIABLE A
!                 ANS(3) - PROPORTION OF THE CASES IN THE HIGHER
!                          CATEGORY OF VARIABLE B
!                 ANS(4) - PROPORTION OF THE CASES IN THE LOWER
!                          CATEGORY OF VARIABLE B
!                 ANS(5) - MEAN OF VARIABLE A FOR THOSE CASES FALLING
!                          INTO THE HIGHER CATEGORY OF VARIABLE B
!                 ANS(6) - MEAN OF VARIABLE A FOR THOSE CASES FALLING
!                          INTO THE LOWER CATEGORY OF VARIABLE B
!                 ANS(7) - BISERIAL CORRELATION COEFFICIENT
!                 ANS(8) - STANDARD ERROR OF BISERIAL CORRELATION
!                          COEFFICIENT
!           IER -  1, IF NO CASES ARE IN THE LOWER CATEGORY OF VARIABLE
!                 B.
!                 -1, IF ALL CASES ARE IN THE LOWER CATEGORY OF
!                 VARIABLE B.
!                 0, OTHERWISE.
!                 IF IER IS NON-ZERO, ANS(I)=10**75,I=5,...,8.
!
!        REMARKS
!           THE VALUES OF THE DICHOTOMIZED VARIABLE, B, MUST BE IN
!           NUMERIC FORM.  THEY CANNOR BE SPECIFIED BY MEANS OF
!           ALPHABETIC OR SPECIAL CHARACTERS.
!
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
!           NDTRI
!
!        METHOD
!           REFER TO P. HORST, 'PSYCHOLOGICAL MEASUREMENT AND
!           PREDICTION', P.95-96 (WADSWORTH, 1966).
!
!     ..................................................................
!
SUBROUTINE biser(N,A,B,Hi,Ans,Ier)
   IMPLICIT NONE

   REAL A,Ans,B,er,fn,Hi,p,q,r,sum,sum2,x,y
   INTEGER i,Ier,N

!
   DIMENSION A(*),B(*),Ans(*)
!
!        COMPUTE MEAN AND STANDARD DEVIATION OF VARIABLE A
!
   Ier = 0
   sum = 0.0
   sum2 = 0.0
   DO i = 1,N
      sum = sum + A(i)
      sum2 = sum2 + A(i)*A(i)
   ENDDO
   fn = N
   Ans(1) = sum/fn
   Ans(2) = (sum2-Ans(1)*sum)/(fn-1.0)
   Ans(2) = sqrt(Ans(2))
!
!        FIND PROPORTIONS OF CASES IN THE HIGHER AND LOWER CATEGORIES
!
   p = 0.0
   sum = 0.0
   sum2 = 0.0
   DO i = 1,N
      IF ( B(i)<Hi ) THEN
         sum2 = sum2 + A(i)
      ELSE
         p = p + 1.0
         sum = sum + A(i)
      ENDIF
   ENDDO
   Ans(4) = 1.0
   Ans(3) = 0.0
   q = fn - p
   IF ( p<=0 ) THEN
      Ier = -1
   ELSE
      Ans(5) = sum/p
      IF ( q<=0 ) THEN
         Ier = 1
         Ans(4) = 0.0
         Ans(3) = 1.0
      ELSE
         Ans(6) = sum2/q
         p = p/fn
         q = 1.0 - p
!
!        FIND ORDINATE OF THE NORMAL DISTRIBUTION CURVE AT THE POINT OF
!        DIVISION BETWEEN SEGMENTS CONTAINING P AND Q PROPORTIONS
!
         CALL ndtri(q,x,y,er)
!
!        COMPUTE THE BISERIAL COEFFICIENT OF CORRELATION
!
         r = ((Ans(5)-Ans(1))/Ans(2))*(p/y)
!
!        COMPUTE THE STANDARD ERROR OF R
!
         Ans(8) = (sqrt(p*q)/y-r*r)/sqrt(fn)
!
!        STORE RESULTS
!
         Ans(3) = p
         Ans(4) = q
         Ans(7) = r
         RETURN
      ENDIF
   ENDIF
   DO i = 5,8
      Ans(i) = huge(0.0)
   ENDDO
!
END SUBROUTINE biser
