!*==besk.f90 processed by SPAG 8.01RF 00:34  2 Mar 2025
!!SPAG Open source Personal, Educational or Academic User  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
!
!     ..................................................................
!
!        SUBROUTINE BESK
!
!           COMPUTE THE K BESSEL FUNCTION FOR A GIVEN ARGUMENT AND ORDER
!
!        USAGE
!           CALL BESK(X,N,BK,IER)
!
!        DESCRIPTION OF PARAMETERS
!           X  -THE ARGUMENT OF THE K BESSEL FUNCTION DESIRED
!           N  -THE ORDER OF THE K BESSEL FUNCTION DESIRED
!           BK -THE RESULTANT K BESSEL FUNCTION
!           IER-RESULTANT ERROR CODE WHERE
!              IER=0  NO ERROR
!              IER=1  N IS NEGATIVE
!              IER=2  X IS ZERO OR NEGATIVE
!              IER=3  X .GT. 170, MACHINE RANGE EXCEEDED
!              IER=4  BK .GT. BIG; where BIG=10**70
!
!        REMARKS
!           N MUST BE GREATER THAN OR EQUAL TO ZERO
!
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
!           NONE
!
!        METHOD
!           COMPUTES ZERO ORDER AND FIRST ORDER BESSEL FUNCTIONS USING
!           SERIES APPROXIMATIONS AND THEN COMPUTES N TH ORDER FUNCTION
!           USING RECURRENCE RELATION.
!           RECURRENCE RELATION AND POLYNOMIAL APPROXIMATION TECHNIQUE
!           AS DESCRIBED BY A.J.M.HITCHCOCK,'POLYNOMIAL APPROXIMATIONS
!           TO BESSEL FUNCTIONS OF ORDER ZERO AND ONE AND TO RELATED
!           FUNCTIONS', M.T.A.C., V.11,1957,PP.86-88, AND G.N. WATSON,
!           'A TREATISE ON THE THEORY OF BESSEL FUNCTIONS', CAMBRIDGE
!           UNIVERSITY PRESS, 1958, P. 62
!
!     ..................................................................
!
SUBROUTINE besk(X,N,Bk,Ier)
IMPLICIT NONE
REAL a , b , Bk , c , fact , g0 , g1 , gj , hj , rj , t , X , x2j
INTEGER Ier , j , l , N
DIMENSION t(12)
real,parameter :: big = huge(0.0) ! was 1.0E70; probably needs calculated
   Bk = .0
   IF ( N<0 ) THEN
      Ier = 1
      RETURN
   ELSEIF ( X<=0 ) THEN
      Ier = 2
      RETURN
   ELSE
      DO WHILE ( X<=170.0 )
         Ier = 0
         IF ( X<=1. ) THEN
            b = X/2.
            a = .5772157 + alog(b)
            c = b*b
            IF ( N==1 ) THEN
               CALL spag_block_2
               RETURN
            ENDIF
!
!     COMPUTE KO USING SERIES EXPANSION
!
            g0 = -a
            x2j = 1.
            fact = 1.
            hj = .0
            DO j = 1 , 6
               rj = 1./float(j)
               x2j = x2j*c
               fact = fact*rj*rj
               hj = hj + rj
               g0 = g0 + x2j*fact*(hj-a)
            ENDDO
            IF ( N/=0 ) THEN
               CALL spag_block_2
               RETURN
            ENDIF
            Bk = g0
            RETURN
         ELSE
            a = exp(-X)
            b = 1./X
            c = sqrt(b)
            t(1) = b
            DO l = 2 , 12
               t(l) = t(l-1)*b
            ENDDO
            IF ( N/=1 ) THEN
!
!     COMPUTE KO USING POLYNOMIAL APPROXIMATION
!
               g0 = a*(1.2533141-.1566642*t(1)+.08811128*t(2)-.09139095*t(3)+.1344596*t(4)-.2299850*t(5)+.3792410*t(6)-.5247277*t(7)&
                  & +.5575368*t(8)-.4262633*t(9)+.2184518*t(10)-.06680977*t(11)+.009189383*t(12))*c
               IF ( N<0 ) CYCLE
               IF ( N==0 ) THEN
                  Bk = g0
                  RETURN
               ENDIF
            ENDIF
!
!     COMPUTE K1 USING POLYNOMIAL APPROXIMATION
!
            g1 = a*(1.2533141+.4699927*t(1)-.1468583*t(2)+.1280427*t(3)-.1736432*t(4)+.2847618*t(5)-.4594342*t(6)+.6283381*t(7)     &
               & -.6632295*t(8)+.5050239*t(9)-.2581304*t(10)+.07880001*t(11)-.01082418*t(12))*c
            IF ( N<1 ) THEN
            ELSEIF ( N==1 ) THEN
               Bk = g1
               RETURN
            ELSE
               CALL spag_block_1
               RETURN
            ENDIF
         ENDIF
      ENDDO
      Ier = 3
      RETURN
   ENDIF
CONTAINS
   SUBROUTINE spag_block_1
!
!     FROM KO,K1 COMPUTE KN USING RECURRENCE RELATION
!
      SPAG_Loop_1_1: DO j = 2 , N
         gj = 2.*(float(j)-1.)*g1/X + g0
         IF ( gj<= BIG ) THEN
            g0 = g1
            g1 = gj
         ELSE
            Ier = 4
            EXIT SPAG_Loop_1_1
         ENDIF
      ENDDO SPAG_Loop_1_1
      Bk = gj
      RETURN
   END SUBROUTINE spag_block_1
   SUBROUTINE spag_block_2
!
!     COMPUTE K1 USING SERIES EXPANSION
!
      x2j = b
      fact = 1.
      hj = 1.
      g1 = 1./X + x2j*(.5+a-hj)
      DO j = 2 , 8
         x2j = x2j*c
         rj = 1./float(j)
         fact = fact*rj*rj
         hj = hj + rj
         g1 = g1 + x2j*fact*(.5+(a-hj)*float(j))
      ENDDO
      IF ( N/=1 ) THEN
         CALL spag_block_1
         RETURN
      ENDIF
      Bk = g1
   END SUBROUTINE spag_block_2
END SUBROUTINE besk
