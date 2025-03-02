!*==canor.f90 processed by SPAG 8.01RF 00:35  2 Mar 2025
!..................................................................
!   SUBROUTINE CANOR
!
!   PURPOSE
!      COMPUTE THE CANONICAL CORRELATIONS BETWEEN TWO SETS OF
!      VARIABLES.  CANOR IS NORMALLY PRECEDED BY A CALL TO SUBROU-
!      TINE CORRE.
!
!   USAGE
!      CALL CANOR (N,MP,MQ,RR,ROOTS,WLAM,CANR,CHISQ,NDF,COEFR,
!                  COEFL,R)
!
!   DESCRIPTION OF PARAMETERS
!      N     - NUMBER OF OBSERVATIONS
!      MP    - NUMBER OF LEFT HAND VARIABLES
!      MQ    - NUMBER OF RIGHT HAND VARIABLES
!      RR    - INPUT MATRIX (ONLY UPPER TRIANGULAR PORTION OF THE
!              SYMMETRIC MATRIX OF M X M, WHERE M = MP + MQ)
!              CONTAINING CORRELATION COEFFICIENTS.  (STORAGE MODE
!              OF 1)
!      ROOTS - OUTPUT VECTOR OF LENGTH MQ CONTAINING EIGENVALUES
!              COMPUTED IN THE NROOT SUBROUTINE.
!      WLAM  - OUTPUT VECTOR OF LENGTH MQ CONTAINING LAMBDA.
!      CANR  - OUTPUT VECTOR OF LENGTH MQ CONTAINING CANONICAL
!              CORRELATIONS.
!      CHISQ - OUTPUT VECTOR OF LENGTH MQ CONTAINING THE
!              VALUES OF CHI-SQUARES.
!      NDF   - OUTPUT VECTOR OF LENGTH MQ CONTAINING THE DEGREES
!              OF FREEDOM ASSOCIATED WITH CHI-SQUARES.
!      COEFR - OUTPUT MATRIX (MQ X MQ) CONTAINING MQ SETS OF
!              RIGHT HAND COEFFICIENTS COLUMNWISE.
!      COEFL - OUTPUT MATRIX (MP X MQ) CONTAINING MQ SETS OF
!              LEFT HAND COEFFICIENTS COLUMNWISE.
!      R     - WORK MATRIX (M X M)
!
!   REMARKS
!      THE NUMBER OF LEFT HAND VARIABLES (MP) SHOULD BE GREATER
!      THAN OR EQUAL TO THE NUMBER OF RIGHT HAND VARIABLES (MQ).
!      THE VALUES OF CANONICAL CORRELATION, LAMBDA, CHI-SQUARE,
!      DEGREES OF FREEDOM, AND CANONICAL COEFFICIENTS ARE COMPUTED
!      ONLY FOR THOSE EIGENVALUES IN ROOTS WHICH ARE GREATER THAN
!      ZERO.
!
!   SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
!      MINV
!      NROOT  (WHICH, IN TURN, CALLS THE SUBROUTINE EIGEN.)
!
!   METHOD
!      REFER TO W. W. COOLEY AND P. R. LOHNES, 'MULTIVARIATE PRO-
!      CEDURES FOR THE BEHAVIORAL SCIENCES', JOHN WILEY AND SONS,
!      1962, CHAPTER 3.
!..................................................................
SUBROUTINE canor(N,Mp,Mq,Rr,Roots,Wlam,Canr,Chisq,Ndf,Coefr,Coefl,R)
IMPLICIT NONE

REAL Canr,Chisq,Coefl,Coefr,det,fmp,fmq,fn,R,Roots,Rr,sum,Wlam
INTEGER i,j,jj,k,l,m,Mp,Mq,N,n1,n2,n3,Ndf
integer,allocatable :: iwork(:)

DIMENSION Rr(*),Roots(*),Wlam(*),Canr(*),Chisq(*),Ndf(*),Coefr(*),Coefl(*),R(*)
!...............................................................
!
!   IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE
!   C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION
!   STATEMENT WHICH FOLLOWS.
!
! DOUBLE PRECISION RR,ROOTS,WLAM,CANR,CHISQ,COEFR,COEFL,R,DET,SUM
!
!   THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS
!   APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS
!   ROUTINE.
!
!   THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO
!   CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT IN STATEMENT
!   165 MUST BE CHANGED TO DSQRT.  ALOG IN STATEMENT 175 MUST BE
!   CHANGED TO DLOG.
!
!   ...............................................................
!
! PARTITION INTERCORRELATIONS AMONG LEFT HAND VARIABLES, BETWEEN
! LEFT AND RIGHT HAND VARIABLES, AND AMONG RIGHT HAND VARIABLES.
!
   m = Mp + Mq
   n1 = 0
   DO i = 1,m
      DO j = 1,m
         IF ( i<j ) THEN
            l = i + (j*j-j)/2
         ELSE
            l = j + (i*i-i)/2
         ENDIF
         n1 = n1 + 1
         R(n1) = Rr(l)
      ENDDO
   ENDDO
   l = Mp
   DO j = 2,Mp
      n1 = m*(j-1)
      DO i = 1,Mp
         l = l + 1
         n1 = n1 + 1
         R(l) = R(n1)
      ENDDO
   ENDDO
   n2 = Mp + 1
   l = 0
   DO j = n2,m
      n1 = m*(j-1)
      DO i = 1,Mp
         l = l + 1
         n1 = n1 + 1
         Coefl(l) = R(n1)
      ENDDO
   ENDDO
   l = 0
   DO j = n2,m
      n1 = m*(j-1) + Mp
      DO i = n2,m
         l = l + 1
         n1 = n1 + 1
         Coefr(l) = R(n1)
      ENDDO
   ENDDO
!
!     SOLVE THE CANONICAL EQUATION
!
   l = Mp*Mp + 1
   k = l + Mp
   if(allocated(iwork))deallocate(iwork)
   ! SYMMETRIC MATRIX OF M X M, WHERE M = MP + MQ)
   allocate( iwork((mp+mq)*(mp+mq)) )
   CALL minv(R,Mp,det,iwork(l:l+Mp-1),iwork(k:k+Mp-1))
!
!        CALCULATE T = INVERSE OF R11 * R12
!
   DO i = 1,Mp
      n2 = 0
      DO j = 1,Mq
         n1 = i - Mp
         Roots(j) = 0.0
         DO k = 1,Mp
            n1 = n1 + Mp
            n2 = n2 + 1
            Roots(j) = Roots(j) + R(n1)*Coefl(n2)
         ENDDO
      ENDDO
      l = i - Mp
      DO j = 1,Mq
         l = l + Mp
         R(l) = Roots(j)
      ENDDO
   ENDDO
!
!        CALCULATE A = R21 * T
!
   l = Mp*Mq
   n3 = l + 1
   DO j = 1,Mq
      n1 = 0
      DO i = 1,Mq
         n2 = Mp*(j-1)
         sum = 0.0
         DO k = 1,Mp
            n1 = n1 + 1
            n2 = n2 + 1
            sum = sum + Coefl(n1)*R(n2)
         ENDDO
         l = l + 1
         R(l) = sum
      ENDDO
   ENDDO
!
!        CALCULATE EIGENVALUES WITH ASSOCIATED EIGENVECTORS OF THE
!        INVERSE OF R22 * A
!
   l = l + 1
   CALL nroot(Mq,R(n3),Coefr,Roots,R(l))
!
!     FOR EACH VALUE OF I = 1, 2, ..., MQ, CALCULATE THE FOLLOWING
!     STATISTICS
!
   SPAG_Loop_1_1: DO i = 1,Mq
!
!        TEST WHETHER EIGENVALUE IS GREATER THAN ZERO
!
      IF ( Roots(i)<=0 ) EXIT SPAG_Loop_1_1
!
!        CANONICAL CORRELATION
!
      Canr(i) = sqrt(Roots(i))
!
!        CHI-SQUARE
!
      Wlam(i) = 1.0
      DO j = i,Mq
         Wlam(i) = Wlam(i)*(1.0-Roots(j))
      ENDDO
      fn = N
      fmp = Mp
      fmq = Mq
      Chisq(i) = -(fn-0.5*(fmp+fmq+1.0))*alog(Wlam(i))
!
!        DEGREES OF FREEDOM FOR CHI-SQUARE
!
      n1 = i - 1
      Ndf(i) = (Mp-n1)*(Mq-n1)
!
!        I-TH SET OF RIGHT HAND COEFFICIENTS
!
      n1 = Mq*(i-1)
      n2 = Mq*(i-1) + l - 1
      DO j = 1,Mq
         n1 = n1 + 1
         n2 = n2 + 1
         Coefr(n1) = R(n2)
      ENDDO
!
!        I-TH SET OF LEFT HAND COEFFICIENTS
!
      DO j = 1,Mp
         n1 = j - Mp
         n2 = Mq*(i-1)
         k = Mp*(i-1) + j
         Coefl(k) = 0.0
         DO jj = 1,Mq
            n1 = n1 + Mp
            n2 = n2 + 1
            Coefl(k) = Coefl(k) + R(n1)*Coefr(n2)
         ENDDO
         Coefl(k) = Coefl(k)/Canr(i)
      ENDDO
   ENDDO SPAG_Loop_1_1
END SUBROUTINE canor
