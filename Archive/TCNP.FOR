C                                                                       TCNP  10
C     ..................................................................TCNP  20
C                                                                       TCNP  30
C        SUBROUTINE TCNP                                                TCNP  40
C                                                                       TCNP  50
C        PURPOSE                                                        TCNP  60
C           A SERIES EXPANSION IN CHEBYSHEV POLYNOMIALS WITH INDEPENDENTTCNP  70
C           VARIABLE X IS TRANSFORMED TO A POLYNOMIAL WITH INDEPENDENT  TCNP  80
C           VARIABLE Z, WHERE X=A*Z+B.                                  TCNP  90
C                                                                       TCNP 100
C        USAGE                                                          TCNP 110
C           CALL TCNP(A,B,POL,N,C,WORK)                                 TCNP 120
C                                                                       TCNP 130
C        DESCRIPTION OF PARAMETERS                                      TCNP 140
C           A     - FACTOR OF LINEAR TERM IN GIVEN LINEAR TRANSFORMATIONTCNP 150
C           B     - CONSTANT TERM IN GIVEN LINEAR TRANSFORMATION        TCNP 160
C           POL   - COEFFICIENT VECTOR OF POLYNOMIAL (RESULTANT VALUE)  TCNP 170
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           TCNP 180
C           N     - DIMENSION OF COEFFICIENT VECTORS POL AND C          TCNP 190
C           C     - GIVEN COEFFICIENT VECTOR OF EXPANSION               TCNP 200
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           TCNP 210
C                   POL AND C MAY BE IDENTICALLY LOCATED                TCNP 220
C           WORK  - WORKING STORAGE OF DIMENSION 2*N                    TCNP 230
C                                                                       TCNP 240
C        REMARKS                                                        TCNP 250
C           COEFFICIENT VECTOR C REMAINS UNCHANGED IF NOT COINCIDING    TCNP 260
C           WITH COEFFICIENT VECTOR POL.                                TCNP 270
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1.                TCNP 280
C           THE LINEAR TRANSFORMATION X=A*Z+B OR Z=(1/A)(X-B) TRANSFORMSTCNP 290
C           THE RANGE (-1,+1) IN X TO THE RANGE (ZL,ZR) IN Z, WHERE     TCNP 300
C           ZL=-(1+B)/A AND ZR=(1-B)/A.                                 TCNP 310
C           FOR GIVEN ZL, ZR WE HAVE A=2/(ZR-ZL) AND B=-(ZR+ZL)/(ZR-ZL) TCNP 320
C                                                                       TCNP 330
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  TCNP 340
C           NONE                                                        TCNP 350
C                                                                       TCNP 360
C        METHOD                                                         TCNP 370
C           THE TRANSFORMATION IS BASED ON THE RECURRENCE EQUATION      TCNP 380
C           FOR CHEBYSHEV POLYNOMIALS T(N,X)                            TCNP 390
C           T(N+1,X)=2*X*T(N,X)-T(N-1,X),                               TCNP 400
C           WHERE THE FIRST TERM IN BRACKETS IS THE INDEX,              TCNP 410
C           THE SECOND IS THE ARGUMENT.                                 TCNP 420
C           STARTING VALUES ARE T(0,X)=1, T(1,X)=X.                     TCNP 430
C           THE TRANSFORMATION IS IMPLICITLY DEFINED BY MEANS OF        TCNP 440
C           X = A*Z+B TOGETHER WITH                                     TCNP 450
C           SUM(POL(I)*Z**(I-1), SUMMED OVER I FROM 1 TO N)             TCNP 460
C           =SUM(C(I)*T(I-1,X), SUMMED OVER I FROM 1 TO N).             TCNP 470
C                                                                       TCNP 480
C     ..................................................................TCNP 490
C                                                                       TCNP 500
      SUBROUTINE TCNP(A,B,POL,N,C,WORK)                                 TCNP 510
C                                                                       TCNP 520
      DIMENSION POL(1),C(1),WORK(1)                                     TCNP 530
C                                                                       TCNP 540
C        TEST OF DIMENSION                                              TCNP 550
      IF(N-1)2,1,3                                                      TCNP 560
C                                                                       TCNP 570
C        DIMENSION LESS THAN 2                                          TCNP 580
    1 POL(1)=C(1)                                                       TCNP 590
    2 RETURN                                                            TCNP 600
C                                                                       TCNP 610
    3 POL(1)=C(1)+C(2)*B                                                TCNP 620
      POL(2)=C(2)*A                                                     TCNP 630
      IF(N-2)2,2,4                                                      TCNP 640
C                                                                       TCNP 650
C        INITIALIZATION                                                 TCNP 660
    4 WORK(1)=1.                                                        TCNP 670
      WORK(2)=B                                                         TCNP 680
      WORK(3)=0.                                                        TCNP 690
      WORK(4)=A                                                         TCNP 700
      XD=A+A                                                            TCNP 710
      X0=B+B                                                            TCNP 720
C                                                                       TCNP 730
C        CALCULATE COEFFICIENT VECTOR OF NEXT CHEBYSHEV POLYNOMIAL      TCNP 740
C        AND ADD MULTIPLE OF THIS VECTOR TO POLYNOMIAL POL              TCNP 750
      DO 6 J=3,N                                                        TCNP 760
      P=0.                                                              TCNP 770
C                                                                       TCNP 780
      DO 5 K=2,J                                                        TCNP 790
      H=P-WORK(2*K-3)+X0*WORK(2*K-2)                                    TCNP 800
      P=WORK(2*K-2)                                                     TCNP 810
      WORK(2*K-2)=H                                                     TCNP 820
      WORK(2*K-3)=P                                                     TCNP 830
      POL(K-1)=POL(K-1)+H*C(J)                                          TCNP 840
    5 P=XD*P                                                            TCNP 850
      WORK(2*J-1)=0.                                                    TCNP 860
      WORK(2*J)=P                                                       TCNP 870
    6 POL(J)=C(J)*P                                                     TCNP 880
      RETURN                                                            TCNP 890
      END                                                               TCNP 900
