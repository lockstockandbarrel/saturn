C                                                                       DTCN  10
C     ..................................................................DTCN  20
C                                                                       DTCN  30
C        SUBROUTINE DTCNP                                               DTCN  40
C                                                                       DTCN  50
C        PURPOSE                                                        DTCN  60
C           A SERIES EXPANSION IN CHEBYSHEV POLYNOMIALS WITH INDEPENDENTDTCN  70
C           VARIABLE X IS TRANSFORMED TO A POLYNOMIAL WITH INDEPENDENT  DTCN  80
C           VARIABLE Z, WHERE X=A*Z+B.                                  DTCN  90
C                                                                       DTCN 100
C        USAGE                                                          DTCN 110
C           CALL DTCNP(A,B,POL,N,C,WORK)                                DTCN 120
C                                                                       DTCN 130
C        DESCRIPTION OF PARAMETERS                                      DTCN 140
C           A     - FACTOR OF LINEAR TERM IN GIVEN LINEAR TRANSFORMATIONDTCN 150
C                   DOUBLE PRECISION VARIABLE                           DTCN 160
C           B     - CONSTANT TERM IN GIVEN LINEAR TRANSFORMATION        DTCN 170
C                   DOUBLE PRECISION VARIABLE                           DTCN 180
C           POL   - COEFFICIENT VECTOR OF POLYNOMIAL (RESULTANT VALUE)  DTCN 190
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DTCN 200
C                   DOUBLE PRECISION VECTOR                             DTCN 210
C           N     - DIMENSION OF COEFFICIENT VECTORS POL AND C          DTCN 220
C           C     - GIVEN COEFFICIENT VECTOR OF EXPANSION               DTCN 230
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DTCN 240
C                   POL AND C MAY BE IDENTICALLY LOCATED                DTCN 250
C                   DOUBLE PRECISION VECTOR                             DTCN 260
C           WORK  - WORKING STORAGE OF DIMENSION 2*N                    DTCN 270
C                   DOUBLE PRECISION ARRAY                              DTCN 280
C                                                                       DTCN 290
C        REMARKS                                                        DTCN 300
C           COEFFICIENT VECTOR C REMAINS UNCHANGED IF NOT COINCIDING    DTCN 310
C           WITH COEFFICIENT VECTOR POL.                                DTCN 320
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1.                DTCN 330
C           THE LINEAR TRANSFORMATION X=A*Z+B OR Z=(1/A)(X-B) TRANSFORMSDTCN 340
C           THE RANGE (-1,+1) IN X TO THE RANGE (ZL,ZR) IN Z, WHERE     DTCN 350
C           ZL=-(1+B)/A AND ZR=(1-B)/A.                                 DTCN 360
C           FOR GIVEN ZL, ZR WE HAVE A=2/(ZR-ZL) AND B=-(ZR+ZL)/(ZR-ZL) DTCN 370
C                                                                       DTCN 380
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DTCN 390
C           NONE                                                        DTCN 400
C                                                                       DTCN 410
C        METHOD                                                         DTCN 420
C           THE TRANSFORMATION IS BASED ON THE RECURRENCE EQUATION      DTCN 430
C           FOR CHEBYSHEV POLYNOMIALS T(N,X)                            DTCN 440
C           T(N+1,X)=2*X*T(N,X)-T(N-1,X),                               DTCN 450
C           WHERE THE FIRST TERM IN BRACKETS IS THE INDEX,              DTCN 460
C           THE SECOND IS THE ARGUMENT.                                 DTCN 470
C           STARTING VALUES ARE T(0,X)=1, T(1,X)=X.                     DTCN 480
C           THE TRANSFORMATION IS IMPLICITLY DEFINED BY MEANS OF        DTCN 490
C           X = A*Z+B TOGETHER WITH                                     DTCN 500
C           SUM(POL(I)*Z**(I-1), SUMMED OVER I FROM 1 TO N)             DTCN 510
C           =SUM(C(I)*T(I-1,X), SUMMED OVER I FROM 1 TO N).             DTCN 520
C                                                                       DTCN 530
C     ..................................................................DTCN 540
C                                                                       DTCN 550
      SUBROUTINE DTCNP(A,B,POL,N,C,WORK)                                DTCN 560
C                                                                       DTCN 570
      DIMENSION POL(1),C(1),WORK(1)                                     DTCN 580
      DOUBLE PRECISION A,B,POL,C,WORK,H,P,XD,X0                         DTCN 590
C                                                                       DTCN 600
C        TEST OF DIMENSION                                              DTCN 610
      IF(N-1)2,1,3                                                      DTCN 620
C                                                                       DTCN 630
C        DIMENSION LESS THAN 2                                          DTCN 640
    1 POL(1)=C(1)                                                       DTCN 650
    2 RETURN                                                            DTCN 660
C                                                                       DTCN 670
    3 POL(1)=C(1)+C(2)*B                                                DTCN 680
      POL(2)=C(2)*A                                                     DTCN 690
      IF(N-2)2,2,4                                                      DTCN 700
C                                                                       DTCN 710
C        INITIALIZATION                                                 DTCN 720
    4 WORK(1)=1.D0                                                      DTCN 730
      WORK(2)=B                                                         DTCN 740
      WORK(3)=0.D0                                                      DTCN 750
      WORK(4)=A                                                         DTCN 760
      XD=A+A                                                            DTCN 770
      X0=B+B                                                            DTCN 780
C                                                                       DTCN 790
C        CALCULATE COEFFICIENT VECTOR OF NEXT CHEBYSHEV POLYNOMIAL      DTCN 800
C        AND ADD MULTIPLE OF THIS VECTOR TO POLYNOMIAL POL              DTCN 810
      DO 6 J=3,N                                                        DTCN 820
      P=0.D0                                                            DTCN 830
C                                                                       DTCN 840
      DO 5 K=2,J                                                        DTCN 850
      H=P-WORK(2*K-3)+X0*WORK(2*K-2)                                    DTCN 860
      P=WORK(2*K-2)                                                     DTCN 870
      WORK(2*K-2)=H                                                     DTCN 880
      WORK(2*K-3)=P                                                     DTCN 890
      POL(K-1)=POL(K-1)+H*C(J)                                          DTCN 900
    5 P=XD*P                                                            DTCN 910
      WORK(2*J-1)=0.D0                                                  DTCN 920
      WORK(2*J)=P                                                       DTCN 930
    6 POL(J)=C(J)*P                                                     DTCN 940
      RETURN                                                            DTCN 950
      END                                                               DTCN 960
