C                                                                       DTLE  10
C     ..................................................................DTLE  20
C                                                                       DTLE  30
C        SUBROUTINE DTLEP                                               DTLE  40
C                                                                       DTLE  50
C        PURPOSE                                                        DTLE  60
C           A SERIES EXPANSION IN LEGENDRE POLYNOMIALS WITH INDEPENDENT DTLE  70
C           VARIABLE X IS TRANSFORMED TO A POLYNOMIAL WITH INDEPENDENT  DTLE  80
C           VARIABLE Z, WHERE X=A*Z+B                                   DTLE  90
C                                                                       DTLE 100
C        USAGE                                                          DTLE 110
C           CALL DTLEP(A,B,POL,N,C,WORK)                                DTLE 120
C                                                                       DTLE 130
C        DESCRIPTION OF PARAMETERS                                      DTLE 140
C           A     - FACTOR OF LINEAR TERM IN GIVEN LINEAR TRANSFORMATIONDTLE 150
C                   DOUBLE PRECISION VARIABLE                           DTLE 160
C           B     - CONSTANT TERM IN GIVEN LINEAR TRANSFORMATION        DTLE 170
C                   DOUBLE PRECISION VARIABLE                           DTLE 180
C           POL   - COEFFICIENT VECTOR OF POLYNOMIAL (RESULTANT VALUE)  DTLE 190
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DTLE 200
C                   DOUBLE PRECISION VECTOR                             DTLE 210
C           N     - DIMENSION OF COEFFICIENT VECTORS POL AND C          DTLE 220
C           C     - GIVEN COEFFICIENT VECTOR OF EXPANSION               DTLE 230
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DTLE 240
C                   POL AND C MAY BE IDENTICALLY LOCATED                DTLE 250
C                   DOUBLE PRECISION VECTOR                             DTLE 260
C           WORK  - WORKING STORAGE OF DIMENSION 2*N                    DTLE 270
C                   DOUBLE PRECISION ARRAY                              DTLE 280
C                                                                       DTLE 290
C        REMARKS                                                        DTLE 300
C           COEFFICIENT VECTOR C REMAINS UNCHANGED IF NOT COINCIDING    DTLE 310
C           WITH COEFFICIENT VECTOR POL.                                DTLE 320
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1.                DTLE 330
C           THE LINEAR TRANSFORMATION X=A*Z+B OR Z=(1/A)(X-B) TRANSFORMSDTLE 340
C           THE RANGE (-1,+1) IN X TO THE RANGE (ZL,ZR) IN Z, WHERE     DTLE 350
C           ZL=-(1+B)/A AND ZR=(1-B)/A.                                 DTLE 360
C           FOR GIVEN ZL, ZR WE HAVE A=2/(ZR-ZL) AND B=-(ZR+ZL)/(ZR-ZL) DTLE 370
C                                                                       DTLE 380
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DTLE 390
C           NONE                                                        DTLE 400
C                                                                       DTLE 410
C        METHOD                                                         DTLE 420
C           THE TRANSFORMATION IS BASED ON THE RECURRENCE EQUATION      DTLE 430
C           FOR LEGENDRE POLYNOMIALS P(N,X)                             DTLE 440
C           P(N+1,X)=2*X*P(N,X)-P(N-1,X)-(X*P(N,X)-P(N-1,X))/(N+1),     DTLE 450
C           WHERE THE FIRST TERM IN BRACKETS IS THE INDEX,              DTLE 460
C           THE SECOND IS THE ARGUMENT.                                 DTLE 470
C           STARTING VALUES ARE P(0,X)=1, P(1,X)=X.                     DTLE 480
C           THE TRANSFORMATION IS IMPLICITLY DEFINED BY MEANS OF        DTLE 490
C           X=A*Z+B TOGETHER WITH                                       DTLE 500
C           SUM(POL(I)*Z**(I-1), SUMMED OVER I FROM 1 TO N)             DTLE 510
C           =SUM(C(I)*P(I-1,X), SUMMED OVER I FROM 1 TO N).             DTLE 520
C                                                                       DTLE 530
C     ..................................................................DTLE 540
C                                                                       DTLE 550
      SUBROUTINE DTLEP(A,B,POL,N,C,WORK)                                DTLE 560
C                                                                       DTLE 570
      DIMENSION POL(1),C(1),WORK(1)                                     DTLE 580
      DOUBLE PRECISION A,B,POL,C,WORK,H,P,Q,Q1,FI                       DTLE 590
C                                                                       DTLE 600
C        TEST OF DIMENSION                                              DTLE 610
      IF(N-1)2,1,3                                                      DTLE 620
C                                                                       DTLE 630
C        DIMENSION LESS THAN 2                                          DTLE 640
    1 POL(1)=C(1)                                                       DTLE 650
    2 RETURN                                                            DTLE 660
C                                                                       DTLE 670
    3 POL(1)=C(1)+B*C(2)                                                DTLE 680
      POL(2)=A*C(2)                                                     DTLE 690
      IF(N-2)2,2,4                                                      DTLE 700
C                                                                       DTLE 710
C        INITIALIZATION                                                 DTLE 720
    4 WORK(1)=1.D0                                                      DTLE 730
      WORK(2)=B                                                         DTLE 740
      WORK(3)=0.D0                                                      DTLE 750
      WORK(4)=A                                                         DTLE 760
      FI=1.D0                                                           DTLE 770
C                                                                       DTLE 780
C        CALCULATE COEFFICIENT VECTOR OF NEXT LEGENDRE POLYNOMIAL       DTLE 790
C        AND ADD MULTIPLE OF THIS VECTOR TO POLYNOMIAL POL              DTLE 800
      DO 6 J=3,N                                                        DTLE 810
      FI=FI+1.D0                                                        DTLE 820
      Q=1.D0/FI-1.D0                                                    DTLE 830
      Q1=1.D0-Q                                                         DTLE 840
      P=0.D0                                                            DTLE 850
C                                                                       DTLE 860
      DO 5 K=2,J                                                        DTLE 870
      H=(A*P+B*WORK(2*K-2))*Q1+Q*WORK(2*K-3)                            DTLE 880
      P=WORK(2*K-2)                                                     DTLE 890
      WORK(2*K-2)=H                                                     DTLE 900
      WORK(2*K-3)=P                                                     DTLE 910
    5 POL(K-1)=POL(K-1)+H*C(J)                                          DTLE 920
      WORK(2*J-1)=0.D0                                                  DTLE 930
      WORK(2*J)=A*P*Q1                                                  DTLE 940
    6 POL(J)=C(J)*WORK(2*J)                                             DTLE 950
      RETURN                                                            DTLE 960
      END                                                               DTLE 970
