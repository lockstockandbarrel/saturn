C                                                                       DTLA  10
C     ..................................................................DTLA  20
C                                                                       DTLA  30
C        SUBROUTINE DTLAP                                               DTLA  40
C                                                                       DTLA  50
C        PURPOSE                                                        DTLA  60
C           A SERIES EXPANSION IN LAGUERRE POLYNOMIALS WITH INDEPENDENT DTLA  70
C           VARIABLE X IS TRANSFORMED TO A POLYNOMIAL WITH INDEPENDENT  DTLA  80
C           VARIABLE Z, WHERE X=A*Z+B                                   DTLA  90
C                                                                       DTLA 100
C        USAGE                                                          DTLA 110
C           CALL DTLAP(A,B,POL,N,C,WORK)                                DTLA 120
C                                                                       DTLA 130
C        DESCRIPTION OF PARAMETERS                                      DTLA 140
C           A     - FACTOR OF LINEAR TERM IN GIVEN LINEAR TRANSFORMATIONDTLA 150
C                   DOUBLE PRECISION VARIABLE                           DTLA 160
C           B     - CONSTANT TERM IN GIVEN LINEAR TRANSFORMATION        DTLA 170
C                   DOUBLE PRECISION VARIABLE                           DTLA 180
C           POL   - COEFFICIENT VECTOR OF POLYNOMIAL (RESULTANT VALUE)  DTLA 190
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DTLA 200
C                   DOUBLE PRECISION VECTOR                             DTLA 210
C           N     - DIMENSION OF COEFFICIENT VECTORS POL AND C          DTLA 220
C           C     - GIVEN COEFFICIENT VECTOR OF EXPANSION               DTLA 230
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DTLA 240
C                   POL AND C MAY BE IDENTICALLY LOCATED                DTLA 250
C                   DOUBLE PRECISION VECTOR                             DTLA 260
C           WORK  - WORKING STORAGE OF DIMENSION 2*N                    DTLA 270
C                   DOUBLE PRECISION ARRAY                              DTLA 280
C                                                                       DTLA 290
C        REMARKS                                                        DTLA 300
C           COEFFICIENT VECTOR C REMAINS UNCHANGED IF NOT COINCIDING    DTLA 310
C           WITH COEFFICIENT VECTOR POL.                                DTLA 320
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1.                DTLA 330
C           THE LINEAR TRANSFORMATION X=A*Z+B OR Z=(1/A)(X-B) TRANSFORMSDTLA 340
C           THE RANGE (0,C) IN X TO THE RANGE (ZL,ZR) IN Z, WHERE       DTLA 350
C           ZL=-B/A AND ZR=(C-B)/A.                                     DTLA 360
C           FOR GIVEN ZL, ZR AND C WE HAVE A=C/(ZR-ZL) AND              DTLA 370
C           B=-C*ZL/(ZR-ZL)                                             DTLA 380
C                                                                       DTLA 390
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DTLA 400
C           NONE                                                        DTLA 410
C                                                                       DTLA 420
C        METHOD                                                         DTLA 430
C           THE TRANSFORMATION IS BASED ON THE RECURRENCE EQUATION      DTLA 440
C           FOR LAGUERRE POLYNOMIALS L(N,X)                             DTLA 450
C           L(N+1,X)=2*L(N,X)-L(N-1,X)-((1+X)*L(N,X)-L(N-1,X))/(N+1),   DTLA 460
C           WHERE THE FIRST TERM IN BRACKETS IS THE INDEX,              DTLA 470
C           THE SECOND IS THE ARGUMENT.                                 DTLA 480
C           STARTING VALUES ARE L(0,X)=1, L(1,X)=1-X.                   DTLA 490
C           THE TRANSFORMATION IS IMPLICITLY DEFINED BY MEANS OF        DTLA 500
C           X=A*Z+B TOGETHER WITH                                       DTLA 510
C           SUM(POL(I)*Z**(I-1), SUMMED OVER I FROM 1 TO N)             DTLA 520
C           =SUM(C(I)*L(I-1,X), SUMMED OVER I FROM 1 TO N).             DTLA 530
C                                                                       DTLA 540
C     ..................................................................DTLA 550
C                                                                       DTLA 560
      SUBROUTINE DTLAP(A,B,POL,N,C,WORK)                                DTLA 570
C                                                                       DTLA 580
      DIMENSION POL(1),C(1),WORK(1)                                     DTLA 590
      DOUBLE PRECISION A,B,POL,C,WORK,H,P,Q,Q1,Q2,FI                    DTLA 600
C                                                                       DTLA 610
C        TEST OF DIMENSION                                              DTLA 620
      IF(N-1)2,1,3                                                      DTLA 630
C                                                                       DTLA 640
C        DIMENSION LESS THAN 2                                          DTLA 650
    1 POL(1)=C(1)                                                       DTLA 660
    2 RETURN                                                            DTLA 670
C                                                                       DTLA 680
    3 POL(1)=C(1)+C(2)-B*C(2)                                           DTLA 690
      POL(2)=-C(2)*A                                                    DTLA 700
      IF(N-2)2,2,4                                                      DTLA 710
C                                                                       DTLA 720
C        INITIALIZATION                                                 DTLA 730
    4 WORK(1)=1.D0                                                      DTLA 740
      WORK(2)=1.D0-B                                                    DTLA 750
      WORK(3)=0.D0                                                      DTLA 760
      WORK(4)=-A                                                        DTLA 770
      FI=1.D0                                                           DTLA 780
C                                                                       DTLA 790
C        CALCULATE COEFFICIENT VECTOR OF NEXT LAGUERRE POLYNOMIAL       DTLA 800
C        AND ADD MULTIPLE OF THIS VECTOR TO POLYNOMIAL POL              DTLA 810
      DO 6 J=3,N                                                        DTLA 820
      FI=FI+1.D0                                                        DTLA 830
      Q=1.D0/FI                                                         DTLA 840
      Q1=Q-1.D0                                                         DTLA 850
      Q2=1.D0-Q1-B*Q                                                    DTLA 860
      Q=Q*A                                                             DTLA 870
      P=0.D0                                                            DTLA 880
C                                                                       DTLA 890
      DO 5 K=2,J                                                        DTLA 900
      H=-P*Q+WORK(2*K-2)*Q2+WORK(2*K-3)*Q1                              DTLA 910
      P=WORK(2*K-2)                                                     DTLA 920
      WORK(2*K-2)=H                                                     DTLA 930
      WORK(2*K-3)=P                                                     DTLA 940
    5 POL(K-1)=POL(K-1)+H*C(J)                                          DTLA 950
      WORK(2*J-1)=0.D0                                                  DTLA 960
      WORK(2*J)=-Q*P                                                    DTLA 970
    6 POL(J)=C(J)*WORK(2*J)                                             DTLA 980
      RETURN                                                            DTLA 990
      END                                                               DTLA1000
