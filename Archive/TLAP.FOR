C                                                                       TLAP  10
C     ..................................................................TLAP  20
C                                                                       TLAP  30
C        SUBROUTINE TLAP                                                TLAP  40
C                                                                       TLAP  50
C        PURPOSE                                                        TLAP  60
C           A SERIES EXPANSION IN LAGUERRE POLYNOMIALS WITH INDEPENDENT TLAP  70
C           VARIABLE X IS TRANSFORMED TO A POLYNOMIAL WITH INDEPENDENT  TLAP  80
C           VARIABLE Z, WHERE X=A*Z+B                                   TLAP  90
C                                                                       TLAP 100
C        USAGE                                                          TLAP 110
C           CALL TLAP(A,B,POL,N,C,WORK)                                 TLAP 120
C                                                                       TLAP 130
C        DESCRIPTION OF PARAMETERS                                      TLAP 140
C           A     - FACTOR OF LINEAR TERM IN GIVEN LINEAR TRANSFORMATIONTLAP 150
C           B     - CONSTANT TERM IN GIVEN LINEAR TRANSFORMATION        TLAP 160
C           POL   - COEFFICIENT VECTOR OF POLYNOMIAL (RESULTANT VALUE)  TLAP 170
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           TLAP 180
C           N     - DIMENSION OF COEFFICIENT VECTORS POL AND C          TLAP 190
C           C     - GIVEN COEFFICIENT VECTOR OF EXPANSION               TLAP 200
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           TLAP 210
C                   POL AND C MAY BE IDENTICALLY LOCATED                TLAP 220
C           WORK  - WORKING STORAGE OF DIMENSION 2*N                    TLAP 230
C                                                                       TLAP 240
C        REMARKS                                                        TLAP 250
C           COEFFICIENT VECTOR C REMAINS UNCHANGED IF NOT COINCIDING    TLAP 260
C           WITH COEFFICIENT VECTOR POL.                                TLAP 270
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1.                TLAP 280
C           THE LINEAR TRANSFORMATION X=A*Z+B OR Z=(1/A)(X-B) TRANSFORMSTLAP 290
C           THE RANGE (0,C) IN X TO THE RANGE (ZL,ZR) IN Z, WHERE       TLAP 300
C           ZL=-B/A AND ZR=(C-B)/A.                                     TLAP 310
C           FOR GIVEN ZL, ZR AND C WE HAVE A=C/(ZR-ZL) AND              TLAP 320
C           B=-C*ZL/(ZR-ZL)                                             TLAP 330
C                                                                       TLAP 340
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  TLAP 350
C           NONE                                                        TLAP 360
C                                                                       TLAP 370
C        METHOD                                                         TLAP 380
C           THE TRANSFORMATION IS BASED ON THE RECURRENCE EQUATION      TLAP 390
C           FOR LAGUERRE POLYNOMIALS L(N,X)                             TLAP 400
C           L(N+1,X)=2*L(N,X)-L(N-1,X)-((1+X)*L(N,X)-L(N-1,X))/(N+1),   TLAP 410
C           WHERE THE FIRST TERM IN BRACKETS IS THE INDEX,              TLAP 420
C           THE SECOND IS THE ARGUMENT.                                 TLAP 430
C           STARTING VALUES ARE L(0,X)=1, L(1,X)=1-X.                   TLAP 440
C           THE TRANSFORMATION IS IMPLICITLY DEFINED BY MEANS OF        TLAP 450
C           X=A*Z+B TOGETHER WITH                                       TLAP 460
C           SUM(POL(I)*Z**(I-1), SUMMED OVER I FROM 1 TO N)             TLAP 470
C           =SUM(C(I)*L(I-1,X), SUMMED OVER I FROM 1 TO N).             TLAP 480
C                                                                       TLAP 490
C     ..................................................................TLAP 500
C                                                                       TLAP 510
      SUBROUTINE TLAP(A,B,POL,N,C,WORK)                                 TLAP 520
C                                                                       TLAP 530
      DIMENSION POL(1),C(1),WORK(1)                                     TLAP 540
C                                                                       TLAP 550
C        TEST OF DIMENSION                                              TLAP 560
      IF(N-1)2,1,3                                                      TLAP 570
C                                                                       TLAP 580
C        DIMENSION LESS THAN 2                                          TLAP 590
    1 POL(1)=C(1)                                                       TLAP 600
    2 RETURN                                                            TLAP 610
C                                                                       TLAP 620
    3 POL(1)=C(1)+C(2)-B*C(2)                                           TLAP 630
      POL(2)=-C(2)*A                                                    TLAP 640
      IF(N-2)2,2,4                                                      TLAP 650
C                                                                       TLAP 660
C        INITIALIZATION                                                 TLAP 670
    4 WORK(1)=1.                                                        TLAP 680
      WORK(2)=1.D0-B                                                    TLAP 690
      WORK(3)=0.                                                        TLAP 700
      WORK(4)=-A                                                        TLAP 710
      FI=1.                                                             TLAP 720
C                                                                       TLAP 730
C        CALCULATE COEFFICIENT VECTOR OF NEXT LAGUERRE POLYNOMIAL       TLAP 740
C        AND ADD MULTIPLE OF THIS VECTOR TO POLYNOMIAL POL              TLAP 750
      DO 6 J=3,N                                                        TLAP 760
      FI=FI+1.                                                          TLAP 770
      Q=1./FI                                                           TLAP 780
      Q1=Q-1.                                                           TLAP 790
      Q2=1.-Q1-B*Q                                                      TLAP 800
      Q=Q*A                                                             TLAP 810
      P=0.                                                              TLAP 820
C                                                                       TLAP 830
      DO 5 K=2,J                                                        TLAP 840
      H=-P*Q+WORK(2*K-2)*Q2+WORK(2*K-3)*Q1                              TLAP 850
      P=WORK(2*K-2)                                                     TLAP 860
      WORK(2*K-2)=H                                                     TLAP 870
      WORK(2*K-3)=P                                                     TLAP 880
    5 POL(K-1)=POL(K-1)+H*C(J)                                          TLAP 890
      WORK(2*J-1)=0.                                                    TLAP 900
      WORK(2*J)=-Q*P                                                    TLAP 910
    6 POL(J)=C(J)*WORK(2*J)                                             TLAP 920
      RETURN                                                            TLAP 930
      END                                                               TLAP 940
