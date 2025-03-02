C                                                                       DTHE  10
C     ..................................................................DTHE  20
C                                                                       DTHE  30
C        SUBROUTINE DTHEP                                               DTHE  40
C                                                                       DTHE  50
C        PURPOSE                                                        DTHE  60
C           A SERIES EXPANSION IN HERMITE POLYNOMIALS WITH INDEPENDENT  DTHE  70
C           VARIABLE X IS TRANSFORMED TO A POLYNOMIAL WITH INDEPENDENT  DTHE  80
C           VARIABLE Z, WHERE X=A*Z+B                                   DTHE  90
C                                                                       DTHE 100
C        USAGE                                                          DTHE 110
C           CALL DTHEP(A,B,POL,N,C,WORK)                                DTHE 120
C                                                                       DTHE 130
C        DESCRIPTION OF PARAMETERS                                      DTHE 140
C           A     - FACTOR OF LINEAR TERM IN GIVEN LINEAR TRANSFORMATIONDTHE 150
C                   DOUBLE PRECISION VARIABLE                           DTHE 160
C           B     - CONSTANT TERM IN GIVEN LINEAR TRANSFORMATION        DTHE 170
C                   DOUBLE PRECISION VARIABLE                           DTHE 180
C           POL   - COEFFICIENT VECTOR OF POLYNOMIAL (RESULTANT VALUE)  DTHE 190
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DTHE 200
C                   DOUBLE PRECISION VECTOR                             DTHE 210
C           N     - DIMENSION OF COEFFICIENT VECTOR POL AND C           DTHE 220
C           C     - COEFFICIENT VECTOR OF GIVEN EXPANSION               DTHE 230
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DTHE 240
C                   POL AND C MAY BE IDENTICALLY LOCATED                DTHE 250
C                   DOUBLE PRECISION VECTOR                             DTHE 260
C           WORK  - WORKING STORAGE OF DIMENSION 2*N                    DTHE 270
C                   DOUBLE PRECISION ARRAY                              DTHE 280
C                                                                       DTHE 290
C        REMARKS                                                        DTHE 300
C           COEFFICIENT VECTOR C REMAINS UNCHANGED IF NOT COINCIDING    DTHE 310
C           WITH COEFFICIENT VECTOR POL.                                DTHE 320
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1.                DTHE 330
C           THE LINEAR TRANSFORMATION X=A*Z+B OR Z=(1/A)(X-B) TRANSFORMSDTHE 340
C           THE RANGE (-C,C) IN X TO THE RANGE (ZL,ZR) IN Z WHERE       DTHE 350
C           ZL=-(C+B)/A AND ZR=(C-B)/A.                                 DTHE 360
C           FOR GIVEN ZL, ZR AND C WE HAVE A=2C/(ZR-ZL) AND             DTHE 370
C           B=-C(ZR+ZL)/(ZR-ZL)                                         DTHE 380
C                                                                       DTHE 390
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DTHE 400
C           NONE                                                        DTHE 410
C                                                                       DTHE 420
C        METHOD                                                         DTHE 430
C           THE TRANSFORMATION IS BASED ON THE RECURRENCE EQUATION      DTHE 440
C           FOR HERMITE POLYNOMIALS H(N,X)                              DTHE 450
C           H(N+1,X)=2*(X*H(N,X)-N*H(N-1,X)),                           DTHE 460
C           WHERE THE FIRST TERM IN BRACKETS IS THE INDEX               DTHE 470
C           THE SECOND IS THE ARGUMENT.                                 DTHE 480
C           STARTING VALUES ARE H(0,X)=1,H(1,X)=2*X.                    DTHE 490
C           THE TRANSFORMATION IS IMPLICITLY DEFINED BY MEANS OF        DTHE 500
C           X=A*Z+B TOGETHER WITH                                       DTHE 510
C           SUM(POL(I)*Z**(I-1), SUMMED OVER I FROM 1 TO N)             DTHE 520
C           =SUM(C(I)*H(I-1,X), SUMMED OVER I FROM 1 TO N).             DTHE 530
C                                                                       DTHE 540
C     ..................................................................DTHE 550
C                                                                       DTHE 560
      SUBROUTINE DTHEP(A,B,POL,N,C,WORK)                                DTHE 570
C                                                                       DTHE 580
      DIMENSION POL(1),C(1),WORK(1)                                     DTHE 590
      DOUBLE PRECISION A,B,POL,C,WORK,H,P,FI,XD,X0                      DTHE 600
C                                                                       DTHE 610
C        TEST OF DIMENSION                                              DTHE 620
      IF(N-1)2,1,3                                                      DTHE 630
C                                                                       DTHE 640
C        DIMENSION LESS THAN 2                                          DTHE 650
    1 POL(1)=C(1)                                                       DTHE 660
    2 RETURN                                                            DTHE 670
C                                                                       DTHE 680
    3 XD=A+A                                                            DTHE 690
      X0=B+B                                                            DTHE 700
      POL(1)=C(1)+C(2)*X0                                               DTHE 710
      POL(2)=C(2)*XD                                                    DTHE 720
      IF(N-2)2,2,4                                                      DTHE 730
C                                                                       DTHE 740
C        INITIALIZATION                                                 DTHE 750
    4 WORK(1)=1.D0                                                      DTHE 760
      WORK(2)=X0                                                        DTHE 770
      WORK(3)=0.D0                                                      DTHE 780
      WORK(4)=XD                                                        DTHE 790
      FI=2.D0                                                           DTHE 800
C                                                                       DTHE 810
C        CALCULATE COEFFICIENT VECTOR OF NEXT HERMITE POLYNOMIAL        DTHE 820
C        AND ADD MULTIPLE OF THIS VECTOR TO POLYNOMIAL POL              DTHE 830
      DO 6 J=3,N                                                        DTHE 840
      P=0.D0                                                            DTHE 850
C                                                                       DTHE 860
      DO 5 K=2,J                                                        DTHE 870
      H=P*XD+WORK(2*K-2)*X0-FI*WORK(2*K-3)                              DTHE 880
      P=WORK(2*K-2)                                                     DTHE 890
      WORK(2*K-2)=H                                                     DTHE 900
      WORK(2*K-3)=P                                                     DTHE 910
    5 POL(K-1)=POL(K-1)+H*C(J)                                          DTHE 920
      WORK(2*J-1)=0.D0                                                  DTHE 930
      WORK(2*J)=P*XD                                                    DTHE 940
      FI=FI+2.D0                                                        DTHE 950
    6 POL(J)=C(J)*WORK(2*J)                                             DTHE 960
      RETURN                                                            DTHE 970
      END                                                               DTHE 980
