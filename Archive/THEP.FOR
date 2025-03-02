C                                                                       THEP  10
C     ..................................................................THEP  20
C                                                                       THEP  30
C        SUBROUTINE THEP                                                THEP  40
C                                                                       THEP  50
C        PURPOSE                                                        THEP  60
C           A SERIES EXPANSION IN HERMITE POLYNOMIALS WITH INDEPENDENT  THEP  70
C           VARIABLE X IS TRANSFORMED TO A POLYNOMIAL WITH INDEPENDENT  THEP  80
C           VARIABLE Z, WHERE X=A*Z+B                                   THEP  90
C                                                                       THEP 100
C        USAGE                                                          THEP 110
C           CALL THEP(A,B,POL,N,C,WORK)                                 THEP 120
C                                                                       THEP 130
C        DESCRIPTION OF PARAMETERS                                      THEP 140
C           A     - FACTOR OF LINEAR TERM IN GIVEN LINEAR TRANSFORMATIONTHEP 150
C           B     - CONSTANT TERM IN GIVEN LINEAR TRANSFORMATION        THEP 160
C           POL   - COEFFICIENT VECTOR OF POLYNOMIAL (RESULTANT VALUE)  THEP 170
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           THEP 180
C           N     - DIMENSION OF COEFFICIENT VECTOR POL AND C           THEP 190
C           C     - COEFFICIENT VECTOR OF GIVEN EXPANSION               THEP 200
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           THEP 210
C                   POL AND C MAY BE IDENTICALLY LOCATED                THEP 220
C           WORK  - WORKING STORAGE OF DIMENSION 2*N                    THEP 230
C                                                                       THEP 240
C        REMARKS                                                        THEP 250
C           COEFFICIENT VECTOR C REMAINS UNCHANGED IF NOT COINCIDING    THEP 260
C           WITH COEFFICIENT VECTOR POL.                                THEP 270
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1.                THEP 280
C           THE LINEAR TRANSFORMATION X=A*Z+B OR Z=(1/A)(X-B) TRANSFORMSTHEP 290
C           THE RANGE (-C,C) IN X TO THE RANGE (ZL,ZR) IN Z WHERE       THEP 300
C           ZL=-(C+B)/A AND ZR=(C-B)/A.                                 THEP 310
C           FOR GIVEN ZL, ZR AND C WE HAVE A=2C/(ZR-ZL) AND             THEP 320
C           B=-C(ZR+ZL)/(ZR-ZL)                                         THEP 330
C                                                                       THEP 340
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  THEP 350
C           NONE                                                        THEP 360
C                                                                       THEP 370
C        METHOD                                                         THEP 380
C           THE TRANSFORMATION IS BASED ON THE RECURRENCE EQUATION      THEP 390
C           FOR HERMITE POLYNOMIALS H(N,X)                              THEP 400
C           H(N+1,X)=2*(X*H(N,X)-N*H(N-1,X)),                           THEP 410
C           WHERE THE FIRST TERM IN BRACKETS IS THE INDEX               THEP 420
C           THE SECOND IS THE ARGUMENT.                                 THEP 430
C           STARTING VALUES ARE H(0,X)=1,H(1,X)=2*X.                    THEP 440
C           THE TRANSFORMATION IS IMPLICITLY DEFINED BY MEANS OF        THEP 450
C           X=A*Z+B TOGETHER WITH                                       THEP 460
C           SUM(POL(I)*Z**(I-1), SUMMED OVER I FROM 1 TO N)             THEP 470
C           =SUM(C(I)*H(I-1,X), SUMMED OVER I FROM 1 TO N).             THEP 480
C                                                                       THEP 490
C     ..................................................................THEP 500
C                                                                       THEP 510
      SUBROUTINE THEP(A,B,POL,N,C,WORK)                                 THEP 520
C                                                                       THEP 530
      DIMENSION POL(1),C(1),WORK(1)                                     THEP 540
C                                                                       THEP 550
C        TEST OF DIMENSION                                              THEP 560
      IF(N-1)2,1,3                                                      THEP 570
C                                                                       THEP 580
C        DIMENSION LESS THAN 2                                          THEP 590
    1 POL(1)=C(1)                                                       THEP 600
    2 RETURN                                                            THEP 610
C                                                                       THEP 620
    3 XD=A+A                                                            THEP 630
      X0=B+B                                                            THEP 640
      POL(1)=C(1)+C(2)*X0                                               THEP 650
      POL(2)=C(2)*XD                                                    THEP 660
      IF(N-2)2,2,4                                                      THEP 670
C                                                                       THEP 680
C        INITIALIZATION                                                 THEP 690
    4 WORK(1)=1.                                                        THEP 700
      WORK(2)=X0                                                        THEP 710
      WORK(3)=0.                                                        THEP 720
      WORK(4)=XD                                                        THEP 730
      FI=2.                                                             THEP 740
C                                                                       THEP 750
C        CALCULATE COEFFICIENT VECTOR OF NEXT HERMITE POLYNOMIAL        THEP 760
C        AND ADD MULTIPLE OF THIS VECTOR TO POLYNOMIAL POL              THEP 770
      DO 6 J=3,N                                                        THEP 780
      P=0.                                                              THEP 790
C                                                                       THEP 800
      DO 5 K=2,J                                                        THEP 810
      H=P*XD+WORK(2*K-2)*X0-FI*WORK(2*K-3)                              THEP 820
      P=WORK(2*K-2)                                                     THEP 830
      WORK(2*K-2)=H                                                     THEP 840
      WORK(2*K-3)=P                                                     THEP 850
    5 POL(K-1)=POL(K-1)+H*C(J)                                          THEP 860
      WORK(2*J-1)=0.                                                    THEP 870
      WORK(2*J)=P*XD                                                    THEP 880
      FI=FI+2.                                                          THEP 890
    6 POL(J)=C(J)*WORK(2*J)                                             THEP 900
      RETURN                                                            THEP 910
      END                                                               THEP 920
