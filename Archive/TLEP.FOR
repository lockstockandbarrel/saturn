C                                                                       TLEP  10
C     ..................................................................TLEP  20
C                                                                       TLEP  30
C        SUBROUTINE TLEP                                                TLEP  40
C                                                                       TLEP  50
C        PURPOSE                                                        TLEP  60
C           A SERIES EXPANSION IN LEGENDRE POLYNOMIALS WITH INDEPENDENT TLEP  70
C           VARIABLE X IS TRANSFORMED TO A POLYNOMIAL WITH INDEPENDENT  TLEP  80
C           VARIABLE Z, WHERE X=A*Z+B                                   TLEP  90
C                                                                       TLEP 100
C        USAGE                                                          TLEP 110
C           CALL TLEP(A,B,POL,N,C,WORK)                                 TLEP 120
C                                                                       TLEP 130
C        DESCRIPTION OF PARAMETERS                                      TLEP 140
C           A     - FACTOR OF LINEAR TERM IN GIVEN LINEAR TRANSFORMATIONTLEP 150
C           B     - CONSTANT TERM IN GIVEN LINEAR TRANSFORMATION        TLEP 160
C           POL   - COEFFICIENT VECTOR OF POLYNOMIAL (RESULTANT VALUE)  TLEP 170
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           TLEP 180
C           N     - DIMENSION OF COEFFICIENT VECTORS POL AND C          TLEP 190
C           C     - GIVEN COEFFICIENT VECTOR OF EXPANSION               TLEP 200
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           TLEP 210
C                   POL AND C MAY BE IDENTICALLY LOCATED                TLEP 220
C           WORK  - WORKING STORAGE OF DIMENSION 2*N                    TLEP 230
C                                                                       TLEP 240
C        REMARKS                                                        TLEP 250
C           COEFFICIENT VECTOR C REMAINS UNCHANGED IF NOT COINCIDING    TLEP 260
C           WITH COEFFICIENT VECTOR POL.                                TLEP 270
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1.                TLEP 280
C           THE LINEAR TRANSFORMATION X=A*Z+B OR Z=(1/A)(X-B) TRANSFORMSTLEP 290
C           THE RANGE (-1,+1) IN X TO THE RANGE (ZL,ZR) IN Z, WHERE     TLEP 300
C           ZL=-(1+B)/A AND ZR=(1-B)/A.                                 TLEP 310
C           FOR GIVEN ZL, ZR WE HAVE A=2/(ZR-ZL) AND B=-(ZR+ZL)/(ZR-ZL) TLEP 320
C                                                                       TLEP 330
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  TLEP 340
C           NONE                                                        TLEP 350
C                                                                       TLEP 360
C        METHOD                                                         TLEP 370
C           THE TRANSFORMATION IS BASED ON THE RECURRENCE EQUATION      TLEP 380
C           FOR LEGENDRE POLYNOMIALS P(N,X)                             TLEP 390
C           P(N+1,X)=2*X*P(N,X)-P(N-1,X)-(X*P(N,X)-P(N-1,X))/(N+1),     TLEP 400
C           WHERE THE FIRST TERM IN BRACKETS IS THE INDEX,              TLEP 410
C           THE SECOND IS THE ARGUMENT.                                 TLEP 420
C           STARTING VALUES ARE P(0,X)=1, P(1,X)=X.                     TLEP 430
C           THE TRANSFORMATION IS IMPLICITLY DEFINED BY MEANS OF        TLEP 440
C           X=A*Z+B TOGETHER WITH                                       TLEP 450
C           SUM(POL(I)*Z**(I-1), SUMMED OVER I FROM 1 TO N)             TLEP 460
C           =SUM(C(I)*P(I-1,X), SUMMED OVER I FROM 1 TO N).             TLEP 470
C                                                                       TLEP 480
C     ..................................................................TLEP 490
C                                                                       TLEP 500
      SUBROUTINE TLEP(A,B,POL,N,C,WORK)                                 TLEP 510
C                                                                       TLEP 520
      DIMENSION POL(1),C(1),WORK(1)                                     TLEP 530
C                                                                       TLEP 540
C        TEST OF DIMENSION                                              TLEP 550
      IF(N-1)2,1,3                                                      TLEP 560
C                                                                       TLEP 570
C        DIMENSION LESS THAN 2                                          TLEP 580
    1 POL(1)=C(1)                                                       TLEP 590
    2 RETURN                                                            TLEP 600
C                                                                       TLEP 610
    3 POL(1)=C(1)+B*C(2)                                                TLEP 620
      POL(2)=A*C(2)                                                     TLEP 630
      IF(N-2)2,2,4                                                      TLEP 640
C                                                                       TLEP 650
C        INITIALIZATION                                                 TLEP 660
    4 WORK(1)=1.                                                        TLEP 670
      WORK(2)=B                                                         TLEP 680
      WORK(3)=0.                                                        TLEP 690
      WORK(4)=A                                                         TLEP 700
      FI=1.                                                             TLEP 710
C                                                                       TLEP 720
C        CALCULATE COEFFICIENT VECTOR OF NEXT LEGENDRE POLYNOMIAL       TLEP 730
C        AND ADD MULTIPLE OF THIS VECTOR TO POLYNOMIAL POL              TLEP 740
      DO 6 J=3,N                                                        TLEP 750
      FI=FI+1.                                                          TLEP 760
      Q=1./FI-1.                                                        TLEP 770
      Q1=1.-Q                                                           TLEP 780
      P=0.                                                              TLEP 790
C                                                                       TLEP 800
      DO 5 K=2,J                                                        TLEP 810
      H=(A*P+B*WORK(2*K-2))*Q1+Q*WORK(2*K-3)                            TLEP 820
      P=WORK(2*K-2)                                                     TLEP 830
      WORK(2*K-2)=H                                                     TLEP 840
      WORK(2*K-3)=P                                                     TLEP 850
    5 POL(K-1)=POL(K-1)+H*C(J)                                          TLEP 860
      WORK(2*J-1)=0.                                                    TLEP 870
      WORK(2*J)=A*P*Q1                                                  TLEP 880
    6 POL(J)=C(J)*WORK(2*J)                                             TLEP 890
      RETURN                                                            TLEP 900
      END                                                               TLEP 910
