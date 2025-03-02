C                                                                       DTCS  10
C     ..................................................................DTCS  20
C                                                                       DTCS  30
C        SUBROUTINE DTCSP                                               DTCS  40
C                                                                       DTCS  50
C        PURPOSE                                                        DTCS  60
C           A SERIES EXPANSION IN SHIFTED CHEBYSHEV POLYNOMIALS WITH    DTCS  70
C           INDEPENDENT VARIABLE X IS TRANSFORMED TO A POLYNOMIAL WITH  DTCS  80
C           INDEPENDENT VARIABLE Z, WHERE X=A*Z+B.                      DTCS  90
C                                                                       DTCS 100
C        USAGE                                                          DTCS 110
C           CALL DTCSP(A,B,POL,N,C,WORK)                                DTCS 120
C                                                                       DTCS 130
C        DESCRIPTION OF PARAMETERS                                      DTCS 140
C           A     - FACTOR OF LINEAR TERM IN GIVEN LINEAR TRANSFORMATIONDTCS 150
C                   DOUBLE PRECISION VARIABLE                           DTCS 160
C           B     - CONSTANT TERM IN GIVEN LINEAR TRANSFORMATION        DTCS 170
C                   DOUBLE PRECISION VARIABLE                           DTCS 180
C           POL   - COEFFICIENT VECTOR OF POLYNOMIAL (RESULTANT VALUE)  DTCS 190
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DTCS 200
C                   DOUBLE PRECISION VECTOR                             DTCS 210
C           N     - DIMENSION OF COEFFICIENT VECTORS POL AND C          DTCS 220
C           C     - GIVEN COEFFICIENT VECTOR OF EXPANSION               DTCS 230
C                   POL AND C MAY BE IDENTICALLY LOCATED                DTCS 240
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           DTCS 250
C                   DOUBLE PRECISION VECTOR                             DTCS 260
C           WORK  - WORKING STORAGE OF DIMENSION 2*N                    DTCS 270
C                   DOUBLE PRECISION ARRAY                              DTCS 280
C                                                                       DTCS 290
C        REMARKS                                                        DTCS 300
C           COEFFICIENT VECTOR C REMAINS UNCHANGED IF NOT COINCIDING    DTCS 310
C           WITH COEFFICIENT VECTOR POL.                                DTCS 320
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1.                DTCS 330
C           THE LINEAR TRANSFORMATION X=A*Z+B OR Z=(1/A)(X-B) TRANSFORMSDTCS 340
C           THE RANGE (0,1) IN X TO THE RANGE (ZL,ZR) IN Z, WHERE       DTCS 350
C           ZL=-B/A AND ZR=(1-B)/A.                                     DTCS 360
C           FOR GIVEN ZL, ZR WE HAVE A=1/(ZR-ZL) AND B=-ZL/(ZR-ZL).     DTCS 370
C                                                                       DTCS 380
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DTCS 390
C           NONE                                                        DTCS 400
C                                                                       DTCS 410
C        METHOD                                                         DTCS 420
C           THE TRANSFORMATION IS BASED ON THE RECURRENCE EQUATION FOR  DTCS 430
C           SHIFTED CHEBYSHEV POLYNOMIALS TS(N,X)                       DTCS 440
C           TS(N+1,X)=(4*X-2)*TS(N,X)-TS(N-1,X),                        DTCS 450
C           WHERE THE FIRST TERM IN BRACKETS IS THE INDEX,              DTCS 460
C           THE SECOND IS THE ARGUMENT.                                 DTCS 470
C           STARTING VALUES ARE TS(0,X)=1, TS(1,X)=2*X-1.               DTCS 480
C           THE TRANSFORMATION IS IMPLICITLY DEFINED BY MEANS OF        DTCS 490
C           X=A*Z+B TOGETHER WITH                                       DTCS 500
C           SUM(POL(I)*Z**(I-1), SUMMED OVER I FROM 1 TO N)             DTCS 510
C           =SUM(C(I)*TS(I-1,X), SUMMED OVER I FROM 1 TO N).            DTCS 520
C                                                                       DTCS 530
C     ..................................................................DTCS 540
C                                                                       DTCS 550
      SUBROUTINE DTCSP(A,B,POL,N,C,WORK)                                DTCS 560
C                                                                       DTCS 570
      DIMENSION POL(1),C(1),WORK(1)                                     DTCS 580
      DOUBLE PRECISION A,B,POL,C,WORK,H,P,XD,X0                         DTCS 590
C                                                                       DTCS 600
C        TEST OF DIMENSION                                              DTCS 610
      IF(N-1)2,1,3                                                      DTCS 620
C                                                                       DTCS 630
C        DIMENSION LESS THAN 2                                          DTCS 640
    1 POL(1)=C(1)                                                       DTCS 650
    2 RETURN                                                            DTCS 660
C                                                                       DTCS 670
    3 XD=A+A                                                            DTCS 680
      X0=B+B-1.D0                                                       DTCS 690
      POL(1)=C(1)+C(2)*X0                                               DTCS 700
      POL(2)=C(2)*XD                                                    DTCS 710
      IF(N-2)2,2,4                                                      DTCS 720
C                                                                       DTCS 730
C        INITIALIZATION                                                 DTCS 740
    4 WORK(1)=1.D0                                                      DTCS 750
      WORK(2)=X0                                                        DTCS 760
      WORK(3)=0.D0                                                      DTCS 770
      WORK(4)=XD                                                        DTCS 780
      XD=XD+XD                                                          DTCS 790
      X0=X0+X0                                                          DTCS 800
C                                                                       DTCS 810
C        CALCULATE COEFFICIENT VECTOR OF NEXT SHIFTED CHEBYSHEV         DTCS 820
C        POLYNOMIAL AND ADD MULTIPLE OF THIS VECTOR TO POLYNOMIAL POL   DTCS 830
      DO 6 J=3,N                                                        DTCS 840
      P=0.D0                                                            DTCS 850
C                                                                       DTCS 860
      DO 5 K=2,J                                                        DTCS 870
      H=P-WORK(2*K-3)+X0*WORK(2*K-2)                                    DTCS 880
      P=WORK(2*K-2)                                                     DTCS 890
      WORK(2*K-2)=H                                                     DTCS 900
      WORK(2*K-3)=P                                                     DTCS 910
      POL(K-1)=POL(K-1)+H*C(J)                                          DTCS 920
    5 P=XD*P                                                            DTCS 930
      WORK(2*J-1)=0.D0                                                  DTCS 940
      WORK(2*J)=P                                                       DTCS 950
    6 POL(J)=C(J)*P                                                     DTCS 960
      RETURN                                                            DTCS 970
      END                                                               DTCS 980
