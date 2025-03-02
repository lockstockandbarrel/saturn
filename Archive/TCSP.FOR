C                                                                       TCSP  10
C     ..................................................................TCSP  20
C                                                                       TCSP  30
C        SUBROUTINE TCSP                                                TCSP  40
C                                                                       TCSP  50
C        PURPOSE                                                        TCSP  60
C           A SERIES EXPANSION IN SHIFTED CHEBYSHEV POLYNOMIALS WITH    TCSP  70
C           INDEPENDENT VARIABLE X IS TRANSFORMED TO A POLYNOMIAL WITH  TCSP  80
C           INDEPENDENT VARIABLE Z, WHERE X=A*Z+B.                      TCSP  90
C                                                                       TCSP 100
C        USAGE                                                          TCSP 110
C           CALL TCSP(A,B,POL,N,C,WORK)                                 TCSP 120
C                                                                       TCSP 130
C        DESCRIPTION OF PARAMETERS                                      TCSP 140
C           A     - FACTOR OF LINEAR TERM IN GIVEN LINEAR TRANSFORMATIONTCSP 150
C           B     - CONSTANT TERM IN GIVEN LINEAR TRANSFORMATION        TCSP 160
C           POL   - COEFFICIENT VECTOR OF POLYNOMIAL (RESULTANT VALUE)  TCSP 170
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           TCSP 180
C           N     - DIMENSION OF COEFFICIENT VECTORS POL AND C          TCSP 190
C           C     - GIVEN COEFFICIENT VECTOR OF EXPANSION               TCSP 200
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH           TCSP 210
C                   POL AND C MAY BE IDENTICALLY LOCATED                TCSP 220
C           WORK  - WORKING STORAGE OF DIMENSION 2*N                    TCSP 230
C                                                                       TCSP 240
C        REMARKS                                                        TCSP 250
C           COEFFICIENT VECTOR C REMAINS UNCHANGED IF NOT COINCIDING    TCSP 260
C           WITH COEFFICIENT VECTOR POL.                                TCSP 270
C           OPERATION IS BYPASSED IN CASE N LESS THAN 1.                TCSP 280
C           THE LINEAR TRANSFORMATION X=A*Z+B OR Z=(1/A)(X-B) TRANSFORMSTCSP 290
C           THE RANGE (0,1) IN X TO THE RANGE (ZL,ZR) IN Z, WHERE       TCSP 300
C           ZL=-B/A AND ZR=(1-B)/A.                                     TCSP 310
C           FOR GIVEN ZL, ZR WE HAVE A=1/(ZR-ZL) AND B=-ZL/(ZR-ZL).     TCSP 320
C                                                                       TCSP 330
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  TCSP 340
C           NONE                                                        TCSP 350
C                                                                       TCSP 360
C        METHOD                                                         TCSP 370
C           THE TRANSFORMATION IS BASED ON THE RECURRENCE EQUATION FOR  TCSP 380
C           SHIFTED CHEBYSHEV POLYNOMIALS TS(N,X)                       TCSP 390
C           TS(N+1,X)=(4*X-2)*TS(N,X)-TS(N-1,X),                        TCSP 400
C           WHERE THE FIRST TERM IN BRACKETS IS THE INDEX,              TCSP 410
C           THE SECOND IS THE ARGUMENT.                                 TCSP 420
C           STARTING VALUES ARE TS(0,X)=1, TS(1,X)=2*X-1.               TCSP 430
C           THE TRANSFORMATION IS IMPLICITLY DEFINED BY MEANS OF        TCSP 440
C           X=A*Z+B TOGETHER WITH                                       TCSP 450
C           SUM(POL(I)*Z**(I-1), SUMMED OVER I FROM 1 TO N)             TCSP 460
C           =SUM(C(I)*TS(I-1,X), SUMMED OVER I FROM 1 TO N).            TCSP 470
C                                                                       TCSP 480
C     ..................................................................TCSP 490
C                                                                       TCSP 500
      SUBROUTINE TCSP(A,B,POL,N,C,WORK)                                 TCSP 510
C                                                                       TCSP 520
      DIMENSION POL(1),C(1),WORK(1)                                     TCSP 530
C                                                                       TCSP 540
C        TEST OF DIMENSION                                              TCSP 550
      IF(N-1)2,1,3                                                      TCSP 560
C                                                                       TCSP 570
C        DIMENSION LESS THAN 2                                          TCSP 580
    1 POL(1)=C(1)                                                       TCSP 590
    2 RETURN                                                            TCSP 600
C                                                                       TCSP 610
    3 XD=A+A                                                            TCSP 620
      X0=B+B-1.                                                         TCSP 630
      POL(1)=C(1)+C(2)*X0                                               TCSP 640
      POL(2)=C(2)*XD                                                    TCSP 650
      IF(N-2)2,2,4                                                      TCSP 660
C                                                                       TCSP 670
C        INITIALIZATION                                                 TCSP 680
    4 WORK(1)=1.                                                        TCSP 690
      WORK(2)=X0                                                        TCSP 700
      WORK(3)=0.                                                        TCSP 710
      WORK(4)=XD                                                        TCSP 720
      XD=XD+XD                                                          TCSP 730
      X0=X0+X0                                                          TCSP 740
C                                                                       TCSP 750
C        CALCULATE COEFFICIENT VECTOR OF NEXT SHIFTED CHEBYSHEV         TCSP 760
C        POLYNOMIAL AND ADD MULTIPLE OF THIS VECTOR TO POLYNOMIAL POL   TCSP 770
      DO 6 J=3,N                                                        TCSP 780
      P=0.                                                              TCSP 790
C                                                                       TCSP 800
      DO 5 K=2,J                                                        TCSP 810
      H=P-WORK(2*K-3)+X0*WORK(2*K-2)                                    TCSP 820
      P=WORK(2*K-2)                                                     TCSP 830
      WORK(2*K-2)=H                                                     TCSP 840
      WORK(2*K-3)=P                                                     TCSP 850
      POL(K-1)=POL(K-1)+H*C(J)                                          TCSP 860
    5 P=XD*P                                                            TCSP 870
      WORK(2*J-1)=0.                                                    TCSP 880
      WORK(2*J)=P                                                       TCSP 890
    6 POL(J)=C(J)*P                                                     TCSP 900
      RETURN                                                            TCSP 910
      END                                                               TCSP 920
