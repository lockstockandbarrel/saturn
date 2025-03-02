C                                                                       BESY  10
C     ..................................................................BESY  20
C                                                                       BESY  30
C        SUBROUTINE BESY                                                BESY  40
C                                                                       BESY  50
C        PURPOSE                                                        BESY  60
C           COMPUTE THE Y BESSEL FUNCTION FOR A GIVEN ARGUMENT AND ORDERBESY  70
C                                                                       BESY  80
C        USAGE                                                          BESY  90
C           CALL BESY(X,N,BY,IER)                                       BESY 100
C                                                                       BESY 110
C        DESCRIPTION OF PARAMETERS                                      BESY 120
C           X  -THE ARGUMENT OF THE Y BESSEL FUNCTION DESIRED           BESY 130
C           N  -THE ORDER OF THE Y BESSEL FUNCTION DESIRED              BESY 140
C           BY -THE RESULTANT Y BESSEL FUNCTION                         BESY 150
C           IER-RESULTANT ERROR CODE WHERE                              BESY 160
C              IER=0  NO ERROR                                          BESY 170
C              IER=1  N IS NEGATIVE                                     BESY 180
C              IER=2  X IS NEGATIVE OR ZERO                             BESY 190
C              IER=3  BY HAS EXCEEDED MAGNITUDE OF 10**70               BESY 200
C                                                                       BESY 210
C        REMARKS                                                        BESY 220
C           VERY SMALL VALUES OF X MAY CAUSE THE RANGE OF THE LIBRARY   BESY 230
C           FUNCTION ALOG TO BE EXCEEDED                                BESY 240
C           X MUST BE GREATER THAN ZERO                                 BESY 250
C           N MUST BE GREATER THAN OR EQUAL TO ZERO                     BESY 260
C                                                                       BESY 270
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  BESY 280
C           NONE                                                        BESY 290
C                                                                       BESY 300
C        METHOD                                                         BESY 310
C           RECURRENCE RELATION AND POLYNOMIAL APPROXIMATION TECHNIQUE  BESY 320
C           AS DESCRIBED BY A.J.M.HITCHCOCK,'POLYNOMIAL APPROXIMATIONS  BESY 330
C           TO BESSEL FUNCTIONS OF ORDER ZERO AND ONE AND TO RELATED    BESY 340
C           FUNCTIONS', M.T.A.C., V.11,1957,PP.86-88, AND G.N. WATSON,  BESY 350
C           'A TREATISE ON THE THEORY OF BESSEL FUNCTIONS', CAMBRIDGE   BESY 360
C           UNIVERSITY PRESS, 1958, P. 62                               BESY 370
C                                                                       BESY 380
C     ..................................................................BESY 390
C                                                                       BESY 400
      SUBROUTINE BESY(X,N,BY,IER)                                       BESY 410
C                                                                       BESY 420
C     CHECK FOR ERRORS IN N AND X                                       BESY 430
C                                                                       BESY 440
      IF(N)180,10,10                                                    BESY 450
   10 IER=0                                                             BESY 460
      IF(X)190,190,20                                                   BESY 470
C                                                                       BESY 480
C     BRANCH IF X LESS THAN OR EQUAL 4                                  BESY 490
C                                                                       BESY 500
   20 IF(X-4.0)40,40,30                                                 BESY 510
C                                                                       BESY 520
C       COMPUTE Y0 AND Y1 FOR X GREATER THAN 4                          BESY 530
C                                                                       BESY 540
   30 T1=4.0/X                                                          BESY 550
      T2=T1*T1                                                          BESY 560
      P0=((((-.0000037043*T2+.0000173565)*T2-.0000487613)*T2            BESY 570
     1  +.00017343)*T2-.001753062)*T2+.3989423                          BESY 580
      Q0=((((.0000032312*T2-.0000142078)*T2+.0000342468)*T2             BESY 590
     1  -.0000869791)*T2+.0004564324)*T2-.01246694                      BESY 600
      P1=((((.0000042414*T2-.0000200920)*T2+.0000580759)*T2             BESY 610
     1  -.000223203)*T2+.002921826)*T2+.3989423                         BESY 620
      Q1=((((-.0000036594*T2+.00001622)*T2-.0000398708)*T2              BESY 630
     1  +.0001064741)*T2-.0006390400)*T2+.03740084                      BESY 640
      A=2.0/SQRT(X)                                                     BESY 650
      B=A*T1                                                            BESY 660
      C=X-.7853982                                                      BESY 670
      Y0=A*P0*SIN(C)+B*Q0*COS(C)                                        BESY 680
      Y1=-A*P1*COS(C)+B*Q1*SIN(C)                                       BESY 690
      GO TO 90                                                          BESY 700
C                                                                       BESY 710
C       COMPUTE Y0 AND Y1 FOR X LESS THAN OR EQUAL TO 4                 BESY 720
C                                                                       BESY 730
   40 XX=X/2.                                                           BESY 740
      X2=XX*XX                                                          BESY 750
      T=ALOG(XX)+.5772157                                               BESY 760
      SUM=0.                                                            BESY 770
      TERM=T                                                            BESY 780
      Y0=T                                                              BESY 790
      DO 70 L=1,15                                                      BESY 800
      IF(L-1)50,60,50                                                   BESY 810
   50 SUM=SUM+1./FLOAT(L-1)                                             BESY 820
   60 FL=L                                                              BESY 830
      TS=T-SUM                                                          BESY 840
      TERM=(TERM*(-X2)/FL**2)*(1.-1./(FL*TS))                           BESY 850
   70 Y0=Y0+TERM                                                        BESY 860
      TERM = XX*(T-.5)                                                  BESY 870
      SUM=0.                                                            BESY 880
      Y1=TERM                                                           BESY 890
      DO 80 L=2,16                                                      BESY 900
      SUM=SUM+1./FLOAT(L-1)                                             BESY 910
      FL=L                                                              BESY 920
      FL1=FL-1.                                                         BESY 930
      TS=T-SUM                                                          BESY 940
      TERM=(TERM*(-X2)/(FL1*FL))*((TS-.5/FL)/(TS+.5/FL1))               BESY 950
   80 Y1=Y1+TERM                                                        BESY 960
      PI2=.6366198                                                      BESY 970
      Y0=PI2*Y0                                                         BESY 980
      Y1=-PI2/X+PI2*Y1                                                  BESY 990
C                                                                       BESY1000
C     CHECK IF ONLY Y0 OR Y1 IS DESIRED                                 BESY1010
C                                                                       BESY1020
   90 IF(N-1)100,100,130                                                BESY1030
C                                                                       BESY1040
C     RETURN EITHER Y0 OR Y1 AS REQUIRED                                BESY1050
C                                                                       BESY1060
  100 IF(N)110,120,110                                                  BESY1070
  110 BY=Y1                                                             BESY1080
      GO TO 170                                                         BESY1090
  120 BY=Y0                                                             BESY1100
      GO TO 170                                                         BESY1110
C                                                                       BESY1120
C    PERFORM RECURRENCE OPERATIONS TO FIND YN(X)                        BESY1130
C                                                                       BESY1140
  130 YA=Y0                                                             BESY1150
      YB=Y1                                                             BESY1160
      K=1                                                               BESY1170
  140 T=FLOAT(2*K)/X                                                    BESY1180
      YC=T*YB-YA                                                        BESY1190
      IF(ABS(YC)-1.0E70)145,145,141                                     BESY1200
  141 IER=3                                                             BESY1210
      RETURN                                                            BESY1220
  145 K=K+1                                                             BESY1230
      IF(K-N)150,160,150                                                BESY1240
  150 YA=YB                                                             BESY1250
      YB=YC                                                             BESY1260
      GO TO 140                                                         BESY1270
  160 BY=YC                                                             BESY1280
  170 RETURN                                                            BESY1290
  180 IER=1                                                             BESY1300
      RETURN                                                            BESY1310
  190 IER=2                                                             BESY1320
      RETURN                                                            BESY1330
      END                                                               BESY1340
