C                                                                       BESI  10
C     ..................................................................BESI  20
C                                                                       BESI  30
C        SUBROUTINE BESI                                                BESI  40
C                                                                       BESI  50
C        PURPOSE                                                        BESI  60
C           COMPUTE THE I BESSEL FUNCTION FOR A GIVEN ARGUMENT AND ORDERBESI  70
C                                                                       BESI  80
C        USAGE                                                          BESI  90
C           CALL BESI(X,N,BI,IER)                                       BESI 100
C                                                                       BESI 110
C        DESCRIPTION OF PARAMETERS                                      BESI 120
C           X  -THE ARGUMENT OF THE I BESSEL FUNCTION DESIRED           BESI 130
C           N  -THE ORDER OF THE I BESSEL FUNCTION DESIRED              BESI 140
C           BI -THE RESULTANT I BESSEL FUNCTION                         BESI 150
C           IER-RESULTANT ERROR CODE WHERE                              BESI 160
C              IER=0 NO ERROR                                           BESI 170
C              IER=1 N IS NEGATIVE                                      BESI 180
C              IER=2 X IS NEGATIVE                                      BESI 190
C              IER=3 UNDERFLOW, BI .LT. 1.E-69, BI SET TO 0.0           BESI 200
C              IER=4 OVERFLOW, X .GT. 170 WHERE X .GT. N                BESI 210
C                                                                       BESI 220
C        REMARKS                                                        BESI 230
C           N AND X MUST BE .GE. ZERO                                   BESI 240
C                                                                       BESI 250
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  BESI 260
C           NONE                                                        BESI 270
C                                                                       BESI 280
C        METHOD                                                         BESI 290
C           COMPUTES I BESSEL FUNCTION USING SERIES OR ASYMPTOTIC       BESI 300
C           APPROXIMATION DEPENDING ON RANGE OF ARGUMENTS.              BESI 310
C                                                                       BESI 320
C     ..................................................................BESI 330
C                                                                       BESI 340
      SUBROUTINE BESI(X,N, BI,IER)                                      BESI 350
C                                                                       BESI 360
C     CHECK FOR ERRORS IN N AND X AND EXIT IF ANY ARE PRESENT           BESI 370
C                                                                       BESI 380
      IER=0                                                             BESI 390
      BI=1.0                                                            BESI 400
      IF(N)150,15,10                                                    BESI 410
   10 IF(X)160,20,20                                                    BESI 420
   15 IF(X)160,17,20                                                    BESI 430
   17 RETURN                                                            BESI 440
C                                                                       BESI 450
C     DEFINE TOLERANCE                                                  BESI 460
C                                                                       BESI 470
   20 TOL=1.E-6                                                         BESI 480
C                                                                       BESI 490
C     IF ARGUMENT GT 12 AND GT N, USE ASYMPTOTIC FORM                   BESI 500
C                                                                       BESI 510
      IF(X-12.)40,40,30                                                 BESI 520
   30 IF(X-FLOAT(N))40,40,110                                           BESI 530
C                                                                       BESI 540
C     COMPUTE FIRST TERM OF SERIES AND SET INITIAL VALUE OF THE SUM     BESI 550
C                                                                       BESI 560
   40 XX=X/2.                                                           BESI 570
   50 TERM=1.0                                                          BESI 580
      IF(N) 70,70,55                                                    BESI 590
   55 DO 60 I=1,N                                                       BESI 600
      FI=I                                                              BESI 610
      IF(ABS(TERM)-1.E-68)56,60,60                                      BESI 620
   56 IER=3                                                             BESI 630
      BI=0.0                                                            BESI 640
      RETURN                                                            BESI 650
   60 TERM=TERM*XX/FI                                                   BESI 660
   70 BI=TERM                                                           BESI 670
      XX=XX*XX                                                          BESI 680
C                                                                       BESI 690
C     COMPUTE TERMS, STOPPING WHEN ABS(TERM) LE ABS(SUM OF TERMS)       BESI 700
C     TIMES TOLERANCE                                                   BESI 710
C                                                                       BESI 720
      DO 90 K=1,1000                                                    BESI 730
      IF(ABS(TERM)-ABS(BI*TOL))100,100,80                               BESI 740
   80 FK=K*(N+K)                                                        BESI 750
      TERM=TERM*(XX/FK)                                                 BESI 760
   90 BI=BI+TERM                                                        BESI 770
C                                                                       BESI 780
C     RETURN BI AS ANSWER                                               BESI 790
C                                                                       BESI 800
  100 RETURN                                                            BESI 810
C                                                                       BESI 820
C     X GT 12 AND X GT N, SO USE ASYMPTOTIC APPROXIMATION               BESI 830
C                                                                       BESI 840
  110 FN=4*N*N                                                          BESI 850
      IF(X-170.0)115,111,111                                            BESI 860
  111 IER=4                                                             BESI 870
      RETURN                                                            BESI 880
  115 XX=1./(8.*X)                                                      BESI 890
      TERM=1.                                                           BESI 900
      BI=1.                                                             BESI 910
      DO 130 K=1,30                                                     BESI 920
      IF(ABS(TERM)-ABS(TOL*BI))140,140,120                              BESI 930
  120 FK=(2*K-1)**2                                                     BESI 940
      TERM=TERM*XX*(FK-FN)/FLOAT(K)                                     BESI 950
  130 BI=BI+TERM                                                        BESI 960
C                                                                       BESI 970
C     SIGNIFICANCE LOST AFTER 30 TERMS, TRY SERIES                      BESI 980
C                                                                       BESI 990
      GO TO 40                                                          BESI1000
  140 PI=3.141592653                                                    BESI1010
      BI=BI*EXP(X)/SQRT(2.*PI*X)                                        BESI1020
      GO TO 100                                                         BESI1030
  150 IER=1                                                             BESI1040
      GO TO 100                                                         BESI1050
  160 IER=2                                                             BESI1060
      GO TO 100                                                         BESI1070
      END                                                               BESI1080
