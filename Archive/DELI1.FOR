C                                                                       DEL1  10
C     ..................................................................DEL1  20
C                                                                       DEL1  30
C        SUBROUTINE DELI1                                               DEL1  40
C                                                                       DEL1  50
C        PURPOSE                                                        DEL1  60
C           COMPUTES THE ELLIPTIC INTEGRAL OF FIRST KIND                DEL1  70
C                                                                       DEL1  80
C        USAGE                                                          DEL1  90
C           CALL DELI1(RES,X,CK)                                        DEL1 100
C                                                                       DEL1 110
C        DESCRIPTION OF PARAMETERS                                      DEL1 120
C           RES   - RESULT VALUE IN DOUBLE PRECISION                    DEL1 130
C           X     - UPPER INTEGRATION BOUND (ARGUMENT OF ELLIPTIC       DEL1 140
C                   INTEGRAL OF FIRST KIND) IN DOUBLE PRECISION         DEL1 150
C           CK    - COMPLEMENTARY MODULUS IN DOUBLE PRECISION           DEL1 160
C                                                                       DEL1 170
C        REMARKS                                                        DEL1 180
C           DOUBLE PRECISION MODULUS K = DSQRT(1.D0-CK*CK).             DEL1 190
C                                                                       DEL1 200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DEL1 210
C           NONE                                                        DEL1 220
C                                                                       DEL1 230
C        METHOD                                                         DEL1 240
C           DEFINITION                                                  DEL1 250
C           RES=INTEGRAL(1/SQRT((1+T*T)*(1+(CK*T)**2)), SUMMED          DEL1 260
C           OVER T FROM 0 TO X).                                        DEL1 270
C           EQUIVALENT ARE THE DEFINITIONS                              DEL1 280
C           RES=INTEGRAL(1/(COS(T)*SQRT(1+(CK*TAN(T))**2)), SUMMED      DEL1 290
C           OVER T FROM 0 TO ATAN(X)),                                  DEL1 300
C           RES=INTEGRAL(1/SQRT(1-(K*SIN(T))**2), SUMMED OVER           DEL1 310
C           T FROM 0 TO ATAN(X)).                                       DEL1 320
C           EVALUATION                                                  DEL1 330
C           LANDENS TRANSFORMATION IS USED FOR CALCULATION.             DEL1 340
C           REFERENCE                                                   DEL1 350
C           R. BULIRSCH, NUMERICAL CALCULATION OF ELLIPTIC INTEGRALS ANDDEL1 360
C                  ELLIPTIC FUNCTIONS.                                  DEL1 370
C                  HANDBOOK SERIES OF SPECIAL FUNCTIONS                 DEL1 380
C                  NUMERISCHE MATHEMATIK VOL. 7, 1965, PP. 78-90.       DEL1 390
C                                                                       DEL1 400
C     ..................................................................DEL1 410
C                                                                       DEL1 420
      SUBROUTINE DELI1(RES,X,CK)                                        DEL1 430
C                                                                       DEL1 440
      DOUBLE PRECISION RES,X,CK,ANGLE,GEO,ARI,PIM,SQGEO,AARI,TEST       DEL1 450
C                                                                       DEL1 460
      IF(X)2,1,2                                                        DEL1 470
    1 RES=0.D0                                                          DEL1 480
      RETURN                                                            DEL1 490
C                                                                       DEL1 500
    2 IF(CK)4,3,4                                                       DEL1 510
    3 RES=DLOG(DABS(X)+DSQRT(1.D0+X*X))                                 DEL1 520
      GOTO 13                                                           DEL1 530
C                                                                       DEL1 540
    4 ANGLE=DABS(1.D0/X)                                                DEL1 550
      GEO=DABS(CK)                                                      DEL1 560
      ARI=1.D0                                                          DEL1 570
      PIM=0.D0                                                          DEL1 580
    5 SQGEO=ARI*GEO                                                     DEL1 590
      AARI=ARI                                                          DEL1 600
      ARI=GEO+ARI                                                       DEL1 610
      ANGLE=-SQGEO/ANGLE+ANGLE                                          DEL1 620
      SQGEO=DSQRT(SQGEO)                                                DEL1 630
      IF(ANGLE)7,6,7                                                    DEL1 640
C                                                                       DEL1 650
C        REPLACE 0 BY SMALL VALUE                                       DEL1 660
C                                                                       DEL1 670
    6 ANGLE=SQGEO*1.D-17                                                DEL1 680
    7 TEST=AARI*1.D-9                                                   DEL1 690
      IF(DABS(AARI-GEO)-TEST)10,10,8                                    DEL1 700
    8 GEO=SQGEO+SQGEO                                                   DEL1 710
      PIM=PIM+PIM                                                       DEL1 720
      IF(ANGLE)9,5,5                                                    DEL1 730
    9 PIM=PIM+3.1415926535897932                                        DEL1 740
      GOTO 5                                                            DEL1 750
   10 IF(ANGLE)11,12,12                                                 DEL1 760
   11 PIM=PIM+3.1415926535897932                                        DEL1 770
   12 RES=(DATAN(ARI/ANGLE)+PIM)/ARI                                    DEL1 780
   13 IF(X)14,15,15                                                     DEL1 790
   14 RES=-RES                                                          DEL1 800
   15 RETURN                                                            DEL1 810
      END                                                               DEL1 820
