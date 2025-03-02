C                                                                       ELI1  10
C     ..................................................................ELI1  20
C                                                                       ELI1  30
C        SUBROUTINE ELI1                                                ELI1  40
C                                                                       ELI1  50
C        PURPOSE                                                        ELI1  60
C           COMPUTES THE ELLIPTIC INTEGRAL OF FIRST KIND                ELI1  70
C                                                                       ELI1  80
C        USAGE                                                          ELI1  90
C           CALL ELI1(RES,X,CK)                                         ELI1 100
C                                                                       ELI1 110
C        DESCRIPTION OF PARAMETERS                                      ELI1 120
C           RES   - RESULT VALUE                                        ELI1 130
C           X     - UPPER INTEGRATION BOUND (ARGUMENT OF ELLIPTIC       ELI1 140
C                   INTEGRAL OF FIRST KIND)                             ELI1 150
C           CK    - COMPLEMENTARY MODULUS                               ELI1 160
C                                                                       ELI1 170
C        REMARKS                                                        ELI1 180
C           MODULUS K = SQRT(1.-CK*CK).                                 ELI1 190
C                                                                       ELI1 200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  ELI1 210
C           NONE                                                        ELI1 220
C                                                                       ELI1 230
C        METHOD                                                         ELI1 240
C           DEFINITION                                                  ELI1 250
C           RES=INTEGRAL(1/SQRT((1+T*T)*(1+(CK*T)**2)), SUMMED          ELI1 260
C           OVER T FROM 0 TO X).                                        ELI1 270
C           EQUIVALENT ARE THE DEFINITIONS                              ELI1 280
C           RES=INTEGRAL(1/(COS(T)*SQRT(1+(CK*TAN(T))**2)), SUMMED      ELI1 290
C           OVER T FROM 0 TO ATAN(X)),                                  ELI1 300
C           RES=INTEGRAL(1/SQRT(1-(K*SIN(T))**2), SUMMED OVER           ELI1 310
C           T FROM 0 TO ATAN(X)).                                       ELI1 320
C           EVALUATION                                                  ELI1 330
C           LANDENS TRANSFORMATION IS USED FOR CALCULATION.             ELI1 340
C           REFERENCE                                                   ELI1 350
C           R. BULIRSCH, NUMERICAL CALCULATION OF ELLIPTIC INTEGRALS ANDELI1 360
C                  ELLIPTIC FUNCTIONS.                                  ELI1 370
C                  HANDBOOK SERIES OF SPECIAL FUNCTIONS                 ELI1 380
C                  NUMERISCHE MATHEMATIK VOL. 7, 1965, PP. 78-90.       ELI1 390
C                                                                       ELI1 400
C     ..................................................................ELI1 410
C                                                                       ELI1 420
      SUBROUTINE ELI1(RES,X,CK)                                         ELI1 430
C                                                                       ELI1 440
      IF(X)2,1,2                                                        ELI1 450
    1 RES=0.                                                            ELI1 460
      RETURN                                                            ELI1 470
    2 IF(CK)4,3,4                                                       ELI1 480
    3 RES=ALOG(ABS(X)+SQRT(1.+X*X))                                     ELI1 490
      GOTO 13                                                           ELI1 500
    4 ANGLE=ABS(1./X)                                                   ELI1 510
      GEO=ABS(CK)                                                       ELI1 520
      ARI=1.                                                            ELI1 530
      PIM=0.                                                            ELI1 540
    5 SQGEO=ARI*GEO                                                     ELI1 550
      AARI=ARI                                                          ELI1 560
      ARI=GEO+ARI                                                       ELI1 570
      ANGLE=-SQGEO/ANGLE+ANGLE                                          ELI1 580
      SQGEO=SQRT(SQGEO)                                                 ELI1 590
      IF(ANGLE)7,6,7                                                    ELI1 600
C     REPLACE 0 BY SMALL VALUE                                          ELI1 610
    6 ANGLE=SQGEO*1.E-8                                                 ELI1 620
    7 TEST=AARI*1.E-4                                                   ELI1 630
      IF(ABS(AARI-GEO)-TEST)10,10,8                                     ELI1 640
    8 GEO=SQGEO+SQGEO                                                   ELI1 650
      PIM=PIM+PIM                                                       ELI1 660
      IF(ANGLE)9,5,5                                                    ELI1 670
    9 PIM=PIM+3.1415927                                                 ELI1 680
      GOTO 5                                                            ELI1 690
   10 IF(ANGLE)11,12,12                                                 ELI1 700
   11 PIM=PIM+3.1415927                                                 ELI1 710
   12 RES=(ATAN(ARI/ANGLE)+PIM)/ARI                                     ELI1 720
   13 IF(X)14,15,15                                                     ELI1 730
   14 RES=-RES                                                          ELI1 740
   15 RETURN                                                            ELI1 750
      END                                                               ELI1 760
