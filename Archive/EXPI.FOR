C                                                                       EXPI  10
C     ..................................................................EXPI  20
C                                                                       EXPI  30
C        SUBROUTINE EXPI                                                EXPI  40
C                                                                       EXPI  50
C        PURPOSE                                                        EXPI  60
C           COMPUTES THE EXPONENTIAL INTEGRAL -EI(-X)                   EXPI  70
C                                                                       EXPI  80
C        USAGE                                                          EXPI  90
C           CALL EXPI(X,RES)                                            EXPI 100
C                                                                       EXPI 110
C        DESCRIPTION OF PARAMETERS                                      EXPI 120
C           X     - ARGUMENT OF EXPONENTIAL INTEGRAL                    EXPI 130
C           RES   - RESULT VALUE                                        EXPI 140
C           AUX   - RESULTANT AUXILIARY VALUE                           EXPI 150
C                                                                       EXPI 160
C        REMARKS                                                        EXPI 170
C           X GT 170 (X LT -174) MAY CAUSE UNDERFLOW (OVERFLOW)         EXPI 180
C           WITH THE EXPONENTIAL FUNCTION                               EXPI 190
C           FOR X = 0 THE RESULT VALUE IS SET TO -1.E75                 EXPI 200
C                                                                       EXPI 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  EXPI 220
C           NONE                                                        EXPI 230
C                                                                       EXPI 240
C        METHOD                                                         EXPI 250
C           DEFINITION                                                  EXPI 260
C           RES=INTEGRAL(EXP(-T)/T, SUMMED OVER T FROM X TO INFINITY).  EXPI 270
C           EVALUATION                                                  EXPI 280
C           THREE DIFFERENT RATIONAL APPROXIMATIONS ARE USED IN THE     EXPI 290
C           RANGES 1 LE X, X LE -9 AND -9 LT X LE -3 RESPECTIVELY,      EXPI 300
C           A POLYNOMIAL APPROXIMATION IS USED IN -3 LT X LT 1.         EXPI 310
C                                                                       EXPI 320
C     ..................................................................EXPI 330
C                                                                       EXPI 340
      SUBROUTINE EXPI(X,RES,AUX)                                        EXPI 350
      IF(X-1.)2,1,1                                                     EXPI 360
    1 Y=1./X                                                            EXPI 370
      AUX=1.-Y*(((Y+3.377358E0)*Y+2.052156E0)*Y+2.709479E-1)/((((Y*     EXPI 380
     11.072553E0+5.716943E0)*Y+6.945239E0)*Y+2.593888E0)*Y+2.709496E-1) EXPI 390
      RES=AUX*Y*EXP(-X)                                                 EXPI 400
      RETURN                                                            EXPI 410
    2 IF(X+3.)6,6,3                                                     EXPI 420
    3 AUX=(((((((7.122452E-7*X-1.766345E-6)*X+2.928433E-5)*X-2.335379E-4EXPI 430
     1)*X+1.664156E-3)*X-1.041576E-2)*X+5.555682E-2)*X-2.500001E-1)*X   EXPI 440
     2+9.999999E-1                                                      EXPI 450
      RES=-1.E75                                                        EXPI 460
      IF(X)4,5,4                                                        EXPI 470
    4 RES=X*AUX-ALOG(ABS(X))-5.772157E-1                                EXPI 480
    5 RETURN                                                            EXPI 490
    6 IF(X+9.)8,8,7                                                     EXPI 500
    7 AUX=1.-((((5.176245E-2*X+3.061037E0)*X+3.243665E1)*X+2.244234E2)*XEXPI 510
     1+2.486697E2)/((((X+3.995161E0)*X+3.893944E1)*X+2.263818E1)*X      EXPI 520
     2+1.807837E2)                                                      EXPI 530
      GOTO 9                                                            EXPI 540
    8 Y=9./X                                                            EXPI 550
      AUX=1.-Y*(((Y+7.659824E-1)*Y-7.271015E-1)*Y-1.080693E0)/((((Y     EXPI 560
     1*2.518750E0+1.122927E1)*Y+5.921405E0)*Y-8.666702E0)*Y-9.724216E0) EXPI 570
    9 RES=AUX*EXP(-X)/X                                                 EXPI 580
      RETURN                                                            EXPI 590
      END                                                               EXPI 600
