C                                                                       CEL1  10
C     ..................................................................CEL1  20
C                                                                       CEL1  30
C        SUBROUTINE CEL1                                                CEL1  40
C                                                                       CEL1  50
C        PURPOSE                                                        CEL1  60
C           CALCULATE COMPLETE ELLIPTIC INTEGRAL OF FIRST KIND          CEL1  70
C                                                                       CEL1  80
C        USAGE                                                          CEL1  90
C           CALL CEL1(RES,AK,IER)                                       CEL1 100
C                                                                       CEL1 110
C        DESCRIPTION OF PARAMETERS                                      CEL1 120
C           RES   - RESULT VALUE                                        CEL1 130
C           AK    - MODULUS (INPUT)                                     CEL1 140
C           IER   - RESULTANT ERROR CODE WHERE                          CEL1 150
C                   IER=0  NO ERROR                                     CEL1 160
C                   IER=1  AK NOT IN RANGE -1 TO +1                     CEL1 170
C                                                                       CEL1 180
C        REMARKS                                                        CEL1 190
C           THE RESULT IS SET TO 1.E75 IF ABS(AK) GE 1                  CEL1 200
C           FOR MODULUS AK AND COMPLEMENTARY MODULUS CK,                CEL1 210
C           EQUATION AK*AK+CK*CK=1.0 IS USED.                           CEL1 220
C           AK MUST BE IN THE RANGE -1 TO +1                            CEL1 230
C                                                                       CEL1 240
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CEL1 250
C           NONE                                                        CEL1 260
C                                                                       CEL1 270
C        METHOD                                                         CEL1 280
C           DEFINITION                                                  CEL1 290
C           CEL1(AK)=INTEGRAL(1/SQRT((1+T*T)*(1+(CK*T)**2)), SUMMED     CEL1 300
C           OVER T FROM 0 TO INFINITY).                                 CEL1 310
C           EQUIVALENT ARE THE DEFINITIONS                              CEL1 320
C           CEL1(AK)=INTEGRAL(1/(COS(T)SQRT(1+(CK*TAN(T))**2)),SUMMED   CEL1 330
C           OVER T FROM 0 TO PI/2),                                     CEL1 340
C           CEL1(AK)=INTEGRAL(1/SQRT(1-(AK*SIN(T))**2),SUMMED OVER T    CEL1 350
C           FROM 0 TO PI/2), WHERE K=SQRT(1.-CK*CK).                    CEL1 360
C           EVALUATION                                                  CEL1 370
C           LANDENS TRANSFORMATION IS USED FOR CALCULATION.             CEL1 380
C           REFERENCE                                                   CEL1 390
C           R.BULIRSCH, 'NUMERICAL CALCULATION OF ELLIPTIC INTEGRALS    CEL1 400
C           AND ELLIPTIC FUNCTIONS', HANDBOOK SERIES SPECIAL FUNCTIONS, CEL1 410
C           NUMERISCHE MATHEMATIK VOL. 7, 1965, PP. 78-90.              CEL1 420
C                                                                       CEL1 430
C     ..................................................................CEL1 440
C                                                                       CEL1 450
      SUBROUTINE CEL1(RES,AK,IER)                                       CEL1 460
      IER=0                                                             CEL1 470
      ARI=2.                                                            CEL1 480
      GEO=(0.5-AK)+0.5                                                  CEL1 490
      GEO=GEO+GEO*AK                                                    CEL1 500
      RES=0.5                                                           CEL1 510
      IF(GEO)1,2,4                                                      CEL1 520
    1 IER=1                                                             CEL1 530
    2 RES=1.E75                                                         CEL1 540
      RETURN                                                            CEL1 550
    3 GEO=GEO*AARI                                                      CEL1 560
    4 GEO=SQRT(GEO)                                                     CEL1 570
      GEO=GEO+GEO                                                       CEL1 580
      AARI=ARI                                                          CEL1 590
      ARI=ARI+GEO                                                       CEL1 600
      RES=RES+RES                                                       CEL1 610
      IF(GEO/AARI-0.9999)3,5,5                                          CEL1 620
    5 RES=RES/ARI*6.283185E0                                            CEL1 630
      RETURN                                                            CEL1 640
      END                                                               CEL1 650
