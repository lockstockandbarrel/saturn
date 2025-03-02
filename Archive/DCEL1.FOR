C                                                                       DCE1  10
C     ..................................................................DCE1  20
C                                                                       DCE1  30
C        SUBROUTINE DCEL1                                               DCE1  40
C                                                                       DCE1  50
C        PURPOSE                                                        DCE1  60
C           CALCULATE COMPLETE ELLIPTIC INTEGRAL OF FIRST KIND          DCE1  70
C                                                                       DCE1  80
C        USAGE                                                          DCE1  90
C           CALL DCEL1(RES,AK,IER)                                      DCE1 100
C                                                                       DCE1 110
C        DESCRIPTION OF PARAMETERS                                      DCE1 120
C           RES   - RESULT VALUE IN DOUBLE PRECISION                    DCE1 130
C           AK    - MODULUS (INPUT) IN DOUBLE PRECISION                 DCE1 140
C           IER   - RESULTANT ERROR CODE WHERE                          DCE1 150
C                   IER=0  NO ERROR                                     DCE1 160
C                   IER=1  AK NOT IN RANGE -1 TO +1                     DCE1 170
C                                                                       DCE1 180
C        REMARKS                                                        DCE1 190
C           THE RESULT IS SET TO 1.E75 IF ABS(AK) GE 1                  DCE1 200
C           FOR MODULUS AK AND COMPLEMENTARY MODULUS CK,                DCE1 210
C           EQUATION AK*AK+CK*CK=1.D0 IS USED.                          DCE1 220
C           AK MUST BE IN THE RANGE -1 TO +1                            DCE1 230
C                                                                       DCE1 240
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DCE1 250
C           NONE                                                        DCE1 260
C                                                                       DCE1 270
C        METHOD                                                         DCE1 280
C           DEFINITION                                                  DCE1 290
C           CEL1(AK)=INTEGRAL(1/SQRT((1+T*T)*(1+(CK*T)**2)), SUMMED     DCE1 300
C           OVER T FROM 0 TO INFINITY).                                 DCE1 310
C           EQUIVALENT ARE THE DEFINITIONS                              DCE1 320
C           CEL1(AK)=INTEGRAL(1/(COS(T)SQRT(1+(CK*TAN(T))**2)),SUMMED   DCE1 330
C           OVER T FROM 0 TO PI/2),                                     DCE1 340
C           CEL1(AK)=INTEGRAL(1/SQRT(1-(AK*SIN(T))**2),SUMMED OVER T    DCE1 350
C           FROM 0 TO PI/2), WHERE K=SQRT(1.-CK*CK).                    DCE1 360
C           EVALUATION                                                  DCE1 370
C           LANDENS TRANSFORMATION IS USED FOR CALCULATION.             DCE1 380
C           REFERENCE                                                   DCE1 390
C           R.BULIRSCH, 'NUMERICAL CALCULATION OF ELLIPTIC INTEGRALS    DCE1 400
C           AND ELLIPTIC FUNCTIONS', HANDBOOK SERIES SPECIAL FUNCTIONS, DCE1 410
C           NUMERISCHE MATHEMATIK VOL. 7, 1965, PP. 78-90.              DCE1 420
C                                                                       DCE1 430
C     ..................................................................DCE1 440
C                                                                       DCE1 450
      SUBROUTINE DCEL1(RES,AK,IER)                                      DCE1 460
      DOUBLE PRECISION RES,AK,GEO,ARI,AARI                              DCE1 470
      IER=0                                                             DCE1 480
      ARI=2.D0                                                          DCE1 490
      GEO=(0.5D0-AK)+0.5D0                                              DCE1 500
      GEO=GEO+GEO*AK                                                    DCE1 510
      RES=0.5D0                                                         DCE1 520
      IF(GEO)1,2,4                                                      DCE1 530
    1 IER=1                                                             DCE1 540
    2 RES=1.D75                                                         DCE1 550
      RETURN                                                            DCE1 560
    3 GEO=GEO*AARI                                                      DCE1 570
    4 GEO=DSQRT(GEO)                                                    DCE1 580
      GEO=GEO+GEO                                                       DCE1 590
      AARI=ARI                                                          DCE1 600
      ARI=ARI+GEO                                                       DCE1 610
      RES=RES+RES                                                       DCE1 620
      IF(GEO/AARI-0.999999995D0)3,5,5                                   DCE1 630
    5 RES=RES/ARI*6.2831853071795865D0                                  DCE1 640
      RETURN                                                            DCE1 650
      END                                                               DCE1 660
