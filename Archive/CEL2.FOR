C                                                                       CEL2  10
C     ..................................................................CEL2  20
C                                                                       CEL2  30
C        SUBROUTINE CEL2                                                CEL2  40
C                                                                       CEL2  50
C        PURPOSE                                                        CEL2  60
C           COMPUTES THE GENERALIZED COMPLETE ELLIPTIC INTEGRAL OF      CEL2  70
C           SECOND KIND.                                                CEL2  80
C                                                                       CEL2  90
C        USAGE                                                          CEL2 100
C           CALL CEL2(RES,AK,A,B,IER)                                   CEL2 110
C                                                                       CEL2 120
C        DESCRIPTION OF PARAMETERS                                      CEL2 130
C           RES   - RESULT VALUE                                        CEL2 140
C           AK    - MODULUS (INPUT)                                     CEL2 150
C           A     - CONSTANT TERM IN NUMERATOR                          CEL2 160
C           B     - FACTOR OF QUADRATIC TERM IN NUMERATOR               CEL2 170
C           IER   - RESULTANT ERROR CODE WHERE                          CEL2 180
C                   IER=0  NO ERROR                                     CEL2 190
C                   IER=1  AK NOT IN RANGE -1 TO +1                     CEL2 200
C                                                                       CEL2 210
C        REMARKS                                                        CEL2 220
C           FOR ABS(AK) GE 1 THE RESULT IS SET TO 1.E75 IF B IS         CEL2 230
C           POSITIVE, TO -1.E75 IF B IS NEGATIVE.                       CEL2 240
C           SPECIAL CASES ARE                                           CEL2 250
C           K(K) OBTAINED WITH A = 1, B = 1                             CEL2 260
C           E(K) OBTAINED WITH A = 1, B = CK*CK WHERE CK IS             CEL2 270
C           COMPLEMENTARY MODULUS.                                      CEL2 280
C           B(K) OBTAINED WITH A = 1, B = 0                             CEL2 290
C           D(K) OBTAINED WITH A = 0, B = 1                             CEL2 300
C           WHERE K, E, B, D DEFINE SPECIAL CASES OF THE GENERALIZED    CEL2 310
C           COMPLETE ELLIPTIC INTEGRAL OF SECOND KIND IN THE USUAL      CEL2 320
C           NOTATION, AND THE ARGUMENT K OF THESE FUNCTIONS MEANS       CEL2 330
C           THE MODULUS.                                                CEL2 340
C                                                                       CEL2 350
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CEL2 360
C           NONE                                                        CEL2 370
C                                                                       CEL2 380
C        METHOD                                                         CEL2 390
C           DEFINITION                                                  CEL2 400
C           RES=INTEGRAL((A+B*T*T)/(SQRT((1+T*T)*(1+(CK*T)**2))*(1+T*T))CEL2 410
C           SUMMED OVER T FROM 0 TO INFINITY).                          CEL2 420
C           EVALUATION                                                  CEL2 430
C           LANDENS TRANSFORMATION IS USED FOR CALCULATION.             CEL2 440
C           REFERENCE                                                   CEL2 450
C           R.BULIRSCH, 'NUMERICAL CALCULATION OF ELLIPTIC INTEGRALS    CEL2 460
C           AND ELLIPTIC FUNCTIONS', HANDBOOK SERIES SPECIAL FUNCTIONS, CEL2 470
C           NUMERISCHE MATHEMATIK VOL. 7, 1965, PP. 78-90.              CEL2 480
C                                                                       CEL2 490
C     ..................................................................CEL2 500
C                                                                       CEL2 510
      SUBROUTINE CEL2(RES,AK,A,B,IER)                                   CEL2 520
      IER=0                                                             CEL2 530
      ARI=2.                                                            CEL2 540
      GEO=(0.5-AK)+0.5                                                  CEL2 550
      GEO=GEO+GEO*AK                                                    CEL2 560
      RES=A                                                             CEL2 570
      A1=A+B                                                            CEL2 580
      B0=B+B                                                            CEL2 590
      IF(GEO)1,2,6                                                      CEL2 600
    1 IER=1                                                             CEL2 610
    2 IF(B)3,8,4                                                        CEL2 620
    3 RES=-1.E75                                                        CEL2 630
      RETURN                                                            CEL2 640
    4 RES=1.E75                                                         CEL2 650
      RETURN                                                            CEL2 660
    5 GEO=GEO*AARI                                                      CEL2 670
    6 GEO=SQRT(GEO)                                                     CEL2 680
      GEO=GEO+GEO                                                       CEL2 690
      AARI=ARI                                                          CEL2 700
      ARI=ARI+GEO                                                       CEL2 710
      B0=B0+RES*GEO                                                     CEL2 720
      RES=A1                                                            CEL2 730
      B0=B0+B0                                                          CEL2 740
      A1=B0/ARI+A1                                                      CEL2 750
      IF(GEO/AARI-0.9999)5,7,7                                          CEL2 760
    7 RES=A1/ARI                                                        CEL2 770
      RES=RES+0.5707963E0*RES                                           CEL2 780
    8 RETURN                                                            CEL2 790
      END                                                               CEL2 800
