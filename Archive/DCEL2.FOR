C                                                                       DCE2  10
C     ..................................................................DCE2  20
C                                                                       DCE2  30
C        SUBROUTINE DCEL2                                               DCE2  40
C                                                                       DCE2  50
C        PURPOSE                                                        DCE2  60
C           COMPUTES THE GENERALIZED COMPLETE ELLIPTIC INTEGRAL OF      DCE2  70
C           SECOND KIND.                                                DCE2  80
C                                                                       DCE2  90
C        USAGE                                                          DCE2 100
C           CALL DCEL2(RES,AK,A,B,IER)                                  DCE2 110
C                                                                       DCE2 120
C        DESCRIPTION OF PARAMETERS                                      DCE2 130
C           RES   - RESULT VALUE IN DOUBLE PRECISION                    DCE2 140
C           AK    - MODULUS (INPUT) IN DOUBLE PRECISION                 DCE2 150
C           A     - DOUBLE PRECISION CONSTANT TERM IN NUMERATOR         DCE2 160
C           B     - DOUBLE PRECISION FACTOR OF QUADRATIC TERM           DCE2 170
C                   IN NUMERATOR                                        DCE2 180
C           IER   - RESULTANT ERROR CODE WHERE                          DCE2 190
C                   IER=0  NO ERROR                                     DCE2 200
C                   IER=1  AK NOT IN RANGE -1 TO +1                     DCE2 210
C                                                                       DCE2 220
C        REMARKS                                                        DCE2 230
C           FOR ABS(AK) GE 1 THE RESULT IS SET TO 1.E75 IF B IS         DCE2 240
C           POSITIVE, TO -1.D75 IF B IS NEGATIVE.                       DCE2 250
C           SPECIAL CASES ARE                                           DCE2 260
C           K(K) OBTAINED WITH A = 1, B = 1                             DCE2 270
C           E(K) OBTAINED WITH A = 1, B = CK*CK WHERE CK IS             DCE2 280
C           COMPLEMENTARY MODULUS.                                      DCE2 290
C           B(K) OBTAINED WITH A = 1, B = 0                             DCE2 300
C           D(K) OBTAINED WITH A = 0, B = 1                             DCE2 310
C           WHERE K, E, B, D DEFINE SPECIAL CASES OF THE GENERALIZED    DCE2 320
C           COMPLETE ELLIPTIC INTEGRAL OF SECOND KIND IN THE USUAL      DCE2 330
C           NOTATION, AND THE ARGUMENT K OF THESE FUNCTIONS MEANS       DCE2 340
C           THE MODULUS.                                                DCE2 350
C                                                                       DCE2 360
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DCE2 370
C           NONE                                                        DCE2 380
C                                                                       DCE2 390
C        METHOD                                                         DCE2 400
C           DEFINITION                                                  DCE2 410
C           RES=INTEGRAL((A+B*T*T)/(SQRT((1+T*T)*(1+(CK*T)**2))*(1+T*T))DCE2 420
C           SUMMED OVER T FROM 0 TO INFINITY).                          DCE2 430
C           EVALUATION                                                  DCE2 440
C           LANDENS TRANSFORMATION IS USED FOR CALCULATION.             DCE2 450
C           REFERENCE                                                   DCE2 460
C           R.BULIRSCH, 'NUMERICAL CALCULATION OF ELLIPTIC INTEGRALS    DCE2 470
C           AND ELLIPTIC FUNCTIONS', HANDBOOK SERIES SPECIAL FUNCTIONS, DCE2 480
C           NUMERISCHE MATHEMATIK VOL. 7, 1965, PP. 78-90.              DCE2 490
C                                                                       DCE2 500
C     ..................................................................DCE2 510
C                                                                       DCE2 520
      SUBROUTINE DCEL2(RES,AK,A,B,IER)                                  DCE2 530
      DOUBLE PRECISION RES,AK,A,B,GEO,ARI,AARI,B0,A1                    DCE2 540
      IER=0                                                             DCE2 550
      ARI=2.D0                                                          DCE2 560
      GEO=(0.5D0-AK)+0.5D0                                              DCE2 570
      GEO=GEO+GEO*AK                                                    DCE2 580
      RES=A                                                             DCE2 590
      A1=A+B                                                            DCE2 600
      B0=B+B                                                            DCE2 610
      IF(GEO)1,2,6                                                      DCE2 620
    1 IER=1                                                             DCE2 630
    2 IF(B)3,8,4                                                        DCE2 640
    3 RES=-1.D75                                                        DCE2 650
      RETURN                                                            DCE2 660
    4 RES=1.D75                                                         DCE2 670
      RETURN                                                            DCE2 680
    5 GEO=GEO*AARI                                                      DCE2 690
    6 GEO=DSQRT(GEO)                                                    DCE2 700
      GEO=GEO+GEO                                                       DCE2 710
      AARI=ARI                                                          DCE2 720
      ARI=ARI+GEO                                                       DCE2 730
      B0=B0+RES*GEO                                                     DCE2 740
      RES=A1                                                            DCE2 750
      B0=B0+B0                                                          DCE2 760
      A1=B0/ARI+A1                                                      DCE2 770
      IF(GEO/AARI-0.999999995D0)5,7,7                                   DCE2 780
    7 RES=A1/ARI                                                        DCE2 790
      RES=RES+0.57079632679489662D0*RES                                 DCE2 800
    8 RETURN                                                            DCE2 810
      END                                                               DCE2 820
