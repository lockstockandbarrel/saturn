C                                                                       PGCD  10
C     ..................................................................PGCD  20
C                                                                       PGCD  30
C        SUBROUTINE PGCD                                                PGCD  40
C                                                                       PGCD  50
C        PURPOSE                                                        PGCD  60
C           DETERMINE GREATEST COMMON DIVISOR OF TWO POLYNOMIALS        PGCD  70
C                                                                       PGCD  80
C        USAGE                                                          PGCD  90
C           CALL PGCD(X,IDIMX,Y,IDIMY,WORK,EPS,IER)                     PGCD 100
C                                                                       PGCD 110
C        DESCRIPTION OF PARAMETERS                                      PGCD 120
C           X     -  VECTOR OF COEFFICIENTS FOR FIRST POLYNOMIAL,       PGCD 130
C                    ORDERED FROM SMALLEST TO LARGEST POWER             PGCD 140
C           IDIMX -  DIMENSION OF X                                     PGCD 150
C           Y     -  VECTOR OF COEFFICIENTS FOR SECOND POLYNOMIAL,      PGCD 160
C                    ORDERED FROM SMALLEST TO LARGEST POWER.            PGCD 170
C                    THIS IS REPLACED BY GREATEST COMMON DIVISOR        PGCD 180
C           IDIMY -  DIMENSION OF Y                                     PGCD 190
C           WORK  -  WORKING STORAGE ARRAY                              PGCD 200
C           EPS   -  TOLERANCE VALUE BELOW WHICH COEFFICIENT IS         PGCD 210
C                    ELIMINATED DURING NORMALIZATION                    PGCD 220
C           IER   -  RESULTANT ERROR CODE WHERE                         PGCD 230
C                    IER=0  NO ERROR                                    PGCD 240
C                    IER=1  X OR Y IS ZERO POLYNOMIAL                   PGCD 250
C                                                                       PGCD 260
C        REMARKS                                                        PGCD 270
C           IDIMX MUST BE GREATER THAN IDIMY                            PGCD 280
C           IDIMY=1 ON RETURN MEANS X AND Y ARE PRIME, THE GCD IS A     PGCD 290
C           CONSTANT. IDIMX IS DESTROYED DURING COMPUTATION.            PGCD 300
C                                                                       PGCD 310
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PGCD 320
C           PDIV                                                        PGCD 330
C           PNORM                                                       PGCD 340
C                                                                       PGCD 350
C        METHOD                                                         PGCD 360
C           GREATEST COMMON DIVISOR OF TWO POLYNOMIALS X AND Y IS       PGCD 370
C           DETERMINED BY MEANS OF EUCLIDEAN ALGORITHM. COEFFICIENT     PGCD 380
C           VECTORS X AND Y ARE DESTROYED AND GREATEST COMMON           PGCD 390
C           DIVISOR IS GENERATED IN Y.                                  PGCD 400
C                                                                       PGCD 410
C     ..................................................................PGCD 420
C                                                                       PGCD 430
      SUBROUTINE PGCD(X,IDIMX,Y,IDIMY,WORK,EPS,IER)                     PGCD 440
      DIMENSION X(1),Y(1),WORK(1)                                       PGCD 450
C                                                                       PGCD 460
C     DIMENSION REQUIRED FOR VECTOR NAMED  WORK  IS   IDIMX-IDIMY+1     PGCD 470
C                                                                       PGCD 480
    1 CALL PDIV(WORK,NDIM,X,IDIMX,Y,IDIMY,EPS,IER)                      PGCD 490
      IF(IER) 5,2,5                                                     PGCD 500
    2 IF(IDIMX) 5,5,3                                                   PGCD 510
C                                                                       PGCD 520
C     INTERCHANGE X AND Y                                               PGCD 530
C                                                                       PGCD 540
    3 DO 4 J=1,IDIMY                                                    PGCD 550
      WORK(1)=X(J)                                                      PGCD 560
      X(J)=Y(J)                                                         PGCD 570
    4 Y(J)=WORK(1)                                                      PGCD 580
      NDIM=IDIMX                                                        PGCD 590
      IDIMX=IDIMY                                                       PGCD 600
      IDIMY=NDIM                                                        PGCD 610
      GO TO 1                                                           PGCD 620
    5 RETURN                                                            PGCD 630
      END                                                               PGCD 640
