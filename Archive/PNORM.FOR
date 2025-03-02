C                                                                       PNOR  10
C     ..................................................................PNOR  20
C                                                                       PNOR  30
C        SUBROUTINE PNORM                                               PNOR  40
C                                                                       PNOR  50
C        PURPOSE                                                        PNOR  60
C           NORMALIZE COEFFICIENT VECTOR OF A POLYNOMIAL                PNOR  70
C                                                                       PNOR  80
C        USAGE                                                          PNOR  90
C           CALL PNORM(X,IDIMX,EPS)                                     PNOR 100
C                                                                       PNOR 110
C        DESCRIPTION OF PARAMETERS                                      PNOR 120
C           X      - VECTOR OF ORIGINAL COEFFICIENTS, ORDERED FROM      PNOR 130
C                    SMALLEST TO LARGEST POWER. IT REMAINS UNCHANGED    PNOR 140
C           IDIMX  - DIMENSION OF X. IT IS REPLACED BY FINAL DIMENSION  PNOR 150
C           EPS    - TOLERANCE BELOW WHICH COEFFICIENT IS ELIMINATED    PNOR 160
C                                                                       PNOR 170
C        REMARKS                                                        PNOR 180
C           IF ALL COEFFICIENTS ARE LESS THAN EPS, RESULT IS A ZERO     PNOR 190
C           POLYNOMIAL WITH IDIMX=0 BUT VECTOR X REMAINS INTACT         PNOR 200
C                                                                       PNOR 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PNOR 220
C           NONE                                                        PNOR 230
C                                                                       PNOR 240
C        METHOD                                                         PNOR 250
C           DIMENSION OF VECTOR X IS REDUCED BY ONE FOR EACH TRAILING   PNOR 260
C           COEFFICIENT WITH AN ABSOLUTE VALUE LESS THAN OR EQUAL TO EPSPNOR 270
C                                                                       PNOR 280
C     ..................................................................PNOR 290
C                                                                       PNOR 300
      SUBROUTINE PNORM(X,IDIMX,EPS)                                     PNOR 310
      DIMENSION X(1)                                                    PNOR 320
C                                                                       PNOR 330
    1 IF(IDIMX) 4,4,2                                                   PNOR 340
    2 IF(ABS(X(IDIMX))-EPS) 3,3,4                                       PNOR 350
    3 IDIMX=IDIMX-1                                                     PNOR 360
      GO TO 1                                                           PNOR 370
    4 RETURN                                                            PNOR 380
      END                                                               PNOR 390
