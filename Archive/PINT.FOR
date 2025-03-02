C                                                                       PINT  10
C     ..................................................................PINT  20
C                                                                       PINT  30
C        SUBROUTINE PINT                                                PINT  40
C                                                                       PINT  50
C        PURPOSE                                                        PINT  60
C           FIND INTEGRAL OF A POLYNOMIAL WITH CONSTANT OF INTEGRATION  PINT  70
C           EQUAL TO ZERO                                               PINT  80
C                                                                       PINT  90
C        USAGE                                                          PINT 100
C           CALL PINT(Y,IDIMY,X,IDIMX)                                  PINT 110
C                                                                       PINT 120
C        DESCRIPTION OF PARAMETERS                                      PINT 130
C           Y     - VECTOR OF COEFFICIENTS FOR INTEGRAL, ORDERED FROM   PINT 140
C                   SMALLEST TO LARGEST POWER                           PINT 150
C           IDIMY - DIMENSION OF Y (EQUAL TO IDIMX+1)                   PINT 160
C           X     - VECTOR OF COEFFICIENTS FOR ORIGINAL POLYNOMIAL,     PINT 170
C                   ORDERED FROM SMALLEST TO LARGEST POWER              PINT 180
C           IDIMX - DIMENSION OF X                                      PINT 190
C                                                                       PINT 200
C        REMARKS                                                        PINT 210
C           NONE                                                        PINT 220
C                                                                       PINT 230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PINT 240
C           NONE                                                        PINT 250
C                                                                       PINT 260
C        METHOD                                                         PINT 270
C           DIMENSION OF Y IS SET AT DIMENSION OF X PLUS ONE, AND THE   PINT 280
C           CONSTANT TERM IS SET TO ZERO. INTEGRAL IS THEN CALCULATED   PINT 290
C           BY DIVIDING COEFFICIENTS BY THEIR RESPECTIVE EXPONENTS.     PINT 300
C                                                                       PINT 310
C     ..................................................................PINT 320
C                                                                       PINT 330
      SUBROUTINE PINT(Y,IDIMY,X,IDIMX)                                  PINT 340
      DIMENSION X(1),Y(1)                                               PINT 350
C                                                                       PINT 360
      IDIMY=IDIMX+1                                                     PINT 370
      Y(1)=0.                                                           PINT 380
      IF(IDIMX)1,1,2                                                    PINT 390
    1 RETURN                                                            PINT 400
    2 EXPT=1.                                                           PINT 410
      DO 3 I=2,IDIMY                                                    PINT 420
      Y(I)=X(I-1)/EXPT                                                  PINT 430
    3 EXPT=EXPT+1.                                                      PINT 440
      GO TO 1                                                           PINT 450
      END                                                               PINT 460
