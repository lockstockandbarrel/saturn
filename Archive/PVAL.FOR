C                                                                       PVAL  10
C     ..................................................................PVAL  20
C                                                                       PVAL  30
C        SUBROUTINE PVAL                                                PVAL  40
C                                                                       PVAL  50
C        PURPOSE                                                        PVAL  60
C           EVALUATE A POLYNOMIAL FOR A GIVEN VALUE OF THE VARIABLE     PVAL  70
C                                                                       PVAL  80
C        USAGE                                                          PVAL  90
C           CALL PVAL(RES,ARG,X,IDIMX)                                  PVAL 100
C                                                                       PVAL 110
C        DESCRIPTION OF PARAMETERS                                      PVAL 120
C           RES    - RESULTANT VALUE OF POLYNOMIAL                      PVAL 130
C           ARG    - GIVEN VALUE OF THE VARIABLE                        PVAL 140
C           X      - VECTOR OF COEFFICIENTS, ORDERED FROM SMALLEST TO   PVAL 150
C                    LARGEST POWER                                      PVAL 160
C           IDIMX  - DIMENSION OF X                                     PVAL 170
C                                                                       PVAL 180
C        REMARKS                                                        PVAL 190
C           NONE                                                        PVAL 200
C                                                                       PVAL 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PVAL 220
C           NONE                                                        PVAL 230
C                                                                       PVAL 240
C        METHOD                                                         PVAL 250
C           EVALUATION IS DONE BY MEANS OF NESTED MULTIPLICATION        PVAL 260
C                                                                       PVAL 270
C     ..................................................................PVAL 280
C                                                                       PVAL 290
      SUBROUTINE PVAL(RES,ARG,X,IDIMX)                                  PVAL 300
      DIMENSION X(1)                                                    PVAL 310
C                                                                       PVAL 320
      RES=0.                                                            PVAL 330
      J=IDIMX                                                           PVAL 340
    1 IF(J)3,3,2                                                        PVAL 350
    2 RES=RES*ARG+X(J)                                                  PVAL 360
      J=J-1                                                             PVAL 370
      GO TO 1                                                           PVAL 380
    3 RETURN                                                            PVAL 390
      END                                                               PVAL 400
