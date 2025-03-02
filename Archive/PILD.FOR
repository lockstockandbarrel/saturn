C                                                                       PILD  10
C     ..................................................................PILD  20
C                                                                       PILD  30
C        SUBROUTINE PILD                                                PILD  40
C                                                                       PILD  50
C        PURPOSE                                                        PILD  60
C           EVALUATE POLYNOMIAL AND ITS FIRST DERIVATIVE FOR A GIVEN    PILD  70
C           ARGUMENT                                                    PILD  80
C                                                                       PILD  90
C        USAGE                                                          PILD 100
C           CALL PILD(POLY,DVAL,ARGUM,X,IDIMX)                          PILD 110
C                                                                       PILD 120
C        DESCRIPTION OF PARAMETERS                                      PILD 130
C           POLY  - VALUE OF POLYNOMIAL                                 PILD 140
C           DVAL  - DERIVATIVE                                          PILD 150
C           ARGUM - ARGUMENT                                            PILD 160
C           X     - VECTOR OF COEFFICIENTS FOR POLYNOMIAL, ORDERED      PILD 170
C                   FROM SMALLEST TO LARGEST POWER                      PILD 180
C           IDIMX - DIMENSION OF X                                      PILD 190
C                                                                       PILD 200
C        REMARKS                                                        PILD 210
C           NONE                                                        PILD 220
C                                                                       PILD 230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PILD 240
C           PQSD                                                        PILD 250
C                                                                       PILD 260
C        METHOD                                                         PILD 270
C           EVALUATION IS DONE BY MEANS OF SUBROUTINE PQSD (QUADRATIC   PILD 280
C           SYNTHETIC DIVISION)                                         PILD 290
C                                                                       PILD 300
C     ..................................................................PILD 310
C                                                                       PILD 320
      SUBROUTINE PILD (POLY,DVAL,ARGUM,X,IDIMX)                         PILD 330
      DIMENSION X(1)                                                    PILD 340
C                                                                       PILD 350
      P=ARGUM+ARGUM                                                     PILD 360
      Q=-ARGUM*ARGUM                                                    PILD 370
C                                                                       PILD 380
      CALL PQSD (DVAL,POLY,P,Q,X,IDIMX)                                 PILD 390
C                                                                       PILD 400
      POLY=ARGUM*DVAL+POLY                                              PILD 410
C                                                                       PILD 420
      RETURN                                                            PILD 430
      END                                                               PILD 440
