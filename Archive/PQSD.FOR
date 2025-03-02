C                                                                       PQSD  10
C     ..................................................................PQSD  20
C                                                                       PQSD  30
C        SUBROUTINE PQSD                                                PQSD  40
C                                                                       PQSD  50
C        PURPOSE                                                        PQSD  60
C           PERFORM QUADRATIC SYNTHETIC DIVISION                        PQSD  70
C                                                                       PQSD  80
C        USAGE                                                          PQSD  90
C           CALL PQSD(A,B,P,Q,X,IDIMX)                                  PQSD 100
C                                                                       PQSD 110
C        DESCRIPTION OF PARAMETERS                                      PQSD 120
C           A     - COEFFICIENT OF Z IN REMAINDER (CALCULATED)          PQSD 130
C           B     - CONSTANT TERM IN REMAINDER (CALCULATED)             PQSD 140
C           P     - COEFFICIENT OF Z IN QUADRATIC POLYNOMIAL            PQSD 150
C           Q     - CONSTANT TERM IN QUADRATIC POLYNOMIAL               PQSD 160
C           X     - COEFFICIENT VECTOR FOR GIVEN POLYNOMIAL, ORDERED    PQSD 170
C                   FROM SMALLEST TO LARGEST POWER                      PQSD 180
C           IDIMX - DIMENSION OF X                                      PQSD 190
C                                                                       PQSD 200
C        REMARKS                                                        PQSD 210
C           NONE                                                        PQSD 220
C                                                                       PQSD 230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PQSD 240
C           NONE                                                        PQSD 250
C                                                                       PQSD 260
C        METHOD                                                         PQSD 270
C           POLYNOMIAL IS DIVIDED BY THE QUADRATIC Z**2-P*Z-Q GIVING    PQSD 280
C           THE LINEAR REMAINDER A*Z+B                                  PQSD 290
C                                                                       PQSD 300
C     ..................................................................PQSD 310
C                                                                       PQSD 320
      SUBROUTINE PQSD(A,B,P,Q,X,IDIMX)                                  PQSD 330
      DIMENSION X(1)                                                    PQSD 340
C                                                                       PQSD 350
      A=0.                                                              PQSD 360
      B=0.                                                              PQSD 370
      J=IDIMX                                                           PQSD 380
    1 IF(J)3,3,2                                                        PQSD 390
    2 Z=P*A+B                                                           PQSD 400
      B=Q*A+X(J)                                                        PQSD 410
      A=Z                                                               PQSD 420
      J=J-1                                                             PQSD 430
      GO TO 1                                                           PQSD 440
    3 RETURN                                                            PQSD 450
      END                                                               PQSD 460
