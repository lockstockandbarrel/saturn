C                                                                       PCLD  10
C     ..................................................................PCLD  20
C                                                                       PCLD  30
C        SUBROUTINE PCLD                                                PCLD  40
C                                                                       PCLD  50
C        PURPOSE                                                        PCLD  60
C           SHIFT OF ORIGIN (COMPLETE LINEAR SYNTHETIC DIVISION)        PCLD  70
C                                                                       PCLD  80
C        USAGE                                                          PCLD  90
C           CALL PCLD(X,IDIMX,U)                                        PCLD 100
C                                                                       PCLD 110
C        DESCRIPTION OF PARAMETERS                                      PCLD 120
C           X     - VECTOR OF COEFFICIENTS, ORDERED FROM SMALLEST TO    PCLD 130
C                   LARGEST POWER. IT IS REPLACED BY VECTOR OF          PCLD 140
C                   TRANSFORMED COEFFICIENTS.                           PCLD 150
C           IDIMX - DIMENSION OF X                                      PCLD 160
C           U     - SHIFT PARAMETER                                     PCLD 170
C                                                                       PCLD 180
C        REMARKS                                                        PCLD 190
C           NONE                                                        PCLD 200
C                                                                       PCLD 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PCLD 220
C           NONE                                                        PCLD 230
C                                                                       PCLD 240
C        METHOD                                                         PCLD 250
C           COEFFICIENT VECTOR X(I) OF POLYNOMIAL P(Z) IS TRANSFORMED   PCLD 260
C           SUCH THAT Q(Z)=P(Z-U) WHERE Q(Z) DENOTES THE POLYNOMIAL     PCLD 270
C           WITH TRANSFORMED COEFFICIENT VECTOR.                        PCLD 280
C                                                                       PCLD 290
C     ..................................................................PCLD 300
C                                                                       PCLD 310
      SUBROUTINE PCLD (X,IDIMX,U)                                       PCLD 320
      DIMENSION X(1)                                                    PCLD 330
C                                                                       PCLD 340
      K=1                                                               PCLD 350
    1 J=IDIMX                                                           PCLD 360
    2 IF (J-K) 4,4,3                                                    PCLD 370
    3 X(J-1)=X(J-1)+U*X(J)                                              PCLD 380
      J=J-1                                                             PCLD 390
      GO TO 2                                                           PCLD 400
    4 K=K+1                                                             PCLD 410
      IF (IDIMX-K) 5,5,1                                                PCLD 420
    5 RETURN                                                            PCLD 430
      END                                                               PCLD 440
