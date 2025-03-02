C                                                                       PDER  10
C     ..................................................................PDER  20
C                                                                       PDER  30
C        SUBROUTINE PDER                                                PDER  40
C                                                                       PDER  50
C        PURPOSE                                                        PDER  60
C           FIND DERIVATIVE OF A POLYNOMIAL                             PDER  70
C                                                                       PDER  80
C        USAGE                                                          PDER  90
C           CALL PDER(Y,IDIMY,X,IDIMX)                                  PDER 100
C                                                                       PDER 110
C        DESCRIPTION OF PARAMETERS                                      PDER 120
C           Y     - VECTOR OF COEFFICIENTS FOR DERIVATIVE, ORDERED FROM PDER 130
C                   SMALLEST TO LARGEST POWER                           PDER 140
C           IDIMY - DIMENSION OF Y (EQUAL TO IDIMX-1)                   PDER 150
C           X     - VECTOR OF COEFFICIENTS FOR ORIGINAL POLYNOMIAL,     PDER 160
C                   ORDERED FROM SMALLEST TO LARGEST POWER              PDER 170
C           IDIMX - DIMENSION OF X                                      PDER 180
C                                                                       PDER 190
C        REMARKS                                                        PDER 200
C           NONE                                                        PDER 210
C                                                                       PDER 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PDER 230
C           NONE                                                        PDER 240
C                                                                       PDER 250
C        METHOD                                                         PDER 260
C           DIMENSION OF Y IS SET AT DIMENSION OF X LESS ONE. DERIVATIVEPDER 270
C           IS THEN CALCULATED BY MULTIPLYING COEFFICIENTS BY THEIR     PDER 280
C           RESPECTIVE EXPONENTS.                                       PDER 290
C                                                                       PDER 300
C     ..................................................................PDER 310
C                                                                       PDER 320
      SUBROUTINE PDER(Y,IDIMY,X,IDIMX)                                  PDER 330
      DIMENSION X(1),Y(1)                                               PDER 340
C                                                                       PDER 350
C     TEST OF DIMENSION                                                 PDER 360
      IF (IDIMX-1) 3,3,1                                                PDER 370
    1 IDIMY=IDIMX-1                                                     PDER 380
      EXPT=0.                                                           PDER 390
      DO 2 I=1,IDIMY                                                    PDER 400
      EXPT=EXPT+1.                                                      PDER 410
    2 Y(I)=X(I+1)*EXPT                                                  PDER 420
      GO TO 4                                                           PDER 430
    3 IDIMY=0                                                           PDER 440
    4 RETURN                                                            PDER 450
      END                                                               PDER 460
