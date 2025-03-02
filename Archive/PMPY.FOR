C                                                                       PMPY  10
C     ..................................................................PMPY  20
C                                                                       PMPY  30
C        SUBROUTINE PMPY                                                PMPY  40
C                                                                       PMPY  50
C        PURPOSE                                                        PMPY  60
C           MULTIPLY TWO POLYNOMIALS                                    PMPY  70
C                                                                       PMPY  80
C        USAGE                                                          PMPY  90
C           CALL PMPY(Z,IDIMZ,X,IDIMX,Y,IDIMY)                          PMPY 100
C                                                                       PMPY 110
C        DESCRIPTION OF PARAMETERS                                      PMPY 120
C           Z     - VECTOR OF RESULTANT COEFFICIENTS, ORDERED FROM      PMPY 130
C                   SMALLEST TO LARGEST POWER                           PMPY 140
C           IDIMZ - DIMENSION OF Z (CALCULATED)                         PMPY 150
C           X     - VECTOR OF COEFFICIENTS FOR FIRST POLYNOMIAL, ORDEREDPMPY 160
C                   FROM SMALLEST TO LARGEST POWER                      PMPY 170
C           IDIMX - DIMENSION OF X (DEGREE IS IDIMX-1)                  PMPY 180
C           Y     - VECTOR OF COEFFICIENTS FOR SECOND POLYNOMIAL,       PMPY 190
C                   ORDERED FROM SMALLEST TO LARGEST POWER              PMPY 200
C           IDIMY - DIMENSION OF Y (DEGREE IS IDIMY-1)                  PMPY 210
C                                                                       PMPY 220
C        REMARKS                                                        PMPY 230
C           Z CANNOT BE IN THE SAME LOCATION AS X                       PMPY 240
C           Z CANNOT BE IN THE SAME LOCATION AS Y                       PMPY 250
C                                                                       PMPY 260
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PMPY 270
C           NONE                                                        PMPY 280
C                                                                       PMPY 290
C        METHOD                                                         PMPY 300
C           DIMENSION OF Z IS CALCULATED AS IDIMX+IDIMY-1               PMPY 310
C           THE COEFFICIENTS OF Z ARE CALCULATED AS SUM OF PRODUCTS     PMPY 320
C           OF COEFFICIENTS OF X AND Y , WHOSE EXPONENTS ADD UP TO THE  PMPY 330
C           CORRESPONDING EXPONENT OF Z.                                PMPY 340
C                                                                       PMPY 350
C     ..................................................................PMPY 360
C                                                                       PMPY 370
      SUBROUTINE PMPY(Z,IDIMZ,X,IDIMX,Y,IDIMY)                          PMPY 380
      DIMENSION Z(1),X(1),Y(1)                                          PMPY 390
C                                                                       PMPY 400
      IF(IDIMX*IDIMY)10,10,20                                           PMPY 410
   10 IDIMZ=0                                                           PMPY 420
      GO TO 50                                                          PMPY 430
   20 IDIMZ=IDIMX+IDIMY-1                                               PMPY 440
      DO 30 I=1,IDIMZ                                                   PMPY 450
   30 Z(I)=0.                                                           PMPY 460
      DO 40 I=1,IDIMX                                                   PMPY 470
      DO 40 J=1,IDIMY                                                   PMPY 480
      K=I+J-1                                                           PMPY 490
   40 Z(K)=X(I)*Y(J)+Z(K)                                               PMPY 500
   50 RETURN                                                            PMPY 510
      END                                                               PMPY 520
