C                                                                       PADD  10
C     ..................................................................PADD  20
C                                                                       PADD  30
C        SUBROUTINE PADD                                                PADD  40
C                                                                       PADD  50
C        PURPOSE                                                        PADD  60
C           ADD TWO POLYNOMIALS                                         PADD  70
C                                                                       PADD  80
C        USAGE                                                          PADD  90
C           CALL PADD(Z,IDIMZ,X,IDIMX,Y,IDIMY)                          PADD 100
C                                                                       PADD 110
C        DESCRIPTION OF PARAMETERS                                      PADD 120
C           Z     - VECTOR OF RESULTANT COEFFICIENTS, ORDERED FROM      PADD 130
C                   SMALLEST TO LARGEST POWER                           PADD 140
C           IDIMZ - DIMENSION OF Z (CALCULATED)                         PADD 150
C           X     - VECTOR OF COEFFICIENTS FOR FIRST POLYNOMIAL, ORDEREDPADD 160
C                   FROM SMALLEST TO LARGEST POWER                      PADD 170
C           IDIMX - DIMENSION OF X (DEGREE IS IDIMX-1)                  PADD 180
C           Y     - VECTOR OF COEFFICIENTS FOR SECOND POLYNOMIAL,       PADD 190
C                   ORDERED FROM SMALLEST TO LARGEST POWER              PADD 200
C           IDIMY - DIMENSION OF Y (DEGREE IS IDIMY-1)                  PADD 210
C                                                                       PADD 220
C        REMARKS                                                        PADD 230
C           VECTOR Z MAY BE IN SAME LOCATION AS EITHER VECTOR X OR      PADD 240
C           VECTOR Y ONLY IF THE DIMENSION OF THAT VECTOR IS NOT LESS   PADD 250
C           THAN THE OTHER INPUT VECTOR                                 PADD 260
C           THE RESULTANT POLYNOMIAL MAY HAVE TRAILING ZERO COEFFICIENTSPADD 270
C                                                                       PADD 280
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PADD 290
C           NONE                                                        PADD 300
C                                                                       PADD 310
C        METHOD                                                         PADD 320
C           DIMENSION OF RESULTANT VECTOR IDIMZ IS CALCULATED AS THE    PADD 330
C           LARGER OF THE TWO INPUT VECTOR DIMENSIONS. CORRESPONDING    PADD 340
C           COEFFICIENTS ARE THEN ADDED TO FORM Z.                      PADD 350
C                                                                       PADD 360
C     ..................................................................PADD 370
C                                                                       PADD 380
      SUBROUTINE PADD(Z,IDIMZ,X,IDIMX,Y,IDIMY)                          PADD 390
      DIMENSION Z(1),X(1),Y(1)                                          PADD 400
C                                                                       PADD 410
C     TEST DIMENSIONS OF SUMMANDS                                       PADD 420
C                                                                       PADD 430
      NDIM=IDIMX                                                        PADD 440
      IF (IDIMX-IDIMY) 10,20,20                                         PADD 450
   10 NDIM=IDIMY                                                        PADD 460
   20 IF(NDIM) 90,90,30                                                 PADD 470
   30 DO 80 I=1,NDIM                                                    PADD 480
      IF(I-IDIMX) 40,40,60                                              PADD 490
   40 IF(I-IDIMY) 50,50,70                                              PADD 500
   50 Z(I)=X(I)+Y(I)                                                    PADD 510
      GO TO 80                                                          PADD 520
   60 Z(I)=Y(I)                                                         PADD 530
      GO TO 80                                                          PADD 540
   70 Z(I)=X(I)                                                         PADD 550
   80 CONTINUE                                                          PADD 560
   90 IDIMZ=NDIM                                                        PADD 570
      RETURN                                                            PADD 580
      END                                                               PADD 590
