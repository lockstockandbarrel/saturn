C                                                                       PSUB  10
C     ..................................................................PSUB  20
C                                                                       PSUB  30
C        SUBROUTINE PSUB                                                PSUB  40
C                                                                       PSUB  50
C        PURPOSE                                                        PSUB  60
C           SUBTRACT ONE POLYNOMIAL FROM ANOTHER                        PSUB  70
C                                                                       PSUB  80
C        USAGE                                                          PSUB  90
C           CALL PSUB(Z,IDIMZ,X,IDIMX,Y,IDIMY)                          PSUB 100
C                                                                       PSUB 110
C        DESCRIPTION OF PARAMETERS                                      PSUB 120
C           Z     - VECTOR OF RESULTANT COEFFICIENTS, ORDERED FROM      PSUB 130
C                   SMALLEST TO LARGEST POWER                           PSUB 140
C           IDIMZ - DIMENSION OF Z (CALCULATED)                         PSUB 150
C           X     - VECTOR OF COEFFICIENTS FOR FIRST POLYNOMIAL, ORDEREDPSUB 160
C                   FROM SMALLEST TO LARGEST POWER                      PSUB 170
C           IDIMX - DIMENSION OF X (DEGREE IS IDIMX-1)                  PSUB 180
C           Y     - VECTOR OF COEFFICIENTS FOR SECOND POLYNOMIAL,       PSUB 190
C                   ORDERED FROM SMALLEST TO LARGEST POWER              PSUB 200
C           IDIMY - DIMENSION OF Y (DEGREE IS IDIMY-1)                  PSUB 210
C                                                                       PSUB 220
C        REMARKS                                                        PSUB 230
C           VECTOR Z MAY BE IN SAME LOCATION AS EITHER VECTOR X OR      PSUB 240
C           VECTOR Y ONLY IF THE DIMENSION OF THAT VECTOR IS NOT LESS   PSUB 250
C           THAN THE OTHER INPUT VECTOR                                 PSUB 260
C           THE RESULTANT POLYNOMIAL MAY HAVE TRAILING ZERO COEFFICIENTSPSUB 270
C                                                                       PSUB 280
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PSUB 290
C           NONE                                                        PSUB 300
C                                                                       PSUB 310
C        METHOD                                                         PSUB 320
C           DIMENSION OF RESULTANT VECTOR IDIMZ IS CALCULATED AS THE    PSUB 330
C           LARGER OF THE TWO INPUT VECTOR DIMENSIONS. COEFFICIENTS IN  PSUB 340
C           VECTOR Y ARE THEN SUBTRACTED FROM CORRESPONDING COEFFICIENTSPSUB 350
C           IN VECTOR X.                                                PSUB 360
C                                                                       PSUB 370
C     ..................................................................PSUB 380
C                                                                       PSUB 390
      SUBROUTINE PSUB(Z,IDIMZ,X,IDIMX,Y,IDIMY)                          PSUB 400
      DIMENSION Z(1),X(1),Y(1)                                          PSUB 410
C                                                                       PSUB 420
C     TEST DIMENSIONS OF SUMMANDS                                       PSUB 430
C                                                                       PSUB 440
      NDIM=IDIMX                                                        PSUB 450
      IF (IDIMX-IDIMY) 10,20,20                                         PSUB 460
   10 NDIM=IDIMY                                                        PSUB 470
   20 IF (NDIM) 90,90,30                                                PSUB 480
   30 DO 80 I=1,NDIM                                                    PSUB 490
      IF (I-IDIMX) 40,40,60                                             PSUB 500
   40 IF (I-IDIMY) 50,50,70                                             PSUB 510
   50 Z(I)=X(I)-Y(I)                                                    PSUB 520
      GO TO 80                                                          PSUB 530
   60 Z(I)=-Y(I)                                                        PSUB 540
      GO TO 80                                                          PSUB 550
   70 Z(I)=X(I)                                                         PSUB 560
   80 CONTINUE                                                          PSUB 570
   90 IDIMZ=NDIM                                                        PSUB 580
      RETURN                                                            PSUB 590
      END                                                               PSUB 600
