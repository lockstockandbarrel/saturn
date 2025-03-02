C                                                                       PDDM  10
C     ..................................................................PDDM  20
C                                                                       PDDM  30
C        SUBROUTINE PADDM                                               PDDM  40
C                                                                       PDDM  50
C        PURPOSE                                                        PDDM  60
C           ADD COEFFICIENTS OF ONE POLYNOMIAL TO THE PRODUCT OF A      PDDM  70
C           FACTOR BY COEFFICIENTS OF ANOTHER POLYNOMIAL                PDDM  80
C                                                                       PDDM  90
C        USAGE                                                          PDDM 100
C           CALL PADDM(Z,IDIMZ,X,IDIMX,FACT,Y,IDIMY)                    PDDM 110
C                                                                       PDDM 120
C        DESCRIPTION OF PARAMETERS                                      PDDM 130
C           Z     - VECTOR OF RESULTANT COEFFICIENTS, ORDERED FROM      PDDM 140
C                   SMALLEST TO LARGEST POWER                           PDDM 150
C           IDIMZ - DIMENSION OF Z (CALCULATED)                         PDDM 160
C           X     - VECTOR OF COEFFICIENTS FOR FIRST POLYNOMIAL, ORDEREDPDDM 170
C                   FROM SMALLEST TO LARGEST POWER                      PDDM 180
C           IDIMX - DIMENSION OF X (DEGREE IS IDIMX-1)                  PDDM 190
C           FACT  - FACTOR TO BE MULTIPLIED BY VECTOR Y                 PDDM 200
C           Y     - VECTOR OF COEFFICIENTS FOR SECOND POLYNOMIAL,       PDDM 210
C                   ORDERED FROM SMALLEST TO LARGEST POWER              PDDM 220
C           IDIMY - DIMENSION OF Y (DEGREE IS IDIMY-1)                  PDDM 230
C                                                                       PDDM 240
C        REMARKS                                                        PDDM 250
C           VECTOR Z MAY BE IN SAME LOCATION AS EITHER VECTOR X OR      PDDM 260
C           VECTOR Y ONLY IF THE DIMENSION OF THAT VECTOR IS NOT LESS   PDDM 270
C           THAN THE OTHER INPUT VECTOR                                 PDDM 280
C           THE RESULTANT POLYNOMIAL MAY HAVE TRAILING ZERO COEFFICIENTSPDDM 290
C                                                                       PDDM 300
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PDDM 310
C           NONE                                                        PDDM 320
C                                                                       PDDM 330
C        METHOD                                                         PDDM 340
C           DIMENSION OF RESULTANT VECTOR IDIMZ IS CALCULATED AS THE    PDDM 350
C           LARGER OF THE TWO INPUT VECTOR DIMENSIONS. COEFFICIENT IN   PDDM 360
C           VECTOR X IS THEN ADDED TO COEFFICIENT IN VECTOR Y MULTIPLIEDPDDM 370
C           BY FACTOR TO FORM Z.                                        PDDM 380
C                                                                       PDDM 390
C     ..................................................................PDDM 400
C                                                                       PDDM 410
      SUBROUTINE PADDM(Z,IDIMZ,X,IDIMX,FACT,Y,IDIMY)                    PDDM 420
      DIMENSION Z(1),X(1),Y(1)                                          PDDM 430
C                                                                       PDDM 440
C     TEST DIMENSIONS OF SUMMANDS                                       PDDM 450
C                                                                       PDDM 460
      NDIM=IDIMX                                                        PDDM 470
      IF(IDIMX-IDIMY) 10,20,20                                          PDDM 480
   10 NDIM=IDIMY                                                        PDDM 490
   20 IF(NDIM) 90,90,30                                                 PDDM 500
   30 DO 80 I=1,NDIM                                                    PDDM 510
      IF(I-IDIMX) 40,40,60                                              PDDM 520
   40 IF(I-IDIMY) 50,50,70                                              PDDM 530
   50 Z(I)=FACT*Y(I)+X(I)                                               PDDM 540
      GO TO 80                                                          PDDM 550
   60 Z(I)=FACT*Y(I)                                                    PDDM 560
      GO TO 80                                                          PDDM 570
   70 Z(I)=X(I)                                                         PDDM 580
   80 CONTINUE                                                          PDDM 590
   90 IDIMZ=NDIM                                                        PDDM 600
      RETURN                                                            PDDM 610
      END                                                               PDDM 620
