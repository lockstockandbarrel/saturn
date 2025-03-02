C                                                                       PVSU  10
C     ..................................................................PVSU  20
C                                                                       PVSU  30
C        SUBROUTINE PVSUB                                               PVSU  40
C                                                                       PVSU  50
C        PURPOSE                                                        PVSU  60
C           SUBSTITUTE VARIABLE OF A POLYNOMIAL BY ANOTHER POLYNOMIAL   PVSU  70
C                                                                       PVSU  80
C        USAGE                                                          PVSU  90
C           CALL PVSUB(Z,IDIMZ,X,IDIMX,Y,IDIMY,WORK1,WORK2)             PVSU 100
C                                                                       PVSU 110
C        DESCRIPTION OF PARAMETERS                                      PVSU 120
C           Z     - VECTOR OF COEFFICIENTS FOR RESULTANT POLYNOMIAL,    PVSU 130
C                   ORDERED FROM SMALLEST TO LARGEST POWER              PVSU 140
C           IDIMZ - DIMENSION OF Z                                      PVSU 150
C           X     - VECTOR OF COEFFICIENTS FOR ORIGINAL POLYNOMIAL,     PVSU 160
C                   ORDERED FROM SMALLEST TO LARGEST POWER              PVSU 170
C           IDIMX - DIMENSION OF X                                      PVSU 180
C           Y     - VECTOR OF COEFFICIENTS FOR POLYNOMIAL WHICH IS      PVSU 190
C                   SUBSTITUTED FOR VARIABLE, ORDERED FROM SMALLEST TO  PVSU 200
C                   LARGEST POWER                                       PVSU 210
C           IDIMY - DIMENSION OF Y                                      PVSU 220
C           WORK1 - WORKING STORAGE ARRAY (SAME DIMENSION AS Z)         PVSU 230
C           WORK2 - WORKING STORAGE ARRAY (SAME DIMENSION AS Z)         PVSU 240
C                                                                       PVSU 250
C        REMARKS                                                        PVSU 260
C           NONE                                                        PVSU 270
C                                                                       PVSU 280
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PVSU 290
C           PMPY                                                        PVSU 300
C           PADDM                                                       PVSU 310
C           PCLA                                                        PVSU 320
C                                                                       PVSU 330
C        METHOD                                                         PVSU 340
C           VARIABLE OF POLYNOMIAL X IS SUBSTITUTED BY POLYNOMIAL Y     PVSU 350
C           TO FORM POLYNOMIAL Z. DIMENSION OF NEW POLYNOMIAL IS        PVSU 360
C           (IDIMX-1)*(IDIMY-1)+1. SUBROUTINE REQUIRES TWO WORK AREAS   PVSU 370
C                                                                       PVSU 380
C     ..................................................................PVSU 390
C                                                                       PVSU 400
      SUBROUTINE PVSUB(Z,IDIMZ,X,IDIMX,Y,IDIMY,WORK1,WORK2)             PVSU 410
      DIMENSION Z(1),X(1),Y(1),WORK1(1),WORK2(1)                        PVSU 420
C                                                                       PVSU 430
C     TEST OF DIMENSIONS                                                PVSU 440
C                                                                       PVSU 450
      IF (IDIMX-1) 1,3,3                                                PVSU 460
    1 IDIMZ=0                                                           PVSU 470
    2 RETURN                                                            PVSU 480
C                                                                       PVSU 490
    3 IDIMZ=1                                                           PVSU 500
      Z(1)=X(1)                                                         PVSU 510
      IF (IDIMY*IDIMX-IDIMY) 2,2,4                                      PVSU 520
    4 IW1=1                                                             PVSU 530
      WORK1(1)=1.                                                       PVSU 540
C                                                                       PVSU 550
      DO 5 I=2,IDIMX                                                    PVSU 560
      CALL PMPY(WORK2,IW2,Y,IDIMY,WORK1,IW1)                            PVSU 570
      CALL PCLA(WORK1,IW1,WORK2,IW2)                                    PVSU 580
      FACT=X(I)                                                         PVSU 590
      CALL PADDM(Z,IDIMR,Z,IDIMZ,FACT,WORK1,IW1)                        PVSU 600
      IDIMZ=IDIMR                                                       PVSU 605
    5 CONTINUE                                                          PVSU 610
      GO TO 2                                                           PVSU 620
      END                                                               PVSU 630
