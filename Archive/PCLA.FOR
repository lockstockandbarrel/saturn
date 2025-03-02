C                                                                       PCLA  10
C     ..................................................................PCLA  20
C                                                                       PCLA  30
C        SUBROUTINE PCLA                                                PCLA  40
C                                                                       PCLA  50
C        PURPOSE                                                        PCLA  60
C           MOVE POLYNOMIAL X TO Y                                      PCLA  70
C                                                                       PCLA  80
C        USAGE                                                          PCLA  90
C           CALL PCLA(Y,IDIMY,X,IDIMX)                                  PCLA 100
C                                                                       PCLA 110
C        DESCRIPTION OF PARAMETERS                                      PCLA 120
C           Y     - VECTOR OF RESULTANT COEFFICIENTS, ORDERED FROM      PCLA 130
C                   SMALLEST TO LARGEST POWER                           PCLA 140
C           IDIMY - DIMENSION OF Y                                      PCLA 150
C           X     - VECTOR OF COEFFICIENTS FOR POLYNOMIAL, ORDERED      PCLA 160
C                   FROM SMALLEST TO LARGEST POWER                      PCLA 170
C           IDIMX - DIMENSION OF X                                      PCLA 180
C                                                                       PCLA 190
C        REMARKS                                                        PCLA 200
C           NONE                                                        PCLA 210
C                                                                       PCLA 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PCLA 230
C           NONE                                                        PCLA 240
C                                                                       PCLA 250
C        METHOD                                                         PCLA 260
C           IDIMY IS REPLACED BY IDIMX AND VECTOR X IS MOVED TO Y       PCLA 270
C                                                                       PCLA 280
C     ..................................................................PCLA 290
C                                                                       PCLA 300
      SUBROUTINE PCLA (Y,IDIMY,X,IDIMX)                                 PCLA 310
      DIMENSION X(1),Y(1)                                               PCLA 320
C                                                                       PCLA 330
      IDIMY=IDIMX                                                       PCLA 340
      IF(IDIMX) 30,30,10                                                PCLA 350
   10 DO 20 I=1,IDIMX                                                   PCLA 360
   20 Y(I)=X(I)                                                         PCLA 370
   30 RETURN                                                            PCLA 380
      END                                                               PCLA 390
