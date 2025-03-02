C                                                                       ABSN  10
C     ..................................................................ABSN  20
C                                                                       ABSN  30
C        SUBROUTINE ABSNT                                               ABSN  40
C                                                                       ABSN  50
C        PURPOSE                                                        ABSN  60
C           TEST MISSING OR ZERO VALUES FOR EACH OBSERVATION IN         ABSN  70
C           MATRIX A.                                                   ABSN  80
C                                                                       ABSN  90
C        USAGE                                                          ABSN 100
C           CALL ABSNT (A,S,NO,NV)                                      ABSN 110
C                                                                       ABSN 120
C        DESCRIPTION OF PARAMETERS                                      ABSN 130
C           A  - OBSERVATION MATRIX, NO BY NV                           ABSN 140
C           S  - OUTPUT VECTOR OF LENGTH NO INDICATING THE FOLLOWING    ABSN 150
C                CODES FOR EACH OBSERVATION.                            ABSN 160
C                1  THERE IS NOT A MISSING OR ZERO VALUE.               ABSN 170
C                0  AT LEAST ONE VALUE IS MISSING OR ZERO.              ABSN 180
C           NO - NUMBER OF OBSERVATIONS. NO MUST BE > OR = TO 1.        ABSN 190
C           NV - NUMBER OF VARIABLES FOR EACH OBSERVATION. NV MUST BE   ABSN 200
C                GREATER THAN OR EQUAL TO 1.                            ABSN 201
C                                                                       ABSN 210
C        REMARKS                                                        ABSN 220
C           NONE                                                        ABSN 230
C                                                                       ABSN 240
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  ABSN 250
C           NONE                                                        ABSN 260
C                                                                       ABSN 270
C        METHOD                                                         ABSN 280
C           A TEST IS MADE FOR EACH ROW (OBSERVATION) OF THE MATRIX A.  ABSN 290
C           IF THERE IS NOT A MISSING OR ZERO VALUE, 1 IS PLACED IN     ABSN 300
C           S(J). IF AT LEAST ONE VALUE IS MISSING OR ZERO, 0 IS PLACED ABSN 310
C           IN S(J).                                                    ABSN 320
C                                                                       ABSN 330
C     ..................................................................ABSN 340
C                                                                       ABSN 350
      SUBROUTINE ABSNT(A,S,NO,NV)                                       ABSN 360
      DIMENSION A(1),S(1)                                               ABSN 370
C                                                                       ABSN 380
      DO 20 J=1,NO                                                      ABSN 390
      IJ=J-NO                                                           ABSN 400
      S(J)=1.0                                                          ABSN 410
      DO 10 I=1,NV                                                      ABSN 420
      IJ=IJ+NO                                                          ABSN 430
      IF(A(IJ)) 10,5,10                                                 ABSN 440
    5 S(J)=0                                                            ABSN 450
      GO TO 20                                                          ABSN 460
   10 CONTINUE                                                          ABSN 470
   20 CONTINUE                                                          ABSN 480
      RETURN                                                            ABSN 490
      END                                                               ABSN 500
