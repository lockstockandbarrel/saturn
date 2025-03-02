C                                                                       RECP  10
C     ..................................................................RECP  20
C                                                                       RECP  30
C        FUNCTION RECP                                                  RECP  40
C                                                                       RECP  50
C        PURPOSE                                                        RECP  60
C           CALCULATE RECIPROCAL OF AN ELEMENT. THIS IS A FORTRAN       RECP  70
C           FUNCTION SUBPROGRAM WHICH MAY BE USED AS AN ARGUMENT BY     RECP  80
C           SUBROUTINE MFUN.                                            RECP  90
C                                                                       RECP 100
C        USAGE                                                          RECP 110
C           RECP(E)                                                     RECP 120
C                                                                       RECP 130
C        DESCRIPTION OF PARAMETERS                                      RECP 140
C           E - MATRIX ELEMENT                                          RECP 150
C                                                                       RECP 160
C        REMARKS                                                        RECP 170
C           RECIPROCAL OF ZERO IS TAKEN TO BE 1.0E75                    RECP 180
C                                                                       RECP 190
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RECP 200
C           NONE                                                        RECP 210
C                                                                       RECP 220
C        METHOD                                                         RECP 230
C           RECIPROCAL OF ELEMENT E IS PLACED IN RECP                   RECP 240
C                                                                       RECP 250
C     ..................................................................RECP 260
C                                                                       RECP 270
      FUNCTION RECP(E)                                                  RECP 280
C                                                                       RECP 290
      BIG=1.0E75                                                        RECP 300
C                                                                       RECP 310
C        TEST ELEMENT FOR ZERO                                          RECP 320
C                                                                       RECP 330
      IF(E) 1,2,1                                                       RECP 340
C                                                                       RECP 350
C        IF NON-ZERO, CALCULATE RECIPROCAL                              RECP 360
C                                                                       RECP 370
    1 RECP=1.0/E                                                        RECP 380
      RETURN                                                            RECP 390
C                                                                       RECP 400
C        IF ZERO, SET EQUAL TO INFINITY                                 RECP 410
C                                                                       RECP 420
    2 RECP=SIGN(BIG,E)                                                  RECP 430
      RETURN                                                            RECP 440
      END                                                               RECP 450
