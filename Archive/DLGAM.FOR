C                                                                       DLGA  10
C     ..................................................................DLGA  20
C                                                                       DLGA  30
C        SUBROUTINE DLGAM                                               DLGA  40
C                                                                       DLGA  50
C        PURPOSE                                                        DLGA  60
C           COMPUTES THE DOUBLE PRECISION NATURAL LOGARITHM OF THE      DLGA  70
C           GAMMA FUNCTION OF A GIVEN DOUBLE PRECISION ARGUMENT.        DLGA  80
C                                                                       DLGA  90
C        USAGE                                                          DLGA 100
C           CALL DLGAM(XX,DLNG,IER)                                     DLGA 110
C                                                                       DLGA 120
C        DESCRIPTION OF PARAMETERS                                      DLGA 130
C           XX   - THE DOUBLE PRECISION ARGUMENT FOR THE LOG GAMMA      DLGA 140
C                  FUNCTION.                                            DLGA 150
C           DLNG - THE RESULTANT DOUBLE PRECISION LOG GAMMA FUNCTION    DLGA 160
C                  VALUE.                                               DLGA 170
C           IER  - RESULTANT ERROR CODE WHERE                           DLGA 180
C                  IER= 0----NO ERROR.                                  DLGA 190
C                  IER=-1----XX IS WITHIN 10**(-9) OF BEING ZERO OR XX  DLGA 200
C                            IS NEGATIVE.  DLNG IS SET TO -1.OD75.      DLGA 210
C                  IER=+1----XX IS GREATER THAN 10**70. DLNG IS SET TO  DLGA 220
C                            +1.OD75.                                   DLGA 230
C                                                                       DLGA 240
C        REMARKS                                                        DLGA 250
C           NONE                                                        DLGA 260
C                                                                       DLGA 270
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DLGA 280
C           NONE                                                        DLGA 290
C                                                                       DLGA 300
C        METHOD                                                         DLGA 310
C           THE EULER-MCLAURIN EXPANSION TO THE SEVENTH DERIVATIVE TERM DLGA 320
C           IS USED, AS GIVEN BY M. ABRAMOWITZ AND I.A. STEGUN,         DLGA 330
C           'HANDBOOK OF MATHEMATICAL FUNCTIONS', U. S. DEPARTMENT OF   DLGA 340
C           COMMERCE, NATIONAL BUREAU OF STANDARDS APPLIED MATHEMATICS  DLGA 350
C           SERIES, 1966, EQUATION 6.1.41.                              DLGA 360
C                                                                       DLGA 370
C     ..................................................................DLGA 380
C                                                                       DLGA 390
      SUBROUTINE DLGAM(XX,DLNG,IER)                                     DLGA 400
      DOUBLE PRECISION XX,ZZ,TERM,RZ2,DLNG                              DLGA 410
      IER=0                                                             DLGA 420
      ZZ=XX                                                             DLGA 430
      IF(XX-1.D10) 2,2,1                                                DLGA 440
    1 IF(XX-1.D70) 8,9,9                                                DLGA 450
C                                                                       DLGA 460
C        SEE IF XX IS NEAR ZERO OR NEGATIVE                             DLGA 470
C                                                                       DLGA 480
    2 IF(XX-1.D-9) 3,3,4                                                DLGA 490
    3 IER=-1                                                            DLGA 500
      DLNG=-1.D75                                                       DLGA 510
      GO TO 10                                                          DLGA 520
C                                                                       DLGA 530
C        XX GREATER THAN ZERO AND LESS THAN OR EQUAL TO 1.D+10          DLGA 540
C                                                                       DLGA 550
    4 TERM=1.D0                                                         DLGA 560
    5 IF(ZZ-18.D0) 6,6,7                                                DLGA 570
    6 TERM=TERM*ZZ                                                      DLGA 580
      ZZ=ZZ+1.D0                                                        DLGA 590
      GO TO 5                                                           DLGA 600
    7 RZ2=1.D0/ZZ**2                                                    DLGA 610
      DLNG =(ZZ-0.5D0)*DLOG(ZZ)-ZZ +0.9189385332046727 -DLOG(TERM)+     DLGA 620
     1(1.D0/ZZ)*(.8333333333333333D-1 -(RZ2*(.2777777777777777D-2 +(RZ2*DLGA 630
     2(.7936507936507936D-3 -(RZ2*(.5952380952380952D-3)))))))          DLGA 640
      GO TO 10                                                          DLGA 650
C                                                                       DLGA 660
C        XX GREATER THAN 1.D+10 AND LESS THAN 1.D+70                    DLGA 670
C                                                                       DLGA 680
    8 DLNG=ZZ*(DLOG(ZZ)-1.D0)                                           DLGA 690
      GO TO 10                                                          DLGA 700
C                                                                       DLGA 710
C        XX GREATER THAN OR EQUAL TO 1.D+70                             DLGA 720
C                                                                       DLGA 730
    9 IER=+1                                                            DLGA 740
      DLNG=1.D75                                                        DLGA 750
   10 RETURN                                                            DLGA 760
      END                                                               DLGA 770
