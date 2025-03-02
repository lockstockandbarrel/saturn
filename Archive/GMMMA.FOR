C                                                                       GMMM  10
C     ..................................................................GMMM  20
C                                                                       GMMM  30
C        SUBROUTINE GMMMA                                               GMMM  40
C                                                                       GMMM  50
C        PURPOSE                                                        GMMM  60
C           COMPUTES THE GAMMA FUNCTION FOR A GIVEN ARGUMENT            GMMM  70
C                                                                       GMMM  80
C        USAGE                                                          GMMM  90
C           CALL GMMMA(XX,GX,IER)                                       GMMM 100
C                                                                       GMMM 110
C        DESCRIPTION OF PARAMETERS                                      GMMM 120
C           XX -THE ARGUMENT FOR THE GAMMA FUNCTION                     GMMM 130
C           GX -THE RESULTANT GAMMA FUNCTION VALUE                      GMMM 140
C           IER-RESULTANT ERROR CODE WHERE                              GMMM 150
C               IER=0  NO ERROR                                         GMMM 160
C               IER=1  XX IS WITHIN .000001 OF BEING A NEGATIVE INTEGER GMMM 170
C               IER=2  XX GT 57, OVERFLOW, GX SET TO 1.0E75             GMMM 180
C                                                                       GMMM 190
C        REMARKS                                                        GMMM 200
C           NONE                                                        GMMM 210
C                                                                       GMMM 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  GMMM 230
C           NONE                                                        GMMM 240
C                                                                       GMMM 250
C        METHOD                                                         GMMM 260
C           THE RECURSION RELATION AND POLYNOMIAL APPROXIMATION         GMMM 270
C           BY C.HASTINGS,JR., 'APPROXIMATIONS FOR DIGITAL COMPUTERS',  GMMM 280
C           PRINCETON UNIVERSITY PRESS, 1955                            GMMM 290
C                                                                       GMMM 300
C     ..................................................................GMMM 310
C                                                                       GMMM 320
      SUBROUTINE GMMMA(XX,GX,IER)                                       GMMM 330
      IF(XX-57.)6,6,4                                                   GMMM 340
    4 IER=2                                                             GMMM 350
      GX=1.E75                                                          GMMM 360
      RETURN                                                            GMMM 370
    6 X=XX                                                              GMMM 380
      ERR=1.0E-6                                                        GMMM 390
      IER=0                                                             GMMM 400
      GX=1.0                                                            GMMM 410
      IF(X-2.0)50,50,15                                                 GMMM 420
   10 IF(X-2.0)110,110,15                                               GMMM 430
   15 X=X-1.0                                                           GMMM 440
      GX=GX*X                                                           GMMM 450
      GO TO 10                                                          GMMM 460
   50 IF(X-1.0)60,120,110                                               GMMM 470
C                                                                       GMMM 480
C        SEE IF X IS NEAR NEGATIVE INTEGER OR ZERO                      GMMM 490
C                                                                       GMMM 500
   60 IF(X-ERR)62,62,80                                                 GMMM 510
   62 Y=FLOAT(INT(X))-X                                                 GMMM 520
      IF(ABS(Y)-ERR)130,130,64                                          GMMM 530
   64 IF(1.0-Y-ERR)130,130,70                                           GMMM 540
C                                                                       GMMM 550
C        X NOT NEAR A NEGATIVE INTEGER OR ZERO                          GMMM 560
C                                                                       GMMM 570
   70 IF(X-1.0)80,80,110                                                GMMM 580
   80 GX=GX/X                                                           GMMM 590
      X=X+1.0                                                           GMMM 600
      GO TO 70                                                          GMMM 610
  110 Y=X-1.0                                                           GMMM 620
      GY=1.0+Y*(-0.5771017+Y*(+0.9858540+Y*(-0.8764218+Y*(+0.8328212+   GMMM 630
     1Y*(-0.5684729+Y*(+0.2548205+Y*(-0.05149930)))))))                 GMMM 640
      GX=GX*GY                                                          GMMM 650
  120 RETURN                                                            GMMM 660
  130 IER=1                                                             GMMM 670
      RETURN                                                            GMMM 680
      END                                                               GMMM 690
