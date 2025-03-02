C                                                                       NTRI  10
C.......................................................................NTRI  20
C                                                                       NTRI  30
C        SUBROUTINE NDTRI                                               NTRI  40
C                                                                       NTRI  50
C        PURPOSE                                                        NTRI  60
C          COMPUTES X = P**(-1)(Y), THE ARGUMENT X SUCH THAT Y= P(X) =  NTRI  70
C          THE PROBABILITY THAT THE RANDOM VARIABLE U, DISTRIBUTED      NTRI  80
C          NORMALLY(0,1), IS LESS THAN OR EQUAL TO X.  F(X), THE        NTRI  90
C          ORDINATE OF THE NORMAL DENSITY, AT X, IS ALSO COMPUTED.      NTRI 100
C                                                                       NTRI 110
C        USAGE                                                          NTRI 120
C          CALL NDTRI(P,X,D,IER)                                        NTRI 130
C                                                                       NTRI 140
C        DESCRIPTION OF PARAMETERS                                      NTRI 150
C          P   - INPUT PROBABILITY.                                     NTRI 160
C          X   - OUTPUT ARGUMENT SUCH THAT P = Y = THE PROBABILITY THAT NTRI 170
C                   U, THE RANDOM VARIABLE, IS LESS THAN OR EQUAL TO X. NTRI 180
C          D   - OUTPUT DENSITY, F(X).                                  NTRI 190
C          IER - OUTPUT ERROR CODE                                      NTRI 200
C                = -1 IF P IS NOT IN THE INTERVAL (0,1), INCLUSIVE.     NTRI 210
C                  X=D=.99999E+74 IN THIS CASE                          NTRI 211
C                = 0 IF THERE IS NO ERROR.  SEE REMARKS, BELOW.         NTRI 220
C                                                                       NTRI 230
C        REMARKS                                                        NTRI 240
C          MAXIMUM ERROR IS 0.00045.                                    NTRI 250
C          IF P = 0, X IS SET TO -(10)**74.  D IS SET TO 0.             NTRI 260
C          IF P = 1, X IS SET TO  (10)**74.  D IS SET TO 0.             NTRI 270
C                                                                       NTRI 280
C        SUBROUTINES AND SUBPROGRAMS REQUIRED                           NTRI 290
C          NONE                                                         NTRI 300
C                                                                       NTRI 310
C        METHOD                                                         NTRI 320
C          BASED ON APPROXIMATIONS IN C. HASTINGS, APPROXIMATIONS FOR   NTRI 330
C          DIGITAL COMPUTERS, PRINCETON UNIV. PRESS, PRINCETON, N.J.,   NTRI 340
C          1955.  SEE EQUATION 26.2.23, HANDBOOK OF MATHEMATICAL        NTRI 350
C          FUNCTIONS, ABRAMOWITZ AND STEGUN, DOVER PUBLICATIONS, INC.,  NTRI 360
C          NEW YORK.                                                    NTRI 370
C                                                                       NTRI 380
C.......................................................................NTRI 390
C                                                                       NTRI 400
      SUBROUTINE NDTRI(P,X,D,IE)                                        NTRI 410
C                                                                       NTRI 420
      IE=0                                                              NTRI 430
      X=.99999E+74                                                      NTRI 431
      D=X                                                               NTRI 432
      IF(P)1,4,2                                                        NTRI 440
    1 IE=-1                                                             NTRI 450
      GO TO 12                                                          NTRI 460
    2 IF (P-1.0)7,5,1                                                   NTRI 470
    4 X=-.999999E+74                                                    NTRI 480
    5 D=0.0                                                             NTRI 490
      GO TO 12                                                          NTRI 500
C                                                                       NTRI 510
C                                                                       NTRI 520
    7 D=P                                                               NTRI 530
      IF(D-0.5)9,9,8                                                    NTRI 540
    8 D=1.0-D                                                           NTRI 550
    9 T2=ALOG(1.0/(D*D))                                                NTRI 560
      T=SQRT(T2)                                                        NTRI 570
      X=T-(2.515517+0.802853*T+0.010328*T2)/(1.0+1.432788*T+0.189269*T2 NTRI 580
     1  +0.001308*T*T2)                                                 NTRI 590
      IF(P-0.5)10,10,11                                                 NTRI 600
   10 X=-X                                                              NTRI 610
   11 D=0.3989423*EXP(-X*X/2.0)                                         NTRI 620
   12 RETURN                                                            NTRI 630
      END                                                               NTRI 640
