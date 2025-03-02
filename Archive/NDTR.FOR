C                                                                       NDTR  10
C.......................................................................NDTR  20
C                                                                       NDTR  30
C        SUBROUTINE NDTR                                                NDTR  40
C                                                                       NDTR  50
C        PURPOSE                                                        NDTR  60
C           COMPUTES Y = P(X) = PROBABILITY THAT THE RANDOM VARIABLE  U,NDTR  70
C           DISTRIBUTED NORMALLY(0,1), IS LESS THAN OR EQUAL TO X.      NDTR  80
C           F(X), THE ORDINATE OF THE NORMAL DENSITY AT X, IS ALSO      NDTR  90
C           COMPUTED.                                                   NDTR 100
C                                                                       NDTR 110
C        USAGE                                                          NDTR 120
C           CALL NDTR(X,P,D)                                            NDTR 130
C                                                                       NDTR 140
C        DESCRIPTION OF PARAMETERS                                      NDTR 150
C           X--INPUT SCALAR FOR WHICH P(X) IS COMPUTED.                 NDTR 160
C           P--OUTPUT PROBABILITY.                                      NDTR 170
C           D--OUTPUT DENSITY.                                          NDTR 180
C                                                                       NDTR 190
C        REMARKS                                                        NDTR 200
C           MAXIMUM ERROR IS 0.0000007.                                 NDTR 210
C                                                                       NDTR 220
C        SUBROUTINES AND SUBPROGRAMS REQUIRED                           NDTR 230
C           NONE                                                        NDTR 240
C                                                                       NDTR 250
C        METHOD                                                         NDTR 260
C           BASED ON APPROXIMATIONS IN C. HASTINGS, APPROXIMATIONS FOR  NDTR 270
C           DIGITAL COMPUTERS, PRINCETON UNIV. PRESS, PRINCETON, N.J.,  NDTR 280
C           1955.  SEE EQUATION 26.2.17, HANDBOOK OF MATHEMATICAL       NDTR 290
C           FUNCTIONS, ABRAMOWITZ AND STEGUN, DOVER PUBLICATIONS, INC., NDTR 300
C           NEW YORK.                                                   NDTR 310
C                                                                       NDTR 320
C.......................................................................NDTR 330
C                                                                       NDTR 340
      SUBROUTINE NDTR(X,P,D)                                            NDTR 350
C                                                                       NDTR 360
      AX=ABS(X)                                                         NDTR 370
      T=1.0/(1.0+.2316419*AX)                                           NDTR 380
      D=0.3989423*EXP(-X*X/2.0)                                         NDTR 390
      P = 1.0 - D*T*((((1.330274*T - 1.821256)*T + 1.781478)*T -        NDTR 400
     1  0.3565638)*T + 0.3193815)                                       NDTR 410
      IF(X)1,2,2                                                        NDTR 420
    1 P=1.0-P                                                           NDTR 430
    2 RETURN                                                            NDTR 440
      END                                                               NDTR 450
