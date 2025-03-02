C                                                                       DPCS  10
C     ..................................................................DPCS  20
C                                                                       DPCS  30
C        SUBROUTINE DPECS                                               DPCS  40
C                                                                       DPCS  50
C        PURPOSE                                                        DPCS  60
C           ECONOMIZATION OF A POLYNOMIAL FOR UNSYMMETRIC RANGE         DPCS  70
C                                                                       DPCS  80
C        USAGE                                                          DPCS  90
C           CALL DPECS(P,N,BOUND,EPS,TOL,WORK)                          DPCS 100
C                                                                       DPCS 110
C        DESCRIPTION OF PARAMETERS                                      DPCS 120
C           P     - DOUBLE PRECISION COEFFICIENT VECTOR OF GIVEN        DPCS 130
C                   POLYNOMIAL                                          DPCS 140
C           N     - DIMENSION OF COEFFICIENT VECTOR P                   DPCS 150
C           BOUND - SINGLE PRECISION RIGHT HAND BOUNDARY OF INTERVAL    DPCS 160
C           EPS   - SINGLE PRECISION INITIAL ERROR BOUND                DPCS 170
C           TOL   - SINGLE PRECISION TOLERANCE FOR ERROR                DPCS 180
C           WORK  - DOUBLE PRECISION WORKING STORAGE OF DIMENSION N     DPCS 190
C                                                                       DPCS 200
C        REMARKS                                                        DPCS 210
C           THE INITIAL COEFFICIENT VECTOR P IS REPLACED BY THE         DPCS 220
C           ECONOMIZED VECTOR.                                          DPCS 230
C           THE INITIAL ERROR BOUND EPS IS REPLACED BY A FINAL          DPCS 240
C           ERROR BOUND.                                                DPCS 250
C           N IS REPLACED BY THE DIMENSION OF THE REDUCED POLYNOMIAL.   DPCS 260
C           IN CASE OF AN ARBITRARY INTERVAL (XL,XR) IT IS NECESSARY    DPCS 270
C           FIRST TO CALCULATE THE EXPANSION OF THE GIVEN POLYNOMIAL    DPCS 280
C           WITH ARGUMENT X IN POWERS OF T = (X-XL).                    DPCS 290
C           THIS IS ACCOMPLISHED THROUGH SUBROUTINE DPCLD.              DPCS 300
C           OPERATION IS BYPASSED IN CASE OF N LESS THAN 1.             DPCS 310
C                                                                       DPCS 320
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DPCS 330
C           NONE                                                        DPCS 340
C                                                                       DPCS 350
C        METHOD                                                         DPCS 360
C           SUBROUTINE DPECS TAKES AN (N-1)ST DEGREE POLYNOMIAL         DPCS 370
C           APPROXIMATION TO A FUNCTION F(X) VALID WITHIN A TOLERANCE   DPCS 380
C           EPS OVER THE INTERVAL (0,BOUND) AND REDUCES IT IF POSSIBLE  DPCS 390
C           TO A POLYNOMIAL OF LOWER DEGREE VALID WITHIN TOLERANCE      DPCS 400
C           TOL.                                                        DPCS 410
C           THE COEFFICIENT VECTOR OF THE N-TH SHIFTED CHEBYSHEV        DPCS 420
C           POLYNOMIAL IS CALCULATED FROM THE RECURSION FORMULA         DPCS 430
C           A(K) = -A(K+1)*K*L*(2*K-1)/(2*(N+K-1)*(N-K+1)).             DPCS 440
C           REFERENCE                                                   DPCS 450
C           K. A. BRONS, ALGORITHM 37, TELESCOPE 1, CACM VOL. 4, 1961,  DPCS 460
C           NO. 3, PP. 151.                                             DPCS 470
C                                                                       DPCS 480
C     ..................................................................DPCS 490
C                                                                       DPCS 500
      SUBROUTINE DPECS(P,N,BOUND,EPS,TOL,WORK)                          DPCS 510
C                                                                       DPCS 520
      DIMENSION P(1),WORK(1)                                            DPCS 530
      DOUBLE PRECISION P,WORK                                           DPCS 540
C                                                                       DPCS 550
      FL=BOUND*0.5                                                      DPCS 560
C                                                                       DPCS 570
C        TEST OF DIMENSION                                              DPCS 580
C                                                                       DPCS 590
    1 IF(N-1)2,3,6                                                      DPCS 600
    2 RETURN                                                            DPCS 610
C                                                                       DPCS 620
    3 IF(EPS+ABS(SNGL(P(1)))-TOL)4,4,5                                  DPCS 630
    4 N=0                                                               DPCS 640
      EPS=EPS+ABS(SNGL(P(1)))                                           DPCS 650
    5 RETURN                                                            DPCS 660
C                                                                       DPCS 670
C        CALCULATE EXPANSION OF CHEBYSHEV POLYNOMIAL                    DPCS 680
C                                                                       DPCS 690
    6 NEND=N-1                                                          DPCS 700
      WORK(N)=-P(N)                                                     DPCS 710
      DO 7 J=1,NEND                                                     DPCS 720
      K=N-J                                                             DPCS 730
      FN=(NEND-1+K)*(N-K)                                               DPCS 740
      FK=K*(K+K-1)                                                      DPCS 750
    7 WORK(K)=-WORK(K+1)*DBLE(FK)*DBLE(FL)/DBLE(FN)                     DPCS 760
C                                                                       DPCS 770
C        TEST FOR FEASIBILITY OF REDUCTION                              DPCS 780
C                                                                       DPCS 790
      FN=DABS(WORK(1))                                                  DPCS 800
      IF(EPS+FN-TOL)8,8,5                                               DPCS 810
C                                                                       DPCS 820
C        REDUCE POLYNOMIAL                                              DPCS 830
C                                                                       DPCS 840
    8 EPS=EPS+FN                                                        DPCS 850
      N=NEND                                                            DPCS 860
      DO 9 J=1,NEND                                                     DPCS 870
    9 P(J)=P(J)+WORK(J)                                                 DPCS 880
      GOTO 1                                                            DPCS 890
      END                                                               DPCS 900
