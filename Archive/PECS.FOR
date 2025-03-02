C                                                                       PECS  10
C     ..................................................................PECS  20
C                                                                       PECS  30
C        SUBROUTINE PECS                                                PECS  40
C                                                                       PECS  50
C        PURPOSE                                                        PECS  60
C           ECONOMIZATION OF A POLYNOMIAL FOR UNSYMMETRIC RANGE         PECS  70
C                                                                       PECS  80
C        USAGE                                                          PECS  90
C           CALL PECS (P,N,BOUND,EPS,TOL,WORK)                          PECS 100
C                                                                       PECS 110
C        DESCRIPTION OF PARAMETERS                                      PECS 120
C           P     - COEFFICIENT VECTOR OF GIVEN POLYNOMIAL              PECS 130
C           N     - DIMENSION OF COEFFICIENT VECTOR                     PECS 140
C           BOUND - RIGHT HAND BOUNDARY OF INTERVAL                     PECS 150
C           EPS   - INITIAL ERROR BOUND                                 PECS 160
C           TOL   - TOLERANCE FOR ERROR                                 PECS 170
C           WORK  - WORKING STORAGE OF DIMENSION N                      PECS 180
C                                                                       PECS 190
C        REMARKS                                                        PECS 200
C           THE INITIAL COEFFICIENT VECTOR P IS REPLACED BY THE         PECS 210
C           ECONOMIZED VECTOR.                                          PECS 220
C           THE INITIAL ERROR BOUND EPS IS REPLACED BY A FINAL          PECS 230
C           ERROR BOUND.                                                PECS 240
C           N IS REPLACED BY THE DIMENSION OF THE REDUCED POLYNOMIAL.   PECS 250
C           IN CASE OF AN ARBITRARY INTERVAL (XL,XR) IT IS NECESSARY    PECS 260
C           FIRST TO CALCULATE THE EXPANSION OF THE GIVEN POLYNOMIAL    PECS 270
C           WITH ARGUMENT X IN POWERS OF T = (X-XL).                    PECS 280
C           THIS IS ACCOMPLISHED THROUGH SUBROUTINE PCLD.               PECS 290
C           OPERATION IS BYPASSED IN CASE OF N LESS THAN 1.             PECS 300
C                                                                       PECS 310
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PECS 320
C           NONE                                                        PECS 330
C                                                                       PECS 340
C        METHOD                                                         PECS 350
C           SUBROUTINE PECS TAKES AN (N-1)ST DEGREE POLYNOMIAL          PECS 360
C           APPROXIMATION TO A FUNCTION F(X) VALID WITHIN A TOLERANCE   PECS 370
C           EPS OVER THE INTERVAL (0,BOUND) AND REDUCES IT IF POSSIBLE  PECS 380
C           TO A POLYNOMIAL OF LOWER DEGREE VALID WITHIN TOLERANCE      PECS 390
C           TOL.                                                        PECS 400
C           THE COEFFICIENT VECTOR OF THE N-TH SHIFTED CHEBYSHEV        PECS 410
C           POLYNOMIAL IS CALCULATED FROM THE RECURSION FORMULA         PECS 420
C           A(K) = -A(K+1)*K*L*(2*K-1)/(2*(N+K-1)*(N-K+1)).             PECS 430
C           REFERENCE                                                   PECS 440
C           K. A. BRONS, ALGORITHM 37, TELESCOPE 1, CACM VOL. 4, 1961,  PECS 450
C           NO. 3, PP. 151.                                             PECS 460
C                                                                       PECS 470
C     ..................................................................PECS 480
C                                                                       PECS 490
      SUBROUTINE PECS(P,N,BOUND,EPS,TOL,WORK)                           PECS 500
C                                                                       PECS 510
      DIMENSION P(1),WORK(1)                                            PECS 520
      FL=BOUND*0.5                                                      PECS 530
C                                                                       PECS 540
C     TEST OF DIMENSION                                                 PECS 550
C                                                                       PECS 560
    1 IF(N-1)2,3,6                                                      PECS 570
    2 RETURN                                                            PECS 580
    3 IF(EPS+ABS(P(1))-TOL)4,4,5                                        PECS 590
    4 N=0                                                               PECS 600
      EPS=EPS+ABS(P(1))                                                 PECS 610
    5 RETURN                                                            PECS 620
C                                                                       PECS 630
C     CALCULATE EXPANSION OF CHEBYSHEV POLYNOMIAL                       PECS 640
C                                                                       PECS 650
    6 NEND=N-1                                                          PECS 660
      WORK(N)=-P(N)                                                     PECS 670
      DO 7 J=1,NEND                                                     PECS 680
      K=N-J                                                             PECS 690
      FN=(NEND-1+K)*(N-K)                                               PECS 700
      FK=K*(K+K-1)                                                      PECS 710
    7 WORK(K)=-WORK(K+1)*FK*FL/FN                                       PECS 720
C                                                                       PECS 730
C        TEST FOR FEASIBILITY OF REDUCTION                              PECS 740
C                                                                       PECS 750
      FN=ABS(WORK(1))                                                   PECS 760
      IF(EPS+FN-TOL)8,8,5                                               PECS 770
C                                                                       PECS 780
C     REDUCE POLYNOMIAL                                                 PECS 790
C                                                                       PECS 800
    8 EPS=EPS+FN                                                        PECS 810
      N=NEND                                                            PECS 820
      DO 9 J=1,NEND                                                     PECS 830
    9 P(J)=P(J)+WORK(J)                                                 PECS 840
      GOTO 1                                                            PECS 850
      END                                                               PECS 860
