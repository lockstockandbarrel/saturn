C                                                                       PECN  10
C     ..................................................................PECN  20
C                                                                       PECN  30
C        SUBROUTINE PECN                                                PECN  40
C                                                                       PECN  50
C        PURPOSE                                                        PECN  60
C           ECONOMIZE A POLYNOMIAL FOR SYMMETRIC RANGE                  PECN  70
C                                                                       PECN  80
C        USAGE                                                          PECN  90
C           CALL PECN (P,N,BOUND,EPS,TOL,WORK)                          PECN 100
C                                                                       PECN 110
C        DESCRIPTION OF PARAMETERS                                      PECN 120
C           P     - COEFFICIENT VECTOR OF GIVEN POLYNOMIAL              PECN 130
C                   ON RETURN P CONTAINS THE ECONOMIZED POLYNOMIAL      PECN 140
C           N     - DIMENSION OF COEFFICIENT VECTOR P                   PECN 150
C                   ON RETURN N CONTAINS DIMENSION OF ECONOMIZED        PECN 160
C                   POLYNOMIAL                                          PECN 170
C           BOUND - RIGHT HAND BOUNDARY OF RANGE                        PECN 180
C           EPS   - INITIAL ERROR BOUND                                 PECN 190
C                   ON RETURN EPS CONTAINS AN ERROR BOUND FOR THE       PECN 200
C                   ECONOMIZED POLYNOMIAL                               PECN 210
C           TOL   - TOLERANCE FOR ERROR                                 PECN 220
C                   FINAL VALUE OF EPS MUST BE LESS THAN TOL            PECN 230
C           WORK  - WORKING STORAGE OF DIMENSION N (STARTING VALUE      PECN 240
C                   OF N RATHER THAN FINAL VALUE)                       PECN 250
C                                                                       PECN 260
C        REMARKS                                                        PECN 270
C           THE OPERATION IS BYPASSED IN CASE OF N LESS THAN 1.         PECN 280
C           IN CASE OF AN ARBITRARY INTERVAL (XL,XR) IT IS NECESSARY    PECN 290
C           FIRST TO CALCULATE THE EXPANSION OF THE GIVEN POLYNOMIAL    PECN 300
C           WITH ARGUMENT X IN POWERS OF T = (X-(XR-XL)/2).             PECN 310
C           THIS IS ACCOMPLISHED THROUGH SUBROUTINE PCLD.               PECN 320
C                                                                       PECN 330
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PECN 340
C           NONE                                                        PECN 350
C                                                                       PECN 360
C        METHOD                                                         PECN 370
C           SUBROUTINE PECN TAKES AN (N-1)ST DEGREE POLYNOMIAL          PECN 380
C           APPROXIMATION TO A FUNCTION F(X) VALID WITHIN A TOLERANCE   PECN 390
C           EPS OVER THE INTERVAL (-BOUND,BOUND) AND REDUCES IT IF      PECN 400
C           POSSIBLE TO A POLYNOMIAL OF LOWER DEGREE VALID WITHIN       PECN 410
C           THE GIVEN TOLERANCE TOL.                                    PECN 420
C           THE INITIAL COEFFICIENT VECTOR P IS REPLACED BY THE FINAL   PECN 430
C           VECTOR. THE INITIAL ERROR BOUND EPS IS REPLACED BY A FINAL  PECN 440
C           ERROR BOUND.                                                PECN 450
C           N IS REPLACED BY THE DIMENSION OF THE REDUCED POLYNOMIAL.   PECN 460
C           THE COEFFICIENT VECTOR OF THE N-TH CHEBYSHEV POLYNOMIAL     PECN 470
C           IS CALCULATED FROM THE RECURSION FORMULA                    PECN 480
C           A(K-1)=-A(K+1)*K*L*L*(K-1)/((N+K-2)*(N-K+2))                PECN 490
C           REFERENCE                                                   PECN 500
C           K. A. BRONS, ALGORITHM 38, TELESCOPE 2, CACM VOL. 4, 1961,  PECN 510
C           NO. 3, PP. 151-152.                                         PECN 520
C                                                                       PECN 530
C     ..................................................................PECN 540
C                                                                       PECN 550
      SUBROUTINE PECN(P,N,BOUND,EPS,TOL,WORK)                           PECN 560
C                                                                       PECN 570
      DIMENSION P(1),WORK(1)                                            PECN 580
      FL=BOUND*BOUND                                                    PECN 590
C                                                                       PECN 600
C     TEST OF DIMENSION                                                 PECN 610
C                                                                       PECN 620
    1 IF(N-1)2,3,6                                                      PECN 630
    2 RETURN                                                            PECN 640
    3 IF(EPS+ABS(P(1))-TOL)4,4,5                                        PECN 650
    4 N=0                                                               PECN 660
      EPS=EPS+ABS(P(1))                                                 PECN 670
    5 RETURN                                                            PECN 680
C                                                                       PECN 690
C     CALCULATE EXPANSION OF CHEBYSHEV POLYNOMIAL                       PECN 700
C                                                                       PECN 710
    6 NEND=N-2                                                          PECN 720
      WORK(N)=-P(N)                                                     PECN 730
      DO 7 J=1,NEND,2                                                   PECN 740
      K=N-J                                                             PECN 750
      FN=(NEND-1+K)*(NEND+3-K)                                          PECN 760
      FK=K*(K-1)                                                        PECN 770
    7 WORK(K-1)=-WORK(K+1)*FK*FL/FN                                     PECN 780
C                                                                       PECN 790
C     TEST FOR FEASIBILITY OF REDUCTION                                 PECN 800
C                                                                       PECN 810
      IF(K-2)8,8,9                                                      PECN 820
    8 FN=ABS(WORK(1))                                                   PECN 830
      GOTO 10                                                           PECN 840
    9 FN=N-1                                                            PECN 850
      FN=ABS(WORK(2)/FN)                                                PECN 860
   10 IF(EPS+FN-TOL)11,11,5                                             PECN 870
C                                                                       PECN 880
C     REDUCE POLYNOMIAL                                                 PECN 890
C                                                                       PECN 900
   11 EPS=EPS+FN                                                        PECN 910
      N=N-1                                                             PECN 920
      DO 12 J=K,N,2                                                     PECN 930
   12 P(J-1)=P(J-1)+WORK(J-1)                                           PECN 940
      GOTO 1                                                            PECN 950
      END                                                               PECN 960
