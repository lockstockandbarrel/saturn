C                                                                       DPCN  10
C     ..................................................................DPCN  20
C                                                                       DPCN  30
C        SUBROUTINE DPECN                                               DPCN  40
C                                                                       DPCN  50
C        PURPOSE                                                        DPCN  60
C           ECONOMIZE A POLYNOMIAL FOR SYMMETRIC RANGE                  DPCN  70
C                                                                       DPCN  80
C        USAGE                                                          DPCN  90
C           CALL DPECN(P,N,BOUND,EPS,TOL,WORK)                          DPCN 100
C                                                                       DPCN 110
C        DESCRIPTION OF PARAMETERS                                      DPCN 120
C           P     - DOUBLE PRECISION COEFFICIENT VECTOR OF GIVEN        DPCN 130
C                   POLYNOMIAL                                          DPCN 140
C                   ON RETURN P CONTAINS THE ECONOMIZED POLYNOMIAL      DPCN 150
C           N     - DIMENSION OF COEFFICIENT VECTOR P                   DPCN 160
C                   ON RETURN N CONTAINS DIMENSION OF ECONOMIZED        DPCN 170
C                   POLYNOMIAL                                          DPCN 180
C           BOUND - SINGLE PRECISION RIGHT HAND BOUNDARY OF RANGE       DPCN 190
C           EPS   - SINGLE PRECISION INITIAL ERROR BOUND                DPCN 200
C                   ON RETURN EPS CONTAINS AN ERROR BOUND FOR THE       DPCN 210
C                   ECONOMIZED POLYNOMIAL                               DPCN 220
C           TOL   - SINGLE PRECISION TOLERANCE FOR ERROR                DPCN 230
C                   FINAL VALUE OF EPS MUST BE LESS THAN TOL            DPCN 240
C           WORK  - DOUBLE PRECISION WORKING STORAGE OF DIMENSION N     DPCN 250
C                   (STARTING VALUE OF N RATHER THAN FINAL VALUE)       DPCN 260
C                                                                       DPCN 270
C        REMARKS                                                        DPCN 280
C           THE OPERATION IS BYPASSED IN CASE OF N LESS THAN 1.         DPCN 290
C           IN CASE OF AN ARBITRARY INTERVAL (XL,XR) IT IS NECESSARY    DPCN 300
C           FIRST TO CALCULATE THE EXPANSION OF THE GIVEN POLYNOMIAL    DPCN 310
C           WITH ARGUMENT X IN POWERS OF T = (X-(XR-XL)/2).             DPCN 320
C           THIS IS ACCOMPLISHED THROUGH SUBROUTINE DPCLD.              DPCN 330
C                                                                       DPCN 340
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DPCN 350
C           NONE                                                        DPCN 360
C                                                                       DPCN 370
C        METHOD                                                         DPCN 380
C           SUBROUTINE DPECN TAKES AN (N-1)ST DEGREE POLYNOMIAL         DPCN 390
C           APPROXIMATION TO A FUNCTION F(X) VALID WITHIN A TOLERANCE   DPCN 400
C           EPS OVER THE INTERVAL (-BOUND,BOUND) AND REDUCES IT IF      DPCN 410
C           POSSIBLE TO A POLYNOMIAL OF LOWER DEGREE VALID WITHIN       DPCN 420
C           THE GIVEN TOLERANCE TOL.                                    DPCN 430
C           THE INITIAL COEFFICIENT VECTOR P IS REPLACED BY THE FINAL   DPCN 440
C           VECTOR. THE INITIAL ERROR BOUND EPS IS REPLACED BY A FINAL  DPCN 450
C           ERROR BOUND.                                                DPCN 460
C           N IS REPLACED BY THE DIMENSION OF THE REDUCED POLYNOMIAL.   DPCN 470
C           THE COEFFICIENT VECTOR OF THE N-TH CHEBYSHEV POLYNOMIAL     DPCN 480
C           IS CALCULATED FROM THE RECURSION FORMULA                    DPCN 490
C           A(K-1)=-A(K+1)*K*L*L*(K-1)/((N+K-2)*(N-K+2))                DPCN 500
C           REFERENCE                                                   DPCN 510
C           K. A. BRONS, ALGORITHM 38, TELESCOPE 2, CACM VOL. 4, 1961,  DPCN 520
C           NO. 3, PP. 151-152.                                         DPCN 530
C                                                                       DPCN 540
C     ..................................................................DPCN 550
C                                                                       DPCN 560
      SUBROUTINE DPECN(P,N,BOUND,EPS,TOL,WORK)                          DPCN 570
C                                                                       DPCN 580
      DIMENSION P(1),WORK(1)                                            DPCN 590
      DOUBLE PRECISION P,WORK                                           DPCN 600
C                                                                       DPCN 610
      FL=BOUND*BOUND                                                    DPCN 620
C                                                                       DPCN 630
C        TEST OF DIMENSION                                              DPCN 640
C                                                                       DPCN 650
    1 IF(N-1)2,3,6                                                      DPCN 660
    2 RETURN                                                            DPCN 670
C                                                                       DPCN 680
    3 IF(EPS+ABS(SNGL(P(1)))-TOL)4,4,5                                  DPCN 690
    4 N=0                                                               DPCN 700
      EPS=EPS+ABS(SNGL(P(1)))                                           DPCN 710
    5 RETURN                                                            DPCN 720
C                                                                       DPCN 730
C        CALCULATE EXPANSION OF CHEBYSHEV POLYNOMIAL                    DPCN 740
C                                                                       DPCN 750
    6 NEND=N-2                                                          DPCN 760
      WORK(N)=-P(N)                                                     DPCN 770
      DO 7 J=1,NEND,2                                                   DPCN 780
      K=N-J                                                             DPCN 790
      FN=(NEND-1+K)*(NEND+3-K)                                          DPCN 800
      FK=K*(K-1)                                                        DPCN 810
    7 WORK(K-1)=-WORK(K+1)*DBLE(FK*FL/FN)                               DPCN 820
C                                                                       DPCN 830
C        TEST FOR FEASIBILITY OF REDUCTION                              DPCN 840
C                                                                       DPCN 850
      IF(K-2)8,8,9                                                      DPCN 860
    8 FN=DABS(WORK(1))                                                  DPCN 870
      GOTO 10                                                           DPCN 880
    9 FN=N-1                                                            DPCN 890
      FN=ABS(SNGL(WORK(2))/FN)                                          DPCN 900
   10 IF(EPS+FN-TOL)11,11,5                                             DPCN 910
C                                                                       DPCN 920
C        REDUCE POLYNOMIAL                                              DPCN 930
C                                                                       DPCN 940
   11 EPS=EPS+FN                                                        DPCN 950
      N=N-1                                                             DPCN 960
      DO 12 J=K,N,2                                                     DPCN 970
   12 P(J-1)=P(J-1)+WORK(J-1)                                           DPCN 980
      GOTO 1                                                            DPCN 990
      END                                                               DPCN1000
