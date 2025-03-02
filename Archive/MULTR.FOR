C                                                                       MULT  10
C     ..................................................................MULT  20
C                                                                       MULT  30
C        SUBROUTINE MULTR                                               MULT  40
C                                                                       MULT  50
C        PURPOSE                                                        MULT  60
C           PERFORM A MULTIPLE LINEAR REGRESSION ANALYSIS FOR A         MULT  70
C           DEPENDENT VARIABLE AND A SET OF INDEPENDENT VARIABLES.  THISMULT  80
C           SUBROUTINE IS NORMALLY USED IN THE PERFORMANCE OF MULTIPLE  MULT  90
C           AND POLYNOMIAL REGRESSION ANALYSES.                         MULT 100
C                                                                       MULT 110
C        USAGE                                                          MULT 120
C           CALL MULTR (N,K,XBAR,STD,D,RX,RY,ISAVE,B,SB,T,ANS)          MULT 130
C                                                                       MULT 140
C        DESCRIPTION OF PARAMETERS                                      MULT 150
C           N     - NUMBER OF OBSERVATIONS.                             MULT 160
C           K     - NUMBER OF INDEPENDENT VARIABLES IN THIS REGRESSION. MULT 170
C           XBAR  - INPUT VECTOR OF LENGTH M CONTAINING MEANS OF ALL    MULT 180
C                   VARIABLES. M IS NUMBER OF VARIABLES IN OBSERVATIONS.MULT 190
C           STD   - INPUT VECTOR OF LENGTH M CONTAINING STANDARD DEVI-  MULT 200
C                   ATIONS OF ALL VARIABLES.                            MULT 210
C           D     - INPUT VECTOR OF LENGTH M CONTAINING THE DIAGONAL OF MULT 220
C                   THE MATRIX OF SUMS OF CROSS-PRODUCTS OF DEVIATIONS  MULT 230
C                   FROM MEANS FOR ALL VARIABLES.                       MULT 240
C           RX    - INPUT MATRIX (K X K) CONTAINING THE INVERSE OF      MULT 250
C                   INTERCORRELATIONS AMONG INDEPENDENT VARIABLES.      MULT 260
C           RY    - INPUT VECTOR OF LENGTH K CONTAINING INTERCORRELA-   MULT 270
C                   TIONS OF INDEPENDENT VARIABLES WITH DEPENDENT       MULT 280
C                   VARIABLE.                                           MULT 290
C           ISAVE - INPUT VECTOR OF LENGTH K+1 CONTAINING SUBSCRIPTS OF MULT 300
C                   INDEPENDENT VARIABLES IN ASCENDING ORDER.  THE      MULT 310
C                   SUBSCRIPT OF THE DEPENDENT VARIABLE IS STORED IN    MULT 320
C                   THE LAST, K+1, POSITION.                            MULT 330
C           B     - OUTPUT VECTOR OF LENGTH K CONTAINING REGRESSION     MULT 340
C                   COEFFICIENTS.                                       MULT 350
C           SB    - OUTPUT VECTOR OF LENGTH K CONTAINING STANDARD       MULT 360
C                   DEVIATIONS OF REGRESSION COEFFICIENTS.              MULT 370
C           T     - OUTPUT VECTOR OF LENGTH K CONTAINING T-VALUES.      MULT 380
C           ANS   - OUTPUT VECTOR OF LENGTH 10 CONTAINING THE FOLLOWING MULT 390
C                   INFORMATION..                                       MULT 400
C                   ANS(1)  INTERCEPT                                   MULT 410
C                   ANS(2)  MULTIPLE CORRELATION COEFFICIENT            MULT 420
C                   ANS(3)  STANDARD ERROR OF ESTIMATE                  MULT 430
C                   ANS(4)  SUM OF SQUARES ATTRIBUTABLE TO REGRES-      MULT 440
C                           SION (SSAR)                                 MULT 450
C                   ANS(5)  DEGREES OF FREEDOM ASSOCIATED WITH SSAR     MULT 460
C                   ANS(6)  MEAN SQUARE OF SSAR                         MULT 470
C                   ANS(7)  SUM OF SQUARES OF DEVIATIONS FROM REGRES-   MULT 480
C                           SION (SSDR)                                 MULT 490
C                   ANS(8)  DEGREES OF FREEDOM ASSOCIATED WITH SSDR     MULT 500
C                   ANS(9)  MEAN SQUARE OF SSDR                         MULT 510
C                   ANS(10) F-VALUE                                     MULT 520
C                                                                       MULT 530
C        REMARKS                                                        MULT 540
C           N MUST BE GREATER THAN K+1.                                 MULT 550
C                                                                       MULT 560
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MULT 570
C           NONE                                                        MULT 580
C                                                                       MULT 590
C        METHOD                                                         MULT 600
C           THE GAUSS-JORDAN METHOD IS USED IN THE SOLUTION OF THE      MULT 610
C           NORMAL EQUATIONS.  REFER TO W. W. COOLEY AND P. R. LOHNES,  MULT 620
C           'MULTIVARIATE PROCEDURES FOR THE BEHAVIORAL SCIENCES',      MULT 630
C           JOHN WILEY AND SONS, 1962, CHAPTER 3, AND B. OSTLE,         MULT 640
C           'STATISTICS IN RESEARCH', THE IOWA STATE COLLEGE PRESS,     MULT 650
C           1954, CHAPTER 8.                                            MULT 660
C                                                                       MULT 670
C     ..................................................................MULT 680
C                                                                       MULT 690
      SUBROUTINE MULTR (N,K,XBAR,STD,D,RX,RY,ISAVE,B,SB,T,ANS)          MULT 700
      DIMENSION XBAR(1),STD(1),D(1),RX(1),RY(1),ISAVE(1),B(1),SB(1),    MULT 710
     1          T(1),ANS(1)                                             MULT 720
C                                                                       MULT 730
C        ...............................................................MULT 740
C                                                                       MULT 750
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  MULT 760
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      MULT 770
C        STATEMENT WHICH FOLLOWS.                                       MULT 780
C                                                                       MULT 790
C     DOUBLE PRECISION XBAR,STD,D,RX,RY,B,SB,T,ANS,RM,BO,SSAR,SSDR,SY,  MULT 800
C    1                 FN,FK,SSARM,SSDRM,F                              MULT 810
C                                                                       MULT 820
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    MULT 830
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      MULT 840
C        ROUTINE.                                                       MULT 850
C                                                                       MULT 860
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO      MULT 870
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT AND ABS IN   MULT 880
C        STATEMENTS 122, 125, AND 135 MUST BE CHANGED TO DSQRT AND DABS.MULT 890
C                                                                       MULT 900
C        ...............................................................MULT 910
C                                                                       MULT 920
      MM=K+1                                                            MULT 930
C                                                                       MULT 940
C        BETA WEIGHTS                                                   MULT 950
C                                                                       MULT 960
      DO 100 J=1,K                                                      MULT 970
  100 B(J)=0.0                                                          MULT 980
      DO 110 J=1,K                                                      MULT 990
      L1=K*(J-1)                                                        MULT1000
      DO 110 I=1,K                                                      MULT1010
      L=L1+I                                                            MULT1020
  110 B(J)=B(J)+RY(I)*RX(L)                                             MULT1030
      RM=0.0                                                            MULT1040
      BO=0.0                                                            MULT1050
      L1=ISAVE(MM)                                                      MULT1060
C                                                                       MULT1070
C        COEFFICIENT OF DETERMINATION                                   MULT1080
C                                                                       MULT1090
      DO 120 I=1,K                                                      MULT1100
      RM=RM+B(I)*RY(I)                                                  MULT1110
C                                                                       MULT1120
C        REGRESSION COEFFICIENTS                                        MULT1130
C                                                                       MULT1140
      L=ISAVE(I)                                                        MULT1150
      B(I)=B(I)*(STD(L1)/STD(L))                                        MULT1160
C                                                                       MULT1170
C        INTERCEPT                                                      MULT1180
C                                                                       MULT1190
  120 BO=BO+B(I)*XBAR(L)                                                MULT1200
      BO=XBAR(L1)-BO                                                    MULT1210
C                                                                       MULT1220
C        SUM OF SQUARES ATTRIBUTABLE TO REGRESSION                      MULT1230
C                                                                       MULT1240
      SSAR=RM*D(L1)                                                     MULT1250
C                                                                       MULT1260
C        MULTIPLE CORRELATION COEFFICIENT                               MULT1270
C                                                                       MULT1280
  122 RM= SQRT( ABS(RM))                                                MULT1290
C                                                                       MULT1300
C        SUM OF SQUARES OF DEVIATIONS FROM REGRESSION                   MULT1310
C                                                                       MULT1320
      SSDR=D(L1)-SSAR                                                   MULT1330
C                                                                       MULT1340
C        VARIANCE OF ESTIMATE                                           MULT1350
C                                                                       MULT1360
      FN=N-K-1                                                          MULT1370
      SY=SSDR/FN                                                        MULT1380
C                                                                       MULT1390
C        STANDARD DEVIATIONS OF REGRESSION COEFFICIENTS                 MULT1400
C                                                                       MULT1410
      DO 130 J=1,K                                                      MULT1420
      L1=K*(J-1)+J                                                      MULT1430
      L=ISAVE(J)                                                        MULT1440
  125 SB(J)= SQRT( ABS((RX(L1)/D(L))*SY))                               MULT1450
C                                                                       MULT1460
C        COMPUTED T-VALUES                                              MULT1470
C                                                                       MULT1480
  130 T(J)=B(J)/SB(J)                                                   MULT1490
C                                                                       MULT1500
C        STANDARD ERROR OF ESTIMATE                                     MULT1510
C                                                                       MULT1520
  135 SY= SQRT( ABS(SY))                                                MULT1530
C                                                                       MULT1540
C        F VALUE                                                        MULT1550
C                                                                       MULT1560
      FK=K                                                              MULT1570
      SSARM=SSAR/FK                                                     MULT1580
      SSDRM=SSDR/FN                                                     MULT1590
      F=SSARM/SSDRM                                                     MULT1600
C                                                                       MULT1610
      ANS(1)=BO                                                         MULT1620
      ANS(2)=RM                                                         MULT1630
      ANS(3)=SY                                                         MULT1640
      ANS(4)=SSAR                                                       MULT1650
      ANS(5)=FK                                                         MULT1660
      ANS(6)=SSARM                                                      MULT1670
      ANS(7)=SSDR                                                       MULT1680
      ANS(8)=FN                                                         MULT1690
      ANS(9)=SSDRM                                                      MULT1700
      ANS(10)=F                                                         MULT1710
      RETURN                                                            MULT1720
      END                                                               MULT1730
