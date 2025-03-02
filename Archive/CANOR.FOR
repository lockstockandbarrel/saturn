C                                                                       CANO  10
C     ..................................................................CANO  20
C                                                                       CANO  30
C        SUBROUTINE CANOR                                               CANO  40
C                                                                       CANO  50
C        PURPOSE                                                        CANO  60
C           COMPUTE THE CANONICAL CORRELATIONS BETWEEN TWO SETS OF      CANO  70
C           VARIABLES.  CANOR IS NORMALLY PRECEDED BY A CALL TO SUBROU- CANO  80
C           TINE CORRE.                                                 CANO  90
C                                                                       CANO 100
C        USAGE                                                          CANO 110
C           CALL CANOR (N,MP,MQ,RR,ROOTS,WLAM,CANR,CHISQ,NDF,COEFR,     CANO 120
C                       COEFL,R)                                        CANO 130
C                                                                       CANO 140
C        DESCRIPTION OF PARAMETERS                                      CANO 150
C           N     - NUMBER OF OBSERVATIONS                              CANO 160
C           MP    - NUMBER OF LEFT HAND VARIABLES                       CANO 170
C           MQ    - NUMBER OF RIGHT HAND VARIABLES                      CANO 180
C           RR    - INPUT MATRIX (ONLY UPPER TRIANGULAR PORTION OF THE  CANO 190
C                   SYMMETRIC MATRIX OF M X M, WHERE M = MP + MQ)       CANO 200
C                   CONTAINING CORRELATION COEFFICIENTS.  (STORAGE MODE CANO 210
C                   OF 1)                                               CANO 220
C           ROOTS - OUTPUT VECTOR OF LENGTH MQ CONTAINING EIGENVALUES   CANO 230
C                   COMPUTED IN THE NROOT SUBROUTINE.                   CANO 240
C           WLAM  - OUTPUT VECTOR OF LENGTH MQ CONTAINING LAMBDA.       CANO 250
C           CANR  - OUTPUT VECTOR OF LENGTH MQ CONTAINING CANONICAL     CANO 260
C                   CORRELATIONS.                                       CANO 270
C           CHISQ - OUTPUT VECTOR OF LENGTH MQ CONTAINING THE           CANO 280
C                   VALUES OF CHI-SQUARES.                              CANO 290
C           NDF   - OUTPUT VECTOR OF LENGTH MQ CONTAINING THE DEGREES   CANO 300
C                   OF FREEDOM ASSOCIATED WITH CHI-SQUARES.             CANO 310
C           COEFR - OUTPUT MATRIX (MQ X MQ) CONTAINING MQ SETS OF       CANO 320
C                   RIGHT HAND COEFFICIENTS COLUMNWISE.                 CANO 330
C           COEFL - OUTPUT MATRIX (MP X MQ) CONTAINING MQ SETS OF       CANO 340
C                   LEFT HAND COEFFICIENTS COLUMNWISE.                  CANO 350
C           R     - WORK MATRIX (M X M)                                 CANO 360
C                                                                       CANO 370
C        REMARKS                                                        CANO 380
C           THE NUMBER OF LEFT HAND VARIABLES (MP) SHOULD BE GREATER    CANO 390
C           THAN OR EQUAL TO THE NUMBER OF RIGHT HAND VARIABLES (MQ).   CANO 400
C           THE VALUES OF CANONICAL CORRELATION, LAMBDA, CHI-SQUARE,    CANO 410
C           DEGREES OF FREEDOM, AND CANONICAL COEFFICIENTS ARE COMPUTED CANO 420
C           ONLY FOR THOSE EIGENVALUES IN ROOTS WHICH ARE GREATER THAN  CANO 430
C           ZERO.                                                       CANO 440
C                                                                       CANO 450
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  CANO 460
C           MINV                                                        CANO 470
C           NROOT  (WHICH, IN TURN, CALLS THE SUBROUTINE EIGEN.)        CANO 480
C                                                                       CANO 490
C        METHOD                                                         CANO 500
C           REFER TO W. W. COOLEY AND P. R. LOHNES, 'MULTIVARIATE PRO-  CANO 510
C           CEDURES FOR THE BEHAVIORAL SCIENCES', JOHN WILEY AND SONS,  CANO 520
C           1962, CHAPTER 3.                                            CANO 530
C                                                                       CANO 540
C     ..................................................................CANO 550
C                                                                       CANO 560
      SUBROUTINE CANOR (N,MP,MQ,RR,ROOTS,WLAM,CANR,CHISQ,NDF,COEFR,     CANO 570
     1                  COEFL,R)                                        CANO 580
      DIMENSION RR(1),ROOTS(1),WLAM(1),CANR(1),CHISQ(1),NDF(1),COEFR(1),CANO 590
     1          COEFL(1),R(1)                                           CANO 600
C                                                                       CANO 610
C        ...............................................................CANO 620
C                                                                       CANO 630
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  CANO 640
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      CANO 650
C        STATEMENT WHICH FOLLOWS.                                       CANO 660
C                                                                       CANO 670
C     DOUBLE PRECISION RR,ROOTS,WLAM,CANR,CHISQ,COEFR,COEFL,R,DET,SUM   CANO 680
C                                                                       CANO 690
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    CANO 700
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      CANO 710
C        ROUTINE.                                                       CANO 720
C                                                                       CANO 730
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO      CANO 740
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  SQRT IN STATEMENT CANO 750
C        165 MUST BE CHANGED TO DSQRT.  ALOG IN STATEMENT 175 MUST BE   CANO 760
C        CHANGED TO DLOG.                                               CANO 770
C                                                                       CANO 780
C        ...............................................................CANO 790
C                                                                       CANO 800
C     PARTITION INTERCORRELATIONS AMONG LEFT HAND VARIABLES, BETWEEN    CANO 810
C     LEFT AND RIGHT HAND VARIABLES, AND AMONG RIGHT HAND VARIABLES.    CANO 820
C                                                                       CANO 830
      M=MP+MQ                                                           CANO 840
      N1=0                                                              CANO 850
      DO 105 I=1,M                                                      CANO 860
      DO 105 J=1,M                                                      CANO 870
      IF(I-J) 102, 103, 103                                             CANO 880
  102 L=I+(J*J-J)/2                                                     CANO 890
      GO TO 104                                                         CANO 900
  103 L=J+(I*I-I)/2                                                     CANO 910
  104 N1=N1+1                                                           CANO 920
  105 R(N1)=RR(L)                                                       CANO 930
      L=MP                                                              CANO 940
      DO 108 J=2,MP                                                     CANO 950
      N1=M*(J-1)                                                        CANO 960
      DO 108 I=1,MP                                                     CANO 970
      L=L+1                                                             CANO 980
      N1=N1+1                                                           CANO 990
  108 R(L)=R(N1)                                                        CANO1000
      N2=MP+1                                                           CANO1010
      L=0                                                               CANO1020
      DO 110 J=N2,M                                                     CANO1030
      N1=M*(J-1)                                                        CANO1040
      DO 110 I=1,MP                                                     CANO1050
      L=L+1                                                             CANO1060
      N1=N1+1                                                           CANO1070
  110 COEFL(L)=R(N1)                                                    CANO1080
      L=0                                                               CANO1090
      DO 120 J=N2,M                                                     CANO1100
      N1=M*(J-1)+MP                                                     CANO1110
      DO 120 I=N2,M                                                     CANO1120
      L=L+1                                                             CANO1130
      N1=N1+1                                                           CANO1140
  120 COEFR(L)=R(N1)                                                    CANO1150
C                                                                       CANO1160
C     SOLVE THE CANONICAL EQUATION                                      CANO1170
C                                                                       CANO1180
      L=MP*MP+1                                                         CANO1190
      K=L+MP                                                            CANO1200
      CALL MINV (R,MP,DET,R(L),R(K))                                    CANO1210
C                                                                       CANO1220
C        CALCULATE T = INVERSE OF R11 * R12                             CANO1230
C                                                                       CANO1240
      DO 140 I=1,MP                                                     CANO1250
      N2=0                                                              CANO1260
      DO 130 J=1,MQ                                                     CANO1270
      N1=I-MP                                                           CANO1280
      ROOTS(J)=0.0                                                      CANO1290
      DO 130 K=1,MP                                                     CANO1300
      N1=N1+MP                                                          CANO1310
      N2=N2+1                                                           CANO1320
  130 ROOTS(J)=ROOTS(J)+R(N1)*COEFL(N2)                                 CANO1330
      L=I-MP                                                            CANO1340
      DO 140 J=1,MQ                                                     CANO1350
      L=L+MP                                                            CANO1360
  140 R(L)=ROOTS(J)                                                     CANO1370
C                                                                       CANO1380
C        CALCULATE A = R21 * T                                          CANO1390
C                                                                       CANO1400
      L=MP*MQ                                                           CANO1410
      N3=L+1                                                            CANO1420
      DO 160 J=1,MQ                                                     CANO1430
      N1=0                                                              CANO1440
      DO 160 I=1,MQ                                                     CANO1450
      N2=MP*(J-1)                                                       CANO1460
      SUM=0.0                                                           CANO1470
      DO 150 K=1,MP                                                     CANO1480
      N1=N1+1                                                           CANO1490
      N2=N2+1                                                           CANO1500
  150 SUM=SUM+COEFL(N1)*R(N2)                                           CANO1510
      L=L+1                                                             CANO1520
  160 R(L)=SUM                                                          CANO1530
C                                                                       CANO1540
C        CALCULATE EIGENVALUES WITH ASSOCIATED EIGENVECTORS OF THE      CANO1550
C        INVERSE OF R22 * A                                             CANO1560
C                                                                       CANO1570
      L=L+1                                                             CANO1580
      CALL NROOT (MQ,R(N3),COEFR,ROOTS,R(L))                            CANO1590
C                                                                       CANO1600
C     FOR EACH VALUE OF I = 1, 2, ..., MQ, CALCULATE THE FOLLOWING      CANO1610
C     STATISTICS                                                        CANO1620
C                                                                       CANO1630
      DO 210 I=1,MQ                                                     CANO1640
C                                                                       CANO1650
C        TEST WHETHER EIGENVALUE IS GREATER THAN ZERO                   CANO1660
C                                                                       CANO1670
      IF(ROOTS(I)) 220, 220, 165                                        CANO1680
C                                                                       CANO1690
C        CANONICAL CORRELATION                                          CANO1700
C                                                                       CANO1710
  165 CANR(I)= SQRT(ROOTS(I))                                           CANO1720
C                                                                       CANO1730
C        CHI-SQUARE                                                     CANO1740
C                                                                       CANO1750
      WLAM(I)=1.0                                                       CANO1760
      DO 170 J=I,MQ                                                     CANO1770
  170 WLAM(I)=WLAM(I)*(1.0-ROOTS(J))                                    CANO1780
      FN=N                                                              CANO1790
      FMP=MP                                                            CANO1800
      FMQ=MQ                                                            CANO1810
  175 CHISQ(I)=-(FN-0.5*(FMP+FMQ+1.0))*ALOG(WLAM(I))                    CANO1820
C                                                                       CANO1830
C        DEGREES OF FREEDOM FOR CHI-SQUARE                              CANO1840
C                                                                       CANO1850
      N1=I-1                                                            CANO1860
      NDF(I)=(MP-N1)*(MQ-N1)                                            CANO1870
C                                                                       CANO1880
C        I-TH SET OF RIGHT HAND COEFFICIENTS                            CANO1890
C                                                                       CANO1900
      N1=MQ*(I-1)                                                       CANO1910
      N2=MQ*(I-1)+L-1                                                   CANO1920
      DO 180 J=1,MQ                                                     CANO1930
      N1=N1+1                                                           CANO1940
      N2=N2+1                                                           CANO1950
  180 COEFR(N1)=R(N2)                                                   CANO1960
C                                                                       CANO1970
C        I-TH SET OF LEFT HAND COEFFICIENTS                             CANO1980
C                                                                       CANO1990
      DO 200 J=1,MP                                                     CANO2000
      N1=J-MP                                                           CANO2010
      N2=MQ*(I-1)                                                       CANO2020
      K=MP*(I-1)+J                                                      CANO2030
      COEFL(K)=0.0                                                      CANO2040
      DO 190 JJ=1,MQ                                                    CANO2050
      N1=N1+MP                                                          CANO2060
      N2=N2+1                                                           CANO2070
  190 COEFL(K)=COEFL(K)+R(N1)*COEFR(N2)                                 CANO2080
  200 COEFL(K)=COEFL(K)/CANR(I)                                         CANO2090
  210 CONTINUE                                                          CANO2100
  220 RETURN                                                            CANO2110
      END                                                               CANO2120
