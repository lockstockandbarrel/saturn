C                                                                       KLM2  10
C     ..................................................................KLM2  20
C                                                                       KLM2  30
C        SUBROUTINE KOLM2                                               KLM2  40
C                                                                       KLM2  50
C        PURPOSE                                                        KLM2  60
C                                                                       KLM2  70
C           TESTS THE DIFFERENCE BETWEEN TWO SAMPLE DISTRIBUTION        KLM2  80
C           FUNCTIONS USING THE KOLMOGOROV-SMIRNOV TEST                 KLM2  90
C                                                                       KLM2 100
C        USAGE                                                          KLM2 110
C           CALL KOLM2(X,Y,N,M,Z,PROB)                                  KLM2 120
C                                                                       KLM2 130
C        DESCRIPTION OF PARAMETERS                                      KLM2 140
C           X    - INPUT VECTOR OF N INDEPENDENT OBSERVATIONS.  ON      KLM2 150
C                  RETURN FROM KOLM2, X HAS BEEN SORTED INTO A          KLM2 160
C                  MONOTONIC NON-DECREASING SEQUENCE.                   KLM2 170
C           Y    - INPUT VECTOR OF M INDEPENDENT OBSERVATIONS.  ON      KLM2 180
C                  RETURN FROM KOLM2, Y HAS BEEN SORTED INTO A          KLM2 190
C                  MONOTONIC NON-DECREASING SEQUENCE.                   KLM2 200
C           N    - NUMBER OF OBSERVATIONS IN X                          KLM2 210
C           M    - NUMBER OF OBSERVATIONS IN Y                          KLM2 220
C           Z    - OUTPUT VARIABLE CONTAINING THE GREATEST VALUE WITH   KLM2 230
C                  RESPECT TO THE SPECTRUM OF X AND Y OF                KLM2 240
C                  SQRT((M*N)/(M+N))*ABS(FN(X)-GM(Y)) WHERE             KLM2 250
C                  FN(X) IS THE EMPIRICAL DISTRIBUTION FUNCTION OF THE  KLM2 260
C                  SET (X) AND GM(Y) IS THE EMPIRICAL DISTRIBUTION      KLM2 270
C                  FUNCTION OF THE SET (Y).                             KLM2 280
C           PROB - OUTPUT VARIABLE CONTAINING THE PROBABILITY OF        KLM2 290
C                  THE STATISTIC BEING GREATER THAN OR EQUAL TO Z IF    KLM2 300
C                  THE HYPOTHESIS THAT X AND Y ARE FROM THE SAME PDF IS KLM2 310
C                  TRUE.  E.G., PROB= 0.05 IMPLIES THAT ONE CAN REJECT  KLM2 320
C                  THE NULL HYPOTHESIS THAT THE SETS X AND Y ARE FROM   KLM2 330
C                  THE SAME DENSITY WITH 5 PER CENT PROBABILITY OF BEINGKLM2 340
C                  INCORRECT.  PROB = 1. - SMIRN(Z).                    KLM2 350
C                                                                       KLM2 360
C        REMARKS                                                        KLM2 370
C           N AND M SHOULD BE GREATER THAN OR EQUAL TO 100.  (SEE THE   KLM2 380
C           MATHEMATICAL DESCRIPTION FOR THIS SUBROUTINE AND FOR THE    KLM2 390
C           SUBROUTINE SMIRN, CONCERNING ASYMPTOTIC FORMULAE).          KLM2 400
C                                                                       KLM2 410
C           DOUBLE PRECISION USAGE---IT IS DOUBTFUL THAT THE USER WILL  KLM2 420
C           WISH TO PERFORM THIS TEST USING DOUBLE PRECISION ACCURACY.  KLM2 430
C           IF ONE WISHES TO COMMUNICATE WITH KOLM2 IN A DOUBLE         KLM2 440
C           PRECISION PROGRAM, HE SHOULD CALL THE FORTRAN SUPPLIED      KLM2 450
C           PROGRAM SNGL(X) PRIOR TO CALLING KOLM2, AND CALL THE        KLM2 460
C           FORTRAN SUPPLIED PROGRAM DBLE(X) AFTER EXITING FROM KOLM2.  KLM2 470
C           (NOTE THAT SUBROUTINE SMIRN DOES HAVE DOUBLE PRECISION      KLM2 480
C           CAPABILITY AS SUPPLIED BY THIS PACKAGE.)                    KLM2 490
C                                                                       KLM2 500
C                                                                       KLM2 510
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  KLM2 520
C           SMIRN                                                       KLM2 530
C                                                                       KLM2 540
C        METHOD                                                         KLM2 550
C           FOR REFERENCE, SEE (1) W. FELLER--ON THE KOLMOGOROV-SMIRNOV KLM2 560
C           LIMIT THEOREMS FOR EMPIRICAL DISTRIBUTIONS--                KLM2 570
C           ANNALS OF MATH. STAT., 19, 1948.  177-189,                  KLM2 580
C           (2) N. SMIRNOV--TABLE FOR ESTIMATING THE GOODNESS OF FIT    KLM2 590
C           OF EMPIRICAL DISTRIBUTIONS--ANNALS OF MATH. STAT., 19,      KLM2 600
C           1948.  279-281.                                             KLM2 610
C           (3) R. VON MISES--MATHEMATICAL THEORY OF PROBABILITY AND    KLM2 620
C           STATISTICS--ACADEMIC PRESS, NEW YORK, 1964.  490-493,       KLM2 630
C           (4) B.V. GNEDENKO--THE THEORY OF PROBABILITY--CHELSEA       KLM2 640
C           PUBLISHING COMPANY, NEW YORK, 1962.  384-401.               KLM2 650
C                                                                       KLM2 660
C     ..................................................................KLM2 670
C                                                                       KLM2 680
      SUBROUTINE KOLM2(X,Y,N,M,Z,PROB)                                  KLM2 690
      DIMENSION X(1),Y(1)                                               KLM2 700
C                                                                       KLM2 710
C        SORT X INTO ASCENDING SEQUENCE                                 KLM2 720
C                                                                       KLM2 730
      DO 5 I=2,N                                                        KLM2 740
      IF(X(I)-X(I-1))1,5,5                                              KLM2 750
    1 TEMP=X(I)                                                         KLM2 760
      IM=I-1                                                            KLM2 770
      DO 3 J=1,IM                                                       KLM2 780
      L=I-J                                                             KLM2 790
      IF(TEMP-X(L))2,4,4                                                KLM2 800
    2 X(L+1)=X(L)                                                       KLM2 810
    3 CONTINUE                                                          KLM2 820
      X(1)=TEMP                                                         KLM2 830
      GO TO 5                                                           KLM2 840
    4 X(L+1)=TEMP                                                       KLM2 850
    5 CONTINUE                                                          KLM2 860
C                                                                       KLM2 870
C        SORT Y INTO ASCENDING SEQUENCE                                 KLM2 880
C                                                                       KLM2 890
      DO 10 I=2,M                                                       KLM2 900
      IF(Y(I)-Y(I-1))6,10,10                                            KLM2 910
    6 TEMP=Y(I)                                                         KLM2 920
      IM=I-1                                                            KLM2 930
      DO 8  J=1,IM                                                      KLM2 940
      L=I-J                                                             KLM2 950
      IF(TEMP-Y(L))7,9,9                                                KLM2 960
    7 Y(L+1)=Y(L)                                                       KLM2 970
    8 CONTINUE                                                          KLM2 980
      Y(1)=TEMP                                                         KLM2 990
      GO TO 10                                                          KLM21000
    9 Y(L+1)=TEMP                                                       KLM21010
   10 CONTINUE                                                          KLM21020
C                                                                       KLM21030
C        CALCULATE D = ABS(FN-GM) OVER THE SPECTRUM OF X AND Y          KLM21040
C                                                                       KLM21050
      XN=FLOAT(N)                                                       KLM21060
      XN1=1./XN                                                         KLM21070
      XM=FLOAT(M)                                                       KLM21080
      XM1=1./XM                                                         KLM21090
      D=0.0                                                             KLM21100
      I=0                                                               KLM21110
      J=0                                                               KLM21120
      K=0                                                               KLM21130
      L=0                                                               KLM21140
   11 IF(X(I+1)-Y(J+1))12,13,18                                         KLM21150
   12 K=1                                                               KLM21160
      GO TO 14                                                          KLM21170
   13 K=0                                                               KLM21180
   14 I=I+1                                                             KLM21190
      IF(I-N)15,21,21                                                   KLM21200
   15 IF(X(I+1)-X(I))14,14,16                                           KLM21210
   16 IF(K)17,18,17                                                     KLM21220
C                                                                       KLM21230
C        CHOOSE THE MAXIMUM DIFFERENCE, D                               KLM21240
C                                                                       KLM21250
   17 D=AMAX1(D,ABS(FLOAT(I)*XN1-FLOAT(J)*XM1))                         KLM21260
      IF(L)22,11,22                                                     KLM21270
   18 J=J+1                                                             KLM21280
      IF(J-M)19,20,20                                                   KLM21290
   19 IF(Y(J+1)-Y(J))18,18,17                                           KLM21300
   20 L=1                                                               KLM21310
      GO TO 17                                                          KLM21320
   21 L=1                                                               KLM21330
      GO TO 16                                                          KLM21340
C                                                                       KLM21350
C        CALCULATE THE STATISTIC Z                                      KLM21360
C                                                                       KLM21370
   22 Z=D*SQRT((XN*XM)/(XN+XM))                                         KLM21380
C                                                                       KLM21390
C        CALCULATE THE PROBABILITY ASSOCIATED WITH Z                    KLM21400
C                                                                       KLM21410
      CALL SMIRN(Z,PROB)                                                KLM21420
      PROB=1.0-PROB                                                     KLM21430
      RETURN                                                            KLM21440
      END                                                               KLM21450
