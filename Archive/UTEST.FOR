C                                                                       UTES  10
C     ..................................................................UTES  20
C                                                                       UTES  30
C        SUBROUTINE UTEST                                               UTES  40
C                                                                       UTES  50
C        PURPOSE                                                        UTES  60
C           TEST WHETHER TWO INDEPENDENT GROUPS ARE FROM THE SAME       UTES  70
C           POPULATION BY MEANS OF MANN-WHITNEY U-TEST                  UTES  80
C                                                                       UTES  90
C        USAGE                                                          UTES 100
C           CALL UTEST(A,R,N1,N2,U,Z,IER)                               UTES 110
C                                                                       UTES 120
C        DESCRIPTION OF PARAMETERS                                      UTES 130
C           A  - INPUT VECTOR OF CASES CONSISTING OF TWO INDEPENDENT    UTES 140
C                GROUPS . SMALLER GROUP PRECEDES LARGER GROUP. LENGTH   UTES 150
C                IS N1+N2.                                              UTES 160
C           R  - OUTPUT VECTOR OF RANKS. SMALLEST VALUE IS RANKED 1,    UTES 170
C                LARGEST IS RANKED N. TIES ARE ASSIGNED AVERAGE OF TIED UTES 180
C                RANKS. LENGTH IS N1+N2.                                UTES 190
C           N1 - NUMBER OF CASES IN SMALLER GROUP                       UTES 200
C           N2 - NUMBER OF CASES IN LARGER GROUP                        UTES 210
C           U  - STATISTIC USED TO TEST HOMOGENEITY OF THE TWO          UTES 220
C                GROUPS (OUTPUT)                                        UTES 230
C           Z  - MEASURE OF SIGNIFICANCE OF U IN TERMS OF NORMAL        UTES 240
C                DISTRIBUTION (OUTPUT)                                  UTES 250
C           IER- 0, IF NO ERROR.                                        UTES 251
C              - 1, IF ALL VALUES OF ONE GROUP ARE TIED.                UTES 252
C                                                                       UTES 260
C        REMARKS                                                        UTES 270
C           Z IS SET TO ZERO IF N2 IS LESS THAN 20                      UTES 280
C                                                                       UTES 290
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  UTES 300
C           RANK                                                        UTES 310
C           TIE                                                         UTES 320
C                                                                       UTES 330
C        METHOD                                                         UTES 340
C           DESCRIBED IN S. SIEGEL, 'NONPARAMETRIC STATISTICS FOR THE   UTES 350
C           BEHAVIORAL SCIENCES', MCGRAW-HILL, NEW YORK, 1956,          UTES 360
C           CHAPTER 6                                                   UTES 370
C                                                                       UTES 380
C     ..................................................................UTES 390
C                                                                       UTES 400
      SUBROUTINE UTEST(A,R,N1,N2,U,Z)                                   UTES 410
      DIMENSION A(1),R(1)                                               UTES 420
C                                                                       UTES 430
C        RANK SCORES FROM BOTH GROUP TOGETHER IN ASCENDING ORDER, AND   UTES 440
C        ASSIGN TIED OBSERVATIONS AVERAGE OF TIED RANKS                 UTES 450
C                                                                       UTES 460
      N=N1+N2                                                           UTES 470
      CALL RANK(A,R,N)                                                  UTES 480
      Z=0.0                                                             UTES 490
C                                                                       UTES 500
C        SUM RANKS IN LARGER GROUP                                      UTES 510
C                                                                       UTES 520
      R2=0.0                                                            UTES 530
      NP=N1+1                                                           UTES 540
      DO 10 I=NP,N                                                      UTES 550
   10 R2=R2+R(I)                                                        UTES 560
C                                                                       UTES 570
C        CALCULATE U                                                    UTES 580
C                                                                       UTES 590
      FNX=N1*N2                                                         UTES 600
      FN=N                                                              UTES 610
      FN2=N2                                                            UTES 620
      UP=FNX+FN2*((FN2+1.0)/2.0)-R2                                     UTES 630
      U=FNX-UP                                                          UTES 640
      IF(UP-U) 20,30,30                                                 UTES 650
   20 U=UP                                                              UTES 660
C                                                                       UTES 670
C        TEST FOR N2 LESS THAN 20                                       UTES 680
C                                                                       UTES 690
   30 IF(N2-20) 80,40,40                                                UTES 700
C                                                                       UTES 710
C        COMPUTE STANDARD DEVIATION                                     UTES 720
C                                                                       UTES 730
   40 KT=1                                                              UTES 740
      CALL TIE(R,N,KT,TS)                                               UTES 750
      IF(TS) 50,60,50                                                   UTES 760
   50 IF (TS-(FN*FN*FN-FN)/12)52,51,52                                  UTES 761
   51 IER=1                                                             UTES 762
      GO TO 80                                                          UTES 763
   52 S=SQRT((FNX/(FN*(FN-1.0)))*(((FN*FN*FN-FN)/12.0)-TS))             UTES 770
      GO TO 70                                                          UTES 780
   60 S=SQRT(FNX*(FN+1.0)/12.0)                                         UTES 790
C                                                                       UTES 800
C        COMPUTE Z                                                      UTES 810
C                                                                       UTES 820
   70 Z=(U-FNX*0.5)/S                                                   UTES 830
   80 RETURN                                                            UTES 840
      END                                                               UTES 850
