C                                                                       APFS  10
C     ..................................................................APFS  20
C                                                                       APFS  30
C        SUBROUTINE APFS                                                APFS  40
C                                                                       APFS  50
C        PURPOSE                                                        APFS  60
C           PERFORM SYMMETRIC FACTORIZATION OF THE MATRIX OF THE NORMAL APFS  70
C           EQUATIONS FOLLOWED BY CALCULATION OF THE LEAST SQUARES FIT  APFS  80
C           OPTIONALLY                                                  APFS  90
C                                                                       APFS 100
C        USAGE                                                          APFS 110
C           CALL APFS(WORK,IP,IRES,IOP,EPS,ETA,IER)                     APFS 120
C                                                                       APFS 130
C        DESCRIPTION OF PARAMETERS                                      APFS 140
C           WORK  - GIVEN SYMMETRIC COEFFICIENT MATRIX, STORED          APFS 150
C                   COMPRESSED, I.E UPPER TRIANGULAR PART COLUMNWISE.   APFS 160
C                   THE GIVEN RIGHT HAND SIDE OCCUPIES THE NEXT IP      APFS 170
C                   LOCATIONS IN WORK. THE VERY LAST COMPONENT OF WORK  APFS 180
C                   CONTAINS THE SQUARE SUM OF FUNCTION VALUES E0       APFS 190
C                   THIS SCHEME OF STORAGE ALLOCATION IS PRODUCED E.G.  APFS 200
C                   BY SUBROUTINE APLL.                                 APFS 210
C                   THE GIVEN MATRIX IS FACTORED IN THE FORM            APFS 220
C                   TRANSPOSE(T)*T AND THE GIVEN RIGHT HAND SIDE IS     APFS 230
C                   DIVIDED BY TRANSPOSE(T).                            APFS 240
C                   THE UPPER TRIANGULAR FACTOR T IS RETURNED IN WORK IFAPFS 250
C                   IOP EQUALS ZERO.                                    APFS 260
C                   IN CASE OF NONZERO IOP THE CALCULATED SOLUTIONS ARE APFS 270
C                   STORED IN THE COLUMNS OF TRIANGULAR ARRAY WORK OF   APFS 280
C                   CORRESPONDING DIMENSION AND E0  IS REPLACED BY THE  APFS 290
C                   SQUARE SUM OF THE ERRORS FOR FIT OF DIMENSION IRES. APFS 300
C                   THE TOTAL DIMENSION OF WORK IS (IP+1)*(IP+2)/2      APFS 310
C           IP    - NUMBER OF FUNDAMENTAL FUNCTIONS USED FOR LEAST      APFS 320
C                   SQUARES FIT                                         APFS 330
C           IRES  - DIMENSION OF CALCULATED LEAST SQUARES FIT.          APFS 340
C                   LET N1, N2, DENOTE THE FOLLOWING NUMBERS            APFS 350
C                   N1 = MAXIMAL DIMENSION FOR WHICH NO LOSS OF         APFS 360
C                        SIGNIFICANCE WAS INDICATED DURING FACTORIZATIONAPFS 370
C                   N2 = SMALLEST DIMENSION FOR WHICH THE SQUARE SUM OF APFS 380
C                        THE ERRORS DOES NOT EXCEED TEST=ABS(ETA*FSQ)   APFS 390
C                   THEN IRES=MINO(IP,N1) IF IOP IS NONNEGATIVE         APFS 400
C                   AND  IRES=MINO(IP,N1,N2) IF IOP IS NEGATIVE         APFS 410
C           IOP   - INPUT PARAMETER FOR SELECTION OF OPERATION          APFS 420
C                   IOP = 0 MEANS TRIANGULAR FACTORIZATION, DIVISION OF APFS 430
C                           THE RIGHT HAND SIDE BY TRANSPOSE(T) AND     APFS 440
C                           CALCULATION OF THE SQUARE SUM OF ERRORS IS  APFS 450
C                           PERFORMED ONLY                              APFS 460
C                   IOP = +1 OR -1 MEANS THE SOLUTION OF DIMENSION IRES APFS 470
C                           IS CALCULATED ADDITIONALLY                  APFS 480
C                   IOP = +2 OR -2 MEANS ALL SOLUTIONS FOR DIMENSION ONEAPFS 490
C                           UP TO IRES ARE CALCULATED ADDITIONALLY      APFS 500
C           EPS   - RELATIVE TOLERANCE FOR TEST ON LOSS OF SIGNIFICANCE.APFS 510
C                   A SENSIBLE VALUE IS BETWEEN 1.E-3 AND 1.E-6         APFS 520
C           ETA   - RELATIVE TOLERANCE FOR TOLERATED SQUARE SUM OF      APFS 530
C                   ERRORS. A REALISTIC VALUE IS BETWEEN 1.E0 AND 1.E-6 APFS 540
C           IER   - RESULTANT ERROR PARAMETER                           APFS 550
C                   IER =-1 MEANS NONPOSITIVE IP                        APFS 560
C                   IER = 0 MEANS NO LOSS OF SIGNIFICANCE DETECTED      APFS 570
C                           AND SPECIFIED TOLERANCE OF ERRORS REACHED   APFS 580
C                   IER = 1 MEANS LOSS OF SIGNIFICANCE DETECTED OR      APFS 590
C                           SPECIFIED TOLERANCE OF ERRORS NOT REACHED   APFS 600
C                                                                       APFS 610
C        REMARKS                                                        APFS 620
C           THE ABSOLUTE TOLERANCE USED INTERNALLY FOR TEST ON LOSS OF  APFS 630
C           SIGNIFICANCE IS TOL=ABS(EPS*WORK(1)).                       APFS 640
C           THE ABSOLUTE TOLERANCE USED INTERNALLY FOR THE SQUARE SUM OFAPFS 650
C           ERRORS IS ABS(ETA*FSQ).                                     APFS 660
C           IOP GREATER THAN 2 HAS THE SAME EFFECT AS IOP = 2.          APFS 670
C           IOP LESS THAN -2 HAS THE SAME EFFECT AS IOP =-2.            APFS 680
C           IRES = 0 MEANS THE ABSOLUTE VALUE OF EPS IS NOT LESS THAN   APFS 690
C           ONE AND/OR WORK(1) IS NOT POSITIVE AND/OR IP IS NOT POSITIVEAPFS 700
C                                                                       APFS 710
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  APFS 720
C           NONE                                                        APFS 730
C                                                                       APFS 740
C        METHOD                                                         APFS 750
C           CALCULATION OF THE LEAST SQUARES FITS IS DONE USING         APFS 760
C           CHOLESKYS SQUARE ROOT METHOD FOR SYMMETRIC FACTORIZATION.   APFS 770
C           THE INCORPORATED TEST ON LOSS OF SIGNIFICANCE MEANS EACH    APFS 780
C           RADICAND MUST BE GREATER THAN THE INTERNAL ABSOLUTE         APFS 790
C           TOLERANCE TOL=ABS(EPS*WORK(1)).                             APFS 800
C           IN CASE OF LOSS OF SIGNIFICANCE IN THE ABOVE SENSE ONLY A   APFS 810
C           SUBSYSTEM OF THE NORMAL EQUATIONS IS SOLVED.                APFS 820
C           IN CASE OF NEGATIVE IOP THE TRIANGULAR FACTORIZATION IS     APFS 830
C           TERMINATED PREMATURELY EITHER IF THE SQUARE SUM OF THE      APFS 840
C           ERRORS DOES NOT EXCEED ETA*FSQ OR IF THERE IS INDICATION    APFS 850
C           FOR LOSS OF SIGNIFICANCE                                    APFS 860
C                                                                       APFS 870
C     ..................................................................APFS 880
C                                                                       APFS 890
      SUBROUTINE APFS(WORK,IP,IRES,IOP,EPS,ETA,IER)                     APFS 900
C                                                                       APFS 910
C                                                                       APFS 920
C        DIMENSIONED DUMMY VARIABLES                                    APFS 930
      DIMENSION WORK(1)                                                 APFS 940
      IRES=0                                                            APFS 950
C                                                                       APFS 960
C        TEST OF SPECIFIED DIMENSION                                    APFS 970
      IF(IP)1,1,2                                                       APFS 980
C                                                                       APFS 990
C        ERROR RETURN IN CASE OF ILLEGAL DIMENSION                      APFS1000
    1 IER=-1                                                            APFS1010
      RETURN                                                            APFS1020
C                                                                       APFS1030
C        INITIALIZE FACTORIZATION PROCESS                               APFS1040
    2 IPIV=0                                                            APFS1050
      IPP1=IP+1                                                         APFS1060
      IER=1                                                             APFS1070
      ITE=IP*IPP1/2                                                     APFS1080
      IEND=ITE+IPP1                                                     APFS1090
      TOL=ABS(EPS*WORK(1))                                              APFS1100
      TEST=ABS(ETA*WORK(IEND))                                          APFS1110
C                                                                       APFS1120
C        START LOOP OVER ALL ROWS OF WORK                               APFS1130
      DO 11 I=1,IP                                                      APFS1140
      IPIV=IPIV+I                                                       APFS1150
      JA=IPIV-IRES                                                      APFS1160
      JE=IPIV-1                                                         APFS1170
C                                                                       APFS1180
C        FORM SCALAR PRODUCT NEEDED TO MODIFY CURRENT ROW ELEMENTS      APFS1190
      JK=IPIV                                                           APFS1200
      DO 9 K=I,IPP1                                                     APFS1210
      SUM=0.                                                            APFS1220
      IF(IRES)5,5,3                                                     APFS1230
    3 JK=JK-IRES                                                        APFS1240
      DO 4 J=JA,JE                                                      APFS1250
      SUM=SUM+WORK(J)*WORK(JK)                                          APFS1260
    4 JK=JK+1                                                           APFS1270
    5 IF(JK-IPIV)6,6,8                                                  APFS1280
C                                                                       APFS1290
C        TEST FOR LOSS OF SIGNIFICANCE                                  APFS1300
    6 SUM=WORK(IPIV)-SUM                                                APFS1310
      IF(SUM-TOL)12,12,7                                                APFS1320
    7 SUM=SQRT(SUM)                                                     APFS1330
      WORK(IPIV)=SUM                                                    APFS1340
      PIV=1./SUM                                                        APFS1350
      GOTO 9                                                            APFS1360
C                                                                       APFS1370
C        UPDATE OFF-DIAGONAL TERMS                                      APFS1380
    8 SUM=(WORK(JK)-SUM)*PIV                                            APFS1390
      WORK(JK)=SUM                                                      APFS1400
    9 JK=JK+K                                                           APFS1410
C                                                                       APFS1420
C        UPDATE SQUARE SUM OF ERRORS                                    APFS1430
      WORK(IEND)=WORK(IEND)-SUM*SUM                                     APFS1440
C                                                                       APFS1450
C        RECORD ADDRESS OF LAST PIVOT ELEMENT                           APFS1460
      IRES=IRES+1                                                       APFS1470
      IADR=IPIV                                                         APFS1480
C                                                                       APFS1490
C        TEST FOR TOLERABLE ERROR IF SPECIFIED                          APFS1500
      IF(IOP)10,11,11                                                   APFS1510
   10 IF(WORK(IEND)-TEST)13,13,11                                       APFS1520
   11 CONTINUE                                                          APFS1530
      IF(IOP)12,22,12                                                   APFS1540
C                                                                       APFS1550
C        PERFORM BACK SUBSTITUTION IF SPECIFIED                         APFS1560
   12 IF(IOP)14,23,14                                                   APFS1570
   13 IER=0                                                             APFS1580
   14 IPIV=IRES                                                         APFS1590
   15 IF(IPIV)23,23,16                                                  APFS1600
   16 SUM=0.                                                            APFS1610
      JA=ITE+IPIV                                                       APFS1620
      JJ=IADR                                                           APFS1630
      JK=IADR                                                           APFS1640
      K=IPIV                                                            APFS1650
      DO 19 I=1,IPIV                                                    APFS1660
      WORK(JK)=(WORK(JA)-SUM)/WORK(JJ)                                  APFS1670
      IF(K-1)20,20,17                                                   APFS1680
   17 JE=JJ-1                                                           APFS1690
      SUM=0.                                                            APFS1700
      DO 18 J=K,IPIV                                                    APFS1710
      SUM=SUM+WORK(JK)*WORK(JE)                                         APFS1720
      JK=JK+1                                                           APFS1730
   18 JE=JE+J                                                           APFS1740
      JK=JE-IPIV                                                        APFS1750
      JA=JA-1                                                           APFS1760
      JJ=JJ-K                                                           APFS1770
   19 K=K-1                                                             APFS1780
   20 IF(IOP/2)21,23,21                                                 APFS1790
   21 IADR=IADR-IPIV                                                    APFS1800
      IPIV=IPIV-1                                                       APFS1810
      GOTO 15                                                           APFS1820
C                                                                       APFS1830
C        NORMAL RETURN                                                  APFS1840
   22 IER=0                                                             APFS1850
   23 RETURN                                                            APFS1860
      END                                                               APFS1870
