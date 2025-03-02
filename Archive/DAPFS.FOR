C                                                                       DAPF  10
C     ..................................................................DAPF  20
C                                                                       DAPF  30
C        SUBROUTINE DAPFS                                               DAPF  40
C                                                                       DAPF  50
C        PURPOSE                                                        DAPF  60
C           PERFORM SYMMETRIC FACTORIZATION OF THE MATRIX OF THE NORMAL DAPF  70
C           EQUATIONS FOLLOWED BY CALCULATION OF THE LEAST SQUARES FIT  DAPF  80
C           OPTIONALLY                                                  DAPF  90
C                                                                       DAPF 100
C        USAGE                                                          DAPF 110
C           CALL DAPFS(WORK,IP,IRES,IOP,EPS,ETA,IER)                    DAPF 120
C                                                                       DAPF 130
C        DESCRIPTION OF PARAMETERS                                      DAPF 140
C           WORK  - GIVEN SYMMETRIC COEFFICIENT MATRIX, STORED          DAPF 150
C                   COMPRESSED, I.E UPPER TRIANGULAR PART COLUMNWISE.   DAPF 160
C                   THE GIVEN RIGHT HAND SIDE OCCUPIES THE NEXT IP      DAPF 170
C                   LOCATIONS IN WORK. THE VERY LAST COMPONENT OF WORK  DAPF 180
C                   CONTAINS THE SQUARE SUM OF FUNCTION VALUES E0       DAPF 190
C                   THIS SCHEME OF STORAGE ALLOCATION IS PRODUCED E.G.  DAPF 200
C                   BY SUBROUTINE APLL.                                 DAPF 210
C                   THE GIVEN MATRIX IS FACTORED IN THE FORM            DAPF 220
C                   TRANSPOSE(T)*T AND THE GIVEN RIGHT HAND SIDE IS     DAPF 230
C                   DIVIDED BY TRANSPOSE(T).                            DAPF 240
C                   THE UPPER TRIANGULAR FACTOR T IS RETURNED IN WORK IFDAPF 250
C                   IOP EQUALS ZERO.                                    DAPF 260
C                   IN CASE OF NONZERO IOP THE CALCULATED SOLUTIONS ARE DAPF 270
C                   STORED IN THE COLUMNS OF TRIANGULAR ARRAY WORK OF   DAPF 280
C                   CORRESPONDING DIMENSION AND E0  IS REPLACED BY THE  DAPF 290
C                   SQUARE SUM OF THE ERRORS FOR FIT OF DIMENSION IRES. DAPF 300
C                   THE TOTAL DIMENSION OF WORK IS (IP+1)*(IP+2)/2      DAPF 310
C           IP    - NUMBER OF FUNDAMENTAL FUNCTIONS USED FOR LEAST      DAPF 320
C                   SQUARES FIT                                         DAPF 330
C           IRES  - DIMENSION OF CALCULATED LEAST SQUARES FIT.          DAPF 340
C                   LET N1, N2, DENOTE THE FOLLOWING NUMBERS            DAPF 350
C                   N1 = MAXIMAL DIMENSION FOR WHICH NO LOSS OF         DAPF 360
C                        SIGNIFICANCE WAS INDICATED DURING FACTORIZATIONDAPF 370
C                   N2 = SMALLEST DIMENSION FOR WHICH THE SQUARE SUM OF DAPF 380
C                        THE ERRORS DOES NOT EXCEED TEST=ABS(ETA*FSQ)   DAPF 390
C                   THEN IRES=MINO(IP,N1) IF IOP IS NONNEGATIVE         DAPF 400
C                   AND  IRES=MINO(IP,N1,N2) IF IOP IS NEGATIVE         DAPF 410
C           IOP   - INPUT PARAMETER FOR SELECTION OF OPERATION          DAPF 420
C                   IOP = 0 MEANS TRIANGULAR FACTORIZATION, DIVISION OF DAPF 430
C                           THE RIGHT HAND SIDE BY TRANSPOSE(T) AND     DAPF 440
C                           CALCULATION OF THE SQUARE SUM OF ERRORS IS  DAPF 450
C                           PERFORMED ONLY                              DAPF 460
C                   IOP = +1 OR -1 MEANS THE SOLUTION OF DIMENSION IRES DAPF 470
C                           IS CALCULATED ADDITIONALLY                  DAPF 480
C                   IOP = +2 OR -2 MEANS ALL SOLUTIONS FOR DIMENSION ONEDAPF 490
C                           UP TO IRES ARE CALCULATED ADDITIONALLY      DAPF 500
C           EPS   - RELATIVE TOLERANCE FOR TEST ON LOSS OF SIGNIFICANCE.DAPF 510
C                   A SENSIBLE VALUE IS BETWEEN 1.E-3 AND 1.E-6         DAPF 520
C           ETA   - RELATIVE TOLERANCE FOR TOLERATED SQUARE SUM OF      DAPF 530
C                   ERRORS. A REALISTIC VALUE IS BETWEEN 1.E0 AND 1.E-6 DAPF 540
C           IER   - RESULTANT ERROR PARAMETER                           DAPF 550
C                   IER =-1 MEANS NONPOSITIVE IP                        DAPF 560
C                   IER = 0 MEANS NO LOSS OF SIGNIFICANCE DETECTED      DAPF 570
C                           AND SPECIFIED TOLERANCE OF ERRORS REACHED   DAPF 580
C                   IER = 1 MEANS LOSS OF SIGNIFICANCE DETECTED OR      DAPF 590
C                           SPECIFIED TOLERANCE OF ERRORS NOT REACHED   DAPF 600
C        NOTE: WORK,SUM AND PIV ARE DOUBLE PRECISION                    DAPF 610
C        REMARKS                                                        DAPF 620
C           THE ABSOLUTE TOLERANCE USED INTERNALLY FOR TEST ON LOSS OF  DAPF 630
C           SIGNIFICANCE IS TOL=ABS(EPS*WORK(1)).                       DAPF 640
C           THE ABSOLUTE TOLERANCE USED INTERNALLY FOR THE SQUARE SUM OFDAPF 650
C           ERRORS IS ABS(ETA*FSQ).                                     DAPF 660
C           IOP GREATER THAN 2 HAS THE SAME EFFECT AS IOP = 2.          DAPF 670
C           IOP LESS THAN -2 HAS THE SAME EFFECT AS IOP =-2.            DAPF 680
C           IRES = 0 MEANS THE ABSOLUTE VALUE OF EPS IS NOT LESS THAN   DAPF 690
C           ONE AND/OR WORK(1) IS NOT POSITIVE AND/OR IP IS NOT POSITIVEDAPF 700
C                                                                       DAPF 710
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DAPF 720
C           NONE                                                        DAPF 730
C                                                                       DAPF 740
C        METHOD                                                         DAPF 750
C           CALCULATION OF THE LEAST SQUARES FITS IS DONE USING         DAPF 760
C           CHOLESKYS SQUARE ROOT METHOD FOR SYMMETRIC FACTORIZATION.   DAPF 770
C           THE INCORPORATED TEST ON LOSS OF SIGNIFICANCE MEANS EACH    DAPF 780
C           RADICAND MUST BE GREATER THAN THE INTERNAL ABSOLUTE         DAPF 790
C           TOLERANCE TOL=ABS(EPS*WORK(1)).                             DAPF 800
C           IN CASE OF LOSS OF SIGNIFICANCE IN THE ABOVE SENSE ONLY A   DAPF 810
C           SUBSYSTEM OF THE NORMAL EQUATIONS IS SOLVED.                DAPF 820
C           IN CASE OF NEGATIVE IOP THE TRIANGULAR FACTORIZATION IS     DAPF 830
C           TERMINATED PREMATURELY EITHER IF THE SQUARE SUM OF THE      DAPF 840
C           ERRORS DOES NOT EXCEED ETA*FSQ OR IF THERE IS INDICATION    DAPF 850
C           FOR LOSS OF SIGNIFICANCE                                    DAPF 860
C                                                                       DAPF 870
C     ..................................................................DAPF 880
C                                                                       DAPF 890
      SUBROUTINE DAPFS(WORK,IP,IRES,IOP,EPS,ETA,IER)                    DAPF 900
C                                                                       DAPF 910
         DOUBLE PRECISION WORK,SUM,PIV                                  DAPF 920
C        DIMENSIONED DUMMY VARIABLES                                    DAPF 930
      DIMENSION WORK(1)                                                 DAPF 940
      IRES=0                                                            DAPF 950
C                                                                       DAPF 960
C        TEST OF SPECIFIED DIMENSION                                    DAPF 970
      IF(IP)1,1,2                                                       DAPF 980
C                                                                       DAPF 990
C        ERROR RETURN IN CASE OF ILLEGAL DIMENSION                      DAPF1000
    1 IER=-1                                                            DAPF1010
      RETURN                                                            DAPF1020
C                                                                       DAPF1030
C        INITIALIZE FACTORIZATION PROCESS                               DAPF1040
    2 IPIV=0                                                            DAPF1050
      IPP1=IP+1                                                         DAPF1060
      IER=1                                                             DAPF1070
      ITE=IP*IPP1/2                                                     DAPF1080
      IEND=ITE+IPP1                                                     DAPF1090
      TOL=ABS(EPS*WORK(1))                                              DAPF1100
      TEST=ABS(ETA*WORK(IEND))                                          DAPF1110
C                                                                       DAPF1120
C        START LOOP OVER ALL ROWS OF WORK                               DAPF1130
      DO 11 I=1,IP                                                      DAPF1140
      IPIV=IPIV+I                                                       DAPF1150
      JA=IPIV-IRES                                                      DAPF1160
      JE=IPIV-1                                                         DAPF1170
C                                                                       DAPF1180
C        FORM SCALAR PRODUCT NEEDED TO MODIFY CURRENT ROW ELEMENTS      DAPF1190
      JK=IPIV                                                           DAPF1200
      DO 9 K=I,IPP1                                                     DAPF1210
      SUM=0.                                                            DAPF1220
      IF(IRES)5,5,3                                                     DAPF1230
    3 JK=JK-IRES                                                        DAPF1240
      DO 4 J=JA,JE                                                      DAPF1250
      SUM=SUM+WORK(J)*WORK(JK)                                          DAPF1260
    4 JK=JK+1                                                           DAPF1270
    5 IF(JK-IPIV)6,6,8                                                  DAPF1280
C                                                                       DAPF1290
C        TEST FOR LOSS OF SIGNIFICANCE                                  DAPF1300
    6 SUM=WORK(IPIV)-SUM                                                DAPF1310
      IF(SUM-TOL)12,12,7                                                DAPF1320
    7 SUM=SQRT(SUM)                                                     DAPF1330
      WORK(IPIV)=SUM                                                    DAPF1340
      PIV=1.0D0/SUM                                                     DAPF1350
      GOTO 9                                                            DAPF1360
C                                                                       DAPF1370
C        UPDATE OFF-DIAGONAL TERMS                                      DAPF1380
    8 SUM=(WORK(JK)-SUM)*PIV                                            DAPF1390
      WORK(JK)=SUM                                                      DAPF1400
    9 JK=JK+K                                                           DAPF1410
C                                                                       DAPF1420
C        UPDATE SQUARE SUM OF ERRORS                                    DAPF1430
      WORK(IEND)=WORK(IEND)-SUM*SUM                                     DAPF1440
C                                                                       DAPF1450
C        RECORD ADDRESS OF LAST PIVOT ELEMENT                           DAPF1460
      IRES=IRES+1                                                       DAPF1470
      IADR=IPIV                                                         DAPF1480
C                                                                       DAPF1490
C        TEST FOR TOLERABLE ERROR IF SPECIFIED                          DAPF1500
      IF(IOP)10,11,11                                                   DAPF1510
   10 IF(WORK(IEND)-TEST)13,13,11                                       DAPF1520
   11 CONTINUE                                                          DAPF1530
      IF(IOP)12,22,12                                                   DAPF1540
C                                                                       DAPF1550
C        PERFORM BACK SUBSTITUTION IF SPECIFIED                         DAPF1560
   12 IF(IOP)14,23,14                                                   DAPF1570
   13 IER=0                                                             DAPF1580
   14 IPIV=IRES                                                         DAPF1590
   15 IF(IPIV)23,23,16                                                  DAPF1600
   16 SUM=0.0D0                                                         DAPF1610
      JA=ITE+IPIV                                                       DAPF1620
      JJ=IADR                                                           DAPF1630
      JK=IADR                                                           DAPF1640
      K=IPIV                                                            DAPF1650
      DO 19 I=1,IPIV                                                    DAPF1660
      WORK(JK)=(WORK(JA)-SUM)/WORK(JJ)                                  DAPF1670
      IF(K-1)20,20,17                                                   DAPF1680
   17 JE=JJ-1                                                           DAPF1690
      SUM=0.0D0                                                         DAPF1700
      DO 18 J=K,IPIV                                                    DAPF1710
      SUM=SUM+WORK(JK)*WORK(JE)                                         DAPF1720
      JK=JK+1                                                           DAPF1730
   18 JE=JE+J                                                           DAPF1740
      JK=JE-IPIV                                                        DAPF1750
      JA=JA-1                                                           DAPF1760
      JJ=JJ-K                                                           DAPF1770
   19 K=K-1                                                             DAPF1780
   20 IF(IOP/2)21,23,21                                                 DAPF1790
   21 IADR=IADR-IPIV                                                    DAPF1800
      IPIV=IPIV-1                                                       DAPF1810
      GOTO 15                                                           DAPF1820
C                                                                       DAPF1830
C        NORMAL RETURN                                                  DAPF1840
   22 IER=0                                                             DAPF1850
   23 RETURN                                                            DAPF1860
      END                                                               DAPF1870
