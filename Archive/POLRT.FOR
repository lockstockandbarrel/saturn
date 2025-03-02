C                                                                       PLRT  10
C     ..................................................................PLRT  20
C                                                                       PLRT  30
C        SUBROUTINE POLRT                                               PLRT  40
C                                                                       PLRT  50
C        PURPOSE                                                        PLRT  60
C           COMPUTES THE REAL AND COMPLEX ROOTS OF A REAL POLYNOMIAL    PLRT  70
C                                                                       PLRT  80
C        USAGE                                                          PLRT  90
C           CALL POLRT(XCOF,COF,M,ROOTR,ROOTI,IER)                      PLRT 100
C                                                                       PLRT 110
C        DESCRIPTION OF PARAMETERS                                      PLRT 120
C           XCOF -VECTOR OF M+1 COEFFICIENTS OF THE POLYNOMIAL          PLRT 130
C                 ORDERED FROM SMALLEST TO LARGEST POWER                PLRT 140
C           COF  -WORKING VECTOR OF LENGTH M+1                          PLRT 150
C           M    -ORDER OF POLYNOMIAL                                   PLRT 160
C           ROOTR-RESULTANT VECTOR OF LENGTH M CONTAINING REAL ROOTS    PLRT 170
C                 OF THE POLYNOMIAL                                     PLRT 180
C           ROOTI-RESULTANT VECTOR OF LENGTH M CONTAINING THE           PLRT 190
C                 CORRESPONDING IMAGINARY ROOTS OF THE POLYNOMIAL       PLRT 200
C           IER  -ERROR CODE WHERE                                      PLRT 210
C                 IER=0  NO ERROR                                       PLRT 220
C                 IER=1  M LESS THAN ONE                                PLRT 230
C                 IER=2  M GREATER THAN 36                              PLRT 240
C                 IER=3  UNABLE TO DETERMINE ROOT WITH 500 INTERATIONS  PLRT 250
C                        ON 5 STARTING VALUES                           PLRT 260
C                 IER=4  HIGH ORDER COEFFICIENT IS ZERO                 PLRT 270
C                                                                       PLRT 280
C        REMARKS                                                        PLRT 290
C           LIMITED TO 36TH ORDER POLYNOMIAL OR LESS.                   PLRT 300
C           FLOATING POINT OVERFLOW MAY OCCUR FOR HIGH ORDER            PLRT 310
C           POLYNOMIALS BUT WILL NOT AFFECT THE ACCURACY OF THE RESULTS.PLRT 320
C                                                                       PLRT 330
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  PLRT 340
C           NONE                                                        PLRT 350
C                                                                       PLRT 360
C        METHOD                                                         PLRT 370
C           NEWTON-RAPHSON ITERATIVE TECHNIQUE.  THE FINAL ITERATIONS   PLRT 380
C           ON EACH ROOT ARE PERFORMED USING THE ORIGINAL POLYNOMIAL    PLRT 390
C           RATHER THAN THE REDUCED POLYNOMIAL TO AVOID ACCUMULATED     PLRT 400
C           ERRORS IN THE REDUCED POLYNOMIAL.                           PLRT 410
C                                                                       PLRT 420
C     ..................................................................PLRT 430
C                                                                       PLRT 440
      SUBROUTINE POLRT(XCOF,COF,M,ROOTR,ROOTI,IER)                      PLRT 450
      DIMENSION XCOF(1),COF(1),ROOTR(1),ROOTI(1)                        PLRT 460
      DOUBLE PRECISION XO,YO,X,Y,XPR,YPR,UX,UY,V,YT,XT,U,XT2,YT2,SUMSQ, PLRT 470
     1 DX,DY,TEMP,ALPHA                                                 PLRT 480
C                                                                       PLRT 490
C        ...............................................................PLRT 500
C                                                                       PLRT 510
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  PLRT 520
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      PLRT 530
C        STATEMENT WHICH FOLLOWS.                                       PLRT 540
C                                                                       PLRT 550
C     DOUBLE PRECISION XCOF,COF,ROOTR,ROOTI                             PLRT 560
C                                                                       PLRT 570
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    PLRT 580
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      PLRT 590
C        ROUTINE.                                                       PLRT 600
C        THE DOUBLE PRECISION VERSION MAY BE MODIFIED BY CHANGING THE   PLRT 610
C        CONSTANT IN STATEMENT 78 TO 1.0D-12 AND IN STATEMENT 122 TO    PLRT 620
C        1.0D-10.  THIS WILL PROVIDE HIGHER PRECISION RESULTS AT THE    PLRT 630
C        COST OF EXECUTION TIME                                         PLRT 640
C                                                                       PLRT 650
C        ...............................................................PLRT 660
C                                                                       PLRT 670
      IFIT=0                                                            PLRT 680
      N=M                                                               PLRT 690
      IER=0                                                             PLRT 700
      IF(XCOF(N+1))10,25,10                                             PLRT 710
   10 IF(N) 15,15,32                                                    PLRT 720
C                                                                       PLRT 730
C        SET ERROR CODE TO 1                                            PLRT 740
C                                                                       PLRT 750
   15 IER=1                                                             PLRT 760
   20 RETURN                                                            PLRT 770
C                                                                       PLRT 780
C        SET ERROR CODE TO 4                                            PLRT 790
C                                                                       PLRT 800
   25 IER=4                                                             PLRT 810
      GO TO 20                                                          PLRT 820
C                                                                       PLRT 830
C        SET ERROR CODE TO 2                                            PLRT 840
C                                                                       PLRT 850
   30 IER=2                                                             PLRT 860
      GO TO 20                                                          PLRT 870
   32 IF(N-36) 35,35,30                                                 PLRT 880
   35 NX=N                                                              PLRT 890
      NXX=N+1                                                           PLRT 900
      N2=1                                                              PLRT 910
      KJ1 = N+1                                                         PLRT 920
      DO 40 L=1,KJ1                                                     PLRT 930
      MT=KJ1-L+1                                                        PLRT 940
   40 COF(MT)=XCOF(L)                                                   PLRT 950
C                                                                       PLRT 960
C        SET INITIAL VALUES                                             PLRT 970
C                                                                       PLRT 980
   45 XO=.00500101                                                      PLRT 990
      YO=0.01000101                                                     PLRT1000
C                                                                       PLRT1010
C        ZERO INITIAL VALUE COUNTER                                     PLRT1020
C                                                                       PLRT1030
      IN=0                                                              PLRT1040
   50 X=XO                                                              PLRT1050
C                                                                       PLRT1060
C        INCREMENT INITIAL VALUES AND COUNTER                           PLRT1070
C                                                                       PLRT1080
      XO=-10.0*YO                                                       PLRT1090
      YO=-10.0*X                                                        PLRT1100
C                                                                       PLRT1110
C        SET X AND Y TO CURRENT VALUE                                   PLRT1120
C                                                                       PLRT1130
      X=XO                                                              PLRT1140
      Y=YO                                                              PLRT1150
      IN=IN+1                                                           PLRT1160
      GO TO 59                                                          PLRT1170
   55 IFIT=1                                                            PLRT1180
      XPR=X                                                             PLRT1190
      YPR=Y                                                             PLRT1200
C                                                                       PLRT1210
C        EVALUATE POLYNOMIAL AND DERIVATIVES                            PLRT1220
C                                                                       PLRT1230
   59 ICT=0                                                             PLRT1240
   60 UX=0.0                                                            PLRT1250
      UY=0.0                                                            PLRT1260
      V =0.0                                                            PLRT1270
      YT=0.0                                                            PLRT1280
      XT=1.0                                                            PLRT1290
      U=COF(N+1)                                                        PLRT1300
      IF(U) 65,130,65                                                   PLRT1310
   65 DO 70 I=1,N                                                       PLRT1320
      L =N-I+1                                                          PLRT1330
      TEMP=COF(L)                                                       PLRT1340
      XT2=X*XT-Y*YT                                                     PLRT1350
      YT2=X*YT+Y*XT                                                     PLRT1360
      U=U+TEMP*XT2                                                      PLRT1370
      V=V+TEMP*YT2                                                      PLRT1380
      FI=I                                                              PLRT1390
      UX=UX+FI*XT*TEMP                                                  PLRT1400
      UY=UY-FI*YT*TEMP                                                  PLRT1410
      XT=XT2                                                            PLRT1420
   70 YT=YT2                                                            PLRT1430
      SUMSQ=UX*UX+UY*UY                                                 PLRT1440
      IF(SUMSQ) 75,110,75                                               PLRT1450
   75 DX=(V*UY-U*UX)/SUMSQ                                              PLRT1460
      X=X+DX                                                            PLRT1470
      DY=-(U*UY+V*UX)/SUMSQ                                             PLRT1480
      Y=Y+DY                                                            PLRT1490
   78 IF(DABS(DY)+DABS(DX)-1.0D-05) 100,80,80                           PLRT1500
C                                                                       PLRT1510
C        STEP ITERATION COUNTER                                         PLRT1520
C                                                                       PLRT1530
   80 ICT=ICT+1                                                         PLRT1540
      IF(ICT-500) 60,85,85                                              PLRT1550
   85 IF(IFIT)100,90,100                                                PLRT1560
   90 IF(IN-5) 50,95,95                                                 PLRT1570
C                                                                       PLRT1580
C        SET ERROR CODE TO 3                                            PLRT1590
C                                                                       PLRT1600
   95 IER=3                                                             PLRT1610
      GO TO 20                                                          PLRT1620
  100 DO 105 L=1,NXX                                                    PLRT1630
      MT=KJ1-L+1                                                        PLRT1640
      TEMP=XCOF(MT)                                                     PLRT1650
      XCOF(MT)=COF(L)                                                   PLRT1660
  105 COF(L)=TEMP                                                       PLRT1670
      ITEMP=N                                                           PLRT1680
      N=NX                                                              PLRT1690
      NX=ITEMP                                                          PLRT1700
      IF(IFIT) 120,55,120                                               PLRT1710
  110 IF(IFIT) 115,50,115                                               PLRT1720
  115 X=XPR                                                             PLRT1730
      Y=YPR                                                             PLRT1740
  120 IFIT=0                                                            PLRT1750
  122 IF(DABS(Y)-1.0D-4*DABS(X)) 135,125,125                            PLRT1760
  125 ALPHA=X+X                                                         PLRT1770
      SUMSQ=X*X+Y*Y                                                     PLRT1780
      N=N-2                                                             PLRT1790
      GO TO 140                                                         PLRT1800
  130 X=0.0                                                             PLRT1810
      NX=NX-1                                                           PLRT1820
      NXX=NXX-1                                                         PLRT1830
  135 Y=0.0                                                             PLRT1840
      SUMSQ=0.0                                                         PLRT1850
      ALPHA=X                                                           PLRT1860
      N=N-1                                                             PLRT1870
  140 COF(2)=COF(2)+ALPHA*COF(1)                                        PLRT1880
  145 DO 150 L=2,N                                                      PLRT1890
  150 COF(L+1)=COF(L+1)+ALPHA*COF(L)-SUMSQ*COF(L-1)                     PLRT1900
  155 ROOTI(N2)=Y                                                       PLRT1910
      ROOTR(N2)=X                                                       PLRT1920
      N2=N2+1                                                           PLRT1930
      IF(SUMSQ) 160,165,160                                             PLRT1940
  160 Y=-Y                                                              PLRT1950
      SUMSQ=0.0                                                         PLRT1960
      GO TO 155                                                         PLRT1970
  165 IF(N) 20,20,45                                                    PLRT1980
      END                                                               PLRT1990
