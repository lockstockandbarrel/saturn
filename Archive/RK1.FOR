C                                                                       RK1   10
C     ..................................................................RK1   20
C                                                                       RK1   30
C        SUBROUTINE RK1                                                 RK1   40
C                                                                       RK1   50
C        PURPOSE                                                        RK1   60
C           INTEGRATES A FIRST ORDER DIFFERENTIAL EQUATION              RK1   70
C           DY/DX=FUN(X,Y) UP TO A SPECIFIED FINAL VALUE                RK1   80
C                                                                       RK1   90
C        USAGE                                                          RK1  100
C           CALL RK1(FUN,HI,XI,YI,XF,YF,ANSX,ANSY,IER)                  RK1  110
C                                                                       RK1  120
C        DESCRIPTION OF PARAMETERS                                      RK1  130
C           FUN -USER-SUPPLIED FUNCTION SUBPROGRAM WITH ARGUMENTS X,Y   RK1  140
C                WHICH GIVES DY/DX                                      RK1  150
C           HI  -THE STEP SIZE                                          RK1  160
C           XI  -INITIAL VALUE OF X                                     RK1  170
C           YI  -INITIAL VALUE OF Y WHERE YI=Y(XI)                      RK1  180
C           XF  -FINAL VALUE OF X                                       RK1  190
C           YF  -FINAL VALUE OF Y                                       RK1  200
C           ANSX-RESULTANT FINAL VALUE OF X                             RK1  210
C           ANSY-RESULTANT FINAL VALUE OF Y                             RK1  220
C                EITHER ANSX WILL EQUAL XF OR ANSY WILL EQUAL YF        RK1  230
C                DEPENDING ON WHICH IS REACHED FIRST                    RK1  240
C           IER -ERROR CODE                                             RK1  250
C                IER=0 NO ERROR                                         RK1  260
C                IER=1 STEP SIZE IS ZERO                                RK1  270
C                                                                       RK1  280
C        REMARKS                                                        RK1  290
C           IF XI IS GREATER THAN XF, ANSX=XI AND ANSY=YI               RK1  300
C           IF H IS ZERO, IER IS SET TO ONE, ANSX IS SET TO XI, AND     RK1  310
C           ANSY IS SET TO ZERO                                         RK1  320
C                                                                       RK1  330
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  RK1  340
C           FUN IS A TWO ARGUMENT FUNCTION SUBPROGRAM FURNISHED BY THE  RK1  350
C           USER.  DY/DX=FUN (X,Y)                                      RK1  360
C           CALLING PROGRAM MUST HAVE FORTRAN EXTERNAL STATEMENT        RK1  370
C           CONTAINING NAMES OF FUNCTION SUBPROGRAMS LISTED IN CALL TO  RK1  380
C           RK1                                                         RK1  390
C                                                                       RK1  400
C        METHOD                                                         RK1  410
C           USES FOURTH ORDER RUNGE-KUTTA INTEGRATION PROCESS ON A      RK1  420
C           RECURSIVE BASIS AS SHOWN IN F.B. HILDEBRAND, 'INTRODUCTION  RK1  430
C           TO NUMERICAL ANALYSIS',MCGRAW-HILL,1956. PROCESS IS         RK1  440
C           TERMINATED AND FINAL VALUE ADJUSTED WHEN EITHER XF OR YF    RK1  450
C           IS REACHED.                                                 RK1  460
C                                                                       RK1  470
C     ..................................................................RK1  480
C                                                                       RK1  490
      SUBROUTINE RK1(FUN,HI,XI,YI,XF,YF,ANSX,ANSY,IER)                  RK1  500
C                                                                       RK1  510
C        ...............................................................RK1  520
C                                                                       RK1  530
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  RK1  540
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      RK1  550
C        STATEMENT WHICH FOLLOWS.                                       RK1  560
C                                                                       RK1  570
C     DOUBLE PRECISION HI,XI,YI,XF,YF,ANSX,ANSY,H,XN,YN,HNEW,XN1,YN1,   RK1  580
C    1                 XX,YY,XNEW,YNEW,H2,T1,T2,T3,T4,FUN               RK1  590
C                                                                       RK1  600
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    RK1  610
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      RK1  620
C        ROUTINE.                                                       RK1  630
C                                                                       RK1  640
C        USER FUNCTION SUBPROGRAM, FUN, MUST BE IN DOUBLE PRECISION.    RK1  650
C                                                                       RK1  660
C        ...............................................................RK1  670
C                                                                       RK1  680
C     IF XF IS LESS THAN OR EQUAL TO XI, RETURN XI,YI AS ANSWER         RK1  690
C                                                                       RK1  700
      IER=0                                                             RK1  705
      IF(XF-XI) 11,11,12                                                RK1  710
   11 ANSX=XI                                                           RK1  720
      ANSY=YI                                                           RK1  730
      RETURN                                                            RK1  740
C                                                                       RK1  750
C     TEST INTERVAL VALUE                                               RK1  760
C                                                                       RK1  770
   12 H=HI                                                              RK1  780
      IF(HI) 16,14,20                                                   RK1  790
   14 IER=1                                                             RK1  800
      ANSX=XI                                                           RK1  810
      ANSY=0.0                                                          RK1  820
      RETURN                                                            RK1  830
   16 H=-HI                                                             RK1  840
C                                                                       RK1  850
C     SET XN=INITIAL X,YN=INITIAL Y                                     RK1  860
C                                                                       RK1  870
   20 XN=XI                                                             RK1  880
      YN=YI                                                             RK1  890
C                                                                       RK1  900
C     INTEGRATE ONE TIME STEP                                           RK1  910
C                                                                       RK1  920
      HNEW=H                                                            RK1  930
      JUMP=1                                                            RK1  940
      GO TO 170                                                         RK1  950
   25 XN1=XX                                                            RK1  960
      YN1=YY                                                            RK1  970
C                                                                       RK1  980
C     COMPARE XN1 (=X(N+1)) TO X FINAL AND BRANCH ACCORDINGLY           RK1  990
C                                                                       RK1 1000
      IF(XN1-XF)50,30,40                                                RK1 1010
C                                                                       RK1 1020
C     XN1=XF, RETURN (XF,YN1) AS ANSWER                                 RK1 1030
C                                                                       RK1 1040
   30 ANSX=XF                                                           RK1 1050
      ANSY=YN1                                                          RK1 1060
      GO TO 160                                                         RK1 1070
C                                                                       RK1 1080
C     XN1 GREATER THAN XF, SET NEW STEP SIZE AND INTEGRATE ONE STEP     RK1 1090
C     RETURN RESULTS OF INTEGRATION AS ANSWER                           RK1 1100
C                                                                       RK1 1110
   40 HNEW=XF-XN                                                        RK1 1120
      JUMP=2                                                            RK1 1130
      GO TO 170                                                         RK1 1140
   45 ANSX=XX                                                           RK1 1150
      ANSY=YY                                                           RK1 1160
      GO TO 160                                                         RK1 1170
C                                                                       RK1 1180
C     XN1 LESS THAN X FINAL, CHECK IF (YN,YN1) SPAN Y FINAL             RK1 1190
C                                                                       RK1 1200
C                                                                       RK1 1210
   50 IF((YN1-YF)*(YF-YN))60,70,110                                     RK1 1220
C                                                                       RK1 1230
C     YN1 AND YN DO NOT SPAN YF. SET (XN,YN) AS (XN1,YN1) AND REPEAT    RK1 1240
C                                                                       RK1 1250
   60 YN=YN1                                                            RK1 1260
      XN=XN1                                                            RK1 1270
      GO TO 170                                                         RK1 1280
C                                                                       RK1 1290
C     EITHER YN OR YN1 =YF. CHECK WHICH AND SET PROPER (X,Y) AS ANSWER  RK1 1300
C                                                                       RK1 1310
   70 IF(YN1-YF)80,100,80                                               RK1 1320
   80 ANSY=YN                                                           RK1 1330
      ANSX=XN                                                           RK1 1340
      GO TO 160                                                         RK1 1350
  100 ANSY=YN1                                                          RK1 1360
      ANSX=XN1                                                          RK1 1370
      GO TO 160                                                         RK1 1380
C                                                                       RK1 1390
C     YN AND YN1 SPAN YF. TRY TO FIND X VALUE ASSOCIATED WITH YF        RK1 1400
C                                                                       RK1 1410
  110 DO 140 I=1,10                                                     RK1 1420
C                                                                       RK1 1430
C     INTERPOLATE TO FIND NEW TIME STEP AND INTEGRATE ONE STEP          RK1 1440
C     TRY TEN INTERPOLATIONS AT MOST                                    RK1 1450
C                                                                       RK1 1460
      HNEW=((YF-YN )/(YN1-YN))*(XN1-XN)                                 RK1 1470
      JUMP=3                                                            RK1 1480
      GO TO 170                                                         RK1 1490
  115 XNEW=XX                                                           RK1 1500
      YNEW=YY                                                           RK1 1510
C                                                                       RK1 1520
C     COMPARE COMPUTED Y VALUE WITH YF AND BRANCH                       RK1 1530
C                                                                       RK1 1540
      IF(YNEW-YF)120,150,130                                            RK1 1550
C                                                                       RK1 1560
C     ADVANCE, YF IS BETWEEN YNEW AND YN1                               RK1 1570
C                                                                       RK1 1580
  120 YN=YNEW                                                           RK1 1590
      XN=XNEW                                                           RK1 1600
      GO TO 140                                                         RK1 1610
C                                                                       RK1 1620
C     ADVANCE, YF IS BETWEEN YN AND YNEW                                RK1 1630
C                                                                       RK1 1640
  130 YN1=YNEW                                                          RK1 1650
      XN1=XNEW                                                          RK1 1660
  140 CONTINUE                                                          RK1 1670
C                                                                       RK1 1680
C     RETURN (XNEW,YF) AS ANSWER                                        RK1 1690
C                                                                       RK1 1700
  150 ANSX=XNEW                                                         RK1 1710
      ANSY=YF                                                           RK1 1720
  160 RETURN                                                            RK1 1730
C                                                                       RK1 1740
  170 H2=HNEW/2.0                                                       RK1 1750
      T1=HNEW*FUN(XN,YN)                                                RK1 1760
      T2=HNEW*FUN(XN+H2,YN+T1/2.0)                                      RK1 1770
      T3=HNEW*FUN(XN+H2,YN+T2/2.0)                                      RK1 1780
      T4=HNEW*FUN(XN+HNEW,YN+T3)                                        RK1 1790
      YY=YN+(T1+2.0*T2+2.0*T3+T4)/6.0                                   RK1 1800
      XX=XN+HNEW                                                        RK1 1810
      GO TO (25,45,115), JUMP                                           RK1 1820
C                                                                       RK1 1830
      END                                                               RK1 1840
