C                                                                       APMM  10
C     ..................................................................APMM  20
C                                                                       APMM  30
C        SUBROUTINE APMM                                                APMM  40
C                                                                       APMM  50
C        PURPOSE                                                        APMM  60
C           APPROXIMATE A FUNCTION TABULATED IN N POINTS BY ANY LINEAR  APMM  70
C           COMBINATION OF M GIVEN CONTINUOUS FUNCTIONS IN THE SENSE    APMM  80
C           OF CHEBYSHEV.                                               APMM  90
C                                                                       APMM 100
C        USAGE                                                          APMM 110
C           CALL APMM(FCT,N,M,TOP,IHE,PIV,T,ITER,IER)                   APMM 120
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT IN THE         APMM 130
C           CALLING PROGRAM.                                            APMM 140
C                                                                       APMM 150
C        DESCRIPTION OF PARAMETERS                                      APMM 160
C           FCT    - NAME OF SUBROUTINE TO BE SUPPLIED BY THE USER.     APMM 170
C                    IT COMPUTES VALUES OF M GIVEN FUNCTIONS FOR        APMM 180
C                    ARGUMENT VALUE X.                                  APMM 190
C                    USAGE                                              APMM 200
C                       CALL FCT(Y,X,K)                                 APMM 210
C                    DESCRIPTION OF PARAMETERS                          APMM 220
C                       Y   - RESULT VECTOR OF DIMENSION M CONTAINING   APMM 230
C                             THE VALUES OF GIVEN CONTINUOUS FUNCTIONS  APMM 240
C                             FOR GIVEN ARGUMENT X                      APMM 250
C                       X   - ARGUMENT VALUE                            APMM 260
C                       K   - AN INTEGER VALUE WHICH IS EQUAL TO M-1    APMM 270
C                    REMARKS                                            APMM 280
C                       IF APPROXIMATION BY NORMAL CHEBYSHEV, SHIFTED   APMM 290
C                       CHEBYSHEV, LEGENDRE, LAGUERRE, HERMITE POLYNO-  APMM 300
C                       MIALS IS DESIRED SUBROUTINES CNP, CSP, LEP,     APMM 310
C                       LAP, HEP, RESPECTIVELY FROM SSP COULD BE USED.  APMM 320
C           N      - NUMBER OF DATA POINTS DEFINING THE FUNCTION WHICH  APMM 330
C                    IS TO BE APPROXIMATED                              APMM 340
C           M      - NUMBER OF GIVEN CONTINUOUS FUNCTIONS FROM WHICH    APMM 350
C                    THE APPROXIMATING FUNCTION IS CONSTRUCTED.         APMM 360
C           TOP    - VECTOR OF DIMENSION 3*N.                           APMM 370
C                    ON ENTRY IT MUST CONTAIN FROM TOP(1) UP TO TOP(N)  APMM 380
C                    THE GIVEN N FUNCTION VALUES AND FROM TOP(N+1) UP   APMM 390
C                    TO TOP(2*N) THE CORRESPONDING NODES                APMM 400
C                    ON RETURN TOP CONTAINS FROM TOP(1) UP TO TOP(N)    APMM 410
C                    THE ERRORS AT THOSE N NODES.                       APMM 420
C                    OTHER VALUES OF TOP ARE SCRATCH.                   APMM 430
C           IHE    - INTEGER VECTOR OF DIMENSION 3*M+4*N+6              APMM 440
C           PIV    - VECTOR OF DIMENSION 3*M+6.                         APMM 450
C                    ON RETURN PIV CONTAINS AT PIV(1) UP TO PIV(M) THE  APMM 460
C                    RESULTING COEFFICIENTS OF LINEAR APPROXIMATION.    APMM 470
C           T      - AUXILIARY VECTOR OF DIMENSION (M+2)*(M+2)          APMM 480
C           ITER   - RESULTANT INTEGER WHICH SPECIFIES THE NUMBER OF    APMM 490
C                    ITERATIONS NEEDED                                  APMM 500
C           IER    - RESULTANT ERROR PARAMETER CODED IN THE FOLLOWING   APMM 510
C                    FORM                                               APMM 520
C                     IER=0  - NO ERROR                                 APMM 530
C                     IER=1  - THE NUMBER OF ITERATIONS HAS REACHED     APMM 540
C                              THE INTERNAL MAXIMUM N+M                 APMM 550
C                     IER=-1 - NO RESULT BECAUSE OF WRONG INPUT PARA-   APMM 560
C                              METER M OR N OR SINCE AT SOME ITERATION  APMM 570
C                              NO SUITABLE PIVOT COULD BE FOUND         APMM 580
C                                                                       APMM 590
C        REMARKS                                                        APMM 600
C           NO ACTION BESIDES ERROR MESSAGE IN CASE M LESS THAN 1 OR    APMM 610
C           N LESS THAN 2.                                              APMM 620
C                                                                       APMM 630
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  APMM 640
C           THE EXTERNAL SUBROUTINE FCT MUST BE FURNISHED BY THE USER.  APMM 650
C                                                                       APMM 660
C        METHOD                                                         APMM 670
C           THE PROBLEM OF APPROXIMATION A TABULATED FUNCTION BY ANY    APMM 680
C           LINEAR COMBINATION OF GIVEN FUNCTIONS IN THE SENSE OF       APMM 690
C           CHEBYSHEV (I.E. TO MINIMIZE THE MAXIMUM ERROR) IS TRANS-    APMM 700
C           FORMED INTO A LINEAR PROGRAMMING PROBLEM. APMM USES A       APMM 710
C           REVISED SIMPLEX METHOD TO SOLVE A CORRESPONDING DUAL        APMM 720
C           PROBLEM. FOR REFERENCE, SEE                                 APMM 730
C           I.BARRODALE/A.YOUNG, ALGORITHMS FOR BEST L-SUB-ONE AND      APMM 740
C           L-SUB-INFINITY, LINEAR APPROXIMATIONS ON A DISCRETE SET,    APMM 750
C           NUMERISCHE MATHEMATIK, VOL.8, ISS.3 (1966), PP.295-306.     APMM 760
C                                                                       APMM 770
C     ..................................................................APMM 780
C                                                                       APMM 790
      SUBROUTINE APMM(FCT,N,M,TOP,IHE,PIV,T,ITER,IER)                   APMM 800
C                                                                       APMM 810
C                                                                       APMM 820
      DIMENSION TOP(1),IHE(1),PIV(1),T(1)                               APMM 830
      DOUBLE PRECISION DSUM                                             APMM 840
C                                                                       APMM 850
C        TEST ON WRONG INPUT PARAMETERS N AND M                         APMM 860
      IER=-1                                                            APMM 870
      IF (N-1) 81,81,1                                                  APMM 880
    1 IF(M) 81,81,2                                                     APMM 890
C                                                                       APMM 900
C        INITIALIZE CHARACTERISTIC VECTORS FOR THE TABLEAU              APMM 910
    2 IER=0                                                             APMM 920
C                                                                       APMM 930
C        PREPARE TOP-ROW TOP                                            APMM 940
      DO 3 I=1,N                                                        APMM 950
      K=I+N                                                             APMM 960
      J=K+N                                                             APMM 970
      TOP(J)=TOP(K)                                                     APMM 980
    3 TOP(K)=-TOP(I)                                                    APMM 990
C                                                                       APMM1000
C        PREPARE INVERSE TRANSFORMATION MATRIX T                        APMM1010
      L=M+2                                                             APMM1020
      LL=L*L                                                            APMM1030
      DO 4 I=1,LL                                                       APMM1040
    4 T(I)=0.                                                           APMM1050
      K=1                                                               APMM1060
      J=L+1                                                             APMM1070
      DO 5 I=1,L                                                        APMM1080
      T(K)=1.                                                           APMM1090
    5 K=K+J                                                             APMM1100
C                                                                       APMM1110
C        PREPARE INDEX-VECTOR IHE                                       APMM1120
      DO 6 I=1,L                                                        APMM1130
      K=I+L                                                             APMM1140
      J=K+L                                                             APMM1150
      IHE(I)=0                                                          APMM1160
      IHE(K)=I                                                          APMM1170
    6 IHE(J)=1-I                                                        APMM1180
      NAN=N+N                                                           APMM1190
      K=L+L+L                                                           APMM1200
      J=K+NAN                                                           APMM1210
      DO 7 I=1,NAN                                                      APMM1220
      K=K+1                                                             APMM1230
      IHE(K)=I                                                          APMM1240
      J=J+1                                                             APMM1250
    7 IHE(J)=I                                                          APMM1260
C                                                                       APMM1270
C        SET COUNTER ITER FOR ITERATION-STEPS                           APMM1280
      ITER=-1                                                           APMM1290
    8 ITER=ITER+1                                                       APMM1300
C                                                                       APMM1310
C        TEST FOR MAXIMUM ITERATION-STEPS                               APMM1320
      IF(N+M-ITER) 9,9,10                                               APMM1330
    9 IER=1                                                             APMM1340
      GO TO 69                                                          APMM1350
C                                                                       APMM1360
C        DETERMINE THE COLUMN WITH THE MOST POSITIVE ELEMENT IN TOP     APMM1370
   10 ISE=0                                                             APMM1380
      IPIV=0                                                            APMM1390
      K=L+L+L                                                           APMM1400
      SAVE=0.                                                           APMM1410
C                                                                       APMM1420
C        START TOP-LOOP                                                 APMM1430
      DO 14 I=1,NAN                                                     APMM1440
      IDO=K+I                                                           APMM1450
      HELP=TOP(I)                                                       APMM1460
      IF(HELP-SAVE) 12,12,11                                            APMM1470
   11 SAVE=HELP                                                         APMM1480
      IPIV=I                                                            APMM1490
   12 IF(IHE(IDO)) 14,13,14                                             APMM1500
   13 ISE=I                                                             APMM1510
   14 CONTINUE                                                          APMM1520
C        END OF TOP-LOOP                                                APMM1530
C                                                                       APMM1540
C        IS OPTIMAL TABLEAU REACHED                                     APMM1550
      IF(IPIV) 69,69,15                                                 APMM1560
C                                                                       APMM1570
C        DETERMINE THE PIVOT-ELEMENT FOR THE COLUMN CHOSEN UPOVE        APMM1580
   15 ILAB=1                                                            APMM1590
      IND=0                                                             APMM1600
      J=ISE                                                             APMM1610
      IF(J) 21,21,34                                                    APMM1620
C                                                                       APMM1630
C        TRANSFER K-TH COLUMN FROM T TO PIV                             APMM1640
   16 K=(K-1)*L                                                         APMM1650
      DO 17 I=1,L                                                       APMM1660
      J=L+I                                                             APMM1670
      K=K+1                                                             APMM1680
   17 PIV(J)=T(K)                                                       APMM1690
C                                                                       APMM1700
C        IS ANOTHER COLUMN NEEDED FOR SEARCH FOR PIVOT-ELEMENT          APMM1710
   18 IF(ISE) 22,22,19                                                  APMM1720
   19 ISE=-ISE                                                          APMM1730
C                                                                       APMM1740
C        TRANSFER COLUMNS IN PIV                                        APMM1750
      J=L+1                                                             APMM1760
      IDO=L+L                                                           APMM1770
      DO 20 I=J,IDO                                                     APMM1780
      K=I+L                                                             APMM1790
   20 PIV(K)=PIV(I)                                                     APMM1800
   21 J=IPIV                                                            APMM1810
      GO TO 34                                                          APMM1820
C                                                                       APMM1830
C        SEARCH PIVOT-ELEMENT PIV(IND)                                  APMM1840
   22 SAVE=1.E38                                                        APMM1850
      IDO=0                                                             APMM1860
      K=L+1                                                             APMM1870
      LL=L+L                                                            APMM1880
      IND=0                                                             APMM1890
C                                                                       APMM1900
C        START PIVOT-LOOP                                               APMM1910
      DO 29 I=K,LL                                                      APMM1920
      J=I+L                                                             APMM1930
      HELP=PIV(I)                                                       APMM1940
      IF(HELP) 29,29,23                                                 APMM1950
   23 HELP=-HELP                                                        APMM1960
      IF(ISE) 26,24,26                                                  APMM1970
   24 IF(IHE(J)) 27,25,27                                               APMM1980
   25 IDO=I                                                             APMM1990
      GO TO 29                                                          APMM2000
   26 HELP=-PIV(J)/HELP                                                 APMM2010
   27 IF(HELP-SAVE) 28,29,29                                            APMM2020
   28 SAVE=HELP                                                         APMM2030
      IND=I                                                             APMM2040
   29 CONTINUE                                                          APMM2050
C        END OF PIVOT-LOOP                                              APMM2060
C                                                                       APMM2070
C        TEST FOR SUITABLE PIVOT-ELEMENT                                APMM2080
      IF(IND) 30,30,32                                                  APMM2090
   30 IF(IDO) 68,68,31                                                  APMM2100
   31 IND=IDO                                                           APMM2110
C        PIVOT-ELEMENT IS STORED IN PIV(IND)                            APMM2120
C                                                                       APMM2130
C        COMPUTE THE RECIPROCAL OF THE PIVOT-ELEMENT REPI               APMM2140
   32 REPI=1./PIV(IND)                                                  APMM2150
      IND=IND-L                                                         APMM2160
C                                                                       APMM2170
C        UPDATE THE TOP-ROW TOP OF THE TABLEAU                          APMM2180
      ILAB=0                                                            APMM2190
      SAVE=-TOP(IPIV)*REPI                                              APMM2200
      TOP(IPIV)=SAVE                                                    APMM2210
C                                                                       APMM2220
C        INITIALIZE J AS COUNTER FOR TOP-LOOP                           APMM2230
      J=NAN                                                             APMM2240
   33 IF(J-IPIV) 34,53,34                                               APMM2250
   34 K=0                                                               APMM2260
C                                                                       APMM2270
C        SEARCH COLUMN IN TRANSFORMATION-MATRIX T                       APMM2280
      DO 36 I=1,L                                                       APMM2290
      IF(IHE(I)-J) 36,35,36                                             APMM2300
   35 K=I                                                               APMM2310
      IF(ILAB) 50,50,16                                                 APMM2320
   36 CONTINUE                                                          APMM2330
C                                                                       APMM2340
C        GENERATE COLUMN USING SUBROUTINE FCT AND TRANSFORMATION-MATRIX APMM2350
      I=L+L+L+NAN+J                                                     APMM2360
      I=IHE(I)-N                                                        APMM2370
      IF(I) 37,37,38                                                    APMM2380
   37 I=I+N                                                             APMM2390
      K=1                                                               APMM2400
   38 I=I+NAN                                                           APMM2410
C                                                                       APMM2420
C        CALL SUBROUTINE FCT                                            APMM2430
      CALL FCT(PIV,TOP(I),M-1)                                          APMM2440
C                                                                       APMM2450
C        PREPARE THE CALLED VECTOR PIV                                  APMM2460
      DSUM=0.D0                                                         APMM2470
      IDO=M                                                             APMM2480
      DO 41 I=1,M                                                       APMM2490
      HELP=PIV(IDO)                                                     APMM2500
      IF(K) 39,39,40                                                    APMM2510
   39 HELP=-HELP                                                        APMM2520
   40 DSUM=DSUM+DBLE(HELP)                                              APMM2530
      PIV(IDO+1)=HELP                                                   APMM2540
   41 IDO=IDO-1                                                         APMM2550
      PIV(L)=-DSUM                                                      APMM2560
      PIV(1)=1.                                                         APMM2570
C                                                                       APMM2580
C        TRANSFORM VECTOR PIV WITH ROWS OF MATRIX T                     APMM2590
      IDO=IND                                                           APMM2600
      IF(ILAB) 44,44,42                                                 APMM2610
   42 K=1                                                               APMM2620
   43 IDO=K                                                             APMM2630
   44 DSUM=0.D0                                                         APMM2640
      HELP=0.                                                           APMM2650
C                                                                       APMM2660
C        START MULTIPLICATION-LOOP                                      APMM2670
      DO 46 I=1,L                                                       APMM2680
      DSUM=DSUM+DBLE(PIV(I)*T(IDO))                                     APMM2690
      TOL=ABS(SNGL(DSUM))                                               APMM2700
      IF(TOL-HELP) 46,46,45                                             APMM2710
   45 HELP=TOL                                                          APMM2720
   46 IDO=IDO+L                                                         APMM2730
C        END OF MULTIPLICATION-LOOP                                     APMM2740
C                                                                       APMM2750
      TOL=1.E-5*HELP                                                    APMM2760
      IF(ABS(SNGL(DSUM))-TOL) 47,47,48                                  APMM2770
   47 DSUM=0.D0                                                         APMM2780
   48 IF(ILAB) 51,51,49                                                 APMM2790
   49 I=K+L                                                             APMM2800
      PIV(I)=DSUM                                                       APMM2810
C                                                                       APMM2820
C        TEST FOR LAST COLUMN-TERM                                      APMM2830
      K=K+1                                                             APMM2840
      IF(K-L) 43,43,18                                                  APMM2850
   50 I=(K-1)*L+IND                                                     APMM2860
      DSUM=T(I)                                                         APMM2870
C                                                                       APMM2880
C        COMPUTE NEW TOP-ELEMENT                                        APMM2890
   51 DSUM=DSUM*DBLE(SAVE)                                              APMM2900
      TOL=1.E-5*ABS(SNGL(DSUM))                                         APMM2910
      TOP(J)=TOP(J)+SNGL(DSUM)                                          APMM2920
      IF(ABS(TOP(J))-TOL) 52,52,53                                      APMM2930
   52 TOP(J)=0.                                                         APMM2940
C                                                                       APMM2950
C        TEST FOR LAST TOP-TERM                                         APMM2960
   53 J=J-1                                                             APMM2970
      IF(J) 54,54,33                                                    APMM2980
C        END OF TOP-LOOP                                                APMM2990
C                                                                       APMM3000
C        TRANSFORM PIVOT-COLUMN                                         APMM3010
   54 I=IND+L                                                           APMM3020
      PIV(I)=-1.                                                        APMM3030
      DO 55 I=1,L                                                       APMM3040
      J=I+L                                                             APMM3050
   55 PIV(I)=-PIV(J)*REPI                                               APMM3060
C                                                                       APMM3070
C        UPDATE TRANSFORMATION-MATRIX T                                 APMM3080
      J=0                                                               APMM3090
      DO 57 I=1,L                                                       APMM3100
      IDO=J+IND                                                         APMM3110
      SAVE=T(IDO)                                                       APMM3120
      T(IDO)=0.                                                         APMM3130
      DO 56 K=1,L                                                       APMM3140
      ISE=K+J                                                           APMM3150
   56 T(ISE)=T(ISE)+SAVE*PIV(K)                                         APMM3160
   57 J=J+L                                                             APMM3170
C                                                                       APMM3180
C        UPDATE INDEX-VECTOR IHE                                        APMM3190
C        INITIALIZE CHARACTERISTICS                                     APMM3200
      J=0                                                               APMM3210
      K=0                                                               APMM3220
      ISE=0                                                             APMM3230
      IDO=0                                                             APMM3240
C                                                                       APMM3250
C        START QUESTION-LOOP                                            APMM3260
      DO 61 I=1,L                                                       APMM3270
      LL=I+L                                                            APMM3280
      ILAB=IHE(LL)                                                      APMM3290
      IF(IHE(I)-IPIV) 59,58,59                                          APMM3300
   58 ISE=I                                                             APMM3310
      J=ILAB                                                            APMM3320
   59 IF(ILAB-IND) 61,60,61                                             APMM3330
   60 IDO=I                                                             APMM3340
      K=IHE(I)                                                          APMM3350
   61 CONTINUE                                                          APMM3360
C        END OF QUESTION-LOOP                                           APMM3370
C                                                                       APMM3380
C        START MODIFICATION                                             APMM3390
      IF(K) 62,62,63                                                    APMM3400
   62 IHE(IDO)=IPIV                                                     APMM3410
      IF(ISE) 67,67,65                                                  APMM3420
   63 IF(IND-J) 64,66,64                                                APMM3430
   64 LL=L+L+L+NAN                                                      APMM3440
      K=K+LL                                                            APMM3450
      I=IPIV+LL                                                         APMM3460
      ILAB=IHE(K)                                                       APMM3470
      IHE(K)=IHE(I)                                                     APMM3480
      IHE(I)=ILAB                                                       APMM3490
      IF(ISE) 67,67,65                                                  APMM3500
   65 IDO=IDO+L                                                         APMM3510
      I=ISE+L                                                           APMM3520
      IHE(IDO)=J                                                        APMM3530
      IHE(I)=IND                                                        APMM3540
   66 IHE(ISE)=0                                                        APMM3550
   67 LL=L+L                                                            APMM3560
      J=LL+IND                                                          APMM3570
      I=LL+L+IPIV                                                       APMM3580
      ILAB=IHE(I)                                                       APMM3590
      IHE(I)=IHE(J)                                                     APMM3600
      IHE(J)=ILAB                                                       APMM3610
C        END OF MODIFICATION                                            APMM3620
C                                                                       APMM3630
      GO TO 8                                                           APMM3640
C                                                                       APMM3650
C        SET ERROR PARAMETER IER=-1 SINCE NO SUITABLE PIVOT IS FOUND    APMM3660
   68 IER=-1                                                            APMM3670
C                                                                       APMM3680
C        EVALUATE FINAL TABLEAU                                         APMM3690
C        COMPUTE SAVE AS MAXIMUM ERROR OF APPROXIMATION AND             APMM3700
C        HELP AS ADDITIVE CONSTANCE FOR RESULTING COEFFICIENTS          APMM3710
   69 SAVE=0.                                                           APMM3720
      HELP=0.                                                           APMM3730
      K=L+L+L                                                           APMM3740
      DO 73 I=1,NAN                                                     APMM3750
      IDO=K+I                                                           APMM3760
      J=IHE(IDO)                                                        APMM3770
      IF(J) 71,70,73                                                    APMM3780
   70 SAVE=-TOP(I)                                                      APMM3790
   71 IF(M+J+1) 73,72,73                                                APMM3800
   72 HELP=TOP(I)                                                       APMM3810
   73 CONTINUE                                                          APMM3820
C                                                                       APMM3830
C        PREPARE T,TOP,PIV                                              APMM3840
      T(1)=SAVE                                                         APMM3850
      IDO=NAN+1                                                         APMM3860
      J=NAN+N                                                           APMM3870
      DO 74 I=IDO,J                                                     APMM3880
   74 TOP(I)=SAVE                                                       APMM3890
      DO 75 I=1,M                                                       APMM3900
   75 PIV(I)=HELP                                                       APMM3910
C                                                                       APMM3920
C        COMPUTE COEFFICIENTS OF RESULTING POLYNOMIAL IN PIV(1) UP TO PIAPMM3930
C        AND CALCULATE ERRORS AT GIVEN NODES IN TOP(1) UP TO TOP(N)     APMM3940
      DO 79 I=1,NAN                                                     APMM3950
      IDO=K+I                                                           APMM3960
      J=IHE(IDO)                                                        APMM3970
      IF(J) 76,79,77                                                    APMM3980
   76 J=-J                                                              APMM3990
      PIV(J)=HELP-TOP(I)                                                APMM4000
      GO TO 79                                                          APMM4010
   77 IF(J-N) 78,78,79                                                  APMM4020
   78 J=J+NAN                                                           APMM4030
      TOP(J)=SAVE+TOP(I)                                                APMM4040
   79 CONTINUE                                                          APMM4050
      DO 80 I=1,N                                                       APMM4060
      IDO=NAN+I                                                         APMM4070
   80 TOP(I)=TOP(IDO)                                                   APMM4080
   81 RETURN                                                            APMM4090
      END                                                               APMM4100
