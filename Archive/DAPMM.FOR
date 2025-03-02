C                                                                       DAPM  10
C     ..................................................................DAPM  20
C                                                                       DAPM  30
C        SUBROUTINE DAPMM                                               DAPM  40
C                                                                       DAPM  50
C        PURPOSE                                                        DAPM  60
C           APPROXIMATE A FUNCTION TABULATED IN N POINTS BY ANY LINEAR  DAPM  70
C           COMBINATION OF M GIVEN CONTINUOUS FUNCTIONS IN THE SENSE    DAPM  80
C           OF CHEBYSHEV.                                               DAPM  90
C                                                                       DAPM 100
C        USAGE                                                          DAPM 110
C           CALL DAPMM(FCT,N,M,TOP,IHE,PIV,T,ITER,IER)                  DAPM 120
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT IN THE         DAPM 130
C           CALLING PROGRAM.                                            DAPM 140
C                                                                       DAPM 150
C        DESCRIPTION OF PARAMETERS                                      DAPM 160
C           FCT    - NAME OF SUBROUTINE TO BE SUPPLIED BY THE USER.     DAPM 170
C                    IT COMPUTES VALUES OF M GIVEN FUNCTIONS FOR        DAPM 180
C                    ARGUMENT VALUE X.                                  DAPM 190
C                    USAGE                                              DAPM 200
C                       CALL FCT(Y,X,K)                                 DAPM 210
C                    DESCRIPTION OF PARAMETERS                          DAPM 220
C                       Y   - RESULT VECTOR OF DIMENSION M CONTAINING   DAPM 230
C                             THE VALUES OF GIVEN CONTINUOUS FUNCTIONS  DAPM 240
C                             FOR GIVEN ARGUMENT X                      DAPM 250
C                       X   - ARGUMENT VALUE                            DAPM 260
C                       K   - AN INTEGER VALUE WHICH IS EQUAL TO M-1    DAPM 270
C                    REMARKS                                            DAPM 280
C                       IF APPROXIMATION BY NORMAL CHEBYSHEV, SHIFTED   DAPM 290
C                       CHEBYSHEV, LEGENDRE, LAGUERRE, HERMITE POLYNO-  DAPM 300
C                       MIALS IS DESIRED SUBROUTINES CNP, CSP, LEP,     DAPM 310
C                       LAP, HEP, RESPECTIVELY FROM SSP COULD BE USED.  DAPM 320
C           N      - NUMBER OF DATA POINTS DEFINING THE FUNCTION WHICH  DAPM 330
C                    IS TO BE APPROXIMATED                              DAPM 340
C           M      - NUMBER OF GIVEN CONTINUOUS FUNCTIONS FROM WHICH    DAPM 350
C                    THE APPROXIMATING FUNCTION IS CONSTRUCTED.         DAPM 360
C           TOP    - VECTOR OF DIMENSION 3*N.                           DAPM 370
C                    ON ENTRY IT MUST CONTAIN FROM TOP(1) UP TO TOP(N)  DAPM 380
C                    THE GIVEN N FUNCTION VALUES AND FROM TOP(N+1) UP   DAPM 390
C                    TO TOP(2*N) THE CORRESPONDING NODES                DAPM 400
C                    ON RETURN TOP CONTAINS FROM TOP(1) UP TO TOP(N)    DAPM 410
C                    THE ERRORS AT THOSE N NODES.                       DAPM 420
C                    OTHER VALUES OF TOP ARE SCRATCH.                   DAPM 430
C           IHE    - INTEGER VECTOR OF DIMENSION 3*M+4*N+6              DAPM 440
C           PIV    - VECTOR OF DIMENSION 3*M+6.                         DAPM 450
C                    ON RETURN PIV CONTAINS AT PIV(1) UP TO PIV(M) THE  DAPM 460
C                    RESULTING COEFFICIENTS OF LINEAR APPROXIMATION.    DAPM 470
C           T      - AUXILIARY VECTOR OF DIMENSION (M+2)*(M+2)          DAPM 480
C           ITER   - RESULTANT INTEGER WHICH SPECIFIES THE NUMBER OF    DAPM 490
C                    ITERATIONS NEEDED                                  DAPM 500
C           IER    - RESULTANT ERROR PARAMETER CODED IN THE FOLLOWING   DAPM 510
C                    FORM                                               DAPM 520
C                     IER=0  - NO ERROR                                 DAPM 530
C                     IER=1  - THE NUMBER OF ITERATIONS HAS REACHED     DAPM 540
C                              THE INTERNAL MAXIMUM N+M                 DAPM 550
C                     IER=-1 - NO RESULT BECAUSE OF WRONG INPUT PARA-   DAPM 560
C                              METER M OR N OR SINCE AT SOME ITERATION  DAPM 570
C                              NO SUITABLE PIVOT COULD BE FOUND         DAPM 580
C        NOTE: DSUM,TOP,PIV,T,SAVE,HELP,RECI,TOL ARE DOUBLE PRECISION   DAPM 590
C        REMARKS                                                        DAPM 600
C           NO ACTION BESIDES ERROR MESSAGE IN CASE M LESS THAN 1 OR    DAPM 610
C           N LESS THAN 2.                                              DAPM 620
C                                                                       DAPM 630
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DAPM 640
C           THE EXTERNAL SUBROUTINE FCT MUST BE FURNISHED BY THE USER.  DAPM 650
C                                                                       DAPM 660
C        METHOD                                                         DAPM 670
C           THE PROBLEM OF APPROXIMATION A TABULATED FUNCTION BY ANY    DAPM 680
C           LINEAR COMBINATION OF GIVEN FUNCTIONS IN THE SENSE OF       DAPM 690
C           CHEBYSHEV (I.E. TO MINIMIZE THE MAXIMUM ERROR) IS TRANS-    DAPM 700
C           FORMED INTO A LINEAR PROGRAMMING PROBLEM. APMM USES A       DAPM 710
C           REVISED SIMPLEX METHOD TO SOLVE A CORRESPONDING DUAL        DAPM 720
C           PROBLEM. FOR REFERENCE, SEE                                 DAPM 730
C           I.BARRODALE/A.YOUNG, ALGORITHMS FOR BEST L-SUB-ONE AND      DAPM 740
C           L-SUB-INFINITY, LINEAR APPROXIMATIONS ON A DISCRETE SET,    DAPM 750
C           NUMERISCHE MATHEMATIK, VOL.8, ISS.3 (1966), PP.295-306.     DAPM 760
C                                                                       DAPM 770
C     ..................................................................DAPM 780
C                                                                       DAPM 790
      SUBROUTINE DAPMM(FCT,N,M,TOP,IHE,PIV,T,ITER,IER)                  DAPM 800
C                                                                       DAPM 810
      DOUBLE PRECISION DSUM,TOP,PIV,T,SAVE,HELP,REPI,TOL                DAPM 820
      DIMENSION TOP(1),IHE(1),PIV(1),T(1)                               DAPM 830
      DOUBLE PRECISION DSUM                                             DAPM 840
C                                                                       DAPM 850
C        TEST ON WRONG INPUT PARAMETERS N AND M                         DAPM 860
      IER=-1                                                            DAPM 870
      IF (N-1) 81,81,1                                                  DAPM 880
    1 IF(M) 81,81,2                                                     DAPM 890
C                                                                       DAPM 900
C        INITIALIZE CHARACTERISTIC VECTORS FOR THE TABLEAU              DAPM 910
    2 IER=0                                                             DAPM 920
C                                                                       DAPM 930
C        PREPARE TOP-ROW TOP                                            DAPM 940
      DO 3 I=1,N                                                        DAPM 950
      K=I+N                                                             DAPM 960
      J=K+N                                                             DAPM 970
      TOP(J)=TOP(K)                                                     DAPM 980
    3 TOP(K)=-TOP(I)                                                    DAPM 990
C                                                                       DAPM1000
C        PREPARE INVERSE TRANSFORMATION MATRIX T                        DAPM1010
      L=M+2                                                             DAPM1020
      LL=L*L                                                            DAPM1030
      DO 4 I=1,LL                                                       DAPM1040
    4 T(I)=0.                                                           DAPM1050
      K=1                                                               DAPM1060
      J=L+1                                                             DAPM1070
      DO 5 I=1,L                                                        DAPM1080
      T(K)=1.                                                           DAPM1090
    5 K=K+J                                                             DAPM1100
C                                                                       DAPM1110
C        PREPARE INDEX-VECTOR IHE                                       DAPM1120
      DO 6 I=1,L                                                        DAPM1130
      K=I+L                                                             DAPM1140
      J=K+L                                                             DAPM1150
      IHE(I)=0                                                          DAPM1160
      IHE(K)=I                                                          DAPM1170
    6 IHE(J)=1-I                                                        DAPM1180
      NAN=N+N                                                           DAPM1190
      K=L+L+L                                                           DAPM1200
      J=K+NAN                                                           DAPM1210
      DO 7 I=1,NAN                                                      DAPM1220
      K=K+1                                                             DAPM1230
      IHE(K)=I                                                          DAPM1240
      J=J+1                                                             DAPM1250
    7 IHE(J)=I                                                          DAPM1260
C                                                                       DAPM1270
C        SET COUNTER ITER FOR ITERATION-STEPS                           DAPM1280
      ITER=-1                                                           DAPM1290
    8 ITER=ITER+1                                                       DAPM1300
C                                                                       DAPM1310
C        TEST FOR MAXIMUM ITERATION-STEPS                               DAPM1320
      IF(N+M-ITER) 9,9,10                                               DAPM1330
    9 IER=1                                                             DAPM1340
      GO TO 69                                                          DAPM1350
C                                                                       DAPM1360
C        DETERMINE THE COLUMN WITH THE MOST POSITIVE ELEMENT IN TOP     DAPM1370
   10 ISE=0                                                             DAPM1380
      IPIV=0                                                            DAPM1390
      K=L+L+L                                                           DAPM1400
      SAVE=0.0D0                                                        DAPM1410
C                                                                       DAPM1420
C        START TOP-LOOP                                                 DAPM1430
      DO 14 I=1,NAN                                                     DAPM1440
      IDO=K+I                                                           DAPM1450
      HELP=TOP(I)                                                       DAPM1460
      IF(HELP-SAVE) 12,12,11                                            DAPM1470
   11 SAVE=HELP                                                         DAPM1480
      IPIV=I                                                            DAPM1490
   12 IF(IHE(IDO)) 14,13,14                                             DAPM1500
   13 ISE=I                                                             DAPM1510
   14 CONTINUE                                                          DAPM1520
C        END OF TOP-LOOP                                                DAPM1530
C                                                                       DAPM1540
C        IS OPTIMAL TABLEAU REACHED                                     DAPM1550
      IF(IPIV) 69,69,15                                                 DAPM1560
C                                                                       DAPM1570
C        DETERMINE THE PIVOT-ELEMENT FOR THE COLUMN CHOSEN UPOVE        DAPM1580
   15 ILAB=1                                                            DAPM1590
      IND=0                                                             DAPM1600
      J=ISE                                                             DAPM1610
      IF(J) 21,21,34                                                    DAPM1620
C                                                                       DAPM1630
C        TRANSFER K-TH COLUMN FROM T TO PIV                             DAPM1640
   16 K=(K-1)*L                                                         DAPM1650
      DO 17 I=1,L                                                       DAPM1660
      J=L+I                                                             DAPM1670
      K=K+1                                                             DAPM1680
   17 PIV(J)=T(K)                                                       DAPM1690
C                                                                       DAPM1700
C        IS ANOTHER COLUMN NEEDED FOR SEARCH FOR PIVOT-ELEMENT          DAPM1710
   18 IF(ISE) 22,22,19                                                  DAPM1720
   19 ISE=-ISE                                                          DAPM1730
C                                                                       DAPM1740
C        TRANSFER COLUMNS IN PIV                                        DAPM1750
      J=L+1                                                             DAPM1760
      IDO=L+L                                                           DAPM1770
      DO 20 I=J,IDO                                                     DAPM1780
      K=I+L                                                             DAPM1790
   20 PIV(K)=PIV(I)                                                     DAPM1800
   21 J=IPIV                                                            DAPM1810
      GO TO 34                                                          DAPM1820
C                                                                       DAPM1830
C        SEARCH PIVOT-ELEMENT PIV(IND)                                  DAPM1840
   22 SAVE=1.D38                                                        DAPM1850
      IDO=0                                                             DAPM1860
      K=L+1                                                             DAPM1870
      LL=L+L                                                            DAPM1880
      IND=0                                                             DAPM1890
C                                                                       DAPM1900
C        START PIVOT-LOOP                                               DAPM1910
      DO 29 I=K,LL                                                      DAPM1920
      J=I+L                                                             DAPM1930
      HELP=PIV(I)                                                       DAPM1940
      IF(HELP) 29,29,23                                                 DAPM1950
   23 HELP=-HELP                                                        DAPM1960
      IF(ISE) 26,24,26                                                  DAPM1970
   24 IF(IHE(J)) 27,25,27                                               DAPM1980
   25 IDO=I                                                             DAPM1990
      GO TO 29                                                          DAPM2000
   26 HELP=-PIV(J)/HELP                                                 DAPM2010
   27 IF(HELP-SAVE) 28,29,29                                            DAPM2020
   28 SAVE=HELP                                                         DAPM2030
      IND=I                                                             DAPM2040
   29 CONTINUE                                                          DAPM2050
C        END OF PIVOT-LOOP                                              DAPM2060
C                                                                       DAPM2070
C        TEST FOR SUITABLE PIVOT-ELEMENT                                DAPM2080
      IF(IND) 30,30,32                                                  DAPM2090
   30 IF(IDO) 68,68,31                                                  DAPM2100
   31 IND=IDO                                                           DAPM2110
C        PIVOT-ELEMENT IS STORED IN PIV(IND)                            DAPM2120
C                                                                       DAPM2130
C        COMPUTE THE RECIPROCAL OF THE PIVOT-ELEMENT REPI               DAPM2140
   32 REPI=1.0D0/PIV(IND)                                               DAPM2150
      IND=IND-L                                                         DAPM2160
C                                                                       DAPM2170
C        UPDATE THE TOP-ROW TOP OF THE TABLEAU                          DAPM2180
      ILAB=0                                                            DAPM2190
      SAVE=-TOP(IPIV)*REPI                                              DAPM2200
      TOP(IPIV)=SAVE                                                    DAPM2210
C                                                                       DAPM2220
C        INITIALIZE J AS COUNTER FOR TOP-LOOP                           DAPM2230
      J=NAN                                                             DAPM2240
   33 IF(J-IPIV) 34,53,34                                               DAPM2250
   34 K=0                                                               DAPM2260
C                                                                       DAPM2270
C        SEARCH COLUMN IN TRANSFORMATION-MATRIX T                       DAPM2280
      DO 36 I=1,L                                                       DAPM2290
      IF(IHE(I)-J) 36,35,36                                             DAPM2300
   35 K=I                                                               DAPM2310
      IF(ILAB) 50,50,16                                                 DAPM2320
   36 CONTINUE                                                          DAPM2330
C                                                                       DAPM2340
C        GENERATE COLUMN USING SUBROUTINE FCT AND TRANSFORMATION-MATRIX DAPM2350
      I=L+L+L+NAN+J                                                     DAPM2360
      I=IHE(I)-N                                                        DAPM2370
      IF(I) 37,37,38                                                    DAPM2380
   37 I=I+N                                                             DAPM2390
      K=1                                                               DAPM2400
   38 I=I+NAN                                                           DAPM2410
C                                                                       DAPM2420
C        CALL SUBROUTINE FCT                                            DAPM2430
      CALL FCT(PIV,TOP(I),M-1)                                          DAPM2440
C                                                                       DAPM2450
C        PREPARE THE CALLED VECTOR PIV                                  DAPM2460
      DSUM=0.D0                                                         DAPM2470
      IDO=M                                                             DAPM2480
      DO 41 I=1,M                                                       DAPM2490
      HELP=PIV(IDO)                                                     DAPM2500
      IF(K) 39,39,40                                                    DAPM2510
   39 HELP=-HELP                                                        DAPM2520
   40 DSUM=DSUM+DBLE(HELP)                                              DAPM2530
      PIV(IDO+1)=HELP                                                   DAPM2540
   41 IDO=IDO-1                                                         DAPM2550
      PIV(L)=-DSUM                                                      DAPM2560
      PIV(1)=1.0D0                                                      DAPM2570
C                                                                       DAPM2580
C        TRANSFORM VECTOR PIV WITH ROWS OF MATRIX T                     DAPM2590
      IDO=IND                                                           DAPM2600
      IF(ILAB) 44,44,42                                                 DAPM2610
   42 K=1                                                               DAPM2620
   43 IDO=K                                                             DAPM2630
   44 DSUM=0.D0                                                         DAPM2640
      HELP=0.0D0                                                        DAPM2650
C                                                                       DAPM2660
C        START MULTIPLICATION-LOOP                                      DAPM2670
      DO 46 I=1,L                                                       DAPM2680
      DSUM=DSUM+DBLE(PIV(I)*T(IDO))                                     DAPM2690
      TOL=DABS(SNGL(DSUM))                                              DAPM2700
      IF(TOL-HELP) 46,46,45                                             DAPM2710
   45 HELP=TOL                                                          DAPM2720
   46 IDO=IDO+L                                                         DAPM2730
C        END OF MULTIPLICATION-LOOP                                     DAPM2740
C                                                                       DAPM2750
      TOL=1.D-5*HELP                                                    DAPM2760
      IF(DABS(SNGL(DSUM))-TOL) 47,47,48                                 DAPM2770
   47 DSUM=0.D0                                                         DAPM2780
   48 IF(ILAB) 51,51,49                                                 DAPM2790
   49 I=K+L                                                             DAPM2800
      PIV(I)=DSUM                                                       DAPM2810
C                                                                       DAPM2820
C        TEST FOR LAST COLUMN-TERM                                      DAPM2830
      K=K+1                                                             DAPM2840
      IF(K-L) 43,43,18                                                  DAPM2850
   50 I=(K-1)*L+IND                                                     DAPM2860
      DSUM=T(I)                                                         DAPM2870
C                                                                       DAPM2880
C        COMPUTE NEW TOP-ELEMENT                                        DAPM2890
   51 DSUM=DSUM*DBLE(SAVE)                                              DAPM2900
      TOL=1.D-5*DABS(SNGL(DSUM))                                        DAPM2910
      TOP(J)=TOP(J)+SNGL(DSUM)                                          DAPM2920
      IF(DABS(TOP(J))-TOL) 52,52,53                                     DAPM2930
   52 TOP(J)=0.D0                                                       DAPM2940
C                                                                       DAPM2950
C        TEST FOR LAST TOP-TERM                                         DAPM2960
   53 J=J-1                                                             DAPM2970
      IF(J) 54,54,33                                                    DAPM2980
C        END OF TOP-LOOP                                                DAPM2990
C                                                                       DAPM3000
C        TRANSFORM PIVOT-COLUMN                                         DAPM3010
   54 I=IND+L                                                           DAPM3020
      PIV(I)=-1.                                                        DAPM3030
      DO 55 I=1,L                                                       DAPM3040
      J=I+L                                                             DAPM3050
   55 PIV(I)=-PIV(J)*REPI                                               DAPM3060
C                                                                       DAPM3070
C        UPDATE TRANSFORMATION-MATRIX T                                 DAPM3080
      J=0                                                               DAPM3090
      DO 57 I=1,L                                                       DAPM3100
      IDO=J+IND                                                         DAPM3110
      SAVE=T(IDO)                                                       DAPM3120
      T(IDO)=0.D0                                                       DAPM3130
      DO 56 K=1,L                                                       DAPM3140
      ISE=K+J                                                           DAPM3150
   56 T(ISE)=T(ISE)+SAVE*PIV(K)                                         DAPM3160
   57 J=J+L                                                             DAPM3170
C                                                                       DAPM3180
C        UPDATE INDEX-VECTOR IHE                                        DAPM3190
C        INITIALIZE CHARACTERISTICS                                     DAPM3200
      J=0                                                               DAPM3210
      K=0                                                               DAPM3220
      ISE=0                                                             DAPM3230
      IDO=0                                                             DAPM3240
C                                                                       DAPM3250
C        START QUESTION-LOOP                                            DAPM3260
      DO 61 I=1,L                                                       DAPM3270
      LL=I+L                                                            DAPM3280
      ILAB=IHE(LL)                                                      DAPM3290
      IF(IHE(I)-IPIV) 59,58,59                                          DAPM3300
   58 ISE=I                                                             DAPM3310
      J=ILAB                                                            DAPM3320
   59 IF(ILAB-IND) 61,60,61                                             DAPM3330
   60 IDO=I                                                             DAPM3340
      K=IHE(I)                                                          DAPM3350
   61 CONTINUE                                                          DAPM3360
C        END OF QUESTION-LOOP                                           DAPM3370
C                                                                       DAPM3380
C        START MODIFICATION                                             DAPM3390
      IF(K) 62,62,63                                                    DAPM3400
   62 IHE(IDO)=IPIV                                                     DAPM3410
      IF(ISE) 67,67,65                                                  DAPM3420
   63 IF(IND-J) 64,66,64                                                DAPM3430
   64 LL=L+L+L+NAN                                                      DAPM3440
      K=K+LL                                                            DAPM3450
      I=IPIV+LL                                                         DAPM3460
      ILAB=IHE(K)                                                       DAPM3470
      IHE(K)=IHE(I)                                                     DAPM3480
      IHE(I)=ILAB                                                       DAPM3490
      IF(ISE) 67,67,65                                                  DAPM3500
   65 IDO=IDO+L                                                         DAPM3510
      I=ISE+L                                                           DAPM3520
      IHE(IDO)=J                                                        DAPM3530
      IHE(I)=IND                                                        DAPM3540
   66 IHE(ISE)=0                                                        DAPM3550
   67 LL=L+L                                                            DAPM3560
      J=LL+IND                                                          DAPM3570
      I=LL+L+IPIV                                                       DAPM3580
      ILAB=IHE(I)                                                       DAPM3590
      IHE(I)=IHE(J)                                                     DAPM3600
      IHE(J)=ILAB                                                       DAPM3610
C        END OF MODIFICATION                                            DAPM3620
C                                                                       DAPM3630
      GO TO 8                                                           DAPM3640
C                                                                       DAPM3650
C        SET ERROR PARAMETER IER=-1 SINCE NO SUITABLE PIVOT IS FOUND    DAPM3660
   68 IER=-1                                                            DAPM3670
C                                                                       DAPM3680
C        EVALUATE FINAL TABLEAU                                         DAPM3690
C        COMPUTE SAVE AS MAXIMUM ERROR OF APPROXIMATION AND             DAPM3700
C        HELP AS ADDITIVE CONSTANCE FOR RESULTING COEFFICIENTS          DAPM3710
   69 SAVE=0.D0                                                         DAPM3720
      HELP=0.D0                                                         DAPM3730
      K=L+L+L                                                           DAPM3740
      DO 73 I=1,NAN                                                     DAPM3750
      IDO=K+I                                                           DAPM3760
      J=IHE(IDO)                                                        DAPM3770
      IF(J) 71,70,73                                                    DAPM3780
   70 SAVE=-TOP(I)                                                      DAPM3790
   71 IF(M+J+1) 73,72,73                                                DAPM3800
   72 HELP=TOP(I)                                                       DAPM3810
   73 CONTINUE                                                          DAPM3820
C                                                                       DAPM3830
C        PREPARE T,TOP,PIV                                              DAPM3840
      T(1)=SAVE                                                         DAPM3850
      IDO=NAN+1                                                         DAPM3860
      J=NAN+N                                                           DAPM3870
      DO 74 I=IDO,J                                                     DAPM3880
   74 TOP(I)=SAVE                                                       DAPM3890
      DO 75 I=1,M                                                       DAPM3900
   75 PIV(I)=HELP                                                       DAPM3910
C                                                                       DAPM3920
C        COMPUTE COEFFICIENTS OF RESULTING POLYNOMIAL IN PIV(1) UP TO PIDAPM3930
C        AND CALCULATE ERRORS AT GIVEN NODES IN TOP(1) UP TO TOP(N)     DAPM3940
      DO 79 I=1,NAN                                                     DAPM3950
      IDO=K+I                                                           DAPM3960
      J=IHE(IDO)                                                        DAPM3970
      IF(J) 76,79,77                                                    DAPM3980
   76 J=-J                                                              DAPM3990
      PIV(J)=HELP-TOP(I)                                                DAPM4000
      GO TO 79                                                          DAPM4010
   77 IF(J-N) 78,78,79                                                  DAPM4020
   78 J=J+NAN                                                           DAPM4030
      TOP(J)=SAVE+TOP(I)                                                DAPM4040
   79 CONTINUE                                                          DAPM4050
      DO 80 I=1,N                                                       DAPM4060
      IDO=NAN+I                                                         DAPM4070
   80 TOP(I)=TOP(IDO)                                                   DAPM4080
   81 RETURN                                                            DAPM4090
      END                                                               DAPM4100
