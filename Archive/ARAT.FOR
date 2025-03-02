C                                                                       ARAT  10
C     ..................................................................ARAT  20
C                                                                       ARAT  30
C        SUBROUTINE ARAT                                                ARAT  40
C                                                                       ARAT  50
C        PURPOSE                                                        ARAT  60
C           CALCULATE BEST RATIONAL APPROXIMATION OF A DISCRETE         ARAT  70
C           FUNCTION IN THE LEAST SQUARES SENSE                         ARAT  80
C                                                                       ARAT  90
C        USAGE                                                          ARAT 100
C           CALL ARAT(DATI,N,WORK,P,IP,IQ,IER)                          ARAT 110
C                                                                       ARAT 120
C        DESCRIPTION OF PARAMETERS                                      ARAT 130
C           DATI  - TWODIMENSIONAL ARRAY WITH 3 COLUMNS AND N ROWS      ARAT 140
C                   THE FIRST COLUMN MUST CONTAIN THE GIVEN ARGUMENTS,  ARAT 150
C                   THE SECOND COLUMN THE GIVEN FUNCTION VALUES AND     ARAT 160
C                   THE THIRD COLUMN THE GIVEN WEIGHTS IF ANY.          ARAT 170
C                   IF NO WEIGHTS ARE TO BE USED THEN THE THIRD         ARAT 180
C                   COLUMN MAY BE DROPPED , EXCEPT THE FIRST ELEMENT    ARAT 190
C                   WHICH MUST CONTAIN A NONPOSITIVE VALUE              ARAT 200
C           N     - NUMBER OF NODES OF THE GIVEN DISCRETE FUNCTION      ARAT 210
C           WORK  - WORKING STORAGE WHICH IS OF DIMENSION               ARAT 220
C                   (IP+IQ)*(IP+IQ+1)+4*N+1 AT LEAST.                   ARAT 230
C                   ON RETURN THE VALUES OF THE NUMERATOR ARE CONTAINED ARAT 240
C                   IN WORK(N+1) UP TO WORK(2*N), WHILE THE VALUES OF   ARAT 250
C                   THE DENOMINATOR ARE STORED IN WORK(2*N+1) UP TO     ARAT 260
C                   WORK(3*N)                                           ARAT 270
C           P     - RESULTANT COEFFICIENT VECTOR OF DENOMINATOR AND     ARAT 280
C                   NUMERATOR. THE DENOMINATOR IS STORED IN FIRST IQ    ARAT 290
C                   LOCATIONS, THE NUMERATOR IN THE FOLLOWING IP        ARAT 300
C                   LOCATIONS.                                          ARAT 310
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH.          ARAT 320
C           IP    - DIMENSION OF THE NUMERATOR   (INPUT VALUE)          ARAT 330
C           IQ    - DIMENSION OF THE DENOMINATOR (INPUT VALUE)          ARAT 340
C           IER   - RESULTANT ERROR PARAMETER                           ARAT 350
C                   IER =-1 MEANS FORMAL ERRORS                         ARAT 360
C                   IER = 0 MEANS NO ERRORS                             ARAT 370
C                   IER = 1,2 MEANS POOR CONVERGENCE OF ITERATION       ARAT 380
C                   IER IS ALSO USED AS INPUT VALUE                     ARAT 390
C                   A NONZERO INPUT VALUE INDICATES AVAILABILITY OF AN  ARAT 400
C                   INITIAL APPROXIMATION STORED IN P                   ARAT 410
C                                                                       ARAT 420
C        REMARKS                                                        ARAT 430
C           THE COEFFICIENT VECTORS OF THE DENOMINATOR AND NUMERATOR    ARAT 440
C           OF THE RATIONAL APPROXIMATION ARE BOTH STORED IN P          ARAT 450
C           STARTING WITH LOW POWERS (DENOMINATOR FIRST).               ARAT 460
C           IP+IQ MUST NOT EXCEED N, ALL THREE VALUES MUST BE POSITIVE. ARAT 470
C           SINCE CHEBYSHEV POLYNOMIALS ARE USED AS FUNDAMENTAL         ARAT 480
C           FUNCTIONS, THE ARGUMENTS SHOULD BE REDUCED TO THE INTERVAL  ARAT 490
C           (-1,1). THIS CAN ALWAYS BE ACCOMPLISHED BY MEANS OF A LINEARARAT 500
C           TRANSFORMATION OF THE ORIGINALLY GIVEN ARGUMENTS.           ARAT 510
C           IF A FIT IN OTHER FUNCTIONS IS REQUIRED, CNP AND CNPS MUST  ARAT 520
C           BE REPLACED BY SUBROUTINES WHICH ARE OF ANALOGOUS DESIGN.   ARAT 530
C                                                                       ARAT 540
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  ARAT 550
C           APLL, APFS, FRAT, CNPS, CNP                                 ARAT 560
C           CNP IS REQUIRED WITHIN FRAT                                 ARAT 570
C                                                                       ARAT 580
C        METHOD                                                         ARAT 590
C           THE ITERATIVE SCHEME USED FOR CALCULATION OF THE            ARAT 600
C           APPROXIMATION IS REPEATED SOLUTION OF THE NORMAL EQUATIONS  ARAT 610
C           WHICH ARE OBTAINED BY LINEARIZATION.                        ARAT 620
C           A REFINED TECHNIQUE OF THIS LINEAR LEAST SQUARES APPROACH   ARAT 630
C           IS USED WHICH GUARANTEES THAT THE DENOMINATOR IS FREE OF    ARAT 640
C           ZEROES WITHIN THE APPROXIMATION INTERVAL.                   ARAT 650
C           FOR REFERENCE SEE                                           ARAT 660
C           D.BRAESS, UEBER DAEMPFUNG BEI MINIMALISIERUNGSVERFAHREN,    ARAT 670
C           COMPUTING(1966), VOL.1, ED.3, PP.264-272.                   ARAT 680
C           D.W.MARQUARDT, AN ALGORITHM FOR LEAST-SQUARES ESTIMATION    ARAT 690
C           OF NONLINEAR PARAMETERS,                                    ARAT 700
C           JSIAM(1963), VOL.11, ED.2, PP.431-441.                      ARAT 710
C                                                                       ARAT 720
C     ..................................................................ARAT 730
C                                                                       ARAT 740
      SUBROUTINE ARAT(DATI,N,WORK,P,IP,IQ,IER)                          ARAT 750
C                                                                       ARAT 760
C                                                                       ARAT 770
      EXTERNAL FRAT                                                     ARAT 780
C                                                                       ARAT 790
C        DIMENSIONED LOCAL VARIABLE                                     ARAT 800
      DIMENSION IERV(3)                                                 ARAT 810
C                                                                       ARAT 820
C        DIMENSIONED DUMMY VARIABLES                                    ARAT 830
      DIMENSION DATI(1),WORK(1),P(1)                                    ARAT 840
C                                                                       ARAT 850
C        INITIALIZE TESTVALUES                                          ARAT 860
      LIMIT=20                                                          ARAT 870
      ETA =1.E-11                                                       ARAT 880
      EPS=1.E-5                                                         ARAT 890
C                                                                       ARAT 900
C        CHECK FOR FORMAL ERRORS                                        ARAT 910
      IF(N)4,4,1                                                        ARAT 920
    1 IF(IP)4,4,2                                                       ARAT 930
    2 IF(IQ)4,4,3                                                       ARAT 940
    3 IPQ=IP+IQ                                                         ARAT 950
      IF(N-IPQ)4,5,5                                                    ARAT 960
C                                                                       ARAT 970
C        ERROR RETURN IN CASE OF FORMAL ERRORS                          ARAT 980
    4 IER=-1                                                            ARAT 990
      RETURN                                                            ARAT1000
C                                                                       ARAT1010
C        INITIALIZE ITERATION PROCESS                                   ARAT1020
    5 KOUNT=0                                                           ARAT1030
      IERV(2)=IP                                                        ARAT1040
      IERV(3)=IQ                                                        ARAT1050
      NDP=N+N+1                                                         ARAT1060
      NNE=NDP+NDP                                                       ARAT1070
      IX=IPQ-1                                                          ARAT1080
      IQP1=IQ+1                                                         ARAT1090
      IRHS=NNE+IPQ*IX/2                                                 ARAT1100
      IEND=IRHS+IX                                                      ARAT1110
C                                                                       ARAT1120
C        TEST FOR AVAILABILITY OF AN INITIAL APPROXIMATION              ARAT1130
      IF(IER)8,6,8                                                      ARAT1140
C                                                                       ARAT1150
C        INITIALIZE NUMERATOR AND DENOMINATOR                           ARAT1160
    6 DO 7 I=2,IPQ                                                      ARAT1170
    7 P(I)=0.                                                           ARAT1180
      P(1)=1.                                                           ARAT1190
C                                                                       ARAT1200
C        CALCULATE VALUES OF NUMERATOR AND DENOMINATOR FOR INITIAL      ARAT1210
C        APPROXIMATION                                                  ARAT1220
    8 DO 9 J=1,N                                                        ARAT1230
      T=DATI(J)                                                         ARAT1240
      I=J+N                                                             ARAT1250
      CALL CNPS(WORK(I),T,P(IQP1),IP)                                   ARAT1260
      K=I+N                                                             ARAT1270
    9 CALL CNPS(WORK(K),T,P,IQ)                                         ARAT1280
C                                                                       ARAT1290
C        SET UP NORMAL EQUATIONS (MAIN LOOP OF ITERATION)               ARAT1300
   10 CALL APLL(FRAT,N,IX,WORK,WORK(IEND+1),DATI,IERV)                  ARAT1310
C                                                                       ARAT1320
C        CHECK FOR ZERO DENOMINATOR                                     ARAT1330
      IF(IERV(1))4,11,4                                                 ARAT1340
   11 INCR=0                                                            ARAT1350
      RELAX=2.                                                          ARAT1360
C                                                                       ARAT1370
C        RESTORE MATRIX IN WORKING STORAGE                              ARAT1380
   12 J=IEND                                                            ARAT1390
      DO 13 I=NNE,IEND                                                  ARAT1400
      J=J+1                                                             ARAT1410
   13 WORK(I)=WORK(J)                                                   ARAT1420
      IF(KOUNT)14,14,15                                                 ARAT1430
C                                                                       ARAT1440
C        SAVE SQUARE SUM OF ERRORS                                      ARAT1450
   14 OSUM=WORK(IEND)                                                   ARAT1460
      DIAG=OSUM*EPS                                                     ARAT1470
      K=IQ                                                              ARAT1480
C                                                                       ARAT1490
C        ADD CONSTANT TO DIAGONAL                                       ARAT1500
      IF(WORK(NNE))17,17,19                                             ARAT1510
   15 IF(INCR)19,19,16                                                  ARAT1520
   16 K=IPQ                                                             ARAT1530
   17 J=NNE-1                                                           ARAT1540
      DO 18 I=1,K                                                       ARAT1550
      WORK(J)=WORK(J)+DIAG                                              ARAT1560
   18 J=J+I                                                             ARAT1570
C                                                                       ARAT1580
C        SOLVE NORMAL EQUATIONS                                         ARAT1590
   19 CALL APFS(WORK(NNE),IX,IRES,1,EPS,ETA,IER)                        ARAT1600
C                                                                       ARAT1610
C        CHECK FOR FAILURE OF EQUATION SOLVER                           ARAT1620
      IF(IRES)4,4,20                                                    ARAT1630
C                                                                       ARAT1640
C        TEST FOR DEFECTIVE NORMALEQUATIONS                             ARAT1650
   20 IF(IRES-IX)21,24,24                                               ARAT1660
   21 IF(INCR)22,22,23                                                  ARAT1670
   22 DIAG=DIAG*0.125                                                   ARAT1680
   23 DIAG=DIAG+DIAG                                                    ARAT1690
      INCR=INCR+1                                                       ARAT1700
C                                                                       ARAT1710
C        START WITH OVER RELAXATION                                     ARAT1720
      RELAX=8.                                                          ARAT1730
      IF(INCR-LIMIT)12,45,45                                            ARAT1740
C                                                                       ARAT1750
C        CALCULATE VALUES OF CHANGE OF NUMERATOR AND DENOMINATOR        ARAT1760
   24 L=NDP                                                             ARAT1770
      J=NNE+IRES*(IRES-1)/2-1                                           ARAT1780
      K=J+IQ                                                            ARAT1790
      WORK(J)=0.                                                        ARAT1800
      IRQ=IQ                                                            ARAT1810
      IRP=IRES-IQ+1                                                     ARAT1820
      IF(IRP)25,26,26                                                   ARAT1830
   25 IRQ=IRES+1                                                        ARAT1840
   26 DO 29 I=1,N                                                       ARAT1850
      T=DATI(I)                                                         ARAT1860
      WORK(I)=0.                                                        ARAT1870
      CALL CNPS(WORK(I),T,WORK(K),IRP)                                  ARAT1880
      M=L+N                                                             ARAT1890
      CALL CNPS(WORK(M),T,WORK(J),IRQ)                                  ARAT1900
      IF(WORK(M)*WORK(L))27,29,29                                       ARAT1910
   27 SUM=WORK(L)/WORK(M)                                               ARAT1920
      IF(RELAX+SUM)29,29,28                                             ARAT1930
   28 RELAX=-SUM                                                        ARAT1940
   29 L=L+1                                                             ARAT1950
C                                                                       ARAT1960
C        MODIFY RELAXATION FACTOR IF NECESSARY                          ARAT1970
      SSOE=OSUM                                                         ARAT1980
      ITER=LIMIT                                                        ARAT1990
   30 SUM=0.                                                            ARAT2000
      RELAX=RELAX*0.5                                                   ARAT2010
      DO 32 I=1,N                                                       ARAT2020
      M=I+N                                                             ARAT2030
      K=M+N                                                             ARAT2040
      L=K+N                                                             ARAT2050
      SAVE=DATI(M)-(WORK(M)+RELAX*WORK(I))/(WORK(K)+RELAX*WORK(L))      ARAT2060
      SAVE=SAVE*SAVE                                                    ARAT2070
      IF(DATI(NDP))32,32,31                                             ARAT2080
   31 SAVE=SAVE*DATI(K)                                                 ARAT2090
   32 SUM=SUM+SAVE                                                      ARAT2100
      IF(ITER)45,33,33                                                  ARAT2110
   33 ITER=ITER-1                                                       ARAT2120
      IF(SUM-OSUM)34,37,35                                              ARAT2130
   34 OSUM=SUM                                                          ARAT2140
      GOTO 30                                                           ARAT2150
C                                                                       ARAT2160
C        TEST FOR IMPROVEMENT                                           ARAT2170
   35 IF(OSUM-SSOE)36,30,30                                             ARAT2180
   36 RELAX=RELAX+RELAX                                                 ARAT2190
   37 T=0.                                                              ARAT2200
      SAVE=0.                                                           ARAT2210
      K=IRES+1                                                          ARAT2220
      DO 38 I=2,K                                                       ARAT2230
      J=J+1                                                             ARAT2240
      T=T+ABS(P(I))                                                     ARAT2250
      P(I)=P(I)+RELAX*WORK(J)                                           ARAT2260
   38 SAVE=SAVE+ABS(P(I))                                               ARAT2270
C                                                                       ARAT2280
C        UPDATE CURRENT VALUES OF NUMERATOR AND DENOMINATOR             ARAT2290
      DO 39 I=1,N                                                       ARAT2300
      J=I+N                                                             ARAT2310
      K=J+N                                                             ARAT2320
      L=K+N                                                             ARAT2330
      WORK(J)=WORK(J)+RELAX*WORK(I)                                     ARAT2340
   39 WORK(K)=WORK(K)+RELAX*WORK(L)                                     ARAT2350
C                                                                       ARAT2360
C        TEST FOR CONVERGENCE                                           ARAT2370
      IF(INCR)40,40,42                                                  ARAT2380
   40 IF(SSOE-OSUM-RELAX*EPS*OSUM)46,46,41                              ARAT2390
   41 IF(ABS(T-SAVE)-RELAX*EPS*SAVE)46,46,42                            ARAT2400
   42 IF(OSUM-ETA*SAVE)46,46,43                                         ARAT2410
   43 KOUNT=KOUNT+1                                                     ARAT2420
      IF(KOUNT-LIMIT)10,44,44                                           ARAT2430
C                                                                       ARAT2440
C        ERROR RETURN IN CASE OF POOR CONVERGENCE                       ARAT2450
   44 IER=2                                                             ARAT2460
      RETURN                                                            ARAT2470
   45 IER=1                                                             ARAT2480
      RETURN                                                            ARAT2490
C                                                                       ARAT2500
C        NORMAL RETURN                                                  ARAT2510
   46 IER=0                                                             ARAT2520
      RETURN                                                            ARAT2530
      END                                                               ARAT2540
