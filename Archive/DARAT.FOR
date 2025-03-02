C                                                                       DARA  10
C     ..................................................................DARA  20
C                                                                       DARA  30
C        SUBROUTINE DARAT                                               DARA  40
C                                                                       DARA  50
C        PURPOSE                                                        DARA  60
C           CALCULATE BEST RATIONAL APPROXIMATION OF A DISCRETE         DARA  70
C           FUNCTION IN THE LEAST SQUARES SENSE                         DARA  80
C                                                                       DARA  90
C        USAGE                                                          DARA 100
C           CALL DARAT(DATI,N,WORK,P,IP,IQ,IER)                         DARA 110
C                                                                       DARA 120
C        DESCRIPTION OF PARAMETERS                                      DARA 130
C           DATI  - TWO DIMENSIONAL ARRAY WITH 3 COLUMNS AND N ROWS     DARA 140
C                   THE FIRST COLUMN MUST CONTAIN THE GIVEN ARGUMENTS,  DARA 150
C                   THE SECOND COLUMN THE GIVEN FUNCTION VALUES AND     DARA 160
C                   THE THIRD COLUMN THE GIVEN WEIGHTS IF ANY.          DARA 170
C                   IF NO WEIGHTS ARE TO BE USED THEN THE THIRD         DARA 180
C                   COLUMN MAY BE DROPPED , EXCEPT THE FIRST ELEMENT    DARA 190
C                   WHICH MUST CONTAIN A NONPOSITIVE VALUE              DARA 200
C           N     - NUMBER OF NODES OF THE GIVEN DISCRETE FUNCTION      DARA 210
C           WORK  - WORKING STORAGE WHICH IS OF DIMENSION               DARA 220
C                   (IP+IQ)*(IP+IQ+1)+4*N+1 AT LEAST.                   DARA 230
C                   ON RETURN THE VALUES OF THE NUMERATOR ARE CONTAINED DARA 240
C                   IN WORK(N+1) UP TO WORK(2*N), WHILE THE VALUES OF   DARA 250
C                   THE DENOMINATOR ARE STORED IN WORK(2*N+1) UP TO     DARA 260
C                   WORK(3*N)                                           DARA 270
C           P     - RESULTANT COEFFICIENT VECTOR OF DENOMINATOR AND     DARA 280
C                   NUMERATOR. THE DENOMINATOR IS STORED IN FIRST IQ    DARA 290
C                   LOCATIONS, THE NUMERATOR IN THE FOLLOWING IP        DARA 300
C                   LOCATIONS.                                          DARA 310
C                   COEFFICIENTS ARE ORDERED FROM LOW TO HIGH.          DARA 320
C           IP    - DIMENSION OF THE NUMERATOR   (INPUT VALUE)          DARA 330
C           IQ    - DIMENSION OF THE DENOMINATOR (INPUT VALUE)          DARA 340
C           IER   - RESULTANT ERROR PARAMETER                           DARA 350
C                   IER =-1 MEANS FORMAL ERRORS                         DARA 360
C                   IER = 0 MEANS NO ERRORS                             DARA 370
C                   IER = 1,2 MEANS POOR CONVERGENCE OF ITERATION       DARA 380
C                   IER IS ALSO USED AS INPUT VALUE                     DARA 390
C                   A NONZERO INPUT VALUE INDICATES AVAILABILITY OF AN  DARA 400
C                   INITIAL APPROXIMATION STORED IN P                   DARA 410
C        NOTE: DATI,WORK,P,T,OSUM,DIAG,RELAX,SUM,SSOE,SAVE (DOUB. PRE.) DARA 420
C        REMARKS                                                        DARA 430
C           THE COEFFICIENT VECTORS OF THE DENOMINATOR AND NUMERATOR    DARA 440
C           OF THE RATIONAL APPROXIMATION ARE BOTH STORED IN P          DARA 450
C           STARTING WITH LOW POWERS (DENOMINATOR FIRST).               DARA 460
C           IP+IQ MUST NOT EXCEED N, ALL THREE VALUES MUST BE POSITIVE. DARA 470
C           SINCE CHEBYSHEV POLYNOMIALS ARE USED AS FUNDAMENTAL         DARA 480
C           FUNCTIONS, THE ARGUMENTS SHOULD BE REDUCED TO THE INTERVAL  DARA 490
C           (-1,1). THIS CAN ALWAYS BE ACCOMPLISHED BY MEANS OF A LINEARDARA 500
C           TRANSFORMATION OF THE ORIGINALLY GIVEN ARGUMENTS.           DARA 510
C           IF A FIT IN OTHER FUNCTIONS IS REQUIRED, CNP AND CNPS MUST  DARA 520
C           BE REPLACED BY SUBROUTINES WHICH ARE OF ANALOGOUS DESIGN.   DARA 530
C                                                                       DARA 540
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DARA 550
C           DAPLL, DAPFS, DFRAT, DCNPS, DCNP                            DARA 560
C           DCNP IS REQUIRED WITHIN DFRAT                               DARA 570
C                                                                       DARA 580
C        METHOD                                                         DARA 590
C           THE ITERATIVE SCHEME USED FOR CALCULATION OF THE            DARA 600
C           APPROXIMATION IS REPEATED SOLUTION OF THE NORMAL EQUATIONS  DARA 610
C           WHICH ARE OBTAINED BY LINEARIZATION.                        DARA 620
C           A REFINED TECHNIQUE OF THIS LINEAR LEAST SQUARES APPROACH   DARA 630
C           IS USED WHICH GUARANTEES THAT THE DENOMINATOR IS FREE OF    DARA 640
C           ZEROES WITHIN THE APPROXIMATION INTERVAL.                   DARA 650
C           FOR REFERENCE SEE                                           DARA 660
C           D.BRAESS, UEBER DAEMPFUNG BEI MINIMALISIERUNGSVERFAHREN,    DARA 670
C           COMPUTING(1966), VOL.1, ED.3, PP.264-272.                   DARA 680
C           D.W.MARQUARDT, AN ALGORITHM FOR LEAST-SQUARES ESTIMATION    DARA 690
C           OF NONLINEAR PARAMETERS,                                    DARA 700
C           JSIAM(1963), VOL.11, ED.2, PP.431-441.                      DARA 710
C                                                                       DARA 720
C     ..................................................................DARA 730
C                                                                       DARA 740
      SUBROUTINE DARAT(DATI,N,WORK,P,IP,IQ,IER)                         DARA 750
C                                                                       DARA 760
C                                                                       DARA 770
      EXTERNAL DFRAT                                                    DARA 780
C                                                                       DARA 790
C        DIMENSIONED LOCAL VARIABLE                                     DARA 800
      DIMENSION IERV(3)                                                 DARA 810
C                                                                       DARA 820
C        DIMENSIONED DUMMY VARIABLES                                    DARA 830
      DIMENSION DATI(1),WORK(1),P(1)                                    DARA 840
      DOUBLE PRECISION DATI,WORK,P,T,OSUM,DIAG,RELAX,SUM,SSOE,SAVE      DARA 850
C        INITIALIZE TESTVALUES                                          DARA 860
      LIMIT=20                                                          DARA 870
      ETA =1.E-29                                                       DARA 880
      EPS=1.E-14                                                        DARA 890
C                                                                       DARA 900
C        CHECK FOR FORMAL ERRORS                                        DARA 910
      IF(N)4,4,1                                                        DARA 920
    1 IF(IP)4,4,2                                                       DARA 930
    2 IF(IQ)4,4,3                                                       DARA 940
    3 IPQ=IP+IQ                                                         DARA 950
      IF(N-IPQ)4,5,5                                                    DARA 960
C                                                                       DARA 970
C        ERROR RETURN IN CASE OF FORMAL ERRORS                          DARA 980
    4 IER=-1                                                            DARA 990
      RETURN                                                            DARA1000
C                                                                       DARA1010
C        INITIALIZE ITERATION PROCESS                                   DARA1020
    5 KOUNT=0                                                           DARA1030
      IERV(2)=IP                                                        DARA1040
      IERV(3)=IQ                                                        DARA1050
      NDP=N+N+1                                                         DARA1060
      NNE=NDP+NDP                                                       DARA1070
      IX=IPQ-1                                                          DARA1080
      IQP1=IQ+1                                                         DARA1090
      IRHS=NNE+IPQ*IX/2                                                 DARA1100
      IEND=IRHS+IX                                                      DARA1110
C                                                                       DARA1120
C        TEST FOR AVAILABILITY OF AN INITIAL APPROXIMATION              DARA1130
      IF(IER)8,6,8                                                      DARA1140
C                                                                       DARA1150
C        INITIALIZE NUMERATOR AND DENOMINATOR                           DARA1160
    6 DO 7 I=2,IPQ                                                      DARA1170
    7 P(I)=0.D0                                                         DARA1180
      P(1)=1.D0                                                         DARA1190
C                                                                       DARA1200
C        CALCULATE VALUES OF NUMERATOR AND DENOMINATOR FOR INITIAL      DARA1210
C        APPROXIMATION                                                  DARA1220
    8 DO 9 J=1,N                                                        DARA1230
      T=DATI(J)                                                         DARA1240
      I=J+N                                                             DARA1250
      CALL DCNPS(WORK(I),T,P(IQP1),IP)                                  DARA1260
      K=I+N                                                             DARA1270
    9 CALL DCNPS(WORK(K),T,P,IQ)                                        DARA1280
C                                                                       DARA1290
C        SET UP NORMAL EQUATIONS (MAIN LOOP OF ITERATION)               DARA1300
   10 CALL DAPLL(DFRAT,N,IX,WORK,WORK(IEND+1),DATI,IERV)                DARA1310
C                                                                       DARA1320
C        CHECK FOR ZERO DENOMINATOR                                     DARA1330
      IF(IERV(1))4,11,4                                                 DARA1340
   11 INCR=0                                                            DARA1350
      RELAX=2.D0                                                        DARA1360
C                                                                       DARA1370
C        RESTORE MATRIX IN WORKING STORAGE                              DARA1380
   12 J=IEND                                                            DARA1390
      DO 13 I=NNE,IEND                                                  DARA1400
      J=J+1                                                             DARA1410
   13 WORK(I)=WORK(J)                                                   DARA1420
      IF(KOUNT)14,14,15                                                 DARA1430
C                                                                       DARA1440
C        SAVE SQUARE SUM OF ERRORS                                      DARA1450
   14 OSUM=WORK(IEND)                                                   DARA1460
      DIAG=OSUM*EPS                                                     DARA1470
      K=IQ                                                              DARA1480
C                                                                       DARA1490
C        ADD CONSTANT TO DIAGONAL                                       DARA1500
      IF(WORK(NNE))17,17,19                                             DARA1510
   15 IF(INCR)19,19,16                                                  DARA1520
   16 K=IPQ                                                             DARA1530
   17 J=NNE-1                                                           DARA1540
      DO 18 I=1,K                                                       DARA1550
      WORK(J)=WORK(J)+DIAG                                              DARA1560
   18 J=J+I                                                             DARA1570
C                                                                       DARA1580
C        SOLVE NORMAL EQUATIONS                                         DARA1590
   19 CALL DAPFS(WORK(NNE),IX,IRES,1,EPS,ETA,IER)                       DARA1600
C                                                                       DARA1610
C        CHECK FOR FAILURE OF EQUATION SOLVER                           DARA1620
      IF(IRES)4,4,20                                                    DARA1630
C                                                                       DARA1640
C        TEST FOR DEFECTIVE NORMALEQUATIONS                             DARA1650
   20 IF(IRES-IX)21,24,24                                               DARA1660
   21 IF(INCR)22,22,23                                                  DARA1670
   22 DIAG=DIAG*0.125                                                   DARA1680
   23 DIAG=DIAG+DIAG                                                    DARA1690
      INCR=INCR+1                                                       DARA1700
C                                                                       DARA1710
C        START WITH OVER RELAXATION                                     DARA1720
      RELAX=8.D0                                                        DARA1730
      IF(INCR-LIMIT)12,45,45                                            DARA1740
C                                                                       DARA1750
C        CALCULATE VALUES OF CHANGE OF NUMERATOR AND DENOMINATOR        DARA1760
   24 L=NDP                                                             DARA1770
      J=NNE+IRES*(IRES-1)/2-1                                           DARA1780
      K=J+IQ                                                            DARA1790
      WORK(J)=0.D0                                                      DARA1800
      IRQ=IQ                                                            DARA1810
      IRP=IRES-IQ+1                                                     DARA1820
      IF(IRP)25,26,26                                                   DARA1830
   25 IRQ=IRES+1                                                        DARA1840
   26 DO 29 I=1,N                                                       DARA1850
      T=DATI(I)                                                         DARA1860
      WORK(I)=0.D0                                                      DARA1870
      CALL DCNPS(WORK(I),T,WORK(K),IRP)                                 DARA1880
      M=L+N                                                             DARA1890
      CALL DCNPS(WORK(M),T,WORK(J),IRQ)                                 DARA1900
      IF(WORK(M)*WORK(L))27,29,29                                       DARA1910
   27 SUM=WORK(L)/WORK(M)                                               DARA1920
      IF(RELAX+SUM)29,29,28                                             DARA1930
   28 RELAX=-SUM                                                        DARA1940
   29 L=L+1                                                             DARA1950
C                                                                       DARA1960
C        MODIFY RELAXATION FACTOR IF NECESSARY                          DARA1970
      SSOE=OSUM                                                         DARA1980
      ITER=LIMIT                                                        DARA1990
   30 SUM=0.D0                                                          DARA2000
      RELAX=RELAX*0.5D0                                                 DARA2010
      DO 32 I=1,N                                                       DARA2020
      M=I+N                                                             DARA2030
      K=M+N                                                             DARA2040
      L=K+N                                                             DARA2050
      SAVE=DATI(M)-(WORK(M)+RELAX*WORK(I))/(WORK(K)+RELAX*WORK(L))      DARA2060
      SAVE=SAVE*SAVE                                                    DARA2070
      IF(DATI(NDP))32,32,31                                             DARA2080
   31 SAVE=SAVE*DATI(K)                                                 DARA2090
   32 SUM=SUM+SAVE                                                      DARA2100
      IF(ITER)45,33,33                                                  DARA2110
   33 ITER=ITER-1                                                       DARA2120
      IF(SUM-OSUM)34,37,35                                              DARA2130
   34 OSUM=SUM                                                          DARA2140
      GOTO 30                                                           DARA2150
C                                                                       DARA2160
C        TEST FOR IMPROVEMENT                                           DARA2170
   35 IF(OSUM-SSOE)36,30,30                                             DARA2180
   36 RELAX=RELAX+RELAX                                                 DARA2190
   37 T=0.D0                                                            DARA2200
      SAVE=0.D0                                                         DARA2210
      K=IRES+1                                                          DARA2220
      DO 38 I=2,K                                                       DARA2230
      J=J+1                                                             DARA2240
      T=T+DABS(P(I))                                                    DARA2250
      P(I)=P(I)+RELAX*WORK(J)                                           DARA2260
   38 SAVE=SAVE+DABS(P(I))                                              DARA2270
C                                                                       DARA2280
C        UPDATE CURRENT VALUES OF NUMERATOR AND DENOMINATOR             DARA2290
      DO 39 I=1,N                                                       DARA2300
      J=I+N                                                             DARA2310
      K=J+N                                                             DARA2320
      L=K+N                                                             DARA2330
      WORK(J)=WORK(J)+RELAX*WORK(I)                                     DARA2340
   39 WORK(K)=WORK(K)+RELAX*WORK(L)                                     DARA2350
C                                                                       DARA2360
C        TEST FOR CONVERGENCE                                           DARA2370
      IF(INCR)40,40,42                                                  DARA2380
   40 IF(SSOE-OSUM-RELAX*EPS*OSUM)46,46,41                              DARA2390
   41 IF(DABS(T-SAVE)-RELAX*EPS*SAVE)46,46,42                           DARA2400
   42 IF(OSUM-ETA*SAVE)46,46,43                                         DARA2410
   43 KOUNT=KOUNT+1                                                     DARA2420
      IF(KOUNT-LIMIT)10,44,44                                           DARA2430
C                                                                       DARA2440
C        ERROR RETURN IN CASE OF POOR CONVERGENCE                       DARA2450
   44 IER=2                                                             DARA2460
      RETURN                                                            DARA2470
   45 IER=1                                                             DARA2480
      RETURN                                                            DARA2490
C                                                                       DARA2500
C        NORMAL RETURN                                                  DARA2510
   46 IER=0                                                             DARA2520
      RETURN                                                            DARA2530
      END                                                               DARA2540
