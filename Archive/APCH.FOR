C                                                                       APCH  10
C     ..................................................................APCH  20
C                                                                       APCH  30
C        SUBROUTINE APCH                                                APCH  40
C                                                                       APCH  50
C        PURPOSE                                                        APCH  60
C           SET UP NORMAL EQUATIONS OF LEAST SQUARES FIT IN TERMS OF    APCH  70
C           CHEBYSHEV POLYNOMIALS FOR A GIVEN DISCRETE FUNCTION         APCH  80
C                                                                       APCH  90
C        USAGE                                                          APCH 100
C           CALL APCH(DATI,N,IP,XD,X0,WORK,IER)                         APCH 110
C                                                                       APCH 120
C        DESCRIPTION OF PARAMETERS                                      APCH 130
C           DATI  - VECTOR OF DIMENSION 3*N (OR DIMENSION 2*N+1)        APCH 140
C                   CONTAINING THE GIVEN ARGUMENTS, FOLLOWED BY THE     APCH 150
C                   FUNCTION VALUES AND N (RESPECTIVELY 1) WEIGHT       APCH 160
C                   VALUES. THE CONTENT OF VECTOR DATI REMAINS          APCH 170
C                   UNCHANGED.                                          APCH 180
C           N     - NUMBER OF GIVEN POINTS                              APCH 190
C           IP    - DIMENSION OF LEAST SQUARES FIT, I.E. NUMBER OF      APCH 200
C                   CHEBYSHEV POLYNOMIALS USED AS FUNDAMENTAL FUNCTIONS APCH 210
C                   IP SHOULD NOT EXCEED N                              APCH 220
C           XD    - RESULTANT MULTIPLICATIVE CONSTANT FOR LINEAR        APCH 230
C                   TRANSFORMATION OF ARGUMENT RANGE                    APCH 240
C           X0    - RESULTANT ADDITIVE CONSTANT FOR LINEAR              APCH 250
C                   TRANSFORMATION OF ARGUMENT RANGE                    APCH 260
C           WORK  - WORKING STORAGE OF DIMENSION (IP+1)*(IP+2)/2        APCH 270
C                   ON RETURN WORK CONTAINS THE SYMMETRIC COEFFICIENT   APCH 280
C                   MATRIX OF THE NORMAL EQUATIONS IN COMPRESSED FORM   APCH 290
C                   FOLLOWED IMMEDIATELY BY RIGHT HAND SIDE             APCH 300
C                   AND SQUARE SUM OF FUNCTION VALUES                   APCH 310
C           IER   - RESULTING ERROR PARAMETER                           APCH 320
C                   IER =-1 MEANS FORMAL ERRORS IN DIMENSION            APCH 330
C                   IER = 0 MEANS NO ERRORS                             APCH 340
C                   IER = 1 MEANS COINCIDING ARGUMENTS                  APCH 350
C                                                                       APCH 360
C        REMARKS                                                        APCH 370
C           NO WEIGHTS ARE USED IF THE VALUE OF DATI(2*N+1) IS          APCH 380
C           NOT POSITIVE.                                               APCH 390
C           EXECUTION OF SUBROUTINE APCH IS A PREPARATORY STEP FOR      APCH 400
C           CALCULATION OF LEAST SQUARES FITS IN CHEBYSHEV POLYNOMIALS  APCH 410
C           IT SHOULD BE FOLLOWED BY EXECUTION OF SUBROUTINE APFS       APCH 420
C                                                                       APCH 430
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  APCH 440
C           NONE                                                        APCH 450
C                                                                       APCH 460
C        METHOD                                                         APCH 470
C           THE LEAST SQUARE FIT IS DETERMINED USING CHEBYSHEV          APCH 480
C           POLYNOMIALS AS FUNDAMENTAL FUNCTION SYSTEM.                 APCH 490
C           THE METHOD IS DISCUSSED IN THE ARTICLE                      APCH 500
C           A.T.BERZTISS, LEAST SQUARES FITTING TO IRREGULARLY SPACED   APCH 510
C           DATA, SIAM REVIEW, VOL.6, ISS.3, 1964, PP. 203-227.         APCH 520
C                                                                       APCH 530
C     ..................................................................APCH 540
C                                                                       APCH 550
      SUBROUTINE APCH(DATI,N,IP,XD,X0,WORK,IER)                         APCH 560
C                                                                       APCH 570
C                                                                       APCH 580
C       DIMENSIONED DUMMY VARIABLES                                     APCH 590
      DIMENSION DATI(1),WORK(1)                                         APCH 600
C                                                                       APCH 610
C        CHECK FOR FORMAL ERRORS IN SPECIFIED DIMENSIONS                APCH 620
      IF(N-1)19,20,1                                                    APCH 630
    1 IF(IP)19,19,2                                                     APCH 640
C                                                                       APCH 650
C        SEARCH SMALLEST AND LARGEST ARGUMENT                           APCH 660
    2 IF(IP-N)3,3,19                                                    APCH 670
    3 XA=DATI(1)                                                        APCH 680
      X0=XA                                                             APCH 690
      XE=0.                                                             APCH 700
      DO 7 I=1,N                                                        APCH 710
      XM=DATI(I)                                                        APCH 720
      IF(XA-XM)5,5,4                                                    APCH 730
    4 XA=XM                                                             APCH 740
    5 IF(X0-XM)6,7,7                                                    APCH 750
    6 X0=XM                                                             APCH 760
    7 CONTINUE                                                          APCH 770
C                                                                       APCH 780
C        INITIALIZE CALCULATION OF NORMAL EQUATIONS                     APCH 790
      XD=X0-XA                                                          APCH 800
      M=(IP*(IP+1))/2                                                   APCH 810
      IEND=M+IP+1                                                       APCH 820
      MT2=IP+IP                                                         APCH 830
      MT2M=MT2-1                                                        APCH 840
C                                                                       APCH 850
C        SET WORKING STORAGE AND RIGHT HAND SIDE TO ZERO                APCH 860
      DO 8 I=1,IP                                                       APCH 870
      J=MT2-I                                                           APCH 880
      WORK(J)=0.                                                        APCH 890
      WORK(I)=0.                                                        APCH 900
      K=M+I                                                             APCH 910
    8 WORK(K)=0.                                                        APCH 920
C                                                                       APCH 930
C        CHECK FOR DEGENERATE ARGUMENT RANGE                            APCH 940
      IF(XD)20,20,9                                                     APCH 950
C                                                                       APCH 960
C        CALCULATE CONSTANTS FOR REDUCTION OF ARGUMENTS                 APCH 970
    9 X0=-(X0+XA)/XD                                                    APCH 980
      XD=2./XD                                                          APCH 990
      SUM=0.                                                            APCH1000
C                                                                       APCH1010
C        START GREAT LOOP OVER ALL GIVEN POINTS                         APCH1020
      DO 15 I=1,N                                                       APCH1030
      T=DATI(I)*XD+X0                                                   APCH1040
      J=I+N                                                             APCH1050
      DF=DATI(J)                                                        APCH1060
C                                                                       APCH1070
C        CALCULATE AND STORE VALUES OF CHEBYSHEV POLYNOMIALS            APCH1080
C        FOR ARGUMENT T                                                 APCH1090
      XA=1.                                                             APCH1100
      XM=T                                                              APCH1110
      IF(DATI(2*N+1))11,11,10                                           APCH1120
   10 J=J+N                                                             APCH1130
      XA=DATI(J)                                                        APCH1140
      XM=T*XA                                                           APCH1150
   11 T=T+T                                                             APCH1160
      SUM=SUM+DF*DF*XA                                                  APCH1170
      DF=DF+DF                                                          APCH1180
      J=1                                                               APCH1190
   12 K=M+J                                                             APCH1200
      WORK(K)=WORK(K)+DF*XA                                             APCH1210
   13 WORK(J)=WORK(J)+XA                                                APCH1220
      IF(J-MT2M)14,15,15                                                APCH1230
   14 J=J+1                                                             APCH1240
      XE=T*XM-XA                                                        APCH1250
      XA=XM                                                             APCH1260
      XM=XE                                                             APCH1270
      IF(J-IP)12,12,13                                                  APCH1280
   15 CONTINUE                                                          APCH1290
      WORK(IEND)=SUM+SUM                                                APCH1300
C                                                                       APCH1310
C        CALCULATE MATRIX OF NORMAL EQUATIONS                           APCH1320
      LL=M                                                              APCH1330
      KK=MT2M                                                           APCH1340
      JJ=1                                                              APCH1350
      K=KK                                                              APCH1360
      DO 18 J=1,M                                                       APCH1370
      WORK(LL)=WORK(K)+WORK(JJ)                                         APCH1380
      LL=LL-1                                                           APCH1390
      IF(K-JJ)16,16,17                                                  APCH1400
   16 KK=KK-2                                                           APCH1410
      K=KK                                                              APCH1420
      JJ=1                                                              APCH1430
      GOTO 18                                                           APCH1440
   17 JJ=JJ+1                                                           APCH1450
      K=K-1                                                             APCH1460
   18 CONTINUE                                                          APCH1470
      IER=0                                                             APCH1480
      RETURN                                                            APCH1490
C                                                                       APCH1500
C        ERROR RETURN IN CASE OF FORMAL ERRORS                          APCH1510
   19 IER=-1                                                            APCH1520
      RETURN                                                            APCH1530
C                                                                       APCH1540
C        ERROR RETURN IN CASE OF COINCIDING ARGUMENTS                   APCH1550
   20 IER=1                                                             APCH1560
      RETURN                                                            APCH1570
      END                                                               APCH1580
