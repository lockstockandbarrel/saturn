C                                                                       FRAT  10
C     ..................................................................FRAT  20
C                                                                       FRAT  30
C        SUBROUTINE FRAT                                                FRAT  40
C                                                                       FRAT  50
C        PURPOSE                                                        FRAT  60
C           FRAT IS USED FOR HANDLING OF DATA AND FUNDAMENTAL FUNCTIONS FRAT  70
C           WITH RATIONAL APPROXIMATION. IT IS A SUBSTANTIAL PART OF    FRAT  80
C           RATIONAL APPROXIMATION AND HAS NO MEANING INDEPENDENTLY     FRAT  90
C                                                                       FRAT 100
C        USAGE                                                          FRAT 110
C           CALL FRAT(I,N,M,P,DATI,WGT,IER)                             FRAT 120
C                                                                       FRAT 130
C        DESCRIPTION OF PARAMETERS                                      FRAT 140
C           I     - SUBSCRIPT OF CURRENT DATA POINT                     FRAT 150
C           N     - NUMBER OF ALL DATA POINTS                           FRAT 160
C           M     - NUMBER OF FUNDAMENTAL FUNCTIONS USED                FRAT 170
C           P     - ARRAY OF DIMENSION M+1 AT LEAST, WHICH CONTAINS     FRAT 180
C                   ON RETURN THE VALUES OF THE M FUNDAMENTAL           FRAT 190
C                   FUNCTIONS, FOLLOWED BY CURRENT FUNCTION VALUE       FRAT 200
C           DATI  - ARRAY CONTAINING GIVEN N ARGUMENTS, FOLLOWED        FRAT 210
C                   BY N FUNCTION VALUES AND FINALLY BY 1 RESPECTIVELY  FRAT 220
C                   N WEIGHT VALUES                                     FRAT 230
C           WGT   - RESULTANT WEIGHT FACTOR USED FOR I-TH TERM          FRAT 240
C           IER   - RESULTANT ERROR PARAMETER, COMBINED WITH INPUT      FRAT 250
C                   VALUES FOR CONTROL                                  FRAT 260
C                   IER(2) MEANS DIMENSION OF NUMERATOR                 FRAT 270
C                   IER(3) MEANS DIMENSION OF DENOMINATOR               FRAT 280
C                   IER(1) IS USED AS RESULTANT ERROR PARAMETER,        FRAT 290
C                   IER(1) = 0 IN CASE OF NO ERRORS                     FRAT 300
C                   IER(1) = 1 OTHERWISE (ZERO DENOMINATOR)             FRAT 310
C                                                                       FRAT 320
C        REMARKS                                                        FRAT 330
C           VECTOR IER IS USED FOR COMMUNICATION BETWEEN ARAT AND FRAT  FRAT 340
C                                                                       FRAT 350
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  FRAT 360
C           CNP                                                         FRAT 370
C                                                                       FRAT 380
C        METHOD                                                         FRAT 390
C           CF. MATHEMATICAL DESCRIPTION OF SUBROUTINE ARAT             FRAT 400
C                                                                       FRAT 410
C     ..................................................................FRAT 420
C                                                                       FRAT 430
      SUBROUTINE FRAT(I,N,M,P,DATI,WGT,IER)                             FRAT 440
C                                                                       FRAT 450
C                                                                       FRAT 460
C        DIMENSIONED DUMMY VARIABLES                                    FRAT 470
      DIMENSION P(1),DATI(1),IER(1)                                     FRAT 480
C                                                                       FRAT 490
C        INITIALIZATION                                                 FRAT 500
      IP=IER(2)                                                         FRAT 510
      IQ=IER(3)                                                         FRAT 520
      IQM1=IQ-1                                                         FRAT 530
      IPQ=IP+IQ                                                         FRAT 540
C                                                                       FRAT 550
C        LOOK UP ARGUMENT, FUNCTION VALUE AND WEIGHT                    FRAT 560
C        LOOK UP NUMERATOR AND DENOMINATOR                              FRAT 570
      T=DATI(I)                                                         FRAT 580
      J=I+N                                                             FRAT 590
      F=DATI(J)                                                         FRAT 600
      FNUM=P(J)                                                         FRAT 610
      J=J+N                                                             FRAT 620
      WGT=1.                                                            FRAT 630
      IF(DATI(2*N+1))2,2,1                                              FRAT 640
    1 WGT=DATI(J)                                                       FRAT 650
    2 FDEN=P(J)                                                         FRAT 660
C                                                                       FRAT 670
C        CALCULATE FUNCTION VALUE USED                                  FRAT 680
      F=F*FDEN-FNUM                                                     FRAT 690
C                                                                       FRAT 700
C        CHECK FOR ZERO DENOMINATOR                                     FRAT 710
      IF(FDEN)4,3,4                                                     FRAT 720
C                                                                       FRAT 730
C        ERROR RETURN IN CASE OF ZERO DENOMINATOR                       FRAT 740
    3 IER(1)=1                                                          FRAT 750
      RETURN                                                            FRAT 760
C                                                                       FRAT 770
C        CALCULATE WEIGHT FACTORS USED                                  FRAT 780
    4 WGT=WGT/(FDEN*FDEN)                                               FRAT 790
      FNUM=-FNUM/FDEN                                                   FRAT 800
C                                                                       FRAT 810
C        CALCULATE FUNDAMENTAL FUNCTIONS                                FRAT 820
      J=IQM1                                                            FRAT 830
      IF(IP-IQ)6,6,5                                                    FRAT 840
    5 J=IP-1                                                            FRAT 850
    6 CALL CNP(P(IQ),T,J)                                               FRAT 860
C                                                                       FRAT 870
C        STORE VALUES OF DENOMINATOR FUNDAMENTAL FUNCTIONS              FRAT 880
    7 IF(IQM1)10,10,8                                                   FRAT 890
    8 DO 9 II=1,IQM1                                                    FRAT 900
      J=II+IQ                                                           FRAT 910
    9 P(II)=P(J)*FNUM                                                   FRAT 920
C                                                                       FRAT 930
C        STORE FUNCTION VALUE                                           FRAT 940
   10 P(IPQ)=F                                                          FRAT 950
C                                                                       FRAT 960
C        NORMAL RETURN                                                  FRAT 970
      IER(1)=0                                                          FRAT 980
      RETURN                                                            FRAT 990
      END                                                               FRAT1000
