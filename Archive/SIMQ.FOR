C                                                                       SIMQ  10
C     ..................................................................SIMQ  20
C                                                                       SIMQ  30
C        SUBROUTINE SIMQ                                                SIMQ  40
C                                                                       SIMQ  50
C        PURPOSE                                                        SIMQ  60
C           OBTAIN SOLUTION OF A SET OF SIMULTANEOUS LINEAR EQUATIONS,  SIMQ  70
C           AX=B                                                        SIMQ  80
C                                                                       SIMQ  90
C        USAGE                                                          SIMQ 100
C           CALL SIMQ(A,B,N,KS)                                         SIMQ 110
C                                                                       SIMQ 120
C        DESCRIPTION OF PARAMETERS                                      SIMQ 130
C           A - MATRIX OF COEFFICIENTS STORED COLUMNWISE.  THESE ARE    SIMQ 140
C               DESTROYED IN THE COMPUTATION.  THE SIZE OF MATRIX A IS  SIMQ 150
C               N BY N.                                                 SIMQ 160
C           B - VECTOR OF ORIGINAL CONSTANTS (LENGTH N). THESE ARE      SIMQ 170
C               REPLACED BY FINAL SOLUTION VALUES, VECTOR X.            SIMQ 180
C           N - NUMBER OF EQUATIONS AND VARIABLES. N MUST BE .GT. ONE.  SIMQ 190
C           KS - OUTPUT DIGIT                                           SIMQ 200
C                0 FOR A NORMAL SOLUTION                                SIMQ 210
C                1 FOR A SINGULAR SET OF EQUATIONS                      SIMQ 220
C                                                                       SIMQ 230
C        REMARKS                                                        SIMQ 240
C           MATRIX A MUST BE GENERAL.                                   SIMQ 250
C           IF MATRIX IS SINGULAR , SOLUTION VALUES ARE MEANINGLESS.    SIMQ 260
C           AN ALTERNATIVE SOLUTION MAY BE OBTAINED BY USING MATRIX     SIMQ 270
C           INVERSION (MINV) AND MATRIX PRODUCT (GMPRD).                SIMQ 280
C                                                                       SIMQ 290
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SIMQ 300
C           NONE                                                        SIMQ 310
C                                                                       SIMQ 320
C        METHOD                                                         SIMQ 330
C           METHOD OF SOLUTION IS BY ELIMINATION USING LARGEST PIVOTAL  SIMQ 340
C           DIVISOR. EACH STAGE OF ELIMINATION CONSISTS OF INTERCHANGINGSIMQ 350
C           ROWS WHEN NECESSARY TO AVOID DIVISION BY ZERO OR SMALL      SIMQ 360
C           ELEMENTS.                                                   SIMQ 370
C           THE FORWARD SOLUTION TO OBTAIN VARIABLE N IS DONE IN        SIMQ 380
C           N STAGES. THE BACK SOLUTION FOR THE OTHER VARIABLES IS      SIMQ 390
C           CALCULATED BY SUCCESSIVE SUBSTITUTIONS. FINAL SOLUTION      SIMQ 400
C           VALUES ARE DEVELOPED IN VECTOR B, WITH VARIABLE 1 IN B(1),  SIMQ 410
C           VARIABLE 2 IN B(2),........, VARIABLE N IN B(N).            SIMQ 420
C           IF NO PIVOT CAN BE FOUND EXCEEDING A TOLERANCE OF 0.0,      SIMQ 430
C           THE MATRIX IS CONSIDERED SINGULAR AND KS IS SET TO 1. THIS  SIMQ 440
C           TOLERANCE CAN BE MODIFIED BY REPLACING THE FIRST STATEMENT. SIMQ 450
C                                                                       SIMQ 460
C     ..................................................................SIMQ 470
C                                                                       SIMQ 480
      SUBROUTINE SIMQ(A,B,N,KS)                                         SIMQ 490
      DIMENSION A(1),B(1)                                               SIMQ 500
C                                                                       SIMQ 510
C        FORWARD SOLUTION                                               SIMQ 520
C                                                                       SIMQ 530
      TOL=0.0                                                           SIMQ 540
      KS=0                                                              SIMQ 550
      JJ=-N                                                             SIMQ 560
      DO 65 J=1,N                                                       SIMQ 570
      JY=J+1                                                            SIMQ 580
      JJ=JJ+N+1                                                         SIMQ 590
      BIGA=0                                                            SIMQ 600
      IT=JJ-J                                                           SIMQ 610
      DO 30 I=J,N                                                       SIMQ 620
C                                                                       SIMQ 630
C        SEARCH FOR MAXIMUM COEFFICIENT IN COLUMN                       SIMQ 640
C                                                                       SIMQ 650
      IJ=IT+I                                                           SIMQ 660
      IF(ABS(BIGA)-ABS(A(IJ))) 20,30,30                                 SIMQ 670
   20 BIGA=A(IJ)                                                        SIMQ 680
      IMAX=I                                                            SIMQ 690
   30 CONTINUE                                                          SIMQ 700
C                                                                       SIMQ 710
C        TEST FOR PIVOT LESS THAN TOLERANCE (SINGULAR MATRIX)           SIMQ 720
C                                                                       SIMQ 730
      IF(ABS(BIGA)-TOL) 35,35,40                                        SIMQ 740
   35 KS=1                                                              SIMQ 750
      RETURN                                                            SIMQ 760
C                                                                       SIMQ 770
C        INTERCHANGE ROWS IF NECESSARY                                  SIMQ 780
C                                                                       SIMQ 790
   40 I1=J+N*(J-2)                                                      SIMQ 800
      IT=IMAX-J                                                         SIMQ 810
      DO 50 K=J,N                                                       SIMQ 820
      I1=I1+N                                                           SIMQ 830
      I2=I1+IT                                                          SIMQ 840
      SAVE=A(I1)                                                        SIMQ 850
      A(I1)=A(I2)                                                       SIMQ 860
      A(I2)=SAVE                                                        SIMQ 870
C                                                                       SIMQ 880
C        DIVIDE EQUATION BY LEADING COEFFICIENT                         SIMQ 890
C                                                                       SIMQ 900
   50 A(I1)=A(I1)/BIGA                                                  SIMQ 910
      SAVE=B(IMAX)                                                      SIMQ 920
      B(IMAX)=B(J)                                                      SIMQ 930
      B(J)=SAVE/BIGA                                                    SIMQ 940
C                                                                       SIMQ 950
C        ELIMINATE NEXT VARIABLE                                        SIMQ 960
C                                                                       SIMQ 970
      IF(J-N) 55,70,55                                                  SIMQ 980
   55 IQS=N*(J-1)                                                       SIMQ 990
      DO 65 IX=JY,N                                                     SIMQ1000
      IXJ=IQS+IX                                                        SIMQ1010
      IT=J-IX                                                           SIMQ1020
      DO 60 JX=JY,N                                                     SIMQ1030
      IXJX=N*(JX-1)+IX                                                  SIMQ1040
      JJX=IXJX+IT                                                       SIMQ1050
   60 A(IXJX)=A(IXJX)-(A(IXJ)*A(JJX))                                   SIMQ1060
   65 B(IX)=B(IX)-(B(J)*A(IXJ))                                         SIMQ1070
C                                                                       SIMQ1080
C        BACK SOLUTION                                                  SIMQ1090
C                                                                       SIMQ1100
   70 NY=N-1                                                            SIMQ1110
      IT=N*N                                                            SIMQ1120
      DO 80 J=1,NY                                                      SIMQ1130
      IA=IT-J                                                           SIMQ1140
      IB=N-J                                                            SIMQ1150
      IC=N                                                              SIMQ1160
      DO 80 K=1,J                                                       SIMQ1170
      B(IB)=B(IB)-A(IA)*B(IC)                                           SIMQ1180
      IA=IA-N                                                           SIMQ1190
   80 IC=IC-1                                                           SIMQ1200
      RETURN                                                            SIMQ1210
      END                                                               SIMQ1220
