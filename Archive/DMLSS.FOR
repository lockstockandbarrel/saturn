C                                                                       DMLS  10
C     ..................................................................DMLS  20
C                                                                       DMLS  30
C        SUBROUTINE DMLSS                                               DMLS  40
C                                                                       DMLS  50
C        PURPOSE                                                        DMLS  60
C           SUBROUTINE DMLSS IS THE SECOND STEP IN THE PROCEDURE FOR    DMLS  70
C           CALCULATING THE LEAST SQUARES SOLUTION OF MINIMAL LENGTH    DMLS  80
C           OF A SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS WITH SYMMETRIC DMLS  90
C           POSITIVE SEMI-DEFINITE COEFFICIENT MATRIX.                  DMLS 100
C                                                                       DMLS 110
C        USAGE                                                          DMLS 120
C           CALL DMLSS(A,N,IRANK,TRAC,INC,RHS,IER)                      DMLS 130
C                                                                       DMLS 140
C        DESCRIPTION OF PARAMETERS                                      DMLS 150
C           A     - COEFFICIENT MATRIX IN FACTORED FORM AS GENERATED    DMLS 160
C                   BY SUBROUTINE MFSS FROM INITIALLY GIVEN SYMMETRIC   DMLS 170
C                   COEFFICIENT MATRIX A STORED IN N*(N+1)/2 LOCATIONS  DMLS 180
C                   A REMAINS UNCHANGED                                 DMLS 190
C                   A MUST BE OF DOUBLE PRECISION                       DMLS 200
C           N     - DIMENSION OF COEFFICIENT MATRIX                     DMLS 210
C           IRANK - RANK OF COEFFICIENT MATRIX, CALCULATED BY MEANS OF  DMLS 220
C                   SUBROUTINE DMFSS                                    DMLS 230
C           TRAC  - VECTOR OF DIMENSION N CONTAINING THE                DMLS 240
C                   SUBSCRIPTS OF PIVOT ROWS AND COLUMNS, I.E. THE      DMLS 250
C                   PRODUCT REPRESENTATION IN TRANSPOSITIONS OF THE     DMLS 260
C                   PERMUTATION WHICH WAS APPLIED TO ROWS AND COLUMNS   DMLS 270
C                   OF A IN THE FACTORIZATION PROCESS                   DMLS 280
C                   TRAC IS A RESULTANT ARRAY OF SUBROUTINE MFSS        DMLS 290
C                   TRAC MUST BE OF DOUBLE PRECISION                    DMLS 300
C           INC   - INPUT VARIABLE WHICH SHOULD CONTAIN THE VALUE ZERO  DMLS 310
C                   IF THE SYSTEM OF SIMULTANEOUS EQUATIONS IS KNOWN    DMLS 320
C                   TO BE COMPATIBLE AND A NONZERO VALUE OTHERWISE      DMLS 330
C           RHS   - VECTOR OF DIMENSION N CONTAINING THE RIGHT HAND SIDEDMLS 340
C                   ON RETURN RHS CONTAINS THE MINIMAL LENGTH SOLUTION  DMLS 350
C                   RHS MUST BE OF DOUBLE PRECISION                     DMLS 360
C           IER   - RESULTANT ERROR PARAMETER                           DMLS 370
C                   IER = 0 MEANS NO ERRORS                             DMLS 380
C                   IER =-1 MEANS N AND/OR IRANK IS NOT POSITIVE AND/OR DMLS 390
C                           IRANK IS GREATER THAN N                     DMLS 400
C                   IER = 1 MEANS THE FACTORIZATION CONTAINED IN A HAS  DMLS 410
C                           ZERO DIVISORS AND/OR TRAC CONTAINS          DMLS 420
C                           VALUES OUTSIDE THE FEASIBLE RANGE 1 UP TO N DMLS 430
C                                                                       DMLS 440
C        REMARKS                                                        DMLS 450
C           THE MINIMAL LENGTH SOLUTION IS PRODUCED IN THE STORAGE      DMLS 460
C           LOCATIONS OCCUPIED BY THE RIGHT HAND SIDE.                  DMLS 470
C           SUBROUTINE DMLSS DOES TAKE CARE OF THE PERMUTATION          DMLS 480
C           WHICH WAS APPLIED TO ROWS AND COLUMNS OF A.                 DMLS 490
C           OPERATION IS BYPASSED IN CASE OF A NON POSITIVE VALUE       DMLS 500
C           OF IRANK                                                    DMLS 510
C                                                                       DMLS 520
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DMLS 530
C           NONE                                                        DMLS 540
C                                                                       DMLS 550
C        METHOD                                                         DMLS 560
C           LET T, U, TU BE THE COMPONENTS OF THE FACTORIZATION OF A,   DMLS 570
C           AND LET THE RIGHT HAND SIDE BE PARTITIONED INTO A FIRST     DMLS 580
C           PART X1 OF DIMENSION IRANK AND A SECOND PART X2 OF DIMENSIONDMLS 590
C           N-IRANK. THEN THE FOLLOWING OPERATIONS ARE APPLIED IN       DMLS 600
C           SEQUENCE                                                    DMLS 610
C           (1) INTERCHANGE RIGHT HAND SIDE                             DMLS 620
C           (2) X1 = X1 + U * X2                                        DMLS 630
C           (3) X2 =-TRANSPOSE(U) * X1                                  DMLS 640
C           (4) X2 = INVERSE(TU) * INVERSE(TRANSPOSE(TU)) * X2          DMLS 650
C           (5) X1 = X1 + U * X2                                        DMLS 660
C           (6) X1 = INVERSE(T) * INVERSE(TRANSPOSE(T)) * X1            DMLS 670
C           (7) X2 =-TRANSPOSE(U) * X1                                  DMLS 680
C           (8) X2 = INVERSE(TU) * INVERSE(TRANSPOSE(TU)) * X2          DMLS 690
C           (9) X1 = X1 + U * X2                                        DMLS 700
C           (10)X2 = TRANSPOSE(U) * X1                                  DMLS 710
C           (11) REINTERCHANGE CALCULATED SOLUTION                      DMLS 720
C           IF THE SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS IS SPECIFIED DMLS 730
C           TO BE COMPATIBLE THEN STEPS (2), (3), (4) AND (5) ARE       DMLS 740
C           CANCELLED.                                                  DMLS 750
C           IF THE COEFFICIENT MATRIX HAS RANK N, THEN THE ONLY STEPS   DMLS 760
C           PERFORMED ARE (1), (6) AND (11).                            DMLS 770
C                                                                       DMLS 780
C     ..................................................................DMLS 790
C                                                                       DMLS 800
      SUBROUTINE DMLSS(A,N,IRANK,TRAC,INC,RHS,IER)                      DMLS 810
C                                                                       DMLS 820
C                                                                       DMLS 830
C        DIMENSIONED DUMMY VARIABLES                                    DMLS 840
      DIMENSION A(1),TRAC(1),RHS(1)                                     DMLS 850
      DOUBLE PRECISION SUM,A,RHS,TRAC,HOLD                              DMLS 860
C                                                                       DMLS 870
C        TEST OF SPECIFIED DIMENSIONS                                   DMLS 880
      IDEF=N-IRANK                                                      DMLS 890
      IF(N)33,33,1                                                      DMLS 900
    1 IF(IRANK)33,33,2                                                  DMLS 910
    2 IF(IDEF)33,3,3                                                    DMLS 920
C                                                                       DMLS 930
C        CALCULATE AUXILIARY VALUES                                     DMLS 940
    3 ITE=IRANK*(IRANK+1)/2                                             DMLS 950
      IX2=IRANK+1                                                       DMLS 960
      NP1=N+1                                                           DMLS 970
      IER=0                                                             DMLS 980
C                                                                       DMLS 990
C        INTERCHANGE RIGHT HAND SIDE                                    DMLS1000
      JJ=1                                                              DMLS1010
      II=1                                                              DMLS1020
    4 DO 6 I=1,N                                                        DMLS1030
      J=TRAC(II)                                                        DMLS1040
      IF(J)31,31,5                                                      DMLS1050
    5 HOLD=RHS(II)                                                      DMLS1060
      RHS(II)=RHS(J)                                                    DMLS1070
      RHS(J)=HOLD                                                       DMLS1080
    6 II=II+JJ                                                          DMLS1090
      IF(JJ)32,7,7                                                      DMLS1100
C                                                                       DMLS1110
C        PERFORM STEP 2 IF NECESSARY                                    DMLS1120
    7 ISW=1                                                             DMLS1130
      IF(INC*IDEF)8,28,8                                                DMLS1140
C                                                                       DMLS1150
C        CALCULATE X1 = X1 + U * X2                                     DMLS1160
    8 ISTA=ITE                                                          DMLS1170
      DO 10 I=1,IRANK                                                   DMLS1180
      ISTA=ISTA+1                                                       DMLS1190
      JJ=ISTA                                                           DMLS1200
      SUM=0.D0                                                          DMLS1210
      DO 9 J=IX2,N                                                      DMLS1220
      SUM=SUM+A(JJ)*RHS(J)                                              DMLS1230
    9 JJ=JJ+J                                                           DMLS1240
   10 RHS(I)=RHS(I)+SUM                                                 DMLS1250
      GOTO(11,28,11),ISW                                                DMLS1260
C                                                                       DMLS1270
C        CALCULATE X2 = TRANSPOSE(U) * X1                               DMLS1280
   11 ISTA=ITE                                                          DMLS1290
      DO 15 I=IX2,N                                                     DMLS1300
      JJ=ISTA                                                           DMLS1310
      SUM=0.D0                                                          DMLS1320
      DO 12 J=1,IRANK                                                   DMLS1330
      JJ=JJ+1                                                           DMLS1340
   12 SUM=SUM+A(JJ)*RHS(J)                                              DMLS1350
      GOTO(13,13,14),ISW                                                DMLS1360
   13 SUM=-SUM                                                          DMLS1370
   14 RHS(I)=SUM                                                        DMLS1380
   15 ISTA=ISTA+I                                                       DMLS1390
      GOTO(16,29,30),ISW                                                DMLS1400
C                                                                       DMLS1410
C        INITIALIZE STEP (4) OR STEP (8)                                DMLS1420
   16 ISTA=IX2                                                          DMLS1430
      IEND=N                                                            DMLS1440
      JJ=ITE+ISTA                                                       DMLS1450
C                                                                       DMLS1460
C        DIVISION OF X1 BY TRANSPOSE OF TRIANGULAR MATRIX               DMLS1470
   17 SUM=0.D0                                                          DMLS1480
      DO 20 I=ISTA,IEND                                                 DMLS1490
      IF(A(JJ))18,31,18                                                 DMLS1500
   18 RHS(I)=(RHS(I)-SUM)/A(JJ)                                         DMLS1510
      IF(I-IEND)19,21,21                                                DMLS1520
   19 JJ=JJ+ISTA                                                        DMLS1530
      SUM=0.D0                                                          DMLS1540
      DO 20 J=ISTA,I                                                    DMLS1550
      SUM=SUM+A(JJ)*RHS(J)                                              DMLS1560
   20 JJ=JJ+1                                                           DMLS1570
C                                                                       DMLS1580
C        DIVISION OF X1 BY TRIANGULAR MATRIX                            DMLS1590
   21 SUM=0.D0                                                          DMLS1600
      II=IEND                                                           DMLS1610
      DO 24 I=ISTA,IEND                                                 DMLS1620
      RHS(II)=(RHS(II)-SUM)/A(JJ)                                       DMLS1630
      IF(II-ISTA)25,25,22                                               DMLS1640
   22 KK=JJ-1                                                           DMLS1650
      SUM=0.D0                                                          DMLS1660
      DO 23 J=II,IEND                                                   DMLS1670
      SUM=SUM+A(KK)*RHS(J)                                              DMLS1680
   23 KK=KK+J                                                           DMLS1690
      JJ=JJ-II                                                          DMLS1700
   24 II=II-1                                                           DMLS1710
   25 IF(IDEF)26,30,26                                                  DMLS1720
   26 GOTO(27,11,8),ISW                                                 DMLS1730
C                                                                       DMLS1740
C        PERFORM STEP (5)                                               DMLS1750
   27 ISW=2                                                             DMLS1760
      GOTO 8                                                            DMLS1770
C                                                                       DMLS1780
C        PERFORM STEP (6)                                               DMLS1790
   28 ISTA=1                                                            DMLS1800
      IEND=IRANK                                                        DMLS1810
      JJ=1                                                              DMLS1820
      ISW=2                                                             DMLS1830
      GOTO 17                                                           DMLS1840
C                                                                       DMLS1850
C        PERFORM STEP (8)                                               DMLS1860
   29 ISW=3                                                             DMLS1870
      GOTO 16                                                           DMLS1880
C                                                                       DMLS1890
C        REINTERCHANGE CALCULATED SOLUTION                              DMLS1900
   30 II=N                                                              DMLS1910
      JJ=-1                                                             DMLS1920
      GOTO 4                                                            DMLS1930
C                                                                       DMLS1940
C        ERROR RETURN IN CASE OF ZERO DIVISOR                           DMLS1950
   31 IER=1                                                             DMLS1960
   32 RETURN                                                            DMLS1970
C                                                                       DMLS1980
C        ERROR RETURN IN CASE OF ILLEGAL DIMENSION                      DMLS1990
   33 IER=-1                                                            DMLS2000
      RETURN                                                            DMLS2010
      END                                                               DMLS2020
