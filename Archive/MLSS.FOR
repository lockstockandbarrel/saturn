C                                                                       MLSS  10
C     ..................................................................MLSS  20
C                                                                       MLSS  30
C        SUBROUTINE MLSS                                                MLSS  40
C                                                                       MLSS  50
C        PURPOSE                                                        MLSS  60
C           SUBROUTINE MLSS IS THE SECOND STEP IN THE PROCEDURE FOR     MLSS  70
C           CALCULATING THE LEAST SQUARES SOLUTION OF MINIMAL LENGTH    MLSS  80
C           OF A SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS WITH SYMMETRIC MLSS  90
C           POSITIVE SEMI-DEFINITE COEFFICIENT MATRIX.                  MLSS 100
C                                                                       MLSS 110
C        USAGE                                                          MLSS 120
C           CALL MLSS(A,N,IRANK,TRAC,INC,RHS,IER)                       MLSS 130
C                                                                       MLSS 140
C        DESCRIPTION OF PARAMETERS                                      MLSS 150
C           A     - COEFFICIENT MATRIX IN FACTORED FORM AS GENERATED    MLSS 160
C                   BY SUBROUTINE MFSS FROM INITIALLY GIVEN SYMMETRIC   MLSS 170
C                   COEFFICIENT MATRIX A STORED IN N*(N+1)/2 LOCATIONS  MLSS 180
C                   A REMAINS UNCHANGED                                 MLSS 190
C           N     - DIMENSION OF COEFFICIENT MATRIX                     MLSS 200
C           IRANK - RANK OF COEFFICIENT MATRIX, CALCULATED BY MEANS OF  MLSS 210
C                   SUBROUTINE MFSS                                     MLSS 220
C           TRAC  - VECTOR OF DIMENSION N CONTAINING THE                MLSS 230
C                   SUBSCRIPTS OF PIVOT ROWS AND COLUMNS, I.E. THE      MLSS 240
C                   PRODUCT REPRESENTATION IN TRANSPOSITIONS OF THE     MLSS 250
C                   PERMUTATION WHICH WAS APPLIED TO ROWS AND COLUMNS   MLSS 260
C                   OF A IN THE FACTORIZATION PROCESS                   MLSS 270
C                   TRAC IS A RESULTANT ARRAY OF SUBROUTINE MFSS        MLSS 280
C           INC   - INPUT VARIABLE WHICH SHOULD CONTAIN THE VALUE ZERO  MLSS 290
C                   IF THE SYSTEM OF SIMULTANEOUS EQUATIONS IS KNOWN    MLSS 300
C                   TO BE COMPATIBLE AND A NONZERO VALUE OTHERWISE      MLSS 310
C           RHS   - VECTOR OF DIMENSION N CONTAINING THE RIGHT HAND SIDEMLSS 320
C                   ON RETURN RHS CONTAINS THE MINIMAL LENGTH SOLUTION  MLSS 330
C           IER   - RESULTANT ERROR PARAMETER                           MLSS 340
C                   IER = 0 MEANS NO ERRORS                             MLSS 350
C                   IER =-1 MEANS N AND/OR IRANK IS NOT POSITIVE AND/OR MLSS 360
C                           IRANK IS GREATER THAN N                     MLSS 370
C                   IER = 1 MEANS THE FACTORIZATION CONTAINED IN A HAS  MLSS 380
C                           ZERO DIVISORS AND/OR TRAC CONTAINS          MLSS 390
C                           VALUES OUTSIDE THE FEASIBLE RANGE 1 UP TO N MLSS 400
C                                                                       MLSS 410
C        REMARKS                                                        MLSS 420
C           THE MINIMAL LENGTH SOLUTION IS PRODUCED IN THE STORAGE      MLSS 430
C           LOCATIONS OCCUPIED BY THE RIGHT HAND SIDE.                  MLSS 440
C           SUBROUTINE MLSS DOES TAKE CARE OF THE PERMUTATION           MLSS 450
C           WHICH WAS APPLIED TO ROWS AND COLUMNS OF A.                 MLSS 460
C           OPERATION IS BYPASSED IN CASE OF A NON POSITIVE VALUE       MLSS 470
C           OF IRANK                                                    MLSS 480
C                                                                       MLSS 490
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MLSS 500
C           NONE                                                        MLSS 510
C                                                                       MLSS 520
C        METHOD                                                         MLSS 530
C           LET T, U, TU BE THE COMPONENTS OF THE FACTORIZATION OF A,   MLSS 540
C           AND LET THE RIGHT HAND SIDE BE PARTITIONED INTO A FIRST     MLSS 550
C           PART X1 OF DIMENSION IRANK AND A SECOND PART X2 OF DIMENSIONMLSS 560
C           N-IRANK. THEN THE FOLLOWING OPERATIONS ARE APPLIED IN       MLSS 570
C           SEQUENCE                                                    MLSS 580
C           (1) INTERCHANGE RIGHT HAND SIDE                             MLSS 590
C           (2) X1 = X1 + U * X2                                        MLSS 600
C           (3) X2 =-TRANSPOSE(U) * X1                                  MLSS 610
C           (4) X2 = INVERSE(TU) * INVERSE(TRANSPOSE(TU)) * X2          MLSS 620
C           (5) X1 = X1 + U * X2                                        MLSS 630
C           (6) X1 = INVERSE(T) * INVERSE(TRANSPOSE(T)) * X1            MLSS 640
C           (7) X2 =-TRANSPOSE(U) * X1                                  MLSS 650
C           (8) X2 = INVERSE(TU) * INVERSE(TRANSPOSE(TU)) * X2          MLSS 660
C           (9) X1 = X1 + U * X2                                        MLSS 670
C           (10)X2 = TRANSPOSE(U) * X1                                  MLSS 680
C           (11) REINTERCHANGE CALCULATED SOLUTION                      MLSS 690
C           IF THE SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS IS SPECIFIED MLSS 700
C           TO BE COMPATIBLE THEN STEPS (2), (3), (4) AND (5) ARE       MLSS 710
C           CANCELLED.                                                  MLSS 720
C           IF THE COEFFICIENT MATRIX HAS RANK N, THEN THE ONLY STEPS   MLSS 730
C           PERFORMED ARE (1), (6) AND (11).                            MLSS 740
C                                                                       MLSS 750
C     ..................................................................MLSS 760
C                                                                       MLSS 770
      SUBROUTINE MLSS(A,N,IRANK,TRAC,INC,RHS,IER)                       MLSS 780
C                                                                       MLSS 790
C                                                                       MLSS 800
C        DIMENSIONED DUMMY VARIABLES                                    MLSS 810
      DIMENSION A(1),TRAC(1),RHS(1)                                     MLSS 820
      DOUBLE PRECISION SUM                                              MLSS 830
C                                                                       MLSS 840
C        TEST OF SPECIFIED DIMENSIONS                                   MLSS 850
      IDEF=N-IRANK                                                      MLSS 860
      IF(N)33,33,1                                                      MLSS 870
    1 IF(IRANK)33,33,2                                                  MLSS 880
    2 IF(IDEF)33,3,3                                                    MLSS 890
C                                                                       MLSS 900
C        CALCULATE AUXILIARY VALUES                                     MLSS 910
    3 ITE=IRANK*(IRANK+1)/2                                             MLSS 920
      IX2=IRANK+1                                                       MLSS 930
      NP1=N+1                                                           MLSS 940
      IER=0                                                             MLSS 950
C                                                                       MLSS 960
C        INTERCHANGE RIGHT HAND SIDE                                    MLSS 970
      JJ=1                                                              MLSS 980
      II=1                                                              MLSS 990
    4 DO 6 I=1,N                                                        MLSS1000
      J=TRAC(II)                                                        MLSS1010
      IF(J)31,31,5                                                      MLSS1020
    5 HOLD=RHS(II)                                                      MLSS1030
      RHS(II)=RHS(J)                                                    MLSS1040
      RHS(J)=HOLD                                                       MLSS1050
    6 II=II+JJ                                                          MLSS1060
      IF(JJ)32,7,7                                                      MLSS1070
C                                                                       MLSS1080
C        PERFORM STEP 2 IF NECESSARY                                    MLSS1090
    7 ISW=1                                                             MLSS1100
      IF(INC*IDEF)8,28,8                                                MLSS1110
C                                                                       MLSS1120
C        CALCULATE X1 = X1 + U * X2                                     MLSS1130
    8 ISTA=ITE                                                          MLSS1140
      DO 10 I=1,IRANK                                                   MLSS1150
      ISTA=ISTA+1                                                       MLSS1160
      JJ=ISTA                                                           MLSS1170
      SUM=0.D0                                                          MLSS1180
      DO 9 J=IX2,N                                                      MLSS1190
      SUM=SUM+A(JJ)*RHS(J)                                              MLSS1200
    9 JJ=JJ+J                                                           MLSS1210
   10 RHS(I)=RHS(I)+SUM                                                 MLSS1220
      GOTO(11,28,11),ISW                                                MLSS1230
C                                                                       MLSS1240
C        CALCULATE X2 = TRANSPOSE(U) * X1                               MLSS1250
   11 ISTA=ITE                                                          MLSS1260
      DO 15 I=IX2,N                                                     MLSS1270
      JJ=ISTA                                                           MLSS1280
      SUM=0.D0                                                          MLSS1290
      DO 12 J=1,IRANK                                                   MLSS1300
      JJ=JJ+1                                                           MLSS1310
   12 SUM=SUM+A(JJ)*RHS(J)                                              MLSS1320
      GOTO(13,13,14),ISW                                                MLSS1330
   13 SUM=-SUM                                                          MLSS1340
   14 RHS(I)=SUM                                                        MLSS1350
   15 ISTA=ISTA+I                                                       MLSS1360
      GOTO(16,29,30),ISW                                                MLSS1370
C                                                                       MLSS1380
C        INITIALIZE STEP (4) OR STEP (8)                                MLSS1390
   16 ISTA=IX2                                                          MLSS1400
      IEND=N                                                            MLSS1410
      JJ=ITE+ISTA                                                       MLSS1420
C                                                                       MLSS1430
C        DIVISION OF X1 BY TRANSPOSE OF TRIANGULAR MATRIX               MLSS1440
   17 SUM=0.D0                                                          MLSS1450
      DO 20 I=ISTA,IEND                                                 MLSS1460
      IF(A(JJ))18,31,18                                                 MLSS1470
   18 RHS(I)=(RHS(I)-SUM)/A(JJ)                                         MLSS1480
      IF(I-IEND)19,21,21                                                MLSS1490
   19 JJ=JJ+ISTA                                                        MLSS1500
      SUM=0.D0                                                          MLSS1510
      DO 20 J=ISTA,I                                                    MLSS1520
      SUM=SUM+A(JJ)*RHS(J)                                              MLSS1530
   20 JJ=JJ+1                                                           MLSS1540
C                                                                       MLSS1550
C        DIVISION OF X1 BY TRIANGULAR MATRIX                            MLSS1560
   21 SUM=0.D0                                                          MLSS1570
      II=IEND                                                           MLSS1580
      DO 24 I=ISTA,IEND                                                 MLSS1590
      RHS(II)=(RHS(II)-SUM)/A(JJ)                                       MLSS1600
      IF(II-ISTA)25,25,22                                               MLSS1610
   22 KK=JJ-1                                                           MLSS1620
      SUM=0.D0                                                          MLSS1630
      DO 23 J=II,IEND                                                   MLSS1640
      SUM=SUM+A(KK)*RHS(J)                                              MLSS1650
   23 KK=KK+J                                                           MLSS1660
      JJ=JJ-II                                                          MLSS1670
   24 II=II-1                                                           MLSS1680
   25 IF(IDEF)26,30,26                                                  MLSS1690
   26 GOTO(27,11,8),ISW                                                 MLSS1700
C                                                                       MLSS1710
C        PERFORM STEP (5)                                               MLSS1720
   27 ISW=2                                                             MLSS1730
      GOTO 8                                                            MLSS1740
C                                                                       MLSS1750
C        PERFORM STEP (6)                                               MLSS1760
   28 ISTA=1                                                            MLSS1770
      IEND=IRANK                                                        MLSS1780
      JJ=1                                                              MLSS1790
      ISW=2                                                             MLSS1800
      GOTO 17                                                           MLSS1810
C                                                                       MLSS1820
C        PERFORM STEP (8)                                               MLSS1830
   29 ISW=3                                                             MLSS1840
      GOTO 16                                                           MLSS1850
C                                                                       MLSS1860
C        REINTERCHANGE CALCULATED SOLUTION                              MLSS1870
   30 II=N                                                              MLSS1880
      JJ=-1                                                             MLSS1890
      GOTO 4                                                            MLSS1900
C                                                                       MLSS1910
C        ERROR RETURN IN CASE OF ZERO DIVISOR                           MLSS1920
   31 IER=1                                                             MLSS1930
   32 RETURN                                                            MLSS1940
C                                                                       MLSS1950
C        ERROR RETURN IN CASE OF ILLEGAL DIMENSION                      MLSS1960
   33 IER=-1                                                            MLSS1970
      RETURN                                                            MLSS1980
      END                                                               MLSS1990
