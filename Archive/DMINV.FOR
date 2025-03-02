C                                                                       DMNV0010
C     ..................................................................DMNV0020
C                                                                       DMNV0030
C        SUBROUTINE DMINV          NOTE  THIS IS SSP(MINV)+ CHANGES     DMNV0040
C                                        TO SEQ..   40 THRU 80          DMNV0050
C        PURPOSE                                   100, 330, 420        DMNV0060
C           INVERT A MATRIX                        680                  DMNV0070
C                                        SEE TEXT  380 THRU 500         DMNV0080
C        USAGE                                                          DMNV0090
C           CALL DMINV(A,N,D,L,M)                                       DMNV0100
C                                                                       DMNV0110
C        DESCRIPTION OF PARAMETERS                                      DMNV0120
C           A - INPUT MATRIX, DESTROYED IN COMPUTATION AND REPLACED BY  DMNV0130
C               RESULTANT INVERSE.                                      DMNV0140
C           N - ORDER OF MATRIX A                                       DMNV0150
C           D - RESULTANT DETERMINANT                                   DMNV0160
C           L - WORK VECTOR OF LENGTH N                                 DMNV0170
C           M - WORK VECTOR OF LENGTH N                                 DMNV0180
C                                                                       DMNV0190
C        REMARKS                                                        DMNV0200
C           MATRIX A MUST BE A GENERAL MATRIX                           DMNV0210
C                                                                       DMNV0220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DMNV0230
C           NONE                                                        DMNV0240
C                                                                       DMNV0250
C        METHOD                                                         DMNV0260
C           THE STANDARD GAUSS-JORDAN METHOD IS USED. THE DETERMINANT   DMNV0270
C           IS ALSO CALCULATED. A DETERMINANT OF ZERO INDICATES THAT    DMNV0280
C           THE MATRIX IS SINGULAR.                                     DMNV0290
C                                                                       DMNV0300
C     ..................................................................DMNV0310
C                                                                       DMNV0320
      SUBROUTINE DMINV(A,N,D,L,M)                                       DMNV0330
      DIMENSION A(1),L(1),M(1)                                          DMNV0340
C                                                                       DMNV0350
C        ...............................................................DMNV0360
C                                                                       DMNV0370
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  DMNV0380
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      DMNV0390
C        STATEMENT WHICH FOLLOWS.                                       DMNV0400
C                                                                       DMNV0410
      DOUBLE PRECISION A,D,BIGA,HOLD                                    DMNV0420
C                                                                       DMNV0430
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    DMNV0440
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      DMNV0450
C        ROUTINE.                                                       DMNV0460
C                                                                       DMNV0470
C        THE DOUBLE PRECISION VERSION OF THIS SUBROUTINE MUST ALSO      DMNV0480
C        CONTAIN DOUBLE PRECISION FORTRAN FUNCTIONS.  ABS IN STATEMENT  DMNV0490
C        10 MUST BE CHANGED TO DABS.                                    DMNV0500
C                                                                       DMNV0510
C        ...............................................................DMNV0520
C                                                                       DMNV0530
C        SEARCH FOR LARGEST ELEMENT                                     DMNV0540
C                                                                       DMNV0550
      D=1.0                                                             DMNV0560
      NK=-N                                                             DMNV0570
      DO 80 K=1,N                                                       DMNV0580
      NK=NK+N                                                           DMNV0590
      L(K)=K                                                            DMNV0600
      M(K)=K                                                            DMNV0610
      KK=NK+K                                                           DMNV0620
      BIGA=A(KK)                                                        DMNV0630
      DO 20 J=K,N                                                       DMNV0640
      IZ=N*(J-1)                                                        DMNV0650
      DO 20 I=K,N                                                       DMNV0660
      IJ=IZ+I                                                           DMNV0670
   10 IF(DABS(BIGA)-DABS(A(IJ))) 15,20,20                               DMNV0680
   15 BIGA=A(IJ)                                                        DMNV0690
      L(K)=I                                                            DMNV0700
      M(K)=J                                                            DMNV0710
   20 CONTINUE                                                          DMNV0720
C                                                                       DMNV0730
C        INTERCHANGE ROWS                                               DMNV0740
C                                                                       DMNV0750
      J=L(K)                                                            DMNV0760
      IF(J-K) 35,35,25                                                  DMNV0770
   25 KI=K-N                                                            DMNV0780
      DO 30 I=1,N                                                       DMNV0790
      KI=KI+N                                                           DMNV0800
      HOLD=-A(KI)                                                       DMNV0810
      JI=KI-K+J                                                         DMNV0820
      A(KI)=A(JI)                                                       DMNV0830
   30 A(JI) =HOLD                                                       DMNV0840
C                                                                       DMNV0850
C        INTERCHANGE COLUMNS                                            DMNV0860
C                                                                       DMNV0870
   35 I=M(K)                                                            DMNV0880
      IF(I-K) 45,45,38                                                  DMNV0890
   38 JP=N*(I-1)                                                        DMNV0900
      DO 40 J=1,N                                                       DMNV0910
      JK=NK+J                                                           DMNV0920
      JI=JP+J                                                           DMNV0930
      HOLD=-A(JK)                                                       DMNV0940
      A(JK)=A(JI)                                                       DMNV0950
   40 A(JI) =HOLD                                                       DMNV0960
C                                                                       DMNV0970
C        DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS        DMNV0980
C        CONTAINED IN BIGA)                                             DMNV0990
C                                                                       DMNV1000
   45 IF(BIGA) 48,46,48                                                 DMNV1010
   46 D=0.0                                                             DMNV1020
      RETURN                                                            DMNV1030
   48 DO 55 I=1,N                                                       DMNV1040
      IF(I-K) 50,55,50                                                  DMNV1050
   50 IK=NK+I                                                           DMNV1060
      A(IK)=A(IK)/(-BIGA)                                               DMNV1070
   55 CONTINUE                                                          DMNV1080
C                                                                       DMNV1090
C        REDUCE MATRIX                                                  DMNV1100
C                                                                       DMNV1110
      DO 65 I=1,N                                                       DMNV1120
      IK=NK+I                                                           DMNV1130
      HOLD=A(IK)                                                        DMNV1140
      IJ=I-N                                                            DMNV1150
      DO 65 J=1,N                                                       DMNV1160
      IJ=IJ+N                                                           DMNV1170
      IF(I-K) 60,65,60                                                  DMNV1180
   60 IF(J-K) 62,65,62                                                  DMNV1190
   62 KJ=IJ-I+K                                                         DMNV1200
      A(IJ)=HOLD*A(KJ)+A(IJ)                                            DMNV1210
   65 CONTINUE                                                          DMNV1220
C                                                                       DMNV1230
C        DIVIDE ROW BY PIVOT                                            DMNV1240
C                                                                       DMNV1250
      KJ=K-N                                                            DMNV1260
      DO 75 J=1,N                                                       DMNV1270
      KJ=KJ+N                                                           DMNV1280
      IF(J-K) 70,75,70                                                  DMNV1290
   70 A(KJ)=A(KJ)/BIGA                                                  DMNV1300
   75 CONTINUE                                                          DMNV1310
C                                                                       DMNV1320
C        PRODUCT OF PIVOTS                                              DMNV1330
C                                                                       DMNV1340
      D=D*BIGA                                                          DMNV1350
C                                                                       DMNV1360
C        REPLACE PIVOT BY RECIPROCAL                                    DMNV1370
C                                                                       DMNV1380
      A(KK)=1.0/BIGA                                                    DMNV1390
   80 CONTINUE                                                          DMNV1400
C                                                                       DMNV1410
C        FINAL ROW AND COLUMN INTERCHANGE                               DMNV1420
C                                                                       DMNV1430
      K=N                                                               DMNV1440
  100 K=(K-1)                                                           DMNV1450
      IF(K) 150,150,105                                                 DMNV1460
  105 I=L(K)                                                            DMNV1470
      IF(I-K) 120,120,108                                               DMNV1480
  108 JQ=N*(K-1)                                                        DMNV1490
      JR=N*(I-1)                                                        DMNV1500
      DO 110 J=1,N                                                      DMNV1510
      JK=JQ+J                                                           DMNV1520
      HOLD=A(JK)                                                        DMNV1530
      JI=JR+J                                                           DMNV1540
      A(JK)=-A(JI)                                                      DMNV1550
  110 A(JI) =HOLD                                                       DMNV1560
  120 J=M(K)                                                            DMNV1570
      IF(J-K) 100,100,125                                               DMNV1580
  125 KI=K-N                                                            DMNV1590
      DO 130 I=1,N                                                      DMNV1600
      KI=KI+N                                                           DMNV1610
      HOLD=A(KI)                                                        DMNV1620
      JI=KI-K+J                                                         DMNV1630
      A(KI)=-A(JI)                                                      DMNV1640
  130 A(JI) =HOLD                                                       DMNV1650
      GO TO 100                                                         DMNV1660
  150 RETURN                                                            DMNV1670
      END                                                               DMNV1680
