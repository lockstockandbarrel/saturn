C                                                                       DMTD  10
C     ..................................................................DMTD  20
C                                                                       DMTD  30
C        SUBROUTINE DMTDS                                               DMTD  40
C                                                                       DMTD  50
C        PURPOSE                                                        DMTD  60
C           MULTIPLY A GENERAL MATRIX A ON THE LEFT OR RIGHT BY         DMTD  70
C           INVERSE(T),INVERSE(TRANSPOSE(T)) OR INVERSE(TRANSPOSE(T*T)) DMTD  80
C           THE TRIANGULAR MATRIX T IS STORED COLUMNWISE IN COMPRESSED  DMTD  90
C           FORM, I.E. UPPER TRIANGULAR PART ONLY.                      DMTD 100
C                                                                       DMTD 110
C        USAGE                                                          DMTD 120
C           CALL DMTDS(A,M,N,T,IOP,IER)                                 DMTD 130
C                                                                       DMTD 140
C        DESCRIPTION OF PARAMETERS                                      DMTD 150
C           A     - GIVEN GENERAL MATRIX WITH  M ROWS AND N COLUMNS.    DMTD 160
C                   A MUST BE OF DOUBLE PRECISION                       DMTD 170
C           M     - NUMBER OF ROWS OF MATRIX A                          DMTD 180
C           N     - NUMBER OF COLUMNS OF MATRIX A                       DMTD 190
C           T     - GIVEN TRIANGULAR MATRIX STORED COLUMNWISE UPPER     DMTD 200
C                   TRIANGULAR PART ONLY. ITS NUMBER OF ROWS AND        DMTD 210
C                   COLUMNS K IS IMPLIED BY COMPATIBILITY.              DMTD 220
C                   K = M IF IOP IS POSITIVE,                           DMTD 230
C                   K = N IF IOP IS NEGATIVE.                           DMTD 240
C                   T OCCUPIES K*(K+1)/2 STORAGE POSITIONS.             DMTD 250
C                   T MUST BE OF DOUBLE PRECISION                       DMTD 260
C           IOP   - INPUT VARIABLE FOR SELECTION OF OPERATION           DMTD 270
C                   IOP = 1 - A IS REPLACED BY INVERSE(T)*A             DMTD 280
C                   IOP =-1 - A IS REPLACED BY A*INVERSE(T)             DMTD 290
C                   IOP = 2 - A IS REPLACED BY INVERSE(TRANSPOSE(T))*A  DMTD 300
C                   IOP =-2 - A IS REPLACED BY A*INVERSE(TRANSPOSE(T))  DMTD 310
C                   IOP = 3 - A IS REPLACED BY INVERSE(TRANSPOSE(T)*T)*ADMTD 320
C                   IOP =-3 - A IS REPLACED BY A*INVERSE(TRANSPOSE(T)*T)DMTD 330
C           IER   - RESULTING ERROR PARAMETER                           DMTD 340
C                   IER =-1 MEANS M AND N ARE NOT BOTH POSITIVE         DMTD 350
C                                 AND/OR IOP IS ILLEGAL                 DMTD 360
C                   IER = 0 MEANS OPERATION WAS SUCCESSFUL              DMTD 370
C                   IER = 1 MEANS TRIANGULAR MATRIX T IS SINGULAR       DMTD 380
C                                                                       DMTD 390
C        REMARKS                                                        DMTD 400
C           SUBROUTINE DMTDS MAY BE USED TO CALCULATE THE SOLUTION OF   DMTD 410
C           A SYSTEM OF EQUATIONS WITH SYMMETRIC POSITIVE DEFINITE      DMTD 420
C           COEFFICIENT MATRIX. THE FIRST STEP TOWARDS THE SOLUTION     DMTD 430
C           IS TRIANGULAR FACTORIZATION BY MEANS OF DMFSD, THE SECOND   DMTD 440
C           STEP IS APPLICATION OF DMTDS.                               DMTD 450
C           SUBROUTINES DMFSD AND DMTDS MAY BE USED IN ORDER TO         DMTD 460
C           CACULATE THE PRODUCT TRANSPOSE(A)*INVERSE(B)*A WITH GIVEN   DMTD 470
C           SYMMETRIC POSITIVE DEFINITE B AND GIVEN A IN THREE STEPS    DMTD 480
C           1) TRIANGULAR FACTORIZATION OF B (B=TRANSPOSE(T)*T)         DMTD 490
C           2) MULTIPLICATION OF A ON THE LEFT BY INVERSE(TRANSPOSE(T)) DMTD 500
C              A IS REPLACED BY C=INVERSE(TRANSPOSE(T))*A               DMTD 510
C           3) CALCULATION OF THE RESULT FORMING TRANSPOSE(C)*C         DMTD 520
C                                                                       DMTD 530
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DMTD 540
C           NONE                                                        DMTD 550
C                                                                       DMTD 560
C        METHOD                                                         DMTD 570
C           CALCULATION OF X = INVERSE(T)*A IS DONE USING BACKWARD      DMTD 580
C           SUBSTITUTION TO OBTAIN X FROM T*X = A.                      DMTD 590
C           CALCULATION OF Y = INVERSE(TRANSPOSE(T))*A IS DONE USING    DMTD 600
C           FORWARD SUBSTITUTION TO OBTAIN Y FROM TRANSPOSE(T)*Y = A.   DMTD 610
C           CALCULATION OF Z = INVERSE(TRANSPOSE(T)*T)*A IS DONE        DMTD 620
C           SOLVING FIRST TRANSPOSE(T)*Y = A AND THEN T*Z = Y, IE.      DMTD 630
C           USING THE ABOVE TWO STEPS IN REVERSE ORDER                  DMTD 640
C                                                                       DMTD 650
C     ..................................................................DMTD 660
C                                                                       DMTD 670
      SUBROUTINE DMTDS(A,M,N,T,IOP,IER)                                 DMTD 680
C                                                                       DMTD 690
C                                                                       DMTD 700
      DIMENSION A(1),T(1)                                               DMTD 710
      DOUBLE PRECISION DSUM,A,T                                         DMTD 720
C                                                                       DMTD 730
C        TEST OF DIMENSION                                              DMTD 740
      IF(M)2,2,1                                                        DMTD 750
    1 IF(N)2,2,4                                                        DMTD 760
C                                                                       DMTD 770
C        ERROR RETURN IN CASE OF ILLEGAL DIMENSIONS                     DMTD 780
    2 IER=-1                                                            DMTD 790
      RETURN                                                            DMTD 800
C                                                                       DMTD 810
C        ERROR RETURN IN CASE OF SINGULAR MATRIX T                      DMTD 820
    3 IER=1                                                             DMTD 830
      RETURN                                                            DMTD 840
C                                                                       DMTD 850
C        INITIALIZE DIVISION PROCESS                                    DMTD 860
    4 MN=M*N                                                            DMTD 870
      MM=M*(M+1)/2                                                      DMTD 880
      MM1=M-1                                                           DMTD 890
      IER=0                                                             DMTD 900
      ICS=M                                                             DMTD 910
      IRS=1                                                             DMTD 920
      IMEND=M                                                           DMTD 930
C                                                                       DMTD 940
C        TEST SPECIFIED OPERATION                                       DMTD 950
      IF(IOP)5,2,6                                                      DMTD 960
    5 MM=N*(N+1)/2                                                      DMTD 970
      MM1=N-1                                                           DMTD 980
      IRS=M                                                             DMTD 990
      ICS=1                                                             DMTD1000
      IMEND=MN-M+1                                                      DMTD1010
      MN=M                                                              DMTD1020
    6 IOPE=MOD(IOP+3,3)                                                 DMTD1030
      IF(IABS(IOP)-3)7,7,2                                              DMTD1040
    7 IF(IOPE-1)8,18,8                                                  DMTD1050
C                                                                       DMTD1060
C        INITIALIZE SOLUTION OF TRANSPOSE(T)*X = A                      DMTD1070
    8 MEND=1                                                            DMTD1080
      LLD=IRS                                                           DMTD1090
      MSTA=1                                                            DMTD1100
      MDEL=1                                                            DMTD1110
      MX=1                                                              DMTD1120
      LD=1                                                              DMTD1130
      LX=0                                                              DMTD1140
C                                                                       DMTD1150
C        TEST FOR NONZERO DIAGONAL TERM IN T                            DMTD1160
    9 IF(T(MSTA))10,3,10                                                DMTD1170
   10 DO 11 I=MEND,MN,ICS                                               DMTD1180
   11 A(I)=A(I)/T(MSTA)                                                 DMTD1190
C                                                                       DMTD1200
C        IS M EQUAL 1                                                   DMTD1210
      IF(MM1)2,15,12                                                    DMTD1220
   12 DO 14 J=1,MM1                                                     DMTD1230
      MSTA=MSTA+MDEL                                                    DMTD1240
      MDEL=MDEL+MX                                                      DMTD1250
      DO 14 I=MEND,MN,ICS                                               DMTD1260
      DSUM=0.D0                                                         DMTD1270
      L=MSTA                                                            DMTD1280
      LDX=LD                                                            DMTD1290
      LL=I                                                              DMTD1300
      DO 13 K=1,J                                                       DMTD1310
      DSUM=DSUM-T(L)*A(LL)                                              DMTD1320
      LL=LL+LLD                                                         DMTD1330
      L=L+LDX                                                           DMTD1340
   13 LDX=LDX+LX                                                        DMTD1350
      IF(T(L))14,3,14                                                   DMTD1360
   14 A(LL)=(DSUM+A(LL))/T(L)                                           DMTD1370
C                                                                       DMTD1380
C        TEST END OF OPERATION                                          DMTD1390
   15 IF(IER)16,17,16                                                   DMTD1400
   16 IER=0                                                             DMTD1410
      RETURN                                                            DMTD1420
   17 IF(IOPE)18,18,16                                                  DMTD1430
C                                                                       DMTD1440
C        INITIALIZE SOLUTION OF T*X = A                                 DMTD1450
   18 IER=1                                                             DMTD1460
      MEND=IMEND                                                        DMTD1470
      MN=M*N                                                            DMTD1480
      LLD=-IRS                                                          DMTD1490
      MSTA=MM                                                           DMTD1500
      MDEL=-1                                                           DMTD1510
      MX=0                                                              DMTD1520
      LD=-MM1                                                           DMTD1530
      LX=1                                                              DMTD1540
      GOTO 9                                                            DMTD1550
      END                                                               DMTD1560
