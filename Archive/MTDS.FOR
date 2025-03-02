C                                                                       MTDS  10
C     ..................................................................MTDS  20
C                                                                       MTDS  30
C        SUBROUTINE MTDS                                                MTDS  40
C                                                                       MTDS  50
C        PURPOSE                                                        MTDS  60
C           MULTIPLY A GENERAL MATRIX A ON THE LEFT OR RIGHT BY         MTDS  70
C           INVERSE(T),INVERSE(TRANSPOSE(T)) OR INVERSE(TRANSPOSE(T*T)) MTDS  80
C           THE TRIANGULAR MATRIX T IS STORED COLUMNWISE IN COMPRESSED  MTDS  90
C           FORM, I.E. UPPER TRIANGULAR PART ONLY.                      MTDS 100
C                                                                       MTDS 110
C        USAGE                                                          MTDS 120
C           CALL MTDS(A,M,N,T,IOP,IER)                                  MTDS 130
C                                                                       MTDS 140
C        DESCRIPTION OF PARAMETERS                                      MTDS 150
C           A     - GIVEN GENERAL MATRIX WHITH M ROWS AND N COLUMNS.    MTDS 160
C           M     - NUMBER OF ROWS OF MATRIX A                          MTDS 170
C           N     - NUMBER OF COLUMNS OF MATRIX A                       MTDS 180
C           T     - GIVEN TRIANGULAR MATRIX STORED COLUMNWISE UPPER     MTDS 190
C                   TRIANGULAR PART ONLY. ITS NUMBER OF ROWS AND        MTDS 200
C                   COLUMNS K IS IMPLIED BY COMPATIBILITY.              MTDS 210
C                   K = M IF IOP IS POSITIVE,                           MTDS 220
C                   K = N IF IOP IS NEGATIVE.                           MTDS 230
C                   T OCCUPIES K*(K+1)/2 STORAGE POSITIONS.             MTDS 240
C           IOP   - INPUT VARIABLE FOR SELECTION OF OPERATION           MTDS 250
C                   IOP = 1 - A IS REPLACED BY INVERSE(T)*A             MTDS 260
C                   IOP =-1 - A IS REPLACED BY A*INVERSE(T)             MTDS 270
C                   IOP = 2 - A IS REPLACED BY INVERSE(TRANSPOSE(T))*A  MTDS 280
C                   IOP =-2 - A IS REPLACED BY A*INVERSE(TRANSPOSE(T))  MTDS 290
C                   IOP = 3 - A IS REPLACED BY INVERSE(TRANSPOSE(T)*T)*AMTDS 300
C                   IOP =-3 - A IS REPLACED BY A*INVERSE(TRANSPOSE(T)*T)MTDS 310
C           IER   - RESULTING ERROR PARAMETER                           MTDS 320
C                   IER =-1 MEANS M AND N ARE NOT BOTH POSITIVE         MTDS 330
C                                 AND/OR IOP IS ILLEGAL                 MTDS 340
C                   IER = 0 MEANS OPERATION WAS SUCCESSFUL              MTDS 350
C                   IER = 1 MEANS TRIANGULAR MATRIX T IS SINGULAR       MTDS 360
C                                                                       MTDS 370
C        REMARKS                                                        MTDS 380
C           SUBROUTINE MTDS MAY BE USED TO CALCULATE THE SOLUTION OF    MTDS 390
C           A SYSTEM OF EQUATIONS WITH SYMMETRIC POSITIVE DEFINITE      MTDS 400
C           COEFFICIENT MATRIX. THE FIRST STEP TOWARDS THE SOLUTION     MTDS 410
C           IS TRIANGULAR FACTORIZATION BY MEANS OF MFSD, THE SECOND    MTDS 420
C           STEP IS APPLICATION OF MTDS.                                MTDS 430
C           SUBROUTINES MFSD AND MTDS MAY BE USED IN ORDER TO CALCULATE MTDS 440
C           THE PRODUCT TRANSPOSE(A)*INVERSE(B)*A WITH GIVEN SYMMETRIC  MTDS 450
C           POSITIVE DEFINITE B AND GIVEN A EFFICIENTLY IN THREE STEPS  MTDS 460
C           1) TRIANGULAR FACTORIZATION OF B (B=TRANSPOSE(T)*T)         MTDS 470
C           2) MULTIPLICATION OF A ON THE LEFT BY INVERSE(TRANSPOSE(T)) MTDS 480
C              A IS REPLACED BY C=INVERSE(TRANSPOSE(T))*A               MTDS 490
C           3) CALCULATION OF THE RESULT FORMING TRANSPOSE(C)*C         MTDS 500
C                                                                       MTDS 510
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  MTDS 520
C           NONE                                                        MTDS 530
C                                                                       MTDS 540
C        METHOD                                                         MTDS 550
C           CALCULATION OF X = INVERSE(T)*A IS DONE USING BACKWARD      MTDS 560
C           SUBSTITUTION TO OBTAIN X FROM T*X = A.                      MTDS 570
C           CALCULATION OF Y = INVERSE(TRANSPOSE(T))*A IS DONE USING    MTDS 580
C           FORWARD SUBSTITUTION TO OBTAIN Y FROM TRANSPOSE(T)*Y = A.   MTDS 590
C           CALCULATION OF Z = INVERSE(TRANSPOSE(T)*T)*A IS DONE        MTDS 600
C           SOLVING FIRST TRANSPOSE(T)*Y = A AND THEN T*Z = Y, IE.      MTDS 610
C           USING THE ABOVE TWO STEPS IN REVERSE ORDER                  MTDS 620
C                                                                       MTDS 630
C     ..................................................................MTDS 640
C                                                                       MTDS 650
      SUBROUTINE MTDS(A,M,N,T,IOP,IER)                                  MTDS 660
C                                                                       MTDS 670
C                                                                       MTDS 680
      DIMENSION A(1),T(1)                                               MTDS 690
      DOUBLE PRECISION DSUM                                             MTDS 700
C                                                                       MTDS 710
C        TEST OF DIMENSION                                              MTDS 720
      IF(M)2,2,1                                                        MTDS 730
    1 IF(N)2,2,4                                                        MTDS 740
C                                                                       MTDS 750
C        ERROR RETURN IN CASE OF ILLEGAL DIMENSIONS                     MTDS 760
    2 IER=-1                                                            MTDS 770
      RETURN                                                            MTDS 780
C                                                                       MTDS 790
C        ERROR RETURN IN CASE OF SINGULAR MATRIX T                      MTDS 800
    3 IER=1                                                             MTDS 810
      RETURN                                                            MTDS 820
C                                                                       MTDS 830
C        INITIALIZE DIVISION PROCESS                                    MTDS 840
    4 MN=M*N                                                            MTDS 850
      MM=M*(M+1)/2                                                      MTDS 860
      MM1=M-1                                                           MTDS 870
      IER=0                                                             MTDS 880
      ICS=M                                                             MTDS 890
      IRS=1                                                             MTDS 900
      IMEND=M                                                           MTDS 910
C                                                                       MTDS 920
C        TEST SPECIFIED OPERATION                                       MTDS 930
      IF(IOP)5,2,6                                                      MTDS 940
    5 MM=N*(N+1)/2                                                      MTDS 950
      MM1=N-1                                                           MTDS 960
      IRS=M                                                             MTDS 970
      ICS=1                                                             MTDS 980
      IMEND=MN-M+1                                                      MTDS 990
      MN=M                                                              MTDS1000
    6 IOPE=MOD(IOP+3,3)                                                 MTDS1010
      IF(IABS(IOP)-3)7,7,2                                              MTDS1020
    7 IF(IOPE-1)8,18,8                                                  MTDS1030
C                                                                       MTDS1040
C        INITIALIZE SOLUTION OF TRANSPOSE(T)*X = A                      MTDS1050
    8 MEND=1                                                            MTDS1060
      LLD=IRS                                                           MTDS1070
      MSTA=1                                                            MTDS1080
      MDEL=1                                                            MTDS1090
      MX=1                                                              MTDS1100
      LD=1                                                              MTDS1110
      LX=0                                                              MTDS1120
C                                                                       MTDS1130
C        TEST FOR NONZERO DIAGONAL TERM IN T                            MTDS1140
    9 IF(T(MSTA))10,3,10                                                MTDS1150
   10 DO 11 I=MEND,MN,ICS                                               MTDS1160
   11 A(I)=A(I)/DBLE(T(MSTA))                                           MTDS1170
C                                                                       MTDS1180
C        IS M EQUAL 1                                                   MTDS1190
      IF(MM1)2,15,12                                                    MTDS1200
   12 DO 14 J=1,MM1                                                     MTDS1210
      MSTA=MSTA+MDEL                                                    MTDS1220
      MDEL=MDEL+MX                                                      MTDS1230
      DO 14 I=MEND,MN,ICS                                               MTDS1240
      DSUM=0.D0                                                         MTDS1250
      L=MSTA                                                            MTDS1260
      LDX=LD                                                            MTDS1270
      LL=I                                                              MTDS1280
      DO 13 K=1,J                                                       MTDS1290
      DSUM=DSUM-T(L)*A(LL)                                              MTDS1300
      LL=LL+LLD                                                         MTDS1310
      L=L+LDX                                                           MTDS1320
   13 LDX=LDX+LX                                                        MTDS1330
      IF(T(L))14,3,14                                                   MTDS1340
   14 A(LL)=(DSUM+A(LL))/T(L)                                           MTDS1350
C                                                                       MTDS1360
C        TEST END OF OPERATION                                          MTDS1370
   15 IF(IER)16,17,16                                                   MTDS1380
   16 IER=0                                                             MTDS1390
      RETURN                                                            MTDS1400
   17 IF(IOPE)18,18,16                                                  MTDS1410
C                                                                       MTDS1420
C        INITIALIZE SOLUTION OF T*X = A                                 MTDS1430
   18 IER=1                                                             MTDS1440
      MEND=IMEND                                                        MTDS1450
      MN=M*N                                                            MTDS1460
      LLD=-IRS                                                          MTDS1470
      MSTA=MM                                                           MTDS1480
      MDEL=-1                                                           MTDS1490
      MX=0                                                              MTDS1500
      LD=-MM1                                                           MTDS1510
      LX=1                                                              MTDS1520
      GOTO 9                                                            MTDS1530
      END                                                               MTDS1540
