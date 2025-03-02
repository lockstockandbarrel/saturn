C                                                                       FRIF  10
C     ..................................................................FRIF  20
C                                                                       FRIF  30
C        SUBROUTINE FORIF                                               FRIF  40
C                                                                       FRIF  50
C        PURPOSE                                                        FRIF  60
C           FOURIER ANALYSIS OF A GIVEN PERIODIC FUNCTION IN THE        FRIF  70
C           RANGE 0-2PI                                                 FRIF  80
C           COMPUTES THE COEFFICIENTS OF THE DESIRED NUMBER OF TERMS    FRIF  90
C           IN THE FOURIER SERIES F(X)=A(0)+SUM(A(K)COS KX+B(K)SIN KX)  FRIF 100
C           WHERE K=1,2,...,M TO APPROXIMATE THE COMPUTED VALUES OF A   FRIF 110
C           GIVEN FUNCTION SUBPROGRAM                                   FRIF 120
C                                                                       FRIF 130
C        USAGE                                                          FRIF 140
C           CALL FORIF(FUN,N,M,A,B,IER)                                 FRIF 150
C                                                                       FRIF 160
C        DESCRIPTION OF PARAMETERS                                      FRIF 170
C           FUN-NAME OF FUNCTION SUBPROGRAM TO BE USED FOR COMPUTING    FRIF 180
C               DATA POINTS                                             FRIF 190
C           N  -DEFINES THE INTERVAL SUCH THAT 2N+1 POINTS ARE TAKEN    FRIF 200
C               OVER THE INTERVAL (0,2PI). THE SPACING IS THUS 2PI/2N+1 FRIF 210
C           M  -THE MAXIMUM ORDER OF THE HARMONICS TO BE FITTED         FRIF 220
C           A  -RESULTANT VECTOR OF FOURIER COSINE COEFFICIENTS OF      FRIF 230
C               LENGTH M+1                                              FRIF 240
C               A SUB 0, A SUB 1,..., A SUB M                           FRIF 250
C           B  -RESULTANT VECTOR OF FOURIER SINE COEFFICIENTS OF        FRIF 260
C               LENGTH M+1                                              FRIF 270
C               B SUB 0, B SUB 1,..., B SUB M                           FRIF 280
C           IER-RESULTANT ERROR CODE WHERE                              FRIF 290
C               IER=0  NO ERROR                                         FRIF 300
C               IER=1  N NOT GREATER OR EQUAL TO M                      FRIF 310
C               IER=2  M LESS THAN 0                                    FRIF 320
C                                                                       FRIF 330
C        REMARKS                                                        FRIF 340
C           M MUST BE GREATER THAN OR EQUAL TO ZERO                     FRIF 350
C           N MUST BE GREATER THAN OR EQUAL TO M                        FRIF 360
C           THE FIRST ELEMENT IN VECTOR B IS ZERO IN ALL CASES          FRIF 370
C                                                                       FRIF 380
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  FRIF 390
C           FUN-NAME OF USER FUNCTION SUBPROGRAM USED FOR COMPUTING     FRIF 400
C               DATA POINTS                                             FRIF 410
C           CALLING PROGRAM MUST HAVE FORTRAN EXTERNAL STATEMENT        FRIF 420
C           CONTAINING NAMES OF FUNCTION SUBPROGRAMS LISTED IN CALL TO  FRIF 430
C           FORIF                                                       FRIF 440
C                                                                       FRIF 450
C        METHOD                                                         FRIF 460
C           USES RECURSIVE TECHNIQUE DESCRIBED IN A. RALSTON, H. WILF,  FRIF 470
C           'MATHEMATICAL METHODS FOR DIGITAL COMPUTERS', JOHN WILEY    FRIF 480
C           AND SONS, NEW YORK, 1960, CHAPTER 24. THE METHOD OF         FRIF 490
C           INDEXING THROUGH THE PROCEDURE HAS BEEN MODIFIED TO         FRIF 500
C           SIMPLIFY THE COMPUTATION.                                   FRIF 510
C                                                                       FRIF 520
C     ..................................................................FRIF 530
C                                                                       FRIF 540
      SUBROUTINE FORIF(FUN,N,M,A,B,IER)                                 FRIF 550
      DIMENSION A(1),B(1)                                               FRIF 560
C                                                                       FRIF 570
C        CHECK FOR PARAMETER ERRORS                                     FRIF 580
C                                                                       FRIF 590
      IER=0                                                             FRIF 600
   20 IF(M) 30,40,40                                                    FRIF 610
   30 IER=2                                                             FRIF 620
      RETURN                                                            FRIF 630
   40 IF(M-N) 60,60,50                                                  FRIF 640
   50 IER=1                                                             FRIF 650
      RETURN                                                            FRIF 660
C                                                                       FRIF 670
C        COMPUTE AND PRESET CONSTANTS                                   FRIF 680
C                                                                       FRIF 690
   60 AN=N                                                              FRIF 700
      COEF=2.0/(2.0*AN+1.0)                                             FRIF 710
      CONST=3.141593*COEF                                               FRIF 720
      S1=SIN(CONST)                                                     FRIF 730
      C1=COS(CONST)                                                     FRIF 740
      C=1.0                                                             FRIF 750
      S=0.0                                                             FRIF 760
      J=1                                                               FRIF 770
      FUNZ=FUN(0.0)                                                     FRIF 780
   70 U2=0.0                                                            FRIF 790
      U1=0.0                                                            FRIF 800
      AI=2*N                                                            FRIF 810
C                                                                       FRIF 820
C        FORM FOURIER COEFFICIENTS RECURSIVELY                          FRIF 830
C                                                                       FRIF 840
   75 X=AI*CONST                                                        FRIF 850
      U0=FUN(X)+2.0*C*U1-U2                                             FRIF 860
      U2=U1                                                             FRIF 870
      U1=U0                                                             FRIF 880
      AI=AI-1.0                                                         FRIF 890
      IF(AI) 80,80,75                                                   FRIF 900
   80 A(J)=COEF*(FUNZ+C*U1-U2)                                          FRIF 910
      B(J)=COEF*S*U1                                                    FRIF 920
      IF(J-(M+1)) 90,100,100                                            FRIF 930
   90 Q=C1*C-S1*S                                                       FRIF 940
      S=C1*S+S1*C                                                       FRIF 950
      C=Q                                                               FRIF 960
      J=J+1                                                             FRIF 970
      GO TO 70                                                          FRIF 980
  100 A(1)=A(1)*0.5                                                     FRIF 990
      RETURN                                                            FRIF1000
      END                                                               FRIF1010
