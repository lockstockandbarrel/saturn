C                                                                       FRIT  10
C     ..................................................................FRIT  20
C                                                                       FRIT  30
C        SUBROUTINE FORIT                                               FRIT  40
C                                                                       FRIT  50
C        PURPOSE                                                        FRIT  60
C           FOURIER ANALYSIS OF A PERIODICALLY TABULATED FUNCTION.      FRIT  70
C           COMPUTES THE COEFFICIENTS OF THE DESIRED NUMBER OF TERMS    FRIT  80
C           IN THE FOURIER SERIES F(X)=A(0)+SUM(A(K)COS KX+B(K)SIN KX)  FRIT  90
C           WHERE K=1,2,...,M TO APPROXIMATE A GIVEN SET OF             FRIT 100
C           PERIODICALLY TABULATED VALUES OF A FUNCTION.                FRIT 110
C                                                                       FRIT 120
C        USAGE                                                          FRIT 130
C           CALL FORIT(FNT,N,M,A,B,IER)                                 FRIT 140
C                                                                       FRIT 150
C        DESCRIPTION OF PARAMETERS                                      FRIT 160
C           FNT-VECTOR OF TABULATED FUNCTION VALUES OF LENGTH 2N+1      FRIT 170
C           N  -DEFINES THE INTERVAL SUCH THAT 2N+1 POINTS ARE TAKEN    FRIT 180
C               OVER THE INTERVAL (0,2PI). THE SPACING IS THUS 2PI/2N+1 FRIT 190
C           M  -MAXIMUM ORDER OF HARMONICS TO BE FITTED                 FRIT 200
C           A  -RESULTANT VECTOR OF FOURIER COSINE COEFFICIENTS OF      FRIT 210
C               LENGTH M+1                                              FRIT 220
C               A SUB 0, A SUB 1,..., A SUB M                           FRIT 230
C           B  -RESULTANT VECTOR OF FOURIER SINE COEFFICIENTS OF        FRIT 240
C               LENGTH M+1                                              FRIT 250
C               B SUB 0, B SUB 1,..., B SUB M                           FRIT 260
C           IER-RESULTANT ERROR CODE WHERE                              FRIT 270
C               IER=0  NO ERROR                                         FRIT 280
C               IER=1  N NOT GREATER OR EQUAL TO M                      FRIT 290
C               IER=2  M LESS THAN 0                                    FRIT 300
C                                                                       FRIT 310
C        REMARKS                                                        FRIT 320
C           M MUST BE GREATER THAN OR EQUAL TO ZERO                     FRIT 330
C           N MUST BE GREATER THAN OR EQUAL TO M                        FRIT 340
C           THE FIRST ELEMENT OF VECTOR B IS ZERO IN ALL CASES          FRIT 350
C                                                                       FRIT 360
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  FRIT 370
C           NONE                                                        FRIT 380
C                                                                       FRIT 390
C        METHOD                                                         FRIT 400
C           USES RECURSIVE TECHNIQUE DESCRIBED IN A. RALSTON, H. WILF,  FRIT 410
C           'MATHEMATICAL METHODS FOR DIGITAL COMPUTERS', JOHN WILEY    FRIT 420
C           AND SONS, NEW YORK, 1960, CHAPTER 24. THE METHOD OF INDEXINGFRIT 430
C           THROUGH THE PROCEDURE HAS BEEN MODIFIED TO SIMPLIFY THE     FRIT 440
C           COMPUTATION.                                                FRIT 450
C                                                                       FRIT 460
C     ..................................................................FRIT 470
C                                                                       FRIT 480
      SUBROUTINE FORIT(FNT,N,M,A,B,IER)                                 FRIT 490
      DIMENSION A(1),B(1),FNT(1)                                        FRIT 500
C                                                                       FRIT 510
C        CHECK FOR PARAMETER ERRORS                                     FRIT 520
C                                                                       FRIT 530
      IER=0                                                             FRIT 540
   20 IF(M) 30,40,40                                                    FRIT 550
   30 IER=2                                                             FRIT 560
      RETURN                                                            FRIT 570
   40 IF(M-N) 60,60,50                                                  FRIT 580
   50 IER=1                                                             FRIT 590
      RETURN                                                            FRIT 600
C                                                                       FRIT 610
C        COMPUTE AND PRESET CONSTANTS                                   FRIT 620
C                                                                       FRIT 630
   60 AN=N                                                              FRIT 640
      COEF=2.0/(2.0*AN+1.0)                                             FRIT 650
      CONST=3.141593*COEF                                               FRIT 660
      S1=SIN(CONST)                                                     FRIT 670
      C1=COS(CONST)                                                     FRIT 680
      C=1.0                                                             FRIT 690
      S=0.0                                                             FRIT 700
      J=1                                                               FRIT 710
      FNTZ=FNT(1)                                                       FRIT 720
   70 U2=0.0                                                            FRIT 730
      U1=0.0                                                            FRIT 740
      I=2*N+1                                                           FRIT 750
C                                                                       FRIT 760
C        FORM FOURIER COEFFICIENTS RECURSIVELY                          FRIT 770
C                                                                       FRIT 780
   75 U0=FNT(I)+2.0*C*U1-U2                                             FRIT 790
      U2=U1                                                             FRIT 800
      U1=U0                                                             FRIT 810
      I=I-1                                                             FRIT 820
      IF(I-1) 80,80,75                                                  FRIT 830
   80 A(J)=COEF*(FNTZ+C*U1-U2)                                          FRIT 840
      B(J)=COEF*S*U1                                                    FRIT 850
      IF(J-(M+1)) 90,100,100                                            FRIT 860
   90 Q=C1*C-S1*S                                                       FRIT 870
      S=C1*S+S1*C                                                       FRIT 880
      C=Q                                                               FRIT 890
      J=J+1                                                             FRIT 900
      GO TO 70                                                          FRIT 910
  100 A(1)=A(1)*0.5                                                     FRIT 920
      RETURN                                                            FRIT 930
      END                                                               FRIT 940
