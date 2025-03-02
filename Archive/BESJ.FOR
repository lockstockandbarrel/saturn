C                                                                       BESJ  10
C     ..................................................................BESJ  20
C                                                                       BESJ  30
C        SUBROUTINE BESJ                                                BESJ  40
C                                                                       BESJ  50
C        PURPOSE                                                        BESJ  60
C           COMPUTE THE J BESSEL FUNCTION FOR A GIVEN ARGUMENT AND ORDERBESJ  70
C                                                                       BESJ  80
C        USAGE                                                          BESJ  90
C           CALL BESJ(X,N,BJ,D,IER)                                     BESJ 100
C                                                                       BESJ 110
C        DESCRIPTION OF PARAMETERS                                      BESJ 120
C           X  -THE ARGUMENT OF THE J BESSEL FUNCTION DESIRED           BESJ 130
C           N  -THE ORDER OF THE J BESSEL FUNCTION DESIRED              BESJ 140
C           BJ -THE RESULTANT J BESSEL FUNCTION                         BESJ 150
C           D  -REQUIRED ACCURACY                                       BESJ 160
C           IER-RESULTANT ERROR CODE WHERE                              BESJ 170
C              IER=0  NO ERROR                                          BESJ 180
C              IER=1  N IS NEGATIVE                                     BESJ 190
C              IER=2  X IS NEGATIVE OR ZERO                             BESJ 200
C              IER=3  REQUIRED ACCURACY NOT OBTAINED                    BESJ 210
C              IER=4  RANGE OF N COMPARED TO X NOT CORRECT (SEE REMARKS)BESJ 220
C                                                                       BESJ 230
C        REMARKS                                                        BESJ 240
C           N MUST BE GREATER THAN OR EQUAL TO ZERO, BUT IT MUST BE     BESJ 250
C           LESS THAN                                                   BESJ 260
C              20+10*X-X** 2/3   FOR X LESS THAN OR EQUAL TO 15         BESJ 270
C              90+X/2           FOR X GREATER THAN 15                   BESJ 280
C                                                                       BESJ 290
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  BESJ 300
C           NONE                                                        BESJ 310
C                                                                       BESJ 320
C        METHOD                                                         BESJ 330
C           RECURRENCE RELATION TECHNIQUE DESCRIBED BY H. GOLDSTEIN AND BESJ 340
C           R.M. THALER,'RECURRENCE TECHNIQUES FOR THE CALCULATION OF   BESJ 350
C           BESSEL FUNCTIONS',M.T.A.C.,V.13,PP.102-108 AND I.A. STEGUN  BESJ 360
C           AND M. ABRAMOWITZ,'GENERATION OF BESSEL FUNCTIONS ON HIGH   BESJ 370
C           SPEED COMPUTERS',M.T.A.C.,V.11,1957,PP.255-257              BESJ 380
C                                                                       BESJ 390
C     ..................................................................BESJ 400
C                                                                       BESJ 410
      SUBROUTINE BESJ(X,N,BJ,D,IER)                                     BESJ 420
C                                                                       BESJ 430
      BJ=.0                                                             BESJ 440
      IF(N)10,20,20                                                     BESJ 450
   10 IER=1                                                             BESJ 460
      RETURN                                                            BESJ 470
   20 IF(X)30,30,31                                                     BESJ 480
   30 IER=2                                                             BESJ 490
      RETURN                                                            BESJ 500
   31 IF(X-15.)32,32,34                                                 BESJ 510
   32 NTEST=20.+10.*X-X** 2/3                                           BESJ 520
      GO TO 36                                                          BESJ 530
   34 NTEST=90.+X/2.                                                    BESJ 540
   36 IF(N-NTEST)40,38,38                                               BESJ 550
   38 IER=4                                                             BESJ 560
      RETURN                                                            BESJ 570
   40 IER=0                                                             BESJ 580
      N1=N+1                                                            BESJ 590
      BPREV=.0                                                          BESJ 600
C                                                                       BESJ 610
C     COMPUTE STARTING VALUE OF M                                       BESJ 620
C                                                                       BESJ 630
      IF(X-5.)50,60,60                                                  BESJ 640
   50 MA=X+6.                                                           BESJ 650
      GO TO 70                                                          BESJ 660
   60 MA=1.4*X+60./X                                                    BESJ 670
   70 MB=N+IFIX(X)/4+2                                                  BESJ 680
      MZERO=MAX0(MA,MB)                                                 BESJ 690
C                                                                       BESJ 700
C     SET UPPER LIMIT OF M                                              BESJ 710
C                                                                       BESJ 720
      MMAX=NTEST                                                        BESJ 730
  100 DO 190 M=MZERO,MMAX,3                                             BESJ 740
C                                                                       BESJ 750
C     SET F(M),F(M-1)                                                   BESJ 760
C                                                                       BESJ 770
      FM1=1.0E-28                                                       BESJ 780
      FM=.0                                                             BESJ 790
      ALPHA=.0                                                          BESJ 800
      IF(M-(M/2)*2)120,110,120                                          BESJ 810
  110 JT=-1                                                             BESJ 820
      GO TO 130                                                         BESJ 830
  120 JT=1                                                              BESJ 840
  130 M2=M-2                                                            BESJ 850
      DO 160 K=1,M2                                                     BESJ 860
      MK=M-K                                                            BESJ 870
      BMK=2.*FLOAT(MK)*FM1/X-FM                                         BESJ 880
      FM=FM1                                                            BESJ 890
      FM1=BMK                                                           BESJ 900
      IF(MK-N-1)150,140,150                                             BESJ 910
  140 BJ=BMK                                                            BESJ 920
  150 JT=-JT                                                            BESJ 930
      S=1+JT                                                            BESJ 940
  160 ALPHA=ALPHA+BMK*S                                                 BESJ 950
      BMK=2.*FM1/X-FM                                                   BESJ 960
      IF(N)180,170,180                                                  BESJ 970
  170 BJ=BMK                                                            BESJ 980
  180 ALPHA=ALPHA+BMK                                                   BESJ 990
      BJ=BJ/ALPHA                                                       BESJ1000
      IF(ABS(BJ-BPREV)-ABS(D*BJ))200,200,190                            BESJ1010
  190 BPREV=BJ                                                          BESJ1020
      IER=3                                                             BESJ1030
  200 RETURN                                                            BESJ1040
      END                                                               BESJ1050
