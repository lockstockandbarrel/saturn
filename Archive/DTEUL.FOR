C                                                                       DTEU  10
C     ..................................................................DTEU  20
C                                                                       DTEU  30
C        SUBROUTINE DTEUL                                               DTEU  40
C                                                                       DTEU  50
C        PURPOSE                                                        DTEU  60
C           COMPUTE THE SUM OF FCT(K) FOR K FROM ONE UP TO INFINITY.    DTEU  70
C                                                                       DTEU  80
C        USAGE                                                          DTEU  90
C           CALL DTEUL(FCT,SUM,MAX,EPS,IER)                             DTEU 100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT.               DTEU 110
C                                                                       DTEU 120
C        DESCRIPTION OF PARAMETERS                                      DTEU 130
C           FCT    - NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION      DTEU 140
C                    SUBPROGRAM USED. IT COMPUTES THE K-TH TERM OF THE  DTEU 150
C                    SERIES TO ANY GIVEN INDEX K.                       DTEU 160
C           SUM    - RESULTANT VALUE IN DOUBLE PRECISION CONTAINING ON  DTEU 170
C                    RETURN THE SUM OF THE GIVEN SERIES.                DTEU 180
C           MAX    - INPUT VALUE, WHICH SPECIFIES THE MAXIMAL NUMBER    DTEU 190
C                    OF TERMS OF THE SERIES THAT ARE RESPECTED.         DTEU 200
C           EPS    - SINGLE PRECISION INPUT VALUE, WHICH SPECIFIES THE  DTEU 210
C                    UPPER BOUND OF THE RELATIVE ERROR.                 DTEU 220
C                    SUMMATION IS STOPPED AS SOON AS FIVE TIMES IN      DTEU 230
C                    SUCCESSION THE ABSOLUTE VALUE OF THE TERMS OF THE  DTEU 240
C                    TRANSFORMED SERIES ARE FOUND TO BE LESS THAN       DTEU 250
C                    EPS*(ABSOLUTE VALUE OF CURRENT SUM).               DTEU 260
C           IER    - RESULTANT ERROR PARAMETER CODED IN THE FOLLOWING   DTEU 270
C                    FORM                                               DTEU 280
C                     IER=0  - NO ERROR                                 DTEU 290
C                     IER=1  - REQUIRED ACCURACY NOT REACHED WITH       DTEU 300
C                              MAXIMAL NUMBER OF TERMS                  DTEU 310
C                     IER=-1 - THE INTEGER MAX IS LESS THAN ONE.        DTEU 320
C                                                                       DTEU 330
C        REMARKS                                                        DTEU 340
C           NO ACTION BESIDES ERROR MESSAGE IN CASE MAX LESS THAN ONE.  DTEU 350
C                                                                       DTEU 360
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DTEU 370
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(K) MUST BE FURNISHED   DTEU 380
C           BY THE USER.                                                DTEU 390
C                                                                       DTEU 400
C        METHOD                                                         DTEU 410
C           EVALUATION IS DONE BY MEANS OF A SUITABLY REFINED EULER     DTEU 420
C           TRANSFORMATION. FOR REFERENCE, SEE                          DTEU 430
C           F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,         DTEU 440
C           MCGRAW/HILL, NEW YORK/TORONTO/LONDON, 1956, PP.155-160, AND DTEU 450
C           P. NAUR, REPORT ON THE ALGORITHMIC LANGUAGE ALGOL 60,       DTEU 460
C           CACM, VOL.3, ISS.5 (1960), PP.311.                          DTEU 470
C                                                                       DTEU 480
C     ..................................................................DTEU 490
C                                                                       DTEU 500
      SUBROUTINE DTEUL (FCT,SUM,MAX,EPS,IER)                            DTEU 510
C                                                                       DTEU 520
      DIMENSION Y(15)                                                   DTEU 530
      DOUBLE PRECISION FCT,SUM,Y,AMN,AMP                                DTEU 540
C                                                                       DTEU 550
C        TEST ON WRONG INPUT PARAMETER MAX                              DTEU 560
C                                                                       DTEU 570
      IF(MAX)1,1,2                                                      DTEU 580
    1 IER=-1                                                            DTEU 590
      GOTO 12                                                           DTEU 600
C                                                                       DTEU 610
C        INITIALIZE EULER TRANSFORMATION                                DTEU 620
C                                                                       DTEU 630
    2 IER=1                                                             DTEU 640
      I=1                                                               DTEU 650
      M=1                                                               DTEU 660
      N=1                                                               DTEU 670
      Y(1)=FCT(N)                                                       DTEU 680
      SUM=Y(1)*.5D0                                                     DTEU 690
C                                                                       DTEU 700
C        START EULER-LOOP                                               DTEU 710
C                                                                       DTEU 720
    3 J=0                                                               DTEU 730
    4 I=I+1                                                             DTEU 740
      IF(I-MAX)5,5,12                                                   DTEU 750
    5 N=I                                                               DTEU 760
      AMN=FCT(N)                                                        DTEU 770
      DO 6 K=1,M                                                        DTEU 780
      AMP=(AMN+Y(K))*.5D0                                               DTEU 790
      Y(K)=AMN                                                          DTEU 800
    6 AMN=AMP                                                           DTEU 810
C                                                                       DTEU 820
C        CHECK EULER TRANSFORMATION                                     DTEU 830
C                                                                       DTEU 840
      IF(DABS(AMN)-DABS(Y(M)))7,9,9                                     DTEU 850
    7 IF(M-15)8,9,9                                                     DTEU 860
    8 M=M+1                                                             DTEU 870
      Y(M)=AMN                                                          DTEU 880
      AMN=.5D0*AMN                                                      DTEU 890
C                                                                       DTEU 900
C        UPDATE SUM                                                     DTEU 910
C                                                                       DTEU 920
    9 SUM=SUM+AMN                                                       DTEU 930
      IF(ABS(SNGL(AMN))-EPS*ABS(SNGL(SUM)))10,10,3                      DTEU 940
C                                                                       DTEU 950
C        TEST END OF PROCEDURE                                          DTEU 960
C                                                                       DTEU 970
   10 J=J+1                                                             DTEU 980
      IF(J-5)4,11,11                                                    DTEU 990
   11 IER=0                                                             DTEU1000
   12 RETURN                                                            DTEU1010
      END                                                               DTEU1020
