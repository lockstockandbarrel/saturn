C                                                                       TEUL  10
C     ..................................................................TEUL  20
C                                                                       TEUL  30
C        SUBROUTINE TEUL                                                TEUL  40
C                                                                       TEUL  50
C        PURPOSE                                                        TEUL  60
C           COMPUTE THE SUM OF FCT(K) FOR K FROM ONE UP TO INFINITY.    TEUL  70
C                                                                       TEUL  80
C        USAGE                                                          TEUL  90
C           CALL TEUL(FCT,SUM,MAX,EPS,IER)                              TEUL 100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT.               TEUL 110
C                                                                       TEUL 120
C        DESCRIPTION OF PARAMETERS                                      TEUL 130
C           FCT    - NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.      TEUL 140
C                    IT COMPUTES THE K-TH TERM OF THE SERIES TO ANY     TEUL 150
C                    GIVEN INDEX K.                                     TEUL 160
C           SUM    - RESULTANT VALUE CONTAINING ON RETURN THE SUM OF    TEUL 170
C                    THE GIVEN SERIES.                                  TEUL 180
C           MAX    - INPUT VALUE, WHICH SPECIFIES THE MAXIMAL NUMBER    TEUL 190
C                    OF TERMS OF THE SERIES THAT ARE RESPECTED.         TEUL 200
C           EPS    - INPUT VALUE, WHICH SPECIFIES THE UPPER BOUND OF    TEUL 210
C                    THE RELATIVE ERROR.                                TEUL 220
C                    SUMMATION IS STOPPED AS SOON AS FIVE TIMES IN      TEUL 230
C                    SUCCESSION THE ABSOLUTE VALUE OF THE TERMS OF THE  TEUL 240
C                    TRANSFORMED SERIES ARE FOUND TO BE LESS THAN       TEUL 250
C                    EPS*(ABSOLUTE VALUE OF CURRENT SUM).               TEUL 260
C           IER    - RESULTANT ERROR PARAMETER CODED IN THE FOLLOWING   TEUL 270
C                    FORM                                               TEUL 280
C                     IER=0  - NO ERROR                                 TEUL 290
C                     IER=1  - REQUIRED ACCURACY NOT REACHED WITH       TEUL 300
C                              MAXIMAL NUMBER OF TERMS                  TEUL 310
C                     IER=-1 - THE INTEGER MAX IS LESS THAN ONE.        TEUL 320
C                                                                       TEUL 330
C        REMARKS                                                        TEUL 340
C           NO ACTION BESIDES ERROR MESSAGE IN CASE MAX LESS THAN ONE.  TEUL 350
C                                                                       TEUL 360
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  TEUL 370
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(K) MUST BE FURNISHED   TEUL 380
C           BY THE USER.                                                TEUL 390
C                                                                       TEUL 400
C        METHOD                                                         TEUL 410
C           EVALUATION IS DONE BY MEANS OF A SUITABLY REFINED EULER     TEUL 420
C           TRANSFORMATION. FOR REFERENCE, SEE                          TEUL 430
C           F.B.HILDEBRAND, INTRODUCTION TO NUMERICAL ANALYSIS,         TEUL 440
C           MCGRAW/HILL, NEW YORK/TORONTO/LONDON, 1956, PP.155-160, AND TEUL 450
C           P. NAUR, REPORT ON THE ALGORITHMIC LANGUAGE ALGOL 60,       TEUL 460
C           CACM, VOL.3, ISS.5 (1960), PP.311.                          TEUL 470
C                                                                       TEUL 480
C     ..................................................................TEUL 490
C                                                                       TEUL 500
      SUBROUTINE TEUL (FCT,SUM,MAX,EPS,IER)                             TEUL 510
C                                                                       TEUL 520
      DIMENSION Y(15)                                                   TEUL 530
C                                                                       TEUL 540
C        TEST ON WRONG INPUT PARAMETER MAX                              TEUL 550
C                                                                       TEUL 560
      IF(MAX)1,1,2                                                      TEUL 570
    1 IER=-1                                                            TEUL 580
      GOTO 12                                                           TEUL 590
C                                                                       TEUL 600
C        INITIALIZE EULER TRANSFORMATION                                TEUL 610
C                                                                       TEUL 620
    2 IER=1                                                             TEUL 630
      I=1                                                               TEUL 640
      M=1                                                               TEUL 650
      N=1                                                               TEUL 660
      Y(1)=FCT(N)                                                       TEUL 670
      SUM=Y(1)*.5                                                       TEUL 680
C                                                                       TEUL 690
C        START EULER-LOOP                                               TEUL 700
C                                                                       TEUL 710
    3 J=0                                                               TEUL 720
    4 I=I+1                                                             TEUL 730
      IF(I-MAX)5,5,12                                                   TEUL 740
    5 N=I                                                               TEUL 750
      AMN=FCT(N)                                                        TEUL 760
      DO 6 K=1,M                                                        TEUL 770
      AMP=(AMN+Y(K))*.5                                                 TEUL 780
      Y(K)=AMN                                                          TEUL 790
    6 AMN=AMP                                                           TEUL 800
C                                                                       TEUL 810
C        CHECK EULER TRANSFORMATION                                     TEUL 820
C                                                                       TEUL 830
      IF(ABS(AMN)-ABS(Y(M)))7,9,9                                       TEUL 840
    7 IF(M-15)8,9,9                                                     TEUL 850
    8 M=M+1                                                             TEUL 860
      Y(M)=AMN                                                          TEUL 870
      AMN=.5*AMN                                                        TEUL 880
C                                                                       TEUL 890
C        UPDATE SUM                                                     TEUL 900
C                                                                       TEUL 910
    9 SUM=SUM+AMN                                                       TEUL 920
      IF(ABS(AMN)-EPS*ABS(SUM))10,10,3                                  TEUL 930
C                                                                       TEUL 940
C        TEST END OF PROCEDURE                                          TEUL 950
C                                                                       TEUL 960
   10 J=J+1                                                             TEUL 970
      IF(J-5)4,11,11                                                    TEUL 980
   11 IER=0                                                             TEUL 990
   12 RETURN                                                            TEUL1000
      END                                                               TEUL1010
