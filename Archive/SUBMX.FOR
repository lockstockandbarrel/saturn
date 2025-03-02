C                                                                       SUBM  10
C     ..................................................................SUBM  20
C                                                                       SUBM  30
C        SUBROUTINE SUBMX                                               SUBM  40
C                                                                       SUBM  50
C        PURPOSE                                                        SUBM  60
C           BASED ON VECTOR S DERIVED FROM SUBROUTINE SUBST OR ABSNT,   SUBM  70
C           THIS SUBROUTINE COPIES FROM A LARGER MATRIX OF OBSERVATION  SUBM  80
C           DATA A SUBSET MATRIX OF THOSE OBSERVATIONS WHICH HAVE       SUBM  90
C           SATISFIED CERTAIN CONDITION.  THIS SUBROUTINE IS NORMALLY   SUBM 100
C           USED PRIOR TO STATISTICAL ANALYSES (E.G., MULTIPLE REGRES-  SUBM 110
C           SION, FACTOR ANALYSIS).                                     SUBM 120
C                                                                       SUBM 130
C        USAGE                                                          SUBM 140
C           CALL SUBMX (A,D,S,NO,NV,N)                                  SUBM 150
C                                                                       SUBM 160
C        DESCRIPTION OF PARAMETERS                                      SUBM 170
C           A  - INPUT MATRIX OF OBSERVATIONS, NO BY NV.                SUBM 180
C           D  - OUTPUT MATRIX OF OBSERVATIONS, N BY NV.                SUBM 190
C           S -  INPUT VECTOR OF LENGTH NO CONTAINING THE CODES DERIVED SUBM 200
C                FROM SUBROUTINE SUBST OR ABSNT.                        SUBM 210
C           NO - NUMBER OF OBSERVATIONS. NO MUST BE > OR = TO 1.        SUBM 220
C           NV - NUMBER OF VARIABLES. NV MUST BE > OR = TO 1.           SUBM 230
C           N  - OUTPUT VARIABLE CONTAINING THE NUMBER OF NON-ZERO CODESSUBM 240
C                IN VECTOR S.                                           SUBM 250
C                                                                       SUBM 260
C        REMARKS                                                        SUBM 270
C           MATRIX D CAN BE IN THE SAME LOCATION AS MATRIX A.           SUBM 280
C                                                                       SUBM 290
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SUBM 300
C           NONE                                                        SUBM 310
C                                                                       SUBM 320
C        METHOD                                                         SUBM 330
C           IF S(I) CONTAINS A NON-ZERO CODE, I-TH OBSERVATION IS       SUBM 340
C           COPIED FROM THE INPUT MATRIX TO THE OUTPUT MATRIX.          SUBM 350
C                                                                       SUBM 360
C     ..................................................................SUBM 370
C                                                                       SUBM 380
      SUBROUTINE SUBMX (A,D,S,NO,NV,N)                                  SUBM 390
      DIMENSION A(1),D(1),S(1)                                          SUBM 400
C                                                                       SUBM 410
      L=0                                                               SUBM 420
      LL=0                                                              SUBM 430
      DO 20 J=1,NV                                                      SUBM 440
      DO 15 I=1,NO                                                      SUBM 450
      L=L+1                                                             SUBM 460
      IF(S(I)) 15, 15, 10                                               SUBM 470
   10 LL=LL+1                                                           SUBM 480
      D(LL)=A(L)                                                        SUBM 490
   15 CONTINUE                                                          SUBM 500
   20 CONTINUE                                                          SUBM 510
C                                                                       SUBM 520
C        COUNT NON-ZERO CODES IN VECTOR S                               SUBM 530
C                                                                       SUBM 540
      N=0                                                               SUBM 550
      DO 30 I=1,NO                                                      SUBM 560
      IF(S(I)) 30, 30, 25                                               SUBM 570
   25 N=N+1                                                             SUBM 580
   30 CONTINUE                                                          SUBM 590
C                                                                       SUBM 600
      RETURN                                                            SUBM 610
      END                                                               SUBM 620
