C                                                                       ORDE  10
C     ..................................................................ORDE  20
C                                                                       ORDE  30
C        SUBROUTINE ORDER                                               ORDE  40
C                                                                       ORDE  50
C        PURPOSE                                                        ORDE  60
C           CONSTRUCT FROM A LARGER MATRIX OF CORRELATION COEFFICIENTS  ORDE  70
C           A SUBSET MATRIX OF INTERCORRELATIONS AMONG INDEPENDENT      ORDE  80
C           VARIABLES AND A VECTOR OF INTERCORRELATIONS OF INDEPENDENT  ORDE  90
C           VARIABLES WITH DEPENDENT VARIABLE.  THIS SUBROUTINE IS      ORDE 100
C           NORMALLY USED IN THE PERFORMANCE OF MULTIPLE AND POLYNOMIAL ORDE 110
C           REGRESSION ANALYSES.                                        ORDE 120
C                                                                       ORDE 130
C        USAGE                                                          ORDE 140
C           CALL ORDER (M,R,NDEP,K,ISAVE,RX,RY)                         ORDE 150
C                                                                       ORDE 160
C        DESCRIPTION OF PARAMETERS                                      ORDE 170
C           M     - NUMBER OF VARIABLES AND ORDER OF MATRIX R.          ORDE 180
C           R     - INPUT MATRIX CONTAINING CORRELATION COEFFICIENTS.   ORDE 190
C                   THIS SUBROUTINE EXPECTS ONLY UPPER TRIANGULAR       ORDE 200
C                   PORTION OF THE SYMMETRIC MATRIX TO BE STORED (BY    ORDE 210
C                   COLUMN) IN R.  (STORAGE MODE OF 1)                  ORDE 220
C           NDEP  - THE SUBSCRIPT NUMBER OF THE DEPENDENT VARIABLE.     ORDE 230
C           K     - NUMBER OF INDEPENDENT VARIABLES TO BE INCLUDED      ORDE 240
C                   IN THE FORTHCOMING REGRESSION. K MUST BE GREATER    ORDE 250
C                   THAN OR EQUAL TO 1.                                 ORDE 251
C           ISAVE - INPUT VECTOR OF LENGTH K+1 CONTAINING, IN ASCENDING ORDE 260
C                   ORDER, THE SUBSCRIPT NUMBERS OF K INDEPENDENT       ORDE 270
C                   VARIABLES TO BE INCLUDED IN THE FORTHCOMING REGRES- ORDE 280
C                   SION.                                               ORDE 290
C                   UPON RETURNING TO THE CALLING ROUTINE, THIS VECTOR  ORDE 300
C                   CONTAINS, IN ADDITION, THE SUBSCRIPT NUMBER OF      ORDE 310
C                   THE DEPENDENT VARIABLE IN K+1 POSITION.             ORDE 320
C           RX    - OUTPUT MATRIX (K X K) CONTAINING INTERCORRELATIONS  ORDE 330
C                   AMONG INDEPENDENT VARIABLES TO BE USED IN FORTH-    ORDE 340
C                   COMING REGRESSION.                                  ORDE 350
C           RY    - OUTPUT VECTOR OF LENGTH K CONTAINING INTERCORRELA-  ORDE 360
C                   TIONS OF INDEPENDENT VARIABLES WITH DEPENDENT       ORDE 370
C                   VARIABLES.                                          ORDE 380
C                                                                       ORDE 390
C        REMARKS                                                        ORDE 400
C           NONE                                                        ORDE 410
C                                                                       ORDE 420
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  ORDE 430
C           NONE                                                        ORDE 440
C                                                                       ORDE 450
C        METHOD                                                         ORDE 460
C           FROM THE SUBSCRIPT NUMBERS OF THE VARIABLES TO BE INCLUDED  ORDE 470
C           IN THE FORTHCOMING REGRESSION, THE SUBROUTINE CONSTRUCTS THEORDE 480
C           MATRIX RX AND THE VECTOR RY.                                ORDE 490
C                                                                       ORDE 500
C     ..................................................................ORDE 510
C                                                                       ORDE 520
      SUBROUTINE ORDER (M,R,NDEP,K,ISAVE,RX,RY)                         ORDE 530
      DIMENSION R(1),ISAVE(1),RX(1),RY(1)                               ORDE 540
C                                                                       ORDE 550
C        ...............................................................ORDE 560
C                                                                       ORDE 570
C        IF A DOUBLE PRECISION VERSION OF THIS ROUTINE IS DESIRED, THE  ORDE 580
C        C IN COLUMN 1 SHOULD BE REMOVED FROM THE DOUBLE PRECISION      ORDE 590
C        STATEMENT WHICH FOLLOWS.                                       ORDE 600
C                                                                       ORDE 610
C     DOUBLE PRECISION R,RX,RY                                          ORDE 620
C                                                                       ORDE 630
C        THE C MUST ALSO BE REMOVED FROM DOUBLE PRECISION STATEMENTS    ORDE 640
C        APPEARING IN OTHER ROUTINES USED IN CONJUNCTION WITH THIS      ORDE 650
C        ROUTINE.                                                       ORDE 660
C                                                                       ORDE 670
C        ...............................................................ORDE 680
C                                                                       ORDE 690
C     COPY INTERCORRELATIONS OF INDEPENDENT VARIABLES                   ORDE 700
C     WITH DEPENDENT VARIABLE                                           ORDE 710
C                                                                       ORDE 720
      MM=0                                                              ORDE 730
      DO 130 J=1,K                                                      ORDE 740
      L2=ISAVE(J)                                                       ORDE 750
      IF(NDEP-L2) 122, 123, 123                                         ORDE 760
  122 L=NDEP+(L2*L2-L2)/2                                               ORDE 770
      GO TO 125                                                         ORDE 780
  123 L=L2+(NDEP*NDEP-NDEP)/2                                           ORDE 790
  125 RY(J)=R(L)                                                        ORDE 800
C                                                                       ORDE 810
C     COPY A SUBSET MATRIX OF INTERCORRELATIONS AMONG                   ORDE 820
C     INDEPENDENT VARIABLES                                             ORDE 830
C                                                                       ORDE 840
      DO 130 I=1,K                                                      ORDE 850
      L1=ISAVE(I)                                                       ORDE 860
      IF(L1-L2) 127, 128, 128                                           ORDE 870
  127 L=L1+(L2*L2-L2)/2                                                 ORDE 880
      GO TO 129                                                         ORDE 890
  128 L=L2+(L1*L1-L1)/2                                                 ORDE 900
  129 MM=MM+1                                                           ORDE 910
  130 RX(MM)=R(L)                                                       ORDE 920
C                                                                       ORDE 930
C     PLACE THE SUBSCRIPT NUMBER OF THE DEPENDENT                       ORDE 940
C     VARIABLE IN ISAVE(K+1)                                            ORDE 950
C                                                                       ORDE 960
      ISAVE(K+1)=NDEP                                                   ORDE 970
      RETURN                                                            ORDE 980
      END                                                               ORDE 990
