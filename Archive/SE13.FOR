C                                                                       SE13  10
C     ..................................................................SE13  20
C                                                                       SE13  30
C        SUBROUTINE SE13                                                SE13  40
C                                                                       SE13  50
C        PURPOSE                                                        SE13  60
C           TO COMPUTE A VECTOR OF SMOOTHED FUNCTION VALUES GIVEN A     SE13  70
C           VECTOR OF FUNCTION VALUES WHOSE ENTRIES CORRESPOND TO       SE13  80
C           EQUIDISTANTLY SPACED ARGUMENT VALUES.                       SE13  90
C                                                                       SE13 100
C        USAGE                                                          SE13 110
C           CALL SE13(Y,Z,NDIM,IER)                                     SE13 120
C                                                                       SE13 130
C        DESCRIPTION OF PARAMETERS                                      SE13 140
C           Y     -  GIVEN VECTOR OF FUNCTION VALUES (DIMENSION NDIM)   SE13 150
C           Z     -  RESULTING VECTOR OF SMOOTHED FUNCTION VALUES       SE13 160
C                    (DIMENSION NDIM)                                   SE13 170
C           NDIM  -  DIMENSION OF VECTORS Y AND Z                       SE13 180
C           IER   -  RESULTING ERROR PARAMETER                          SE13 190
C                    IER = -1  - NDIM IS LESS THAN 3                    SE13 200
C                    IER =  0  - NO ERROR                               SE13 210
C                                                                       SE13 220
C        REMARKS                                                        SE13 230
C           (1)  IF IER=-1 THERE HAS BEEN NO COMPUTATION.               SE13 240
C           (2)   Z CAN HAVE THE SAME STORAGE ALLOCATION AS Y.  IF Y    SE13 250
C                 IS DISTINCT FROM Z, THEN IT IS NOT DESTROYED.         SE13 260
C                                                                       SE13 270
C        SUBROUTINES AND SUBPROGRAMS REQUIRED                           SE13 280
C           NONE                                                        SE13 290
C                                                                       SE13 300
C        METHOD                                                         SE13 310
C           IF X IS THE (SUPPRESSED) VECTOR OF ARGUMENT VALUES, THEN    SE13 320
C           EXCEPT AT THE ENDPOINTS X(1) AND X(NDIM), EACH SMOOTHED     SE13 330
C           VALUE Z(I) IS OBTAINED BY EVALUATING AT X(I) THE LEAST-     SE13 340
C           SQUARES POLYNOMIAL OF DEGREE 1 RELEVANT TO THE 3 SUCCESSIVE SE13 350
C           POINTS (X(I+K),Y(I+K)) K = -1,0,1.  (SEE HILDEBRAND, F.B.,  SE13 360
C           INTRODUCTION TO NUMERICAL ANALYSIS, MC GRAW-HILL, NEW YORK/ SE13 370
C           TORONTO/LONDON, 1956, PP. 295-302.)                         SE13 380
C                                                                       SE13 390
C     ..................................................................SE13 400
C                                                                       SE13 410
      SUBROUTINE SE13(Y,Z,NDIM,IER)                                     SE13 420
C                                                                       SE13 430
C                                                                       SE13 440
      DIMENSION Y(1),Z(1)                                               SE13 450
C                                                                       SE13 460
C        TEST OF DIMENSION                                              SE13 470
      IF(NDIM-3)3,1,1                                                   SE13 480
C                                                                       SE13 490
C        PREPARE LOOP                                                   SE13 500
    1 B=.1666667*(5.*Y(1)+Y(2)+Y(2)-Y(3))                               SE13 510
      C=.1666667*(5.*Y(NDIM)+Y(NDIM-1)+Y(NDIM-1)-Y(NDIM-2))             SE13 520
C                                                                       SE13 530
C        START LOOP                                                     SE13 540
      DO 2 I=3,NDIM                                                     SE13 550
      A=B                                                               SE13 560
      B=.3333333*(Y(I-2)+Y(I-1)+Y(I))                                   SE13 570
    2 Z(I-2)=A                                                          SE13 580
C        END OF LOOP                                                    SE13 590
C                                                                       SE13 600
C        UPDATE LAST TWO COMPONENTS                                     SE13 610
      Z(NDIM-1)=B                                                       SE13 620
      Z(NDIM)=C                                                         SE13 630
      IER=0                                                             SE13 640
      RETURN                                                            SE13 650
C                                                                       SE13 660
C        ERROR EXIT IN CASE NDIM IS LESS THAN 3                         SE13 670
    3 IER=-1                                                            SE13 680
      RETURN                                                            SE13 690
      END                                                               SE13 700
