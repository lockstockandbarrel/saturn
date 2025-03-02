C                                                                       DE13  10
C     ..................................................................DE13  20
C                                                                       DE13  30
C        SUBROUTINE DSE13                                               DE13  40
C                                                                       DE13  50
C        PURPOSE                                                        DE13  60
C           TO COMPUTE A VECTOR OF SMOOTHED FUNCTION VALUES GIVEN A     DE13  70
C           VECTOR OF FUNCTION VALUES WHOSE ENTRIES CORRESPOND TO       DE13  80
C           EQUIDISTANTLY SPACED ARGUMENT VALUES.                       DE13  90
C                                                                       DE13 100
C        USAGE                                                          DE13 110
C           CALL DSE13(Y,Z,NDIM,IER)                                    DE13 120
C                                                                       DE13 130
C        DESCRIPTION OF PARAMETERS                                      DE13 140
C           Y     -  GIVEN VECTOR OF DOUBLE PRECISION FUNCTION VALUES   DE13 150
C                    (DIMENSION NDIM)                                   DE13 160
C           Z     -  RESULTING VECTOR OF DOUBLE PRECISION SMOOTHED      DE13 170
C                    FUNCTION VALUES (DIMENSION NDIM)                   DE13 180
C           NDIM  -  DIMENSION OF VECTORS Y AND Z                       DE13 190
C           IER   -  RESULTING ERROR PARAMETER                          DE13 200
C                    IER = -1  - NDIM IS LESS THAN 3                    DE13 210
C                    IER =  0  - NO ERROR                               DE13 220
C                                                                       DE13 230
C        REMARKS                                                        DE13 240
C           (1)  IF IER=-1 THERE HAS BEEN NO COMPUTATION.               DE13 250
C           (2)   Z CAN HAVE THE SAME STORAGE ALLOCATION AS Y.  IF Y    DE13 260
C                 IS DISTINCT FROM Z, THEN IT IS NOT DESTROYED.         DE13 270
C                                                                       DE13 280
C        SUBROUTINES AND SUBPROGRAMS REQUIRED                           DE13 290
C           NONE                                                        DE13 300
C                                                                       DE13 310
C        METHOD                                                         DE13 320
C           IF X IS THE (SUPPRESSED) VECTOR OF ARGUMENT VALUES, THEN    DE13 330
C           EXCEPT AT THE ENDPOINTS X(1) AND X(NDIM), EACH SMOOTHED     DE13 340
C           VALUE Z(I) IS OBTAINED BY EVALUATING AT X(I) THE LEAST-     DE13 350
C           SQUARES POLYNOMIAL OF DEGREE 1 RELEVANT TO THE 3 SUCCESSIVE DE13 360
C           POINTS (X(I+K),Y(I+K)) K = -1,0,1.  (SEE HILDEBRAND, F.B.,  DE13 370
C           INTRODUCTION TO NUMERICAL ANALYSIS, MC GRAW-HILL, NEW YORK/ DE13 380
C           TORONTO/LONDON, 1956, PP. 295-302.)                         DE13 390
C                                                                       DE13 400
C     ..................................................................DE13 410
C                                                                       DE13 420
      SUBROUTINE DSE13(Y,Z,NDIM,IER)                                    DE13 430
C                                                                       DE13 440
      DIMENSION Y(1),Z(1)                                               DE13 450
      DOUBLE PRECISION Y,Z,A,B,C                                        DE13 460
C                                                                       DE13 470
C        TEST OF DIMENSION                                              DE13 480
      IF(NDIM-3)3,1,1                                                   DE13 490
C                                                                       DE13 500
C        PREPARE LOOP                                                   DE13 510
    1 B=.16666666666666667D0*(5.D0*Y(1)+Y(2)+Y(2)-Y(3))                 DE13 520
      C=.16666666666666667*(5.D0*Y(NDIM)+Y(NDIM-1)+Y(NDIM-1)-Y(NDIM-2)) DE13 530
C                                                                       DE13 540
C        START LOOP                                                     DE13 550
      DO 2 I=3,NDIM                                                     DE13 560
      A=B                                                               DE13 570
      B=.33333333333333333D0*(Y(I-2)+Y(I-1)+Y(I))                       DE13 580
    2 Z(I-2)=A                                                          DE13 590
C        END OF LOOP                                                    DE13 600
C                                                                       DE13 610
C        UPDATE LAST TWO COMPONENTS                                     DE13 620
      Z(NDIM-1)=B                                                       DE13 630
      Z(NDIM)=C                                                         DE13 640
      IER=0                                                             DE13 650
      RETURN                                                            DE13 660
C                                                                       DE13 670
C        ERROR EXIT IN CASE NDIM IS LESS THAN 3                         DE13 680
    3 IER=-1                                                            DE13 690
      RETURN                                                            DE13 700
      END                                                               DE13 710
