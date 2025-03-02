C                                                                       DE15  10
C     ..................................................................DE15  20
C                                                                       DE15  30
C        SUBROUTINE DSE15                                               DE15  40
C                                                                       DE15  50
C        PURPOSE                                                        DE15  60
C           TO COMPUTE A VECTOR OF SMOOTHED FUNCTION VALUES GIVEN A     DE15  70
C           VECTOR OF FUNCTION VALUES WHOSE ENTRIES CORRESPOND TO       DE15  80
C           EQUIDISTANTLY SPACED ARGUMENT VALUES.                       DE15  90
C                                                                       DE15 100
C        USAGE                                                          DE15 110
C           CALL DSE15(Y,Z,NDIM,IER)                                    DE15 120
C                                                                       DE15 130
C        DESCRIPTION OF PARAMETERS                                      DE15 140
C           Y     -  GIVEN VECTOR OF DOUBLE PRECISION FUNCTION VALUES   DE15 150
C                    (DIMENSION NDIM)                                   DE15 160
C           Z     -  RESULTING VECTOR OF DOUBLE PRECISION SMOOTHED      DE15 170
C                    FUNCTION VALUES (DIMENSION NDIM)                   DE15 180
C           NDIM  -  DIMENSION OF VECTORS Y AND Z                       DE15 190
C           IER   -  RESULTING ERROR PARAMETER                          DE15 200
C                    IER = -1  - NDIM IS LESS THAN 5                    DE15 210
C                    IER =  0  - NO ERROR                               DE15 220
C                                                                       DE15 230
C        REMARKS                                                        DE15 240
C           (1)  IF IER=-1 THERE HAS BEEN NO COMPUTATION.               DE15 250
C           (2)   Z CAN HAVE THE SAME STORAGE ALLOCATION AS Y.  IF Y IS DE15 260
C                 DISTINCT FROM Z, THEN IT IS NOT DESTROYED.            DE15 270
C                                                                       DE15 280
C        SUBROUTINE AND FUNCTION SUBPROGRAMS REQUIRED                   DE15 290
C           NONE                                                        DE15 300
C                                                                       DE15 310
C        METHOD                                                         DE15 320
C           IF X IS THE (SUPPRESSED) VECTOR OF ARGUMENT VALUES, THEN    DE15 330
C           EXCEPT AT THE POINTS X(1),X(2),X(NDIM-1) AND X(NDIM), EACH  DE15 340
C           SMOOTHED VALUE Z(I) IS OBTAINED BY EVALUATING AT X(I) THE   DE15 350
C           LEAST-SQUARES POLYNOMIAL OF DEGREE 1 RELEVANT TO THE 5      DE15 360
C           SUCCESSIVE POINTS (X(I+K),Y(I+K)) K = -2,-1,...,2.  (SEE    DE15 370
C           HILDEBRAND, F.B., INTRODUCTION TO NUMERICAL ANALYSIS,       DE15 380
C           MC GRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP. 295-302.)  DE15 390
C                                                                       DE15 400
C     ..................................................................DE15 410
C                                                                       DE15 420
      SUBROUTINE DSE15(Y,Z,NDIM,IER)                                    DE15 430
C                                                                       DE15 440
C                                                                       DE15 450
      DIMENSION Y(1),Z(1)                                               DE15 460
      DOUBLE PRECISION Y,Z,A,B,C                                        DE15 470
C                                                                       DE15 480
C        TEST OF DIMENSION                                              DE15 490
      IF(NDIM-5)3,1,1                                                   DE15 500
C                                                                       DE15 510
C        PREPARE LOOP                                                   DE15 520
    1 A=Y(1)+Y(1)                                                       DE15 530
      C=Y(2)+Y(2)                                                       DE15 540
      B=.2D0*(A+Y(1)+C+Y(3)-Y(5))                                       DE15 550
      C=.1D0*(A+A+C+Y(2)+Y(3)+Y(3)+Y(4))                                DE15 560
C                                                                       DE15 570
C        START LOOP                                                     DE15 580
      DO 2 I=5,NDIM                                                     DE15 590
      A=B                                                               DE15 600
      B=C                                                               DE15 610
      C=.2D0*(Y(I-4)+Y(I-3)+Y(I-2)+Y(I-1)+Y(I))                         DE15 620
    2 Z(I-4)=A                                                          DE15 630
C        END OF LOOP                                                    DE15 640
C                                                                       DE15 650
C        UPDATE LAST FOUR COMPONENTS                                    DE15 660
      A=Y(NDIM)+Y(NDIM)                                                 DE15 670
     0A=.1D0*(A+A+Y(NDIM-1)+Y(NDIM-1)+Y(NDIM-1)+Y(NDIM-2)+Y(NDIM-2)     DE15 680
     1      +Y(NDIM-3))                                                 DE15 690
      Z(NDIM-3)=B                                                       DE15 700
      Z(NDIM-2)=C                                                       DE15 710
      Z(NDIM-1)=A                                                       DE15 720
      Z(NDIM)=A+A-C                                                     DE15 730
      IER=0                                                             DE15 740
      RETURN                                                            DE15 750
C                                                                       DE15 760
C        ERROR EXIT IN CASE NDIM IS LESS THAN 5                         DE15 770
    3 IER=-1                                                            DE15 780
      RETURN                                                            DE15 790
      END                                                               DE15 800
