C                                                                       SE15  10
C     ..................................................................SE15  20
C                                                                       SE15  30
C        SUBROUTINE SE15                                                SE15  40
C                                                                       SE15  50
C        PURPOSE                                                        SE15  60
C           TO COMPUTE A VECTOR OF SMOOTHED FUNCTION VALUES GIVEN A     SE15  70
C           VECTOR OF FUNCTION VALUES WHOSE ENTRIES CORRESPOND TO       SE15  80
C           EQUIDISTANTLY SPACED ARGUMENT VALUES.                       SE15  90
C                                                                       SE15 100
C        USAGE                                                          SE15 110
C           CALL SE15(Y,Z,NDIM,IER)                                     SE15 120
C                                                                       SE15 130
C        DESCRIPTION OF PARAMETERS                                      SE15 140
C           Y     -  GIVEN VECTOR OF FUNCTION VALUES (DIMENSION NDIM)   SE15 150
C           Z     -  RESULTING VECTOR OF SMOOTHED FUNCTION VALUES       SE15 160
C                    (DIMENSION NDIM)                                   SE15 170
C           NDIM  -  DIMENSION OF VECTORS Y AND Z                       SE15 180
C           IER   -  RESULTING ERROR PARAMETER                          SE15 190
C                    IER = -1  - NDIM IS LESS THAN 5                    SE15 200
C                    IER =  0  - NO ERROR                               SE15 210
C                                                                       SE15 220
C        REMARKS                                                        SE15 230
C           (1)  IF IER=-1 THERE HAS BEEN NO COMPUTATION.               SE15 240
C           (2)   Z CAN HAVE THE SAME STORAGE ALLOCATION AS Y.  IF Y IS SE15 250
C                 DISTINCT FROM Z, THEN IT IS NOT DESTROYED.            SE15 260
C                                                                       SE15 270
C        SUBROUTINE AND FUNCTION SUBPROGRAMS REQUIRED                   SE15 280
C           NONE                                                        SE15 290
C                                                                       SE15 300
C        METHOD                                                         SE15 310
C           IF X IS THE (SUPPRESSED) VECTOR OF ARGUMENT VALUES, THEN    SE15 320
C           EXCEPT AT THE POINTS X(1),X(2),X(NDIM-1) AND X(NDIM), EACH  SE15 330
C           SMOOTHED VALUE Z(I) IS OBTAINED BY EVALUATING AT X(I) THE   SE15 340
C           LEAST-SQUARES POLYNOMIAL OF DEGREE 1 RELEVANT TO THE 5      SE15 350
C           SUCCESSIVE POINTS (X(I+K),Y(I+K)) K = -2,-1,...,2.  (SEE    SE15 360
C           HILDEBRAND, F.B., INTRODUCTION TO NUMERICAL ANALYSIS,       SE15 370
C           MC GRAW-HILL, NEW YORK/TORONTO/LONDON, 1956, PP. 295-302.)  SE15 380
C                                                                       SE15 390
C     ..................................................................SE15 400
C                                                                       SE15 410
      SUBROUTINE SE15(Y,Z,NDIM,IER)                                     SE15 420
C                                                                       SE15 430
C                                                                       SE15 440
      DIMENSION Y(1),Z(1)                                               SE15 450
C                                                                       SE15 460
C        TEST OF DIMENSION                                              SE15 470
      IF(NDIM-5)3,1,1                                                   SE15 480
C                                                                       SE15 490
C        PREPARE LOOP                                                   SE15 500
    1 A=Y(1)+Y(1)                                                       SE15 510
      C=Y(2)+Y(2)                                                       SE15 520
      B=.2*(A+Y(1)+C+Y(3)-Y(5))                                         SE15 530
      C=.1*(A+A+C+Y(2)+Y(3)+Y(3)+Y(4))                                  SE15 540
C                                                                       SE15 550
C        START LOOP                                                     SE15 560
      DO 2 I=5,NDIM                                                     SE15 570
      A=B                                                               SE15 580
      B=C                                                               SE15 590
      C=.2*(Y(I-4)+Y(I-3)+Y(I-2)+Y(I-1)+Y(I))                           SE15 600
    2 Z(I-4)=A                                                          SE15 610
C        END OF LOOP                                                    SE15 620
C                                                                       SE15 630
C        UPDATE LAST FOUR COMPONENTS                                    SE15 640
      A=Y(NDIM)+Y(NDIM)                                                 SE15 650
     0A=.1*(A+A+Y(NDIM-1)+Y(NDIM-1)+Y(NDIM-1)+Y(NDIM-2)+Y(NDIM-2)       SE15 660
     1      +Y(NDIM-3))                                                 SE15 670
      Z(NDIM-3)=B                                                       SE15 680
      Z(NDIM-2)=C                                                       SE15 690
      Z(NDIM-1)=A                                                       SE15 700
      Z(NDIM)=A+A-C                                                     SE15 710
      IER=0                                                             SE15 720
      RETURN                                                            SE15 730
C                                                                       SE15 740
C        ERROR EXIT IN CASE NDIM IS LESS THAN 5                         SE15 750
    3 IER=-1                                                            SE15 760
      RETURN                                                            SE15 770
      END                                                               SE15 780
