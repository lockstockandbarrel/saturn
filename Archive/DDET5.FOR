C                                                                       DDT5  10
C     ..................................................................DDT5  20
C                                                                       DDT5  30
C        SUBROUTINE DDET5                                               DDT5  40
C                                                                       DDT5  50
C        PURPOSE                                                        DDT5  60
C           TO COMPUTE A VECTOR OF DERIVATIVE VALUES GIVEN A VECTOR OF  DDT5  70
C           FUNCTION VALUES WHOSE ENTRIES CORRESPOND TO EQUIDISTANTLY   DDT5  80
C           SPACED ARGUMENT VALUES.                                     DDT5  90
C                                                                       DDT5 100
C        USAGE                                                          DDT5 110
C           CALL DDET5(H,Y,Z,NDIM,IER)                                  DDT5 120
C                                                                       DDT5 130
C        DESCRIPTION OF PARAMETERS                                      DDT5 140
C           H     -  DOUBLE PRECISION CONSTANT DIFFERENCE BETWEEN       DDT5 150
C                    SUCCESSIVE ARGUMENT VALUES (H IS POSITIVE IF THE   DDT5 160
C                    ARGUMENT VALUES INCREASE AND NEGATIVE OTHERWISE)   DDT5 170
C           Y     -  GIVEN VECTOR OF DOUBLE PRECISION FUNCTION VALUES   DDT5 180
C                    (DIMENSION NDIM)                                   DDT5 190
C           Z     -  RESULTING VECTOR OF DOUBLE PRECISION DERIVATIVE    DDT5 200
C                    VALUES (DIMENSION NDIM)                            DDT5 210
C           NDIM  -  DIMENSION OF VECTORS Y AND Z                       DDT5 220
C           IER   -  RESULTING ERROR PARAMETER                          DDT5 230
C                    IER = -1  - NDIM IS LESS THAN 5                    DDT5 240
C                    IER =  0  - NO ERROR                               DDT5 250
C                    IER =  1  - H = 0                                  DDT5 260
C                                                                       DDT5 270
C        REMARKS                                                        DDT5 280
C           (1)   IF IER = -1,1, THEN THERE IS NO COMPUTATION.          DDT5 290
C           (2)   Z CAN HAVE THE SAME STORAGE ALLOCATION AS Y.  IF Y IS DDT5 300
C                 DISTINCT FROM Z, THEN IT IS NOT DESTROYED.            DDT5 310
C                                                                       DDT5 320
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DDT5 330
C           NONE                                                        DDT5 340
C                                                                       DDT5 350
C        METHOD                                                         DDT5 360
C           IF X IS THE (SUPPRESSED) VECTOR OF ARGUMENT VALUES, THEN    DDT5 370
C           EXCEPT AT THE POINTS X(1),X(2),X(NDIM-1) AND X(NDIM), Z(I)  DDT5 380
C           IS THE DERIVATIVE AT X(I) OF THE LAGRANGIAN INTERPOLATION   DDT5 390
C           POLYNOMIAL OF DEGREE 4 RELEVANT TO THE 5 SUCCESSIVE POINTS  DDT5 400
C           (X(I+K),Y(I+K)) K = -2,-1,...,2.  (SEE HILDEBRAND, F.B.,    DDT5 410
C           INTRODUCTION TO NUMERICAL ANALYSIS, MC GRAW-HILL, NEW YORK/ DDT5 420
C           TORONTO/LONDON, 1956, PP. 82-84.)                           DDT5 430
C                                                                       DDT5 440
C     ..................................................................DDT5 450
C                                                                       DDT5 460
      SUBROUTINE DDET5(H,Y,Z,NDIM,IER)                                  DDT5 470
C                                                                       DDT5 480
C                                                                       DDT5 490
      DIMENSION Y(1),Z(1)                                               DDT5 500
      DOUBLE PRECISION H,Y,Z,HH,YY,A,B,C                                DDT5 510
C                                                                       DDT5 520
C        TEST OF DIMENSION                                              DDT5 530
      IF(NDIM-5)4,1,1                                                   DDT5 540
C                                                                       DDT5 550
C        TEST OF STEPSIZE                                               DDT5 560
    1 IF(H)2,5,2                                                        DDT5 570
C                                                                       DDT5 580
C        PREPARE DIFFERENTIATION LOOP                                   DDT5 590
    2 HH=.08333333333333333D0/H                                         DDT5 600
      YY=Y(NDIM-4)                                                      DDT5 610
      B=HH*(-25.D0*Y(1)+48.D0*Y(2)-36.D0*Y(3)+16.D0*Y(4)-3.D0*Y(5))     DDT5 620
      C=HH*(-3.D0*Y(1)-10.D0*Y(2)+18.D0*Y(3)-6.D0*Y(4)+Y(5))            DDT5 630
C                                                                       DDT5 640
C        START DIFFERENTIATION LOOP                                     DDT5 650
      DO 3 I=5,NDIM                                                     DDT5 660
      A=B                                                               DDT5 670
      B=C                                                               DDT5 680
      C=HH*(Y(I-4)-Y(I)+8.D0*(Y(I-1)-Y(I-3)))                           DDT5 690
    3 Z(I-4)=A                                                          DDT5 700
C        END OF DIFFERENTIATION LOOP                                    DDT5 710
C                                                                       DDT5 720
C        NORMAL EXIT                                                    DDT5 730
      IER=0                                                             DDT5 740
     0A=HH*(-YY+6.D0*Y(NDIM-3)-18.D0*Y(NDIM-2)+10.D0*Y(NDIM-1)          DDT5 750
     1      +3.D0*Y(NDIM))                                              DDT5 760
     0Z(NDIM)=HH*(3.D0*YY-16.D0*Y(NDIM-3)+36.D0*Y(NDIM-2)               DDT5 770
     1            -48.D0*Y(NDIM-1)+25.D0*Y(NDIM))                       DDT5 780
      Z(NDIM-1)=A                                                       DDT5 790
      Z(NDIM-2)=C                                                       DDT5 800
      Z(NDIM-3)=B                                                       DDT5 810
      RETURN                                                            DDT5 820
C                                                                       DDT5 830
C        ERROR EXIT IN CASE NDIM IS LESS THAN 5                         DDT5 840
    4 IER=-1                                                            DDT5 850
      RETURN                                                            DDT5 860
C                                                                       DDT5 870
C        ERROR EXIT IN CASE OF ZERO STEPSIZE                            DDT5 880
    5 IER=1                                                             DDT5 890
      RETURN                                                            DDT5 900
      END                                                               DDT5 910
