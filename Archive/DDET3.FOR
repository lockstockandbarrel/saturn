C                                                                       DDT3  10
C     ..................................................................DDT3  20
C                                                                       DDT3  30
C        SUBROUTINE DDET3                                               DDT3  40
C                                                                       DDT3  50
C        PURPOSE                                                        DDT3  60
C           TO COMPUTE A VECTOR OF DERIVATIVE VALUES GIVEN A VECTOR OF  DDT3  70
C           FUNCTION VALUES WHOSE ENTRIES CORRESPOND TO EQUIDISTANTLY   DDT3  80
C           SPACED ARGUMENT VALUES.                                     DDT3  90
C                                                                       DDT3 100
C        USAGE                                                          DDT3 110
C           CALL DDET3(H,Y,Z,NDIM,IER)                                  DDT3 120
C                                                                       DDT3 130
C        DESCRIPTION OF PARAMETERS                                      DDT3 140
C           H     -  DOUBLE PRECISION CONSTANT DIFFERENCE BETWEEN       DDT3 150
C                    SUCCESSIVE ARGUMENT VALUES (H IS POSITIVE IF THE   DDT3 160
C                    ARGUMENT VALUES INCREASE AND NEGATIVE OTHERWISE)   DDT3 170
C           Y     -  GIVEN VECTOR OF DOUBLE PRECISION FUNCTION VALUES   DDT3 180
C                    (DIMENSION NDIM)                                   DDT3 190
C           Z     -  RESULTING VECTOR OF DOUBLE PRECISION DERIVATIVE    DDT3 200
C                    VALUES (DIMENSION NDIM)                            DDT3 210
C           NDIM  -  DIMENSION OF VECTORS Y AND Z                       DDT3 220
C           IER   -  RESULTING ERROR PARAMETER                          DDT3 230
C                    IER = -1  - NDIM IS LESS THAN 3                    DDT3 240
C                    IER =  0  - NO ERROR                               DDT3 250
C                    IER =  1  - H = 0                                  DDT3 260
C                                                                       DDT3 270
C        REMARKS                                                        DDT3 280
C           (1)   IF IER = -1,1, THEN THERE IS NO COMPUTATION.          DDT3 290
C           (2)   Z CAN HAVE THE SAME STORAGE ALLOCATION AS Y. IF Y IS  DDT3 300
C                 DISTINCT FROM Z, THEN IT IS NOT DESTROYED.            DDT3 310
C                                                                       DDT3 320
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DDT3 330
C           NONE                                                        DDT3 340
C                                                                       DDT3 350
C        METHOD                                                         DDT3 360
C           IF X IS THE (SUPPRESSED) VECTOR OF ARGUMENT VALUES, THEN    DDT3 370
C           EXCEPT AT THE ENDPOINTS X(1) AND X(NDIM), Z(I) IS THE       DDT3 380
C           DERIVATIVE AT X(I) OF THE LAGRANGIAN INTERPOLATION          DDT3 390
C           POLYNOMIAL OF DEGREE 2 RELEVANT TO THE 3 SUCCESSIVE POINTS  DDT3 400
C           (X(I+K),Y(I+K)) K = -1,0,1.  (SEE HILDEBRAND, F.B.,         DDT3 410
C           INTRODUCTION TO NUMERICAL ANALYSIS, MC-GRAW-HILL, NEW YORK/ DDT3 420
C           TORONTO/LONDON, 1956, PP.82-84.)                            DDT3 430
C                                                                       DDT3 440
C     ..................................................................DDT3 450
C                                                                       DDT3 460
      SUBROUTINE DDET3(H,Y,Z,NDIM,IER)                                  DDT3 470
C                                                                       DDT3 480
C                                                                       DDT3 490
      DIMENSION Y(1),Z(1)                                               DDT3 500
      DOUBLE PRECISION H,Y,Z,HH,YY,A,B                                  DDT3 510
C                                                                       DDT3 520
C        TEST OF DIMENSION                                              DDT3 530
      IF(NDIM-3)4,1,1                                                   DDT3 540
C                                                                       DDT3 550
C        TEST OF STEPSIZE                                               DDT3 560
    1 IF(H)2,5,2                                                        DDT3 570
C                                                                       DDT3 580
C        PREPARE DIFFERENTIATION LOOP                                   DDT3 590
    2 HH=.5D0/H                                                         DDT3 600
      YY=Y(NDIM-2)                                                      DDT3 610
      B=Y(2)+Y(2)                                                       DDT3 620
      B=HH*(B+B-Y(3)-Y(1)-Y(1)-Y(1))                                    DDT3 630
C                                                                       DDT3 640
C        START DIFFERENTIATION LOOP                                     DDT3 650
      DO 3 I=3,NDIM                                                     DDT3 660
      A=B                                                               DDT3 670
      B=HH*(Y(I)-Y(I-2))                                                DDT3 680
    3 Z(I-2)=A                                                          DDT3 690
C        END OF DIFFERENTIATION LOOP                                    DDT3 700
C                                                                       DDT3 710
C        NORMAL EXIT                                                    DDT3 720
      IER=0                                                             DDT3 730
      A=Y(NDIM-1)+Y(NDIM-1)                                             DDT3 740
      Z(NDIM)=HH*(Y(NDIM)+Y(NDIM)+Y(NDIM)-A-A+YY)                       DDT3 750
      Z(NDIM-1)=B                                                       DDT3 760
      RETURN                                                            DDT3 770
C                                                                       DDT3 780
C        ERROR EXIT IN CASE NDIM IS LESS THAN 3                         DDT3 790
    4 IER=-1                                                            DDT3 800
      RETURN                                                            DDT3 810
C                                                                       DDT3 820
C        ERROR EXIT IN CASE OF ZERO STEPSIZE                            DDT3 830
    5 IER=1                                                             DDT3 840
      RETURN                                                            DDT3 850
      END                                                               DDT3 860
