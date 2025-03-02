C                                                                       DG32  10
C     ..................................................................DG32  20
C                                                                       DG32  30
C        SUBROUTINE DQG32                                               DG32  40
C                                                                       DG32  50
C        PURPOSE                                                        DG32  60
C           TO COMPUTE INTEGRAL(FCT(X), SUMMED OVER X FROM XL TO XU)    DG32  70
C                                                                       DG32  80
C        USAGE                                                          DG32  90
C           CALL DQG32 (XL,XU,FCT,Y)                                    DG32 100
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DG32 110
C                                                                       DG32 120
C        DESCRIPTION OF PARAMETERS                                      DG32 130
C           XL     - DOUBLE PRECISION LOWER BOUND OF THE INTERVAL.      DG32 140
C           XU     - DOUBLE PRECISION UPPER BOUND OF THE INTERVAL.      DG32 150
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DG32 160
C                    SUBPROGRAM USED.                                   DG32 170
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DG32 180
C                                                                       DG32 190
C        REMARKS                                                        DG32 200
C           NONE                                                        DG32 210
C                                                                       DG32 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DG32 230
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DG32 240
C           MUST BE FURNISHED BY THE USER.                              DG32 250
C                                                                       DG32 260
C        METHOD                                                         DG32 270
C           EVALUATION IS DONE BY MEANS OF 32-POINT GAUSS QUADRATURE    DG32 280
C           FORMULA, WHICH INTEGRATES POLYNOMIALS UP TO DEGREE 63       DG32 290
C           EXACTLY. FOR REFERENCE, SEE                                 DG32 300
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           DG32 310
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.100-111 AND 337-340.   DG32 320
C                                                                       DG32 330
C     ..................................................................DG32 340
C                                                                       DG32 350
      SUBROUTINE DQG32(XL,XU,FCT,Y)                                     DG32 360
C                                                                       DG32 370
C                                                                       DG32 380
      DOUBLE PRECISION XL,XU,Y,A,B,C,FCT                                DG32 390
C                                                                       DG32 400
      A=.5D0*(XU+XL)                                                    DG32 410
      B=XU-XL                                                           DG32 420
      C=.49863193092474078D0*B                                          DG32 430
      Y=.35093050047350483D-2*(FCT(A+C)+FCT(A-C))                       DG32 440
      C=.49280575577263417D0*B                                          DG32 450
      Y=Y+.8137197365452835D-2*(FCT(A+C)+FCT(A-C))                      DG32 460
      C=.48238112779375322D0*B                                          DG32 470
      Y=Y+.12696032654631030D-1*(FCT(A+C)+FCT(A-C))                     DG32 480
      C=.46745303796886984D0*B                                          DG32 490
      Y=Y+.17136931456510717D-1*(FCT(A+C)+FCT(A-C))                     DG32 500
      C=.44816057788302606D0*B                                          DG32 510
      Y=Y+.21417949011113340D-1*(FCT(A+C)+FCT(A-C))                     DG32 520
      C=.42468380686628499D0*B                                          DG32 530
      Y=Y+.25499029631188088D-1*(FCT(A+C)+FCT(A-C))                     DG32 540
      C=.39724189798397120D0*B                                          DG32 550
      Y=Y+.29342046739267774D-1*(FCT(A+C)+FCT(A-C))                     DG32 560
      C=.36609105937014484D0*B                                          DG32 570
      Y=Y+.32911111388180923D-1*(FCT(A+C)+FCT(A-C))                     DG32 580
      C=.33152213346510760D0*B                                          DG32 590
      Y=Y+.36172897054424253D-1*(FCT(A+C)+FCT(A-C))                     DG32 600
      C=.29385787862038116D0*B                                          DG32 610
      Y=Y+.39096947893535153D-1*(FCT(A+C)+FCT(A-C))                     DG32 620
      C=.25344995446611470D0*B                                          DG32 630
      Y=Y+.41655962113473378D-1*(FCT(A+C)+FCT(A-C))                     DG32 640
      C=.21067563806531767D0*B                                          DG32 650
      Y=Y+.43826046502201906D-1*(FCT(A+C)+FCT(A-C))                     DG32 660
      C=.16593430114106382D0*B                                          DG32 670
      Y=Y+.45586939347881942D-1*(FCT(A+C)+FCT(A-C))                     DG32 680
      C=.11964368112606854D0*B                                          DG32 690
      Y=Y+.46922199540402283D-1*(FCT(A+C)+FCT(A-C))                     DG32 700
      C=.7223598079139825D-1*B                                          DG32 710
      Y=Y+.47819360039637430D-1*(FCT(A+C)+FCT(A-C))                     DG32 720
      C=.24153832843869158D-1*B                                         DG32 730
      Y=B*(Y+.48270044257363900D-1*(FCT(A+C)+FCT(A-C)))                 DG32 740
      RETURN                                                            DG32 750
      END                                                               DG32 760
