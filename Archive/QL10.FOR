C                                                                       QL10  10
C     ..................................................................QL10  20
C                                                                       QL10  30
C        SUBROUTINE QL10                                                QL10  40
C                                                                       QL10  50
C        PURPOSE                                                        QL10  60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X FROM 0    QL10  70
C                               TO INFINITY).                           QL10  80
C                                                                       QL10  90
C        USAGE                                                          QL10 100
C           CALL QL10(FCT,Y)                                            QL10 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                QL10 120
C                                                                       QL10 130
C        DESCRIPTION OF PARAMETERS                                      QL10 140
C           FCT    - THE NAME OF AN EXTERNAL FUNCTION SUBPROGRAM USED.  QL10 150
C           Y      - THE RESULTING INTEGRAL VALUE.                      QL10 160
C                                                                       QL10 170
C        REMARKS                                                        QL10 180
C           NONE                                                        QL10 190
C                                                                       QL10 200
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  QL10 210
C           THE EXTERNAL FUNCTION SUBPROGRAM FCT(X) MUST BE FURNISHED   QL10 220
C           BY THE USER.                                                QL10 230
C                                                                       QL10 240
C        METHOD                                                         QL10 250
C           EVALUATION IS DONE BY MEANS OF 10-POINT GAUSSIAN-LAGUERRE   QL10 260
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER       QL10 270
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 19.                     QL10 280
C           FOR REFERENCE, SEE                                          QL10 290
C           V.I.KRYLOV, APPROXIMATE CALCULATION OF INTEGRALS,           QL10 300
C           MACMILLAN, NEW YORK/LONDON, 1962, PP.130-132 AND 347-352.   QL10 310
C                                                                       QL10 320
C     ..................................................................QL10 330
C                                                                       QL10 340
      SUBROUTINE QL10(FCT,Y)                                            QL10 350
C                                                                       QL10 360
C                                                                       QL10 370
      X=29.92070                                                        QL10 380
      Y=.9911827E-12*FCT(X)                                             QL10 390
      X=21.99659                                                        QL10 400
      Y=Y+.1839565E-8*FCT(X)                                            QL10 410
      X=16.27926                                                        QL10 420
      Y=Y+.4249314E-6*FCT(X)                                            QL10 430
      X=11.84379                                                        QL10 440
      Y=Y+.2825923E-4*FCT(X)                                            QL10 450
      X=8.330153                                                        QL10 460
      Y=Y+.7530084E-3*FCT(X)                                            QL10 470
      X=5.552496                                                        QL10 480
      Y=Y+.009501517*FCT(X)                                             QL10 490
      X=3.401434                                                        QL10 500
      Y=Y+.06208746*FCT(X)                                              QL10 510
      X=1.808343                                                        QL10 520
      Y=Y+.2180683*FCT(X)                                               QL10 530
      X=.7294545                                                        QL10 540
      Y=Y+.4011199*FCT(X)                                               QL10 550
      X=.1377935                                                        QL10 560
      Y=Y+.3084411*FCT(X)                                               QL10 570
      RETURN                                                            QL10 580
      END                                                               QL10 590
