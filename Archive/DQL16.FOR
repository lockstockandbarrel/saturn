C                                                                       DL16  10
C     ..................................................................DL16  20
C                                                                       DL16  30
C        SUBROUTINE DQL16                                               DL16  40
C                                                                       DL16  50
C        PURPOSE                                                        DL16  60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X), SUMMED OVER X           DL16  70
C                               FROM 0 TO INFINITY).                    DL16  80
C                                                                       DL16  90
C        USAGE                                                          DL16 100
C           CALL DQL16 (FCT,Y)                                          DL16 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DL16 120
C                                                                       DL16 130
C        DESCRIPTION OF PARAMETERS                                      DL16 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DL16 150
C                    SUBPROGRAM USED.                                   DL16 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DL16 170
C                                                                       DL16 180
C        REMARKS                                                        DL16 190
C           NONE                                                        DL16 200
C                                                                       DL16 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DL16 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DL16 230
C           MUST BE FURNISHED BY THE USER.                              DL16 240
C                                                                       DL16 250
C        METHOD                                                         DL16 260
C           EVALUATION IS DONE BY MEANS OF 16-POINT GAUSSIAN-LAGUERRE   DL16 270
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY,               DL16 280
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 31.            DL16 290
C           FOR REFERENCE, SEE                                          DL16 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DL16 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DL16 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DL16 330
C           TR00.1100 (MARCH 1964), PP.24-25.                           DL16 340
C                                                                       DL16 350
C     ..................................................................DL16 360
C                                                                       DL16 370
      SUBROUTINE DQL16(FCT,Y)                                           DL16 380
C                                                                       DL16 390
C                                                                       DL16 400
      DOUBLE PRECISION X,Y,FCT                                          DL16 410
C                                                                       DL16 420
      X=.51701160339543318D2                                            DL16 430
      Y=.41614623703728552D-21*FCT(X)                                   DL16 440
      X=.41940452647688333D2                                            DL16 450
      Y=Y+.50504737000355128D-17*FCT(X)                                 DL16 460
      X=.34583398702286626D2                                            DL16 470
      Y=Y+.62979670025178678D-14*FCT(X)                                 DL16 480
      X=.28578729742882140D2                                            DL16 490
      Y=Y+.21270790332241030D-11*FCT(X)                                 DL16 500
      X=.23515905693991909D2                                            DL16 510
      Y=Y+.28623502429738816D-9*FCT(X)                                  DL16 520
      X=.19180156856753135D2                                            DL16 530
      Y=Y+.18810248410796732D-7*FCT(X)                                  DL16 540
      X=.15441527368781617D2                                            DL16 550
      Y=Y+.68283193308711996D-6*FCT(X)                                  DL16 560
      X=.12214223368866159D2                                            DL16 570
      Y=Y+.14844586873981299D-4*FCT(X)                                  DL16 580
      X=.9438314336391939D1                                             DL16 590
      Y=Y+.20427191530827846D-3*FCT(X)                                  DL16 600
      X=.70703385350482341D1                                            DL16 610
      Y=Y+.18490709435263109D-2*FCT(X)                                  DL16 620
      X=.50780186145497679D1                                            DL16 630
      Y=Y+.11299900080339453D-1*FCT(X)                                  DL16 640
      X=.34370866338932066D1                                            DL16 650
      Y=Y+.47328928694125219D-1*FCT(X)                                  DL16 660
      X=.21292836450983806D1                                            DL16 670
      Y=Y+.13629693429637754D0*FCT(X)                                   DL16 680
      X=.11410577748312269D1                                            DL16 690
      Y=Y+.26579577764421415D0*FCT(X)                                   DL16 700
      X=.46269632891508083D0                                            DL16 710
      Y=Y+.33105785495088417D0*FCT(X)                                   DL16 720
      X=.8764941047892784D-1                                            DL16 730
      Y=Y+.20615171495780099D0*FCT(X)                                   DL16 740
      RETURN                                                            DL16 750
      END                                                               DL16 760
