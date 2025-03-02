C                                                                       DA32  10
C     ..................................................................DA32  20
C                                                                       DA32  30
C        SUBROUTINE DQA32                                               DA32  40
C                                                                       DA32  50
C        PURPOSE                                                        DA32  60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X)/SQRT(X), SUMMED OVER X   DA32  70
C                               FROM 0 TO INFINITY).                    DA32  80
C                                                                       DA32  90
C        USAGE                                                          DA32 100
C           CALL DQA32 (FCT,Y)                                          DA32 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DA32 120
C                                                                       DA32 130
C        DESCRIPTION OF PARAMETERS                                      DA32 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DA32 150
C                    SUBPROGRAM USED.                                   DA32 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DA32 170
C                                                                       DA32 180
C        REMARKS                                                        DA32 190
C           NONE                                                        DA32 200
C                                                                       DA32 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DA32 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DA32 230
C           MUST BE FURNISHED BY THE USER.                              DA32 240
C                                                                       DA32 250
C        METHOD                                                         DA32 260
C           EVALUATION IS DONE BY MEANS OF 32-POINT GENERALIZED GAUSS-  DA32 270
C           LAGUERRE QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY       DA32 280
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 63.            DA32 290
C           FOR REFERENCE, SEE                                          DA32 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DA32 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DA32 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DA32 330
C           TR00.1100 (MARCH 1964), PP.15-16.                           DA32 340
C                                                                       DA32 350
C     ..................................................................DA32 360
C                                                                       DA32 370
      SUBROUTINE DQA32(FCT,Y)                                           DA32 380
C                                                                       DA32 390
C                                                                       DA32 400
      DOUBLE PRECISION X,Y,FCT                                          DA32 410
C                                                                       DA32 420
      X=.11079926894707576D3                                            DA32 430
      Y=.11071413071713886D-47*FCT(X)                                   DA32 440
      X=.9791671642606276D2                                             DA32 450
      Y=Y+.33594959802163184D-42*FCT(X)                                 DA32 460
      X=.8795611994313352D2                                             DA32 470
      Y=Y+.68422760225114810D-38*FCT(X)                                 DA32 480
      X=.7933909652982320D2                                             DA32 490
      Y=Y+.31147912492595276D-34*FCT(X)                                 DA32 500
      X=.71868499359551422D2                                            DA32 510
      Y=Y+.50993217982259985D-31*FCT(X)                                 DA32 520
      X=.65184426376135782D2                                            DA32 530
      Y=Y+.38582071909299337D-28*FCT(X)                                 DA32 540
      X=.59129027934391951D2                                            DA32 550
      Y=Y+.15723595577851821D-25*FCT(X)                                 DA32 560
      X=.53697231826148512D2                                            DA32 570
      Y=Y+.38234137666012957D-23*FCT(X)                                 DA32 580
      X=.48514583867416048D2                                            DA32 590
      Y=Y+.59657255685597023D-21*FCT(X)                                 DA32 600
      X=.43825336369903902D2                                            DA32 610
      Y=Y+.63045091330075628D-19*FCT(X)                                 DA32 620
      X=.39488797123368127D2                                            DA32 630
      Y=Y+.47037694213516382D-17*FCT(X)                                 DA32 640
      X=.35469961396173283D2                                            DA32 650
      Y=Y+.25601367826448761D-15*FCT(X)                                 DA32 660
      X=.31742543790616606D2                                            DA32 670
      Y=Y+.10437247453181695D-13*FCT(X)                                 DA32 680
      X=.28284583194970531D2                                            DA32 690
      Y=Y+.32566814614194407D-12*FCT(X)                                 DA32 700
      X=.25077856544198053D2                                            DA32 710
      Y=Y+.7918355533895448D-11*FCT(X)                                  DA32 720
      X=.22107070382206007D2                                            DA32 730
      Y=Y+.15230434500290903D-9*FCT(X)                                  DA32 740
      X=.19359271087268714D2                                            DA32 750
      Y=Y+.23472334846430987D-8*FCT(X)                                  DA32 760
      X=.16823405362953694D2                                            DA32 770
      Y=Y+.29302506329522187D-7*FCT(X)                                  DA32 780
      X=.14489986690780274D2                                            DA32 790
      Y=Y+.29910658734544941D-6*FCT(X)                                  DA32 800
      X=.12350838217714770D2                                            DA32 810
      Y=Y+.25166805020623692D-5*FCT(X)                                  DA32 820
      X=.10398891905552624D2                                            DA32 830
      Y=Y+.17576998461700718D-4*FCT(X)                                  DA32 840
      X=.8628029857405929D1                                             DA32 850
      Y=Y+.10251858271572549D-3*FCT(X)                                  DA32 860
      X=.70329577982838936D1                                            DA32 870
      Y=Y+.50196739702612497D-3*FCT(X)                                  DA32 880
      X=.56091034574961513D1                                            DA32 890
      Y=Y+.20726581990151553D-2*FCT(X)                                  DA32 900
      X=.43525345293301410D1                                            DA32 910
      Y=Y+.7245173957068918D-2*FCT(X)                                   DA32 920
      X=.32598922564569419D1                                            DA32 930
      Y=Y+.21512081019758274D-1*FCT(X)                                  DA32 940
      X=.23283376682103970D1                                            DA32 950
      Y=Y+.54406257907377837D-1*FCT(X)                                  DA32 960
      X=.15555082314789380D1                                            DA32 970
      Y=Y+.11747996392819887D0*FCT(X)                                   DA32 980
      X=.9394832145007343D0                                             DA32 990
      Y=Y+.21699669861237368D0*FCT(X)                                   DA321000
      X=.47875647727748885D0                                            DA321010
      Y=Y+.34337168469816740D0*FCT(X)                                   DA321020
      X=.17221572414539558D0                                            DA321030
      Y=Y+.46598957212535609D0*FCT(X)                                   DA321040
      X=.19127510968446856D-1                                           DA321050
      Y=Y+.54275484988260796D0*FCT(X)                                   DA321060
      RETURN                                                            DA321070
      END                                                               DA321080
