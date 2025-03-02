C                                                                       DA24  10
C     ..................................................................DA24  20
C                                                                       DA24  30
C        SUBROUTINE DQA24                                               DA24  40
C                                                                       DA24  50
C        PURPOSE                                                        DA24  60
C           TO COMPUTE INTEGRAL(EXP(-X)*FCT(X)/SQRT(X), SUMMED OVER X   DA24  70
C                               FROM 0 TO INFINITY).                    DA24  80
C                                                                       DA24  90
C        USAGE                                                          DA24 100
C           CALL DQA24 (FCT,Y)                                          DA24 110
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT                DA24 120
C                                                                       DA24 130
C        DESCRIPTION OF PARAMETERS                                      DA24 140
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION  DA24 150
C                    SUBPROGRAM USED.                                   DA24 160
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.     DA24 170
C                                                                       DA24 180
C        REMARKS                                                        DA24 190
C           NONE                                                        DA24 200
C                                                                       DA24 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DA24 220
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)    DA24 230
C           MUST BE FURNISHED BY THE USER.                              DA24 240
C                                                                       DA24 250
C        METHOD                                                         DA24 260
C           EVALUATION IS DONE BY MEANS OF 24-POINT GENERALIZED GAUSS-  DA24 270
C           LAGUERRE QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY       DA24 280
C           WHENEVER FCT(X) IS A POLYNOMIAL UP TO DEGREE 47.            DA24 290
C           FOR REFERENCE, SEE                                          DA24 300
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF    DA24 310
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED     DA24 320
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT       DA24 330
C           TR00.1100 (MARCH 1964), PP.15-16.                           DA24 340
C                                                                       DA24 350
C     ..................................................................DA24 360
C                                                                       DA24 370
      SUBROUTINE DQA24(FCT,Y)                                           DA24 380
C                                                                       DA24 390
C                                                                       DA24 400
      DOUBLE PRECISION X,Y,FCT                                          DA24 410
C                                                                       DA24 420
      X=.8055628081995041D2                                             DA24 430
      Y=.15871102921547994D-34*FCT(X)                                   DA24 440
      X=.69068601975304369D2                                            DA24 450
      Y=Y+.11969225386627757D-29*FCT(X)                                 DA24 460
      X=.60206666963057223D2                                            DA24 470
      Y=Y+.7370072160301340D-26*FCT(X)                                  DA24 480
      X=.52795432527283630D2                                            DA24 490
      Y=Y+.11129154937804570D-22*FCT(X)                                 DA24 500
      X=.46376979557540133D2                                            DA24 510
      Y=Y+.63767746470102769D-20*FCT(X)                                 DA24 520
      X=.40711598185543107D2                                            DA24 530
      Y=Y+.17460319202373353D-17*FCT(X)                                 DA24 540
      X=.35653703516328212D2                                            DA24 550
      Y=Y+.26303192453168170D-15*FCT(X)                                 DA24 560
      X=.31106464709046565D2                                            DA24 570
      Y=Y+.23951797309583587D-13*FCT(X)                                 DA24 580
      X=.27001406056472356D2                                            DA24 590
      Y=Y+.14093865163091778D-11*FCT(X)                                 DA24 600
      X=.23287932824879917D2                                            DA24 610
      Y=Y+.56305930756763382D-10*FCT(X)                                 DA24 620
      X=.19927425875242462D2                                            DA24 630
      Y=Y+.15860934990330765D-8*FCT(X)                                  DA24 640
      X=.16889671928527108D2                                            DA24 650
      Y=Y+.32450282717915397D-7*FCT(X)                                  DA24 660
      X=.14150586187285759D2                                            DA24 670
      Y=Y+.49373179873395010D-6*FCT(X)                                  DA24 680
      X=.11690695926056073D2                                            DA24 690
      Y=Y+.56945173834696962D-5*FCT(X)                                  DA24 700
      X=.9494095330026488D1                                             DA24 710
      Y=Y+.50571980554969778D-4*FCT(X)                                  DA24 720
      X=.7547704680023454D1                                             DA24 730
      Y=Y+.35030086360234566D-3*FCT(X)                                  DA24 740
      X=.58407332713236080D1                                            DA24 750
      Y=Y+.19127846396388306D-2*FCT(X)                                  DA24 760
      X=.43642830769353062D1                                            DA24 770
      Y=Y+.8306009823955105D-2*FCT(X)                                   DA24 780
      X=.31110524551477130D1                                            DA24 790
      Y=Y+.28889923149962199D-1*FCT(X)                                  DA24 800
      X=.20751129098523806D1                                            DA24 810
      Y=Y+.8095935396920770D-1*FCT(X)                                   DA24 820
      X=.12517406323627464D1                                            DA24 830
      Y=Y+.18364459415857036D0*FCT(X)                                   DA24 840
      X=.63729027873266879D0                                            DA24 850
      Y=Y+.33840894389128221D0*FCT(X)                                   DA24 860
      X=.22910231649262433D0                                            DA24 870
      Y=Y+.50792308532951820D0*FCT(X)                                   DA24 880
      X=.25437996585689359D-1                                           DA24 890
      Y=Y+.62200206075592616D0*FCT(X)                                   DA24 900
      RETURN                                                            DA24 910
      END                                                               DA24 920
