C                                                                       SICI  10
C     ..................................................................SICI  20
C                                                                       SICI  30
C        SUBROUTINE SICI                                                SICI  40
C                                                                       SICI  50
C        PURPOSE                                                        SICI  60
C           COMPUTES THE SINE AND COSINE INTEGRAL                       SICI  70
C                                                                       SICI  80
C        USAGE                                                          SICI  90
C           CALL SICI(SI,CI,X)                                          SICI 100
C                                                                       SICI 110
C        DESCRIPTION OF PARAMETERS                                      SICI 120
C           SI    - THE RESULTANT VALUE SI(X)                           SICI 130
C           CI    - THE RESULTANT VALUE CI(X)                           SICI 140
C           X     - THE ARGUMENT OF SI(X) AND CI(X)                     SICI 150
C                                                                       SICI 160
C        REMARKS                                                        SICI 170
C           THE ARGUMENT VALUE REMAINS UNCHANGED                        SICI 180
C                                                                       SICI 190
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  SICI 200
C           NONE                                                        SICI 210
C                                                                       SICI 220
C        METHOD                                                         SICI 230
C           DEFINITION                                                  SICI 240
C           SI(X)=INTEGRAL(SIN(T)/T)                                    SICI 250
C           CI(X)=INTEGRAL(COS(T)/T)                                    SICI 260
C           EVALUATION                                                  SICI 270
C           REDUCTION OF RANGE USING SYMMETRY.                          SICI 280
C           DIFFERENT APPROXIMATIONS ARE USED FOR ABS(X) GREATER        SICI 290
C           THAN 4 AND FOR ABS(X) LESS THAN 4.                          SICI 300
C           REFERENCE                                                   SICI 310
C           LUKE AND WIMP, 'POLYNOMIAL APPROXIMATIONS TO INTEGRAL       SICI 320
C           TRANSFORMS',  MATHEMATICAL TABLES AND OTHER AIDS TO         SICI 330
C           COMPUTATION, VOL. 15, 1961, ISSUE 74, PP. 174-178.          SICI 340
C                                                                       SICI 350
C     ..................................................................SICI 360
C                                                                       SICI 370
      SUBROUTINE SICI(SI,CI,X)                                          SICI 380
      Z=ABS(X)                                                          SICI 390
      IF(Z-4.)1,1,4                                                     SICI 400
    1 Y=(4.-Z)*(4.+Z)                                                   SICI 410
      SI=-1.570797E0                                                    SICI 420
      IF(Z)3,2,3                                                        SICI 430
    2 CI=-1.E75                                                         SICI 440
      RETURN                                                            SICI 450
    3 SI=X*(((((1.753141E-9*Y+1.568988E-7)*Y+1.374168E-5)*Y+6.939889E-4)SICI 460
     1*Y+1.964882E-2)*Y+4.395509E-1+SI/X)                               SICI 470
      CI=((5.772156E-1+ALOG(Z))/Z-Z*(((((1.386985E-10*Y+1.584996E-8)*Y  SICI 480
     1+1.725752E-6)*Y+1.185999E-4)*Y+4.990920E-3)*Y+1.315308E-1))*Z     SICI 490
      RETURN                                                            SICI 500
    4 SI=SIN(Z)                                                         SICI 510
      Y=COS(Z)                                                          SICI 520
      Z=4./Z                                                            SICI 530
      U=((((((((4.048069E-3*Z-2.279143E-2)*Z+5.515070E-2)*Z-7.261642E-2)SICI 540
     1*Z+4.987716E-2)*Z-3.332519E-3)*Z-2.314617E-2)*Z-1.134958E-5)*Z    SICI 550
     2+6.250011E-2)*Z+2.583989E-10                                      SICI 560
      V=(((((((((-5.108699E-3*Z+2.819179E-2)*Z-6.537283E-2)*Z           SICI 570
     1+7.902034E-2)*Z-4.400416E-2)*Z-7.945556E-3)*Z+2.601293E-2)*Z      SICI 580
     2-3.764000E-4)*Z-3.122418E-2)*Z-6.646441E-7)*Z+2.500000E-1         SICI 590
      CI=Z*(SI*V-Y*U)                                                   SICI 600
      SI=-Z*(SI*U+Y*V)                                                  SICI 610
      IF(X)5,6,6                                                        SICI 620
    5 SI=3.141593E0-SI                                                  SICI 630
    6 RETURN                                                            SICI 640
      END                                                               SICI 650
