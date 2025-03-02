C                                                                       GMAD  10
C     ..................................................................GMAD  20
C                                                                       GMAD  30
C        SUBROUTINE GMADD                                               GMAD  40
C                                                                       GMAD  50
C        PURPOSE                                                        GMAD  60
C           ADD TWO GENERAL MATRICES TO FORM RESULTANT GENERAL MATRIX   GMAD  70
C                                                                       GMAD  80
C        USAGE                                                          GMAD  90
C           CALL GMADD(A,B,R,N,M)                                       GMAD 100
C                                                                       GMAD 110
C        DESCRIPTION OF PARAMETERS                                      GMAD 120
C           A - NAME OF FIRST INPUT MATRIX                              GMAD 130
C           B - NAME OF SECOND INPUT MATRIX                             GMAD 140
C           R - NAME OF OUTPUT MATRIX                                   GMAD 150
C           N - NUMBER OF ROWS IN A,B,R                                 GMAD 160
C           M - NUMBER OF COLUMNS IN A,B,R                              GMAD 170
C                                                                       GMAD 180
C        REMARKS                                                        GMAD 190
C           ALL MATRICES MUST BE STORED AS GENERAL MATRICES             GMAD 200
C                                                                       GMAD 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  GMAD 220
C           NONE                                                        GMAD 230
C                                                                       GMAD 240
C        METHOD                                                         GMAD 250
C           ADDITION IS PERFORMED ELEMENT BY ELEMENT                    GMAD 260
C                                                                       GMAD 270
C     ..................................................................GMAD 280
C                                                                       GMAD 290
      SUBROUTINE GMADD(A,B,R,N,M)                                       GMAD 300
      DIMENSION A(1),B(1),R(1)                                          GMAD 310
C                                                                       GMAD 320
C        CALCULATE NUMBER OF ELEMENTS                                   GMAD 330
C                                                                       GMAD 340
      NM=N*M                                                            GMAD 350
C                                                                       GMAD 360
C        ADD MATRICES                                                   GMAD 370
C                                                                       GMAD 380
      DO 10 I=1,NM                                                      GMAD 390
   10 R(I)=A(I)+B(I)                                                    GMAD 400
      RETURN                                                            GMAD 410
      END                                                               GMAD 420
