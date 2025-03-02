C                                                                       GMSU  10
C     ..................................................................GMSU  20
C                                                                       GMSU  30
C        SUBROUTINE GMSUB                                               GMSU  40
C                                                                       GMSU  50
C        PURPOSE                                                        GMSU  60
C           SUBTRACT ONE GENERAL MATRIX FROM ANOTHER TO FORM RESULTANT  GMSU  70
C           MATRIX                                                      GMSU  80
C                                                                       GMSU  90
C        USAGE                                                          GMSU 100
C           CALL GMSUB(A,B,R,N,M)                                       GMSU 110
C                                                                       GMSU 120
C        DESCRIPTION OF PARAMETERS                                      GMSU 130
C           A - NAME OF FIRST INPUT MATRIX                              GMSU 140
C           B - NAME OF SECOND INPUT MATRIX                             GMSU 150
C           R - NAME OF OUTPUT MATRIX                                   GMSU 160
C           N - NUMBER OF ROWS IN A,B,R                                 GMSU 170
C           M - NUMBER OF COLUMNS IN A,B,R                              GMSU 180
C                                                                       GMSU 190
C        REMARKS                                                        GMSU 200
C           ALL MATRICES MUST BE STORED AS GENERAL MATRICES             GMSU 210
C                                                                       GMSU 220
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  GMSU 230
C           NONE                                                        GMSU 240
C                                                                       GMSU 250
C        METHOD                                                         GMSU 260
C           MATRIX B ELEMENTS ARE SUBTRACTED FROM CORRESPONDING MATRIX AGMSU 270
C           ELEMENTS                                                    GMSU 280
C                                                                       GMSU 290
C     ..................................................................GMSU 300
C                                                                       GMSU 310
      SUBROUTINE GMSUB(A,B,R,N,M)                                       GMSU 320
      DIMENSION A(1),B(1),R(1)                                          GMSU 330
C                                                                       GMSU 340
C        CALCULATE NUMBER OF ELEMENTS                                   GMSU 350
C                                                                       GMSU 360
      NM=N*M                                                            GMSU 370
C                                                                       GMSU 380
C        SUBTRACT MATRICES                                              GMSU 390
C                                                                       GMSU 400
      DO 10 I=1,NM                                                      GMSU 410
   10 R(I)=A(I)-B(I)                                                    GMSU 420
      RETURN                                                            GMSU 430
      END                                                               GMSU 440
