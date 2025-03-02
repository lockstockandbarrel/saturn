C                                                                       DCPY  10
C     ..................................................................DCPY  20
C                                                                       DCPY  30
C        SUBROUTINE DCPY                                                DCPY  40
C                                                                       DCPY  50
C        PURPOSE                                                        DCPY  60
C           COPY DIAGONAL ELEMENTS OF A MATRIX INTO A VECTOR            DCPY  70
C                                                                       DCPY  80
C        USAGE                                                          DCPY  90
C           CALL DCPY (A,R,N,MS)                                        DCPY 100
C                                                                       DCPY 110
C        DESCRIPTION OF PARAMETERS                                      DCPY 120
C           A - NAME OF INPUT MATRIX                                    DCPY 130
C           R - NAME OF OUTPUT VECTOR OF LENGTH N                       DCPY 140
C           N - NUMBER OF ROWS AND COLUMNS IN MATRIX A                  DCPY 150
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX A         DCPY 160
C                  0 - GENERAL                                          DCPY 170
C                  1 - SYMMETRIC                                        DCPY 180
C                  2 - DIAGONAL                                         DCPY 190
C                                                                       DCPY 200
C        REMARKS                                                        DCPY 210
C           INPUT MATRIX MUST BE A SQUARE MATRIX                        DCPY 220
C                                                                       DCPY 230
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  DCPY 240
C           LOC                                                         DCPY 250
C                                                                       DCPY 260
C        METHOD                                                         DCPY 270
C           ELEMENTS ON DIAGONAL OF MATRIX ARE MOVED TO CORRESPONDING   DCPY 280
C           POSITIONS OF VECTOR R                                       DCPY 290
C                                                                       DCPY 300
C     ..................................................................DCPY 310
C                                                                       DCPY 320
      SUBROUTINE DCPY(A,R,N,MS)                                         DCPY 330
      DIMENSION A(1),R(1)                                               DCPY 340
C                                                                       DCPY 350
      DO 3 J=1,N                                                        DCPY 360
C                                                                       DCPY 370
C        LOCATE DIAGONAL ELEMENT FOR ANY MATRIX STORAGE MODE            DCPY 380
C                                                                       DCPY 390
      CALL LOC(J,J,IJ,N,N,MS)                                           DCPY 400
C                                                                       DCPY 410
C        MOVE DIAGONAL ELEMENT TO VECTOR R                              DCPY 420
C                                                                       DCPY 430
    3 R(J)=A(IJ)                                                        DCPY 440
      RETURN                                                            DCPY 450
      END                                                               DCPY 460
