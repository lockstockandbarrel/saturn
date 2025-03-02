C                                                                       GMTR  10
C     ..................................................................GMTR  20
C                                                                       GMTR  30
C        SUBROUTINE GMTRA                                               GMTR  40
C                                                                       GMTR  50
C        PURPOSE                                                        GMTR  60
C           TRANSPOSE A GENERAL MATRIX                                  GMTR  70
C                                                                       GMTR  80
C        USAGE                                                          GMTR  90
C           CALL GMTRA(A,R,N,M)                                         GMTR 100
C                                                                       GMTR 110
C        DESCRIPTION OF PARAMETERS                                      GMTR 120
C           A - NAME OF MATRIX TO BE TRANSPOSED                         GMTR 130
C           R - NAME OF RESULTANT MATRIX                                GMTR 140
C           N - NUMBER OF ROWS OF A AND COLUMNS OF R                    GMTR 150
C           M - NUMBER OF COLUMNS OF A AND ROWS OF R                    GMTR 160
C                                                                       GMTR 170
C        REMARKS                                                        GMTR 180
C           MATRIX R CANNOT BE IN THE SAME LOCATION AS MATRIX A         GMTR 190
C           MATRICES A AND R MUST BE STORED AS GENERAL MATRICES         GMTR 200
C                                                                       GMTR 210
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  GMTR 220
C           NONE                                                        GMTR 230
C                                                                       GMTR 240
C        METHOD                                                         GMTR 250
C           TRANSPOSE N BY M MATRIX A TO FORM M BY N MATRIX R           GMTR 260
C                                                                       GMTR 270
C     ..................................................................GMTR 280
C                                                                       GMTR 290
      SUBROUTINE GMTRA(A,R,N,M)                                         GMTR 300
      DIMENSION A(1),R(1)                                               GMTR 310
C                                                                       GMTR 320
      IR=0                                                              GMTR 330
      DO 10 I=1,N                                                       GMTR 340
      IJ=I-N                                                            GMTR 350
      DO 10 J=1,M                                                       GMTR 360
      IJ=IJ+N                                                           GMTR 370
      IR=IR+1                                                           GMTR 380
   10 R(IR)=A(IJ)                                                       GMTR 390
      RETURN                                                            GMTR 400
      END                                                               GMTR 410
