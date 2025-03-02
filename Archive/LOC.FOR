C                                                                       LOC   10
C     ..................................................................LOC   20
C                                                                       LOC   30
C        SUBROUTINE LOC                                                 LOC   40
C                                                                       LOC   50
C        PURPOSE                                                        LOC   60
C           COMPUTE A VECTOR SUBSCRIPT FOR AN ELEMENT IN A MATRIX OF    LOC   70
C           SPECIFIED STORAGE MODE                                      LOC   80
C                                                                       LOC   90
C        USAGE                                                          LOC  100
C           CALL LOC (I,J,IR,N,M,MS)                                    LOC  110
C                                                                       LOC  120
C        DESCRIPTION OF PARAMETERS                                      LOC  130
C           I   - ROW NUMBER OF ELEMENT                                 LOC  140
C           J   - COLUMN NUMBER  OF ELEMENT                             LOC  150
C           IR  - RESULTANT VECTOR SUBSCRIPT                            LOC  160
C           N   - NUMBER OF ROWS IN MATRIX                              LOC  170
C           M   - NUMBER OF COLUMNS IN MATRIX                           LOC  180
C           MS  - ONE DIGIT NUMBER FOR STORAGE MODE OF MATRIX           LOC  190
C                  0 - GENERAL                                          LOC  200
C                  1 - SYMMETRIC                                        LOC  210
C                  2 - DIAGONAL                                         LOC  220
C                                                                       LOC  230
C        REMARKS                                                        LOC  240
C           NONE                                                        LOC  250
C                                                                       LOC  260
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  LOC  270
C           NONE                                                        LOC  280
C                                                                       LOC  290
C        METHOD                                                         LOC  300
C           MS=0   SUBSCRIPT IS COMPUTED FOR A MATRIX WITH N*M ELEMENTS LOC  310
C                  IN STORAGE (GENERAL MATRIX)                          LOC  320
C           MS=1   SUBSCRIPT IS COMPUTED FOR A MATRIX WITH N*(N+1)/2 IN LOC  330
C                  STORAGE (UPPER TRIANGLE OF SYMMETRIC MATRIX). IF     LOC  340
C                  ELEMENT IS IN LOWER TRIANGULAR PORTION, SUBSCRIPT IS LOC  350
C                  CORRESPONDING ELEMENT IN UPPER TRIANGLE.             LOC  360
C           MS=2   SUBSCRIPT IS COMPUTED FOR A MATRIX WITH N ELEMENTS   LOC  370
C                  IN STORAGE (DIAGONAL ELEMENTS OF DIAGONAL MATRIX).   LOC  380
C                  IF ELEMENT IS NOT ON DIAGONAL (AND THEREFORE NOT IN  LOC  390
C                  STORAGE), IR IS SET TO ZERO.                         LOC  400
C                                                                       LOC  410
C     ..................................................................LOC  420
C                                                                       LOC  430
      SUBROUTINE LOC(I,J,IR,N,M,MS)                                     LOC  440
C                                                                       LOC  450
      IX=I                                                              LOC  460
      JX=J                                                              LOC  470
      IF(MS-1) 10,20,30                                                 LOC  480
   10 IRX=N*(JX-1)+IX                                                   LOC  490
      GO TO 36                                                          LOC  500
   20 IF(IX-JX) 22,24,24                                                LOC  510
   22 IRX=IX+(JX*JX-JX)/2                                               LOC  520
      GO TO 36                                                          LOC  530
   24 IRX=JX+(IX*IX-IX)/2                                               LOC  540
      GO TO 36                                                          LOC  550
   30 IRX=0                                                             LOC  560
      IF(IX-JX) 36,32,36                                                LOC  570
   32 IRX=IX                                                            LOC  580
   36 IR=IRX                                                            LOC  590
      RETURN                                                            LOC  600
      END                                                               LOC  610
