      SUBROUTINE FFT(DATA,NN,ISIGN)                                     FFT00010
C *** This is the Danielson and Lanczos implementation of the fast      FFT00020
C *** Fourier transform as described in Numerical Recipes, Press et     FFT00030
C *** al in section 12.2.  It has been tested by comparing with         FFT00040
C *** THE ORIGINAL COOLEY-TUKEY TRANSFORM, which is a fortran 4         FFT00050
C *** implementation of the same code.                                  FFT00060
C ***  TRANSFORM(K)=SUM(DATA(J)*EXP(ISIGN*                              FFT00070
C ***  2PI*SQRT(-1)*(J-1)*(K-1)/NN)). SUMMED OVER ALL J                 FFT00080
C *** AND K FROM 1 TO NN. DATA IS IN A ONE-DIMENSIONAL                  FFT00090
C *** COMPLEX ARRAY (I.E.,THE REAL AND IMAGINARY                        FFT00100
C *** PARTS ARE ADJACENT IN STORAGE ,SUCH AS FORTRAN IV                 FFT00110
C *** PLACES THEM) WHOSE LENGTH NN=2**K, K.GE.0 (IF NECESSARY           FFT00120
C *** APPEND ZEROES TO THE DATA). ISIGN IS +1 OR -1. IF A -1            FFT00130
C *** TRANSFORM IS FOLLOWED BY A +1 ONE (OR A +1 BY A -1) THE           FFT00140
C *** ORIGINAL DATA REAPPEAR, MULTIPLIED BY NN. TRANSFORM               FFT00150
C *** VALUES ARE RETURNED IN ARRAY DATA, REPLACING THE INPUT.           FFT00160
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA                                  FFT00170
      DIMENSION DATA(2*NN)                                              FFT00180
      N=2*NN                                                            FFT00190
      J=1                                                               FFT00200
      DO I=1,N,2                                                        FFT00210
        IF(J.GT.I)THEN                                                  FFT00220
          TEMPR=DATA(J)                                                 FFT00230
          TEMPI=DATA(J+1)                                               FFT00240
          DATA(J)=DATA(I)                                               FFT00250
          DATA(J+1)=DATA(I+1)                                           FFT00260
          DATA(I)=TEMPR                                                 FFT00270
          DATA(I+1)=TEMPI                                               FFT00280
        ENDIF                                                           FFT00290
        M=N/2                                                           FFT00300
1       IF((M.GE.2).AND.(J.GT.M))THEN                                   FFT00310
          J=J-M                                                         FFT00320
          M=M/2                                                         FFT00330
          GOTO 1                                                        FFT00340
        ENDIF                                                           FFT00350
        J=J+M                                                           FFT00360
      ENDDO                                                             FFT00370
C *** Here begins the Danielson-Lanczos section (outer loop executed    FFT00380
C *** Log2 (NN) times                                                   FFT00390
      MMAX=2                                                            FFT00400
2     IF(N.GT.MMAX)THEN                                                 FFT00410
        ISTEP=2*MMAX                                                    FFT00420
        THETA=6.28318530717959D0/(ISIGN*MMAX)                           FFT00430
        WPR=-2*DSIN(0.5D0*THETA)**2                                     FFT00440
        WPI=DSIN(THETA)                                                 FFT00450
        WR=1                                                            FFT00460
        WI=0                                                            FFT00470
        DO M=1,MMAX,2                                                   FFT00480
          DO I=M,N,ISTEP                                                FFT00490
            J=I+MMAX                                                    FFT00500
            TEMPR=WR*DATA(J)-WI*DATA(J+1)                               FFT00510
            TEMPI=WR*DATA(J+1)+WI*DATA(J)                               FFT00520
            DATA(J)=DATA(I)-TEMPR                                       FFT00530
            DATA(J+1)=DATA(I+1)-TEMPI                                   FFT00540
            DATA(I)=DATA(I)+TEMPR                                       FFT00550
            DATA(I+1)=DATA(I+1)+TEMPI                                   FFT00560
          ENDDO                                                         FFT00570
          WTEMP=WR                                                      FFT00580
          WR=WR*WPR-WI*WPI+WR                                           FFT00590
          WI=WI*WPR+WTEMP*WPI+WI                                        FFT00600
        ENDDO                                                           FFT00610
        MMAX=ISTEP                                                      FFT00620
        GOTO 2                                                          FFT00630
      ENDIF                                                             FFT00640
      RETURN                                                            FFT00650
      END                                                               FFT00660
