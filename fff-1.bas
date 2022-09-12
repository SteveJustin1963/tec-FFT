5000 //COMPLEX DFT BY CORRELATION 
5010 //Upon entry, N% contains the number of points in the DFT, and 
5020 //XR[ ] and XI[ ] contain the real and imaginary parts of the time domain. 
5030 //Upon return, REX[ ] and IMX[ ] contain the frequency domain data 
5040 //All signals run from 0 to N%-1. 
5050 
5060 PI = 3.14159265 //Set constants 
5070 
5080 FOR K% = 0 TO N%-1    //Zero REX[ ] and IMX[ ], so they can be used 
5090   REX[K%] = 0           //as accumulators during the correlation 
5100   IMX[K%] = 0 
5110 NEXT K%
5120 
5130 FOR K% = 0 TO N%-1     //Loop for each value in frequency domain 
5140   FOR I% = 0 TO N%-1     //Correlate with the complex sinusoid, SR & SI 

5150 
5160     SR = COS(2*PI*K%I%/N%)   //Calculate complex sinusoid 
5170     SI = -S1N(2*PI*K%I%/N%) 
5180     REX[K%] = REX[K%] + XR[I%]*SR - XI[I%]*SI 
5190     IMX[K%] = IMX[K%] + XR[I%]*SI + XI[I%]*SR 
5200 
5210   NEXT I% 
5220 NEXT K% 
5230 
5240 RETURN 


//


SUBROUTINE FFT(X,M) 
COMPLEX X(4096),U,S,T 
PI=3.14159265 
N=2**M 
DO 20 L=1,M 
LE=2**(M+1-L) 
LE2=LFJ2 
U=(1.0,0.0) 
S—CMPLX(COS(PI/FLOAT(LE2)),-SIN(PI/FLOAT(LE2))) 
DO 20 J=1,LE2 
DO 10 I=J,N,LE 
IP=I+LE2 
T=X(I)+X(IP) 
X(IP)=(X(1)-X(1P))*U 
10 X(I)=T 
20 U=U*S 
ND2=N/2 
NM1=N-1 
J-1 
DO 50 I=1,NM1 
IF(I.GE.J) GO TO 30 
T=X(J) 
X(J)=X(I) 
X(I)=T 
30 K=ND2 
40 IF(K.GE J) GO TO 50 
J=J-K 
K=K/2 
GO TO 40 
50 J=J+K 
RETURN 
END 

========
convert to basic 

5 DIM X(4096),U,S,T
10 PI=3.14159265
20 N=2^M
30 FOR L=1 TO M
40 LE=2^(M+1-L)
50 LE2=LE/2
60 U=(1.0,0.0)
70 S=CMPLX(COS(PI/FLOAT(LE2)),-SIN(PI/FLOAT(LE2)))
80 FOR J=1 TO LE2
90 FOR I=J TO N STEP LE
100 IP=I+LE2
110 T=X(I)+X(IP)
120 X(IP)=(X(I)-X(IP))*U
130 X(I)=T
140 NEXT I
150 NEXT J
160 NEXT L
170 ND2=N/2
180 NM1=N-1
190 J=1
200 FOR I=1 TO NM1
210 IF(I>=J) GOTO 240
220 T=X(J)
230 X(J)=X(I)
240 K=ND2
250 IF(K>=J) GOTO 290
260 J=J-K
270 K=K/2
280 GOTO 250
290 J=J+K
300 NEXT I
RETURN
====================











//TABLE 12-3 The Fast Fourier Transform in FORTRAN. Data are passed to this subroutine in the variables X( ) and M The integer, M, is the base two logarithm of the length of the DFT, i.e., M = 8 for a 256 point DFT, M = 12 for a 4096 point DFT, etc. The complex away, X( ), holds the time domain data upon entering the DFT. Upon return from this subroutine, X() is overwritten with the frequency domain data. Take note: this subroutine requires that the input and output signals run fromX(1) through X(N), rather than the customary X(0) through X(N-1). 


1000 THE FAST FOURIER TRANSFORM 
1010 Upon entry, N% contains the number of points in the DFT, REX[ ] and 
1020 IMX[ ] contain the real and imaginary parts of the input. Upon return,
1030 REX[ ] and IMX[ ] contain the DFT output. All signals run from 0 to N%-1. 
1040 
1050 PI = 3.14159265		//Set constants 
1060 NM1% = N%-1 
1070 ND2% = N%/2 
1080 M% = CINT(LOG(N%)/LOG(2)) 
1090 J% = ND2% 
1100 
1110 FOR I% =1 TO N%-2 		//Bit reversal sorting 
1120 IF I% >= J% THEN GOTO 1190 
1130 TR = REX[J%] 
1140 TI = IMX[J%] 
1150 REX[J%] = REX[I%] 
1160 IMX[J%] = IMX[I%] 
1170 REX[I%] = TR 
1180 IMX[I%] = TI 
1190 K% = ND2% 
1200 IF K% > J% THEN GOTO 1240 
1210 J% = J%-K% 
1220 K% = K%/2 
1230 GOTO 1200 
1240 J% = J%+K% 
1250 NEXT I% 
1260
1270 FOR L% =1 TO M% 	//Loop for each stage 
1280 LE% = CINT(2^L%) 
1290 LE2% = LE%/2 
1300 UR = 1 
1310 UI = 0 
1320 SR = COS(PI/LE2%) 	//Calculate sine & cosine values 
1330 SI = -SIN(PI/LE2%) 
1340 FOR J% = 1 TO LE2% //Loop for each sub DFT
1350   JM1% = J%-1 
1360   FOR I% = JM1% TO NM1% STEP LE% 	//Loop for each butterfly
1370     IP% = I%+LE2% 
1380     TR = REX[IP%]*UR - IMX[IP%)*UI //Butterfly calculation 
1390     TI = REX[IP%]*UI + IMX[IP%]*UR 
1400     REX[IP%] = REX[I%]-TR 
1410     IMX[IP%] = IMX[I%]-TI 
1420     REX[I%] = REX[I%]+TR 
1430     IMX[I%] = IMX[I%]+TI 
1440   NEXT I% 
1450   TR = UR 
1460   UR = TR*SR - UI*SI 
1470   UI = TR*SI + UI*SR 
1480 NEXT J% 
1490 NEXT L% 
1500 
1510 RETURN 

//'INVERSE FAST FOURIER TRANSFORM SUBROUTINE 
//'Upon entry, N% contains the number of points in the IDFT, REX[ ] and 
//'IMX[ ] contain the real and imaginary parts of the complex frequency domain. 
//'Upon return, REX[ ] and LMX[ ] contain the complex time domain. 
//'All signals run from 0 to N%-1. 

2060 FOR K% = 0 TO N%-1 	//'Change the sign of IMX[ ] 
2070 IMX[K%] = -IMX[K%] 
2080 NEXT K% 
2100 GOSUB 1000 		//'Calculate forward FFT (Table 12-3) 
2120 FOR I% = 0 TO N%-1 	//'Divide the time domain by N% and 'change the sign of IlvfX[ ] 
2130 REX[I%] = REX[I%]/N% 
2140 IMX[I%] = -IMX[I%]/N% 
2150 NEXT 1% 
2160  
2170 RETURN 
 
==================


