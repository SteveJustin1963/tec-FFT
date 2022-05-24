void bfft(float data[], int nn, int isign)

/* This is the Danielson and Lanczos implementation of the fast
   Fourier transform as described in Numerical Recipes, Press et
   al in section 12.2.  It has been tested by comparing with
   THE ORIGINAL COOLEY-TUKEY TRANSFORM, which is a fortran 4
   implementation of the same code. Bob Coldwell converted the
   Fortran to C.  The data is assumed to have real values
   in data[0],data[2],data[4],... and imaginary value in
   data[1],data[3],... This is not the same as Press's four1 which starts
   with a real in data[1]
   transform[k]=sum(data(j)*exp(isign*2*pi*sqrt(-1)*j*k/nn))
   The sum is over all j from 0 to nn-1
   nn must be a power of 2.  The dimension of the data is 2*nn, that is
   it has nn real values and nn imaginary values.  The transform is
   placed in data repalcing the original values. */
{double wr,wi,wpr,wpi,wtemp,theta;
 float tempr,tempi;
 int m,n,i,j,istep,mmax;
 n=2*nn;
 j=1;
 for (i=1;i<n;i+=2)
   {if(j > i)
     {tempr=data[j-1];
      tempi=data[j];
      data[j-1]=data[i-1];
      data[j]=data[i];
      data[i-1]=tempr;
      data[i]=tempi;}
    m=n/2;
    was1:
    if(m >= 2 && j > m)
      {j=j-m;
       m=m/2;
       goto was1;}
    j=j+m;}
/* Here begins the Danielson-Lanczos section (outer loop executed    FFT00380
   Log2 (NN) times */
  mmax=2;
  while (n > mmax)
    {istep=2*mmax;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j-1]-wi*data[j];
				tempi=wr*data[j]+wi*data[j-1];
				data[j-1]=data[i-1]-tempr;
				data[j]=data[i]-tempi;
				data[i-1] += tempr;
				data[i] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
