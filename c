
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXSIZE 1000

typedef struct {
    double real;
    double imag;
} complex;

void FFT(complex buf[], int n, int step)
{
    if (step < n) {
        FFT(buf, n, step * 2);
        FFT(buf + step, n, step * 2);
 
        for (int i = 0; i < n; i += 2 * step) {
            complex t = {
                .real = buf[i].real + buf[i + step].real,
                .imag = buf[i].imag + buf[i + step].imag
            };
            buf[i / 2] = t;
 
            complex s = {
                .real = buf[i].real - buf[i + step].real,
                .imag = buf[i].imag - buf[i + step].imag
            };
            buf[(i + n)/ 2] = s;
        }
    }
}
 
void show(complex buf[], int n)
{
    for (int i = 0; i < n; i++) {
        printf("%g %g\n", buf[i].real, buf[i].imag);
    }
}
 
int main()
{
    complex buf[MAXSIZE] = {0};
 
    for (int i = 0; i < MAXSIZE; i++) {
        buf[i].real = sin(2 * M_PI * i / MAXSIZE);
    }
 
    show(buf, MAXSIZE);
    FFT(buf, MAXSIZE, 1);
    show(buf, MAXSIZE);
 
    return 0;
}
