//
////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace std;

const int N = 1024;
const float PI = 3.1416;


class wavelet
{
    public:
 	float real[N];
 	float imag[N];
 	
}W;

class cwt
{
    public:
 	float real[N];
 	float imag[N];
 	
}C;


class fftdata
{
    public:
 	float real[N];
 	float imag[N];
 	
}F;



//class tag Array()
//{
// 	for(int i=0;i<10;i++)
// {
//  X.a[i]=i+1;
//        cout<<X.a[i]<<" "<<endl;
// }
// cout<<endl;
// return X;
//}

inline void swap (float &a, float &b)
{
    float t;
    t = a;
    a = b;
    b = t;
}


wavelet  Wavelet_Morlet (float scale, float frequency, int N)
{
    // ���ٸ���Ҷ��任
    //float Morlet [N];
    int t;
    
    
    for (t=0; t<N; t++)
    
    {
    	float ts = t/scale;
    	W.real[t] = pow(PI,-0.25)*cos(frequency*ts)*exp((-1/2)*pow(ts,2));
    	W.imag[t] = pow(PI,-0.25)*sin(frequency*ts)*exp((-1/2)*pow(ts,2));
	}
	
	return W;
}




/*
float CWT_imag(float xreal [], float ximag [], float scale, float frequency, int N, float dt)
{
	
	float creal [], float cimag [], wavefft_real [], float wavefft_imag [];
	
	wavefft_real[] = FFT_real(Wavelet_Morlet(), N);
	wavefft_imag[] = FFT_imag(Wavelet_Morlet()), N;
	
	
	for (i=0, i<N, i++)
	{
		creal[i] = wavefft_real[i]*xreal[i] + wavefft_imag[i]*ximag[i];
		cimag[i] = wavefft_real[i]*ximag[i] - wavefft_imag[i]*xreal[i];
	}
		
	return IFFT_imag (creal, cimag, N);
}
*/


void bitrp (float xreal [], float ximag [], int n)
{
    // λ��ת�û� Bit-reversal Permutation
    int i, j, a, b, p;

    for (i = 1, p = 0; i < n; i *= 2)
        {
        p ++;
        }
    for (i = 0; i < n; i ++)
        {
        a = i;
        b = 0;
        for (j = 0; j < p; j ++)
            {
            b = (b << 1) + (a & 1);    // b = b * 2 + a % 2;
            a >>= 1;        // a = a / 2;
            }
        if ( b > i)
            {
            swap (xreal [i], xreal [b]);
            swap (ximag [i], ximag [b]);
            }
        }
}

fftdata FFT(float xreal [], float ximag [], int n)
{
    // ���ٸ���Ҷ�任�������� x �任���Ա����� x �У�xreal, ximag �ֱ��� x ��ʵ�����鲿
    float wreal [N / 2], wimag [N / 2], treal, timag, ureal, uimag, arg;
    int m, k, j, t, index1, index2;
    
    bitrp (xreal, ximag, n);

    // ���� 1 ��ǰ n / 2 �� n �η����Ĺ���� W'j = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
    arg = - 2 * PI / n;
    treal = cos (arg);
    timag = sin (arg);
    wreal [0] = 1.0;
    wimag [0] = 0.0;
    for (j = 1; j < n / 2; j ++)
        {
        wreal [j] = wreal [j - 1] * treal - wimag [j - 1] * timag;
        wimag [j] = wreal [j - 1] * timag + wimag [j - 1] * treal;
        }

    for (m = 2; m <= n; m *= 2)
        {
        for (k = 0; k < n; k += m)
            {
            for (j = 0; j < m / 2; j ++)
                {
                index1 = k + j;
                index2 = index1 + m / 2;
                t = n * j / m;    // ��ת���� w ��ʵ���� wreal [] �е��±�Ϊ t
                treal = wreal [t] * xreal [index2] - wimag [t] * ximag [index2];
                timag = wreal [t] * ximag [index2] + wimag [t] * xreal [index2];
                ureal = xreal [index1];
                uimag = ximag [index1];
                
                xreal [index1] = ureal + treal;
                F.real [index1] = xreal [index1];
                ximag [index1] = uimag + timag;
                F.imag [index1] = ximag [index1];
                
                xreal [index2] = ureal - treal;
                F.real [index2] = xreal [index2];
                ximag [index2] = uimag - timag;
                F.imag [index2] = ximag [index2];
                }
            }
        }
        
	return F;
}

fftdata  IFFT (float xreal [], float ximag [], int n)
{
    // ���ٸ���Ҷ��任
    float wreal [N / 2], wimag [N / 2], treal, timag, ureal, uimag, arg;
    int m, k, j, t, index1, index2;
    
    bitrp (xreal, ximag, n);

    // ���� 1 ��ǰ n / 2 �� n �η��� Wj = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
    arg = 2 * PI / n;
    treal = cos (arg);
    timag = sin (arg);
    wreal [0] = 1.0;
    wimag [0] = 0.0;
    for (j = 1; j < n / 2; j ++)
        {
        wreal [j] = wreal [j - 1] * treal - wimag [j - 1] * timag;
        wimag [j] = wreal [j - 1] * timag + wimag [j - 1] * treal;
        }

    for (m = 2; m <= n; m *= 2)
        {
        for (k = 0; k < n; k += m)
            {
            for (j = 0; j < m / 2; j ++)
                {
                index1 = k + j;
                index2 = index1 + m / 2;
                t = n * j / m;    // ��ת���� w ��ʵ���� wreal [] �е��±�Ϊ t
                treal = wreal [t] * xreal [index2] - wimag [t] * ximag [index2];
                timag = wreal [t] * ximag [index2] + wimag [t] * xreal [index2];
                ureal = xreal [index1];
                uimag = ximag [index1];
                xreal [index1] = ureal + treal;
                ximag [index1] = uimag + timag;
                xreal [index2] = ureal - treal;
                ximag [index2] = uimag - timag;
                }
            }
        }

    for (j=0; j < n; j ++)
        {
        xreal [j] /= n;
        ximag [j] /= n;
        
        F.real [j] = xreal [j];
        F.imag [j] = ximag [j];
        }
        
        
        return F;
}



cwt WaveletTransform(float xreal [], float ximag [], float scale, float frequency, int N, float dt)
{
	int i;
	//float creal [], float cimag [], wavefft_real [], float wavefft_imag [];
	
	//wavefft_real[] = FFT_real(Wavelet_Morlet(), N);
	//wavefft_imag[] = FFT_imag(Wavelet_Morlet()), N;
	
	fftdata InputFFT = FFT(xreal, ximag, N);
	
	wavelet InputWave = Wavelet_Morlet(scale, frequency, N);
	
	fftdata InputWaveFFT = FFT(InputWave.real, InputWave.imag, N);
	
	float Normal = sqrt(2*PI*scale/dt);
	
	for (i=0; i<N; i++)
	{
		C.real[i] = Normal*(InputWaveFFT.real[i]*InputFFT.real[i] + InputWaveFFT.imag[i]*InputFFT.imag[i]);
		C.imag[i] = Normal*(InputWaveFFT.real[i]*InputFFT.imag[i] - InputWaveFFT.imag[i]*InputFFT.real[i]);
	}
	
	fftdata IFFTresult = IFFT(C.real, C.imag, N);
	
	for (i=0; i<N; i++)
	{
		C.real[i] = IFFTresult.real[i];
		C.imag[i] = IFFTresult.imag[i];
	}
		
	return C;

}


void CWT_test ()
{
    
	float s = 5;
	float f = 10;
	float dt = 1;
	
	char inputfile [] = "input.txt";    // ���ļ� input.txt �ж���ԭʼ����
    char outputfile [] = "output.txt";    // �����������ļ� output.txt ��
    float xreal [N] = {}, ximag [N] = {};
    int n, i;
    FILE *input, *output;

    if (!(input = fopen (inputfile, "r")))
        {
        printf ("Cannot open file. ");
        exit (1);
        }
    if (!(output = fopen (outputfile, "w")))
        {
        printf ("Cannot open file. ");
        exit (1);
        }
        
    i = 0;
    while ((fscanf (input, "%f%f", xreal + i, ximag + i)) != EOF)
        {
        i ++;
        }
    n = i;    // Ҫ�� n Ϊ 2 ��������
    while (i > 1)
        {
        if (i % 2)
            {
            fprintf (output, "%d is not a power of 2! ", n);
            exit (1);
            }
        i /= 2;
        }

    cwt Result = WaveletTransform(xreal, ximag, s, f, n, dt);
    
    
    
    fprintf (output, "FFT:    i	    real	imag ");
    for (i = 0; i < n; i ++)
        {
        fprintf (output, " \n %4d    %8.4f    %8.4f", i, Result.real [i], Result.imag [i]);
        }
    //fprintf (output, "================================= ");


    if ( fclose (input))
        {
        printf ("File close error. ");
        exit (1);
        }
    if ( fclose (output))
        {
        printf ("File close error. ");
        exit (1);
        }
}

int main ()
{
    CWT_test ();

    return 0;
}

