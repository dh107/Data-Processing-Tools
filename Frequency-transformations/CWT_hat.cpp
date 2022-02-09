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

    int t;
    
    
    for (t=0; t<N; t++)
    
    {
    	float ts = t/scale;
    	//float ts = t/(2*pow(2,(0.125*scale)));
    	W.real[t] = pow(PI,-0.25)*cos(frequency*ts)*exp((-1/2)*pow(ts,2));
    	W.imag[t] = pow(PI,-0.25)*sin(frequency*ts)*exp((-1/2)*pow(ts,2));
	}
	
	return W;
}


wavelet  Wavelet_Hat (float scale, int n, float dt)
{

    int t;
    float s0, s; 
    
    
    for (t=0; t<n; t++)
    
    {
    	float ww = (2*PI*t)/(n*dt);
    	s0 = 2*dt;
		s =  s0*pow(2,(0.125*scale));
    	W.real[t] = (-1*sqrt(8))*(1/(sqrt(3)))*pow(PI,0.25)*pow(s,2.5)*pow(ww,2)*exp((-1/2)*pow((s*ww),2));
    	W.imag[t] = 0.0;
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
    // 位反转置换 Bit-reversal Permutation
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
    // 快速傅立叶变换，将复数 x 变换后仍保存在 x 中，xreal, ximag 分别是 x 的实部和虚部
    float wreal [N / 2], wimag [N / 2], treal, timag, ureal, uimag, arg;
    int m, k, j, t, index1, index2;
    
    bitrp (xreal, ximag, n);

    // 计算 1 的前 n / 2 个 n 次方根的共轭复数 W'j = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
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
                t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
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

        F.real [j] = xreal [j];
        F.imag [j] = ximag [j];
        }
        
	return F;
}

fftdata  IFFT (float xreal [], float ximag [], int n)
{
    // 快速傅立叶逆变换
    float wreal [N / 2], wimag [N / 2], treal, timag, ureal, uimag, arg;
    int m, k, j, t, index1, index2;
    
    bitrp (xreal, ximag, n);

    // 计算 1 的前 n / 2 个 n 次方根 Wj = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
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
                t = n * j / m;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
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



cwt WaveletTransform(float xreal [], float ximag [], float scale, float frequency, int n, float dt)
{
	int i;
    float s0, s; 
    
    s0 = 2*dt;
	s =  s0*pow(2,(0.125*scale));
	
	fftdata InputFFT = FFT(xreal, ximag, n);
	
	wavelet InputWaveFFT = Wavelet_Hat(scale, n, dt);
	//wavelet InputWave = Wavelet_Morlet(scale, frequency, N);
	//fftdata InputWaveFFT = FFT(InputWave.real, InputWave.imag, n);
	
	float Normal = sqrt(2*PI*s/dt);
	
	for (i=0; i<n; i++)
	{
		C.real[i] = 1*(InputWaveFFT.real[i]*InputFFT.real[i] + InputWaveFFT.imag[i]*InputFFT.imag[i]);
		C.imag[i] = 1*(InputWaveFFT.real[i]*InputFFT.imag[i] - InputWaveFFT.imag[i]*InputFFT.real[i]);
	}
	
	fftdata IFFTresult = IFFT(C.real, C.imag, n);
	
	for (i=0; i<n; i++)
	{
		C.real[i] = Normal*IFFTresult.real[i];
		C.imag[i] = Normal*IFFTresult.imag[i];
	}
		
	return C;

}


void CWT_test ()
{
    
	//float s[5] = {0.04, 0.08, 0.12, 0.16, 0.20} ;
	int M=10;
	float s[M] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0} ;
	float f = 10;
	float dt = 1.0;
	
	char inputfile [] = "input4.txt";    // 从文件 input.txt 中读入原始数据
    char outputfile [] = "output.txt";    // 将结果输出到文件 output.txt 中
    float xreal [N] = {}, ximag [N] = {};
    int n, i, j;
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
    n = i;    // 要求 n 为 2 的整数幂
    while (i > 1)
        {
        if (i % 2)
            {
            fprintf (output, "%d is not a power of 2! ", n);
            exit (1);
            }
        i /= 2;
        }

	for (j=0; j<M; j++)
	{
    		cwt Result = WaveletTransform(xreal, ximag, s[j], f, n, dt);
    
    
    
    //fprintf (output, "FFT:    i	    real	imag ");
    for (i = 0; i < n; i ++)
        {
        fprintf (output, " \n %4d   %8.2f   %8.4f    %8.4f", i, s[j], Result.real [i], Result.imag [i]); 
        }
    //fprintf (output, "================================= ");
    
   

	}

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

