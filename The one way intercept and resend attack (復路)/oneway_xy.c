#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <float.h>

// 複素数型の定義
typedef struct{
	double real; // 実部
	double imag; // 虚部
} complex;

// 複素数の和をとる関数
complex sum(complex a, complex b){
	complex ans;

	ans.real = a.real + b.real;
	ans.imag = a.imag + b.imag;

	return(ans);
}

// 複素数の内積を行う関数
complex pro(complex a, complex b){
	complex ans;

	ans.real = a.real*b.real - a.imag*b.imag;
	ans.imag = a.real*b.imag + a.imag*b.real;

	return(ans);
}

// 複素数に定数倍を行う関数
complex mul(double a, complex b){
	complex ans;

	ans.real = a*b.real;
	ans.imag = a*b.imag;

	return(ans);
}

// 共役複素数の積を行う関数
double Pccn(complex a){
    double ans;

    ans = a.real*a.real + a.imag*a.imag;
    
    return (ans);
}

int main (void){
    complex onep = {1, 1}; // (1+i)
    complex onen = {1, -1}; // (1-i)
    for (int i = 0; i < 360; i++){
        double c = i * M_PI / 180; // ラジアン表記
        complex ecp = {cos(c), sin(c)}; // e^(iγ)
        complex ecn = {cos(c), -sin(c)}; // e^(-iγ)
        complex onepecp = pro(onep,ecp); // (1+i)e^(iγ)
        complex onenecp = pro(onen, ecp); // (1-i)e^(iγ)
        complex onepecn = pro(onep,ecn); // (1+i)e^(-iγ)
        complex onenecn = pro(onen, ecn); // (1-i)e^(-iγ)
        //printf("Now %d\n",i);
        #pragma omp parallel for
        for (int k = 0; k<360; k++){
            double d = k * M_PI / 180; // ラジアン表記
            complex tmp;

            /////////////// r1 (復路) ///////////////
            complex r1m0x0, r1m0x1, r1m1x0, r1m1x1;
            r1m0x0.real = (2+cos(d)+2*cos(c)*sin(d)+sin(c)*sin(d))/8;
            r1m0x0.imag = (cos(d)-sin(c)*sin(d))/8;
            r1m0x1.real = (cos(d)+sin(c)*sin(d))/8;
            r1m0x1.imag = (-cos(d)+sin(c)*sin(d))/8;
            r1m1x0.real = (2-cos(d)-2*cos(c)*sin(d)-sin(c)*sin(d))/8;
            r1m1x0.imag = (-cos(d)+sin(c)*sin(d))/8;
            r1m1x1.real = (-cos(d)-sin(c)*sin(d))/8;
            r1m1x1.imag = (cos(d)-sin(c)*sin(d))/8;

            /////////////// r2 (復路) ///////////////
            complex r2m0x0, r2m0x1, r2m1x0, r2m1x1;
            r2m0x0.real = (cos(d)-sin(c)*sin(d))/8;
            r2m0x0.imag = (-cos(d)-sin(c)*sin(d))/8;
            r2m0x1.real = (2+cos(d)-2*cos(c)*sin(d)-sin(c)*sin(d))/8;
            r2m0x1.imag = (cos(d)+sin(c)*sin(d))/8;
            r2m1x0.real = (-cos(d)+sin(c)*sin(d))/8;
            r2m1x0.imag = (cos(d)+sin(c)*sin(d))/8;
            r2m1x1.real = (2-cos(d)+2*cos(c)*sin(d)+sin(c)*sin(d))/8;
            r2m1x1.imag = (-cos(d)-sin(c)*sin(d))/8;

            /////////////// r3 (復路) ///////////////
            complex r3m0x0, r3m0x1, r3m1x0, r3m1x1;
            r3m0x0.real = (2-cos(d)+2*cos(c)*sin(d)-sin(c)*sin(d))/8;
            r3m0x0.imag = (-cos(d)+sin(c)*sin(d))/8;
            r3m0x1.real = (-cos(d)-sin(c)*sin(d))/8;
            r3m0x1.imag = (cos(d)-sin(c)*sin(d))/8;
            r3m1x0.real = (2+cos(d)-2*cos(c)*sin(d)+sin(c)*sin(d))/8;
            r3m1x0.imag = (+cos(d)-sin(c)*sin(d))/8;
            r3m1x1.real = (cos(d)+sin(c)*sin(d))/8;
            r3m1x1.imag = (-cos(d)+sin(c)*sin(d))/8;

            /////////////// r4 (復路) ///////////////
            complex r4m0x0, r4m0x1, r4m1x0, r4m1x1;
            r4m0x0.real = (-cos(d)+sin(c)*sin(d))/8;
            r4m0x0.imag = (cos(d)+sin(c)*sin(d))/8;
            r4m0x1.real = (2-cos(d)-2*cos(c)*sin(d)+sin(c)*sin(d))/8;
            r4m0x1.imag = (-cos(d)-sin(c)*sin(d))/8;
            r4m1x0.real = (cos(d)-sin(c)*sin(d))/8;
            r4m1x0.imag = (-cos(d)-sin(c)*sin(d))/8;
            r4m1x1.real = (2+cos(d)+2*cos(c)*sin(d)-sin(c)*sin(d))/8;
            r4m1x1.imag = (cos(d)+sin(c)*sin(d))/8;

            /////////////// r1 (復路) ///////////////
            complex r1m0y0, r1m0y1, r1m1y0, r1m1y1;
            r1m0y0.real = (2+cos(d)+cos(c)*sin(d)+2*sin(c)*sin(d))/8;
            r1m0y0.imag = (-cos(d)+cos(c)*sin(d))/8;
            r1m0y1.real = (cos(d)+cos(c)*sin(d))/8;
            r1m0y1.imag = (cos(d)-cos(c)*sin(d))/8;
            r1m1y0.real = (2-cos(d)-cos(c)*sin(d)-2*sin(c)*sin(d))/8;
            r1m1y0.imag = (cos(d)-cos(c)*sin(d))/8;
            r1m1y1.real = (-cos(d)-cos(c)*sin(d))/8;
            r1m1y1.imag = (-cos(d)+cos(c)*sin(d))/8;

            /////////////// r2 (復路) ///////////////
            complex r2m0y0, r2m0y1, r2m1y0, r2m1y1;
            r2m0y0.real = (cos(d)-cos(c)*sin(d))/8;
            r2m0y0.imag = (cos(d)+cos(c)*sin(d))/8;
            r2m0y1.real = (2+cos(d)-cos(c)*sin(d)-2*sin(c)*sin(d))/8;
            r2m0y1.imag = (-cos(d)-cos(c)*sin(d))/8;
            r2m1y0.real = (-cos(d)+cos(c)*sin(d))/8;
            r2m1y0.imag = (-cos(d)-cos(c)*sin(d))/8;
            r2m1y1.real = (2-cos(d)+cos(c)*sin(d)+2*sin(c)*sin(d))/8;
            r2m1y1.imag = (cos(d)+cos(c)*sin(d))/8;

            /////////////// r3 (復路) ///////////////
            complex r3m0y0, r3m0y1, r3m1y0, r3m1y1;
            r3m0y0.real = (-cos(d)+cos(c)*sin(d))/8;
            r3m0y0.imag = (-cos(d)-cos(c)*sin(d))/8;
            r3m0y1.real = (2-cos(d)+cos(c)*sin(d)-2*sin(c)*sin(d))/8;
            r3m0y1.imag = (cos(d)+cos(c)*sin(d))/8;
            r3m1y0.real = (cos(d)-cos(c)*sin(d))/8;
            r3m1y0.imag = (cos(d)+cos(c)*sin(d))/8;
            r3m1y1.real = (2+cos(d)-cos(c)*sin(d)+2*sin(c)*sin(d))/8;
            r3m1y1.imag = (-cos(d)-cos(c)*sin(d))/8;

            /////////////// r4 (復路) ///////////////
            complex r4m0y0, r4m0y1, r4m1y0, r4m1y1;
            r4m0y0.real = (2-cos(d)-cos(c)*sin(d)+2*sin(c)*sin(d))/8;
            r4m0y0.imag = (cos(d)-cos(c)*sin(d))/8;
            r4m0y1.real = (-cos(d)-cos(c)*sin(d))/8;
            r4m0y1.imag = (-cos(d)+cos(c)*sin(d))/8;
            r4m1y0.real = (2+cos(d)+cos(c)*sin(d)-2*sin(c)*sin(d))/8;
            r4m1y0.imag = (-cos(d)+cos(c)*sin(d))/8;
            r4m1y1.real = (cos(d)+cos(c)*sin(d))/8;
            r4m1y1.imag = (cos(d)-cos(c)*sin(d))/8;

            // 関数定義部分
            double f5,f6,f7,f8,l5,l6,l7,l8;
            f5 = Pccn(r3m0x0) + Pccn(r3m1x0);
            f6 = Pccn(r3m0x1) + Pccn(r3m1x1);
            f7 = Pccn(r4m0x0) + Pccn(r4m1x0);
            f8 = Pccn(r4m0x1) + Pccn(r4m1x1);
            l5 = Pccn(r3m0y0) + Pccn(r3m1y0);
            l6 = Pccn(r3m0y1) + Pccn(r3m1y1);
            l7 = Pccn(r4m0y0) + Pccn(r4m1y0);
            l8 = Pccn(r4m0y1) + Pccn(r4m1y1);

            double f1,f2,f3,f4,l1,l2,l3,l4;
            f1 = Pccn(r1m0x0) + Pccn(r1m1x0);
            f2 = Pccn(r1m0x1) + Pccn(r1m1x1);
            f3 = Pccn(r2m0x0) + Pccn(r2m1x0);
            f4 = Pccn(r2m0x1) + Pccn(r2m1x1);
            l1 = Pccn(r1m0y0) + Pccn(r1m1y0);
            l2 = Pccn(r1m0y1) + Pccn(r1m1y1);
            l3 = Pccn(r2m0y0) + Pccn(r2m1y0);
            l4 = Pccn(r2m0y1) + Pccn(r2m1y1);

            double R = f5+f6+f7+f8+l5+l6+l7+l8;
            double R3 = (f5+f6+l5+l6)/R;
            double R4 = (f7+f8+l7+l8)/R;

            double PAB = R3*(f5/(f5+f6) + l6/(l5+l6))/2 + R4*(f8/(f7+f8) + l7/(l7+l8))/2;

            // 関数定義部分
            double g1,g2,g3,g4,m1,m2,m3,m4;
            g1 = Pccn(r1m0x0) + Pccn(r1m0x1);
            g2 = Pccn(r1m1x0) + Pccn(r1m1x1);
            g3 = Pccn(r2m0x0) + Pccn(r2m0x1);
            g4 = Pccn(r2m1x0) + Pccn(r2m1x1);
            m1 = Pccn(r1m0y0) + Pccn(r1m0y1);
            m2 = Pccn(r1m1y0) + Pccn(r1m1y1);
            m3 = Pccn(r2m0y0) + Pccn(r2m0y1);
            m4 = Pccn(r2m1y0) + Pccn(r2m1y1);

            R = g1+g2+g3+g4+m1+m2+m3+m4;
            double R1 = (g1+g2+m1+m2)/R;
            double R2 = (g3+g4+m3+m4)/R;

            double PE = R1*(g1/(g1+g2) +m1/(m1+m2))/2 + R2*(g4/(g3+g4) + m4/(m3+m4))/2;

            double R12 = g1+g2+g3+g4+m1+m2+m3+m4;
            double A0E0 = (g1+m1)/R12;
            double A0E1 = (g2+m2)/R12;
            double A1E0 = (g3+m3)/R12;
            double A1E1 = (g4+m4)/R12;
            double IAE = A0E0*log2(A0E0/((A0E0+A0E1)*(A0E0+A1E0)))+A0E1*log2(A0E1/((A0E0+A0E1)*(A0E1+A1E1)))+A1E0*log2(A1E0/((A1E0+A1E1)*(A0E0+A1E0)))+A1E1*log2(A1E1/((A1E0+A1E1)*(A0E1+A1E1)));

            R = f1+f2+f3+f4+l1+l2+l3+l4;
            double A0B0 = (f1+l1)/R;
            double A0B1 = (f2+l2)/R;
            double A1B0 = (f3+l3)/R;
            double A1B1 = (f4+l4)/R;
            double IAB = A0B0*log2(A0B0/((A0B0+A0B1)*(A0B0+A1B0)))+A0B1*log2(A0B1/((A0B0+A0B1)*(A0B1+A1B1)))+A1B0*log2(A1B0/((A1B0+A1B1)*(A0B0+A1B0)))+A1B1*log2(A1B1/((A1B0+A1B1)*(A0B1+A1B1)));

            double r1 = (m1+m2+g1+g2);
            double r2 = (m3+m4+g3+g4);
            double r3 = (l5+l6+f5+f6);
            double r4 = (l7+l8+f7+f8);
            

            //if (PAB >= 0.899999){
            //    printf("a = %d, b = %d, PAB=%lf, PE=%lf, IAE=%lf, IAB=%lf, %lf,%lf,%lf,%lf\n", i, k, PAB,PE, IAE, IAB, r1,r2,r3,r4);
            //}
                    
            // 探索判定部
            //if (PE >= 0.92426){
            //    printf("a = %d, b = %d, PAB=%lf, PE=%lf, IAE=%lf, IAB=%lf, %lf,%lf,%lf,%lf\n", i, k, PAB,PE, IAE, IAB, r1,r2,r3,r4);
            //}
            // 探索判定部
            //if (IAE > 0.612979){
            //    printf("a = %d, b = %d, PAB=%lf, PE=%lf, IAE=%lf, IAB=%lf, %lf,%lf,%lf,%lf\n", i, k, PAB,PE, IAE, IAB, r1,r2,r3,r4);
            //}
        }
    }
    return 0;
}
