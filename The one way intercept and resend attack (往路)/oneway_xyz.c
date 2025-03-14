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
    FILE *file;

    file = fopen("FinalData_xyz.csv", "w");

    if (file == NULL) {
        printf("ファイルを開けませんでした。\n");
        return 1;
    }
    fprintf(file, "a,b,PAB,PE,IAE,IAB\n");
    complex onep = {1, 1}; // (1+i)
    complex onen = {1, -1}; // (1-i)
    //printf("開始！\n");
    // i, j, k,lは角度。a,bはラジアン。
    //#pragma omp parallel for
    for (int i = 0; i < 360; i++){
        double a = i * M_PI / 180; // ラジアン表記
        complex eap = {cos(a), sin(a)}; // e^(iα)
        complex ean = {cos(a), -sin(a)}; // e^(-iα)
        complex onepeap = pro(onep,eap); // (1+i)e^(iα)
        complex oneneap = pro(onen, eap); // (1-i)e^(iα)
        complex onepean = pro(onep,ean); // (1+i)e^(-iα)
        complex onenean = pro(onen, ean); // (1-i)e^(-iα)
        #pragma omp parallel for
        for (int k = 0; k<360; k++){
            double b = k * M_PI / 180; // ラジアン表記
            complex tmp;
            
            // 要素全部
            /////////////// r1 (往路) ///////////////
            complex r1x0g0, r1x0g1, r1x1g0, r1x1g1;
            r1x0g0.real = (2+cos(b)+2*cos(a)*sin(b)+sin(a)*sin(b))/8;
            r1x0g0.imag = (-cos(b)+sin(a)*sin(b))/8;
            r1x0g1.real = (2-cos(b)-2*cos(a)*sin(b)-sin(a)*sin(b))/8;
            r1x0g1.imag = (cos(b)-sin(a)*sin(b))/8;
            r1x1g0.real = (cos(b)+sin(a)*sin(b))/8;
            r1x1g0.imag = (cos(b)-sin(a)*sin(b))/8;
            r1x1g1.real = (-cos(b)-sin(a)*sin(b))/8;
            r1x1g1.imag = (-cos(b)+sin(a)*sin(b))/8;

            /////////////// r2 (往路) ///////////////
            complex r2x0g0, r2x0g1, r2x1g0, r2x1g1;
            r2x0g0.real = (cos(b)-sin(a)*sin(b))/8;
            r2x0g0.imag = (cos(b)+sin(a)*sin(b))/8;
            r2x0g1.real = (-cos(b)+sin(a)*sin(b))/8;
            r2x0g1.imag = (-cos(b)-sin(a)*sin(b))/8;
            r2x1g0.real = (2+cos(b)-2*cos(a)*sin(b)-sin(a)*sin(b))/8;
            r2x1g0.imag = (-cos(b)-sin(a)*sin(b))/8;
            r2x1g1.real = (2-cos(b)+2*cos(a)*sin(b)+sin(a)*sin(b))/8;
            r2x1g1.imag = (cos(b)+sin(a)*sin(b))/8;

            /////////////// r3 (往路) ///////////////
            complex r3x0g0, r3x0g1, r3x1g0, r3x1g1;
            r3x0g0.real = (2-cos(b)+2*cos(a)*sin(b)-sin(a)*sin(b))/8;
            r3x0g0.imag = (cos(b)-sin(a)*sin(b))/8;
            r3x0g1.real = (2+cos(b)-2*cos(a)*sin(b)+sin(a)*sin(b))/8;
            r3x0g1.imag = (-cos(b)+sin(a)*sin(b))/8;
            r3x1g0.real = (-cos(b)-sin(a)*sin(b))/8;
            r3x1g0.imag = (-cos(b)+sin(a)*sin(b))/8;
            r3x1g1.real = (cos(b)+sin(a)*sin(b))/8;
            r3x1g1.imag = (cos(b)-sin(a)*sin(b))/8;

            /////////////// r4 (往路) ///////////////
            complex r4x0g0, r4x0g1, r4x1g0, r4x1g1;
            r4x0g0.real = (-cos(b)+sin(a)*sin(b))/8;
            r4x0g0.imag = (-cos(b)-sin(a)*sin(b))/8;
            r4x0g1.real = (cos(b)-sin(a)*sin(b))/8;
            r4x0g1.imag = (cos(b)+sin(a)*sin(b))/8;
            r4x1g0.real = (2-cos(b)-2*cos(a)*sin(b)+sin(a)*sin(b))/8;
            r4x1g0.imag = (cos(b)+sin(a)*sin(b))/8;
            r4x1g1.real = (2+cos(b)+2*cos(a)*sin(b)-sin(a)*sin(b))/8;
            r4x1g1.imag = (-cos(b)-sin(a)*sin(b))/8;

            /////////////// r1 (往路) ///////////////
            complex r1z0g0, r1z0g1, r1z1g0, r1z1g1;
            r1z0g0 = mul(sin(b)/8, onepean);
            r1z0g0.real += (1+cos(b))/4; 
            r1z0g1 = mul(-sin(b)/8,onepean);
            r1z0g1.real += (1-cos(b))/4;
            r1z1g0 = mul(sin(b)/8, oneneap);
            r1z1g1 = mul(-sin(b)/8, oneneap);
            
            /////////////// r2 (往路) ///////////////
            complex r2z0g0, r2z0g1, r2z1g0, r2z1g1;
            r2z0g0 = mul(-sin(b)/8, onepean);
            r2z0g0.real += (1+cos(b))/4; 
            r2z0g1 = mul(sin(b)/8,onepean);
            r2z0g1.real += (1-cos(b))/4;
            r2z1g0 = mul(-sin(b)/8, oneneap);
            r2z1g1 = mul(sin(b)/8, oneneap);

            /////////////// r3 (往路) ///////////////
            complex r3z0g0, r3z0g1, r3z1g0, r3z1g1;
            r3z0g0 = mul(sin(b)/8, onenean);
            r3z0g1 = mul(-sin(b)/8, onenean);
            r3z1g0 = mul(sin(b)/8, onepeap);
            r3z1g0.real += (1-cos(b))/4; 
            r3z1g1 = mul(-sin(b)/8,onepeap);
            r3z1g1.real += (1+cos(b))/4;

            /////////////// r4 (往路) ///////////////
            complex r4z0g0, r4z0g1, r4z1g0, r4z1g1;
            r4z0g0 = mul(-sin(b)/8, onenean);
            r4z0g1 = mul(sin(b)/8, onenean);
            r4z1g0 = mul(-sin(b)/8, onepeap);
            r4z1g0.real += (1-cos(b))/4; 
            r4z1g1 = mul(sin(b)/8,onepeap);
            r4z1g1.real += (1+cos(b))/4;

            /////////////// r1 (往路) ///////////////
            complex r1y0g0, r1y0g1, r1y1g0, r1y1g1;
            r1y0g0.real = (2+cos(b)+cos(a)*sin(b)+2*sin(a)*sin(b))/8;
            r1y0g0.imag = (cos(b)-cos(a)*sin(b))/8;
            r1y0g1.real = (2-cos(b)-cos(a)*sin(b)-2*sin(a)*sin(b))/8;
            r1y0g1.imag = (-cos(b)+cos(a)*sin(b))/8;
            r1y1g0.real = (cos(b)+cos(a)*sin(b))/8;
            r1y1g0.imag = (-cos(b)+cos(a)*sin(b))/8;
            r1y1g1.real = (-cos(b)-cos(a)*sin(b))/8;
            r1y1g1.imag = (cos(b)-cos(a)*sin(b))/8;

            /////////////// r2 (往路) ///////////////
            complex r2y0g0, r2y0g1, r2y1g0, r2y1g1;
            r2y0g0.real = (cos(b)-cos(a)*sin(b))/8;
            r2y0g0.imag = (-cos(b)-cos(a)*sin(b))/8;
            r2y0g1.real = (-cos(b)+cos(a)*sin(b))/8;
            r2y0g1.imag = (cos(b)+cos(a)*sin(b))/8;
            r2y1g0.real = (2+cos(b)-cos(a)*sin(b)-2*sin(a)*sin(b))/8;
            r2y1g0.imag = (cos(b)+cos(a)*sin(b))/8;
            r2y1g1.real = (2-cos(b)+cos(a)*sin(b)+2*sin(a)*sin(b))/8;
            r2y1g1.imag = (-cos(b)-cos(a)*sin(b))/8;


            /////////////// r3 (往路) ///////////////
            complex r3y0g0, r3y0g1, r3y1g0, r3y1g1;
            r3y0g0.real = (-cos(b)+cos(a)*sin(b))/8;
            r3y0g0.imag = (cos(b)+cos(a)*sin(b))/8;
            r3y0g1.real = (cos(b)-cos(a)*sin(b))/8;
            r3y0g1.imag = (-cos(b)-cos(a)*sin(b))/8;
            r3y1g0.real = (2-cos(b)+cos(a)*sin(b)-2*sin(a)*sin(b))/8;
            r3y1g0.imag = (-cos(b)-cos(a)*sin(b))/8;
            r3y1g1.real = (2+cos(b)-cos(a)*sin(b)+2*sin(a)*sin(b))/8;
            r3y1g1.imag = (cos(b)+cos(a)*sin(b))/8;

            /////////////// r4 (往路) ///////////////
            complex r4y0g0, r4y0g1, r4y1g0, r4y1g1;
            r4y0g0.real = (2-cos(b)-cos(a)*sin(b)+2*sin(a)*sin(b))/8;
            r4y0g0.imag = (-cos(b)+cos(a)*sin(b))/8;
            r4y0g1.real = (2+cos(b)+cos(a)*sin(b)-2*sin(a)*sin(b))/8;
            r4y0g1.imag = (cos(b)-cos(a)*sin(b))/8;
            r4y1g0.real = (-cos(b)-cos(a)*sin(b))/8;
            r4y1g0.imag = (cos(b)-cos(a)*sin(b))/8;
            r4y1g1.real = (cos(b)+cos(a)*sin(b))/8;
            r4y1g1.imag = (-cos(b)+cos(a)*sin(b))/8;


            // 関数定義部分
            double f1,f2,f3,f4,f5,f6,f7,f8;
            f1 = Pccn(r1x0g0) + Pccn(r1x0g1);
            f2 = Pccn(r1x1g0) + Pccn(r1x1g1);
            f3 = Pccn(r2x0g0) + Pccn(r2x0g1);
            f4 = Pccn(r2x1g0) + Pccn(r2x1g1);
            f5 = Pccn(r3x0g0) + Pccn(r3x0g1);
            f6 = Pccn(r3x1g0) + Pccn(r3x1g1);
            f7 = Pccn(r4x0g0) + Pccn(r4x0g1);
            f8 = Pccn(r4x1g0) + Pccn(r4x1g1);

            double l1,l2,l3,l4,l5,l6,l7,l8;
            l1 = Pccn(r1y0g0) + Pccn(r1y0g1);
            l2 = Pccn(r1y1g0) + Pccn(r1y1g1);
            l3 = Pccn(r2y0g0) + Pccn(r2y0g1);
            l4 = Pccn(r2y1g0) + Pccn(r2y1g1);
            l5 = Pccn(r3y0g0) + Pccn(r3y0g1);
            l6 = Pccn(r3y1g0) + Pccn(r3y1g1);
            l7 = Pccn(r4y0g0) + Pccn(r4y0g1);
            l8 = Pccn(r4y1g0) + Pccn(r4y1g1);

            double u1,u2,u3,u4,u5,u6,u7,u8;
            u1 = Pccn(r1z0g0) + Pccn(r1z0g1);
            u2 = Pccn(r1z1g0) + Pccn(r1z1g1);
            u3 = Pccn(r2z0g0) + Pccn(r2z0g1);
            u4 = Pccn(r2z1g0) + Pccn(r2z1g1);
            u5 = Pccn(r3z0g0) + Pccn(r3z0g1);
            u6 = Pccn(r3z1g0) + Pccn(r3z1g1);
            u7 = Pccn(r4z0g0) + Pccn(r4z0g1);
            u8 = Pccn(r4z1g0) + Pccn(r4z1g1);

            // 関数定義部分
            double g1,g2,g3,g4,g7,g8,m1,m2,m3,m4,m5,m6,v1,v2,v5,v6,v7,v8;
            g1 = Pccn(r1x0g0) + Pccn(r1x1g0);
            g2 = Pccn(r1x0g1) + Pccn(r1x1g1);
            g3 = Pccn(r2x0g0) + Pccn(r2x1g0);
            g4 = Pccn(r2x0g1) + Pccn(r2x1g1);
            g7 = Pccn(r4x0g0) + Pccn(r4x1g0);
            g8 = Pccn(r4x0g1) + Pccn(r4x1g1);
            m1 = Pccn(r1y0g0) + Pccn(r1y1g0);
            m2 = Pccn(r1y0g1) + Pccn(r1y1g1);
            m3 = Pccn(r2y0g0) + Pccn(r2y1g0);
            m4 = Pccn(r2y0g1) + Pccn(r2y1g1);
            m5 = Pccn(r3y0g0) + Pccn(r3y1g0);
            m6 = Pccn(r3y0g1) + Pccn(r3y1g1);
            v1 = Pccn(r1z0g0) + Pccn(r1z1g0);
            v2 = Pccn(r1z0g1) + Pccn(r1z1g1);
            v5 = Pccn(r3z0g0) + Pccn(r3z1g0);
            v6 = Pccn(r3z0g1) + Pccn(r3z1g1);
            v7 = Pccn(r4z0g0) + Pccn(r4z1g0);
            v8 = Pccn(r4z0g1) + Pccn(r4z1g1);

            // xy 
            double R = f5+f6+f7+f8+l5+l6+l7+l8;
            double R3 = (f5+f6+l5+l6)/R;
            double R4 = (f7+f8+l7+l8)/R;
            double PABxy = R3*(f5/(f5+f6) + l6/(l5+l6))/2 + R4*(f8/(f7+f8) + l7/(l7+l8))/2;

            R = g1+g2+g3+g4+m1+m2+m3+m4;
            double R1 = (g1+g2+m1+m2)/R;
            double R2 = (g3+g4+m3+m4)/R;
            double PExy = R1*(g1/(g1+g2) + m1/(m1+m2))/2 + R2*(g4/(g3+g4) + m4/(m3+m4))/2;

            R = g1+g2+g3+g4+m1+m2+m3+m4;
            double A0E0 = (g1+m1)/R;
            double A0E1 = (g2+m2)/R;
            double A1E0 = (g3+m3)/R;
            double A1E1 = (g4+m4)/R;
            double IAExy = A0E0*log2(A0E0/((A0E0+A0E1)*(A0E0+A1E0)))+A0E1*log2(A0E1/((A0E0+A0E1)*(A0E1+A1E1)))+A1E0*log2(A1E0/((A1E0+A1E1)*(A0E0+A1E0)))+A1E1*log2(A1E1/((A1E0+A1E1)*(A0E1+A1E1)));

            R = f1+f2+f3+f4+l1+l2+l3+l4;
            double A0B0 = (f1+l1)/R;
            double A0B1 = (f2+l2)/R;
            double A1B0 = (f3+l3)/R;
            double A1B1 = (f4+l4)/R;
            double IABxy = A0B0*log2(A0B0/((A0B0+A0B1)*(A0B0+A1B0)))+A0B1*log2(A0B1/((A0B0+A0B1)*(A0B1+A1B1)))+A1B0*log2(A1B0/((A1B0+A1B1)*(A0B0+A1B0)))+A1B1*log2(A1B1/((A1B0+A1B1)*(A0B1+A1B1)));

            // xz 
            R = g1+g2+g7+g8+v1+v2+v7+v8;
            R1 = (g1+g2+v1+v2)/R;
            R4 = (g7+g8+v7+v8)/R;
            double PExz = R1*(g1/(g1+g2) + v1/(v1+v2))/2 + R4*(g8/(g7+g8) + v8/(v7+v8))/2;

            R = f3+f4+f5+f6+u3+u4+u5+u6;
            R2 = (f3+f4+u3+u4)/R;
            R3 = (f5+f6+u5+u6)/R;
            double PABxz = R2*(f4/(f3+f4) + u3/(u3+u4))/2 + R3*(f5/(f5+f6) + u6/(u5+u6))/2;

            R = g1+g2+g7+g8+v1+v2+v7+v8;
            A0E0 = (g1+v1)/R;
            A0E1 = (g2+v2)/R;
            A1E0 = (g7+v7)/R;
            A1E1 = (g8+v8)/R;
            double IAExz = A0E0*log2(A0E0/((A0E0+A0E1)*(A0E0+A1E0)))+A0E1*log2(A0E1/((A0E0+A0E1)*(A0E1+A1E1)))+A1E0*log2(A1E0/((A1E0+A1E1)*(A0E0+A1E0)))+A1E1*log2(A1E1/((A1E0+A1E1)*(A0E1+A1E1)));

            R = f1+f2+f7+f8+u1+u2+u7+u8;
            A0B0 = (f1+u1)/R;
            A0B1 = (f2+u2)/R;
            A1B0 = (f7+u7)/R;
            A1B1 = (f8+u8)/R;
            double IABxz = A0B0*log2(A0B0/((A0B0+A0B1)*(A0B0+A1B0)))+A0B1*log2(A0B1/((A0B0+A0B1)*(A0B1+A1B1)))+A1B0*log2(A1B0/((A1B0+A1B1)*(A0B0+A1B0)))+A1B1*log2(A1B1/((A1B0+A1B1)*(A0B1+A1B1)));
            
            // yz 
            R = l3+l4+l7+l8+u3+u4+u7+u8;
            R2 = (l3+l4+u3+u4)/R;
            R4 = (l7+l8+u7+u8)/R;
            double PAByz = R2*(l4/(l3+l4) + u3/(u3+u4))/2 + R4*(l7/(l7+l8) + u8/(u7+u8))/2;

            R = m1+m2+m5+m6+v1+v2+v5+v6;
            R1 = (m1+m2+v1+v2)/R;
            R3 = (m5+m6+v5+v6)/R;
            double PEyz = R1*(v1/(v1+v2) +m1/(m1+m2))/2 + R3*(v6/(v5+v6) + m6/(m5+m6))/2;

            R = m1+v1+m2+v2+m5+v5+m6+v6;
            A0E0 = (m1+v1)/R;
            A0E1 = (m2+v2)/R;
            A1E0 = (m5+v5)/R;
            A1E1 = (m6+v6)/R;
            double IAEyz = A0E0*log2(A0E0/((A0E0+A0E1)*(A0E0+A1E0)))+A0E1*log2(A0E1/((A0E0+A0E1)*(A0E1+A1E1)))+A1E0*log2(A1E0/((A1E0+A1E1)*(A0E0+A1E0)))+A1E1*log2(A1E1/((A1E0+A1E1)*(A0E1+A1E1)));

            R = l1+l2+l5+l6+u1+u2+u5+u6;
            A0B0 = (l1+u1)/R;
            A0B1 = (l2+u2)/R;
            A1B0 = (l5+u5)/R;
            A1B1 = (l6+u6)/R;
            double IAByz = A0B0*log2(A0B0/((A0B0+A0B1)*(A0B0+A1B0)))+A0B1*log2(A0B1/((A0B0+A0B1)*(A0B1+A1B1)))+A1B0*log2(A1B0/((A1B0+A1B1)*(A0B0+A1B0)))+A1B1*log2(A1B1/((A1B0+A1B1)*(A0B1+A1B1)));

            double PAB = (PABxz + PABxy + PAByz)/3;
            double PE  = (PExz  + PExy  + PEyz )/3;
            double IAE = (IAExz + IAExy + IAEyz)/3;
            double IAB = (IABxz + IABxy + IAByz)/3;

            // xz 
            double r11 = g1+g2+v1+v2;
            double r12 = f3+f4+u3+u4;
            double r13 = f5+f6+u5+u6;
            double r14 = g7+g8+v7+v8;
            
            // yz 
            double r21 = m1+m2+v1+v2;
            double r22 = l3+l4+u3+u4;
            double r23 = m5+m6+v5+v6;
            double r24 = l7+l8+u7+u8;
            // xy 
            double r31 = m1+m2+g1+g2;
            double r32 = m3+m4+g3+g4;
            double r33 = l5+l6+f5+f6;
            double r34 = l7+l8+f7+f8;

            if(fabs(r11 - r12) <= DBL_EPSILON * fmax(1, fmax(fabs(r11), fabs(r12))) && fabs(r11 - r13) <= DBL_EPSILON * fmax(1, fmax(fabs(r11), fabs(r13))) && fabs(r11 - r14) <= DBL_EPSILON * fmax(1, fmax(fabs(r11), fabs(r14))) && fabs(r12 - r13) <= DBL_EPSILON * fmax(1, fmax(fabs(r12), fabs(r13))) && fabs(r12 - r14) <= DBL_EPSILON * fmax(1, fmax(fabs(r12), fabs(r14))) && fabs(r13 - r14) <= DBL_EPSILON * fmax(1, fmax(fabs(r13), fabs(r14)))){
                if(fabs(r21 - r22) <= DBL_EPSILON * fmax(1, fmax(fabs(r21), fabs(r22))) && fabs(r21 - r23) <= DBL_EPSILON * fmax(1, fmax(fabs(r21), fabs(r23))) && fabs(r21 - r24) <= DBL_EPSILON * fmax(1, fmax(fabs(r21), fabs(r24))) && fabs(r22 - r23) <= DBL_EPSILON * fmax(1, fmax(fabs(r22), fabs(r23))) && fabs(r22 - r24) <= DBL_EPSILON * fmax(1, fmax(fabs(r22), fabs(r24))) && fabs(r23 - r24) <= DBL_EPSILON * fmax(1, fmax(fabs(r23), fabs(r24)))){
                    if(fabs(r31 - r32) <= DBL_EPSILON * fmax(1, fmax(fabs(r31), fabs(r32))) && fabs(r31 - r33) <= DBL_EPSILON * fmax(1, fmax(fabs(r31), fabs(r33))) && fabs(r31 - r34) <= DBL_EPSILON * fmax(1, fmax(fabs(r31), fabs(r34))) && fabs(r32 - r33) <= DBL_EPSILON * fmax(1, fmax(fabs(r32), fabs(r33))) && fabs(r32 - r34) <= DBL_EPSILON * fmax(1, fmax(fabs(r32), fabs(r34))) && fabs(r33 - r34) <= DBL_EPSILON * fmax(1, fmax(fabs(r33), fabs(r34)))){
                        fprintf(file, "%d,%d,%.32lf,%.32lf,%.32lf, %.32lf\n", i,k,PAB,PE,IAE,IAB);
                    }
                }
            }
        }
    }
    fclose(file);
    return 0;
}
