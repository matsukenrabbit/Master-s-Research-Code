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

    file = fopen("FD_xyzalt(0~359).csv", "w");

    if (file == NULL) {
        printf("ファイルを開けませんでした。\n");
        return 1;
    }
    fprintf(file, "a,b,PAB,PE,IAE,IAB\n");
    complex onep = {1, 1}; // (1+i)
    complex onen = {1, -1}; // (1-i)
    //printf("開始！\n");
    // i, j, k,lは角度。a,b,c,dはラジアン。
    //#pragma omp parallel for
    for (int i = 0; i < 360; i++){
        double c = i * M_PI / 180; // ラジアン表記
        complex ecp = {cos(c), sin(c)}; // e^(iγ)
        complex ecn = {cos(c), -sin(c)}; // e^(-iγ)
        complex onepecp = pro(onep,ecp); // (1+i)e^(iγ)
        complex onenecp = pro(onen, ecp); // (1-i)e^(iγ)
        complex onepecn = pro(onep,ecn); // (1+i)e^(-iγ)
        complex onenecn = pro(onen, ecn); // (1-i)e^(-iγ)
        #pragma omp parallel for
        for (int k = 0; k<360; k++){
            double d = k * M_PI / 180; // ラジアン表記
            complex tmp;

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

            /////////////// r1 (復路) ///////////////
            complex r1m0z0, r1m0z1, r1m1z0, r1m1z1;
            r1m0z0 = mul(sin(d)/8, onenecp);
            r1m0z0.real += (1+cos(d))/4;
            r1m0z1 = mul(sin(d)/8, onepecn);
            r1m1z0 = mul(-sin(d)/8, onenecp);
            r1m1z0.real += (1-cos(d))/4;
            r1m1z1 = mul(-sin(d)/8, onepecn);

            /////////////// r2 (復路) ///////////////
            complex r2m0z0, r2m0z1, r2m1z0, r2m1z1;
            r2m0z0 = mul(-sin(d)/8, onenecp);
            r2m0z0.real += (1+cos(d))/4;
            r2m0z1 = mul(-sin(d)/8, onepecn);
            r2m1z0 = mul(sin(d)/8, onenecp);
            r2m1z0.real += (1-cos(d))/4;
            r2m1z1 = mul(sin(d)/8, onepecn);

            /////////////// r3 (復路) ///////////////
            complex r3m0z0, r3m0z1, r3m1z0, r3m1z1;
            r3m0z0 = mul(sin(d)/8, onepecp);
            r3m0z1 = mul(sin(d)/8, onenecn);
            r3m0z1.real += (1-cos(d))/4;
            r3m1z0 = mul(-sin(d)/8, onepecp);
            r3m1z1 = mul(-sin(d)/8, onenecn);
            r3m1z1.real += (1+cos(d))/4;

            /////////////// r4 (復路) ///////////////
            complex r4m0z0, r4m0z1, r4m1z0, r4m1z1;
            r4m0z0 = mul(-sin(d)/8, onepecp);
            r4m0z1 = mul(-sin(d)/8, onenecn);
            r4m0z1.real += (1-cos(d))/4;
            r4m1z0 = mul(sin(d)/8, onepecp);
            r4m1z1 = mul(sin(d)/8, onenecn);
            r4m1z1.real += (1+cos(d))/4;

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

            // 関数定義部分
            double f1,f2,f3,f4,f5,f6,f7,f8;
            f1 = Pccn(r1m0x0) + Pccn(r1m1x0);
            f2 = Pccn(r1m0x1) + Pccn(r1m1x1);
            f3 = Pccn(r2m0x0) + Pccn(r2m1x0);
            f4 = Pccn(r2m0x1) + Pccn(r2m1x1);
            f5 = Pccn(r3m0x0) + Pccn(r3m1x0);
            f6 = Pccn(r3m0x1) + Pccn(r3m1x1);
            f7 = Pccn(r4m0x0) + Pccn(r4m1x0);
            f8 = Pccn(r4m0x1) + Pccn(r4m1x1);

            double l1,l2,l3,l4,l5,l6,l7,l8;
            l1 = Pccn(r1m0y0) + Pccn(r1m1y0);
            l2 = Pccn(r1m0y1) + Pccn(r1m1y1);
            l3 = Pccn(r2m0y0) + Pccn(r2m1y0);
            l4 = Pccn(r2m0y1) + Pccn(r2m1y1);
            l5 = Pccn(r3m0y0) + Pccn(r3m1y0);
            l6 = Pccn(r3m0y1) + Pccn(r3m1y1);
            l7 = Pccn(r4m0y0) + Pccn(r4m1y0);
            l8 = Pccn(r4m0y1) + Pccn(r4m1y1);

            double u1,u2,u3,u4,u5,u6,u7,u8;
            u1 = Pccn(r1m0z0) + Pccn(r1m1z0);
            u2 = Pccn(r1m0z1) + Pccn(r1m1z1);
            u3 = Pccn(r2m0z0) + Pccn(r2m1z0);
            u4 = Pccn(r2m0z1) + Pccn(r2m1z1);
            u5 = Pccn(r3m0z0) + Pccn(r3m1z0);
            u6 = Pccn(r3m0z1) + Pccn(r3m1z1);
            u7 = Pccn(r4m0z0) + Pccn(r4m1z0);
            u8 = Pccn(r4m0z1) + Pccn(r4m1z1);

            // 関数定義部分
            double g1,g2,g3,g4,g7,g8,m1,m2,m3,m4,m5,m6,v1,v2,v5,v6,v7,v8;
            g1 = Pccn(r1m0x0) + Pccn(r1m0x1);
            g2 = Pccn(r1m1x0) + Pccn(r1m1x1);
            g3 = Pccn(r2m0x0) + Pccn(r2m0x1);
            g4 = Pccn(r2m1x0) + Pccn(r2m1x1);
            g7 = Pccn(r4m0x0) + Pccn(r4m0x1);
            g8 = Pccn(r4m1x0) + Pccn(r4m1x1);
            m1 = Pccn(r1m0y0) + Pccn(r1m0y1);
            m2 = Pccn(r1m1y0) + Pccn(r1m1y1);
            m3 = Pccn(r2m0y0) + Pccn(r2m0y1);
            m4 = Pccn(r2m1y0) + Pccn(r2m1y1);
            m5 = Pccn(r3m0y0) + Pccn(r3m0y1);
            m6 = Pccn(r3m1y0) + Pccn(r3m1y1);
            v1 = Pccn(r1m0z0) + Pccn(r1m0z1);
            v2 = Pccn(r1m1z0) + Pccn(r1m1z1);
            v5 = Pccn(r3m0z0) + Pccn(r3m0z1);
            v6 = Pccn(r3m1z0) + Pccn(r3m1z1);
            v7 = Pccn(r4m0z0) + Pccn(r4m0z1);
            v8 = Pccn(r4m1z0) + Pccn(r4m1z1);

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

            // xz OK
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
