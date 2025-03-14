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

    file = fopen("FD_yzalt(0~359).csv", "w");

    if (file == NULL) {
        printf("ファイルを開けませんでした。\n");
        return 1;
    }
    fprintf(file, "c,d,PAB,PE,IAE,IAB\n");
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

            // 関数定義部分
            double l3,l4,l7,l8,u3,u4,u7,u8;
            l3 = Pccn(r2m0y0) + Pccn(r2m1y0);
            l4 = Pccn(r2m0y1) + Pccn(r2m1y1);
            l7 = Pccn(r4m0y0) + Pccn(r4m1y0);
            l8 = Pccn(r4m0y1) + Pccn(r4m1y1);
            u3 = Pccn(r2m0z0) + Pccn(r2m1z0);
            u4 = Pccn(r2m0z1) + Pccn(r2m1z1);
            u7 = Pccn(r4m0z0) + Pccn(r4m1z0);
            u8 = Pccn(r4m0z1) + Pccn(r4m1z1);

            double l1,l2,l5,l6,u1,u2,u5,u6;
            l1 = Pccn(r1m0y0) + Pccn(r1m1y0);
            l2 = Pccn(r1m0y1) + Pccn(r1m1y1);
            l5 = Pccn(r3m0y0) + Pccn(r3m1y0);
            l6 = Pccn(r3m0y1) + Pccn(r3m1y1);
            u1 = Pccn(r1m0z0) + Pccn(r1m1z0);
            u2 = Pccn(r1m0z1) + Pccn(r1m1z1);
            u5 = Pccn(r3m0z0) + Pccn(r3m1z0);
            u6 = Pccn(r3m0z1) + Pccn(r3m1z1);
            

            double R = l3+l4+l7+l8+u3+u4+u7+u8;
            double R2 = (l3+l4+u3+u4)/R;
            double R4 = (l7+l8+u7+u8)/R;

            double PAB = R2*(l4/(l3+l4) + u3/(u3+u4))/2 + R4*(l7/(l7+l8) + u8/(u7+u8))/2;

            R = l1+l2+l5+l6+u1+u2+u5+u6;
            double A0B0 = (l1+u1)/R;
            double A0B1 = (l2+u2)/R;
            double A1B0 = (l5+u5)/R;
            double A1B1 = (l6+u6)/R;
            double IAB = A0B0*log2(A0B0/((A0B0+A0B1)*(A0B0+A1B0)))+A0B1*log2(A0B1/((A0B0+A0B1)*(A0B1+A1B1)))+A1B0*log2(A1B0/((A1B0+A1B1)*(A0B0+A1B0)))+A1B1*log2(A1B1/((A1B0+A1B1)*(A0B1+A1B1)));

            // 関数定義部分
            double m1,m2,m5,m6,v1,v2,v5,v6;
            m1 = Pccn(r1m0y0) + Pccn(r1m0y1);
            m2 = Pccn(r1m1y0) + Pccn(r1m1y1);
            m5 = Pccn(r3m0y0) + Pccn(r3m0y1);
            m6 = Pccn(r3m1y0) + Pccn(r3m1y1);
            v1 = Pccn(r1m0z0) + Pccn(r1m0z1);
            v2 = Pccn(r1m1z0) + Pccn(r1m1z1);
            v5 = Pccn(r3m0z0) + Pccn(r3m0z1);
            v6 = Pccn(r3m1z0) + Pccn(r3m1z1);


            R = m1+m2+m5+m6+v1+v2+v5+v6;
            double R1 = (m1+m2+v1+v2)/R;
            double R3 = (m5+m6+v5+v6)/R;

            double PE = R1*(v1/(v1+v2) +m1/(m1+m2))/2 + R3*(v6/(v5+v6) + m6/(m5+m6))/2;

            double A0E0 = (v1+m1)/R;
            double A0E1 = (v2+m2)/R;
            double A1E0 = (v5+m5)/R;
            double A1E1 = (v6+m6)/R;
            double IAE = A0E0*log2(A0E0/((A0E0+A0E1)*(A0E0+A1E0)))+A0E1*log2(A0E1/((A0E0+A0E1)*(A0E1+A1E1)))+A1E0*log2(A1E0/((A1E0+A1E1)*(A0E0+A1E0)))+A1E1*log2(A1E1/((A1E0+A1E1)*(A0E1+A1E1)));

            double r1 = (m1+m2+v1+v2);
            double r2 = (l3+l4+u3+u4);
            double r3 = (m5+m6+v5+v6);
            double r4 = (l7+l8+u7+u8);
            

            //if (PAB >= 0.89998){
            //    printf("a = %d, b = %d, PAB=%lf, PE=%lf, IAE=%lf, IAB=%lf, %lf, %lf, %lf, %lf\n", i, k, PAB, PE, IAE, IAB,r1,r2,r3,r4);
            //}
                    
            // 探索判定部
            //if (PE >= 0.92425){
            //    printf("a = %d, b = %d, PAB=%lf, PE=%lf\n", i, k, PAB,PE,IAE,IAB);
            //}

            // 探索判定部
            //if (IAE > 0.612979){
            //    printf("a = %d, b = %d, PAB=%lf, PE=%lf, IAE=%lf, IAB=%lf, %lf,%lf,%lf,%lf\n", i, k, PAB,PE, IAE, IAB, r1,r2,r3,r4);
            //}


            if(fabs(r1 - r2) <= DBL_EPSILON * fmax(1, fmax(fabs(r1), fabs(r2))) && fabs(r2 - r3) <= DBL_EPSILON * fmax(1, fmax(fabs(r2), fabs(r3))) && fabs(r3 - r4) <= DBL_EPSILON * fmax(1, fmax(fabs(r3), fabs(r4)))){
                if(PAB >= 0.875 || PE >= 0.875)
                    fprintf(file, "%d,%d,%.32lf,%.32lf\n", i,k,PAB,PE,IAE,IAB);
            }
        }
    }
    fclose(file);
    return 0;
}
