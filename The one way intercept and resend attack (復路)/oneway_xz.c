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

    file = fopen("FD_xzalt(0~359).csv", "w");

    if (file == NULL) {
        printf("ファイルを開けませんでした。\n");
        return 1;
    }
    fprintf(file, "a,b,PAB,PE,IAE,IAB\n");
    complex onep = {1, 1}; // (1+i)
    complex onen = {1, -1}; // (1-i)
    //printf("開始！\n");
    // i, j, k,lは角度。c,dはラジアン。
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

            double g1,g2,g7,g8,v1,v2,v7,v8;
            g1 = Pccn(r1m0x0) + Pccn(r1m0x1);
            g2 = Pccn(r1m1x0) + Pccn(r1m1x1);
            g7 = Pccn(r4m0x0) + Pccn(r4m0x1);
            g8 = Pccn(r4m1x0) + Pccn(r4m1x1);
            v1 = Pccn(r1m0z0) + Pccn(r1m0z1);
            v2 = Pccn(r1m1z0) + Pccn(r1m1z1);
            v7 = Pccn(r4m0z0) + Pccn(r4m0z1);
            v8 = Pccn(r4m1z0) + Pccn(r4m1z1);
            
            double R = g1+g2+g7+g8+v1+v2+v7+v8;
            double R1 = (g1+g2+v1+v2)/R;
            double R4 = (g7+g8+v7+v8)/R;
            
            double PE = R1*(g1/(g1+g2) + v1/(v1+v2))/2 + R4*(g8/(g7+g8) + v8/(v7+v8))/2;

            double f3,f4,u3,u4,f5,f6,u5,u6;
            f3 = Pccn(r2m0x0) + Pccn(r2m1x0);
            f4 = Pccn(r2m0x1) + Pccn(r2m1x1);
            f5 = Pccn(r3m0x0) + Pccn(r3m1x0);
            f6 = Pccn(r3m0x1) + Pccn(r3m1x1);
            u3 = Pccn(r2m0z0) + Pccn(r2m1z0);
            u4 = Pccn(r2m0z1) + Pccn(r2m1z1);
            u5 = Pccn(r3m0z0) + Pccn(r3m1z0);
            u6 = Pccn(r3m0z1) + Pccn(r3m1z1);

            double f1,f2,f7,f8,u1,u2,u7,u8;
            f1 = Pccn(r1m0x0) + Pccn(r1m1x0);
            f2 = Pccn(r1m0x1) + Pccn(r1m1x1);
            f7 = Pccn(r4m0x0) + Pccn(r4m1x0);
            f8 = Pccn(r4m0x1) + Pccn(r4m1x1);
            u1 = Pccn(r1m0z0) + Pccn(r1m1z0);
            u2 = Pccn(r1m0z1) + Pccn(r1m1z1);
            u7 = Pccn(r4m0z0) + Pccn(r4m1z0);
            u8 = Pccn(r4m0z1) + Pccn(r4m1z1);
            
            R = f3+f4+f5+f6+u3+u4+u5+u6;
            double R2 = (f3+f4+u3+u4)/R;
            double R3 = (f5+f6+u5+u6)/R;
            
            double PAB = R2*(f4/(f3+f4) + u3/(u3+u4))/2 + R3*(f5/(f5+f6) + u6/(u5+u6))/2;

            // I(A;E)
            R = g1+g2+g7+g8+v1+v2+v7+v8;
            double A0E0 = (g1+v1)/R;
            double A0E1 = (g2+v2)/R;
            double A1E0 = (g7+v7)/R;
            double A1E1 = (g8+v8)/R;
            double IAE = A0E0*log2(A0E0/((A0E0+A0E1)*(A0E0+A1E0)))+A0E1*log2(A0E1/((A0E0+A0E1)*(A0E1+A1E1)))+A1E0*log2(A1E0/((A1E0+A1E1)*(A0E0+A1E0)))+A1E1*log2(A1E1/((A1E0+A1E1)*(A0E1+A1E1)));

            // I(A;B)
            R = f1+f2+f7+f8+u1+u2+u7+u8;
            double A0B0 = (f1+u1)/R;
            double A0B1 = (f2+u2)/R;
            double A1B0 = (f7+u7)/R;
            double A1B1 = (f8+u8)/R;
            double IAB = A0B0*log2(A0B0/((A0B0+A0B1)*(A0B0+A1B0)))+A0B1*log2(A0B1/((A0B0+A0B1)*(A0B1+A1B1)))+A1B0*log2(A1B0/((A1B0+A1B1)*(A0B0+A1B0)))+A1B1*log2(A1B1/((A1B0+A1B1)*(A0B1+A1B1)));

            double r1,r2,r3,r4;
            r1 = g1+g2+v1+v2;
            r2 = f3+f4+u3+u4;
            r3 = f5+f6+u5+u6;
            r4 = g7+g8+v7+v8;

            //if (PAB >= 0.899987){
            //    printf("a = %d, b = %d, PAB=%lf, PE=%lf, IAE=%lf, IAB=%lf, %lf,%lf,%lf,%lf\n", i, k, PAB, PE, IAE, IAB, r1,r2,r3,r4);
            //}
                    
            // 探索判定部
            //if (PE >= 0.924252){
            //    printf("a = %d, b = %d, PAB=%lf, PE=%lf, IAE=%lf, IAB=%lf, %lf,%lf,%lf,%lf\n", i, k, PAB,PE, IAE, IAB, r1,r2,r3,r4);
            //}
            

            //if (0.612981 <IAE){
            //    printf("a = %d, b = %d, PAB=%lf, PE=%lf, IAE=%lf, IAB=%lf, %lf,%lf,%lf,%lf\n", i, k, PAB,PE, IAE, IAB, r1,r2,r3,r4);
            //}

            if(fabs(r1 - r2) <= DBL_EPSILON * fmax(1, fmax(fabs(r1), fabs(r2))) && fabs(r1 - r3) <= DBL_EPSILON * fmax(1, fmax(fabs(r1), fabs(r3))) && fabs(r1 - r4) <= DBL_EPSILON * fmax(1, fmax(fabs(r1), fabs(r4))) && fabs(r2 - r3) <= DBL_EPSILON * fmax(1, fmax(fabs(r2), fabs(r3))) && fabs(r2 - r4) <= DBL_EPSILON * fmax(1, fmax(fabs(r2), fabs(r4))) && fabs(r3 - r4) <= DBL_EPSILON * fmax(1, fmax(fabs(r3), fabs(r4)))){
                fprintf(file, "%d,%d,%.32lf,%.32lf,%.32lf,%.32lf\n", i,k,PAB,PE,IAE,IAB);
            }

        }
    }
    fclose(file);
    return 0;
}
