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

    file = fopen("FinalData_xy.csv", "w");

    if (file == NULL) {
        printf("ファイルを開けませんでした。\n");
        return 1;
    }
    fprintf(file, "a,b,PAB,PE,IAE,IAB\n");
    complex onep = {1, 1}; // (1+i)
    complex onen = {1, -1}; // (1-i)
    //printf("開始！\n");
    // i,kは角度。a,bはラジアン。
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

            // PAB定義部分
            double f5,f6,f7,f8,l5,l6,l7,l8;
            f5 = Pccn(r3x0g0) + Pccn(r3x0g1);
            f6 = Pccn(r3x1g0) + Pccn(r3x1g1);
            f7 = Pccn(r4x0g0) + Pccn(r4x0g1);
            f8 = Pccn(r4x1g0) + Pccn(r4x1g1);
            l5 = Pccn(r3y0g0) + Pccn(r3y0g1);
            l6 = Pccn(r3y1g0) + Pccn(r3y1g1);
            l7 = Pccn(r4y0g0) + Pccn(r4y0g1);
            l8 = Pccn(r4y1g0) + Pccn(r4y1g1); 

	    double R = f5+f6+f7+f8+l5+l6+l7+l8;
            double R3 = (f5+f6+l5+l6)/R;
            double R4 = (f7+f8+l7+l8)/R;
	    double PAB = R3*(f5/(f5+f6)+l6/(l5+l6))/2+R4*(f8/(f7+f8)+l7/(l7+l8))/2;

	    // PE定義部分
            double g1,g2,g3,g4,m1,m2,m3,m4;
            g1 = Pccn(r1x0g0) + Pccn(r1x1g0);
            g2 = Pccn(r1x0g1) + Pccn(r1x1g1);
            g3 = Pccn(r2x0g0) + Pccn(r2x1g0);
            g4 = Pccn(r2x0g1) + Pccn(r2x1g1);
            m1 = Pccn(r1y0g0) + Pccn(r1y1g0);
            m2 = Pccn(r1y0g1) + Pccn(r1y1g1);
            m3 = Pccn(r2y0g0) + Pccn(r2y1g0);
            m4 = Pccn(r2y0g1) + Pccn(r2y1g1);

            R = g1+g2+g3+g4+m1+m2+m3+m4;
            double R1 = (g1+g2+m1+m2)/R;
            double R2 = (g3+g4+m3+m4)/R;
	    double PE = R1*(g1/(g1+g2)+m1/(m1+m2))/2+R2*(g4/(g3+g4)+m4/(m3+m4))/2;

            // Aliceがr_aを得る確率
	    double r1 = m1+m2+g1+g2;
            double r2 = m3+m4+g3+g4;
            double r3 = l5+l6+f5+f6;
            double r4 = l7+l8+f7+f8;
		
            // IAEとIABの導出
	    double f1,f2,f3,f4,l1,l2,l3,l4;
            f1 = Pccn(r1x0g0) + Pccn(r1x0g1);
            f2 = Pccn(r1x1g0) + Pccn(r1x1g1);
            f3 = Pccn(r2x0g0) + Pccn(r2x0g1);
            f4 = Pccn(r2x1g0) + Pccn(r2x1g1);
            l1 = Pccn(r1y0g0) + Pccn(r1y0g1);
            l2 = Pccn(r1y1g0) + Pccn(r1y1g1);
            l3 = Pccn(r2y0g0) + Pccn(r2y0g1);
            l4 = Pccn(r2y1g0) + Pccn(r2y1g1);
            
            double R12 = g1+g2+g3+g4+m1+m2+m3+m4;
            double A0E0 = (g1+m1)/R12;
            double A0E1 = (g2+m2)/R12;
            double A1E0 = (g3+m3)/R12;
            double A1E1 = (g4+m4)/R12;
            double IAE = A0E0*log2(A0E0/((A0E0+A0E1)*(A0E0+A1E0)))+A0E1*log2(A0E1/((A0E0+A0E1)*(A0E1+A1E1)))+A1E0*log2(A1E0/((A1E0+A1E1)*(A0E0+A1E0)))+A1E1*log2(A1E1/((A1E0+A1E1)*(A0E1+A1E1)));

            R12 = f1+f2+f3+f4+l1+l2+l3+l4;
            double A0B0 = (f1+l1)/R12;
            double A0B1 = (f2+l2)/R12;
            double A1B0 = (f3+l3)/R12;
            double A1B1 = (f4+l4)/R12;
            double IAB = A0B0*log2(A0B0/((A0B0+A0B1)*(A0B0+A1B0)))+A0B1*log2(A0B1/((A0B0+A0B1)*(A0B1+A1B1)))+A1B0*log2(A1B0/((A1B0+A1B1)*(A0B0+A1B0)))+A1B1*log2(A1B1/((A1B0+A1B1)*(A0B1+A1B1)));

            //if (PAB >= 0.899999){
            //    printf("a = %d, b = %d, PAB=%lf, PE=%lf, IAE=%lf,IAB=%lf, %lf, %lf, %lf, %lf\n", i, k, PAB, PE, IAE,IAB,r1,r2,r3,r4);
            //}
                    
            // 探索判定部
            //if (PE >= 0.92426){
            //    printf("a = %d, b = %d, PAB=%lf, PE=%lf, IAE=%lf, IAB=%lf, %lf,%lf,%lf,%lf\n", i, k, PAB, PE, IAE, IAB, r1,r2,r3,r4);
            //}

            //if (IAE > 0.61302){
            //    printf("a = %d, b = %d, PAB=%lf, PE=%lf, IAE=%lf, IAB=%lf, %lf,%lf,%lf,%lf\n", i, k, PAB, PE, IAE, IAB,r1,r2,r3,r4);
            //}

            if(fabs(r1 - r2) <= DBL_EPSILON * fmax(1, fmax(fabs(r1), fabs(r2))) && fabs(r1 - r3) <= DBL_EPSILON * fmax(1, fmax(fabs(r1), fabs(r3))) && fabs(r1 - r4) <= DBL_EPSILON * fmax(1, fmax(fabs(r1), fabs(r4))) && fabs(r2 - r3) <= DBL_EPSILON * fmax(1, fmax(fabs(r2), fabs(r3))) && fabs(r2 - r4) <= DBL_EPSILON * fmax(1, fmax(fabs(r2), fabs(r4))) && fabs(r3 - r4) <= DBL_EPSILON * fmax(1, fmax(fabs(r3), fabs(r4)))){
                fprintf(file, "%d,%d,%.32lf,%.32lf,%.32lf,%.32lf\n", i,k,PAB,PE,IAE,IAB);
            }
        }
    }
    fclose(file);
    return 0;
}
