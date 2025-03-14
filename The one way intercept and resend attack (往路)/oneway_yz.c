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

    file = fopen("FinalData_yz.csv", "w");

    if (file == NULL) {
        printf("ファイルを開けませんでした。\n");
        return 1;
    }
    fprintf(file, "a,b,PAB,PE,IAE,IAB\n");
    complex onep = {1, 1}; // (1+i)
    complex onen = {1, -1}; // (1-i)
    // i, j, k,lは角度。a,b,c,dはラジアン。
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

            // 関数定義部分
            double l3,l4,l7,l8,u3,u4,u7,u8;
            l3 = Pccn(r2y0g0) + Pccn(r2y0g1);
            l4 = Pccn(r2y1g0) + Pccn(r2y1g1);
            l7 = Pccn(r4y0g0) + Pccn(r4y0g1);
            l8 = Pccn(r4y1g0) + Pccn(r4y1g1);
            u3 = Pccn(r2z0g0) + Pccn(r2z0g1);
            u4 = Pccn(r2z1g0) + Pccn(r2z1g1);
            u7 = Pccn(r4z0g0) + Pccn(r4z0g1);
            u8 = Pccn(r4z1g0) + Pccn(r4z1g1);

            double l1,l2,l5,l6,u1,u2,u5,u6;
            l1 = Pccn(r1y0g0) + Pccn(r1y0g1);
            l2 = Pccn(r1y1g0) + Pccn(r1y1g1);
            l5 = Pccn(r3y0g0) + Pccn(r3y0g1);
            l6 = Pccn(r3y1g0) + Pccn(r3y1g1);
            u1 = Pccn(r1z0g0) + Pccn(r1z0g1);
            u2 = Pccn(r1z1g0) + Pccn(r1z1g1);
            u5 = Pccn(r3z0g0) + Pccn(r3z0g1);
            u6 = Pccn(r3z1g0) + Pccn(r3z1g1);

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

            double m1,m2,m5,m6,v1,v2,v5,v6;
            m1 = Pccn(r1y0g0) + Pccn(r1y1g0);
            m2 = Pccn(r1y0g1) + Pccn(r1y1g1);
            m5 = Pccn(r3y0g0) + Pccn(r3y1g0);
            m6 = Pccn(r3y0g1) + Pccn(r3y1g1);
            v1 = Pccn(r1z0g0) + Pccn(r1z1g0);
            v2 = Pccn(r1z0g1) + Pccn(r1z1g1);
            v5 = Pccn(r3z0g0) + Pccn(r3z1g0);
            v6 = Pccn(r3z0g1) + Pccn(r3z1g1);

            R = m1+m2+m5+m6+v1+v2+v5+v6;
            double R1 = (m1+m2+v1+v2)/R;
            double R3 = (m5+m6+v5+v6)/R;
            double PE = R1*(v1/(v1+v2) +m1/(m1+m2))/2 + R3*(v6/(v5+v6) + m6/(m5+m6))/2;

	          // AliceとEveの相互情報量
            double A0E0 = (v1+m1)/R;
            double A0E1 = (v2+m2)/R;
            double A1E0 = (v5+m5)/R;
            double A1E1 = (v6+m6)/R;
            double IAE = A0E0*log2(A0E0/((A0E0+A0E1)*(A0E0+A1E0)))+A0E1*log2(A0E1/((A0E0+A0E1)*(A0E1+A1E1)))+A1E0*log2(A1E0/((A1E0+A1E1)*(A0E0+A1E0)))+A1E1*log2(A1E1/((A1E0+A1E1)*(A0E1+A1E1)));


            double r1 = (m1+m2+v1+v2);
            double r2 = (l3+l4+u3+u4);
            double r3 = (m5+m6+v5+v6);
            double r4 = (l7+l8+u7+u8);
            

            //if (PAB >= 0.89999){
            //    printf("a = %d, b = %d, PAB=%lf, PE=%lf, IAE=%lf, IAB=%lf, %lf,%lf,%lf,%lf\n", i, k, PAB, PE, IAE, IAB,r1,r2,r3,r4);
            //}
                    
            // 探索判定部
            ///if (PE >= 0.92426){
            //    printf("a = %d, b = %d, PAB=%lf, PE=%lf, IAE=%lf, IAB=%lf, %lf, %lf, %lf, %lf\n", i, k, PAB, PE, IAE, IAB, r1,r2,r3,r4);
            //}

            //if (IAE > 0.612981){
            //    printf("a = %d, b = %d, PAB=%lf, PE=%lf, IAE=%lf, IAB=%lf, %lf, %lf, %lf, %lf\n", i, k, PAB, PE, IAE, IAB,r1,r2,r3,r4);
            //}

            if(fabs(r1 - r2) <= DBL_EPSILON * fmax(1, fmax(fabs(r1), fabs(r2))) && fabs(r1 - r3) <= DBL_EPSILON * fmax(1, fmax(fabs(r1), fabs(r3))) && fabs(r1 - r4) <= DBL_EPSILON * fmax(1, fmax(fabs(r1), fabs(r4))) && fabs(r2 - r3) <= DBL_EPSILON * fmax(1, fmax(fabs(r2), fabs(r3))) && fabs(r2 - r4) <= DBL_EPSILON * fmax(1, fmax(fabs(r2), fabs(r4))) && fabs(r3 - r4) <= DBL_EPSILON * fmax(1, fmax(fabs(r3), fabs(r4)))){
                fprintf(file, "%d,%d,%.32lf,%.32lf,%.32lf,%.32lf\n", i,k,PAB,PE,IAE,IAB);
            }
        }
    }
    fclose(file);
    return 0;
}
