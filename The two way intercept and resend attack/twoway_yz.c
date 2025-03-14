#define _USE_MATH_DEFINES
#include <stdio.h>
#include <math.h>
#include <omp.h>

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
    FILE *file;

    file = fopen("FinalData_twoway_yz.csv", "w");

    if (file == NULL) {
        printf("ファイルを開けませんでした。\n");
        return 1;
    }
    fprintf(file, "a,b,c,d,PAB,PE,IAE,IAB\n");
    printf("開始！\n");
    // i, j, k,lは角度。a,b,c,dはラジアン。
    for (int i = 0; i < 360; i++){
        double a = i * M_PI / 180; // ラジアン表記
        complex eap = {cos(a), sin(a)}; // e^(iα)
        complex ean = {cos(a), -sin(a)}; // e^(-iα)
        complex onepeap = pro(onep,eap); // (1+i)e^(iα)
        complex oneneap = pro(onen, eap); // (1-i)e^(iα)
        complex onepean = pro(onep,ean); // (1+i)e^(-iα)
        complex onenean = pro(onen, ean); // (1-i)e^(-iα)
        printf("Now %d\n",i);
        #pragma omp parallel for
        for (int k = 0; k<360; k++){
            double c = k * M_PI / 180; // ラジアン表記
            complex ecp = {cos(c), sin(c)}; // e^(iγ)
            complex ecn = {cos(c), -sin(c)}; // e^(-iγ)
            complex onepecp = pro(onep,ecp); // (1+i)e^(iγ)
            complex onenecp = pro(onen,ecp); // (1-i)e^(iγ)
            complex onepecn = pro(onep,ecn); // (1+i)e^(-iγ)
            complex onenecn = pro(onen,ecn); // (1-i)e^(-iγ)
            complex eapcp = {cos(a+c), sin(a+c)}; // e^(iα+iγ)
            complex eapcn = {cos(a-c), sin(a-c)}; // e^(iα-iγ)
            complex eancp = {cos(a-c), -sin(a-c)}; // e^(-iα+iγ)
            complex eancn = {cos(a+c), -sin(a+c)}; // e^(-iα-iγ)
            complex ecpi = {-sin(c), cos(c)}; // ie^(iγ)
            complex ecni = {sin(c), cos(c)}; // ie^(-iγ)
            complex eapi = {-sin(a), cos(a)}; // ie^(iα)
            complex eani = {sin(a), cos(a)}; // ie^(-iα)
            complex onepeapcp = pro(onep,eapcp); // (1+i)e^(iα+iγ)
            complex oneneapcp = pro(onen,eapcp); // (1-i)e^(iα+iγ)
            complex onepeapcn = pro(onep,eapcn); // (1+i)e^(iα-iγ)
            complex oneneapcn = pro(onen,eapcn); // (1-i)e^(iα-iγ)
            complex onepeancp = pro(onep,eancp); // (1+i)e^(-iα+iγ)
            complex oneneancp = pro(onen,eancp); // (1-i)e^(-iα+iγ)
            complex onepeancn = pro(onep,eancn); // (1+i)e^(-iα-iγ)
            complex oneneancn = pro(onen,eancn); // (1-i)e^(-iα-iγ)
            for (int j = 0; j<360; j++){
                double b = j * M_PI / 180; // ラジアン表記
                for  (int l = 0; l<360; l++){
                    double d = l * M_PI / 180; // ラジアン表記
                    complex tmp;

                    /////////////// r1 ///////////////
                    // 5 bd
                    tmp = sum(mul((1+cos(b))*sin(d)/16, onenecp), mul(sin(b)*(1+cos(d))/16, onepean));//OK
                    complex r1m0z0g0 = {(1+cos(b))*(1+cos(d))/8, 0};//OK
                    r1m0z0g0 = sum(r1m0z0g0, tmp);//OK
                    // 5 bd
                    tmp = sum(mul((1-cos(b))*sin(d)/16, onenecp), mul(-sin(b)*(1+cos(d))/16, onepean));//OK
                    complex r1m0z0g1 = {(1-cos(b))*(1+cos(d))/8, 0};//OK
                    r1m0z0g1 = sum(r1m0z0g1, tmp);//OK
                    // 5 bd
                    tmp = sum(mul((1-cos(b))*sin(d)/16, onepecn), mul(sin(b)*(1-cos(d))/16, oneneap));//OK
                    complex r1m0z1g0= mul(sin(b)*sin(d)/8, eapcn);//OK
                    r1m0z1g0 = sum(r1m0z1g0, tmp);//OK
                    // 5 bd
                    tmp = sum(mul((1+cos(b))*sin(d)/16, onepecn), mul(-sin(b)*(1-cos(d))/16, oneneap));//OK
                    complex r1m0z1g1= mul(-sin(b)*sin(d)/8, eapcn);//OK
                    r1m0z1g1 = sum(r1m0z1g1, tmp);//OK
                    // 5 bd
                    tmp = sum(mul(-(1+cos(b))*sin(d)/16, onenecp), mul(sin(b)*(1-cos(d))/16, onepean));//OK
                    complex r1m1z0g0 = {(1+cos(b))*(1-cos(d))/8, 0};//OK
                    r1m1z0g0 = sum(r1m1z0g0, tmp);
                    // 5 bd
                    tmp = sum(mul(-(1-cos(b))*sin(d)/16, onenecp), mul(-sin(b)*(1-cos(d))/16, onepean)); //OK
                    complex r1m1z0g1 = {(1-cos(b))*(1-cos(d))/8, 0};//OK
                    r1m1z0g1 = sum(r1m1z0g1, tmp);//OK
                    // 5 bd
                    tmp = sum(mul(-(1-cos(b))*sin(d)/16, onepecn), mul(sin(b)*(1+cos(d))/16, oneneap)); //OK
                    complex r1m1z1g0 = mul(-sin(b)*sin(d)/8, eapcn);//OK
                    r1m1z1g0 = sum(r1m1z1g0, tmp);//OK
                    // 5 bd
                    tmp = sum(mul(-(1+cos(b))*sin(d)/16, onepecn), mul(-sin(b)*(1+cos(d))/16, oneneap)); // OK
                    complex r1m1z1g1 = mul(sin(b)*sin(d)/8, eapcn); // OK
                    r1m1z1g1 = sum(r1m1z1g1, tmp); // OK
                    /////////////// r2 ///////////////
                    // 5 bd
                    tmp = sum(mul(-(1+cos(b))*sin(d)/16, onenecp), mul(-sin(b)*(1+cos(d))/16, onepean));//OK
                    complex r2m0z0g0 = {(1+cos(b))*(1+cos(d))/8, 0};//OK
                    r2m0z0g0 = sum(r2m0z0g0, tmp);//OK
                    // 5 bd
                    tmp = sum(mul(-(1-cos(b))*sin(d)/16, onenecp), mul(sin(b)*(1+cos(d))/16, onepean));//OK
                    complex r2m0z0g1 = {(1-cos(b))*(1+cos(d))/8, 0};//OK
                    r2m0z0g1 = sum(r2m0z0g1, tmp);//OK
                    // 5 bd
                    tmp = sum(mul(-(1-cos(b))*sin(d)/16, onepecn), mul(-sin(b)*(1-cos(d))/16, oneneap));//OK
                    complex r2m0z1g0= mul(sin(b)*sin(d)/8, eapcn);//OK
                    r2m0z1g0 = sum(r2m0z1g0, tmp);//OK
                    // 5 bd
                    tmp = sum(mul(-(1+cos(b))*sin(d)/16, onepecn), mul(sin(b)*(1-cos(d))/16, oneneap));//OK
                    complex r2m0z1g1= mul(-sin(b)*sin(d)/8, eapcn);//OK
                    r2m0z1g1 = sum(r2m0z1g1, tmp);//OK
                    // 5 bd
                    tmp = sum(mul((1+cos(b))*sin(d)/16, onenecp), mul(-sin(b)*(1-cos(d))/16, onepean));//OK
                    complex r2m1z0g0 = {(1+cos(b))*(1-cos(d))/8, 0};//OK
                    r2m1z0g0 = sum(r2m1z0g0, tmp);
                    // 5 bd
                    tmp = sum(mul((1-cos(b))*sin(d)/16, onenecp), mul(sin(b)*(1-cos(d))/16, onepean)); //OK
                    complex r2m1z0g1 = {(1-cos(b))*(1-cos(d))/8, 0};//OK
                    r2m1z0g1 = sum(r2m1z0g1, tmp);//OK
                    // 5 bd
                    tmp = sum(mul((1-cos(b))*sin(d)/16, onepecn), mul(-sin(b)*(1+cos(d))/16, oneneap)); //OK
                    complex r2m1z1g0= mul(-sin(b)*sin(d)/8, eapcn);//OK
                    r2m1z1g0 = sum(r2m1z1g0, tmp);//OK
                    // 5 bd
                    tmp = sum(mul((1+cos(b))*sin(d)/16, onepecn), mul(sin(b)*(1+cos(d))/16, oneneap)); // OK
                    complex r2m1z1g1 = mul(sin(b)*sin(d)/8, eapcn); // OK
                    r2m1z1g1 = sum(r2m1z1g1, tmp); // OK
                    /////////////// r3 ///////////////
                    // 5 bd
                    tmp = sum(mul((1+cos(b))*sin(d)/16, onepecp), mul(sin(b)*(1+cos(d))/16, onenean));// OK
                    complex r3m0z0g0 = mul(sin(b)*sin(d)/8, eancp); // OK
                    r3m0z0g0 = sum(r3m0z0g0, tmp);
                    // 5 bd
                    tmp = sum(mul((1-cos(b))*sin(d)/16, onepecp), mul(-sin(b)*(1+cos(d))/16, onenean));// OK
                    complex r3m0z0g1 = mul(-sin(b)*sin(d)/8, eancp); // OK
                    r3m0z0g1 = sum(r3m0z0g1, tmp);
                    // 5 bd
                    tmp = sum(mul((1-cos(b))*sin(d)/16, onenecn), mul(sin(b)*(1-cos(d))/16, onepeap));// OK
                    complex r3m0z1g0 = {(1-cos(b))*(1-cos(d))/8, 0};// OK
                    r3m0z1g0 = sum(r3m0z1g0, tmp);
                    // 5 bd
                    tmp = sum(mul((1+cos(b))*sin(d)/16, onenecn), mul(-sin(b)*(1-cos(d))/16, onepeap));// OK
                    complex r3m0z1g1 = {(1+cos(b))*(1-cos(d))/8,0};// OK
                    r3m0z1g1 = sum(r3m0z1g1, tmp);
                    // 5 bd
                    tmp = sum(mul(-(1+cos(b))*sin(d)/16, onepecp), mul(sin(b)*(1-cos(d))/16, onenean));// OK
                    complex r3m1z0g0 = mul(-sin(b)*sin(d)/8, eancp);// OK
                    r3m1z0g0 = sum(r3m1z0g0, tmp);
                    // 5 bd
                    tmp = sum(mul(-(1-cos(b))*sin(d)/16, onepecp), mul(-sin(b)*(1-cos(d))/16, onenean));// OK
                    complex r3m1z0g1 = mul(sin(b)*sin(d)/8, eancp);// OK
                    r3m1z0g1 = sum(r3m1z0g1, tmp);
                    // 5 bd
                    tmp = sum(mul(-(1-cos(b))*sin(d)/16, onenecn), mul(sin(b)*(1+cos(d))/16, onepeap));//OK
                    complex r3m1z1g0 = {(1-cos(b))*(1+cos(d))/8, 0};//OK
                    r3m1z1g0 = sum(r3m1z1g0, tmp);
                    // 5 bd
                    tmp = sum(mul(-(1+cos(b))*sin(d)/16, onenecn), mul(-sin(b)*(1+cos(d))/16, onepeap));//OK
                    complex r3m1z1g1 = {(1+cos(b))*(1+cos(d))/8, 0};
                    r3m1z1g1 = sum(r3m1z1g1, tmp);
                    /////////////// r4 ///////////////
                    // 5 bd
                    tmp = sum(mul(-(1+cos(b))*sin(d)/16, onepecp), mul(-sin(b)*(1+cos(d))/16, onenean));// OK
                    complex r4m0z0g0 = mul(sin(b)*sin(d)/8, eancp); // OK
                    r4m0z0g0 = sum(r4m0z0g0, tmp);
                    // 5 bd
                    tmp = sum(mul(-(1-cos(b))*sin(d)/16, onepecp), mul(sin(b)*(1+cos(d))/16, onenean));// OK
                    complex r4m0z0g1 = mul(-sin(b)*sin(d)/8, eancp); // OK
                    r4m0z0g1 = sum(r4m0z0g1, tmp);
                    // 5 bd
                    tmp = sum(mul(-(1-cos(b))*sin(d)/16, onenecn), mul(-sin(b)*(1-cos(d))/16, onepeap));// OK
                    complex r4m0z1g0 = {(1-cos(b))*(1-cos(d))/8, 0};// OK
                    r4m0z1g0 = sum(r4m0z1g0, tmp);
                    // 5 bd
                    tmp = sum(mul(-(1+cos(b))*sin(d)/16, onenecn), mul(sin(b)*(1-cos(d))/16, onepeap));// OK
                    complex r4m0z1g1 = {(1+cos(b))*(1-cos(d))/8, 0};// OK
                    r4m0z1g1 = sum(r4m0z1g1, tmp);
                    // 5 bd
                    tmp = sum(mul((1+cos(b))*sin(d)/16, onepecp), mul(-sin(b)*(1-cos(d))/16, onenean));// OK
                    complex r4m1z0g0 = mul(-sin(b)*sin(d)/8, eancp);// OK
                    r4m1z0g0 = sum(r4m1z0g0, tmp);
                    // 5 bd
                    tmp = sum(mul((1-cos(b))*sin(d)/16, onepecp), mul(sin(b)*(1-cos(d))/16, onenean));// OK
                    complex r4m1z0g1 = mul(sin(b)*sin(d)/8, eancp);// OK
                    r4m1z0g1 = sum(r4m1z0g1, tmp);
                    // 5 bd
                    tmp = sum(mul((1-cos(b))*sin(d)/16, onenecn), mul(-sin(b)*(1+cos(d))/16, onepeap));// OK
                    complex r4m1z1g0 = {(1-cos(b))*(1+cos(d))/8, 0};// OK
                    r4m1z1g0 = sum(r4m1z1g0, tmp);
                    // 5 bd
                    tmp = sum(mul((1+cos(b))*sin(d)/16, onenecn), mul(sin(b)*(1+cos(d))/16, onepeap));// OK
                    complex r4m1z1g1 = {(1+cos(b))*(1+cos(d))/8, 0};// OK
                    r4m1z1g1 = sum(r4m1z1g1, tmp);

                    /////////////// r1 ///////////////
                    // 5 bd
                    complex r201_r202   = sum(mul((1-cos(b))*sin(d)/32, onepecn), mul((1-cos(b))*(1+cos(d))/32, onen)); // OK
                    complex r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneancn), mul(sin(b)*(1+cos(d))/32, onepean)); // OK
                    complex r205_r206   = sum(mul(sin(b)*(1-cos(d))/32, oneneap), mul(-sin(b)*sin(d)/32,onepeapcp));// OK
                    complex r207_r208   = sum(mul((1+cos(b))*(1-cos(d))/32,onep), mul((1+cos(b))*sin(d)/32,onenecp));// OK
                    complex r209_r2010  = sum(mul(sin(b)*sin(d)/16, eapcn), mul(-sin(b)*(1+cos(d))/16,eapi));// OK
                    complex r2011         = mul((1+cos(b))*sin(d)/16, ecni); // OK
                    complex r1m0y0g0 = {(1+cos(b))*(1+cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r1m0y0g0 = sum(sum(r1m0y0g0, tmp), sum(r209_r2010, r2011));//OK

                    // 5 bd
                    r201_r202   = sum(mul((1+cos(b))*sin(d)/32, onepecn), mul((1+cos(b))*(1+cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneancn), mul(-sin(b)*(1+cos(d))/32, onepean)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1-cos(d))/32, oneneap), mul(sin(b)*sin(d)/32,onepeapcp));// OK
                    r207_r208   = sum(mul((1-cos(b))*(1-cos(d))/32,onep), mul((1-cos(b))*sin(d)/32,onenecp));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eapcn), mul(sin(b)*(1+cos(d))/16,eapi));// OK
                    r2011         = mul((1-cos(b))*sin(d)/16, ecni); // OK
                    complex r1m0y0g1 = {(1-cos(b))*(1+cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r1m0y0g1 = sum(sum(r1m0y0g1, tmp), sum(r209_r2010, r2011));//OK

                    // 5 bd
                    r201_r202   = sum(mul((1-cos(b))*sin(d)/32, onepecn), mul(-(1-cos(b))*(1+cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneancn), mul(sin(b)*(1+cos(d))/32, onepean)); // OK
                    r205_r206   = sum(mul(sin(b)*(1-cos(d))/32, oneneap), mul(sin(b)*sin(d)/32,onepeapcp));// OK
                    r207_r208   = sum(mul(-(1+cos(b))*(1-cos(d))/32,onep), mul((1+cos(b))*sin(d)/32,onenecp));// OK
                    r209_r2010  = sum(mul(sin(b)*sin(d)/16, eapcn), mul(sin(b)*(1+cos(d))/16,eapi));// OK
                    r2011         = mul(-(1+cos(b))*sin(d)/16, ecni); // OK
                    complex r1m0y1g0 = {(1+cos(b))*(1+cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r1m0y1g0 = sum(sum(r1m0y1g0, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul((1+cos(b))*sin(d)/32, onepecn), mul(-(1+cos(b))*(1+cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneancn), mul(-sin(b)*(1+cos(d))/32, onepean)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1-cos(d))/32, oneneap), mul(-sin(b)*sin(d)/32,onepeapcp));// OK
                    r207_r208   = sum(mul(-(1-cos(b))*(1-cos(d))/32,onep), mul((1-cos(b))*sin(d)/32,onenecp));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eapcn), mul(-sin(b)*(1+cos(d))/16,eapi));// OK
                    r2011         = mul(-(1-cos(b))*sin(d)/16, ecni); // OK
                    complex r1m0y1g1 = {(1-cos(b))*(1+cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r1m0y1g1 = sum(sum(r1m0y1g1, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul(-(1-cos(b))*sin(d)/32, onepecn), mul((1-cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneancn), mul(sin(b)*(1-cos(d))/32, onepean)); // OK
                    r205_r206   = sum(mul(sin(b)*(1+cos(d))/32, oneneap), mul(sin(b)*sin(d)/32,onepeapcp));// OK
                    r207_r208   = sum(mul((1+cos(b))*(1+cos(d))/32,onep), mul(-(1+cos(b))*sin(d)/32,onenecp));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eapcn), mul(-sin(b)*(1-cos(d))/16,eapi));// OK
                    r2011         = mul(-(1+cos(b))*sin(d)/16, ecni); // OK
                    complex r1m1y0g0 = {(1+cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r1m1y0g0 = sum(sum(r1m1y0g0, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul(-(1+cos(b))*sin(d)/32, onepecn), mul((1+cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneancn), mul(-sin(b)*(1-cos(d))/32, onepean)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1+cos(d))/32, oneneap), mul(-sin(b)*sin(d)/32,onepeapcp));// OK
                    r207_r208   = sum(mul((1-cos(b))*(1+cos(d))/32,onep), mul(-(1-cos(b))*sin(d)/32,onenecp));// OK
                    r209_r2010  = sum(mul(sin(b)*sin(d)/16, eapcn), mul(sin(b)*(1-cos(d))/16,eapi));// OK
                    r2011         = mul(-(1-cos(b))*sin(d)/16, ecni); // OK
                    complex r1m1y0g1 = {(1-cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r1m1y0g1 = sum(sum(r1m1y0g1, tmp), sum(r209_r2010, r2011));//OK

                    // 5 bd
                    r201_r202   = sum(mul(-(1-cos(b))*sin(d)/32, onepecn), mul(-(1-cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneancn), mul(sin(b)*(1-cos(d))/32, onepean)); // OK
                    r205_r206   = sum(mul(sin(b)*(1+cos(d))/32, oneneap), mul(-sin(b)*sin(d)/32,onepeapcp));// OK
                    r207_r208   = sum(mul(-(1+cos(b))*(1+cos(d))/32,onep), mul(-(1+cos(b))*sin(d)/32,onenecp));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eapcn), mul(sin(b)*(1-cos(d))/16,eapi));// OK
                    r2011         = mul((1+cos(b))*sin(d)/16, ecni); // OK
                    complex r1m1y1g0 = {(1+cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r1m1y1g0 = sum(sum(r1m1y1g0, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul(-(1+cos(b))*sin(d)/32, onepecn), mul(-(1+cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneancn), mul(-sin(b)*(1-cos(d))/32, onepean)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1+cos(d))/32, oneneap), mul(sin(b)*sin(d)/32,onepeapcp));// OK
                    r207_r208   = sum(mul(-(1-cos(b))*(1+cos(d))/32,onep), mul(-(1-cos(b))*sin(d)/32,onenecp));// OK
                    r209_r2010  = sum(mul(sin(b)*sin(d)/16, eapcn), mul(-sin(b)*(1-cos(d))/16,eapi));// OK
                    r2011         = mul((1-cos(b))*sin(d)/16, ecni); // OK
                    complex r1m1y1g1 = {(1-cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r1m1y1g1 = sum(sum(r1m1y1g1, tmp), sum(r209_r2010, r2011));// OK

                    /////////////// r2 ///////////////
                    // 5 bd
                    r201_r202   = sum(mul(-(1-cos(b))*sin(d)/32, onepecn), mul(-(1-cos(b))*(1+cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneancn), mul(-sin(b)*(1+cos(d))/32, onepean)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1-cos(d))/32, oneneap), mul(sin(b)*sin(d)/32,onepeapcp));// OK
                    r207_r208   = sum(mul(-(1+cos(b))*(1-cos(d))/32,onep), mul(-(1+cos(b))*sin(d)/32,onenecp));// OK
                    r209_r2010  = sum(mul(sin(b)*sin(d)/16, eapcn), mul(-sin(b)*(1+cos(d))/16,eapi));// OK
                    r2011         = mul((1+cos(b))*sin(d)/16, ecni); // OK
                    complex r2m0y0g0 = {(1+cos(b))*(1+cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r2m0y0g0 = sum(sum(r2m0y0g0, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul(-(1+cos(b))*sin(d)/32, onepecn), mul(-(1+cos(b))*(1+cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneancn), mul(sin(b)*(1+cos(d))/32, onepean)); // OK
                    r205_r206   = sum(mul(sin(b)*(1-cos(d))/32, oneneap), mul(-sin(b)*sin(d)/32,onepeapcp));// OK
                    r207_r208   = sum(mul(-(1-cos(b))*(1-cos(d))/32,onep), mul(-(1-cos(b))*sin(d)/32,onenecp));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eapcn), mul(sin(b)*(1+cos(d))/16,eapi));// OK
                    r2011         = mul((1-cos(b))*sin(d)/16, ecni); // OK
                    complex r2m0y0g1 = {(1-cos(b))*(1+cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r2m0y0g1 = sum(sum(r2m0y0g1, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul(-(1-cos(b))*sin(d)/32, onepecn), mul((1-cos(b))*(1+cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneancn), mul(-sin(b)*(1+cos(d))/32, onepean)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1-cos(d))/32, oneneap), mul(-sin(b)*sin(d)/32,onepeapcp));// OK
                    r207_r208   = sum(mul((1+cos(b))*(1-cos(d))/32,onep), mul(-(1+cos(b))*sin(d)/32,onenecp));// OK
                    r209_r2010  = sum(mul(sin(b)*sin(d)/16, eapcn), mul(sin(b)*(1+cos(d))/16,eapi));// OK
                    r2011         = mul(-(1+cos(b))*sin(d)/16, ecni); // OK
                    complex r2m0y1g0 = {(1+cos(b))*(1+cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r2m0y1g0 = sum(sum(r2m0y1g0, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul(-(1+cos(b))*sin(d)/32, onepecn), mul((1+cos(b))*(1+cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneancn), mul(sin(b)*(1+cos(d))/32, onepean)); // OK
                    r205_r206   = sum(mul(sin(b)*(1-cos(d))/32, oneneap), mul(sin(b)*sin(d)/32,onepeapcp));// OK
                    r207_r208   = sum(mul((1-cos(b))*(1-cos(d))/32,onep), mul(-(1-cos(b))*sin(d)/32,onenecp));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eapcn), mul(-sin(b)*(1+cos(d))/16,eapi));// OK
                    r2011         = mul(-(1-cos(b))*sin(d)/16, ecni); // OK
                    complex r2m0y1g1 = {(1-cos(b))*(1+cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r2m0y1g1 = sum(sum(r2m0y1g1, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul((1-cos(b))*sin(d)/32, onepecn), mul(-(1-cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneancn), mul(-sin(b)*(1-cos(d))/32, onepean)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1+cos(d))/32, oneneap), mul(-sin(b)*sin(d)/32,onepeapcp));// OK
                    r207_r208   = sum(mul(-(1+cos(b))*(1+cos(d))/32,onep), mul((1+cos(b))*sin(d)/32,onenecp));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eapcn), mul(-sin(b)*(1-cos(d))/16,eapi));// OK
                    r2011         = mul(-(1+cos(b))*sin(d)/16, ecni); // OK
                    complex r2m1y0g0 = {(1+cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r2m1y0g0 = sum(sum(r2m1y0g0, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul((1+cos(b))*sin(d)/32, onepecn), mul(-(1+cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneancn), mul(sin(b)*(1-cos(d))/32, onepean)); // OK
                    r205_r206   = sum(mul(sin(b)*(1+cos(d))/32, oneneap), mul(sin(b)*sin(d)/32,onepeapcp));// OK
                    r207_r208   = sum(mul(-(1-cos(b))*(1+cos(d))/32,onep), mul((1-cos(b))*sin(d)/32,onenecp));// OK
                    r209_r2010  = sum(mul(sin(b)*sin(d)/16, eapcn), mul(sin(b)*(1-cos(d))/16,eapi));// OK
                    r2011         = mul(-(1-cos(b))*sin(d)/16, ecni); // OK
                    complex r2m1y0g1 = {(1-cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r2m1y0g1 = sum(sum(r2m1y0g1, tmp), sum(r209_r2010, r2011));//OK

                    // 5 bd
                    r201_r202   = sum(mul((1-cos(b))*sin(d)/32, onepecn), mul((1-cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneancn), mul(-sin(b)*(1-cos(d))/32, onepean)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1+cos(d))/32, oneneap), mul(sin(b)*sin(d)/32,onepeapcp));// OK
                    r207_r208   = sum(mul((1+cos(b))*(1+cos(d))/32,onep), mul((1+cos(b))*sin(d)/32,onenecp));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eapcn), mul(sin(b)*(1-cos(d))/16,eapi));// OK
                    r2011         = mul((1+cos(b))*sin(d)/16, ecni); // OK
                    complex r2m1y1g0 = {(1+cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r2m1y1g0 = sum(sum(r2m1y1g0, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul((1+cos(b))*sin(d)/32, onepecn), mul((1+cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneancn), mul(sin(b)*(1-cos(d))/32, onepean)); // OK
                    r205_r206   = sum(mul(sin(b)*(1+cos(d))/32, oneneap), mul(-sin(b)*sin(d)/32,onepeapcp));// OK
                    r207_r208   = sum(mul((1-cos(b))*(1+cos(d))/32,onep), mul((1-cos(b))*sin(d)/32,onenecp));// OK
                    r209_r2010  = sum(mul(sin(b)*sin(d)/16, eapcn), mul(-sin(b)*(1-cos(d))/16,eapi));// OK
                    r2011         = mul((1-cos(b))*sin(d)/16, ecni); // OK
                    complex r2m1y1g1 = {(1-cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r2m1y1g1 = sum(sum(r2m1y1g1, tmp), sum(r209_r2010, r2011));// OK

                    /////////////// r3 ///////////////
                    // 5 bd
                    r201_r202   = sum(mul((1+cos(b))*sin(d)/32, onepecp), mul(-(1+cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneapcp), mul(sin(b)*(1-cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(sin(b)*(1+cos(d))/32, onenean), mul(sin(b)*sin(d)/32,onepeancn));// OK 
                    r207_r208   = sum(mul(-(1-cos(b))*(1+cos(d))/32,onep), mul((1-cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(sin(b)*sin(d)/16, eancp), mul(sin(b)*(1-cos(d))/16,eani));// OK
                    r2011         = mul(-(1-cos(b))*sin(d)/16, ecpi); // OK
                    complex r3m0y0g0 = {(1-cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r3m0y0g0 = sum(sum(r3m0y0g0, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul((1-cos(b))*sin(d)/32, onepecp), mul(-(1-cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneapcp), mul(-sin(b)*(1-cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1+cos(d))/32, onenean), mul(-sin(b)*sin(d)/32,onepeancn));// OK
                    r207_r208   = sum(mul(-(1+cos(b))*(1+cos(d))/32,onep), mul((1+cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eancp), mul(-sin(b)*(1-cos(d))/16,eani));// OK
                    r2011         = mul(-(1+cos(b))*sin(d)/16, ecpi); // OK
                    complex r3m0y0g1 = {(1+cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r3m0y0g1 = sum(sum(r3m0y0g1, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul((1+cos(b))*sin(d)/32, onepecp), mul((1+cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneapcp), mul(sin(b)*(1-cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(sin(b)*(1+cos(d))/32, onenean), mul(-sin(b)*sin(d)/32,onepeancn));// OK
                    r207_r208   = sum(mul((1-cos(b))*(1+cos(d))/32,onep), mul((1-cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(sin(b)*sin(d)/16, eancp), mul(-sin(b)*(1-cos(d))/16,eani));// OK
                    r2011         = mul((1-cos(b))*sin(d)/16, ecpi); // OK
                    complex r3m0y1g0 = {(1-cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r3m0y1g0 = sum(sum(r3m0y1g0, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul((1-cos(b))*sin(d)/32, onepecp), mul((1-cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneapcp), mul(-sin(b)*(1-cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1+cos(d))/32, onenean), mul(sin(b)*sin(d)/32,onepeancn));// OK
                    r207_r208   = sum(mul((1+cos(b))*(1+cos(d))/32,onep), mul((1+cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eancp), mul(sin(b)*(1-cos(d))/16,eani));// OK
                    r2011         = mul((1+cos(b))*sin(d)/16, ecpi); // OK
                    complex r3m0y1g1 = {(1+cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r3m0y1g1 = sum(sum(r3m0y1g1, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul(-(1+cos(b))*sin(d)/32, onepecp), mul(-(1+cos(b))*(1+cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneapcp), mul(sin(b)*(1+cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(sin(b)*(1-cos(d))/32, onenean), mul(-sin(b)*sin(d)/32,onepeancn));// OK
                    r207_r208   = sum(mul(-(1-cos(b))*(1-cos(d))/32,onep), mul(-(1-cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eancp), mul(sin(b)*(1+cos(d))/16,eani));// OK
                    r2011         = mul((1-cos(b))*sin(d)/16, ecpi); // OK
                    complex r3m1y0g0 = {(1-cos(b))*(1+cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r3m1y0g0 = sum(sum(r3m1y0g0, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul(-(1-cos(b))*sin(d)/32, onepecp), mul(-(1-cos(b))*(1+cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneapcp), mul(-sin(b)*(1+cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1-cos(d))/32, onenean), mul(sin(b)*sin(d)/32,onepeancn));// OK
                    r207_r208   = sum(mul(-(1+cos(b))*(1-cos(d))/32,onep), mul(-(1+cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(sin(b)*sin(d)/16, eancp), mul(-sin(b)*(1+cos(d))/16,eani));// OK
                    r2011         = mul((1+cos(b))*sin(d)/16, ecpi); // OK
                    complex r3m1y0g1 = {(1+cos(b))*(1+cos(d))/16, 0}; // 
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r3m1y0g1 = sum(sum(r3m1y0g1, tmp), sum(r209_r2010, r2011));//OK

                    // 5 bd
                    r201_r202   = sum(mul(-(1+cos(b))*sin(d)/32, onepecp), mul((1+cos(b))*(1+cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneapcp), mul(sin(b)*(1+cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(sin(b)*(1-cos(d))/32, onenean), mul(sin(b)*sin(d)/32,onepeancn));// OK
                    r207_r208   = sum(mul((1-cos(b))*(1-cos(d))/32,onep), mul(-(1-cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eancp), mul(-sin(b)*(1+cos(d))/16,eani));// OK
                    r2011         = mul(-(1-cos(b))*sin(d)/16, ecpi); // OK
                    complex r3m1y1g0 = {(1-cos(b))*(1+cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r3m1y1g0 = sum(sum(r3m1y1g0, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul(-(1-cos(b))*sin(d)/32, onepecp), mul((1-cos(b))*(1+cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneapcp), mul(-sin(b)*(1+cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1-cos(d))/32, onenean), mul(-sin(b)*sin(d)/32,onepeancn));// OK
                    r207_r208   = sum(mul((1+cos(b))*(1-cos(d))/32,onep), mul(-(1+cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(sin(b)*sin(d)/16, eancp), mul(sin(b)*(1+cos(d))/16,eani));// OK
                    r2011         = mul(-(1+cos(b))*sin(d)/16, ecpi); // OK
                    complex r3m1y1g1 = {(1+cos(b))*(1+cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r3m1y1g1 = sum(sum(r3m1y1g1, tmp), sum(r209_r2010, r2011));// OK

                    /////////////// r4 ///////////////
                    // 5 bd
                    r201_r202   = sum(mul(-(1+cos(b))*sin(d)/32, onepecp), mul((1+cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneapcp), mul(-sin(b)*(1-cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1+cos(d))/32, onenean), mul(-sin(b)*sin(d)/32,onepeancn));// OK
                    r207_r208   = sum(mul((1-cos(b))*(1+cos(d))/32,onep), mul(-(1-cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(sin(b)*sin(d)/16, eancp), mul(sin(b)*(1-cos(d))/16,eani));// OK
                    r2011         = mul(-(1-cos(b))*sin(d)/16, ecpi); // OK
                    complex r4m0y0g0 = {(1-cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r4m0y0g0 = sum(sum(r4m0y0g0, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul(-(1-cos(b))*sin(d)/32, onepecp), mul((1-cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneapcp), mul(sin(b)*(1-cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(sin(b)*(1+cos(d))/32, onenean), mul(sin(b)*sin(d)/32,onepeancn));// OK
                    r207_r208   = sum(mul((1+cos(b))*(1+cos(d))/32,onep), mul(-(1+cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eancp), mul(-sin(b)*(1-cos(d))/16,eani));// OK
                    r2011         = mul(-(1+cos(b))*sin(d)/16, ecpi); // OK
                    complex r4m0y0g1 = {(1+cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r4m0y0g1 = sum(sum(r4m0y0g1, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul(-(1+cos(b))*sin(d)/32, onepecp), mul(-(1+cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneapcp), mul(-sin(b)*(1-cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1+cos(d))/32, onenean), mul(sin(b)*sin(d)/32,onepeancn));// OK
                    r207_r208   = sum(mul(-(1-cos(b))*(1+cos(d))/32,onep), mul(-(1-cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(sin(b)*sin(d)/16, eancp), mul(-sin(b)*(1-cos(d))/16,eani));// OK
                    r2011         = mul((1-cos(b))*sin(d)/16, ecpi); // OK
                    complex r4m0y1g0 = {(1-cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r4m0y1g0 = sum(sum(r4m0y1g0, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul(-(1-cos(b))*sin(d)/32, onepecp), mul(-(1-cos(b))*(1-cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneapcp), mul(sin(b)*(1-cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(sin(b)*(1+cos(d))/32, onenean), mul(-sin(b)*sin(d)/32,onepeancn));// OK
                    r207_r208   = sum(mul(-(1+cos(b))*(1+cos(d))/32,onep), mul(-(1+cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eancp), mul(sin(b)*(1-cos(d))/16,eani));// OK
                    r2011         = mul((1+cos(b))*sin(d)/16, ecpi); // OK
                    complex r4m0y1g1 = {(1+cos(b))*(1-cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r4m0y1g1 = sum(sum(r4m0y1g1, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul((1+cos(b))*sin(d)/32, onepecp), mul((1+cos(b))*(1+cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneapcp), mul(-sin(b)*(1+cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1-cos(d))/32, onenean), mul(sin(b)*sin(d)/32,onepeancn));// OK
                    r207_r208   = sum(mul((1-cos(b))*(1-cos(d))/32,onep), mul((1-cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eancp), mul(sin(b)*(1+cos(d))/16,eani));// OK
                    r2011         = mul((1-cos(b))*sin(d)/16, ecpi); // OK
                    complex r4m1y0g0 = {(1-cos(b))*(1+cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r4m1y0g0 = sum(sum(r4m1y0g0, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul((1-cos(b))*sin(d)/32, onepecp), mul((1-cos(b))*(1+cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneapcp), mul(sin(b)*(1+cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(sin(b)*(1-cos(d))/32, onenean), mul(-sin(b)*sin(d)/32,onepeancn));// OK
                    r207_r208   = sum(mul((1+cos(b))*(1-cos(d))/32,onep), mul((1+cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(sin(b)*sin(d)/16, eancp), mul(-sin(b)*(1+cos(d))/16,eani));// OK
                    r2011         = mul((1+cos(b))*sin(d)/16, ecpi); // OK
                    complex r4m1y0g1 = {(1+cos(b))*(1+cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r4m1y0g1 = sum(sum(r4m1y0g1, tmp), sum(r209_r2010, r2011));//OK

                    // 5 bd
                    r201_r202   = sum(mul((1+cos(b))*sin(d)/32, onepecp), mul(-(1+cos(b))*(1+cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(-sin(b)*sin(d)/32, oneneapcp), mul(-sin(b)*(1+cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(-sin(b)*(1-cos(d))/32, onenean), mul(-sin(b)*sin(d)/32,onepeancn));// OK
                    r207_r208   = sum(mul(-(1-cos(b))*(1-cos(d))/32,onep), mul((1-cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(-sin(b)*sin(d)/16, eancp), mul(-sin(b)*(1+cos(d))/16,eani));// OK
                    r2011         = mul(-(1-cos(b))*sin(d)/16, ecpi); // OK
                    complex r4m1y1g0 = {(1-cos(b))*(1+cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r4m1y1g0 = sum(sum(r4m1y1g0, tmp), sum(r209_r2010, r2011));// OK

                    // 5 bd
                    r201_r202   = sum(mul((1-cos(b))*sin(d)/32, onepecp), mul(-(1-cos(b))*(1+cos(d))/32, onen)); // OK
                    r203_r204   = sum(mul(sin(b)*sin(d)/32, oneneapcp), mul(sin(b)*(1+cos(d))/32, onepeap)); // OK
                    r205_r206   = sum(mul(sin(b)*(1-cos(d))/32, onenean), mul(sin(b)*sin(d)/32,onepeancn));// OK
                    r207_r208   = sum(mul(-(1+cos(b))*(1-cos(d))/32,onep), mul((1+cos(b))*sin(d)/32,onenecn));// OK
                    r209_r2010  = sum(mul(sin(b)*sin(d)/16, eancp), mul(sin(b)*(1+cos(d))/16,eani));// OK
                    r2011         = mul(-(1+cos(b))*sin(d)/16, ecpi); // OK
                    complex r4m1y1g1 = {(1+cos(b))*(1+cos(d))/16, 0}; // OK
                    tmp = sum(sum(r201_r202, r203_r204), sum(r205_r206, r207_r208));// OK
                    r4m1y1g1 = sum(sum(r4m1y1g1, tmp), sum(r209_r2010, r2011));// OK

                    
                    // 関数定義部分
                    double l3,l4,l7,l8,u3,u4,u7,u8;
                    l3 = Pccn(r2m0y0g0) + Pccn(r2m0y0g1) + Pccn(r2m1y0g0) + Pccn(r2m1y0g1);
                    l4 = Pccn(r2m0y1g0) + Pccn(r2m0y1g1) + Pccn(r2m1y1g0) + Pccn(r2m1y1g1);
                    l7 = Pccn(r4m0y0g0) + Pccn(r4m0y0g1) + Pccn(r4m1y0g0) + Pccn(r4m1y0g1);
                    l8 = Pccn(r4m0y1g0) + Pccn(r4m0y1g1) + Pccn(r4m1y1g0) + Pccn(r4m1y1g1);
                    u3 = Pccn(r2m0z0g0) + Pccn(r2m0z0g1) + Pccn(r2m1z0g0) + Pccn(r2m1z0g1);
                    u4 = Pccn(r2m0z1g0) + Pccn(r2m0z1g1) + Pccn(r2m1z1g0) + Pccn(r2m1z1g1);
                    u7 = Pccn(r4m0z0g0) + Pccn(r4m0z0g1) + Pccn(r4m1z0g0) + Pccn(r4m1z0g1);
                    u8 = Pccn(r4m0z1g0) + Pccn(r4m0z1g1) + Pccn(r4m1z1g0) + Pccn(r4m1z1g1);

                    // 関数定義部分
                    double l1,l2,l5,l6,u1,u2,u5,u6;
                    l1 = Pccn(r1m0y0g0) + Pccn(r1m0y0g1) + Pccn(r1m1y0g0) + Pccn(r1m1y0g1);
                    l2 = Pccn(r1m0y1g0) + Pccn(r1m0y1g1) + Pccn(r1m1y1g0) + Pccn(r1m1y1g1);
                    l5 = Pccn(r3m0y0g0) + Pccn(r3m0y0g1) + Pccn(r3m1y0g0) + Pccn(r3m1y0g1);
                    l6 = Pccn(r3m0y1g0) + Pccn(r3m0y1g1) + Pccn(r3m1y1g0) + Pccn(r3m1y1g1);
                    u1 = Pccn(r1m0z0g0) + Pccn(r1m0z0g1) + Pccn(r1m1z0g0) + Pccn(r1m1z0g1);
                    u2 = Pccn(r1m0z1g0) + Pccn(r1m0z1g1) + Pccn(r1m1z1g0) + Pccn(r1m1z1g1);
                    u5 = Pccn(r3m0z0g0) + Pccn(r3m0z0g1) + Pccn(r3m1z0g0) + Pccn(r3m1z0g1);
                    u6 = Pccn(r3m0z1g0) + Pccn(r3m0z1g1) + Pccn(r3m1z1g0) + Pccn(r3m1z1g1);

                    double R = l3+l4+l7+l8+u3+u4+u7+u8;
                    double R2 = (l3+l4+u3+u4)/R;
                    double R4 = (l7+l8+u7+u8)/R;

                    double PAB = R2*(l4/(l3+l4) + u3/(u3+u4))/2 + R4*(l7/(l7+l8) + u8/(u7+u8))/2;

                    double m1,m2,m5,m6,v1,v2,v5,v6;
                    m1 = Pccn(r1m0y0g0) + Pccn(r1m0y1g0) + Pccn(r1m1y0g1) + Pccn(r1m1y1g1);
                    m2 = Pccn(r1m0y0g1) + Pccn(r1m0y1g1) + Pccn(r1m1y0g0) + Pccn(r1m1y1g0);
                    m5 = Pccn(r3m0y0g0) + Pccn(r3m0y1g0) + Pccn(r3m1y0g1) + Pccn(r3m1y1g1);
                    m6 = Pccn(r3m0y0g1) + Pccn(r3m0y1g1) + Pccn(r3m1y0g0) + Pccn(r3m1y1g0);
                    v1 = Pccn(r1m0z0g0) + Pccn(r1m0z1g0) + Pccn(r1m1z0g1) + Pccn(r1m1z1g1);
                    v2 = Pccn(r1m0z0g1) + Pccn(r1m0z1g1) + Pccn(r1m1z0g0) + Pccn(r1m1z1g0);
                    v5 = Pccn(r3m0z0g0) + Pccn(r3m0z1g0) + Pccn(r3m1z0g1) + Pccn(r3m1z1g1);
                    v6 = Pccn(r3m0z0g1) + Pccn(r3m0z1g1) + Pccn(r3m1z0g0) + Pccn(r3m1z1g0);

                    double Rd = m1+m2+m5+m6+v1+v2+v5+v6;
                    double R1 = (m1+m2+v1+v2)/Rd;
                    double R3 = (m5+m6+v5+v6)/Rd;
                    double PE = R1*(v1/(v1+v2) +m1/(m1+m2))/2 + R3*(v6/(v5+v6) + m6/(m5+m6))/2;

                    R = l1+l2+l5+l6+u1+u2+u5+u6;
                    double A0B0 = (l1 + u1)/R;
                    double A0B1 = (l2 + u2)/R;
                    double A1B0 = (l5 + u5)/R;
                    double A1B1 = (l6 + u6)/R;
                    double IAB = A0B0*log2(A0B0/((A0B0+A0B1)*(A0B0+A1B0)))+A0B1*log2(A0B1/((A0B0+A0B1)*(A0B1+A1B1)))+A1B0*log2(A1B0/((A1B0+A1B1)*(A0B0+A1B0)))+A1B1*log2(A1B1/((A1B0+A1B1)*(A0B1+A1B1)));
                    
                    // AliceとEveの相互情報量
                    R = v1+v2+v5+v6+m1+m2+m5+m6;
                    double A0E0 = (v1 + m1)/R;
                    double A0E1 = (v2 + m2)/R;
                    double A1E0 = (v5 + m5)/R;
                    double A1E1 = (v6 + m6)/R;
                    double IAE = A0E0*log2(A0E0/((A0E0+A0E1)*(A0E0+A1E0)))+A0E1*log2(A0E1/((A0E0+A0E1)*(A0E1+A1E1)))+A1E0*log2(A1E0/((A1E0+A1E1)*(A0E0+A1E0)))+A1E1*log2(A1E1/((A1E0+A1E1)*(A0E1+A1E1)));

                    // 探索判定部
                    if (PE > 0.853 || PAB > 0.899 || IAE > 0.39898){
                        fprintf(file,"%d,%d,%d,%d,%lf,%lf,%lf,%lf\n", i, j, k, l, PAB, PE, IAE, IAB);
                    }
                }
            }
        }
    }
    fclose(file);
    return 0;
}
// 2025 / 2 
