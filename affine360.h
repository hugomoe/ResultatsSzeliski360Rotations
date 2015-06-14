#ifndef _AFFINE360_H
#define _AFFINE360_H



//parameter for filters
    #define TAPS 4
    #define PREC 10
//parameters for raised cosine-weighted sinc filter
    #define PERIOD 1.
    //#define BETA 0.36
//parameter for B-Spline
    //#define ORDER 3
//constants
#define PREC2 20
#define PI 3.14159265358979323

bool eq(double a,double b);
void transpo_opt(float *img,double *a,int wh[2]);

float filter_fun(float x);

float filter_h(float *img,int w,int h,float *H,int N,int prec,double a0,double a1,int i,int j,int l);
float filter_v(float *img,int w,int h,float *H,int N,int prec,double a0,double a1,int i,int j,int l);

int apply_rh(float *img1,float *img2,int w,int h,double s,double a0,double a1);
int apply_rv(float *img1,float *img2,int w,int h,double s,double a0,double a1);

int apply_affinite(float *img,float *img_f,int w,int h,int w_f,int h_f,double *affinity);

#endif // _AFFINE360_H
