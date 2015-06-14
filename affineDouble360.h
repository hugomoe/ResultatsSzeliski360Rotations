#ifndef _AFFINEDOUBLE360_H
#define _AFFINEDOUBLE360_H



//parameter for filters
    #define TAPS 4
    #define PREC 10
//parameters for raised cosine-weighted sinc filter
    #define PERIOD 1.
    //#define BETA 0.36
//parameter for B-Spline
    //#define ORDER 5
//constants
#define PREC2 20
#define PI 3.14159265358979323

bool eq(double a,double b);
void transpo_opt(double *img,double *a,int wh[2]);

double filter_fun(double x);

double filter_h(double *img,int w,int h,double *H,int N,int prec,double a0,double a1,int i,int j,int l);
double filter_v(double *img,int w,int h,double *H,int N,int prec,double a0,double a1,int i,int j,int l);

int apply_rh(double *img1,double *img2,int w,int h,double s,double a0,double a1);
int apply_rv(double *img1,double *img2,int w,int h,double s,double a0,double a1);

int apply_affinite(double *img,double *img_f,int w,int h,int w_f,int h_f,double *affinity);

#endif // _AFFINEDOUBLE360_H
